package cnv.manage;

import gwas.PhenoPrep;
import gwas.Qc;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import seq.manage.VCFOps;
import seq.manage.VCOps;
import stats.Rscript.RScatter;
import stats.Rscript.SCATTER_TYPE;
import stats.Rscript.SlopeLine;
import stats.StatsCrossTabs;
import stats.StatsCrossTabs.STAT_TYPE;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.Sort;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerSet;
import cnv.filesys.Pedigree;
import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;

/**
 * Run a gwas of original (Array) genotypes for mitochondrial cn for estimates from particular pc-estimates/pedigree/covars <br>
 * Currently hard-codes the recommended analysis (PC15,PC150 for natural order PCs and PC15, PCLast for Stepwise)<br>
 * Requires Plink2
 */
public class MitoGWAS {
	private static final String SUB_DIR = "_eval/typed/";

	/**
	 * @param proj
	 * @param pedFile
	 *            ped file with dna, full path
	 * @param pcFile
	 *            pc file ( must have valid eval data structure), full path
	 * @param covarFile
	 *            covariates to use, full path
	 * @param outputDir
	 *            where the analysis will happen, full path
	 */
	public static void analyze(Project proj, String pedFile, String pcFile, String covarFile, String outputDir, int numthreads) {
		Pedigree ped = new Pedigree(proj, pedFile, false);

		String fullOut = outputDir + ext.rootOf(proj.getPropertyFilename()) + "/";
		new File(fullOut).mkdirs();

		String root = fullOut + ext.rootOf(proj.getPropertyFilename());
		String newPed = root + ".pedigree.dat";

		ped.writeToFile(newPed, true);
		proj.PEDIGREE_FILENAME.setValue(newPed);

		String plinkPed = root + ".ped";
		String plinkMap = root + ".map";
		if (!Files.exists(plinkMap) || !Files.exists(plinkPed)) {
			MarkerSet markerSet = proj.getMarkerSet();
			String[] markerNames = markerSet.getMarkerNames();
			byte[] chr = markerSet.getChrs();
			ArrayList<String> markersToAnalyze = new ArrayList<String>();
			int numSkipped = 0;
			for (int i = 0; i < markerNames.length; i++) {
				if (proj.getArrayType().isCNOnly(markerNames[i]) || chr[i] < 1) {
					numSkipped++;
				} else {
					markersToAnalyze.add(markerNames[i]);
				}
			}
			proj.getLog().reportTimeInfo(numSkipped + " copy number or chr0 only probes  were removed, " + markersToAnalyze.size() + " remaining");
			String exportList = root + "_markers.txt";
			Files.writeList(Array.toStringArray(markersToAnalyze), exportList);
			String blankCluster = root + ".blankCluster.ser";
			proj.getLog().reportTimeWarning("Using blank cluster filter file to ensure AB lookup");
			new ClusterFilterCollection().serialize(blankCluster);
			proj.GC_THRESHOLD.setValue((double) 0);
			proj.getLog().reportTimeWarning("Setting gc threshold to 0");
			PlinkData.saveGenvisisToPlinkPedSet(proj, root, "", blankCluster, exportList, numthreads);
		} else {
			proj.getLog().reportTimeInfo(plinkMap + " and " + plinkPed + "exist, skipping");
		}

		double maf = (double) 5 / ped.getDnas().length;
		double geno = 0.1;
		double mind = 0.1;
		proj.getLog().reportTimeInfo("using maf of " + maf + " (MAF= " + maf + "  of " + ped.getDnas().length);

		ArrayList<String> plinkConverCommand = new ArrayList<String>();
		plinkConverCommand.add("plink2");
		plinkConverCommand.add("--file");
		plinkConverCommand.add(root);
		plinkConverCommand.add("--make-bed");
		plinkConverCommand.add("--out");
		root = root + "_maf_" + ext.roundToSignificantFigures(maf, 5) + "_geno_" + geno + "_mind_" + mind;
		plinkConverCommand.add(root);
		plinkConverCommand.add("--maf");
		plinkConverCommand.add(maf + "");
		plinkConverCommand.add("--geno");
		plinkConverCommand.add(geno + "");
		plinkConverCommand.add("--mind");
		plinkConverCommand.add("" + mind);

		String[] in = new String[] { plinkPed, plinkMap };
		String fam = root + ".fam";
		String bim = root + ".bim";
		String bed = root + ".bed";

		Files.copyFileUsingFileChannels(fam, fullOut + "plink.fam", proj.getLog());
		Files.copyFileUsingFileChannels(bim, fullOut + "plink.bim", proj.getLog());
		Files.copyFileUsingFileChannels(bed, fullOut + "plink.bed", proj.getLog());
		// Qc.fullGamut(fullOut, true, proj.getLog());

		String[] out = new String[] { fam, bim, bed };
		CmdLine.runCommandWithFileChecks(Array.toStringArray(plinkConverCommand), "", in, out, true, false, false, proj.getLog());
		ArrayList<PlinkAssoc> plinkCommands = new ArrayList<MitoGWAS.PlinkAssoc>();
		String subDir = ext.rootOf(pcFile, false) + SUB_DIR;
		String natty = subDir + "WITHOUT_BUILDERS_NATURAL_WITHOUT_INDEPS_finalSummary.estimates.txt.gz";
		String[] nattyTitles = new String[] { "PC15", "PC150" };
		plinkCommands.addAll(generateCommands(proj, root, natty, nattyTitles, covarFile, ped));

		String stepWise = subDir + "WITHOUT_BUILDERS_STEPWISE_RANK_R2_WITHOUT_INDEPS_finalSummary.estimates.txt.gz";
		String[] header = Files.getHeaderOfFile(stepWise, proj.getLog());
		String[] stepWiseTitles = new String[] { "PC15", header[header.length - 1] };

		plinkCommands.addAll(generateCommands(proj, root, stepWise, stepWiseTitles, covarFile, ped));

		PlinkAssocProducer producer = new PlinkAssocProducer(plinkCommands, proj.getLog());
		WorkerTrain<Boolean> train = new WorkerTrain<Boolean>(producer, numthreads, 2, proj.getLog());
		while (train.hasNext()) {
			train.next();
		}

		summarize(proj, root, plinkCommands);
	}

	private static void summarize(Project proj, String root, ArrayList<PlinkAssoc> plinkCommands) {
		System.exit(1);
		proj.getLog().reportTimeInfo("Computing correlation matrix of results");
		double[][] emp1s = new double[plinkCommands.size()][];
		int[][] keys = new int[plinkCommands.size()][];
		String[] empTitles = new String[plinkCommands.size()];
		String outQQ = root + "summaryPvals";
		String pvalDB = outQQ + ".txt";
		String pvalQQ = outQQ + ".qq.txt";
		// /String plot = outQQ + ".jpeg";
		String outTabs = root + "stabs.correl.txt";
		for (int j = 0; j < plinkCommands.size(); j++) {
			String results = plinkCommands.get(j).getOutputs()[0];
			empTitles[j] = ext.removeDirectoryInfo(results);
		}
		empTitles = Array.untag(empTitles, true, true);
		String[] pvalFiles = new String[plinkCommands.size()];
		for (int i = 0; i < pvalFiles.length; i++) {
			pvalFiles[i] = ext.parseDirectoryOfFile(pvalDB) + empTitles[i] + ".pvalsQQ.txt";
		}
		String[] empLogP = Array.tagOn(empTitles, "p_", null);
		if (!Files.exists(pvalDB) || !Files.exists(pvalQQ)) {
			for (int j = 0; j < plinkCommands.size(); j++) {
				String results = plinkCommands.get(j).getOutputs()[1];
				proj.getLog().reportTimeInfo("Loading " + results);
				ProjectDataParserBuilder builderPermResults = new ProjectDataParserBuilder();
				builderPermResults.treatAllNumeric(false);
				builderPermResults.sampleBased(false);
				builderPermResults.hasHeader(true);
				builderPermResults.dataKeyColumnName("SNP");
				builderPermResults.requireAll(false);
				builderPermResults.separator("[\\s]+");
				builderPermResults.setInvalidNumericToNaN(true);
				builderPermResults.firstEntryOnly(true);
				builderPermResults.numericDataTitles(new String[] { "P" });
				ExtProjectDataParser permParser;
				try {
					permParser = builderPermResults.build(proj, results);
					permParser.determineIndicesFromTitles();
					permParser.loadData();
					emp1s[j] = permParser.getNumericDataForTitle("P");
					Files.writeList(Array.toStringArray(Array.removeNaN(emp1s[j])), pvalFiles[j]);
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					return;
				}

			}
			StatsCrossTabs sTabs = new StatsCrossTabs(emp1s, null, null, empTitles, STAT_TYPE.SPEARMAN_CORREL, true, proj.getLog());
			sTabs.computeTable();
			sTabs.dumpTables(outTabs);
			MarkerSet markerSet = proj.getMarkerSet();
			boolean[] valids = Array.booleanArray(emp1s[0].length, true);
			VCFFileReader annoReader = new VCFFileReader(proj.BLAST_ANNOTATION_FILENAME.getValue(), true);
			String[] annotations = VCFOps.getAnnotationKeys(proj.BLAST_ANNOTATION_FILENAME.getValue(), proj.getLog())[0];

			CloseableIterator<VariantContext> vcIter = annoReader.iterator();
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(pvalDB));
				writer.println("SNP\tCHR\tBP\tREF\tALT" + Array.toStr(empTitles) + "\t" + Array.toStr(annotations));
				int numInvalid = 0;
				for (int i = 0; i < emp1s[0].length; i++) {
					if (!vcIter.hasNext()) {
						writer.close();
						annoReader.close();
						throw new IllegalArgumentException("Invalid marker annotation " + proj.BLAST_ANNOTATION_FILENAME.getValue());
					}
					VariantContext markVC = vcIter.next();
					if (!markVC.getID().equals(markerSet.getMarkerNames()[i])) {
						writer.close();
						annoReader.close();
						throw new IllegalArgumentException("Mismatched file/project order for " + proj.BLAST_ANNOTATION_FILENAME.getValue());
					}

					boolean valid = true;
					StringBuilder builder = new StringBuilder();
					for (int j = 0; j < emp1s.length; j++) {
						if (valid) {
							valid = Double.isFinite(emp1s[j][i]) && markerSet.getChrs()[i] > 0;
							if (!valid) {
								numInvalid++;
								valids[i] = false;
							}
							if (j == 0) {
								builder.append(markerSet.getMarkerNames()[i] + "\t" + markerSet.getChrs()[i] + "\t" + markerSet.getPositions()[i] + "\t" + markVC.getReference().getDisplayString() + "\t" + markVC.getAlternateAlleles().toString() + emp1s[j][i]);
							} else {
								builder.append("\t" + emp1s[j][i]);
							}
						}
					}
					builder.append(VCOps.getAnnotationsFor(annotations, markVC, "."));
					if (valid) {
						writer.println(builder.toString());
					}

				}
				if (numInvalid > 0) {
					proj.getLog().reportTimeWarning(numInvalid + " markers had an NaN or chr 0 in one of the gwas runs ( or was not included in the analysis), removed");
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + pvalDB);
				proj.getLog().reportException(e);
			}
			annoReader.close();

			try {

				for (int i = 0; i < emp1s.length; i++) {
					emp1s[i] = Array.subArray(emp1s[i], valids);
					keys[i] = Sort.trickSort(emp1s[i]);
				}
				PrintWriter writer = new PrintWriter(new FileWriter(pvalQQ));
				writer.println("RANK\t" + Array.toStr(empLogP));
				for (int i = 0; i < emp1s[0].length; i++) {
					StringBuilder builder = new StringBuilder();
					double ranke = (double) (i + 1);
					double rankP = -1 * Math.log10((double) ranke / emp1s[0].length);
					builder.append(rankP);
					for (int j = 0; j < emp1s.length; j++) {
						double plog = -1 * Math.log10(emp1s[j][keys[j][i]]);
						builder.append("\t" + plog);
					}
					writer.println(builder.toString());
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + pvalQQ);
				proj.getLog().reportException(e);
			}
		}
		proj.QQ_FILENAMES.setValue(pvalFiles);
		proj.saveProperties();
		String out = ext.addToRoot(pvalQQ, ".qqplot");
		RScatter rScatter = new RScatter(pvalQQ, out + ".rscript", ext.removeDirectoryInfo(out), out + ".jpeg", "RANK", empLogP, null, SCATTER_TYPE.POINT, proj.getLog());
		SlopeLine slopeLine = new SlopeLine(0, 1);
		rScatter.setSlopeLines(new SlopeLine[] { slopeLine });
		rScatter.setOverWriteExisting(true);
		rScatter.execute();

		// ArrayList<RScatter> rs = new ArrayList<RScatter>();
		// RScatters rsScatters = new RScatters(rs.toArray(new RScatter[rs.size()]), pvalQQ + "finalRscript", pvalQQ + "final.pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, proj.getLog());

		// String[] script = generateManhatQQScript(pvalDB, empTitles, plot);
		// String rscript = plot + ".rscript";
		// Files.writeList(script, rscript);
		// CmdLine.runCommandWithFileChecks(new String[] { "/panfs/roc/itascasoft/R/3.1.1/bin/Rscript", rscript }, "", null, new String[] { plot }, true, true, false, proj.getLog());
	}

	// private static String[] generateMadnhatQQScript(String db, String[] pvalColumns, String output) {
	// ArrayList<String> command = new ArrayList<String>();
	// ArrayList<String> order = new ArrayList<String>();
	//
	// command.add("library(qqman)");
	//
	// String main = "data";
	// command.add(main + "=read.table(\"" + db + "\", header=TRUE)");
	// //
	// for (int i = 0; i < pvalColumns.length; i++) {
	// command.add("jpeg(file=\"" + ext.parseDirectoryOfFile(output) + pvalColumns[i] + ".jpeg" + "\",height=2000,width=2500)");
	// // ,onefile = TRUE
	// command.add("op <- par(mfrow=c(1,2))");
	//
	// order.add(pvalColumns[i]);
	// String man = pvalColumns[i] + "man =manhattan(" + main + ", p =\"" + pvalColumns[i] + "\")";
	// command.add(man);
	// command.add(pvalColumns[i] + "q =qq(" + main + "$" + pvalColumns[i] + ")");
	// command.add(pvalColumns[i] + "man");
	// command.add("title(main = \"" + pvalColumns[i] + "\")");
	// command.add(pvalColumns[i] + "q");
	// command.add("title(main = \"" + pvalColumns[i] + "\")");
	// command.add("par(op)");
	// command.add("dev.off()");
	// }
	//
	// Files.writeList(Array.toStringArray(order), output + ".order.txt");
	// return Array.toStringArray(command);
	// }

	private static ArrayList<PlinkAssoc> generateCommands(Project proj, String root, String mtPhenoFile, String[] titles, String covarFile, Pedigree ped) {
		ArrayList<PlinkAssoc> plinkCommands = new ArrayList<PlinkAssoc>();

		ProjectDataParserBuilder builderCurrent = new ProjectDataParserBuilder();
		builderCurrent.numericDataTitles(titles);
		builderCurrent.sampleBased(true);
		builderCurrent.dataKeyColumnName("DNA");
		builderCurrent.treatAllNumeric(false);
		builderCurrent.requireAll(true);
		ExtProjectDataParser parser;
		try {
			parser = builderCurrent.build(proj, mtPhenoFile);
			parser.determineIndicesFromTitles();
			parser.loadData();

			ProjectDataParserBuilder builderCovar = new ProjectDataParserBuilder();
			builderCovar.sampleBased(true);
			builderCovar.hasHeader(true);
			builderCovar.dataKeyColumnName("DNA");
			builderCovar.requireAll(false);
			builderCovar.numericDataTitles(Array.subArray(Files.getHeaderOfFile(covarFile, proj.getLog()), 1));
			ExtProjectDataParser covarParser = builderCovar.build(proj, covarFile);
			covarParser.determineIndicesFromTitles();
			covarParser.loadData();
			boolean[] hasVarianceWithinPed = Array.booleanArray(covarParser.getNumericDataTitles().length, false);
			for (int i = 0; i < hasVarianceWithinPed.length; i++) {
				Hashtable<String, String> varHash = new Hashtable<String, String>();
				for (int j = 0; j < ped.getDnas().length; j++) {
					int sampIndex = ext.indexOfStr(ped.getDnas()[j], proj.getSamples());
					String val = covarParser.getNumericDataForTitle(covarParser.getNumericDataTitles()[i])[sampIndex] + "";
					varHash.put(val, val);
					if (varHash.size() > 1) {
						hasVarianceWithinPed[i] = true;
						break;
					}
				}
				if (varHash.size() < 2) {
					proj.getLog().reportTimeWarning(covarParser.getNumericDataTitles()[i] + " had no variance in ped samples, removing");
				}
			}

			for (int i = 0; i < titles.length; i++) {
				String outCurrent = root + ext.rootOf(mtPhenoFile) + "_" + titles[i] + ".txt";
				for (int j = 0; j < covarParser.getNumericDataTitles().length; j++) {
					if (hasVarianceWithinPed[j]) {
						outCurrent = ext.addToRoot(outCurrent, "_" + covarParser.getNumericDataTitles()[j]);
					}
				}
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(outCurrent));
					writer.print("FID\tIID\t" + titles[i]);
					for (int j = 0; j < covarParser.getNumericDataTitles().length; j++) {
						if (hasVarianceWithinPed[j]) {
							writer.print("\t" + covarParser.getNumericDataTitles()[j]);
						}
					}
					writer.println();
					for (int j = 0; j < ped.getDnas().length; j++) {
						int sampIndex = ext.indexOfStr(ped.getDnas()[j], proj.getSamples());
						writer.print(ped.getFID(j) + "\t" + ped.getIID(j) + "\t" + parser.getNumericDataForTitle(titles[i])[sampIndex]);
						for (int k = 0; k < covarParser.getNumericDataTitles().length; k++) {
							if (hasVarianceWithinPed[k]) {
								writer.print("\t" + covarParser.getNumericDataForTitle(covarParser.getNumericDataTitles()[k])[sampIndex]);
							}
						}
						writer.println();
					}
					writer.close();

					String covarTitles = Array.toStr(Array.subArray(covarParser.getNumericDataTitles(), hasVarianceWithinPed), ",");
					String outPrepReg = ext.addToRoot(outCurrent, ".prepped");

					PlinkAssoc regCommand = prepareAssoc(root, false, outCurrent, titles[i], covarTitles, outPrepReg, proj.getLog());
					plinkCommands.add(regCommand);
					String outPrepInv = ext.addToRoot(outCurrent, ".prepped.inv");
					PlinkAssoc invCommand = prepareAssoc(root, true, outCurrent, titles[i], covarTitles, outPrepInv, proj.getLog());
					plinkCommands.add(invCommand);

				} catch (Exception e) {
					proj.getLog().reportError("Error writing to " + outCurrent);
					proj.getLog().reportException(e);
				}
			}

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return plinkCommands;

	}

	private static class PlinkAssocProducer implements Producer<Boolean> {
		private ArrayList<PlinkAssoc> assocs;
		private int index;

		public PlinkAssocProducer(ArrayList<PlinkAssoc> assocs, Logger log) {
			super();
			this.assocs = assocs;
			this.index = 0;

		}

		@Override
		public boolean hasNext() {

			// TODO Auto-generated method stub
			return index < assocs.size();
		}

		@Override
		public Callable<Boolean> next() {
			PlinkAssoc toReturn = assocs.get(index);
			index++;
			return toReturn;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}
	}

	private static class PlinkAssoc implements Callable<Boolean> {
		private ArrayList<String> command;
		private String[] inputs;
		private String[] outputs;
		private Logger log;

		public PlinkAssoc(ArrayList<String> command, String[] inputs, String[] outputs, Logger log) {
			super();
			this.command = command;
			this.inputs = inputs;
			this.outputs = outputs;
			this.log = log;
		}

		public String[] getOutputs() {
			return outputs;
		}

		@Override
		public Boolean call() throws Exception {
			return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, outputs, true, false, false, log);
		}
	}

	private static PlinkAssoc prepareAssoc(String root, boolean inverse, String inputDb, String pheno, String covars, String output, Logger log) {
		String processed = ext.addToRoot(output, "_pheno");
		if (!Files.exists(processed)) {
			PhenoPrep.parse("", inputDb, "IID", pheno, null, 3.0, false, false, true, false, inverse, covars, root + ".fam", true, true, false, false, false, false, null, output, true, false, false, false, false, false, log);
		}
		ArrayList<String> plink = new ArrayList<String>();
		plink.add("plink2");

		plink.add("--linear");
		plink.add("perm");

		plink.add("--bfile");
		plink.add(root);
		plink.add("--covar-name");
		plink.add(covars);
		plink.add("--covar");
		plink.add(inputDb);
		plink.add("--pheno");
		plink.add(processed);
		plink.add("--out");
		plink.add(ext.rootOf(output, false));
		plink.add("--threads");
		plink.add(2 + "");

		String[] inputs = new String[] { processed, inputDb };
		String[] outputs = new String[] { ext.rootOf(output, false) + ".assoc.linear.perm", ext.rootOf(output, false) + ".assoc.linear" };
		PlinkAssoc assoc = new PlinkAssoc(plink, inputs, outputs, log);
		return assoc;
	}

	public static void test() {
		Project proj = new Project("/home/pankrat2/lanej/projects/gedi_gwas.properties", false);
		String ped = "/panfs/roc/groups/5/pankrat2/shared/gedi_gwas/mitoAnalyze/gedi_gwas.ped";
		String covar = "/panfs/roc/groups/5/pankrat2/shared/gedi_gwas/mitoAnalyze/gedi_gwas.covar";
		int numthreads = 12;
		String outputDir = "/scratch.global/lanej/mitoAnalyze/";
		String pcFile = proj.PROJECT_DIRECTORY.getValue() + "gedi_gwasALL_1000PCs_OHW_40_ws15_gc_corrected_recomp.PCs.extrapolated.txt";
		analyze(proj, ped, pcFile, covar, outputDir, numthreads);
	}

	public static void main(String[] args) {

		int numArgs = args.length;
		String filename = null;
		String ped = null;
		String outputDir = null;
		int numthreads = 24;
		String pcFile = null;
		String covar = null;
		String usage = "\n" + "one.JL.MitoAnalyze requires 0-1 arguments\n";
		usage += "   (1) filename (i.e. proj= (no default))\n" + "";
		usage += "   (2) ped (i.e. ped= (no default))\n" + "";
		usage += "   (3) output directory (i.e. out= (no default))\n" + "";
		usage += "   (4) numthreads (i.e. numthreads=" + numthreads + " ( default))\n" + "";
		usage += "   (5) PC file (i.e. pc= (no default))\n" + "";
		usage += "   (6) covarFile (i.e. cov= (no default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			}
			if (args[i].startsWith("-t")) {
				System.out.println("Running test suite");
				test();
				System.exit(1);
			}
			else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ped=")) {
				ped = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outputDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pc=")) {
				pcFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cov=")) {
				covar = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("numthreads=")) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Project proj = new Project(filename, false);
			analyze(proj, ped, pcFile, covar, outputDir, numthreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
