package cnv.manage;

import gwas.PhenoPrep;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import common.Array;
import common.CmdLine;
import common.Files;
import common.ext;
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
	public static void analyze(Project proj, String pedFile, String pcFile, String covarFile, String outputDir) {
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
			PlinkData.saveGenvisisToPlinkPedSet(proj, root, "", null, null);
		} else {
			proj.getLog().reportTimeInfo(plinkMap + " and " + plinkPed + "exist, skipping");
		}

		ArrayList<String> plinkConverCommand = new ArrayList<String>();
		plinkConverCommand.add("plink2");
		plinkConverCommand.add("--file");
		plinkConverCommand.add(root);
		plinkConverCommand.add("--make-bed");
		plinkConverCommand.add("--out");
		plinkConverCommand.add(root);
		String[] in = new String[] { plinkPed, plinkMap };
		String fam = root + ".fam";
		String bim = root + ".bim";
		String bed = root + ".bed";

		String[] out = new String[] { fam, bim, bed };
		CmdLine.runCommandWithFileChecks(Array.toStringArray(plinkConverCommand), "", in, out, true, false, false, proj.getLog());

		String subDir = ext.rootOf(pcFile, false) + SUB_DIR;
		String natty = subDir + "WITHOUT_BUILDERS_NATURAL_WITHOUT_INDEPS_finalSummary.estimates.txt.gz";
		String[] nattyTitles = new String[] { "PC15", "PC150" };
		runGwas(proj, root, natty, nattyTitles, covarFile, ped);

		String stepWise = subDir + "WITHOUT_BUILDERS_STEPWISE_RANK_R2_WITHOUT_INDEPS_finalSummary.estimates.txt.gz";
		String[] header = Files.getHeaderOfFile(stepWise, proj.getLog());
		String[] stepWiseTitles = new String[] { "PC15", header[header.length - 1] };
		runGwas(proj, root, stepWise, stepWiseTitles, covarFile, ped);

	}

	private static void runGwas(Project proj, String root, String mtPhenoFile, String[] titles, String covarFile, Pedigree ped) {
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
			for (int i = 0; i < titles.length; i++) {
				String outCurrent = root + ext.rootOf(mtPhenoFile) + "_" + titles[i] + ".txt";
				for (int j = 0; j < covarParser.getNumericDataTitles().length; j++) {
					outCurrent = ext.addToRoot(outCurrent, covarParser.getNumericDataTitles()[j]);
				}
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(outCurrent));
					writer.print("FID\tIID\t" + titles[i]);
					for (int j = 0; j < covarParser.getNumericDataTitles().length; j++) {
						writer.print("\t" + covarParser.getNumericDataTitles()[j]);
					}
					writer.println();
					for (int j = 0; j < ped.getDnas().length; j++) {
						int sampIndex = ext.indexOfStr(ped.getDnas()[j], proj.getSamples());
						writer.print(ped.getFID(j) + "\t" + ped.getIID(j) + "\t" + parser.getNumericDataForTitle(titles[i])[sampIndex]);
						for (int k = 0; k < covarParser.getNumericDataTitles().length; k++) {
							writer.print("\t" + covarParser.getNumericDataForTitle(covarParser.getNumericDataTitles()[k])[sampIndex]);
						}
						writer.println();
					}
					writer.close();
					String outPrepReg = ext.addToRoot(outCurrent, ".prepped");
					PhenoPrep.parse("", outCurrent, "IID", titles[i], null, 3.0, false, false, true, false, false, Array.toStr(covarParser.getNumericDataTitles(), ","), root + ".fam", true, true, false, false, false, false, null, outPrepReg, true, false, false, false, false, false, proj.getLog());
					String outPrepInv = ext.addToRoot(outCurrent, ".prepped.inv");
					PhenoPrep.parse("", outCurrent, "IID", titles[i], null, 3.0, false, false, true, false, true, Array.toStr(covarParser.getNumericDataTitles(), ","), root + ".fam", true, true, false, false, false, false, null, outPrepInv, true, false, false, false, false, false, proj.getLog());

				} catch (Exception e) {
					proj.getLog().reportError("Error writing to " + outCurrent);
					proj.getLog().reportException(e);
				}
			}

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static void test() {
		Project proj = new Project("/home/pankrat2/lanej/projects/gedi_gwas.properties", false);
		String ped = "/panfs/roc/groups/5/pankrat2/shared/gedi_gwas/mitoAnalyze/gedi_gwas.ped";
		String covar = "/panfs/roc/groups/5/pankrat2/shared/gedi_gwas/mitoAnalyze/gedi_gwas.covar";

		String outputDir = "/scratch.global/lanej/mitoAnalyze/";
		String pcFile = proj.PROJECT_DIRECTORY.getValue() + "gedi_gwasALL_1000PCs_OHW_40_ws15_gc_corrected_recomp.PCs.extrapolated.txt";
		analyze(proj, ped, pcFile, covar, outputDir);
	}

	public static void main(String[] args) {
		test();
		// int numArgs = args.length;
		// String filename = null;
		// String logfile = null;
		//
		// String usage = "\n" + "one.JL.MitoAnalyze requires 0-1 arguments\n";
		// usage += "   (1) filename (i.e. proj=" + filename + " (default))\n" + "";
		//
		// for (int i = 0; i < args.length; i++) {
		// if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
		// System.err.println(usage);
		// System.exit(1);
		// } else if (args[i].startsWith("file=")) {
		// filename = args[i].split("=")[1];
		// numArgs--;
		// } else if (args[i].startsWith("log=")) {
		// logfile = args[i].split("=")[1];
		// numArgs--;
		// } else {
		// System.err.println("Error - invalid argument: " + args[i]);
		// }
		// }
		// if (numArgs != 0) {
		// System.err.println(usage);
		// System.exit(1);
		// }
		// try {
		// log = new Logger(logfile);
		// parse(filename, log);
		// } catch (Exception e) {
		// e.printStackTrace();
		// }
	}

}
