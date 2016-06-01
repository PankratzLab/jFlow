package one.JL;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.Callable;

import common.Array;
import common.CmdLine;
import common.ExcelConverter;
import common.Files;
import common.HashVec;
import common.Logger;
import common.WorkerHive;
import common.ext;
import cnv.analysis.FilterCalls;
import filesys.CNVariant;
import filesys.GeneData;
import filesys.GeneTrack;
import filesys.LocusSet;
import filesys.Segment;

public class PlinkCNV {

	private static void run(String cnvaFile, String sampFile) {

		final String dir = ext.parseDirectoryOfFile(cnvaFile);
		String merge = dir + "penncnv.cnv";
		if (!Files.exists(merge)) {
			FilterCalls.mergeCNVs(cnvaFile, merge, FilterCalls.DEFAULT_CLEAN_FACTOR, null);
		}
		final Logger log = new Logger(dir + "cnv.log");
		String[] filters = Files.list(dir, ".crf", false);

		WorkerHive<PlinkResult> hive = new WorkerHive<PlinkResult>(6, 10, log);

		for (int filt = 0; filt < filters.length; filt++) {
			String[][] phenos = HashVec.loadFileToStringMatrix(dir + "pheno.dat", false, null, false);
			final String filtFile = dir + filters[filt];
			String out = dir + ext.rootOf(filters[filt]) + ".cnv";

			if (!Files.exists(out)) {
				FilterCalls.fromParameters(dir + filters[filt], new Logger());
			}
			String sampCNVs = ext.addToRoot(out, "." + ext.rootOf(sampFile));
			HashSet<String> sampSet = HashVec.loadFileToHashSet(sampFile, false);
			if (!Files.exists(sampCNVs)) {
				LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(out, log);

				try {
					int numTotal = 0;
					int numSamps = 0;
					PrintWriter writer = new PrintWriter(new FileWriter(sampCNVs));
					writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
					for (int i = 0; i < cnvs.getLoci().length; i++) {
						numTotal++;
						if (sampSet.contains(cnvs.getLoci()[i].getIndividualID())) {
							writer.println(cnvs.getLoci()[i].toPlinkFormat());
							numSamps++;
						}
					}
					writer.close();
					log.reportTimeInfo(numTotal + " - > " + numSamps);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}

			for (int i = 2; i < phenos[0].length; i++) {

				final String pheno = phenos[0][i];
				final String key = phenos[0][i] + "_" + ext.rootOf(sampFile) + "_" + ext.rootOf(filters[filt]);
				final String opDir = dir + key + "/";
				boolean quant = false;

				new File(opDir).mkdirs();
				String phenoCNV = opDir + pheno + ".cnv";
				if (!Files.exists(phenoCNV)) {
					Files.copyFile(sampCNVs, phenoCNV);
				}
				ArrayList<String> fam = new ArrayList<String>();
				ArrayList<String> quants = new ArrayList<String>();
				HashSet<String> types = new HashSet<String>();
				for (int j = 1; j < phenos.length; j++) {
					if (sampSet.contains(phenos[j][0])) {
						fam.add(phenos[j][0] + "\t" + phenos[j][0] + "\t" + 0 + "\t" + 0 + "\t" + phenos[j][1] + "\t" + phenos[j][i]);
						quants.add(phenos[j][0] + "\t" + phenos[j][0] + "\t" + phenos[j][i]);
					}
					try {
						double val = Double.parseDouble(phenos[j][i]);
						if (val != 1 && val != 2)
							types.add(val + "");
					} catch (NumberFormatException nfe) {

					}
				}
				quant = types.size() > 0;
				if (quant) {
					log.reportTimeInfo("Detected phenotype " + pheno + " is quantatative");
				} else {
					log.reportTimeInfo("Detected phenotype " + pheno + " is discrete");
				}

				Files.writeArrayList(fam, opDir + pheno + ".fam");
				Files.writeArrayList(quants, opDir + pheno + ".qPheno.txt");

				String map = opDir + pheno + ".cnv.map";
				if (!Files.exists(map)) {
					CmdLine.run(dir + "plink --cnv-list " + pheno + ".cnv --cnv-make-map --out " + pheno, opDir);
				}

				if (quant) {
					Callable<PlinkResult> c1 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {
							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--pheno");
							cmd.add(opDir + pheno + ".qPheno.txt");
							cmd.add("--mperm");
							cmd.add("10000");
							cmd.add("--out");
							cmd.add(opDir + pheno);
							String out = opDir + pheno + ".cnv.qt.summary.mperm";
							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, false, true, log);
							return new PlinkResult(filtFile, key, out, Array.toStr(Array.toStringArray(cmd), " "), complete, true);
						}
					};
					hive.addCallable(c1);

				} else {
					Callable<PlinkResult> c1 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {
							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--cnv-indiv-perm");
							cmd.add("--mperm");
							cmd.add("10000");
							cmd.add("--out");
							cmd.add(opDir + pheno);
							String out = opDir + pheno + ".cnv.summary.mperm";
							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, false, true, log);

							return new PlinkResult(filtFile, key, out, Array.toStr(Array.toStringArray(cmd), " "), complete, false);
						}
					};
					Callable<PlinkResult> c2 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {
							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--mperm");
							cmd.add("10000");
							cmd.add("--out");
							cmd.add(opDir + pheno + "_position");
							String out = opDir + pheno + "_position.cnv.summary.mperm";

							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, false, true, log);
							return new PlinkResult(filtFile, key + "_pos", out, Array.toStr(Array.toStringArray(cmd), " "), complete, false);

						}
					};
					Callable<PlinkResult> c3 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {
							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--mperm");
							cmd.add("10000");
							cmd.add("--cnv-test-window");
							cmd.add("200");
							cmd.add("--out");
							cmd.add(opDir + pheno + "_window");
							String out = opDir + pheno + "_window.cnv.summary.mperm";

							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, false, true, log);
							return new PlinkResult(filtFile, key + "_window", out, Array.toStr(Array.toStringArray(cmd), " "), complete, false);

						}
					};
					hive.addCallable(c1);
					hive.addCallable(c2);
					hive.addCallable(c3);
				}
			}

		}
		hive.execute(true);
		ArrayList<PlinkResult> results = hive.getResults();
		ArrayList<String> filesToCombine = new ArrayList<String>();
		GeneTrack geneTrack = GeneTrack.load("N:/statgen/NCBI/RefSeq_hg18.gtrack", false);

		String finalDir = dir + ext.rootOf(sampFile) + "_results/";
		filesToCombine.add(finalDir + "key.txt");
		ArrayList<String> key = new ArrayList<String>();
		ArrayList<String> combo = new ArrayList<String>();
		String comboOut = finalDir + "combo.txt";
		String allSigGenes = finalDir + "sigGenes.txt";
		HashSet<String> allSigGenesList = new HashSet<String>();
		String allSigGenes5k = finalDir + "sigGenes5k.txt";
		HashSet<String> allSigGenesList5k = new HashSet<String>();

		String allSigGenesEMP2 = finalDir + "sigGenesEMP2_0_5.txt";
		HashSet<String> allSigGenesListEMP2 = new HashSet<String>();

		filesToCombine.add(comboOut);
		filesToCombine.add(allSigGenesEMP2);
		filesToCombine.add(allSigGenes);
		filesToCombine.add(allSigGenes5k);

		combo.add("Type\tCHR\tSNP\tEMP1\tEMP2\tUCSC\tGENE\tGENE_5k\tLink_buffer\tAFF\tUNAFF\tNCNV\tM0\tM1\tCMD");
		new File(finalDir).mkdirs();
		int index = 1;
		for (PlinkResult result : results) {
			if (result.complete) {
				String out = finalDir + result.key.substring(0, 26) + "_" + index + ".txt";
				String[][] data = HashVec.loadFileToStringMatrix(result.file, false, null, "[\\s]+", false, 1000, true);
				log.reportTimeInfo("Loading paired " + ext.rootOf(result.file, false));
				String[][] dataCount = HashVec.loadFileToStringMatrix(ext.rootOf(result.file, false), false, null, "[\\s]+", false, 1000, true);
				String allSigGenesSpecific = finalDir + result.key + "_Genes" + index + ".txt";
				HashSet<String> allSigGenesListSpecific = new HashSet<String>();

				ArrayList<String> toReport = new ArrayList<String>();
				toReport.add("Type\t" + Array.toStr(data[0]) + "\tUCSC\tGENE\tGENE_5k\tLink_buffer\tAFF\tUNAFF\tNCNV\tM0\tM1\tCMD");

				for (int i = 1; i < data.length; i++) {
					if (data[i].length > 3 && Double.parseDouble(data[i][2]) < 0.05) {
						int loc = Integer.parseInt(data[i][1].replaceAll("p" + data[i][0] + "-", ""));
						Segment seg = new Segment(Byte.parseByte(data[i][0]), loc, loc);
						GeneData[] geneDatas = geneTrack.getOverlappingGenes(seg);
						GeneData[] geneDatasBuff = geneTrack.getOverlappingGenes(seg.getBufferedSegment(5000));
						if (!data[i][1].equals(dataCount[i][1])) {
							throw new IllegalArgumentException("mismatched files");
						}
						StringBuilder builder = new StringBuilder(result.key + "\t" + Array.toStr(data[i]) + "\t" + seg.getUCSClocation() + "\t");
						for (int j = 0; j < geneDatas.length; j++) {
							builder.append((j == 0 ? "" : ":") + geneDatas[j].getGeneName());
							allSigGenesList.add(geneDatas[j].getGeneName() + "\t" + geneDatas[j].getUCSClocation() + "\t" + result.key);
							allSigGenesListSpecific.add(geneDatas[j].getGeneName() + "\t" + geneDatas[j].getUCSClocation() + "\t" + result.key);
							if (Double.parseDouble(data[i][3]) < 1) {
								allSigGenesListEMP2.add(geneDatas[j].getGeneName() + "\t" + geneDatas[j].getUCSClocation() + "\t" + result.key);
							}
						}
						builder.append("\t");
						for (int j = 0; j < geneDatasBuff.length; j++) {
							builder.append((j == 0 ? "" : ":") + geneDatasBuff[j].getGeneName());
							allSigGenesList5k.add(geneDatasBuff[j].getGeneName() + "\t" + geneDatasBuff[j].getUCSClocation() + "\t" + result.key);
						}
						builder.append("\t" + seg.getBufferedSegment(5000).getUCSCLink("hg18"));
						if (result.quant) {
							builder.append("\tNA\tNA" + "\t" + dataCount[i][3] + "\t" + dataCount[i][4] + "\t" + dataCount[i][5]);

						} else {
							builder.append("\t" + dataCount[i][3] + "\t" + dataCount[i][4] + "\tNA\tNA\tNA");
						}
						builder.append("\t" + result.command);
						toReport.add(builder.toString());
						combo.add(builder.toString());
					}
				}
				if (toReport.size() > 1) {

					key.add(result.key + "\t" + index);
					index++;
					Files.writeArrayList(toReport, out);
					Files.writeHashSet(allSigGenesListSpecific, allSigGenesSpecific);
					filesToCombine.add(out);
					// filesToCombine.add(result.filtFile);
				}
			}
		}

		Files.writeArrayList(key, finalDir + "key.txt");
		Files.writeArrayList(combo, comboOut);
		Files.writeHashSet(allSigGenesList, allSigGenes);
		Files.writeHashSet(allSigGenesList5k, allSigGenes5k);
		Files.writeHashSet(allSigGenesListEMP2, allSigGenesEMP2);

		ExcelConverter excelConverter = new ExcelConverter(filesToCombine, finalDir + "p_0_05.xlsx", log);
		excelConverter.convert(true);
	}

	private static class PlinkResult {
		private String filtFile;
		private String file;
		private String key;
		private String command;
		private boolean quant;
		private boolean complete;

		public PlinkResult(String filtFile, String key, String file, String command, boolean complete, boolean quant) {
			super();
			this.filtFile = filtFile;
			this.key = key;
			this.file = file;
			this.command = command;
			this.complete = complete;
			this.quant = quant;
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "C:/data/ARIC/shadowCNVs/";
		String cnvFile = dir + "combinedMF.cnv";
		String[] sampFiles = new String[] { dir + "whites.txt", dir + "all.txt", dir + "casesOnly.txt" };

		String usage = "\n" +
				"one.JL.ARICCNV requires 0-1 arguments\n" +
				"   (1) cnvs (i.e. cnvs=" + cnvFile + " (default))\n" +
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				cnvFile = args[i].split("=")[1];
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
			for (int i = 0; i < sampFiles.length; i++) {
				run(cnvFile, sampFiles[i]);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
