package org.genvisis.imputation;

import java.io.File;
import java.util.ArrayList;

import org.genvisis.cnv.manage.Resources;
import org.genvisis.common.CmdLine;
import org.genvisis.common.CmdLine.Command;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class Imputation {

	public static class Shapeit {
		private String shapeitBin;
		private int numThreads;

		/**
		 * @param shapeitBin filepath of shapeit binary
		 * @param numThreads number of threads to use for shapeit calls
		 */
		public Shapeit(String shapeitBin, int numThreads) {
			super();
			this.shapeitBin = shapeitBin;
			this.numThreads = numThreads;
		}

		private String convertHapsToVCFCommands(String hapFile, String outFile) {
			StringBuilder commands = new StringBuilder();
			commands.append(shapeitBin + " -convert --input-haps " + hapFile + " --output-vcf " + outFile
											+ " -T " + numThreads);
			commands.append("\n");
			commands.append("rm " + hapFile + "*");
			commands.append("\n");
			return commands.toString();
		}

		public String generateShapeitCommands(String vcfFile, String geneticMapFile, String outFile) {
			String phasedTempRoot = ext.parseDirectoryOfFile(outFile) + "TEMP_" + ext.rootOf(outFile);
			StringBuilder commands = new StringBuilder();
			commands.append(shapeitBin + " -V " + vcfFile + " -M " + geneticMapFile + " -O "
											+ phasedTempRoot + " -T " + numThreads);
			commands.append("\n");
			commands.append(convertHapsToVCFCommands(phasedTempRoot, outFile));
			commands.append("\n");
			return commands.toString();
		}

		public String generateShapeitCommands(String plinkFileroot, boolean pedFile,
																					String geneticMapFile, String outFile) {
			String phasedTempRoot = ext.parseDirectoryOfFile(outFile) + "TEMP_" + ext.rootOf(outFile);
			StringBuilder commands = new StringBuilder();
			commands.append(shapeitBin + " " + (pedFile ? "-P " : "-B ") + plinkFileroot + " -M "
											+ geneticMapFile + " -O " + phasedTempRoot + " -T " + numThreads);
			commands.append("\n");
			commands.append(convertHapsToVCFCommands(phasedTempRoot, outFile));
			commands.append("\n");
			return commands.toString();
		}

	}

	public static class Minimac3 {
		private String minimac3Bin;
		private int numThreads;

		/**
		 * @param minimac3Bin filepath of Minimac3 bin (should be Minimac3-omp if multi-threading
		 * @param numThreads number of threads to use for Minimac3 calls
		 */
		public Minimac3(String minimac3Bin, int numThreads) {
			super();
			this.minimac3Bin = minimac3Bin;
			this.numThreads = numThreads;
		}

		public String generateMinimacCommands(String vcfFile, String outPrefix, String refPanelFile) {
			return minimac3Bin + " --refHaps " + refPanelFile + " --haps " + vcfFile + " --prefix "
						 + outPrefix + " --cpus " + numThreads + "\n";
		}
	}

//	public static void imputeOneChr(String chr, String plinkFileroot, String outDir, int numThreads,
//																	String shapeitBin, String minimac3Bin, String geneticMap,
//																	String refPanel, Logger log) {
//		String chrPlinkroot = outDir + ext.removeDirectoryInfo(plinkFileroot) + "_chr" + chr;
//		String root = ext.rootOf(plinkFileroot, true);
//		String phasedFile = outDir + root + "_phased.vcf.gz";
//		String outRoot = outDir + root + "_imputed";
//
//		ArrayList<CmdLine.Command> commands = new ArrayList<CmdLine.Command>();
//		commands.add(new Command(new String[] {"plink2", "--noweb", "--bfile", plinkFileroot, "--chr",
//																					 "--make-bed", "--out", chrPlinkroot},
//														 new String[] {plinkFileroot + ".bed", plinkFileroot + ".bim",
//																					 plinkFileroot + ".fam"},
//														 new String[] {chrPlinkroot + ".bed", chrPlinkroot + ".bim",
//																					 chrPlinkroot + ".fam"},
//														 ""));
//		commands.add(new Command(new String[] {shapeitBin, "-B", plinkFileroot, "-M", geneticMap, "-O",
//																					 chrPlinkroot},
//														 new String[] {plinkFileroot + ".bed", plinkFileroot + ".bim",
//																					 plinkFileroot + ".fam"},
//														 new String[] {chrPlinkroot + ".bed", chrPlinkroot + ".bim",
//																					 chrPlinkroot + ".fam"},
//														 ""));
//
//		if (!Files.exists(chrPlinkroot)) {
//			log.report("Generating PLINK files for chromosome " + chrs[i]);
//			// CmdLine.runDefaults("plink2 --noweb --bfile " + plinkFileroot + " --chr " + chrs[i] + "
//			// --make-bed" + " --out " + chrPlinkroot, "");
//			commands.append("plink2 --noweb --bfile " + plinkFileroot + " --chr " + chrs[i]
//											+ " --make-bed" + " --out " + chrPlinkroot);
//			commands.append("\n");
//		}
//
//		// if (!Files.exists(chrVcf)) {
//		// log.report("Generating VCF for chrosome " + chrs[i]);
//		// commands.append("plink2 --noweb --bfile " + plinkFileroot + " --recode vcf --out " +
//		// chrPlinkroot);
//		// commands.append("\n");
//		// }
//
//		if (!Files.exists(phasedFile)) {
//			log.report("Generating SHAPEIT commands to phase data for chromosome " + chrs[i]);
//			commands.append(shapeit.generateShapeitCommands(chrPlinkroot, false, geneticMaps[i],
//																											phasedFile));
//		}
//
//		if (!Files.exists(outRoot + ".vcf.gz")) {
//			log.report("Generating Minimac commands to impute data for chromosome " + chrs[i]);
//			commands.append(minimac3.generateMinimacCommands(phasedFile, outRoot, refPanels[i]));
//		}
//
//		if (commands.length() > 0) {
//			Files.qsub(pbs, commands.toString(), 60000, 48, numThreads);
//			qsubCommands.add("qsub " + pbs);
//		}
//	}

	/**
	 * 
	 * @param chrs chromosomes to impute
	 * @param vcfFiles corresponding vcf to each chromosome in chrs
	 * @param outDir directory to generate outputs
	 * @param numThreads
	 * @param log
	 * @return
	 */
//	public static boolean generateImputationScript(String[] chrs, String plinkFileroot, String outDir,
//																								 int numThreads, Logger log) {
//		String shapeitBin;
//		String minimac3Bin;
//		String[] geneticMaps;
//		String[] refPanels;
//		Shapeit shapeit;
//		Minimac3 minimac3;
//		ArrayList<String> qsubCommands;
//		String qsubScriptName;
//
//		if (Files.isWindows()) {
//			log.reportError("Imputation pipeline cannot be run on windows, SHAPEIT and Minimac3 are UNIX programs");
//			return false;
//		}
//		shapeitBin = Resources.BIN_RESOURCE_TYPE.SHAPEIT.getResource().getResource(log);
//		if (shapeitBin == null) {
//			log.reportError("SHAPEIT could not be located or downloaded");
//			return false;
//		}
//		minimac3Bin = Resources.BIN_RESOURCE_TYPE.MINIMAC3.getResource().getResource(log);
//		if (minimac3Bin == null) {
//			log.reportError("Minimac3 could not be located or downloaded");
//			return false;
//		}
//		geneticMaps = new String[chrs.length];
//		refPanels = new String[chrs.length];
//		for (int i = 0; i < chrs.length; i++) {
//			String chr = chrs[i];
//			geneticMaps[i] = Resources.GENOME_CHROMOSOME_RESOURCE_TYPE.GENETIC_MAP.getResource(GENOME_BUILD.HG19,
//																																												 chr)
//																																						.getResource(log);
//			if (geneticMaps[i] == null) {
//				log.reportError("Genetic Map for chr" + chr + " could not be located or downloaded");
//				return false;
//			}
//			refPanels[i] = Resources.GENOME_CHROMOSOME_RESOURCE_TYPE.G1K_PHASE3v5_REF_PANEL.getResource(GENOME_BUILD.HG19,
//																																																	chr)
//																																										 .getResource(log);
//			if (refPanels[i] == null) {
//				log.reportError("Reference panel for chr" + chr + " could not be located or downloaded");
//				return false;
//			}
//		}
//
//		shapeit = new Shapeit(shapeitBin, numThreads);
//		minimac3 = new Minimac3(minimac3Bin, numThreads);
//
//		qsubCommands = new ArrayList<String>();
//
//		new File(outDir).mkdirs();
//
//		for (int i = 0; i < chrs.length; i++) {
//			String chrPlinkroot = outDir + ext.removeDirectoryInfo(plinkFileroot) + "_chr" + chrs[i];
//			// String chrVcf = chrPlinkroot + ".vcf";
//			String root = ext.rootOf(chrPlinkroot, true);
//			String phasedFile = outDir + root + "_phased.vcf.gz";
//			String outRoot = outDir + root + "_imputed";
//			String pbs = outDir + root + "_imputation.pbs";
//
//			StringBuilder commands = new StringBuilder();
//
//
//			if (!Files.exists(chrPlinkroot)) {
//				log.report("Generating PLINK files for chromosome " + chrs[i]);
//				// CmdLine.runDefaults("plink2 --noweb --bfile " + plinkFileroot + " --chr " + chrs[i] + "
//				// --make-bed" + " --out " + chrPlinkroot, "");
//				commands.append("plink2 --noweb --bfile " + plinkFileroot + " --chr " + chrs[i]
//												+ " --make-bed" + " --out " + chrPlinkroot);
//				commands.append("\n");
//			}
//
//			// if (!Files.exists(chrVcf)) {
//			// log.report("Generating VCF for chrosome " + chrs[i]);
//			// commands.append("plink2 --noweb --bfile " + plinkFileroot + " --recode vcf --out " +
//			// chrPlinkroot);
//			// commands.append("\n");
//			// }
//
//			if (!Files.exists(phasedFile)) {
//				log.report("Generating SHAPEIT commands to phase data for chromosome " + chrs[i]);
//				commands.append(shapeit.generateShapeitCommands(chrPlinkroot, false, geneticMaps[i],
//																												phasedFile));
//			}
//
//			if (!Files.exists(outRoot + ".vcf.gz")) {
//				log.report("Generating Minimac commands to impute data for chromosome " + chrs[i]);
//				commands.append(minimac3.generateMinimacCommands(phasedFile, outRoot, refPanels[i]));
//			}
//
//			if (commands.length() > 0) {
//				Files.qsub(pbs, commands.toString(), 60000, 48, numThreads);
//				qsubCommands.add("qsub " + pbs);
//			}
//		}
//		if (qsubCommands.isEmpty()) {
//			return false;
//		}
//		qsubCommands.add(0, "cd " + outDir);
//		qsubScriptName = outDir + "submitImputationqsubs.sh";
//		Files.writeIterable(qsubCommands, qsubScriptName);
//		Files.chmod(qsubScriptName);
//
//		log.report("Pipeline ready, run " + qsubScriptName + " to submit imputation jobs");
//
//		return true;
//	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
