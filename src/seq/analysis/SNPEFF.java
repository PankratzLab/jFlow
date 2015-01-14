package seq.analysis;

import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.ext;

//TODO counts
public class SNPEFF {
	public static final String SNP_EFF_COMMAND = "snpEff=";
	public static final String SNP_EFF_BUILD_COMMAND = "build=";

	public static final String SNP_EFF_NO_ANNO_COMMAND = "-noAnno";

	public static final String SNP_EFF = "snpEff.jar";
	public static final String EFF = "eff";
	public static final String COUNT = "count";
	public static final String COUNT_OUT = "" + COUNT;

	public static final String T = "-t";
	public static final String V = "-v";
	public static final String JAVA = "java";
	public static final String JAR = "-jar";

	public static final String ONLY_CODING = "-onlyCoding";
	public static final String I = "-i";
	public static final String INTERVAL = "-interval";

	public static final String O = "-o";
	public static final String GATK = "gatk";
	public static final String[] BUILDS = { "hg19" };
	public static final String CARROT = ">";

	private String snpEffLocation;
	private boolean fail, verbose, overWriteExistingOutput;
	private Logger log;

	public SNPEFF(String snpEffLocation, boolean verbose, boolean overWriteExistingOutput, Logger log) {
		super();

		this.snpEffLocation = snpEffLocation;
		this.log = log;
		this.verbose = verbose;
		this.overWriteExistingOutput = overWriteExistingOutput;
		this.fail = !verify();
	}

	public boolean isFail() {
		return fail;
	}

	private boolean verify() {
		boolean verify = true;
		if (!Files.exists(snpEffLocation)) {
			verify = false;
			log.reportError("Warning - could not find SNP EFF location " + snpEffLocation);
		}
		if (!Files.exists(snpEffLocation + SNP_EFF)) {
			verify = false;
			log.reportError("Warning - could not find the SNP EFF jar file " + snpEffLocation + SNP_EFF);
		}
		return verify;
	}

	public String getSnpEffLocation() {
		return snpEffLocation;
	}

	public void setSnpEffLocation(String snpEffLocation) {
		this.snpEffLocation = snpEffLocation;
	}

	public boolean runSnpEffCountOnBamDirectory(String inputDirectory, String output, String build, String match, String bedFile, int numThreads) {
		String[] inputBams = Files.toFullPaths(Files.list(inputDirectory, match, false), inputDirectory);
		if (inputBams != null && inputBams.length > 0) {
			return runSnpEFFCount(inputBams, output, build, bedFile, numThreads);
		} else {
			log.reportError("Error - no files ending with " + match + " were found in dirctory " + inputDirectory);
			return false;
		}
	}

	public boolean runSnpEFFCount(String[] inputBams, String output, String build, String bedFile, int numThreads) {
		boolean progress = true;
		String[] command = new String[] { JAVA, JAR, snpEffLocation + SNP_EFF, COUNT };
		String[] inputs = inputBams;
		if (bedFile != null) {
			inputs = Array.concatAll(inputBams, new String[] { bedFile });
			command = Array.concatAll(command, new String[] { INTERVAL, bedFile });
		}
		command = Array.concatAll(command, new String[] { V, build });

		// command = Array.concatAll(command, new String[] { V, build });

		command = Array.concatAll(command, inputBams);
		command = Array.concatAll(command, new String[] { CARROT, output });
		String batFile = ext.addToRoot(output, ".bat");
		Files.write(Array.toStr(command, " "), batFile);
		Files.chmod(batFile);
		progress = CmdLine.runCommandWithFileChecks(new String[] { batFile }, "", inputs, new String[] { output }, verbose, overWriteExistingOutput, false, log);
		return progress;
	}

	public SnpEffResult annotateAVCF(String vcfFile, String build) {
		SnpEffResult snpEffResult = new SnpEffResult(vcfFile, log);
		snpEffResult.parse();
		boolean progress = runSnpEFF(snpEffResult, build);
		snpEffResult.setFail(!progress);
		return snpEffResult;
	}

	private boolean runSnpEFF(SnpEffResult snpEffResult, String build) {
		boolean progress = true;
		String[] inputFiles = new String[] { snpEffResult.getInputVCF() };
		String[] outputFiles = new String[] { snpEffResult.getOutputSnpEffVCF() };
		String[] command = new String[] { JAVA, JAR, snpEffLocation + SNP_EFF, V, O, GATK, build, snpEffResult.getInputVCF(), CARROT, snpEffResult.getOutputSnpEffVCF() };
		String batFile = ext.addToRoot(snpEffResult.getInputVCF(), ".bat");
		Files.write(Array.toStr(command, " "), batFile);
		Files.chmod(batFile);

		progress = CmdLine.runCommandWithFileChecks(new String[] { batFile }, "", inputFiles, outputFiles, verbose, overWriteExistingOutput, false, log);
		return progress;

	}

	public static class SnpEffResult {
		public static final String OUT_EFF = ".eff";
		public static final String OUT_GATK = ".gatk";

		private String vcfFile;
		private String outputSnpEffVCF;
		private String outputGatkSnpEffVCF;
		private Logger log;
		private boolean fail;

		public SnpEffResult(String inputVCF, Logger log) {
			super();
			this.vcfFile = inputVCF;
			this.log = log;
		}

		public void parse() {
			this.outputSnpEffVCF = ext.addToRoot(vcfFile, OUT_EFF);
			this.outputGatkSnpEffVCF = ext.addToRoot(outputSnpEffVCF, OUT_GATK);
		}

		public String getInputVCF() {
			return vcfFile;
		}

		public void setInputVCF(String inputVCF) {
			this.vcfFile = inputVCF;
		}

		public String getOutputSnpEffVCF() {
			return outputSnpEffVCF;
		}

		public void setOutputSnpEffVCF(String outputSnpEffVCF) {
			this.outputSnpEffVCF = outputSnpEffVCF;
		}

		public String getOutputGatkSnpEffVCF() {
			return outputGatkSnpEffVCF;
		}

		public void setOutputGatkSnpEffVCF(String outputGatkSnpEffVCF) {
			this.outputGatkSnpEffVCF = outputGatkSnpEffVCF;
		}

		public boolean isFail() {
			return fail;
		}

		public Logger getLog() {
			return log;
		}

		public void setLog(Logger log) {
			this.log = log;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

	}

	public static void testCounts(String snpEffLocation, String inputDirectory, String build, String match, String bedFile, int numThreads, Logger log) {
		SNPEFF snpeff = new SNPEFF(snpEffLocation, true, true, log);
		snpeff.runSnpEffCountOnBamDirectory(inputDirectory, inputDirectory + "SnpEff.counts", build, match, bedFile, numThreads);
	}

	public static void main(String[] args) {
		String snpEffLocation = "/home/pankrat2/public/bin/snpEff/";
		String inputDirectory = "/home/tsaim/shared/Project_Tsai_Project_021/bam/";
		String bedFile = "/home/tsaim/lane0212/bin/ref/S04380219_Regions.bed";
		String build = BUILDS[0];
		String match = ".bam";
		int numThreads = 2;
		Logger log = new Logger(inputDirectory + "snpEff.log");
		testCounts(snpEffLocation, inputDirectory, build, match, bedFile, numThreads, log);

		// String usage = "\n" + "seq.analysis.SNPEFF requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";
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
