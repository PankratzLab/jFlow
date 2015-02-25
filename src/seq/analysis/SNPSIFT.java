package seq.analysis;

import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.PSF;
import common.ext;

public class SNPSIFT {
	public static final String SNP_SIFT_LOCATION_COMMAND = "snpSift=";
	public static final String SNP_SIFT_JAR = "SnpSift.jar";
	public static final String DB_NSFP = "dbnsfp";
	public static final String TSTV = "tstv";
	public static final String HOM = "hom";

	public static final String V = "-v";

	private String snpSiftLocation;
	private boolean verbose, fail, overWriteExisting;
	private Logger log;

	public SNPSIFT(String snpSiftLocation, boolean verbose, boolean overWriteExisting, Logger log) {
		super();
		this.snpSiftLocation = snpSiftLocation;
		this.verbose = verbose;
		this.overWriteExisting = overWriteExisting;
		this.log = log;
		this.fail = !verify();
	}

	public String getSnpSiftLocation() {
		return snpSiftLocation;
	}

	public SnpSiftResult annotateDbnsfp(String inputVCF, Logger log) {
		boolean progress = true;
		SnpSiftResult snSiftResult = new SnpSiftResult(inputVCF, verbose, log);
		snSiftResult.parse();
		progress = dbnsfpAVCF(snSiftResult.getInputVCF(), snSiftResult.getOutputVCF(), log);
		snSiftResult.setFail(!progress);
		return snSiftResult;
	}

	public SnpSiftResult tsTv(String inputVCF, Logger log) {
		boolean progress = true;
		SnpSiftResult snSiftResult = new SnpSiftResult(inputVCF, verbose, log);
		snSiftResult.parse();
		progress = tsTvVcf(snSiftResult.getInputVCF(), snSiftResult.getOutputTsTv(), log);
		snSiftResult.setFail(!progress);
		return snSiftResult;
	}

	private boolean dbnsfpAVCF(String inputVCF, String outputVCF, Logger log) {
		boolean progress = true;
		if (!fail) {
			String[] command = PSF.Java.buildJavaJar(snpSiftLocation + SNP_SIFT_JAR);
			String[] args = new String[] { DB_NSFP, V, inputVCF, PSF.Ext.CARROT, outputVCF };
			command = Array.concatAll(command, args);
			String batFile = outputVCF + ".bat";
			command = CmdLine.prepareBatchForCommandLine(command, batFile, true, log);
			progress = CmdLine.runCommandWithFileChecks(command, "", new String[] { inputVCF, batFile }, new String[] { outputVCF }, verbose, overWriteExisting, false, log);
		} else {
			progress = false;
		}
		return progress;
	}

	private boolean tsTvVcf(String inputVCF, String outputTsTv, Logger log) {
		boolean progress = true;
		if (!fail) {
			String[] command = PSF.Java.buildJavaJar(snpSiftLocation + SNP_SIFT_JAR);
			String[] args = new String[] { TSTV, HOM, inputVCF, PSF.Ext.CARROT, outputTsTv };
			command = Array.concatAll(command, args);
			String batFile = outputTsTv + ".bat";
			command = CmdLine.prepareBatchForCommandLine(command, batFile, true, log);
			progress = CmdLine.runCommandWithFileChecks(command, "", new String[] { inputVCF, batFile }, new String[] { outputTsTv }, verbose, overWriteExisting, false, log);
		} else {
			progress = false;
		}
		return progress;
	}

	private boolean verify() {
		boolean verify = true;
		if (!Files.exists(snpSiftLocation)) {
			verify = false;
			log.reportError("Warning - could not find SNP EFF location " + snpSiftLocation);
		}
		if (!Files.exists(snpSiftLocation + SNP_SIFT_JAR)) {
			verify = false;
			log.reportError("Warning - could not find the SNP EFF jar file " + snpSiftLocation + SNP_SIFT_JAR);
		}
		return verify;
	}

	public static class SnpSiftResult {
		public static final String[] TYPES = { ".dbnsfp", ".tstv" };
		private String inputVCF, outputVCF, outputTsTv;
		private boolean verbose, fail;
		private Logger log;

		public SnpSiftResult(String inputVCF, boolean verbose, Logger log) {
			super();
			this.inputVCF = inputVCF;
			this.verbose = verbose;
			this.log = log;

		}

		public void parse() {
			this.outputVCF = ext.addToRoot(inputVCF, TYPES[0]);
			this.outputTsTv = ext.rootOf(inputVCF, false) + TYPES[1];
		}

		public boolean isVerbose() {
			return verbose;
		}

		public Logger getLog() {
			return log;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public String getInputVCF() {
			return inputVCF;
		}

		public String getOutputVCF() {
			return outputVCF;
		}

		public String getOutputTsTv() {
			return outputTsTv;
		}
	}

}
