package seq.analysis;

import seq.manage.VCFOps;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.PSF;

/**
 * Class for automating the annotation of a vcf using ANNOVAR
 *
 */
public class ANNOVAR {
	public static final String ANNOVAR_COMMAND = "annovar=";
	public static final String TABLE_ANNOVAR = "table_annovar.pl";
	private static final String PROTOCOL = "-protocol";
	private static final String DEFAULT_PROTOCOLS = "refGene,cytoBand,genomicSuperDups,esp6500si_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,popfreq_max_20150413,popfreq_all_20150413,esp6500siv2_all,esp6500siv2_aa,esp6500siv2_ea,cosmic70,dbnsfp30a";
	private static final String OPERATION = "-operation";
	private static final String DEFAULT_OPERATIONS = "g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f";
	private static final String REMOVE = "-remove";
	private static final String DEFUALT_ANNOVAR_DB = "humandb/";
	private static final String BUILD_VERSION = "-buildver";
	private static final String OUT = "-out";
	private static final String NA_STRING = "-nastring";
	private static final String DEFAULT_NA_STRING = ".";
	private static final String VCF_INPUT = "-vcfinput";
	private static final String MULTI_ANNO = "multianno.vcf";
	private static final String THREAD ="--thread";

	private String annovarLocation;
	private String annodvarDB;
	private boolean fail, verbose, overWriteExistingOutput;
	private Logger log;

	public ANNOVAR(String annovarLocation, boolean verbose, boolean overWriteExistingOutput, Logger log) {
		super();
		this.annovarLocation = annovarLocation;
		this.verbose = verbose;
		this.overWriteExistingOutput = overWriteExistingOutput;
		this.log = log;
		this.fail = !verify();
	}

	private boolean verify() {
		boolean verify = true;
		if (annodvarDB == null) {
			// we don't care for now
		}
		if (!Files.exists(annovarLocation)) {
			verify = false;
			log.reportError("Warning - could not find Annovar directory  " + annovarLocation);
		}
		if (!Files.exists(annovarLocation + TABLE_ANNOVAR)) {
			verify = false;
			log.reportError("Warning - could not find " + annovarLocation + TABLE_ANNOVAR);
		}
		if (!Files.exists(annovarLocation + DEFUALT_ANNOVAR_DB)) {
			verify = false;
			log.reportError("Warning - could not find " + annovarLocation + DEFUALT_ANNOVAR_DB);
		}
		return verify;
	}

	public String getAnnovarLocation() {
		return annovarLocation;
	}

	public boolean isFail() {
		return fail;
	}

	public AnnovarResults AnnovarAVCF(String inputVCF, String build,  int numthreads,Logger log) {
		AnnovarResults annovarResults = new AnnovarResults(inputVCF, build);
		annovarResults.parse();
		boolean progress = AnnovarAVCF(annovarResults.getInputVCF(), annovarResults.getOutput(), annovarResults.getOutputVCF(), annovarResults.getBuild(), numthreads, log);
		annovarResults.setFail(!progress);
		if (!progress) {
			log.reportTimeError("The annovar has failed :( \n This is most likely caused by not having the required database files. Please see \"http://www.openbioinformatics.org/annovar/annovar_startup.html\" for commands to download the necessary files");
		}
		return annovarResults;
	}

	private boolean AnnovarAVCF(String inputVCF, String outputBase, String outputVCF, String build, int numthreads, Logger log) {
		boolean progress = !fail;
		if (progress) {
			String[] command = new String[] { PSF.Cmd.PERL, annovarLocation + TABLE_ANNOVAR, inputVCF, annovarLocation + DEFUALT_ANNOVAR_DB, BUILD_VERSION, build, OUT, outputBase, THREAD, numthreads + "", REMOVE, PROTOCOL, DEFAULT_PROTOCOLS, OPERATION, DEFAULT_OPERATIONS, NA_STRING, DEFAULT_NA_STRING, VCF_INPUT };
			progress = CmdLine.runCommandWithFileChecks(command, "", new String[] { inputVCF }, new String[] { outputVCF }, verbose, overWriteExistingOutput, false, log);
		}
		return progress;
	}

	public static class AnnovarResults {
		private String inputVCF;
		private String build;
		private String output;
		private String outputVCF;
		private boolean fail;
		private Logger log;

		public AnnovarResults(String inputVCF, String build) {
			super();
			this.inputVCF = inputVCF;
			this.build = build;
			this.fail = false;
		}

		public void parse() {
			this.output = VCFOps.getAppropriateRoot(inputVCF, false);
			this.outputVCF = output + "." + build + "_" + MULTI_ANNO;
		}

		public String getInputVCF() {
			return inputVCF;
		}

		public String getBuild() {
			return build;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public String getOutput() {
			return output;
		}

		public String getOutputVCF() {
			return outputVCF;
		}

		public Logger getLog() {
			return log;
		}

	}

}
