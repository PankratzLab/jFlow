package seq.analysis;

import java.util.ArrayList;

import seq.manage.BamOps;
import stats.Rscript;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.ext;

public class ExomeDepth {
	private static final String EXOME_COUNTS = "ExomeCount.dafr";
	private static final String MY_COUNTS = "my.counts";

	private String[] allReferenceBAMFiles, allReferenceBAIFiles, analysisBamFiles;
	private String outputDir;
	private String outputRoot;
	private boolean fail;

	private Logger log;

	public ExomeDepth(String[] allReferenceBamFiles, String[] analysisBamFiles, String outputDir, String outputRoot, Logger log) {
		super();
		this.allReferenceBAMFiles = allReferenceBamFiles;
		this.fail = gatherBai();
		this.analysisBamFiles = analysisBamFiles;
		this.outputDir = outputDir;
		this.outputRoot = outputRoot;
		this.log = log;
		if (!fail) {
			fail = Array.countIf(ext.indexLargeFactors(analysisBamFiles, allReferenceBamFiles, true, log, true, false), -1) > 0;
			if (fail) {
				log.reportTimeError("Could not detect all analysis .bam files in the complete reference set");
			}
		}
	}

	private boolean gatherBai() {

		boolean verify = true;
		this.allReferenceBAIFiles = new String[allReferenceBAMFiles.length];
		for (int i = 0; i < allReferenceBAMFiles.length; i++) {
			String bai = ext.rootOf(allReferenceBAMFiles[i], false) + BamOps.BAI_EXT;
			if (!Files.exists(allReferenceBAMFiles[i]) || !Files.exists(bai)) {
				log.reportTimeError("Could not find " + allReferenceBAMFiles[i] + " with corresponding .bai file" + bai);
			} else {
				allReferenceBAIFiles[i] = bai;
			}
		}
		return verify;
	}

	private String generateBamBaiScript() {
		String script = addBaseLoadScript("");
		String[] bamBaiV = generateBamBaiRVectors();
		script += "BAMFILES <- " + bamBaiV[0] + "\n";
		script += "BAIFILES <- " + bamBaiV[1] + "\n";
		return script;
	}

	private boolean generateCountFile() {
		boolean created = true;
		String scriptFile = outputDir + outputRoot + "Exome_Counts.Rscript";
		String script = generateCountsScripts();
		CmdLine.prepareBatchForCommandLine(new String[] { script }, scriptFile, true, log);
		created = CmdLine.runCommandWithFileChecks(new String[] { "Rscript", scriptFile }, "", allReferenceBAMFiles, new String[] { getCountFile() }, true, false, false, log);

		return created;
	}

	private String generateCountsScripts() {
		String script = generateBamBaiScript();
		script += MY_COUNTS + " <- getBamCounts(bed.frame = exons.hg19 , bam.files=BAMFILES, include.chr=TRUE, index.files=BAIFILES )\n";
		script += EXOME_COUNTS + " <- as(" + MY_COUNTS + "[, colnames(" + MY_COUNTS + ")], 'data.frame')\n";
		script += "save(" + EXOME_COUNTS + ",file=\"" + getCountFile() + "\")\n";
		return script;
	}

	private String getCountFile() {
		return outputDir + outputRoot + EXOME_COUNTS + ".Rda";
	}

	private String loadCountFileScript() {
		return "load(\"" + getCountFile() + "\")";
	}

	private String addBaseLoadScript(String script) {
		script += "library(ExomeDepth)\n";
		script += "data(exons.hg19)\n";
		return script;
	}

	private String[] generateBamBaiRVectors() {
		String[] bamBaiV = new String[2];
		bamBaiV[0] = Rscript.generateRVector(allReferenceBAMFiles);
		bamBaiV[1] = Rscript.generateRVector(allReferenceBAIFiles);
		return bamBaiV;
	}

	public static void batchExomeDepth(String bams, String outputDir, String outputRoot, int numBatches, int numthreads, Logger log) {
		String[] allReferenceBamFiles = Files.isDirectory(bams) ? Files.list(bams, BamOps.BAM_EXT, false) : HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, true);
		ArrayList<String[]> batches = Array.splitUpArray(allReferenceBamFiles, numBatches, log);
		for (int i = 0; i < batches.size(); i++) {
			ExomeDepth depth = new ExomeDepth(allReferenceBamFiles, batches.get(i), outputDir, outputRoot, log);
			if (!Files.exists(depth.getCountFile())) {
				log.reportTimeWarning("Did not find " + depth.getCountFile() + ", generating it now (takes a long time)");
				depth.generateCountFile();
				System.exit(1);
			}
			String script = depth.generateBamBaiScript();

			System.out.println(script);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String bams = "bams/";
		String outputDir = null;
		String outputRoot = "ExomeDepth";

		int numBatches = 1;
		int numthreads = 1;
		String logfile = null;
		Logger log;

		String usage = "\n" + "seq.analysis.ExomeDepth requires 0-1 arguments\n";
		usage += "   (1) full path to a directory of or file of bams (i.e. bams=" + bams + " (default))\n" + "";
		usage += "   (2) number of batches to run (i.e. numBatches=" + bams + " (default))\n" + "";
		usage += PSF.Ext.getOutputDirCommand(3, outputDir);
		usage += "   (4) output root command (i.e. root=" + outputRoot + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(5, numthreads);

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("bams=")) {
				bams = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("numBatches=")) {
				numBatches = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
				outputDir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("root")) {
				outputRoot = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
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

			log = new Logger(logfile);
			batchExomeDepth(bams, outputDir, outputRoot, numBatches, numthreads, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
