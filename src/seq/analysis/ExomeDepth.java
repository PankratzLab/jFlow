package seq.analysis;

import java.util.ArrayList;
import java.util.concurrent.Callable;

import seq.manage.BamOps;
import stats.Rscript;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;

public class ExomeDepth {
	private static final String EXOME_COUNTS_DAFR = "ExomeCount.dafr";
	private static final String MY_COUNTS_VAR = "my.counts";
	private static final String MY_TEST_VAR = "my.test";
	private static final String MY_REF_SAMPLES_VAR = "my.ref.samples";
	private static final String MY_REF_SET_VAR = "my.reference.set";
	private static final String MY_REF_SET_SELECTED_VAR = "my.reference.selected";
	private static final String MY_MATRIX_VAR = "my.matrix";
	private static final String MY_ALL_EXONS_VAR = "all.exons";
	private static final String MY_CHOICE_VAR = "my.choice";

	private String[] allReferenceBAMFiles, allReferenceBAIFiles, analysisBamFiles;
	private String outputDir;
	private String outputRoot;
	private boolean fail;

	private Logger log;

	public ExomeDepth(String[] allReferenceBamFiles, String[] analysisBamFiles, String outputDir, String outputRoot, Logger log) {
		super();
		this.allReferenceBAMFiles = allReferenceBamFiles;
		this.log = log;

		this.fail = gatherBai();
		this.analysisBamFiles = analysisBamFiles;
		this.outputDir = outputDir;
		this.outputRoot = outputRoot;
		if (!fail) {
			fail = Array.countIf(ext.indexLargeFactors(analysisBamFiles, allReferenceBamFiles, true, log, true, false), -1) > 0;
			if (fail) {
				log.reportTimeError("Could not detect all analysis .bam files in the complete reference set");
			}
		}
	}

	public String[] getAnalysisBamFiles() {
		return analysisBamFiles;
	}

	public String[] getAllReferenceBAMFiles() {
		return allReferenceBAMFiles;
	}

	public String getOutputDir() {
		return outputDir;
	}

	public String getOutputRoot() {
		return outputRoot;
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

	/**
	 * @param eAnalysis
	 *            the script property is modified
	 * @return
	 */
	private ExomeDepthAnalysis generateCallingScript(final ExomeDepthAnalysis eAnalysis) {
		String script = addBaseLoadScript("");
		script += loadCountFileScript() + "\n";
		script += MY_TEST_VAR + " <- " + EXOME_COUNTS_DAFR + "$" + ext.removeDirectoryInfo(eAnalysis.getInputBam()) + "\n";
		String[] tmpRef = new String[allReferenceBAMFiles.length - 1];
		int index = 0;
		for (int i = 0; i < allReferenceBAMFiles.length; i++) {
			if (!allReferenceBAMFiles[i].equals(eAnalysis.getInputBam())) {
				tmpRef[index] = ext.removeDirectoryInfo(allReferenceBAMFiles[i]);
				index++;
			}
		}

		script += MY_REF_SAMPLES_VAR + " <- " + Rscript.generateRVector(tmpRef) + "\n";
		script += MY_REF_SET_VAR + " <- " + "as.matrix(" + EXOME_COUNTS_DAFR + "[, " + MY_REF_SAMPLES_VAR + "])\n";

		script += MY_CHOICE_VAR + " <- " + "select.reference.set (test.counts = " + MY_TEST_VAR + ",reference.counts = " + MY_REF_SET_VAR;
		script += ",bin.length = (" + EXOME_COUNTS_DAFR + "$end - " + EXOME_COUNTS_DAFR + "$start)/1000,n.bins.reduced = 10000)\n";
		script += "print(" + MY_CHOICE_VAR + "[[1]])\n";

		script += MY_MATRIX_VAR + " <- as.matrix( " + EXOME_COUNTS_DAFR + "[," + MY_CHOICE_VAR + "$reference.choice, drop = FALSE])\n";
		script += MY_REF_SET_SELECTED_VAR + "<- apply(X =" + MY_MATRIX_VAR + ",MAR = 1,FUN = sum)\n";

		script += MY_ALL_EXONS_VAR + "<- new ('ExomeDepth', test = " + MY_TEST_VAR + " , reference = " + MY_REF_SET_SELECTED_VAR + ",";
		script += "formula = 'cbind(test,reference) ~ 1')\n";

		script += MY_ALL_EXONS_VAR + "<- CallCNVs(x = " + MY_ALL_EXONS_VAR + " ,";
		script += "chromosome = " + EXOME_COUNTS_DAFR + "$space,";
		script += "start = " + EXOME_COUNTS_DAFR + "$start,";
		script += "end = " + EXOME_COUNTS_DAFR + "$end,";
		script += "name = " + EXOME_COUNTS_DAFR + "$names)\n";
		script += "write.table(" + MY_ALL_EXONS_VAR + "@CNV.calls, " + "\"" + eAnalysis.getExomeDepthOutput() + "\", sep=\"\\t\")";
		eAnalysis.setScript(script);
		Files.write(script, eAnalysis.getrScriptFile());
		return eAnalysis;

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
		script += MY_COUNTS_VAR + " <- getBamCounts(bed.frame = exons.hg19 , bam.files=BAMFILES, include.chr=TRUE, index.files=BAIFILES )\n";
		script += EXOME_COUNTS_DAFR + " <- as(" + MY_COUNTS_VAR + "[, colnames(" + MY_COUNTS_VAR + ")], 'data.frame')\n";
		script += EXOME_COUNTS_DAFR + "$chromosome <- gsub(as.charachter(" + EXOME_COUNTS_DAFR + "$space),pattern = 'chr', replacement = '')\n";
		script += "save(" + EXOME_COUNTS_DAFR + ",file=\"" + getCountFile() + "\")\n";
		return script;
	}

	private String getCountFile() {
		return outputDir + outputRoot + EXOME_COUNTS_DAFR + ".Rda";
	}

	private String loadCountFileScript() {
		return "load(\"" + getCountFile() + "\")";
	}

	private String addBaseLoadScript(String script) {
		script += "library(ExomeDepth)\n";
		script += "data(exons.hg19)\n";
		script += "data(Conrad.hg19)\n";
		return script;
	}

	private String[] generateBamBaiRVectors() {
		String[] bamBaiV = new String[2];
		bamBaiV[0] = Rscript.generateRVector(allReferenceBAMFiles);
		bamBaiV[1] = Rscript.generateRVector(allReferenceBAIFiles);
		return bamBaiV;
	}

	private static class ExomeDepthAnalysisProducer implements Producer<ExomeDepthAnalysis> {
		private ExomeDepth exomeDepth;
		private int index;
		private Logger log;

		public ExomeDepthAnalysisProducer(ExomeDepth exomeDepth, Logger log) {
			super();
			this.exomeDepth = exomeDepth;
			this.index = 0;
			this.log = log;
		}

		@Override
		public boolean hasNext() {

			// TODO Auto-generated method stub
			return index < exomeDepth.getAnalysisBamFiles().length;
		}

		@Override
		public Callable<ExomeDepthAnalysis> next() {
			final ExomeDepthAnalysis eAnalysis = new ExomeDepthAnalysis(exomeDepth.getAnalysisBamFiles()[index], exomeDepth.getOutputDir(), exomeDepth.getOutputRoot());
			exomeDepth.generateCallingScript(eAnalysis);
			Callable<ExomeDepthAnalysis> callable = new Callable<ExomeDepth.ExomeDepthAnalysis>() {

				@Override
				public ExomeDepthAnalysis call() throws Exception {
					log.reportTimeInfo("Running ExomeDepth on " + exomeDepth.getAnalysisBamFiles()[index]);
					CmdLine.runCommandWithFileChecks(new String[] { "Rscript", eAnalysis.getrScriptFile() }, "", exomeDepth.getAllReferenceBAMFiles(), new String[] { eAnalysis.getExomeDepthOutput() }, true, false, false, log);
					return eAnalysis;
				}
			};
			index++;
			return callable;
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static class ExomeDepthAnalysis {
		private String inputBam;
		private String outputDir;
		private String outputRoot;
		private String exomeDepthOutput;
		private String script;
		private String rScriptFile;
		private boolean fail;
		private Logger log;

		private ExomeDepthAnalysis(String inputBam, String outputDir, String outputRoot) {
			super();
			this.inputBam = inputBam;
			this.outputDir = outputDir;
			this.outputRoot = outputRoot;
			this.exomeDepthOutput = outputDir + outputRoot + ext.rootOf(inputBam) + ".exCNV";
			this.rScriptFile = outputDir + outputRoot + ext.rootOf(inputBam) + ".RScript";
		}

		public String getrScriptFile() {
			return rScriptFile;
		}

		public String getScript() {
			return script;
		}

		public void setScript(String script) {
			this.script = script;
		}

		public String getInputBam() {
			return inputBam;
		}

		public String getExomeDepthOutput() {
			return exomeDepthOutput;
		}

	}

	public static void runExomeDepth(String bams, String outputDir, String outputRoot, int numBatches, int numthreads, int wallTimeInHours, int memoryInMb, Logger log) {
		String[] allReferenceBamFiles = Files.isDirectory(bams) ? Files.listFullPaths(bams, BamOps.BAM_EXT, false) : HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, true);
		outputDir = outputDir == null ? ext.parseDirectoryOfFile(bams) : outputDir;

		log.reportTimeInfo("found " + allReferenceBamFiles.length + " bam files in " + bams);
		if (numBatches > 0) {
			log.reportTimeInfo("number of batches set to " + numBatches + ", preparing for batched run...");
			ArrayList<String[]> batches = Array.splitUpArray(allReferenceBamFiles, numBatches, log);
			for (int i = 0; i < batches.size(); i++) {
				ExomeDepth exomeDepth = new ExomeDepth(allReferenceBamFiles, batches.get(i), outputDir, outputRoot, log);
				if (!Files.exists(exomeDepth.getCountFile())) {
					log.reportTimeWarning("Did not find " + exomeDepth.getCountFile() + ", generating it now (takes a long time)");
					exomeDepth.generateCountFile();
				} else {
					log.reportTimeWarning("Using existing count file " + exomeDepth.getCountFile());
				}
				String tmpBams = outputDir + outputRoot + "ExomeDepth_" + i + ".bams";
				Files.writeList(batches.get(i), tmpBams);

				String qsub = outputDir + outputRoot + "ExomeDepth_" + i + ".pbs";
				String command = "";
				command += Array.toStr(PSF.Java.buildJavaCPXMX(PSF.Java.GENVISIS, "seq.analysis.ExomeDepth", memoryInMb), " ");
				command += " bams=" + tmpBams;
				command += " numBatches=0";
				command += " " + PSF.Ext.OUTPUT_DIR_COMMAND + outputDir;
				command += " root=" + outputRoot;
				command += " " + PSF.Ext.NUM_THREADS_COMMAND + numthreads;
				Files.qsub(qsub, command, memoryInMb, wallTimeInHours, numthreads);
				String script = exomeDepth.generateBamBaiScript();
				System.out.println(script);
			}
		} else {
			ExomeDepth exomeDepth = new ExomeDepth(allReferenceBamFiles, allReferenceBamFiles, outputDir, outputRoot, log);
			if (!Files.exists(exomeDepth.getCountFile())) {
				log.reportTimeError("Did not find " + exomeDepth.getCountFile() + ", please start the analysis in batch mode first");
				return;
			}
			ExomeDepthAnalysisProducer producer = new ExomeDepthAnalysisProducer(exomeDepth, log);
			WorkerTrain<ExomeDepthAnalysis> train = new WorkerTrain<ExomeDepth.ExomeDepthAnalysis>(producer, numthreads, numthreads, log);
			while (train.hasNext()) {
				train.next();
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String bams = "bams/";
		String outputDir = null;
		String outputRoot = "ExomeDepth";
		int memoryInMb = 62000;
		int wallTimeInHours = 47;
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
		usage += PSF.Ext.getMemoryMbCommand(6, memoryInMb);
		usage += PSF.Ext.getWallTimeCommand(7, wallTimeInHours);

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
			} else if (args[i].startsWith(PSF.Ext.WALLTIME_HRS)) {
				wallTimeInHours = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.MEMORY_MB)) {
				memoryInMb = ext.parseIntArg(args[i]);
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
			runExomeDepth(bams, outputDir, outputRoot, numBatches, numthreads, wallTimeInHours, memoryInMb, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
