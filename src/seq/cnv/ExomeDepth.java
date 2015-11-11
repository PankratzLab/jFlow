package seq.cnv;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import cnv.var.CNVariant;
import cnv.var.LocusSet;
import cnv.var.LocusSet.TO_STRING_TYPE;
import seq.cnv.CNVExtraInfo.EXTRA_INFO_TYPE;
import seq.manage.BamOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import stats.Rscript;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.Positions;
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
	private static final String[] RESULT_PARSE = new String[] { "type", "nexons", "start", "end", "chromosome", "BF", "reads.expected", "reads.observed", "reads.ratio" };

	private String[] allReferenceBAMFiles, allReferenceBAIFiles, allSampleNames, analysisBamFiles;
	private String outputDir;
	private String outputRoot;
	private boolean fail;
	private Hashtable<String, HashSet<String>> sampleSpecificExclude;// for excluding relatives from the reference set
	private HashSet<String> globalExclude; // never use these for a reference

	private Logger log;

	public ExomeDepth(String[] allReferenceBamFiles, String[] analysisBamFiles, String outputDir, String outputRoot, Logger log) {
		super();
		this.allReferenceBAMFiles = allReferenceBamFiles;
		this.log = log;
		this.sampleSpecificExclude = new Hashtable<String, HashSet<String>>();
		this.globalExclude = new HashSet<String>();
		this.allSampleNames = BamOps.getSampleNames(allReferenceBamFiles);
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

	public String[] getAllSampleNames() {
		return allSampleNames;
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

	public void parseVpop(VcfPopulation vpop) {
		for (int i = 0; i < allSampleNames.length; i++) {
			if (vpop.getSuperPop().get(VcfPopulation.EXCLUDE).contains(allSampleNames[i])) {
				globalExclude.add(allSampleNames[i]);
			}
			String[] sampSpecificExclude = vpop.getPopulationForInd(allSampleNames[i], RETRIEVE_TYPE.SUB);
			for (int j = 0; j < sampSpecificExclude.length; j++) {
				Set<String> curSet = vpop.getSubPop().get(sampSpecificExclude[j]);
				for (String samp : curSet) {
					if (!samp.equals(allSampleNames[i])) {
						sampleSpecificExclude.get(allSampleNames[i]).add(samp);
					}
				}
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
				sampleSpecificExclude.put(allSampleNames[i], new HashSet<String>());
			}
		}
		return verify;
	}

	public Hashtable<String, HashSet<String>> getExcludeFromReference() {
		return sampleSpecificExclude;
	}

	/**
	 * @param eAnalysis
	 *            the script property is modified
	 * @return
	 */
	private ExomeDepthAnalysis generateCallingScript(final ExomeDepthAnalysis eAnalysis) {
		String script = addBaseLoadScript("");
		script += loadCountFileScript() + "\n";
		script += MY_TEST_VAR + " <- " + EXOME_COUNTS_DAFR + "$" + Rscript.makeRSafe(ext.removeDirectoryInfo(eAnalysis.getInputBam())) + "\n";
		ArrayList<String> tmpRef = new ArrayList<String>();
		for (int i = 0; i < allReferenceBAMFiles.length; i++) {
			if (!allReferenceBAMFiles[i].equals(eAnalysis.getInputBam()) && !eAnalysis.getExcludeFromRef().contains(allSampleNames[i]) && !globalExclude.contains(allSampleNames[i])) {
				tmpRef.add(Rscript.makeRSafe(ext.removeDirectoryInfo(allReferenceBAMFiles[i])));
			}
		}

		script += MY_REF_SAMPLES_VAR + " <- " + Rscript.generateRVector(Array.toStringArray(tmpRef), true) + "\n";
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
		script += "write.table(" + MY_ALL_EXONS_VAR + "@CNV.calls, " + "\"" + eAnalysis.getExomeDepthOutput() + "\", sep=\"\\t\", row.names = FALSE , quote=FALSE)\n";
		script += MY_ALL_EXONS_VAR + "<- AnnotateExtra(x = " + MY_ALL_EXONS_VAR + ",";
		script += "reference.annotation = Conrad.hg19.common.CNVs, min.overlap = 0.5, column.name = 'Conrad.hg19')\n";
		script += "write.table(" + MY_ALL_EXONS_VAR + "@CNV.calls, " + "\"" + eAnalysis.getAnnoExomeDepthOutput() + "\", sep=\"\\t\", row.names = FALSE , quote=FALSE)\n";
		script += "save(" + MY_ALL_EXONS_VAR + ",file=\"" + eAnalysis.getRDafrExomeDepthOutput() + "\")\n";
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

	public boolean generateCountFile() {
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
		script += EXOME_COUNTS_DAFR + "$chromosome <- gsub(as.character(" + EXOME_COUNTS_DAFR + "$space),pattern = 'chr', replacement = '')\n";
		script += EXOME_COUNTS_DAFR + "$space <- gsub(as.character(" + EXOME_COUNTS_DAFR + "$space),pattern = 'chr', replacement = '')\n";
		script += "save(" + EXOME_COUNTS_DAFR + ",file=\"" + getCountFile() + "\")\n";
		return script;
	}

	public String getCountFile() {
		return outputDir + outputRoot + EXOME_COUNTS_DAFR + ".Rda";
	}

	private String loadCountFileScript() {
		return "load(\"" + getCountFile() + "\")";
	}

	private static String addBaseLoadScript(String script) {
		script += "library(ExomeDepth)\n";
		script += "data(exons.hg19)\n";
		script += "data(Conrad.hg19)\n";
		return script;
	}

	private String[] generateBamBaiRVectors() {
		String[] bamBaiV = new String[2];
		bamBaiV[0] = Rscript.generateRVector(allReferenceBAMFiles, true);
		bamBaiV[1] = Rscript.generateRVector(allReferenceBAIFiles, true);
		return bamBaiV;
	}

	static class ExomeDepthAnalysisProducer implements Producer<ExomeDepthAnalysis> {
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
			return index < exomeDepth.getAnalysisBamFiles().length;
		}

		@Override
		public Callable<ExomeDepthAnalysis> next() {
			final ExomeDepthAnalysis eAnalysis = new ExomeDepthAnalysis(exomeDepth.getAnalysisBamFiles()[index], exomeDepth.getOutputDir(), exomeDepth.getOutputRoot(), log);
			eAnalysis.setExcludeFromRef(exomeDepth.getExcludeFromReference().get(exomeDepth.getAllSampleNames()[index]));
			exomeDepth.generateCallingScript(eAnalysis);
			Callable<ExomeDepthAnalysis> callable = new Callable<ExomeDepth.ExomeDepthAnalysis>() {

				@Override
				public ExomeDepthAnalysis call() throws Exception {
					log.reportTimeInfo("Running ExomeDepth on " + eAnalysis.getSampleName());
					if (eAnalysis.getExcludeFromRef().size() > 0) {
						log.reportTimeInfo("Excluding sample specific " + eAnalysis.getExcludeFromRef().toString() + " samples from reference set for " + eAnalysis.getSampleName());
					}
					CmdLine.runCommandWithFileChecks(new String[] { "Rscript", eAnalysis.getrScriptFile() }, "", exomeDepth.getAllReferenceBAMFiles(), new String[] { eAnalysis.getExomeDepthOutput(), eAnalysis.getAnnoExomeDepthOutput() }, true, false, false, log);
					eAnalysis.plotCNVs(0.5);
					return eAnalysis;
				}
			};
			index++;
			return callable;
		}

		@Override
		public void remove() {
		}

		@Override
		public void shutdown() {
		}
	}

	static class ExomeDepthAnalysis {
		private String sampleName;
		private String inputBam;
		private String exomeDepthOutput;
		private String rDafrexomeDepthOutput;
		private String exomeDepthPDFOutput;
		private HashSet<String> excludeFromRef;
		private String script;
		private String rScriptFile;
		//private boolean fail;
		private Logger log;

		private ExomeDepthAnalysis(String inputBam, String outputDir, String outputRoot, Logger log) {
			super();
			this.inputBam = inputBam;
			this.exomeDepthOutput = outputDir + outputRoot + ext.rootOf(inputBam) + ".exCNV";
			this.rScriptFile = outputDir + outputRoot + ext.rootOf(inputBam) + ".RScript";
			this.rDafrexomeDepthOutput = outputDir + outputRoot + ext.rootOf(inputBam) + ".dafr";
			this.exomeDepthPDFOutput = outputDir + outputRoot + ext.rootOf(inputBam) + "cnvs.pdf";
			this.sampleName = BamOps.getSampleName(inputBam);
			this.excludeFromRef = new HashSet<String>();
			this.log = log;
		}

		public HashSet<String> getExcludeFromRef() {
			return excludeFromRef;
		}

		public String getSampleName() {
			return sampleName;
		}

		public String getRDafrExomeDepthOutput() {
			return rDafrexomeDepthOutput;
		}

		public String getrScriptFile() {
			return rScriptFile;
		}

		public void setExcludeFromRef(HashSet<String> excludeFromRef) {
			this.excludeFromRef = excludeFromRef;
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

		public String getAnnoExomeDepthOutput() {
			return ext.addToRoot(exomeDepthOutput, ".anno");
		}

		public boolean plotCNVs(double bufferPercent) {
			String script = "";
			script += addBaseLoadScript(script);
			script += "pdf(file = \"" + exomeDepthPDFOutput + "\")\n";
			script += getCNVPlotScript(bufferPercent);
			script += "dev.off()";
			String scriptFile = exomeDepthPDFOutput + ".Rscript";
			CmdLine.prepareBatchForCommandLine(new String[] { script }, scriptFile, true, log);
			boolean created = CmdLine.runCommandWithFileChecks(new String[] { "Rscript", scriptFile }, "", new String[] { rDafrexomeDepthOutput }, new String[] { exomeDepthPDFOutput }, true, true, false, log);
			return created;

		}

		public String getCNVPlotScript(double bufferPercent) {
			SeqCNVariant[] cnvs = processCNVs();
			ArrayList<SeqCNVariant> tmp = CNVariant.sortByQuality(cnvs, 1);
			cnvs = tmp.toArray(new SeqCNVariant[tmp.size()]);
			String script = "";
			script += "load(\"" + rDafrexomeDepthOutput + "\")\n";
			for (int i = 0; i < cnvs.length; i++) {
				script = getPlotFor(bufferPercent, cnvs[i], script);
			}
			return script;
		}

		private static String getPlotFor(double bufferPercent, SeqCNVariant cnv, String script) {
			script += "plot(all.exons , sequence = \"" + cnv.getChr() + "\",";
			int buffer = (int) (bufferPercent * cnv.getSize());
			String[] curBoundary = new String[] { cnv.getStart() + " - " + buffer, cnv.getStop() + " + " + buffer };
			script += "xlim = " + Rscript.generateRVector(curBoundary, false) + ",";
			script += "main = '" + cnv.getIndividualID();
			// + sampleName + "_" + i;
			for (int i = 0; i < cnv.getcExtraInfos().length; i++) {
				script += "\\n" + cnv.getcExtraInfos()[i].getsExtra() + "=" + cnv.getcExtraInfos()[i].getdExtra();
			}
			script += "\\nscore=" + cnv.getScore() + "' , cex.lab = 0.8,";
			script += "with.gene = FALSE)\n";
			return script;
		}

		public SeqCNVariant[] processCNVs() {

			ArrayList<SeqCNVariant> cnvs = new ArrayList<SeqCNVariant>();
			if (!Files.exists(getExomeDepthOutput())) {
				log.reportFileNotFound(getExomeDepthOutput());
			} else {
				String[] header = Files.getHeaderOfFile(getExomeDepthOutput(), log);
				int[] indices = ext.indexFactors(RESULT_PARSE, header, true, false);
				if (Array.countIf(indices, -1) > 0) {
					log.reportTimeError("Did not find complete header " + Array.toStr(RESULT_PARSE) + " in " + getExomeDepthOutput());
				} else {
					try {
						BufferedReader reader = Files.getAppropriateReader(getExomeDepthOutput());
						reader.readLine();
						while (reader.ready()) {
							String[] line = reader.readLine().trim().split("\t");
							try {
								int cn = -99;
								if (line[indices[0]].equals("deletion")) {
									cn = 1;
								} else if (line[indices[0]].equals("duplication")) {
									cn = 3;
								} else {
									log.reportTimeError("Invalid copy number type on line " + Array.toStr(line));
									return null;
								}
								int nexons = Integer.parseInt(line[indices[1]]);
								int start = Integer.parseInt(line[indices[2]]);
								int stop = Integer.parseInt(line[indices[3]]);
								byte chr = Positions.chromosomeNumber(line[indices[4]]);
								double score = Double.parseDouble(line[indices[5]]);
								ExomeDepthEI[] eis = new ExomeDepthEI[3];
								eis[0] = new ExomeDepthEI(EXTRA_INFO_TYPE.EXOME_DEPTH, RESULT_PARSE[6], Double.parseDouble(line[indices[6]]));
								eis[1] = new ExomeDepthEI(EXTRA_INFO_TYPE.EXOME_DEPTH, RESULT_PARSE[7], Double.parseDouble(line[indices[7]]));
								eis[2] = new ExomeDepthEI(EXTRA_INFO_TYPE.EXOME_DEPTH, RESULT_PARSE[8], Double.parseDouble(line[indices[8]]));

								SeqCNVariant cnVariant = new SeqCNVariant(sampleName, sampleName, chr, start, stop, cn, score, nexons, 99, eis);
								cnvs.add(cnVariant);
							} catch (NumberFormatException nfe) {
								log.reportTimeError("Invalid number on line " + Array.toStr(line));
								return null;
							}
						}
						reader.close();
					} catch (FileNotFoundException e) {
						e.printStackTrace();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}
			return cnvs.toArray(new SeqCNVariant[cnvs.size()]);
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
				command += Array.toStr(PSF.Java.buildJavaCPXMX(PSF.Java.GENVISIS, "seq.cnv.ExomeDepth", memoryInMb), " ");
				command += " bams=" + tmpBams;
				command += " numBatches=0";
				command += " " + PSF.Ext.OUTPUT_DIR_COMMAND + outputDir;
				command += " root=" + outputRoot;
				command += " " + PSF.Ext.NUM_THREADS_COMMAND + numthreads;
				Files.qsub(qsub, command, memoryInMb, wallTimeInHours, numthreads);
			}
		} else {
			ExomeDepth exomeDepth = new ExomeDepth(allReferenceBamFiles, allReferenceBamFiles, outputDir, outputRoot, log);
			callCNVs(exomeDepth, outputDir, outputRoot, numthreads, log);
		}
	}

	public static ExomeDepthAnalysis[] callCNVs(ExomeDepth exomeDepth, String outputDir, String outputRoot, int numthreads, Logger log) {
		if (!Files.exists(exomeDepth.getCountFile())) {
			log.reportTimeError("Did not find " + exomeDepth.getCountFile() + ", please start the analysis in batch mode first");
			return null;
		}
		ExomeDepthAnalysis[] eAnalysis = new ExomeDepthAnalysis[exomeDepth.getAnalysisBamFiles().length];

		ExomeDepthAnalysisProducer producer = new ExomeDepthAnalysisProducer(exomeDepth, log);
		WorkerTrain<ExomeDepthAnalysis> train = new WorkerTrain<ExomeDepth.ExomeDepthAnalysis>(producer, numthreads, numthreads, log);
		int index = 0;
		while (train.hasNext()) {
			eAnalysis[index] = train.next();
			index++;
		}
		ArrayList<SeqCNVariant> allTmp = new ArrayList<SeqCNVariant>();
		for (int i = 0; i < eAnalysis.length; i++) {
			SeqCNVariant[] tmp = eAnalysis[i].processCNVs();
			for (int j = 0; j < tmp.length; j++) {
				allTmp.add(tmp[j]);
			}
		}
		allTmp = CNVariant.sortByQuality(allTmp.toArray(new SeqCNVariant[allTmp.size()]), 1);
		LocusSet<SeqCNVariant> set = new LocusSet<SeqCNVariant>(allTmp.toArray(new SeqCNVariant[allTmp.size()]), true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;
		};
		set.writeRegions(outputDir + outputRoot + ".all.cnvs", TO_STRING_TYPE.REGULAR, true, log);
		return eAnalysis;
	}

	private static class ExomeDepthEI extends CNVExtraInfo {
		private static final long serialVersionUID = 1L;

		public ExomeDepthEI(EXTRA_INFO_TYPE type, String title, double extra) {
			super(EXTRA_INFO_TYPE.EXOME_DEPTH);
			this.dExtra = extra;
			this.sExtra = title;
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
