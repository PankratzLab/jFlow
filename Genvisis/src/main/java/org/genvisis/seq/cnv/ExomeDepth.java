package org.genvisis.seq.cnv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.LocusSet.TO_STRING_TYPE;
import org.genvisis.seq.cnv.CNVExtraInfo.EXTRA_INFO_TYPE;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.genvisis.seq.manage.VCFTumorNormalOps;
import org.genvisis.seq.manage.VCFTumorNormalOps.TNSample;
import org.genvisis.stats.Rscript;

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
	private static final String[] RESULT_PARSE = new String[] {	"type", "nexons", "start", "end",
																															"chromosome", "BF", "reads.expected",
																															"reads.observed", "reads.ratio"};

	private final String[] allReferenceBAMFiles;
	private String[] allReferenceBAIFiles;
	private final String[] allSampleNames;
	private final String[] analysisBamFiles;
	private final String outputDir;
	private final String outputRoot;
	private boolean fail;
	private final Hashtable<String, HashSet<String>> sampleSpecificExclude;// for excluding relatives
																																					// from the reference set
	private final HashSet<String> globalExclude; // never use these for a reference
	private final String rLoc;
	private final Logger log;

	public ExomeDepth(String[] allReferenceBamFiles, String[] analysisBamFiles, String outputDir,
										String outputRoot, String rLoc, Logger log) {
		super();
		allReferenceBAMFiles = allReferenceBamFiles;
		this.log = log;
		sampleSpecificExclude = new Hashtable<String, HashSet<String>>();
		globalExclude = new HashSet<String>();
		allSampleNames = BamOps.getSampleNames(allReferenceBamFiles, log);
		fail = gatherBai();

		this.analysisBamFiles = analysisBamFiles;
		this.outputDir = outputDir;
		new File(outputDir).mkdirs();
		this.outputRoot = outputRoot;
		if (!fail) {
			fail = Array.countIf(	ext.indexLargeFactors(analysisBamFiles, allReferenceBamFiles, true, log,
																									true, false),
														-1) > 0;
			if (fail) {
				log.reportError("Could not detect all analysis .bam files in the complete reference set");
			}
		}
		this.rLoc = rLoc;
	}

	public String[] getAnalysisBamFiles() {
		return analysisBamFiles;
	}

	public String getRscriptCall() {
		return rLoc == null ? "Rscript" : rLoc + "Rscript";
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
		if (vpop.getType() == POPULATION_TYPE.TUMOR_NORMAL) {
			log.reportTimeInfo("Setting up exome depth for somatic cnv calling");
			TNSample[] tnsaSamples = VCFTumorNormalOps.matchSamples(allReferenceBAMFiles, vpop, log);
			for (int i = 0; i < tnsaSamples.length; i++) {
				for (int j = 0; j < tnsaSamples.length; j++) {
					if (i != j) {
						sampleSpecificExclude	.get(tnsaSamples[i].getTumorSample())
																	.add(tnsaSamples[j].getTumorSample());
						sampleSpecificExclude	.get(tnsaSamples[i].getTumorSample())
																	.add(tnsaSamples[j].getNormalSample());
					}
				}
			}

		} else {
			for (int i = 0; i < allSampleNames.length; i++) {
				if (vpop.getSuperPop().get(VcfPopulation.EXCLUDE).contains(allSampleNames[i])) {
					globalExclude.add(allSampleNames[i]);
				}
				String[] sampSpecificExclude = vpop.getPopulationForInd(allSampleNames[i],
																																RETRIEVE_TYPE.SUB);
				for (String element : sampSpecificExclude) {
					Set<String> curSet = vpop.getSubPop().get(element);
					for (String samp : curSet) {
						if (!samp.equals(allSampleNames[i])) {
							sampleSpecificExclude.get(allSampleNames[i]).add(samp);
						}
					}
				}
			}
		}

	}

	private boolean gatherBai() {

		boolean verify = true;
		allReferenceBAIFiles = new String[allReferenceBAMFiles.length];
		for (int i = 0; i < allReferenceBAMFiles.length; i++) {
			String bai = ext.rootOf(allReferenceBAMFiles[i], false) + BamOps.BAI_EXT;
			if (!Files.exists(allReferenceBAMFiles[i]) || !Files.exists(bai)) {
				log.reportError("Could not find "	+ allReferenceBAMFiles[i]
														+ " with corresponding .bai file" + bai);
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
	 * @param eAnalysis the script property is modified
	 * @return
	 */
	private ExomeDepthAnalysis generateCallingScript(final ExomeDepthAnalysis eAnalysis) {
		String script = addBaseLoadScript("");
		script += loadCountFileScript() + "\n";
		script += MY_TEST_VAR	+ " <- " + EXOME_COUNTS_DAFR + "$"
							+ Rscript.makeRSafe(ext.removeDirectoryInfo(eAnalysis.getInputBam())) + "\n";
		ArrayList<String> tmpRef = new ArrayList<String>();
		for (int i = 0; i < allReferenceBAMFiles.length; i++) {
			if (!allReferenceBAMFiles[i].equals(eAnalysis.getInputBam())
						&& !eAnalysis.getExcludeFromRef().contains(allSampleNames[i])
					&& !globalExclude.contains(allSampleNames[i])) {
				tmpRef.add(Rscript.makeRSafe(ext.removeDirectoryInfo(allReferenceBAMFiles[i])));
			}
		}

		script += MY_REF_SAMPLES_VAR	+ " <- "
							+ Rscript.generateRVector(Array.toStringArray(tmpRef), true) + "\n";
		script += MY_REF_SET_VAR	+ " <- " + "as.matrix(" + EXOME_COUNTS_DAFR + "[, "
							+ MY_REF_SAMPLES_VAR + "])\n";
		if (tmpRef.size() == 0) {
			throw new IllegalArgumentException("0 size reference set, cannot call cnvs with exomeDepth");
		}

		if (tmpRef.size() > 1) {
			script += MY_CHOICE_VAR	+ " <- " + "select.reference.set (test.counts = " + MY_TEST_VAR
								+ ",reference.counts = " + MY_REF_SET_VAR;
			script += ",bin.length = ("	+ EXOME_COUNTS_DAFR + "$end - " + EXOME_COUNTS_DAFR
								+ "$start)/1000,n.bins.reduced = 10000)\n";
			script += "print(" + MY_CHOICE_VAR + "[[1]])\n";

			script += MY_MATRIX_VAR	+ " <- as.matrix( " + EXOME_COUNTS_DAFR + "[," + MY_CHOICE_VAR
								+ "$reference.choice, drop = FALSE])\n";
			script += MY_REF_SET_SELECTED_VAR + "<- apply(X =" + MY_MATRIX_VAR + ",MAR = 1,FUN = sum)\n";
		} else {// causes R error otherwise, and no need to select ref
			script += MY_REF_SET_SELECTED_VAR	+ " <- " + EXOME_COUNTS_DAFR + "$"
								+ Rscript.makeRSafe(ext.removeDirectoryInfo(tmpRef.get(0))) + "\n";
		}
		script += MY_ALL_EXONS_VAR	+ "<- new ('ExomeDepth', test = " + MY_TEST_VAR + " , reference = "
							+ MY_REF_SET_SELECTED_VAR + ",";
		script += "formula = 'cbind(test,reference) ~ 1')\n";

		script += MY_ALL_EXONS_VAR + "<- CallCNVs(x = " + MY_ALL_EXONS_VAR + " ,";
		script += "chromosome = " + EXOME_COUNTS_DAFR + "$space,";
		script += "start = " + EXOME_COUNTS_DAFR + "$start,";
		script += "end = " + EXOME_COUNTS_DAFR + "$end,";
		script += "name = " + EXOME_COUNTS_DAFR + "$names)\n";
		script += "write.table("	+ MY_ALL_EXONS_VAR + "@CNV.calls, " + "\""
							+ eAnalysis.getExomeDepthOutput()
							+ "\", sep=\"\\t\", row.names = FALSE , quote=FALSE)\n";
		script += MY_ALL_EXONS_VAR + "<- AnnotateExtra(x = " + MY_ALL_EXONS_VAR + ",";
		script +=
						"reference.annotation = Conrad.hg19.common.CNVs, min.overlap = 0.5, column.name = 'Conrad.hg19')\n";
		script += "write.table("	+ MY_ALL_EXONS_VAR + "@CNV.calls, " + "\""
							+ eAnalysis.getAnnoExomeDepthOutput()
							+ "\", sep=\"\\t\", row.names = FALSE , quote=FALSE)\n";
		script += "save("	+ MY_ALL_EXONS_VAR + ",file=\"" + eAnalysis.getRDafrExomeDepthOutput()
							+ "\")\n";
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

		CmdLine.prepareBatchForCommandLine(new String[] {script}, scriptFile, true, log);
		created = CmdLine.runCommandWithFileChecks(	new String[] {getRscriptCall(), scriptFile}, "",
																								allReferenceBAMFiles, new String[] {getCountFile()},
																								true, false, false, log);
		return created;
	}

	private String generateCountsScripts() {
		String script = generateBamBaiScript();
		script += MY_COUNTS_VAR
							+ " <- getBamCounts(bed.frame = exons.hg19 , bam.files=BAMFILES, include.chr=TRUE, index.files=BAIFILES )\n";
		script += EXOME_COUNTS_DAFR	+ " <- as(" + MY_COUNTS_VAR + "[, colnames(" + MY_COUNTS_VAR
							+ ")], 'data.frame')\n";
		script += EXOME_COUNTS_DAFR	+ "$chromosome <- gsub(as.character(" + EXOME_COUNTS_DAFR
							+ "$space),pattern = 'chr', replacement = '')\n";
		script += EXOME_COUNTS_DAFR	+ "$space <- gsub(as.character(" + EXOME_COUNTS_DAFR
							+ "$space),pattern = 'chr', replacement = '')\n";
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

	static class ExomeDepthAnalysisProducer extends AbstractProducer<ExomeDepthAnalysis> {
		private final ExomeDepth exomeDepth;
		private int index;
		private final Logger log;

		public ExomeDepthAnalysisProducer(ExomeDepth exomeDepth, Logger log) {
			super();
			this.exomeDepth = exomeDepth;
			index = 0;
			this.log = log;
		}

		@Override
		public boolean hasNext() {
			return index < exomeDepth.getAnalysisBamFiles().length;
		}

		@Override
		public Callable<ExomeDepthAnalysis> next() {
			final ExomeDepthAnalysis eAnalysis = new ExomeDepthAnalysis(exomeDepth.getAnalysisBamFiles()[index],
																																	exomeDepth.getOutputDir(),
																																	exomeDepth.getOutputRoot(),
																																	exomeDepth.getRscriptCall(), log);
			eAnalysis.setExcludeFromRef(exomeDepth.getExcludeFromReference()
																						.get(exomeDepth.getAllSampleNames()[index]));
			exomeDepth.generateCallingScript(eAnalysis);
			Callable<ExomeDepthAnalysis> callable = new Callable<ExomeDepth.ExomeDepthAnalysis>() {

				@Override
				public ExomeDepthAnalysis call() throws Exception {
					log.reportTimeInfo("Running ExomeDepth on " + eAnalysis.getSampleName());
					if (eAnalysis.getExcludeFromRef().size() > 0) {
						log.reportTimeInfo("Excluding sample specific "
																	+ eAnalysis.getExcludeFromRef().toString()
																+ " samples from reference set for " + eAnalysis.getSampleName());
					}
					CmdLine.runCommandWithFileChecks(	new String[] {eAnalysis.getrScriptCall(),
																													eAnalysis.getrScriptFile()},
																						"", exomeDepth.getAllReferenceBAMFiles(),
																						new String[] {eAnalysis.getExomeDepthOutput(),
																													eAnalysis.getAnnoExomeDepthOutput()},
																						true, false, false, log);
					eAnalysis.plotCNVs(0.5);
					eAnalysis.dumpRawData();
					return eAnalysis;
				}
			};
			index++;
			return callable;
		}
	}

	static class ExomeDepthAnalysis {
		private final String sampleName;
		private final String inputBam;
		private final String exomeDepthOutput;
		private final String rDafrexomeDepthOutput;
		private final String exomeDepthPDFOutput;
		private final String exomeDepthRawDataOutput;
		private HashSet<String> excludeFromRef;
		private String script;
		private final String rScriptFile;
		private final String rScriptCall;
		// private boolean fail;
		private final Logger log;

		private ExomeDepthAnalysis(	String inputBam, String outputDir, String outputRoot,
																String rScriptCall, Logger log) {
			super();
			this.inputBam = inputBam;
			exomeDepthOutput = outputDir + outputRoot + ext.rootOf(inputBam) + ".exCNV";
			rScriptFile = outputDir + outputRoot + ext.rootOf(inputBam) + ".RScript";
			rDafrexomeDepthOutput = outputDir + outputRoot + ext.rootOf(inputBam) + ".dafr";
			exomeDepthPDFOutput = outputDir + outputRoot + ext.rootOf(inputBam) + "cnvs.pdf";
			exomeDepthRawDataOutput = outputDir + outputRoot + ext.rootOf(inputBam) + "rawData.txt";
			sampleName = BamOps.getSampleName(inputBam);
			excludeFromRef = new HashSet<String>();
			this.rScriptCall = rScriptCall;
			this.log = log;
		}

		public String getExomeDepthRawDataOutput() {
			return exomeDepthRawDataOutput;
		}

		public String getrScriptCall() {
			return rScriptCall;
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
			CmdLine.prepareBatchForCommandLine(new String[] {script}, scriptFile, true, log);
			boolean created = CmdLine.runCommandWithFileChecks(	new String[] {"Rscript", scriptFile}, "",
																													new String[] {rDafrexomeDepthOutput},
																													new String[] {exomeDepthPDFOutput}, true,
																													false, false, log);
			return created;
		}

		public String getCNVPlotScript(double bufferPercent) {
			SeqCNVariant[] cnvs = processCNVs();
			CNVariant.sortInPlaceByQuality(cnvs, true);
			String script = "";
			script += "load(\"" + rDafrexomeDepthOutput + "\")\n";
			for (SeqCNVariant cnv : cnvs) {
				script = getPlotFor(bufferPercent, cnv, script);
			}
			return script;
		}

		public boolean dumpRawData() {
			String script = "";
			script += addBaseLoadScript(script);
			script += getRawDataDumpScript();
			String scriptFile = exomeDepthRawDataOutput + ".rawData.Rscript";
			CmdLine.prepareBatchForCommandLine(new String[] {script}, scriptFile, true, log);
			boolean created = CmdLine.runCommandWithFileChecks(	new String[] {"Rscript", scriptFile}, "",
																													new String[] {rDafrexomeDepthOutput},
																													new String[] {exomeDepthRawDataOutput},
																													true, false, false, log);
			return created;
		}

		public String getRawDataDumpScript() {
			String script = "";
			script += "load(\"" + rDafrexomeDepthOutput + "\")\n";
			script += "expected = all.exons@expected\n";
			script += "test <- all.exons@test\n";
			script += "reference <- all.exons@reference\n";
			script += "freq = test/ (reference + test)\n";
			script += "ratio <-  freq/ expected\n";
			script += "anno <- all.exons@annotations\n";
			script +=
							"exomeObject = data.frame(anno$chromosome,anno$start,anno$end, test,expected,reference,ratio)\n";
			script += "write.table(exomeObject, \""	+ exomeDepthRawDataOutput
								+ "\",row.names = FALSE,quote=FALSE )\n";
			return script;

		}

		private static String getPlotFor(double bufferPercent, SeqCNVariant cnv, String script) {
			script += "plot(all.exons , sequence = \"" + cnv.getChr() + "\",";
			int buffer = (int) (bufferPercent * cnv.getSize());
			String[] curBoundary = new String[] {cnv.getStart()	+ " - " + buffer,
																						cnv.getStop() + " + " + buffer};
			script += "xlim = " + Rscript.generateRVector(curBoundary, false) + ",";
			script += "main = '" + cnv.getIndividualID();
			// + sampleName + "_" + i;
			for (int i = 0; i < cnv.getcExtraInfos().length; i++) {
				script += "\\n"	+ cnv.getcExtraInfos()[i].getsExtra() + "="
									+ cnv.getcExtraInfos()[i].getdExtra();
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
					log.reportError("Did not find complete header "	+ Array.toStr(RESULT_PARSE) + " in "
															+ getExomeDepthOutput());
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
									log.reportError("Invalid copy number type on line " + Array.toStr(line));
									return null;
								}
								int nexons = Integer.parseInt(line[indices[1]]);
								int start = Integer.parseInt(line[indices[2]]);
								int stop = Integer.parseInt(line[indices[3]]);
								byte chr = Positions.chromosomeNumber(line[indices[4]]);
								double score = Double.parseDouble(line[indices[5]]);
								ExomeDepthEI[] eis = new ExomeDepthEI[3];
								eis[0] = new ExomeDepthEI(EXTRA_INFO_TYPE.EXOME_DEPTH, RESULT_PARSE[6],
																					Double.parseDouble(line[indices[6]]));
								eis[1] = new ExomeDepthEI(EXTRA_INFO_TYPE.EXOME_DEPTH, RESULT_PARSE[7],
																					Double.parseDouble(line[indices[7]]));
								eis[2] = new ExomeDepthEI(EXTRA_INFO_TYPE.EXOME_DEPTH, RESULT_PARSE[8],
																					Double.parseDouble(line[indices[8]]));

								SeqCNVariant cnVariant = new SeqCNVariant(sampleName, sampleName, chr, start, stop,
																													cn, score, nexons, 99, eis);
								cnvs.add(cnVariant);
							} catch (NumberFormatException nfe) {
								log.reportError("Invalid number on line " + Array.toStr(line));
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

	public static void runExomeDepth(	String bams, String outputDir, String outputRoot, String Rloc,
																		int numBatches, int numthreads, int wallTimeInHours,
																		int memoryInMb, Logger log) {
		String[] allReferenceBamFiles = Files.isDirectory(bams)
																															? Files.listFullPaths(bams,
																																									BamOps.BAM_EXT,
																																									false)
																														: HashVec.loadFileToStringArray(bams,
																																														false,
																																														new int[] {0},
																																														true);
		outputDir = outputDir == null ? ext.parseDirectoryOfFile(bams) : outputDir;

		log.reportTimeInfo("found " + allReferenceBamFiles.length + " bam files in " + bams);
		if (numBatches > 0) {
			log.reportTimeInfo("number of batches set to "	+ numBatches
													+ ", preparing for batched run...");
			List<String[]> batches = Array.splitUpArray(allReferenceBamFiles, numBatches, log);
			for (int i = 0; i < batches.size(); i++) {
				ExomeDepth exomeDepth = new ExomeDepth(	allReferenceBamFiles, batches.get(i), outputDir,
																								outputRoot, Rloc, log);
				if (!Files.exists(exomeDepth.getCountFile())) {
					log.reportTimeWarning("Did not find "	+ exomeDepth.getCountFile()
																+ ", generating it now (takes a long time)");
					exomeDepth.generateCountFile();
				} else {
					log.reportTimeWarning("Using existing count file " + exomeDepth.getCountFile());
				}
				String tmpBams = outputDir + outputRoot + "ExomeDepth_" + i + ".bams";
				Files.writeArray(batches.get(i), tmpBams);

				String qsub = outputDir + outputRoot + "ExomeDepth_" + i + ".pbs";
				String command = "";
				command += Array.toStr(	PSF.Java.buildJavaCPXMX(PSF.Java.GENVISIS, "seq.cnv.ExomeDepth",
																												memoryInMb),
																" ");
				command += " bams=" + tmpBams;
				command += " numBatches=0";
				command += " " + PSF.Ext.OUTPUT_DIR_COMMAND + outputDir;
				command += " root=" + outputRoot;
				command += " " + PSF.Ext.NUM_THREADS_COMMAND + numthreads;
				if (Rloc != null) {
					command += " rDir=" + Rloc;
				}
				Files.qsub(qsub, command, memoryInMb, wallTimeInHours, numthreads);
			}
		} else {
			ExomeDepth exomeDepth = new ExomeDepth(	allReferenceBamFiles, allReferenceBamFiles, outputDir,
																							outputRoot, Rloc, log);
			callCNVs(exomeDepth, outputDir, outputRoot, numthreads, log);
		}
	}

	public static ExomeDepthAnalysis[] callCNVs(ExomeDepth exomeDepth, String outputDir,
																							String outputRoot, int numthreads, Logger log) {
		if (!Files.exists(exomeDepth.getCountFile())) {
			log.reportError("Did not find " + exomeDepth.getCountFile());
			log.reportError("This is most likely caused by a failure of running Rscript");

			return null;
		}
		ExomeDepthAnalysis[] eAnalysis = new ExomeDepthAnalysis[exomeDepth.getAnalysisBamFiles().length];

		ExomeDepthAnalysisProducer producer = new ExomeDepthAnalysisProducer(exomeDepth, log);
		WorkerTrain<ExomeDepthAnalysis> train = new WorkerTrain<ExomeDepth.ExomeDepthAnalysis>(	producer,
																																														numthreads,
																																														numthreads,
																																														log);
		int index = 0;
		while (train.hasNext()) {
			eAnalysis[index] = train.next();
			index++;
		}
		ArrayList<SeqCNVariant> allTmp = new ArrayList<SeqCNVariant>();
		for (ExomeDepthAnalysis eAnalysi : eAnalysis) {
			SeqCNVariant[] tmp = eAnalysi.processCNVs();
			for (SeqCNVariant element : tmp) {
				allTmp.add(element);
			}
		}
		SeqCNVariant[] cnvs = allTmp.toArray(new SeqCNVariant[allTmp.size()]);
		CNVariant.sortInPlaceByQuality(cnvs, true);
		LocusSet<SeqCNVariant> set = new LocusSet<SeqCNVariant>(cnvs, true, log) {

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
			dExtra = extra;
			sExtra = title;
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
		String Rloc = null;

		String usage = "\n" + "seq.analysis.ExomeDepth requires 0-1 arguments\n";
		usage += "   (1) full path to a directory of or file of bams (i.e. bams="	+ bams
							+ " (default))\n" + "";
		usage += "   (2) number of batches to run (i.e. numBatches=" + bams + " (default))\n" + "";
		usage += PSF.Ext.getOutputDirCommand(3, outputDir);
		usage += "   (4) output root command (i.e. root=" + outputRoot + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(5, numthreads);
		usage += PSF.Ext.getMemoryMbCommand(6, memoryInMb);
		usage += PSF.Ext.getWallTimeCommand(7, wallTimeInHours);
		usage += "   (8) alternative R location (i.e. rDir= (no default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("bams=")) {
				bams = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("rDir=")) {
				Rloc = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("numBatches=")) {
				numBatches = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.WALLTIME_HRS)) {
				wallTimeInHours = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.MEMORY_MB)) {
				memoryInMb = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
				outputDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("root")) {
				outputRoot = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(logfile);
			runExomeDepth(bams, outputDir, outputRoot, Rloc, numBatches, numthreads, wallTimeInHours,
										memoryInMb, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
