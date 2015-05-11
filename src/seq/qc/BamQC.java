package seq.qc;

import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Positions;
import common.ext;
import filesys.Segment;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import seq.analysis.SNPEFF;
import stats.Histogram;

public class BamQC {
	public static final String[] QC_HEADER = { "Input File", "numTotal", "numUnique", "numDuplicated", "numUnMapped", "numOnTarget", "NumSecondaryAlignments", "NumInvalidAlignments", "PercentDuplicated", "PercentUnMapped", "PercentOnTarget", "AverageInsertSize", "AverageOnTargetInsertSize", "MappingQualityFilter", "PhreadScoreFilter", "Total Base Pairs Targeted" };
	public static final String[] HIST_HEADER = { "Bin" };
	public static int NUM_GC_HISTOGRAMS = 2;
	public static int NO_QC_GC_HISTOGRAM = 0;
	public static int QC_GC_HISTOGRAM = 1;
	public static int NUM_INSERT_HISTOGRAMS = 2;
	public static int NO_QC_INSERT_HISTOGRAM = 0;
	public static int QC_INSERT_HISTOGRAM = 1;
	public static final String READ_DEPTH_AT = "Percent Coverage At";
	// private static final int DEFUALT_GC_SIG_FIGS = 2;
	// private static final double DEFUALT_MIN_GC = 0.0D;
	// private static final double DEFUALT_MAX_GC = 1.0D;
	// private static final double NORM_UNITS = 1000000.0D;
	// private static final int DEFUALT_INSERT_SIG_FIGS = 0;
	// private static final double DEFUALT_MIN_INSERT = 0.0D;
	// private static final double DEFUALT_MAX_INSERT = 2000.0D;
	private String inputSamOrBamFile;
	private String libraryReadDepthResultsFile;
	private int numDuplicated;
	private int numUnMapped;
	private int numOnTarget;
	private int numTotal;
	private int numUnique;
	private int totalTargetedBasePairs;
	private int numSecondaryAlignments;
	private int numInvalidAligments;
	private double percentDuplicated;
	private double percentUnMapped;
	private double percentOnTarget;
	private double[] percentCoverageAtDepth;
	private double[] averageInsertSize;
	private Histogram.DynamicHistogram[] gHistograms;
	private Histogram.DynamicHistogram[] insertHistograms;
	private FilterNGS filterNGS;
	public static final String COMMAND_DIR = "dir=";
	public static final String COMMAND_COMMON_EXT = "commonExt=";
	public static final String COMMAND_INPUT_FILE = "fileOfinputSamOrBams=";
	public static final String COMMAND_TARGET_LIBRARY = "target=";
	public static final String COMMAND_OUTPUT = "output=";
	public static final String COMMAND_SKIP_NUM_LINES = "skip=";
	public static final String COMMAND_NUM_THREADS = "numThreads=";
	public static final String COMMAND_LOG_FILE = "log=";
	public static final String COMMAND_MAPPING_SCORE = "mapQ=";
	public static final String COMMAND_PHREAD_SCORE = "phread=";
	public static final String COMMAND_READ_DEPTH = "readDepth=";
	public static final String COMMAND_BAITS_LIBRARY = "baits=";
	public static final String COMMAND_OUTPUT_DIRECTORY = "outDir=";
	public static final String COMMAND_BAITS_AS_TARGET = "-baits";
	public static final String COMMAND_NORM_READS = "normReads=";

	public BamQC(String inputSamOrBamFile, String outputDir, FilterNGS filterNGS) {
		this.inputSamOrBamFile = inputSamOrBamFile;
		this.libraryReadDepthResultsFile = Files.getSerializedFileName(outputDir, inputSamOrBamFile);
		this.numTotal = 0;
		this.numUnique = 0;
		this.numDuplicated = 0;
		this.numUnMapped = 0;
		this.numOnTarget = 0;
		this.numSecondaryAlignments = 0;
		this.numInvalidAligments = 0;
		this.averageInsertSize = new double[NUM_INSERT_HISTOGRAMS];
		this.totalTargetedBasePairs = 0;
		this.percentDuplicated = 0;
		this.percentUnMapped = 0;
		this.percentOnTarget = 0;
		this.percentCoverageAtDepth = new double[filterNGS.getReadDepthFilter().length];
		this.gHistograms = Histogram.DynamicHistogram.initHistograms(NUM_GC_HISTOGRAMS, 0, 1.0, 2);
		this.insertHistograms = Histogram.DynamicHistogram.initHistograms(NUM_INSERT_HISTOGRAMS, 0.0, 2000.0, 0);
		this.filterNGS = filterNGS;
	}

	public void addToGHistoGram(double gcContent, int which) {
		this.gHistograms[which].addDataPointToHistogram(gcContent);
	}

	public void addToInsertHistoGram(double insertSize, int which) {
		this.insertHistograms[which].addDataPointToHistogram(insertSize);
	}

	public Histogram.DynamicHistogram getGcHistogram(int which) {
		return this.gHistograms[which];
	}

	public Histogram.DynamicHistogram getInsertHistogram(int which) {
		return this.insertHistograms[which];
	}

	public void computePercents() {
		this.percentDuplicated = computePercent(this.numDuplicated, this.numTotal);
		this.percentUnMapped = computePercent(this.numUnMapped, this.numUnique);
		this.percentOnTarget = computePercent(this.numOnTarget, this.numUnique);
		this.averageInsertSize[NO_QC_INSERT_HISTOGRAM] /= Array.sum(getInsertHistogram(NO_QC_INSERT_HISTOGRAM).getCounts());
		this.averageInsertSize[QC_INSERT_HISTOGRAM] /= Array.sum(getInsertHistogram(QC_INSERT_HISTOGRAM).getCounts());
	}

	public void setPercentCoverageAtDepth(double[] percentCoverageAtDepth) {
		this.percentCoverageAtDepth = percentCoverageAtDepth;
	}

	public void addNumTotal() {
		this.numTotal += 1;
	}

	public double addInsertSize(SAMRecord samRecord, int which) {
		double insert = (0.0D / 0.0D);
		if (samRecord.getProperPairFlag()) {
			insert = Math.abs(samRecord.getInferredInsertSize());
			this.averageInsertSize[which] += insert;
			addToInsertHistoGram(insert, which);
		}
		return insert;
	}

	public String getLibraryReadDepthResultsFile() {
		return this.libraryReadDepthResultsFile;
	}

	public void setLibraryReadDepthResultsFile(String libraryReadDepthResultsFile) {
		this.libraryReadDepthResultsFile = libraryReadDepthResultsFile;
	}

	public void setTotalTargetedBasePairs(int totalTargetedBasePairs) {
		this.totalTargetedBasePairs = totalTargetedBasePairs;
	}

	public void addNumUnique(SAMRecord samRecord) {
		if (!samRecord.getDuplicateReadFlag()) {
			this.numUnique += 1;
		}
	}

	public void addNumUnMapped(SAMRecord samRecord) {
		if (samRecord.getMateUnmappedFlag()) {
			this.numUnMapped += 1;
		}
	}

	public void addNumDuplicate(SAMRecord samRecord) {
		if (samRecord.getDuplicateReadFlag()) {
			this.numDuplicated += 1;
		}
	}

	public void addNumOnTarget() {
		this.numOnTarget += 1;
	}

	public String getInputSamOrBamFile() {
		return this.inputSamOrBamFile;
	}

	public int getNumDuplicated() {
		return this.numDuplicated;
	}

	public int getNumUnMapped() {
		return this.numUnMapped;
	}

	public int getNumOnTarget() {
		return this.numOnTarget;
	}

	public int getNumTotalReads() {
		return this.numTotal;
	}

	public double getPercentDuplicated() {
		return this.percentDuplicated;
	}

	public double getPercentUnMapped() {
		return this.percentUnMapped;
	}

	public double getPercentOnTarget() {
		return this.percentOnTarget;
	}

	public int getNumSecondaryAlignments() {
		return this.numSecondaryAlignments;
	}

	public int getNumInvalidAligments() {
		return this.numInvalidAligments;
	}

	public boolean secondaryAlignment(SAMRecord samRecord) {
		if (samRecord.isSecondaryOrSupplementary()) {
			this.numSecondaryAlignments += 1;
			return true;
		}
		return false;
	}

	public boolean invalidAlignment(SAMRecord samRecord) {
		if (samRecord.isValid() != null) {
			this.numInvalidAligments += 1;
			return true;
		}
		return false;
	}

	public String getSummary() {
		String summary = "";
		summary = summary + this.inputSamOrBamFile + "\t" + this.numTotal + "\t" + this.numUnique + "\t" + this.numDuplicated + "\t" + this.numUnMapped + "\t" + this.numOnTarget + "\t" + this.numSecondaryAlignments + "\t" + this.numInvalidAligments + "\t" + this.percentDuplicated + "\t" + this.percentUnMapped + "\t" + this.percentOnTarget + "\t" + this.averageInsertSize[NO_QC_INSERT_HISTOGRAM] + "\t" + this.averageInsertSize[QC_INSERT_HISTOGRAM] + "\t" + this.filterNGS.getMappingQualityFilter() + "\t" + this.filterNGS.getPhreadScoreFilter() + "\t" + this.totalTargetedBasePairs + "\t" + Array.toStr(this.percentCoverageAtDepth);
		return summary;
	}

	public static BamQC[] qcBams(String dir, String outputDir, String commonExt, String fileOfinputSamOrBams, String targetLibraryFile, String baitLibraryFile, int skipNumLines, FilterNGS filterNGS, int numThreads, String output, String snpEffLocation, boolean baitsAsTarget, double normalizeDepthsTo, Logger log) {
		long time = System.currentTimeMillis();
		log.report(ext.getTime() + " Info - beginning qc");
		BamQC[] bamQCs = null;
		if ((dir == null) && (fileOfinputSamOrBams == null)) {
			log.reportError("Error - must provide a directory to search or a file listing files to qc");
		} else {
			String[] inputbams = determineInputFiles(dir, commonExt, fileOfinputSamOrBams, log);
			if ((inputbams == null) || (inputbams.length < 1)) {
				log.reportError("Error - did not find any files to qc");
			} else {
				LibraryNGS libraryNGS = null;
				LibraryNGS.BaitsLibrary lBaitsLibrary = null;
				if ((baitsAsTarget) && ((baitLibraryFile == null) || (!Files.exists(baitLibraryFile)))) {
					log.reportError("Error - baits a library was flagged, but could not find the baits file" + (baitLibraryFile == null ? "" : baitLibraryFile));
					return null;
				}
				if ((targetLibraryFile == null) && (!baitsAsTarget)) {
					log.report("Warning - a target library file was not provided, on target percentages will not be computed");
				} else {
					if (baitLibraryFile != null) {
						lBaitsLibrary = LibraryNGS.BaitsLibrary.loadBaitLibrary(baitLibraryFile, log);
					}
					if (baitsAsTarget) {
						libraryNGS = new LibraryNGS(lBaitsLibrary.getBaits(), log);
					} else {
						libraryNGS = LibraryNGS.getLibraryNGS(targetLibraryFile, skipNumLines, log);
					}
					if (lBaitsLibrary != null) {
						log.report(ext.getTime() + " Info - mapping baits library...");
						libraryNGS.mapBaits(lBaitsLibrary, baitsAsTarget);
					}
				}
				bamQCs = qcBams(inputbams, outputDir, libraryNGS, filterNGS, normalizeDepthsTo, numThreads, log);
				summarize(bamQCs, outputDir, filterNGS, output, log);
				String librarySummary = outputDir + ext.addToRoot(output, baitsAsTarget ? ".libraryBaitsResults.summary" : ".libraryResults.summary");
				LibraryNGS.summarizeLibraries(libraryNGS, getLibraryReadDepthResultsFiles(bamQCs), librarySummary, filterNGS, log);
				SNPEFF snpeff = new SNPEFF(snpEffLocation, true, true, log);
				snpeff.runSnpEFFCount(inputbams, outputDir + ext.addToRoot(output, new StringBuilder(String.valueOf(baitsAsTarget ? ".libraryBaitsResults.summary" : ".libraryResults.summary")).append("count").toString()), SNPEFF.BUILDS[0], baitsAsTarget ? null : targetLibraryFile, numThreads);
			}
		}
		log.report(ext.getTime() + " Info - finished qc in " + ext.getTimeElapsed(time));
		return bamQCs;
	}

	public static BamQC[] qcBams(String[] inputbams, String outputDir, LibraryNGS libraryNGS, FilterNGS filterNGS, double normalizeDepthsTo, int numThreads, Logger log) {
		BamQC[] bamQCs = new BamQC[inputbams.length];
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		Hashtable<String, Future<BamQC>> tmpResults = new Hashtable<String, Future<BamQC>>();
		for (int i = 0; i < bamQCs.length; i++) {
			tmpResults.put(i + "", executor.submit(new WorkerBamQC(inputbams[i], outputDir, libraryNGS, filterNGS, normalizeDepthsTo, log)));
		}
		for (int i = 0; i < bamQCs.length; i++) {
			try {
				try {
					bamQCs[i] = tmpResults.get(i + "").get();
				} catch (NullPointerException npe) {
					log.reportTimeError("Could not get qc for " + inputbams[i]);
				}
			} catch (InterruptedException e) {
				log.reportError("Error - in file " + inputbams[i]);
				log.reportException(e);
			} catch (ExecutionException e) {
				log.reportError("Error - in file " + inputbams[i]);
				log.reportException(e);
			}
		}
		executor.shutdown();
		try {
			executor.awaitTermination(10L, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			log.reportException(e);
		}
		return bamQCs;
	}

	private static String[] getLibraryReadDepthResultsFiles(BamQC[] bamQCs) {
		String[] libraryReadDepthResultsFiles = new String[bamQCs.length];
		for (int i = 0; i < bamQCs.length; i++) {
			libraryReadDepthResultsFiles[i] = bamQCs[i].getLibraryReadDepthResultsFile();
		}
		return libraryReadDepthResultsFiles;
	}

	public static BamQC qcBam(String inputSamOrBamFile, String outputDir, LibraryNGS libraryNGS, FilterNGS filterNGS, double normalizeDepthsTo, Logger log) {

		LibraryNGS.ReadDepth readDepth = null;
		if (libraryNGS != null) {
			readDepth = new LibraryNGS.ReadDepth(libraryNGS, filterNGS, log);
		}
		BamQC bamQC = initBamQC(inputSamOrBamFile, outputDir, filterNGS);
		System.out.println(bamQC.getLibraryReadDepthResultsFile());
		// if (!Files.exists(bamQC.getLibraryReadDepthResultsFile())) {

		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
		samReaderFactory.validationStringency(ValidationStringency.LENIENT);
		SamReader reader = samReaderFactory.open(new File(inputSamOrBamFile));
		log.report(ext.getTime() + " Info - beginning processing of " + inputSamOrBamFile);
		int numTotalReads = 0;
		for (SAMRecord samRecord : reader) {
			if (!bamQC.invalidAlignment(samRecord)) {
				bamQC.addInsertSize(samRecord, NO_QC_INSERT_HISTOGRAM);
				if ((samRecord.getReadUnmappedFlag()) || (!bamQC.secondaryAlignment(samRecord))) {
					double currentGC = getGCContent(samRecord);
					bamQC.addToGHistoGram(currentGC, NO_QC_GC_HISTOGRAM);
					numTotalReads++;
					if (passesMapQ(samRecord, filterNGS)) {
						bamQC.addNumTotal();
						bamQC.addNumDuplicate(samRecord);
						bamQC.addNumUnique(samRecord);
						bamQC.addNumUnMapped(samRecord);
						if ((!samRecord.getReferenceName().equals("*")) && (libraryNGS != null) && (readDepth != null)) {
							Segment segment = new Segment(Positions.chromosomeNumber(samRecord.getReferenceName(), log), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd());
							if (segment.getChr() >= 0) {
								int[] libraryIndices = libraryNGS.indicesInLibrary(segment);
								if ((!samRecord.getDuplicateReadFlag()) && (!samRecord.getReadUnmappedFlag()) && (libraryIndices != null)) {
									bamQC.addInsertSize(samRecord, QC_INSERT_HISTOGRAM);
									bamQC.addToGHistoGram(currentGC, QC_GC_HISTOGRAM);
									for (int i = 0; i < libraryIndices.length; i++) {
										readDepth.addCounts(libraryIndices[i], libraryNGS.getTargetSegments()[libraryIndices[i]], samRecord, filterNGS);
									}
									bamQC.addNumOnTarget();
								}
							}
						}
						if ((numTotalReads != 0) && (numTotalReads % 1000000 == 0)) {
							float usedMemory = (float) (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
							float freeMemory = (float) Runtime.getRuntime().maxMemory() - usedMemory;
							float maxMemory = (float) Runtime.getRuntime().maxMemory();
							log.report(ext.getTime() + "Info - processed " + numTotalReads + " total reads from " + inputSamOrBamFile + "\nFree memory: " + Math.round(freeMemory / maxMemory * 100.0D) + "%");
						}
					}
				}
			}
		}
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		bamQC.computePercents();
		if ((libraryNGS != null) && (readDepth != null)) {
			double normalizeFactor = 1.0D;
			if (normalizeDepthsTo != 0.0D) {
				normalizeFactor = normalizeDepthsTo * 1000000.0D / bamQC.getNumOnTarget();
				log.report("Info -normalizing depth results from " + bamQC.getInputSamOrBamFile() + " with " + bamQC.getNumOnTarget() + " reads to " + normalizeDepthsTo * 1000000.0D + " reads");
				log.report("Info -normalizing depth results from " + bamQC.getInputSamOrBamFile() + " using a factor of (depth)*" + normalizeFactor);
			}
			LibraryNGS.LibraryReadDepthResults lDepthResults = readDepth.getDepthResults(libraryNGS, filterNGS, normalizeFactor);
			bamQC.setPercentCoverageAtDepth(lDepthResults.getTotalPercentCoveredAtDepth());
			bamQC.setTotalTargetedBasePairs(lDepthResults.getTotalBasePairsTargeted());
			lDepthResults.dump(ext.rootOf(bamQC.getLibraryReadDepthResultsFile(), false) + ".txt", log);
			lDepthResults.serialize(bamQC.getLibraryReadDepthResultsFile());
		} else {
			double[] nans = new double[filterNGS.getReadDepthFilter().length];
			Arrays.fill(nans, Double.NaN);
			bamQC.setPercentCoverageAtDepth(nans);
			bamQC.setTotalTargetedBasePairs(0);
		}
		log.report(ext.getTime() + " Summary for: " + inputSamOrBamFile + "\n");
		log.report(ext.getTime() + " Summary:" + bamQC.getSummary());
		if (bamQC.getNumInvalidAligments() > 0) {
			log.reportError("Warning - " + bamQC.getNumInvalidAligments() + " read(s) were invalid and were skipped in file " + inputSamOrBamFile);
		}
		if (bamQC.getNumSecondaryAlignments() > 0) {
			log.report("Info - " + bamQC.getNumSecondaryAlignments() + " secondary alignment(s) in " + inputSamOrBamFile + " were skipped for computing QC metrics");
		}
		// } else {
		// log.reportTimeInfo("loading pre-computed library results");
		// LibraryNGS.LibraryReadDepthResults lDepthResults = LibraryNGS.LibraryReadDepthResults.load(bamQC.getLibraryReadDepthResultsFile(), false);
		// bamQC.setPercentCoverageAtDepth(lDepthResults.getTotalPercentCoveredAtDepth());
		// bamQC.setTotalTargetedBasePairs(lDepthResults.getTotalBasePairsTargeted());
		// }
		return bamQC;
	}

	private static double getGCContent(SAMRecord samRecord) {
		String sr = samRecord.getReadString();
		int gscs = 0;
		for (int i = 0; i < sr.length(); i++) {
			String base = sr.charAt(i) + "";
			if ((base.equals("G")) || (base.equals("C"))) {
				gscs++;
			}
		}
		return gscs / sr.length();
	}

	private static class WorkerBamQC implements Callable<BamQC> {
		private String inputSamOrBamFile;
		private String outputDir;
		private LibraryNGS libraryNGS;
		private FilterNGS filterNGS;
		private double normalizeDepthsTo;
		private Logger log;

		public WorkerBamQC(String inputSamOrBamFile, String outputDir, LibraryNGS libraryNGS, FilterNGS filterNGS, double normalizeDepthsTo, Logger log) {
			this.inputSamOrBamFile = inputSamOrBamFile;
			this.outputDir = outputDir;
			this.libraryNGS = libraryNGS;
			this.filterNGS = filterNGS;
			this.normalizeDepthsTo = normalizeDepthsTo;
			this.log = log;
		}

		public BamQC call() {
			return BamQC.qcBam(this.inputSamOrBamFile, this.outputDir, this.libraryNGS, this.filterNGS, this.normalizeDepthsTo, this.log);
		}
	}

	private static String[] determineInputFiles(String dir, String commonExt, String fileOfinputSamOrBams, Logger log) {
		String[] inputbams = null;
		if (fileOfinputSamOrBams != null) {
			if (Files.exists(fileOfinputSamOrBams)) {
				log.report(ext.getTime() + " Info - loading files to qc from " + fileOfinputSamOrBams);
				inputbams = HashVec.loadFileToStringArray(fileOfinputSamOrBams, false, new int[] { 0 }, false);
				if (!Files.exists("", inputbams)) {
					log.reportError("Error - found files that do not exist in " + fileOfinputSamOrBams);
					inputbams = null;
				}
			} else {
				log.reportError("Error - could not find file " + fileOfinputSamOrBams);
			}
		} else {
			log.report(ext.getTime() + " Info - loading all files from dir " + dir + " with extension " + commonExt);
			inputbams = Files.toFullPaths(Files.list(dir, "", commonExt, true, false), dir);
		}
		return inputbams;
	}

	private static double computePercent(int top, int bottom) {
		return bottom > 0 ? (double) top / bottom : 0;
	}

	private static BamQC initBamQC(String inputSamOrBamFile, String outputDir, FilterNGS filterNGS) {
		if (filterNGS == null) {
			return new BamQC(inputSamOrBamFile, outputDir, new FilterNGS(0.0, 0.0, new int[1]));
		}
		BamQC bamQC = new BamQC(inputSamOrBamFile, outputDir, filterNGS);
		return bamQC;
	}

	private static void summarize(BamQC[] bamQCs, String outputDir, FilterNGS filterNGS, String output, Logger log) {
		if (outputDir == null) {
			for (int i = 0; i < bamQCs.length; i++) {
				if (bamQCs[i] != null) {
					output = ext.parseDirectoryOfFile(bamQCs[i].getInputSamOrBamFile()) + ext.removeDirectoryInfo(output);
					break;
				}
			}
		} else {
			output = outputDir + output;
		}
		summarizeQC(bamQCs, filterNGS, output, log);
		summarizeGCContent(bamQCs, ext.addToRoot(output, "raw"), NO_QC_GC_HISTOGRAM, log);
		summarizeGCContent(bamQCs, ext.addToRoot(output, "onTarget"), QC_GC_HISTOGRAM, log);
		summarizeInsertSize(bamQCs, ext.addToRoot(output, "raw"), NO_QC_INSERT_HISTOGRAM, log);
		summarizeInsertSize(bamQCs, ext.addToRoot(output, "onTarget"), QC_INSERT_HISTOGRAM, log);
	}

	private static void summarizeQC(BamQC[] bamQCs, FilterNGS filterNGS, String output, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.print(Array.toStr(QC_HEADER));
			for (int i = 0; i < filterNGS.getReadDepthFilter().length; i++) {
				writer.print("\tPercent Coverage At " + filterNGS.getReadDepthFilter()[i]);
			}
			writer.println();
			for (int i = 0; i < bamQCs.length; i++) {
				if (bamQCs[i] != null) {
					writer.println(bamQCs[i].getSummary());
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}
	}

	private static void summarizeGCContent(BamQC[] bamQCs, String output, int which, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(ext.rootOf(output, false) + ".gcs"));
			double[] bins = null;
			writer.print(Array.toStr(HIST_HEADER));
			for (int i = 0; i < bamQCs.length; i++) {
				if (bamQCs[i] != null) {
					writer.print("\t" + ext.rootOf(bamQCs[i].getInputSamOrBamFile()));
					if (bins == null) {
						bins = bamQCs[i].getGcHistogram(which).getBins();
					}
				}
			}
			writer.println();
			for (int i = 0; i < bins.length; i++) {
				writer.print(bins[i]);
				for (int j = 0; j < bamQCs.length; j++) {
					if (bamQCs[j] != null) {
						writer.print("\t" + bamQCs[j].getGcHistogram(which).getCounts()[i]);
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}
	}

	private static void summarizeInsertSize(BamQC[] bamQCs, String output, int which, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(ext.rootOf(output, false) + ".insert"));
			double[] bins = null;
			writer.print(Array.toStr(HIST_HEADER));
			for (int i = 0; i < bamQCs.length; i++) {
				if (bamQCs[i] != null) {
					writer.print("\t" + ext.rootOf(bamQCs[i].getInputSamOrBamFile()));
					if (bins == null) {
						bins = bamQCs[i].getInsertHistogram(which).getBins();
					}
				}
			}
			writer.println();
			for (int i = 0; i < bins.length; i++) {
				writer.print(bins[i]);
				for (int j = 0; j < bamQCs.length; j++) {
					if (bamQCs[j] != null) {
						writer.print("\t" + bamQCs[j].getInsertHistogram(which).getCounts()[i]);
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}
	}

	private static boolean passesMapQ(SAMRecord samRecord, FilterNGS filterNGS) {
		boolean passesMapQ = false;
		if (filterNGS.getMappingQualityFilter() == 0.0) {
			passesMapQ = true;
		} else if (((filterNGS.getMappingQualityFilter() == 0.0) || (samRecord.getMappingQuality() != 255)) && (samRecord.getMappingQuality() >= filterNGS.getMappingQualityFilter())) {
			passesMapQ = true;
		}
		return passesMapQ;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = null;
		String commonExt = ".bam";
		String fileOfinputSamOrBams = null;
		String targetLibraryFile = null;
		String baitLibraryFile = null;
		String output = "bamQC.txt";
		String outputDir = null;
		int skipNumLines = 2;
		int numThreads = 4;
		double mappingQuality = 0;
		double phreadScore = 0.0;
		int[] readDepth = { 0, 1, 2, 3, 4, 10, 20, 30, 40 };
		boolean baitsAsTarget = false;
		double normalizeDepthsTo = 0.0D;
		String snpEffLocation = null;
		String logfile = "bamQC.log";

		String usage = "\nseq.BamQC requires 0-1 arguments\n";
		usage = usage + "    NOTE: all bam files are assumed to be sorted:\n";
		usage = usage + "   (1) directory containing bam or sam files to qc (i.e. dir=" + dir + " (no default))\n";
		usage = usage + "   (2) full path to a file listing bam or sam files (one per line, in first column, no header) to qc (i.e. fileOfinputSamOrBams=" + fileOfinputSamOrBams + " (no default))\n";
		usage = usage + "   (3) full path to a target library file to compute on target percentage (i.e. target=" + targetLibraryFile + " (no default))\n";
		usage = usage + "   (3) full path to a bait library file to map to the target library (i.e. baits=" + baitLibraryFile + " (no default))\n";
		usage = usage + "   (3) full path to output directory (i.e. outDir=" + outputDir + " (defualts to directory of first file qc -ed))\n";

		usage = usage + "   (4) output summary file name (i.e. output=" + output + " (default - defaults to directory of first file qc -ed))\n";
		usage = usage + "   (5) file extension to search for in a directory (i.e. commonExt=" + commonExt + " (default))\n";
		usage = usage + "   (6) if using a target library file, number of lines to skip before targets are listed (i.e. skip=" + skipNumLines + " (default))\n";
		usage = usage + "   (7) number of threads  (i.e. numThreads=" + numThreads + " (default))\n";
		usage = usage + "   (8) log file  (i.e. log=" + logfile + " (default))\n";
		usage = usage + "   (9) mapping quality score filter to compute qc metrics with (i.e. mapQ=" + mappingQuality + " (default))\n";
		usage = usage + "   (10) phread score filter to compute qc metrics with (i.e. phread=" + phreadScore + " (default))\n";
		usage = usage + "   (11) read depth filter (comma-delimited if multiple) to compute qc metrics with (i.e. readDepth=" + Array.toStr(Array.toStringArray(readDepth), ",") + " (default))\n";
		usage = usage + "   (12) compute QC against the baits library (i.e. -baits" + baitsAsTarget + " (default))\n";
		usage = usage + "   (13) normalize coverage to number of Reads (in millions)  (i.e. normReads=" + normalizeDepthsTo + " (default, no normalizing))\n";
		usage = usage + "   (14) full path to the SNP EFF directory (i.e. snpEff= ( no default))\n";
		for (int i = 0; i < args.length; i++) {
			if ((args[i].equals("-h")) || (args[i].equals("-help")) || (args[i].equals("/h")) || (args[i].equals("/help"))) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("fileOfinputSamOrBams=")) {
				fileOfinputSamOrBams = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("target=")) {
				targetLibraryFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("-baits")) {
				baitsAsTarget = true;
				numArgs--;
			} else if (args[i].startsWith("outDir=")) {
				outputDir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("output=")) {
				output = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("commonExt=")) {
				commonExt = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("baits=")) {
				baitLibraryFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("skip=")) {
				skipNumLines = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("numThreads=")) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("mapQ=")) {
				mappingQuality = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("phread=")) {
				phreadScore = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("phread=")) {
				phreadScore = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("normReads=")) {
				normalizeDepthsTo = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("snpEff=")) {
				snpEffLocation = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("readDepth=")) {
				readDepth = Array.toIntArray(ext.parseStringArg(args[i], "").split(","));
				phreadScore = ext.parseDoubleArg(args[i]);
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
			Logger log = new Logger(logfile);
			FilterNGS filterNGS = new FilterNGS(mappingQuality, phreadScore, readDepth);
			qcBams(dir, outputDir, commonExt, fileOfinputSamOrBams, targetLibraryFile, baitLibraryFile, skipNumLines, filterNGS, numThreads, output, snpEffLocation, baitsAsTarget, normalizeDepthsTo, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
