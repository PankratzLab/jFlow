package seq.qc;

import common.Array;
import common.Files;
import common.Logger;
import common.ext;
import filesys.Segment;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;


public class LibraryNGS implements Serializable {
	private static final long serialVersionUID = 1L;
	private Segment[] targetSegments;
	private int[] baitsPerTarget;
	private double[] baitGCContent;
	private Logger log;

	public LibraryNGS(Segment[] targetSegments, Logger log) {
		this.targetSegments = targetSegments;
		this.log = log;
	}

	public void mapBaits(BaitsLibrary baitsLibrary, boolean baitsAsTarget) {
		baitsLibrary.mapToLibrary(this, baitsAsTarget);
	}

	public int[] getBaitsPerTarget() {
		return this.baitsPerTarget;
	}

	public void setBaitsPerTarget(int[] baitsPerTarget) {
		this.baitsPerTarget = baitsPerTarget;
	}

	public double[] getBaitGCContent() {
		return this.baitGCContent;
	}

	public void setBaitGCContent(double[] baitGCContent) {
		this.baitGCContent = baitGCContent;
	}

	public void setTargetSegments(Segment[] targetSegments) {
		this.targetSegments = targetSegments;
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static LibraryNGS load(String filename, boolean jar) {
		return (LibraryNGS) Files.readSerial(filename, jar, true);
	}

	public int[] indicesInLibrary(Segment segment) {
		return Segment.binarySearchForAllOverLappingIndices(segment, this.targetSegments);
	}

	public Segment getTargetSegmentsAt(int libraryIndex) {
		return getTargetSegments()[libraryIndex];
	}

	public Segment[] getTargetSegments() {
		return this.targetSegments;
	}

	public boolean inLibrary(Segment segment) {
		return indicesInLibrary(segment) != null;
	}

	public Logger getLog() {
		return this.log;
	}

	public static LibraryNGS getLibraryNGS(String libraryNGSFile, int skipNumLines, Logger log) {
		LibraryNGS libraryNGS = null;
		String serLibrary = ext.rootOf(libraryNGSFile, false) + ".ser";
		if (!Files.exists(serLibrary)) {
			libraryNGS = new LibraryNGS(Segment.loadRegions(libraryNGSFile, 0, 1, 2, skipNumLines, true, true, true,0), log);
			libraryNGS.serialize(serLibrary);
		} else {
			libraryNGS = load(serLibrary, false);
		}
		return libraryNGS;
	}

	public static class ReadDepth {
		public static final String[] BASES = { "A", "T", "C", "G", "N", "Other" };
		public static final String[] GC = { "G", "C" };
		private int[][] depths;
		private int[][] GCs;
		private int[] numReads;
		private double[] phreadsAvg;
		private double[] mapQsAvg;
		private double[] averageInsertSize;
		//private Histogram.DynamicHistogram[] dyHistograms;
		private Logger log;

		public ReadDepth(LibraryNGS libraryNGS, FilterNGS filterNGS, Logger log) {
			this.GCs = initInt(libraryNGS.getTargetSegments(), log);
			this.depths = initInt(libraryNGS.getTargetSegments(), log);
			this.numReads = new int[libraryNGS.getTargetSegments().length];
			this.phreadsAvg = new double[libraryNGS.getTargetSegments().length];
			this.mapQsAvg = new double[libraryNGS.getTargetSegments().length];
			this.averageInsertSize = new double[libraryNGS.getTargetSegments().length];
			//this.dyHistograms = new Histogram.DynamicHistogram[filterNGS.getReadDepthFilter().length];
			this.log = log;
		}

		public void computeAverageQualities() {
			for (int i = 0; i < this.depths.length; i++) {
				int numBases = Array.sum(this.depths[i]);
				if (numBases > 0) {
					this.phreadsAvg[i] /= numBases;
					this.mapQsAvg[i] /= this.numReads[i];
					this.averageInsertSize[i] /= this.numReads[i];
				}
			}
		}

		public double getAvgPhreadAt(int i) {
			return this.phreadsAvg[i];
		}

		public double getAvgMapQAt(int i) {
			return this.mapQsAvg[i];
		}

		public double getAvgInsertSizeAt(int i) {
			return this.averageInsertSize[i];
		}

		public void normalize(double normalizeFactor) {
			if ((normalizeFactor != 1.0D) && (normalizeFactor > 0.0D)) {
				for (int i = 0; i < this.depths.length; i++) {
					for (int j = 0; j < this.depths[i].length; j++) {
						this.depths[i][j] = ((int) Math.round(this.depths[i][j] * normalizeFactor));
						this.GCs[i][j] = ((int) Math.round(this.GCs[i][j] * normalizeFactor));
					}
				}
			}
		}

		public double getMedianCoverage() {
			double medianCoverage = (0.0D / 0.0D);
			int numBases = 0;
			ArrayList<Double> allCounts = new ArrayList<Double>(75000000);
			for (int i = 0; i < this.depths.length; i++) {
				for (int j = 0; j < this.depths[i].length; j++) {
					numBases++;
					allCounts.add(Double.valueOf(this.depths[i][j]));
				}
			}
			if (numBases > 0) {
				medianCoverage = Array.median(Array.toDoubleArray(allCounts));
			}
			return medianCoverage;
		}

		public void addCounts(int libraryIndex, Segment segment, SAMRecord samRecord, FilterNGS filterNGS) {
			List<AlignmentBlock> alList = samRecord.getAlignmentBlocks();
			double[] phreads = LibraryNGS.parseBytesToDouble(samRecord.getBaseQualities());
			String[] bases = LibraryNGS.parseToString(samRecord.getReadString());
			this.numReads[libraryIndex] += 1;
			if (samRecord.getProperPairFlag()) {
				this.averageInsertSize[libraryIndex] += Math.abs(samRecord.getInferredInsertSize());
			}
			if (samRecord.getMappingQuality() != 255) {
				this.mapQsAvg[libraryIndex] += samRecord.getMappingQuality();
			}
			for (int i = 0; i < alList.size(); i++) {
				AlignmentBlock alignmentBlock = (AlignmentBlock) alList.get(i);
				BlockParser blockParser = new BlockParser(alignmentBlock, segment, this.log);
				blockParser.parseIndices();
				if (blockParser.getTargetStartArrayIndex() >= 0) {
					int seqIndex = blockParser.getSequenceStartArrayIndex();
					for (int j = blockParser.getTargetStartArrayIndex(); j < blockParser.getTargetStopArrayIndex(); j++) {
						if ((filterNGS.getPhreadScoreFilter() == 0.0D) || (phreads[seqIndex] >= filterNGS.getPhreadScoreFilter())) {
							this.depths[libraryIndex][j] += 1;
							this.phreadsAvg[libraryIndex] += phreads[seqIndex];
							if ((bases[seqIndex].equals(GC[0])) || (bases[seqIndex].equals(GC[1]))) {
								this.GCs[libraryIndex][j] += 1;
							}
						}
						seqIndex++;
					}
				}
			}
		}

		private static class BlockParser {
			private AlignmentBlock alignmentBlock;
			private Segment targetSegment;
			private int sequenceStartArrayIndex;
			private int targetStartArrayIndex;
			private int targetStopArrayIndex;
			private Logger log;

			public BlockParser(AlignmentBlock alignmentBlock, Segment targetSegment, Logger log) {
				this.alignmentBlock = alignmentBlock;
				this.targetSegment = targetSegment;
				this.log = log;
			}

			public void parseIndices() {
				int alignmentRefStart = this.alignmentBlock.getReferenceStart();
				int alignmentRefStop = alignmentRefStart + this.alignmentBlock.getLength() - 1;
				this.sequenceStartArrayIndex = (this.alignmentBlock.getReadStart() - 1);
				if (alligns(alignmentRefStart, alignmentRefStop)) {
					if ((alignmentRefStart >= this.targetSegment.getStart()) && (alignmentRefStop <= this.targetSegment.getStop())) {
						this.targetStartArrayIndex = (alignmentRefStart - this.targetSegment.getStart());
						this.targetStopArrayIndex = (alignmentRefStart - this.targetSegment.getStart() + this.alignmentBlock.getLength());
					} else if ((this.targetSegment.getStart() >= alignmentRefStart) && (this.targetSegment.getStop() <= alignmentRefStop)) {
						this.targetStartArrayIndex = 0;
						this.targetStopArrayIndex = this.targetSegment.getSize();
					} else if ((alignmentRefStart >= this.targetSegment.getStart()) && (alignmentRefStop >= this.targetSegment.getStop())) {
						this.targetStartArrayIndex = (alignmentRefStart - this.targetSegment.getStart());
						this.targetStopArrayIndex = this.targetSegment.getSize();
					} else if ((alignmentRefStart <= this.targetSegment.getStart()) && (alignmentRefStop <= this.targetSegment.getStop())) {
						this.targetStartArrayIndex = 0;
						this.targetStopArrayIndex = (alignmentRefStart + this.alignmentBlock.getLength() - this.targetSegment.getStart());
					} else {
						this.log.reportError("Internal Error - could not determine block alignments for " + alignmentRefStart + "\t" + alignmentRefStop + "\t" + this.targetSegment.getUCSClocation());
						this.targetStartArrayIndex = -1;
						this.targetStopArrayIndex = -1;
					}
					if ((this.targetStartArrayIndex < 0) || (this.targetStopArrayIndex < 0)) {
						this.log.reportError("Internal Error - could not determine block alignments for " + alignmentRefStart + "\t" + alignmentRefStop + "\t" + this.targetSegment.getUCSClocation());
						this.targetStartArrayIndex = -1;
						this.targetStopArrayIndex = -1;
					}
				} else {
					this.targetStartArrayIndex = -2;
					this.targetStopArrayIndex = -2;
				}
			}

			public boolean alligns(int alignmentRefStart, int alignmentRefStop) {
				return new Segment(this.targetSegment.getChr(), alignmentRefStart, alignmentRefStop).overlaps(this.targetSegment);
			}

//			public Segment getTargetSegment() {
//				return this.targetSegment;
//			}

			public int getSequenceStartArrayIndex() {
				return this.sequenceStartArrayIndex;
			}

			public int getTargetStartArrayIndex() {
				return this.targetStartArrayIndex;
			}

			public int getTargetStopArrayIndex() {
				return this.targetStopArrayIndex;
			}
		}

		public LibraryNGS.LibraryReadDepthResults getDepthResults(LibraryNGS libraryNGS, FilterNGS filterNGS, double normalizeFactor) {
			LibraryNGS.LibraryReadDepthResults readDepthResults = new LibraryNGS.LibraryReadDepthResults(this, libraryNGS.getTargetSegments(), filterNGS);
			readDepthResults.populateSummary(libraryNGS, this, normalizeFactor, filterNGS);
			return readDepthResults;
		}

		public int[][] getDepths() {
			return this.depths;
		}

		public int[] getNumReads() {
			return this.numReads;
		}

		public int[][] getGCs() {
			return this.GCs;
		}

		public Logger getLog() {
			return this.log;
		}

		private static int[][] initInt(Segment[] targetSegments, Logger log) {
			int[][] depths = new int[targetSegments.length][];
			for (int i = 0; i < depths.length; i++) {
				depths[i] = new int[targetSegments[i].getSize()];
			}
			return depths;
		}
	}

	public static class LibraryReadDepthResults implements Serializable {
		public static final String[] SummaryHeader = { "UCSC", "numBasePairsTargeted", "numBasePairsSequenced", "numGCSequenced", "numReads", "averageCoverage", "averageGC", "averagePerBaseBias", "averagePhreadScore", "averageMapQ", "averageInsertSize", "numBaitsPerTarget", "AverageBaitGC" };
		public static final String[] ADD_HEADER = { "Percent_Covered_at_depth", "Percent_GC_at_depth", "num_GC_at_depth" };
		public static final String LIBRARY_READ_DEPTH_EXTENSION = ".libraryResults.ser";
		public static final String SUMMARY = ".libraryResults.summary";
		public static final String SUMMARY_BAITS = ".libraryBaitsResults.summary";
		private static final long serialVersionUID = 1L;
		private LibraryNGS.TargetReadDepthResults[] targetReadDepthResults;
		private double[] totalPercentCoveredAtDepth;
		private double[] totalPercentGCAtDepth;
		private Segment[] targetSegments;
		private FilterNGS filterNGS;
		private int totalBasePairsTargeted;
		private int totalBasePairsSequenced;

		public LibraryReadDepthResults(LibraryNGS.ReadDepth readDepth, Segment[] targetSegments, FilterNGS filterNGS) {
			this.targetReadDepthResults = new LibraryNGS.TargetReadDepthResults[readDepth.getDepths().length];
			this.targetSegments = targetSegments;
			this.filterNGS = filterNGS;
			this.totalBasePairsTargeted = 0;
			this.totalBasePairsSequenced = 0;
		}

		public LibraryReadDepthResults(Segment[] targetSegments, FilterNGS filterNGS) {
			this.targetReadDepthResults = new LibraryNGS.TargetReadDepthResults[targetSegments.length];
			this.targetSegments = targetSegments;
			this.filterNGS = filterNGS;
			this.totalBasePairsTargeted = 0;
			this.totalBasePairsSequenced = 0;
		}

		public LibraryReadDepthResults(LibraryNGS.TargetReadDepthResults[] targetReadDepthResults, double[] totalPercentCoveredAtDepth, double[] totalPercentGCAtDepth, Segment[] targetSegments, FilterNGS filterNGS, int totalBasePairsTargeted, int totalBasePairsSequenced) {
			this.targetReadDepthResults = targetReadDepthResults;
			this.totalPercentCoveredAtDepth = totalPercentCoveredAtDepth;
			this.totalPercentGCAtDepth = totalPercentGCAtDepth;
			this.targetSegments = targetSegments;
			this.filterNGS = filterNGS;
			this.totalBasePairsTargeted = totalBasePairsTargeted;
			this.totalBasePairsSequenced = totalBasePairsSequenced;
		}

		public String getHeader() {
			String header = "";
			header = header + Array.toStr(SummaryHeader);
			for (int i = 0; i < ADD_HEADER.length; i++) {
				for (int j = 0; j < this.filterNGS.getReadDepthFilter().length; j++) {
					header = header + "\t" + ADD_HEADER[i] + "_" + this.filterNGS.getReadDepthFilter()[j];
				}
			}
			return header;
		}

		public double[] getTotalPercentCoveredAtDepth() {
			return this.totalPercentCoveredAtDepth;
		}

		public double[] getTotalPercentGCAtDepth() {
			return totalPercentGCAtDepth;
		}

		public LibraryNGS.TargetReadDepthResults[] getTargetReadDepthResults() {
			return this.targetReadDepthResults;
		}

		public int getTotalBasePairsTargeted() {
			return this.totalBasePairsTargeted;
		}

		public void serialize(String filename) {
			Files.writeSerial(this, filename);
		}

		public static LibraryReadDepthResults load(String filename, boolean jar) {
			return (LibraryReadDepthResults) Files.readSerial(filename, jar, true);
		}

		public int getTotalNumReads() {
			int totalNumReads = 0;
			for (int i = 0; i < this.targetReadDepthResults.length; i++) {
				totalNumReads += this.targetReadDepthResults[i].getNumReadsSequenced();
			}
			return totalNumReads;
		}

		public String getSummaryFor(int libraryIndex) {
			String summary = "";
			summary = summary + this.targetSegments[libraryIndex].getUCSClocation();
			summary = summary + "\t" + this.targetReadDepthResults[libraryIndex].getSummmary();
			return summary;
		}

		public void dump(String filename, Logger log) {
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(filename));
				writer.println(getHeader());
				for (int i = 0; i < this.targetReadDepthResults.length; i++) {
					writer.println(getSummaryFor(i));
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + filename);
				log.reportException(e);
			}
		}

		public void setTotalPercentCoveredAtDepth(double[] totalPercentCoveredAtDepth) {
			this.totalPercentCoveredAtDepth = totalPercentCoveredAtDepth;
		}

		public void setTotalPercentGCAtDepth(double[] totalPercentGCAtDepth) {
			this.totalPercentGCAtDepth = totalPercentGCAtDepth;
		}

		public void populateSummary(LibraryNGS libraryNGS, LibraryNGS.ReadDepth readDepth, double normalizeFactor, FilterNGS filterNGS) {
			this.totalPercentCoveredAtDepth = new double[filterNGS.getReadDepthFilter().length];
			this.totalPercentGCAtDepth = new double[filterNGS.getReadDepthFilter().length];
			double avgCoverage = readDepth.getMedianCoverage();
			readDepth.computeAverageQualities();
			readDepth.normalize(normalizeFactor);
			for (int i = 0; i < readDepth.getDepths().length; i++) {
				this.totalBasePairsTargeted += readDepth.getDepths()[i].length;
				if ((libraryNGS.getBaitsPerTarget() != null) && (libraryNGS.getBaitGCContent() != null)) {
					this.targetReadDepthResults[i] = new LibraryNGS.TargetReadDepthResults(readDepth.getNumReads()[i], readDepth.getDepths()[i].length, libraryNGS.getBaitsPerTarget()[i], libraryNGS.getBaitGCContent()[i], filterNGS);
				} else {
					this.targetReadDepthResults[i] = new LibraryNGS.TargetReadDepthResults(readDepth.getNumReads()[i], readDepth.getDepths()[i].length, 0, 0.0D, filterNGS);
				}
				this.targetReadDepthResults[i].setNumSamples(1);
				this.targetReadDepthResults[i].setAverageCoverage(Array.mean(readDepth.getDepths()[i]));
				for (int j = 0; j < readDepth.getDepths()[i].length; j++) {
					this.targetReadDepthResults[i].parseCount(readDepth.getDepths()[i][j], readDepth.getGCs()[i][j]);
					this.totalBasePairsSequenced += readDepth.getDepths()[i][j];
				}
				if (this.targetReadDepthResults[i].getNumBPSequencedInTarget() > 0) {
					this.targetReadDepthResults[i].setAverageGCContent(this.targetReadDepthResults[i].getNumGCSequencedInTarget() / this.targetReadDepthResults[i].getNumBPSequencedInTarget());
				}
				this.targetReadDepthResults[i].computePerBaseBias(avgCoverage);
				this.targetReadDepthResults[i].setAvgPhreadScore(readDepth.getAvgPhreadAt(i));
				this.targetReadDepthResults[i].setAvgMapQ(readDepth.getAvgMapQAt(i));
				this.targetReadDepthResults[i].setAverageInsertSize(readDepth.getAvgInsertSizeAt(i));
			}
			summarizePercents(filterNGS);
		}

		public void summarizePercents(FilterNGS filterNGS) {
			for (int i = 0; i < this.targetReadDepthResults.length; i++) {
				this.targetReadDepthResults[i].computePercents();
				for (int j = 0; j < filterNGS.getReadDepthFilter().length; j++) {
					this.totalPercentCoveredAtDepth[j] += this.targetReadDepthResults[i].getnumBPSequencedAtDepth(j);
					this.totalPercentGCAtDepth[j] += this.targetReadDepthResults[i].getPercentGCSequencedAtDepth(j);
				}
			}
			for (int i = 0; i < filterNGS.getReadDepthFilter().length; i++) {
				this.totalPercentCoveredAtDepth[i] /= this.totalBasePairsTargeted;
				this.totalPercentGCAtDepth[i] /= this.totalBasePairsSequenced;
			}
		}
	}

	public static void summarizeLibraries(LibraryNGS originalLibraryNGS, String[] libraryReadDepthResultFiles, String output, FilterNGS filterNGS, Logger log) {
		LibraryReadDepthResults summaryResults = new LibraryReadDepthResults(originalLibraryNGS.getTargetSegments(), filterNGS);
		summaryResults.setTotalPercentCoveredAtDepth(new double[filterNGS.getReadDepthFilter().length]);
		summaryResults.setTotalPercentGCAtDepth(new double[filterNGS.getReadDepthFilter().length]);
		for (int i = 0; i < libraryReadDepthResultFiles.length; i++) {
			LibraryReadDepthResults curReadDepthResults = LibraryReadDepthResults.load(libraryReadDepthResultFiles[i], false);
			for (int j = 0; j < summaryResults.getTargetReadDepthResults().length; j++) {
				if (summaryResults.getTargetReadDepthResults()[j] == null) {
					summaryResults.getTargetReadDepthResults()[j] = new TargetReadDepthResults(0, curReadDepthResults.getTargetReadDepthResults()[j].getNumBPInTarget(), originalLibraryNGS.getBaitsPerTarget()[j], originalLibraryNGS.getBaitGCContent()[j], filterNGS);
					summaryResults.getTargetReadDepthResults()[j].setNumSamples(libraryReadDepthResultFiles.length);
				}
				summaryResults.getTargetReadDepthResults()[j].addFromAnother(curReadDepthResults.getTargetReadDepthResults()[j]);
			}
		}
		summaryResults.summarizePercents(filterNGS);
		summaryResults.dump(output, log);
	}

	public static class TargetReadDepthResults implements Serializable {
		private static final long serialVersionUID = 1L;
		private FilterNGS filterNGS;
		private int numBPSequencedInTarget;
		private int numBPInTarget;
		private int numReadsSequenced;
		private int numBaits;
		private int numSamples;
		private int numGCSequencedInTarget;
		private double baitGC;
		private int[] numBPSequencedAtDepth;
		private int[] numGCSequencedAtDepth;
		private double[] percentGCSequencedAtDepth;
		private double[] percentAtDepth;
		private double[] percentGCAtDepth;
		private double[] averagedNumGCAtDepth;
		private double averageCoverage;
		private double averageGCContent;
		private double averagePerBaseBias;
		private double averagePhreadScore;
		private double averageMapQ;
		private double averageInsertSize;

		public TargetReadDepthResults(FilterNGS filterNGS, int numBPSequencedInTarget, int numBPInTarget, int numReadsSequenced, int numBaits, int numSamples, int numGCSequencedInTarget, double baitGC, int[] numBPSequencedAtDepth, int[] numGCSequencedAtDepth, double[] percentGCSequencedAtDepth, double[] percentAtDepth, double[] percentGCAtDepth, double[] averagedNumGCAtDepth, double averageCoverage, double averageGCContent, double averagePerBaseBias, double averagePhreadScore, double averageMapQ, double averageInsertSize) {
			this.filterNGS = filterNGS;
			this.numBPSequencedInTarget = numBPSequencedInTarget;
			this.numBPInTarget = numBPInTarget;
			this.numReadsSequenced = numReadsSequenced;
			this.numBaits = numBaits;
			this.numSamples = numSamples;
			this.numGCSequencedInTarget = numGCSequencedInTarget;
			this.baitGC = baitGC;
			this.numBPSequencedAtDepth = numBPSequencedAtDepth;
			this.numGCSequencedAtDepth = numGCSequencedAtDepth;
			this.percentGCSequencedAtDepth = percentGCSequencedAtDepth;
			this.percentAtDepth = percentAtDepth;
			this.percentGCAtDepth = percentGCAtDepth;
			this.averagedNumGCAtDepth = averagedNumGCAtDepth;
			this.averageCoverage = averageCoverage;
			this.averageGCContent = averageGCContent;
			this.averagePerBaseBias = averagePerBaseBias;
			this.averagePhreadScore = averagePhreadScore;
			this.averageMapQ = averageMapQ;
			this.averageInsertSize = averageInsertSize;
		}

		public TargetReadDepthResults(int numReadsSequenced, int numBPInTarget, int numBaits, double baitGC, FilterNGS filterNGS) {
			this.filterNGS = filterNGS;
			this.numBPSequencedInTarget = 0;
			this.numGCSequencedInTarget = 0;
			this.averageCoverage = 0.0D;
			this.averageGCContent = 0.0D;
			this.averagePerBaseBias = 0.0D;
			this.averagePhreadScore = 0.0D;
			this.averageMapQ = 0.0D;
			this.averageInsertSize = 0.0D;
			this.numBPInTarget = numBPInTarget;
			this.numBaits = numBaits;
			this.baitGC = baitGC;
			this.numBPSequencedAtDepth = new int[filterNGS.getReadDepthFilter().length];
			this.numGCSequencedAtDepth = new int[filterNGS.getReadDepthFilter().length];
			this.percentGCSequencedAtDepth = new double[filterNGS.getReadDepthFilter().length];
			this.percentAtDepth = new double[filterNGS.getReadDepthFilter().length];
			this.percentGCAtDepth = new double[filterNGS.getReadDepthFilter().length];
			this.averagedNumGCAtDepth = new double[filterNGS.getReadDepthFilter().length];
			this.numReadsSequenced = numReadsSequenced;
		}

		public String getSummmary() {
			String summary = "";
			summary = summary + this.numBPInTarget;
			summary = summary + "\t" + this.numBPSequencedInTarget;
			summary = summary + "\t" + this.numGCSequencedInTarget;
			summary = summary + "\t" + this.numReadsSequenced;
			summary = summary + "\t" + this.averageCoverage;
			summary = summary + "\t" + this.averageGCContent;
			summary = summary + "\t" + this.averagePerBaseBias;
			summary = summary + "\t" + this.averagePhreadScore;
			summary = summary + "\t" + this.averageMapQ;
			summary = summary + "\t" + this.averageInsertSize;
			summary = summary + "\t" + this.numBaits;
			summary = summary + "\t" + this.baitGC;
			summary = summary + "\t" + Array.toStr(this.percentAtDepth);
			summary = summary + "\t" + Array.toStr(this.percentGCAtDepth);
			summary = summary + "\t" + Array.toStr(this.averagedNumGCAtDepth);
			return summary;
		}

		public void setNumGCSequencedInTarget(int numGCSequencedInTarget) {
			this.numGCSequencedInTarget = numGCSequencedInTarget;
		}

		public void addNumBPSequencedInTarget(int num) {
			this.numBPSequencedInTarget += num;
		}

		public void addGCSequencedInTarget(int num) {
			this.numGCSequencedInTarget += num;
		}

		public int getnumBPSequencedAtDepth(int dIndex) {
			return this.numBPSequencedAtDepth[dIndex];
		}

		public double getAverageGCContent() {
			return this.averageGCContent;
		}

		public void setAverageGCContent(double averageGCContent) {
			this.averageGCContent = averageGCContent;
		}

		public void setAverageCoverage(double averageCoverage) {
			this.averageCoverage = averageCoverage;
		}

		public double getPercentGCSequencedAtDepth(int dIndex) {
			return this.percentGCSequencedAtDepth[dIndex];
		}

		public int getNumGCSequencedAtDepth(int dIndex) {
			return this.numGCSequencedAtDepth[dIndex];
		}

		public double getAverageInsertSize() {
			return this.averageInsertSize;
		}

		public void setAverageInsertSize(double averageInsertSize) {
			this.averageInsertSize = averageInsertSize;
		}

		public int getNumBPSequencedInTarget() {
			return this.numBPSequencedInTarget;
		}

		public void setNumSamples(int numSamples) {
			this.numSamples = numSamples;
		}

		public int getNumGCSequencedInTarget() {
			return this.numGCSequencedInTarget;
		}

		public double getAverageCoverage() {
			return this.averageCoverage;
		}

		public double getAveragePerBaseBias() {
			return this.averagePerBaseBias;
		}

		public double getAvgPhreadScore() {
			return this.averagePhreadScore;
		}

		public void setAvgPhreadScore(double avgPhreadScore) {
			this.averagePhreadScore = avgPhreadScore;
		}

		public double getAvgMapQ() {
			return this.averageMapQ;
		}

		public void setAvgMapQ(double avgMapQ) {
			this.averageMapQ = avgMapQ;
		}

		public void computePercents() {
			double divideDepthBy = this.numBPInTarget * this.numSamples;
			this.averageCoverage = (this.numBPInTarget > 0 ? this.averageCoverage / this.numSamples : 0.0D);
			this.averageGCContent = (this.numBPInTarget > 0 ? this.averageGCContent / this.numSamples : 0.0D);
			this.averagePerBaseBias = (this.numBPInTarget > 0 ? this.averagePerBaseBias / this.numSamples : 0.0D);
			this.averageMapQ = (this.numBPInTarget > 0 ? this.averageMapQ / this.numSamples : 0.0D);
			this.averagePhreadScore = (this.numBPInTarget > 0 ? this.averagePhreadScore / this.numSamples : 0.0D);
			this.averageInsertSize = (this.numBPInTarget > 0 ? this.averageInsertSize / this.numSamples : 0.0D);
			if (Double.isNaN(this.averageGCContent)) {
				System.out.println(this.numBPInTarget + "\t" + this.averageGCContent + "\t" + this.numSamples + "\t" + this.numBPSequencedInTarget);
			}
			for (int i = 0; i < this.percentAtDepth.length; i++) {
				this.percentAtDepth[i] = (this.numBPInTarget > 0 ? this.numBPSequencedAtDepth[i] / divideDepthBy : 0.0D);
				this.percentGCAtDepth[i] = (this.numBPInTarget > 0 ? this.percentGCSequencedAtDepth[i] / divideDepthBy : 0.0D);
				this.averagedNumGCAtDepth[i] = (this.numBPInTarget > 0 ? this.numGCSequencedAtDepth[i] / this.numSamples : 0.0D);
			}
		}

		public void computePerBaseBias(double libraryAverageCoverage) {
			if (libraryAverageCoverage > 0.0D) {
				this.averagePerBaseBias = (this.averageCoverage / libraryAverageCoverage);
			} else {
				this.averagePerBaseBias = (0.0D / 0.0D);
			}
		}

		public void addFromAnother(TargetReadDepthResults targetReadDepthResults) {
			this.averageCoverage += targetReadDepthResults.getAverageCoverage();
			this.averageGCContent += targetReadDepthResults.getAverageGCContent();
			this.averagePerBaseBias += targetReadDepthResults.getAveragePerBaseBias();
			this.averagePhreadScore += targetReadDepthResults.getAvgPhreadScore();
			this.averageMapQ += targetReadDepthResults.getAvgMapQ();
			this.numBPSequencedInTarget += targetReadDepthResults.getNumBPSequencedInTarget();
			this.numGCSequencedInTarget += targetReadDepthResults.getNumGCSequencedInTarget();
			this.numReadsSequenced += targetReadDepthResults.getNumReadsSequenced();
			this.averageInsertSize += targetReadDepthResults.getAverageInsertSize();
			for (int i = 0; i < this.filterNGS.getReadDepthFilter().length; i++) {
				this.numBPSequencedAtDepth[i] += targetReadDepthResults.getnumBPSequencedAtDepth(i);
				this.numGCSequencedAtDepth[i] += targetReadDepthResults.getNumGCSequencedAtDepth(i);
				this.percentGCSequencedAtDepth[i] += targetReadDepthResults.getPercentGCSequencedAtDepth(i);
			}
		}

		public int getNumReadsSequenced() {
			return this.numReadsSequenced;
		}

		public int getNumBPInTarget() {
			return this.numBPInTarget;
		}

		public void parseCount(int depth, int numGCs) {
			addNumBPSequencedInTarget(depth);
			addGCSequencedInTarget(numGCs);
			for (int i = 0; i < this.filterNGS.getReadDepthFilter().length; i++) {
				if (depth >= this.filterNGS.getReadDepthFilter()[i]) {
					this.numBPSequencedAtDepth[i] += 1;
					this.numGCSequencedAtDepth[i] += numGCs;
					if (depth > 0) {
						this.percentGCSequencedAtDepth[i] += numGCs / depth;
					}
					if (numGCs > depth) {
						System.err.println("Error internal accounting error, number of GCs to big");
					}
				}
			}
		}
	}

	private static double[] parseBytesToDouble(byte[] bytes) {
		double[] d = new double[bytes.length];
		for (int i = 0; i < d.length; i++) {
			d[i] = bytes[i];
		}
		return d;
	}

	private static String[] parseToString(String s) {
		String[] sa = new String[s.length()];
		for (int i = 0; i < s.length(); i++) {
			sa[i] = s.charAt(i) + "";
		}
		return sa;
	}

	public static class BaitsLibrary implements Serializable {
		private static final long serialVersionUID = 1L;
		public static final String[] BAITS_HEADER = { "TargetID", "ProbeID", "Sequence", "Replication", "Strand", "Coordinates" };
		private String[] targetIDs;
		private String[] ProbeIDs;
		private Segment[] baits;
		private double[] gcContentOfBait;
		//private Logger log;

		public BaitsLibrary(String[] targetIDs, String[] probeIDs, Segment[] baits, double[] gcContentOfBait, Logger log) {
			this.targetIDs = targetIDs;
			this.ProbeIDs = probeIDs;
			this.baits = baits;
			this.gcContentOfBait = gcContentOfBait;
		//	this.log = log;
		}

		public void mapToLibrary(LibraryNGS libraryNGS, boolean baitsAsTarget) {
			int[] baitsPerTarget = new int[libraryNGS.getTargetSegments().length];
			double[] baitGCContent = new double[libraryNGS.getTargetSegments().length];
			for (int i = 0; i < this.baits.length; i++) {
				int[] libraryIndices = libraryNGS.indicesInLibrary(this.baits[i]);
				if (libraryIndices == null) {
					//this.log.reportError("Error - the bait " + this.baits[i].getUCSClocation() + " could not be mapped to the library");
				} else {
					for (int j = 0; j < libraryIndices.length; j++) {
						baitsPerTarget[libraryIndices[j]] += 1;
						baitGCContent[libraryIndices[j]] += getGcContentOfBait()[i];
					}
				}
			}
			for (int i = 0; i < baitGCContent.length; i++) {
				if (baitsPerTarget[i] > 0) {
					baitGCContent[i] /= baitsPerTarget[i];
				} else {
					baitGCContent[i] = (0.0D / 0.0D);
				}
			}
			libraryNGS.setBaitsPerTarget(baitsPerTarget);
			libraryNGS.setBaitGCContent(baitGCContent);
		}

		public static BaitsLibrary loadSerial(String filename) {
			return (BaitsLibrary) Files.readSerial(filename);
		}

		public void serialize(String filename) {
			Files.writeSerial(this, filename);
		}

		public String[] getTargetIDs() {
			return this.targetIDs;
		}

		public String[] getProbeIDs() {
			return this.ProbeIDs;
		}

		public Segment[] getBaits() {
			return this.baits;
		}

		public double[] getGcContentOfBait() {
			return this.gcContentOfBait;
		}

		public static BaitsLibrary loadBaitLibrary(String fullPathToBaitLibrary, Logger log) {
			if (Files.serializedVersionExists(null, fullPathToBaitLibrary)) {
				log.report(ext.getTime() + " Info - loading baits library from " + Files.getSerializedFileName(null, fullPathToBaitLibrary));
				return loadSerial(Files.getSerializedFileName(null, fullPathToBaitLibrary));
			}
			log.report(ext.getTime() + " Info - loading baits library from " + fullPathToBaitLibrary);
			ArrayList<String> tmpTargetIDs = new ArrayList<String>(9000000);
			ArrayList<String> tmpProbeIDs = new ArrayList<String>(9000000);
			ArrayList<Segment> tmpBaits = new ArrayList<Segment>(9000000);
			ArrayList<Double> tmpGCContent = new ArrayList<Double>(9000000);
			try {
				BufferedReader reader = Files.getAppropriateReader(fullPathToBaitLibrary);
				int[] indices = ext.indexFactors(reader.readLine().trim().split("\t"), BAITS_HEADER, true, false);
				if (Array.countIf(indices, -1) > 0) {
					log.reportError("Error - could not detect proper header in baits file " + fullPathToBaitLibrary);
					log.report("Header must be " + Array.toStr(BAITS_HEADER));
					reader.close();
					reader.close();
					return null;
				}
				while (reader.ready()) {
					String[] line = reader.readLine().split("\t", -1);
					line[indices[5]] = line[indices[5]].trim();
					tmpTargetIDs.add(line[indices[0]]);
					tmpProbeIDs.add(line[indices[1]]);
					try {
						tmpBaits.add(new Segment(line[indices[5]]));
						String[] seq = LibraryNGS.parseToString(line[indices[2]]);
						int Gs = Array.countIf(seq, "G");
						int Cs = Array.countIf(seq, "C");
						double GCcontent = (Gs + Cs) / seq.length;
						tmpGCContent.add(Double.valueOf(GCcontent));
					} catch (NumberFormatException numberFormatException) {
						log.reportError("Error - could not parse line " + Array.toStr(line) + ", skipping");
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + fullPathToBaitLibrary + "\" not found in current directory");
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + fullPathToBaitLibrary + "\"");
			}
			int numBaits = tmpBaits.size();
			String[] targetIDs = new String[numBaits];
			String[] probeIDs = new String[numBaits];
			Segment[] baits = new Segment[numBaits];
			double[] gcContentOfBait = new double[numBaits];

			int[] sortedOrder = Segment.quicksort((Segment[]) tmpBaits.toArray(new Segment[tmpBaits.size()]));
			for (int i = 0; i < sortedOrder.length; i++) {
				targetIDs[i] = ((String) tmpTargetIDs.get(sortedOrder[i]));
				probeIDs[i] = ((String) tmpProbeIDs.get(sortedOrder[i]));
				baits[i] = ((Segment) tmpBaits.get(sortedOrder[i]));
				gcContentOfBait[i] = ((Double) tmpGCContent.get(sortedOrder[i])).doubleValue();
			}
			BaitsLibrary baitsLibrary = new BaitsLibrary(targetIDs, probeIDs, baits, gcContentOfBait, log);
			baitsLibrary.serialize(Files.getSerializedFileName(null, fullPathToBaitLibrary));
			log.report(ext.getTime() + " Info - finished baits library from " + fullPathToBaitLibrary);

			return baitsLibrary;
		}
	}
}
