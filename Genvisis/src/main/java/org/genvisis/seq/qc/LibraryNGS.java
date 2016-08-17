package org.genvisis.seq.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;

import com.google.common.primitives.Doubles;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;


public class LibraryNGS implements Serializable {
  public static class BaitsLibrary implements Serializable {
    private static final long serialVersionUID = 1L;
    public static final String[] BAITS_HEADER =
        {"TargetID", "ProbeID", "Sequence", "Replication", "Strand", "Coordinates"};

    public static BaitsLibrary loadBaitLibrary(String fullPathToBaitLibrary, Logger log) {
      if (SerializedFiles.serializedVersionExists(null, fullPathToBaitLibrary)) {
        log.report(ext.getTime() + " Info - loading baits library from "
            + SerializedFiles.getSerializedFileName(null, fullPathToBaitLibrary));
        return loadSerial(SerializedFiles.getSerializedFileName(null, fullPathToBaitLibrary));
      }
      log.report(ext.getTime() + " Info - loading baits library from " + fullPathToBaitLibrary);
      ArrayList<String> tmpTargetIDs = new ArrayList<String>(9000000);
      ArrayList<String> tmpProbeIDs = new ArrayList<String>(9000000);
      ArrayList<Segment> tmpBaits = new ArrayList<Segment>(9000000);
      ArrayList<Double> tmpGCContent = new ArrayList<Double>(9000000);
      try {
        BufferedReader reader = Files.getAppropriateReader(fullPathToBaitLibrary);
        int[] indices =
            ext.indexFactors(reader.readLine().trim().split("\t"), BAITS_HEADER, true, false);
        if (Array.countIf(indices, -1) > 0) {
          log.reportError(
              "Error - could not detect proper header in baits file " + fullPathToBaitLibrary);
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
        log.reportError(
            "Error: file \"" + fullPathToBaitLibrary + "\" not found in current directory");
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + fullPathToBaitLibrary + "\"");
      }
      int numBaits = tmpBaits.size();
      String[] targetIDs = new String[numBaits];
      String[] probeIDs = new String[numBaits];
      Segment[] baits = new Segment[numBaits];
      double[] gcContentOfBait = new double[numBaits];

      int[] sortedOrder = Segment.quicksort(tmpBaits.toArray(new Segment[tmpBaits.size()]));
      for (int i = 0; i < sortedOrder.length; i++) {
        targetIDs[i] = (tmpTargetIDs.get(sortedOrder[i]));
        probeIDs[i] = (tmpProbeIDs.get(sortedOrder[i]));
        baits[i] = (tmpBaits.get(sortedOrder[i]));
        gcContentOfBait[i] = tmpGCContent.get(sortedOrder[i]).doubleValue();
      }
      BaitsLibrary baitsLibrary =
          new BaitsLibrary(targetIDs, probeIDs, baits, gcContentOfBait, log);
      baitsLibrary.serialize(SerializedFiles.getSerializedFileName(null, fullPathToBaitLibrary));
      log.report(ext.getTime() + " Info - finished baits library from " + fullPathToBaitLibrary);

      return baitsLibrary;
    }

    public static BaitsLibrary loadSerial(String filename) {
      return (BaitsLibrary) SerializedFiles.readSerial(filename);
    }

    private final String[] targetIDs;
    private final String[] ProbeIDs;

    private final Segment[] baits;

    private final double[] gcContentOfBait;
    // private Logger log;

    public BaitsLibrary(String[] targetIDs, String[] probeIDs, Segment[] baits,
        double[] gcContentOfBait, Logger log) {
      this.targetIDs = targetIDs;
      ProbeIDs = probeIDs;
      this.baits = baits;
      this.gcContentOfBait = gcContentOfBait;
      // this.log = log;
    }

    public Segment[] getBaits() {
      return baits;
    }

    public double[] getGcContentOfBait() {
      return gcContentOfBait;
    }

    public String[] getProbeIDs() {
      return ProbeIDs;
    }

    public String[] getTargetIDs() {
      return targetIDs;
    }

    public void mapToLibrary(LibraryNGS libraryNGS, boolean baitsAsTarget) {
      int[] baitsPerTarget = new int[libraryNGS.getTargetSegments().length];
      double[] baitGCContent = new double[libraryNGS.getTargetSegments().length];
      for (int i = 0; i < baits.length; i++) {
        int[] libraryIndices = libraryNGS.indicesInLibrary(baits[i]);
        if (libraryIndices == null) {
          // this.log.reportError("Error - the bait " + this.baits[i].getUCSClocation() + " could
          // not be mapped to the library");
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

    public void serialize(String filename) {
      SerializedFiles.writeSerial(this, filename);
    }
  }
  public static class LibraryReadDepthResults implements Serializable {
    public static final String[] SummaryHeader =
        {"UCSC", "numBasePairsTargeted", "numBasePairsSequenced", "numGCSequenced", "numReads",
            "averageCoverage", "averageGC", "averagePerBaseBias", "averagePhreadScore",
            "averageMapQ", "averageInsertSize", "numBaitsPerTarget", "AverageBaitGC"};
    public static final String[] ADD_HEADER =
        {"Percent_Covered_at_depth", "Percent_GC_at_depth", "num_GC_at_depth"};
    public static final String LIBRARY_READ_DEPTH_EXTENSION = ".libraryResults.ser";
    public static final String SUMMARY = ".libraryResults.summary";
    public static final String SUMMARY_BAITS = ".libraryBaitsResults.summary";
    private static final long serialVersionUID = 1L;

    public static LibraryReadDepthResults load(String filename, boolean jar) {
      return (LibraryReadDepthResults) SerializedFiles.readSerial(filename, jar, true);
    }

    private final LibraryNGS.TargetReadDepthResults[] targetReadDepthResults;
    private double[] totalPercentCoveredAtDepth;
    private double[] totalPercentGCAtDepth;
    private final Segment[] targetSegments;
    private final FilterNGS filterNGS;
    private int totalBasePairsTargeted;

    private int totalBasePairsSequenced;

    public LibraryReadDepthResults(LibraryNGS.ReadDepth readDepth, Segment[] targetSegments,
        FilterNGS filterNGS) {
      targetReadDepthResults = new LibraryNGS.TargetReadDepthResults[readDepth.getDepths().length];
      this.targetSegments = targetSegments;
      this.filterNGS = filterNGS;
      totalBasePairsTargeted = 0;
      totalBasePairsSequenced = 0;
    }

    public LibraryReadDepthResults(LibraryNGS.TargetReadDepthResults[] targetReadDepthResults,
        double[] totalPercentCoveredAtDepth, double[] totalPercentGCAtDepth,
        Segment[] targetSegments, FilterNGS filterNGS, int totalBasePairsTargeted,
        int totalBasePairsSequenced) {
      this.targetReadDepthResults = targetReadDepthResults;
      this.totalPercentCoveredAtDepth = totalPercentCoveredAtDepth;
      this.totalPercentGCAtDepth = totalPercentGCAtDepth;
      this.targetSegments = targetSegments;
      this.filterNGS = filterNGS;
      this.totalBasePairsTargeted = totalBasePairsTargeted;
      this.totalBasePairsSequenced = totalBasePairsSequenced;
    }

    public LibraryReadDepthResults(Segment[] targetSegments, FilterNGS filterNGS) {
      targetReadDepthResults = new LibraryNGS.TargetReadDepthResults[targetSegments.length];
      this.targetSegments = targetSegments;
      this.filterNGS = filterNGS;
      totalBasePairsTargeted = 0;
      totalBasePairsSequenced = 0;
    }

    public void dump(String filename, Logger log) {
      try {
        PrintWriter writer = new PrintWriter(new FileWriter(filename));
        writer.println(getHeader());
        for (int i = 0; i < targetReadDepthResults.length; i++) {
          writer.println(getSummaryFor(i));
        }
        writer.close();
      } catch (Exception e) {
        log.reportError("Error writing to " + filename);
        log.reportException(e);
      }
    }

    public String getHeader() {
      String header = "";
      header = header + Array.toStr(SummaryHeader);
      for (String element : ADD_HEADER) {
        for (int j = 0; j < filterNGS.getReadDepthFilter().length; j++) {
          header = header + "\t" + element + "_" + filterNGS.getReadDepthFilter()[j];
        }
      }
      return header;
    }

    public String getSummaryFor(int libraryIndex) {
      String summary = "";
      summary = summary + targetSegments[libraryIndex].getUCSClocation();
      summary = summary + "\t" + targetReadDepthResults[libraryIndex].getSummmary();
      return summary;
    }

    public LibraryNGS.TargetReadDepthResults[] getTargetReadDepthResults() {
      return targetReadDepthResults;
    }

    public int getTotalBasePairsTargeted() {
      return totalBasePairsTargeted;
    }

    public int getTotalNumReads() {
      int totalNumReads = 0;
      for (TargetReadDepthResults targetReadDepthResult : targetReadDepthResults) {
        totalNumReads += targetReadDepthResult.getNumReadsSequenced();
      }
      return totalNumReads;
    }

    public double[] getTotalPercentCoveredAtDepth() {
      return totalPercentCoveredAtDepth;
    }

    public double[] getTotalPercentGCAtDepth() {
      return totalPercentGCAtDepth;
    }

    public void populateSummary(LibraryNGS libraryNGS, LibraryNGS.ReadDepth readDepth,
        double normalizeFactor, FilterNGS filterNGS) {
      totalPercentCoveredAtDepth = new double[filterNGS.getReadDepthFilter().length];
      totalPercentGCAtDepth = new double[filterNGS.getReadDepthFilter().length];
      double avgCoverage = readDepth.getMedianCoverage();
      readDepth.computeAverageQualities();
      readDepth.normalize(normalizeFactor);
      for (int i = 0; i < readDepth.getDepths().length; i++) {
        totalBasePairsTargeted += readDepth.getDepths()[i].length;
        if ((libraryNGS.getBaitsPerTarget() != null) && (libraryNGS.getBaitGCContent() != null)) {
          targetReadDepthResults[i] = new LibraryNGS.TargetReadDepthResults(
              readDepth.getNumReads()[i], readDepth.getDepths()[i].length,
              libraryNGS.getBaitsPerTarget()[i], libraryNGS.getBaitGCContent()[i], filterNGS);
        } else {
          targetReadDepthResults[i] = new LibraryNGS.TargetReadDepthResults(
              readDepth.getNumReads()[i], readDepth.getDepths()[i].length, 0, 0.0D, filterNGS);
        }
        targetReadDepthResults[i].setNumSamples(1);
        targetReadDepthResults[i].setAverageCoverage(Array.mean(readDepth.getDepths()[i]));
        for (int j = 0; j < readDepth.getDepths()[i].length; j++) {
          targetReadDepthResults[i].parseCount(readDepth.getDepths()[i][j],
              readDepth.getGCs()[i][j]);
          totalBasePairsSequenced += readDepth.getDepths()[i][j];
        }
        if (targetReadDepthResults[i].getNumBPSequencedInTarget() > 0) {
          targetReadDepthResults[i]
              .setAverageGCContent(targetReadDepthResults[i].getNumGCSequencedInTarget()
                  / targetReadDepthResults[i].getNumBPSequencedInTarget());
        }
        targetReadDepthResults[i].computePerBaseBias(avgCoverage);
        targetReadDepthResults[i].setAvgPhreadScore(readDepth.getAvgPhreadAt(i));
        targetReadDepthResults[i].setAvgMapQ(readDepth.getAvgMapQAt(i));
        targetReadDepthResults[i].setAverageInsertSize(readDepth.getAvgInsertSizeAt(i));
      }
      summarizePercents(filterNGS);
    }

    public void serialize(String filename) {
      SerializedFiles.writeSerial(this, filename);
    }

    public void setTotalPercentCoveredAtDepth(double[] totalPercentCoveredAtDepth) {
      this.totalPercentCoveredAtDepth = totalPercentCoveredAtDepth;
    }

    public void setTotalPercentGCAtDepth(double[] totalPercentGCAtDepth) {
      this.totalPercentGCAtDepth = totalPercentGCAtDepth;
    }

    public void summarizePercents(FilterNGS filterNGS) {
      for (TargetReadDepthResults targetReadDepthResult : targetReadDepthResults) {
        targetReadDepthResult.computePercents();
        for (int j = 0; j < filterNGS.getReadDepthFilter().length; j++) {
          totalPercentCoveredAtDepth[j] += targetReadDepthResult.getnumBPSequencedAtDepth(j);
          totalPercentGCAtDepth[j] += targetReadDepthResult.getPercentGCSequencedAtDepth(j);
        }
      }
      for (int i = 0; i < filterNGS.getReadDepthFilter().length; i++) {
        totalPercentCoveredAtDepth[i] /= totalBasePairsTargeted;
        totalPercentGCAtDepth[i] /= totalBasePairsSequenced;
      }
    }
  }
  public static class ReadDepth {
    private static class BlockParser {
      private final AlignmentBlock alignmentBlock;
      private final Segment targetSegment;
      private int sequenceStartArrayIndex;
      private int targetStartArrayIndex;
      private int targetStopArrayIndex;
      private final Logger log;

      public BlockParser(AlignmentBlock alignmentBlock, Segment targetSegment, Logger log) {
        this.alignmentBlock = alignmentBlock;
        this.targetSegment = targetSegment;
        this.log = log;
      }

      public boolean alligns(int alignmentRefStart, int alignmentRefStop) {
        return new Segment(targetSegment.getChr(), alignmentRefStart, alignmentRefStop)
            .overlaps(targetSegment);
      }

      public int getSequenceStartArrayIndex() {
        return sequenceStartArrayIndex;
      }

      // public Segment getTargetSegment() {
      // return this.targetSegment;
      // }

      public int getTargetStartArrayIndex() {
        return targetStartArrayIndex;
      }

      public int getTargetStopArrayIndex() {
        return targetStopArrayIndex;
      }

      public void parseIndices() {
        int alignmentRefStart = alignmentBlock.getReferenceStart();
        int alignmentRefStop = alignmentRefStart + alignmentBlock.getLength() - 1;
        sequenceStartArrayIndex = (alignmentBlock.getReadStart() - 1);
        if (alligns(alignmentRefStart, alignmentRefStop)) {
          if ((alignmentRefStart >= targetSegment.getStart())
              && (alignmentRefStop <= targetSegment.getStop())) {
            targetStartArrayIndex = (alignmentRefStart - targetSegment.getStart());
            targetStopArrayIndex =
                (alignmentRefStart - targetSegment.getStart() + alignmentBlock.getLength());
          } else if ((targetSegment.getStart() >= alignmentRefStart)
              && (targetSegment.getStop() <= alignmentRefStop)) {
            targetStartArrayIndex = 0;
            targetStopArrayIndex = targetSegment.getSize();
          } else if ((alignmentRefStart >= targetSegment.getStart())
              && (alignmentRefStop >= targetSegment.getStop())) {
            targetStartArrayIndex = (alignmentRefStart - targetSegment.getStart());
            targetStopArrayIndex = targetSegment.getSize();
          } else if ((alignmentRefStart <= targetSegment.getStart())
              && (alignmentRefStop <= targetSegment.getStop())) {
            targetStartArrayIndex = 0;
            targetStopArrayIndex =
                (alignmentRefStart + alignmentBlock.getLength() - targetSegment.getStart());
          } else {
            log.reportError(
                "Internal Error - could not determine block alignments for " + alignmentRefStart
                    + "\t" + alignmentRefStop + "\t" + targetSegment.getUCSClocation());
            targetStartArrayIndex = -1;
            targetStopArrayIndex = -1;
          }
          if ((targetStartArrayIndex < 0) || (targetStopArrayIndex < 0)) {
            log.reportError(
                "Internal Error - could not determine block alignments for " + alignmentRefStart
                    + "\t" + alignmentRefStop + "\t" + targetSegment.getUCSClocation());
            targetStartArrayIndex = -1;
            targetStopArrayIndex = -1;
          }
        } else {
          targetStartArrayIndex = -2;
          targetStopArrayIndex = -2;
        }
      }
    }

    public static final String[] BASES = {"A", "T", "C", "G", "N", "Other"};
    public static final String[] GC = {"G", "C"};

    private static int[][] initInt(Segment[] targetSegments, Logger log) {
      int[][] depths = new int[targetSegments.length][];
      for (int i = 0; i < depths.length; i++) {
        depths[i] = new int[targetSegments[i].getSize()];
      }
      return depths;
    }

    private final int[][] depths;
    private final int[][] GCs;
    private final int[] numReads;
    private final double[] phreadsAvg;
    private final double[] mapQsAvg;

    private final double[] averageInsertSize;

    // private Histogram.DynamicHistogram[] dyHistograms;
    private final Logger log;

    public ReadDepth(LibraryNGS libraryNGS, FilterNGS filterNGS, Logger log) {
      GCs = initInt(libraryNGS.getTargetSegments(), log);
      depths = initInt(libraryNGS.getTargetSegments(), log);
      numReads = new int[libraryNGS.getTargetSegments().length];
      phreadsAvg = new double[libraryNGS.getTargetSegments().length];
      mapQsAvg = new double[libraryNGS.getTargetSegments().length];
      averageInsertSize = new double[libraryNGS.getTargetSegments().length];
      // this.dyHistograms = new Histogram.DynamicHistogram[filterNGS.getReadDepthFilter().length];
      this.log = log;
    }

    public void addCounts(int libraryIndex, Segment segment, SAMRecord samRecord,
        FilterNGS filterNGS) {
      List<AlignmentBlock> alList = samRecord.getAlignmentBlocks();
      double[] phreads = LibraryNGS.parseBytesToDouble(samRecord.getBaseQualities());
      String[] bases = LibraryNGS.parseToString(samRecord.getReadString());
      numReads[libraryIndex] += 1;
      if (samRecord.getProperPairFlag()) {
        averageInsertSize[libraryIndex] += Math.abs(samRecord.getInferredInsertSize());
      }
      if (samRecord.getMappingQuality() != 255) {
        mapQsAvg[libraryIndex] += samRecord.getMappingQuality();
      }
      for (int i = 0; i < alList.size(); i++) {
        AlignmentBlock alignmentBlock = alList.get(i);
        BlockParser blockParser = new BlockParser(alignmentBlock, segment, log);
        blockParser.parseIndices();
        if (blockParser.getTargetStartArrayIndex() >= 0) {
          int seqIndex = blockParser.getSequenceStartArrayIndex();
          for (int j = blockParser.getTargetStartArrayIndex(); j < blockParser
              .getTargetStopArrayIndex(); j++) {
            if ((filterNGS.getPhreadScoreFilter() == 0.0D)
                || (phreads[seqIndex] >= filterNGS.getPhreadScoreFilter())) {
              depths[libraryIndex][j] += 1;
              phreadsAvg[libraryIndex] += phreads[seqIndex];
              if ((bases[seqIndex].equals(GC[0])) || (bases[seqIndex].equals(GC[1]))) {
                GCs[libraryIndex][j] += 1;
              }
            }
            seqIndex++;
          }
        }
      }
    }

    public void computeAverageQualities() {
      for (int i = 0; i < depths.length; i++) {
        int numBases = Array.sum(depths[i]);
        if (numBases > 0) {
          phreadsAvg[i] /= numBases;
          mapQsAvg[i] /= numReads[i];
          averageInsertSize[i] /= numReads[i];
        }
      }
    }

    public double getAvgInsertSizeAt(int i) {
      return averageInsertSize[i];
    }

    public double getAvgMapQAt(int i) {
      return mapQsAvg[i];
    }

    public double getAvgPhreadAt(int i) {
      return phreadsAvg[i];
    }

    public LibraryNGS.LibraryReadDepthResults getDepthResults(LibraryNGS libraryNGS,
        FilterNGS filterNGS, double normalizeFactor) {
      LibraryNGS.LibraryReadDepthResults readDepthResults =
          new LibraryNGS.LibraryReadDepthResults(this, libraryNGS.getTargetSegments(), filterNGS);
      readDepthResults.populateSummary(libraryNGS, this, normalizeFactor, filterNGS);
      return readDepthResults;
    }

    public int[][] getDepths() {
      return depths;
    }

    public int[][] getGCs() {
      return GCs;
    }

    public Logger getLog() {
      return log;
    }

    public double getMedianCoverage() {
      double medianCoverage = (0.0D / 0.0D);
      int numBases = 0;
      ArrayList<Double> allCounts = new ArrayList<Double>(75000000);
      for (int[] depth : depths) {
        for (int j = 0; j < depth.length; j++) {
          numBases++;
          allCounts.add(Double.valueOf(depth[j]));
        }
      }
      if (numBases > 0) {
        medianCoverage = Array.median(Doubles.toArray(allCounts));
      }
      return medianCoverage;
    }

    public int[] getNumReads() {
      return numReads;
    }

    public void normalize(double normalizeFactor) {
      if ((normalizeFactor != 1.0D) && (normalizeFactor > 0.0D)) {
        for (int i = 0; i < depths.length; i++) {
          for (int j = 0; j < depths[i].length; j++) {
            depths[i][j] = ((int) Math.round(depths[i][j] * normalizeFactor));
            GCs[i][j] = ((int) Math.round(GCs[i][j] * normalizeFactor));
          }
        }
      }
    }
  }
  public static class TargetReadDepthResults implements Serializable {
    private static final long serialVersionUID = 1L;
    private final FilterNGS filterNGS;
    private int numBPSequencedInTarget;
    private final int numBPInTarget;
    private int numReadsSequenced;
    private final int numBaits;
    private int numSamples;
    private int numGCSequencedInTarget;
    private final double baitGC;
    private final int[] numBPSequencedAtDepth;
    private final int[] numGCSequencedAtDepth;
    private final double[] percentGCSequencedAtDepth;
    private final double[] percentAtDepth;
    private final double[] percentGCAtDepth;
    private final double[] averagedNumGCAtDepth;
    private double averageCoverage;
    private double averageGCContent;
    private double averagePerBaseBias;
    private double averagePhreadScore;
    private double averageMapQ;
    private double averageInsertSize;

    public TargetReadDepthResults(FilterNGS filterNGS, int numBPSequencedInTarget,
        int numBPInTarget, int numReadsSequenced, int numBaits, int numSamples,
        int numGCSequencedInTarget, double baitGC, int[] numBPSequencedAtDepth,
        int[] numGCSequencedAtDepth, double[] percentGCSequencedAtDepth, double[] percentAtDepth,
        double[] percentGCAtDepth, double[] averagedNumGCAtDepth, double averageCoverage,
        double averageGCContent, double averagePerBaseBias, double averagePhreadScore,
        double averageMapQ, double averageInsertSize) {
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

    public TargetReadDepthResults(int numReadsSequenced, int numBPInTarget, int numBaits,
        double baitGC, FilterNGS filterNGS) {
      this.filterNGS = filterNGS;
      numBPSequencedInTarget = 0;
      numGCSequencedInTarget = 0;
      averageCoverage = 0.0D;
      averageGCContent = 0.0D;
      averagePerBaseBias = 0.0D;
      averagePhreadScore = 0.0D;
      averageMapQ = 0.0D;
      averageInsertSize = 0.0D;
      this.numBPInTarget = numBPInTarget;
      this.numBaits = numBaits;
      this.baitGC = baitGC;
      numBPSequencedAtDepth = new int[filterNGS.getReadDepthFilter().length];
      numGCSequencedAtDepth = new int[filterNGS.getReadDepthFilter().length];
      percentGCSequencedAtDepth = new double[filterNGS.getReadDepthFilter().length];
      percentAtDepth = new double[filterNGS.getReadDepthFilter().length];
      percentGCAtDepth = new double[filterNGS.getReadDepthFilter().length];
      averagedNumGCAtDepth = new double[filterNGS.getReadDepthFilter().length];
      this.numReadsSequenced = numReadsSequenced;
    }

    public void addFromAnother(TargetReadDepthResults targetReadDepthResults) {
      averageCoverage += targetReadDepthResults.getAverageCoverage();
      averageGCContent += targetReadDepthResults.getAverageGCContent();
      averagePerBaseBias += targetReadDepthResults.getAveragePerBaseBias();
      averagePhreadScore += targetReadDepthResults.getAvgPhreadScore();
      averageMapQ += targetReadDepthResults.getAvgMapQ();
      numBPSequencedInTarget += targetReadDepthResults.getNumBPSequencedInTarget();
      numGCSequencedInTarget += targetReadDepthResults.getNumGCSequencedInTarget();
      numReadsSequenced += targetReadDepthResults.getNumReadsSequenced();
      averageInsertSize += targetReadDepthResults.getAverageInsertSize();
      for (int i = 0; i < filterNGS.getReadDepthFilter().length; i++) {
        numBPSequencedAtDepth[i] += targetReadDepthResults.getnumBPSequencedAtDepth(i);
        numGCSequencedAtDepth[i] += targetReadDepthResults.getNumGCSequencedAtDepth(i);
        percentGCSequencedAtDepth[i] += targetReadDepthResults.getPercentGCSequencedAtDepth(i);
      }
    }

    public void addGCSequencedInTarget(int num) {
      numGCSequencedInTarget += num;
    }

    public void addNumBPSequencedInTarget(int num) {
      numBPSequencedInTarget += num;
    }

    public void computePerBaseBias(double libraryAverageCoverage) {
      if (libraryAverageCoverage > 0.0D) {
        averagePerBaseBias = (averageCoverage / libraryAverageCoverage);
      } else {
        averagePerBaseBias = (0.0D / 0.0D);
      }
    }

    public void computePercents() {
      double divideDepthBy = numBPInTarget * numSamples;
      averageCoverage = (numBPInTarget > 0 ? averageCoverage / numSamples : 0.0D);
      averageGCContent = (numBPInTarget > 0 ? averageGCContent / numSamples : 0.0D);
      averagePerBaseBias = (numBPInTarget > 0 ? averagePerBaseBias / numSamples : 0.0D);
      averageMapQ = (numBPInTarget > 0 ? averageMapQ / numSamples : 0.0D);
      averagePhreadScore = (numBPInTarget > 0 ? averagePhreadScore / numSamples : 0.0D);
      averageInsertSize = (numBPInTarget > 0 ? averageInsertSize / numSamples : 0.0D);
      if (Double.isNaN(averageGCContent)) {
        System.out.println(numBPInTarget + "\t" + averageGCContent + "\t" + numSamples + "\t"
            + numBPSequencedInTarget);
      }
      for (int i = 0; i < percentAtDepth.length; i++) {
        percentAtDepth[i] = (numBPInTarget > 0 ? numBPSequencedAtDepth[i] / divideDepthBy : 0.0D);
        percentGCAtDepth[i] =
            (numBPInTarget > 0 ? percentGCSequencedAtDepth[i] / divideDepthBy : 0.0D);
        averagedNumGCAtDepth[i] =
            (numBPInTarget > 0 ? numGCSequencedAtDepth[i] / numSamples : 0.0D);
      }
    }

    public double getAverageCoverage() {
      return averageCoverage;
    }

    public double getAverageGCContent() {
      return averageGCContent;
    }

    public double getAverageInsertSize() {
      return averageInsertSize;
    }

    public double getAveragePerBaseBias() {
      return averagePerBaseBias;
    }

    public double getAvgMapQ() {
      return averageMapQ;
    }

    public double getAvgPhreadScore() {
      return averagePhreadScore;
    }

    public int getNumBPInTarget() {
      return numBPInTarget;
    }

    public int getnumBPSequencedAtDepth(int dIndex) {
      return numBPSequencedAtDepth[dIndex];
    }

    public int getNumBPSequencedInTarget() {
      return numBPSequencedInTarget;
    }

    public int getNumGCSequencedAtDepth(int dIndex) {
      return numGCSequencedAtDepth[dIndex];
    }

    public int getNumGCSequencedInTarget() {
      return numGCSequencedInTarget;
    }

    public int getNumReadsSequenced() {
      return numReadsSequenced;
    }

    public double getPercentGCSequencedAtDepth(int dIndex) {
      return percentGCSequencedAtDepth[dIndex];
    }

    public String getSummmary() {
      String summary = "";
      summary = summary + numBPInTarget;
      summary = summary + "\t" + numBPSequencedInTarget;
      summary = summary + "\t" + numGCSequencedInTarget;
      summary = summary + "\t" + numReadsSequenced;
      summary = summary + "\t" + averageCoverage;
      summary = summary + "\t" + averageGCContent;
      summary = summary + "\t" + averagePerBaseBias;
      summary = summary + "\t" + averagePhreadScore;
      summary = summary + "\t" + averageMapQ;
      summary = summary + "\t" + averageInsertSize;
      summary = summary + "\t" + numBaits;
      summary = summary + "\t" + baitGC;
      summary = summary + "\t" + Array.toStr(percentAtDepth);
      summary = summary + "\t" + Array.toStr(percentGCAtDepth);
      summary = summary + "\t" + Array.toStr(averagedNumGCAtDepth);
      return summary;
    }

    public void parseCount(int depth, int numGCs) {
      addNumBPSequencedInTarget(depth);
      addGCSequencedInTarget(numGCs);
      for (int i = 0; i < filterNGS.getReadDepthFilter().length; i++) {
        if (depth >= filterNGS.getReadDepthFilter()[i]) {
          numBPSequencedAtDepth[i] += 1;
          numGCSequencedAtDepth[i] += numGCs;
          if (depth > 0) {
            percentGCSequencedAtDepth[i] += numGCs / depth;
          }
          if (numGCs > depth) {
            System.err.println("Error internal accounting error, number of GCs to big");
          }
        }
      }
    }

    public void setAverageCoverage(double averageCoverage) {
      this.averageCoverage = averageCoverage;
    }

    public void setAverageGCContent(double averageGCContent) {
      this.averageGCContent = averageGCContent;
    }

    public void setAverageInsertSize(double averageInsertSize) {
      this.averageInsertSize = averageInsertSize;
    }

    public void setAvgMapQ(double avgMapQ) {
      averageMapQ = avgMapQ;
    }

    public void setAvgPhreadScore(double avgPhreadScore) {
      averagePhreadScore = avgPhreadScore;
    }

    public void setNumGCSequencedInTarget(int numGCSequencedInTarget) {
      this.numGCSequencedInTarget = numGCSequencedInTarget;
    }

    public void setNumSamples(int numSamples) {
      this.numSamples = numSamples;
    }
  }

  private static final long serialVersionUID = 1L;

  public static LibraryNGS getLibraryNGS(String libraryNGSFile, int skipNumLines, Logger log) {
    LibraryNGS libraryNGS = null;
    String serLibrary = ext.rootOf(libraryNGSFile, false) + ".ser";
    if (!Files.exists(serLibrary)) {
      libraryNGS = new LibraryNGS(
          Segment.loadRegions(libraryNGSFile, 0, 1, 2, skipNumLines, true, true, true, 0), log);
      libraryNGS.serialize(serLibrary);
    } else {
      libraryNGS = load(serLibrary, false);
    }
    return libraryNGS;
  }

  public static LibraryNGS load(String filename, boolean jar) {
    return (LibraryNGS) SerializedFiles.readSerial(filename, jar, true);
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

  public static void summarizeLibraries(LibraryNGS originalLibraryNGS,
      String[] libraryReadDepthResultFiles, String output, FilterNGS filterNGS, Logger log) {
    LibraryReadDepthResults summaryResults =
        new LibraryReadDepthResults(originalLibraryNGS.getTargetSegments(), filterNGS);
    summaryResults.setTotalPercentCoveredAtDepth(new double[filterNGS.getReadDepthFilter().length]);
    summaryResults.setTotalPercentGCAtDepth(new double[filterNGS.getReadDepthFilter().length]);
    for (String libraryReadDepthResultFile : libraryReadDepthResultFiles) {
      LibraryReadDepthResults curReadDepthResults =
          LibraryReadDepthResults.load(libraryReadDepthResultFile, false);
      for (int j = 0; j < summaryResults.getTargetReadDepthResults().length; j++) {
        if (summaryResults.getTargetReadDepthResults()[j] == null) {
          summaryResults.getTargetReadDepthResults()[j] = new TargetReadDepthResults(0,
              curReadDepthResults.getTargetReadDepthResults()[j].getNumBPInTarget(),
              originalLibraryNGS.getBaitsPerTarget()[j], originalLibraryNGS.getBaitGCContent()[j],
              filterNGS);
          summaryResults.getTargetReadDepthResults()[j]
              .setNumSamples(libraryReadDepthResultFiles.length);
        }
        summaryResults.getTargetReadDepthResults()[j]
            .addFromAnother(curReadDepthResults.getTargetReadDepthResults()[j]);
      }
    }
    summaryResults.summarizePercents(filterNGS);
    summaryResults.dump(output, log);
  }

  private Segment[] targetSegments;

  private int[] baitsPerTarget;

  private double[] baitGCContent;

  private final Logger log;

  public LibraryNGS(Segment[] targetSegments, Logger log) {
    this.targetSegments = targetSegments;
    this.log = log;
  }

  public double[] getBaitGCContent() {
    return baitGCContent;
  }

  public int[] getBaitsPerTarget() {
    return baitsPerTarget;
  }

  public Logger getLog() {
    return log;
  }

  public Segment[] getTargetSegments() {
    return targetSegments;
  }

  public Segment getTargetSegmentsAt(int libraryIndex) {
    return getTargetSegments()[libraryIndex];
  }

  public int[] indicesInLibrary(Segment segment) {
    return Segment.binarySearchForAllOverLappingIndices(segment, targetSegments);
  }

  public boolean inLibrary(Segment segment) {
    return indicesInLibrary(segment) != null;
  }

  public void mapBaits(BaitsLibrary baitsLibrary, boolean baitsAsTarget) {
    baitsLibrary.mapToLibrary(this, baitsAsTarget);
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  public void setBaitGCContent(double[] baitGCContent) {
    this.baitGCContent = baitGCContent;
  }

  public void setBaitsPerTarget(int[] baitsPerTarget) {
    this.baitsPerTarget = baitsPerTarget;
  }

  public void setTargetSegments(Segment[] targetSegments) {
    this.targetSegments = targetSegments;
  }
}
