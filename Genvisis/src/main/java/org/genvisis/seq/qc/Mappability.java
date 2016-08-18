package org.genvisis.seq.qc;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BEDFileReader;
import org.genvisis.seq.manage.BEDFileReader.BEDFeatureSeg;
import org.genvisis.seq.manage.BedOps;
import org.genvisis.stats.Histogram.DynamicAveragingHistogram;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.bed.BEDFeature;

public class Mappability<SEGMENT extends Segment> {

  private final LocusSet<SEGMENT> set;
  private ArrayList<MappabilityResult<SEGMENT>> mappabilityResults;
  private final String mappabilityFile;
  private final String callSubsetBed;
  private final Logger log;

  public Mappability(LocusSet<SEGMENT> set, String mappabilityFile, String callSubsetBed,
                     Logger log) {
    super();
    this.set = set;
    this.mappabilityFile = mappabilityFile;
    this.callSubsetBed = callSubsetBed;
    this.mappabilityResults = null;
    this.log = log;
  }

  public ArrayList<MappabilityResult<SEGMENT>> getMappabilityResults() {
    return mappabilityResults;
  }

  public void computeMappability() {
    if (BedOps.verifyBedIndex(mappabilityFile, log) && BedOps.verifyBedIndex(callSubsetBed, log)) {
      this.mappabilityResults = new ArrayList<Mappability.MappabilityResult<SEGMENT>>();
      BEDFileReader mapReader = new BEDFileReader(mappabilityFile, true);
      BEDFileReader callSubsetReader = new BEDFileReader(callSubsetBed, true);
      for (int i = 0; i < set.getLoci().length; i++) {
        if (i % 100 == 0) {
          log.reportTimeInfo(i + " of " + set.getLoci().length);
        }
        mappabilityResults.add(new MappabilityResult<SEGMENT>(mapReader, callSubsetReader,
                                                              set.getLoci()[i], log));
      }
      mapReader.close();
    }
  }

  public Hashtable<String, Integer> generateInternalGeneCounts() {
    Hashtable<String, Integer> geneCounts = new Hashtable<String, Integer>();
    for (int i = 0; i < mappabilityResults.size(); i++) {
      MappabilityResult<SEGMENT> cnMapp = mappabilityResults.get(i);
      for (int j = 0; j < cnMapp.getSubsetNames().length; j++) {
        String geneName = cnMapp.getSubsetNames()[j];
        if (!geneCounts.containsKey(geneName)) {
          geneCounts.put(geneName, 0);
        }
        geneCounts.put(geneName, geneCounts.get(geneName) + 1);
      }

    }
    return geneCounts;
  }

  private Hashtable<String, Integer> generateGeneCounts(LocusSet<GeneData> gLocusSet) {
    Hashtable<String, Integer> geneCounts = new Hashtable<String, Integer>();
    for (int i = 0; i < mappabilityResults.size(); i++) {
      MappabilityResult<SEGMENT> cnMapp = mappabilityResults.get(i);
      GeneData[] overlappingGenes = gLocusSet.getOverLappingLoci(cnMapp.getT());
      if (overlappingGenes != null) {
        for (GeneData overlappingGene : overlappingGenes) {
          String geneName = overlappingGene.getGeneName();
          if (!geneCounts.containsKey(geneName)) {
            geneCounts.put(geneName, 0);
          }
          geneCounts.put(geneName, geneCounts.get(geneName) + 1);
        }

      }
    }
    return geneCounts;
  }

  public static void computeCNVMappability(String mappabilityFile, String cnvFile,
                                           String geneTrackFile, String callSubsetBed, Logger log) {
    BedOps.verifyBedIndex(mappabilityFile, log);
    LocusSet<GeneData> gLocusSet = GeneTrack.load(geneTrackFile, false).convertToLocusSet(log);
    CNVariant[] cnvs = CNVariant.loadPlinkFile(cnvFile, false);
    LocusSet<CNVariant> cLocusSet = new LocusSet<CNVariant>(cnvs, true, log) {

      /**
       * 
       */
      private static final long serialVersionUID = 1L;

    };

    Mappability<CNVariant> cnMappability =
        new Mappability<CNVariant>(cLocusSet, mappabilityFile, callSubsetBed, log);
    cnMappability.computeMappability();
    Hashtable<String, Integer> geneCounts = cnMappability.generateGeneCounts(gLocusSet);
    String outputRoot = ext.rootOf(cnvFile, false) + ".mappability.summary";
    // ArrayList<GeomText> geneLabels = new ArrayList<GeomText>();
    String[] header1 =
        new String[] {"GENE_NAME", "CNV_MAP_SCORE", "NUM_TOTAL_CNVS", "STRING_FACTOR_COUNT"};
    String out1 = outputRoot + ".txt";
    summarize1(log, gLocusSet, cnMappability, geneCounts, header1, out1);
    RScatter rScatterBox =
        new RScatter(out1, out1 + ".rscript", ext.removeDirectoryInfo(out1), out1 + ".jpeg",
                     header1[1], new String[] {header1[2]}, SCATTER_TYPE.POINT, log);
    rScatterBox.setOverWriteExisting(true);
    rScatterBox.execute();

    int maxNumPer = 0;
    for (String gene : geneCounts.keySet()) {
      if (geneCounts.get(gene) > maxNumPer) {
        maxNumPer = geneCounts.get(gene);
      }
    }
    DynamicAveragingHistogram dynamicAveragingHistogramCNVCentered =
        new DynamicAveragingHistogram(0, maxNumPer, 0);
    DynamicAveragingHistogram dynamicAveragingHistogramMapCentered =
        new DynamicAveragingHistogram(0, 1, 2);

    for (int i = 0; i < cnMappability.getMappabilityResults().size(); i++) {
      MappabilityResult<CNVariant> cnMapp = cnMappability.getMappabilityResults().get(i);
      GeneData[] overlappingGenes = gLocusSet.getOverLappingLoci(cnMapp.getT());
      if (overlappingGenes == null) {
      } else {
        for (GeneData overlappingGene : overlappingGenes) {
          String curGene = overlappingGene.getGeneName();
          int counts = geneCounts.get(curGene);
          dynamicAveragingHistogramCNVCentered.addDataPair(counts, cnMapp.getAverageMapScore());
          dynamicAveragingHistogramMapCentered.addDataPair(cnMapp.getAverageMapScore(), counts);
        }
      }
    }
    dynamicAveragingHistogramCNVCentered.average();
    dynamicAveragingHistogramMapCentered.average();
    String[] header2 = new String[] {"CNVS_PER_GENE", "AVG_MAP_SCORE"};
    String out2 = outputRoot + ".hist.cnvCenter.txt";
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(out2));
      writer.println(Array.toStr(header2));
      for (int i = 0; i < dynamicAveragingHistogramCNVCentered.getAverages().length; i++) {
        if (dynamicAveragingHistogramCNVCentered.getCounts()[i] > 0) {
          System.out.println(dynamicAveragingHistogramCNVCentered.getAverages()[i]);
          writer.println(dynamicAveragingHistogramCNVCentered.getBins()[i] + "\t"
                         + dynamicAveragingHistogramCNVCentered.getAverages()[i]);
        }
      }

      writer.close();

    } catch (Exception e) {
      log.reportError("Error writing to " + out2);
      log.reportException(e);
    }

    RScatter rScatterAverageCNVCENTER =
        new RScatter(out2, out2 + ".rscript", ext.removeDirectoryInfo(out2), out2 + ".jpeg",
                     header2[0], new String[] {header2[1]}, SCATTER_TYPE.POINT, log);
    rScatterAverageCNVCENTER.setxLabel("CNVS Per Gene");
    rScatterAverageCNVCENTER.setyLabel("Average CNV mapping score");

    rScatterAverageCNVCENTER.setOverWriteExisting(true);
    rScatterAverageCNVCENTER.execute();

    String[] header3 = new String[] {"MAP_SCORE", "AVG_CNVS_PER_GENE"};
    String out3 = outputRoot + ".hist.scoreCenter.txt";

    try {
      PrintWriter writer = new PrintWriter(new FileWriter(out3));
      writer.println(Array.toStr(header3));
      for (int i = 0; i < dynamicAveragingHistogramMapCentered.getAverages().length; i++) {
        if (dynamicAveragingHistogramMapCentered.getCounts()[i] > 0) {
          writer.println(dynamicAveragingHistogramMapCentered.getBins()[i] + "\t"
                         + dynamicAveragingHistogramMapCentered.getAverages()[i]);
        }
      }

      writer.close();

    } catch (Exception e) {
      log.reportError("Error writing to " + out2);
      log.reportException(e);
    }
    RScatter rScatterAverageSCORECENTER =
        new RScatter(out3, out3 + ".rscript", ext.removeDirectoryInfo(out3), out3 + ".jpeg",
                     header3[0], new String[] {header3[1]}, SCATTER_TYPE.POINT, log);
    rScatterAverageSCORECENTER.setxLabel("CNV mapping score");
    rScatterAverageSCORECENTER.setyLabel("Average CNVs Per gene");

    rScatterAverageSCORECENTER.setOverWriteExisting(true);
    rScatterAverageSCORECENTER.execute();

    String[] header4 =
        Array.concatAll(CNVariant.PLINK_CNV_HEADER, new String[] {"GENE_NAME", "MAP_SCORE"});
    String out4 = outputRoot + ".cnvSummary.txt";
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(out4));
      writer.println(Array.toStr(header4));
      for (int i = 0; i < cnMappability.getMappabilityResults().size(); i++) {
        MappabilityResult<CNVariant> cnMapp = cnMappability.getMappabilityResults().get(i);
        GeneData[] overlappingGenes = gLocusSet.getOverLappingLoci(cnMapp.getT());
        if (overlappingGenes == null) {
          writer.println(cnMapp.getT().toAnalysisString() + "\tNA" + "\t"
                         + cnMapp.getAverageMapScore());
        } else {
          String genes = overlappingGenes[0].getGeneName();
          for (int j = 1; j < overlappingGenes.length; j++) {
            genes += ";" + overlappingGenes[j].getGeneName();

          }
          writer.println(cnMapp.getT().toAnalysisString() + "\t" + genes + "\t"
                         + cnMapp.getAverageMapScore());
        }
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + out4);
      log.reportException(e);
    }

    // Mappability<GeneData> gMappability = new Mappability<GeneData>(gLocusSet, mappabilityFile,
    // log);
    // gMappability.computeMappability();

  }

  private static void summarize1(Logger log, LocusSet<GeneData> gLocusSet,
                                 Mappability<CNVariant> cnMappability,
                                 Hashtable<String, Integer> geneCounts, String[] header1,
                                 String out1) {
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(out1));
      writer.println(Array.toStr(header1));

      for (int i = 0; i < cnMappability.getMappabilityResults().size(); i++) {
        MappabilityResult<CNVariant> cnMapp = cnMappability.getMappabilityResults().get(i);
        GeneData[] overlappingGenes = gLocusSet.getOverLappingLoci(cnMapp.getT());
        if (overlappingGenes == null) {
          log.reportTimeError("Could not find overlapping gene for "
                              + cnMapp.getT().toAnalysisString());
        } else {
          for (GeneData overlappingGene : overlappingGenes) {
            String curGene = overlappingGene.getGeneName();
            int counts = geneCounts.get(curGene);

            String stringCount = (counts < 10 ? "0" : "") + counts + "_CNVS";
            writer.println(curGene + "\t" + cnMapp.getAverageMapScore() + "\t" + counts + "\t"
                           + stringCount);
          }
        }
      }
      writer.close();

      // RScatter rScatterNum = new RScatter(out1, out1 + ".rscript", ext.removeDirectoryInfo(out1),
      // out1 + ".jpeg", header[2], new String[] { header[1] }, SCATTER_TYPE.POINT, log);
      // rScatterNum.setOverWriteExisting(true);
      // //rScatterNum.setgTexts(geneLabels.toArray(new GeomText[geneLabels.size()]));
      // rScatterNum.execute();
    } catch (Exception e) {
      log.reportError("Error writing to " + out1);
      log.reportException(e);
    }
  }

  public static class MappabilityResult<SEGMENT extends Segment> {
    private final Logger log;
    private double cumulativeMapScore;
    private int numBases;
    // private int calledOnCount = 0;
    private double averageMapScore;
    private String[] subsetNames;
    private final SEGMENT t;

    public double getAverageMapScore() {
      return averageMapScore;
    }

    public MappabilityResult(BEDFileReader mapReader, BEDFileReader callSubsetReader, SEGMENT t,
                             Logger log) {
      super();
      this.t = t;
      this.log = log;
      this.subsetNames = null;
      computeMappability(mapReader, callSubsetReader);
    }

    public SEGMENT getT() {
      return t;
    }

    public String[] getSubsetNames() {
      return subsetNames;
    }

    private void computeMappability(BEDFileReader reader, BEDFileReader callSubsetReader) {
      CloseableIterator<BEDFeature> iterator =
          reader.query(Positions.getChromosomeUCSC(t.getChr(), true), t.getStart(), t.getStop());
      // LocusSet<Segment> callSegs
      // =LocusSet.loadSegmentSetFromFile("C:/bin/ExomeDepth/exons.hg19.sort.bed", 0, 1, 2, 0, true,
      // true, 0, log);
      // callSegs.writeRegions(ext.addToRoot("C:/bin/ExomeDepth/exons.hg19.sort.bed", ".chr"),
      // TO_STRING_TYPE.REGULAR, false, log);
      LocusSet<BEDFeatureSeg> callSegs = callSubsetReader.loadSegsFor(t, log);

      //
      HashSet<String> subsetNamesAl = new HashSet<String>();
      this.numBases = 0;
      double currentScore = -1;
      // int calledOnCount = 0;
      while (iterator.hasNext()) {
        BEDFeature bedFeature = iterator.next();
        Segment bedSeg = BedOps.getSegment(bedFeature, log);
        double mapScore = Double.NaN;
        try {
          mapScore = Double.parseDouble(bedFeature.getName());

        } catch (NumberFormatException nfe) {
          String error = "Could not convert " + bedFeature.getName() + " to mappability score";
          log.reportTimeError(error);
          throw new IllegalArgumentException(error);
        }
        if (mapScore < 0) {
          String error = "Could not convert " + bedFeature.getName()
                         + " to mappability score, negative mapScore";
          log.reportTimeError(error);
          throw new IllegalArgumentException(error);
        }
        Segment union = bedSeg.getIntersection(t, log);

        BEDFeatureSeg[] called = callSegs.getOverLappingLoci(union);
        if (called != null) {
          for (BEDFeatureSeg element : called) {
            if (element.getBedFeature().getName() == null
                || element.getBedFeature().getName() == "") {
              String error =
                  "The call subset file is expected to have the gene name in the fourth column";
              log.reportTimeError(error);
              throw new IllegalArgumentException(error);
            }
            subsetNamesAl.add(element.getBedFeature().getName().split("_")[0]);

            double tmpMap = mapScore;
            Segment calledUnion = union.getIntersection(element, log);
            numBases += calledUnion.getSize();
            tmpMap *= calledUnion.getSize();
            // calledOnCount++;
            if (currentScore > 0) {
              currentScore += tmpMap;
            } else {
              currentScore = tmpMap;
            }
          }

        }
      }
      // this.calledOnCount = calledOnCount;
      this.cumulativeMapScore = currentScore;
      this.averageMapScore = cumulativeMapScore / numBases;
      this.subsetNames = subsetNamesAl.toArray(new String[subsetNamesAl.size()]);
      if (averageMapScore > 1) {
        System.out.println(averageMapScore + "\t" + cumulativeMapScore + "\t" + numBases);
        String error = "Detected an average mapping score greater than 1";
        log.reportTimeError(error);
        throw new IllegalArgumentException(error);
      }
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String mappabilityFile = "Mappability.bed";
    String cnvFile = "cnvs.cnv";
    String geneTrackFile = "RefSeq_hg19.gtrack";
    String callSubsetBed = "exons.hg19.bed";
    String usage = "\n" + "one.JL.Mappability requires 0-1 arguments\n";
    usage += "   (1) mappability file (i.e. mapFile=" + mappabilityFile + " (default))\n" + "";
    usage += "   (2) cnv file (i.e. cnvs=" + cnvFile + " (default))\n" + "";
    usage += "   (3) geneTrackFile  (i.e. genes=" + cnvFile + " (default))\n" + "";
    usage += "   (4) call subsetBed  (i.e. callSubset=" + callSubsetBed + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("mapFile=")) {
        mappabilityFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("cnvFile=")) {
        cnvFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("genes=")) {
        geneTrackFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("callSubset=")) {
        callSubsetBed = arg.split("=")[1];
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
      Logger log = new Logger(ext.rootOf(cnvFile, false) + ".mappability.log");
      computeCNVMappability(mappabilityFile, cnvFile, geneTrackFile, callSubsetBed, log);
    } catch (Exception e) {

      e.printStackTrace();
    }
  }

}
