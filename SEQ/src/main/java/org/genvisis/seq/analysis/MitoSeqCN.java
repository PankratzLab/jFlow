package org.genvisis.seq.analysis;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.concurrent.Callable;
import org.genvisis.cnv.LocusSet;
import org.genvisis.cnv.manage.ReferenceGenome;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.manage.BEDFileReader;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.BamOps.BamIndexStats;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerTrain;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;
import org.pankratzlab.core.CLI;
import org.pankratzlab.shared.filesys.Segment;
import org.genvisis.seq.manage.BedOps;
import org.genvisis.seq.manage.SamRecordOps;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

/**
 * @author lane0212 Inspired by Kendall
 */
public class MitoSeqCN {

  private MitoSeqCN() {

  }

  /**
   * @param fileOfBams bam files to estimate mtDNA CN for
   * @param outDir output directory for results
   * @param captureBed defining targeted capture regions
   * @param genomeBuild reference genome build
   * @param aName assembly name
   * @param aType assembly type
   * @param numthreads
   * @param log
   * @return the name of the output file
   */
  public static String run(String fileOfBams, String outDir, String captureBed,
                           String referenceGenomeFasta, ASSEMBLY_NAME aName, ASSAY_TYPE aType,
                           int numthreads, Logger log) {
    new File(outDir).mkdirs();

    String output = outDir + ext.rootOf(fileOfBams) + aType + "_" + aName + "_mtDNACN.summary.txt";

    if (!Files.exists(output)) {
      String[] bams = HashVec.loadFileToStringArray(fileOfBams, false, new int[] {0}, true);
      log.reportTimeInfo("Detected " + bams.length + " bam files");
      ReferenceGenome referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);

      LocusSet<Segment> genomeBinsMinusBinsCaputure = referenceGenome.getBins(20000).autosomal(true,
                                                                                               log);
      if (aType == ASSAY_TYPE.WXS) {
        if (captureBed == null || !Files.exists(captureBed)) {// required for WXS
          throw new IllegalArgumentException("A valid capture bed file must be provided for "
                                             + ASSAY_TYPE.WXS);
        }
        BedOps.verifyBedIndex(captureBed, log);
        BEDFileReader readerCapture = new BEDFileReader(captureBed, false);

        genomeBinsMinusBinsCaputure = genomeBinsMinusBinsCaputure.removeThese(readerCapture.loadAll(log)
                                                                                           .getStrictSegmentSet(),
                                                                              21000)
                                                                 .autosomal(true, log);
        readerCapture.close();

      } else {
        log.reportTimeWarning("Not using capture target subset for WGS ");

      }
      log.reportTimeInfo(genomeBinsMinusBinsCaputure.getBpCovered()
                         + " bp covered by reference bin regions");
      if (!referenceGenome.requestContig(aName.getMitoContig()).hasContig()
          || !referenceGenome.requestContig(aName.getxContig()).hasContig()
          || !referenceGenome.requestContig(aName.getyContig()).hasContig()) {
        throw new IllegalArgumentException("Required contig for " + aName + " is missing ( "
                                           + aName.getMitoContig() + " ," + aName.getxContig()
                                           + ", " + aName.getyContig() + " from "
                                           + referenceGenome.getReferenceFasta());
      } else {
        int mitoLength = referenceGenome.getContigLength(aName.getMitoContig());
        log.reportTimeInfo("Mitochondrial genome length = " + mitoLength);

        MitoCNProducer producer = new MitoCNProducer(bams, referenceGenome,
                                                     genomeBinsMinusBinsCaputure, outDir, aName,
                                                     log);

        ArrayList<MitoCNResult> results = new ArrayList<>();

        try (WorkerTrain<MitoCNResult> train = new WorkerTrain<>(producer, numthreads, numthreads,
                                                                 log);
             PrintWriter writer = Files.openAppropriateWriter(output)) {
          writer.println(ArrayUtils.toStr(MitoCNResult.header));
          while (train.hasNext()) {
            MitoCNResult result = train.next();
            if (result != null) {
              log.reportTimeInfo(ArrayUtils.toStr(result.getResult()));
              writer.println(ArrayUtils.toStr(result.getResult()));
              results.add(result);
            }

          }
        } catch (Exception e) {
          log.reportError("Error writing to " + output);
          log.reportException(e);
        }
      }
    } else {
      log.reportTimeWarning(output + " exists, skipping mtDNA CN estimation");
    }
    return output;
  }

  /**
   * Stores mtDNA CN estimation results for NGS data
   */
  public static class MitoCNResult {

    private static final String[] header = new String[] {"Sample", "NumMitoReads",
                                                         "TotalAlignedReads", "TotalUnAlignedReads",
                                                         "XReads", "YReads",
                                                         "AutosomalAlignedReads",
                                                         "NormalizationReads", "MitoLen", "normLen",
                                                         "MTBamFile", "MTBamFileTrim",
                                                         "EstimatedReadLength", "mtDNACNEstimate"};
    private final String sample;
    private final int numMitoReads;
    private final int numXReads;
    private final int numYReads;
    private int autosomalOnTargetReads;
    private final int normalizationReads;
    private final int mitoLen;
    private final int estimatedReadLength;
    private final long normLen;
    private final BamIndexStats bamIndexStats;
    private final String outBam;

    private MitoCNResult(String sample, int numMitoReads, int numXReads, int numYReads,
                         int offTargetReads, int mitoLen, long offTLen, BamIndexStats bamIndexStats,
                         int estimatedReadLength, String outBam) {
      super();
      this.sample = sample;
      this.numMitoReads = numMitoReads;
      this.numXReads = numXReads;
      this.numYReads = numYReads;
      this.mitoLen = mitoLen;
      this.normLen = offTLen;
      this.normalizationReads = offTargetReads;
      this.bamIndexStats = bamIndexStats;
      this.outBam = outBam;
      this.estimatedReadLength = estimatedReadLength;
      autosomalOnTargetReads = bamIndexStats.getAlignedRecordCount();
      autosomalOnTargetReads -= numMitoReads;
      autosomalOnTargetReads -= numXReads;
      autosomalOnTargetReads -= numYReads;

    }

    /**
     * @param numMitoReads number of mitochondrial reads counted
     * @param normalizationReads number of normalization reads to use
     * @param estimatedReadLength read length
     * @param normLen length of normalization region scanned
     * @param mitoLen length of mitochondrial genome
     * @return
     */
    private static double computemtDNACNEstimate(long numMitoReads, long normalizationReads,
                                                 long estimatedReadLength, long normLen,
                                                 long mitoLen) {
      double mtCov = (double) (numMitoReads * estimatedReadLength) / mitoLen;
      double autoCov = (double) (normalizationReads * estimatedReadLength) / normLen;

      return (2 * mtCov) / autoCov;

    }

    private String[] getResult() {

      ArrayList<String> result = new ArrayList<>();
      result.add(sample);
      result.add(Integer.toString(numMitoReads));
      result.add(Integer.toString(bamIndexStats.getAlignedRecordCount()));
      result.add(Integer.toString(bamIndexStats.getUnalignedRecordCount()));
      result.add(Integer.toString(numXReads));
      result.add(Integer.toString(numYReads));
      result.add(Integer.toString(autosomalOnTargetReads));
      result.add(Integer.toString(normalizationReads));
      result.add(Integer.toString(mitoLen));
      result.add(Long.toString(normLen));
      result.add(outBam);
      result.add(ext.rootOf(ext.rootOf(outBam)));
      result.add(Integer.toString(estimatedReadLength));
      result.add(Double.toString(computemtDNACNEstimate(numMitoReads, normalizationReads,
                                                        estimatedReadLength, normLen, mitoLen)));

      return ArrayUtils.toStringArray(result);
    }
  }

  private static class MitoCNWorker implements Callable<MitoCNResult> {

    private final String bam;
    private final String outDir;
    private final int mitoLength;
    private final int xLength;
    private final int yLength;
    private final LocusSet<Segment> genomeBinsMinusBinsCaputure;
    private final ASSEMBLY_NAME params;
    private final Logger log;

    private MitoCNWorker(String bam, LocusSet<Segment> genomeBinsMinusBinsCaputure, String outDir,
                         int mitoLength, int xLength, int yLength, ASSEMBLY_NAME params,
                         Logger log) {
      super();
      this.bam = bam;
      this.genomeBinsMinusBinsCaputure = genomeBinsMinusBinsCaputure;
      this.outDir = outDir;
      this.mitoLength = mitoLength;
      this.xLength = xLength;
      this.yLength = yLength;
      this.params = params;
      this.log = log;

    }

    private static boolean passesFilter(SAMRecord samRecord) {
      return !samRecord.getReadUnmappedFlag() && !samRecord.getDuplicateReadFlag();
    }

    @Override
    public MitoCNResult call() throws Exception {

      try {
        String sample = BamOps.getSampleName(bam, log);
        log.reportTimeInfo("Processing sample " + sample);
        String outputMTBam = outDir + ext.addToRoot(ext.removeDirectoryInfo(bam), ".chrM");

        SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.LENIENT);

        SAMFileWriter sAMFileWriter = new SAMFileWriterFactory().setCreateIndex(true)
                                                                .makeSAMOrBAMWriter(reader.getFileHeader(),
                                                                                    true,
                                                                                    new File(outputMTBam));

        ArrayList<Segment> toSearch = new ArrayList<>();
        toSearch.add(new Segment(params.getMitoContig(), 0, mitoLength + 1));
        toSearch.add(new Segment(params.getxContig(), 0, xLength + 1));
        toSearch.add(new Segment(params.getyContig(), 0, yLength + 1));
        for (Segment segment : toSearch) {
          log.reportTimeInfo("Will search : " + segment.getUCSClocation());
        }
        QueryInterval[] queryInterestIntervals = BamOps.convertSegsToQI(toSearch.toArray(new Segment[toSearch.size()]),
                                                                        reader.getFileHeader(), 0,
                                                                        true, params.addChr(), log);
        SAMRecordIterator sIterator = reader.query(queryInterestIntervals, false);
        int numMitoReads = 0;
        int numXReads = 0;
        int numYReads = 0;
        int numOffTarget = 0;
        while (sIterator.hasNext()) {
          SAMRecord samRecord = sIterator.next();

          if (passesFilter(samRecord)) {
            if (samRecord.getContig().equals(params.getMitoContig())) {
              sAMFileWriter.addAlignment(samRecord);
              numMitoReads++;
            } else if (samRecord.getContig().equals(params.getxContig())) {
              numXReads++;
            } else if (samRecord.getContig().equals(params.getyContig())) {
              numYReads++;
            } else {
              reader.close();
              throw new IllegalArgumentException("Invalid contig " + samRecord.getContig());
            }
          }
        }
        sIterator.close();
        sAMFileWriter.close();

        QueryInterval[] offTargetIntervalse = BamOps.convertSegsToQI(genomeBinsMinusBinsCaputure.getLoci(),
                                                                     reader.getFileHeader(), 0,
                                                                     true, params.addChr(), log);
        sIterator = reader.query(offTargetIntervalse, false);
        while (sIterator.hasNext()) {
          SAMRecord samRecord = sIterator.next();
          if (passesFilter(samRecord)) {
            numOffTarget++;

            if (numOffTarget % 1000000 == 0) {
              log.reportTimeInfo("Processing normalization-reads for sample " + sample + " , found "
                                 + numOffTarget + ", currently on "
                                 + SamRecordOps.getDisplayLoc(samRecord));
            }
          }

        }

        BamIndexStats bamIndexStats = BamOps.getBamIndexStats(reader);
        reader.close();
        int estimatedReadLength = (int) BamOps.estimateReadSize(bam, log).mean();
        return new MitoCNResult(sample, numMitoReads, numXReads, numYReads, numOffTarget,
                                mitoLength, genomeBinsMinusBinsCaputure.getBpCovered(),
                                bamIndexStats, estimatedReadLength, outputMTBam);
      } catch (Exception e) {
        log.reportError("Could not process " + bam);
        log.reportException(e);
        return null;
      }
    }
  }

  private static class MitoCNProducer extends AbstractProducer<MitoCNResult> {

    private final String[] bams;
    private final String outDir;
    private int index;
    private final int mitoLength;
    private final int xLength;
    private final int yLength;
    private final LocusSet<Segment> genomeBinsMinusBinsCaputure;
    private final ASSEMBLY_NAME params;
    private final Logger log;

    private MitoCNProducer(String[] bams, ReferenceGenome referenceGenome,
                           LocusSet<Segment> genomeBinsMinusBinsCaputure, String outDir,
                           ASSEMBLY_NAME params, Logger log) {
      super();
      this.bams = bams;
      this.outDir = outDir;
      mitoLength = referenceGenome.getContigLength(params.getMitoContig());
      xLength = referenceGenome.getContigLength(params.getxContig());
      yLength = referenceGenome.getContigLength(params.getyContig());
      this.genomeBinsMinusBinsCaputure = genomeBinsMinusBinsCaputure;
      index = 0;
      this.params = params;
      this.log = log;
    }

    @Override
    public boolean hasNext() {
      return index < bams.length;
    }

    @Override
    public Callable<MitoCNResult> next() {
      String currentBam = bams[index];

      index++;
      return new MitoCNWorker(currentBam, genomeBinsMinusBinsCaputure, outDir, mitoLength, xLength,
                              yLength, params, log);
    }
  }

  public static void main(String[] args) {

    CLI c = new CLI(MitoSeqCN.class);
    c.addArg("bams", "file of .bams with one (full path) .bam file per line ", true);
    c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "mitoWES/");
    c.addArgWithDefault(CLI.ARG_REFERENCE_GENOME,
                        CLI.DESC_REFERENCE_GENOME + " . The .dict file (such as hg19.dict) must also exist and the contigs must match the .bams analyzed",
                        "hg19.fa");
    c.addArgWithDefault("assayType",
                        "The assay type being analyzed (Options are "
                                     + ArrayUtils.toStr(ASSAY_TYPE.values(), ","),
                        ASSAY_TYPE.WXS.toString());
    c.addArgWithDefault("assemblyName",
                        "The assembly name being analyzed (Options are "
                                        + ArrayUtils.toStr(ASSEMBLY_NAME.values(), ",") + ")",
                        ASSEMBLY_NAME.HG19.toString());

    c.addArg("captureBed", "If using " + ASSAY_TYPE.WXS + ", a capture bed is required",
             "AgilentCaptureRegions.txt");
    c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, "8");

    c.parseWithExit(args);

    ASSAY_TYPE aType = ASSAY_TYPE.valueOf(c.get("assayType"));

    run(c.get("bams"), c.get(CLI.ARG_OUTDIR), c.get("captureBed"), c.get(CLI.ARG_REFERENCE_GENOME),
        ASSEMBLY_NAME.valueOf(c.get("assemblyName")), aType, c.getI(CLI.ARG_THREADS), new Logger());

  }
}
