package org.genvisis.seq.manage;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;

import com.google.common.primitives.Doubles;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReader.Indexing;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SamReaderFactory.Option;
import htsjdk.samtools.ValidationStringency;

/**
 * Class for common bamFile manips
 *
 */
public class BamOps {

  public static final String BAM_EXT = ".bam";
  public static final String BAI_EXT = ".bai";
  /**
   * Default number of reads to scan in order to estimate read length for a sample
   */
  public static final int NUM_READ_ESTIMATOR = 100000;

  private BamOps() {

  }

  /**
   * This method will check if an appropriate .bai index file exists for a given .bam, and create it
   * if not
   *
   * @param bamFile bam file to verify
   * @param log
   * @return
   */
  public static boolean verifyIndex(String bamFile, Logger log) {
    ArrayList<Option> options = new ArrayList<SamReaderFactory.Option>();
    options.add(Option.INCLUDE_SOURCE_IN_RECORDS);
    String index = getAssociatedBamIndex(bamFile);
    if (!Files.exists(index)) {
      log.reportTimeInfo("Attempting to generate index " + index);
      htsjdk.samtools.BAMIndexer.createIndex(BamOps.getDefaultReader(bamFile,
                                                                     ValidationStringency.STRICT,
                                                                     options),
                                             new File(index));
    }
    return Files.exists(index);
  }

  /**
   * @param bamFile
   * @return the associated
   */
  public static String getAssociatedBamIndex(String bamFile) {
    return ext.rootOf(bamFile, false) + BAI_EXT;
  }

  /**
   * @param bamOrSam .bam or .sam file
   * @param stringency Stringency validation for the records
   * @return new reader
   */
  public static SamReader getDefaultReader(String bamOrSam, ValidationStringency stringency) {

    return getDefaultReader(bamOrSam, stringency, new ArrayList<SamReaderFactory.Option>());
  }

  /**
   * @param bamFile the input bam file
   * @param outputBam the output bam file
   * @param set the {@link LocusSet} that will serve as an inclusion or exclusion filter
   * @param include if true, reads that overlap segments in the set will be included in the output.
   *        If false, reads that do not overlap the set will be included in the output
   * @param log
   * @return
   */
  public static boolean subsetBam(String bamFile, String outputBam, LocusSet<Segment> set,
                                  boolean include, Logger log) {
    SamReader reader = getDefaultReader(bamFile, ValidationStringency.STRICT);
    SAMFileWriter sAMFileWriter =
                                new SAMFileWriterFactory().setCreateIndex(true)
                                                          .makeSAMOrBAMWriter(reader.getFileHeader(),
                                                                              true, new File(outputBam));

    for (SAMRecord samRecord : reader) {
      Segment seg = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
      int[] indices = set.getOverlappingIndices(seg);
      if ((indices == null || indices.length == 0) && !include) {
        sAMFileWriter.addAlignment(samRecord);
      } else if ((indices != null && indices.length > 0) && include) {
        sAMFileWriter.addAlignment(samRecord);
      } else {
        throw new IllegalStateException("Misguided logic introduced");
      }
    }
    try {
      reader.close();
      sAMFileWriter.close();
    } catch (IOException e) {
      log.reportException(e);
      return false;
    }
    return true;
  }

  /**
   * @param bamOrSam .bam or .sam file
   * @param stringency Stringency validation for the records
   * @param options {@link SamReaderFactory.Option} to apply to the reader
   * @return new reader
   */
  public static SamReader getDefaultReader(String bamOrSam, ValidationStringency stringency,
                                           List<Option> options) {
    SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
    samReaderFactory.validationStringency(stringency);
    for (Option option : options) {
      samReaderFactory.enable(option);
    }
    return samReaderFactory.open(new File(bamOrSam));
  }

  /**
   * @param segs Genvisis type segments to search
   * @param sFileHeader an {@link SAMFileHeader}
   * @param bpBuffer bp buffer to be added to the segments
   * @param optimize perform an extra step to make a more efficient query
   * @param chr When querying, add "chr" to the contig
   * @param log
   * @return array of {@link QueryInterval} that can be queried by a bamfile reader
   */
  public static QueryInterval[] convertSegsToQI(Segment[] segs, SAMFileHeader sFileHeader,
                                                int bpBuffer, boolean optimize, boolean chr,
                                                Logger log) {
    QueryInterval[] qIntervals = new QueryInterval[segs.length];
    segs = Segment.sortSegments(segs);
    for (int i = 0; i < qIntervals.length; i++) {
      String sequenceName = Positions.getChromosomeUCSC(segs[i].getChr(), chr);
      int referenceIndex = sFileHeader.getSequenceIndex(sequenceName);
      if (referenceIndex < 0) {
        referenceIndex = sFileHeader.getSequenceIndex(sequenceName + "T");// MT
        if (referenceIndex < 0) {
          log.reportError("Error - could not find " + sequenceName
                          + " in the sequence dictionary, halting");
          return null;
        }
      }
      qIntervals[i] = new QueryInterval(referenceIndex, segs[i].getStart() - bpBuffer,
                                        segs[i].getStop() + bpBuffer);
    }
    if (optimize) {
      qIntervals = QueryInterval.optimizeIntervals(qIntervals);
    }
    return qIntervals;
  }

  public static Segment[] converQItoSegs(QueryInterval[] qIntervals, SAMFileHeader sFileHeader,
                                         Logger log) {
    Segment[] segs = new Segment[qIntervals.length];
    for (int i = 0; i < segs.length; i++) {
      segs[i] = new Segment(
                            Positions.chromosomeNumber(sFileHeader.getSequence(qIntervals[i].referenceIndex)
                                                                  .getSequenceName()),
                            qIntervals[i].start, qIntervals[i].end);
    }
    return segs;
  }

  /**
   * @param bamfile get the {@link SAMFileHeader} for this bam file
   * @param log
   * @return
   */
  public static SAMFileHeader getHeader(String bamfile, Logger log) {
    SamReader reader = getDefaultReader(bamfile, ValidationStringency.STRICT);
    SAMFileHeader samFileHeader = reader.getFileHeader();
    try {
      reader.close();
    } catch (IOException e) {
      log.reportException(e);
    }
    return samFileHeader;

  }

  /**
   * @param bams get all bar codes from the bams in this array
   * @param log
   * @return
   */
  public static String[] getAllBarCodes(String[] bams, Logger log) {
    HashSet<String> unique = new HashSet<String>();
    for (String bam : bams) {
      unique.addAll(getBarcodesFor(bam, log));
    }
    return unique.toArray(new String[unique.size()]);
  }

  /**
   * @param bam Get a bar code for this specific bam file
   * @param log
   * @return
   */
  public static ArrayList<String> getBarcodesFor(String bam, Logger log) {
    ArrayList<String> barcodes = new ArrayList<String>();
    SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
    samReaderFactory.validationStringency(ValidationStringency.LENIENT);
    SamReader reader = samReaderFactory.open(new File(bam));
    List<SAMReadGroupRecord> rgs = reader.getFileHeader().getReadGroups();
    HashSet<String> barcodesUnique = new HashSet<String>();
    for (SAMReadGroupRecord samReadGroupRecord : rgs) {
      String[] id = samReadGroupRecord.getId().split("_");
      String[] tmpCodes = id[id.length - 3].split("-");
      if (tmpCodes.length != 2) {
        throw new IllegalArgumentException("Could not parse barcodes for RG " + samReadGroupRecord
                                           + " in bam file " + bam);

      } else {
        for (String tmpCode : tmpCodes) {
          if (tmpCode.replaceAll("A", "").replaceAll("C", "").replaceAll("T", "")
                     .replaceAll("G", "").length() != 0) {
            throw new IllegalArgumentException("Invalid barcode " + tmpCode);
          } else {
            barcodesUnique.add(tmpCode);
          }
        }
      }

    }
    barcodes.addAll(barcodesUnique);
    return barcodes;
  }

  /**
   * @param bamFiles
   * @param log
   * @return sample names for all the bam files of interest
   */
  public static String[] getSampleNames(String[] bamFiles, Logger log) {
    String[] sampleNames = new String[bamFiles.length];
    for (int i = 0; i < sampleNames.length; i++) {
      sampleNames[i] = getSampleName(bamFiles[i], log);
    }
    return sampleNames;
  }

  /**
   * @author lane0212 Stores some simple counts that can be quickly retrieved from the index file
   */
  public static class BamIndexStats {
    private final int alignedRecordCount;
    private final int unalignedRecordCount;

    public BamIndexStats(int alignedRecordCount, int unalignedRecordCount) {
      super();
      this.alignedRecordCount = alignedRecordCount;
      this.unalignedRecordCount = unalignedRecordCount;
    }

    public int getAlignedRecordCount() {
      return alignedRecordCount;
    }

    public int getUnalignedRecordCount() {
      return unalignedRecordCount;
    }
  }

  public static BamIndexStats getBamIndexStats(String bamFile) {
    SamReader reader = BamOps.getDefaultReader(bamFile, ValidationStringency.STRICT);
    return getBamIndexStats(reader);
  }

  public static BamIndexStats getBamIndexStats(SamReader reader) {
    BAMIndexMetaData[] result = getIndexMetaData(reader);
    int alignedRecordCount = 0;
    int unalignedRecordCount = 0;
    for (BAMIndexMetaData element : result) {
      alignedRecordCount += element.getAlignedRecordCount();
      unalignedRecordCount += element.getUnalignedRecordCount();
    }

    return new BamIndexStats(alignedRecordCount, unalignedRecordCount);

  }

  public static BAMIndexMetaData[] getIndexMetaData(SamReader reader) {
    Indexing index = reader.indexing();
    BAMIndex bamIndex = index.getIndex();
    List<SAMSequenceRecord> records = reader.getFileHeader().getSequenceDictionary().getSequences();
    BAMIndexMetaData[] result = new BAMIndexMetaData[records.size()];
    for (int i = 0; i < result.length; i++) {
      result[i] = bamIndex.getMetaData(i);
    }
    return result;

  }

  /**
   * @param bamFile estimate the length of reads comprising this bam file
   * @param log
   * @return the estimated read size
   */
  public static int estimateReadSize(String bamFile, Logger log) {
    return estimateReadSize(bamFile, NUM_READ_ESTIMATOR, log);
  }



  /**
   * @author Kitty Stores information regarding the insert size estimate
   */
  public static class InsertSizeEstimate {
    private double avgInsertSize;
    private double stDevInsertSize;

    /**
     * @param avgInsertSize
     * @param stDevInsertSize
     */
    private InsertSizeEstimate(double avgInsertSize, double stDevInsertSize) {
      super();
      this.avgInsertSize = avgInsertSize;
      this.stDevInsertSize = stDevInsertSize;
    }

    public double getAvgInsertSize() {
      return avgInsertSize;
    }

    public double getStDevInsertSize() {
      return stDevInsertSize;
    }

    public String getSummary() {
      return "AvgInsertSize\t" + avgInsertSize + System.lineSeparator() + "StDevInsertSize\t"
             + stDevInsertSize;

    }

    public static InsertSizeEstimate load(String file) {
      String[] data = HashVec.loadFileToStringArray(file, false, new int[] {1}, false);
      return new InsertSizeEstimate(Double.parseDouble(data[0]), Double.parseDouble(data[1]));
    }



  }

  /**
   * @param bamFile estimate the insert size of reads comprising this bam file
   * @param numReads the number of reads to compute the estimate from (reasonable performance with
   *        100000)
   * @param log
   * @return the estimated insert size
   */
  public static InsertSizeEstimate estimateInsertSize(String bamFile, int numReads, Logger log) {
    SamReader reader = getDefaultReader(bamFile, ValidationStringency.STRICT);

    SAMRecordIterator iterator = reader.iterator();
    ArrayList<Double> insertSizes = new ArrayList<Double>();
    int readsScanned = 0;
    while (iterator.hasNext()) {
      SAMRecord samRecord = iterator.next();
      if (samRecord.getProperPairFlag() && !samRecord.getReadUnmappedFlag()
          && samRecord.getCigar().getCigarElements().size() == 1) {
        insertSizes.add((double) samRecord.getInferredInsertSize());
        readsScanned++;

      }
      if (readsScanned > numReads) {
        break;
      }
    }
    try {
      reader.close();
    } catch (IOException e) {
      log.reportException(e);

    }
    if (readsScanned > 0) {
      double[] finals = Doubles.toArray(insertSizes);
      double averageInsertSize = Array.mean(finals);
      double stDevInsertSize = Array.stdev(finals);
      return new InsertSizeEstimate(averageInsertSize, stDevInsertSize);

    }
    return new InsertSizeEstimate(0, 0);
  }

  /**
   * @param bamFile estimate the length of reads comprising this bam file
   * @param numReads the number of reads to compute the estimate from (reasonable performance with
   *        100000)
   * @param log
   * @return the estimated read size
   */
  public static int estimateReadSize(String bamFile, int numReads, Logger log) {
    SamReader reader = getDefaultReader(bamFile, ValidationStringency.STRICT);
    int readsize = 0;
    SAMRecordIterator iterator = reader.iterator();

    int readsScanned = 0;
    while (iterator.hasNext()) {
      SAMRecord samRecord = iterator.next();
      if (!samRecord.getReadUnmappedFlag() && samRecord.getCigar().getCigarElements().size() == 1) {

        readsize += samRecord.getReadLength();
        readsScanned++;

      }
      if (readsScanned > numReads) {
        break;
      }
    }
    try {
      reader.close();
    } catch (IOException e) {
      log.reportException(e);

    }
    if (readsScanned > 0) {
      double avg = (double) readsize / readsScanned;
      return (int) avg;
    }
    return -1;
  }

  @Deprecated
  public static String getSampleName(String bamFile) {
    return getSampleName(bamFile, new Logger());
  }

  /**
   * @param bamFile the bam file to extract sample name from
   * @param log your friendly logger
   * @return
   */
  public static String getSampleName(String bamFile, Logger log) {
    SamReader reader = getDefaultReader(bamFile, ValidationStringency.STRICT);
    String sample = reader.getFileHeader().getReadGroups().get(0).getSample();
    try {
      reader.close();
    } catch (IOException e) {
      log.reportException(e);
    }
    return sample;
  }

  private static class SampleNameExtractor implements Callable<SampleNameExtractor> {
    private final String bamFile;
    private String sampleName;

    public SampleNameExtractor(String bamFile) {
      super();
      this.bamFile = bamFile;
    }

    @Override
    public SampleNameExtractor call() throws Exception {
      sampleName = getSampleName(bamFile);
      return this;
    }

    public String getBamFile() {
      return bamFile;
    }

    public String getName() {
      return sampleName;
    }

  }

  private static class SampleNameProducer extends AbstractProducer<SampleNameExtractor> {

    private final String[] bamFiles;
    private int index;

    public SampleNameProducer(String[] bamFiles) {
      super();
      this.bamFiles = bamFiles;
      index = 0;
    }

    @Override
    public boolean hasNext() {
      return index < bamFiles.length;
    }

    @Override
    public Callable<SampleNameExtractor> next() {
      SampleNameExtractor ex = new SampleNameExtractor(bamFiles[index]);
      index++;
      return ex;
    }
  }

  /**
   * Designed to take the sample names from
   * {@link VCFOps#getSamplesInFile(htsjdk.variant.vcf.VCFFileReader)} and match to an array of bams
   *
   * @param samples samples to match
   * @param variantSets variant sets that may be appended to the vcf sample names
   * @param bams the bam files
   * @param numThreads
   * @param log
   * @return Hashtable of the sample -> bam file mapping
   */
  public static Map<String, String> matchToVcfSamplesToBamFiles(String[] samples,
                                                                Set<String> variantSets,
                                                                String[] bams, int numThreads,
                                                                Logger log) {
    HashMap<String, String> matched = new HashMap<String, String>();
    HashMap<String, String> bamSamples = new HashMap<String, String>();
    SampleNameProducer producer = new SampleNameProducer(bams);
    WorkerTrain<SampleNameExtractor> train = new WorkerTrain<SampleNameExtractor>(producer,
                                                                                  numThreads, 10,
                                                                                  log);

    while (train.hasNext()) {
      SampleNameExtractor ex = train.next();
      String bamSamp = ex.getName();
      if (bamSamples.containsKey(bamSamp)) {
        throw new IllegalArgumentException("Bams must be sample unique");
      } else {
        bamSamples.put(bamSamp, ex.getBamFile());
        if (variantSets != null) {
          for (String set : variantSets) {
            bamSamples.put(bamSamp + set, ex.getBamFile());
          }
        }
      }
    }

    for (int i = 0; i < samples.length; i++) {
      if (!bamSamples.containsKey(samples[i])) {
        log.reportTimeWarning("Did not find matching bam file for " + samples[i]);
      } else {
        if (matched.containsKey(samples[i])) {
          throw new IllegalArgumentException("Multiple bam files matched sample " + samples[i]
                                             + ", perhaps because of variant sets?");
        }
        matched.put(samples[i], bamSamples.get(samples[i]));
      }
    }

    log.reportTimeInfo("Found matching bam files for" + matched.size() + " of " + samples.length
                       + " samples");

    return matched;

  }

  // TODO
  public static void removeSegments(String inputBam, String outputBam, String bedFile, int buffer,
                                    Logger log) {


  }

}
