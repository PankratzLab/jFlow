package org.genvisis.seq;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.manage.SeqOps;
import org.genvisis.seq.manage.SeqOps.GC_COMP_METHOD;
import org.genvisis.seq.manage.VCOps;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.ArrayUtils.BYTE_DECODE_FORMAT;
import org.pankratzlab.common.filesys.LocusSet;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.filesys.Segment;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;

public class ReferenceGenome {

  public static final GenomeBuild DEFAULT_BUILD =GenomeBuild.HG19;
  private final String referenceFasta;
  private IndexedFastaSequenceFile indexedFastaSequenceFile;
  private byte[] currentSeq;
  private byte[] inMemoryContig;
  private int defaultBuffer;
  private ReferenceSequence referenceSequence;
  private final Logger log;

  public ReferenceGenome(String referenceFasta, Logger log) {
    super();
    this.referenceFasta = referenceFasta;
    this.log = log;
    indexedFastaSequenceFile = null;
    referenceSequence = null;
    defaultBuffer = 0;
    try {
      indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(referenceFasta));
    } catch (FileNotFoundException e) {
      log.reportFileNotFound(referenceFasta);
      e.printStackTrace();
    }
  }

  public SAMSequenceDictionary getDictionary() {
    return indexedFastaSequenceFile.getSequenceDictionary();
  }

  /**
   * @param bpBinSize the bin size (non-sliding) to break the reference genome into;
   * @return
   */
  public LocusSet<Segment> getBins(int bpBinSize) {
    ArrayList<Segment> bins = new ArrayList<>();
    SAMSequenceDictionary samSequenceDictionary = indexedFastaSequenceFile.getSequenceDictionary();
    System.out.println(samSequenceDictionary.getSequences().size() + " contigs detected");
    for (SAMSequenceRecord samSequenceRecord : samSequenceDictionary.getSequences()) {
      int length = samSequenceRecord.getSequenceLength();
      byte chr = Positions.chromosomeNumber(samSequenceRecord.getSequenceName());
      int currentStart = 1;
      int currentStop = Math.min(length, bpBinSize);
      int binsAdded = 0;
      while (currentStop < length) {
        Segment seg = new Segment(chr, currentStart, currentStop);
        currentStart += bpBinSize + 1;
        currentStop += bpBinSize + 1;
        currentStop = Math.min(length, currentStop);
        if (seg.getSize() != bpBinSize) {
          log.reportTimeWarning("bin " + seg.getUCSClocation() + " size (" + seg.getSize()
                                + ") does not equal the bin size of " + bpBinSize);
        }
        bins.add(seg);
        binsAdded++;
      }
      if (currentStop == length) {
        Segment seg = new Segment(chr, currentStart, currentStop);
        bins.add(seg);
        if (seg.getSize() != bpBinSize) {
          log.reportTimeWarning("bin " + seg.getUCSClocation() + " size (" + seg.getSize()
                                + ") does not equal the bin size of " + bpBinSize);
        }
      } else {
        String error = "BUG: End of bin for " + samSequenceRecord.getSequenceName()
                       + " did not end at " + length;
        log.reportError(error);
        throw new IllegalStateException(error);
      }
      log.reportTimeInfo(samSequenceRecord.getSequenceName() + " -> " + length + "bp; "
                         + (binsAdded + 1) + " " + bpBinSize + "bp bins");
    }

    LocusSet<Segment> binsToReturn = new LocusSet<Segment>(bins.toArray(new Segment[bins.size()]),
                                                           true, log) {

      /**
       * 
       */
      private static final long serialVersionUID = 1L;

    };
    log.reportTimeInfo(referenceFasta + " broken up to " + bins.size() + " bins ");
    return binsToReturn;
  }

  public IndexedFastaSequenceFile getIndexedFastaSequenceFile() {
    return indexedFastaSequenceFile;
  }

  public String getReferenceFasta() {
    return referenceFasta;
  }

  public ReferenceSequence getReferenceSequence() {
    return referenceSequence;
  }

  public int getDefaultBuffer() {
    return defaultBuffer;
  }

  public void setDefaultBuffer(int defaultBuffer) {
    this.defaultBuffer = defaultBuffer;
  }

  public String[][] getSequencesFor(Segment[] segs, boolean memoryMode) {
    return getSequencesFor(segs, -1, memoryMode);
  }

  /**
   * @param segs these should be sorted in chromosomal order if speed is desired
   * @return
   */
  public String[][] getSequencesFor(Segment[] segs, int reportEvery, boolean memoryMode) {
    String[][] seqs = new String[segs.length][];
    for (int i = 0; i < seqs.length; i++) {
      if (reportEvery > 0 && i % reportEvery == 0) {
        log.reportTimeInfo((i + 1) + " segments queried of " + segs.length + " total from "
                           + referenceFasta);
        log.memoryFree();
        log.memoryTotal();
      }
      seqs[i] = getSequenceFor(segs[i], memoryMode);
    }
    return seqs;
  }

  /**
   * Stores results of contig presence/absence, as well as if formatting was needed (remove "chr"
   * etc)
   */
  public static class ContigRequest {

    private final String requestedContig;
    private final boolean hasContig;
    private final String contigPresent;

    /**
     * @param requestedContig the requested contig
     * @param hasContig whether this contig is present as request, or trimmed
     * @param contigPresent the contig actually present
     */
    private ContigRequest(String requestedContig, boolean hasContig, String contigPresent) {
      super();
      this.requestedContig = requestedContig;
      this.hasContig = hasContig;
      this.contigPresent = contigPresent;
    }

    /**
     * @return the requestedContig
     */
    public String getRequestedContig() {
      return requestedContig;
    }

    /**
     * @return the hasContig
     */
    public boolean hasContig() {
      return hasContig;
    }

    /**
     * @return the contigPresent
     */
    public String getContigPresent() {
      return contigPresent;
    }

  }

  public ContigRequest requestContig(String contig) {
    if (indexedFastaSequenceFile.getSequenceDictionary() == null) {
      log.reportError("Could not find sequence dictionary for " + referenceFasta + " ("
                      + ext.rootOf(referenceFasta) + ".dict)"
                      + ", will not be able to check if contig " + contig + " is present");
      return new ContigRequest(contig, false, null);
    }
    if (indexedFastaSequenceFile.getSequenceDictionary().getSequence(contig) != null) {
      return new ContigRequest(contig, true, contig);
    }
    String contigTrim = contig.replaceAll("chr", "");
    if (indexedFastaSequenceFile.getSequenceDictionary().getSequence(contigTrim) != null) {
      return new ContigRequest(contig, true, contigTrim);
    }
    return new ContigRequest(contig, false, null);
  }

  public String[] getSequenceFor(Segment segment) {
    return getSequenceFor(segment, ASSEMBLY_NAME.HG19, false);
  }

  public String[] getSequenceFor(Segment segment, boolean memoryMode) {
    return getSequenceFor(segment, ASSEMBLY_NAME.HG19, memoryMode);
  }

  public int getContigLength(Segment seg) {
    return getContigLength(Positions.getChromosomeUCSC(seg.getChr(), true, true));
  }

  public int getContigLength(String contig) {
    ContigRequest cr = requestContig(contig);
    if (cr.hasContig()) {
      return indexedFastaSequenceFile.getSequenceDictionary().getSequence(cr.getContigPresent())
                                     .getSequenceLength();
    } else {
      return -1;
    }
  }

  /**
   * @param segment
   * @param memoryMode this stores an entire contig in memory, which is faster for many large, in
   *          order queries.
   * @return will return null if {@link Positions#getChromosomeUCSC(int, boolean, boolean)} returns
   *         a contig not in the files {@link SAMSequenceDictionary}
   */
  public String[] getSequenceFor(Segment segment, ASSEMBLY_NAME aName, boolean memoryMode) {
    String requestedContig = Positions.getChromosomeUCSC(segment.getChr(), aName.addChr(), true);
    if (segment.getChr() == 26) {
      requestedContig = aName.getMitoContig();
    }

    ContigRequest cr = requestContig(requestedContig);

    if (cr.hasContig()) {
      requestedContig = cr.getContigPresent();
      int seqLength = indexedFastaSequenceFile.getSequenceDictionary().getSequence(requestedContig)
                                              .getSequenceLength();
      int start = segment.getStart() - defaultBuffer;
      if (start < 0) {
        start = 0;
      }
      int stop = segment.getStop() + defaultBuffer;
      if (stop > seqLength) {
        stop = seqLength - 1;
      }
      String[] requestedSeq = null;
      if (memoryMode) {
        if (referenceSequence == null || !referenceSequence.getName().equals(requestedContig)) {
          log.reportTimeInfo("loading reference sequence for " + requestedContig);
          referenceSequence = indexedFastaSequenceFile.getSequence(requestedContig);
          inMemoryContig = referenceSequence.getBases();
          log.reportTimeInfo("reference sequence in memory for " + requestedContig);

        } else {
          // log.reportTimeInfo("Memory works");
        }
        try {
          int startIndex = Math.max(0, start - 1);
          int stopIndex = Math.min(inMemoryContig.length - 1, stop);
          requestedSeq = new String[stopIndex - startIndex];
          for (int i = 0; i < requestedSeq.length; i++) {
            requestedSeq[i] = new String(new byte[] {inMemoryContig[startIndex + i]},
                                         ext.UTF_8).toUpperCase();
          }
        } catch (Exception e) {
          log.reportError("Invalid query " + segment.getUCSClocation() + "; buffer " + defaultBuffer
                          + "; current contig " + referenceSequence.getName());
        }
      } else {
        ReferenceSequence subReferenceSequence = indexedFastaSequenceFile.getSubsequenceAt(requestedContig,
                                                                                           start,
                                                                                           stop);
        requestedSeq = ArrayUtils.decodeByteArray(subReferenceSequence.getBases(),
                                                  BYTE_DECODE_FORMAT.UPPER_CASE, log);
      }
      return requestedSeq;

    } else {
      // log.reportTimeError("Requested contig " + requestedContig + " was not in the sequence
      // dictionary for " + referenceFasta);
      return null;
    }

  }

  /**
   * This method is slow, it loads the complete chromosome into memory and then subsets.
   */
  @Deprecated
  public String[] getSequenceForOld(Segment segment) {
    String requestedContig = Positions.getChromosomeUCSC(segment.getChr(), true);
    if (referenceSequence == null || currentSeq == null
        || !referenceSequence.getName().equals(requestedContig)) {
      referenceSequence = indexedFastaSequenceFile.getSequence(requestedContig);
      currentSeq = referenceSequence.getBases();
    }
    int start = segment.getStart() - 1 - defaultBuffer;
    int stop = segment.getStop() + defaultBuffer;
    if (start < 0) {
      log.reportTimeWarning("Buffer of " + defaultBuffer
                            + " adjusts base pair extraction to index less than 0, adjusting start to index 0 (bp 1)");
      start = 0;
    }
    if (stop >= currentSeq.length) {
      stop = currentSeq.length - 1;
      log.reportTimeWarning("Buffer of " + defaultBuffer
                            + " ,adjusts base pair extraction to index greater than sequence length,  adjusting stop to index "
                            + (currentSeq.length - 1));
    }
    byte[] subTmp = null;
    try {
      subTmp = Arrays.copyOfRange(currentSeq, start, stop);
    } catch (Exception e) {
      log.reportError("Could not extract bases:");
      log.reportError("Segment: " + segment.getUCSClocation());
      log.reportError("Ref: " + referenceSequence.getName());
      log.reportError("Start : " + start);
      log.reportError("Stop : " + stop);
      return new String[] {};
    }
    if (!referenceSequence.getName().equals(requestedContig)) {
      log.reportError("Mismatched request");
      log.reportError("Segment: " + segment.getUCSClocation());
      log.reportError("Ref: " + referenceSequence.getName());
      log.reportError("Start : " + start);
      log.reportError("Stop : " + stop);
      return new String[] {};
    }

    String[] requestedSeq = new String[subTmp.length];
    try {
      for (int i = 0; i < requestedSeq.length; i++) {
        requestedSeq[i] = new String(new byte[] {subTmp[i]}, ext.UTF_8).toUpperCase();
      }
    } catch (UnsupportedEncodingException e) {
      log.reportError("Could not convert reference byte to string");
      e.printStackTrace();
    }
    return requestedSeq;
  }

  public double getGCContentFor(Segment seg) {
    return getGCContentFor(seg, false);
  }

  public double getGCContentFor(Segment seg, boolean memoryMode) {
    String[] seq = getSequenceFor(seg, memoryMode);
    if (seq != null) {
      return SeqOps.getProportionGC(seq, GC_COMP_METHOD.GCTA_ONLY);

    } else {
      return Double.NaN;
    }
  }

  public static double getPercentGC(String[] seq) {
    return SeqOps.getProportionGC(seq, GC_COMP_METHOD.GCTA_ONLY) * 100;
  }

  public double getGCContentFor(VariantContext vc) {
    // defaultBuffer=0;
    // System.out.println("REF\t"+Array.toStr(getSequenceFor(VCOps.getSegment(vc))));
    // System.out.println("VARIANT_A\t"+vc.getAlleles().toString());
    return getGCContentFor(VCOps.getSegment(vc), false);
  }

}
