package org.genvisis.seq.manage;

import java.io.Closeable;
import java.io.IOException;
import java.util.ArrayList;

import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author lane0212<br>
 *         Trying to make an indexed bed file reader, mimicking {@link VCFFileReader}
 */
public class BEDFileReader implements Closeable, Iterable<BEDFeature> {
  private final FeatureReader<BEDFeature> reader;

  public BEDFileReader(final String file, final boolean requireIndex) {
    reader = AbstractFeatureReader.getFeatureReader(file, new BEDCodec(), requireIndex);
  }

  /** Queries for records within the region specified. */
  public CloseableIterator<BEDFeature> query(final String chrom, final int start, final int end) {
    try {
      return reader.query(chrom, start, end);
    } catch (final IOException ioe) {
      throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
    }
  }

  /** Queries for records within the region specified. */
  public CloseableIterator<BEDFeature> query(Segment seg) {
    return query(Positions.getChromosomeUCSC(seg.getChr(), true), seg.getStart(), seg.getStop());
  }

  public LocusSet<BEDFeatureSeg> loadSegsFor(Segment segment, Logger log) {
    CloseableIterator<BEDFeature> iterator = query(segment);
    ArrayList<BEDFeatureSeg> bedSegs = new ArrayList<BEDFeatureSeg>();
    while (iterator.hasNext()) {
      BEDFeature bedFeature = iterator.next();
      bedSegs.add(new BEDFeatureSeg(bedFeature, log));

    }
    LocusSet<BEDFeatureSeg> segSet =
                                   new LocusSet<BEDFeatureSeg>(bedSegs.toArray(new BEDFeatureSeg[bedSegs.size()]),
                                                               true, log) {

                                     /**
                                      * 
                                      */
                                     private static final long serialVersionUID = 1L;

                                   };

    return segSet;
  }

  public LocusSet<BEDFeatureSeg> loadAll(Logger log) {
    ArrayList<BEDFeatureSeg> bedSegs = new ArrayList<BEDFeatureSeg>();
    CloseableIterator<BEDFeature> iterator;
    try {
      iterator = reader.iterator();
      while (iterator.hasNext()) {
        BEDFeature bedFeature = iterator.next();
        bedSegs.add(new BEDFeatureSeg(bedFeature, log));
      }
    } catch (IOException e) {
      log.reportException(e);
      e.printStackTrace();
      return null;
    }

    LocusSet<BEDFeatureSeg> segSet =
                                   new LocusSet<BEDFeatureSeg>(bedSegs.toArray(new BEDFeatureSeg[bedSegs.size()]),
                                                               true, log) {

                                     /**
                                      * 
                                      */
                                     private static final long serialVersionUID = 1L;

                                   };

    return segSet;

  }

  public static class BEDFeatureSeg extends Segment {
    private final BEDFeature bedFeature;
    /**
     *
     */
    private static final long serialVersionUID = 1L;

    public BEDFeatureSeg(BEDFeature bedFeature, Logger log) {
      super(bedFeature.getContig(), bedFeature.getStart(), bedFeature.getEnd());
      this.bedFeature = bedFeature;
    }

    public BEDFeature getBedFeature() {
      return bedFeature;
    }

  }

  public BEDFeature[] loadBEDFeaturesFor(Segment segment, Logger log) {
    CloseableIterator<BEDFeature> iterator = query(segment);
    ArrayList<BEDFeature> bedSegs = new ArrayList<BEDFeature>();
    while (iterator.hasNext()) {
      bedSegs.add(iterator.next());
    }

    return bedSegs.toArray(new BEDFeature[bedSegs.size()]);
  }

  @Override
  public void close() {
    try {
      reader.close();
    } catch (final IOException ioe) {
      throw new TribbleException("Could not close a bed context feature reader.", ioe);
    }
  }

  /** Returns an iterator over all records in this VCF/BCF file. */
  @Override
  public CloseableIterator<BEDFeature> iterator() {
    try {
      return reader.iterator();
    } catch (final IOException ioe) {
      throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
    }
  }

}
