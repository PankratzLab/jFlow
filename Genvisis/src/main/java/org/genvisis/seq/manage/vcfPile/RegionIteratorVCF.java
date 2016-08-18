package org.genvisis.seq.manage.vcfPile;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCOps;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author lane0212 Iterate by region over a vcf;
 *
 */
public class RegionIteratorVCF<T extends Segment> implements Iterator<VariantContext[]> {
  private String vcfFile;// force the file name to prevent reader being passed across threads
  private final VCFFileReader reader;
  private final LocusSet<T> set;
  private int setIndex;
  private final Logger log;

  public RegionIteratorVCF(String vcfFile, LocusSet<T> segs, Logger log) {
    super();
    this.reader = new VCFFileReader(new File(vcfFile), true);
    this.setIndex = 0;
    this.set = segs;
    this.log = log;
  }

  public Logger getLog() {
    return log;
  }

  public String getVcfFile() {
    return vcfFile;
  }

  @Override
  public boolean hasNext() {
    boolean hasNext = setIndex < set.getLoci().length;
    if (!hasNext) {
      reader.close();
    }
    return hasNext;
  }

  public T getCurrentIndex() {
    return set.getLoci()[setIndex];
  }

  @Override
  public VariantContext[] next() {
    T region = set.getLoci()[setIndex];
    setIndex++;
    CloseableIterator<VariantContext> rIter =
        reader.query(Positions.getChromosomeUCSC(region.getChr(), true), region.getStart(),
                     region.getStop());
    ArrayList<VariantContext> tmp = new ArrayList<VariantContext>();
    while (rIter.hasNext()) {
      VariantContext vc = rIter.next();// apparently an approximate query
      if (VCOps.getSegment(vc).overlaps(region)) {
        tmp.add(vc);
      }
    }
    return tmp.toArray(new VariantContext[tmp.size()]);
  }

  @Override
  public void remove() {
    // TODO Auto-generated method stub
  }

}
