package org.genvisis.seq.manage.vcfPile;

import java.io.Serializable;
import java.util.HashSet;

import org.genvisis.common.Logger;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

class PiledVcfRegion<T extends Segment> implements Serializable {
  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  private final T region;
  private final double avgGC;
  private final String[] vcfSamples;
  private final int[] numHomVar;
  private final int[] numHet;
  private final int[] numHomRef;
  private final int[] totalCalledVar;
  private final int[] noCall;
  private final double[] avgGQ;
  private final double[] avgDP;

  public PiledVcfRegion(T region, ReferenceGenome referenceGenome, String[] vcfSamples) {
    super();
    this.region = region;
    this.avgGC = referenceGenome.getGCContentFor(region);
    this.vcfSamples = vcfSamples;
    this.numHomVar = new int[vcfSamples.length];
    this.numHet = new int[vcfSamples.length];
    this.numHomRef = new int[vcfSamples.length];
    this.totalCalledVar = new int[vcfSamples.length];
    this.noCall = new int[vcfSamples.length];
    this.avgGQ = new double[vcfSamples.length];
    this.avgDP = new double[vcfSamples.length];
  }

  public void addVariantContext(VariantContext vc, Logger log) {
    if (!VCOps.getSegment(vc).overlaps(region)) {
      String error = "Attempting to add a variant does not overlap the region of interest";
      error += "\nCurrent region = " + region.getUCSClocation();
      error += "\nVariant region = " + VCOps.getSegment(vc).getUCSClocation();
      log.reportTimeError(error);
      throw new IllegalStateException(error);
    } else {
      for (int i = 0; i < vcfSamples.length; i++) {
        HashSet<String> tmp = new HashSet<String>();
        tmp.add(vcfSamples[i]);
        VariantContext vcTmp = VCOps.getSubset(vc, tmp, VC_SUBSET_TYPE.SUBSET_STRICT);
        Genotype g = vcTmp.getGenotype(0);

        if (g.isHomRef()) {
          numHomRef[i]++;
        } else if (g.isHet()) {
          numHet[i]++;
        } else if (g.isHomVar()) {
          numHomVar[i]++;
        } else if (g.isNoCall()) {
          noCall[i]++;
        } else {
          String error = "Invalid variant type with genotype " + g.toString();
          log.reportTimeError(error);
          throw new IllegalStateException(error);
        }
        if (!g.isNoCall()) {
          totalCalledVar[i]++;
        }
        avgGQ[i] += g.getGQ();
        avgDP[i] += g.getDP();
      }
    }
  }

  public double[] getAvgDP() {
    return avgDP;
  }

  public double getAvgGC() {
    return avgGC;
  }

  public double[] getAvgGQ() {
    return avgGQ;
  }

  public T getRegion() {
    return region;
  }

  public int[] getTotalCalledVar() {
    return totalCalledVar;
  }

}
