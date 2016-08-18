package org.genvisis.seq.qc;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.genvisis.common.Logger;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilterPass;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantFilterSample {
  private final String filterName;
  private final VariantContextFilter vCFiltAltGeno;
  private final VariantContextFilter vCFiltHomRefGeno;
  private final VariantContextFilter vCFiltAllGeno;
  private final VariantContextFilter vCFiltIndel;
  private final VariantContextFilter vCFiltWholeVariant;
  private final Logger log;

  public VariantFilterSamplePass filter(VariantContext vc, FILTER_METHOD method) {
    VariantFilterSamplePass vcPass = new VariantFilterSamplePass();
    Set<String> samples = vc.getSampleNames();
    HashSet<String> samplesToSetMissing = new HashSet<String>();
    VariantContextFilterPass tmp = null;
    if (vCFiltWholeVariant != null) {
      tmp = vCFiltWholeVariant.filter(vc);
      vcPass.getVcfps().add(tmp);
    }
    if (tmp == null || tmp.passed()) {
      HashSet<String> passingSamples = new HashSet<String>();
      int samplesTested = 0;
      for (String sample : samples) {
        HashSet<String> curSamp = new HashSet<String>();
        curSamp.add(sample);
        VariantContext vcSamp = VCOps.getSubset(vc, curSamp, VC_SUBSET_TYPE.SUBSET_STRICT, false);
        Genotype g = vcSamp.getGenotype(0);

        boolean pass = true;
        if (pass && vCFiltAllGeno != null) {
          tmp = vCFiltAllGeno.filter(vcSamp);
          vcPass.getVcfps().add(tmp);
          pass = tmp.passed();
          if (!pass) {
            samplesTested++;
          }
        }
        if (pass && vCFiltIndel != null && vcSamp.isIndel()) {
          tmp = vCFiltIndel.filter(vcSamp);
          vcPass.getVcfps().add(tmp);
          pass = tmp.passed();
          if (!pass) {
            samplesTested++;
          }
        }
        if (pass && vCFiltAltGeno != null && !g.isHomRef() && !g.isNoCall()) {// alt
          tmp = vCFiltAltGeno.filter(vcSamp);
          vcPass.getVcfps().add(tmp);
          pass = tmp.passed();
          if (!pass) {
            samplesTested++;
          }
        }
        if (pass && vCFiltHomRefGeno != null && g.isHomRef()) {
          tmp = vCFiltHomRefGeno.filter(vcSamp);
          vcPass.getVcfps().add(tmp);
          pass = tmp.passed();
          if (!pass) {
            samplesTested++;
          }
        }

        if (pass) {
          samplesTested++;

          passingSamples.add(sample);
        } else {

          samplesToSetMissing.add(sample);
        }
      }
      if (samplesTested != vc.getSampleNames().size()) {
        throw new IllegalStateException("Not all samples passed through the filters, internal error ,tested "
                                        + samplesTested + " and should have tested "
                                        + vc.getSampleNames().size());
      }
      switch (method) {
        case REMOVE_FROM_CONTEXT:// keep passing only
          vcPass.setPassingContext(VCOps.getSubset(vc, passingSamples, VC_SUBSET_TYPE.SUBSET_STRICT,
                                                   false));
          break;
        case SET_TO_MISSING:// set non passing to missing
          vcPass.setPassingContext(VCOps.setTheseSamplesToMissing(vc, samplesToSetMissing));
          break;
        default:
          log.reportTimeError("Invalid method " + method);
          return null;
      }
    } else {
      switch (method) {
        case REMOVE_FROM_CONTEXT:// remove all
          vcPass.setPassingContext(VCOps.getSubset(vc, new HashSet<String>(),
                                                   VC_SUBSET_TYPE.SUBSET_STRICT, false));
          break;
        case SET_TO_MISSING:// set all to missing
          samplesToSetMissing.addAll(samples);
          vcPass.setPassingContext(VCOps.setTheseSamplesToMissing(vc, samplesToSetMissing));
          break;
        default:
          log.reportTimeError("Invalid method " + method);
          return null;
      }
    }
    return vcPass;
  }

  public String getFilterName() {
    return filterName;
  }

  public enum FILTER_METHOD {
                             /**
                              * Samples that do not pass the filters will be removed from the final
                              * {@link VariantContext}
                              */
                             REMOVE_FROM_CONTEXT,
                             /**
                              * Samples that do not pass the filters will be set to missing in the
                              * final {@link VariantContext}
                              */
                             SET_TO_MISSING;
  }

  public static class VariantFilterSamplePass {
    private VariantContext passingContext;
    private ArrayList<VariantContextFilterPass> vcfps;

    public VariantFilterSamplePass() {
      super();
      vcfps = new ArrayList<VariantContextFilterPass>();
    }

    public VariantContext getPassingContext() {
      return passingContext;
    }

    public void setPassingContext(VariantContext passingContext) {
      this.passingContext = passingContext;
    }

    public ArrayList<VariantContextFilterPass> getVcfps() {
      return vcfps;
    }

    public void setVcfps(ArrayList<VariantContextFilterPass> vcfps) {
      this.vcfps = vcfps;
    }

  }

  public static class Builder {
    private VariantContextFilter vCFiltAltGeno = null;
    private VariantContextFilter vCFiltHomRefGeno = null;
    private VariantContextFilter vCFiltAllGeno = null;
    private VariantContextFilter vCFiltIndel = null;
    private VariantContextFilter vCFiltWholeVariant = null;

    /**
     * @param vCFiltAltGeno this filter will be applied to alternate calls only
     * @return
     */
    public Builder vCFiltAltGeno(VariantContextFilter vCFiltAltGeno) {
      this.vCFiltAltGeno = vCFiltAltGeno;
      return this;
    }

    /**
     * @param vCFiltHomRefGeno this filter will be applied to reference (Hom Ref) calls only
     * @return
     */
    public Builder vCFiltHomRefGeno(VariantContextFilter vCFiltHomRefGeno) {
      this.vCFiltHomRefGeno = vCFiltHomRefGeno;
      return this;
    }

    /**
     * @param vCFiltAllGeno this filter will be applied to all genotypes
     * @return
     */
    public Builder vCFiltAllGeno(VariantContextFilter vCFiltAllGeno) {
      this.vCFiltAllGeno = vCFiltAllGeno;
      return this;
    }

    /**
     * @param vCFiltIndel this filter will be applied to Indels only
     * @return
     */
    public Builder vCFiltIndel(VariantContextFilter vCFiltIndel) {
      this.vCFiltIndel = vCFiltIndel;
      return this;
    }

    /**
     * @param vCFiltWholeVariant this filter will be applied to the whole variant (Annotations,
     *        biallelic, etc)
     * @return
     */
    public Builder vCFiltWholeVariant(VariantContextFilter vCFiltWholeVariant) {
      this.vCFiltWholeVariant = vCFiltWholeVariant;
      return this;
    }

    public VariantFilterSample build(String filterName, Logger log) {
      return new VariantFilterSample(this, filterName, log);
    }
  }

  private VariantFilterSample(Builder builder, String filterName, Logger log) {
    this.filterName = filterName;
    vCFiltAltGeno = builder.vCFiltAltGeno;
    vCFiltHomRefGeno = builder.vCFiltHomRefGeno;
    vCFiltAllGeno = builder.vCFiltAllGeno;
    vCFiltIndel = builder.vCFiltIndel;
    vCFiltWholeVariant = builder.vCFiltWholeVariant;
    this.log = log;
  }
}
