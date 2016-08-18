package org.genvisis.cnv.annotation;

import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.filesys.Segment;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class LocusAnnotation {
  private final String locusName;
  private final Segment seg;
  private Allele ref;
  private Allele[] alts;
  private final Logger log;
  private AnnotationData[] annotations;

  public void setRef(Allele ref) {
    this.ref = ref;
  }

  public void setAlts(Allele[] alts) {
    this.alts = alts;
  }

  /**
   * @param skipDefaultValue if true, annotations that are set to their default value will not be
   *        added to the {@link VariantContext}
   * @return
   */
  public VariantContext annotateLocus(boolean skipDefaultValue) {
    ArrayList<Allele> alleles = new ArrayList<Allele>();
    Allele refA = Allele.create(ref, false);
    alleles.add(refA);
    int maxAlleleLength = refA.getBases().length - 1;
    boolean isIndel = false;
    for (Allele alt : alts) {
      Allele altA = Allele.create(alt, false);
      alleles.add(altA);
      if (!isIndel) {
        isIndel = alt.getBases().length != refA.getBases().length;
      }
    }

    VariantContextBuilder vBuilder = new VariantContextBuilder();
    vBuilder.chr(Positions.getChromosomeUCSC(seg.getChr(), true));
    if (seg.getChr() == 0) {
      vBuilder.start(seg.getStart() <= 0 ? 1 : seg.getStart());
      vBuilder.stop(seg.getStop() <= 0 ? 1 + maxAlleleLength : seg.getStop() + maxAlleleLength);

    } else {
      vBuilder.start(seg.getStart());
      vBuilder.stop(seg.getStop() + maxAlleleLength);
    }
    vBuilder.id(locusName);
    vBuilder.alleles(alleles);
    if (annotations != null) {
      for (int i = 0; i < annotations.length; i++) {
        if (!skipDefaultValue
            || !annotations[i].getDefaultValue().equals(annotations[i].getData())) {
          annotations[i].addAnnotation(vBuilder);
        }
      }
    }
    VariantContext vc = null;
    try {
      vc = vBuilder.make();
    } catch (Exception e) {
      log.reportException(e);
      log.reportTimeError("Could not create VC at " + seg.getUCSClocation() + " for " + locusName);
      log.reportTimeError("Ref " + ref.getDisplayString());
      for (Allele alt : alts) {
        log.reportTimeError("Alt " + alt.getDisplayString());

      }
      for (AnnotationData annotation : annotations) {
        log.reportTimeError("Annotation: " + annotation.getData());
      }
    }
    return vc;
  }

  public Segment getSeg() {
    return seg;
  }

  public String getLocusName() {

    return locusName;
  }

  public void setAnnotations(AnnotationData[] annotations) {
    this.annotations = annotations;
  }

  public AnnotationData[] getAnnotations() {
    return annotations;
  }

  public void addAnnotation(AnnotationData annotationData) {
    annotations = Array.concatAll(annotations, new AnnotationData[] {annotationData});
  }

  public static class Builder {
    private Allele ref = Allele.create("A", true);
    private Allele[] alts = new Allele[] {Allele.create("N", false)};
    private AnnotationData[] annotations = null;
    private Logger log = new Logger();

    public Builder ref(Allele ref) {
      this.ref = ref;
      return this;
    }

    public Builder log(Logger log) {
      this.log = log;
      return this;
    }

    public Builder alts(Allele[] alts) {
      this.alts = alts;
      return this;
    }

    public Builder annotations(AnnotationData[] annotations) {
      this.annotations = annotations;
      return this;
    }

    public LocusAnnotation build(String locusName, Segment seg) {
      return new LocusAnnotation(this, locusName, seg);
    }
  }

  public LocusAnnotation(Builder builder, String locusName, Segment seg) {
    this.locusName = locusName;
    this.seg = seg;
    annotations = builder.annotations;
    ref = builder.ref;
    alts = builder.alts;
    log = builder.log;
  }
}
