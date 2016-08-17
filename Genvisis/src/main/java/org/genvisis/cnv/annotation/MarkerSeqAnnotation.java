package org.genvisis.cnv.annotation;

import java.util.List;

import org.genvisis.cnv.annotation.BlastAnnotationTypes.TOP_BOT;
import org.genvisis.cnv.util.CNVHelper;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCOps;

import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class MarkerSeqAnnotation extends AnnotationData {

  private static final String DEFAULT_NAME = "PROBE_DESIGN";
  private static final String DESCRIPTION =
      "The probe sequence for the A probe, probe sequence for the B probe (if different than A) interrogation position,strand by design,TOB_BOTTOM SNP designation,TOP_BOTTOM Reference Designation, A allele, and B allele ";

  public static MarkerSeqAnnotation getDefault() {
    return new MarkerSeqAnnotation();
  }

  private String seqA;
  private String seqB;
  private Strand strand;
  private int interrogationPosition;
  private Segment seg;
  private TOP_BOT topBotProbe;
  private TOP_BOT topBotRef;
  private Allele A;
  private Allele B;
  private Allele ref;
  private Allele[] alts;
  private boolean indel;
  private boolean aInsertion;
  private boolean bInsertion;
  private boolean aDeletion;

  private boolean bDeletion;

  public MarkerSeqAnnotation() {
    super(VCFHeaderLineType.String, null, 1, DEFAULT_NAME, DESCRIPTION, DEFUALT_VALUE,
        DEFUALT_VALUE);
  }

  public Allele getA() {
    return A;
  }

  public Allele[] getAlts() {
    return alts;
  }

  public Allele getB() {
    return B;
  }

  public int getInterrogationPosition() {
    return interrogationPosition;
  }

  public Allele getRef() {
    return ref;
  }

  public Segment getSeg() {
    return seg;
  }

  public String getSeqB() {
    return seqB;
  }

  public String getSequence() {
    return seqA;
  }

  public Strand getStrand() {
    return strand;
  }

  public TOP_BOT getTopBotProbe() {
    return topBotProbe;
  }

  public TOP_BOT getTopBotRef() {
    return topBotRef;
  }

  public boolean goLeft() {
    return topBotProbe == topBotRef;
  }

  public boolean isaDeletion() {
    return aDeletion;
  }

  public boolean isaInsertion() {
    return aInsertion;
  }

  public boolean isbDeletion() {
    return bDeletion;
  }

  public boolean isbInsertion() {
    return bInsertion;

  }

  public boolean isIndel() {
    return indel;
  }

  @Override
  public void parseAnnotation(VariantContext vc, Logger log) {
    if (vc.hasAttribute(getName())) {
      setData(vc.getAttributeAsString(getName(), DEFAULT_NAME));
      if (!ext.isMissingValue(getData())) {
        List<String> data = getDataAsList();
        seqA = data.get(0);
        seqB = data.get(1);

        interrogationPosition = -1;
        try {
          interrogationPosition = Integer.parseInt(data.get(2));

        } catch (NumberFormatException nfe) {

        }
        strand = Strand.toStrand(data.get(3));
        seg = VCOps.getSegment(vc);
        topBotProbe = TOP_BOT.valueOf(data.get(4));
        topBotRef = TOP_BOT.valueOf(data.get(5));
        ref = vc.getReference();
        List<Allele> alleles = vc.getAlternateAlleles();
        alts = new Allele[alleles.size()];
        for (int i = 0; i < alts.length; i++) {
          alts[i] = alleles.get(i);
        }
        A = Allele.create(data.get(6), ref.basesMatch(data.get(6)));
        B = Allele.create(data.get(7), ref.basesMatch(data.get(7)));
        if (A.isReference() && B.isReference()) {
          throw new IllegalArgumentException("A and B alleles cannot both be reference");
        }
        indel = vc.isIndel();
        if (indel) {
          aInsertion = A.getBases().length > ref.getBases().length;
          bInsertion = B.getBases().length > ref.getBases().length;
          aDeletion = A.getBases().length < ref.getBases().length;
          bDeletion = B.getBases().length < ref.getBases().length;

          if ((aDeletion && aInsertion) || (bDeletion && bInsertion)) {
            throw new IllegalStateException("Allele cannot be both insertion and deletion");
          }
        }
      }
    }
  }

  public void setDesignData(String seqA, String segB, int interrogationPosition, Strand strand,
      TOP_BOT topBotProbe, TOP_BOT topBotRef, Allele A, Allele B) {
    this.seqA = seqA;
    seqB = segB;
    this.interrogationPosition = interrogationPosition;
    this.strand = strand;
    this.topBotProbe = topBotProbe;
    this.topBotRef = topBotRef;
    this.A = A;
    this.B = B;
    // this.seg=seg; populate on load only
    setData(seqA + DEFUALT_DELIMITER + segB + DEFUALT_DELIMITER + interrogationPosition
        + DEFUALT_DELIMITER + CNVHelper.decode(strand) + DEFUALT_DELIMITER + topBotProbe
        + DEFUALT_DELIMITER + topBotRef + DEFUALT_DELIMITER + A.getDisplayString()
        + DEFUALT_DELIMITER + B.getDisplayString());
  }

  public void setSegment(Segment seg1) {
    seg = seg1;
  }

  public void setStrand(Strand strand2) {
    strand = strand2;
  }

}
