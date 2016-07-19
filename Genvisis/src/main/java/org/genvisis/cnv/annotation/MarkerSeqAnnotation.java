package org.genvisis.cnv.annotation;

import java.util.List;

import org.genvisis.cnv.annotation.BlastAnnotationTypes.TOP_BOT;
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
	private static final String DESCRIPTION = "The probe sequence for the A probe, probe sequence for the B probe (if different than A) interrogation position,strand by design,TOB_BOTTOM SNP designation,TOP_BOTTOM Reference Designation, A allele, and B allele ";
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
		super(VCFHeaderLineType.String, null, 1, DEFAULT_NAME, DESCRIPTION, DEFUALT_VALUE, DEFUALT_VALUE);
	}

	public static MarkerSeqAnnotation getDefault() {
		return new MarkerSeqAnnotation();
	}

	public void setDesignData(String seqA, String segB, int interrogationPosition, Strand strand, TOP_BOT topBotProbe, TOP_BOT topBotRef, Allele A, Allele B) {
		this.seqA = seqA;
		this.seqB = segB;
		this.interrogationPosition = interrogationPosition;
		this.strand = strand;
		this.topBotProbe = topBotProbe;
		this.topBotRef = topBotRef;
		this.A = A;
		this.B = B;
		// this.seg=seg; populate on load only
		setData(seqA + DEFUALT_DELIMITER + segB + DEFUALT_DELIMITER + interrogationPosition + DEFUALT_DELIMITER + strand.getEncoding() + DEFUALT_DELIMITER + topBotProbe + DEFUALT_DELIMITER + topBotRef + DEFUALT_DELIMITER + A.getDisplayString() + DEFUALT_DELIMITER + B.getDisplayString());
	}

	public Allele getA() {
		return A;
	}

	public Allele getB() {
		return B;
	}

	public Allele getRef() {
		return ref;
	}

	public Allele[] getAlts() {
		return alts;
	}

	public String getSeqB() {
		return seqB;
	}

	@Override
	public void parseAnnotation(VariantContext vc, Logger log) {
		if (vc.hasAttribute(getName())) {
			setData(vc.getAttributeAsString(getName(), DEFAULT_NAME));
			if (!ext.isMissingValue(getData())) {
				List<String> data = getDataAsList();
				this.seqA = data.get(0);
				this.seqB = data.get(1);

				this.interrogationPosition = -1;
				try {
					interrogationPosition = Integer.parseInt(data.get(2));

				} catch (NumberFormatException nfe) {

				}
				this.strand = Strand.toStrand(data.get(3));
				this.seg = VCOps.getSegment(vc);
				this.topBotProbe = TOP_BOT.valueOf(data.get(4));
				this.topBotRef = TOP_BOT.valueOf(data.get(5));
				this.ref = vc.getReference();
				List<Allele> alleles = vc.getAlternateAlleles();
				this.alts = new Allele[alleles.size()];
				for (int i = 0; i < alts.length; i++) {
					alts[i] = alleles.get(i);
				}
				this.A = Allele.create(data.get(6), ref.basesMatch(data.get(6)));
				this.B = Allele.create(data.get(7), ref.basesMatch(data.get(7)));
				if (A.isReference() && B.isReference()) {
					throw new IllegalArgumentException("A and B alleles cannot both be reference");
				}
				this.indel = vc.isIndel();
				if (indel) {
					this.aInsertion = A.getBases().length > ref.getBases().length;
					this.bInsertion = B.getBases().length > ref.getBases().length;
					this.aDeletion = A.getBases().length < ref.getBases().length;
					this.bDeletion = B.getBases().length < ref.getBases().length;

					if ((aDeletion && aInsertion) || (bDeletion && bInsertion)) {
						throw new IllegalStateException("Allele cannot be both insertion and deletion");
					}
				}
			}
		}
	}

	public boolean isIndel() {
		return indel;
	}

	public boolean isaInsertion() {
		return aInsertion;
	}

	public boolean isbInsertion() {
		return bInsertion;

	}
	public boolean isaDeletion() {
		return aDeletion;
	}

	public boolean isbDeletion() {
		return bDeletion;
	}

	public String getSequence() {
		return seqA;
	}

	public Strand getStrand() {
		return strand;
	}

	public void setStrand(Strand strand2) {
        this.strand = strand2;
    }

    public int getInterrogationPosition() {
		return interrogationPosition;
	}

	public Segment getSeg() {
		return seg;
	}

	public void setSegment(Segment seg1) {
	    seg = seg1;
	}

	public boolean goLeft() {
		return topBotProbe == topBotRef;
	}

	public TOP_BOT getTopBotProbe() {
		return topBotProbe;
	}

	public TOP_BOT getTopBotRef() {
		return topBotRef;
	}

}
