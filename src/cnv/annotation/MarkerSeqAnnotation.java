package cnv.annotation;

import java.util.List;

import cnv.annotation.BlastAnnotationTypes.TOP_BOT;
import seq.manage.VCOps;
import common.Logger;
import filesys.Segment;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class MarkerSeqAnnotation extends AnnotationData {

	private static final String DEFAULT_NAME = "PROBE_DESIGN";
	private static final String DESCRIPTION = "The probe sequence, interrogation position,strand by design,TOB_BOTTOM SNP designation,TOP_BOTTOM Reference Designation, A allele, and B allele ";
	private String sequence;
	private Strand strand;
	private int interrogationPosition;
	private Segment seg;
	private TOP_BOT topBotProbe;
	private TOP_BOT topBotRef;
	private Allele A;
	private Allele B;
	private Allele ref;
	private Allele[] alts;

	public MarkerSeqAnnotation() {
		super(VCFHeaderLineType.String, null, 1, DEFAULT_NAME, DESCRIPTION, DEFUALT_VALUE, DEFUALT_VALUE);
	}

	public static MarkerSeqAnnotation getDefault() {
		return new MarkerSeqAnnotation();
	}

	public void setDesignData(String sequence, int interrogationPosition, Strand strand, TOP_BOT topBotProbe, TOP_BOT topBotRef, Allele A, Allele B) {
		this.sequence = sequence;
		this.interrogationPosition = interrogationPosition;
		this.strand = strand;
		this.topBotProbe = topBotProbe;
		this.topBotRef = topBotRef;
		this.A = A;
		this.B = B;
		// this.seg=seg; populate on load only
		setData(sequence + DEFUALT_DELIMITER + interrogationPosition + DEFUALT_DELIMITER + strand.getEncoding() + DEFUALT_DELIMITER + topBotProbe + DEFUALT_DELIMITER + topBotRef + DEFUALT_DELIMITER + A.getDisplayString() + DEFUALT_DELIMITER + B.getDisplayString());
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

	@Override
	public void parseAnnotation(VariantContext vc, Logger log) {
		if (vc.hasAttribute(getName())) {
			setData(vc.getAttributeAsString(getName(), DEFAULT_NAME));
			List<String> data = getDataAsList();
			this.sequence = data.get(0);
			this.interrogationPosition = -1;
			try {
				interrogationPosition = Integer.parseInt(data.get(1));

			} catch (NumberFormatException nfe) {

			}
			this.strand = Strand.toStrand(data.get(2));
			this.seg = VCOps.getSegment(vc);
			this.topBotProbe = TOP_BOT.valueOf(data.get(3));
			this.topBotRef = TOP_BOT.valueOf(data.get(4));
			this.ref = vc.getReference();
			List<Allele> alleles = vc.getAlternateAlleles();
			this.alts = new Allele[alleles.size()];
			for (int i = 0; i < alts.length; i++) {
				alts[i] = alleles.get(i);
			}
			this.A = Allele.create(data.get(5), ref.basesMatch(data.get(5)));
			this.B = Allele.create(data.get(6), ref.basesMatch(data.get(6)));
			if (A.isReference() && B.isReference()) {
				throw new IllegalArgumentException("A and B alleles cannot both be reference");
			}
		}
	}

	public String getSequence() {
		return sequence;
	}

	public Strand getStrand() {
		return strand;
	}

	public int getInterrogationPosition() {
		return interrogationPosition;
	}

	public Segment getSeg() {
		return seg;
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
