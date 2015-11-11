package cnv.annotation;

import java.util.List;

import cnv.annotation.BlastAnnotationTypes.TOP_BOT;
import seq.manage.VCOps;
import common.Logger;
import filesys.Segment;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class MarkerSeqAnnotation extends AnnotationData {

	private static final String DEFAULT_NAME = "PROBE_DESIGN";
	private static final String DESCRIPTION = "The probe sequence, interrogation position,TOB_BOTTOM SNP designation,TOP_BOTTOM Reference Designation and strand by design";
	private String sequence;
	private Strand strand;
	private int interrogationPosition;
	private Segment seg;
	private TOP_BOT topBotProbe;
	private TOP_BOT topBotRef;

	public MarkerSeqAnnotation() {
		super(VCFHeaderLineType.String, null, 1, DEFAULT_NAME, DESCRIPTION, DEFUALT_VALUE, DEFUALT_VALUE);
	}

	public static MarkerSeqAnnotation getDefault() {
		return new MarkerSeqAnnotation();
	}

	public void setDesignData(String sequence, int interrogationPosition, Strand strand, TOP_BOT topBotProbe, TOP_BOT topBotRef) {
		this.sequence = sequence;
		this.interrogationPosition = interrogationPosition;
		this.strand = strand;
		this.topBotProbe=topBotProbe;
		this.topBotRef =topBotRef;
		//this.seg=seg; populate on load only
		setData(sequence + DEFUALT_DELIMITER + interrogationPosition + DEFUALT_DELIMITER + strand.getEncoding() + DEFUALT_DELIMITER + topBotProbe + DEFUALT_DELIMITER + topBotRef);
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
