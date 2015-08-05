package cnv.annotation;

import java.util.List;

import common.Logger;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class MarkerSeqAnnotation extends AnnotationData {

	private static final String DEFAULT_NAME = "PROBE_DESIGN";
	private static final String DESCRIPTION = "The probe sequence and interrogation position by design";
	private String sequence;
	private Strand strand;
	private int interrogationPosition;

	public MarkerSeqAnnotation() {
		super(VCFHeaderLineType.String, null, 1, DEFAULT_NAME, DESCRIPTION, DEFUALT_VALUE, DEFUALT_VALUE);
	}

	public static MarkerSeqAnnotation getDefault() {
		return new MarkerSeqAnnotation();
	}

	public void setDesignData(String sequence, int interrogationPosition, Strand strand) {
		this.sequence = sequence;
		this.interrogationPosition = interrogationPosition;
		this.strand = strand;
		setData(sequence + DEFUALT_DELIMITER + interrogationPosition + DEFUALT_DELIMITER + strand);
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
			this.strand = Strand.valueOf(data.get(2));
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

}
