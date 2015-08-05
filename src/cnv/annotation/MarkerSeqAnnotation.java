package cnv.annotation;

import java.util.List;

import common.Logger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class MarkerSeqAnnotation extends AnnotationData {

	private static final String DEFAULT_NAME = "PROBE_DESIGN";
	private static final String DESCRIPTION = "The probe sequence and interrogation position by design";
	private String sequence;
	private int interrogationPosition;

	public MarkerSeqAnnotation() {
		super(VCFHeaderLineType.String, null, 1, DEFAULT_NAME, DESCRIPTION, DEFUALT_VALUE, DEFUALT_VALUE);
	}

	public static MarkerSeqAnnotation getDefault() {
		return new MarkerSeqAnnotation();
	}

	public void setDesignData(String sequence, int interrogationPosition) {
		this.sequence = sequence;
		this.interrogationPosition = interrogationPosition;
		setData(sequence + DEFUALT_DELIMITER + interrogationPosition);
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
		}
	}

	public String getSequence() {
		return sequence;
	}

	public int getInterrogationPosition() {
		return interrogationPosition;
	}

}
