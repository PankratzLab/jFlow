package org.genvisis.cnv.annotation;

import java.util.List;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Logger;
import org.genvisis.stats.Histogram.DynamicHistogram;

import htsjdk.variant.variantcontext.VariantContext;

public class MarkerBlastHistogramAnnotation extends HistogramAnnotation {

	public static final String DEFAULT_NAME = "ALIGNMENT_HISTOGRAM";
	public static final String DEFAULT_DESCRIPTION = "A histogam of blast alignment lengths, alignment lengths include gapped and mismatched alignments. Histogram is reported from the smallest alignment seen to the probe length";

	private List<String> blastAlignmentCounts;

	public MarkerBlastHistogramAnnotation(String name, String description,
																				DynamicHistogram dynamicHistogram) {
		super(name, description, dynamicHistogram);

	}

	public static Annotation getDefaultBlastAnnotation() {
		return new HistogramAnnotation(DEFAULT_NAME, DEFAULT_DESCRIPTION, null) {};
	}

	@Override
	public void parseAnnotation(VariantContext vc, Logger log) {
		// System.out.println(getName());
		if (vc.hasAttribute(getName())) {
			setData(vc.getAttributeAsString(getName(), DEFAULT_NAME));
			blastAlignmentCounts = getDataAsList();
			setFound(true);
			// System.out.println(blastAlignmentCounts.toString());
			// System.out.println(vc.toStringWithoutGenotypes());

		}
	}

	public List<String> getBlastAlignmentCounts() {
		return blastAlignmentCounts;
	}

	public void setBlastAlignmentCounts(List<String> blastAlignmentCounts) {
		this.blastAlignmentCounts = blastAlignmentCounts;
	}

	/**
	 * @param proj
	 * @return an int array representing the counts at each alignment length
	 */
	public int[] formatCountsForProject(Project proj) {
		int[] counts = new int[proj.getArrayType().getProbeLength()];
		if (blastAlignmentCounts != null && blastAlignmentCounts.size() > 0) {
			if (blastAlignmentCounts.size() > 1 || !blastAlignmentCounts.get(0).equals(DEFUALT_VALUE)) {
				if (blastAlignmentCounts.size() > proj.getArrayType().getProbeLength()) {
					proj.getLog()
							.reportTimeError("Aligment counts had more entries ("	+ blastAlignmentCounts.size()
																+ " than the projects probe size ("
																+ proj.getArrayType().getProbeLength() + ")");
				} else {
					int countIndex = counts.length - 1;
					for (int i = blastAlignmentCounts.size() - 1; i >= 0; i--) {
						counts[countIndex] = Integer.parseInt(blastAlignmentCounts.get(i));
						countIndex--;
					}
				}
			}
		}
		return counts;
	}

}
