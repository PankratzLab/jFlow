package org.genvisis.cnv.annotation.markers;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.stats.Histogram.DynamicHistogram;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * TODO, gotta be a better way to track without the string conversions
 *
 *
 */
public class MarkerEvalueHistogramAnnotation extends HistogramAnnotation {

	public static final String DEFAULT_NAME = "BLAST_EVALUE_COUNTS";
	public static final String DEFUALT_DELIM = "/";
	public static final String DEFAULT_DESCRIPTION = "A list of evalues and counts for all blast hits a marker";
	private List<String> blastEvalueCounts;
	private List<String> blastEvalues;
	private final Hashtable<String, Integer> exactHistogram;

	public MarkerEvalueHistogramAnnotation(String name, String description) {
		super(name, description, new DynamicHistogram(0, 1, 2));// the Dynamic histogram is not actualy
																														// used
		exactHistogram = new Hashtable<String, Integer>();
	}

	public void setDataToExactHistogram() {
		String[] numericKeys = HashVec.getNumericKeys(exactHistogram);

		ArrayList<String> data = new ArrayList<String>();
		for (String key : numericKeys) {
			data.add(key);
			data.add(exactHistogram.get(key).toString());
		}
		setData(Array.toStr(data, DEFUALT_DELIMITER));
	}

	public void addExactHistogram(double d) {
		String tmp = d + "";
		if (exactHistogram.containsKey(tmp)) {
			exactHistogram.put(tmp, exactHistogram.get(tmp) + 1);
		} else {
			exactHistogram.put(tmp, 1);
		}
	}

	public static Annotation getDefaultBlastAnnotation() {
		return new HistogramAnnotation(DEFAULT_NAME, DEFAULT_DESCRIPTION, null) {};
	}

	@Override
	public void parseAnnotation(VariantContext vc, Logger log) {
		if (vc.hasAttribute(getName())) {
			setData(vc.getAttributeAsString(getName(), DEFAULT_NAME));
			List<String> tmp = getDataAsList();
			if (!tmp.get(0).equals(".")) {// skip no alignments
				// System.out.println(tmp.toString());
				ArrayList<String> tmpCounts = new ArrayList<String>();
				ArrayList<String> tmpEvals = new ArrayList<String>();
				for (int i = 0; i < tmp.size(); i += 2) {
					tmpEvals.add(tmp.get(i));
					tmpCounts.add(tmp.get(i + 1));
				}
				blastEvalueCounts = tmpCounts;
				blastEvalues = tmpEvals;
			}
			setFound(true);

		}
	}

	/**
	 * Not really a histogram since the values are not binned, but thats cuz I didn't know the best
	 * way to bin
	 *
	 */
	public static class EvalueHistogram {
		private final double[] evalues;
		private final int[] counts;

		public EvalueHistogram(double[] evalues, int[] counts) {
			super();
			this.evalues = evalues;
			this.counts = counts;
		}

		public int[] getCounts() {
			return counts;
		}

		public double[] getEvalues() {
			return evalues;
		}

	}

	public EvalueHistogram formatHistogram() {
		if (blastEvalueCounts == null	|| blastEvalues == null
				|| blastEvalueCounts.size() != blastEvalues.size()) {
			throw new IllegalStateException("Invalid sizes for counts and evalue arrays");

		} else {
			double[] evalues = new double[blastEvalues.size()];
			int[] counts = new int[blastEvalueCounts.size()];
			for (int i = 0; i < counts.length; i++) {
				counts[i] = Integer.parseInt(blastEvalueCounts.get(i));
				evalues[i] = Double.parseDouble(blastEvalues.get(i));
			}
			return new EvalueHistogram(evalues, counts);
		}

	}

	public List<String> getBlastEvalues() {
		return blastEvalues;
	}

	public List<String> getBlastEvalueCounts() {
		return blastEvalueCounts;
	}

}
