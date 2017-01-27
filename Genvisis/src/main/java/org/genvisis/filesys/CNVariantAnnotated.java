/**
 * 
 */
package org.genvisis.filesys;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.genvisis.cnv.annotation.segments.SegmentAnnotationKeys;
import org.genvisis.cnv.annotation.segments.SegmentAnotation;
import org.genvisis.common.Array;


public class CNVariantAnnotated extends CNVariant {
	/**
	 * TODO, load from file methods
	 */
	private static final long serialVersionUID = 1L;

	private final SegmentAnotation segmentAnotation;


	/**
	 * @param cnv base cnv
	 * @param segmentAnotation the {@link SegmentAnotation} annotation
	 */
	public CNVariantAnnotated(CNVariant cnv, SegmentAnotation segmentAnotation) {
		super(cnv);
		this.segmentAnotation = segmentAnotation;
	}



	@Override
	public String toString() {
		return getIndividualID();
	}

	@Override
	public String toAnalysisString() {
		StringBuilder builder = new StringBuilder();
		builder.append(toPlinkFormat());
		for (SegmentAnnotationKeys key : SegmentAnnotationKeys.values()) {

			if (segmentAnotation.getAttributes().containsKey(key.toString())) {
				builder.append("\t" + segmentAnotation.getAttributes().get(key.toString()));
			} else {
				builder.append("\t" + key.getMissingValue());
			}
		}
		return builder.toString();
	}

	@Override
	public String[] getHeader() {
		List<String> header = new ArrayList<String>();
		for (String head : PLINK_CNV_HEADER) {
			header.add(head);
		}
		for (SegmentAnnotationKeys key : SegmentAnnotationKeys.values()) {
			header.add(key.toString());
		}
		return Array.toStringArray(header);
	}



	public static class TallyResult {
		private HashSet<String> allCN;
		private HashSet<String> dupCN;
		private HashSet<String> delCN;
		private HashSet<CNVariant> allLocs;

		public TallyResult() {
			super();
			this.allCN = new HashSet<String>();
			this.dupCN = new HashSet<String>();
			this.delCN = new HashSet<String>();
			this.allLocs = new HashSet<CNVariant>();
		}

		public HashSet<String> getAllCN() {
			return allCN;
		}



		public HashSet<CNVariant> getAllLocs() {
			return allLocs;
		}

		public HashSet<String> getDupCN() {
			return dupCN;
		}

		public HashSet<String> getDelCN() {
			return delCN;
		}



	}

	/**
	 * @param cList list of {@link CNVariantAnnotated} to tally
	 * @param key the annotation key to tally
	 * @return count map
	 */
	public static Map<String, TallyResult> tallyAnnotation(	List<CNVariantAnnotated> cList,
																													SegmentAnnotationKeys key) {

		Map<String, TallyResult> tally = new HashMap<String, TallyResult>();
		for (CNVariantAnnotated cnVarAn : cList) {
			if (cnVarAn.segmentAnotation.getAttributes().containsKey(key.toString())) {
				List<String> anns = cnVarAn.segmentAnotation.getAttributes().get(key.toString());
				for (String ann : anns) {
					if (!tally.containsKey(ann)) {
						tally.put(ann, new TallyResult());


					}

					tally.get(ann).getAllCN().add(cnVarAn.getIndividualID());
					tally.get(ann).getAllLocs().add(cnVarAn);
					if (cnVarAn.getCN() < 2) {
						tally.get(ann).getDelCN().add(cnVarAn.getIndividualID());

					} else if (cnVarAn.getCN() > 2) {
						tally.get(ann).getDupCN().add(cnVarAn.getIndividualID());

					}
				}
			}
		}
		return tally;
	}

}
