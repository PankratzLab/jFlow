/**
 * 
 */
package org.genvisis.cnv.analysis.pod;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.genvisis.cnv.analysis.pod.InformativeBAF.BAF_STRATEGY;
import org.genvisis.cnv.analysis.pod.InformativeBAF.InformativeResult;
import org.genvisis.cnv.analysis.pod.PODAnalysis.PODResults.Builder;
import org.genvisis.cnv.analysis.pod.PODGenotype.BAF_EFFECT;
import org.genvisis.cnv.analysis.pod.PODGenotype.Genotype;
import org.genvisis.cnv.analysis.pod.PODGenotype.POD;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.ArrayUtils;
import org.genvisis.filesys.Segment;


/**
 * 
 * Methods to determine parent of origin effects from BAFs/Genotypes
 */
public class PODAnalysis {


	private enum COMPLETION_STATUS {
		FULL, MO_ONLY, FA_ONLY
	}

	private static class Trio {
		private SubSample off;
		private SubSample mo;
		private SubSample fa;
		private BAF_EFFECT[] offEffects;

		private Trio(SubSample off, SubSample mo, SubSample fa, BAF_EFFECT[] offEffects) {
			super();
			this.off = off;
			this.mo = mo;
			this.fa = fa;
			this.offEffects = offEffects;
		}
	}

	private static class SubSample {
		private byte[] genos;

		private SubSample(byte[] genos, int[] subsetTo) {
			super();
			this.genos = ArrayUtils.subArray(genos, subsetTo);
		}

		private byte[] getGenos() {
			return genos;
		}
	}

	/**
	 * Results from parent of origin detection
	 *
	 */
	public static class PODResults {
		private String offDNA;
		private String moDNA;
		private String faDNA;
		private List<Segment> searchSpace;
		private int numtotal;
		private int numHet;
		private int numBAFInformative;
		private int numSNPInformative;
		private int numPODMaternal;
		private int numPODPaternal;
		private int numNONE;
		private int numMaternalAlleleMatch;
		private int numMaternalExactMatch;

		private int numPaternalAlleleMatch;
		private int numPaternalExactMatch;

		private COMPLETION_STATUS completionStatus;

		private PODResults(Builder builder) {
			this.offDNA = builder.offDNA;
			this.moDNA = builder.moDNA;
			this.faDNA = builder.faDNA;
			this.searchSpace = builder.searchSpace;
			this.numtotal = builder.numtotal;
			this.numHet = builder.numHet;
			this.numBAFInformative = builder.numBAFInformative;
			this.numSNPInformative = builder.numSNPInformative;
			this.numPODMaternal = builder.numPODMaternal;
			this.numPODPaternal = builder.numPODPaternal;
			this.numNONE = builder.numNONE;
			this.numMaternalAlleleMatch = builder.numMaternalAlleleMatch;
			this.numMaternalExactMatch = builder.numMaternalExactMatch;
			this.numPaternalAlleleMatch = builder.numPaternalAlleleMatch;
			this.numPaternalExactMatch = builder.numPaternalExactMatch;
			this.completionStatus = builder.completionStatus;
		}

		private double getProportionPOD(POD pod) {
			if (numSNPInformative > 0) {
				switch (pod) {
					case MATERNAL:
						return (double) numPODMaternal / numSNPInformative;
					case PATERNAL:
						return (double) numPODPaternal / numSNPInformative;
					case NONE:
					default:
						throw new IllegalArgumentException();

				}
			} else {
				return Double.NaN;
			}
		}

		String getSearchSpaceString() {
			StringBuilder stringBuilder = new StringBuilder(searchSpace.get(0).getUCSClocation());
			for (int i = 1; i < searchSpace.size(); i++) {
				stringBuilder.append(";" + searchSpace.get(i).getUCSClocation());
			}
			return stringBuilder.toString();
		}

		public static String[] getHeader() {
			return new String[] {"OFFSPRING_DNA", "MO_DNA", "FA_DNA", "SEARCH_SPACE", "NUM_TOTAL",
													 "NUM_HET",
													 "numBAFInformative",
													 "numSNPInformative", "numPODMaternal", "numPODPaternal", "numNONE",
													 "proportionMaternal", "proportionPaternal",
													 "numMaternalAlleleMatch", "numMaternalExactMatch",
													 "numPaternalAlleleMatch", "numPaternalExactMatch",
													 "completionStatus"

			};
		}

		@Override
		public String toString() {
			StringJoiner joiner = new StringJoiner("\t");
			joiner.add(offDNA);
			joiner.add(moDNA);
			joiner.add(faDNA);
			joiner.add(getSearchSpaceString());
			joiner.add(Integer.toString(numtotal));
			joiner.add(Integer.toString(numHet));
			joiner.add(Integer.toString(numBAFInformative));
			joiner.add(Integer.toString(numSNPInformative));
			joiner.add(Integer.toString(numPODMaternal));
			joiner.add(Integer.toString(numPODPaternal));
			joiner.add(Integer.toString(numNONE));
			joiner.add(Double.toString(getProportionPOD(POD.MATERNAL)));
			joiner.add(Double.toString(getProportionPOD(POD.PATERNAL)));
			joiner.add(Integer.toString(numMaternalAlleleMatch));
			joiner.add(Integer.toString(numMaternalExactMatch));
			joiner.add(Integer.toString(numPaternalAlleleMatch));
			joiner.add(Integer.toString(numPaternalExactMatch));
			joiner.add(completionStatus.toString());
			return joiner.toString();
		}


		/**
		 * Creates builder to build {@link PODResults}.
		 * 
		 * @return created builder
		 */
		private static Builder builder() {
			return new Builder();
		}

		/**
		 * Builder to build {@link PODResults}.
		 */
		static class Builder {
			private String offDNA;
			private String moDNA = ".";
			private String faDNA = ".";
			private List<Segment> searchSpace = new ArrayList<>();
			private int numtotal = -1;
			private int numHet = -1;
			private int numBAFInformative = -1;
			private int numSNPInformative = -1;
			private int numPODMaternal = -1;
			private int numPODPaternal = -1;
			private int numNONE = -1;
			private int numMaternalAlleleMatch = -1;
			private int numMaternalExactMatch = -1;
			private int numPaternalAlleleMatch = -1;
			private int numPaternalExactMatch = -1;
			private COMPLETION_STATUS completionStatus;

			public Builder() {
				// dont need init params
			}

			private void withOffDNA(String offDNA) {
				this.offDNA = offDNA;
			}

			private void withMoDNA(String moDNA) {
				this.moDNA = moDNA;
			}

			private void withFaDNA(String faDNA) {
				this.faDNA = faDNA;
			}

			private void withSearchSpace(List<Segment> searchSpace) {
				this.searchSpace = searchSpace;
			}

			private void withNumtotal(int numtotal) {
				this.numtotal = numtotal;
			}

			private void withNumHet(int numHet) {
				this.numHet = numHet;
			}

			private void withNumBAFInformative(int numBAFInformative) {
				this.numBAFInformative = numBAFInformative;
			}

			private void withNumSNPInformative(int numSNPInformative) {
				this.numSNPInformative = numSNPInformative;
			}

			private void withNumPODMaternal(int numPODMaternal) {
				this.numPODMaternal = numPODMaternal;
			}

			private void withPODNumPaternal(int numPaternal) {
				this.numPODPaternal = numPaternal;
			}

			private void withNumNONE(int numNONE) {
				this.numNONE = numNONE;
			}

			private void withNumMaternalAlleleMatch(int numMaternalAlleleMatch) {
				this.numMaternalAlleleMatch = numMaternalAlleleMatch;
			}

			private void withNumMaternalExactMatch(int numMaternalExactMatch) {
				this.numMaternalExactMatch = numMaternalExactMatch;
			}

			private void withNumPaternalAlleleMatch(int numPaternalAlleleMatch) {
				this.numPaternalAlleleMatch = numPaternalAlleleMatch;
			}

			private void withNumPaternalExactMatch(int numPaternalExactMatch) {
				this.numPaternalExactMatch = numPaternalExactMatch;
			}

			private void withCompletionStatus(COMPLETION_STATUS completionStatus) {
				this.completionStatus = completionStatus;
			}

			private PODResults build() {
				return new PODResults(this);
			}
		}
	}



	/**
	 * How segments will be searched
	 *
	 */
	public enum SEARCH_SPACE_TYPE {
		/**
		 * Include all markers in this searchSpace, i.e union of all segments
		 */
		INCLUSIVE,
		/**
		 * Run for each individual segment
		 */
		INDIVIDUAL;
	}

	/**
	 * @param proj
	 * @param markerDetailSet
	 * @param segmentSearchSpaces
	 * @param offDNA
	 * @param moDNA
	 * @param faDNA
	 * @param type how the segmentSearchSpaces will be searched
	 * @return List of {@link PODResults} for every {@link Segment} in segmentSearchSpaces according
	 *         to the {@link SEARCH_SPACE_TYPE}
	 */
	public static List<PODResults> analyze(Project proj, MarkerDetailSet markerDetailSet,
																				 List<Segment> segmentSearchSpaces,

																				 String offDNA, String moDNA, String faDNA,
																				 SEARCH_SPACE_TYPE type) {

		proj.getLog().reportTimeInfo("POD for " + offDNA + " with type " + type);

		Sample off = proj.getFullSampleFromRandomAccessFile(offDNA);
		Sample mo = moDNA == null ? null : proj.getFullSampleFromRandomAccessFile(moDNA);
		Sample fa = faDNA == null ? null : proj.getFullSampleFromRandomAccessFile(faDNA);
		List<PODResults> results = new ArrayList<>();

		double[] bafs = ArrayUtils.toDoubleArray(off.getBAFs());
		InformativeResult informativeResult = InformativeBAF.getInformativeIndices(bafs,
																																							 off.getAB_Genotypes(),
																																							 BAF_STRATEGY.HET_ONLY,
																																							 InformativeBAF.CHEBYSHEV);
		if (type == SEARCH_SPACE_TYPE.INDIVIDUAL) {
			for (Segment searchSpace : segmentSearchSpaces) {

				PODResults podResults = computeForSegment(proj, markerDetailSet, offDNA, moDNA, faDNA, off,
																									mo, fa,
																									searchSpace, informativeResult);
				results.add(podResults);
			}
		} else if (type == SEARCH_SPACE_TYPE.INCLUSIVE) {
			List<Integer> indices = new ArrayList<>();
			for (Segment searchSpace : segmentSearchSpaces) {
				@SuppressWarnings("deprecation")
				int[] tmp = markerDetailSet.getIndicesOfMarkersIn(searchSpace, null,
																													proj.getLog());
				for (int i = 0; i < tmp.length; i++) {
					indices.add(tmp[i]);
				}
			}
			int[] markerIndicesToUse = indices.stream().mapToInt(i -> i).toArray();
			Builder builder = computeForIndices(offDNA, moDNA, faDNA, off, mo, fa,
																					markerIndicesToUse, informativeResult);
			builder.withSearchSpace(segmentSearchSpaces);
			results.add(builder.build());

		}
		return results;
	}

	private static PODResults computeForSegment(Project proj, MarkerDetailSet markerDetailSet,
																							String offDNA, String moDNA, String faDNA, Sample off,
																							Sample mo, Sample fa,
																							Segment searchSpace,
																							InformativeResult informativeResult) {
		@SuppressWarnings("deprecation")
		int[] markerIndicesToUse = markerDetailSet.getIndicesOfMarkersIn(searchSpace, null,
																																		 proj.getLog());

		Builder builder = computeForIndices(offDNA, moDNA, faDNA, off, mo, fa,
																				markerIndicesToUse, informativeResult);
		builder.withSearchSpace(Arrays.asList(new Segment[] {searchSpace}));

		return builder.build();
	}

	private static Builder computeForIndices(String offDNA, String moDNA, String faDNA, Sample off,
																					 Sample mo, Sample fa,
																					 int[] markerIndicesToUse,
																					 InformativeResult informativeResult) {
		Builder builder = PODResults.builder();
		builder.withOffDNA(offDNA);

		builder.withNumHet(ArrayUtils.countIf(ArrayUtils.toIntArray(ArrayUtils.subArray(off.getAB_Genotypes(),
																																										markerIndicesToUse)),
																					1));
		builder.withNumtotal(markerIndicesToUse.length);
		if (moDNA != null && faDNA != null) {
			bafGenotypePOD(markerIndicesToUse, builder, off, mo, fa, informativeResult);
		}
		SubSample offSub = new SubSample(off.getAB_Genotypes(),
																		 markerIndicesToUse);
		if (moDNA != null) {
			SubSample moSub = new SubSample(mo.getAB_Genotypes(),
																			markerIndicesToUse);
			GenoCompResult genoCompResult = compGenos(offSub, moSub);
			builder.withNumMaternalExactMatch(genoCompResult.numExactMatch);
			builder.withNumMaternalAlleleMatch(genoCompResult.numAlleleMatch);
			if (faDNA == null) {
				builder.withCompletionStatus(COMPLETION_STATUS.MO_ONLY);
			}
			builder.withMoDNA(moDNA);

		}
		if (faDNA != null) {
			SubSample faSub = new SubSample(fa.getAB_Genotypes(),
																			markerIndicesToUse);
			GenoCompResult genoCompResult = compGenos(offSub, faSub);
			builder.withNumPaternalExactMatch(genoCompResult.numExactMatch);
			builder.withNumPaternalAlleleMatch(genoCompResult.numAlleleMatch);
			if (moDNA == null) {
				builder.withCompletionStatus(COMPLETION_STATUS.FA_ONLY);
			}
			builder.withFaDNA(faDNA);
		}
		return builder;
	}

	/**
	 * Determine POD ala https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3680018/ using BAFs and
	 * genotypes
	 */
	static InformativeResult bafGenotypePOD(int[] markerIndicesToUse, Builder builder,
																					Sample off,
																					Sample mo, Sample fa,
																					InformativeResult informativeResult) {
		builder.withCompletionStatus(COMPLETION_STATUS.FULL);

		List<Integer> informativeSpecific = new ArrayList<>();
		List<BAF_EFFECT> effectsSpecific = new ArrayList<>();

		for (int i = 0; i < markerIndicesToUse.length; i++) {

			if (informativeResult.getInformatives().contains(markerIndicesToUse[i])) {
				informativeSpecific.add(markerIndicesToUse[i]);
				effectsSpecific.add(informativeResult.getEffects()
																						 .get(informativeResult.getInformatives()
																																	 .indexOf(markerIndicesToUse[i])));
			}

		}
		builder.withNumBAFInformative(informativeSpecific.size());
		int[] subset = informativeSpecific.stream().mapToInt(i -> i).toArray();


		SubSample offSub = new SubSample(off.getAB_Genotypes(),
																		 subset);
		SubSample moSub = new SubSample(mo.getAB_Genotypes(),
																		subset);
		SubSample faSub = new SubSample(fa.getAB_Genotypes(),
																		subset);

		Trio trio = new Trio(offSub, moSub, faSub,
												 effectsSpecific.toArray(new BAF_EFFECT[effectsSpecific.size()]));

		List<POD> pods = PODAnalysis.analyze(trio);
		informativeResult.setPods(pods);
		Map<POD, Long> result = pods.stream().collect(
																									Collectors.groupingBy(
																																				Function.identity(),
																																				Collectors.counting()));

		for (POD pod : POD.values()) {
			int count = (int) (result.containsKey(pod) ? result.get(pod) : 0);
			switch (pod) {
				case MATERNAL:
					builder.withNumPODMaternal(count);
					break;
				case NONE:
					builder.withNumNONE(count);
					break;
				case PATERNAL:
					builder.withPODNumPaternal(count);
					break;
				default:
					break;

			}
		}
		builder.withNumSNPInformative(builder.numPODMaternal + builder.numPODPaternal);
		return informativeResult;
	}

	static class GenoCompResult {
		int numExactMatch;
		int numAlleleMatch;

		GenoCompResult(int numExactMatch, int numAlleleMatch) {
			super();
			this.numExactMatch = numExactMatch;
			this.numAlleleMatch = numAlleleMatch;
		}

		int getNumExactMatch() {
			return numExactMatch;
		}

		int getNumAlleleMatch() {
			return numAlleleMatch;
		}

	}



	private static GenoCompResult compGenos(SubSample off, SubSample p) {
		int numExactMatch = 0;
		int numAlleleMatch = 0;
		for (int i = 0; i < off.getGenos().length; i++) {
			GenoCompResult tmp = Genotype.getSharedAlleleCount(Genotype.fromByte(off.getGenos()[i]),
																												 Genotype.fromByte(p.getGenos()[i]));
			numExactMatch += tmp.numExactMatch;
			numAlleleMatch += tmp.numAlleleMatch;
		}

		return new GenoCompResult(numExactMatch, numAlleleMatch);
	}



	private static List<POD> analyze(Trio trio) {

		List<POD> pods = new ArrayList<>();
		for (int i = 0; i < trio.off.genos.length; i++) {
			pods.add(PODGenotype.getPodEffect(Genotype.fromByte(trio.mo.genos[i]),
																				Genotype.fromByte(trio.fa.genos[i]),
																				trio.offEffects[i]));
		}
		return pods;
	}

}
