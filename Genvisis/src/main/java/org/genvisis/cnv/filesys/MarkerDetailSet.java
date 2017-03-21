package org.genvisis.cnv.filesys;

import java.io.PrintWriter;
import java.io.Serializable;
import java.lang.ref.Reference;
import java.lang.ref.SoftReference;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.cnv.annotation.markers.AnnotationFileLoader.QUERY_TYPE;
import org.genvisis.cnv.annotation.markers.AnnotationParser;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.markers.MarkerAnnotationLoader;
import org.genvisis.cnv.annotation.markers.MarkerBlastAnnotation;
import org.genvisis.cnv.annotation.markers.MarkerSeqAnnotation;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker.GenomicPosition;
import org.genvisis.cnv.manage.TextExport;
import org.genvisis.cnv.util.Java6Helper;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.filesys.Segment;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;
import com.google.common.collect.SetMultimap;
import com.google.common.collect.Sets;

import htsjdk.variant.variantcontext.Allele;

public class MarkerDetailSet implements MarkerSetInfo, Serializable, TextExport {

	public static class Marker implements Serializable, Comparable<Marker> {

		public enum RefAllele {
			A, B
		}

		public static class GenomicPosition implements Serializable, Comparable<GenomicPosition> {

			private static final long serialVersionUID = 1L;

			private final byte chr;
			private final int position;

			/**
			 * @param chr
			 * @param position
			 */
			public GenomicPosition(byte chr, int position) {
				super();
				this.chr = chr;
				this.position = position;
			}

			public byte getChr() {
				return chr;
			}

			public int getPosition() {
				return position;
			}

			@Override
			public int hashCode() {
				final int prime = 31;
				int result = 1;
				result = prime * result + chr;
				result = prime * result + position;
				return result;
			}

			@Override
			public boolean equals(Object obj) {
				if (this == obj)
					return true;
				if (obj == null)
					return false;
				if (getClass() != obj.getClass())
					return false;
				GenomicPosition other = (GenomicPosition) obj;
				if (chr != other.chr)
					return false;
				if (position != other.position)
					return false;
				return true;
			}

			@Override
			public int compareTo(GenomicPosition o) {
				int cmp = Java6Helper.compare(chr, o.chr);
				if (cmp != 0)
					return cmp;
				cmp = Java6Helper.compare(position, o.position);
				return cmp;
			}



		}

		private static final long serialVersionUID = 2L;

		private final String name;
		private final GenomicPosition genomicPosition;
		private final char a;
		private final char b;
		private final RefAllele refAllele;


		/**
		 * @param name
		 * @param chr
		 * @param position
		 * @param a
		 * @param b
		 * @param refAllele
		 */
		private Marker(String name, byte chr, int position, char a, char b, RefAllele refAllele) {
			super();
			this.name = name;
			this.genomicPosition = new GenomicPosition(chr, position);
			this.a = a;
			this.b = b;
			this.refAllele = refAllele;
		}

		public String getName() {
			return name;
		}

		public byte getChr() {
			return genomicPosition.getChr();
		}

		public int getPosition() {
			return genomicPosition.getPosition();
		}

		public GenomicPosition getGenomicPosition() {
			return genomicPosition;
		}

		public char getA() {
			return a;
		}

		public char getB() {
			return b;
		}

		public RefAllele getRefAllele() {
			return refAllele;
		}

		public char[] getAB() {
			return new char[] {a, b};
		}

		public char getRef() {
			if (refAllele == null) {
				return ABLookup.MISSING_ALLELE;
			}
			switch (refAllele) {
				case A:
					return a;
				case B:
					return b;
				default:
					throw new IllegalStateException("Undefined refAllele instance:" + refAllele.toString());

			}
		}

		public char getAlt() {
			if (refAllele == null) {
				return ABLookup.MISSING_ALLELE;
			}
			switch (refAllele) {
				case A:
					return b;
				case B:
					return a;
				default:
					throw new IllegalStateException("Undefined refAllele instance:" + refAllele.toString());

			}
		}



		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + a;
			result = prime * result + b;
			result = prime * result + ((genomicPosition == null) ? 0 : genomicPosition.hashCode());
			result = prime * result + ((name == null) ? 0 : name.hashCode());
			result = prime * result + ((refAllele == null) ? 0 : refAllele.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Marker other = (Marker) obj;
			if (a != other.a)
				return false;
			if (b != other.b)
				return false;
			if (genomicPosition == null) {
				if (other.genomicPosition != null)
					return false;
			} else if (!genomicPosition.equals(other.genomicPosition))
				return false;
			if (name == null) {
				if (other.name != null)
					return false;
			} else if (!name.equals(other.name))
				return false;
			if (refAllele != other.refAllele)
				return false;
			return true;
		}

		@Override
		public int compareTo(Marker o) {
			int cmp = genomicPosition.compareTo(o.genomicPosition);
			if (cmp != 0)
				return cmp;
			cmp = name.compareTo(o.name);
			if (cmp != 0)
				return cmp;
			cmp = Java6Helper.compare(a, o.a);
			if (cmp != 0)
				return cmp;
			cmp = Java6Helper.compare(b, o.b);
			return cmp;
		}
	}

	public static final long serialVersionUID = 1L;

	private final ImmutableList<Marker> markers;
	private final int hashCode;
	private final long markerSetFingerprint;

	private transient Reference<Map<String, Marker>> markerNameMapRef = null;
	private transient Reference<ListMultimap<Byte, Marker>> chrMapRef = null;
	private transient Reference<SetMultimap<GenomicPosition, Marker>> genomicPositionMapRef = null;
	private transient Reference<Map<Marker, Integer>> markerIndexMapRef = null;

	public MarkerDetailSet(MarkerSetInfo markerSet) {
		this(markerSet.getMarkerNames(), markerSet.getChrs(), markerSet.getPositions(),
				 markerSet.getFingerprint());
	}

	public MarkerDetailSet(MarkerSetInfo markerSet, char[][] abAlleles) {
		this(markerSet.getMarkerNames(), markerSet.getChrs(), markerSet.getPositions(), abAlleles,
				 markerSet.getFingerprint());
	}

	public MarkerDetailSet(String[] markerNames, byte[] chrs, int[] positions) {
		this(markerNames, chrs, positions, null);
	}

	public MarkerDetailSet(String[] markerNames, byte[] chrs, int[] positions,
												 long markerSetFingerprint) {
		this(markerNames, chrs, positions, null, markerSetFingerprint);
	}

	@SuppressWarnings("deprecation")
	public MarkerDetailSet(String[] markerNames, byte[] chrs, int[] positions, char[][] abAlleles) {
		this(markerNames, chrs, positions, abAlleles, MarkerSet.fingerprint(markerNames));
	}

	public MarkerDetailSet(String[] markerNames, byte[] chrs, int[] positions, char[][] abAlleles,
												 long markerSetFingerprint) {
		super();
		int numMarkers = markerNames.length;
		if (numMarkers != chrs.length || numMarkers != positions.length
				|| (abAlleles != null && numMarkers != abAlleles.length)) {
			throw new IllegalArgumentException(this.getClass().getName()
																				 + " cannot be constructed with mismatched list of Markers and positions or AB Alleles");
		}
		ImmutableList.Builder<Marker> markersBuilder = ImmutableList.builder();

		for (int i = 0; i < markerNames.length; i++) {
			char a;
			char b;
			if (abAlleles != null) {
				a = abAlleles[i][0];
				b = abAlleles[i][1];
			} else {
				a = 'A';
				b = 'B';
			}
			markersBuilder.add(new Marker(markerNames[i], chrs[i], positions[i], a, b, null));
		}
		markers = markersBuilder.build();
		this.markerSetFingerprint = markerSetFingerprint;
		hashCode = generateHashCode();
	}

	@SuppressWarnings("deprecation")
	public MarkerDetailSet(Iterable<Marker> markers) {
		this.markers = ImmutableList.copyOf(markers);
		this.markerSetFingerprint = MarkerSet.fingerprint(getMarkerNames());
		this.hashCode = generateHashCode();
	}

	private static BlastAnnotation closestChrMatch(byte naiveChr, int naivePosition,
																								 Iterable<BlastAnnotation> matches) {
		// If more than one match, try choosing one with chr that matches naiveChr (chr
		// won't be affected by build) and closest position to naivePosition
		BlastAnnotation bestMatch = null;
		for (BlastAnnotation annotation : matches) {
			if (annotation.getRefLoc().getChr() == naiveChr) {
				if (bestMatch == null) {
					bestMatch = annotation;
				} else {
					if (Math.abs(naivePosition
											 - annotation.getRefLoc().getStart()) < Math.abs(
																																			 naivePosition
																																			 - bestMatch.getRefLoc()
																																									.getStart())) {
						bestMatch = annotation;
					}
				}
			}
		}
		return bestMatch;
	}

	private int generateHashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((markers == null) ? 0 : markers.hashCode());
		return result;
	}

	private Map<String, Marker> generateMarkerNameMap() {
		ImmutableMap.Builder<String, Marker> markerNameMapBuilder = ImmutableMap.builder();
		for (Marker marker : markers) {
			markerNameMapBuilder.put(marker.getName(), marker);
		}
		return markerNameMapBuilder.build();
	}

	private ListMultimap<Byte, Marker> generateChrMap() {
		ImmutableListMultimap.Builder<Byte, Marker> chrMapBuilder = ImmutableListMultimap.builder();
		for (Marker marker : markers) {
			chrMapBuilder.put(marker.getChr(), marker);
		}
		chrMapBuilder.orderKeysBy(Ordering.natural());
		chrMapBuilder.orderValuesBy(Ordering.natural());
		return chrMapBuilder.build();
	}

	private SetMultimap<GenomicPosition, Marker> generateGenomicPositionMap() {
		ImmutableSetMultimap.Builder<GenomicPosition, Marker> genomicPositionMapBuilder = ImmutableSetMultimap.builder();
		for (Marker marker : markers) {
			genomicPositionMapBuilder.put(marker.getGenomicPosition(), marker);
		}
		genomicPositionMapBuilder.orderKeysBy(Ordering.natural());
		return genomicPositionMapBuilder.build();
	}

	private Map<Marker, Integer> generateMarkerIndexMap() {
		ImmutableMap.Builder<Marker, Integer> indexMapBuilder = ImmutableMap.builder();
		for (int i = 0; i < markers.size(); i++) {
			indexMapBuilder.put(markers.get(i), i);
		}
		return indexMapBuilder.build();
	}

	public static MarkerDetailSet parseFromBLASTAnnotation(MarkerSetInfo naiveMarkerSet,
																												 String blastAnnotation,
																												 Logger log) {
		String[] markerNames = naiveMarkerSet.getMarkerNames();
		List<Marker> markers = Lists.newArrayListWithCapacity(markerNames.length);
		if (!Files.exists(blastAnnotation)) {
			log.reportTimeWarning("Could not find " + blastAnnotation
														+ ", cannot generate BLAST Annotation based Marker Set");
			return null;
		}
		Map<String, MarkerBlastAnnotation> masterMarkerList = MarkerBlastAnnotation.initForMarkers(markerNames);
		MarkerAnnotationLoader annotationLoader = new MarkerAnnotationLoader(null, blastAnnotation,
																																				 new MarkerDetailSet(naiveMarkerSet),
																																				 true, log);
		annotationLoader.setReportEvery(10000);
		List<Map<String, ? extends AnnotationParser>> parsers = Lists.newArrayList();
		parsers.add(masterMarkerList);
		annotationLoader.fillAnnotations(null, parsers, QUERY_TYPE.ONE_TO_ONE);

		int missingABcount = 0;
		int missingPositionCount = 0;
		int ambiguousPositionCount = 0;
		for (int i = 0; i < markerNames.length; i++) {

			MarkerBlastAnnotation markerBlastAnnotation = masterMarkerList.get(markerNames[i]);
			int interrogationPosition = markerBlastAnnotation.getMarkerSeqAnnotation()
																											 .getInterrogationPosition();
			int probeLength = markerBlastAnnotation.getMarkerSeqAnnotation().getSequence().length();
			int positionOffset = interrogationPosition - probeLength + 1;
			char a = ABLookup.MISSING_ALLELE;
			char b = ABLookup.MISSING_ALLELE;
			Marker.RefAllele refAllele = null;
			try {
				MarkerSeqAnnotation markerSeqAnnotation = markerBlastAnnotation.getMarkerSeqAnnotation();
				char[] ab = ABLookup.parseABFromMarkerSeqAnnotation(markerSeqAnnotation);
				a = ab[0];
				b = ab[1];
				Allele ref = markerSeqAnnotation.getRef();
				if (ref.basesMatch(String.valueOf(a)))
					refAllele = Marker.RefAllele.A;
				else if (ref.basesMatch(String.valueOf(b)))
					refAllele = Marker.RefAllele.B;
			} catch (NullPointerException npe) {
				missingABcount++;
			}
			boolean ambiguousPosition = false;
			List<BlastAnnotation> perfectMatches = markerBlastAnnotation.getAnnotationsFor(BLAST_ANNOTATION_TYPES.PERFECT_MATCH,
																																										 log);
			BlastAnnotation bestMatch = null;
			if (perfectMatches.size() == 1) {
				bestMatch = perfectMatches.get(0);
			} else if (!perfectMatches.isEmpty()) {
				// TODO implement liftOver here to check which match is actually in naiveMarkerSet
				ambiguousPosition = true;
				bestMatch = closestChrMatch(naiveMarkerSet.getChrs()[i], naiveMarkerSet.getPositions()[i],
																		perfectMatches);
				if (bestMatch == null) {
					log.reportTimeWarning("Arbitrarily selecting first of " + perfectMatches.size()
																+ " perfect matches for " + markerNames[i]
																+ ", none of which are on the same chromosome as indicated by the array");
					bestMatch = perfectMatches.get(0);
				}
			} else {
				// Otherwise, choose lowest e-value match from non-perfect matches
				List<BlastAnnotation> onTargetMatches = markerBlastAnnotation.getAnnotationsFor(BLAST_ANNOTATION_TYPES.ON_T_ALIGNMENTS_NON_PERFECT,
																																												log);
				List<BlastAnnotation> offTargetMatches = markerBlastAnnotation.getAnnotationsFor(BLAST_ANNOTATION_TYPES.OFF_T_ALIGNMENTS,
																																												 log);
				Set<BlastAnnotation> nonPerfectMatches = Sets.newHashSet();
				nonPerfectMatches.addAll(onTargetMatches);
				nonPerfectMatches.addAll(offTargetMatches);
				List<BlastAnnotation> bestMatches = Lists.newArrayList();
				double bestMatchEVal = Double.MAX_VALUE;
				for (BlastAnnotation annotation : nonPerfectMatches) {
					double eVal = annotation.geteValue();
					if (eVal < bestMatchEVal) {
						bestMatches = Lists.newArrayList(annotation);
						bestMatchEVal = eVal;
					} else if (eVal == bestMatchEVal) {
						bestMatches.add(annotation);
					}
				}
				if (bestMatches.size() == 1) {
					bestMatch = bestMatches.get(0);
				} else if (!bestMatches.isEmpty()) {
					//
					ambiguousPosition = true;
					bestMatch = closestChrMatch(naiveMarkerSet.getChrs()[i], naiveMarkerSet.getPositions()[i],
																			bestMatches);
					if (bestMatch == null) {
						log.reportTimeWarning("Arbitrarily selecting first of " + bestMatches.size()
																	+ " best matches with identical E-values for " + markerNames[i]
																	+ ", none of which are on the same chromosome as indicated by the array");
						bestMatch = bestMatches.get(0);
					}
				}
			}
			byte chr;
			int position;
			if (bestMatch == null) {
				missingPositionCount++;
				chr = 0;
				position = 0;
			} else {
				Segment seg = bestMatch.getRefLoc();
				chr = seg.getChr();
				position = bestMatch.getEffectiveInterrogationPosition(positionOffset, log);
			}
			if (ambiguousPosition) {
				ambiguousPositionCount++;
			}
			markers.add(new Marker(markerNames[i], chr, position, a, b, refAllele));
		}
		if (missingABcount > 0) {
			log.reportError("Warning - there " + (missingABcount > 1 ? "were " : "was ")
											+ missingABcount + " marker" + (missingABcount > 1 ? "s" : "")
											+ " without an AB value");
		}
		if (ambiguousPositionCount > 0) {
			log.reportError("Warning - there " + (ambiguousPositionCount > 1 ? "were " : "was ")
											+ ambiguousPositionCount + " marker" + (ambiguousPositionCount > 1 ? "s" : "")
											+ " with ambiguous BLAST matches");
		}
		if (missingPositionCount > 0) {
			log.reportError("Warning - there " + (missingPositionCount > 1 ? "were " : "was ")
											+ missingPositionCount + " marker" + (missingPositionCount > 1 ? "s" : "")
											+ " missing a BLAST result and association position");
		}
		return new MarkerDetailSet(markers);


	}

	public static MarkerDetailSet load(String filename, boolean jar) {
		return (MarkerDetailSet) SerializedFiles.readSerial(filename, jar, false);
	}

	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}

	/**
	 * 
	 * @return a List of {@link Marker}s in defined order
	 */
	public List<Marker> getMarkers() {
		return markers;
	}


	public Map<String, Marker> getMarkerNameMap() {
		Map<String, Marker> markerNameMap = markerNameMapRef == null ? null : markerNameMapRef.get();
		if (markerNameMap == null) {
			markerNameMap = generateMarkerNameMap();
			markerNameMapRef = new SoftReference<Map<String, Marker>>(markerNameMap);
		}
		return markerNameMap;
	}

	/**
	 * 
	 * @return a {@link ListMultimap} from chromosome (sorted) to {@link Marker}s (sorted by position)
	 */
	public ListMultimap<Byte, Marker> getChrMap() {
		ListMultimap<Byte, Marker> chrMap = chrMapRef == null ? null : chrMapRef.get();
		if (chrMap == null) {
			chrMap = generateChrMap();
			chrMapRef = new SoftReference<ListMultimap<Byte, Marker>>(chrMap);
		}
		return chrMap;
	}

	/**
	 * 
	 * @return a {@link SetMultimap} from {@link GenomicPosition} to {@link Marker}s (sorted by
	 *         GenomicPosition)
	 */
	public SetMultimap<GenomicPosition, Marker> getGenomicPositionMap() {
		SetMultimap<GenomicPosition, Marker> genomicPositionMap = genomicPositionMapRef == null ? null
																																														: genomicPositionMapRef.get();
		if (genomicPositionMap == null) {
			genomicPositionMap = generateGenomicPositionMap();
			genomicPositionMapRef = new SoftReference<SetMultimap<GenomicPosition, Marker>>(genomicPositionMap);
		}
		return genomicPositionMap;
	}


	public Collection<Marker> getMarkersSortedByChrPos() {
		return getGenomicPositionMap().values();
	}

	public Map<Marker, Integer> getMarkerIndexMap() {
		Map<Marker, Integer> markerIndexMap = markerIndexMapRef == null ? null
																																		: markerIndexMapRef.get();
		if (markerIndexMap == null) {
			markerIndexMap = generateMarkerIndexMap();
			markerIndexMapRef = new SoftReference<Map<Marker, Integer>>(markerIndexMap);
		}
		return markerIndexMap;
	}

	@Override
	public void exportToText(Project proj, String filename) {
		PrintWriter writer;

		try {
			writer = Files.getAppropriateWriter(filename);
			for (Marker marker : markers) {
				writer.println(Joiner.on('\t').join(marker.getName(), marker.getChr(), marker.getPosition(),
																						marker.getA(), marker.getB()));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing " + filename);
			e.printStackTrace();
		}
	}

	public LinkedHashSet<Marker> getMarkersInSeg(Segment seg) {
		List<Marker> chrMarkers = getChrMap().get(seg.getChr());
		Iterator<Marker> markerIter = chrMarkers.iterator();
		Marker curMarker = null;
		while (markerIter.hasNext()) {
			curMarker = markerIter.next();
			if (curMarker.getPosition() >= seg.getStart())
				break;
		}
		LinkedHashSet<Marker> segMarkers = Sets.newLinkedHashSet();
		if (curMarker != null) {
			while (curMarker.getPosition() < seg.getStop()) {
				segMarkers.add(curMarker);
				if (markerIter.hasNext())
					curMarker = markerIter.next();
				else
					break;
			}
		}
		return segMarkers;
	}

	@Override
	public int hashCode() {
		return hashCode;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		MarkerDetailSet other = (MarkerDetailSet) obj;
		if (markers == null) {
			if (other.markers != null)
				return false;
		} else if (!markers.equals(other.markers))
			return false;
		return true;
	}

	@Override
	public long getFingerprint() {
		return markerSetFingerprint;
	}

	@SuppressWarnings("deprecation")
	@Override
	public boolean checkFingerprint(Sample samp) {
		return MarkerSet.checkFingerprints(this, samp);
	}

	@Override
	@Deprecated
	public int[] getIndicesOfMarkersIn(Segment seg, int[][] indicesByChr, Logger log) {
		Collection<Marker> segMarkers = getMarkersInSeg(seg);
		Iterator<Marker> segMarkersIter = segMarkers.iterator();
		Map<Marker, Integer> markerIndexMap = getMarkerIndexMap();
		int[] segMarkerIndices = new int[segMarkers.size()];
		for (int i = 0; i < segMarkerIndices.length; i++) {
			segMarkerIndices[i] = markerIndexMap.get(segMarkersIter.next());
		}
		return segMarkerIndices;
	}


	@Override
	@Deprecated
	public String[] getMarkersIn(Segment seg, int[][] indicesByChr) {
		Collection<Marker> segMarkers = getMarkersInSeg(seg);
		Iterator<Marker> segMarkersIter = segMarkers.iterator();
		String[] segMarkerNames = new String[segMarkers.size()];
		for (int i = 0; i < segMarkerNames.length; i++) {
			segMarkerNames[i] = segMarkersIter.next().getName();
		}
		return segMarkerNames;
	}

	@Override
	@Deprecated
	public String[] getMarkerNames() {
		String[] markerNames = new String[markers.size()];
		for (int i = 0; i < markers.size(); i++) {
			markerNames[i] = markers.get(i).getName();
		}
		return markerNames;
	}

	@Override
	@Deprecated
	public byte[] getChrs() {
		byte[] chrs = new byte[markers.size()];
		for (int i = 0; i < markers.size(); i++) {
			chrs[i] = markers.get(i).getChr();
		}
		return chrs;
	}

	@Override
	@Deprecated
	public int[] getPositions() {
		int[] positions = new int[markers.size()];
		for (int i = 0; i < markers.size(); i++) {
			positions[i] = markers.get(i).getPosition();
		}
		return positions;
	}

	@Override
	@Deprecated
	public int[][] getPositionsByChr() {
		ListMultimap<Byte, Marker> chrMap = getChrMap();
		int[][] positionsByChr = new int[MarkerSet.CHR_INDICES][];

		for (byte chr : chrMap.keySet()) {
			if (chr >= 0 && chr < MarkerSet.CHR_INDICES) {
				List<Marker> markers = chrMap.get(chr);
				positionsByChr[chr] = new int[markers.size()];
				for (int i = 0; i < markers.size(); i++) {
					positionsByChr[chr][i] = markers.get(i).getPosition();
				}
			}
		}
		return positionsByChr;
	}

	@Override
	@Deprecated
	public int[][] getIndicesByChr() {
		ListMultimap<Byte, Marker> chrMap = getChrMap();
		Map<Marker, Integer> markerIndexMap = getMarkerIndexMap();
		int[][] indicesByChr = new int[MarkerSet.CHR_INDICES][0];

		for (byte chr : chrMap.keySet()) {
			if (chr >= 0 && chr < MarkerSet.CHR_INDICES) {
				List<Marker> markers = chrMap.get(chr);
				indicesByChr[chr] = new int[markers.size()];
				for (int i = 0; i < markers.size(); i++) {
					indicesByChr[chr][i] = markerIndexMap.get(markers.get(i));
				}
			}
		}
		return indicesByChr;
	}

	@Deprecated
	public char[][] getABAlleles() {
		char[][] abAlleles = new char[markers.size()][2];
		for (int i = 0; i < markers.size(); i++) {
			abAlleles[i] = markers.get(i).getAB();
		}
		return abAlleles;
	}

	@Deprecated
	public Map<String, Integer> getMarkerIndices() {
		Map<Marker, Integer> markerIndices = getMarkerIndexMap();
		Map<String, Integer> indices = Maps.newHashMapWithExpectedSize(markerIndices.size());
		for (Map.Entry<Marker, Integer> entry : markerIndices.entrySet()) {
			indices.put(entry.getKey().toString(), entry.getValue());
		}
		return indices;
	}
}
