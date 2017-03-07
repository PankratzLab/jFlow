package org.genvisis.cnv.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.cnv.annotation.markers.AnnotationFileLoader.QUERY_TYPE;
import org.genvisis.cnv.annotation.markers.AnnotationParser;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.markers.MarkerAnnotationLoader;
import org.genvisis.cnv.annotation.markers.MarkerBlastAnnotation;
import org.genvisis.cnv.manage.TextExport;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;

public class MarkerSet implements Serializable, TextExport {
	public static final long serialVersionUID = 1L;

	private final long fingerprint;
	private final String[] markerNames;
	private final byte[] chrs;
	private final int[] positions;
	private char[][] abAlleles;

	public MarkerSet(String[] markerNames, byte[] chrs, int[] positions) {
		this(markerNames, chrs, positions, null);
	}

	public MarkerSet(String[] markerNames, byte[] chrs, int[] positions, char[][] abAlleles) {
		int numMarkers = markerNames.length;
		if (numMarkers != chrs.length || numMarkers != positions.length
				|| (abAlleles != null && numMarkers != abAlleles.length)) {
			throw new IllegalArgumentException(this.getClass().getName()
																				 + " cannot be constructed with mismatched list of Markers and positions or AB Alleles");
		}
		this.markerNames = markerNames.clone();
		this.chrs = chrs.clone();
		this.positions = positions.clone();
		this.abAlleles = abAlleles == null ? null : abAlleles.clone();
		fingerprint = fingerprint(markerNames);
	}

	public String[] getMarkerNames() {
		return markerNames.clone();
	}

	public byte[] getChrs() {
		return chrs.clone();
	}

	public int[] getPositions() {
		return positions.clone();
	}

	public char[][] getABAlleles() {
		return abAlleles == null ? null : abAlleles.clone();
	}

	public Map<String, Integer> getMarkerIndices() {
		String[] markerNames = getMarkerNames();
		Map<String, Integer> indices = Maps.newHashMap();
		for (int i = 0; i < markerNames.length; i++) {
			indices.put(markerNames[i], i);
		}
		return indices;
	}

	public int[][] getPositionsByChr() {
		IntVector iv;
		byte chr;
		int[][] positionsByChr;
		boolean done;

		positionsByChr = new int[27][0];

		chr = 0;
		iv = new IntVector(20000);
		done = false;
		for (int i = 0; !done; i++) {
			if (i == chrs.length || chrs[i] != chr) {
				positionsByChr[chr] = Ints.toArray(iv);
				chr = i == chrs.length ? 0 : chrs[i];
				iv = new IntVector(20000);
			}
			if (i == chrs.length) {
				done = true;
			} else {
				iv.add(positions[i]);
			}
		}

		return positionsByChr;
	}

	public int[][] getIndicesByChr() {
		IntVector iv;
		byte chr;
		int[][] indicesByChr;
		boolean done;

		indicesByChr = new int[27][0];

		chr = 0;
		iv = new IntVector(20000);
		done = false;
		for (int i = 0; !done; i++) {
			if (i == chrs.length || chrs[i] != chr) {
				indicesByChr[chr] = Ints.toArray(iv);
				chr = i == chrs.length ? 0 : chrs[i];
				iv = new IntVector(20000);
			}
			if (i == chrs.length) {
				done = true;
			} else {
				iv.add(i);
			}
		}

		return indicesByChr;
	}

	public long getFingerprint() {
		return fingerprint;
	}

	@Override
	public void exportToText(Project proj, String filename) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (int i = 0; i < markerNames.length; i++) {
				writer.println(markerNames[i] + "\t" + chrs[i] + "\t" + positions[i]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing " + filename);
			e.printStackTrace();
		}
	}

	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}

	public static MarkerSet load(String filename, boolean jar) {
		return (MarkerSet) SerializedFiles.readSerial(filename, jar, true);
	}

	public static MarkerSet parseFromBLASTAnnotation(MarkerSet naiveMarkerSet, String blastAnnotation,
																									 Logger log) {
		String[] markerNames = naiveMarkerSet.getMarkerNames();
		byte[] chrs = new byte[markerNames.length];
		int[] positions = new int[markerNames.length];
		char[][] abLookup = new char[markerNames.length][];
		if (!Files.exists(blastAnnotation)) {
			log.reportTimeWarning("Could not find " + blastAnnotation
														+ ", cannot generate BLAST Annotation based Marker Set");
			return null;
		}
		Map<String, MarkerBlastAnnotation> masterMarkerList = MarkerBlastAnnotation.initForMarkers(markerNames);
		MarkerAnnotationLoader annotationLoader = new MarkerAnnotationLoader(null, blastAnnotation,
																																				 naiveMarkerSet, true, log);
		annotationLoader.setReportEvery(10000);
		List<Map<String, ? extends AnnotationParser>> parsers = Lists.newArrayList();
		parsers.add(masterMarkerList);
		annotationLoader.fillAnnotations(null, parsers, QUERY_TYPE.ONE_TO_ONE);

		int missingABcount = 0;
		int missingPositionCount = 0;
		int ambiguousPositionCount = 0;
		for (int i = 0; i < markerNames.length; i++) {

			MarkerBlastAnnotation markerBlastAnnotation = masterMarkerList.get(markerNames[i]);
			try {
				abLookup[i] = ABLookup.parseABFromMarkerSeqAnnotation(markerBlastAnnotation.getMarkerSeqAnnotation());
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
			if (bestMatch == null) {
				missingPositionCount++;
				chrs[i] = 0;
				positions[i] = 0;
			} else {
				Segment seg = bestMatch.getRefLoc();
				chrs[i] = seg.getChr();
				positions[i] = seg.getStart();
			}
			if (ambiguousPosition) {
				ambiguousPositionCount++;
			}
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
		return new MarkerSet(markerNames, chrs, positions, abLookup);


	}

	public static BlastAnnotation closestChrMatch(byte naiveChr, int naivePosition,
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

	public static long fingerprint(String[] names) {
		long sum, trav;

		sum = 0;
		for (int i = 0; i < names.length; i++) {
			trav = 0;
			for (int j = 0; j < names[i].length(); j++) {
				trav += names[i].charAt(j);
			}
			sum += trav * (i + 1);
		}

		return sum;
	}

	public int[] getIndicesOfMarkersIn(Segment seg, int[][] indicesByChr, Logger log) {
		return ext.indexLargeFactors(getMarkersIn(seg, indicesByChr), markerNames, true, log, true,
																 false);
	}


	public String[] getMarkersIn(Segment seg, int[][] indicesByChr) {
		int index = seg.getChr();
		ArrayList<String> markersIn = new ArrayList<String>();
		int[] indices = indicesByChr == null ? getIndicesByChr()[index] : indicesByChr[index];
		for (int indice : indices) {
			int bp = positions[indice];

			if (bp >= seg.getStart() && bp <= seg.getStop()) {
				markersIn.add(markerNames[indice]);
			}
			if (bp > seg.getStop()) {
				break;
			}
		}
		return ArrayUtils.toStringArray(markersIn);
	}

	public void checkFingerprint(Sample samp) {
		if (samp.getFingerprint() != fingerprint) {
			System.err.println("Error - Sample has a different fingerprint (" + samp.getFingerprint()
												 + ") than the MarkerSet (" + fingerprint + ")");
		}
	}

	public static void convert(String filename) {
		BufferedReader reader;
		String[] line;
		int count;
		String[] markerNames;
		byte[] chrs;
		int[] positions;

		try {
			reader = new BufferedReader(new FileReader(filename));
			count = 0;
			while (reader.ready()) {
				reader.readLine();
				count++;
			}
			reader.close();

			markerNames = new String[count];
			chrs = new byte[count];
			positions = new int[count];

			reader = new BufferedReader(new FileReader(filename));
			for (int i = 0; i < count; i++) {
				line = reader.readLine().trim().split("[\\s]+");
				markerNames[i] = line[0];
				chrs[i] = Positions.chromosomeNumber(line[1]);
				positions[i] = Integer.parseInt(line[2]);
			}
			reader.close();

			new MarkerSet(markerNames, chrs, positions, null).serialize(filename + ".ser");
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public static MarkerData[] loadFromList(Project proj, String[] markerNames) {
		// TODO remove completely
		return null;
	}

	public static byte[] translateABtoForwardGenotypes(byte[] abGenotypes, char[][] abLookup) {
		byte[] result;
		String geno;

		result = ArrayUtils.byteArray(abGenotypes.length, (byte) -3);
		for (int i = 0; i < abGenotypes.length; i++) {
			switch (abGenotypes[i]) {
				case 0:
					geno = abLookup[i][0] + "" + abLookup[i][0];
					break;
				case 1:
					geno = abLookup[i][0] + "" + abLookup[i][1];
					break;
				case 2:
					geno = abLookup[i][1] + "" + abLookup[i][1];
					break;
				case -1:
					geno = Sample.ALLELE_PAIRS[0];
					break;
				default:
					System.err.println("Error - invalid AB genotype: " + abGenotypes[i]);
					geno = null;
			}
			// System.out.println(geno);
			for (byte j = 0; j < Sample.ALLELE_PAIRS.length; j++) {
				if (geno.equals(Sample.ALLELE_PAIRS[j])) {
					result[i] = j;
				}
			}
		}

		return result;
	}

	/**
	 * Same functionality as a {@link MarkerSet} but indicesByChr are explicitly computed and stored
	 * since that can be expensive<br>
	 * Is nice if a function needs pos, chr, names, and indicesByChr across threads etc;
	 */
	public static class PreparedMarkerSet extends MarkerSet {

		/**
		 *
		 */
		private static final long serialVersionUID = 1L;

		private final int[][] indicesByChr;
		private final int[][] positionsByChr;

		public static PreparedMarkerSet getPreparedMarkerSet(MarkerSet markerSet) {
			if (markerSet == null) {
				return null;
			} else {
				return new PreparedMarkerSet(markerSet);
			}
		}

		private PreparedMarkerSet(MarkerSet markerSet) {
			super(markerSet.markerNames, markerSet.chrs, markerSet.positions, null);
			indicesByChr = super.getIndicesByChr();
			positionsByChr = super.getPositionsByChr();

		}

		@Override
		public int[][] getIndicesByChr() {
			return indicesByChr;
		}

		@Override
		public int[][] getPositionsByChr() {
			return positionsByChr;
		}

	}
}
