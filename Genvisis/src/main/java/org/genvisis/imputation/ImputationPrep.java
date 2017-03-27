package org.genvisis.imputation;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.poi.util.IOUtils;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker.GenomicPosition;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ProgressMonitor;
import org.genvisis.common.ext;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class ImputationPrep {

	private static final Set<Character> VALID_ALLELES = ImmutableSet.of('A', 'T', 'G', 'C', 'I', 'D');
	private static final String[][] REF_COLS;

	static {
		String[] markerNameAliases = Arrays.copyOf(Aliases.MARKER_NAMES,
																							 Aliases.MARKER_NAMES.length + 1);
		markerNameAliases[markerNameAliases.length - 1] = "ID";

		REF_COLS = new String[][] {Aliases.CHRS, Aliases.POSITIONS, markerNameAliases,
															 Aliases.REF_ALLELES, Aliases.ALT_ALLELES};
	}
	private static final Map<Character, Character> PALINDROMIC_PAIRS = ImmutableMap.of('A', 'T',
																																										 'T', 'A',
																																										 'G', 'C',
																																										 'C', 'G');

	private final Project proj;
	private final String referenceFile;
	private final Logger log;

	private final Set<Marker> matchingMarkers;

	private static class ReferencePosition {
		private byte chr;
		private int position;
		private String id;
		private char ref;
		private char alt;

		/**
		 * @param chr
		 * @param position
		 * @param id
		 * @param ref
		 * @param alt
		 */
		public ReferencePosition(byte chr, int position, String id, char ref, char alt) {
			super();
			this.chr = chr;
			this.position = position;
			this.id = id;
			this.ref = ref;
			this.alt = alt;
		}

		public byte getChr() {
			return chr;
		}

		public int getPosition() {
			return position;
		}

		public String getId() {
			return id;
		}

		public char getRef() {
			return ref;
		}

		public char getAlt() {
			return alt;
		}


	}


	/**
	 * @param proj Project to run for. Must have a {@link Project#BLAST_ANNOTATION_FILENAME} or an
	 *        {@link MarkerDetailSet} that otherwise includes Ref/Alt alleles
	 * @param referenceFile The reference set file (preferably gzipped). Must include chromosome,
	 *        position, ref allele, and alt allele
	 */
	public ImputationPrep(Project proj, String referenceFile) {
		super();
		this.proj = proj;
		this.referenceFile = referenceFile;
		this.log = proj.getLog();
		matchingMarkers = filterMarkers();
	}

	public Set<Marker> getMatchingMarkers() {
		return matchingMarkers;
	}

	@SuppressWarnings("resource") // eclipse doesn't recognize reader.close() calls in catch blocks, apparently
	private Map<Byte, Map<Integer, Set<ReferencePosition>>> readRefFile() {
		Set<GenomicPosition> markerSetPositions = proj.getMarkerSet().getGenomicPositionMap().keySet();
		log.report("Parsing reference panel file");
		long time = System.currentTimeMillis();
		Map<Byte, Map<Integer, Set<ReferencePosition>>> referencePositionsBuild = Maps.newHashMap();
		String taskName = "parseRefFile";
		proj.getProgressMonitor().beginDeterminateTask(taskName, "Parsing reference panel file",
																									 Files.countLines(referenceFile, 0),
																									 ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		BufferedReader reader = null;
		try {
			reader = Files.getAppropriateReader(referenceFile);
			String header = reader.readLine();
			if (header == null) {
				log.reportError("Reference file is empty");
				throw new IllegalArgumentException("Reference file is empty");
			}
			String delim = ext.determineDelimiter(header);
			int[] cols = ext.indexFactors(REF_COLS, header.split(delim), false, true, false, false);
			for (int index : cols) {
				if (index < 0) {
					throw new IllegalArgumentException("Invalid Reference File header");
				}
			}
			String temp;
			while ((temp = reader.readLine()) != null) {
				String[] refLine = ArrayUtils.subArray(temp.split(delim), cols);
				proj.getProgressMonitor().updateTask(taskName);
				byte chr = Positions.chromosomeNumber(refLine[0], log);
				int position;
				try {
					position = Integer.parseInt(refLine[1]);
				} catch (NumberFormatException e) {
					throw new IllegalStateException("Imputation reference file (" + referenceFile
																					+ ") contains a non-integer position: " + refLine[1]);
				}
				if (markerSetPositions.contains(new GenomicPosition(chr, position))) {
					String id = refLine[2];
					char ref = refLine[3].length() == 1 ? refLine[3].toUpperCase().charAt(0)
																							: ABLookup.MISSING_ALLELE;
					char alt = refLine[4].length() == 1 ? refLine[4].toUpperCase().charAt(0)
																							: ABLookup.MISSING_ALLELE;
					Map<Integer, Set<ReferencePosition>> posMap = referencePositionsBuild.get(chr);
					if (posMap == null) {
						posMap = Maps.newHashMap();
						referencePositionsBuild.put(chr, posMap);
					}
					Set<ReferencePosition> refPosSet = posMap.get(position);
					if (refPosSet == null) {
						refPosSet = Sets.newHashSet();
						posMap.put(position, refPosSet);
					}
					refPosSet.add(new ReferencePosition(chr, position, id, ref, alt));
				}
			}
		} catch (IOException ioe) {
			log.reportIOException(referenceFile);
		} finally {
			IOUtils.closeQuietly(reader);
		}
		proj.getProgressMonitor().endTask(taskName);
		log.report("Finished parsing Reference File in " + ext.getTimeElapsed(time));
		return referencePositionsBuild;
	}

	private Set<Marker> filterMarkers() {
		MarkerDetailSet markerSet = proj.getMarkerSet();
		List<Marker> markers = markerSet.getMarkers();
		Map<Byte, Map<Integer, Set<ReferencePosition>>> referencePositions = readRefFile();
		ImmutableSet.Builder<Marker> matchingMarkersBuilder = ImmutableSet.builder();
		int mismatchPos = 0;
		int invalidAlleles = 0;
		int alleleFlips = 0;
		int strandFlips = 0;
		int strandAlelleFlips = 0;
		int mismatchAlleles = 0;

		String taskName = "Matching Project Markers to Reference Panel";
		proj.getProgressMonitor().beginDeterminateTask(taskName, taskName,
																									 markers.size(),
																									 ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		List<Marker> noMatches = Lists.newArrayList();
		for (Marker marker : markers) {
			proj.getProgressMonitor().updateTask(taskName);
			Map<Integer, Set<ReferencePosition>> chrMap = referencePositions.get(marker.getChr());
			if (chrMap == null) {
				log.reportError("Warning - Chr " + marker.getChr() + " missing from reference file");
				referencePositions.put(marker.getChr(),
															 new HashMap<Integer, Set<ReferencePosition>>());
				chrMap = referencePositions.get(marker.getChr());
			}
			Set<ReferencePosition> refMatches = chrMap.get(marker.getPosition());
			if (refMatches == null) {
				mismatchPos++;
				noMatches.add(marker);
			} else {
				// Copy the matches in order to pick off bad matches without removing from actual structure
				refMatches = Sets.newHashSet(refMatches);
				char ref = marker.getRef();
				char alt = marker.getAlt();
				// Find the best matched ReferencePosition from all of the matches
				if (!validateAlleles(ref, alt, refMatches)) {
					invalidAlleles++;
				} else if (matches(ref, alt, refMatches)) {
					matchingMarkersBuilder.add(marker);
					// Don't accept anything less than a proper match but count for reporting purposes
				} else if (matches(alt, ref, refMatches)) {
					alleleFlips++;
				} else {
					Character aFlip = PALINDROMIC_PAIRS.get(ref);
					Character bFlip = PALINDROMIC_PAIRS.get(alt);
					if (matches(aFlip, bFlip, refMatches)) {
						strandFlips++;
					} else if (matches(bFlip, aFlip, refMatches)) {
						strandAlelleFlips++;
					} else {
						mismatchAlleles++;
					}
				}
			}
		}

		proj.getProgressMonitor().endTask(taskName);

		Set<Marker> matchingMarkersSet = matchingMarkersBuilder.build();

		log.report(invalidAlleles
							 + " positions had invalid allele codes in the project or reference set");
		log.report(mismatchPos + " positions did not match to a position in the reference set");
		log.report(mismatchAlleles + " positions had mismatched alleles and were excluded");
		log.report(alleleFlips + " positions would have matched after flipping the alleles");
		log.report(strandFlips + " positions would have matched after flipping the strand");
		log.report(strandAlelleFlips
							 + " positions would have matched after flipping the alleles and strand");
		log.report("A total of " + (markers.size() - matchingMarkersSet.size())
							 + " positions were excluded");
		log.report("");
		log.report(matchingMarkersSet.size()
							 + " positions matched the reference set and will be used for imputation");
		return matchingMarkersSet;
	}

	private static boolean validateAlleles(char a, char b, Collection<ReferencePosition> refMatches) {
		if (!validAlleles(a, b)) {
			return false;
		}
		for (Iterator<ReferencePosition> i = refMatches.iterator(); i.hasNext();) {
			ReferencePosition refPos = i.next();
			char ref = refPos.getRef();
			char alt = refPos.getAlt();
			if (!validAlleles(ref, alt)) {
				i.remove();
			}
		}
		return !refMatches.isEmpty();
	}

	private static boolean validAlleles(char... alleles) {
		for (char allele : alleles) {
			if (!VALID_ALLELES.contains(allele)) {
				return false;
			}
		}
		return true;
	}

	private static boolean matches(Character a, Character b,
																 Collection<ReferencePosition> refMatches) {
		if (a != null && b != null) {
			for (ReferencePosition refPos : refMatches) {
				if (a == refPos.getRef() && b == refPos.getAlt()) {
					return true;
				}
			}
		}
		return false;
	}
}
