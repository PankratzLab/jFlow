package org.genvisis.gwas.windows;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.util.Java6Helper;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.SciStringComparator;
import org.genvisis.common.ext;

import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;

/**
 * A utility class for computing hit windows: continuous regions of markers containing one or more
 * of those markers having a significant p-value. The p-value threshold for extending the window
 * (suggestive) is typically lower than the p-value for establishing a window (index).
 * <p>
 * The two primary uses of this class are {@link #findHitWindows(PFileStructure)} which computes the
 * windows, and {@link #generateHitWindows}, which writes the computed windows after computing them.
 * </p>
 */
public class MultiHitWindows {

	// -- command line args --

	public static final String PVAL_FILE = "pVals";
	public static final String MARKER_COLUMNS = "markerCols";
	public static final String CHR_COLUMNS = "chrCols";
	public static final String POS_COLUMNS = "posCols";
	public static final String PVAL_COLUMNS = "pvalCols";
	public static final String MARKER_SET = "markerSet";
	public static final String BUILD = "build";
	public static final String INDEX_PVAL = "iPVal";
	public static final String SUGGESTIVE_PVAL = "sPVal";
	public static final String WINDOW = "win";

	private MultiHitWindows() {
		// prevent instantiation of utility class
	}

	// -- Helper classes --

	// FIXME unify with HitWindow. These two are extremely close... I think this class probably should
	// just be a static utility method that generates the HitWindow. The blocker was tracking p-values
	// in MultiHit. Could add state to MultiHit but that would be sad.
	/**
	 * Statistics about a range of successful hits.
	 */
	private static class MultiHitWindow implements Comparable<MultiHitWindow> {
		/*
		 * Header for printing (matches #stats() output)
		 */
		private static final String[] HEADER = new String[] {"MARKER", "CHR", "START", "END",
																												 "NUM_INDEX",
																												 "NUM_SUGGESTIVE", "NUM_TOTAL", "P-VAL"};

		private final int chr;
		private final long start;
		private final long end;
		private final int numSuggestive;
		private final int numIndex;
		private final int numTotal;
		private final double lowestPval;
		private final List<MultiHit> indexHits;

		public MultiHitWindow(List<MultiHit> markers, String pValLabel, WindowThreshold hitParams) {
			Collections.sort(markers);
			chr = markers.get(0).getChr();
			start = markers.get(0).getPos();
			end = markers.get(markers.size() - 1).getPos();
			int suggestiveMarkers = 0;
			int indexMarkers = 0;
			numTotal = markers.size();
			List<MultiHit> index = new ArrayList<MultiHit>();

			double lowestP = Double.MAX_VALUE;
			for (MultiHit marker : markers) {
				double currentP = marker.pVal(pValLabel);
				// Check if the current marker has the lowest p-value in the window
				if (currentP < lowestP) {
					lowestP = currentP;
					index.clear();
				}
				// Maintain a list of all markers with the lowest p-value in this window
				if (Math.abs(currentP - lowestP) < 0.0000000001) {
					index.add(marker);
				}
				// Count the number of index and suggestive p-value hits
				if (currentP < hitParams.getIndexPval()) {
					indexMarkers++;
				} else if (currentP < hitParams.getSugPval()) {
					suggestiveMarkers++;
				}
			}

			this.numIndex = indexMarkers;
			this.numSuggestive = suggestiveMarkers;
			indexHits = Collections.unmodifiableList(index);
			lowestPval = lowestP;
		}

		/**
		 * @return A list of tab-delimited columns. Each list entry is one line to be printed (in the
		 *         case of multiple index SNPs)
		 */
		public List<String> stats() {
			List<String> stats = new ArrayList<String>();
			for (MultiHit index : indexHits) {
				StringBuilder sb = new StringBuilder();
				sb.append(index.getName()).append("\t");
				sb.append(chr).append("\t");
				sb.append(start).append("\t");
				sb.append(end).append("\t");
				sb.append(numIndex).append("\t");
				sb.append(numSuggestive).append("\t");
				sb.append(numTotal).append("\t");
				sb.append(lowestPval);
				stats.add(sb.toString());
			}
			return stats;
		}

		@Override
		public int compareTo(MultiHitWindow o) {
			int c = Java6Helper.compare(chr, o.chr);
			if (c == 0) {
				c = Java6Helper.compare(start, o.start);
			}
			if (c == 0) {
				c = Java6Helper.compare(end, o.end);
			}
			return c;
		}

	}

	// -- Utility methods --

	/**
	 * Skip centromere splitting
	 *
	 * @see #generateHitWindows(String, PFileStructure, WindowThreshold, int, String, Logger)
	 */
	public static void generateHitWindows(String outSuffix, PFileStructure pFile,
																				WindowThreshold hitParams) {
		generateHitWindows(outSuffix, pFile, hitParams, -1);
	}

	/**
	 * Use the given {@link Project} to find logging and centromere splitting info.
	 *
	 * @see #generateHitWindows(String, PFileStructure, WindowThreshold, int, String, Logger)
	 */
	public static void generateHitWindows(String outSuffix, PFileStructure pFile,
																				WindowThreshold hitParams, Project proj) {
		generateHitWindows(outSuffix, pFile, hitParams,
											 proj.GENOME_BUILD_VERSION.getValue().getBuildInt(),
											 proj.PLINK_DIR_FILEROOTS.getValueString() + "plink.bim",
											 proj.getLog());
	}

	/**
	 * Split centromeres using defaults for the specified genome build
	 *
	 * @see #generateHitWindows(String, PFileStructure, WindowThreshold, int, String, Logger)
	 */
	public static void generateHitWindows(String outSuffix, PFileStructure pFile,
																				WindowThreshold hitParams,
																				int genomeBuild) {
		generateHitWindows(outSuffix, pFile, hitParams, genomeBuild, null);
	}

	/**
	 * Split centromeres with the given markerset file.
	 *
	 * @see #generateHitWindows(String, PFileStructure, WindowThreshold, int, String, Logger)
	 */
	public static void generateHitWindows(String outSuffix, PFileStructure pFile,
																				WindowThreshold hitParams,
																				int genomeBuild, String markerSetFile) {
		generateHitWindows(outSuffix, pFile, hitParams, genomeBuild, markerSetFile, new Logger());
	}

	/**
	 * Find hit windows for the given {@link PFileStructure} and write them to disk.
	 *
	 * @see #findHitWindows(PFileStructure, WindowThreshold, int, String, Logger)
	 */
	public static void generateHitWindows(String outSuffix, PFileStructure pFile,
																				WindowThreshold hitParams,
																				int genomeBuild, String markerSetFile, Logger log) {

		Multimap<String, MultiHitWindow> hitWindows = findHitWindows(pFile, hitParams,
																																 genomeBuild, markerSetFile, log);

		if (hitWindows == null) {
			log.reportError("No hit windows detected. Aborting");
			return;
		}

		for (Entry<String, Collection<MultiHitWindow>> e : hitWindows.asMap().entrySet()) {
			log.report("Writing hit windows for p-value column: " + e.getKey());
			writeWindows(outSuffix, e.getKey(), e.getValue());
		}
	}

	/**
	 * Use default parameters when reading from the given p-value file.
	 *
	 * @see #findHitWindows(PFileStructure, WindowThreshold, int, String, Logger)
	 */
	public static Multimap<String, MultiHitWindow> findHitWindows(String pValFile,
																																String... pValColumns) {
		return findHitWindows(new PFileStructure(pValFile).pvals(pValColumns));
	}

	/**
	 * Use default {@link WindowThreshold}
	 *
	 * @see #findHitWindows(PFileStructure, WindowThreshold, int, String, Logger)
	 */
	public static Multimap<String, MultiHitWindow> findHitWindows(PFileStructure pFile) {
		return findHitWindows(pFile, new WindowThreshold());
	}

	/**
	 * Don't split centromeres
	 *
	 * @see #findHitWindows(PFileStructure, WindowThreshold, int, String, Logger)
	 */
	public static Multimap<String, MultiHitWindow> findHitWindows(PFileStructure pFile,
																																WindowThreshold hitParams) {
		return findHitWindows(pFile, hitParams, -1);
	}

	/**
	 * Use default centromere splitting for the given genome build
	 *
	 * @see #findHitWindows(PFileStructure, WindowThreshold, int, String, Logger)
	 */
	public static Multimap<String, MultiHitWindow> findHitWindows(PFileStructure pFile,
																																WindowThreshold hitParams,
																																int genomeBuild) {
		return findHitWindows(pFile, hitParams, genomeBuild, null);
	}

	/**
	 * Use default logger
	 *
	 * @see #findHitWindows(PFileStructure, WindowThreshold, int, String, Logger)
	 */
	public static Multimap<String, MultiHitWindow> findHitWindows(PFileStructure pFile,
																																WindowThreshold hitParams,
																																int genomeBuild,
																																String markerSetFile) {
		return findHitWindows(pFile, hitParams, genomeBuild, markerSetFile, new Logger());
	}

	/**
	 * @param pFile Descriptor of the p-value file to read
	 * @param hitParams Hit window detection parameters
	 * @param genomeBuild Build number for centromere splitting
	 * @param markerSetFile For centromere splitting
	 * @param log Where to log output
	 * @return Map of p-value columns to the list of hit windows for that p-value
	 */
	public static Multimap<String, MultiHitWindow> findHitWindows(PFileStructure pFile,
																																WindowThreshold hitParams,
																																int genomeBuild,
																																String markerSetFile,
																																Logger log) {

		// 1. parse MarkerPvals
		List<MultiHit> markers = parseMarkers(pFile, log);

		if (markers == null) {
			// No markers found
			return null;
		}

		// 2. Split MarkerPvals list to comparable components
		int[][] centromereBoundaries = Positions.determineCentromereBoundariesFromMarkerSet(markerSetFile,
																																												genomeBuild,
																																												log);
		Multimap<String, MultiHit> continuousMarkers = splitMarkers(markers, centromereBoundaries);

		// 3. Generate hit windows for each p-value column
		Multimap<String, MultiHitWindow> windows = TreeMultimap.<String, MultiHitWindow>create();

		for (int i = 0; i < pFile.getPvals().length; i++) {
			String pValLabel = pFile.getPvals()[i];
			List<MultiHitWindow> wins = findWindows(continuousMarkers, pValLabel, hitParams);
			windows.putAll(pValLabel, wins);
		}

		return windows;
	}

	/**
	 * @return List of all markers, according to the given {@link PFileStructure}
	 */
	private static List<MultiHit> parseMarkers(PFileStructure pFile, Logger log) {
		List<MultiHit> markers = new ArrayList<MultiHit>();

		try {
			BufferedReader reader = Files.getAppropriateReader(pFile.getFile());
			String temp = reader.readLine();
			String delimiter = ext.determineDelimiter(temp);
			String[] line = temp.trim().split(delimiter);

			int markerIndex = ext.indexFactors(new String[][] {pFile.getMarkers()}, line, false, false,
																				 true, true, log, false)[0];
			int chrIndex = ext.indexFactors(new String[][] {pFile.getChrs()}, line, false, false, true,
																			true, log, false)[0];
			int posIndex = ext.indexFactors(new String[][] {pFile.getPos()}, line, false, false, true,
																			true, log, false)[0];
			int[] pvalIndices = ext.indexFactors(pFile.getPvals(), line, false, log, false, false);

			if (markerIndex == -1 || chrIndex == -1 || posIndex == -1) {
				log.reportError("Aborting - did not find identifying column headers (marker, chr, position) in file: "
												+ pFile.getFile());
				return null;
			} else if (ArrayUtils.max(pvalIndices) == -1) {
				log.reportError("Aborting - no p-value columns found in file: " + pFile.getFile());
				return null;
			}

			// Report and prune any missing p-val columns
			boolean[] use = new boolean[pvalIndices.length];
			for (int i = 0; i < pvalIndices.length; i++) {
				if (pvalIndices[i] == -1) {
					log.reportTimeWarning("Column: " + pFile.getPvals()[i] + " not found in file: "
																+ pFile.getFile());
				} else {
					use[i] = true;
				}
			}
			pvalIndices = ArrayUtils.subArray(pvalIndices, use);
			pFile.pvals(ArrayUtils.subArray(pFile.getPvals(), use));

			String name;
			byte chr;
			int pos;

			// Parse the name, position and pvals for each marker
			while (reader.ready()) {
				line = reader.readLine().split(delimiter);
				name = line[markerIndex];
				chr = Byte.parseByte(line[chrIndex]);
				pos = Integer.parseInt(line[posIndex]);
				Map<String, Double> pValMap = new HashMap<String, Double>();
				for (int i = 0; i < pvalIndices.length; i++) {
					String label = pFile.getPvals()[i];
					String p = line[pvalIndices[i]];
					Double d = ext.isMissingValue(p) ? 999 : Double.parseDouble(p);
					pValMap.put(label, d);
				}

				markers.add(new MultiHit(name, chr, pos, pValMap));
			}

			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + pFile.getFile() + "\" not found in current directory");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + pFile.getFile() + "\"");
		}
		return markers;
	}

	/**
	 * Helper method to identify the the hit windows in each continuous set of markers, for the given
	 * p-value column
	 */
	private static List<MultiHitWindow> findWindows(Multimap<String, MultiHit> continuousMarkers,
																									String pValLabel, WindowThreshold hitParams) {
		List<MultiHitWindow> windows = new ArrayList<MultiHitWindow>();

		List<String> segments = new ArrayList<String>(continuousMarkers.keySet());
		Collections.sort(segments, new SciStringComparator());
		for (String segment : segments) {
			List<MultiHit> markers = new ArrayList<MultiHit>(continuousMarkers.get(segment));
			windows.addAll(findWindows(markers, pValLabel, hitParams));
		}

		return windows;
	}

	/**
	 * Helper method to find the hit windows in the given set of markers, for the given p-value column
	 */
	private static List<MultiHitWindow> findWindows(List<MultiHit> markers,
																									String pValLabel, WindowThreshold hitParams) {
		List<MultiHitWindow> windows = new ArrayList<MultiHitWindow>();

		for (int i = 0; i < markers.size(); i++) {
			if (markers.get(i).pVal(pValLabel) < hitParams.getIndexPval()) {

				int start = findIndex(markers, i, i - 1, -1, pValLabel, hitParams);
				int end = findIndex(markers, i, i + 1, 1, pValLabel, hitParams);
				windows.add(new MultiHitWindow(markers.subList(start, end + 1), pValLabel, hitParams));
				i = end;
			}
		}

		return windows;
	}

	/**
	 * Recursive helper to find the index of the last marker in the given window range
	 */
	private static int findIndex(List<MultiHit> markers, int lastGoodIndex, int currentIndex,
															 int searchStep, String pValLabel, WindowThreshold hitParams) {
		// boundary check
		if (currentIndex > 0 && currentIndex < markers.size()) {
			MultiHit currentMarker = markers.get(currentIndex);
			if (currentMarker.pVal(pValLabel) < hitParams.getSugPval()) {
				// current marker is good
				return findIndex(markers, currentIndex, currentIndex + searchStep, searchStep, pValLabel,
												 hitParams);
			} else if (Math.abs(currentMarker.getPos()
													- markers.get(lastGoodIndex).getPos()) <= hitParams.getWindow()) {
				// keep looking but don't update the last good index
				return findIndex(markers, lastGoodIndex, currentIndex + searchStep, searchStep, pValLabel,
												 hitParams);
			}
		}
		return lastGoodIndex;
	}

	/**
	 * Write the given hit windows
	 */
	private static void writeWindows(String outSuffix, String pvalName,
																	 Collection<MultiHitWindow> hitWindows) {
		String delim = "\t";
		PrintWriter writer = Files.getAppropriateWriter(pvalName + "_" + outSuffix);
		writer.println(ArrayUtils.toStr(MultiHitWindow.HEADER, delim));

		for (MultiHitWindow hw : hitWindows) {
			for (String line : hw.stats()) {
				writer.println(line);
			}
		}

		writer.flush();
		writer.close();
	}

	/**
	 * The purpose of this method is to bin a set of markers to facilitate future processing of them.
	 * Each key in the returned map will map to a list of markers on a) the same chromosome and b) the
	 * same side of the centromere.
	 *
	 * @param markers Genome-wide list of {@link MultiHit}
	 * @param centromereBoundaries (optional) list collection of centromere start and end points,
	 *        indexed by chromosome
	 * @return A map of continuous chromosomal segments to the {@link MultiHit} in those segments.
	 *         Each individual set will be sorted.
	 */
	private static Multimap<String, MultiHit> splitMarkers(List<MultiHit> markers,
																												 int[][] centromereBoundaries) {
		Multimap<String, MultiHit> continuousMarkers = TreeMultimap.<String, MultiHit>create();
		// TODO: an improvement would be to read over all problematic regions and break groups at each
		// point

		for (MultiHit mp : markers) {
			int chr = mp.getChr();
			if (centromereBoundaries == null) {
				continuousMarkers.put(Integer.toString(chr), mp);
			} else if (mp.getPos() < centromereBoundaries[chr][0]) {
				continuousMarkers.put(chr + "a", mp);
			} else if (mp.getPos() > centromereBoundaries[chr][1]) {
				continuousMarkers.put(chr + "b", mp);
			}
		}

		return continuousMarkers;
	}

	/**
	 * Command-line entry point
	 */
	public static void main(String... args) {
		CLI c = new CLI(MultiHitWindows.class);

		c.addArg(PVAL_FILE,
						 "Input file containing one or more p-values for a set of marker names and positions",
						 "plink_results_pVals_Additive.xln", true);
		c.addArgWithDefault(CLI.ARG_OUTFILE, "Suffix for all output files", "hits.out");
		c.addArgWithDefault(INDEX_PVAL, "p-value threshold for defining a window",
												Double.toString(WindowThreshold.DEFAULT_INDEX_PVAL),
												CLI.Arg.NUMBER);
		c.addArgWithDefault(SUGGESTIVE_PVAL, "p-value threshold for extending a window",
												Double.toString(WindowThreshold.DEFAULT_SUGGESTIVE_PVAL),
												CLI.Arg.NUMBER);
		c.addArgWithDefault(WINDOW, "window search space (in bp)",
												Integer.toString(WindowThreshold.DEFAULT_WINDOW), CLI.Arg.NUMBER);
		c.addArgWithDefault(MARKER_COLUMNS,
												"Comma-separated list of column aliases containing marker names",
												ArrayUtils.toStr(Aliases.MARKER_NAMES, ","));
		c.addArgWithDefault(CHR_COLUMNS,
												"Comma-separated list column aliases containing marker positions (on chromosome)",
												ArrayUtils.toStr(Aliases.CHRS, ","));
		c.addArgWithDefault(POS_COLUMNS,
												"Comma-separated list of column aliases containing marker positions (in base pairs)",
												ArrayUtils.toStr(Aliases.POSITIONS, ","));
		c.addArg(PVAL_COLUMNS,
						 "Comma-separated list of columns containing p-values. NB: hits will be computed for ALL matching p-value columns.",
						 "tdt_P,emim_P", true);
		c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ + " - used for centromere determination", false);
		c.addArg(MARKER_SET, "Path to marker set file - for centromere determination if no project",
						 false);
		c.addArg(BUILD, "GRC genome build number - for centromere determination if no project", false,
						 CLI.Arg.NUMBER);

		c.parseWithExit(args);

		final String outSuffix = c.get(CLI.ARG_OUTFILE);
		PFileStructure pFile = new PFileStructure(c.get(PVAL_FILE));
		pFile.markers(c.get(MARKER_COLUMNS).split(","))
				 .chrs(c.get(CHR_COLUMNS).split(","))
				 .pos(c.get(POS_COLUMNS).split(","))
				 .pvals(c.get(PVAL_COLUMNS).split(","));

		WindowThreshold hitParams = new WindowThreshold().index(c.getD(INDEX_PVAL))
																										 .suggestive(c.getD(SUGGESTIVE_PVAL))
																										 .window(c.getI(WINDOW));

		if (c.has(CLI.ARG_PROJ)) {
			generateHitWindows(outSuffix, pFile, hitParams, new Project(c.get(CLI.ARG_PROJ), false));
		} else if (c.has(MARKER_SET)) {
			generateHitWindows(outSuffix, pFile, hitParams, -1, c.get(MARKER_SET));
		} else if (c.has(BUILD)) {
			generateHitWindows(outSuffix, pFile, hitParams, c.getI(BUILD));
		} else {
			generateHitWindows(outSuffix, pFile, hitParams);
		}
	}
}
