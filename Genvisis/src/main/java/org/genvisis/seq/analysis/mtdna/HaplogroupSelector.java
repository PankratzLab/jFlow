/**
 * 
 */
package org.genvisis.seq.analysis.mtdna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import org.genvisis.CLI;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.mtdna.HaploTrie.HaplogroupMatchResult;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;

/**
 * @author Kitty
 *
 */
public class HaplogroupSelector {
	/**
	 * 
	 */
	private static final String HAPLO_SELECT_RESULTS_X = "haploSelectResults_X";
	public static final String PHE_EXT = ".phe";
	public static final String KEEP_EXT = ".keep";

	private HaplogroupSelector() {

	}


	/**
	 * @param haplogrepFile
	 * @param selectFor select haplogroups for this group
	 * @param selectFrom select haplogroups from these groups
	 * @param vpopFile
	 * @param excludeFile remove from selection
	 * @param outDir
	 * @param xfactor number of unique samples to select per case
	 * @param minMatch minimum string matching for returned haplogroup
	 * @return
	 */
	public static String run(	String haplogrepFile, String selectFor, String[] selectFrom,
														String vpopFile, String excludeFile, String outDir, int xfactor,
														int minMatch) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "hapSelector.log");
		log.reportTimeInfo("Selecting haplogroup matches for "+ selectFor + " from "
												+ Array.toStr(selectFrom, ", and ") + " using ..."
												+ ext.removeDirectoryInfo(haplogrepFile));
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
		HashSet<String> samplesToChooseFrom = new HashSet<String>();
		HashSet<String> samplesToChooseFor = new HashSet<String>();

		for (String pop : selectFrom) {
			samplesToChooseFrom.addAll(vpop.getSuperPop().get(pop));
		}
		samplesToChooseFor.addAll(vpop.getSuperPop().get(selectFor));
		if (excludeFile != null) {
			HashSet<String> excludes = HashVec.loadFileToHashSet(excludeFile, false);
			samplesToChooseFrom.removeAll(excludes);
			samplesToChooseFor.removeAll(excludes);
		}
		log.reportTimeInfo("Selecting from "+ samplesToChooseFrom.size() + " controls for "
												+ samplesToChooseFor.size() + " cases");
		HaplogroupStruct haplogroupStruct = loadHaplogroups(haplogrepFile, samplesToChooseFor,
																												samplesToChooseFrom, log);
		HaploTrie haploTrie = populateTrie(haplogroupStruct);
		ArrayList<HaploMatch> haploMatchs = searchForMatches(	haplogroupStruct, haploTrie, xfactor,
																													minMatch);
		String[] header = HaploMatch.BASE_HEADER;
		for (int i = 0; i < xfactor; i++) {
			header = Array.concatAll(header, Array.tagOn(HaploMatch.BASE_HEADER, "Control_r1_", null));
		}
		StringBuilder builder = new StringBuilder(Array.toStr(header));
		StringBuilder caseControlBuilder = new StringBuilder();
		caseControlBuilder.append("##phe1,Integer,-9,\"1/2 = ARIC/CUSHING\"\n");
		caseControlBuilder.append("#ID\tphe1");
		ArrayList<String> samplesToKeep = new ArrayList<String>();
		ArrayList<String> unMatched = new ArrayList<String>();
		unMatched.add("UnMatchedSample\tHaplogroup");

		for (HaploMatch haploMatch : haploMatchs) {
			if (haploMatch.controlsMatched.size() == xfactor) {
				samplesToKeep.add(haploMatch.sample);
				caseControlBuilder.append("\n" + haploMatch.sample + "\t2");
				for (HaploMatch controlSelected : haploMatch.controlsMatched) {
					samplesToKeep.add(controlSelected.sample);
					caseControlBuilder.append("\n" + controlSelected.sample + "\t1");
				}
			} else {
				unMatched.add(haploMatch.sample+ "\t" + haploMatch.haplogroup + "\t"
											+ haploMatch.getResults());
			}
			builder.append("\n" + haploMatch.getResults());

		}
		String outFile = outDir + HAPLO_SELECT_RESULTS_X + xfactor + "_min" + minMatch + ".txt";
		String outFilePhe = outDir + HAPLO_SELECT_RESULTS_X + xfactor + "_min" + minMatch + PHE_EXT;
		String outFileKeeps = outDir + HAPLO_SELECT_RESULTS_X + xfactor + "_min" + minMatch + KEEP_EXT;
		String outFileUnmatched = outDir+ HAPLO_SELECT_RESULTS_X + xfactor + "_min" + minMatch
															+ ".unmatched";

		Files.write(builder.toString(), outFile);
		Files.write(caseControlBuilder.toString(), outFilePhe);
		Files.writeIterable(samplesToKeep, outFileKeeps);
		Files.writeIterable(unMatched, outFileUnmatched);

		return outFilePhe;
	}


	private static class HaploMatch implements Comparable<HaploMatch> {
		private static final String[] BASE_HEADER = new String[] {"Sample", "haplogroup"};

		private String sample;
		private String haplogroup;
		private int haploStringLength;
		private ArrayList<HaploMatch> controlsMatched;

		private HaploMatch(String sample, String haplogroup) {
			super();
			this.sample = sample;
			this.haplogroup = haplogroup;
			this.haploStringLength = haplogroup.length();
			this.controlsMatched = new ArrayList<HaplogroupSelector.HaploMatch>();
		}

		private String getResults() {
			StringBuilder builder = new StringBuilder();
			builder.append(sample);
			builder.append("\t" + haplogroup);
			for (HaploMatch coHaploMatch : controlsMatched) {
				builder.append("\t" + coHaploMatch.sample + "\t" + coHaploMatch.haplogroup);

			}
			return builder.toString();

		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(HaploMatch other) {
			return -1 * Integer.compare(this.haploStringLength, other.haploStringLength);
		}



	}

	private static ArrayList<HaploMatch> searchForMatches(HaplogroupStruct haplogroupStruct,
																												HaploTrie haploTrie, int xFactor,
																												int minMatch) {
		ArrayList<HaploMatch> haploMatches = new ArrayList<HaploMatch>();
		for (String haplogroup : haplogroupStruct.haploCase.keySet()) {
			for (String acase : haplogroupStruct.haploCase.get(haplogroup)) {
				haploMatches.add(new HaploMatch(acase, haplogroup));
			}
		}
		Collections.sort(haploMatches);// try to sort specific to general
		HashSet<String> taken = new HashSet<String>();
		for (HaploMatch haploMatch : haploMatches) {
			for (int i = 0; i < xFactor; i++) {

				HaplogroupMatchResult haplogroupMatchResult =
																										haploTrie.getBestMatchHaplotypeSamples(haploMatch.haplogroup);
				String sample = haplogroupMatchResult.getSample();
				if (taken.contains(sample)) {
					throw new IllegalStateException("Non-unique sample selected");
				}
				if (haploMatch.haplogroup.startsWith(haplogroupStruct.fullSampleMap	.get(sample)
																																						.substring(	0,
																																												minMatch))) {
					haploMatch.controlsMatched.add(new HaploMatch(sample,
																												haplogroupStruct.fullSampleMap.get(sample)));
					taken.add(sample);

				}
			}
		}
		return haploMatches;


	}


	private static HaploTrie populateTrie(

																				HaplogroupStruct haplogroupStruct) {
		HaploTrie haploTrie = new HaploTrie();
		for (String haplogroup : haplogroupStruct.haploControl.keySet()) {
			haploTrie.insert(haplogroup, haplogroupStruct.haploControl.get(haplogroup));
		}

		return haploTrie;

	}


	/**
	 * Stores relevant haplogroup info
	 *
	 */
	private static class HaplogroupStruct {
		private HashMap<String, HashSet<String>> haploControl;
		private HashMap<String, HashSet<String>> haploCase;
		private HashMap<String, String> fullSampleMap;

		private HaplogroupStruct(	HashMap<String, HashSet<String>> haploControl,
															HashMap<String, HashSet<String>> haploCase,
															HashMap<String, String> fullSampleMap) {
			super();
			this.haploControl = haploControl;
			this.haploCase = haploCase;
			this.fullSampleMap = fullSampleMap;
		}

	}

	private static HaplogroupStruct loadHaplogroups(String haplogrepFile, HashSet<String> choosefor,
																									HashSet<String> chooseFrom, Logger log) {
		HashMap<String, HashSet<String>> haploControl = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> haploCase = new HashMap<String, HashSet<String>>();
		HashMap<String, String> fullSampleMap = new HashMap<String, String>();

		int hapIndex = ext.indexOfStr("Haplogroup", Files.getHeaderOfFile(haplogrepFile, log));

		try {
			BufferedReader reader = Files.getAppropriateReader(haplogrepFile);

			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				String haplogroup = line[hapIndex];
				if (choosefor.contains(line[0])) {
					add(haploCase, line, haplogroup);// add sample
				} else if (chooseFrom.contains(line[0])) {
					add(haploControl, line, haplogroup);
				}
				fullSampleMap.put(line[0], haplogroup);

			}
			reader.close();
		} catch (FileNotFoundException e) {
			log.reportException(e);
		} catch (IOException e) {
			log.reportException(e);
		}


		return new HaplogroupStruct(haploControl, haploCase, fullSampleMap);
	}


	private static void add(HashMap<String, HashSet<String>> haploCase, String[] line,
													String haplogroup) {
		if (!haploCase.containsKey(haplogroup)) {
			haploCase.put(haplogroup, new HashSet<String>());
		}
		haploCase.get(haplogroup).add(line[0]);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		CLI c = new CLI(HaplogroupSelector.class);
		String haplogrepFile = "haplogrepResults.txt";
		String vpop = "a.vpop";
		String outDir = "haplogroups/";
		String selectFor = "Cases";
		String selectFrom = "Controls";
		String excludeFile = null;
		c.addArgWithDefault("haplogrepResults", "Haplogrep result file", haplogrepFile);
		c.addArgWithDefault("vpop", "VPOP file", vpop);
		c.addArgWithDefault("selectFor", "Select matches for ", selectFor);
		c.addArgWithDefault("selectFrom", "Select matches from (comma delimited) ", selectFrom);
		c.addArgWithDefault("excludeFile", "Samples in this file will not be used", excludeFile);

		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, outDir);

		c.parseWithExit(args);
		run(c.get("haplogrepResults"), c.get("selectFor"), c.get("selectFrom").split(","),
				c.get("vpop"), c.get("excludeFile"), c.get(CLI.ARG_OUTDIR), 5, 2);

	}

}
