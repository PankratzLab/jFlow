package org.genvisis.one.JL.mica;

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
import org.genvisis.common.Logger;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.analysis.PhaserNGS;
import org.genvisis.seq.manage.BamOps;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

/**
 * Special k-mer counter for mica
 *
 */
public class SpecialK {

	private static final String OTHER = "OTHER";
	private static final String TOTAL = "TOTAL";

	private static final String STR = "GCT";
	private static final String START = "TGCTGTT";
	private static final String END = "ATTTTT";
	private static final String G_INSERT = "G";

	private static final String UPSTREAM_INSERTION = "CAAGTCCCTTTTTTTTCAGG";
	private static final String DOWNSTREAM_VARIANT = "TTTTCTACGTCTGTTGTTGT";

	private static final String AND = "&&";

	private static class Kmers {
		private HashMap<String, Integer> countMap;
		private HashMap<String, String> key;

		public Kmers(HashMap<String, Integer> countMap, HashMap<String, String> key) {
			super();
			this.countMap = countMap;
			this.key = key;
		}

	}

	private static Kmers getCountMap() {
		HashMap<String, Integer> countMap = new HashMap<>();
		HashMap<String, String> key = new HashMap<>();

		for (int i = 4; i < 11; i++) {
			StringBuilder repeat = new StringBuilder();
			StringBuilder repeatRef = new StringBuilder();
			boolean addRef = false;
			repeat.append(START);
			repeatRef.append(START);

			for (int j = 0; j < i; j++) {
				repeat.append(STR);
				repeatRef.append(STR);

				if (i == 5 && j == 1) {
					repeatRef.append(G_INSERT);
					addRef = true;
				}
			}
			repeat.append(END);
			repeatRef.append(END);
			if (addRef) {
				key.put(repeatRef.toString(), "A5.1");
				key.put(repeatRef.toString() + AND + UPSTREAM_INSERTION, "A5.1_UPSTREAM_INSERTION");
				key.put(repeatRef.toString() + AND + DOWNSTREAM_VARIANT, "A5.1_DOWNSTREAM_VARIANT");
				key.put(repeatRef.toString() + AND + UPSTREAM_INSERTION + AND + DOWNSTREAM_VARIANT,
						"A5.1_UPSTREAM_INSERTION_DOWNSTREAM_VARIANT");

				countMap.put(repeatRef.toString(), 0);
				countMap.put(repeatRef.toString() + AND + UPSTREAM_INSERTION, 0);
				countMap.put(repeatRef.toString() + AND + DOWNSTREAM_VARIANT, 0);
				countMap.put(repeatRef.toString() + AND + UPSTREAM_INSERTION + AND + DOWNSTREAM_VARIANT, 0);
			}
			countMap.put(repeat.toString(), 0);
			countMap.put(repeat.toString() + AND + UPSTREAM_INSERTION, 0);
			countMap.put(repeat.toString() + AND + DOWNSTREAM_VARIANT, 0);
			countMap.put(repeat.toString() + AND + UPSTREAM_INSERTION + AND + DOWNSTREAM_VARIANT, 0);

			key.put(repeat.toString(), "A" + i);
			key.put(repeat.toString() + AND + UPSTREAM_INSERTION, "A" + i + "_UPSTREAM_INSERTION");
			key.put(repeat.toString() + AND + DOWNSTREAM_VARIANT, "A" + i + "_DOWNSTREAM_VARIANT");
			key.put(repeat.toString() + AND + UPSTREAM_INSERTION + AND + DOWNSTREAM_VARIANT,
					"A" + i + "_UPSTREAM_INSERTION_DOWNSTREAM_VARIANT");
		}

		countMap.put(TOTAL, 0);
		key.put(TOTAL, TOTAL);
		countMap.put(OTHER, 0);
		key.put(OTHER, OTHER);

		return new Kmers(countMap, key);
	}

	private static void run(String bamDir, String outDir, String mergeFile) {
		new File(outDir).mkdirs();
		HashMap<String, String> toMerge = loadMergeFile(mergeFile);
		Logger log = new Logger(outDir + "sk.log");
		String[] bams = Files.listFullPaths(bamDir, ".bam", false);
		log.reportTimeInfo("Found " + bams.length + " bams");
		Segment loc = new Segment("chr6:31,380,082-31,380,305");
		HashSet<String> upSDs = new HashSet<>();

		upSDs.add(DOWNSTREAM_VARIANT);
		upSDs.add(UPSTREAM_INSERTION);

		Kmers kmers = getCountMap();
		ArrayList<String> allRepeats = new ArrayList<>();
		allRepeats.addAll(kmers.countMap.keySet());
		StringBuilder results = new StringBuilder();
		StringBuilder resultsBasic = new StringBuilder();

		Collections.sort(allRepeats);
		results.append("SAMPLE");
		resultsBasic.append("SAMPLE");

		for (String repeat : allRepeats) {
			results.append("\t" + kmers.key.get(repeat) + "|" + repeat);
			resultsBasic.append("\t" + kmers.key.get(repeat) + "|" + repeat);
		}
		for (int i = 0; i < 3; i++) {
			results.append("\t" + toMerge.get("SAMPLE"));
		}
		results.append("\tOtherCounts\n");
		resultsBasic.append("\tOtherCounts\n");

		for (int i = 0; i < bams.length; i++) {
			// if (i > 10) {
			// break;
			// }
			int countOther = 0;
			HashMap<String, Integer> countMapAll = getCountMap().countMap;
			SamReader reader = BamOps.getDefaultReader(bams[i], ValidationStringency.STRICT);
			QueryInterval[] queryInterestIntervals = BamOps.convertSegsToQI(new Segment[] { loc },
					reader.getFileHeader(), 0, true, true, log);
			SAMRecordIterator sIterator = reader.query(queryInterestIntervals, false);
			log.reportTimeInfo("Analyzing " + bams[i]);
			while (sIterator.hasNext()) {
				SAMRecord samRecord = sIterator.next();
				if (samRecord.getAlignmentStart() < loc.getStart() && samRecord.getAlignmentEnd() > loc.getStop()) {
					String seq = samRecord.getReadString();
					boolean found = false;
					for (String repeats : countMapAll.keySet()) {
						String[] combos = repeats.split(AND);
						boolean foundAll = true;
						for (String repeat : combos) {
							if (!seq.toUpperCase().contains(repeat) && !repeat.equals(OTHER)) {
								foundAll = false;
							}
						}
						if (foundAll) {
							countMapAll.put(repeats, countMapAll.get(repeats) + 1);
							found = true;
						}
					}
					if (!found) {
						countMapAll.put(OTHER, countMapAll.get(OTHER) + 1);
					}
					countMapAll.put(TOTAL, countMapAll.get(TOTAL) + 1);
				}
			}
			sIterator.close();
			String samp = BamOps.getSampleName(bams[i], log);
			results.append(samp);
			resultsBasic.append(samp);

			for (String repeat : allRepeats) {
				results.append("\t" + countMapAll.get(repeat));
				resultsBasic.append("\t" + countMapAll.get(repeat));

			}
			results.append("\t" + toMerge.get(samp) + "\t" + countOther + "\n");
			resultsBasic.append("\t" + countOther + "\n");

			try {
				reader.close();
			} catch (IOException e) {
				log.reportException(e);
			}
		}
		Files.write(results.toString(), outDir + "repeatCounts.txt");
		Files.write(resultsBasic.toString(), outDir + "repeatBasicCounts.txt");
	}

	private static HashMap<String, String> loadMergeFile(String mergeFile) {
		HashMap<String, String> merge = new HashMap<String, String>();
		try {
			BufferedReader reader = Files.getAppropriateReader(mergeFile);
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				if (merge.containsKey(line[4])) {
					merge.put(line[4], merge.get(line[4]) + "\t" + Array.toStr(line));
				} else {
					merge.put(line[4], Array.toStr(line));

				}
			}
			reader.close();
		} catch (FileNotFoundException e) {

		} catch (IOException e) {

		}

		return merge;
	}

	public static void main(String[] args) {
		CLI c = new CLI(PhaserNGS.class);

		c.addArgWithDefault("bams", "directory of bams", null);
		c.addArgWithDefault("outDir", "output directory", null);
		c.addArgWithDefault("mergeFile", "mergeFile", null);

		c.parseWithExit(args);
		run(c.get("bams"), c.get("outDir"), c.get("mergeFile"));
	}

}
