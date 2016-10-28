package org.genvisis.one.JL.mica;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

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

	private static final String STR = "GCT";
	private static final String START = "TGCTGTT";
	private static final String END = "ATTTTT";

	private static final String G_INSERT = "G";

	private static HashMap<String, Integer> getCountMap() {
		HashMap<String, Integer> countMap = new HashMap<String, Integer>();
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
				countMap.put(repeatRef.toString(), 0);
			}
			countMap.put(repeat.toString(), 0);
		}
		return countMap;
	}

	private static void run(String bamDir, String outDir, String mergeFile) {
		new File(outDir).mkdirs();
		HashMap<String, String> toMerge = loadMergeFile(mergeFile);
		Logger log = new Logger(outDir + "sk.log");
		String[] bams = Files.listFullPaths(bamDir, ".bam", false);
		log.reportTimeInfo("Found " + bams.length + " bams");
		Segment loc = new Segment("chr6:31,380,133-31,380,188");
		ArrayList<String> allRepeats = new ArrayList<>();
		allRepeats.addAll(getCountMap().keySet());
		StringBuilder results = new StringBuilder();
		StringBuilder resultsBasic = new StringBuilder();

		Collections.sort(allRepeats);
		results.append("SAMPLE");
		for (String repeat : allRepeats) {
			results.append("\t" + repeat);
			resultsBasic.append("\t" + repeat);

		}
		for (int i = 0; i < 3; i++) {
			results.append("\t" + toMerge.get("SAMPLE"));
		}
		results.append("\n");
		resultsBasic.append("\n");

		for (int i = 0; i < bams.length; i++) {
			HashMap<String, Integer> countMap = getCountMap();
			SamReader reader = BamOps.getDefaultReader(bams[i], ValidationStringency.STRICT);
			QueryInterval[] queryInterestIntervals = BamOps.convertSegsToQI(new Segment[] { loc },
					reader.getFileHeader(), 0, true, true, log);
			SAMRecordIterator sIterator = reader.query(queryInterestIntervals, false);
			log.reportTimeInfo("Analyzing " + bams[i]);
			while (sIterator.hasNext()) {
				SAMRecord samRecord = sIterator.next();
				if (samRecord.getAlignmentStart() < loc.getStart() && samRecord.getAlignmentEnd() > loc.getStop()) {
					String seq = samRecord.getReadString();
					for (String repeat : countMap.keySet()) {
						if (seq.toUpperCase().contains(repeat)) {
							countMap.put(repeat, countMap.get(repeat) + 1);
						}
					}
				}
			}
			sIterator.close();
			String samp = BamOps.getSampleName(bams[i], log);
			results.append(samp);
			for (String repeat : allRepeats) {
				results.append("\t" + countMap.get(repeat));
				resultsBasic.append("\t" + countMap.get(repeat));

			}
			results.append("\t" + toMerge.get(samp) + "\n");
			resultsBasic.append("\n");

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
