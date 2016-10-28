package org.genvisis.one.JL.mica;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.genvisis.CLI;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.analysis.PhaserNGS;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.StrandOps;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.annotation.Strand;

/**
 * Special k-mer counter for mica
 *
 */
public class SpecialK {

	private static final String STR = "GCT";
	private static final String START = "TT";
	private static final String END = "AT";

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
				countMap.put(StrandOps.flipsIfNeeded(repeatRef.toString(), Strand.NEGATIVE, false), 0);

			}
			countMap.put(repeat.toString(), 0);
			countMap.put(StrandOps.flipsIfNeeded(repeat.toString(), Strand.NEGATIVE, false), 0);

		}
		return countMap;
	}

	private static void run(String bamDir, String outDir) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "sk.log");
		String[] bams = Files.listFullPaths(bamDir, ".bam", false);
		log.reportTimeInfo("Found " + bams.length + " bams");
		Segment loc = new Segment("chr6:31,380,133-31,380,188");
		ArrayList<String> allRepeats = new ArrayList<>();
		allRepeats.addAll(getCountMap().keySet());
		StringBuilder results = new StringBuilder();
		Collections.sort(allRepeats);
		results.append("SampleName");
		for (String repeat : allRepeats) {
			results.append("\t" + repeat);
		}
		results.append("\n");

		for (int i = 0; i < bams.length; i++) {
			HashMap<String, Integer> countMap = getCountMap();
			SamReader reader = BamOps.getDefaultReader(bams[i], ValidationStringency.STRICT);
			QueryInterval[] queryInterestIntervals = BamOps.convertSegsToQI(new Segment[] { loc },
					reader.getFileHeader(), 0, true, true, log);
			SAMRecordIterator sIterator = reader.query(queryInterestIntervals, false);
			log.reportTimeInfo("Analyzing " + bams[i]);
			while (sIterator.hasNext()) {
				SAMRecord samRecord = sIterator.next();
				if (samRecord.getAlignmentStart() < loc.getStart() && samRecord.getAlignmentEnd() >loc.getStop()) {
					String seq = samRecord.getReadString();
					for (String repeat : countMap.keySet()) {
						if (seq.toUpperCase().contains(repeat)) {
							countMap.put(repeat, countMap.get(repeat) + 1);
						}
					}
				}
			}
			sIterator.close();
			results.append(BamOps.getSampleName(bams[i], log));
			for (String repeat : allRepeats) {
				results.append("\t" + countMap.get(repeat));
			}
			results.append("\n");

			try {
				reader.close();
			} catch (IOException e) {
				log.reportException(e);
			}

		}
		Files.write(results.toString(), outDir + "repeatCounts.txt");
	}

	public static void main(String[] args) {
		CLI c = new CLI(PhaserNGS.class);

		c.addArgWithDefault("bams", "directory of bams", null);
		c.addArgWithDefault("outDir", "output directory", null);
		c.parseWithExit(args);
		run(c.get("bams"), c.get("outDir"));
	}

}
