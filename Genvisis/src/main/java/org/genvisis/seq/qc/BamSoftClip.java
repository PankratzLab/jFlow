package org.genvisis.seq.qc;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.stats.Histogram.DynamicAveragingHistogram;

/**
 * Checkin relationship of soft clip count to insert size
 *
 */
public class BamSoftClip {

	private static class Summary {
		private Hashtable<String, Integer> seqCount;
		private HashSet<String> rgs;

		public Summary(Hashtable<String, Integer> seqCount, HashSet<String> rgs) {
			super();
			this.seqCount = seqCount;
			this.rgs = rgs;
		}

		public Hashtable<String, Integer> getSeqCount() {
			return seqCount;
		}
		

	}

	private static void clipVInsert(String[] bams, String outputDir, Logger log) {
		ArrayList<Summary> summaries = new ArrayList<Summary>();
		HashSet<String> allSofts = new HashSet<String>();
		for (int i = 0; i < bams.length; i++) {

			SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
			samReaderFactory.validationStringency(ValidationStringency.LENIENT);
			SamReader reader = samReaderFactory.open(new File(bams[i]));
			log.report(ext.getTime() + " Info - beginning processing of " + bams[i]);
			int numTotalReads = 0;
			DynamicAveragingHistogram dynamicAveragingHistogram = new DynamicAveragingHistogram(-700, 700, 0);
			Hashtable<String, Integer> seqCount = new Hashtable<String, Integer>();
			HashSet<String> rgs = new HashSet<String>();
			for (SAMRecord samRecord : reader) {
				numTotalReads++;
				if (numTotalReads % 3000000 == 0) {
					log.reportTimeInfo(numTotalReads + " total reads scanned");
					break;
				}
				if (!samRecord.getReadUnmappedFlag() && samRecord.getReadPairedFlag() && !samRecord.getMateUnmappedFlag() && !samRecord.getDuplicateReadFlag()) {
					int insert = samRecord.getInferredInsertSize();
					int numSoft = 0;
					rgs.add(samRecord.getReadGroup().getId());
					Cigar cigar = samRecord.getCigar();
					String[] bases = Array.decodeByteArray(samRecord.getReadBases(), log);
					int curStart = 0;
					int readIndex = 0;
					for (CigarElement cigarElement : cigar.getCigarElements()) {
						if (cigarElement.getOperator().consumesReadBases()) {
							readIndex += cigarElement.getLength();
						}
						if (cigarElement.getOperator() == CigarOperator.S) {
							numSoft += cigarElement.getLength();
							String softy = Array.toStr(Array.subArray(bases, curStart, readIndex), "");

							if (seqCount.containsKey(softy)) {
								seqCount.put(softy, seqCount.get(softy) + 1);
							} else {
								seqCount.put(softy, 1);
							}
						}
						if (cigarElement.getOperator().consumesReadBases()) {
							curStart += cigarElement.getLength();
						}
					}

					dynamicAveragingHistogram.addDataPair(insert, numSoft);
				}
			}
			dynamicAveragingHistogram.average();
			String out = outputDir + BamOps.getSampleName(bams[i]) + "clipSum.txt";
			ArrayList<String> outPrint = new ArrayList<String>();
			outPrint.add("InsertSize\tCount\tAvgSoftClipped");
			for (int j = 0; j < dynamicAveragingHistogram.getCounts().length; j++) {
				outPrint.add(dynamicAveragingHistogram.getBins()[j] + "\t" + dynamicAveragingHistogram.getCounts()[j] + "\t" + dynamicAveragingHistogram.getAverages()[j]);
			}
			String outCounts = outputDir + BamOps.getSampleName(bams[i]) + "clipSumCounts.txt";
			ArrayList<String> outPrintCount = new ArrayList<String>();
			outPrint.add("clippedBases\tCount");
			ArrayList<String> remove = new ArrayList<String>();
			for (String clip : seqCount.keySet()) {
				if (seqCount.get(clip) > 5 && clip.length() > 3) {
					allSofts.add(clip);
					outPrintCount.add(clip + "\t" + seqCount.get(clip));
				} else {
					remove.add(clip);
				}
			}
			for (String r : remove) {
				seqCount.remove(r);
			}
			Summary summary = new Summary(seqCount, rgs);
			summaries.add(summary);
			Files.writeList(Array.toStringArray(outPrintCount), outCounts);
			Files.writeList(Array.toStringArray(outPrint), out);
		}
		String outFinal = outputDir + "summaryCounts.txt";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outFinal));
			writer.println("Clipped\t" + Array.toStr(Array.tagOn(bams, "Count_", null)));
			for (String clip : allSofts) {
				writer.print(clip);
				for (Summary summary : summaries) {
					if (summary.getSeqCount().containsKey(clip)) {
						writer.print("\t" + summary.getSeqCount().get(clip));
					} else {
						writer.print("\t0");
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outFinal);
			log.reportException(e);
		}
	}

	public static void main(String[] args) {
		String bamFile = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/mutect/tnMatch.txt";
		String[] bams = HashVec.loadFileToStringArray(bamFile, false, new int[] { 1 }, true);
		String outputDir = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/mutect/clipMe/";
		new File(outputDir).mkdirs();
		Logger log = new Logger(outputDir + "clip.log");
		clipVInsert(bams, outputDir, log);
		// jcp seq.qc.BamSoftClip
	}

}
