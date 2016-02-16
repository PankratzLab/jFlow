package seq.qc;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;

import scala.collection.generic.BitOperations.Int;
import seq.manage.BamOps;
import seq.manage.ReferenceGenome;
import stats.Histogram.DynamicAveragingHistogram;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

/**
 * Checkin relationship of soft clip count to insert size
 *
 */
public class BamSoftClip {

	private static void clipVInsert(String[] bams, String outputDir, Logger log) {
		for (int i = 0; i < bams.length; i++) {

			SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
			samReaderFactory.validationStringency(ValidationStringency.LENIENT);
			SamReader reader = samReaderFactory.open(new File(bams[i]));
			log.report(ext.getTime() + " Info - beginning processing of " + bams[i]);
			int numTotalReads = 0;
			DynamicAveragingHistogram dynamicAveragingHistogram = new DynamicAveragingHistogram(-700, 700, 0);
			Hashtable<String, Integer> seqCount = new Hashtable<String, Integer>();
			for (SAMRecord samRecord : reader) {
				numTotalReads++;
				if (numTotalReads % 1000000 == 0) {
					log.reportTimeInfo(numTotalReads + " total reads scanned");
					System.exit(1);
					break;
				}
				if (!samRecord.getReadUnmappedFlag() && samRecord.getReadPairedFlag() && !samRecord.getMateUnmappedFlag() && !samRecord.getDuplicateReadFlag()) {
					int insert = samRecord.getInferredInsertSize();
					int numSoft = 0;
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
							System.out.println(Array.toStr(bases));
							System.out.println(cigar);

							String softy = Array.toStr(Array.subArray(bases, curStart, readIndex), "");
							System.out.println(softy);
							try {
								Thread.sleep(100);
							} catch (InterruptedException ie) {
							}

							if (seqCount.containsKey(softy)) {
								seqCount.put(softy, seqCount.get(softy));
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
			Files.writeList(Array.toStringArray(outPrint), out);

		}
	}

	public static void main(String[] args) {
		String bamFile = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/mutect/tnMatch.txt";
		String[] bams = HashVec.loadFileToStringArray(bamFile, false, new int[] { 1 }, true);
		String outputDir = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/mutect/clipMe/";
		new File(outputDir).mkdirs();
		Logger log = new Logger(outputDir + "clip.log");
		clipVInsert(bams, outputDir, log);
	}

}
