package org.genvisis.seq.manage;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BamExtractor {
	public static final String BP_BUFFER_COMMAND = "bpBuffer=";
	private final Segment[] segmentsToExtract;
	private final String bamFile;
	private final boolean verbose;
	private boolean fail;
	private final boolean overWriteExisting;
	private final int bpBuffer;
	private String[] bed;
	private final Logger log;

	public BamExtractor(Segment[] segmentsToExtract, String bamFile, int bpBuffer, boolean verbose,
											boolean overWriteExisting, Logger log) {
		this.segmentsToExtract = segmentsToExtract;
		this.bamFile = bamFile;
		this.bpBuffer = bpBuffer;
		this.verbose = verbose;
		this.overWriteExisting = overWriteExisting;
		bed = null;
		this.log = log;
		fail = (!verify());
	}

	private boolean verify() {
		if (!Files.exists(bamFile)) {
			log.reportError("Error - could not find file " + bamFile);
			return false;
		}
		return true;
	}

	public String[] getBed() {
		return bed;
	}

	public void setBed(String[] bed) {
		this.bed = bed;
	}

	public BamExtractor extractTo(String outputFile) {
		if ((!fail) && ((!Files.exists(outputFile)) || (overWriteExisting))) {
			SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
			samReaderFactory.validationStringency(ValidationStringency.LENIENT);
			SamReader reader = samReaderFactory.open(new File(bamFile));

			reader.indexing();
			if (!reader.hasIndex()) {
				log.reportError("Error - the bam file "+ bamFile
												+ " must have a \".bai\" index file associated with it, halting");
			} else {
				SAMFileHeader samFileHeader = reader.getFileHeader();
				QueryInterval[] qIntervals = convertSegsToQI(	segmentsToExtract, samFileHeader, bpBuffer,
																											log);
				qIntervals = QueryInterval.optimizeIntervals(qIntervals);
				if (qIntervals != null) {
					if (verbose) {
						log.reportTimeInfo("Attempting to extract "+ qIntervals.length
																+ " intervals with a bp buffer of " + bpBuffer);
						log.reportTimeInfo(bamFile + " " + ">" + " " + outputFile);
					}
					bed = getBedIntervals(qIntervals, samFileHeader);
					SAMFileWriter sAMFileWriter =
																			new SAMFileWriterFactory().setCreateIndex(true)
																																.makeSAMOrBAMWriter(reader.getFileHeader(),
																																										true, new File(outputFile));
					dumpIntervals(reader, sAMFileWriter, ext.removeDirectoryInfo(outputFile), qIntervals,
												log);
					sAMFileWriter.close();
					log.reportTimeInfo("Finished extracting"+ qIntervals.length
															+ " intervals with a bp buffer of " + bpBuffer);
					log.reportTimeInfo(bamFile + " " + ">" + " " + outputFile);
					try {
						reader.close();
					} catch (IOException e) {
						log.reportException(e);
						e.printStackTrace();
					}
				} else {
					fail = true;
				}
			}
		}
		return this;
	}

	private static String[] getBedIntervals(QueryInterval[] qIntervals, SAMFileHeader samFileHeader) {
		String[] ucsc = new String[qIntervals.length];
		for (int i = 0; i < ucsc.length; i++) {
			String chr = samFileHeader.getSequence(qIntervals[i].referenceIndex).getSequenceName();
			int start = qIntervals[i].start;
			int stop = qIntervals[i].end;
			ucsc[i] = (chr + "\t" + start + "\t" + stop);
		}
		return ucsc;
	}

	public static void extractAll(BamSample bamSample, String outputDirectory, int bpBuffer,
																boolean verbose, boolean overWriteExisting, int numThreads,
																Logger log) {
		WorkerHive<BamExtractor> hive = new WorkerHive<BamExtractor>(numThreads, 10, log);

		for (int i = 0; i < bamSample.getSamples().length; i++) {
			String curBam = bamSample.getBamForSampleAt(i);
			if (curBam != null) {
				String outputBam = outputDirectory + ext.removeDirectoryInfo(curBam);
				outputBam = ext.addToRoot(outputBam, ".mini");
				hive.addCallable(new WorkerExtractor(	bamSample.getSegmentsToExtract(), curBam, verbose,
																							overWriteExisting, bpBuffer, outputBam, log));
			}

		}
		hive.execute(true);
		Files.writeArray(hive.getResults().get(0).getBed(), outputDirectory + "regions.bed");

	}

	public static class WorkerExtractor implements Callable<BamExtractor> {
		private final BamExtractor bamExtractor;
		private final String outputBam;

		public WorkerExtractor(	Segment[] segmentsToExtract, String bamFile, boolean verbose,
														boolean overWriteExisting, int bpBuffer, String outputBam, Logger log) {
			bamExtractor = new BamExtractor(segmentsToExtract, bamFile, bpBuffer, verbose,
																			overWriteExisting, log);
			this.outputBam = outputBam;
		}

		@Override
		public BamExtractor call() {
			return bamExtractor.extractTo(outputBam);
		}
	}

	private void dumpIntervals(	SamReader reader, SAMFileWriter sAMFileWriter, String rootFile,
															QueryInterval[] qIntervals, Logger log) {
		SAMRecordIterator sIterator = reader.query(qIntervals, false);
		int count = 0;
		while (sIterator.hasNext()) {
			SAMRecord samRecord = sIterator.next();
			if ((count != 0) && (count % 10000 == 0)) {
				log.reportTimeInfo("Wrote "+ count + " reads (currently on " + samRecord.getReferenceName()
														+ " for file " + rootFile);
			}
			count++;
			sAMFileWriter.addAlignment(samRecord);
		}
	}

	private static QueryInterval[] convertSegsToQI(	Segment[] segs, SAMFileHeader sFileHeader,
																									int bpBuffer, Logger log) {
		QueryInterval[] qIntervals = new QueryInterval[segs.length];
		segs = Segment.sortSegments(segs);
		for (int i = 0; i < qIntervals.length; i++) {
			String sequenceName = Positions.getChromosomeUCSC(segs[i].getChr(), true);
			int referenceIndex = sFileHeader.getSequenceIndex(sequenceName);
			if (referenceIndex < 0) {
				log.reportError("Error - could not find "+ sequenceName
												+ " in the sequence dictionary, halting");
				return null;
			}
			qIntervals[i] = new QueryInterval(referenceIndex, segs[i].getStart() - bpBuffer,
																				segs[i].getStop() + bpBuffer);
		}
		return qIntervals;
	}

	public static void test(String bamFile) {
		Segment seg = new Segment((byte) 2, 0, 90000000);
		Segment[] segs = {seg};
		BamExtractor bamExtractor = new BamExtractor(segs, bamFile, 1000, true, true, new Logger());
		bamExtractor.extractTo(ext.addToRoot(bamFile, ".extracted"));
	}

	public static class BamSample {
		private final String[] bamFiles;
		private final Hashtable<String, String> bamSampleMap;
		private String[] samples;
		private final ArrayList<Segment> segsToExtract;
		private final Logger log;
		private final boolean verbose;
		private boolean fail;

		public BamSample(String[] bamFiles, Logger log, boolean verbose) {
			this.bamFiles = bamFiles;
			bamSampleMap = new Hashtable<String, String>();
			samples = null;
			segsToExtract = new ArrayList<Segment>();
			this.log = log;
			this.verbose = verbose;
		}

		public Hashtable<String, String> getBamSampleMap() {
			return bamSampleMap;
		}

		public String[] getBamFiles() {
			return bamFiles;
		}

		public String[] getSamples() {
			return samples;
		}

		public void addSegmentToExtract(Segment seg) {
			segsToExtract.add(seg);
		}

		public Segment[] getSegmentsToExtract() {
			return segsToExtract.toArray(new Segment[segsToExtract.size()]);
		}

		public String getBamForSampleAt(int index) {
			if (bamSampleMap.containsKey(samples[index])) {
				return bamSampleMap.get(samples[index]);
			} else {
				return null;
			}
		}

		public boolean verify(String[] samples, String[] varset) {
			boolean verified = true;
			for (int i = 0; i < samples.length; i++) {

				if (varset == null && !bamSampleMap.containsKey(samples[i])) {
					verified = false;
					log.reportTimeError("Could not find a matching bam file for sample " + samples[i]);
				} else if (varset == null) {
					verified = true;

				} else {
					boolean hasOne = false;

					for (int j = 0; j < varset.length; j++) {
						if (samples[i].endsWith(varset[j])) {
							if (!bamSampleMap.containsKey(samples[i].replaceAll(varset[j], ""))) {
							} else {
								hasOne = true;
							}
						}
					}
					if (!hasOne) {
						log.reportTimeError("Could not find a matching bam file for sample "+ samples[i]
																+ " in the following var sets" + Array.toStr(varset, ","));

					}
					verified = hasOne;
				}
			}
			return verified;
		}

		public boolean isFail() {
			return fail;
		}

		public void dumpToIGVMap(String correspondingVCF, String[] varSets) {
			String output = correspondingVCF + ".mapping";

			String dumper = "";
			for (String curSample : samples) {
				if (bamSampleMap.containsKey(curSample)) {
					if (varSets != null) {
						for (String varSet : varSets) {
							dumper = dumper+ curSample + varSet + "\t" + "./"
												+ ext.removeDirectoryInfo(bamSampleMap.get(curSample) + "\n");
						}
					} else {
						dumper = dumper+ curSample + "\t" + "./"
											+ ext.removeDirectoryInfo(bamSampleMap.get(curSample) + "\n");
					}

				} else {
					log.reportTimeError("did not find sample " + curSample + " in the sample map");
				}
			}
			Files.write(dumper, output);
		}

		public void generateMap() {
			ArrayList<String> tmpSamples = new ArrayList<String>();
			for (String bamFile2 : bamFiles) {
				if (!fail) {
					if (Files.exists(bamFile2)) {
						SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
						samReaderFactory.validationStringency(ValidationStringency.LENIENT);
						SamReader reader = samReaderFactory.open(new File(bamFile2));
						List<SAMReadGroupRecord> sGroupRecords = reader.getFileHeader().getReadGroups();
						for (int j = 0; j < sGroupRecords.size(); j++) {
							bamSampleMap.put(sGroupRecords.get(j).getSample(), bamFile2);
							tmpSamples.add(sGroupRecords.get(j).getSample());
							// System.out.println((sGroupRecords.get(j)).getSample());
						}
						try {
							reader.close();
						} catch (IOException e) {
							log.reportException(e);
							fail = true;
						}
					} else {
						fail = true;
						log.reportError("Error - could not find bam file " + bamFile2);
					}
				}
			}
			samples = Array.unique(Array.toStringArray(tmpSamples));
			if (samples.length != bamFiles.length) {
				log.reportTimeError("Detected more samples than bam files, currently this only supports single sample bam files");
				fail = true;
			}
			if (verbose) {
				log.reportTimeInfo("Found "+ samples.length + " samples in " + bamFiles.length
														+ " .bam files");
			}
		}
	}

	public static void main(String[] args) {
		// int numArgs = args.length;
		String filename =
										"D:/data/Project_Tsai_Project_021/testBamExtract/rrd_lane_HapMap_Control_CAGAGAGG-CTCTCTAT.merge.sorted.dedup.realigned.bam";
		test(filename);
	}
}
