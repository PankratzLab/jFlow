package seq.manage;

import common.Array;
import common.Files;
import common.Logger;
import common.Positions;
import common.WorkerHive;
import common.ext;
import filesys.Segment;
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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.concurrent.Callable;


final class BamExtractor {
	public static final String BP_BUFFER_COMMAND = "bpBuffer=";
	private Segment[] segmentsToExtract;
	private final String bamFile;
	private boolean verbose;
	private boolean fail;
	private boolean overWriteExisting;
	private int bpBuffer;
	private String[] bed;
	private final Logger log;

	public BamExtractor(Segment[] segmentsToExtract, String bamFile, int bpBuffer, boolean verbose, boolean overWriteExisting, Logger log) {
		this.segmentsToExtract = segmentsToExtract;
		this.bamFile = bamFile;
		this.bpBuffer = bpBuffer;
		this.verbose = verbose;
		this.overWriteExisting = overWriteExisting;
		this.bed = null;
		this.log = log;
		this.fail = (!verify());
	}

	private boolean verify() {
		if (!Files.exists(this.bamFile)) {
			this.log.reportError("Error - could not find file " + this.bamFile);
			return false;
		}
		return true;
	}

	public String[] getBed() {
		return this.bed;
	}

	public void setBed(String[] bed) {
		this.bed = bed;
	}

	public BamExtractor extractTo(String outputFile) {
		if ((!this.fail) && ((!Files.exists(outputFile)) || (this.overWriteExisting))) {
			SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
			samReaderFactory.validationStringency(ValidationStringency.LENIENT);
			SamReader reader = samReaderFactory.open(new File(this.bamFile));

			reader.indexing();
			if (!reader.hasIndex()) {
				this.log.reportError("Error - the bam file " + this.bamFile + " must have a \".bai\" index file associated with it, halting");
			} else {
				SAMFileHeader samFileHeader = reader.getFileHeader();
				QueryInterval[] qIntervals = convertSegsToQI(this.segmentsToExtract, samFileHeader, this.bpBuffer, this.log);
				qIntervals = QueryInterval.optimizeIntervals(qIntervals);
				if (qIntervals != null) {
					if (this.verbose) {
						this.log.reportTimeInfo("Attempting to extract " + qIntervals.length + " intervals with a bp buffer of " + this.bpBuffer);
						this.log.reportTimeInfo(this.bamFile + " " + ">" + " " + outputFile);
					}
					this.bed = getBedIntervals(qIntervals, samFileHeader);
					SAMFileWriter sAMFileWriter = new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reader.getFileHeader(), true, new File(outputFile));
					dumpIntervals(reader, sAMFileWriter, ext.removeDirectoryInfo(outputFile), qIntervals, this.log);
					sAMFileWriter.close();
					this.log.reportTimeInfo("Finished extracting" + qIntervals.length + " intervals with a bp buffer of " + this.bpBuffer);
					this.log.reportTimeInfo(this.bamFile + " " + ">" + " " + outputFile);
					try {
						reader.close();
					} catch (IOException e) {
						this.log.reportException(e);
						e.printStackTrace();
					}
				} else {
					this.fail = true;
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

	public static void extractAll(BamSample bamSample, String outputDirectory, int bpBuffer, boolean verbose, boolean overWriteExisting, int numThreads, Logger log) {
		WorkerHive<BamExtractor> hive = new WorkerHive<BamExtractor>(numThreads, 10, log);

		for (int i = 0; i < bamSample.getSamples().length; i++) {
			String curBam = bamSample.getBamForSampleAt(i);
			if (curBam != null) {
				String outputBam = outputDirectory + ext.removeDirectoryInfo(curBam);
				outputBam = ext.addToRoot(outputBam, ".mini");
				hive.addCallable(new WorkerExtractor(bamSample.getSegmentsToExtract(), curBam, verbose, overWriteExisting, bpBuffer, outputBam, log));
			}

		}
		hive.execute(true);
		Files.writeList(hive.getResults().get(0).getBed(), outputDirectory + "regions.bed");

	}

	private static class WorkerExtractor implements Callable<BamExtractor> {
		private BamExtractor bamExtractor;
		private String outputBam;

		public WorkerExtractor(Segment[] segmentsToExtract, String bamFile, boolean verbose, boolean overWriteExisting, int bpBuffer, String outputBam, Logger log) {
			this.bamExtractor = new BamExtractor(segmentsToExtract, bamFile, bpBuffer, verbose, overWriteExisting, log);
			this.outputBam = outputBam;
		}

		public BamExtractor call() {
			return this.bamExtractor.extractTo(this.outputBam);
		}
	}

	private void dumpIntervals(SamReader reader, SAMFileWriter sAMFileWriter, String rootFile, QueryInterval[] qIntervals, Logger log) {
		SAMRecordIterator sIterator = reader.query(qIntervals, false);
		int count = 0;
		while (sIterator.hasNext()) {
			SAMRecord samRecord = (SAMRecord) sIterator.next();
			if ((count != 0) && (count % 10000 == 0)) {
				log.reportTimeInfo("Wrote " + count + " reads (currently on " + samRecord.getReferenceName() + " for file " + rootFile);
			}
			count++;
			sAMFileWriter.addAlignment(samRecord);
		}
	}

	private static QueryInterval[] convertSegsToQI(Segment[] segs, SAMFileHeader sFileHeader, int bpBuffer, Logger log) {
		QueryInterval[] qIntervals = new QueryInterval[segs.length];
		segs = Segment.sortSegments(segs);
		for (int i = 0; i < qIntervals.length; i++) {
			String sequenceName = Positions.getChromosomeUCSC(segs[i].getChr(), true);
			int referenceIndex = sFileHeader.getSequenceIndex(sequenceName);
			if (referenceIndex < 0) {
				log.reportError("Error - could not find " + sequenceName + " in the sequence dictionary, halting");
				return null;
			}
			qIntervals[i] = new QueryInterval(referenceIndex, segs[i].getStart() - bpBuffer, segs[i].getStop() + bpBuffer);
		}
		return qIntervals;
	}

	public static void test(String bamFile) {
		Segment seg = new Segment((byte) 2, 0, 90000000);
		Segment[] segs = { seg };
		BamExtractor bamExtractor = new BamExtractor(segs, bamFile, 1000, true, true, new Logger());
		bamExtractor.extractTo(ext.addToRoot(bamFile, ".extracted"));
	}

	public static class BamSample {
		private String[] bamFiles;
		private Hashtable<String, String> bamSampleMap;
		private String[] samples;
		private ArrayList<Segment> segsToExtract;
		private Logger log;
		private boolean verbose;
		private boolean fail;

		public BamSample(String[] bamFiles, Logger log, boolean verbose) {
			this.bamFiles = bamFiles;
			this.bamSampleMap = new Hashtable<String, String>();
			this.samples = null;
			this.segsToExtract = new ArrayList<Segment>();
			this.log = log;
			this.verbose = verbose;
		}

		public Hashtable<String, String> getBamSampleMap() {
			return this.bamSampleMap;
		}

		public String[] getSamples() {
			return this.samples;
		}

		public void addSegmentToExtract(Segment seg) {
			this.segsToExtract.add(seg);
		}

		public Segment[] getSegmentsToExtract() {
			return (Segment[]) this.segsToExtract.toArray(new Segment[this.segsToExtract.size()]);
		}

		public String getBamForSampleAt(int index) {
			if (bamSampleMap.containsKey(samples[index])) {
				return (String) this.bamSampleMap.get(this.samples[index]);
			} else {
				return null;
			}
		}

		public boolean verify(String[] samples) {
			boolean verified = true;
			for (int i = 0; i < samples.length; i++) {
				if (!this.bamSampleMap.containsKey(samples[i])) {
					verified = false;
					this.log.reportTimeError("Could not find a matching bam file for sample " + samples[i]);
				}
			}
			return verified;
		}

		public boolean isFail() {
			return this.fail;
		}

		public void dumpToIGVMap(String correspondingVCF) {
			String output = correspondingVCF + ".mapping";
			String dumper = "";
			for (int i = 0; i < this.samples.length; i++) {
				String curSample = this.samples[i];
				if (this.bamSampleMap.containsKey(curSample)) {
					dumper = dumper + (i == 0 ? "" : "\n") + curSample + "\t" + "./"+ext.removeDirectoryInfo(this.bamSampleMap.get(curSample));
				} else {
					this.log.reportTimeError("did not find sample " + curSample + " in the sample map");
				}
			}
			Files.write(dumper, output);
		}

		public void generateMap() {
			ArrayList<String> tmpSamples = new ArrayList<String>();
			for (int i = 0; i < this.bamFiles.length; i++) {
				if (!this.fail) {
					if (Files.exists(this.bamFiles[i])) {
						SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
						samReaderFactory.validationStringency(ValidationStringency.LENIENT);
						SamReader reader = samReaderFactory.open(new File(this.bamFiles[i]));
						List<SAMReadGroupRecord> sGroupRecords = reader.getFileHeader().getReadGroups();
						for (int j = 0; j < sGroupRecords.size(); j++) {
							this.bamSampleMap.put(((SAMReadGroupRecord) sGroupRecords.get(j)).getSample(), this.bamFiles[i]);
							tmpSamples.add(((SAMReadGroupRecord) sGroupRecords.get(j)).getSample());
							//System.out.println((sGroupRecords.get(j)).getSample());
						}
						try {
							reader.close();
						} catch (IOException e) {
							this.log.reportException(e);
							this.fail = true;
						}
					} else {
						this.fail = true;
						this.log.reportError("Error - could not find bam file " + this.bamFiles[i]);
					}
				}
			}
			this.samples = Array.unique(Array.toStringArray(tmpSamples));
			if (this.samples.length != this.bamFiles.length) {
				this.log.reportTimeError("Detected more samples than bam files, currently this only supports single sample bam files");
				this.fail = true;
			}
			if (this.verbose) {
				this.log.reportTimeInfo("Found " + this.samples.length + " samples in " + this.bamFiles.length + " .bam files");
			}
		}
	}

	public static void main(String[] args) {
		//int numArgs = args.length;
		String filename = "D:/data/Project_Tsai_Project_021/testBamExtract/rrd_lane_HapMap_Control_CAGAGAGG-CTCTCTAT.merge.sorted.dedup.realigned.bam";
		test(filename);
	}
}
