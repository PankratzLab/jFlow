package seq.qc;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.concurrent.Callable;

import common.Array;
import common.Files;
import common.Logger;
import common.PSF;
import common.Positions;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;

/**
 * Class to check determine the sex of an individual based on read counts from X and Y chromosomes.
 *
 */
public class SexCheck {
	private static final String[] SEX_CHECK_HEADER = { "DNA", "Bam_File", "CountX", "CountY", "PropX", "PropY" };
	private static final byte X = 23;
	private static final byte Y = 24;
	private static final int MAPQ_FILTER = 60;

	private Logger log;
	private WorkerTrain<SexCheckResults> train;
	private SexCheckProducer producer;

	/**
	 * @param bamFiles
	 *            files to check
	 * @param numThreads
	 *            number of threads for da checking
	 * @param log
	 */
	public SexCheck(String[] bamFiles, int numThreads, Logger log) {
		super();
		this.producer = new SexCheckProducer(bamFiles, log);
		this.train = new WorkerTrain<SexCheckResults>(producer, numThreads, 1, log);
		this.log = log;
	}

	public void checkSex(String fullPathTooutput) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fullPathTooutput));
			writer.println(Array.toStr(SEX_CHECK_HEADER));
			while (train.hasNext()) {
				SexCheckResults sexCheckResults = train.next();
				writer.println(sexCheckResults.getSample() + "\t" + sexCheckResults.getBamFile() + "\t" + sexCheckResults.getNumXReads() + "\t" + sexCheckResults.getNumYReads() + "\t" + sexCheckResults.getPropX() + "\t" + sexCheckResults.getPropY());
				writer.flush();
			}

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + fullPathTooutput);
			log.reportException(e);
		}
	}

	/**
	 * Dishes up callables of {@link SexCheckResults}
	 *
	 */
	private static class SexCheckProducer implements Producer<SexCheckResults> {
		private String[] bamFiles;
		private Logger log;
		private int index;

		public SexCheckProducer(String[] bamFiles, Logger log) {
			super();
			this.bamFiles = bamFiles;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < bamFiles.length;
		}

		@Override
		public Callable<SexCheckResults> next() {
			SexCheckWorker sexCheckWorker = new SexCheckWorker(bamFiles[index], log);
			index++;
			return sexCheckWorker;
		}

		@Override
		public void shutdown() {

		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

	}

	/**
	 * Computes the counts for sex chromosomes
	 *
	 */
	private static class SexCheckWorker implements Callable<SexCheckResults> {
		private String bamFile;
		private Logger log;

		public SexCheckWorker(String bamFile, Logger log) {
			super();
			this.bamFile = bamFile;
			this.log = log;
		}

		@Override
		public SexCheckResults call() throws Exception {
			SexCheckResults sexCheckResults = getSexCheckResults(bamFile, log);
			return sexCheckResults;
		}

	}

	private static SexCheckResults getSexCheckResults(String bamFile, Logger log) {
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
		SamReader reader = samReaderFactory.open(new File(bamFile));
		reader.indexing();
		SexCheckResults sexCheckResults = new SexCheckResults(0, 0, bamFile, bamFile);
		if (!reader.hasIndex()) {
			log.reportError("Error - the bam file " + bamFile + " must have a \".bai\" index file associated with it, halting");
			try {
				reader.close();
			} catch (IOException e) {
				log.reportIOException(bamFile);
				e.printStackTrace();
			}
		} else {
			String seqX = Positions.getChromosomeUCSC((byte) X, true);
			String seqY = Positions.getChromosomeUCSC((byte) Y, true);
			SAMFileHeader header = reader.getFileHeader();
			String sample = header.getReadGroups().get(0).getSample();
			log.reportTimeInfo("Computing sex checks for sample " + sample + " in file " + bamFile);
			int refX = header.getSequenceIndex(seqX);
			int refY = header.getSequenceIndex(seqY);

			QueryInterval qX = new QueryInterval(refX, 0, -1);
			QueryInterval qY = new QueryInterval(refY, 0, -1);
			sexCheckResults = getCountsForSexChr(new QueryInterval[] { qX, qY }, reader, bamFile, sample, log);
			try {
				reader.close();
			} catch (IOException e) {
				log.reportIOException(bamFile);
				e.printStackTrace();
			}
		}

		return sexCheckResults;
	}

	private static SexCheckResults getCountsForSexChr(QueryInterval[] qInterval, SamReader reader, String bamFile, String sample, Logger log) {
		SAMRecordIterator sIterator = reader.query(qInterval, false);
		int numXReads = 0;
		int numYReads = 0;
		int totalReads = 0;
		int goodReads = 0;
		while (sIterator.hasNext()) {
			totalReads++;
			SAMRecord samRecord = sIterator.next();
			if (goodRead(samRecord)) {
				goodReads++;
				byte chr = Positions.chromosomeNumber(samRecord.getReferenceName(), log);
				if (chr == X) {
					numXReads++;
				} else if (chr == Y) {
					numYReads++;
				} else {
					log.reportTimeError("Invalid chromosome found in " + bamFile + ", halting");
					return new SexCheckResults(0, 0, bamFile, sample);
				}
			}
			if (totalReads % 1000000 == 0) {
				log.reportTimeInfo("Sample: " + sample + " Read " + totalReads + " from chrs " + X + " and " + Y + ", " + goodReads + " passed standard filter, " + numXReads + " chr " + X + ", " + numYReads + " chr " + Y);
			}
		}
		log.reportTimeInfo("Sample: " + sample + " Finished with " + totalReads + " from chrs " + X + " and " + Y + ", " + goodReads + " passed standard filter, " + numXReads + " chr " + X + ", " + numYReads + " chr " + Y);

		return new SexCheckResults(numXReads, numYReads, bamFile, sample);
	}

	/**
	 * 
	 * @return true if the read is valid, not a duplicate, is primary, has a proper pair, and passes the MAPQ (trying to avoid pseudoautosomal)
	 */
	private static boolean goodRead(SAMRecord samRecord) {
		if (samRecord.isValid() != null) {
			return false;
		}
		if (samRecord.getDuplicateReadFlag()) {
			return false;
		}
		if (samRecord.isSecondaryOrSupplementary()) {
			return false;
		}
		if (!samRecord.getProperPairFlag()) {
			return false;
		}
		if (samRecord.getMappingQuality() == 255 || samRecord.getMappingQuality() < MAPQ_FILTER) {
			return false;
		}

		return true;
	}

	private static class SexCheckResults {
		private int numXReads;
		private int numYReads;
		private String bamFile;
		private String sample;

		private SexCheckResults(int numXReads, int numYReads, String bamFile, String sample) {
			super();
			this.numXReads = numXReads;
			this.numYReads = numYReads;
			this.bamFile = bamFile;
			this.sample = sample;
		}

		private String getSample() {
			return sample;
		}

		private int getNumXReads() {
			return numXReads;
		}

		private int getNumYReads() {
			return numYReads;
		}

		private String getBamFile() {
			return bamFile;
		}

		private double getPropX() {
			double total = (double) numXReads + numYReads;
			if (total > 0) {
				return (double) numXReads / total;
			} else {
				return 0;
			}
		}

		private double getPropY() {
			double total = (double) numXReads + numYReads;
			if (total > 0) {
				return (double) numYReads / total;
			} else {
				return 0;
			}
		}

	}

	/**
	 * @param dir
	 *            containing bam files
	 * @param fullPathTooutput
	 *            where the checks will be reported
	 * @param numThreads
	 * @param log
	 */
	public static void checkSex(String dir, String fullPathTooutput, int numThreads, Logger log) {
		String[] bamFiles = Files.listFullPaths(dir, ".bam", false);
		if (bamFiles.length < 1) {
			log.reportTimeError("Did not find any bam files in directory " + dir);
		} else {
			SexCheck sexCheck = new SexCheck(bamFiles, numThreads, log);
			sexCheck.checkSex(fullPathTooutput);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String directory = "bams/";
		String output = "bams/sexCheck.txt";
		int numThreads = 4;
		Logger log;

		String usage = "\n" + "seq.qc.SexCheck requires 0-1 arguments\n";
		usage += "   (1) directory of bam files to check (i.e. dir=" + directory + " (default))\n" + "";
		usage += "   (2) full path to output file(i.e. out=" + output + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(3, numThreads);

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				directory = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				output = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			new File(ext.parseDirectoryOfFile(output)).mkdirs();
			log = new Logger(ext.rootOf(output, false) + ".log");
			checkSex(directory, output, numThreads, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
