package seq.analysis;

import filesys.LocusSet;
import filesys.Segment;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import seq.manage.BEDFileReader;
import seq.manage.BamOps;
import seq.manage.BedOps;
import seq.manage.ReferenceGenome;
import seq.manage.BamOps.BamIndexStats;
import common.Array;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;

/**
 * @author lane0212 Inspired by Kendall
 */
public class MitoSeqCN { 
	// TODO, need to adopt to WGS 
	public static void run(String fileOfBams, String outDir, String captureBed, String referenceGenomeFasta,boolean chr, int numthreads) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "mtDNACN.log");
		String[] bams = HashVec.loadFileToStringArray(fileOfBams, false, new int[] { 0 }, true);
		log.reportTimeInfo("Detected " + bams.length + " bam files");
		ReferenceGenome referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);
		System.out.println(captureBed);
		BedOps.verifyBedIndex(captureBed, log);
		BEDFileReader readerCapture = new BEDFileReader(captureBed, false);
		LocusSet<Segment> genomeBinsMinusBinsCaputure = referenceGenome.getBins(20000).removeThese(readerCapture.loadAll(log).getStrictSegmentSet(), 21000).autosomal(true, log);// essentially remove the neighbor
		readerCapture.close();
		log.reportTimeInfo(genomeBinsMinusBinsCaputure.getBpCovered() + " bp covered by reference bins int the anti-on-target regions");

		if (!referenceGenome.hasContig("chrM") || !referenceGenome.hasContig("chrX") || !referenceGenome.hasContig("chrY")) {
			throw new IllegalArgumentException("Required contig chrM,chrX,or chrY missing from " + referenceGenomeFasta);
		} else {
			int mitoLength = referenceGenome.getContigLength("chrM");
			log.reportTimeInfo("Mitochondrial genome length = " + mitoLength);

			MitoCNProducer producer = new MitoCNProducer(bams, referenceGenome, genomeBinsMinusBinsCaputure, outDir, chr, log);
			WorkerTrain<MitoCNResult> train = new WorkerTrain<MitoSeqCN.MitoCNResult>(producer, numthreads, numthreads, log);
			String output = outDir + ext.rootOf(fileOfBams) + "_mtDNACN.summary.txt";
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println(Array.toStr(MitoCNResult.header));
				while (train.hasNext()) {
					MitoCNResult result = train.next();
					if (result != null) {
						log.reportTimeInfo(Array.toStr(result.getResult()));
						writer.println(Array.toStr(result.getResult()));
					}

				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + output);
				log.reportException(e);
			}

		}

	}

	private static class MitoCNResult {
		private static final String[] header = new String[] { "Sample", "NumMitoReads", "TotalAlignedReads", "XReads", "YReads", "AutosomalOnTargetAlignedReads", "OffTargetReads", "MitoLen", "OffTLen" };
		private String sample;
		private int numMitoReads, numXReads, numYReads, autosomalOnTargetReads, offTargetReads, mitoLen;
		private long offTLen;
		private BamIndexStats bamIndexStats;
		private String outBam;

		public MitoCNResult(String sample, int numMitoReads, int numXReads, int numYReads, int offTargetReads, int mitoLen, long offTLen, BamIndexStats bamIndexStats, String outBam) {
			super();
			this.sample = sample;
			this.numMitoReads = numMitoReads;
			this.numXReads = numXReads;
			this.numYReads = numYReads;
			this.mitoLen = mitoLen;
			this.offTLen = offTLen;
			this.offTargetReads = offTargetReads;
			this.bamIndexStats = bamIndexStats;
			this.outBam = outBam;
			this.autosomalOnTargetReads = bamIndexStats.getAlignedRecordCount();
			autosomalOnTargetReads -= numMitoReads;
			autosomalOnTargetReads -= numXReads;
			autosomalOnTargetReads -= numYReads;

		}

		private String[] getResult() {

			ArrayList<String> result = new ArrayList<String>();
			result.add(sample);
			result.add(numMitoReads + "");
			result.add(bamIndexStats.getAlignedRecordCount() + "");
			result.add(numXReads + "");
			result.add(numYReads + "");
			result.add(autosomalOnTargetReads + "");
			result.add(offTargetReads + "");
			result.add(mitoLen + "");
			result.add(offTLen + "");

			return Array.toStringArray(result);

		}

	}

	private static class MitoCNWorker implements Callable<MitoCNResult> {
		private String bam;
		private String outDir;
		private int mitoLength, xLength, yLength;
		private LocusSet<Segment> genomeBinsMinusBinsCaputure;
		private boolean chr;
		private Logger log;

		public MitoCNWorker(String bam, LocusSet<Segment> genomeBinsMinusBinsCaputure, String outDir, int mitoLength, int xLength, int yLength, boolean chr, Logger log) {
			super();
			this.bam = bam;
			this.genomeBinsMinusBinsCaputure = genomeBinsMinusBinsCaputure;
			this.outDir = outDir;
			this.mitoLength = mitoLength;
			this.xLength = xLength;
			this.yLength = yLength;
			this.chr = chr;
			this.log = log;

		}

		@Override
		public MitoCNResult call() throws Exception {

			try {
				String sample = BamOps.getSampleName(bam);
				log.reportTimeInfo("Processing sample " + sample);
				String outputMTBam = outDir + ext.addToRoot(ext.removeDirectoryInfo(bam), ".chrM");

				SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);

				SAMFileWriter sAMFileWriter = new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reader.getFileHeader(), true, new File(outputMTBam));

				ArrayList<Segment> toSearch = new ArrayList<Segment>();
				toSearch.add(new Segment("chrMT", 0, mitoLength + 1));
				toSearch.add(new Segment("chrX", 0, xLength + 1));
				toSearch.add(new Segment("chrY", 0, yLength + 1));
				for (Segment segment : toSearch) {
					log.reportTimeInfo("Will search : " + segment.getUCSClocation());
				}
				QueryInterval[] queryInterestIntervals = BamOps.convertSegsToQI(toSearch.toArray(new Segment[toSearch.size()]), reader.getFileHeader(), 0, true, chr, log);
				SAMRecordIterator sIterator = reader.query(queryInterestIntervals, false);
				int numMitoReads = 0;
				int numXReads = 0;
				int numYReads = 0;
				int numOffTarget = 0;
				while (sIterator.hasNext()) {
					SAMRecord samRecord = sIterator.next();

					if (!samRecord.getReadUnmappedFlag()) {
						if (samRecord.getContig().equals("chrM")) {
							sAMFileWriter.addAlignment(samRecord);
							numMitoReads++;
						} else if (samRecord.getContig().equals("chrX")) {
							numXReads++;
						} else if (samRecord.getContig().equals("chrY")) {
							numYReads++;
						} else {
							throw new IllegalArgumentException("Invalid contig " + samRecord.getContig());
						}
					}
				}
				sIterator.close();
				sAMFileWriter.close();

				QueryInterval[] offTargetIntervalse = BamOps.convertSegsToQI(genomeBinsMinusBinsCaputure.getLoci(), reader.getFileHeader(), 0, true, true, log);
				sIterator = reader.query(offTargetIntervalse, false);
				while (sIterator.hasNext()) {
					SAMRecord samRecord = sIterator.next();
					if (!samRecord.getReadUnmappedFlag()) {
						numOffTarget++;

						if (numOffTarget % 1000000 == 0) {
							log.reportTimeInfo("Processing off target reads for sample " + sample + " , found " + numOffTarget);
						}
					}

				}

				BamIndexStats bamIndexStats = BamOps.getBamIndexStats(reader);

				return new MitoCNResult(sample, numMitoReads, numXReads, numYReads, numOffTarget, mitoLength, genomeBinsMinusBinsCaputure.getBpCovered(), bamIndexStats, outputMTBam);
			} catch (Exception e) {
				log.reportTimeError("Could not process " + bam);
				log.reportException(e);
				return null;
			}
		}
	}

	private static class MitoCNProducer implements Producer<MitoCNResult> {
		private String[] bams;
		private String outDir;
		private int index;
		private int mitoLength, xLength, yLength;
		private LocusSet<Segment> genomeBinsMinusBinsCaputure;
		private boolean chr;
		private Logger log;

		public MitoCNProducer(String[] bams, ReferenceGenome referenceGenome, LocusSet<Segment> genomeBinsMinusBinsCaputure, String outDir,boolean chr, Logger log) {
			super();
			this.bams = bams;
			this.outDir = outDir;
			this.mitoLength = referenceGenome.getContigLength("chrM");
			this.xLength = referenceGenome.getContigLength("chrX");
			this.yLength = referenceGenome.getContigLength("chrY");
			this.genomeBinsMinusBinsCaputure = genomeBinsMinusBinsCaputure;
			this.index = 0;
			this.chr =chr;
			this.log = log;
		}

		@Override
		public boolean hasNext() {
			return index < bams.length;
		}

		@Override
		public Callable<MitoCNResult> next() {
			String currentBam = bams[index];

			index++;
			return new MitoCNWorker(currentBam, genomeBinsMinusBinsCaputure, outDir, mitoLength, xLength, yLength, chr, log);
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String fileOfBams = "fileOfBams.txt";
		String outDir = "mitoWES/";
		int numthreads = 24;
		String referenceGenome = "hg19_canonical.fa";
		String captureBed = "AgilentCaptureRegions.txt";
		String logfile = null;

		String usage = "\n" +
				"seq.analysis.mitoSeqCN requires 0-1 arguments\n" +
				"   (1) file of Bams (i.e. bams=" + fileOfBams + " (default))\n" +
				"   (2) output directory (i.e. outDir=" + outDir + " (default))\n" +
				"   (3) number of threads (i.e. " + PSF.Ext.getNumThreadsCommand(3, numthreads) + "\n" +
				"   (4) reference genome (i.e. ref=" + outDir + " (default))\n" +
				"   (5)  capture Regions (i.e. captureBed=" + captureBed + " (default))\n" +

				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("bams=")) {
				fileOfBams = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outDir=")) {
				outDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("ref=")) {
				referenceGenome = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("captureBed=")) {
				captureBed = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
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

			run(fileOfBams, outDir, captureBed, referenceGenome, true, numthreads);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
