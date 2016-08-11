package org.genvisis.seq.analysis;

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
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BEDFileReader;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.BedOps;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.BamOps.BamIndexStats;

/**
 * @author lane0212 Inspired by Kendall
 */
public class MitoSeqCN {

	/**
	 * Used to determine mitochondrial sequences, and non-autosomal X and Ys,
	 *
	 */
	public enum BUILD_PARAMS {
		GRCH37("MT", "X", "Y"), HG19("chrMT", "chrX", "chrY");

		private String mitoContig;
		private String xContig;
		private String yContig;

		private BUILD_PARAMS(String mitoContig, String xContig, String yContig) {
			this.mitoContig = mitoContig;
			this.xContig = xContig;
			this.yContig = yContig;
		}

	}

	/**
	 * @param fileOfBams
	 *            bam files to estimate mtDNA CN for
	 * @param outDir
	 *            output directory for results
	 * @param captureBed
	 *            defining targeted capture regions
	 * @param referenceGenomeFasta
	 *            reference genomve
	 * @param params
	 * @param numthreads
	 */
	public static List<MitoCNResult> run(String fileOfBams, String outDir, String captureBed,
			String referenceGenomeFasta, BUILD_PARAMS params, int numthreads) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "mtDNACN.log");
		String[] bams = HashVec.loadFileToStringArray(fileOfBams, false, new int[] { 0 }, true);
		log.reportTimeInfo("Detected " + bams.length + " bam files");
		ReferenceGenome referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);
		BedOps.verifyBedIndex(captureBed, log);
		LocusSet<Segment> genomeBinsMinusBinsCaputure = referenceGenome.getBins(20000).autosomal(true, log);
		if (captureBed != null) {// Should only be used for
			BEDFileReader readerCapture = new BEDFileReader(captureBed, false);

			genomeBinsMinusBinsCaputure = genomeBinsMinusBinsCaputure
					.removeThese(readerCapture.loadAll(log).getStrictSegmentSet(), 21000).autosomal(true, log);
			readerCapture.close();
		} else {
			log.reportTimeWarning("No capture targets defined, assuming this is WGS");
		}
		log.reportTimeInfo(genomeBinsMinusBinsCaputure.getBpCovered() + " bp covered by reference bin regions");
		if (!referenceGenome.hasContig(params.mitoContig) || !referenceGenome.hasContig(params.xContig)
				|| !referenceGenome.hasContig(params.yContig)) {
			throw new IllegalArgumentException("Required contig for " + params + " is missing ( " + params.mitoContig
					+ " ," + params.xContig + ", " + params.yContig + " from " + referenceGenomeFasta);
		} else {
			int mitoLength = referenceGenome.getContigLength(params.mitoContig);
			log.reportTimeInfo("Mitochondrial genome length = " + mitoLength);

			String output = outDir + ext.rootOf(fileOfBams) + "_mtDNACN.summary.txt";

			MitoCNProducer producer = new MitoCNProducer(bams, referenceGenome, genomeBinsMinusBinsCaputure, outDir,
					params, log);
			WorkerTrain<MitoCNResult> train = new WorkerTrain<MitoSeqCN.MitoCNResult>(producer, numthreads, numthreads,
					log);
			ArrayList<MitoCNResult> results = new ArrayList<MitoSeqCN.MitoCNResult>();
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println(Array.toStr(MitoCNResult.header));
				while (train.hasNext()) {
					MitoCNResult result = train.next();
					if (result != null) {
						log.reportTimeInfo(Array.toStr(result.getResult()));
						writer.println(Array.toStr(result.getResult()));
						results.add(result);
					}

				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + output);
				log.reportException(e);
			}
			return results;

		}

	}

	/**
	 * Stores mtDNA CN estimation results for NGS data
	 *
	 */
	public static class MitoCNResult {
		private static final String[] header = new String[] { "Sample", "NumMitoReads", "TotalAlignedReads", "XReads",
				"YReads", "AutosomalOnTargetAlignedReads", "OffTargetReads", "MitoLen", "OffTLen", "MTBamFile" };
		private String sample;
		private int numMitoReads;
		private int numXReads;
		private int numYReads;
		private int autosomalOnTargetReads;
		private int offTargetReads;
		private int mitoLen;
		private long offTLen;
		private BamIndexStats bamIndexStats;
		private String outBam;

		private MitoCNResult(String sample, int numMitoReads, int numXReads, int numYReads, int offTargetReads,
				int mitoLen, long offTLen, BamIndexStats bamIndexStats, String outBam) {
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
			result.add(Integer.toString(numMitoReads));
			result.add(Integer.toString(bamIndexStats.getAlignedRecordCount()));
			result.add(Integer.toString(numXReads));
			result.add(Integer.toString(numYReads));
			result.add(Integer.toString(autosomalOnTargetReads));
			result.add(Integer.toString(offTargetReads));
			result.add(Integer.toString(mitoLen));
			result.add(Long.toString(offTLen));
			result.add(outBam);

			return Array.toStringArray(result);

		}

	}

	private static class MitoCNWorker implements Callable<MitoCNResult> {
		private String bam;
		private String outDir;
		private int mitoLength;
		private int xLength;
		private int yLength;
		private LocusSet<Segment> genomeBinsMinusBinsCaputure;
		private BUILD_PARAMS params;
		private Logger log;

		private MitoCNWorker(String bam, LocusSet<Segment> genomeBinsMinusBinsCaputure, String outDir, int mitoLength,
				int xLength, int yLength, BUILD_PARAMS params, Logger log) {
			super();
			this.bam = bam;
			this.genomeBinsMinusBinsCaputure = genomeBinsMinusBinsCaputure;
			this.outDir = outDir;
			this.mitoLength = mitoLength;
			this.xLength = xLength;
			this.yLength = yLength;
			this.params = params;
			this.log = log;

		}

		@Override
		public MitoCNResult call() throws Exception {

			try {
				String sample = BamOps.getSampleName(bam);
				log.reportTimeInfo("Processing sample " + sample);
				String outputMTBam = outDir + ext.addToRoot(ext.removeDirectoryInfo(bam), ".chrM");

				SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);

				SAMFileWriter sAMFileWriter = new SAMFileWriterFactory().setCreateIndex(true)
						.makeSAMOrBAMWriter(reader.getFileHeader(), true, new File(outputMTBam));

				ArrayList<Segment> toSearch = new ArrayList<Segment>();
				toSearch.add(new Segment(params.mitoContig, 0, mitoLength + 1));
				toSearch.add(new Segment(params.xContig, 0, xLength + 1));
				toSearch.add(new Segment(params.yContig, 0, yLength + 1));
				for (Segment segment : toSearch) {
					log.reportTimeInfo("Will search : " + segment.getUCSClocation());
				}
				QueryInterval[] queryInterestIntervals = BamOps.convertSegsToQI(
						toSearch.toArray(new Segment[toSearch.size()]), reader.getFileHeader(), 0, true,
						params == BUILD_PARAMS.HG19, log);
				SAMRecordIterator sIterator = reader.query(queryInterestIntervals, false);
				int numMitoReads = 0;
				int numXReads = 0;
				int numYReads = 0;
				int numOffTarget = 0;
				while (sIterator.hasNext()) {
					SAMRecord samRecord = sIterator.next();

					if (!samRecord.getReadUnmappedFlag() && !samRecord.getDuplicateReadFlag()) {
						if (samRecord.getContig().equals(params.mitoContig)) {
							sAMFileWriter.addAlignment(samRecord);
							numMitoReads++;
						} else if (samRecord.getContig().equals(params.xContig)) {
							numXReads++;
						} else if (samRecord.getContig().equals(params.yContig)) {
							numYReads++;
						} else {
							throw new IllegalArgumentException("Invalid contig " + samRecord.getContig());
						}
					}
				}
				sIterator.close();
				sAMFileWriter.close();

				QueryInterval[] offTargetIntervalse = BamOps.convertSegsToQI(genomeBinsMinusBinsCaputure.getLoci(),
						reader.getFileHeader(), 0, true, params == BUILD_PARAMS.HG19, log);
				sIterator = reader.query(offTargetIntervalse, false);
				while (sIterator.hasNext()) {
					SAMRecord samRecord = sIterator.next();
					if (!samRecord.getReadUnmappedFlag()) {
						numOffTarget++;

						if (numOffTarget % 1000000 == 0) {
							log.reportTimeInfo("Processing off normalization-reads for sample " + sample + " , found "
									+ numOffTarget);
						}
					}

				}

				BamIndexStats bamIndexStats = BamOps.getBamIndexStats(reader);
				reader.close();
				return new MitoCNResult(sample, numMitoReads, numXReads, numYReads, numOffTarget, mitoLength,
						genomeBinsMinusBinsCaputure.getBpCovered(), bamIndexStats, outputMTBam);
			} catch (Exception e) {
				log.reportTimeError("Could not process " + bam);
				log.reportException(e);
				return null;
			}
		}
	}

	private static class MitoCNProducer extends AbstractProducer<MitoCNResult> {
		private String[] bams;
		private String outDir;
		private int index;
		private int mitoLength;
		private int xLength;
		private int yLength;
		private LocusSet<Segment> genomeBinsMinusBinsCaputure;
		private BUILD_PARAMS params;
		private Logger log;

		private MitoCNProducer(String[] bams, ReferenceGenome referenceGenome,
				LocusSet<Segment> genomeBinsMinusBinsCaputure, String outDir, BUILD_PARAMS params, Logger log) {
			super();
			this.bams = bams;
			this.outDir = outDir;
			this.mitoLength = referenceGenome.getContigLength(params.mitoContig);
			this.xLength = referenceGenome.getContigLength(params.xContig);
			this.yLength = referenceGenome.getContigLength(params.yContig);
			this.genomeBinsMinusBinsCaputure = genomeBinsMinusBinsCaputure;
			this.index = 0;
			this.params = params;
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
			return new MitoCNWorker(currentBam, genomeBinsMinusBinsCaputure, outDir, mitoLength, xLength, yLength,
					params, log);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String fileOfBams = "fileOfBams.txt";
		String outDir = "mitoWES/";
		int numthreads = 24;
		String referenceGenome = "hg19_canonical.fa";
		String captureBed = "AgilentCaptureRegions.txt";

		String usage = "\n" + "seq.analysis.mitoSeqCN requires 0-1 arguments\n" + "   (1) file of Bams (i.e. bams="
				+ fileOfBams + " (default))\n" + "   (2) output directory (i.e. outDir=" + outDir + " (default))\n"
				+ "   (3) number of threads (i.e. " + PSF.Ext.getNumThreadsCommand(3, numthreads) + "\n"
				+ "   (4) reference genome (i.e. ref=" + outDir + " (default))\n"
				+ "   (5)  capture Regions (i.e. captureBed=" + captureBed + " (default))\n" +

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
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {

			run(fileOfBams, outDir, captureBed, referenceGenome, BUILD_PARAMS.HG19, numthreads);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
