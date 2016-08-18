package org.genvisis.seq.qc.contamination;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.concurrent.Callable;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.BamPileUp.PILE_TYPE;
import org.genvisis.seq.manage.BamPileUp.bamPileWorker;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.SAM_FILTER_TYPE;
import org.genvisis.stats.Histogram.DynamicHistogram;

/**
 * Uses an mpileup approach to determine contamination levels of a bam file
 *
 */
public class BamContamination {

	/**
	 * Producer that defaults to contamination detection
	 *
	 */
	public static class BamContaminationProducer extends AbstractProducer<DynamicHistogram> {
		private Segment[] q;
		private FilterNGS filterNGS;
		private ReferenceGenome referenceGenome;
		private String[] bamFiles;
		private int index;
		private Logger log;

		public BamContaminationProducer(Segment[] q, FilterNGS filterNGS, ReferenceGenome referenceGenome, String[] bamFiles, Logger log) {
			super();
			this.q = q;
			this.filterNGS = filterNGS;
			this.referenceGenome = referenceGenome;
			this.bamFiles = bamFiles;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return index < bamFiles.length;
		}

		@Override
		public Callable<DynamicHistogram> next() {
			bamPileWorker worker = new bamPileWorker(bamFiles[index], q, filterNGS, referenceGenome, 1, PILE_TYPE.CONTAMINATION, SAM_FILTER_TYPE.GENOTYPE, log);
			index++;
			return worker;
		}
	}

	public static void runContam(String bams, String referenceGenomeFasta, String segFile, String pfbFile, FilterNGS filterNGS, int numthreads, Logger log) {
		String[] bamFiles = null;
		if (Files.isDirectory(bams)) {
			bamFiles = Files.listFullPaths(bams, ".bam", false);
		} else {
			bamFiles = HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, true);
		}
		// bamFiles = Array.subArray(bamFiles, 0, numthreads);
		Segment[] q = segFile == null ? null : Segment.loadRegions(segFile, 0, 1, 2, 0, true, true, true, 0);

		ReferenceGenome referenceGenome = referenceGenomeFasta == null ? null : new ReferenceGenome(referenceGenomeFasta, log);
		BamContaminationProducer producer = new BamContaminationProducer(q, filterNGS, referenceGenome, bamFiles, log);
		log.reportTimeInfo("Detected " + bamFiles.length + " bam files in " + bams);
		DynamicHistogram[] hists = new DynamicHistogram[bamFiles.length];
		WorkerTrain<DynamicHistogram> train = new WorkerTrain<DynamicHistogram>(producer, numthreads, numthreads, log);
		int index = 0;
		while (train.hasNext()) {
			hists[index] = train.next();
			index++;
		}
		String outputHist = bams + "contamSummary.txt";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outputHist));
			writer.print("PROP_REF_BIN");
			for (int i = 0; i < hists.length; i++) {
				writer.print("\t" + ext.rootOf(bamFiles[i]) + "_COUNTS");
			}
			writer.println();
			double[] bins = hists[0].getBins();
			for (int i = 0; i < bins.length; i++) {
				writer.print(bins[i]);
				for (int j = 0; j < hists.length; j++) {
					writer.print("\t" + hists[j].getCounts()[i]);
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outputHist);
			log.reportException(e);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String segFile = null;
		String referenceGenomeFasta = "hg19_canonical.fa";
		String bams = null;
		String pfbFile = null;
		double minPhred = 30;
		double minMapQ = 30;

		int minDepth = 20;
		int minAltDepth = -1;

		int numthreads = 4;
		String logfile = null;
		Logger log;
		String usage = "\n" + "seq.manage.BamPileUp requires 1-2 arguments\n";
		usage += "   (1) full path to a directory of *.bam files or a file listing bam files in the first column, to pile-up (i.e. bams=" + bams + " (no default))\n" + "";
		usage += "   (2) full path to a reference fasta  (i.e. ref= (no default))\n" + "";
		usage += "   (3) full path to a file of segments to subset the pile up  (i.e. segs= (no default))\n" + "";
		usage += "   (4) minimum phred score for a base pair to be piled  (i.e. minPhred=" + minPhred + " (default))\n" + "";
		usage += "   (5) minimum mapping quality score for a read to be piled  (i.e. minMapQ=" + minMapQ + " (default))\n" + "";
		usage += "   (6) minimum total depth for a position to be reported  (i.e. minDepth=" + minDepth + " (default))\n" + "";
		usage += "   (7) minimum alternate allele depth for a position to be reported  (i.e. minAltDepth=" + minAltDepth + " (default, no minimum))\n" + "";
		usage += "   (8) file listing pfbs for known variants  (i.e. pfb= (no default))\n" + "";

		usage += PSF.Ext.getNumThreadsCommand(9, numthreads);

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("bams=")) {
				bams = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("ref=")) {
				referenceGenomeFasta = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("segs=")) {
				segFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("pfb=")) {
				pfbFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("minPhred=")) {
				minPhred = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("minDepth=")) {
				minDepth = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("minAltDepth=")) {
				minAltDepth = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("minMapQ=")) {
				minMapQ = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
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
			logfile = Files.exists(bams) ? bams + "contam.log" : logfile;
			log = new Logger(logfile);
			FilterNGS filterNGS = new FilterNGS(minMapQ, minPhred, new int[] { minDepth, minAltDepth });
			runContam(bams, referenceGenomeFasta, segFile, pfbFile, filterNGS, numthreads, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
