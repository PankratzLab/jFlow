package seq.cnv;

import java.util.ArrayList;
import java.util.concurrent.Callable;

import cnv.var.CNVariant;
import cnv.var.LocusSet;
import seq.manage.BEDFileReader;
import seq.manage.BEDFileReader.BEDFeatureSeg;
import seq.manage.BamPile;
import seq.manage.BamPileUp;
import seq.manage.BamSegPileUp;
import seq.manage.BedOps;
import seq.manage.ReferenceGenome;
import seq.manage.BamPileUp.PILE_TYPE;
import seq.qc.FilterNGS;
import seq.qc.FilterNGS.SAM_FILTER_TYPE;
import seq.qc.Mappability;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.WorkerTrain;
import common.ext;
import common.ArraySpecialList.ArrayIntList;
import common.WorkerTrain.Producer;
import filesys.Segment;

/**
 * Currently geared toward providing qcMetrics for ExomeDepthCalls from raw bam data
 *
 */
public class CnvBamQC {

	public static void qcCNVs(String bams, String cnvFile, String mappabilityFile, String callSubsetBed, String referenceGenomeFasta, int numThreads, Logger log) {
		ReferenceGenome referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);

		LocusSet<CNVariant> cnLocusSet = CNVariant.loadLocSet(cnvFile, log);
		Mappability<CNVariant> mappability = new Mappability<CNVariant>(cnLocusSet, mappabilityFile, callSubsetBed, log);
		String[] bamFiles = null;
		if (Files.isDirectory(bams)) {
			bamFiles = Files.listFullPaths(bams, ".bam", false);
		} else {
			bamFiles = HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, true);
		}
		log.reportTimeInfo("Found " + bamFiles.length + " bam files to search");
		CallSplit callSplit = new CallSplit(cnLocusSet, callSubsetBed, log);
		callSplit.matchAndSplit();
		PileupProducer producer = new PileupProducer(bamFiles, referenceGenomeFasta, log, callSplit);
		WorkerTrain<BamPile[]> train = new WorkerTrain<BamPile[]>(producer, numThreads, 2, log);
		while (train.hasNext()) {
			train.next();
		}
	}

	private static class CallSplit {
		private LocusSet<CNVariant> cnLocusSet;
		private BEDFileReader callSubsetBedReader;
		private Segment[] segsToSearch;
		private Logger log;
		private ArrayIntList[] matched;
		private boolean[] hadProblem;

		public CallSplit(LocusSet<CNVariant> cnLocusSet, String callSubsetBed, Logger log) {
			super();
			this.cnLocusSet = cnLocusSet;
			if (BedOps.verifyBedIndex(callSubsetBed, log)) {
				this.callSubsetBedReader = new BEDFileReader(callSubsetBed, true);

			}
			this.hadProblem = Array.booleanArray(cnLocusSet.getLoci().length, false);
			this.log = log;
			this.matched = new ArrayIntList[cnLocusSet.getLoci().length];
		}

		private void matchAndSplit() {
			int numProblems = 0;
			ArrayList<Segment> tmpSplit = new ArrayList<Segment>();
			int currentIndex = 0;
			for (int i = 0; i < cnLocusSet.getLoci().length; i++) {
				matched[i] = new ArrayIntList(100);
				LocusSet<BEDFeatureSeg> segs = callSubsetBedReader.loadSegsFor(cnLocusSet.getLoci()[i], log);
				Segment[] overlaps = segs.getOverLappingLoci(cnLocusSet.getLoci()[i]);
				for (int j = 0; j < overlaps.length; j++) {
					tmpSplit.add(overlaps[j].getUnion(cnLocusSet.getLoci()[i], log));
					matched[i].add(currentIndex);
					currentIndex++;
				}
				if (matched[i].size() != cnLocusSet.getLoci()[i].getNumMarkers()) {
					hadProblem[i] = true;
					numProblems++;
					// log.reportTimeWarning("found " + matched[i].size() + " overlapping regions and should have found " + cnLocusSet.getLoci()[i].getNumMarkers() + " for " + cnLocusSet.getLoci()[i].toAnalysisString());
				}
			}
			if (numProblems > 0) {
				log.reportTimeWarning(numProblems + " could not be perfectly matched to an identical number of exons");
			}
			this.segsToSearch = tmpSplit.toArray(new Segment[tmpSplit.size()]);
		}

		public Segment[] getSegsToSearch() {
			return segsToSearch;
		}

	}

	private static class PileupProducer implements Producer<BamPile[]> {
		private int index;
		private String[] bamFiles;
		private String referenceGenomeFasta;
		private Logger log;
		private CallSplit callSplit;

		public PileupProducer(String[] bamFiles, String referenceGenomeFasta, Logger log, CallSplit callSplit) {
			super();
			this.bamFiles = bamFiles;
			this.referenceGenomeFasta = referenceGenomeFasta;
			this.log = log;
			this.callSplit = callSplit;
		}

		@Override
		public boolean hasNext() {
			return index < bamFiles.length;
		}

		@Override
		public Callable<BamPile[]> next() {
			PileUpWorker worker = new PileUpWorker(bamFiles[index], referenceGenomeFasta, log, callSplit);
			index++;
			return worker;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static class PileUpWorker implements Callable<BamPile[]> {
		private String bamFile;
		private String referenceGenomeFasta;
		private Logger log;
		private CallSplit callSplit;

		public PileUpWorker(String bamFile, String referenceGenomeFasta, Logger log, CallSplit callSplit) {
			super();
			this.bamFile = bamFile;
			this.referenceGenomeFasta = referenceGenomeFasta;
			this.log = log;
			this.callSplit = callSplit;
		}

		@Override
		public BamPile[] call() throws Exception {
			BamSegPileUp bamSegPileUp = new BamSegPileUp(bamFile, new ReferenceGenome(referenceGenomeFasta, log), new Segment[] { callSplit.getSegsToSearch()[0] }, log);
			ArrayList<BamPile> bamPiles = new ArrayList<BamPile>();
			while (bamSegPileUp.hasNext()) {
				BamPile bamPile = bamSegPileUp.next();
				System.out.println(bamPile.getBin().getUCSClocation() + Array.toStr(bamPile.getCounts()));
			}
			System.exit(1);
			return bamPiles.toArray(new BamPile[bamPiles.size()]);
		}
	}

	public static void main(String[] args) {
		String referenceGenomeFasta = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
		String mappabilityFile = "/home/pankrat2/public/bin/ref/mapability/hg19/wgEncodeCrgMapabilityAlign100mer.bedGraph";
		String cnvFile = "/home/tsaim/shared/Project_Tsai_Project_021/ExomeDepth/results/ExomeDepth.all.cnvs";
		String callSubsetBed = "/panfs/roc/groups/5/pankrat2/public/bin/ExomeDepth/exons.hg19.sort.chr.bed";
		String bams = "/home/tsaim/shared/Project_Tsai_Project_021/bam/";
		Logger log = new Logger(ext.rootOf(cnvFile, false) + ".qc.log");
		int numThreads = 1;
		qcCNVs(bams, cnvFile, mappabilityFile, callSubsetBed, referenceGenomeFasta, numThreads, log);
	}

	public static void maina(String[] args) {
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
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
