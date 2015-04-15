package seq.manage;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.SamRecordFilter;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Callable;

import seq.qc.FilterNGS;
import stats.Histogram.DynamicHistogram;
import common.Array;
import common.Files;
import common.Logger;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;
import filesys.Segment;

public class BamPileUp {
	public enum PILE_TYPE {
		/**
		 * Will only report positions with alternate alleles, and demands a reference genome
		 */
		CONTAMINATION, /**
		 * Pileup to all positions passing filters supplied
		 */
		REGULAR;
	}

	private String bam;
	private FilterNGS filterNGS;
	private ReferenceGenome referenceGenome;
	private Segment[] intervals;
	private int binSize;
	private PILE_TYPE type;
	private Logger log;

	public BamPileUp(String bam, ReferenceGenome referenceGenome, int refBinSize, FilterNGS filterNGS, Segment[] intervals, PILE_TYPE type, Logger log) {
		super();
		this.bam = bam;
		this.binSize = refBinSize;
		this.log = log;
		this.filterNGS = filterNGS;
		this.intervals = intervals;
		this.referenceGenome = referenceGenome;
		this.type = type;
	}

	public DynamicHistogram pileUp() {
		if (type == PILE_TYPE.CONTAMINATION && referenceGenome == null) {
			String error = "A reference genome must be provided to detect reference alleles, cannot detect contamination in bam file " + bam;
			log.reportTimeError(error);
			return null;
		} else if (referenceGenome == null) {
			log.reportTimeWarning("A reference genome was not provided, reference and alternate alleles will not be recorded");
		}

		if (type == PILE_TYPE.CONTAMINATION && binSize != 1) {
			String error = "A bin size of one must be used for contamination detection on  " + bam;
			log.reportTimeError(error);
			return null;
		} else {
			log.reportTimeInfo("Using a bin size of " + binSize);
		}
		WorkerTrain<TmpBamPile> train = new WorkerTrain<BamPileUp.TmpBamPile>(null, 2, 200, log);// an extra thread should be about a half hour speed up per sample
		train.setAutoShutDown(false);
		TmpBamPileProducer tmpBamPileProducer = null;
		SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);
		DynamicHistogram histogram = new DynamicHistogram(0, 1, 2);
		QueryInterval[] queryIntervals = null;
		if (intervals != null) {
			log.reportTimeInfo("Optimizing " + intervals.length + " queries for pile up");
			queryIntervals = BamOps.convertSegsToQI(intervals, reader.getFileHeader(), 0, true, log);
			log.reportTimeInfo("Finished Optimizing " + intervals.length + " queries to " + queryIntervals.length + " intervals for pile up");
		}

		SAMRecordIterator sIterator = queryIntervals == null ? reader.iterator() : reader.query(queryIntervals, false);
		AggregateFilter filter = initializeFilters(filterNGS);
		Segment currentSegment = new Segment((byte) 0, 0, 0);
		ArrayList<BamPile> bamPiles = new ArrayList<BamPileUp.BamPile>();

		int totalReads = 0;
		int readsPiled = 0;
		int positionsPiled = 0;
		long time = System.currentTimeMillis();
		long totalTime = System.currentTimeMillis();

		String output = ext.rootOf(bam, false) + ".bamPile.txt";
		PrintWriter writer = Files.getAppropriateWriter(output);
		writer.println(BamPile.getHeader());

		while (sIterator.hasNext()) {
			totalReads++;
			SAMRecord samRecord = sIterator.next();
			if (!filter.filterOut(samRecord)) {
				readsPiled++;
				if (readsPiled % 100000 == 0) {
					log.reportTimeInfo("~" + readsPiled + " of " + totalReads + " total reads piled to " + positionsPiled + " positions (" + ext.getTimeElapsed(time) + ") " + SamRecordOps.getDisplayLoc(samRecord));
					time = System.currentTimeMillis();
				}
				BamPileInitializer bamPileInitializer = new BamPileInitializer(binSize, currentSegment, samRecord, log);
				while (bamPileInitializer.hasNext()) {
					bamPiles.add(bamPileInitializer.next());
				}
				currentSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
				tmpBamPileProducer = new TmpBamPileProducer(samRecord, currentSegment, bamPiles.toArray(new BamPile[bamPiles.size()]), filterNGS, log);
				train.setProducer(tmpBamPileProducer);
				bamPiles = new ArrayList<BamPileUp.BamPile>(bamPiles.size());
				while (train.hasNext()) {
					TmpBamPile tmpBamPile = train.next();
					if (tmpBamPile.overlapsCurrentRecord()) {
						bamPiles.add(tmpBamPile.getBamPile());// stays in this round
					} else {
						BamPile bamPile = tmpBamPile.getBamPile();
						if (referenceGenome != null || type == PILE_TYPE.CONTAMINATION) {
							bamPile.setReference(referenceGenome);
						}
						if (filterNGS.getReadDepthFilter() == null || bamPile.getTotalDepth(false, false) > filterNGS.getReadDepthFilter()[0]) {
							int altAlleleDepth = filterNGS.getReadDepthFilter().length > 1 ? filterNGS.getReadDepthFilter()[1] : -1;
							if (type == PILE_TYPE.REGULAR || (bamPile.hasAltAllele(log) && bamPile.hasOnlyOneAlt(log) && bamPile.getNumAlt(log) > altAlleleDepth && bamPile.getNumRef(log) > altAlleleDepth)) {
								positionsPiled++;
								histogram.addDataPointToHistogram(bamPile.getPropRef(log));
								writer.println(bamPile.getOuput(log));
							}
						}
					}
				}
			}
		}
		sIterator.close();
		train.shutdown();
		log.reportTimeInfo(readsPiled + " reads piled in total (" + ext.getTimeElapsed(totalTime) + ")");

		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		writer.close();
		return histogram;
	}

	private static AggregateFilter initializeFilters(FilterNGS filterNGS) {
		ArrayList<SamRecordFilter> filters = filterNGS.getStandardSAMRecordFilters();
		filters.add(filterNGS.getSamRecordMapQFilter(filterNGS.getMappingQualityFilter()));
		AggregateFilter filter = new AggregateFilter(filters);
		return filter;
	}

	/**
	 * Feeds up new bamPiles for the next base pairs determined by the next alignment
	 *
	 */
	private static class BamPileInitializer implements Iterator<BamPile> {
		private int binSize;
		private Segment previousBin;
		private Segment nextBin;
		private Logger log;
		private Segment samRecordSeg;
		private SAMRecord samRecord;

		public BamPileInitializer(int binSize, Segment previousBin, SAMRecord samRecord, Logger log) {
			super();
			this.binSize = binSize;
			this.previousBin = previousBin;
			this.samRecord = samRecord;
			this.samRecordSeg = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
			this.nextBin = samRecordSeg.overlaps(previousBin) ? getNextBin(previousBin, binSize) : scanToNext();
		}

		private Segment scanToNext() {
			if (samRecordSeg.getChr() > 0) {
				if (previousBin.getChr() != samRecordSeg.getChr()) {
					previousBin = new Segment(samRecordSeg.getChr(), 1, 1 + binSize);
				}
				while (!samRecordSeg.overlaps(previousBin) && previousBin.getStop() <= samRecordSeg.getStop()) {
					previousBin = new Segment(previousBin.getChr(), previousBin.getStop() + 1, previousBin.getStop() + binSize);
				}
			} else {
				log.reportTimeError("Could not find valid segment for " + samRecord.toString());
			}
			return previousBin;
		}

		private static Segment getNextBin(Segment seg, int binSize) {
			return new Segment(seg.getChr(), seg.getStop() + 1, seg.getStop() + binSize);
		}

		@Override
		public boolean hasNext() {
			return nextBin.overlaps(samRecordSeg) && nextBin.getStop() <= samRecordSeg.getStop();
		}

		@Override
		public BamPile next() {
			BamPile bamPile = new BamPile(nextBin);
			nextBin = getNextBin(nextBin, binSize);
			return bamPile;
			// TODO Auto-generated method stub
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

	}

	private static class TmpBamPileProducer implements Producer<TmpBamPile> {
		private SAMRecord samRecord;
		private Segment samRecordSegment;
		private BamPile[] bamPiles;
		private FilterNGS filterNGS;
		private Logger log;
		private int index = 0;

		private TmpBamPileProducer(SAMRecord samRecord, Segment samRecordSegment, BamPile[] bamPiles, FilterNGS filterNGS, Logger log) {
			super();
			this.samRecord = samRecord;
			this.bamPiles = bamPiles;
			this.filterNGS = filterNGS;
			this.samRecordSegment = samRecordSegment;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < bamPiles.length;
		}

		@Override
		public Callable<TmpBamPile> next() {
			TmpBamPile tmpBamPile = new TmpBamPile(samRecord, samRecordSegment, bamPiles[index], filterNGS, log);
			index++;
			return tmpBamPile;
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static class TmpBamPile implements Callable<TmpBamPile> {
		private SAMRecord samRecord;
		private BamPile bamPile;
		private FilterNGS filterNGS;
		private boolean overlaps;
		private Segment samRecordSegment;
		private Logger log;

		private TmpBamPile(SAMRecord samRecord, Segment samRecordSegment, BamPile bamPile, FilterNGS filterNGS, Logger log) {
			super();
			this.samRecord = samRecord;
			this.bamPile = bamPile;
			this.filterNGS = filterNGS;
			this.samRecordSegment = samRecordSegment;
			this.log = log;
		}

		@Override
		public TmpBamPile call() throws Exception {
			this.overlaps = samRecordSegment.overlaps(bamPile.getBin());
			if (overlaps) {
				bamPile.addRecord(samRecord, filterNGS.getPhreadScoreFilter(), log);
			} else {
				bamPile.summarize();
			}
			return this;
		}

		private BamPile getBamPile() {
			return bamPile;
		}

		private boolean overlapsCurrentRecord() {
			return overlaps;
		}

	}

	/**
	 * Stores the actual pileUp for each position found
	 *
	 */
	private static class BamPile {
		private static final String[] BASE_HEADER = new String[] { "UCSC", "REF", "NUM_REF", "NUM_ALT", "PROP_REF", "PROP_ALT" };
		private static final String[] COUNT_HEADER = new String[] { "NUMA", "NUMG", "NUMC", "NUMT", "NUMN", "NUMDEL", "NUMINS" };
		private static final String[] EXT_HEADER = new String[] { "MAPQ", "PHRED" };

		private Segment bin;
		private int[] counts;
		private String refAllele;

		private double[] avgMapQ;
		private double[] avgPhread;

		private static String getHeader() {
			String header = Array.toStr(BASE_HEADER);
			for (int i = 0; i < COUNT_HEADER.length; i++) {
				header += "\t" + COUNT_HEADER[i];
			}

			for (int i = 0; i < COUNT_HEADER.length; i++) {
				for (int j = 0; j < EXT_HEADER.length; j++) {
					header += "\t" + COUNT_HEADER[i] + "_" + EXT_HEADER[j];
				}
			}

			return header;
		}

		private BamPile(Segment bin) {
			super();
			this.bin = bin;
			this.counts = new int[7];// A,G,C,T,N,Del,Ins
			this.avgMapQ = new double[7];
			this.avgPhread = new double[7];
			this.refAllele = "NA";

		}

		private int getTotalDepth(boolean includeIndels, boolean includeNs) {
			boolean[] indelmask = new boolean[7];
			Arrays.fill(indelmask, true);
			if (!includeIndels) {
				indelmask[6] = false;
				indelmask[5] = false;
			}
			if (!includeNs) {
				indelmask[4] = false;
			}

			return Array.sum(Array.subArray(counts, indelmask));

		}

		private void setReference(ReferenceGenome referenceGenome) {

			String[] ref = referenceGenome.getSequenceFor(bin);
			refAllele = ref[0];
			if (ref.length > 1) {
				System.out.println(Array.toStr(ref));
			}

			if (ref.length > 1) {
				for (int i = 1; i < ref.length; i++) {
					refAllele += "/" + ref[i];
				}
			}
		}

		private String getOuput(Logger log) {
			double numRef = (double) getNumRef(log);
			double numAlt = (double) getNumAlt(log);
			String out = "";
			out += bin.getUCSClocation();
			out += "\t" + refAllele;
			out += "\t" + numRef;
			out += "\t" + numAlt;
			out += "\t" + getPropRef(log);
			out += "\t" + getPropAlt(log);
			out += "\t" + counts[0];
			out += "\t" + counts[1];
			out += "\t" + counts[2];
			out += "\t" + counts[3];
			out += "\t" + counts[4];
			out += "\t" + counts[5];
			out += "\t" + counts[6];
			out += "\t" + Array.toStr(avgMapQ);
			out += "\t" + Array.toStr(avgPhread);
			return out;
		}

		private double getPropRef(Logger log) {
			double numRef = (double) getNumRef(log);
			double numAlt = (double) getNumAlt(log);
			double total = numRef + numAlt;
			return numRef / total;
		}

		private double getPropAlt(Logger log) {
			double numRef = (double) getNumRef(log);
			double numAlt = (double) getNumAlt(log);
			double total = numRef + numAlt;
			return numAlt / total;
		}

		private int[] getAltCounts(Logger log) {
			boolean[] referenceMask = new boolean[7];
			Arrays.fill(referenceMask, true);
			referenceMask[6] = false;
			referenceMask[5] = false;

			if (refAllele.equals("A")) {
				referenceMask[0] = false;
			} else if (refAllele.equals("G")) {
				referenceMask[1] = false;
			} else if (refAllele.equals("C")) {
				referenceMask[2] = false;
			} else if (refAllele.equals("T")) {
				referenceMask[3] = false;
			} else if (refAllele.equals("N")) {
				referenceMask[4] = false;
			}
			if (Array.booleanArraySum(referenceMask) != 4) {
				log.reportTimeError("Invalid number of alternate allele possibilities, found " + Array.booleanArraySum(referenceMask) + " with ref allele" + refAllele);

			}
			return Array.subArray(counts, referenceMask);
		}

		private int getNumAlt(Logger log) {
			return Array.sum(getAltCounts(log));
		}

		private boolean hasOnlyOneAlt(Logger log) {
			return Array.countIf(getAltCounts(log), 0) == 3;// everthing else is 0
		}

		private int getNumRef(Logger log) {
			boolean[] referenceMask = new boolean[7];
			Arrays.fill(referenceMask, false);
			referenceMask[6] = false;
			referenceMask[5] = false;

			if (refAllele.equals("A")) {
				referenceMask[0] = true;
			} else if (refAllele.equals("G")) {
				referenceMask[1] = true;
			} else if (refAllele.equals("C")) {
				referenceMask[2] = true;
			} else if (refAllele.equals("T")) {
				referenceMask[3] = true;
			} else if (refAllele.equals("N")) {
				referenceMask[4] = true;
			}
			if (Array.booleanArraySum(referenceMask) != 1) {
				log.reportTimeError("Invalid number of reference allele possibilities");
			}
			return (Array.sum(Array.subArray(counts, referenceMask)));
		}

		private boolean hasAltAllele(Logger log) {
			return getNumAlt(log) > 0;
		}

		private void summarize() {
			for (int i = 0; i < counts.length; i++) {
				if (counts[i] > 0) {
					avgMapQ[i] = (double) avgMapQ[i] / counts[i];
					avgPhread[i] = (double) avgPhread[i] / counts[i];
				}
			}
		}

		private Segment getBin() {
			return bin;
		}

		private void addRecord(SAMRecord samRecord, double phredFilter, Logger log) {
			Segment samRecordSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
			Segment toPile = bin.getUnion(samRecordSegment, log);
			int mapQ = samRecord.getMappingQuality();
			if (mapQ == 255) {
				String error = "Detected invalid mapping quality (255)";
				log.reportTimeError(error);
				throw new IllegalArgumentException(error);
			} else {

				// System.out.println(toPile.getSize());
				String r = samRecord.getReadString();
				double[] p = SamRecordOps.getReadPhred(samRecord);
				// double[] p = new double[r.length()];

				int curRefBase = samRecord.getAlignmentStart();
				int curReadBase = 0;
				List<CigarElement> cigarEls = samRecord.getCigar().getCigarElements();
				for (CigarElement cigarEl : cigarEls) {
					if (curRefBase > toPile.getStop()) {
						break;

					}
					CigarOperator op = cigarEl.getOperator();
					for (int i = 0; i < cigarEl.getLength(); i++) {
						String base = null;
						if (curRefBase > toPile.getStop()) {
							break;

						}
						if (curRefBase >= toPile.getStart() && curRefBase <= toPile.getStop()) {
							if (p[curReadBase] > phredFilter) {
								base = r.charAt(curReadBase) + "";
							}
						}

						if (op.consumesReadBases() && op.consumesReferenceBases()) {
							if (base != null) {
								addRegBase(base, mapQ, p[curReadBase], log);
							}
							curRefBase++;
							curReadBase++;
						} else if (op.consumesReadBases()) {
							if (base != null) {
								addIns(mapQ, p[curReadBase]);
							}
							curReadBase++;
						} else if (op.consumesReferenceBases()) {
							if (base != null) {
								addDel(mapQ, p[curReadBase]);
							}
							curRefBase++;
						}
					}
				}
			}
		}

		private void addDel(int mapQ, double p) {
			counts[5]++;
			avgPhread[5] += p;
			avgMapQ[5] += mapQ;
		}

		private void addIns(int mapQ, double p) {
			counts[6]++;
			avgPhread[6] += p;
			avgMapQ[6] += mapQ;
		}

		private void addRegBase(String b, int mapQ, double p, Logger log) {
			if (b.equals("A")) {
				counts[0]++;
				avgPhread[0] += p;
				avgMapQ[0] += mapQ;
			} else if (b.equals("G")) {
				counts[1]++;
				avgPhread[1] += p;
				avgMapQ[1] += mapQ;
			} else if (b.equals("C")) {
				counts[2]++;
				avgPhread[2] += p;
				avgMapQ[2] += mapQ;
			} else if (b.equals("T")) {
				counts[3]++;
				avgPhread[3] += p;
				avgMapQ[3] += mapQ;
			} else if (b.equals("N")) {
				counts[4]++;
				avgPhread[4] += p;
				avgMapQ[4] += mapQ;
			} else {
				String error = "Invalid base " + b + " for regular base tracking";
				log.reportTimeError(error);
				throw new IllegalArgumentException(error);
			}
		}

	}

	public static class bamPileWorker implements Callable<DynamicHistogram> {
		private String bamFile;
		private Segment[] q;
		private FilterNGS filterNGS;
		private ReferenceGenome referenceGenome;
		private int binSize;
		private PILE_TYPE type;
		private Logger log;

		public bamPileWorker(String bamFile, Segment[] q, FilterNGS filterNGS, ReferenceGenome referenceGenome, int binSize, PILE_TYPE type, Logger log) {
			super();
			this.bamFile = bamFile;
			this.q = q;
			this.filterNGS = filterNGS;
			this.referenceGenome = referenceGenome;
			this.binSize = binSize;
			this.log = log;
			this.type = type;
		}

		@Override
		public DynamicHistogram call() throws Exception {
			BamPileUp pileUp = new BamPileUp(bamFile, referenceGenome, binSize, filterNGS, q, type, log);
			DynamicHistogram histogram = pileUp.pileUp();
			return histogram;
		}
	}
	
	
	
	public static class BamPileUpSummary{
		
		private DynamicHistogram proRefHistogram;
		private Logger log;
	}
}

//
// public static class BamPileProducer implements Producer<Boolean> {
// private Segment[] q;
// private FilterNGS filterNGS;
// private ReferenceGenome referenceGenome;
// private String[] bamFiles;
// private int index;
// private Logger log;
//
// public BamPileProducer(Segment[] q, FilterNGS filterNGS, ReferenceGenome referenceGenome, String[] bamFiles, Logger log) {
// super();
// this.q = q;
// this.filterNGS = filterNGS;
// this.referenceGenome = referenceGenome;
// this.bamFiles = bamFiles;
// this.log = log;
// this.index = 0;
// }
//
// @Override
// public boolean hasNext() {
// // TODO Auto-generated method stub
// return index < bamFiles.length;
// }
//
// @Override
// public Callable<Boolean> next() {
// bamPileWorker worker = new bamPileWorker(bamFiles[index], q, filterNGS, referenceGenome, log);
// index++;
// return worker;
// }
//
// @Override
// public void remove() {
// // TODO Auto-generated method stub
//
// }
//
// @Override
// public void shutdown() {
// // TODO Auto-generated method stub
//
// }
// }
//
// public static class bamPileWorker implements Callable<Boolean> {
// private String bamFile;
// private Segment[] q;
// private FilterNGS filterNGS;
// private ReferenceGenome referenceGenome;
// private Logger log;
//
// public bamPileWorker(String bamFile, Segment[] q, FilterNGS filterNGS, ReferenceGenome referenceGenome, Logger log) {
// super();
// this.bamFile = bamFile;
// this.q = q;
// this.filterNGS = filterNGS;
// this.referenceGenome = referenceGenome;
// this.log = log;
// }
//
// @Override
// public Boolean call() throws Exception {
// BamPileUp pileUp = new BamPileUp(bamFile, referenceGenome, 1, filterNGS, q, log);
// pileUp.pileUp();
// // TODO Auto-generated method stub
// return true;
// }
//
// }
//
// public static void testDir(Logger log) {
// String segFile = "C:/bin/Agilent/captureLibraries/SureSelectHumanAllExonV5UTRs/AgilentCaptureRegions_chr1.txt";
// String ref = "C:/bin/ref/hg19_canonical.fa";
// String testbamDir = "D:/data/Project_Tsai_Project_021/testPileUp/";
// String[] bamFiles = Files.listFullPaths(testbamDir, ".bam", false);
// Segment[] q = segFile == null ? null : Segment.loadRegions(segFile, 0, 1, 2, 0, true, true, true, 100);
// FilterNGS filterNGS = new FilterNGS(30, 30, new int[] { 15 });
// ReferenceGenome referenceGenome = ref == null ? null : new ReferenceGenome(ref, log);
// BamPileProducer producer = new BamPileProducer(q, filterNGS, referenceGenome, bamFiles, log);
// log.reportTimeInfo("Detected " + bamFiles.length + " bam files in directory " + testbamDir);
// WorkerTrain<Boolean> train = new WorkerTrain<Boolean>(producer, 4, 4, log);
//
// while (train.hasNext()) {
// train.next();
// }
// }
//
// public static void main(String[] args) {
// int numArgs = args.length;
// String segFile = null;
// String referenceGenomeFasta = "hg19_canonical.fa";
// String bamDir = null;
// double minPhred = 30;
// double minMapQ = 30;
//
// String logfile = null;
// Logger log;
//
// String usage = "\n" + "seq.manage.BamPileUp requires 1-2 arguments\n";
// usage += "   (1) full path to a directory of *.bam files to pile-up (i.e. bamDir=" + bamDir + " (no default))\n" + "";
// usage += "   OPTIONAL:";
// usage += "   (1) full path to a reference fasta  (i.e. ref= (no default))\n" + "";
// usage += "   (2) full path to a file of segments to subset the pile up  (i.e. segs= (no default))\n" + "";
//
// for (int i = 0; i < args.length; i++) {
// if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
// System.err.println(usage);
// System.exit(1);
// } else if (args[i].startsWith("bamDir=")) {
// bamDir = ext.parseStringArg(args[i], null);
// numArgs--;
// } else if (args[i].startsWith("ref=")) {
// referenceGenomeFasta = ext.parseStringArg(args[i], null);
// numArgs--;
// } else if (args[i].startsWith("segs=")) {
// segFile = ext.parseStringArg(args[i], null);
// numArgs--;
// } else if (args[i].startsWith("log=")) {
// logfile = args[i].split("=")[1];
// numArgs--;
// } else {
// System.err.println("Error - invalid argument: " + args[i]);
// }
// }
// if (numArgs != 0) {
// System.err.println(usage);
// System.exit(1);
// }
// try {
// log = new Logger(logfile);
//
// testDir(log);
// } catch (Exception e) {
// e.printStackTrace();
// }
// }
// }

// SamRecordFilter dup = new DuplicateReadFilter();
// SamRecordFilter al = new AlignedFilter(true);
// filters.add(new AlignedFilter(true));
// filters.add(new SecondaryAlignmentFilter());
// filters.add(getValidRecordFilter());
// filters.add(getValidReferenceFilter());
// filters.add(getProperlyPairedFilter());
// System.out.println(samRecord.getDuplicateReadFlag());
// System.out.println(samRecord.getMappingQuality());

// if(al.filterOut(samRecord)||dup.filterOut(samRecord)||new FilterNGS().getSamRecordMapQFilter(minMapQ).filterOut(samRecord)||new SecondaryAlignmentFilter().filterOut(samRecord)){
// //System.out.println(samRecord.getDuplicateReadFlag());
// //System.out.println(samRecord.getReadUnmappedFlag());
//
// }else{
// System.out.println(samRecord.getReadUnmappedFlag());
//
// //System.out.println(samRecord.getDuplicateReadFlag());
//
// }
// for (int i = 0; i < filters.size(); i++) {
//
// //
// // if(filter.filterOut(samRecord)){
// // System.out.println(samRecord.getDuplicateReadFlag());
// // System.out.println(samRecord.getMappingQuality());
// // System.out.println(i);
// // }
// }
// long time = System.currentTimeMillis();

// // log.reportTimeInfo(ext.getTimeElapsed(time) + "Filter time");
// time = System.currentTimeMillis();

// System.out.println(samRecord.toString() + "\t" + samRecord.getReferenceName() + "\t" + samRecord.getAlignmentStart() + "\t" + samRecord.getAlignmentEnd());
// log.reportTimeInfo(ext.getTimeElapsed(time) + "\t init time");

// try {
// Thread.sleep(100);
// } catch (InterruptedException ie) {
// }
// System.out.println(bamPiles.size() + "\tSize of pile");
// if(bamPiles.size()<80){
// System.exit(1);
// }
// if(currentSegment.getStop()<samRecord.getAlignmentEnd()){
// log.reportTimeError("Unsorted\t"+currentSegment.getStart()+"\t"+samRecord.getAlignmentStart());
// //System.out.println(bamPiles.size() + "\tSize of pile\t" + bamPiles.get(i).getBin().getUCSClocation() + "\t" + SamRecordOps.getReferenceSegmentForRecord(samRecord, log).getUCSClocation());
//
// System.exit(1);
// }
// System.out.println(samRecord.getMappingQuality());

// if (bamPiles.size() < 101) {
// System.out.println(bamPiles.size() + "\tSize of pile\t" + bamPiles.get(i).getBin().getUCSClocation() + "\t" + SamRecordOps.getReferenceSegmentForRecord(samRecord, log).getUCSClocation());
// try {
// Thread.sleep(100);
// } catch (InterruptedException ie) {
// }
// }

// int start = 0;
// int stop = r.length();
// for (int i = 1; i <= p.length; i++) {
// if (samRecord.getReferencePositionAtReadPosition(i) == toPile.getStart()) {
// start = i - 1;// 1 based - 0 base
// break;
// }
// }
// for (int i = start; i <= p.length; i++) {
// if (samRecord.getReferencePositionAtReadPosition(i) == toPile.getStop()) {
// stop = i - 1;// 1 based - 0 base
// break;
// }
// }

// String[] readBases = SamRecordOps.getReadBases(samRecord);
//

// samRecord.getAlignmentBlocks().get(1).
// samRecord.getCigar().getCigarElements().get(0).getOperator().
// samRecord.get

// TODO here we go
// samRecord.g

// log.reportTimeInfo(ext.getTimeElapsed(time2));

// boolean[] toRemove = new boolean[bamPiles.size()];

//
//
// Segment samRecordSegment = SamRecordOps.getReferenceSegmentForRecord(samRecord, log);
//
// for (int i = 0; i < bamPiles.size(); i++) {
// if (bamPiles.get(i).getBin().overlaps(samRecordSegment)) {
// bamPiles.get(i).addRecord(samRecord, filterNGS.getPhreadScoreFilter(), log);
// toRemove[i] = false;
// } else {
// toRemove[i] = true;
// bamPiles.get(i).summarize();
// }
// }
//
//
// for (int i = 0; i < toRemove.length; i++) {
// if (!toRemove[i]) {
// bamPilesTmp.add(bamPiles.get(i));
// } else {
// writer.println(bamPiles.get(i).getOuput());
// }
// }
// bamPiles = bamPilesTmp;
