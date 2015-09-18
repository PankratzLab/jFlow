package seq.cnv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import cnv.var.CNVariant;
import cnv.var.LocusSet;
import seq.manage.BEDFileReader;
import seq.manage.BEDFileReader.BEDFeatureSeg;
import seq.manage.BamPile;
import seq.manage.BamSegPileUp.PileupProducer;
import seq.manage.BedOps;
import seq.qc.Mappability;
import stats.Histogram.DynamicAveragingHistogram;
import stats.Rscript.COLUMNS_MULTIPLOT;
import stats.Rscript.PLOT_DEVICE;
import stats.Rscript.RScatter;
import stats.Rscript.RScatters;
import stats.Rscript.SCATTER_TYPE;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.WorkerTrain;
import common.ext;
import common.ArraySpecialList.ArrayIntList;
import common.WorkerTrain.Producer;
import filesys.Segment;

/**
 * Currently geared toward providing qcMetrics for ExomeDepth Calls from raw bam data
 *
 */
public class CnvBamQC {
	private static final String[] QC_HEADER = new String[] { "Ref_Map_score", "Samp_MapQ", "Population_MapQ", "Samp_Depth", "Population_Depth", "Diff_Depth", "STATE" };
	private String[] bams;
	private CallSplit callSplit;
	private BamPile[][] bamPiles;
	private Mappability<CNVariant> mappability;
	private Logger log;

	public CnvBamQC(String[] bams, CallSplit callSplit, BamPile[][] bamPiles, Mappability<CNVariant> mappability, Logger log) {
		super();
		this.bams = bams;
		this.callSplit = callSplit;
		this.bamPiles = bamPiles;
		this.mappability = mappability;
		this.log = log;
	}

	public void summarize(String output) {
		mappability.computeMappability();
		HistogramQC histogramQC = new HistogramQC();

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER) + "\t" + Array.toStr(QC_HEADER));
			LocusSet<CNVariant> cnLocusSet = callSplit.getCnLocusSet();

			for (int i = 0; i < cnLocusSet.getLoci().length; i++) {
				CNVariant currentCnVariant = cnLocusSet.getLoci()[i];
				double mapScore = mappability.getMappabilityResults().get(i).getAverageMapScore();
				writer.print(currentCnVariant.toPlinkFormat());
				QcResult popQcResult = getPopQcResult(currentCnVariant, i);
				writer.print("\t" + mapScore);
				writer.print("\t" + popQcResult.getSampMapQ());
				writer.print("\t" + popQcResult.getPopAvgMapQ());
				writer.print("\t" + popQcResult.getSampDepth());
				writer.print("\t" + popQcResult.getPopAvgDepth());
				writer.print("\t" + (popQcResult.getSampDepth() - popQcResult.getPopAvgDepth()));

				// writer.print(mappability.getMappabilityResults().get(i).getAverageMapScore());

				if (currentCnVariant.getCN() > 2) {
					histogramQC.getMapQPopDup().addDataPair(mapScore, popQcResult.getPopAvgMapQ());
					writer.print("\tDUP");
				} else if (currentCnVariant.getCN() < 2) {
					histogramQC.getMapQPopDel().addDataPair(mapScore, popQcResult.getPopAvgMapQ());
					writer.print("\tDEL");

				} else {
					log.reportTimeError("Invalid cnv " + currentCnVariant.toPlinkFormat());
					writer.close();
					return;
				}
				writer.println();

				histogramQC.getMapQPopAll().addDataPair(mapScore, popQcResult.getPopAvgMapQ());
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}
		plotResults(output, histogramQC, log);
	}

	private static RScatter[] plotResults(String summary, HistogramQC histogramQC, Logger log) {
		String[] plotFiles = Array.concatAll(new String[] { summary }, splitCn(summary, log));
		String rscatterAll = ext.parseDirectoryOfFile(summary) + "summary.qc";
		ArrayList<RScatter> rScatters = new ArrayList<RScatter>();
		ArrayList<RScatter> rScatterHists = histogramQC.plotAndDumpMapQPop(ext.addToRoot(summary, ".histogram"), log);
		rScatters.addAll(rScatterHists);
		for (int i = 0; i < plotFiles.length; i++) {

			String plotBox = ext.removeDirectoryInfo(plotFiles[i]) + "_boxCN";
			String outBox = ext.rootOf(plotFiles[i], false) + "_boxCN";
			String[] ysBox = new String[] { QC_HEADER[3], QC_HEADER[4], QC_HEADER[5] };
			RScatter rsScatterBox = new RScatter(plotFiles[i], outBox + ".rscript", plotBox, outBox + ".jpeg", QC_HEADER[6], ysBox, SCATTER_TYPE.BOX, log);
			rsScatterBox.setyLabel("X coverage");
			rsScatterBox.setOverWriteExisting(true);
			rsScatterBox.execute();
			rScatters.add(rsScatterBox);

			String plotMap = ext.removeDirectoryInfo(plotFiles[i]) + "_mapQ";
			String outMap = ext.rootOf(plotFiles[i], false) + "_mapQ";
			String[] ysMap = new String[] { QC_HEADER[1], QC_HEADER[2] };
			RScatter rsScatterMap = new RScatter(plotFiles[i], outMap + ".rscript", plotMap, outMap + ".jpeg", QC_HEADER[0], ysMap, SCATTER_TYPE.POINT, log);
			rsScatterMap.setOverWriteExisting(true);
			rsScatterMap.setxLabel(QC_HEADER[0]);
			rsScatterMap.setyLabel("Phred-scaled MapQ");
			rsScatterMap.setTitle(ext.rootOf(plotFiles[i]) + " Mapping Quality Metrics");
			rsScatterMap.execute();
			rScatters.add(rsScatterMap);

			String plotDepth = ext.removeDirectoryInfo(plotFiles[i]) + "_depth";
			String outDepth = ext.rootOf(plotFiles[i], false) + "_depth";
			System.out.println(plotMap + "\t" + plotDepth);
			String[] ys = new String[] { QC_HEADER[5] };
			RScatter rsScatterDepth = new RScatter(plotFiles[i], outDepth + ".rscript", plotDepth, outDepth + ".jpeg", QC_HEADER[4], ys, SCATTER_TYPE.POINT, log);
			rsScatterDepth.setOverWriteExisting(true);
			rsScatterDepth.setxLabel(QC_HEADER[4]);
			rsScatterDepth.setyLabel(QC_HEADER[5]);
			rsScatterDepth.setTitle(ext.rootOf(plotFiles[i]) + " Depth Quality Metrics");
			rsScatterDepth.execute();
			rScatters.add(rsScatterDepth);
		}

		RScatters rsScatters = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]), rscatterAll + ".rscript", rscatterAll + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, log);

		rsScatters.execute();
		return null;

	}

	/**
	 * @param summary
	 *            Split this summary file by {@link CNVariant} deletion vs duplication
	 * @param log
	 * @return
	 */
	private static String[] splitCn(String summary, Logger log) {
		ArrayList<String> splits = new ArrayList<String>();
		splits.add(ext.addToRoot(summary, ".DUP"));
		splits.add(ext.addToRoot(summary, ".DEL"));
		PrintWriter[] writers = Files.getAppropriateWriters(Array.toStringArray(splits));
		for (int i = 0; i < writers.length; i++) {
			writers[i].println(Array.toStr(CNVariant.PLINK_CNV_HEADER) + "\t" + Array.toStr(QC_HEADER));
		}
		try {
			BufferedReader reader = Files.getAppropriateReader(summary);
			reader.readLine();
			CNVariant[] cnvs = CNVariant.toCNVariantArray(CNVariant.loadPlinkFile(summary, null, false));
			int index = 0;
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("[\\s]+");
				if (cnvs[index].getCN() > 2) {
					writers[0].println(Array.toStr(line));
				} else if (cnvs[index].getCN() < 2) {
					writers[1].println(Array.toStr(line));
				} else {
					log.reportTimeError("Invalid cnv " + cnvs[index].toPlinkFormat());
					return null;
				}
				index++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + summary + "\" not found in current directory");
			return null;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + summary + "\"");
			return null;
		}
		Files.closeAllWriters(writers);
		return Array.toStringArray(splits);

	}

	private QcResult getPopQcResult(CNVariant currentCnVariant, int index) {
		int sampleIndex = getIndex(currentCnVariant, bams, log);
		double sampMapQ = 0;
		double sampDepth = 0;
		double popAvgMapQ = 0;
		double popAvgDepth = 0;

		int[] matchedIndices = Array.toIntArray(callSplit.getMatched()[index]);
		for (int i = 0; i < matchedIndices.length; i++) {
			for (int j = 0; j < bamPiles.length; j++) {
				if (j == sampleIndex) {
					sampMapQ += bamPiles[j][matchedIndices[i]].getOverallAvgMapQ();
					sampDepth += bamPiles[j][matchedIndices[i]].getOverallAvgDepth();
				} else {
					// TODO exclude any other overlapping
					// System.out.println(bamPiles[j].length + "\t" + bamPiles.length + "\t" + j + "\t" + matchedIndices.length + "\t" + matchedIndices[i]);
					popAvgMapQ += bamPiles[j][matchedIndices[i]].getOverallAvgMapQ();
					popAvgDepth += bamPiles[j][matchedIndices[i]].getOverallAvgDepth();

				}

			}
		}
		sampMapQ /= matchedIndices.length;
		sampDepth /= matchedIndices.length;
		popAvgMapQ /= (matchedIndices.length * (bamPiles.length - 1));
		popAvgDepth /= (matchedIndices.length * (bamPiles.length - 1));

		return new QcResult(sampMapQ, sampDepth, popAvgMapQ, popAvgDepth);
	}

	/**
	 * Store histogram of QC metrics across all cnvs
	 *
	 */
	private static class HistogramQC {
		private DynamicAveragingHistogram mapQPopDel;
		private DynamicAveragingHistogram mapQPopDup;
		private DynamicAveragingHistogram mapQPopAll;

		public HistogramQC() {
			super();
			this.mapQPopDel = new DynamicAveragingHistogram(0, 1.1, 2);
			this.mapQPopDup = new DynamicAveragingHistogram(0, 1.1, 2);
			this.mapQPopAll = new DynamicAveragingHistogram(0, 1.1, 2);

		}

		public DynamicAveragingHistogram getMapQPopDel() {
			return mapQPopDel;
		}

		public DynamicAveragingHistogram getMapQPopDup() {
			return mapQPopDup;
		}

		public DynamicAveragingHistogram getMapQPopAll() {
			return mapQPopAll;
		}

		public ArrayList<RScatter> plotAndDumpMapQPop(String output, Logger log) {
			ArrayList<RScatter> rScatters = new ArrayList<RScatter>();
			mapQPopDel.average();
			mapQPopDup.average();
			mapQPopAll.average();
			DynamicAveragingHistogram[] popsMapQ = new DynamicAveragingHistogram[] { mapQPopAll, mapQPopDel, mapQPopDup };
			String[] titles = new String[] { QC_HEADER[2] + "_all_CN", QC_HEADER[2] + "_Deletions", QC_HEADER[2] + "_Duplications" };
			DynamicAveragingHistogram.dumpToSameFile(popsMapQ, titles, output, log);
			System.out.println(output);

			String rootOutAverage = ext.addToRoot(output, ".average");
			RScatter rScatterAverage = new RScatter(output, rootOutAverage + ".rscript", ext.removeDirectoryInfo(rootOutAverage), rootOutAverage + ".pdf", "Bin", Array.tagOn(titles, null, DynamicAveragingHistogram.DUMP_AVG), SCATTER_TYPE.POINT, log);
			rScatterAverage.setTitle("MapQ average");
			rScatterAverage.setxLabel(QC_HEADER[0]);
			rScatterAverage.setyLabel("Average " + QC_HEADER[2]);
			rScatterAverage.setOverWriteExisting(true);
			rScatters.add(rScatterAverage);

			String rootOutCount = ext.addToRoot(output, ".count");

			RScatter rScatterCount = new RScatter(output, rootOutCount + ".rscript", ext.removeDirectoryInfo(rootOutCount), rootOutCount + ".pdf", "Bin", Array.tagOn(titles, null, DynamicAveragingHistogram.DUMP_COUNT), SCATTER_TYPE.POINT, log);
			rScatterCount.setTitle("Count by " + QC_HEADER[0]);
			rScatterCount.setxLabel(QC_HEADER[0]);
			rScatterCount.setyLabel("Count");
			rScatterCount.setOverWriteExisting(true);
			rScatters.add(rScatterCount);
			String rootOutProp = ext.addToRoot(output, ".prop");

			RScatter rScatterProp = new RScatter(output, rootOutProp + ".rscript", ext.removeDirectoryInfo(rootOutProp), rootOutProp + ".pdf", "Bin", Array.tagOn(titles, null, DynamicAveragingHistogram.DUMP_PROP), SCATTER_TYPE.POINT, log);
			rScatterProp.setTitle("Proportion by " + QC_HEADER[0]);

			rScatterProp.setxLabel(QC_HEADER[0]);
			rScatterProp.setyLabel("Proportion");
			rScatterProp.setOverWriteExisting(true);
			rScatters.add(rScatterProp);

			return rScatters;
		}
	}

	/**
	 * Store qc metrics across a {@link CNVariant}
	 *
	 */
	private static class QcResult {
		private double sampMapQ;
		private double sampDepth;
		private double popAvgMapQ;
		private double popAvgDepth;

		public QcResult(double sampMapQ, double sampDepth, double popAvgMapQ, double popAvgDepth) {
			super();
			this.sampMapQ = sampMapQ;
			this.sampDepth = sampDepth;
			this.popAvgMapQ = popAvgMapQ;
			this.popAvgDepth = popAvgDepth;
		}

		public double getSampMapQ() {
			return sampMapQ;
		}

		public double getSampDepth() {
			return sampDepth;
		}

		public double getPopAvgMapQ() {
			return popAvgMapQ;
		}

		public double getPopAvgDepth() {
			return popAvgDepth;
		}

	}

	private static int getIndex(CNVariant cnVariant, String[] bams, Logger log) {
		int index = -1;
		for (int i = 0; i < bams.length; i++) {
			if (ext.rootOf(bams[i]).contains(cnVariant.getIndividualID())) {
				if (index > 0) {
					String error = "Multiple bam files found for " + cnVariant.getIndividualID();
					log.reportTimeError(error);
					throw new IllegalStateException(error);

				} else {
					index = i;
				}
			}
		}
		if (index < 0) {
			String error = "Could not find corresponding bam file for sample " + cnVariant.getIndividualID();
			log.reportTimeError(error);
			throw new IllegalStateException(error);
		}
		return index;
	}

	public static void qcCNVs(String bams, String cnvFile, String serDir, String mappabilityFile, String callSubsetBed, String referenceGenomeFasta, int numThreads, Logger log) {
		// ReferenceGenome referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);

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

		String serToReport = serDir + ext.rootOf(cnvFile) + "_QC/";
		PileupProducer producer = new PileupProducer(bamFiles, serToReport, null, null, callSplit.getSegsToSearch(), log);
		WorkerTrain<BamPile[]> train = new WorkerTrain<BamPile[]>(producer, numThreads, 2, log);
		BamPile[][] bamPiles = new BamPile[bamFiles.length][];
		int index = 0;
		while (train.hasNext()) {
			BamPile[] bPilesTmp = train.next();
			bamPiles[index] = bPilesTmp;
			index++;
		}
		CnvBamQC cnvBamQC = new CnvBamQC(bamFiles, callSplit, bamPiles, mappability, log);
		String summary = serToReport + "qc_summary.txt";
		cnvBamQC.summarize(summary);

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

		public LocusSet<CNVariant> getCnLocusSet() {
			return cnLocusSet;
		}

		public ArrayIntList[] getMatched() {
			return matched;
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
					tmpSplit.add(overlaps[j].getIntersection(cnLocusSet.getLoci()[i], log));
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

	public static void main(String[] args) {
		String referenceGenomeFasta = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
		String mappabilityFile = "/home/pankrat2/public/bin/ref/mapability/hg19/wgEncodeCrgMapabilityAlign100mer.bedGraph";
		String cnvFile = "/home/tsaim/shared/Project_Tsai_Project_021/ExomeDepth/results/ExomeDepth.all.cnvs";
		String callSubsetBed = "/panfs/roc/groups/5/pankrat2/public/bin/ExomeDepth/exons.hg19.sort.chr.bed";
		String bams = "/home/tsaim/shared/Project_Tsai_Project_021/bam/";
		String serDir = "/home/tsaim/shared/Project_Tsai_Project_021/ExomeDepth/QC/";
		Logger log = new Logger(ext.rootOf(cnvFile, false) + ".qc.log");
		int numThreads = 24;
		qcCNVs(bams, cnvFile, serDir, mappabilityFile, callSubsetBed, referenceGenomeFasta, numThreads, log);
	}

	// public static void maina(String[] args) {
	// int numArgs = args.length;
	// String segFile = null;
	// String referenceGenomeFasta = "hg19_canonical.fa";
	// String bams = null;
	// String pfbFile = null;
	// double minPhred = 30;
	// double minMapQ = 30;
	//
	// int minDepth = 20;
	// int minAltDepth = -1;
	//
	// int numthreads = 4;
	// String logfile = null;
	// Logger log;
	// String usage = "\n" + "seq.manage.BamPileUp requires 1-2 arguments\n";
	// usage += "   (1) full path to a directory of *.bam files or a file listing bam files in the first column, to pile-up (i.e. bams=" + bams + " (no default))\n" + "";
	// usage += "   (2) full path to a reference fasta  (i.e. ref= (no default))\n" + "";
	// usage += "   (3) full path to a file of segments to subset the pile up  (i.e. segs= (no default))\n" + "";
	// usage += "   (4) minimum phred score for a base pair to be piled  (i.e. minPhred=" + minPhred + " (default))\n" + "";
	// usage += "   (5) minimum mapping quality score for a read to be piled  (i.e. minMapQ=" + minMapQ + " (default))\n" + "";
	// usage += "   (6) minimum total depth for a position to be reported  (i.e. minDepth=" + minDepth + " (default))\n" + "";
	// usage += "   (7) minimum alternate allele depth for a position to be reported  (i.e. minAltDepth=" + minAltDepth + " (default, no minimum))\n" + "";
	// usage += "   (8) file listing pfbs for known variants  (i.e. pfb= (no default))\n" + "";
	//
	// usage += PSF.Ext.getNumThreadsCommand(9, numthreads);
	//
	// for (int i = 0; i < args.length; i++) {
	// if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
	// System.err.println(usage);
	// System.exit(1);
	// } else if (args[i].startsWith("bams=")) {
	// bams = ext.parseStringArg(args[i], null);
	// numArgs--;
	// } else if (args[i].startsWith("ref=")) {
	// referenceGenomeFasta = ext.parseStringArg(args[i], null);
	// numArgs--;
	// } else if (args[i].startsWith("segs=")) {
	// segFile = ext.parseStringArg(args[i], null);
	// numArgs--;
	// } else if (args[i].startsWith("pfb=")) {
	// pfbFile = ext.parseStringArg(args[i], null);
	// numArgs--;
	// } else if (args[i].startsWith("minPhred=")) {
	// minPhred = ext.parseDoubleArg(args[i]);
	// numArgs--;
	// } else if (args[i].startsWith("minDepth=")) {
	// minDepth = ext.parseIntArg(args[i]);
	// numArgs--;
	// } else if (args[i].startsWith("minAltDepth=")) {
	// minAltDepth = ext.parseIntArg(args[i]);
	// numArgs--;
	// } else if (args[i].startsWith("minMapQ=")) {
	// minMapQ = ext.parseDoubleArg(args[i]);
	// numArgs--;
	// } else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
	// numthreads = ext.parseIntArg(args[i]);
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
	// logfile = Files.exists(bams) ? bams + "contam.log" : logfile;
	// log = new Logger(logfile);
	// FilterNGS filterNGS = new FilterNGS(minMapQ, minPhred, new int[] { minDepth, minAltDepth });
	// } catch (Exception e) {
	// e.printStackTrace();
	// }
	// }
}
