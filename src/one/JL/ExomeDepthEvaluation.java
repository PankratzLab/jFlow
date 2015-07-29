package one.JL;

import java.util.ArrayList;

import stats.Histogram.DynamicHistogram;
import stats.Rscript.COLUMNS_MULTIPLOT;
import stats.Rscript.PLOT_DEVICE;
import stats.Rscript.RScatter;
import stats.Rscript.RScatters;
import stats.Rscript.SCATTER_TYPE;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.var.CNVariant;
import cnv.var.CNVariant.CONSENSUS_TYPE;
import cnv.var.CNVariant.MatchResults;

/**
 * @author lane0212 Comparing cnvs called by PennCnv and Exome Depth in osteo
 *
 */
public class ExomeDepthEvaluation {
	private Project proj;
	private Logger log;
	private String exomeDepthCNVFile;
	private String pennCNVFile;

	public ExomeDepthEvaluation(Project proj, Logger log, String exomeDepthCNVFile, String pennCNVFile) {
		super();
		this.proj = proj;
		this.log = log;
		this.exomeDepthCNVFile = exomeDepthCNVFile;
		this.pennCNVFile = pennCNVFile;
	}

	public void compare() {
		String serFile = ext.parseDirectoryOfFile(exomeDepthCNVFile) + ext.rootOf(exomeDepthCNVFile) + "_v_" + ext.rootOf(pennCNVFile) + ".ser";
		MatchResults matchResults = null;
		MarkerSet markerSet = proj.getMarkerSet();
		if (!Files.exists(serFile)) {
			matchResults = CNVariant.findSignificantConsensus(exomeDepthCNVFile, pennCNVFile, CONSENSUS_TYPE.CN_AWARE);
			matchResults.writeSerial(serFile);
		} else {
			matchResults = MatchResults.readSerial(serFile, log);
		}
		HistogramComparison[] comparisons = new HistogramComparison[] { new HistogramComparison(), new HistogramComparison() };// with sig overlap, without sig overlap
		log.reportTimeInfo(matchResults.getMatched1().size() + " cnvs matching");
		log.reportTimeInfo(matchResults.getUnmatched1().size() + " cnvs not matching");

		for (CNVariant matchingExomeD : matchResults.getMatched1()) {
			comparisons[0].addCNV(matchingExomeD);
		}
		for (CNVariant notMatchingExomeD : matchResults.getUnmatched1()) {

			comparisons[1].addCNV(notMatchingExomeD);
		}

		DynamicHistogram[] exomeDHits = comparisons[0].getAll();
		DynamicHistogram[] exomeDMisses = comparisons[1].getAll();
		ArrayList<RScatter> rScatters = new ArrayList<RScatter>();
		for (int i = 0; i < exomeDMisses.length; i++) {
			String rootCount = plot(exomeDHits, exomeDMisses, rScatters, i, false);
			String rootprop = plot(exomeDHits, exomeDMisses, rScatters, i, true);

		}
		String finalOut = ext.rootOf(serFile, false);
		RScatters rScatters2 = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]), finalOut + ".rscript", finalOut + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, log);
		rScatters2.execute();

		// DynamicHistogram.dumpToSameFile(histograms, titles, output, proportion, log)

	}

	private String plot(DynamicHistogram[] exomeDHits, DynamicHistogram[] exomeDMisses, ArrayList<RScatter> rScatters, int i, boolean prop) {
		String output = ext.parseDirectoryOfFile(exomeDepthCNVFile) + HistogramComparison.XLABELS[i] + (prop ? ".prop" : ".count") + ".txt";
		DynamicHistogram.dumpToSameFile(new DynamicHistogram[] { exomeDHits[i], exomeDMisses[i] }, HistogramComparison.TITLES[i], output, prop, log);
		String root = ext.rootOf(output, false);
		RScatter rsScatterCount = new RScatter(output, root + ".rescript", ext.removeDirectoryInfo(output), root + ".pdf", "Bin", HistogramComparison.TITLES[i], SCATTER_TYPE.POINT, log);
		rsScatterCount.setxLabel(HistogramComparison.XLABELS[i]);
		rsScatterCount.setyLabel((prop ? "Proportion" : "Count"));
		rsScatterCount.setTitle("n=" + Array.sum(exomeDHits[i].getCounts()) + " matching pennCNV; n=" + Array.sum(exomeDMisses[i].getCounts()) + " not matching");
		rsScatterCount.execute();
		rScatters.add(rsScatterCount);
		return root;
	}

	private static class HistogramComparison {
		private static final String[] XLABELS = new String[] { "SCORE", "SITES_exons", "log_LENGTH" };
		private static final String[][] TITLES = new String[][] { { "ExomeDepthMatch_Score", "ExomeDepthMiss_Score" }, { "ExomeDepthMatch_Sites", "ExomeDepthMiss_Sites" }, { "ExomeDepthMatch_Length", "ExomeDepthMiss_Length" } };
		private DynamicHistogram scoreHistogram;
		private DynamicHistogram sitesHistogram;
		private DynamicHistogram lengthHistogram;

		public HistogramComparison() {
			this.scoreHistogram = new DynamicHistogram(0, 200, 0);
			this.sitesHistogram = new DynamicHistogram(0, 100, 0);
			this.lengthHistogram = new DynamicHistogram(0, 20, 0);
		}

		public void addCNV(CNVariant cnVariant) {
			scoreHistogram.addDataPointToHistogram(cnVariant.getScore());
			sitesHistogram.addDataPointToHistogram(cnVariant.getNumMarkers());
			lengthHistogram.addDataPointToHistogram(Math.log(cnVariant.getSize()));
		}

		public DynamicHistogram[] getAll() {
			return new DynamicHistogram[] { scoreHistogram, sitesHistogram, lengthHistogram };
		}

	}

	private static class ExomeDepthComparison {
		private int totalCNVs;
		private int totalOverlap;
		private int totalSigOverlap;
		private int totalPerfectOverlap;

		public int getTotalCNVs() {
			return totalCNVs;
		}

		public void setTotalCNVs(int totalCNVs) {
			this.totalCNVs = totalCNVs;
		}

		public int getTotalOverlap() {
			return totalOverlap;
		}

		public void setTotalOverlap(int totalOverlap) {
			this.totalOverlap = totalOverlap;
		}

		public int getTotalSigOverlap() {
			return totalSigOverlap;
		}

		public void setTotalSigOverlap(int totalSigOverlap) {
			this.totalSigOverlap = totalSigOverlap;
		}

		public int getTotalPerfectOverlap() {
			return totalPerfectOverlap;
		}

		public void setTotalPerfectOverlap(int totalPerfectOverlap) {
			this.totalPerfectOverlap = totalPerfectOverlap;
		}

	}

	public static void compareIt(Project proj, String exomeDepthCNVFile, String pennCNVFile) {
		ExomeDepthEvaluation exomeDepthOsteo = new ExomeDepthEvaluation(proj, proj.getLog(), exomeDepthCNVFile, pennCNVFile);
		exomeDepthOsteo.compare();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String exomeDepthCNVFile = null;
		String pennCNVFile = null;
		String logfile = null;
		Logger log;

		String usage = "\n" + "one.JL.ExomeDepthOsteo requires 0-1 arguments\n";
		usage += "   (1) filename (i.e. proj=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("exomeD=")) {
				exomeDepthCNVFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("penn=")) {
				pennCNVFile = args[i].split("=")[1];
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
			log = new Logger(logfile);
			Project proj = new Project(filename, false);
			compareIt(proj, exomeDepthCNVFile, pennCNVFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
//
// CNVariantHash exomeDCnVariantHash = CNVariantHash.load(exomeDepthCNVFile, 1, false);
// CNVariantHash pennCnVariantHash = CNVariantHash.load(pennCNVFile, 1, false);
// for (String cnvInd : exomeDCnVariantHash.getHashes().keySet()) {
// Hashtable<String, CNVariant[]> indExomeDepth = exomeDCnVariantHash.getDataFor(cnvInd);
// if (pennCnVariantHash.getDataFor(cnvInd).size() > 0 && indExomeDepth.size() > 0) {
// Hashtable<String, CNVariant[]> indPenn = pennCnVariantHash.getDataFor(cnvInd);
// for (String key : indExomeDepth.keySet()) {
// // CNVariant.findSignificantConsensus(file1, file2);
// }
//
// } else {
// log.reportTimeWarning("No penncnv data for" + cnvInd);
//	