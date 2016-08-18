package org.genvisis.one.JL;

import java.util.ArrayList;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariantHash;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.CNVariant.CONSENSUS_TYPE;
import org.genvisis.filesys.CNVariant.MatchResults;
import org.genvisis.filesys.CNVariant.OVERLAP_TYPE;
import org.genvisis.filesys.LocusSet.TO_STRING_TYPE;
import org.genvisis.stats.Histogram.DynamicHistogram;
import org.genvisis.stats.Rscript.COLUMNS_MULTIPLOT;
import org.genvisis.stats.Rscript.PLOT_DEVICE;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.RScatters;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

/**
 * @author lane0212 Comparing cnvs called by PennCnv and Exome Depth in osteo
 *
 */
public class ExomeDepthEvaluation {
	private Project proj;
	private Logger log;
	private String exomeDepthCNVFile;
	private String pennCNVFile;
	private int probeCoverage;
	private boolean removeExcludeArray;

	public ExomeDepthEvaluation(Project proj, Logger log, String exomeDepthCNVFile, String pennCNVFile, int probeCoverage, boolean removeExcludeArray) {
		super();
		this.proj = proj;
		this.log = log;
		this.exomeDepthCNVFile = exomeDepthCNVFile;
		this.pennCNVFile = pennCNVFile;
		this.probeCoverage = probeCoverage;
		this.removeExcludeArray = removeExcludeArray;
	}

	public void compare() {
		MatchResults matchResults = null;
		MarkerSet markerSet = proj.getMarkerSet();
		int[][] indicesByChr = markerSet.getIndicesByChr();
		SampleData sampleData = proj.getSampleData(0, false);
		for (int i = 0; i < CNVariant.CONSENSUS_TYPE.values().length; i++) {
			CONSENSUS_TYPE cType = CNVariant.CONSENSUS_TYPE.values()[i];
			for (int j = 0; j < CNVariant.OVERLAP_TYPE.values().length; j++) {
				OVERLAP_TYPE otype = CNVariant.OVERLAP_TYPE.values()[j];
				String serFile = ext.parseDirectoryOfFile(exomeDepthCNVFile) + ext.rootOf(exomeDepthCNVFile) + "_v_" + ext.rootOf(pennCNVFile) + cType + "_" + otype + ".ser";
				if (removeExcludeArray) {
					serFile = ext.addToRoot(serFile, "._EX_");
				}
				if (!Files.exists(serFile)) {
					matchResults = CNVariant.findSignificantConsensus(exomeDepthCNVFile, pennCNVFile, ext.rootOf(serFile, false) + ".consensus.cnv", cType, otype);
					matchResults.writeSerial(serFile);
				} else {
					matchResults = MatchResults.readSerial(serFile, log);
				}
				HistogramComparison[] comparisons = new HistogramComparison[] { new HistogramComparison(), new HistogramComparison() };// with sig overlap, without sig overlap
				log.reportTimeInfo(matchResults.getMatched1().size() + " cnvs matching");
				log.reportTimeInfo(matchResults.getUnmatched1().size() + " cnvs not matching");

				for (CNVariant matchingExomeD : matchResults.getMatched1()) {
					if (!removeExcludeArray || !sampleData.individualShouldBeExcluded(matchingExomeD.getFamilyID())) {
						comparisons[0].addCNV(matchingExomeD, markerSet, indicesByChr, -1);
					}
				}
				for (CNVariant notMatchingExomeD : matchResults.getUnmatched1()) {
					if (!removeExcludeArray || !sampleData.individualShouldBeExcluded(notMatchingExomeD.getFamilyID())) {
						comparisons[1].addCNV(notMatchingExomeD, markerSet, indicesByChr, probeCoverage);
					}
				}

				DynamicHistogram[] exomeDHits = comparisons[0].getAll();

				DynamicHistogram[] exomeDMisses = comparisons[1].getAll();
				if (cType == CONSENSUS_TYPE.CN_AWARE && otype == OVERLAP_TYPE.OVERLAP_LOC_AND_INDIVIDUAL) {
					String cnvHitsFile = ext.parseDirectoryOfFile(serFile) + cType + "_" + otype + "probeCoverage_" + probeCoverage + "_hits.cnv";
					String cnvMissFile = ext.parseDirectoryOfFile(serFile) + cType + "_" + otype + "_probeCoverage_" + probeCoverage + "_miss.cnv";

					if (removeExcludeArray) {
						cnvHitsFile = ext.addToRoot(cnvHitsFile, "._EX_");
						cnvMissFile = ext.addToRoot(cnvMissFile, "._EX_");
					}

					cnvHitsFile = ext.parseDirectoryOfFile(serFile) + "ExomeDepth_Hits.cnv";
					cnvMissFile = ext.parseDirectoryOfFile(serFile) + "ExomeDepth_Miss.cnv";

					LocusSet<CNVariant> hitSet = new LocusSet<CNVariant>(comparisons[0].getCnvsToStore(Sort.DESCENDING), false, log) {

						/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					};
					hitSet.writeRegions(cnvHitsFile, TO_STRING_TYPE.REGULAR, true, log);
					LocusSet<CNVariant> missSet = new LocusSet<CNVariant>(comparisons[1].getCnvsToStore(Sort.DESCENDING), false, log) {

						/**
						 * 
						 */
						private static final long serialVersionUID = 1L;

					};
					ArrayList<String> listHit = new ArrayList<String>();
					ArrayList<String> listMiss = new ArrayList<String>();

					for (int k = 0; k < missSet.getLoci().length; k++) {
						listMiss.add(missSet.getLoci()[k].getFamilyID() + "\t" + missSet.getLoci()[k].getUCSClocation());
					}
					for (int k = 0; k < hitSet.getLoci().length; k++) {
						listHit.add(hitSet.getLoci()[k].getFamilyID() + "\t" + hitSet.getLoci()[k].getUCSClocation());
					}
					String listHitFile = ext.rootOf(serFile) + ".list.hits.txt";
					String listMissFile = ext.rootOf(serFile) + ".list.miss.txt";

					proj.INDIVIDUAL_CNV_LIST_FILENAMES.setValue(new String[] { listHitFile, listMissFile });

					Files.writeList(listHit.toArray(new String[listHit.size()]), proj.INDIVIDUAL_CNV_LIST_FILENAMES.getValue()[0]);
					Files.writeList(listMiss.toArray(new String[listMiss.size()]), proj.INDIVIDUAL_CNV_LIST_FILENAMES.getValue()[1]);

					missSet.writeRegions(cnvMissFile, TO_STRING_TYPE.REGULAR, true, log);
					proj.CNV_FILENAMES.setValue(new String[] { pennCNVFile, cnvMissFile, cnvHitsFile });
					proj.saveProperties();
					CNVariantHash.load(cnvHitsFile, 1, false, log);
					CNVariantHash.load(cnvMissFile, 1, false, log);
				}
				ArrayList<RScatter> rScatters = new ArrayList<RScatter>();
				for (int k = 0; k < exomeDMisses.length; k++) {
					plot(exomeDHits, exomeDMisses, rScatters, k, false, cType, otype);
					plot(exomeDHits, exomeDMisses, rScatters, k, true, cType, otype);

				}
				String finalOut = ext.rootOf(serFile, false) + "_probeCoverage_" + probeCoverage;
				RScatters rScatters2 = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]), finalOut + ".rscript", finalOut + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, log);
				rScatters2.execute();
			}
		}
	}

	private String plot(DynamicHistogram[] exomeDHits, DynamicHistogram[] exomeDMisses, ArrayList<RScatter> rScatters, int i, boolean prop, CNVariant.CONSENSUS_TYPE cType, CNVariant.OVERLAP_TYPE oType) {
		String output = ext.parseDirectoryOfFile(exomeDepthCNVFile) + HistogramComparison.XLABELS[i] + (prop ? ".prop" : ".count") + cType + "_" + oType + ".txt";
		DynamicHistogram.dumpToSameFile(new DynamicHistogram[] { exomeDHits[i], exomeDMisses[i] }, HistogramComparison.TITLES[i], output, prop, log);
		String root = ext.rootOf(output, false);
		RScatter rsScatterCount = new RScatter(output, root + ".rescript", ext.removeDirectoryInfo(output), root + ".pdf", "Bin", HistogramComparison.TITLES[i], SCATTER_TYPE.POINT, log);
		rsScatterCount.setxLabel(HistogramComparison.XLABELS[i]);
		rsScatterCount.setyLabel((prop ? "Proportion" : "Count"));
		rsScatterCount.setTitle("n=" + Array.sum(exomeDHits[i].getCounts()) + " match; n=" + Array.sum(exomeDMisses[i].getCounts()) + " not match\n " + cType + ";" + oType);
		// rsScatterCount.execute();
		rScatters.add(rsScatterCount);
		return root;
	}

	private static class HistogramComparison {
		private static final String[] XLABELS = new String[] { "SCORE", "SITES_exons", "log_LENGTH", "ArrayProbeCoverage" };
		private static final String[][] TITLES = new String[][] { { "ExomeDepthMatch_Score", "ExomeDepthMiss_Score" }, { "ExomeDepthMatch_Sites", "ExomeDepthMiss_Sites" }, { "ExomeDepthMatch_Length", "ExomeDepthMiss_Length" }, { "ExomeDepthMatch_ProbeCoverage", "ExomeDepthMiss_ProbeCoverage" } };
		private DynamicHistogram scoreHistogram;
		private DynamicHistogram sitesHistogram;
		private DynamicHistogram lengthHistogram;
		private DynamicHistogram probeCoverageHistogram;
		private ArrayList<CNVariant> cnvsToStore;

		public HistogramComparison() {
			this.scoreHistogram = new DynamicHistogram(0, 200, 0);
			this.sitesHistogram = new DynamicHistogram(0, 100, 0);
			this.lengthHistogram = new DynamicHistogram(0, 20, 0);
			this.probeCoverageHistogram = new DynamicHistogram(0, 100, 0);
			this.cnvsToStore = new ArrayList<CNVariant>();
		}

		public void addCNV(CNVariant cnVariant, MarkerSet markerSet, int[][] indicesByChr, int probesRequired) {
			int numprobes = markerSet.getMarkersIn(cnVariant.getBufferedSegment(50), indicesByChr).length;
			if (numprobes > probesRequired) {
				scoreHistogram.addDataPointToHistogram(cnVariant.getScore());
				sitesHistogram.addDataPointToHistogram(cnVariant.getNumMarkers());
				lengthHistogram.addDataPointToHistogram(Math.log(cnVariant.getSize()));
				probeCoverageHistogram.addDataPointToHistogram(numprobes);
				cnvsToStore.add(cnVariant);
			}
		}

		public CNVariant[] getCnvsToStore(int direction) {
			ArrayList<CNVariant> sorted = CNVariant.sortByQuality(cnvsToStore.toArray(new CNVariant[cnvsToStore.size()]), direction);
			return sorted.toArray(new CNVariant[sorted.size()]);
		}

		public DynamicHistogram[] getAll() {
			return new DynamicHistogram[] { scoreHistogram, sitesHistogram, lengthHistogram, probeCoverageHistogram };
		}

	}

	public static void compareIt(Project proj, String exomeDepthCNVFile, String pennCNVFile, int probeCoverage, boolean removeExcludeArray) {
		ExomeDepthEvaluation exomeDepthOsteo = new ExomeDepthEvaluation(proj, proj.getLog(), exomeDepthCNVFile, pennCNVFile, probeCoverage, removeExcludeArray);
		exomeDepthOsteo.compare();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String exomeDepthCNVFile = null;
		String pennCNVFile = null;
		int probeCoverage = -1;
		boolean removeExcludeArray = false;
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
			} else if (args[i].startsWith("probeCoverage=")) {
				probeCoverage = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-remove")) {
				removeExcludeArray = true;
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
			Project proj = new Project(filename, false);
			compareIt(proj, exomeDepthCNVFile, pennCNVFile, probeCoverage, removeExcludeArray);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
//
// private static class ExomeDepthComparison {
// private int totalCNVs;
// private int totalOverlap;
// private int totalSigOverlap;
// private int totalPerfectOverlap;
//
// public int getTotalCNVs() {
// return totalCNVs;
// }
//
// public void setTotalCNVs(int totalCNVs) {
// this.totalCNVs = totalCNVs;
// }
//
// public int getTotalOverlap() {
// return totalOverlap;
// }
//
// public void setTotalOverlap(int totalOverlap) {
// this.totalOverlap = totalOverlap;
// }
//
// public int getTotalSigOverlap() {
// return totalSigOverlap;
// }
//
// public void setTotalSigOverlap(int totalSigOverlap) {
// this.totalSigOverlap = totalSigOverlap;
// }
//
// public int getTotalPerfectOverlap() {
// return totalPerfectOverlap;
// }
//
// public void setTotalPerfectOverlap(int totalPerfectOverlap) {
// this.totalPerfectOverlap = totalPerfectOverlap;
// }
//
// }

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