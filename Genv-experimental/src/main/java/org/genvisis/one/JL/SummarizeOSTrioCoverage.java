package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.concurrent.Callable;

import javax.jms.IllegalStateException;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

public class SummarizeOSTrioCoverage {
	private static final int[] covTargets = new int[] { 5, 10, 20, 30, 40 };

	// /home/spectorl/lanej0/OS_seq/bamQC
	public static void summarize(String vpopFile, String indir, String outDir, int numThreads,
			Logger log) throws IllegalStateException {
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY,
				log);
		vpop.report();
		String[] files = Files.listFullPaths(indir,
				".sorted.dedup.realigned.recal.txt", false);
		log.reportTimeInfo("Found " + files.length + " files to summarize");
		Hashtable<String, String> map = new Hashtable<String, String>();
		for (int i = 0; i < files.length; i++) {
			String key = ext.removeDirectoryInfo(files[i]).split("_")[0];
			map.put(key, files[i]);
		}
		String outputFinal = outDir + "trioCoverageSummary.txt";
		// if (!Files.exists(outputFinal)) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outputFinal));

			ArrayList<String> famsNotFound = new ArrayList<String>();
			SumFamProducer producer = new SumFamProducer(vpop, outDir, log, map);
			WorkerTrain<FamSum> train = new WorkerTrain<SummarizeOSTrioCoverage.FamSum>(producer, numThreads, 2, log);

			int num = 0;
			ArrayList<FamSum> haveData = new ArrayList<FamSum>();
			while (train.hasNext()) {

				FamSum famSum = train.next();
				if (num == 0) {
					writer.println(Array.toStr(famSum.getFamCoverageResults().developHeader()));
				}
				num++;
				if (famSum.getFamCoverageResults().getNumTotalRegions() > 0) {
					writer.println(Array.toStr(famSum.getFamCoverageResults().getData()));
					haveData.add(famSum);
				}

			}
			Files.writeArray(Array.toStringArray(famsNotFound), outDir + "famsNotFound.txt");
			writer.close();

		} catch (Exception e) {
			log.reportError("Error writing to " + outputFinal);
			log.reportException(e);
		}
		// }

		// String scripA = outputFinal + ".avgC.rscript";
		// String plotA = ext.removeDirectoryInfo(scripA);
		// String outA = outputFinal + ".avgC.jpeg";
		// RScatter rscatterAvg = new RScatter(outputFinal, scripA, plotA, outA, "MO_AVG_COVERAGE", new String[] { "FA_AVG_COVERAGE" }, "OFF_AVG_COVERAGE", SCATTER_TYPE.POINT, log);
		// rscatterAvg.setOverWriteExisting(true);
		// rscatterAvg.execute();

		// String scripAIN = outputFinal + ".avgCI.rscript";
		// String plotAIn = ext.removeDirectoryInfo(scripAIN);
		// String outAIN = outputFinal + ".avgCI.jpeg";
		// RScatter rscatterAvgIN = new RScatter(outputFinal, scripAIN, plotAIn, outAIN, "TRIO_AVG_COVERAGE", new String[] { "INTERSECT_ProportionCoveredAt10X", "INTERSECT_ProportionCoveredAt20X", "INTERSECT_ProportionCoveredAt30X" }, null, SCATTER_TYPE.POINT, log);
		// rscatterAvgIN.setrSafeAltYColumnNames(new String[] { "10X", "20X", "30X" });
		// rscatterAvgIN.setOverWriteExisting(true);
		// rscatterAvgIN.setxLabel("Average coverage of trio");
		// rscatterAvgIN.setyLabel("Proportion of targets covered in all members of trio");
		// rscatterAvgIN.setyRange(new double[] { 0, 1 });
		// rscatterAvgIN.setFontsize(10);
		// rscatterAvgIN.execute();

		String scripAINMIN = outputFinal + ".avgCIMin.rscript";
		String plotAInMIN = ext.removeDirectoryInfo(scripAINMIN);
		String outAINMIN = outputFinal + ".avgCIMin.jpeg";
		RScatter rscatterAvgINMIN = new RScatter(outputFinal, scripAINMIN, plotAInMIN, outAINMIN, "MIN_AVG_COVERAGE", new String[] { "INTERSECT_ProportionCoveredAt5X", "INTERSECT_ProportionCoveredAt10X", "INTERSECT_ProportionCoveredAt20X", "INTERSECT_ProportionCoveredAt30X", "INTERSECT_ProportionCoveredAt40X" }, null, SCATTER_TYPE.POINT, log);
		rscatterAvgINMIN.setrSafeAltYColumnNames(new String[] { "5X Coverage", "10X Coverage", "20X Coverage", "30X Coverage", "40X Coverage" });
		rscatterAvgINMIN.setOverWriteExisting(true);
		rscatterAvgINMIN.setxLabel("Minimum average coverage of trio (n=" + Files.countLines(outputFinal, 1) + ")");
		rscatterAvgINMIN.setyLabel("Proportion of targets covered in all members of trio");
		rscatterAvgINMIN.setyRange(new double[] { 0, 1 });
		rscatterAvgINMIN.setxRange(new double[] { 15, 55 });

		rscatterAvgINMIN.setFontsize(10);
		rscatterAvgINMIN.execute();

		String scripBOX = outputFinal + ".avgCIMinbox.rscript";
		String plotABOX = ext.removeDirectoryInfo(scripBOX);
		String outAINBOX = outputFinal + ".avgCIMinbox.jpeg";
		RScatter rscatterAvgBOX = new RScatter(outputFinal, scripBOX, plotABOX, outAINBOX, "MIN_AVG_COVERAGE", new String[] { "INTERSECT_ProportionCoveredAt5X", "INTERSECT_ProportionCoveredAt10X", "INTERSECT_ProportionCoveredAt20X", "INTERSECT_ProportionCoveredAt30X", "INTERSECT_ProportionCoveredAt40X" }, null, SCATTER_TYPE.BOX_NO_MELT, log);
		rscatterAvgBOX.setrSafeAltYColumnNames(new String[] { "5X Coverage", "10X Coverage", "20X Coverage", "30X Coverage", "40X Coverage" });
		rscatterAvgBOX.setOverWriteExisting(true);
		rscatterAvgBOX.setxLabel("Minimum average coverage of trio (n=" + Files.countLines(outputFinal, 1) + ")");
		rscatterAvgBOX.setyLabel("Proportion of targets covered in all members of trio");
		rscatterAvgBOX.setyRange(new double[] { 0, 1 });
		rscatterAvgBOX.setFontsize(10);
		rscatterAvgBOX.execute();

		// RScatters rsScatters = new RScatters(rsArrayList.toArray(new RScatter[rsArrayList.size()]), outputFinal + ".rscript", outputFinal + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, log);
		// rsScatters.execute();

	}

	private static class SumFamProducer extends AbstractProducer<FamSum> {

		private VcfPopulation vpop;
		private Logger log;
		private Hashtable<String, String> map;
		private Iterator<String> iter;
		private String outDir;

		public SumFamProducer(VcfPopulation vpop, String outDir, Logger log, Hashtable<String, String> map) {
			super();
			this.vpop = vpop;
			this.log = log;
			this.outDir = outDir;
			this.map = map;
			this.iter = vpop.getSubPop().keySet().iterator();
		}

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return iter.hasNext();
		}

		@Override
		public Callable<FamSum> next() {
			String fam = iter.next();
			Set<String> famInds = vpop.getSubPop().get(fam);
			FamSum famSum = new FamSum(outDir, log, map, fam, famInds);
			return famSum;
		}
	}

	private static class FamSum implements Callable<FamSum> {
		private String outDir;
		private Logger log;
		private Hashtable<String, String> map;
		private String fam;
		private Set<String> famInds;
		private FamCoverageResults famCoverageResults;

		public FamSum(String outDir, Logger log, Hashtable<String, String> map, String fam, Set<String> famInds) {
			super();
			this.outDir = outDir;
			this.log = log;
			this.map = map;
			this.fam = fam;
			this.famInds = famInds;
		}

		public FamCoverageResults getFamCoverageResults() {
			return famCoverageResults;
		}

		// public void setFamCoverageResults(FamCoverageResults famCoverageResults) {
		// this.famCoverageResults = famCoverageResults;
		// }

		@Override
		public FamSum call() throws Exception {
			try {
				this.famCoverageResults = sumFam(outDir, log, map, fam, famInds);
			} catch (IllegalStateException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return this;
		}

	}

	private static class FamCoverageResults {
		// private String covFile;
		private String fam;
		private String off;
		private String mo;
		private String fa;
		private int numTotalRegions;
		private int[] numCoveredAtOff;
		private int[] numCoveredAtMo;
		private int[] numCoveredAtFa;
		private int[] numCoveredAtIntersect;
		// private int[] coverageTargets;
		private double[] totalAvgCoverage;

		public FamCoverageResults(String covFile, String fam, String off, String mo, String fa, int numTotalRegions, int[] numCoveredAtOff, int[] numCoveredAtMo, int[] numCoveredAtFa, int[] numCoveredAtIntersect, int[] coverageTargets, double[] totalAvgCoverage) {
			super();
			this.fam = fam;
			this.off = off;
			this.mo = mo;
			this.fa = fa;
			this.numTotalRegions = numTotalRegions;
			this.numCoveredAtOff = numCoveredAtOff;
			this.numCoveredAtMo = numCoveredAtMo;
			this.numCoveredAtFa = numCoveredAtFa;
			this.numCoveredAtIntersect = numCoveredAtIntersect;
			// this.coverageTargets = coverageTargets;
			this.totalAvgCoverage = totalAvgCoverage;
		}

		public int getNumTotalRegions() {
			return numTotalRegions;
		}

		private String[] developHeader() {
			ArrayList<String> h = new ArrayList<String>();
			h.add("Family");
			h.add("OFF");
			h.add("MO");
			h.add("FA");
			h.add("TotalTargets");
			h.add("OFF_AVG_COVERAGE");
			h.add("MO_AVG_COVERAGE");
			h.add("FA_AVG_COVERAGE");
			h.add("TRIO_AVG_COVERAGE");
			h.add("MIN_AVG_COVERAGE");
			h.add("MAX_AVG_COVERAGE");

			for (int i = 0; i < covTargets.length; i++) {
				h.add("OFF_ProportionCoveredAt" + covTargets[i] + "X");
				h.add("MO_ProportionCoveredAt" + covTargets[i] + "X");
				h.add("FA_ProportionCoveredAt" + covTargets[i] + "X");
				h.add("INTERSECT_ProportionCoveredAt" + covTargets[i] + "X");
			}
			return Array.toStringArray(h);
		}

		private String[] getData() {
			ArrayList<String> h = new ArrayList<String>();
			h.add(fam);
			h.add(off);
			h.add(mo);
			h.add(fa);
			h.add(numTotalRegions + "");
			h.add(totalAvgCoverage[0] + "");
			h.add(totalAvgCoverage[1] + "");
			h.add(totalAvgCoverage[2] + "");
			h.add(Array.mean(totalAvgCoverage) + "");
			h.add(Array.min(totalAvgCoverage) + "");
			h.add(Array.max(totalAvgCoverage) + "");

			for (int i = 0; i < covTargets.length; i++) {
				h.add((numTotalRegions == 0 ? 0 : (double) numCoveredAtOff[i] / numTotalRegions) + "");
				h.add((numTotalRegions == 0 ? 0 : (double) numCoveredAtMo[i] / numTotalRegions) + "");
				h.add((numTotalRegions == 0 ? 0 : (double) numCoveredAtFa[i] / numTotalRegions) + "");
				h.add((numTotalRegions == 0 ? 0 : (double) numCoveredAtIntersect[i] / numTotalRegions) + "");
			}
			return Array.toStringArray(h);
		}

	}

	private static FamCoverageResults sumFam(String outDir, Logger log, Hashtable<String, String> map, String fam, Set<String> famInds) throws IllegalStateException {
		String off = null;
		String offFile = null;
		String mo = null;
		String moFile = null;
		String fa = null;
		String faFile = null;
		boolean useFam = false;
		int numTotalRegions = 0;
		int[] coverageTargets = covTargets;
		String coverageTrio = outDir + fam + ".covSummary.txt";

		int[] numCoveredAtOff = new int[coverageTargets.length];
		int[] numCoveredAtMo = new int[coverageTargets.length];
		int[] numCoveredAtFa = new int[coverageTargets.length];
		int[] numCoveredAtIntersect = new int[coverageTargets.length];
		double[] avgCoverage = new double[3];
		for (String ind : famInds) {
			ind = ind.replaceAll(".variant.*", "");
			if (ind.endsWith("C")) {
				if (off != null) {
					throw new IllegalStateException("Too many kids");
				}
				off = ind;
				if (map.containsKey(off)) {
					offFile = map.get(off);
					useFam = true;
				} else {
					log.reportTimeInfo("No file for ind " + ind);
				}
			}

			else if (ind.endsWith("D")) {
				if (fa != null) {
					throw new IllegalStateException("Too many pops");
				}
				fa = ind;
				if (map.containsKey(fa)) {
					faFile = map.get(fa);
				}
			} else if (ind.endsWith("M")) {
				if (mo != null) {
					throw new IllegalStateException("Too many mom");
				}
				mo = ind;
				if (map.containsKey(mo)) {
					moFile = map.get(mo);
				}
			} else {
				throw new IllegalStateException("Unknown");

			}
		}
		if (useFam) {
			log.reportTimeInfo("Paste for " + fam);
			String[] tags = new String[] { off, mo, fa };
			String avgC = "averageCoverage";
			int col = ext.indexOfStr(avgC,
					Files.getHeaderOfFile(offFile, log));
			if (!Files.exists(coverageTrio)) {
				Files.paste(new String[] { offFile, moFile, faFile },
						coverageTrio, new int[] { 0, col }, 0, tags, null, log);
			}

			try {
				BufferedReader reader = Files.getAppropriateReader(coverageTrio);
				String[] header = new String[] { avgC + "_" + off, avgC + "_" + mo, avgC + "_" + fa };
				int[] indices = ext.indexFactors(header, Files.getHeaderOfFile(coverageTrio, log), true, true);
				reader.readLine();

				while (reader.ready()) {
					numTotalRegions++;
					String[] line = reader.readLine().trim().split("\t");
					double offC = Double.parseDouble(line[indices[0]]);
					double moC = Double.parseDouble(line[indices[1]]);
					double faC = Double.parseDouble(line[indices[2]]);
					avgCoverage[0] += offC;
					avgCoverage[1] += moC;
					avgCoverage[2] += faC;

					for (int i = 0; i < coverageTargets.length; i++) {
						if (offC >= coverageTargets[i]) {
							numCoveredAtOff[i]++;
						}
						if (moC >= coverageTargets[i]) {
							numCoveredAtMo[i]++;

						}
						if (faC >= coverageTargets[i]) {
							numCoveredAtFa[i]++;

						}
						if (faC >= coverageTargets[i] && moC >= coverageTargets[i] && offC >= coverageTargets[i]) {
							numCoveredAtIntersect[i]++;
						}
					}
				}
				reader.close();
				log.reportTimeInfo("Found " + numTotalRegions + " total regions");
				for (int i = 0; i < avgCoverage.length; i++) {
					avgCoverage[i] = (double) avgCoverage[i] / numTotalRegions;
				}
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + coverageTrio + "\" not found in current directory");
				return null;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + coverageTrio + "\"");
				return null;
			}

		}
		return new FamCoverageResults(coverageTrio, fam, off, mo, fa, numTotalRegions, numCoveredAtOff, numCoveredAtMo, numCoveredAtFa, numCoveredAtIntersect, coverageTargets, avgCoverage);

	}

	public static void main(String[] args) {
		String vpop = "/panfs/roc/groups/12/spectorl/lanej0/OS_seq/trioCoverage/dn.vpop";
		String inDir = "/panfs/roc/groups/12/spectorl/lanej0/OS_seq/bamQC/";
		String outDir = "/panfs/roc/groups/12/spectorl/lanej0/OS_seq/trioCoverage/";
		Logger log = new Logger(outDir + "tc.log");
		int numThreads = 24;
		try {
			summarize(vpop, inDir, outDir, numThreads, log);
		} catch (IllegalStateException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
