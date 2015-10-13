package cnv.analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;

import stats.Rscript.RScatter;
import stats.Rscript.SCATTER_TYPE;
import common.Array;
import common.Sort;
import common.WorkerHive;
import common.ext;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.hmm.CNVCaller;
import cnv.hmm.PFB;
import cnv.manage.ExtProjectDataParser;
import cnv.manage.ExtProjectDataParser.Builder;
import cnv.qc.LrrSd;
import cnv.qc.SampleQC;

/**
 * @author lane0212 Class to quantify mosiacism at a particular locus, based on procedure from http://www.ncbi.nlm.nih.gov/pubmed/22277120
 */
public class MosiacismQuant {
	private static final double MIN_BAF = 0.15;
	private static final double MAX_BAF = 0.85;

	private static SampleMosiac prep(Project proj, String sampleName, String qcMetric, int numControls) {
		SampleMosiacBase[] controls = new SampleMosiacBase[numControls];
		int sampleIndex = ext.indexOfStr(sampleName, proj.getSamples());
		SampleQC sampleQC = SampleQC.loadSampleQC(proj, LrrSd.SAMPLE_COLUMN, new String[] { qcMetric }, true);
		double[] data = sampleQC.getDataFor(qcMetric);
		int[] sorted = Sort.quicksort(data);
		int sampleSortIndex = ext.indexOfStr(sampleIndex + "", Array.toStringArray(sorted));
		int minIndex = Math.max(0, sampleSortIndex - numControls);
		int maxIndex = Math.min(proj.getSamples().length, sampleSortIndex + numControls);
		ArrayList<Double> controlData = new ArrayList<Double>();
		ArrayList<String> controlDataNames = new ArrayList<String>();

		double sampData = data[sorted[sampleSortIndex]];
		for (int i = minIndex; i < maxIndex; i++) {
			if (i != sampleSortIndex) {
				controlData.add(data[sorted[i]]);
				controlDataNames.add(proj.getSamples()[sorted[i]]);
			}
		}
		double[] distanceToSample = Array.InplaceAbs(Array.InplaceSub(Array.toDoubleArray(controlData), sampData));
		int[] minDists = Sort.quicksort(distanceToSample);
		for (int i = 0; i < controls.length; i++) {
			int index = minDists[i];
			SampleMosiacBase control = new SampleMosiacBase(proj, controlDataNames.get(index), qcMetric, controlData.get(index), distanceToSample[index]);
			controls[i] = control;
		}
		return new SampleMosiac(proj, sampleName, qcMetric, sampData, controls);
	}

	private static class SampleMosiac extends SampleMosiacBase {
		private Project proj;
		private String sampleName;

		private SampleMosiacBase[] controls;

		public SampleMosiac(Project proj, String sampleName, String qcMetric, double qcValue, SampleMosiacBase[] controls) {
			super(proj, sampleName, qcMetric, qcValue, qcValue);

			this.controls = controls;
		}

		public SampleMosiacBase[] getControls() {
			return controls;
		}

		public void load(int numThreads) {
			this.samp = proj.getFullSampleFromRandomAccessFile(sampleName);
			WorkerHive<SampleMosiacBase> hive = new WorkerHive<MosiacismQuant.SampleMosiacBase>(numThreads, 10, proj.getLog());
			for (int i = 0; i < controls.length; i++) {
				hive.addCallable(controls[i]);
			}
			hive.execute(true);
			ArrayList<SampleMosiacBase> results = hive.getResults();
			this.controls = results.toArray(new SampleMosiacBase[results.size()]);
		}

	}

	private static class SampleMosiacBase implements Callable<SampleMosiacBase> {
		private Project proj;
		protected Sample samp;

		protected String sampleName;
		private String qcMetric;
		private double qcValue;
		private double qcDistance;

		public SampleMosiacBase(Project proj, String sampleName, String qcMetric, double qcValue, double qcDistance) {
			super();
			this.proj = proj;
			this.sampleName = sampleName;
			this.qcMetric = qcMetric;
			this.qcValue = qcValue;
			this.qcDistance = qcDistance;
		}

		public String getSampleName() {
			return sampleName;
		}

		public String getQcMetric() {
			return qcMetric;
		}

		public double getQcValue() {
			return qcValue;
		}

		public double getQcDistance() {
			return qcDistance;
		}

		private void load() {
			this.samp = proj.getFullSampleFromRandomAccessFile(sampleName);
		}

		@Override
		public SampleMosiacBase call() throws Exception {
			load();
			return this;
		}
	}

	// select case marker/indices
	// select controls
	// select control markers
	// adjust case bafs
	// adjust control bafs
	// create case cdf
	// create control cdf

	public static void test(Project proj) {
		SampleQC sampleQC = SampleQC.loadSampleQC(proj, true);

		Sample samp = proj.getFullSampleFromRandomAccessFile(proj.getSamples()[ext.indexOfStr("7355066051_R03C01", proj.getSamples())]);
		SampleMosiac sampleMosiac = prep(proj, samp.getSampleName(), "BAF1585_SD", 5);
		for (int i = 0; i < sampleMosiac.getControls().length; i++) {
			System.out.println(sampleMosiac.getSampleName() + "\t" + sampleMosiac.getQcValue() + "\t" + sampleMosiac.getControls()[i].getSampleName() + "\t" + sampleMosiac.getControls()[i].getQcValue());
		}
		sampleMosiac.load(5);

		System.exit(1);
		MarkerSet markerSet = proj.getMarkerSet();
		int[][] indices = markerSet.getIndicesByChr();
		String testDir = proj.PROJECT_DIRECTORY.getValue() + "TestMosaic/";
		new File(testDir).mkdirs();
		PFB pfb = PFB.loadPFB(proj);
		// double[] pfbCDF = Array.getValuesBetween(Array.removeNaN(pfb.getPfbs()), MIN_BAF, MAX_BAF);
		// ArrayList<Double> pArrayList = new ArrayList<Double>();
		double[] bafsPop = loadBaf1585Mean(proj);
		for (int i = 0; i < indices.length; i++) {
			// for (int j = 0; j < pfbCDF.length; j++) {
			// pArrayList.add(pfbCDF[i]);
			// }
			// ArrayList<Double> pArrayListChr = new ArrayList<Double>();
			// pArrayListChr.addAll(pArrayList);
			// Collections.shuffle(pArrayListChr);
			// ArrayList<Double> pArrayListChrSuf = new ArrayList<Double>();

			String out = testDir + "chr" + i + "out.txt";
			double[] bafChr = Array.subArray(Array.toDoubleArray(samp.getBAFs()), indices[i]);
			double[] bafsPopChr = Array.subArray(bafsPop, indices[i]);

			// List<Double> shf = pArrayListChr.subList(0, bafs.length);
			// pArrayListChrSuf.addAll(shf);
			// double[] pArrayListChrA = Array.subArray(Array.toDoubleArray(pArrayListChrSuf), 0, bafs.length);

			// int[] popIndex = Sort.quicksort(pArrayListChrA);

			double[] bafsAdjust = CNVCaller.adjustBaf(bafChr, MIN_BAF, MAX_BAF, proj.getLog());
			double[] bafsAdjustRanged = Array.getValuesBetween(bafsAdjust, MIN_BAF, MAX_BAF);

			double[] bafsPopChrAdjust = CNVCaller.adjustBaf(Array.removeNaN(bafsPopChr), MIN_BAF, MAX_BAF, proj.getLog());
			int[] popSortAdjust = Sort.quicksort(bafsPopChrAdjust);
			ArrayList<Double> pArrayList = new ArrayList<Double>();
			boolean low = true;
			int index = 0;
			while (pArrayList.size() < bafsAdjustRanged.length) {
				try {
					if (low) {
						pArrayList.add(bafsPopChrAdjust[popSortAdjust[index]]);
						low = false;
					} else {
						pArrayList.add(bafsPopChrAdjust[popSortAdjust[popSortAdjust.length - 1 - index]]);
						low = true;
						index++;
					}
				} catch (ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException) {
					index = 0;
				}
			}
			double[] pArrayListChrA = Array.toDoubleArray(pArrayList);
			int[] popSort = Sort.quicksort(pArrayListChrA);
			int[] sortIndex = Sort.quicksort(bafsAdjustRanged);
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(out));
				writer.println("Index\tCDF\tValSamp\tValPop\tResid");
				for (int j = 0; j < sortIndex.length; j++) {
					double cdf = (double) j / bafsAdjustRanged.length;
					writer.println(j + "\t" + cdf + "\t" + bafsAdjustRanged[sortIndex[j]] + "\t" + pArrayListChrA[popSort[j]] + "\t" + (bafsAdjustRanged[sortIndex[j]] - pArrayListChrA[popSort[j]]));
				}
				writer.close();
				RScatter rScatter = new RScatter(out, out + ".rscript", ext.removeDirectoryInfo(out), out + ".jpeg", "CDF", new String[] { "ValSamp", "ValPop", "Resid" }, SCATTER_TYPE.POINT, proj.getLog());
				rScatter.setOverWriteExisting(true);
				rScatter.execute();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + out);
				proj.getLog().reportException(e);
			}
		}

	}

	private static double[] loadBaf1585Mean(Project proj) {

		Builder builder = new Builder();
		builder.dataKeyColumnName("MarkerName");
		builder.numericDataTitles(new String[] { "BAF_15_85_MEAN" });
		builder.sampleBased(false);
		builder.requireAll(true);
		builder.treatAllNumeric(false);

		try {
			ExtProjectDataParser extProjectDataParser = builder.build(proj, proj.MARKER_METRICS_FILENAME.getValue());
			extProjectDataParser.determineIndicesFromTitles();
			extProjectDataParser.loadData();
			double[] baf1585s = extProjectDataParser.getNumericDataForTitle("BAF_15_85_MEAN");
			return baf1585s;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			proj.getLog().reportFileNotFound(proj.MARKER_METRICS_FILENAME.getValue());
			return null;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "CircDuCnv.dat";
		String logfile = null;

		String usage = "\n" + "cnv.hmm.CircDuCnv requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
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
			Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);
			test(proj);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
