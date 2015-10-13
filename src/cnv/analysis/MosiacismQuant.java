package cnv.analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.Callable;

import mining.Calcfc;
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
import filesys.Segment;

/**
 * @author lane0212 Class to quantify mosiacism at a particular locus, based on procedure from http://www.ncbi.nlm.nih.gov/pubmed/22277120
 */
public class MosiacismQuant implements Calcfc{
	private static final double MIN_BAF = 0.15;
	private static final double MAX_BAF = 0.85;
	private Project proj;
	private SampleMosiac sampleMosiac;
	private Segment evalSegment;
	private MarkerSet markerSet;
	private int[] indicesByChr;
	
	
	
	public MosiacismQuant(Project proj, SampleMosiac sampleMosiac, Segment evalSegment, MarkerSet markerSet, int[] indicesByChr) {
		super();
		this.proj = proj;
		this.sampleMosiac =  new SampleMosiac(sampleMosiac);
		this.evalSegment = evalSegment;
		this.markerSet = markerSet;
		this.indicesByChr = indicesByChr;
	}

	public void prepareMosiacQuant(int numThreads, double minBaf, double maxBaf) {
		sampleMosiac.load(numThreads);

	}
	
	

	@Override
	public double Compute(int n, int m, double[] x, double[] con) {
		
		
		// TODO Auto-generated method stub
		return 0;
	}
	
	
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
			SampleMosiacBase control = new SampleMosiacBase(proj, controlDataNames.get(index), qcMetric, controlData.get(index), distanceToSample[index], true);
			controls[i] = control;
		}
		return new SampleMosiac(proj, sampleName, qcMetric, sampData, controls);
	}

	private static class SampleMosiac extends SampleMosiacBase {

		private SampleMosiacBase[] controls;
		private CDF smoothedControl;

		public SampleMosiac(SampleMosiac sampleMosiac) {
			super(sampleMosiac);
			this.controls = sampleMosiac.controls;
			this.smoothedControl = sampleMosiac.smoothedControl;
		}
		
		public SampleMosiac(Project proj, String sampleName, String qcMetric, double qcValue, SampleMosiacBase[] controls) {
			super(proj, sampleName, qcMetric, qcValue, qcValue, false);
			this.controls = controls;
		}

		public SampleMosiacBase[] getControls() {
			return controls;
		}

		public void load(int numThreads) {
			load();
			WorkerHive<SampleMosiacBase> hive = new WorkerHive<MosiacismQuant.SampleMosiacBase>(numThreads, 10, getProj().getLog());
			for (int i = 0; i < controls.length; i++) {
				hive.addCallable(controls[i]);
			}
			hive.execute(true);
			ArrayList<SampleMosiacBase> results = hive.getResults();
			this.controls = results.toArray(new SampleMosiacBase[results.size()]);
			for (int i = 0; i < controls.length; i++) {
			}
		}

		public void developCDFs(MarkerSet markerSet, int[][] indicesByChr, Segment seg, double minBaf, double maxBaf) {

			developCDF(markerSet, indicesByChr, seg, minBaf, maxBaf, -1);

			double[] smoothed = new double[getCdf().getVals().length];
			for (int i = 0; i < controls.length; i++) {
				controls[i].developCDF(markerSet, indicesByChr, seg, minBaf, maxBaf, getCdf().getVals().length);
				if (i > 0) {// smoothed random
					CDF tmp = controls[i].getCdf();
					for (int j = 0; j < smoothed.length; j++) {
						smoothed[j] += tmp.getVals()[tmp.getOrder()[j]];
					}
				}
			}
			CDF tmpRandomControl = controls[0].getCdf();
			for (int i = 0; i < smoothed.length; i++) {
				double tmp = smoothed[i] + tmpRandomControl.getVals()[tmpRandomControl.getOrder()[i]];
				smoothed[i] = tmp / controls.length;
			}
			smoothed = CNVCaller.adjustBaf(smoothed, minBaf, maxBaf, proj.getLog());
			this.smoothedControl = new CDF(smoothed);
		}

		private void plotCDFs(String output) {
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println("CDF\tSample\tSmoothedControl");
				for (int i = 0; i < smoothedControl.getVals().length; i++) {
					double cdf = (double) i / smoothedControl.getVals().length;
					writer.println(cdf + "\t" + getCdf().getVals()[getCdf().getOrder()[i]] + "\t" + smoothedControl.getVals()[smoothedControl.getOrder()[i]]);
				}
				writer.close();
				RScatter rScatter = new RScatter(output, output + ".rscript", ext.removeDirectoryInfo(output), output + ".jpeg", "CDF", new String[] { "Sample", "SmoothedControl" }, SCATTER_TYPE.POINT, proj.getLog());
				rScatter.setOverWriteExisting(true);
				rScatter.execute();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + output);
				proj.getLog().reportException(e);
			}
		}
	}

	private static class SampleMosiacBase implements Callable<SampleMosiacBase> {
		protected Project proj;
		protected Sample samp;
		protected CDF cdf;
		protected String sampleName;
		private String qcMetric;
		private double qcValue;
		private double qcDistance;
		private boolean control;

		public SampleMosiacBase(SampleMosiacBase s) {
			this.proj = s.proj;
			this.samp = s.samp;
			this.cdf = s.cdf;
			this.sampleName = s.sampleName;
			this.qcMetric = s.qcMetric;
			this.qcValue = s.qcValue;
			this.qcDistance = s.qcDistance;
			this.control = s.control;
		}
		
		public SampleMosiacBase(Project proj, String sampleName, String qcMetric, double qcValue, double qcDistance, boolean control) {
			super();
			this.proj = proj;
			this.sampleName = sampleName;
			this.qcMetric = qcMetric;
			this.qcValue = qcValue;
			this.qcDistance = qcDistance;
			this.control = control;
		}

		protected String getSampleName() {
			return sampleName;
		}

		protected Project getProj() {
			return proj;
		}

		protected CDF getCdf() {
			return cdf;
		}

		protected String getQcMetric() {
			return qcMetric;
		}

		protected double getQcValue() {
			return qcValue;
		}

		protected double getQcDistance() {
			return qcDistance;
		}

		protected void load() {
			this.samp = proj.getFullSampleFromRandomAccessFile(sampleName);
		}

		protected void developCDF(MarkerSet markerSet, int[][] indicesByChr, Segment seg, double minBaf, double maxBaf, int numControlForce) {
			double[] bafs = selectBafs(markerSet, indicesByChr, seg, minBaf, maxBaf, numControlForce);
			bafs = CNVCaller.adjustBaf(bafs, minBaf, maxBaf, proj.getLog());
			this.cdf = new CDF(bafs);
		}

		private double[] selectBafs(MarkerSet markerSet, int[][] indicesByChr, Segment seg, double minBaf, double maxBaf, int numControlForce) {
			String[] markersInSeg = markerSet.getMarkersIn(seg, indicesByChr);
			int[] indicesInSeg = ext.indexLargeFactors(markersInSeg, markerSet.getMarkerNames(), true, proj.getLog(), true, false);
			double[] bafs = Array.toDoubleArray(samp.getBAFs());
			ArrayList<Integer> bafIndicesToUse = new ArrayList<Integer>();

			for (int i = 0; i < indicesInSeg.length; i++) {
				int index = indicesInSeg[i];
				if (useBaf(minBaf, maxBaf, bafs[index])) {
					bafIndicesToUse.add(index);
				}
			}
			if (!control) {
				return Array.subArray(bafs, Array.toIntArray(bafIndicesToUse));
			} else {
				if (bafIndicesToUse.size() > numControlForce) {
					Collections.shuffle(bafIndicesToUse);
					ArrayList<Integer> bafIndicesToUseTmp = new ArrayList<Integer>();
					for (int i = 0; i < numControlForce; i++) {
						bafIndicesToUseTmp.add(bafIndicesToUse.get(i));
					}
					return Array.subArray(bafs, Array.toIntArray(bafIndicesToUseTmp));
				} else {
					boolean forward = true;
					boolean canGoForward = true;
					boolean canGoBackward = true;
					int nextIndexForward = bafIndicesToUse.size() - 1;
					int nextIndexBackward = bafIndicesToUse.get(0);

					while (bafIndicesToUse.size() < numControlForce && (canGoBackward || canGoForward)) {
						if (forward) {
							if (nextIndexForward < markerSet.getMarkerNames().length && useBaf(minBaf, maxBaf, bafs[nextIndexForward])) {
								bafIndicesToUse.add(nextIndexForward);
							} else {
								canGoForward = false;
							}
							nextIndexForward++;
							forward = false;
						} else {
							if (nextIndexBackward >= 0 && useBaf(minBaf, maxBaf, bafs[nextIndexBackward])) {
								bafIndicesToUse.add(nextIndexBackward);
							} else {
								canGoBackward = false;
							}
							nextIndexBackward--;
							forward = true;
						}
					}
					if (bafIndicesToUse.size() != numControlForce) {
						proj.getLog().reportTimeWarning("Not enough baf values between " + minBaf + " and " + maxBaf + " for sample " + sampleName + ", imputting with random re-sampling");
						Collections.shuffle(bafIndicesToUse);
						ArrayList<Integer> bafIndicesToUseTmp = new ArrayList<Integer>();
						bafIndicesToUseTmp.addAll(bafIndicesToUse);
						int add = 0;
						while (bafIndicesToUseTmp.size() < numControlForce) {
							bafIndicesToUseTmp.add(bafIndicesToUse.get(add));
							add++;
						}
						return Array.subArray(bafs, Array.toIntArray(bafIndicesToUseTmp));
					} else {
						return Array.subArray(bafs, Array.toIntArray(bafIndicesToUse));
					}

				}
			}
		}

		private boolean useBaf(double minBaf, double maxBaf, double baf) {
			return Double.isFinite(baf) && baf > minBaf && baf < maxBaf;
		}

		@Override
		public SampleMosiacBase call() throws Exception {
			load();
			return this;
		}
	}

	private static class CDF {
		private double[] vals;
		private int[] order;

		public CDF(double[] vals) {
			super();
			this.vals = vals;
			this.order = Sort.quicksort(vals);

		}

		public double[] getVals() {
			return vals;
		}

		public int[] getOrder() {
			return order;
		}

	}

	// select case marker/indices check
	// select controls check
	// select control markers check
	// adjust case bafs check
	// adjust control bafs check
	// create case cdf check
	// create control cdf check
	// Non-linear fitting

	public static void test(Project proj) {

		Sample samp = proj.getFullSampleFromRandomAccessFile(proj.getSamples()[ext.indexOfStr("7355066051_R03C01", proj.getSamples())]);
		SampleMosiac sampleMosiac = prep(proj, samp.getSampleName(), "BAF1585_SD", 5);
		for (int i = 0; i < sampleMosiac.getControls().length; i++) {
			System.out.println(sampleMosiac.getSampleName() + "\t" + sampleMosiac.getQcValue() + "\t" + sampleMosiac.getControls()[i].getSampleName() + "\t" + sampleMosiac.getControls()[i].getQcValue());
		}
		sampleMosiac.load(5);
		MarkerSet markerSet = proj.getMarkerSet();
		int[][] indices = markerSet.getIndicesByChr();

		Segment segTest = new Segment((byte) 17, 0, Integer.MAX_VALUE);
		sampleMosiac.developCDFs(markerSet, indices, segTest, MIN_BAF, MAX_BAF);
		String testDir = proj.PROJECT_DIRECTORY.getValue() + "TestMosaic/";
		new File(testDir).mkdirs();
		String out = testDir + ext.replaceWithLinuxSafeCharacters(segTest.getUCSClocation(), true);
		sampleMosiac.plotCDFs(out);
	
		System.exit(1);

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
