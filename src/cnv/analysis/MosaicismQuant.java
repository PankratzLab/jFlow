package cnv.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.Callable;

import org.apache.commons.math3.optim.nonlinear.scalar.gradient.NonLinearConjugateGradientOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.NonLinearConjugateGradientOptimizer.Formula;

import mining.Calcfc;
import mining.Cobyla;
import mining.CobylaExitStatus;
import stats.Rscript.RScatter;
import stats.Rscript.SCATTER_TYPE;
import common.Array;
import common.Sort;
import common.WorkerHive;
import common.ext;
import cnv.analysis.MosaicismQuant.ComputeParams.Builder;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.hmm.CNVCaller;
import cnv.qc.LrrSd;
import cnv.qc.SampleQC;
import filesys.Segment;

/**
 * @author lane0212 Class to quantify mosiacism at a particular locus, based on procedure from http://www.ncbi.nlm.nih.gov/pubmed/22277120
 */
public class MosaicismQuant implements Calcfc  {
	private static final double MIN_BAF = 0.15;
	private static final double MAX_BAF = 0.85;
	private Project proj;
	private SampleMosiac sampleMosiac;
	private Segment evalSegment;
	private MarkerSet markerSet;
	private int[][] indicesByChr;
	private MOSAIC_TYPE type;
	private ComputeParams computeParams;

	public MosaicismQuant(Project proj, SampleMosiac sampleMosiac, MOSAIC_TYPE type, ComputeParams computeParams, Segment evalSegment, MarkerSet markerSet, int[][] indicesByChr) {
		super();
		this.proj = proj;
		this.sampleMosiac = new SampleMosiac(sampleMosiac);
		this.evalSegment = evalSegment;
		this.markerSet = markerSet;
		this.indicesByChr = indicesByChr;
		this.type = type;
		this.computeParams = computeParams;
	}

	private enum MOSAIC_TYPE {
		MONOSOMY_DISOMY, TRISOMY_DISOMY;
	}

	@Override
	public double Compute(int n, int m, double[] x, double[] con) {
		System.out.println(Array.toStr(con));
		if (x[0] - computeParams.getMinF() < 0) {
			con[0] = -100000;
		}
		if (computeParams.getMaxF() - x[0] < 0) {
			con[1] = -100000;
		}
		if (x[1] - computeParams.getMinShift() < 0) {
			con[2] = -100000;
		}
		if (computeParams.getMaxShift() - x[1] < 0) {
			con[3] = -100000;
		}
		if (x[2] - computeParams.getMinNullD() < 0) {
			con[4] = 0;
		}
		if (computeParams.getMaxNullD() - x[2] < 0) {
			con[5] = 0;

		}
		return getValue(x);
	}

	private double getValue(double[] x) {
		double f = x[0];
		double shift = x[1];
		// shift = -96.3918991185974;
		// f = 0.229963195;
		double delta = getDelta(f, type);

		double nullD = x[2];
		double nullDHalf = getHalfNull(nullD);
		int count = sampleMosiac.getCdf().getVals().length;
		int halfCount = getHalfCount(count);
		double halfSigma = getHalfSigma(halfCount, shift);
		System.out.println(halfCount);
		System.out.println(halfSigma);
		System.out.println(count);
		System.out.println(shift);
		System.out.println(shift);
		double[] currentCase = sampleMosiac.getCdf().getValsInOrder();
		double[] currentControl = sampleMosiac.getSmoothedControl().getValsInOrder();
		double diffAverage = Array.mean(currentCase) - Array.mean(currentControl);
		System.out.println(Array.mean(currentControl) + "AV\t" + Array.stdev(currentControl));

		double[] shiftControl = getShift(currentControl, diffAverage, delta, halfSigma);
		System.out.println(Array.mean(shiftControl) + "AV\t" + Array.stdev(shiftControl));

		int[] shiftControlRank = Sort.quicksort(shiftControl);
		double[] resids = getResid(shiftControl, shiftControlRank, currentCase);
		System.out.println(diffAverage);
		System.out.println(getSumResids(resids));
		System.out.println(Array.toStr(x));
		// System.exit(1);
		return getSumResids(resids);
	}

	public void prepareMosiacQuant(int numThreads, double minBaf, double maxBaf) {
		sampleMosiac.load(numThreads);
		sampleMosiac.developCDFs(markerSet, indicesByChr, evalSegment, minBaf, maxBaf);

	}

	private static double getSumResids(double[] resids) {
		double sum = 0;
		for (int i = 0; i < resids.length; i++) {
			sum += Math.sqrt(Math.pow(resids[i], 2));
		}
		return sum;

	}

	private static double[] getResid(double[] shiftControl, int[] shiftControlRank, double[] currentCase) {
		double[] resid = new double[shiftControl.length];
		for (int i = 0; i < resid.length; i++) {
			resid[i] = currentCase[i] - shiftControl[shiftControlRank[i]];
		}
		return resid;
	}

	private static double[] getShift(double[] currentControl, double diffAverage, double delta, double halfSigma) {
		double[] shift = new double[currentControl.length];
		for (int i = 0; i < currentControl.length; i++) {
			double val = diffAverage + currentControl[i];
			if (i - halfSigma > 0) {
				val += delta * -1;
			} else {
				val += delta;
			}
			shift[i] = val;
		}
		return shift;
	}

	private static double getHalfDelta(double f) {
		return 1 * (0.5 - (1 / (2 - f)));
	}

	private static double getHalfSigma(int halfCount, double shift) {
		return (double) halfCount - shift;
	}

	private static int getHalfCount(int count) {
		return (int) (count + 0.5) / 2;
	}

	private static int getHalfNull(double nullD) {
		return (int) (nullD + 0.5) / 2;
	}

	private static double getDelta(double f, MOSAIC_TYPE type) {
		double delta = Double.NaN;
		switch (type) {
		case MONOSOMY_DISOMY:
			delta = -1 * (0.5 - (1 / (2 - f)));
			break;
		case TRISOMY_DISOMY:
			delta = 0.5 - (1 / (2 + f));
			break;
		default:
			break;

		}
		return delta;
	}

	public SampleMosiac getSampleMosiac() {
		return sampleMosiac;
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

		public CDF getSmoothedControl() {
			return smoothedControl;
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
			WorkerHive<SampleMosiacBase> hive = new WorkerHive<MosaicismQuant.SampleMosiacBase>(numThreads, 10, getProj().getLog());
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

		public double[] getValsInOrder() {
			double[] orderedVals = new double[vals.length];
			for (int i = 0; i < orderedVals.length; i++) {
				orderedVals[i] = vals[order[i]];
			}
			return orderedVals;
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

		SampleMosiac sampleMosiac = prep(proj, "7355066051_R03C01", "BAF1585_SD", 5);
		MarkerSet markerSet = proj.getMarkerSet();
		int[][] indices = markerSet.getIndicesByChr();

		Segment segTest = new Segment("chr17:42,963,198-78,940,173");
		Builder builder = new Builder();
		ComputeParams params = builder.build();
		MosaicismQuant mosiacismQuant = new MosaicismQuant(proj, sampleMosiac, MOSAIC_TYPE.MONOSOMY_DISOMY, params, segTest, markerSet, indices);
		mosiacismQuant.prepareMosiacQuant(5, MIN_BAF, MAX_BAF);

		String testDir = proj.PROJECT_DIRECTORY.getValue() + "TestMosaic/";
		new File(testDir).mkdirs();
		String out = testDir + ext.replaceWithLinuxSafeCharacters(segTest.getUCSClocation(), true);
		mosiacismQuant.getSampleMosiac().plotCDFs(out);

		System.out.println(mosiacismQuant.Compute(1, 1, params.getX(), params.getCon()));
		
		CobylaExitStatus cobylaExitStatus = Cobyla.FindMinimum(mosiacismQuant, params.getX().length, params.getCon().length, params.getX(), 100, .00000001, 3, 5000);
		System.out.println(cobylaExitStatus);
		System.exit(1);

	}

	public static class ComputeParams {

		private double f ;
		private double minF ;
		private double maxF ;
		private double shift;
		private double minShift ;
		private double maxShift ;
		private double nullD ;
		private double minNullD ;
		private double maxNullD;

		private double[] getX() {
			return new double[] { f == 0 ? f + .1 : f, shift, nullD };
		}

		private double[] getCon() {
			return new double[] { minF, maxF, minShift, maxShift, minNullD, maxNullD };
		}

		public double getMinF() {
			return minF;
		}

		public double getMaxF() {
			return maxF;
		}

		public double getMinShift() {
			return minShift;
		}

		public double getMaxShift() {
			return maxShift;
		}

		public double getMinNullD() {
			return minNullD;
		}

		public double getMaxNullD() {
			return maxNullD;
		}

		public static class Builder {
			private double f = 0;
			private double minF = 0;
			private double maxF = 1;
			private double shift = 0;
			private double minShift = 0;
			private double maxShift =1;
			private double nullD = 0;
			private double minNullD = 0;
			private double maxNullD = 100;

			public Builder f(double f) {
				this.f = f;
				return this;
			}

			public Builder minF(double minF) {
				this.minF = minF;
				return this;
			}

			public Builder maxF(double maxF) {
				this.maxF = maxF;
				return this;
			}

			public Builder shift(double shift) {
				this.shift = shift;
				return this;
			}

			public Builder minShift(double minShift) {
				this.minShift = minShift;
				return this;
			}

			public Builder maxShift(double maxShift) {
				this.maxShift = maxShift;
				return this;
			}

			public Builder nullD(double nullD) {
				this.nullD = nullD;
				return this;
			}

			public Builder minNullD(double minNullD) {
				this.minNullD = minNullD;
				return this;
			}

			public Builder maxNullD(double maxNullD) {
				this.maxNullD = maxNullD;
				return this;
			}

			public ComputeParams build() {
				return new ComputeParams(this);
			}
		}

		private ComputeParams(Builder builder) {
			this.f = builder.f;
			this.minF = builder.minF;
			this.maxF = builder.maxF;
			this.shift = builder.shift;
			this.minShift = builder.minShift;
			this.maxShift = builder.maxShift;
			this.nullD = builder.nullD;
			this.minNullD = builder.minNullD;
			this.maxNullD = builder.maxNullD;
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
