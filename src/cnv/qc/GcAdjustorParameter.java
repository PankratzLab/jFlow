package cnv.qc;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import common.Files;
import common.Logger;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;
import stats.CrossValidation;
import cnv.analysis.CentroidCompute;
import cnv.analysis.CentroidCompute.Builder;
import cnv.filesys.Centroids;
import cnv.filesys.MarkerSet.PreparedMarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import cnv.qc.GcAdjustor.GcModel;

/**
 * Stores betas to re-create a gc-correction for a particular sample
 *
 */
public class GcAdjustorParameter implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String sample;
	private double[] betas;
	private GC_CORRECTION_METHOD correctionMethod;

	public GcAdjustorParameter(String sample, double[] betas, GC_CORRECTION_METHOD correctionMethod) {
		super();
		this.betas = betas;
		this.sample = sample;
		this.correctionMethod = correctionMethod;
		if (betas == null || betas.length != 2) {
			throw new IllegalArgumentException("betas must be length 2");
		}
	}

	public GC_CORRECTION_METHOD getCorrectionMethod() {
		return correctionMethod;
	}

	private double adjust(GC_CORRECTION_METHOD method, String sampleName, double intensity, double gc) {
		double[] X = new double[betas.length];
		double predicted = 0;
		X[0] = 1;// constant term
		X[1] = gc;
		predicted = betas[0];
		predicted += betas[1] * X[1];
		return intensity - predicted;
	}

	/**
	 * A bare bones adjustment
	 */
	public double[] quickAdjust(GC_CORRECTION_METHOD method, String sampleName, double[] intensities, double[] gcs) {
		double[] ret = new double[intensities.length];
		verify(method, sampleName, null);
		for (int i = 0; i < ret.length; i++) {
			if (!Double.isNaN(intensities[i]) && !Double.isNaN(gcs[i])) {
				ret[i] = adjust(method, sampleName, intensities[i], gcs[i]);
			} else {
				ret[i] = Double.NaN;
			}
		}
		return ret;
	}

	/**
	 * @param val
	 *            intensity values, like LRR
	 * @param gc
	 *            gc content from {@link GcModel}<br>
	 * 
	 * 
	 *            For compatability with {@link GcAdjustor} methods
	 */
	public CrossValidation adjust(GC_CORRECTION_METHOD method, String sampleName, double[] intensity, double[][] gc, boolean verbose, Logger log) {
		verify(method, sampleName, gc);
		CrossValidation cv = new CrossValidation(null, null, intensity, gc, verbose, false, log);
		cv.setBetas(betas);
		cv.computePredictedValues();
		cv.computeResiduals();
		return cv;
	}

	private void verify(GC_CORRECTION_METHOD method, String sampleName, double[][] gc) {
		if (!sample.equals(sampleName)) {
			throw new IllegalArgumentException("Mismatched samples, object is storing " + sample + " and " + sampleName + " was provided");
		}
		if (gc != null && betas.length != gc[0].length + 1) {
			throw new IllegalArgumentException("Independant variable sizes do not match beta length");
		}
		if (method != correctionMethod) {
			throw new IllegalArgumentException("Mismatched correction methods, obect is storing " + correctionMethod + " and " + method + " was provided");
		}
	}

	private static class GcAdjustorWorker implements Callable<GcAdjustorParameter[][]> {
		private Project proj;
		private PreparedMarkerSet markerSet;
		private String sample;
		private GcModel gcmodel;
		private Centroids[] centroids;
		private GC_CORRECTION_METHOD[] correction_METHODs;
		private boolean debugMode;

		public GcAdjustorWorker(Project proj, PreparedMarkerSet markerSet, String sample, GcModel gcmodel, Centroids[] centroids, GC_CORRECTION_METHOD[] correction_METHODs, boolean debugMode) {
			super();
			this.proj = proj;
			this.sample = sample;
			this.gcmodel = gcmodel;
			this.centroids = centroids;
			this.correction_METHODs = correction_METHODs;
			this.markerSet = markerSet;
			this.debugMode = debugMode;
		}

		@Override
		public GcAdjustorParameter[][] call() throws Exception {
			Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
			GcAdjustorParameter[][] params = new GcAdjustorParameter[centroids == null ? 1 : 1 + centroids.length][correction_METHODs.length];
			float[] intensites = samp.getLRRs();
			ArrayList<GcAdjustorParameter> parameters = getParams(intensites);
			params[0] = parameters.toArray(new GcAdjustorParameter[parameters.size()]);
			if (centroids != null) {
				for (int i = 0; i < centroids.length; i++) {
					float[] centIntensites = samp.getLRRs(centroids[i].getCentroids());
					ArrayList<GcAdjustorParameter> centParameters = getParams(centIntensites);
					params[i + 1] = centParameters.toArray(new GcAdjustorParameter[centParameters.size()]);
				}
			}
			return params;
		}

		private ArrayList<GcAdjustorParameter> getParams(float[] intensites) {
			ArrayList<GcAdjustorParameter> parameters = new ArrayList<GcAdjustorParameter>();
			for (int i = 0; i < correction_METHODs.length; i++) {
				//long time = System.currentTimeMillis();
				GcAdjustor gcAdjustor = GcAdjustor.getComputedAdjustor(proj, markerSet, intensites, gcmodel, correction_METHODs[i], false, false, false);
				if (debugMode) {
					//proj.getLog().reportTimeInfo(correction_METHODs[i] + "\toriginal way " + ext.getTimeElapsed(time));
				}
				GcAdjustorParameter gcAdjustorParameters = new GcAdjustorParameter(sample, gcAdjustor.getCrossValidation().getBetas(), correction_METHODs[i]);
				parameters.add(gcAdjustorParameters);
			}
			return parameters;
		}
	}

	private static class GcAdjustorProducer implements Producer<GcAdjustorParameter[][]> {
		private Project proj;
		private String[] samples;
		private GcModel gcmodel;
		private boolean debugMode;
		private Centroids[] centroids;
		private GC_CORRECTION_METHOD[] correction_METHODs;
		private PreparedMarkerSet markerSet;
		private int index;

		public GcAdjustorProducer(Project proj, String[] samples, GcModel gcmodel, boolean debugMode, Centroids[] centroids, GC_CORRECTION_METHOD[] correction_METHODs) {
			super();
			this.proj = proj;
			this.markerSet = new PreparedMarkerSet(proj.getMarkerSet());
			this.samples = samples;
			this.gcmodel = gcmodel;
			this.debugMode = debugMode;
			this.centroids = centroids;
			this.correction_METHODs = correction_METHODs;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<GcAdjustorParameter[][]> next() {
			GcAdjustorWorker worker = new GcAdjustorWorker(proj, markerSet, samples[index], gcmodel, centroids, correction_METHODs, debugMode);
			index++;
			return worker;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}
	}

	public static GcAdjustorParameter[][][] generateAdjustmentParameters(Project proj, GC_CORRECTION_METHOD[] methods, int numThreads) {
		String customCent = proj.DATA_DIRECTORY.getValue() + "custom_recomp.cent";
		proj.CUSTOM_CENTROIDS_FILENAME.setValue(customCent);
		if (!Files.exists(customCent)) {
			CentroidCompute.computeAndDumpCentroids(proj, null, customCent, new Builder(), numThreads, 2);
		}
		return generateAdjustmentParameters(proj, new String[] { customCent }, methods, numThreads);
	}

	public static GcAdjustorParameter[][][] generateAdjustmentParameters(Project proj, String[] centroidFiles, GC_CORRECTION_METHOD[] methods, int numThreads) {
		Logger log = proj.getLog();
		if (!Files.exists(proj.GC_MODEL_FILENAME.getValue())) {
			log.reportTimeError(proj.GC_MODEL_FILENAME.getValue() + " did not exist, cannot generate gc correction parameters");
			return null;
		} else {
			Centroids[] centroids = null;
			if (centroidFiles != null) {
				centroids = new Centroids[centroidFiles.length];
				for (int i = 0; i < centroids.length; i++) {
					centroids[i] = Centroids.load(centroidFiles[i], false);
				}
			}
			GcModel gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(false, false), false, log);
			GcAdjustorProducer producer = new GcAdjustorProducer(proj, proj.getSamples(), gcModel, true, centroids, methods);
			WorkerTrain<GcAdjustorParameter[][]> train = new WorkerTrain<GcAdjustorParameter[][]>(producer, numThreads, 2, proj.getLog());
			int index = 0;

			GcAdjustorParameter[][][] finalParams = new GcAdjustorParameter[centroids == null ? 1 : 1 + centroids.length][methods.length][proj.getSamples().length];
			while (train.hasNext()) {
				GcAdjustorParameter[][] tmp = train.next();// cents,methods
				for (int i = 0; i < tmp.length; i++) {
					for (int j = 0; j < tmp[i].length; j++) {
						finalParams[i][j][index] = tmp[i][j];
					}
				}
				index++;
				if (index % 100 == 0) {
					proj.getLog().reportTimeInfo("Generated gc correction parameters for " + index + " of " + proj.getSamples().length + " samples");
				}
			}
			for (int i = 0; i < finalParams.length; i++) {
				for (int j = 0; j < finalParams[i].length; j++) {
					String output = proj.DATA_DIRECTORY.getValue() + "gc_parameters." + methods[j] + ".ser";
					if (i != 0) {
						output = ext.addToRoot(output, "." + ext.rootOf(centroidFiles[i - 1]));
						GcAdjustorParameters tmp = new GcAdjustorParameters(finalParams[i][j], centroids[i - 1], methods[j], proj.getSampleList().getFingerprint(), proj.getMarkerSet().getFingerprint());
						tmp.writeSerial(output);
					} else {
						GcAdjustorParameters tmp = new GcAdjustorParameters(finalParams[i][j], null, methods[j], proj.getSampleList().getFingerprint(), proj.getMarkerSet().getFingerprint());
						tmp.writeSerial(output);
					}
					proj.GC_CORRECTION_PARAMETERS_FILENAMES.addValue(output);
				}
			}
			proj.saveProperties();
			return finalParams;
		}
	}

	/**
	 * Storage of adjustment parameters, also stores any the centroids that were used
	 *
	 */
	public static class GcAdjustorParameters implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private GcAdjustorParameter[] gcAdjustorParameters;
		private GC_CORRECTION_METHOD correction_METHOD;
		private Centroids centroids;
		private long sampleFingerprint;
		private long markerFingerprint;

		public GcAdjustorParameters(GcAdjustorParameter[] gcAdjustorParameters, Centroids centroids, GC_CORRECTION_METHOD correction_METHOD, long sampleFingerprint, long markerFingerprint) {
			super();
			this.gcAdjustorParameters = gcAdjustorParameters;
			this.centroids = centroids;
			this.sampleFingerprint = sampleFingerprint;
			this.markerFingerprint = markerFingerprint;
			this.correction_METHOD = correction_METHOD;
			if (centroids != null && centroids.getFingerprint() != markerFingerprint) {
				throw new IllegalArgumentException("Mismatched fingerprints for centroids");
			}
			for (int i = 0; i < gcAdjustorParameters.length; i++) {
				if (gcAdjustorParameters[i].getCorrectionMethod() != correction_METHOD) {
					throw new IllegalArgumentException("Mismatched correction methods ");
				}
			}
		}

		public void writeSerial(String filename) {
			Files.writeSerial(this, filename, true);
		}

		public static GcAdjustorParameters readSerial(String filename, Logger log) {
			return (GcAdjustorParameters) Files.readSerial(filename, false, log, false, true);
		}

	}

	public static void main(String[] args) {
		Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);
		generateAdjustmentParameters(proj, GC_CORRECTION_METHOD.values(), 5);
	}
	// double[] tmp = gcAdjustor.getCorrectedIntensities();
	//
	// double[][] tmpGC = new double[gcAdjustor.getFullGcs().length][1];
	// for (int j = 0; j < tmpGC.length; j++) {
	// tmpGC[j][0] = (double) gcAdjustor.getFullGcs()[j];
	// }
	//
	// time = System.currentTimeMillis();
	// GcAdjustor gcAdjustor2 = GcAdjustor.getComputedAdjustor(proj, sample, gcAdjustorParameters, markerSet, intensites, gcmodel, correction_METHODs[i], true, true, false);
	// proj.getLog().reportTimeInfo(correction_METHODs[i] + "\tnew way " + ext.getTimeElapsed(time));
	//
	// time = System.currentTimeMillis();
	// double[] gcs = gcmodel.getGcsFor(markerSet.getMarkerNames());
	// double[] tmp3 = gcAdjustorParameters.quickAdjust(correction_METHODs[i], sample, Array.toDoubleArray(intensites), gcs);
	// proj.getLog().reportTimeInfo(correction_METHODs[i] + "\tfast track way " + ext.getTimeElapsed(time));
	//
	// double[] tmp2 = gcAdjustor2.getCorrectedIntensities();
	// for (int j = 0; j < tmp2.length; j++) {
	// if (!Double.isNaN(tmp3[j]) && !Double.isNaN(tmp[j]) && !Double.isNaN(tmp2[j]) && tmp[j] != tmp2[j] && tmp[j] != tmp3[j]) {
	// System.out.println(tmp[j] + "\t" + tmp2[j]);
	// System.exit(1);
	// }
	// }
}
