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
import cnv.filesys.Centroids;
import cnv.filesys.MarkerSet.PreparedMarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import cnv.qc.GcAdjustor.GcModel;

public class GcAdjustorParameters implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String sample;
	private double[] betas;
	private GC_CORRECTION_METHOD correctionMethod;

	public GcAdjustorParameters(String sample, double[] betas, GC_CORRECTION_METHOD correctionMethod) {
		super();
		this.betas = betas;
		this.sample = sample;
		this.correctionMethod = correctionMethod;
		if (betas == null || betas.length != 2) {
			throw new IllegalArgumentException("betas must be length 2");
		}
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

	public static void writeSerial(GcAdjustorParameters[] gcParameters, String filename) {
		Files.writeSerial(gcParameters, filename, true);
	}

	public static GcAdjustorParameters[] readSerial(String filename, Logger log) {
		return (GcAdjustorParameters[]) Files.readSerial(filename, false, log, false, true);
	}

	private static class GcAdjustorWorker implements Callable<GcAdjustorParameters[]> {
		private Project proj;
		private PreparedMarkerSet markerSet;
		private String sample;
		private GcModel gcmodel;
		private Centroids centroids;
		private GC_CORRECTION_METHOD[] correction_METHODs;
		private boolean debugMode;

		public GcAdjustorWorker(Project proj, PreparedMarkerSet markerSet, String sample, GcModel gcmodel, Centroids centroids, GC_CORRECTION_METHOD[] correction_METHODs, boolean debugMode) {
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
		public GcAdjustorParameters[] call() throws Exception {
			Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
			float[] intensites = null;
			if (centroids != null) {
				intensites = samp.getLRRs(centroids.getCentroids());
			} else {
				intensites = samp.getLRRs();
			}
			ArrayList<GcAdjustorParameters> parameters = new ArrayList<GcAdjustorParameters>();
			for (int i = 0; i < correction_METHODs.length; i++) {
				long time = System.currentTimeMillis();

				GcAdjustor gcAdjustor = GcAdjustor.getComputedAdjustor(proj, markerSet, intensites, gcmodel, correction_METHODs[i], false, false, false);
				proj.getLog().reportTimeInfo(correction_METHODs[i] + "\toriginal way " + ext.getTimeElapsed(time));

				GcAdjustorParameters gcAdjustorParameters = new GcAdjustorParameters(sample, gcAdjustor.getCrossValidation().getBetas(), correction_METHODs[i]);
				parameters.add(gcAdjustorParameters);
				if (debugMode) {
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
			}
			return parameters.toArray(new GcAdjustorParameters[parameters.size()]);
		}
	}

	private static class GcAdjustorProducer implements Producer<GcAdjustorParameters[]> {
		private Project proj;
		private String[] samples;
		private GcModel gcmodel;
		private boolean debugMode;
		private Centroids centroids;
		private GC_CORRECTION_METHOD[] correction_METHODs;
		private PreparedMarkerSet markerSet;
		private int index;

		public GcAdjustorProducer(Project proj, String[] samples, GcModel gcmodel, boolean debugMode, Centroids centroids, GC_CORRECTION_METHOD[] correction_METHODs) {
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
		public Callable<GcAdjustorParameters[]> next() {
			GcAdjustorWorker worker = new GcAdjustorWorker(proj, markerSet, samples[index], gcmodel, centroids, correction_METHODs, debugMode);
			index++;
			return worker;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	public static void generateAdjustmentParameters(Project proj, String centroids, GC_CORRECTION_METHOD[] methods, int numThreads) {
		Logger log = proj.getLog();
		if (!Files.exists(proj.GC_MODEL_FILENAME.getValue())) {
			log.reportTimeError(proj.GC_MODEL_FILENAME.getValue() + " did not exist, cannot generate gc correction parameters");
		} else {
			Centroids cents = null;
			if (centroids != null) {
				cents = Centroids.load(centroids, false);
			}
			GcModel gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(false, false), false, log);
			GcAdjustorProducer producer = new GcAdjustorProducer(proj, proj.getSamples(), gcModel, true, cents, methods);
			WorkerTrain<GcAdjustorParameters[]> train = new WorkerTrain<GcAdjustorParameters[]>(producer, numThreads, 2, proj.getLog());
			int index = 0;
			
			GcAdjustorParameters[][] finalParams = new GcAdjustorParameters[methods.length][proj.getSamples().length];
			while (train.hasNext()) {
				GcAdjustorParameters[] tmp = train.next();
				for (int i = 0; i < tmp.length; i++) {
					finalParams[i][index] = tmp[i];
				}
				index++;
				if (index % 100 == 0) {
					proj.getLog().reportTimeInfo("Generated gc correction parameters for " + index + " of " + proj.getSamples().length + " samples");
				}
			}
			for (int i = 0; i < methods.length; i++) {
				String output = proj.DATA_DIRECTORY.getValue() + "gc_parameters." + methods[i] + ".ser";
				proj.GC_CORRECTION_PARAMETERS_FILENAMES.addValue(output);
				writeSerial(finalParams[i], output);
			}
		}
	}

	public static void main(String[] args) {
		Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);
		generateAdjustmentParameters(proj, null, GC_CORRECTION_METHOD.values(), 5);
	}

}
