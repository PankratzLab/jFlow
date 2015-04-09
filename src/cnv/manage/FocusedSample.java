package cnv.manage;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import cnv.filesys.Project;
import cnv.filesys.Sample;
import common.Array;
import common.Files;
import common.Logger;
import common.WorkerHive;
import common.ext;

/**
 *Subsets samples
 * 
 */
public class FocusedSample {

	private Sample focusedSample;
	private boolean overwriteExisting;
	private Hashtable<String, Float> outliers;

	// private boolean fail;

	public FocusedSample(int[] focusedIndices, Sample sample, long newFingerPrint, boolean overwriteExisting) {
		super();
		this.overwriteExisting = overwriteExisting;
		this.focusedSample = getFocusedSample(sample, focusedIndices, newFingerPrint);
		this.outliers = new Hashtable<String, Float>();
		// this.fail =false;
	}

	public boolean saveToRandomAccessFile(String filename, String sampleName, Logger log) {
		boolean saved = true;
		if (!overwriteExisting && Files.exists(filename)) {
			log.reportTimeWarning(filename + " already exists, will not overwrite");
			saved = false;
		} else {
			focusedSample.saveToRandomAccessFile(filename, outliers, focusedSample.getSampleName());
		}
		return saved;
	}

	public Hashtable<String, Float> getOutliers() {
		return outliers;
	}

	private Sample getFocusedSample(Sample sample, int[] focusedIndices, long newFingerPrint) {
		byte[] abGenotypes = sample.getAB_Genotypes() == null ? null : Array.subArray(sample.getAB_Genotypes(), focusedIndices);
		byte[] forwardGenotypes = sample.getForwardGenotypes() == null ? null : Array.subArray(sample.getForwardGenotypes(), focusedIndices);

		float[] xs = sample.getXs() == null ? null : Array.subArray(sample.getXs(), focusedIndices);
		float[] ys = sample.getYs() == null ? null : Array.subArray(sample.getYs(), focusedIndices);

		float[] lrrs = sample.getLRRs() == null ? null : Array.subArray(sample.getLRRs(), focusedIndices);
		float[] bafs = sample.getBAFs() == null ? null : Array.subArray(sample.getBAFs(), focusedIndices);

		float[] gcs = sample.getGCs() == null ? null : Array.subArray(sample.getGCs(), focusedIndices);
		boolean canXYBeNegative = sample.getCanXYBeNegative();
		Sample focusedSample = new Sample(sample.getSampleName(), newFingerPrint, gcs, xs, ys, bafs, lrrs, forwardGenotypes, abGenotypes, canXYBeNegative);
		return focusedSample;
	}

	private static WorkerSubset[] getWorkers(Project original, Project newFocus, int[] focusedIndices, boolean[] samplesToUse, long newFingerPrint, boolean overwriteExisting) {
		String[] samples = samplesToUse == null ? original.getSamples() : Array.subArray(original.getSamples(), samplesToUse);
		WorkerSubset[] workerSubsets = new WorkerSubset[samples.length];
		for (int i = 0; i < samples.length; i++) {
			workerSubsets[i] = new WorkerSubset(samples[i], original, newFocus, focusedIndices, newFingerPrint, overwriteExisting);
		}
		return workerSubsets;
	}

	private static class WorkerSubset implements Callable<FocusedSample> {
		private String sample;
		private Project original;
		private Project newFocus;
		private int[] focusedIndices;
		private long newFingerPrint;
		private boolean overwriteExisting;

		public WorkerSubset(String sample, Project original, Project newFocus, int[] focusedIndices, long newFingerPrint, boolean overwriteExisting) {
			super();
			this.sample = sample;
			this.original = original;
			this.newFocus = newFocus;
			this.focusedIndices = focusedIndices;
			this.newFingerPrint = newFingerPrint;
			this.overwriteExisting = overwriteExisting;
		}

		@Override
		public FocusedSample call() throws Exception {
			Sample samp = original.getFullSampleFromRandomAccessFile(sample);
			FocusedSample focusedSample = null;
			if (samp != null) {
				String newSampleFileName = newFocus.getDir(newFocus.SAMPLE_DIRECTORY, true, true) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION;
				focusedSample = new FocusedSample(focusedIndices, samp, newFingerPrint, overwriteExisting);
				focusedSample.saveToRandomAccessFile(newSampleFileName, samp.getSampleName(), original.getLog());
			} else {
				original.getLog().reportTimeError("Could not load sample " + sample);
			}
			// TODO Auto-generated method stub
			return focusedSample;
		}

	}

	public static boolean focusAProject(Project original, Project newFocus, String[] markersToUse, boolean[] samplesToUse, long newFingerPrint, int numThreads, boolean overwriteExisting, Logger log) {
		boolean focused = true;
		if (original.getProjectDir().equals(newFocus.getProjectDir())) {
			log.reportTimeError("The focused project must have a different project directory than the original, halting");
			focused = false;
		} else if (original.getDir(original.SAMPLE_DIRECTORY).equals(newFocus.getDir(newFocus.SAMPLE_DIRECTORY))) {
			log.reportTimeError("The focused project must have a different sample directory than the original, halting");
			focused = false;
		} else if (markersToUse == null || samplesToUse == null) {
			log.reportTimeError("Please provide subsets... if you want they can be all markers or samples");
			focused = false;
		} else {
			log.reportTimeInfo("Markers to export = " + markersToUse.length);
			log.reportTimeInfo("Samples to export = " + Array.booleanArraySum(samplesToUse));
			int[] markerToUseIndices = ext.indexLargeFactors(markersToUse, original.getMarkerNames(), true, log, true, false);
			WorkerHive<FocusedSample> hive = new WorkerHive<FocusedSample>(numThreads, 10, log);
			WorkerSubset[] workerSubsets = getWorkers(original, newFocus, markerToUseIndices, samplesToUse, newFingerPrint, overwriteExisting);
			hive.addCallables(workerSubsets);
			hive.setReportEvery(100);
			hive.execute(true);
			FocusedSample[] focusedSamples = hive.getResults().toArray(new FocusedSample[hive.getResults().size()]);
			writeOutliers(newFocus, focusedSamples);
			Array.toByteArray(new ArrayList<Byte>());
			return focused;
		}
		return focused;
	}

	private static void writeOutliers(Project newFocus, FocusedSample[] focusedSamples) {
		Hashtable<String, Float> outliers = new Hashtable<String, Float>();
		for (int i = 0; i < focusedSamples.length; i++) {
			outliers.putAll(focusedSamples[i].getOutliers());
		}
		if (outliers.size() > 0) {
			String outlierFileName = newFocus.getDir(newFocus.SAMPLE_DIRECTORY, true) + "outliers.ser";
			newFocus.getLog().reportTimeInfo("Writing outliers to " + outlierFileName);
			Files.writeSerial(outliers, outlierFileName);
		}
	}

	public static void test() {
		Project proj = new Project(null, false);
		Project proj2 = new Project(null, false);
		proj2.setProperty(proj.PROJECT_DIRECTORY, proj.getProjectDir() + "subSample/");
		String[] markers = proj.getTargetMarkers();
		focusAProject(proj, proj2, markers, null, proj.getMarkerSet().getFingerprint(), 8, true, proj.getLog());
		proj2.getMarkerLookup();
	}

	public static void main(String[] args) {
		test();
	}

}
