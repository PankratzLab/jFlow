package cnv.manage;

import java.io.File;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import common.Array;
import common.Files;
import common.PSF;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.qc.GcAdjustor;
import cnv.qc.LrrSd;
import cnv.qc.MarkerMetrics;
import cnv.qc.GcAdjustor.GcModel;

/**
 * @author lane0212 Correct a project for gc content
 */
public class GcCorrection {
	private Project projOriginal;
	private GcModel gcModel;
	private Project projCorrected;

	public GcCorrection(Project proj, GcModel gcModel) {
		super();
		this.projOriginal = proj;
		this.gcModel = gcModel;
		this.projCorrected = prepareNewProject(proj);
	}

	public void correct(int numThreads) {
		GCProducer producer = new GCProducer(projOriginal, projCorrected, gcModel);
		WorkerTrain<GcCorrectedSample> train = new WorkerTrain<GcCorrectedSample>(producer, numThreads, 2, projOriginal.getLog());
		Hashtable<String, Float> outliers = new Hashtable<String, Float>();
		int numSamples = projOriginal.getSamples().length;
		int index = 0;
		while (train.hasNext()) {
			index++;
			if (index % 50 == 0) {
				projOriginal.getLog().reportTimeInfo(index + " of " + numSamples + " have been corrected");
			}
			GcCorrectedSample gcCorrectedSample = train.next();
			outliers.putAll(gcCorrectedSample.getOutliers());
		}
		Files.writeSerial(outliers, projCorrected.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
		TransposeData.transposeData(projCorrected, 2000000000, false);
		// LrrSd.init(projCorrected, null, null, numThreads);
		// MarkerMetrics.fullQC(projCorrected, null, null, numThreads);
	}

	private static class GCProducer implements Producer<GcCorrectedSample> {
		private Project projOriginal;
		private Project projCorrected;
		private GcModel gcmodel;
		private String[] samples;
		private int index;

		private GCProducer(Project projOriginal, Project projCorrected, GcModel gcmodel) {
			super();
			this.projOriginal = projOriginal;
			this.projCorrected = projCorrected;
			this.gcmodel = gcmodel;
			this.samples = projOriginal.getSamples();
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<GcCorrectedSample> next() {
			GcWorker worker = new GcWorker(projOriginal, projCorrected, samples[index], gcmodel);

			index++;
			return worker;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static class GcWorker implements Callable<GcCorrectedSample> {
		private Project projOriginal;
		private Project projCorrected;

		private String sample;
		private GcModel gcmodel;

		public GcWorker(Project projOriginal, Project projCorrected, String sample, GcModel gcmodel) {
			super();
			this.projOriginal = projOriginal;
			this.projCorrected = projCorrected;
			this.gcmodel = gcmodel;
			this.sample = sample;
		}

		@Override
		public GcCorrectedSample call() throws Exception {
			Sample curSample = projOriginal.getFullSampleFromRandomAccessFile(sample);

			String newSampleFile = projCorrected.SAMPLE_DIRECTORY.getValue(true, false) + curSample.getSampleName() + Sample.SAMPLE_DATA_FILE_EXTENSION;

			Sample correctedSamp = null;

			Hashtable<String, Float> outliers = null;
			if (!Files.exists(newSampleFile)) {
				GcAdjustor gcAdjustor = GcAdjustor.getComputedAdjustor(projOriginal, curSample, gcmodel, false, true, true, false);
				outliers = new Hashtable<String, Float>();
				correctedSamp = new Sample(curSample.getSampleName(), curSample.getFingerprint(), curSample.getGCs(), curSample.getXs(), curSample.getYs(), curSample.getBAFs(), Array.toFloatArray(gcAdjustor.getCorrectedIntensities()), curSample.getForwardGenotypes(), curSample.getAB_Genotypes(), curSample.getCanXYBeNegative());
				correctedSamp.saveToRandomAccessFile(newSampleFile, outliers, curSample.getSampleName());
			} else {
				correctedSamp = projCorrected.getFullSampleFromRandomAccessFile(sample);
				outliers = Sample.loadOutOfRangeValuesFromRandomAccessFile(newSampleFile);
			}
			return new GcCorrectedSample(outliers, correctedSamp);
		}

	}

	private static Project prepareNewProject(Project projOriginal) {
		String newProjectFile = ext.addToRoot(projOriginal.PROJECT_PROPERTIES_FILENAME.getValue(), ".gc_corrected");
		Files.copyFileUsingFileChannels(projOriginal.PROJECT_PROPERTIES_FILENAME.getValue(), newProjectFile, projOriginal.getLog());
		Project projCorrected = new Project(newProjectFile, false);
		String newDir = projOriginal.PROJECT_DIRECTORY.getValue() + "gc_corrected/";
		projOriginal.getLog().reportTimeInfo("Preparing project " + newProjectFile + " in " + newDir);
		new File(newDir).mkdirs();
		projCorrected.PROJECT_DIRECTORY.setValue(newDir);
		projCorrected.PROJECT_NAME.setValue(projOriginal.PROJECT_NAME.getValue() + "_gc_corrected");
		projCorrected.saveProperties();
		projOriginal.copyBasicFiles(projCorrected, false);
		return projCorrected;
	}

	private static class GcCorrectedSample {
		private Hashtable<String, Float> outliers;
		private Sample sample;

		public GcCorrectedSample(Hashtable<String, Float> outliers, Sample sample) {
			super();
			this.outliers = outliers;
			this.sample = sample;
		}

		public Hashtable<String, Float> getOutliers() {
			return outliers;
		}
	}

	public static void gcCorrect(Project proj, int numThreads) {
		String gcModelFile = proj.GC_MODEL_FILENAME.getValue();
		if (!Files.exists(gcModelFile)) {
			proj.getLog().reportFileNotFound(gcModelFile);

		} else {
			GcModel gcModel = GcModel.populateFromFile(gcModelFile, true, proj.getLog());
			GcCorrection gcCorrection = new GcCorrection(proj, gcModel);
			gcCorrection.correct(numThreads);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		int numThreads = 8;
		String usage = "\n" + "one.JL.GcCorrection requires 0-1 arguments\n";
		usage += "   (1) project filename to gc correct (i.e. proj=" + filename + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(2, numThreads);
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
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
			proj.PROJECT_PROPERTIES_FILENAME.setValue(filename);
			// PennCNV.gcModel(proj, "N:/statgen/NCBI/hg19.gc5Base.txt.gz", proj.GC_MODEL_FILENAME.getValue() , 100);
			gcCorrect(proj, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
