package cnv.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import stats.Rscript.COLUMNS_MULTIPLOT;
import stats.Rscript.PLOT_DEVICE;
import stats.Rscript.RScatter;
import stats.Rscript.RScatters;
import stats.Rscript.SCATTER_TYPE;
import common.Array;
import common.Files;
import common.PSF;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;
import cnv.analysis.PennCNV;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.qc.GcAdjustor;
import cnv.qc.LrrSd;
import cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
//import cnv.qc.MarkerMetrics;
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
		String outliersSer = projCorrected.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser";

		GCProducer producer = new GCProducer(projOriginal, projCorrected, gcModel, !Files.exists(outliersSer));
		WorkerTrain<GcCorrectedSample> train = new WorkerTrain<GcCorrectedSample>(producer, numThreads, 2, projOriginal.getLog());
		Hashtable<String, Float> outliers = new Hashtable<String, Float>();
		String[] samples = projOriginal.getSamples();
		String firstSampleFile = projCorrected.SAMPLE_DIRECTORY.getValue() + samples[0] + Sample.SAMPLE_DATA_FILE_EXTENSION;
		if (!Files.exists(firstSampleFile)) {
			int numSamples = samples.length;
			int index = 0;
			while (train.hasNext()) {
				index++;
				if (index % 50 == 0) {
					projOriginal.getLog().reportTimeInfo(index + " of " + numSamples + " have been corrected");
				}
				GcCorrectedSample gcCorrectedSample = train.next();
				outliers.putAll(gcCorrectedSample.getOutliers());
			}
			if (!Files.exists(outliersSer)) {
				Files.writeSerial(outliers, outliersSer);
			} else {
				projOriginal.getLog().reportTimeWarning("Did not write outliers, " + outliersSer + " already exists");
			}
			if (!Files.exists(projCorrected.MARKER_DATA_DIRECTORY.getValue() + "markers.0.mdRAF")) {
				TransposeData.transposeData(projCorrected, 2000000000, false);
			}
		} else {
			projOriginal.getLog().reportTimeWarning(firstSampleFile + " exists, skipping correction");
		}

		summarizeMetrics(numThreads);
	}

	private void summarizeMetrics(int numThreads) {
		projOriginal.SAMPLE_QC_FILENAME.setValue(ext.addToRoot(projCorrected.SAMPLE_QC_FILENAME.getValue(), ".priorToCorrection"));
		if (!Files.exists(projCorrected.SAMPLE_QC_FILENAME.getValue())) {
			LrrSd.init(projCorrected, null, null, null, null, numThreads);
		}
		if (!Files.exists(projOriginal.SAMPLE_QC_FILENAME.getValue())) {
			LrrSd.init(projOriginal, null, null, null, null, numThreads);
		}

		String comboQC = ext.addToRoot(projCorrected.SAMPLE_QC_FILENAME.getValue(), ".combinedQC");
		String[] orginalFiles = new String[] { projOriginal.SAMPLE_QC_FILENAME.getValue(), projCorrected.SAMPLE_QC_FILENAME.getValue() };
		String[] titles = new String[] { "UN_CORRECTED", "GC_CORRECTED" };
		String[] fullHeader = Array.concatAll(new String[] { LrrSd.SAMPLE_COLUMN }, LrrSd.NUMERIC_COLUMNS);
		int[] indices = ext.indexFactors(fullHeader, Files.getHeaderOfFile(orginalFiles[0], projOriginal.getLog()), true, false);
		String[][] newColums = Files.paste(orginalFiles, comboQC, indices, 0, titles, new String[] { LrrSd.SAMPLE_COLUMN }, projOriginal.getLog());

		String status = "STATUS";
		String comboBox = comboQC + ".box";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(comboBox));
			writer.println(status + "\t" + Array.toStr(fullHeader));
			for (int i = 0; i < orginalFiles.length; i++) {

				BufferedReader reader = Files.getAppropriateReader(orginalFiles[i]);
				reader.readLine();
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("[\\s]+");
					if (!line[indices[0]].equals(LrrSd.SAMPLE_COLUMN)) {
						writer.println(titles[i] + "\t" + Array.toStr(Array.subArray(line, indices)));
					}
				}
				reader.close();

			}
			writer.close();

		} catch (FileNotFoundException fnfe) {
			projOriginal.getLog().reportError("Error: file(s) \"" + Array.toStr(orginalFiles) + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			projOriginal.getLog().reportError("Error reading file(s) \"" + Array.toStr(orginalFiles) + "\"");
			return;
		} catch (Exception e) {
			projOriginal.getLog().reportError("Error writing to " + comboBox);
			projOriginal.getLog().reportException(e);
		}

		ArrayList<RScatter> rScatters = new ArrayList<RScatter>();
		String gcLookDir = projCorrected.PROJECT_DIRECTORY.getValue() + "gc_analysis/";
		new File(gcLookDir).mkdirs();

		for (int i = 1; i < newColums[0].length; i++) {

			String curFile = gcLookDir + "gc_" + LrrSd.NUMERIC_COLUMNS[i - 1];
			RScatter rScatter = new RScatter(comboQC, curFile + ".rscript", ext.removeDirectoryInfo(curFile), curFile + ".pdf", newColums[1][i], new String[] { newColums[0][i] }, SCATTER_TYPE.POINT, projCorrected.getLog());
			rScatter.setTitle("Original V Corrected " + LrrSd.NUMERIC_COLUMNS[i - 1]);
			rScatter.setxLabel("Corrected " + LrrSd.NUMERIC_COLUMNS[i - 1]);
			rScatter.setyLabel("Original " + LrrSd.NUMERIC_COLUMNS[i - 1]);
			rScatters.add(rScatter);

			String curFileBox = gcLookDir + "gc_" + LrrSd.NUMERIC_COLUMNS[i - 1] + "_box";
			RScatter rScatterBox = new RScatter(comboBox, curFileBox + ".rscript", ext.removeDirectoryInfo(curFileBox), curFileBox + ".pdf", status, new String[] { LrrSd.NUMERIC_COLUMNS[i - 1] }, SCATTER_TYPE.BOX, projCorrected.getLog());
			rScatterBox.setTitle("Original V Corrected " + LrrSd.NUMERIC_COLUMNS[i - 1]);
			rScatterBox.setyLabel(LrrSd.NUMERIC_COLUMNS[i - 1]);
			rScatters.add(rScatterBox);
		}
		String finalRoot = gcLookDir + "Gc_summary";
		RScatters rScatters2 = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]), finalRoot + ".rscript", finalRoot + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, projCorrected.getLog());
		rScatters2.execute();
	}

	private static class GCProducer implements Producer<GcCorrectedSample> {
		private Project projOriginal;
		private Project projCorrected;
		private GcModel gcmodel;
		private String[] samples;
		private int index;
		private boolean loadOutliers;

		private GCProducer(Project projOriginal, Project projCorrected, GcModel gcmodel, boolean loadOutliers) {
			super();
			this.projOriginal = projOriginal;
			this.projCorrected = projCorrected;
			this.gcmodel = gcmodel;
			this.samples = projOriginal.getSamples();
			this.loadOutliers = loadOutliers;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<GcCorrectedSample> next() {
			GcWorker worker = new GcWorker(projOriginal, projCorrected, samples[index], gcmodel, loadOutliers);

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
		private boolean loadOutliers;

		public GcWorker(Project projOriginal, Project projCorrected, String sample, GcModel gcmodel, boolean loadOutliers) {
			super();
			this.projOriginal = projOriginal;
			this.projCorrected = projCorrected;
			this.gcmodel = gcmodel;
			this.sample = sample;
			this.loadOutliers = loadOutliers;

		}

		@Override
		public GcCorrectedSample call() throws Exception {

			String newSampleFile = projCorrected.SAMPLE_DIRECTORY.getValue(true, false) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION;

			Sample correctedSamp = null;

			Hashtable<String, Float> outliers = null;
			if (!Files.exists(newSampleFile)) {
				Sample curSample = projOriginal.getFullSampleFromRandomAccessFile(sample);

				GcAdjustor gcAdjustor = GcAdjustor.getComputedAdjustor(projOriginal, curSample, null, gcmodel, GC_CORRECTION_METHOD.GENVISIS_GC, true, true, false);
				outliers = new Hashtable<String, Float>();
				correctedSamp = new Sample(curSample.getSampleName(), curSample.getFingerprint(), curSample.getGCs(), curSample.getXs(), curSample.getYs(), curSample.getBAFs(), Array.toFloatArray(gcAdjustor.getCorrectedIntensities()), curSample.getForwardGenotypes(), curSample.getAB_Genotypes(), curSample.getCanXYBeNegative());
				correctedSamp.saveToRandomAccessFile(newSampleFile, outliers, curSample.getSampleName());
			} else {
				correctedSamp = projCorrected.getFullSampleFromRandomAccessFile(sample);
				if (loadOutliers) {
					outliers = Sample.loadOutOfRangeValuesFromRandomAccessFile(newSampleFile);
				}
				if (outliers == null) {
					outliers = new Hashtable<String, Float>();
				}
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
		//private Sample sample;

		public GcCorrectedSample(Hashtable<String, Float> outliers, Sample sample) {
			super();
			this.outliers = outliers;
			//this.sample = sample;
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
			if (!Files.exists(proj.GC_MODEL_FILENAME.getValue())) {
				String gcBase = Files.exists("N:/statgen/NCBI/hg19.gc5Base.txt.gz") ? "N:/statgen/NCBI/hg19.gc5Base.txt.gz" : "/home/pankrat2/public/bin/NCBI/hg19.gc5Base.txt.gz";
				PennCNV.gcModel(proj, gcBase, proj.GC_MODEL_FILENAME.getValue(), 100);
			}
			gcCorrect(proj, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
