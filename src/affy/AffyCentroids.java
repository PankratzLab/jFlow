package affy;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Date;

import cnv.filesys.Centroids;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.filesys.SampleList;
import cnv.manage.MarkerDataLoader;
import cnv.var.SampleData;
import stats.Maths;
import common.Array;
import common.Files;
import common.Logger;
import common.Sort;
import common.ext;

public class AffyCentroids implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String[] AFFY_CENTROID_SUFFIXES = { "Name", "AA T Mean", "AA R Mean", "AB T Mean", "AB R Mean", "BB T Mean", "BB R Mean" };

	private float[][][] AffyCentroids; // marker, genotype (0=AA, 1=AB, 2=BB), coordinates (0=Mean Theta, 1=Mean R) (a.k.a. follows the suffix order above)
	private long fingerprint;

	// private static int stopper = 1855448;

	public AffyCentroids(float[][][] AffyCentroids, long fingerprint) {
		this.AffyCentroids = AffyCentroids;
		this.fingerprint = fingerprint;
	}

	public float[][][] getCentroids() {
		return AffyCentroids;
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static AffyCentroids load(String filename, boolean jar) {
		return (AffyCentroids) Files.readSerial(filename, jar, true);
	}

	public static float calcR(float x, float y) {
		return (float) (Math.max(x, 0.0001) + Math.max(y, 0.0001));
	}

	public static float calcTheta(float x, float y) {
		return (float) (Math.atan(Math.max(y, 0.0001) / Math.max(x, 0.0001)) * 2 / Math.PI);
	}

	public static float calcBAF(float theta, float[][] AffyCentroids) {
		if (AffyCentroids[0] != null && theta < AffyCentroids[0][0]) {
			return 0;
		} else if (AffyCentroids[1] != null && theta < AffyCentroids[1][0]) {
			if (AffyCentroids[0] == null) {
				return 0.50f;
			} else {
				return 0.5f * (theta - AffyCentroids[0][0]) / (AffyCentroids[1][0] - AffyCentroids[0][0]);
			}
		} else if (AffyCentroids[2] != null && theta < AffyCentroids[2][0]) {
			if (AffyCentroids[1] == null) {
				return 1.0f;
			} else {
				return 0.5f + 0.5f * (theta - AffyCentroids[1][0]) / (AffyCentroids[2][0] - AffyCentroids[1][0]);
			}
		} else {
			if (AffyCentroids[2] == null) {
				return 0.50f;
			} else {
				return 1;
			}
		}
	}

	public static float calcLRR(float theta, float r, float[][] AffyCentroids, Logger log) {
		float estimatedR;
		// centroid is either all null, or not.
		if (AffyCentroids[0] == null) {
			estimatedR = Float.NaN;
		} else if (theta < AffyCentroids[1][0]) {
			if (AffyCentroids[1][0] - AffyCentroids[0][0] != 0) {
				estimatedR = AffyCentroids[0][1] + (theta - AffyCentroids[0][0]) * (AffyCentroids[1][1] - AffyCentroids[0][1]) / (AffyCentroids[1][0] - AffyCentroids[0][0]);
			} else {
				estimatedR = AffyCentroids[0][1];
			}
			if (estimatedR < 0) {
				estimatedR = r;
				log.report("Warning - estimatedR < 0 ");
			}
		} else {
			if (AffyCentroids[2][0] - AffyCentroids[1][0] != 0) {
				estimatedR = AffyCentroids[1][1] + (theta - AffyCentroids[1][0]) * (AffyCentroids[2][1] - AffyCentroids[1][1]) / (AffyCentroids[2][0] - AffyCentroids[1][0]);
			} else {
				estimatedR = AffyCentroids[1][1];
			}
			if (estimatedR < 0) {
				estimatedR = r;
				log.report("Warning - estimatedR < 0 ");
			}
		}
		return (float) Maths.log2(r / estimatedR);

	}

	public static float[] getAFFYBAF(String[] markerNames, float[][][] affyCents, float[] Xs, float[] Ys, int i, Logger log) {
		float[] AFFYBAFs;
		AFFYBAFs = new float[markerNames.length];
		for (int k = 0; k < markerNames.length; k++) {
			if (markerNames[k].startsWith("CN_")) {
				AFFYBAFs[k] = 0;
			} else if (!markerNames[k].startsWith("CN_")) {
				AFFYBAFs[k] = calcBAF(calcTheta(Xs[k], Ys[k]), affyCents[k]);
			}
		}
		return AFFYBAFs;
	}

	public static float[] getAFFYLRR(String[] markerNames, float[][][] affyCents, float[] Xs, float[] Ys, int i, Logger log) {
		float[] AFFYLRRs;
		AFFYLRRs = new float[Xs.length];
		for (int k = 0; k < markerNames.length; k++) {
			if (markerNames[k].startsWith("CN_")) {
				// if a copy number probeset, take the log2 intensity value - median log2 intensity
				AFFYLRRs[k] = (float) (Maths.log2(Xs[k]) - affyCents[k][0][1]);
			} else if (!markerNames[k].startsWith("CN_")) {
				AFFYLRRs[k] = calcLRR(calcTheta(Xs[k], Ys[k]), calcR(Xs[k], Ys[k]), affyCents[k], log);
			}
		}
		return AFFYLRRs;
	}

	public static void recompute(Project proj, String centroidsFile, Logger log) {
		MarkerSet markerSet;
		Centroids affyCentroids;
		Sample original, sample;
		String[] samples, markerNames;
		float[][][] affyCents;
		float[] Xs;
		float[] Ys;
		float[] AFFYBAFs;
		float[] AFFYLRRs;
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		affyCentroids = Centroids.load(centroidsFile, proj.getJarStatus());
		if (affyCentroids.getFingerprint() != markerSet.getFingerprint()) {
			log.reportError("Error - fingerprint for Centroids file '" + centroidsFile + "' does not match the fingerprint for the current MarkerSet");
			System.exit(1);
		}
		affyCents = affyCentroids.getCentroids();
		samples = proj.getSamples();
		for (int i = 0; i < samples.length; i++) {
			log.report(samples[i]);
			original = proj.getFullSampleFromRandomAccessFile(samples[i]);
			Xs = original.getXs();
			Ys = original.getYs();
			AFFYBAFs = getAFFYBAF(markerNames, affyCents, Xs, Ys, i, log);
			AFFYLRRs = getAFFYLRR(markerNames, affyCents, Xs, Ys, i, log);
			sample = new Sample(original.getSampleName(), original.getFingerprint(), original.getGCs(), original.getXs(), original.getYs(), AFFYBAFs, AFFYLRRs, original.getForwardGenotypes(), original.getAB_Genotypes(), original.getCanXYBeNegative());
			sample.saveToRandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY) + original.getSampleName() + Sample.SAMPLE_DATA_FILE_EXTENSION);
		}
	}

	public static void parseCentroidsFilteredSamples(Project proj, double missingnessThreshold, double confThreshold, Logger log) {
		if (!proj.getSampleData(1, false).hasExcludedIndividuals()) {
			log.reportError("Error - cannot exclude individuals for centroid computations , no factor named 'Exclude/CLASS=Exclude' in Sample Data");
			System.exit(1);
		} else {
			SampleData sampleData = proj.getSampleData(1, false);
			SampleList sampleList = proj.getSampleList();
			String[] samples = sampleList.getSamples();
			boolean[] samplesToBeUsed = new boolean[samples.length];
			int use =0;
			for (int i = 0; i < samples.length; i++) {
				samplesToBeUsed[i] = !sampleData.individualShouldBeExcluded(samples[i]);
				if(samplesToBeUsed[i]){
					use++;
				}
			}
			log.report("Info - generating new cluster centers using " + use + " individuals");
			parseCentroids(proj, samplesToBeUsed, missingnessThreshold, confThreshold, log);
		}
	}

	public static void parseCentroids(Project proj, boolean[] samplesToBeUsed, double missingnessThreshold, double confThreshold, Logger log) {
		String[] samples, markerNames;
		float[][][] centroids;
		int count;
		SampleList sampleList;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		long time;
		MarkerSet markerSet;
		time = new Date().getTime();
		log.report("Computing centroids from intensity means");
		sampleList = proj.getSampleList();
		samples = sampleList.getSamples();
		if (samples.length != samplesToBeUsed.length) {
			log.reportError("Error - mismatched number of samples in project versus sample mask");
			System.exit(1);
		}
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		centroids = new float[markerNames.length][][];
		markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames, log);

		time = new Date().getTime();
		for (int i = 0; i < markerNames.length; i++) {
			markerData = markerDataLoader.requestMarkerData(i);
			centroids[i] = computeCluster(proj, samplesToBeUsed, missingnessThreshold, samples, markerData, i, confThreshold, log);
			markerDataLoader.releaseIndex(i);
		}
		count = 0;
		for (int i = 0; i < centroids.length; i++) {
			if (centroids[i] == null) {
				if (count == 0) {
					log.report("The following marker(s) could not be computed:");
				}
				count++;
			}
		}
		if (count > 0) {
			log.report("Computed mean genotyped centroids for " + (centroids.length - count) + " of " + centroids.length + " markers, " + count + " missing");
		} else {
			log.report("Computed mean genotyped centroids for all " + centroids.length + " markers");
		}
		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.getFilename(Project.GENOTYPE_CENTROIDS_FILENAME));
		log.report("Computation took " + ext.getTimeElapsed(time));
	}

	private static float[][] computeCluster(Project proj, boolean[] samplesToBeUsed, double missingnessThreshold, String[] samples, MarkerData markerData, int i, double confThreshold, Logger log) {
		float[][] centroid;
		float[] thetas;
		float[] rs;
		float[] confs;
		double[] meanThetas;
		double[] meanRs;
		byte[] genotypes;
		int[] counts;
		PrintWriter writer = null;
		String markerName;
		markerName = markerData.getMarkerName();
		meanThetas = new double[5];
		meanRs = new double[5];
		counts = new int[5];
		boolean use = true;

		if (!markerName.startsWith("CN_")) {
			for (int k = 0; k < samples.length; k++) {
				// use if X chromosome and sex is female, dont use if X chromsome and sex is male
				// use if Ychromosome and sex is male, dont use if Ychromsome and sex is female
				use = checkSexMarker(proj, samples[k], markerData);
				if (!use || !samplesToBeUsed[k]) {
					continue;
				}
				confs = markerData.getGCs();
				genotypes = markerData.getAB_Genotypes();
				thetas = markerData.getThetas();
				rs = markerData.getRs();
				// maybe alter the confidence checks
				if (!Float.isNaN(thetas[k]) && !Float.isNaN(rs[k]) && !Float.isNaN(confs[k]) && confs[k] >= confThreshold) {
					meanThetas[0] += thetas[k];
					meanRs[0] += rs[k];
					counts[0]++;
					meanThetas[genotypes[k] + 2] += thetas[k];
					meanRs[genotypes[k] + 2] += rs[k];
					counts[genotypes[k] + 2]++;
				}
			}
			for (int k = 0; k < 5; k++) {
				meanThetas[k] /= counts[k];
				meanRs[k] /= counts[k];
			}
		} else if (markerName.startsWith("CN_")) {
			// checking sex in log2Array...
			float[] log2Xs = log2Array(proj, markerData, samplesToBeUsed, samples);
			// for computing LRR
			meanRs[2] = median(log2Xs);
			meanRs[3] = Array.mean(log2Xs);
			// for vis
			meanRs[4] = median(markerData.getXs());
			meanThetas[2] = 0.5;
			meanThetas[3] = 0.5;
			meanThetas[4] = 0.5;

			counts[2]++;
			counts[3]++;
			counts[4]++;
			counts[0]++;
		}
		centroid = getCentroid(missingnessThreshold, meanThetas, meanRs, counts, writer, markerName, log);
		return centroid;
	}

	private static float[][] getCentroid(double missingnessThreshold, double[] meanThetas, double[] meanRs, int[] counts, PrintWriter writer, String markerName, Logger log) {
		float[][] centroid;
		centroid = new float[3][];
		if (counts[1] >= counts[0] * missingnessThreshold) {
			for (int k = 0; k < 3; k++) {
				centroid[k] = new float[] { (float) meanThetas[0], (float) meanRs[0] };
			}
		} else {
			processCentroid(centroid, meanThetas, meanRs, counts, writer, markerName, log);
		}
		return centroid;
	}

	private static boolean checkSexMarker(Project proj, String sample, MarkerData markerData) {
		int sex = 0;
		byte chr;
		boolean use = true;
		sex = proj.getSampleData(0, false).getSexForIndividual(sample);
		chr = markerData.getChr();
		if ((int) chr == 23 && sex != 2) {
			use = false;
		}
		if ((int) chr == 24 && sex != 1) {
			use = false;
		}
		return use;
	}

	private static void processCentroid(float[][] centroid, double[] meanThetas, double[] meanRs, int[] counts, PrintWriter writer, String markerName, Logger log) {
		if (markerName.startsWith("CN_")) {
			centroid[0] = new float[] { (float) meanThetas[2], (float) meanRs[2] };
			centroid[1] = new float[] { (float) meanThetas[3], (float) meanRs[3] };
			centroid[2] = new float[] { (float) meanThetas[4], (float) meanRs[4] };
			// System.out.println(markerName + "\t" + meanRs[4]);

		} else if (counts[2] == 0 && counts[3] == 0 && counts[4] == 0) {
			nullify(centroid);
		} else if (counts[2] == 0 && counts[3] == 0) {
			nullify(centroid);
		} else if (counts[2] == 0 && counts[4] == 0) {
			nullify(centroid);
		} else if (counts[3] == 0 && counts[4] == 0) {
			nullify(centroid);
		} else if (counts[2] == 0) {
			// Estimating the AA cluster
			centroid[0] = new float[] { (float) (meanThetas[3] - 0.3), getAltRs(centroid, meanRs, counts, 2, 3, 4, log) };
			centroid[1] = new float[] { (float) meanThetas[3], (float) meanRs[3] };
			centroid[2] = new float[] { (float) meanThetas[4], (float) meanRs[4] };

		} else if (counts[3] == 0) {
			// Estimating the AB cluster
			centroid[0] = new float[] { (float) meanThetas[2], (float) meanRs[2] };
			if (counts[2] > 0 && counts[4] > 0) {
				centroid[1] = new float[] { (float) ((meanThetas[2] + meanThetas[4]) / 2), getAltRs(centroid, meanRs, counts, 3, 2, 4, log) };
			} else {
				centroid[1] = new float[] { (float) (0.5), getAltRs(centroid, meanRs, counts, 3, 2, 4, log) };
			}
			centroid[2] = new float[] { (float) meanThetas[4], (float) meanRs[4] };
		} else if (counts[4] == 0) {
			// Estimating the BB cluster
			centroid[0] = new float[] { (float) meanThetas[2], (float) meanRs[2] };
			centroid[1] = new float[] { (float) meanThetas[3], (float) meanRs[3] };
			centroid[2] = new float[] { (float) (meanThetas[3] + 0.3), getAltRs(centroid, meanRs, counts, 4, 3, 2, log) };
		}
		else {
			assignRegularCentroid(centroid, meanThetas, meanRs, counts, log);
		}
		if (centroid[2] != null && (centroid[0][0] > centroid[1][0] || centroid[2][0] < centroid[1][0])) {
			nullify(centroid);
		}
	}

	private static void assignRegularCentroid(float[][] centroid, double[] meanThetas, double[] meanRs, int[] counts, Logger log) {
		for (int k = 0; k < 3; k++) {
			if (counts[k + 2] > 0) {
				// writer.println(markerName+ "\t"+(float)meanThetas[k+2] +"\t" +(float)meanRs[k+2] );
				centroid[k] = new float[] { (float) meanThetas[k + 2], (float) meanRs[k + 2] };
			} else {
				log.reportError("Error assiging regular centroids");
				System.exit(1);
			}
		}
	}

	private static float getAltRs(float[][] centroid, double[] meanRs, int[] counts, int checkIndex, int primaryAlt, int secondaryAlt, Logger log) {
		// AssignAltRs(centroid ,meanRs , counts , 2 , 3 ,4 );
		float altR;
		if (counts[primaryAlt] > 0) {
			altR = (float) meanRs[primaryAlt];
		} else if (counts[secondaryAlt] > 0) {
			altR = (float) meanRs[secondaryAlt];
		} else {
			altR = 0;
			log.reportError("Error retrieving alternate R");
			System.exit(1);
		}
		return altR;
	}

	private static void nullify(float[][] centroid) {
		centroid[0] = null;
		centroid[1] = null;
		centroid[2] = null;

	}

	public static float median(float[] array) {
		return (quant(array, (float) 0.50));
	}

	public static float quant(float[] array, float q) {
		int keys[] = Sort.quicksort(array);
		try {
			if (q > 1 || q < 0) {
				return (0);
			} else {
				double index = (array.length + 1) * q;
				if (index - (int) index == 0) {
					return array[keys[(int) index - 1]];
				} else {
					return q * array[keys[(int) Math.floor(index) - 1]] + (1 - q) * array[keys[(int) Math.ceil(index) - 1]];
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			return -1234567890;
		}
	}

	// only get samps to be used
	public static float[] log2Array(Project proj, MarkerData markerData, boolean[] samplesToBeUsed, String[] samples) {
		float[] xs = markerData.getXs();
		float[] log2Xs;
		ArrayList<Float> log2XsAL = new ArrayList<Float>();
		for (int k = 0; k < xs.length; k++) {
			if (checkSexMarker(proj, samples[k], markerData) && samplesToBeUsed[k]) {
				log2XsAL.add((float) Maths.log2(xs[k]));
			}
		}
		log2Xs = new float[log2XsAL.size()];
		for (int i = 0; i < log2Xs.length; i++) {
			log2Xs[i] = log2XsAL.get(i);
		}
		return log2Xs;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "C:/workspace/Genvisis/projects/ARICGenvisis_CEL_11908.properties";
		boolean fromGenotypes = false;
		boolean sampleFilter = false;
		Project proj;
		String compute = "";
		double callConfidence = 0.99;
		String centFile = "C:/data/ARIC/ARICGenvisis_CEL_11908/data/genotype.cent";
		String usage = "\n" + "affy.AffyCentroids requires 0-1 arguments\n" + "   (1) project (i.e. proj=" + filename + " (default))\n" + "   (2) filename (i.e. file=" + centFile + " (default))\n" + " OR\n" + "   (2) generate centroids from genotypes (i.e. -fromGenotypes (not the default))\n" + " OR\n" + "   (2) file with intensity only flags (i.e. flags=intensityFlags.dat (not the default))\n" + "   (3) centroid file for clustered markers (see " + Project.GENOTYPE_CENTROIDS_FILENAME + " in the Project properties file)\n" + "   (4) centroid file for intensity only markers (see " + Project.GENOTYPE_CENTROIDS_FILENAME + " in the Project properties file)\n" + " OR\n" + "   (2) recompute BAF/LRR and generate new Sample files using these centroids (i.e. compute=genotype.cent (not the default))\n" + "";
		String logfile = null;
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("file=")) {
				centFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("callConf=")) {
				callConfidence = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-fromGenotypes")) {
				fromGenotypes = true;
				numArgs--;
			} else if (args[i].startsWith("-FilteredfromGenotypes")) {
				sampleFilter = true;
				numArgs--;
			}
			else if (args[i].startsWith("compute=")) {
				compute = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		proj = new Project(filename, false);
		if (logfile == null) {
			logfile = "AffyCentroidlog.txt";
		}
		Logger log = new Logger(proj.getProjectDir() + logfile);
		try {
			if (fromGenotypes) {
				parseCentroids(proj, Array.booleanArray(proj.getSamples().length, true), 1, 0.99, log);
			} else if (sampleFilter) {
				parseCentroidsFilteredSamples(proj, 1, callConfidence, log);
			} else if (!compute.equals("")) {
				recompute(proj, compute, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
