package org.genvisis.affy;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Date;

import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.analysis.CentroidCompute.CentroidBuilder;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.Sort;
import org.genvisis.stats.Maths;

public class AffyCentroids implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String[] AFFY_CENTROID_SUFFIXES = {	"Name", "AA T Mean", "AA R Mean",
																													"AB T Mean", "AB R Mean", "BB T Mean",
																													"BB R Mean"};

	private final float[][][] AffyCentroids; // marker, genotype (0=AA, 1=AB, 2=BB), coordinates
																						// (0=Mean Theta, 1=Mean R) (a.k.a. follows the suffix
																						// order above)
	private final long fingerprint;

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
		SerializedFiles.writeSerial(this, filename);
	}

	public static AffyCentroids load(String filename, boolean jar) {
		return (AffyCentroids) SerializedFiles.readSerial(filename, jar, true);
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
				return 0.5f + 0.5f	* (theta - AffyCentroids[1][0])
											/ (AffyCentroids[2][0] - AffyCentroids[1][0]);
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
				estimatedR = AffyCentroids[0][1]
											+ (theta - AffyCentroids[0][0])	* (AffyCentroids[1][1] - AffyCentroids[0][1])
												/ (AffyCentroids[1][0] - AffyCentroids[0][0]);
			} else {
				estimatedR = AffyCentroids[0][1];
			}
			if (estimatedR < 0) {
				estimatedR = r;
				log.report("Warning - estimatedR < 0 ");
			}
		} else {
			if (AffyCentroids[2][0] - AffyCentroids[1][0] != 0) {
				estimatedR = AffyCentroids[1][1]
											+ (theta - AffyCentroids[1][0])	* (AffyCentroids[2][1] - AffyCentroids[1][1])
												/ (AffyCentroids[2][0] - AffyCentroids[1][0]);
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

	public static float[] getAFFYBAF(	String[] markerNames, float[][][] affyCents, float[] Xs,
																		float[] Ys, int i, Logger log) {
		float[] AFFYBAFs;
		AFFYBAFs = new float[markerNames.length];
		for (int k = 0; k < markerNames.length; k++) {
			if (affyCents[k] == null) {
				AFFYBAFs[k] = Float.NaN;
			} else {
				if (markerNames[k].startsWith("CN_")) {
					AFFYBAFs[k] = 0;
				} else if (!markerNames[k].startsWith("CN_")) {
					AFFYBAFs[k] = calcBAF(calcTheta(Xs[k], Ys[k]), affyCents[k]);
				}
			}
		}
		return AFFYBAFs;
	}

	public static float[] getAFFYLRR(	String[] markerNames, float[][][] affyCents, float[] Xs,
																		float[] Ys, int i, Logger log) {
		float[] AFFYLRRs;
		AFFYLRRs = new float[Xs.length];
		for (int k = 0; k < markerNames.length; k++) {
			if (affyCents[k] == null) {
				AFFYLRRs[k] = Float.NaN;
			} else {
				if (markerNames[k].startsWith("CN_")) {
					// if a copy number probeset, take the log2 intensity value - median log2 intensity
					AFFYLRRs[k] = (float) (Maths.log2(Xs[k]) - affyCents[k][0][1]);
				} else if (!markerNames[k].startsWith("CN_")) {
					AFFYLRRs[k] = calcLRR(calcTheta(Xs[k], Ys[k]), calcR(Xs[k], Ys[k]), affyCents[k], log);
				}
			}
		}
		return AFFYLRRs;
	}

	public static void recompute(Project proj, String centroidsFile) {
		MarkerSet markerSet;
		Centroids affyCentroids;
		Logger log;

		log = proj.getLog();
		markerSet = proj.getMarkerSet();
		markerSet.getMarkerNames();
		affyCentroids = Centroids.load(centroidsFile, proj.JAR_STATUS.getValue());

		if (affyCentroids.getFingerprint() != markerSet.getFingerprint()) {
			log.reportError("Error - fingerprint for Centroids file '"	+ centroidsFile
											+ "' does not match the fingerprint for the current MarkerSet");
			System.exit(1);
		}
		affyCentroids.getCentroids();
		Centroids.recompute(proj, centroidsFile);
		//
		// samples = proj.getSamples();
		// Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();
		// for (int i = 0; i < samples.length; i++) {
		// log.report(samples[i]);
		// original = proj.getFullSampleFromRandomAccessFile(samples[i]);
		// Xs = original.getXs();
		// Ys = original.getYs();
		// AFFYBAFs = getAFFYBAF(markerNames, affyCents, Xs, Ys, i, log);
		// AFFYLRRs = getAFFYLRR(markerNames, affyCents, Xs, Ys, i, log);
		// sample = new Sample(original.getSampleName(), original.getFingerprint(), original.getGCs(),
		// original.getXs(), original.getYs(), AFFYBAFs, AFFYLRRs, original.getForwardGenotypes(),
		// original.getAB_Genotypes(), original.getCanXYBeNegative());
		// sample.saveToRandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(false, true) +
		// original.getSampleName() + Sample.SAMPLE_DATA_FILE_EXTENSION, allOutliers,
		// original.getSampleName());
		// }
		// if (allOutliers.size() > 0) {
		// Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
		// }
	}

	public static void parseCentroidsFilteredSamples(	Project proj, double missingnessThreshold,
																										double confThreshold) {
		Logger log = proj.getLog();

		if (!proj.getSampleData(1, false).hasExcludedIndividuals()) {
			log.reportError("Error - cannot exclude individuals for centroid computations , no factor named 'Exclude/CLASS=Exclude' in Sample Data");
			System.exit(1);
		} else {
			SampleData sampleData = proj.getSampleData(1, false);
			SampleList sampleList = proj.getSampleList();
			String[] samples = sampleList.getSamples();
			boolean[] samplesToBeUsed = new boolean[samples.length];
			int use = 0;
			for (int i = 0; i < samples.length; i++) {
				samplesToBeUsed[i] = !sampleData.individualShouldBeExcluded(samples[i]);
				if (samplesToBeUsed[i]) {
					use++;
				}
			}
			log.report("Info - generating new cluster centers using " + use + " individuals");
			parseCentroids(proj, samplesToBeUsed, missingnessThreshold, confThreshold);
		}
	}

	public static void parseCentroids(Project proj, boolean[] samplesToBeUsed,
																		double missingnessThreshold, double confThreshold) {
		String[] samples;
		SampleList sampleList;
		MarkerSet markerSet;
		Logger log;

		log = proj.getLog();
		new Date().getTime();
		log.report("Computing centroids from intensity means");
		sampleList = proj.getSampleList();
		samples = sampleList.getSamples();
		proj.getSampleData(0, false);
		if (samples.length != samplesToBeUsed.length) {
			log.reportError("Error - mismatched number of samples in project versus sample mask");
			System.exit(1);
		}
		markerSet = proj.getMarkerSet();
		markerSet.getMarkerNames();
		CentroidBuilder builder = new CentroidBuilder();
		CentroidCompute.computeAndDumpCentroids(proj,
																						new String[] {proj.CUSTOM_CENTROIDS_FILENAME.getValue()},
																						new CentroidBuilder[] {builder}, 2, 2);
	}

	private static boolean checkSexMarker(Project proj, int sex, MarkerData markerData) {
		byte chr;
		boolean use = true;
		chr = markerData.getChr();
		if (chr == 23 && sex != 2) {
			use = false;
		}
		if (chr == 24 && sex != 1) {
			use = false;
		}
		return use;
	}

	private static void assignRegularCentroid(float[][] centroid, double[] meanThetas,
																						double[] meanRs, int[] counts, Logger log) {
		for (int k = 0; k < 3; k++) {
			if (counts[k + 2] > 0) {
				// writer.println(markerName+ "\t"+(float)meanThetas[k+2] +"\t" +(float)meanRs[k+2] );
				centroid[k] = new float[] {(float) meanThetas[k + 2], (float) meanRs[k + 2]};
			} else {
				log.reportError("Error assiging regular centroids");
				System.exit(1);
			}
		}
	}

	private static float getAltRs(float[][] centroid, double[] meanRs, int[] counts, int checkIndex,
																int primaryAlt, int secondaryAlt, Logger log) {
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
					return q * array[keys[(int) Math.floor(index) - 1]]
									+ (1 - q) * array[keys[(int) Math.ceil(index) - 1]];
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			return -1234567890;
		}
	}

	// only get samps to be used
	public static float[] log2Array(Project proj, MarkerData markerData, boolean[] samplesToBeUsed,
																	String[] samples, SampleData sampleData) {
		float[] xs = markerData.getXs();
		float[] log2Xs;
		ArrayList<Float> log2XsAL = new ArrayList<Float>();
		for (int k = 0; k < xs.length; k++) {
			if (checkSexMarker(proj, sampleData.getSexForIndividual(samples[k]), markerData)
					&& samplesToBeUsed[k]) {
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
		String usage = "\n"	+ "affy.AffyCentroids requires 0-1 arguments\n"
										+ "   (1) project (i.e. proj=" + filename + " (default))\n"
										+ "   (2) filename (i.e. file=" + centFile + " (default))\n" + " OR\n"
										+ "   (2) generate centroids from genotypes (i.e. -fromGenotypes (not the default))\n"
										+ " OR\n"
										+ "   (2) file with intensity only flags (i.e. flags=intensityFlags.dat (not the default))\n"
										+ "   (3) centroid file for clustered markers (see \"GENOTYPE_CENTROIDS_FILENAME\" property in the Project properties file)\n"
										+ "   (4) centroid file for intensity only markers (see \"GENOTYPE_CENTROIDS_FILENAME\" property in the Project properties file)\n"
										+ " OR\n"
										+ "   (2) recompute BAF/LRR and generate new Sample files using these centroids (i.e. compute=genotype.cent (not the default))\n"
										+ "";
		String logfile = null;
		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("file=")) {
				centFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("callConf=")) {
				callConfidence = Double.parseDouble(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("-fromGenotypes")) {
				fromGenotypes = true;
				numArgs--;
			} else if (arg.startsWith("-FilteredfromGenotypes")) {
				sampleFilter = true;
				numArgs--;
			} else if (arg.startsWith("compute=")) {
				compute = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		proj = new Project(filename, logfile, false);
		try {
			if (fromGenotypes) {
				parseCentroids(proj, Array.booleanArray(proj.getSamples().length, true), 1, callConfidence);
			} else if (sampleFilter) {
				parseCentroidsFilteredSamples(proj, 1, callConfidence);
			} else if (!compute.equals("")) {
				recompute(proj, compute);
			}
			// Files.backup(logfile, proj.PROJECT_DIRECTORY.getValue(),
			// proj.PROJECT_DIRECTORY.getValue());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
