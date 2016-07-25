package org.genvisis.cnv.filesys;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.genvisis.cnv.analysis.AnalysisFormats;
import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.manage.TextExport;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.*;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.stats.Maths;

/**
 * @author lane0212
 *
 */
public class Centroids implements Serializable, TextExport {
	public static final long serialVersionUID = 1L;
	public static final String[] ILLUMINA_CENTROID_SUFFIXES = {"Name", "AA T Mean", "AA R Mean", "AB T Mean", "AB R Mean", "BB T Mean", "BB R Mean"};
	
	private float[][][] centroids;  // marker, genotype (0=AA, 1=AB, 2=BB), coordinates (0=Mean Theta, 1=Mean R) (a.k.a. follows the suffix order above)
	private long fingerprint;

	public Centroids(float[][][] centroids, long fingerprint) {
		this.centroids = centroids;
		this.fingerprint = fingerprint;
	}

	public float[][][] getCentroids() {
		return centroids;
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}
	
	public void exportToText(Project proj, String outputFile) {
	    PrintWriter writer;
        float[][][] centroids;
        MarkerSet markerSet;
        String[] markerNames;
        Logger log = proj.getLog();
        
        markerSet = proj.getMarkerSet();
        markerNames = markerSet.getMarkerNames();
        centroids = getCentroids();
        
        if (markerNames.length != centroids.length) {
            log.reportError("Error - mismatched number of markers in centroid object and the project's marker set; aborting");
            return;
        }

        if (markerSet.getFingerprint() != getFingerprint()) {
            log.reportError("Error - mismatched marker fingerprints in centroid object and the project's marker set ; aborting");
            return;
        }
        
        try {
            writer = new PrintWriter(new FileWriter(outputFile));
            writer.println("marker_fingerprint="+getFingerprint());
            writer.println("MarkerName\tAA_Theta_Mean\tAA_R_Mean\tAB_Theta_Mean\tAB_R_Mean\tBB_Theta_Mean\tBB_R_Mean");
            for (int i = 0; i < markerNames.length; i++) {
                writer.print(markerNames[i]);
                for (int j = 0; j < 3; j++) {
                    if (centroids[i][j] == null) {
                        writer.print("\t.\t.");
                    } else {
                        writer.print("\t"+centroids[i][j][0]+"\t"+centroids[i][j][1]);
                    }
                }
                writer.println();
            }
            writer.close();
        } catch (Exception e) {
            log.reportException(e);
        }
	}
	
	public static Centroids load(String filename, boolean jar) {
		return (Centroids)SerializedFiles.readSerial(filename, jar, true);
	}
	
	public static float calcR(float x, float y) {
		return (float)(Math.max(x, 0.0001)+Math.max(y, 0.0001));
	}
	
	public static float calcTheta(float x, float y) {
		return (float)(Math.atan(Math.max(y, 0.0001)/Math.max(x, 0.0001))*2/Math.PI);
	}
	
	public static float calcBAF(float theta, float[][] centroids) {
		if (centroids[0] != null && theta < centroids[0][0]) {
			return 0;
		} else if (centroids[1] != null && theta < centroids[1][0]) {
			if (centroids[0] == null) {
				return 0.50f;
			} else {
				return 0.5f*(theta-centroids[0][0])/(centroids[1][0]-centroids[0][0]);
			}
		} else if (centroids[2] != null && theta < centroids[2][0]) {
			if (centroids[1] == null) {
				return 1.0f;
			} else {
				return 0.5f+0.5f*(theta-centroids[1][0])/(centroids[2][0]-centroids[1][0]);
			}
		} else {
			if (centroids[2] == null) {
				return 0.50f;
			} else {
				return 1;
			}
		}
	}

	public static float calcLRR(float theta, float r, float[][] centroids) {
		float estimatedR;
		
//		if (centroids[2][1] < 0.0000001) {
//			centroids[2][0] = 1;
//		}
		
		if (centroids[0] != null && theta < centroids[0][0] || (centroids[1] == null && centroids[2] == null)) {
			estimatedR = centroids[0][1];
		} else if (centroids[1] != null && theta < centroids[1][0]) {
			if (centroids[0] == null) {
				estimatedR = centroids[1][1];
			} else {
				estimatedR = centroids[0][1]+(theta-centroids[0][0])*(centroids[1][1]-centroids[0][1])/(centroids[1][0]-centroids[0][0]);
			}
		} else if (centroids[2] != null && theta < centroids[2][0]) {
			if (centroids[1] == null) {
				estimatedR = centroids[2][1];
			} else {
				estimatedR = centroids[1][1]+(theta-centroids[1][0])*(centroids[2][1]-centroids[1][1])/(centroids[2][0]-centroids[1][0]);
			}
		} else {
			if (centroids[2] == null) {
				estimatedR = centroids[1][1];
			} else {
				estimatedR = centroids[2][1];
			}
		}

		return (float)Maths.log2(r/estimatedR);
	}
	
	public static void parseIlluminaCentroidsFromCSV(Project proj, String filename) {
		BufferedReader reader;
		String[] line, header;
		Hashtable<String, String> hash;
		int[] indices;
		float[][][] centroids;
		String[] markerNames;
		int index;
		boolean missing;
		MarkerSet markerSet;
		
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		centroids = new float[markerNames.length][][];
		
		hash = new Hashtable<String, String>();
		for (int i = 0; i < markerNames.length; i++) {
			hash.put(markerNames[i], i+"");
		}		
		
		try {
			reader = new BufferedReader(new FileReader(proj.PROJECT_DIRECTORY.getValue()+filename));
			header = reader.readLine().trim().split(",");
			indices = Array.intArray(ILLUMINA_CENTROID_SUFFIXES.length, -1);
			for (int i = 0; i<header.length; i++) {
				index = ext.indexOfEndsWith(header[i], ILLUMINA_CENTROID_SUFFIXES, true);
				if (index >= 0) {
					indices[index] = i;
				}
            }
			if (Array.min(indices) == -1) {
				System.err.println("Error - Need a column header ending with the following suffixes; missing at least one");
				System.err.println("        "+Array.toStr(ILLUMINA_CENTROID_SUFFIXES, "  "));
			}
//			for (int i = 0; i < markerNames.length; i++) {
			while (reader.ready()) {
				line = reader.readLine().trim().split(",");
				if (!hash.containsKey(line[indices[0]])) {
					System.err.println("Error - marker '"+line[indices[0]]+"' was not found in MarkerSet");
					System.exit(1);
				}
				index = Integer.parseInt(hash.get(line[indices[0]]));
				centroids[index] = new float[3][2];
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 2; k++) {
						centroids[index][j][k] = Float.parseFloat(line[indices[1+j*2+k]]);
					}
				}
			}
			missing = false;
			for (int i = 0; i<centroids.length; i++) { // might want to generate an error log or display the number if greater than, say, 10
				if (centroids[i] == null) {
					if (!missing) {
						System.err.println("Error - did not find a centroid for the following markers:");
						missing = true;
					}
					System.err.println("  "+markerNames[i]);
				}
				
            }
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in "+proj.PROJECT_DIRECTORY.getValue());
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
		
//		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.getFilename(proj.ORIGINAL_CENTROIDS_FILENAME));
//		Files.backup(proj.getFilename(proj.CUSTOM_CENTROIDS_FILENAME), "", "");
//		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.getFilename(proj.CUSTOM_CENTROIDS_FILENAME));
		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.ORIGINAL_CENTROIDS_FILENAME.getValue());
		Files.backup(proj.CUSTOM_CENTROIDS_FILENAME.getValue(), "", "");
		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.CUSTOM_CENTROIDS_FILENAME.getValue());
	}

	public static void parseCentroidsFromGenotypes(Project proj, boolean[] samplesToBeUsed, double missingnessThreshold) {
		String[] samples, markerNames;
		float[][][] centroids;
		float[][] centroid;
		int count;
		float[] thetas, rs;
		double[] meanThetas, meanRs;
		byte[] genotypes;
		int[] counts;
		SampleList sampleList;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		long time;
		MarkerSet markerSet;
		
		time = new Date().getTime();
		System.out.println("Computing centroids from genotype means");
		sampleList = proj.getSampleList();
		samples = sampleList.getSamples();
		if (samples.length != samplesToBeUsed.length) {
			System.err.println("Error - mismatched number of samples in project versus sample mask");
			return;
		}
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		centroids = new float[markerNames.length][][];
		
		markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
		time = new Date().getTime();
		for (int i = 0; i < markerNames.length; i++) {
			markerData = markerDataLoader.requestMarkerData(i);
			
			genotypes = markerData.getAbGenotypes();
			thetas = markerData.getThetas();
			rs = markerData.getRs();
			meanThetas = new double[5];
			meanRs = new double[5];
			counts = new int[5];
		
			for (int k = 0; k<samples.length; k++) {
				if (samplesToBeUsed[k] && !Float.isNaN(thetas[k]) && !Float.isNaN(rs[k])) {
					meanThetas[0] += thetas[k];
					meanRs[0] += rs[k];
					counts[0]++;
					meanThetas[genotypes[k]+2] += thetas[k];
					meanRs[genotypes[k]+2] += rs[k];
					counts[genotypes[k]+2]++;
				}
            }
			for (int k = 0; k<5; k++) {
				meanThetas[k] /= counts[k];
				meanRs[k] /= counts[k];
            }
			centroid = new float[3][];
			if (counts[1] >= counts[0]*missingnessThreshold) {
				for (int k = 0; k<3; k++) {
					centroid[k] = new float[] { (float)meanThetas[0], (float)meanRs[0] };  
                }
			} else {
				for (int k = 0; k<3; k++) {
					if (counts[k+2] > 0) {
						centroid[k] = new float[] { (float)meanThetas[k+2], (float)meanRs[k+2] };
					} else {
						centroid[k] = null;
					}
                }
			}
			
			centroids[i] = centroid;
			markerDataLoader.releaseIndex(i);
		}

		count = 0;
		for (int i = 0; i<centroids.length; i++) {
			if (centroids[i] == null) {
				if (count == 0) {
					System.out.println("The following marker(s) could not be computed:");
				}
				System.out.println(markerNames[i]);
				count++;
			}
        }
		if (count > 0) {
			System.out.println("Computed mean genotyped centroids for "+(centroids.length-count)+" of "+centroids.length+" markers, "+count+" missing");
		} else {
			System.out.println("Computed mean genotyped centroids for all "+centroids.length+" markers");
		}
//		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.getFilename(proj.GENOTYPE_CENTROIDS_FILENAME));
		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.GENOTYPE_CENTROIDS_FILENAME.getValue(true, false));
		System.out.println("Computation took "+ext.getTimeElapsed(time));
	}

	public static void generateChimeraCentroids(Project proj, String intensityOnlyFlagFile) {
        float[][][] cents, clusteredCents, unclusteredCents;
        Centroids clustered, unclustered;
        Hashtable<String,String> hash;
        String[] markerNames;
        MarkerSet markerSet;
        boolean problem, jar;
        String flag;
        long time;

        jar = proj.JAR_STATUS.getValue();
        time = new Date().getTime();
        markerSet = proj.getMarkerSet();        
        markerNames = markerSet.getMarkerNames();
        hash = HashVec.loadFileToHashString(proj.PROJECT_DIRECTORY.getValue()+intensityOnlyFlagFile, false);
//        if (!Files.exists(proj.getFilename(proj.GENOTYPE_CENTROIDS_FILENAME), jar)) {
//        	System.err.println("Error - file '"+proj.getFilename(proj.GENOTYPE_CENTROIDS_FILENAME)+"' does not exist in the project's data directory");
        if (!Files.exists(proj.GENOTYPE_CENTROIDS_FILENAME.getValue(), jar)) {
        	System.err.println("Error - file '"+proj.GENOTYPE_CENTROIDS_FILENAME.getValue()+"' does not exist in the project's data directory");
        	return;
        }
        if (!Files.exists(proj.ORIGINAL_CENTROIDS_FILENAME.getValue(), jar)) {
        	System.err.println("Error - file '"+proj.ORIGINAL_CENTROIDS_FILENAME.getValue()+"' does not exist in the project's data directory");
        	return;
        }
       	clustered = Centroids.load(proj.GENOTYPE_CENTROIDS_FILENAME.getValue(), jar);
       	unclustered = Centroids.load(proj.ORIGINAL_CENTROIDS_FILENAME.getValue(), jar);
        if (clustered.getFingerprint() != unclustered.getFingerprint()) {
        	System.err.println("Error - the two centroid files cannot be merged, because they do not have the same fingerprint");
        	return;
        }
        if (clustered.getFingerprint() != markerSet.getFingerprint()) {
        	System.err.println("Error - the centroid files do not match the fingerprint of the current marker set");
        	return;
        }
        
        problem = false;
        cents = new float[markerNames.length][][];
        clusteredCents = clustered.getCentroids();
        unclusteredCents = unclustered.getCentroids();
        for (int i = 0; i<markerNames.length; i++) {
        	flag = hash.get(markerNames[i]);
        	if (flag == null) {
        		System.err.println("Error - no flag for marker '"+markerNames[i]+"'");
        		problem = true;
        	} else if (flag.equals("1")) {
        		cents[i] = unclusteredCents[i];
        	} else if (flag.equals("0")) {
        		cents[i] = clusteredCents[i];
        	} else {
        		System.err.println("Error - invalid flag for marker '"+markerNames[i]+"'; found '"+flag+"', expecting '1' or '0'");
        		problem = true;
        	}
        }
        
        if (problem) {
        	System.err.println("Error - chimera centroids generation failed");
        } else {
//        	new Centroids(cents, markerSet.getFingerprint()).serialize(proj.getFilename(proj.CHIMERA_CENTROIDS_FILENAME));
    		new Centroids(cents, markerSet.getFingerprint()).serialize(proj.CHIMERA_CENTROIDS_FILENAME.getValue());
        	System.out.println("Created chimera centroids in "+ext.getTimeElapsed(time));
        }
	}
	
	public static void recompute(Project proj, String centroidsFile) {
		recompute(proj, centroidsFile, false, 1);
	}

	/**
	 *Recompute thread, applies centroids and saves new sample file
	 *
	 */
	private static class RecomputeWorker implements Callable<Hashtable<String, Float>>{
		private Project proj;
		private String sample; 
		private Centroids centroids;
		private boolean preserveBafs;
		
		
		public RecomputeWorker(Project proj, String sample, Centroids centroids, boolean preserveBafs) {
			super();
			this.proj = proj;
			this.sample = sample;
			this.centroids = centroids;
			this.preserveBafs = preserveBafs;
		}

		@Override
		public Hashtable<String, Float> call() throws Exception {
			Hashtable<String, Float> outliers = new Hashtable<String, Float>();
			Sample original = proj.getFullSampleFromRandomAccessFile(sample);
			Sample sample = new Sample(original.getSampleName(), original.getFingerprint(), original.getGCs(), original.getXs(), original.getYs(), preserveBafs ? original.getBAFs() : original.getBAFs(centroids.getCentroids()), original.getLRRs(centroids.getCentroids()), original.getForwardGenotypes(), original.getAB_Genotypes(), original.getCanXYBeNegative());
			sample.saveToRandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(false, true) + original.getSampleName() + Sample.SAMPLE_FILE_EXTENSION, outliers, sample.getSampleName());
			return outliers;
		}

	}

	/**
	 * Manages centroid application to each sample
	 *
	 */
	private static class RecomputeProducer implements Producer<Hashtable<String, Float>> {
		private Project proj;
		private String[] samples;
		private Centroids centroids;
		private boolean preserveBafs;
		private int index;

		public RecomputeProducer(Project proj, String[] samples, Centroids centroids, boolean preserveBafs) {
			super();
			this.proj = proj;
			this.samples = samples;
			this.centroids = centroids;
			this.preserveBafs = preserveBafs;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<Hashtable<String, Float>> next() {
			String currentSample = samples[index];
			RecomputeWorker worker = new RecomputeWorker(proj, currentSample, centroids, preserveBafs);
			index++;
			return worker;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
		
	}
	
	
	/**
	 * @param proj
	 * @param centroidsFile
	 * @param preserveBafs
	 *            bafs will not be recomputed from the centroids, useful if an atypical baf value is used
	 */
	public static void recompute(Project proj, String centroidsFile, boolean preserveBafs,int numThreads) {
		MarkerSet markerSet;
		Centroids centroids;
		// Sample original, sample;
		String[] samples;
		// float[][][] cents;

		markerSet = proj.getMarkerSet();
		centroids = load(centroidsFile, proj.JAR_STATUS.getValue());
		if (centroids.getFingerprint() != markerSet.getFingerprint()) {
			System.err.println("Error - fingerprint for Centroids file '" + centroidsFile + "' does not match the fingerprint for the current MarkerSet");
		}

		// cents = centroids.getCentroids();
		samples = proj.getSamples();
		Hashtable<String, Float> outliers = new Hashtable<String, Float>();
		RecomputeProducer producer = new RecomputeProducer(proj, samples, centroids, preserveBafs);
		WorkerTrain<Hashtable<String, Float>> train = new WorkerTrain<Hashtable<String, Float>>(producer, numThreads, 10, proj.getLog());
		while (train.hasNext()) {
			Hashtable<String, Float> currentOutliers = train.next();
			outliers.putAll(currentOutliers);
		}
		// for (int i = 0; i < samples.length; i++) {
		// original = proj.getFullSampleFromRandomAccessFile(samples[i]);
		// sample = new Sample(original.getSampleName(), original.getFingerprint(), original.getGCs(), original.getXs(), original.getYs(), preserveBafs ? original.getBAFs() : original.getBAFs(cents), original.getLRRs(cents), original.getForwardGenotypes(), original.getAB_Genotypes(), original.getCanXYBeNegative());
		// sample.saveToRandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(false, true) + original.getSampleName() + Sample.SAMPLE_DATA_FILE_EXTENSION, outliers, sample.getSampleName());
		// }
		if (outliers.size() > 0) {
			if (Files.exists(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser")) {
				Files.copyFile(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser", ext.addToRoot(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser", ext.getTimestampForFilename()));
			}
		}
		SerializedFiles.writeSerial(outliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");

	}

	public static float[][] computeClusterCenters(MarkerData markerData, boolean[] samplesToBeUsed, double missingnessThreshold) {
		float[][] centers;
		float[] xs, ys;
		double[] meanXs, meanYs;
		byte[] genotypes;
		int[] counts;
		
		centers = new float[3][2];
		genotypes = markerData.getAbGenotypes();
		xs = markerData.getXs();
		ys = markerData.getYs();
		meanXs = new double[5];
		meanYs = new double[5];
		counts = new int[5];

		for (int k = 0; k<xs.length; k++) {
			if ((samplesToBeUsed == null || samplesToBeUsed[k]) && !Float.isNaN(xs[k]) && !Float.isNaN(ys[k])) {
				meanXs[0] += xs[k];
				meanYs[0] += ys[k];
				counts[0]++;
				meanXs[genotypes[k]+2] += xs[k];
				meanYs[genotypes[k]+2] += ys[k];
				counts[genotypes[k]+2]++;
			}
        }
		for (int k = 0; k<5; k++) {
			meanXs[k] /= counts[k];
			meanYs[k] /= counts[k];
        }
		if (counts[1] >= counts[0]*missingnessThreshold) {
			for (int k = 0; k<3; k++) {
//				centers[k] = new float[] { (float)meanXs[0], (float)meanYs[0] };  
				centers[k] = null;  
            }
		} else {
			for (int k = 0; k<3; k++) {
				if (counts[k+2] > 0) {
					centers[k] = new float[] { (float)meanXs[k+2], (float)meanYs[k+2] };
				} else {
					centers[k] = null;
				}
            }
		}
		
		return centers;
	}
	
	
	
	/**
	 * 
	 * @param proj
	 * @param centFilename File path FROM THE PROJECT'S DIRECTORY
	 * @param exportFilename File path FROM THE PROJECT'S DIRECTORY
	 */
	public static void exportToText(Project proj, String centFilename, String exportFilename) {

		Centroids centObject;
		String dir;
		
		dir = proj.PROJECT_DIRECTORY.getValue();
		centObject = Centroids.load(dir+centFilename, false);
		
		centObject.exportToText(proj, dir + exportFilename);
	}
	
	@SuppressWarnings("unchecked")
    public static Centroids[] computeSexSpecificCentroids(final Project proj, final boolean[] includeList, String[] pfbFiles, String[] centFiles, final boolean shiftPFBForSex, int threadCount) {
        PrintWriter writerM, writerF;
        MarkerSet markerSet;
        String sampleDataFile;
        String[] allMarkers, samples, header;//, markersToUse;
        byte[] markerChrs;
        final boolean[] inclSampAll;
        final boolean[] inclSampFemales;
        final boolean[] inclSampMales;
        final int markerCount = Array.booleanArraySum(includeList);
        int[] sampleSex;
        final float[][][] rawCentroidsFemale;
        final float[][][] rawCentroidsMale;
        Vector<String>[] markerLists;
        final Logger log = proj.getLog();
        ExecutorService computeHub;
        final ConcurrentLinkedQueue<Integer>[] markerIndexQueues;
        final Hashtable<Integer, String[][]> pfbInfo;
        final Hashtable<Integer, Integer>[] fullToTruncMarkerIndices;
        Hashtable<String, Vector<String>> sexData;
        final MarkerDataLoader[] markerDataLoaders;
        SampleData sampleData;
        
        markerSet = proj.getMarkerSet();
        sampleData = proj.getSampleData(0, false);


        inclSampAll = proj.getSamplesToInclude(null);
        if (!sampleData.hasExcludedIndividuals()) {
            log.report("Warning - there is no 'Exclude' column in SampleData.txt; centroids will be determined using all samples.");
        }
        samples = proj.getSamples();
        sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
        header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
        int sexInd = -1;
        for (int i = 0; i < header.length; i++) {
            if (("CLASS=" + SexChecks.EST_SEX_HEADER).toUpperCase().equals(header[i].toUpperCase())) {
                sexInd = i;
                break;
            }
        }
        if (sexInd == -1) {
            log.reportError("Error - no estimated sex found in sample data file - please run SexChecks with -check argument to generate the required data");
            return null;
        }
        sexData = HashVec.loadFileToHashVec(sampleDataFile, 0, new int[] { sexInd }, "\t", true, false);
        
        inclSampMales = Array.clone(inclSampAll);
        inclSampFemales = Array.clone(inclSampAll);
        sampleSex = new int[inclSampAll.length];
        for (int i = 0; i < samples.length; i++) {
            int sex = sampleData.getSexForIndividual(samples[i]);
            if (sex == -1) {
                sex = Integer.parseInt(sexData.get(samples[i].toUpperCase()).get(0));
            }
            sampleSex[i] = sex;
            if (sex == 1) {
                inclSampFemales[i] = false;
            } else if (sex == 2) {
                inclSampMales[i] = false;
            } else {
                // Leave these for now, but when computing LRRs and BAFs, will need to be crafty....
            }
        }

        if (threadCount == -1) {
            threadCount = Runtime.getRuntime().availableProcessors();
        }
        markerIndexQueues = new ConcurrentLinkedQueue[threadCount];
        markerLists = new Vector[threadCount];
        fullToTruncMarkerIndices = new Hashtable[threadCount];
        markerDataLoaders = new MarkerDataLoader[threadCount];
        for (int i = 0; i < threadCount; i++) {
            markerLists[i] = new Vector<String>();
            markerIndexQueues[i] = new ConcurrentLinkedQueue<Integer>();
            fullToTruncMarkerIndices[i] = new Hashtable<Integer, Integer>();
        }
        
        allMarkers = markerSet.getMarkerNames();
        if (includeList.length != allMarkers.length) {
            log.reportError("Error - mismatched lists; list of markers to include in centroid computation must be the same size as the marker list for the given project"); 
        }
        markerChrs = markerSet.getChrs();
        
        int qInd = 0;
        for (int i = 0; i < markerChrs.length; i++) {
            if (includeList[i]) {
                markerLists[qInd].add(allMarkers[i]);
                fullToTruncMarkerIndices[qInd].put(i, markerLists[qInd].size() - 1);
                markerIndexQueues[qInd].add(Integer.valueOf(i));
                qInd = (qInd + 1) % threadCount;
            }
        }
        
        for (int i = 0; i < threadCount; i++) {
            markerDataLoaders[i] = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, Array.toStringArray(markerLists[i]));
        }
        
        rawCentroidsMale = new float[allMarkers.length][][];
        rawCentroidsFemale = new float[allMarkers.length][][];

        pfbInfo = new Hashtable<Integer, String[][]>();
        log.report("Computing sex-specific centroids for " + markerCount + " sex-specific markers on " + threadCount + " thread(s).");
        
        computeHub = Executors.newFixedThreadPool(threadCount);
        
        
        for (int i = 0; i < threadCount; i++) {
            final int myIndex = i;
            final long myStartTime = System.currentTimeMillis();
            computeHub.execute(new Runnable() {
                @Override
                public void run() {
                    int myMarkerCount = 0;
                    while(!markerIndexQueues[myIndex].isEmpty()) {
                        Integer indexInt = markerIndexQueues[myIndex].poll();
                        if (indexInt == null) continue;
                        int index = indexInt.intValue();
                        
                        if (!includeList[index]) {
                            rawCentroidsMale[index] = new float[][]{{Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}};
                            rawCentroidsFemale[index] = new float[][]{{Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}};
                            continue;
                        }
                        
                        int markerIndex = fullToTruncMarkerIndices[myIndex].get(index);
                        MarkerData markerData = markerDataLoaders[myIndex].requestMarkerData(markerIndex);
                        CentroidCompute centCompM = new CentroidCompute(markerData, 
                                                    null, 
                                                    inclSampMales, 
                                                    false, // NOT intensity only 
                                                    1, // no filtering
                                                    0,  // no filtering
                                                    null,  // no filtering
                                                    true,  // median, not mean
                                                    proj.getLog());
                        
                        CentroidCompute centCompF = new CentroidCompute(markerData, 
                                                    null, 
                                                    inclSampFemales, 
                                                    false, // NOT intensity only 
                                                    1, // no filtering
                                                    0,  // no filtering
                                                    null,  // no filtering
                                                    true,  // median, not mean
                                                    proj.getLog());
                        
                        
                        centCompM.computeCentroid(true);
                        centCompF.computeCentroid(true);
                        
                        rawCentroidsMale[index] = centCompM.getCentroid();
                        rawCentroidsFemale[index] = centCompF.getCentroid();
                        
                        float[] bafCnt = new float[]{0, 0};
                        float[] bafSum = new float[]{0, 0};
                        float[] genCnt = new float[]{0, 0};
                        float[] bafM = centCompM.getRecomputedBAF();
                        byte[] genM = centCompM.getClustGenotypes();
                        float[] bafF = centCompF.getRecomputedBAF();
                        byte[] genF = centCompF.getClustGenotypes();
                        for (int s = 0; s < inclSampAll.length; s++) {
                            if (inclSampMales[s]) {
                                if (!Float.isNaN(bafM[s])) {
                                    bafSum[0] += bafM[s];
                                    bafCnt[0]++;
                                    if (genM[s] >= 0) {
                                        genCnt[0]++;
                                    }
                                }
                            }
                            if (inclSampFemales[s]) {
                                if (!Float.isNaN(bafF[s])) {
                                    bafSum[1] += bafF[s];
                                    bafCnt[1]++;
                                    if (genF[s] >= 0) {
                                        genCnt[1]++;
                                    }
                                }
                            }
                        }
                        
                        pfbInfo.put(index, new String[][]{
                                {markerData.getMarkerName(), "" + (shiftPFBForSex ? markerData.getChr() - 22 : markerData.getChr()), "" + markerData.getPosition(), "" + (genCnt[0] > 0 ? (bafSum[0] / bafCnt[0]) : 2)},
                                {markerData.getMarkerName(), "" + (shiftPFBForSex ? markerData.getChr() - 22 : markerData.getChr()), "" + markerData.getPosition(), "" + (genCnt[1] > 0 ? (bafSum[1] / bafCnt[1]) : 2)}
                        });
                        if (markerIndex > 0 && markerIndex % 10000 == 0) {
                            log.report(ext.getTime() + "\t...sex centroids computed up to marker " + (markerCount - markerIndex) + " of " + markerCount);
                        }
                        
                        markerDataLoaders[myIndex].releaseIndex(markerIndex);
                        centCompM = null;
                        centCompF = null;
                        
                        myMarkerCount++;
                    }
                    
                    System.out.println("Thread " + myIndex + " processed " + myMarkerCount + " markers in " + ext.getTimeElapsed(myStartTime));
                }
            });
        }
        
        computeHub.shutdown();
        try {
            computeHub.awaitTermination(Long.MAX_VALUE, java.util.concurrent.TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            log.report("Centroid computation was interrupted - .pfb and .cent files may not be complete or correct.");
        }
        computeHub = null;
        
        int nullCnt = 0;
        for (int i = 0; i < rawCentroidsFemale.length; i++) {
            if (rawCentroidsFemale[i] == null) {
                nullCnt++;
            }
        }
        System.out.println(nullCnt + " null cent entries");
        
        if (pfbFiles != null) {
            log.report("Writing sex-specific PFB files");
                
            try {
                writerM = new PrintWriter(new FileWriter(pfbFiles[0]));
                writerF = new PrintWriter(new FileWriter(pfbFiles[1]));
                
                writerM.println("Name\tChr\tPosition\tPFB");
                writerF.println("Name\tChr\tPosition\tPFB");
                
                for (int i = 0; i < allMarkers.length; i++) {
                    if (!includeList[i]) continue;
                    String[][] pfbEntry = pfbInfo.get(Integer.valueOf(i));
                    writerM.println(pfbEntry[0][0] + "\t" + pfbEntry[0][1] + "\t" + pfbEntry[0][2] + "\t" + pfbEntry[0][3]);
                    writerF.println(pfbEntry[1][0] + "\t" + pfbEntry[1][1] + "\t" + pfbEntry[1][2] + "\t" + pfbEntry[1][3]);
                }
                
                writerM.flush();
                writerF.flush();

                writerM.close();
                writerF.close();
            } catch (IOException e1) {
                log.reportError("Error - problem occured when writing to new sex-specific .pfb files");
                log.reportException(e1);
            }
            
            writerM = null;
            writerF = null;
        }
        
        Centroids[] centroids = new Centroids[2]; 
        centroids[0] = new Centroids(rawCentroidsMale, markerSet.getFingerprint());
        centroids[1] = new Centroids(rawCentroidsFemale, markerSet.getFingerprint());
        
        if (centFiles != null) {
            log.report("Writing sex-specific Centroid files");
            
            centroids[0].serialize(centFiles[0]);
            Centroids.exportToText(proj, centFiles[0], centFiles[0] + ".txt", allMarkers);
            proj.SEX_CENTROIDS_MALE_FILENAME.setValue(centFiles[0]);
            
            centroids[1].serialize(centFiles[1]);
            Centroids.exportToText(proj, centFiles[1], centFiles[1] + ".txt", allMarkers);
            proj.SEX_CENTROIDS_FEMALE_FILENAME.setValue(centFiles[1]);

            proj.saveProperties(new Project.Property[]{proj.SEX_CENTROIDS_MALE_FILENAME, proj.SEX_CENTROIDS_FEMALE_FILENAME});
        }
        
        return centroids;
    }
	
	public static void exportToText(Project proj, String centFilename, String exportFilename, String[] markerNames) {
		PrintWriter writer;
		Centroids centObject;
		float[][][] centroids;
		String dir;
		
		dir = proj.PROJECT_DIRECTORY.getValue();
		String file = centFilename.startsWith(dir) || centFilename.contains(":") || centFilename.startsWith("/") ? centFilename : dir + centFilename;
		centObject = Centroids.load(file, false);
		centroids = centObject.getCentroids();
		
		if (markerNames.length != centroids.length) {
			System.err.println("Error - mismatched number of markers in the project's marker set and the imported centroids file ("+centFilename+"); aborting");
			return;
		}
		
		if (MarkerSet.fingerprint(markerNames) != centObject.getFingerprint()) {
			System.err.println("Error - mismatched marker fingerprints in the project's marker set and the imported centroids file ("+centFilename+"); aborting");
			return;
		}
		
		String outFile = exportFilename.startsWith(dir) || exportFilename.contains(":") || exportFilename.startsWith("/") ? exportFilename : dir + exportFilename;
		try {
			writer = new PrintWriter(new FileWriter(outFile));
			writer.println("marker_fingerprint="+centObject.getFingerprint());
			writer.println("MarkerName\tAA_Theta_Mean\tAA_R_Mean\tAB_Theta_Mean\tAB_R_Mean\tBB_Theta_Mean\tBB_R_Mean");
			for (int i = 0; i < markerNames.length; i++) {
				writer.print(markerNames[i]);
				for (int j = 0; j < 3; j++) {
					if (centroids[i] == null || centroids[i][j] == null) {
						writer.print("\t.\t.");
					} else {
						writer.print("\t"+centroids[i][j][0]+"\t"+centroids[i][j][1]);
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + outFile);
			e.printStackTrace();
		}
	}
	
	public static void importFromText(Project proj, String importFilename, String centFilename) {

	}
		
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
//		String centFile = "ForNathan_table.csv";
//		String centFile = "SNP Table Myers raw dataset final 022908.csv";
//		String centFile = "Myers_final_042208_ReclusteredCNV_SNP_Table2.csv";
//		String centFile = "Theta_R_mean_dev_550.csv";
		String centFile = "CentroidExample.csv";
		String intensityFlags = "";
//		String intensityFlags = "intesityFlag.txt";
//		String clusteredCentroids = proj.DEFAULT_ORIGINAL_CENTROIDS_FILENAME;
//		String unclusteredCentroids = proj.DEFAULT_GENOTYPE_CENTROIDS_FILENAME;
		boolean fromGenotypes = false;
		boolean sexSpecific = false;
		Project proj;
		String compute = "";
		String importFile = null;
		String exportFile = null;
		int numThreads = 1;

		String usage = "\n"+
			"cnv.filesys.Centroids requires 0-1 arguments\n"+
			"   (1) project properties filename (i.e. proj="+org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
			"   (2) filename (i.e. file=" + centFile + " (default))\n"+
			" OR\n"+
			"   (2) generate centroids from genotypes (i.e. -fromGenotypes (not the default))\n"+
			" OR\n"+
			"   (2) file with intensity only flags (i.e. flags=intensityFlags.dat (not the default))\n"+
			"   (3) centroid file for clustered markers (see \"GENOTYPE_CENTROIDS_FILENAME\" in the Project properties file)\n"+
			"   (4) centroid file for intensity only markers (see \"GENOTYPE_CENTROIDS_FILENAME\" in the Project properties file)\n"+
			" OR\n"+
			"   (2) recompute BAF/LRR and generate new Sample files using these centroids (i.e. compute=genotype.cent (not the default))\n"+
			"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help")
					|| args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("file=")) {
				centFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("import=")) {
				importFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("export=")) {
				exportFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-fromGenotypes")) {
				fromGenotypes = true;
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
			    numThreads = ext.parseIntArg(args[i]);
			    numArgs--;
			} else if (args[i].startsWith("flags=")) {
				intensityFlags = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("compute=")) {
				compute = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-sexSpecific")) {
			    sexSpecific = true;
			    numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		proj = new Project(filename, false);
//		fromGenotypes = true;
////		compute = "genotype.cent";
//		
//		centFile = "data/genotype.cent";
//		exportFile = "data/genotype.cent.xln";
		try {
			if (exportFile != null) {
				exportToText(proj, centFile, exportFile);
			} else if (sexSpecific) {
			    String pennData, sexDir, malePFB, femalePFB, centFilePathM, centFilePathF;
	            pennData = proj.getProperty(proj.PENNCNV_DATA_DIRECTORY);
	            sexDir = pennData + "sexSpecific/";
	            malePFB = sexDir + "males.pfb";
	            femalePFB = sexDir + "females.pfb";
	            centFilePathM = sexDir + "sexSpecific_Male.cent";
	            centFilePathF = sexDir + "sexSpecific_Female.cent";
	            computeSexSpecificCentroids(proj, AnalysisFormats.getChromosomalMarkersOnly(proj), new String[]{malePFB, femalePFB}, new String[]{centFilePathM, centFilePathF}, true, numThreads);
			} else if (importFile != null) {
				importFromText(proj, importFile, centFile);
			} else if (fromGenotypes) {
				parseCentroidsFromGenotypes(proj, Array.booleanArray(proj.getSamples().length, true), 1);
			} else if (!compute.equals("")) {
				recompute(proj, compute, false, numThreads);
			} else if (!intensityFlags.equals("")) {
				generateChimeraCentroids(proj, intensityFlags);
			} else { 
				parseIlluminaCentroidsFromCSV(proj, centFile);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
