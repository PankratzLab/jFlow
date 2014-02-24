package cnv.filesys;

import java.io.*;
import java.util.*;

import cnv.manage.MarkerDataLoader;

import stats.Maths;
import common.*;

public class Centroids implements Serializable {
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
		Files.writeSerial(this, filename);
	}

	public static Centroids load(String filename, boolean jar) {
		return (Centroids)Files.readSerial(filename, jar, true);
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
			reader = new BufferedReader(new FileReader(proj.getProjectDir()+filename));
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
			System.err.println("Error: file \"" + filename + "\" not found in "+proj.getProjectDir());
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
		
		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.getFilename(Project.ORIGINAL_CENTROIDS_FILENAME));
		Files.backup(proj.getFilename(Project.CUSTOM_CENTROIDS_FILENAME), "", "");
		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.getFilename(Project.CUSTOM_CENTROIDS_FILENAME));
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
		
		markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames, new Logger());
		time = new Date().getTime();
		for (int i = 0; i < markerNames.length; i++) {
			markerData = markerDataLoader.requestMarkerData(i);
			
			genotypes = markerData.getAB_Genotypes();
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
		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.getFilename(Project.GENOTYPE_CENTROIDS_FILENAME));
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

        jar = proj.getJarStatus();
        time = new Date().getTime();
        markerSet = proj.getMarkerSet();        
        markerNames = markerSet.getMarkerNames();
        hash = HashVec.loadFileToHashString(proj.getProjectDir()+intensityOnlyFlagFile, false);
        if (!Files.exists(proj.getFilename(Project.GENOTYPE_CENTROIDS_FILENAME), jar)) {
        	System.err.println("Error - file '"+proj.getFilename(Project.GENOTYPE_CENTROIDS_FILENAME)+"' does not exist in the project's data directory");
        	return;
        }
        if (!Files.exists(proj.getFilename(Project.ORIGINAL_CENTROIDS_FILENAME), jar)) {
        	System.err.println("Error - file '"+proj.getFilename(Project.ORIGINAL_CENTROIDS_FILENAME)+"' does not exist in the project's data directory");
        	return;
        }
       	clustered = Centroids.load(proj.getFilename(Project.GENOTYPE_CENTROIDS_FILENAME), jar);
       	unclustered = Centroids.load(proj.getFilename(Project.ORIGINAL_CENTROIDS_FILENAME), jar);
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
    		new Centroids(cents, markerSet.getFingerprint()).serialize(proj.getFilename(Project.CHIMERA_CENTROIDS_FILENAME));
        	System.out.println("Created chimera centroids in "+ext.getTimeElapsed(time));
        }
	}
	
	public static void recompute(Project proj, String centroidsFile) {
		MarkerSet markerSet;
		Centroids centroids;
        Sample original, sample;
        String[] samples;
        float[][][] cents;
        
		markerSet = proj.getMarkerSet();
		centroids = load(centroidsFile, proj.getJarStatus());
		if (centroids.getFingerprint() != markerSet.getFingerprint()) {
			System.err.println("Error - fingerprint for Centroids file '"+centroidsFile+"' does not match the fingerprint for the current MarkerSet");
		}

        cents = centroids.getCentroids(); 
        samples = proj.getSamples();
        for (int i = 0; i<samples.length; i++) {
        	original = proj.getFullSampleFromRandomAccessFile(samples[i]);
        	sample = new Sample(original.getSampleName(), original.getFingerprint(), original.getGCs(), original.getXs(), original.getYs(), original.getBAFs(cents), original.getLRRs(cents), original.getForwardGenotypes(), original.getAB_Genotypes(), original.getCanXYBeNegative());
        	sample.saveToRandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY) + original.getSampleName() + Sample.SAMPLE_DATA_FILE_EXTENSION);
        }
	}

	public static float[][] computeClusterCenters(MarkerData markerData, boolean[] samplesToBeUsed, double missingnessThreshold) {
		float[][] centers;
		float[] xs, ys;
		double[] meanXs, meanYs;
		byte[] genotypes;
		int[] counts;
		
		centers = new float[3][2];
		genotypes = markerData.getAB_Genotypes();
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
	
	public static void exportToText(Project proj, String centFilename, String exportFilename) {
		PrintWriter writer;
		Centroids centObject;
		float[][][] centroids;
		MarkerSet markerSet;
		String[] markerNames;
		String dir;
		
		dir = proj.getProjectDir();
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		centObject = Centroids.load(dir+centFilename, false);
		centroids = centObject.getCentroids();
		
		if (markerNames.length != centroids.length) {
			System.err.println("Error - mismatched number of markers in the project's marker set and the imported centroids file ("+centFilename+"); aborting");
			return;
		}

		if (markerSet.getFingerprint() != centObject.getFingerprint()) {
			System.err.println("Error - mismatched marker fingerprints in the project's marker set and the imported centroids file ("+centFilename+"); aborting");
			return;
		}
		
		try {
			writer = new PrintWriter(new FileWriter(dir+exportFilename));
			writer.println("marker_fingerprint="+centObject.getFingerprint());
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
			System.err.println("Error writing to " + exportFilename);
			e.printStackTrace();
		}
	}

	public static void importFromText(Project proj, String importFilename, String centFilename) {

	}
		
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;
//		String centFile = "ForNathan_table.csv";
//		String centFile = "SNP Table Myers raw dataset final 022908.csv";
//		String centFile = "Myers_final_042208_ReclusteredCNV_SNP_Table2.csv";
//		String centFile = "Theta_R_mean_dev_550.csv";
		String centFile = "CentroidExample.csv";
		String intensityFlags = "";
//		String intensityFlags = "intesityFlag.txt";
//		String clusteredCentroids = Project.DEFAULT_ORIGINAL_CENTROIDS_FILENAME;
//		String unclusteredCentroids = Project.DEFAULT_GENOTYPE_CENTROIDS_FILENAME;
		boolean fromGenotypes = false;
		Project proj;
		String compute = "";
		String importFile = null;
		String exportFile = null;

		String usage = "\n"+
			"cnv.filesys.Centroids requires 0-1 arguments\n"+
			"   (1) project (i.e. proj=" + filename + " (default))\n"+
			"   (2) filename (i.e. file=" + centFile + " (default))\n"+
			" OR\n"+
			"   (2) generate centroids from genotypes (i.e. -fromGenotypes (not the default))\n"+
			" OR\n"+
			"   (2) file with intensity only flags (i.e. flags=intensityFlags.dat (not the default))\n"+
			"   (3) centroid file for clustered markers (see " + Project.GENOTYPE_CENTROIDS_FILENAME + " in the Project properties file)\n"+
			"   (4) centroid file for intensity only markers (see " + Project.GENOTYPE_CENTROIDS_FILENAME + " in the Project properties file)\n"+
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
			} else if (args[i].startsWith("flags=")) {
				intensityFlags = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("compute=")) {
				compute = args[i].split("=")[1];
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
			} else if (importFile != null) {
				importFromText(proj, importFile, centFile);
			} else if (fromGenotypes) {
				parseCentroidsFromGenotypes(proj, Array.booleanArray(proj.getSamples().length, true), 1);
			} else if (!compute.equals("")) {
				recompute(proj, compute);
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
