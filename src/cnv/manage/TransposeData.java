// -Xms1024M -Xmx1024M     or even better: -Xmx15g
package cnv.manage;

import java.io.*;
import java.util.*;

import cnv.filesys.Compression;
import cnv.filesys.Sample;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerDataCollection;
import cnv.filesys.MarkerLookup;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.SampleList;
import cnv.filesys.Centroids;
import common.*;

public class TransposeData {
	public static final int DEFAULT_MARKERS_PER_FILE = 210000;
	public static final String MARKDATA_FILE_EXTENSION = ".scatRaf";
	public static final byte MARKDATA_PARAMETER_TOTAL_LEN = 21;
	public static final byte MARKDATA_NUMSAMPS_START = 0;
	public static final byte MARKDATA_NUMSAMPS_LEN = 4;
	public static final byte MARKDATA_NUMMARKS_START = 4;
	public static final byte MARKDATA_NUMMARKS_LEN = 4;
	public static final byte MARKDATA_NULLSTATUS_START = 8;
	public static final byte MARKDATA_NULLSTATUS_LEN = 1;
	public static final byte MARKDATA_FINGERPRINT_START = 9;
	public static final byte MARKDATA_FINGERPRINT_LEN = 8;
	public static final byte MARKDATA_MARKNAMELEN_START = 17;
	public static final byte MARKDATA_MARKNAMELEN_LEN = 4;
	public static final byte MARKDATA_MARKNAME_START = 21;


	// private Project proj;
	// private int[][] indicesSets;
	// private long fingerprint;
	// private String[] outputFiles;
	//
	// public ExtractPlots(Project proj, int[][] indicesSets, long fingerprint,
	// String[] outputFiles) {
	// this.proj = proj
	// this.indicesSets = indicesSets;
	// this.outputFiles = outputFiles;
	// this.fingerprint = fingerprint;
	//		
	// if (outputFiles.length != indicesSets.length) {
	// System.err.println("Error - malformed inputs to ExtractPlots");
	// System.exit(1);
	// }
	// }

	// public void run() {
	/*
	public static void extractPlots(Project proj, int[][] indicesSets, long fingerprint, boolean compact_XYonly, String[] outputFiles) {
		for (int i = 0; i<indicesSets.length; i++) {
			sampleToSNPmajor(proj, indicesSets[i], fingerprint, compact_XYonly, outputFiles[i]);
		}
	}
	*/
	
	public static boolean containsLRR(Project proj) {
		return proj.getFullSampleFromRandomAccessFile(proj.getSampleList().getSamples()[0]).getLRRs() != null;
//		
//		boolean result;
//		//Extracts one first;
//		extract(proj, proj.getDir(Project.PLOT_DIRECTORY)+"Markers.test.scat", true);
//		//read the result of the extracted
//		MarkerDataCollection temp = MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+"Markers.test.scat", false);
//		//test to see if the result is null
//		if (temp.getCollection()[0].getLRRs()[0]==(float)0) {
//			result = true;
//		} else {
//			result = false;
//		}
//		(new File(proj.getDir(Project.PLOT_DIRECTORY)+"Markers.test.scat")).delete();
//		return result;			
	}
	

	public static void sampleToSNPmajor(Project proj, int[] indices, long fingerprint, boolean compact_XYonly, String outputFile) {
		SampleList list;
		MarkerSet markerSet;
		String[] samples, allMarkers, markerNames;
		Sample samp;
		MarkerData[] collection;
//		float[][] gcs, xs, ys, xRaws, yRaws, thetas, rs, lrrs, bafs;
		float[][] gcs, xs, ys, thetas, rs, lrrs, bafs;
		byte[][] abGenotypes;
//		float[] samp_gcs, samp_xs, samp_ys, samp_xRaws, samp_yRaws, samp_thetas, samp_rs, samp_lrrs, samp_bafs;
		float[] samp_gcs, samp_xs, samp_ys, samp_thetas, samp_rs, samp_lrrs, samp_bafs;
		byte[] samp_abGenotypes, samp_forwardGenotypes;
		String[][] alleleMappings;
		long time;
		byte[] allChrs, chrs;
		int[] allPositions, positions;

		list = proj.getSampleList();
		samples = list.getSamples();
		if (list.getFingerprint()!=fingerprint) {
			System.err.println("Error - SampleList mixup somewhere along the way...");
			System.exit(1);
		}

		markerSet = proj.getMarkerSet();
		allMarkers = markerSet.getMarkerNames();
		allChrs = markerSet.getChrs();
		allPositions = markerSet.getPositions();

		markerNames = new String[indices.length];
		chrs = new byte[indices.length];
		positions = new int[indices.length];
		for (int i = 0; i<markerNames.length; i++) {
			markerNames[i] = allMarkers[indices[i]];
			chrs[i] = allChrs[indices[i]];
			positions[i] = allPositions[indices[i]];
		}

		System.out.println(ext.getTime()+" "+ext.getDate()+" (serialized in memory; "+outputFile+" will have "+markerNames.length+" markers)");
		time = new Date().getTime();
		gcs = new float[markerNames.length][samples.length];
		xs = new float[markerNames.length][samples.length];
		ys = new float[markerNames.length][samples.length];
		abGenotypes = new byte[markerNames.length][samples.length];
		alleleMappings = new String[markerNames.length][4];
		if (compact_XYonly) {
//			xRaws = new float[markerNames.length][];
//			yRaws = new float[markerNames.length][];
			thetas = new float[markerNames.length][];
			rs = new float[markerNames.length][];
			lrrs = new float[markerNames.length][];
			bafs = new float[markerNames.length][];
		} else {
//			xRaws = new float[markerNames.length][samples.length];
//			yRaws = new float[markerNames.length][samples.length];
			thetas = new float[markerNames.length][samples.length];
			rs = new float[markerNames.length][samples.length];
			lrrs = new float[markerNames.length][samples.length];
			bafs = new float[markerNames.length][samples.length];
		}
		for (int i = 0; i<samples.length; i++) {
			samp = proj.getFullSampleFromRandomAccessFile(samples[i]);
			samp_gcs = samp.getGCs();
			samp_xs = samp.getXs();
			samp_ys = samp.getYs();
			samp_abGenotypes = samp.getAB_Genotypes();
			samp_forwardGenotypes = samp.getForwardGenotypes();
//			samp_xRaws = samp.getX_Raws();
//			samp_yRaws = samp.getY_Raws();
			samp_thetas = samp.getThetas();
			samp_rs = samp.getRs();
			samp_lrrs = samp.getLRRs();
			samp_bafs = samp.getBAFs();
			for (int m = 0; m<markerNames.length; m++) {
				gcs[m][i] = samp_gcs[indices[m]];
				xs[m][i] = samp_xs[indices[m]];
				ys[m][i] = samp_ys[indices[m]];
				if (!compact_XYonly) {
//					xRaws[m][i] = samp_xRaws[indices[m]];
//					yRaws[m][i] = samp_yRaws[indices[m]];
					thetas[m][i] = samp_thetas[indices[m]];
					rs[m][i] = samp_rs[indices[m]];
					lrrs[m][i] = samp_lrrs[indices[m]];
					bafs[m][i] = samp_bafs[indices[m]];
				}
				if (samp_abGenotypes == null) {
					abGenotypes[m][i] = -1;
				} else {
					abGenotypes[m][i] = samp_abGenotypes[indices[m]];
					if (alleleMappings[m][samp_abGenotypes[indices[m]]+1]==null) {
						alleleMappings[m][samp_abGenotypes[indices[m]]+1] = Sample.ALLELE_PAIRS[samp_forwardGenotypes[indices[m]]];
					} else if (!alleleMappings[m][samp_abGenotypes[indices[m]]+1].equals(Sample.ALLELE_PAIRS[samp_forwardGenotypes[indices[m]]])) {
						System.err.println("Error - mismatched alleles for marker "+markerNames[m]);
					}
				}
			}
		}
		collection = new MarkerData[markerNames.length];
		for (int i = 0; i<markerNames.length; i++) {
//			collection[i] = new MarkerData(markerNames[i], chrs[i], positions[i], fingerprint, gcs[i], xRaws[i], yRaws[i], xs[i], ys[i], thetas[i], rs[i], bafs[i], lrrs[i], abGenotypes[i], alleleMappings[i]);
			collection[i] = new MarkerData(markerNames[i], chrs[i], positions[i], fingerprint, gcs[i], null, null, xs[i], ys[i], thetas[i], rs[i], bafs[i], lrrs[i], abGenotypes[i], alleleMappings[i]);
		}
		new MarkerDataCollection(outputFile, fingerprint, markerNames, collection).serialize();

		System.out.println(ext.getTime()+" "+ext.getDate()+" (finished in "+ext.getTimeElapsed(time)+")");
	}

	public static String translateName(String str) {
		return ext.replaceAllWith(str, ":", ".");
	}

	// extracts just those markers listed in filename
	public static void extract(Project proj, String filename, boolean compact) {
		String[] targets;
		SampleList sampleList;
		MarkerSet markerSet;
		String[] markerNames;
		long fingerprint;
		int[] keys;
		Hashtable<String,String> hash;
		boolean prob = false;
		String trav;
		long time;

		sampleList = proj.getSampleList();
		if (sampleList==null) {
			sampleList = SampleList.generateSampleList(proj);
		}
		fingerprint = sampleList.getFingerprint();

		targets = HashVec.loadFileToStringArray(filename, false, false, new int[] {0}, true);
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();

		System.out.println("Assuming no MarkerLookup; extracting from .fsamp files");
		time = new Date().getTime();
		hash = new Hashtable<String,String>();
		for (int i = 0; i<targets.length; i++) {
			hash.put(targets[i], "");
		}
		for (int i = 0; i<markerNames.length; i++) {
			if (hash.containsKey(markerNames[i])) {
				hash.put(markerNames[i], i+"");
			}
		}
		System.out.println("Parsed markerIndices in "+ext.getTimeElapsed(time));
		keys = new int[targets.length];
		for (int i = 0; i<targets.length; i++) {
			trav = hash.get(targets[i]);
			if (trav.equals("")) {
				System.err.println("Error - '"+targets[i]+"' not found in the main database");
				prob = true;
			} else {
				keys[i] = Integer.parseInt(trav);
			}
		}
		if (prob) {
			System.exit(1);
		}
//		new File(proj.getDir(Project.PLOT_DIRECTORY)).mkdirs();
		sampleToSNPmajor(proj, keys, fingerprint, compact, proj.getProjectDir()+ext.rootOf(filename)+".scat");
		createMarkerLookup(proj);
	}

	public static void extractForDemo(Project proj, String filename) {
		String[] targets;
		SampleList sampleList;
		MarkerSet markerSet, demoSet;
		MarkerData[] markerData;
		String[] markerNames, files;
		long fingerprint, time;
		byte[] chrs, demoChrs;
		int[] positions, demoPositions;
		Hashtable<String,String> hash;
		boolean prob = false;
		int index;
		float[][][][] allCents, demoCents;
		Centroids cents;
		boolean jar;
		
		time = new Date().getTime();
		sampleList = proj.getSampleList();
		if (sampleList==null) {
			sampleList = SampleList.generateSampleList(proj);
		}
		fingerprint = sampleList.getFingerprint();
		
		targets = HashVec.loadFileToStringArray(filename, false, false, new int[] {0}, true);
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		
		jar = proj.getJarStatus();
		files = Files.list(proj.getDir(Project.DATA_DIRECTORY), ".cent", jar);
		System.out.println("Found "+files.length+" .cent files in "+proj.getDir(Project.DATA_DIRECTORY));
		allCents = new float[files.length][][][];
		for (int i = 0; i<files.length; i++) {
			cents = Centroids.load(proj.getDir(Project.DATA_DIRECTORY)+files[i], jar);
			if (cents.getFingerprint() != markerSet.getFingerprint()) {
				System.err.println("Error - Centroids file '"+files[i]+"' does not match up with the fingerprint of the current marker set");
				System.exit(1);
			}
			allCents[i] = cents.getCentroids();
        }
		
		demoChrs = new byte[targets.length];
		demoPositions = new int[targets.length];
		demoCents = new float[files.length][targets.length][][];
		hash = new Hashtable<String,String>();
		for (int i = 0; i<targets.length; i++) {
			index = ext.indexOfStr(targets[i], markerNames);
			if (index==-1) {
				System.err.println("Error - could not find "+targets[i]+" in the MarkerSet");
				prob = true;
			} else {
				demoChrs[i] = chrs[index];
				demoPositions[i] = positions[index];
				for (int j = 0; j<allCents.length; j++) {
					demoCents[j][i] = allCents[j][index];
                }
				hash.put(targets[i], ext.rootOf(filename)+".scat\t"+i);
			}
		}
		if (prob) {
			System.exit(1);
		}

		markerData = MarkerSet.loadFromList(proj, targets);
		new File(proj.getDir(Project.DEMO_DIRECTORY)+proj.getProperty(Project.DATA_DIRECTORY)).mkdirs();
		new File(proj.getDir(Project.DEMO_DIRECTORY)+proj.getProperty(Project.PLOT_DIRECTORY)).mkdirs();
//		new MarkerDataCollection(proj.getDir(Project.DEMO_DIRECTORY)+proj.getProperty(Project.PLOT_DIRECTORY)+ext.rootOf(filename)+".scat", fingerprint, markerNames, markerData).serialize();
		new MarkerDataCollection(proj.getDir(Project.DEMO_DIRECTORY)+proj.getProperty(Project.PLOT_DIRECTORY)+ext.rootOf(filename)+".scat", fingerprint, targets, markerData).serialize();
		demoSet = new MarkerSet(targets, demoChrs, demoPositions);
		demoSet.serialize(proj.getDir(Project.DEMO_DIRECTORY)+proj.getProperty(Project.MARKERSET_FILENAME));
		for (int i = 0; i<demoCents.length; i++) {
			new Centroids(demoCents[i], demoSet.getFingerprint()).serialize(proj.getDir(Project.DEMO_DIRECTORY)+files[i]);
        }
		sampleList.serialize(proj.getDir(Project.DEMO_DIRECTORY)+proj.getProperty(Project.SAMPLELIST_FILENAME));
		new MarkerLookup(hash).serialize(proj.getDir(Project.DEMO_DIRECTORY)+proj.getProperty(Project.MARKERLOOKUP_FILENAME));
		Files.copyFile(filename, proj.getDir(Project.DEMO_DIRECTORY)+proj.getProperty(Project.DISPLAY_MARKERS_FILENAME));
		Files.copyFile(proj.getFilename(Project.SAMPLE_DATA_FILENAME), proj.getDir(Project.DEMO_DIRECTORY)+proj.getProperty(Project.SAMPLE_DATA_FILENAME));
		System.out.println("Created demo files in "+ext.getTimeElapsed(time));
	}

	/*
	public static void extractAll(Project proj, int numCycles, boolean compact_XYonly) {
		BufferedReader reader;
		PrintWriter writer;
		SetOfSets[] setsOfSets;
		String[] markerNames;
		SampleList sampleList;
		MarkerSet markerSet;
		// Thread[] threads;
		long fingerprint;
		int[][] sets;
		int count;
		// boolean done;
		long heapMaxSize;
		int bytesPerMarker, numMarkersPerFile;
		boolean done;

		done = false;
		while (!done) {
			try {
				try {
					writer = new PrintWriter(new FileWriter(".current"));
					writer.println(proj.getProperty(Project.PROJECT_PROPERTIES_FILENAME));
					writer.println(numCycles);
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing backup file");
				}

				sampleList = proj.getSampleList();
				if (sampleList==null) {
					System.out.println("No SampleList file found; generating one using all the samples in: "+proj.getDir(Project.SAMPLE_DIRECTORY));
					sampleList = SampleList.generateSampleList(proj);
				}
				fingerprint = sampleList.getFingerprint();

				markerSet = proj.getMarkerSet();
				markerNames = markerSet.getMarkerNames();

				if (numCycles==0) {
					bytesPerMarker = 37*sampleList.getSamples().length+161;
					heapMaxSize = Runtime.getRuntime().maxMemory();
					numMarkersPerFile = (int)(((double)heapMaxSize*0.75)/bytesPerMarker);
					System.out.println("Heap size (memory) was set to "+ext.formDeci((double)heapMaxSize/1024/1024/1024, 2)+" Gb; which allows us to do "+numMarkersPerFile+" markers at a time comfortably");
					System.out.println("Suggestion: Use 90% of available memory (i.e. -Xmx15g if you have 16GB of memory available) and use a 64-bit version of java (and the -d64 option) if you plan to use more than 2GB of memory");
					numCycles = (int)Math.ceil((double)markerNames.length/(double)numMarkersPerFile);
				}
				numMarkersPerFile = (int)Math.ceil((double)markerNames.length/(double)numCycles);
				System.out.println("Max memory: "+Runtime.getRuntime().maxMemory());
				Runtime.getRuntime().gc();
				Runtime.getRuntime().runFinalization();
				Runtime.getRuntime().gc();
				System.out.println("Running using "+numCycles+" cycles ("+numMarkersPerFile+" markers per run)");

				count = 0;
				sets = new int[numCycles][];
				for (int i = 0; i<sets.length; i++) {
					sets[i] = new int[Math.min(numMarkersPerFile, markerNames.length-count)];
					for (int j = 0; j<sets[i].length; j++) {
						sets[i][j] = count;
						count++;
					}
				}

				int numBatches = 1;
				setsOfSets = new SetOfSets[numBatches];
				for (int i = 0; i<numBatches; i++) {
					setsOfSets[i] = new SetOfSets();
				}
				for (int i = 0; i<sets.length; i++) {
					setsOfSets[i%numBatches].add(sets[i], proj.getDir(Project.PLOT_DIRECTORY)+"Markers."+ext.formNum(i+1, 4)+".scat");
				}
				new File(proj.getDir(Project.PLOT_DIRECTORY)).mkdirs();

				// threads = new Thread[numBatches];
				// for (int i = 0; i < numBatches; i++) {
				// threads[i] = new Thread(new ExtractPlots(proj,
				// setsOfSets[i].getSets(), fingerprint,
				// setsOfSets[i].getFilenames()));
				// threads[i].start();
				// try {
				// Thread.sleep(100L);
				// } catch(InterruptedException ex) {}
				// }
				// do {
				// done = true;
				// for (int i = 0; i < numBatches; i++) {
				// if (threads[i].isAlive()) {
				// done = false;
				// }
				// }
				// try {
				// Thread.sleep(10000L);
				// } catch(InterruptedException ex) {}
				// } while (!done);

				for (int i = 0; i<numBatches; i++) {
					extractPlots(proj, setsOfSets[i].getSets(), fingerprint, compact_XYonly, setsOfSets[i].getFilenames());
				}

				createMarkerLookup(proj);
				done = true;
			} catch (OutOfMemoryError oome) {

				System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
				Runtime.getRuntime().gc();
				Runtime.getRuntime().runFinalization();
				Runtime.getRuntime().gc();
				System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());

				try {
					reader = new BufferedReader(new FileReader(".current"));
					while (reader.ready()) {
						proj = new Project(reader.readLine(), false);
						numCycles = Integer.parseInt(reader.readLine());
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+".current"+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+".current"+"\"");
					System.exit(2);
				}
				System.err.println(ext.getTime()+"\tAlgorithm failed to run using "+numCycles+" cycles; will try using "+(numCycles+1));
				numCycles++;
			}
		}
	}
	*/

	/*
		public static void extractAll(Project proj, int numCycles, boolean compact_XYonly) {
			BufferedReader reader;
			PrintWriter writer;
			SetOfSets[] setsOfSets;
			String[] markerNames;
			SampleList sampleList;
			MarkerSet markerSet;
			// Thread[] threads;
			long fingerprint;
			int[][] sets;
			int count;
			// boolean done;
			long heapMaxSize;
			int bytesPerMarker, numMarkersPerFile;
			boolean done;
	
			done = false;
			while (!done) {
				try {
					try {
						writer = new PrintWriter(new FileWriter(".current"));
						writer.println(proj.getProperty(Project.PROJECT_PROPERTIES_FILENAME));
						writer.println(numCycles);
						writer.close();
					} catch (Exception e) {
						System.err.println("Error writing backup file");
					}
	
					sampleList = proj.getSampleList();
					if (sampleList==null) {
						System.out.println("No SampleList file found; generating one using all the samples in: "+proj.getDir(Project.SAMPLE_DIRECTORY));
						sampleList = SampleList.generateSampleList(proj);
					}
					fingerprint = sampleList.getFingerprint();
	
					markerSet = proj.getMarkerSet();
					markerNames = markerSet.getMarkerNames();
	
					if (numCycles==0) {
						bytesPerMarker = 37*sampleList.getSamples().length+161;
						heapMaxSize = Runtime.getRuntime().maxMemory();
						numMarkersPerFile = (int)(((double)heapMaxSize*0.75)/bytesPerMarker);
						System.out.println("Heap size (memory) was set to "+ext.formDeci((double)heapMaxSize/1024/1024/1024, 2)+" Gb; which allows us to do "+numMarkersPerFile+" markers at a time comfortably");
						System.out.println("Suggestion: Use 90% of available memory (i.e. -Xmx15g if you have 16GB of memory available) and use a 64-bit version of java (and the -d64 option) if you plan to use more than 2GB of memory");
						numCycles = (int)Math.ceil((double)markerNames.length/(double)numMarkersPerFile);
					}
					numMarkersPerFile = (int)Math.ceil((double)markerNames.length/(double)numCycles);
					System.out.println("Max memory: "+Runtime.getRuntime().maxMemory());
					Runtime.getRuntime().gc();
					Runtime.getRuntime().runFinalization();
					Runtime.getRuntime().gc();
					System.out.println("Running using "+numCycles+" cycles ("+numMarkersPerFile+" markers per run)");
	
					count = 0;
					sets = new int[numCycles][];
					for (int i = 0; i<sets.length; i++) {
						sets[i] = new int[Math.min(numMarkersPerFile, markerNames.length-count)];
						for (int j = 0; j<sets[i].length; j++) {
							sets[i][j] = count;
							count++;
						}
					}
	
					int numBatches = 1;
					setsOfSets = new SetOfSets[numBatches];
					for (int i = 0; i<numBatches; i++) {
						setsOfSets[i] = new SetOfSets();
					}
					for (int i = 0; i<sets.length; i++) {
						setsOfSets[i%numBatches].add(sets[i], proj.getDir(Project.PLOT_DIRECTORY)+"Markers."+ext.formNum(i+1, 4)+".scat");
					}
					new File(proj.getDir(Project.PLOT_DIRECTORY)).mkdirs();
	
					// threads = new Thread[numBatches];
					// for (int i = 0; i < numBatches; i++) {
					// threads[i] = new Thread(new ExtractPlots(proj,
					// setsOfSets[i].getSets(), fingerprint,
					// setsOfSets[i].getFilenames()));
					// threads[i].start();
					// try {
					// Thread.sleep(100L);
					// } catch(InterruptedException ex) {}
					// }
					// do {
					// done = true;
					// for (int i = 0; i < numBatches; i++) {
					// if (threads[i].isAlive()) {
					// done = false;
					// }
					// }
					// try {
					// Thread.sleep(10000L);
					// } catch(InterruptedException ex) {}
					// } while (!done);
	
					for (int i = 0; i<numBatches; i++) {
						extractPlots(proj, setsOfSets[i].getSets(), fingerprint, compact_XYonly, setsOfSets[i].getFilenames());
					}
	
					createMarkerLookup(proj);
					done = true;
				} catch (OutOfMemoryError oome) {
	
					System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
					Runtime.getRuntime().gc();
					Runtime.getRuntime().runFinalization();
					Runtime.getRuntime().gc();
					System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
	
					try {
						reader = new BufferedReader(new FileReader(".current"));
						while (reader.ready()) {
							proj = new Project(reader.readLine(), false);
							numCycles = Integer.parseInt(reader.readLine());
						}
						reader.close();
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \""+".current"+"\" not found in current directory");
						System.exit(1);
					} catch (IOException ioe) {
						System.err.println("Error reading file \""+".current"+"\"");
						System.exit(2);
					}
					System.err.println(ext.getTime()+"\tAlgorithm failed to run using "+numCycles+" cycles; will try using "+(numCycles+1));
					numCycles++;
				}
			}
		}
		*/
	
		/*
		public static void extractAll(Project proj, int numCycles, boolean compact_XYonly) {
			BufferedReader reader;
			PrintWriter writer;
			SetOfSets[] setsOfSets;
			String[] markerNames;
			SampleList sampleList;
			MarkerSet markerSet;
			// Thread[] threads;
			long fingerprint;
			int[][] sets;
			int count;
			// boolean done;
			long heapMaxSize;
			int bytesPerMarker, numMarkersPerFile;
			boolean done;
	
			done = false;
			while (!done) {
				try {
					try {
						writer = new PrintWriter(new FileWriter(".current"));
						writer.println(proj.getProperty(Project.PROJECT_PROPERTIES_FILENAME));
						writer.println(numCycles);
						writer.close();
					} catch (Exception e) {
						System.err.println("Error writing backup file");
					}
	
					sampleList = proj.getSampleList();
					if (sampleList==null) {
						System.out.println("No SampleList file found; generating one using all the samples in: "+proj.getDir(Project.SAMPLE_DIRECTORY));
						sampleList = SampleList.generateSampleList(proj);
					}
					fingerprint = sampleList.getFingerprint();
	
					markerSet = proj.getMarkerSet();
					markerNames = markerSet.getMarkerNames();
	
					if (numCycles==0) {
						bytesPerMarker = 37*sampleList.getSamples().length+161;
						heapMaxSize = Runtime.getRuntime().maxMemory();
						numMarkersPerFile = (int)(((double)heapMaxSize*0.75)/bytesPerMarker);
						System.out.println("Heap size (memory) was set to "+ext.formDeci((double)heapMaxSize/1024/1024/1024, 2)+" Gb; which allows us to do "+numMarkersPerFile+" markers at a time comfortably");
						System.out.println("Suggestion: Use 90% of available memory (i.e. -Xmx15g if you have 16GB of memory available) and use a 64-bit version of java (and the -d64 option) if you plan to use more than 2GB of memory");
						numCycles = (int)Math.ceil((double)markerNames.length/(double)numMarkersPerFile);
					}
					numMarkersPerFile = (int)Math.ceil((double)markerNames.length/(double)numCycles);
					System.out.println("Max memory: "+Runtime.getRuntime().maxMemory());
					Runtime.getRuntime().gc();
					Runtime.getRuntime().runFinalization();
					Runtime.getRuntime().gc();
					System.out.println("Running using "+numCycles+" cycles ("+numMarkersPerFile+" markers per run)");
	
					count = 0;
					sets = new int[numCycles][];
					for (int i = 0; i<sets.length; i++) {
						sets[i] = new int[Math.min(numMarkersPerFile, markerNames.length-count)];
						for (int j = 0; j<sets[i].length; j++) {
							sets[i][j] = count;
							count++;
						}
					}
	
					int numBatches = 1;
					setsOfSets = new SetOfSets[numBatches];
					for (int i = 0; i<numBatches; i++) {
						setsOfSets[i] = new SetOfSets();
					}
					for (int i = 0; i<sets.length; i++) {
						setsOfSets[i%numBatches].add(sets[i], proj.getDir(Project.PLOT_DIRECTORY)+"Markers."+ext.formNum(i+1, 4)+".scat");
					}
					new File(proj.getDir(Project.PLOT_DIRECTORY)).mkdirs();
	
					// threads = new Thread[numBatches];
					// for (int i = 0; i < numBatches; i++) {
					// threads[i] = new Thread(new ExtractPlots(proj,
					// setsOfSets[i].getSets(), fingerprint,
					// setsOfSets[i].getFilenames()));
					// threads[i].start();
					// try {
					// Thread.sleep(100L);
					// } catch(InterruptedException ex) {}
					// }
					// do {
					// done = true;
					// for (int i = 0; i < numBatches; i++) {
					// if (threads[i].isAlive()) {
					// done = false;
					// }
					// }
					// try {
					// Thread.sleep(10000L);
					// } catch(InterruptedException ex) {}
					// } while (!done);
	
					for (int i = 0; i<numBatches; i++) {
						extractPlots(proj, setsOfSets[i].getSets(), fingerprint, compact_XYonly, setsOfSets[i].getFilenames());
					}
	
					createMarkerLookup(proj);
					done = true;
				} catch (OutOfMemoryError oome) {
	
					System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
					Runtime.getRuntime().gc();
					Runtime.getRuntime().runFinalization();
					Runtime.getRuntime().gc();
					System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
	
					try {
						reader = new BufferedReader(new FileReader(".current"));
						while (reader.ready()) {
							proj = new Project(reader.readLine(), false);
							numCycles = Integer.parseInt(reader.readLine());
						}
						reader.close();
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \""+".current"+"\" not found in current directory");
						System.exit(1);
					} catch (IOException ioe) {
						System.err.println("Error reading file \""+".current"+"\"");
						System.exit(2);
					}
					System.err.println(ext.getTime()+"\tAlgorithm failed to run using "+numCycles+" cycles; will try using "+(numCycles+1));
					numCycles++;
				}
			}
		}
		*/
	
		/*
			public static void extractAll(Project proj, int numCycles, boolean compact_XYonly) {
				BufferedReader reader;
				PrintWriter writer;
				SetOfSets[] setsOfSets;
				String[] markerNames;
				SampleList sampleList;
				MarkerSet markerSet;
				// Thread[] threads;
				long fingerprint;
				int[][] sets;
				int count;
				// boolean done;
				long heapMaxSize;
				int bytesPerMarker, numMarkersPerFile;
				boolean done;
		
				done = false;
				while (!done) {
					try {
						try {
							writer = new PrintWriter(new FileWriter(".current"));
							writer.println(proj.getProperty(Project.PROJECT_PROPERTIES_FILENAME));
							writer.println(numCycles);
							writer.close();
						} catch (Exception e) {
							System.err.println("Error writing backup file");
						}
		
						sampleList = proj.getSampleList();
						if (sampleList==null) {
							System.out.println("No SampleList file found; generating one using all the samples in: "+proj.getDir(Project.SAMPLE_DIRECTORY));
							sampleList = SampleList.generateSampleList(proj);
						}
						fingerprint = sampleList.getFingerprint();
		
						markerSet = proj.getMarkerSet();
						markerNames = markerSet.getMarkerNames();
		
						if (numCycles==0) {
							bytesPerMarker = 37*sampleList.getSamples().length+161;
							heapMaxSize = Runtime.getRuntime().maxMemory();
							numMarkersPerFile = (int)(((double)heapMaxSize*0.75)/bytesPerMarker);
							System.out.println("Heap size (memory) was set to "+ext.formDeci((double)heapMaxSize/1024/1024/1024, 2)+" Gb; which allows us to do "+numMarkersPerFile+" markers at a time comfortably");
							System.out.println("Suggestion: Use 90% of available memory (i.e. -Xmx15g if you have 16GB of memory available) and use a 64-bit version of java (and the -d64 option) if you plan to use more than 2GB of memory");
							numCycles = (int)Math.ceil((double)markerNames.length/(double)numMarkersPerFile);
						}
						numMarkersPerFile = (int)Math.ceil((double)markerNames.length/(double)numCycles);
						System.out.println("Max memory: "+Runtime.getRuntime().maxMemory());
						Runtime.getRuntime().gc();
						Runtime.getRuntime().runFinalization();
						Runtime.getRuntime().gc();
						System.out.println("Running using "+numCycles+" cycles ("+numMarkersPerFile+" markers per run)");
		
						count = 0;
						sets = new int[numCycles][];
						for (int i = 0; i<sets.length; i++) {
							sets[i] = new int[Math.min(numMarkersPerFile, markerNames.length-count)];
							for (int j = 0; j<sets[i].length; j++) {
								sets[i][j] = count;
								count++;
							}
						}
		
						int numBatches = 1;
						setsOfSets = new SetOfSets[numBatches];
						for (int i = 0; i<numBatches; i++) {
							setsOfSets[i] = new SetOfSets();
						}
						for (int i = 0; i<sets.length; i++) {
							setsOfSets[i%numBatches].add(sets[i], proj.getDir(Project.PLOT_DIRECTORY)+"Markers."+ext.formNum(i+1, 4)+".scat");
						}
						new File(proj.getDir(Project.PLOT_DIRECTORY)).mkdirs();
		
						// threads = new Thread[numBatches];
						// for (int i = 0; i < numBatches; i++) {
						// threads[i] = new Thread(new ExtractPlots(proj,
						// setsOfSets[i].getSets(), fingerprint,
						// setsOfSets[i].getFilenames()));
						// threads[i].start();
						// try {
						// Thread.sleep(100L);
						// } catch(InterruptedException ex) {}
						// }
						// do {
						// done = true;
						// for (int i = 0; i < numBatches; i++) {
						// if (threads[i].isAlive()) {
						// done = false;
						// }
						// }
						// try {
						// Thread.sleep(10000L);
						// } catch(InterruptedException ex) {}
						// } while (!done);
		
						for (int i = 0; i<numBatches; i++) {
							extractPlots(proj, setsOfSets[i].getSets(), fingerprint, compact_XYonly, setsOfSets[i].getFilenames());
						}
		
						createMarkerLookup(proj);
						done = true;
					} catch (OutOfMemoryError oome) {
		
						System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
						Runtime.getRuntime().gc();
						Runtime.getRuntime().runFinalization();
						Runtime.getRuntime().gc();
						System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
		
						try {
							reader = new BufferedReader(new FileReader(".current"));
							while (reader.ready()) {
								proj = new Project(reader.readLine(), false);
								numCycles = Integer.parseInt(reader.readLine());
							}
							reader.close();
						} catch (FileNotFoundException fnfe) {
							System.err.println("Error: file \""+".current"+"\" not found in current directory");
							System.exit(1);
						} catch (IOException ioe) {
							System.err.println("Error reading file \""+".current"+"\"");
							System.exit(2);
						}
						System.err.println(ext.getTime()+"\tAlgorithm failed to run using "+numCycles+" cycles; will try using "+(numCycles+1));
						numCycles++;
					}
				}
			}
			*/
		
	public static void transposeData(Project proj, long maxMarkerFileSize, boolean keepAllSampleFilesOpen) {
		String[] allSamplesInProj;
		String[] allMarkersInProj;
		int numBytesPerMarker;
		int numMarkersInWriterBuffer;
		int numMarkersLastReadRoundLastBuffer = 0;
        int numMarkersPerBufferChunk;
		int maxNumMarkersPerFile;
//		int[] numMarkersInEachFile;
		int numChunksInWriteBuffer;
		int numBufferChunksPerFile;
		int numMarkerFiles;
		int numMarkersLastReadRound;
		int timesEachSampleFileBeRead;
		String[] markerFilenames;
		int markerFileIndex;
		int markerFileBufferChunkIndex;
		int firstMarkerInCurrentFile;
		int markerIndex;
//		int writeBufferSeekStep;
		RandomAccessFile[] sampleFiles = null;
		RandomAccessFile sampleFile;
        RandomAccessFile markerFile;
        byte[][] markerDataWriteBuffer = null;
        byte[] markerDataWriteBufferParameter;
        byte[] sampleDataReadBuffer;
        int indexInSampleDataReadBuffer;
        int locationInWriteBufferChunk;
		Hashtable<String,String> markerLookupHash = new Hashtable<String,String>();
		byte[][] markersInEachFile;
		String[] markersInEachFile1;
		byte[] markerNamesByteStream;
//		byte[] sampFileMarkerNameReadBuffer;
		int numSampFileOutliers;
		byte[] sampFileOutliersReadBuffer = null;
//		float[] outOfRangeValuesInSampleFile;
		Hashtable<String, Float> outOfRangeValuesInSampleFile = null;
//		int indexOutOfRangeValuesInSampleFile;
		Hashtable<String, Float>[] markFileWriteBufferOutliers;
		Hashtable<String, Float> markFileWriteBufferOutliersAdj;
		byte[] markFileWriteBufferOutliersBytes;
//		int indexMarkFileOutliers;
		Hashtable<String, Float> allOutliers;
		byte nullStatus = 0;
        byte bytesPerSampMark;
        int numMarkersCurrentLoop;
		long timer1, timer2, timer3, timer4, timer5, timer6;

		allSamplesInProj = proj.getSamples();
//		String[] samplesTemp = proj.getSamples();
//		allSamplesInProj = new String[10];
//		for (int i=0; i<10; i++) {
//			allSamplesInProj[i] = samplesTemp[i];
//		}
		allMarkersInProj = proj.getMarkerNames();
		try {
			sampleFile = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSamplesInProj[0] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
			sampleFile.readInt();
			nullStatus = sampleFile.readByte();
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			return;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
		numBytesPerMarker = allSamplesInProj.length * bytesPerSampMark;
//		writeBufferSeekStep = numBytesPerMarker - BYTES_PER_SAMPLE_MARKER;
		if (new File(proj.getProjectDir()).getFreeSpace() <= (allSamplesInProj.length * (long)allMarkersInProj.length * bytesPerSampMark)) {
			System.err.println("Not enough space (available: "+ext.prettyUpSize(new File(proj.getProjectDir()).getFreeSpace(), 1)+") for all the new data to be created (required: "+ext.prettyUpSize(new File(proj.getProjectDir()).getFreeSpace(), 1)+").");
			return;
		}


		if (maxMarkerFileSize == 0) {
			maxMarkerFileSize = Integer.MAX_VALUE;
		}
		numMarkersInWriterBuffer = Math.min((int)((0.65*(double)Runtime.getRuntime().maxMemory())/numBytesPerMarker), allMarkersInProj.length);
//		numMarkersPerBufferChunk = Math.min((int) Integer.MAX_VALUE/numBytesPerMarker, numMarkersInWriterBuffer);
		if (maxMarkerFileSize > Integer.MAX_VALUE) {
			numMarkersPerBufferChunk = (int) Math.min((double) Integer.MAX_VALUE/numBytesPerMarker, (double) numMarkersInWriterBuffer);
		} else {
			numMarkersPerBufferChunk = (int) Math.min((double) maxMarkerFileSize/numBytesPerMarker, (double) numMarkersInWriterBuffer);
		}
		numChunksInWriteBuffer = numMarkersInWriterBuffer/numMarkersPerBufferChunk;
		numMarkersInWriterBuffer = numChunksInWriteBuffer * numMarkersPerBufferChunk;
		timesEachSampleFileBeRead = (int) Math.ceil((double)allMarkersInProj.length / (double)numMarkersInWriterBuffer);
		numMarkersLastReadRound = allMarkersInProj.length % numMarkersInWriterBuffer;

		maxNumMarkersPerFile = (int) Math.min((double)maxMarkerFileSize/(double)numBytesPerMarker, allMarkersInProj.length);
		//TODO remove here
//		numBufferChunksPerFile = maxNumMarkersPerFile / numMarkersPerBufferChunk;
		numBufferChunksPerFile = (int) Math.ceil((double)maxNumMarkersPerFile / (double)numMarkersPerBufferChunk);
		maxNumMarkersPerFile = numMarkersPerBufferChunk * numBufferChunksPerFile;
		numMarkerFiles = (int) Math.ceil((double)allMarkersInProj.length / (double)maxNumMarkersPerFile);
		markerFilenames = new String[numMarkerFiles];
		markersInEachFile = new byte[numMarkerFiles][];
//		numMarkersInEachFile = new int[numMarkerFiles];
		backupOlderRafs(proj.getDir(Project.PLOT_DIRECTORY), MARKDATA_FILE_EXTENSION);
		markerIndex = 0;
		for (int i=0; i<numMarkerFiles; i++) {
			markerFilenames[i]=proj.getDir(Project.PLOT_DIRECTORY)+"markers." + i + MARKDATA_FILE_EXTENSION;
			numMarkersCurrentLoop = Math.min(maxNumMarkersPerFile, allMarkersInProj.length - markerIndex);
			markersInEachFile1 = new String[numMarkersCurrentLoop];
//			numMarkersInEachFile[i] =; 
			for (int j=0; j<numMarkersCurrentLoop; j++) {
				markerLookupHash.put(allMarkersInProj[markerIndex], "markers." + i + MARKDATA_FILE_EXTENSION + "\t" + j);
				markersInEachFile1[j] = allMarkersInProj[markerIndex];
				markerIndex ++;
			}
			try {
				markersInEachFile[i] = Compression.objToBytes(markersInEachFile1);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		new MarkerLookup(markerLookupHash).serialize(proj.getFilename(Project.MARKERLOOKUP_FILENAME));
		System.out.println("Total markers in the project:\t"+allMarkersInProj.length
						  +"\nTotal samples in the project:\t"+allSamplesInProj.length
						  +"\nMaximum memory size available:\t"+Runtime.getRuntime().maxMemory()/1024/1024/1024+"."+(Runtime.getRuntime().maxMemory()/1024/1024/10 - Runtime.getRuntime().maxMemory()/1024/1024/1024*102)+" gb"
						  +"\nNumber of markers in write buffer at a time:\t"+numMarkersInWriterBuffer
						  +"\nNumber of markers per write buffer Chunk:\t"+numMarkersPerBufferChunk
						  +"\nMaximum file size of marker files:\t"+maxMarkerFileSize/1024/1024/1024+"."+((int)(maxMarkerFileSize/1024/1024/10.24) - (int)(maxMarkerFileSize/1024/1024/1024*102.4))+" gb"
						  +"\nNumber of markers per marker file:\t"+maxNumMarkersPerFile
						  +"\nNumber of write buffer chunks per marker file:\t"+numBufferChunksPerFile
						  +"\nNumber of marker files:\t"+numMarkerFiles);

		timer1 = new Date().getTime();
		timer3 = 0;
		firstMarkerInCurrentFile = 0;
		sampleDataReadBuffer = new byte[numMarkersPerBufferChunk * bytesPerSampMark];
		markerDataWriteBuffer = new byte[numChunksInWriteBuffer][numMarkersPerBufferChunk * numBytesPerMarker];
//		markFileOutliersWriteBufferVec = new Vector[numChunksInWriteBuffer];
		markFileWriteBufferOutliers = new Hashtable[numChunksInWriteBuffer];
		allOutliers = new Hashtable<String, Float>();
		try {
			if (keepAllSampleFilesOpen) {
				sampleFiles = new RandomAccessFile[allSamplesInProj.length];
				for (int i=0; i<allSamplesInProj.length; i++) {
					sampleFiles[i] = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSamplesInProj[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
				}
			}

			markerFileIndex = 0;

			markerNamesByteStream = Compression.objToBytes(markersInEachFile[markerFileIndex]);
			markerDataWriteBufferParameter = new byte[MARKDATA_PARAMETER_TOTAL_LEN + markerNamesByteStream.length];
			System.arraycopy(Compression.intToBytes(allSamplesInProj.length), 0, markerDataWriteBufferParameter, MARKDATA_NUMSAMPS_START, MARKDATA_NUMSAMPS_LEN);
			System.arraycopy(Compression.intToBytes(maxNumMarkersPerFile), 0, markerDataWriteBufferParameter, MARKDATA_NUMMARKS_START, MARKDATA_NUMMARKS_LEN);
			markerDataWriteBufferParameter[MARKDATA_NULLSTATUS_START] = nullStatus;
			System.arraycopy(Compression.longToBytes(MarkerSet.fingerprint(allSamplesInProj)), 0, markerDataWriteBufferParameter, MARKDATA_FINGERPRINT_START, MARKDATA_FINGERPRINT_LEN);
			System.arraycopy(Compression.intToBytes(markerNamesByteStream.length), 0, markerDataWriteBufferParameter, MARKDATA_MARKNAMELEN_START, MARKDATA_MARKNAMELEN_LEN);
			System.arraycopy(markerNamesByteStream, 0, markerDataWriteBufferParameter, MARKDATA_MARKNAME_START, markerNamesByteStream.length);
			markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
			markerFile.write(markerDataWriteBufferParameter);

//			markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
////			sampFileMarkerNameReadBuffer = markersInEachFile[markerFileIndex].getBytes();
//			markerFile.writeInt(allSamplesInProj.length);
//			markerFile.writeInt(maxNumMarkersPerFile);
//			markerFile.writeByte(nullStatus);
//			markerFile.writeLong(MarkerSet.fingerprint(allSamplesInProj));
////			markerFile.writeInt(sampFileMarkerNameReadBuffer.length);
////			markerFile.write(sampFileMarkerNameReadBuffer);
////			markerFile.writeInt(markersInEachFile[markerFileIndex].length);
////			markerFile.write(markersInEachFile[markerFileIndex]);
//			markerFile.writeInt(markerNamesByteStream.length);
//			markerFile.write(markerNamesByteStream);

//			markFileOutliersWriteBufferVecAdj = new Vector<Byte>();
			markFileWriteBufferOutliersAdj = new Hashtable<String, Float>();
			markerFileBufferChunkIndex = 0;
			for(int i=0; i<timesEachSampleFileBeRead; i++) {
				timer2 = new Date().getTime();
				timer5 = 0;
				timer6 = 0;
				if ((i+1)==timesEachSampleFileBeRead && numMarkersLastReadRound!=0) {
					numChunksInWriteBuffer = (int) Math.ceil((double)numMarkersLastReadRound / (double)numMarkersPerBufferChunk);
//					markFileOutliersWriteBufferVec = new Vector[numChunksInWriteBuffer];
					markFileWriteBufferOutliers = new Hashtable[numChunksInWriteBuffer];
					if ( numChunksInWriteBuffer>1 ) {
						markerDataWriteBuffer = null;
						markerDataWriteBuffer = new byte[numChunksInWriteBuffer][numMarkersPerBufferChunk * numBytesPerMarker];
						numMarkersLastReadRoundLastBuffer = numMarkersLastReadRound % numMarkersPerBufferChunk;
						if ( numMarkersLastReadRoundLastBuffer!=0 ) {
							markerDataWriteBuffer[numChunksInWriteBuffer-1] = new byte[numMarkersLastReadRoundLastBuffer * numBytesPerMarker];
						}
					} else {
						markerDataWriteBuffer = null;
						markerDataWriteBuffer = new byte[numChunksInWriteBuffer][numMarkersLastReadRound * numBytesPerMarker];
					}
				}
//				for (int j=0; j<markFileOutliersWriteBufferVec.length; j++) {
//					markFileOutliersWriteBufferVec[j] = new Vector<Byte>();
//				}
				for (int j=0; j<markFileWriteBufferOutliers.length; j++) {
					markFileWriteBufferOutliers[j] = new Hashtable<String, Float>();
				}
				for (int j=0; j<allSamplesInProj.length; j++) {
					if (!keepAllSampleFilesOpen) {
						sampleFile = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSamplesInProj[j] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r"); //TODO
					} else {
						sampleFile = sampleFiles[j];
					}
					markerIndex = firstMarkerInCurrentFile;
				    for (int k=0; k < numChunksInWriteBuffer; k++) {
				    	timer4 = new Date().getTime();
//				    	sampleFiles[j].read(sampleDataReadBuffer);
				    	sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + markerIndex * bytesPerSampMark);
				    	sampleFile.read(sampleDataReadBuffer);
				    	// Read in the over range value array
				    	sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + allMarkersInProj.length * bytesPerSampMark);
				    	numSampFileOutliers = sampleFile.readInt();
				    	if (numSampFileOutliers > 0) {
				    		sampFileOutliersReadBuffer = new byte[numSampFileOutliers];
//				    		outOfRangeValuesInSampleFile = new float[numOutOfRangeValuesInSampleFile];
				    		sampleFile.read(sampFileOutliersReadBuffer);
//				    		outOfRangeValuesInSampleFile = new Hashtable<String, Float>();
				    		outOfRangeValuesInSampleFile = (Hashtable<String, Float>) Compression.bytesToObj(sampFileOutliersReadBuffer);
				    	}
				    	indexInSampleDataReadBuffer = 0;
				    	timer5 += (new Date().getTime() - timer4);
				    	timer4 = new Date().getTime();
						locationInWriteBufferChunk = j * bytesPerSampMark;
				    	for (int l=0; l<numMarkersPerBufferChunk && markerIndex<allMarkersInProj.length; l++) {
				    		for (int m=0; m < bytesPerSampMark; m++) {
			    				markerDataWriteBuffer[k][locationInWriteBufferChunk + m] = sampleDataReadBuffer[indexInSampleDataReadBuffer];
				    			if (m==3 && (sampleDataReadBuffer[indexInSampleDataReadBuffer-1] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[0]
				    					     && sampleDataReadBuffer[indexInSampleDataReadBuffer] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[1])) {
				    				allOutliers.put(markerIndex + "\t" + j +"\tx", outOfRangeValuesInSampleFile.get(markerIndex + "\tx"));
				    				markFileWriteBufferOutliers[k].put(markerIndex + "\t" + j + "\tx", outOfRangeValuesInSampleFile.get(markerIndex + "\tx"));
				    			} else if (m==5 && (sampleDataReadBuffer[indexInSampleDataReadBuffer-1] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[0]
			    					                && sampleDataReadBuffer[indexInSampleDataReadBuffer] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[1])) {
				    				allOutliers.put(markerIndex + "\t" + j + "\ty", outOfRangeValuesInSampleFile.get(markerIndex + "\ty"));
				    				markFileWriteBufferOutliers[k].put(markerIndex + "\t" + j + "\ty", outOfRangeValuesInSampleFile.get(markerIndex + "\ty"));
				    			} else if (m==10 && (sampleDataReadBuffer[indexInSampleDataReadBuffer-2] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[0]
				    								 && sampleDataReadBuffer[indexInSampleDataReadBuffer-1] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[1]
				    								 && sampleDataReadBuffer[indexInSampleDataReadBuffer] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[2])) {
				    				allOutliers.put(markerIndex + "\t" + j + "\tlrr", outOfRangeValuesInSampleFile.get(markerIndex + "\tlrr"));
				    				markFileWriteBufferOutliers[k].put(markerIndex + "\t" + j + "\tlrr", outOfRangeValuesInSampleFile.get(markerIndex + "\tlrr"));
				    			}
//				    			if ((m==2 && (0x80 & sampleDataReadBuffer[indexInSampleDataReadBuffer])>0) || (m==4 && (0x80 & sampleDataReadBuffer[indexInSampleDataReadBuffer])>0)) {
//				    				allOutliers.put(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y"), outOfRangeValuesInSampleFile.get(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y")));
//				    				markFileWriteBufferOutliers[k].put(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y"), outOfRangeValuesInSampleFile.get(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y")));
//				    			}
			    				indexInSampleDataReadBuffer ++;
					    	}
				    		locationInWriteBufferChunk += numBytesPerMarker;
							markerIndex ++;
				    	}
				    	timer6 += (new Date().getTime() - timer4);
					}
				    if (!keepAllSampleFilesOpen) {
				    	sampleFile.close();
				    }
				}
				timer2 = new Date().getTime() - timer2;
				System.out.println("Load in sample files - round " + i + ": " + timer2/1000 + " secs = "+timer5/1000+" (file.read) + "+timer6/1000+" (array shuffle)");
				for (int j=0; j<markerDataWriteBuffer.length; j++) {
			    	timer4 = new Date().getTime();
					markerFile.write(markerDataWriteBuffer[j]);
					markFileWriteBufferOutliersAdj.putAll(markFileWriteBufferOutliers[j]);
			    	timer3 += (new Date().getTime() - timer4);
					markerFileBufferChunkIndex ++;
					if (markerFileBufferChunkIndex >= numBufferChunksPerFile) {
						if (markFileWriteBufferOutliersAdj == null || markFileWriteBufferOutliersAdj.size() == 0) {
							markerFile.writeInt(0);
						} else {
							markFileWriteBufferOutliersBytes = Compression.objToBytes(markFileWriteBufferOutliersAdj);
							markerFile.writeInt(markFileWriteBufferOutliersBytes.length);
							markerFile.write(markFileWriteBufferOutliersBytes);
						}
						markerFile.close();
						System.out.println("Write marker file " + markerFileIndex + ":\t"+ timer3/1000 + " sec."); // \toutlier #: " + (markFileOUtliersWriteBufferArray.length/4));
						markerFileBufferChunkIndex = 0;
						markerFileIndex ++;
						markFileWriteBufferOutliersAdj = new Hashtable<String, Float>();
						if ((i+1) != timesEachSampleFileBeRead || (j+1) != markerDataWriteBuffer.length) {
//							markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
////							sampFileMarkerNameReadBuffer = markersInEachFile[markerFileIndex].getBytes();
//							markerFile.writeInt(allSamplesInProj.length);
//							if ((markerFileIndex+1) == numMarkerFiles) {
//								markerFile.writeInt(allMarkersInProj.length % maxNumMarkersPerFile);
//							} else {
//								markerFile.writeInt(maxNumMarkersPerFile);
//							}
//							markerFile.writeByte(nullStatus);
//							markerFile.writeLong(MarkerSet.fingerprint(allSamplesInProj));
////							markerFile.writeInt(sampFileMarkerNameReadBuffer.length);
////							markerFile.write(sampFileMarkerNameReadBuffer);
//							markerFile.writeInt(markersInEachFile[markerFileIndex].length);
//							markerFile.write(markersInEachFile[markerFileIndex]);

							markerNamesByteStream = Compression.objToBytes(markersInEachFile[markerFileIndex]);
							markerDataWriteBufferParameter = new byte[MARKDATA_PARAMETER_TOTAL_LEN + markerNamesByteStream.length];
							System.arraycopy(Compression.intToBytes(allSamplesInProj.length), 0, markerDataWriteBufferParameter, MARKDATA_NUMSAMPS_START, MARKDATA_NUMSAMPS_LEN);
							System.arraycopy(Compression.intToBytes((markerFileIndex+1)==numMarkerFiles? (allMarkersInProj.length % maxNumMarkersPerFile) : maxNumMarkersPerFile), 0, markerDataWriteBufferParameter, MARKDATA_NUMMARKS_START, MARKDATA_NUMMARKS_LEN);
							markerDataWriteBufferParameter[MARKDATA_NULLSTATUS_START] = nullStatus;
							System.arraycopy(Compression.longToBytes(MarkerSet.fingerprint(allSamplesInProj)), 0, markerDataWriteBufferParameter, MARKDATA_FINGERPRINT_START, MARKDATA_FINGERPRINT_LEN);
							System.arraycopy(Compression.intToBytes(markerNamesByteStream.length), 0, markerDataWriteBufferParameter, MARKDATA_MARKNAMELEN_START, MARKDATA_MARKNAMELEN_LEN);
							System.arraycopy(markerNamesByteStream, 0, markerDataWriteBufferParameter, MARKDATA_MARKNAME_START, markerNamesByteStream.length);
							markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
							markerFile.write(markerDataWriteBufferParameter);

//							markFileOutliersWriteBufferVecAdj = new Vector<Byte>();
						}
						timer3 = 0;
					}
				}
				firstMarkerInCurrentFile += numMarkersInWriterBuffer;
			}
			if (markerFileBufferChunkIndex != 0) {
				if (markFileWriteBufferOutliersAdj == null || markFileWriteBufferOutliersAdj.size() == 0) {
					markerFile.writeInt(0);
				} else {
					markFileWriteBufferOutliersBytes = Compression.objToBytes(markFileWriteBufferOutliersAdj);
					markerFile.writeInt(markFileWriteBufferOutliersBytes.length);
					markerFile.write(markFileWriteBufferOutliersBytes);
				}
				markerFile.close();
			}
			if (allOutliers != null && allOutliers.size() != 0) {
				Files.writeSerial(allOutliers, proj.getDir(Project.PLOT_DIRECTORY) + "_outliers." + MARKDATA_FILE_EXTENSION);
			}
			System.out.println("Write marker file " + markerFileIndex + ":\t" + timer3/1000 + " sec."); // \toutlier #: " + (markFileOUtliersWriteBufferArray.length/4));
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		timer1 = (new Date().getTime() - timer1);
		System.out.println("******\nTotal Time: "+ (timer1/3600000) +" hrs "+ (timer1/60000 - 60*(timer1/3600000)) +" mins "+ (timer1/1000 - 60*(timer1/60000)) +" secs.");
	}

	/*
	 * Delete existing files of a specific file name extension.
	 * @param directory	the folder where the existing files are
	 * @param fileNameExtension	the file name extension, including the ".". For example, ".scat", ".xls", ".csv", and etc. 
	 */
	public static void deleteOlderRafs(String directory, final String fileNameExtesion) {
		File[] files;

		// List of files to be backed up.
		files = new File(directory).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(fileNameExtesion);
			}
		});

		// Move files to backup folder
		for (int i = 0; i<files.length; i++) {
			files[i].delete();
		}
		
		if (files.length >= 0) {
			System.out.println("Older version of the data in " + directory + " has been deleted from the hard drive.");
		}
	}

	/*
	 * Move existing files of a specific file name extension to a sub-folder named "Backup.?"
	 * @param directory	the folder where the existing files are
	 * @param fileNameExtension	the file name extension, including the ".". For example, ".scat", ".xls", ".csv", and etc. 
	 */
	public static void backupOlderRafs(String directory, final String fileNameExtesion) {
		File[] files;
		int count;

		// Check existing backup folders.
		count = 0;
		do {
			count++;
		} while (new File(directory+"Backup."+count).exists());

		// List of files to be backed up.
		files = new File(directory).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(fileNameExtesion);
			}
		});

		if (files.length>0) {
			// Create a new backup folder.
			new File(directory+"Backup."+count+"/").mkdirs();
			
			// Move files to backup folder
			for (int i = 0; i<files.length; i++) {
				files[i].renameTo(new File(directory + "Backup." + count + "/" + files[i].getName()));
			}

			System.out.println("Older version of the data in " + directory + " has been moved to " + directory + "Backup." + count + "/ To save disk space, please manually delete these files.");
		}
	}

	/*
	public static void breakUpMarkerCollections(Project proj, int numMarkersInSmallerFiles) {
		boolean jar = false;
		MarkerDataCollection collection;
		MarkerData[] newMarkerData, oldMarkerData = null;
		long time, fingerprint = 0;
		File[] files;
		int count, oldFileIndex, newFileIndex, oldMarkerIndex;
		String[] newMarkerNames, oldMarkerNames = null;

		time = new Date().getTime();
		count = 0;
		do {
			count++;
		} while (new File(proj.getDir(Project.PLOT_DIRECTORY)+"Original."+count).exists());

		files = new File(proj.getDir(Project.PLOT_DIRECTORY)).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".scat");
			}
		});
		new File(proj.getDir(Project.PLOT_DIRECTORY)+"Original."+count+"/").mkdirs();
		for (int i = 0; i<files.length; i++) {
			files[i].renameTo(new File(proj.getDir(Project.PLOT_DIRECTORY)+"Original."+count+"/"+files[i].getName()));
		}
		files = new File(proj.getDir(Project.PLOT_DIRECTORY)+"Original."+count+"/").listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".scat");
			}
		});
		oldFileIndex = newFileIndex = -1;
		oldMarkerIndex = 0;
		boolean done = false;
		try {
			while (!done) {
				newMarkerData = new MarkerData[oldFileIndex==files.length-1&&oldMarkerData.length-oldMarkerIndex<numMarkersInSmallerFiles?oldMarkerData.length-oldMarkerIndex:numMarkersInSmallerFiles];
				newMarkerNames = new String[newMarkerData.length];
				for (int i = 0; i<newMarkerData.length; i++) {
					if (oldFileIndex==-1||oldMarkerIndex==oldMarkerData.length) {
						oldFileIndex++;
						collection = MarkerDataCollection.load(files[oldFileIndex].getCanonicalPath(), jar);
						oldMarkerData = collection.getCollection();
						oldMarkerNames = collection.getMarkerNames();
						fingerprint = collection.getFingerprint();
						oldMarkerIndex = 0;
					}
					newMarkerData[i] = oldMarkerData[oldMarkerIndex];
					newMarkerNames[i] = oldMarkerNames[oldMarkerIndex];
					oldMarkerIndex++;
				}
				new MarkerDataCollection(proj.getDir(Project.PLOT_DIRECTORY)+"Markers."+ext.formNum(++newFileIndex+1, 4)+".scat", fingerprint, newMarkerNames, newMarkerData).serialize();
				done = oldFileIndex==files.length-1&&oldMarkerIndex==oldMarkerData.length;
			}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		createMarkerLookup(proj);
		System.out.println("Reformatted "+(oldFileIndex+1)+" files into "+(newFileIndex+1)+" in "+ext.getTimeElapsed(time));
	}
	*/

	public static void createMarkerLookup(Project proj) {
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] files, markerNames;
		MarkerDataCollection collection;
		long time;

		time = new Date().getTime();
		System.out.println("Creating MarkerLookup file");
		files = new File(proj.getDir(Project.PLOT_DIRECTORY)).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".scat");
			}
		});
		if (files==null) {
			System.err.println("Error - failed to create MarkerLookup -- directory does not exist: "+proj.getDir(Project.PLOT_DIRECTORY));
		} else if (files.length==0) {
			System.err.println("Error - failed to create MarkerLookup -- no .scat files available");
		} else {
			for (int i = 0; i<files.length; i++) {
				collection = MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+files[i], false);
				markerNames = collection.getMarkerNames();
				for (int j = 0; j<markerNames.length; j++) {
					hash.put(markerNames[j], files[i]+"\t"+j);
				}
			}

			new MarkerLookup(hash).serialize(proj.getFilename(Project.MARKERLOOKUP_FILENAME));

			System.out.println("Created MarkerLookup in "+ext.getTimeElapsed(time));
		}
	}

	/**
	 * Create MarkerLookup from MarkerData (RAF and half-precision format)
	 * @param proj
	 */
	public static void createMarkerLookup2(Project proj) {
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] files, markerNames;
		byte[] readBuffer;
		RandomAccessFile currentFile;
		long time;

		time = new Date().getTime();
		System.out.println("Creating MarkerLookup file");
		files = new File(proj.getDir(Project.PLOT_DIRECTORY)).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".scat");
			}
		});
		if (files==null) {
			System.err.println("Error - failed to create MarkerLookup -- directory does not exist: "+proj.getDir(Project.PLOT_DIRECTORY));
		} else if (files.length==0) {
			System.err.println("Error - failed to create MarkerLookup -- no .scat files available");
		} else {
			for (int i = 0; i<files.length; i++) {
				try {
					currentFile = new RandomAccessFile(proj.getDir(Project.PLOT_DIRECTORY)+files[i], "r");
					readBuffer = new byte[currentFile.readInt()];
					currentFile.read(readBuffer);
					markerNames = new String(readBuffer).split("\t");
					for (int j = 0; j<markerNames.length; j++) {
						hash.put(markerNames[j], files[i]+"\t"+j);
					}
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}

			new MarkerLookup(hash).serialize(proj.getFilename(Project.MARKERLOOKUP_FILENAME));

			System.out.println("Created MarkerLookup in "+ext.getTimeElapsed(time));
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		Project proj;
		String filename = Project.DEFAULT_PROJECT;
		// String filename = "osteo_demo.proj";
		String demo = "";
		// String demo = "flagged_results.txt";
		String extract = "";
		boolean extractAll = false;
//		int numCycles = 0;
		int maxFileSize = 0;
//		int each = DEFAULT_MARKERS_PER_FILE;
		int perFile = 0;
		boolean lookup = false;
		boolean compact = false;

		String usage = "\n"+
		"filesys.ExtractPlots requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) extract all markers (i.e. -extractAll ("+(extractAll?"":"not the ")+"default))\n"+
//		"   (3) number of batches to perform extraction in (i.e. cycles="+numCycles+" (default))\n"+
		"   (3) maximum size of each file in bytes (i.e. cycles="+maxFileSize+" (default))\n"+
		"   (4) extract markers in file (i.e. extract=filename.txt (not the default))\n"+
		"   (5) extract markers in file to make demo (i.e. demo=filename.txt (not the default))\n"+
//		"   (6) change numper of markers per file (i.e. per="+perFile+" (default))\n"+
		"   (7) create marker lookup table (i.e. -lookup ("+(lookup?"":"not the ")+"default))\n"+
		"   (8) only parse X/Y (i.e. for Singleton) when extracting (i.e. -compact ("+(compact?"":"not the ")+"default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-extractAll")) {
				extractAll = true;
				numArgs--;
//			} else if (args[i].startsWith("cycles=")) {
			} else if (args[i].startsWith("filesize=")) {
				maxFileSize = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("extract=")) {
				extract = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("demo=")) {
				demo = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("per=")) {
				perFile = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-lookup")) {
				lookup = true;
				numArgs--;
			} else if (args[i].startsWith("-compact")) {
				compact = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		proj = new Project(filename, false);

		if (!proj.getDir(Project.SOURCE_DIRECTORY).equals("")&&!new File(proj.getDir(Project.SOURCE_DIRECTORY)).exists()) {
			System.err.println("Error - the project source location is invalid: "+proj.getDir(Project.SOURCE_DIRECTORY));
			return;
		}

//		extractAll = true;
//		compact = true;
		lookup = true;
		try {
			if (perFile>0) {
//				breakUpMarkerCollections(proj, perFile);
			} else if (!demo.equals("")) {
				extractForDemo(proj, demo);
			} else if (!extract.equals("")) {
				extract(proj, extract, compact);
			} else if (extractAll) {
//				extractAll(proj, numCycles, compact);
				transposeData(proj, maxFileSize, true);
			} else if (lookup) {
				createMarkerLookup(proj);
			} else {
				System.err.println("Error - no selection was made");
				System.out.println(usage);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


	/*
		public static void extractAll(Project proj, int numCycles, boolean compact_XYonly) {
			BufferedReader reader;
			PrintWriter writer;
			SetOfSets[] setsOfSets;
			String[] markerNames;
			SampleList sampleList;
			MarkerSet markerSet;
			// Thread[] threads;
			long fingerprint;
			int[][] sets;
			int count;
			// boolean done;
			long heapMaxSize;
			int bytesPerMarker, numMarkersPerFile;
			boolean done;
	
			done = false;
			while (!done) {
				try {
					try {
						writer = new PrintWriter(new FileWriter(".current"));
						writer.println(proj.getProperty(Project.PROJECT_PROPERTIES_FILENAME));
						writer.println(numCycles);
						writer.close();
					} catch (Exception e) {
						System.err.println("Error writing backup file");
					}
	
					sampleList = proj.getSampleList();
					if (sampleList==null) {
						System.out.println("No SampleList file found; generating one using all the samples in: "+proj.getDir(Project.SAMPLE_DIRECTORY));
						sampleList = SampleList.generateSampleList(proj);
					}
					fingerprint = sampleList.getFingerprint();
	
					markerSet = proj.getMarkerSet();
					markerNames = markerSet.getMarkerNames();
	
					if (numCycles==0) {
						bytesPerMarker = 37*sampleList.getSamples().length+161;
						heapMaxSize = Runtime.getRuntime().maxMemory();
						numMarkersPerFile = (int)(((double)heapMaxSize*0.75)/bytesPerMarker);
						System.out.println("Heap size (memory) was set to "+ext.formDeci((double)heapMaxSize/1024/1024/1024, 2)+" Gb; which allows us to do "+numMarkersPerFile+" markers at a time comfortably");
						System.out.println("Suggestion: Use 90% of available memory (i.e. -Xmx15g if you have 16GB of memory available) and use a 64-bit version of java (and the -d64 option) if you plan to use more than 2GB of memory");
						numCycles = (int)Math.ceil((double)markerNames.length/(double)numMarkersPerFile);
					}
					numMarkersPerFile = (int)Math.ceil((double)markerNames.length/(double)numCycles);
					System.out.println("Max memory: "+Runtime.getRuntime().maxMemory());
					Runtime.getRuntime().gc();
					Runtime.getRuntime().runFinalization();
					Runtime.getRuntime().gc();
					System.out.println("Running using "+numCycles+" cycles ("+numMarkersPerFile+" markers per run)");
	
					count = 0;
					sets = new int[numCycles][];
					for (int i = 0; i<sets.length; i++) {
						sets[i] = new int[Math.min(numMarkersPerFile, markerNames.length-count)];
						for (int j = 0; j<sets[i].length; j++) {
							sets[i][j] = count;
							count++;
						}
					}
	
					int numBatches = 1;
					setsOfSets = new SetOfSets[numBatches];
					for (int i = 0; i<numBatches; i++) {
						setsOfSets[i] = new SetOfSets();
					}
					for (int i = 0; i<sets.length; i++) {
						setsOfSets[i%numBatches].add(sets[i], proj.getDir(Project.PLOT_DIRECTORY)+"Markers."+ext.formNum(i+1, 4)+".scat");
					}
					new File(proj.getDir(Project.PLOT_DIRECTORY)).mkdirs();
	
					// threads = new Thread[numBatches];
					// for (int i = 0; i < numBatches; i++) {
					// threads[i] = new Thread(new ExtractPlots(proj,
					// setsOfSets[i].getSets(), fingerprint,
					// setsOfSets[i].getFilenames()));
					// threads[i].start();
					// try {
					// Thread.sleep(100L);
					// } catch(InterruptedException ex) {}
					// }
					// do {
					// done = true;
					// for (int i = 0; i < numBatches; i++) {
					// if (threads[i].isAlive()) {
					// done = false;
					// }
					// }
					// try {
					// Thread.sleep(10000L);
					// } catch(InterruptedException ex) {}
					// } while (!done);
	
					for (int i = 0; i<numBatches; i++) {
						extractPlots(proj, setsOfSets[i].getSets(), fingerprint, compact_XYonly, setsOfSets[i].getFilenames());
					}
	
					createMarkerLookup(proj);
					done = true;
				} catch (OutOfMemoryError oome) {
	
					System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
					Runtime.getRuntime().gc();
					Runtime.getRuntime().runFinalization();
					Runtime.getRuntime().gc();
					System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
	
					try {
						reader = new BufferedReader(new FileReader(".current"));
						while (reader.ready()) {
							proj = new Project(reader.readLine(), false);
							numCycles = Integer.parseInt(reader.readLine());
						}
						reader.close();
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \""+".current"+"\" not found in current directory");
						System.exit(1);
					} catch (IOException ioe) {
						System.err.println("Error reading file \""+".current"+"\"");
						System.exit(2);
					}
					System.err.println(ext.getTime()+"\tAlgorithm failed to run using "+numCycles+" cycles; will try using "+(numCycles+1));
					numCycles++;
				}
			}
		}
		*/
	
		/*
			public static void extractAll(Project proj, int numCycles, boolean compact_XYonly) {
				BufferedReader reader;
				PrintWriter writer;
				SetOfSets[] setsOfSets;
				String[] markerNames;
				SampleList sampleList;
				MarkerSet markerSet;
				// Thread[] threads;
				long fingerprint;
				int[][] sets;
				int count;
				// boolean done;
				long heapMaxSize;
				int bytesPerMarker, numMarkersPerFile;
				boolean done;
		
				done = false;
				while (!done) {
					try {
						try {
							writer = new PrintWriter(new FileWriter(".current"));
							writer.println(proj.getProperty(Project.PROJECT_PROPERTIES_FILENAME));
							writer.println(numCycles);
							writer.close();
						} catch (Exception e) {
							System.err.println("Error writing backup file");
						}
		
						sampleList = proj.getSampleList();
						if (sampleList==null) {
							System.out.println("No SampleList file found; generating one using all the samples in: "+proj.getDir(Project.SAMPLE_DIRECTORY));
							sampleList = SampleList.generateSampleList(proj);
						}
						fingerprint = sampleList.getFingerprint();
		
						markerSet = proj.getMarkerSet();
						markerNames = markerSet.getMarkerNames();
		
						if (numCycles==0) {
							bytesPerMarker = 37*sampleList.getSamples().length+161;
							heapMaxSize = Runtime.getRuntime().maxMemory();
							numMarkersPerFile = (int)(((double)heapMaxSize*0.75)/bytesPerMarker);
							System.out.println("Heap size (memory) was set to "+ext.formDeci((double)heapMaxSize/1024/1024/1024, 2)+" Gb; which allows us to do "+numMarkersPerFile+" markers at a time comfortably");
							System.out.println("Suggestion: Use 90% of available memory (i.e. -Xmx15g if you have 16GB of memory available) and use a 64-bit version of java (and the -d64 option) if you plan to use more than 2GB of memory");
							numCycles = (int)Math.ceil((double)markerNames.length/(double)numMarkersPerFile);
						}
						numMarkersPerFile = (int)Math.ceil((double)markerNames.length/(double)numCycles);
						System.out.println("Max memory: "+Runtime.getRuntime().maxMemory());
						Runtime.getRuntime().gc();
						Runtime.getRuntime().runFinalization();
						Runtime.getRuntime().gc();
						System.out.println("Running using "+numCycles+" cycles ("+numMarkersPerFile+" markers per run)");
		
						count = 0;
						sets = new int[numCycles][];
						for (int i = 0; i<sets.length; i++) {
							sets[i] = new int[Math.min(numMarkersPerFile, markerNames.length-count)];
							for (int j = 0; j<sets[i].length; j++) {
								sets[i][j] = count;
								count++;
							}
						}
		
						int numBatches = 1;
						setsOfSets = new SetOfSets[numBatches];
						for (int i = 0; i<numBatches; i++) {
							setsOfSets[i] = new SetOfSets();
						}
						for (int i = 0; i<sets.length; i++) {
							setsOfSets[i%numBatches].add(sets[i], proj.getDir(Project.PLOT_DIRECTORY)+"Markers."+ext.formNum(i+1, 4)+".scat");
						}
						new File(proj.getDir(Project.PLOT_DIRECTORY)).mkdirs();
		
						// threads = new Thread[numBatches];
						// for (int i = 0; i < numBatches; i++) {
						// threads[i] = new Thread(new ExtractPlots(proj,
						// setsOfSets[i].getSets(), fingerprint,
						// setsOfSets[i].getFilenames()));
						// threads[i].start();
						// try {
						// Thread.sleep(100L);
						// } catch(InterruptedException ex) {}
						// }
						// do {
						// done = true;
						// for (int i = 0; i < numBatches; i++) {
						// if (threads[i].isAlive()) {
						// done = false;
						// }
						// }
						// try {
						// Thread.sleep(10000L);
						// } catch(InterruptedException ex) {}
						// } while (!done);
		
						for (int i = 0; i<numBatches; i++) {
							extractPlots(proj, setsOfSets[i].getSets(), fingerprint, compact_XYonly, setsOfSets[i].getFilenames());
						}
		
						createMarkerLookup(proj);
						done = true;
					} catch (OutOfMemoryError oome) {
		
						System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
						Runtime.getRuntime().gc();
						Runtime.getRuntime().runFinalization();
						Runtime.getRuntime().gc();
						System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
		
						try {
							reader = new BufferedReader(new FileReader(".current"));
							while (reader.ready()) {
								proj = new Project(reader.readLine(), false);
								numCycles = Integer.parseInt(reader.readLine());
							}
							reader.close();
						} catch (FileNotFoundException fnfe) {
							System.err.println("Error: file \""+".current"+"\" not found in current directory");
							System.exit(1);
						} catch (IOException ioe) {
							System.err.println("Error reading file \""+".current"+"\"");
							System.exit(2);
						}
						System.err.println(ext.getTime()+"\tAlgorithm failed to run using "+numCycles+" cycles; will try using "+(numCycles+1));
						numCycles++;
					}
				}
			}
			*/


	/*
		public static void extractAll(Project proj, int numCycles, boolean compact_XYonly) {
			BufferedReader reader;
			PrintWriter writer;
			SetOfSets[] setsOfSets;
			String[] markerNames;
			SampleList sampleList;
			MarkerSet markerSet;
			// Thread[] threads;
			long fingerprint;
			int[][] sets;
			int count;
			// boolean done;
			long heapMaxSize;
			int bytesPerMarker, numMarkersPerFile;
			boolean done;
	
			done = false;
			while (!done) {
				try {
					try {
						writer = new PrintWriter(new FileWriter(".current"));
						writer.println(proj.getProperty(Project.PROJECT_PROPERTIES_FILENAME));
						writer.println(numCycles);
						writer.close();
					} catch (Exception e) {
						System.err.println("Error writing backup file");
					}
	
					sampleList = proj.getSampleList();
					if (sampleList==null) {
						System.out.println("No SampleList file found; generating one using all the samples in: "+proj.getDir(Project.SAMPLE_DIRECTORY));
						sampleList = SampleList.generateSampleList(proj);
					}
					fingerprint = sampleList.getFingerprint();
	
					markerSet = proj.getMarkerSet();
					markerNames = markerSet.getMarkerNames();
	
					if (numCycles==0) {
						bytesPerMarker = 37*sampleList.getSamples().length+161;
						heapMaxSize = Runtime.getRuntime().maxMemory();
						numMarkersPerFile = (int)(((double)heapMaxSize*0.75)/bytesPerMarker);
						System.out.println("Heap size (memory) was set to "+ext.formDeci((double)heapMaxSize/1024/1024/1024, 2)+" Gb; which allows us to do "+numMarkersPerFile+" markers at a time comfortably");
						System.out.println("Suggestion: Use 90% of available memory (i.e. -Xmx15g if you have 16GB of memory available) and use a 64-bit version of java (and the -d64 option) if you plan to use more than 2GB of memory");
						numCycles = (int)Math.ceil((double)markerNames.length/(double)numMarkersPerFile);
					}
					numMarkersPerFile = (int)Math.ceil((double)markerNames.length/(double)numCycles);
					System.out.println("Max memory: "+Runtime.getRuntime().maxMemory());
					Runtime.getRuntime().gc();
					Runtime.getRuntime().runFinalization();
					Runtime.getRuntime().gc();
					System.out.println("Running using "+numCycles+" cycles ("+numMarkersPerFile+" markers per run)");
	
					count = 0;
					sets = new int[numCycles][];
					for (int i = 0; i<sets.length; i++) {
						sets[i] = new int[Math.min(numMarkersPerFile, markerNames.length-count)];
						for (int j = 0; j<sets[i].length; j++) {
							sets[i][j] = count;
							count++;
						}
					}
	
					int numBatches = 1;
					setsOfSets = new SetOfSets[numBatches];
					for (int i = 0; i<numBatches; i++) {
						setsOfSets[i] = new SetOfSets();
					}
					for (int i = 0; i<sets.length; i++) {
						setsOfSets[i%numBatches].add(sets[i], proj.getDir(Project.PLOT_DIRECTORY)+"Markers."+ext.formNum(i+1, 4)+".scat");
					}
					new File(proj.getDir(Project.PLOT_DIRECTORY)).mkdirs();
	
					// threads = new Thread[numBatches];
					// for (int i = 0; i < numBatches; i++) {
					// threads[i] = new Thread(new ExtractPlots(proj,
					// setsOfSets[i].getSets(), fingerprint,
					// setsOfSets[i].getFilenames()));
					// threads[i].start();
					// try {
					// Thread.sleep(100L);
					// } catch(InterruptedException ex) {}
					// }
					// do {
					// done = true;
					// for (int i = 0; i < numBatches; i++) {
					// if (threads[i].isAlive()) {
					// done = false;
					// }
					// }
					// try {
					// Thread.sleep(10000L);
					// } catch(InterruptedException ex) {}
					// } while (!done);
	
					for (int i = 0; i<numBatches; i++) {
						extractPlots(proj, setsOfSets[i].getSets(), fingerprint, compact_XYonly, setsOfSets[i].getFilenames());
					}
	
					createMarkerLookup(proj);
					done = true;
				} catch (OutOfMemoryError oome) {
	
					System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
					Runtime.getRuntime().gc();
					Runtime.getRuntime().runFinalization();
					Runtime.getRuntime().gc();
					System.out.println("Free memory: "+Runtime.getRuntime().freeMemory());
	
					try {
						reader = new BufferedReader(new FileReader(".current"));
						while (reader.ready()) {
							proj = new Project(reader.readLine(), false);
							numCycles = Integer.parseInt(reader.readLine());
						}
						reader.close();
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \""+".current"+"\" not found in current directory");
						System.exit(1);
					} catch (IOException ioe) {
						System.err.println("Error reading file \""+".current"+"\"");
						System.exit(2);
					}
					System.err.println(ext.getTime()+"\tAlgorithm failed to run using "+numCycles+" cycles; will try using "+(numCycles+1));
					numCycles++;
				}
			}
		}
		*/
}
