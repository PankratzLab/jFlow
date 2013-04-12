// -Xms1024M -Xmx1024M     or even better: -Xmx15g
package cnv.manage;

import java.io.*;
import java.lang.OutOfMemoryError;
import java.util.*;

import cnv.filesys.FullSample;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerDataCollection;
import cnv.filesys.MarkerLookup;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.SampleList;
import cnv.filesys.Centroids;
import cnv.var.SetOfSets;

import common.*;

//public class ExtractPlots implements Runnable {
public class ExtractPlots {
	public static final int DEFAULT_MARKERS_PER_FILE = 210000;

	// private Project proj;
	// private int[][] indicesSets;
	// private long fingerprint;
	// private String[] outputFiles;
	//
	// public ExtractPlots(Project proj, int[][] indicesSets, long fingerprint,
	// String[] outputFiles) {
	// this.proj = proj;
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
	public static void extractPlots(Project proj, int[][] indicesSets, long fingerprint, boolean compact_XYonly, String[] outputFiles) {
		for (int i = 0; i<indicesSets.length; i++) {
			sampleToSNPmajor(proj, indicesSets[i], fingerprint, compact_XYonly, outputFiles[i]);
		}
	}

	public static void sampleToSNPmajor(Project proj, int[] indices, long fingerprint, boolean compact_XYonly, String outputFile) {
		SampleList list;
		MarkerSet markerSet;
		String[] samples, allMarkers, markerNames;
		FullSample samp;
		MarkerData[] collection;
		float[][] gcs, xs, ys, xRaws, yRaws, thetas, rs, lrrs, bafs;
		byte[][] abGenotypes;
		float[] samp_gcs, samp_xs, samp_ys, samp_xRaws, samp_yRaws, samp_thetas, samp_rs, samp_lrrs, samp_bafs;
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
			xRaws = new float[markerNames.length][];
			yRaws = new float[markerNames.length][];
			thetas = new float[markerNames.length][];
			rs = new float[markerNames.length][];
			lrrs = new float[markerNames.length][];
			bafs = new float[markerNames.length][];
		} else {
			xRaws = new float[markerNames.length][samples.length];
			yRaws = new float[markerNames.length][samples.length];
			thetas = new float[markerNames.length][samples.length];
			rs = new float[markerNames.length][samples.length];
			lrrs = new float[markerNames.length][samples.length];
			bafs = new float[markerNames.length][samples.length];
		}
		for (int i = 0; i<samples.length; i++) {
			samp = proj.getFullSample(samples[i]);
			samp_gcs = samp.getGCs();
			samp_xs = samp.getXs();
			samp_ys = samp.getYs();
			samp_abGenotypes = samp.getAB_Genotypes();
			samp_forwardGenotypes = samp.getForwardGenotypes();
			samp_xRaws = samp.getX_Raws();
			samp_yRaws = samp.getY_Raws();
			samp_thetas = samp.getThetas();
			samp_rs = samp.getRs();
			samp_lrrs = samp.getLRRs();
			samp_bafs = samp.getBAFs();
			for (int m = 0; m<markerNames.length; m++) {
				gcs[m][i] = samp_gcs[indices[m]];
				xs[m][i] = samp_xs[indices[m]];
				ys[m][i] = samp_ys[indices[m]];
				if (!compact_XYonly) {
					xRaws[m][i] = samp_xRaws[indices[m]];
					yRaws[m][i] = samp_yRaws[indices[m]];
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
						alleleMappings[m][samp_abGenotypes[indices[m]]+1] = FullSample.ALLELE_PAIRS[samp_forwardGenotypes[indices[m]]];
					} else if (!alleleMappings[m][samp_abGenotypes[indices[m]]+1].equals(FullSample.ALLELE_PAIRS[samp_forwardGenotypes[indices[m]]])) {
						System.err.println("Error - mismatched alleles for marker "+markerNames[m]);
					}
				}
			}
		}
		collection = new MarkerData[markerNames.length];
		for (int i = 0; i<markerNames.length; i++) {
			collection[i] = new MarkerData(markerNames[i], chrs[i], positions[i], fingerprint, gcs[i], xRaws[i], yRaws[i], xs[i], ys[i], thetas[i], rs[i], bafs[i], lrrs[i], abGenotypes[i], alleleMappings[i]);
		}
		new MarkerDataCollection(outputFile, fingerprint, markerNames, collection).serialize();

		System.out.println(ext.getTime()+" "+ext.getDate()+" (finished in "+ext.getTimeElapsed(time)+")");
	}

	public static String translateName(String str) {
		return ext.replaceAllWith(str, ":", ".");
	}

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
					if (sampleList==null) {
						return;
					}
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

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		Project proj;
		String filename = Project.DEFAULT_PROJECT;
		// String filename = "osteo_demo.proj";
		String demo = "";
		// String demo = "flagged_results.txt";
		String extract = "";
		boolean extractAll = false;
		int numCycles = 0;
		// int each = DEFAULT_MARKERS_PER_FILE;
		int perFile = 0;
		boolean lookup = false;
		boolean compact = false;

		String usage = "\n"+
		"filesys.ExtractPlots requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) extract all markers (i.e. -extractAll ("+(extractAll?"":"not the ")+"default))\n"+
		"   (3) number of batches to perform extraction in (i.e. cycles="+numCycles+" (default))\n"+
		"   (4) extract markers in file (i.e. extract=filename.txt (not the default))\n"+
		"   (5) extract markers in file to make demo (i.e. demo=filename.txt (not the default))\n"+
		"   (6) change numper of markers per file (i.e. per="+perFile+" (default))\n"+
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
			} else if (args[i].startsWith("cycles=")) {
				numCycles = Integer.parseInt(args[i].split("=")[1]);
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
		try {
			if (perFile>0) {
				breakUpMarkerCollections(proj, perFile);
			} else if (!demo.equals("")) {
				extractForDemo(proj, demo);
			} else if (!extract.equals("")) {
				extract(proj, extract, compact);
			} else if (extractAll) {
				extractAll(proj, numCycles, compact);
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
}
