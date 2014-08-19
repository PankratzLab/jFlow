package cnv.qc;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;
import common.Logger;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerFreqs;
import cnv.filesys.MarkerLookup;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.manage.MarkerDataLoader;
import cnv.var.SampleData;

//class to compute marker QCs from all markers, currently only supports a simple minor allele frequency computation
//TODO actual loggers
public class CNVMarkerQC implements Runnable {
	Project proj;
	String[] markerNames;
	boolean[] samplesToBeUsed;
	double[] mafs;
	double [] bafs;

	public CNVMarkerQC(Project proj, String[] markerNames, boolean[] samplesToBeUsed) {
		this.proj = proj;
		this.markerNames = markerNames;
		this.samplesToBeUsed = samplesToBeUsed;
		this.mafs = new double[markerNames.length];

	}

	public String[] getMarkerNames() {
		return markerNames;
	}

	public double[] getMafs() {
		return mafs;
	}

	public void run() {
		// TODO is loading in a separate thread within a separate thread a good thing or bad
		// TODO add sex specific calculation, and include/exclude from sampleData

		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
		for (int i = 0; i < markerNames.length; i++) {
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			mafs[i] = markerData.getMAF(samplesToBeUsed, null, null, 0);
			markerDataLoader.releaseIndex(i);
		}
	}

	public static void computeMAFs(Project proj, int threads, boolean excludeSamples, String outputSer) {
		MarkerLookup markerLookup = proj.getMarkerLookup();
		MarkerSet markerSet = proj.getMarkerSet();
		SampleData sampleData = proj.getSampleData(2, false);
		String[] markerNames = markerSet.getMarkerNames();
		Hashtable<String, ArrayList<String>> markerFiles = assignFiles(markerLookup, markerNames);
		ArrayList<ArrayList<String>> cabinet = getcabinet(markerFiles, threads);
		boolean[] samplesToBeUsed = getSamplesToBeUsed(proj, excludeSamples, sampleData);
		CNVMarkerQC[] markerFrequencies = computeFileMAFS(proj, threads, cabinet, samplesToBeUsed);
		double[] mafs = summarize(proj, markerFrequencies, markerNames);
		new MarkerFreqs(mafs, markerSet.getFingerprint()).serialize(proj.getProjectDir() + outputSer);
	}


	private static boolean[] getSamplesToBeUsed(Project proj, boolean excludeSamples, SampleData sampleData) {
		boolean doExclude = shouldExclude(proj, excludeSamples, sampleData);
		if (!doExclude) {
			return null;
		} else {
			String[] listofSamples = sampleData.getListOfSamples();
			boolean[] samplesToBeUsed = new boolean[listofSamples.length];
			for (int i = 0; i < listofSamples.length; i++) {
				samplesToBeUsed[i] = !sampleData.individualShouldBeExcluded(listofSamples[i]);
			}
			return samplesToBeUsed;
		}
	}

	private static boolean shouldExclude(Project proj, boolean excludeSamples, SampleData sampleData) {
		boolean doExclude = false;
		if (!excludeSamples) {
			doExclude = false;
		} else if (sampleData != null && sampleData.hasExcludedIndividuals()) {
			doExclude = true;
		} else {
			proj.getLog().reportError("Error - cannot exclude individuals for MAF computation, no factor named 'Exclude/CLASS=Exclude' in Sample Data");
			System.exit(1);
		}
		return doExclude;
	}

	private static double[] summarize(Project proj, CNVMarkerQC[] markerFrequencies, String[] markerNames) {
		Hashtable<String, Double> mafBuilder = new Hashtable<String, Double>();
		double[] mafs = new double[markerNames.length];
		for (int i = 0; i < markerFrequencies.length; i++) {
			double[] fileMAFs = markerFrequencies[i].getMafs();
			String[] fileMarkerNames = markerFrequencies[i].getMarkerNames();
			for (int k = 0; k < fileMarkerNames.length; k++) {
				mafBuilder.put(fileMarkerNames[k], fileMAFs[k]);
			}
		}
		for (int i = 0; i < markerNames.length; i++) {
			if (mafBuilder.containsKey(markerNames[i])) {
				mafs[i] = mafBuilder.get(markerNames[i]);
			} else {
				proj.getLog().reportError("Error - could not find marker " + markerNames[i] + ", this should not happen");
				System.exit(1);
			}
		}
		return mafs;
	}

	private static CNVMarkerQC[] computeFileMAFS(Project proj, int threads, ArrayList<ArrayList<String>> cabinet, boolean[] samplesToBeUsed) {
		CNVMarkerQC[] markerFrequencies = new CNVMarkerQC[threads];
		Thread[] runningthreads = new Thread[threads];
		for (int i = 0; i < threads; i++) {
			markerFrequencies[i] = new CNVMarkerQC(proj, toStringArray(cabinet.get(i)), samplesToBeUsed);
			runningthreads[i] = new Thread(markerFrequencies[i]);
			runningthreads[i].start();
		}
		checkThreadStatus(threads, runningthreads);
		return markerFrequencies;
	}

	private static void checkThreadStatus(int processors, Thread[] threads) {
		boolean complete;
		complete = false;
		while (!complete) {
			complete = true;
			for (int i = 0; i < processors; i++) {
				if (threads[i].isAlive()) {
					complete = false;
				}
			}
			if (!complete) {
				try {
					Thread.sleep(1000L);
				} catch (InterruptedException ex) {
				}
			}
		}
	}

	private static ArrayList<ArrayList<String>> getcabinet(Hashtable<String, ArrayList<String>> markerFiles, int threads) {
		ArrayList<ArrayList<String>> cabinet = new ArrayList<ArrayList<String>>();
		Set<String> keys = markerFiles.keySet();
		for (int i = 0; i < threads; i++) {
			cabinet.add(new ArrayList<String>());
		}
		int elementAt = 0;
		for (String key : keys) {
			cabinet.get(elementAt % threads).addAll(markerFiles.get(key));
			elementAt++;
		}
		return cabinet;
	}

	private static Hashtable<String, ArrayList<String>> assignFiles(MarkerLookup markerLookup, String[] markerNames) {
		Hashtable<String, ArrayList<String>> markerFiles = new Hashtable<String, ArrayList<String>>();
		for (int i = 0; i < markerNames.length; i++) {
			if (!markerFiles.containsKey(markerLookup.get(markerNames[i]).split("[\\s]+")[0])) {
				ArrayList<String> al = new ArrayList<String>();
				markerFiles.put(markerLookup.get(markerNames[i]).split("[\\s]+")[0], al);
			}
			markerFiles.get(markerLookup.get(markerNames[i]).split("[\\s]+")[0]).add(markerNames[i]);
		}
		return markerFiles;
	}

	private static String[] toStringArray(ArrayList<String> stringList) {
		return stringList.toArray(new String[stringList.size()]);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "C:/workspace/Genvisis/projects/ARICGenvisis_CEL_11908.properties";
		int numThreads = 4;
		boolean excludeSamples = false;
		boolean convert = false;
		String markerFreqSer = "MarkerFreq.ser";
		String markerFreqTxt = "MarkerFreq.txt";
		String usage = "\n" + "cnv.filesys.MarkerFrequencies requires 0-4 arguments\n" + "   (1) project properties filename (i.e. proj=" + filename + " (default))\n" + " OPTIONAL:\n" + "   (2) number of threads to use (i.e. threads=" + numThreads + " (default))\n" + "   (3) name of serialized  MarkerFreq file (i.e ser=" + markerFreqSer + " (default))\n" + "   (4) name of text MarkerFreq file (i.e txt=" + markerFreqTxt + " (default))\n" + "   (5) exclude samples as defined in sampleData (i.e -exclude (not the default ,currently not supported))\n" + "   (6) convert from an existing .txt file (i.e -convert (not the default))\n";
		// String usage = "\n"+
		// "cnv.filesys.MarkerFrequencies requires 0-4 arguments\n"+
		// "   (1) project properties filename (i.e. proj="+filename+" (default))\n"+
		// " OPTIONAL:\n"+
		// "   (2) number of threads to use (i.e. threads="+numThreads+" (default))\n"+
		// "   (3) name of serialized  MarkerFreq file (i.e ser="+markerFreqSer+" (default))\n"+
		// "   (4) name of text MarkerFreq file (i.e txt="+markerFreqTxt+" (default))\n"+
		// "   (5) exclude samples as defined in sampleData (i.e -exclude (not the default ,currently not supported))\n"+
		// "   (6) convert from an existing .txt file (i.e -convert (not the default))\n";
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				return;
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("ser=")) {
				markerFreqSer = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("txt=")) {
				markerFreqTxt = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-exclude")) {
				excludeSamples = true;
				numArgs--;
			} else if (args[i].startsWith("-convert")) {
				convert = true;
				numArgs--;
			}

		}
		if (numArgs != 0) {
			System.err.println(usage);
			return;
		}
		try {
			Project proj = new Project(filename, false);
			proj.setLog(new Logger(proj.getProjectDir() + "CNVMarkerQCLog.txt"));
			if (!convert) {
				computeMAFs(proj, numThreads, excludeSamples, markerFreqSer);
			// currently always output .txt format as well since it is small
				MarkerFreqs.exportToText(filename, markerFreqSer, markerFreqTxt);
			} else {
				MarkerFreqs.convertMarkerFreqsFromTxt(proj, markerFreqTxt, markerFreqSer);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}