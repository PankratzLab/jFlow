package cnv.analysis;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;

import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.manage.MarkerDataLoader;

import common.Array;
import common.HashVec;
import common.Logger;
import common.ext;

public class Zcall {
	public static final String[] BASIC_HEADER = {"Name", "Chr", "Position"};
	public static final String[] AB_GENOTYPES = {"NC", "AA", "AB", "BB"};

	public static void createZcallInputFile(Project proj, String filenameOfSamplesToInclude, String markersToInclude, Logger log) {
		HashSet<String> hash;
		
		boolean[] samplesToInclude;
		String[] samples;

		samples = proj.getSamples();
		samplesToInclude = new boolean[samples.length];
		
		hash = HashVec.loadFileToHashSet(filenameOfSamplesToInclude, false);
		
		for (int i = 0; i < samples.length; i++) {
			samplesToInclude[i] = hash.contains(samples[i]);
		}

		createZcallInputFile(proj, samplesToInclude, markersToInclude, log);
	}
	
	public static void createZcallInputFile(Project proj, boolean[] samplesToInclude, String markersToInclude, Logger log) {
		PrintWriter writer;
		String[] samples;
		float[] xs, ys;
		MarkerData markerData;
        byte[] abGenotypes;
        String markerName;
        ClusterFilterCollection clusterFilterCollection;
        float gcThreshold;
        long time;
        MarkerDataLoader markerDataLoader;
        String[] markerNames;
        String eol;

        if (System.getProperty("os.name").startsWith("Windows")) {
        	eol = "\r\n";
		} else {
			eol = "\n";
		}
        
        samples = proj.getSamples();
        if (samplesToInclude == null) {
        	samplesToInclude = Array.booleanArray(samples.length, true);
        } else if (samplesToInclude.length != samples.length) {
        	log.reportError("Error - length of samplesToInclude does not match length of samples");
        	return;
        }
 		
        clusterFilterCollection = proj.getClusterFilterCollection();

        gcThreshold = Float.parseFloat(proj.getProperty(Project.GC_THRESHOLD));

		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir() + "GenomeStudioData.txt"));
			writer.print(Array.toStr(BASIC_HEADER));
			for (int i = 0; i < samples.length; i++) {
				if (samplesToInclude[i]) {
					writer.print("\t" + samples[i] + ".gtype\t" + samples[i] + ".X\t" + samples[i] + ".Y");
				}
			}
			writer.print(eol);
			
			if (markersToInclude != null) {
				markerNames = HashVec.loadFileToStringArray(markersToInclude, false, new int[] {0}, false);
			} else {
				markerNames = proj.getMarkerNames();
			}
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames, log);

			time = new Date().getTime();
			for (int i = 0; i < markerNames.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);
				if (i % 1000 == 0) {
					System.out.println(ext.getTime()+"\tMarker "+i+" of "+markerNames.length);
				}

				markerName = markerData.getMarkerName();
				
				xs = markerData.getXs();
				ys = markerData.getYs();
				abGenotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerName, gcThreshold);
				
				writer.print(markerName + "\t" + markerData.getChr() + "\t" + markerData.getPosition());

				for (int j = 0; j < samples.length; j++) {
					if (samplesToInclude[j]) {
						writer.print("\t" + AB_GENOTYPES[1 + abGenotypes[j]] + "\t" + xs[j] + "\t" + ys[j]);
					}
				}
				writer.print(eol);
				
				markerDataLoader.releaseIndex(i);
			}
			writer.close();
			log.report("Finished analyzing "+markerNames.length+" in "+ext.getTimeElapsed(time));
		} catch (Exception e) {
			log.reportError("Error writing marker metrics to "+proj.getFilename(Project.MARKER_METRICS_FILENAME, false, false));
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
//		String filename = Project.DEFAULT_PROJECT;
		String filename = "C:/workspace/Genvisis/projects/GEDI_exome.properties";
		String logfile = null;
		Logger log;
		Project proj;
		String samples = "sampleSubset.txt";
		String markersSubset = null;
		
		String usage = "\n" + 
		"cnv.analysis.Zcall requires 0-1 arguments\n" + 
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) filename of subset of markers to include / otherwise all markers (i.e. markers=" + markersSubset + " (default))\n" + 
		"   (3) filename of subset of samples to include / otherwise all samples (i.e. samples=" + samples + " (default))\n" +
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
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
			log = new Logger(logfile);
			proj = new Project(filename, false);
//			createZcallInputFile(proj, proj.getProjectDir() + samples, proj.getDir(Project.DATA_DIRECTORY) + "test.txt", log);
			createZcallInputFile(proj, proj.getProjectDir() + samples, null, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
