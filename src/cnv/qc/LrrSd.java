package cnv.qc;

import java.io.*;

import cnv.filesys.*;
//import java.util.*;
import common.*;

public class LrrSd extends Parallelizable {
	private Project proj;
	private String[] samples;
	private String centroidsFile;
	private int threadNumber;
	private int numThreads;

	public LrrSd(Project proj, String[] samples, String centroidsFile, int threadNumber, int numThreads) {
		this.proj = proj;
		this.samples = samples;
		this.centroidsFile = centroidsFile;
		this.threadNumber = threadNumber;
		this.numThreads = numThreads;
	}
	
	public void run() {
		PrintWriter writer;
		FullSample fsamp;
		float[][][] cents;
		byte[] chrs;
		float[] lrrs;

		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"lrr_sd."+threadNumber));

			if (centroidsFile==null) {
				cents = null;
			} else {
				cents = Centroids.load(centroidsFile, false).getCentroids(); 
			}
			
			for (int i = 0; i<samples.length; i++) {
	        	System.out.println((i+1)+" of "+samples.length);
				fsamp = proj.getFullSample(samples[i]);
				chrs = proj.getMarkerSet().getChrs();
				lrrs = null;
				if (fsamp == null) {
					System.err.println("Error - "+samples[i]+".fsamp not found in samples directory");
					lrrs =  proj.getSample(samples[i]).getLRRs();
				} else {
					lrrs = cents==null?fsamp.getLRRs():fsamp.getLRRs(cents);
					lrrs = Array.subArray(lrrs, 0, Array.indexOfByte(chrs, (byte)23));
				}
				if (lrrs == null) {
					System.err.println("Error - could not find "+samples[i]+".fsamp or "+samples[i]+".samp");
				} else {
					writer.println(samples[i]+"\t"+Array.stdev(lrrs, true));
					writer.flush();
				}
			}
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void finalAction() {
		String[] files; 
		
		files = Array.stringArraySequence(numThreads, proj.getProjectDir()+"lrr_sd.");
		Files.cat(files, proj.getProjectDir()+"lrr_sd.xln", Array.intArray(files.length, 0), new Logger(null));
		for (int i = 0; i<files.length; i++) {
			new File(files[i]).delete();
        }
	}
	
	public static void init(Project proj, String centroidsFile, int numThreads) {
		String[] samples = proj.getSamples();
		String[][] threadSeeds;
		LrrSd[] runables;
		
		threadSeeds = Parallelizable.splitList(samples, numThreads, false);
		runables = new LrrSd[numThreads];
		for (int i = 0; i<numThreads; i++) {
			runables[i] = new LrrSd(proj, threadSeeds[i], centroidsFile, i+1, numThreads);
        }
		
		Parallelizable.launch(runables);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;
		String centroids = null;
		int numThreads = 1;
		Project proj;

		String usage = "\n"+
		"cnv.qc.LrrSd requires 0-4 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) centroids with which to compute LRRs (i.e. cents=genotype.cent (not the default; to be found in data/ directory))\n"+
		"   (3) number of threads to use (i.e. threads="+numThreads+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cents=")) {
				centroids = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			proj = new Project(filename, false);

			init(proj, centroids, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
