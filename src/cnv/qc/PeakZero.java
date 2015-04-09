package cnv.qc;

import java.io.*;
import java.util.*;

import common.*;
import cnv.filesys.*;
import cnv.manage.MarkerDataLoader;
import stats.Histogram;

public class PeakZero {
	public static final String PEAK_ZERO_FILE = "PeakZeroCalculations.xln";
//	public static final String[] DUMP_THESE = {"rs9786160", "rs1050917", "cnv336p1", "cnv1988p3", "cnv9025PP4", "rs9259925", "cnv14145p2", "cnv14229p1"};
//	public static final String[] DUMP_THESE = {"rs28705211", "rs11721", "rs2649588", "cnv22p4"};
//	public static final String[] DUMP_THESE = {"rs12121864", "rs34620016"};
	public static final String[] DUMP_THESE = {"rs277452", "rs2968487"};
//	public static final String[] DUMP_THESE = null;
		
	public static void checkDists(Project proj, String phenoOfSamplesToInclude) {
        PrintWriter writer;
        String trav;
        Hashtable<String,String> hash, drops;
        int count;
		String[] samples;
		boolean[] use;
		float[] lrrs, bafs, lrrArray, bafArray;
		float[] xs, ys, xArray, yArray;
		Histogram lrrHist, bafHist, xHist, yHist;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		String[] markerNames;
		long time;
		
		samples = proj.getSamples();
		if (!phenoOfSamplesToInclude.equals("")) {
			hash = HashVec.loadFileToHashString(proj.getDir(proj.SAMPLE_DATA_FILENAME), "DNA", new String[] {phenoOfSamplesToInclude}, "");
			use = new boolean[samples.length];
			for (int i = 0; i<samples.length; i++) {
				trav = hash.get(samples[i]);
				if (trav == null) {
					System.err.println("Error - '"+samples[i]+"' was not found in a column with the header 'DNA'");
					use[i] = false;
				} else {
					use[i] = !(trav.equals(".") || trav.equals("NA"));
				}
            }
			System.out.println("Distributions will be created from "+Array.booleanArraySum(use)+" of "+samples.length+" possible samples");
		} else {
			use = Array.booleanArray(samples.length, true);
			System.out.println("Distributions will be created using all "+samples.length+" samples");
		}
		count = Array.booleanArraySum(use);
		
		drops = proj.getFilteredHash();
		
		try {
			new File(proj.getDir(proj.RESULTS_DIRECTORY)).mkdirs();
	        writer = new PrintWriter(new FileWriter(proj.getDir(proj.RESULTS_DIRECTORY)+PEAK_ZERO_FILE));
	        writer.println("Marker\tPeakOffset\t#LRR_Maxima\t#BAF_Maxima\t#X_Maxima\t#Y_Maxima\tDropped");

	        lrrArray = new float[count];
	        bafArray = new float[count];
	        xArray = new float[count];
	        yArray = new float[count];

	        time = new Date().getTime();
			markerNames = proj.getMarkerNames();
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
			for (int i = 0; i < markerNames.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);
				if (i % 100 == 0) {
					System.out.println(ext.getTime()+"\tMarker "+i+" of "+markerNames.length);
				}

	        	count = 0;
	        	lrrs = markerData.getLRRs();
	        	bafs = markerData.getBAFs();
	        	xs = markerData.getXs();
	        	ys = markerData.getYs();
	        	for (int k = 0; k<samples.length; k++) {
	        		if (use[k]) {
	        			lrrArray[count] = lrrs[k];
	        			bafArray[count] = bafs[k];
	        			xArray[count] = xs[k];
	        			yArray[count] = ys[k];
	        			count++;
	        		}
                }
	        	lrrHist = new Histogram(lrrArray, -10, 10, 1);
	        	bafHist = new Histogram(bafArray, 0, 1, 2);
	        	xHist = new Histogram(lrrArray, 0, 5, 1);
	        	yHist = new Histogram(bafArray, 0, 5, 1);

	        	if (DUMP_THESE != null && ext.indexOfStr(markerData.getMarkerName(), DUMP_THESE) >= 0) {
	        		lrrHist.dump(proj.getProjectDir()+markerData.getMarkerName()+"_lrr_hist.xln");
	        		bafHist.dump(proj.getProjectDir()+markerData.getMarkerName()+"_baf_hist.xln");
	        	}
	        	
	        	writer.println(markerData.getMarkerName()+"\t"+lrrHist.getMaxBin()+"\t"+lrrHist.getLocalMaxima().length+"\t"+bafHist.getLocalMaxima().length+"\t"+xHist.getLocalMaxima().length+"\t"+yHist.getLocalMaxima().length+"\t"+(drops.containsKey(markerData.getMarkerName())?1:0));
	        	writer.flush();
				markerDataLoader.releaseIndex(i);
            }
				
			System.out.println("Finished "+markerNames.length+" in "+ext.getTimeElapsed(time));

	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+proj.getDir(proj.RESULTS_DIRECTORY)+PEAK_ZERO_FILE);
	        e.printStackTrace();
        }
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		boolean check = true;

		String usage = "\\n"+
		"cnv.qc.PeakZero requires 0-1 arguments\n"+
		"   (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"   (2) check distributions (i.e. -check (not the default))\n"+
		"";


		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-check")) {
				check = true;
				numArgs--;
			}
		}

		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		try {
			if (check) {
				checkDists(new Project(filename, false), "");
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}