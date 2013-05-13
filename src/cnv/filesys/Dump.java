package cnv.filesys;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.Date;

import cnv.manage.TransposeData;

public class Dump {
	public static void dump(String filename) {
        PrintWriter writer;
        
    	try {
	        if (filename.endsWith(Sample.SAMPLE_DATA_FILE_EXTENSION)) {
	        	Sample samp = Sample.loadFromRandomAccessFile(filename, false);
	        	float[] gcs = samp.getGCs();
	        	float[] xs = samp.getXs();
	        	float[] ys = samp.getYs();
	        	float[] thetas = samp.getThetas();
	        	float[] rs = samp.getRs();
	        	float[] lrrs = samp.getLRRs();
	        	float[] bafs = samp.getBAFs();
	        	byte[] forwardGenotypes = samp.getForwardGenotypes();
	        	byte[] abGenotypes = samp.getAB_Genotypes();

	        	writer = new PrintWriter(new FileWriter(filename + "_" + samp.getSampleName() + "_" + samp.getFingerprint() + ".xln"));
	        	writer.println( (xs==null? "" : "X")
	        				  + (ys==null? "" : "\tY")
			        		  + (thetas==null? "" : "\tTheta")
			        		  + (rs==null? "" : "\tR")
			        		  + (bafs==null? "" : "\tBAF")
			        		  + (lrrs==null? "" : "\tLRR")
			        		  + (gcs==null? "" : "\tGC_score")
			        		  + (abGenotypes==null? "" : "\tAB_Genotypes")
			        		  + (forwardGenotypes==null? "" : "\tForward_Genotypes"));
	        	for (int i = 0; i<xs.length; i++) {
	        		writer.println( (xs==null? "" : xs[i])
	        					  + (ys==null? "" : "\t" + ys[i])
	        					  + (thetas==null? "" : "\t" + thetas[i])
	        					  + (rs==null? "" : "\t" + rs[i])
	        					  + (bafs==null? "" : "\t" + bafs[i])
	        					  + (lrrs==null? "" : "\t" + lrrs[i])
	        					  + (gcs==null? "" : "\t" + gcs[i])
	        					  + (abGenotypes==null? "" : "\t" + (abGenotypes[i]==-1? "--" : Sample.AB_PAIRS[abGenotypes[i]]))
	        					  + (forwardGenotypes==null?"" : "\t" + Sample.ALLELE_PAIRS[forwardGenotypes[i]]));
                }
		        writer.close();

	        } else if (filename.endsWith(".bim")) {
	        	MarkerSet set = MarkerSet.load(filename, false);
	        	String[] markerNames = set.getMarkerNames();
	        	byte[] chrs = set.getChrs();
	        	int[] positions = set.getPositions();

	            writer = new PrintWriter(new FileWriter(filename+"_"+set.getFingerprint()+".xln"));
	            writer.println("MarkerName\tChr\tPosition");
	        	for (int i = 0; i<markerNames.length; i++) {
	        		writer.println(markerNames[i]+"\t"+chrs[i]+"\t"+positions[i]);
                }
		        writer.close();

	        } else if (filename.endsWith(MarkerData.MARKER_DATA_FILE_EXTENSION)) {
	        	MarkerData[] mkData;
	        	int indexStartMarker, indexEndMarker;
	        	
	        	indexStartMarker = 0;
	        	indexEndMarker = 0;
	        	mkData = TransposeData.loadFromRAF(filename, indexStartMarker, indexEndMarker);
	        	for (int i=0; i<mkData.length; i++) {
	        		mkData[i].dump(filename + "_dump_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".xln", null);
	        	}
	        }
        } catch (Exception e) {
            System.err.println("Error dumping data from "+filename+" to a textfile");
            e.printStackTrace();
        }
	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    String filename = "sample.fsamp";

	    String usage = "\n"+
	    "cnv.filesys.Dump requires 0-1 arguments\n"+
	    "   (1) filename (i.e. file="+filename+" (default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("file=")) {
			    filename = args[i].split("=")[1];
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
		    dump(filename);
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
