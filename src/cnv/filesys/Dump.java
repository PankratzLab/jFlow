package cnv.filesys;

import java.io.*;

public class Dump {
	public static void dump(String filename) {
        PrintWriter writer;
        
    	try {
	        if (filename.endsWith(".fsamp")) {
	        	FullSample fsamp = FullSample.load(filename, false);
	        	float[] gcs = fsamp.getGCs();
	        	float[] xs = fsamp.getXs();
	        	float[] ys = fsamp.getYs();
	        	float[] xRaws = fsamp.getX_Raws();
	        	float[] yRaws = fsamp.getY_Raws();
	        	float[] thetas = fsamp.getThetas();
	        	float[] rs = fsamp.getRs();
	        	float[] lrrs = fsamp.getLRRs();
	        	float[] bafs = fsamp.getBAFs();
	        	byte[] forwardGenotypes = fsamp.getForwardGenotypes();
	        	byte[] abGenotypes = fsamp.getAB_Genotypes();

	            writer = new PrintWriter(new FileWriter(filename+"_"+fsamp.getSampleName()+"_"+fsamp.getFingerprint()+".xln"));
	        	writer.println(
	        			(xRaws==null?"":"Raw_X\t")+
	        			(yRaws==null?"":"Raw_Y\t")+
	        			(xs==null?"":"X\t")+
	        			(ys==null?"":"Y\t")+
	        			(thetas==null?"":"Theta\t")+
	        			(rs==null?"":"R\t")+
	        			(bafs==null?"":"BAF\t")+
	        			(lrrs==null?"":"LRR\t")+
	        			(gcs==null?"":"GC_score\t")+
	        			(forwardGenotypes==null?"":"Forward_Genotypes\t")+
	        			(abGenotypes==null?"":"AB_Genotypes\t")
	        			);
	        	for (int i = 0; i<xs.length; i++) {
	        		writer.println(
	        				(xRaws==null?"":xRaws[i]+"\t")+
	        				(yRaws==null?"":yRaws[i]+"\t")+
	        				(xs==null?"":xs[i]+"\t")+
	        				(ys==null?"":ys[i]+"\t")+
	        				(thetas==null?"":thetas[i]+"\t")+
	        				(rs==null?"":rs[i]+"\t")+
	        				(bafs==null?"":bafs[i]+"\t")+
	        				(lrrs==null?"":lrrs[i]+"\t")+
	        				(gcs==null?"":gcs[i]+"\t")+
	        				(forwardGenotypes==null?"":FullSample.ALLELE_PAIRS[forwardGenotypes[i]]+"\t")+
	        				(abGenotypes==null?"":(abGenotypes[i]==-1?"--":FullSample.AB_PAIRS[abGenotypes[i]]))
	        				);
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
