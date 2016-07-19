package one;

import java.io.*;
import common.*;
import filesys.*;

public class ThousandGenomeVariation {
	public static final String[] HEADER = {"chr", "loc", "ref", "alleles", "snp.Q", "av.max.map.Q", "depth.cov", "NA12891", "NA12891.Q", "NA12892", "NA12892.Q", "NA12878", "NA12878.Q", "hwe", "maf", "tdt", "display"}; 

	public static void parseVariation(String dir, String filename, String[] regions, String variation) {
		BufferedReader reader;
        PrintWriter writer, writer2;
        String[] line;
        String temp, trav;
        int count;
        Segment[] segs;
        byte chr;
        int pos;
        String combo;

        segs = Segment.getSegments(regions);
        try {
	        reader = new BufferedReader(new FileReader(dir+filename));
	        ext.checkHeader(reader.readLine().trim().split("[\\s]+"), HEADER, false);
	        writer = new PrintWriter(new FileWriter(dir+variation));
	        writer.println("track name=\"DeepTrioVariation\" description=\"1000genomes variation\" visibility=2 itemRgb=\"On\"");
	        writer2 = new PrintWriter(new FileWriter(dir+"region.xln"));
	        while (reader.ready()) {
	        	temp = reader.readLine();
	        	line = temp.trim().split("[\\s]+");
	        	try {
	        		chr = Positions.chromosomeNumber(line[0]);
		        	pos = Integer.parseInt(line[1]);
	        	} catch (Exception e) {
	        		System.err.println("Error - invalid location (chr: "+line[0]+"; position: "+line[1]+")");
	        		chr = -1;
	        		pos = -1;
	        	}
	        	if (Segment.overlapsAny(new Segment(chr, pos, pos), segs)) {
	        		combo = "";
	        		for (int i = 0; i<3; i++) {
	        			count = 0;
	        			trav = line[7+i*2];
	        			for (int j = 0; j<2; j++) {
		        			if (trav.toLowerCase().charAt(j*2) == line[3].toLowerCase().charAt(2)) {
		        				count++;
		        			}
                        }
	        			combo += count;
                    }
	        		writer.println("chr"+chr+"\t"+pos+"\t"+(pos+1)+"\t"+combo+"\t"+"Chromosome"+"^"+chr+"^"+pos+"^"+line[3].toUpperCase().charAt(0)+"^"+line[3].toUpperCase().charAt(2)+"^1"+"\t"+chr+","+pos+",1,"+line[3].toUpperCase().charAt(0)+"/"+line[3].toUpperCase().charAt(2));
	        		writer.flush();
	        		writer2.println(temp);
	        		writer2.flush();
	        	}
	        }
	        reader.close();
            writer.close();
            writer2.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+dir+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+dir+filename+"\"");
	        System.exit(2);
        }
	}
	
	public static void countVariationInEach(String dir, String[] bedfiles) {
		BufferedReader reader;
        int[] counts;
        String trav;
        
        for (int i = 0; i<bedfiles.length; i++) {
            try {
    	        reader = new BufferedReader(new FileReader(dir+bedfiles[i]));
    	        reader.readLine();
    	        counts = new int[3];
    	        while (reader.ready()) {
    	        	trav = reader.readLine().trim().split("[\\s]+")[3];
    	        	for (int j = 0; j<3; j++) {
    	        		if (Integer.parseInt(trav.charAt(j)+"") == 1) {
    	        			counts[j]++;
    	        		}
                    }
    	        }
    	        System.out.println(Array.toStr(counts));
    	        reader.close();
            } catch (FileNotFoundException fnfe) {
    	        System.err.println("Error: file \""+dir+bedfiles[i]+"\" not found in current directory");
    	        System.exit(1);
            } catch (IOException ioe) {
    	        System.err.println("Error reading file \""+dir+bedfiles[i]+"\"");
    	        System.exit(2);
            }

        }
        
	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    String filename = "CEU.trio.dec.with.x.with.rs.calls";
	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\SequencingProjectWithCIDR\\TestInDeepTrio\\";
	    String variation = "2qRegion.out";
	    String[] regions = new String[] {"chr2:221735000-240816800"};
//	    String[] regions = new String[] {"chr1:1-742430"};
	    String[] bedfiles = {"2qRegionPlus1.bed", "Variation_Overlapping_Baits.bed", "Variation_Overlapping_Captured_Exons.bed", "Variation_Overlapping_Captured_Coding_Exons.bed"};

	    String usage = "\n"+"one.ThousandGenomeVariation requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"";

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
	    	if (!new File(dir+variation).exists()) {
	    		parseVariation(dir, filename, regions, variation);
	    	}
	    	countVariationInEach(dir, bedfiles);
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
