package org.genvisis.bioinformatics;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;
import org.genvisis.filesys.*;

public class GeneDensityInRegion {
	public static void computeDensity(String region, int window) {
        PrintWriter writer;
        GeneSet geneSet;
        int[] loc;
        Vector<GeneData> inRegion;
        GeneData[] genes;
        GeneData regionAsGene;
        Vector<Segment> segments;
        int sum;
        int[][] exonBoundaries;
        
        geneSet = GeneSet.load(Aliases.getPathToFileInReferenceDirectory(GeneSet.REFSEQ_DB, true, new Logger()), false);
        genes = geneSet.getSet();
        loc = Positions.parseUCSClocation(region);
        regionAsGene = new GeneData("", new String[0], (byte)loc[0], true, (byte)0, loc[1], loc[2], new int[][] {}, (byte)0, false);
        inRegion = new Vector<GeneData>();
        for (int i = 0; i<genes.length; i++) {
        	if (genes[i].overlaps(regionAsGene)) {
        		inRegion.add(genes[i]);
        	}
        }
        genes = new GeneSet(inRegion).getSet();
        try {
        	writer = new PrintWriter(new FileWriter(ext.replaceAllWith(region, ":", "_")+".xln"));
        	writer.println("Gene\tAssession #'s\tChr\tStart\tStop\tNumExons");
        	for (int i = 0; i<genes.length; i++) {
        		writer.println(genes[i].getGeneName()+"\t"+Array.toStr(genes[i].getNcbiAssessionNumbers(), "|")+"\t"+genes[i].getChr()+"\t"+genes[i].getStart()+"\t"+genes[i].getStop()+"\t"+genes[i].getExonBoundaries().length);
            }
            writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing results");
	        e.printStackTrace();
        }
        
        try {
        	writer = new PrintWriter(new FileWriter(ext.replaceAllWith(region, ":", "_")+".out"));
        	writer.println("Total region: "+ext.prettyUpDistance(regionAsGene.getSize(), 1));
        	writer.println("Number of non-redundant RefSeq genes: "+genes.length);

        	segments = new Vector<Segment>();
        	for (int i = 0; i<genes.length; i++) {
        		segments.add(new Segment(genes[i].getStart(), genes[i].getStop()));
            }
        	Segment.mergeOverlapsAndSort(segments);
        	sum = 0;
        	for (int i = 0; i<segments.size(); i++) {
        		sum += segments.elementAt(i).getSize();
            }
        	writer.println("Genes span: "+ext.prettyUpDistance(sum, 1)+" ("+(int)((double)sum/(double)regionAsGene.getSize()*100)+"%)");

        	segments = new Vector<Segment>();
        	for (int i = 0; i<genes.length; i++) {
        		segments.add(new Segment(genes[i].getStart()-window, genes[i].getStop()+window));
            }
        	Segment.mergeOverlapsAndSort(segments);
        	sum = 0;
        	for (int i = 0; i<segments.size(); i++) {
        		sum += segments.elementAt(i).getSize();
            }
        	writer.println("Including "+ext.prettyUpDistance(window, 0)+" up- and down-stream, the genes span: "+ext.prettyUpDistance(sum, 1)+" ("+(int)((double)sum/(double)regionAsGene.getSize()*100)+"%)");

        	segments = new Vector<Segment>();
        	for (int i = 0; i<genes.length; i++) {
        		exonBoundaries = genes[i].getExonBoundaries();
        		for (int j = 0; j<exonBoundaries.length; j++) {
            		segments.add(new Segment(exonBoundaries[j][0], exonBoundaries[j][1]));
                }
            }
        	Segment.mergeOverlapsAndSort(segments);
        	sum = 0;
        	for (int i = 0; i<segments.size(); i++) {
        		sum += segments.elementAt(i).getSize();
            }
        	writer.println("Number of exons: "+segments.size());
        	writer.println("Exons span: "+ext.prettyUpDistance(sum, 1)+" ("+(int)((double)sum/(double)regionAsGene.getSize()*100)+"%)");

        	for (int i = 0; i<genes.length; i++) {
        		segments.add(new Segment(genes[i].getStart()-window, genes[i].getStart()));
        		segments.add(new Segment(genes[i].getStop(), genes[i].getStop()+window));
            }
        	Segment.mergeOverlapsAndSort(segments);
        	sum = 0;
        	for (int i = 0; i<segments.size(); i++) {
        		sum += segments.elementAt(i).getSize();
            }
        	writer.println("Including "+ext.prettyUpDistance(window, 0)+" both 5' and 3', the total non-overlapping distance is: "+ext.prettyUpDistance(sum, 1)+" ("+(int)((double)sum/(double)regionAsGene.getSize()*100)+"%)");

        	
            writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing results");
	        e.printStackTrace();
        }
        
        Files.more(ext.replaceAllWith(region, ":", "_")+".out");
		
	}

	public static void main(String[] args) {
	    int numArgs = args.length;
	    String region = "chr2:221735000-240816800"; // Original 19.1 Mb
//	    String region = "chr2:227241128-235934686"; // 8.7 Mb
	    int window = 2000;
	    
	    String usage = "\n"+
	    "bioinformatics.GeneDensityInRegion requires 0-1 arguments\n"+
	    "   (1) region (i.e. region="+region+" (default))\n"+
	    "   (2) window in bases (i.e. win="+window+" (default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("region=")) {
			    region = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("window=")) {
			    window = Integer.parseInt(args[i].split("=")[1]);
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
		    computeDensity(region, window);
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
