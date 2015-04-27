package filesys;

import java.io.*;
import java.util.*;
import common.*;

public class GeneSet implements Serializable {
	public static final long serialVersionUID = 1L;
//	public static final String DIRECTORY = "/home/npankrat/NCBI/";
	public static final String DIRECTORY = "N:/statgen/NCBI/";
//	public static final String REFSEQ_FILE = "refGene.txt";
	public static final String REFSEQ_FILE = "refGene_hg19.txt";
	public static final String REFSEQ_DB = "RefSeq.genes";
	public static final String REFSEQ_TRACK = "RefSeq.gtrack";
	public static final String REFSEQ_SEGS = "RefSeq_2k.segs";
	public static final String REFSEQ_EXONS = "RefSeq_exons_2k.segs";
	public static final String KNOWN_FILE = "knownGene.txt";
	public static final String KNOWN_LOOKUP_FILE = "kgXref.txt";

	private GeneData[] set;
	private long fingerprint;
	
	public GeneSet(GeneData[] set) {
		this.set = set;
		this.fingerprint = (long)(Math.random()*Long.MAX_VALUE);
	}

	public GeneSet(Vector<GeneData> setVec) {
		this.set = new GeneData[setVec.size()];
		for (int i = 0; i<setVec.size(); i++) {
			this.set[i] = setVec.elementAt(i);
        }
		
		this.fingerprint = (long)(Math.random()*Long.MAX_VALUE);
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public GeneData[] getSet() {
		return set;
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static GeneSet load(String filename, boolean jar) {
		return (GeneSet)Files.readSerial(filename, jar, true);
	}
	
	public SegmentLists determineGeneSegments(int window) {
		Hashtable<String, Vector<Segment>> hash = new Hashtable<String,Vector<Segment>>();
		Vector<Segment> segs;
		Segment[][] lists;
		int[] chrs;

		for (int i = 0; i<set.length; i++) {
			if (hash.containsKey(set[i].getChr()+"")) {
				segs = hash.get(set[i].getChr()+"");
			} else {
				hash.put(set[i].getChr()+"", segs = new Vector<Segment>());
			}
			segs.add(new Segment(set[i].getChr(), set[i].getStart()-window, set[i].getStop()+window));
        }
		chrs = Array.toIntArray(HashVec.getKeys(hash));
		lists = new Segment[Array.max(chrs)+1][];
		for (int i = 0; i<chrs.length; i++) {
			segs = hash.get(chrs[i]+"");
        	Segment.mergeOverlapsAndSort(segs);
        	lists[chrs[i]] = Segment.toArray(segs);
        }
		
		return new SegmentLists(lists);
	}
	
	public SegmentLists determineExonSegments(int upDownWindow) {
		Hashtable<String, Vector<Segment>> hash = new Hashtable<String,Vector<Segment>>();
		Vector<Segment> segs;
		Segment[][] lists;
		int[] chrs;
		int[][] exonBoundaries;
		GeneData gene;

		for (int i = 0; i<set.length; i++) {
			gene = set[i];
			if (hash.containsKey(gene.getChr()+"")) {
				segs = hash.get(gene.getChr()+"");
			} else {
				hash.put(gene.getChr()+"", segs = new Vector<Segment>());
			}
			segs.add(new Segment(gene.getChr(), gene.getStart()-upDownWindow, gene.getStart()));
			segs.add(new Segment(gene.getChr(), gene.getStop(), gene.getStop()+upDownWindow));
    		exonBoundaries = gene.getExonBoundaries();
    		for (int j = 0; j<exonBoundaries.length; j++) {
    			segs.add(new Segment(exonBoundaries[j][0], exonBoundaries[j][1]));
            }
        }
		chrs = Array.toIntArray(HashVec.getKeys(hash));
		lists = new Segment[Array.max(chrs)+1][];
		for (int i = 0; i<chrs.length; i++) {
			segs = hash.get(chrs[i]+"");
        	Segment.mergeOverlapsAndSort(segs);
        	lists[chrs[i]] = Segment.toArray(segs);
        }
		
		return new SegmentLists(lists);
	}
	
	public static void parseRefSeqGenes() {
		BufferedReader reader;
        Hashtable<String,Vector<GeneData>> hash = new Hashtable<String,Vector<GeneData>>();
        Vector<GeneData> v, overlapping, finalList;
        GeneData gene;
        String[] geneNames;
        Vector<String> assessionNumbers;
        byte strand;
        Vector<Segment> exons;
        boolean newlyAdded;
        byte count;
        int[][] exonBoundaries;
        GeneSet geneSet;
        
        System.out.println("Parsing gene info...");
        try {
	        reader = new BufferedReader(new FileReader(DIRECTORY+REFSEQ_FILE));
	        while (reader.ready()) {
	        	gene = new GeneData(reader.readLine());
	        	if (gene.isFinalized()) {
		        	if (hash.containsKey(gene.getGeneName())) {
		        		v = hash.get(gene.getGeneName());
		        	} else {
		        		hash.put(gene.getGeneName(), v = new Vector<GeneData>());
		        	}
	        		v.add(gene);
	        	}
	        }	        	
	        reader.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+DIRECTORY+REFSEQ_FILE+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+DIRECTORY+REFSEQ_FILE+"\"");
	        System.exit(2);
        }
        
        System.out.println("Collapsing isoforms...");
        geneNames = HashVec.getKeys(hash);
        finalList = new Vector<GeneData>();
        for (int i = 0; i<geneNames.length; i++) {
    		v = hash.get(geneNames[i]);
    		count = 0;
    		while (v.size() > 0) {
        		overlapping = new Vector<GeneData>();
    			overlapping.add(v.remove(0));
				strand = overlapping.elementAt(0).getStrand();
    			newlyAdded = true;
    			while (newlyAdded) {
        			newlyAdded = false;
        			for (int j = 0; j<v.size(); j++) {
        				for (int k = 0; k<overlapping.size(); k++) {
            				if (v.elementAt(j).overlaps(overlapping.elementAt(k))) {
            	    			overlapping.add(v.remove(j));
            	    			j--;
            	    			k = overlapping.size();
            	    			newlyAdded = true;
            				}
                        }
                    }
    			}
    			count++; // genes can be in more than one "finished" location, and even be on different chromosomes (though half the time it's an X/Y pairing)
    			
    			assessionNumbers = new Vector<String>();
    			exons = new Vector<Segment>();
    			for (int j = 0; j<overlapping.size(); j++) {
    				gene = overlapping.elementAt(j);
    				assessionNumbers.add(gene.getNcbiAssessionNumbers()[0]);
    				if (gene.getStrand() != strand) {
    					strand = GeneData.BOTH_STRANDS; // genes can be listed as being on both strands
    				}
    				exonBoundaries = gene.getExonBoundaries();
    				
    				// different isoforms have different combinations of the same superset of exons
    				for (int k = 0; k<exonBoundaries.length; k++) {
    					Segment.addIfAbsent(new Segment(exonBoundaries[k][0], exonBoundaries[k][1]), exons);
                    }
                }
    			
    			// different isoforms can splice at different spots in the same exon
    			Segment.mergeOverlapsAndSort(exons);
    			
    			exonBoundaries = Segment.convertListToSortedBoundaries(exons);
    			if (geneNames[i].equals("MIR4444-2")) { // an example of a multi-chr gene
					for (int j2 = 0; j2 < overlapping.size(); j2++) {
						System.out.println(i+"\t"+overlapping.get(j2).getUCSClocation()+"\t"+count);
					}
				}
        		finalList.add(new GeneData(geneNames[i],
        				Array.toStringArray(assessionNumbers),
        				overlapping.elementAt(0).getChr(),
        				true,
        				strand, 
        				exonBoundaries[0][0],
        				exonBoundaries[exons.size()-1][1],
        				exonBoundaries,
        				(count==1&&v.size()==0?0:count)
        		));
    		}
        }
        
        geneSet = new GeneSet(finalList);
        geneSet.serialize(DIRECTORY+REFSEQ_DB);
        System.out.println("Parsing gene segments...");
        geneSet.determineGeneSegments(2000).serialize(DIRECTORY+REFSEQ_SEGS);
        System.out.println("Parsing exon segments...");
        geneSet.determineExonSegments(2000).serialize(DIRECTORY+REFSEQ_EXONS);
        System.out.println("Parsing gene track...");
        new GeneTrack(DIRECTORY+REFSEQ_DB).serialize(DIRECTORY+REFSEQ_TRACK);
	}

	public static void parseUCSCknownGenes() {
	}
	
	public static GeneSet loadRefSeqGenes() {
		return load(DIRECTORY+REFSEQ_DB, false);
	}

	public static void main(String[] args) {
	    int numArgs = args.length;
	    boolean refseq = true;
	    boolean known = false;

	    String usage = "\n"+"filesys.GeneSet requires 0-1 arguments\n"+
	    "   (1) parse RefSeq genes (i.e. -refseq (not the default))\n"+
	    "   (2) parse UCSC known genes (i.e. -known (not the default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("-refseq")) {
			    refseq = true;
			    numArgs--;
		    } else if (args[i].startsWith("-known")) {
			    known = true;
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
	    	if (refseq) {
	    		parseRefSeqGenes();
	    	}
	    	if (known) {
	    		parseUCSCknownGenes();
	    	}
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
