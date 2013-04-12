package cnv.analysis;

import java.io.*;
import java.util.Hashtable;

import filesys.*;
import common.*;
import cnv.filesys.*;
import cnv.park.MakeUCSCtrack;
import cnv.var.CNVariant;

public class FilterCalls {
	public static final int DEFAULT_MIN_SIZE_KB = 0;
	public static final int DEFAULT_MIN_NUM_SNPS = 1;
	public static final double DEFAULT_MIN_SCORE = 10.0;
	public static final boolean DEFAULT_FILTER_REGIONS_FLAG = false;
	public static final String DEFAULT_PROBLEMATIC_REGIONS = "data/problematicRegions.dat";
//	public static final String DEFAULT_COMMON_CNP_REFERENCE = "data/polymorphic_CNPs.txt";
//	public static final String DEFAULT_COMMON_CNP_REFERENCE = "data/common_CNPs.txt";
	public static final String DEFAULT_COMMON_CNP_REFERENCE = "data/all_CNPs.txt";
	public static final String[] DEFAULT_REGION_DIRECTORIES = {"C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\", "/home/npankrat/", "P:\\"};
	public static final int COMMON_IN = 1;
	public static final int COMMON_OUT = 2;
	public static final int COMMON_IGNORED = 3;
	public static final int DEFAULT_COMMON_IN_OUT_OR_IGNORED = COMMON_IGNORED;
	public static final boolean DEFAULT_BREAK_CENTROMERE = false;
	
	public static void filter(String dir, String in, String out, int delSize, int dupSize, int number, double score, boolean filterRegions, int commonInOutOrIgnore, String individualsToKeepFile, boolean breakupCentromeres) {
		String[] individualsToKeepList;
		if (individualsToKeepFile != null && !new File(individualsToKeepFile).exists()) {
			System.err.println("Error - could not find \""+individualsToKeepFile+"\" in directory; will not be able to filter by indiviudals");
			individualsToKeepFile = null;
		}
		individualsToKeepList = individualsToKeepFile==null?null:HashVec.loadFileToStringArray(individualsToKeepFile, false, false, new int[] {0,1}, true, false, "\t"); 

		filter(dir, in, out, delSize, dupSize, number, score, filterRegions, commonInOutOrIgnore, individualsToKeepList, breakupCentromeres);
	}
	
	public static void filter(String dir, String in, String out, int delSize, int dupSize, int number, double score, boolean filterRegions, int commonInOutOrIgnore, String[] individualsToKeepList, boolean breakupCentromeres) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		CNVariant cnv;
		Segment[] problemRegions, centromereMidpoints, commonReference;
		Hashtable<String, String> indHash;
		int count, countCentromeric;

		problemRegions = filterRegions?Segment.loadUCSCregions(Files.firstDirectoryThatExists(DEFAULT_REGION_DIRECTORIES, true, true)+DEFAULT_PROBLEMATIC_REGIONS, false):new Segment[0];
		centromereMidpoints = loadCentromereMidpoints(Positions.CENTROMERE_MIDPOINT_SEGS);
		commonReference = commonInOutOrIgnore!=COMMON_IGNORED?Segment.loadUCSCregions(Files.firstDirectoryThatExists(DEFAULT_REGION_DIRECTORIES, true, true)+DEFAULT_COMMON_CNP_REFERENCE, false):new Segment[0];
		
		indHash = individualsToKeepList==null?null:HashVec.loadToHashNull(individualsToKeepList);

		try {
			reader = new BufferedReader(new FileReader(dir+in));
			writer = new PrintWriter(new FileWriter(dir+out));
			System.out.println("Writing to '"+dir+out+"'");
			writer.println(reader.readLine());
			count = 0;
			countCentromeric = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				cnv = new CNVariant(line);
				if (((cnv.getCN() < 2 && cnv.getSize()>=delSize*1000) || (cnv.getCN() > 2 && cnv.getSize()>=dupSize*1000))
						&&cnv.getNumMarkers()>=number&&cnv.getScore()>score&&!inOneOfTheseRegions(cnv, problemRegions)) {
					if (spansCentromereMidpoint(cnv, centromereMidpoints)) {
						if (breakupCentromeres) {
							line[3] = cnv.getStart()+"";
							line[4] = Positions.CENTROMERE_BOUNDARIES_FROM_SNPS[cnv.getChr()][0]+"";
							writer.println(Array.toStr(line));
							line[3] = cnv.getStart()+"";
							line[4] = Positions.CENTROMERE_BOUNDARIES_FROM_SNPS[cnv.getChr()][0]+"";
							writer.println(Array.toStr(line));
						}
						countCentromeric++;
//						System.err.println("Warning - a CNV for "+cnv.getFamilyID()+","+cnv.getIndividualID()+" spans a centromere ("+cnv.getUCSClocation()+") with "+cnv.getNumMarkers()+" markers");
					} else {
						if ( (commonInOutOrIgnore==COMMON_IGNORED
								||(commonInOutOrIgnore==COMMON_IN&&inOneOfTheseRegions(cnv, commonReference))
								||(commonInOutOrIgnore==COMMON_OUT&&!inOneOfTheseRegions(cnv, commonReference)) )
							&& (indHash == null
								|| indHash.containsKey(line[0]+"\t"+line[1])) ) {
							writer.println(Array.toStr(line));
						}
					}
					if (cnv.getSize()>10000000||cnv.getNumMarkers()>500) {
//						System.err.println("Warning - "+cnv.getFamilyID()+","+cnv.getIndividualID()+" has a gigantic CNV spanning "+ext.prettyUpDistance(cnv.getSize(), 0)+" and "+cnv.getNumMarkers()+" markers ("+cnv.getUCSClocation()+")");
						count++;
					}
				}
			}
			System.err.println("Identified "+count+" gigantic CNVs");
			System.err.println("Identified "+countCentromeric+" CNVs that spanned centromeres; these were "+(breakupCentromeres?"broken up into two CNVs, one on each side of the centromere":"ignored"));
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+in+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+in+"\"");
			System.exit(2);
		}
	}
	
	public static void filterOnSegments(String dir, String filein, String fileout, String segmentFile, boolean excludeInsteadOfInclude) {
		BufferedReader reader;
        PrintWriter writer;
        Segment[][] genesByChr;
        CNVariant cnv;
        SegmentList segList;

        if (segmentFile.endsWith(".segs")) {
            segList = SegmentList.load(segmentFile, false);
        } else if (new File(segmentFile+".segs").exists()) {
        	segList = SegmentList.load(segmentFile+".segs", false);
        } else {
        	segList = SegmentList.parseSegmentList(segmentFile, false);
        	segList.serialize(segmentFile+".segs");
        }

        genesByChr = segList.getLists();

        try {
	        reader = new BufferedReader(new FileReader(dir+filein));
	        writer = new PrintWriter(new FileWriter(dir+fileout));
	        writer.println(reader.readLine());
	        while (reader.ready()) {
	        	cnv = new CNVariant(reader.readLine().trim().split("[\\s]+"));
	        	if (genesByChr[cnv.getChr()] != null && Segment.overlapsAny(new Segment((byte)cnv.getChr(), cnv.getStart(), cnv.getStop()), genesByChr[cnv.getChr()])) {
	        		writer.println(cnv.toPlinkFormat());
	        	}
	        }	        
	        reader.close();
            writer.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+dir+filein+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+dir+filein+"\"");
	        System.exit(2);
        }	
	}	

	public static void filterBasedOnNumberOfCNVsAtLocus(Project proj, String filein, String fileout, int totalRequired, int delRequired, int dupRequired, int delLimitedTo, int dupLimitedTo) {
		PrintWriter writer;
		MarkerSet markerSet;
		int[][] positions;
		int[][][] counts;
		int firstSNP, lastSNP, indel;
		CNVariant[] cnvs;
		int index;
		boolean[][] acceptableSNPs;
		boolean accepted;
		int dels, dups;
		
		markerSet = proj.getMarkerSet();
		positions = markerSet.getPositionsByChr();
		counts = new int[positions.length][][];
		acceptableSNPs = new boolean[positions.length][];
		for (int i = 0; i<positions.length; i++) {
			counts[i] = new int[positions[i].length][2];
			acceptableSNPs[i] = new boolean[positions[i].length];
        }
		
		cnvs = CNVariant.loadPlinkFile(filein, false);

		for (int i = 0; i<cnvs.length; i++) {
			firstSNP = Array.binarySearch(positions[cnvs[i].getChr()], cnvs[i].getStart(), true);
			lastSNP = Array.binarySearch(positions[cnvs[i].getChr()], cnvs[i].getStop(), true);
			indel = cnvs[i].getCN()<2?0:1;
			for (int j = firstSNP; j<=lastSNP; j++) {
				counts[cnvs[i].getChr()][j][indel]++;
            }
        }
		
		for (int i = 0; i<positions.length; i++) {
			for (int j = 0; j<positions[i].length; j++) {
				dels = counts[i][j][0];
				dups = counts[i][j][1];
				acceptableSNPs[i][j] = dels + dups >= totalRequired && dels >= delRequired && dups >= dupRequired && dels <= delLimitedTo && dups <= dupLimitedTo;  
            }
        }
		
		try {
	        writer = new PrintWriter(new FileWriter(fileout));
			for (int i = 0; i<cnvs.length; i++) {
				firstSNP = Array.binarySearch(positions[cnvs[i].getChr()], cnvs[i].getStart(), true);
				lastSNP = Array.binarySearch(positions[cnvs[i].getChr()], cnvs[i].getStop(), true);
				indel = cnvs[i].getCN()<2?0:1;

				index = firstSNP;
				accepted = false;
				while (!accepted && index <= lastSNP) {
					if (acceptableSNPs[cnvs[i].getChr()][index]) {
						accepted = true;
					}
					index++;
				}
				
				if (accepted) {
	        		writer.println(cnvs[i].toPlinkFormat());
				}
	        }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+fileout);
	        e.printStackTrace();
        }
	}
	
	public static boolean inOneOfTheseRegions(CNVariant cnv, Segment[] regions) {
		for (int i = 0; i<regions.length; i++) {
			if (cnv.significantOverlap(regions[i])) {
				return true;
			}
		}
		return false;
	}

	public static CNVariant[] loadCentromereMidpoints(String[] midpointPositions) {
		CNVariant[] midpoints = new CNVariant[midpointPositions.length];

		for (int i = 0; i<midpoints.length; i++) {
			midpoints[i] = new CNVariant(midpointPositions[i]);
		}

		return midpoints;
	}

	public static boolean spansCentromereMidpoint(CNVariant cnv, Segment[] midpoints) {
		for (int i = 0; i<midpoints.length; i++) {
			if (cnv.overlaps(midpoints[i])) {
				return true;
			}
		}
		return false;
	}
	
	public static String getFilename(String root, int delSize, int dupSize, int number, double score, int commonInOutOrIgnore) {
		return root+"_"+(delSize == dupSize?delSize+"kb":delSize+","+dupSize+"kb")+"_"+number+"SNP_"+score+"_"+(commonInOutOrIgnore==COMMON_IN?"isCNP":(commonInOutOrIgnore==COMMON_OUT?"notCNP":"CNPstatusIgnored"))+".cnv";

	}
	
	public static void union(String firstCNVfile, String secondCNVfile, String outputfile) {
		PrintWriter writer;
		CNVariant[] list1, list2;
		int count;
		boolean unique;
		
		list1 = CNVariant.loadPlinkFile(firstCNVfile, false);
		list2 = CNVariant.loadPlinkFile(secondCNVfile, false);
		
		try {
	        writer = new PrintWriter(new FileWriter(outputfile));
			for (int i = 0; i<list1.length; i++) {
				writer.println(list1[i].toPlinkFormat());
	        }
			for (int i = 0; i<list2.length; i++) {
				count = 0;
				unique = true;
				while (unique && count < list1.length) {
					if (list1[count].equalsIncludingIndividual(list2[i])) {
						unique = false;
					}
					count++;
				}
				if (unique) {
					writer.println(list2[i].toPlinkFormat());
				}
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+outputfile);
	        e.printStackTrace();
        }
	}
	
	public static void stdFilters(String dir, String filename, boolean makeUCSCtracks, String pedfile) {
		String root;

		root = ext.rootOf(filename);
		FilterCalls.filter(dir, filename, root+"_allAbove10.0_unfiltered.cnv", 1, 1, 1, 10, false, COMMON_IGNORED, pedfile, true);
		FilterCalls.filter(dir, filename, root+"_allAbove10.0_filtered.cnv", 1, 1, 1, 10, true, COMMON_IGNORED, pedfile, true);
		FilterCalls.filter(dir, filename, root+"_ConservativeCalls.cnv", 100, 100, 20, 10, true, COMMON_IGNORED, pedfile, true);
		if (makeUCSCtracks) {
			MakeUCSCtrack.makeTrack(dir+root+"_allAbove10.0_unfiltered.cnv.cnv");
			MakeUCSCtrack.makeTrack(dir+root+"_ConservativeCalls.cnv");
		}
		FilterCalls.filterOnSegments(dir, root+"_allAbove10.0_filtered.cnv", root+"_allAbove10.0_filtered_inGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, false);
		FilterCalls.filterOnSegments(dir, root+"_allAbove10.0_filtered.cnv", root+"_allAbove10.0_filtered_inExons.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_EXONS, false);

//		FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new Project(Project.DEFAULT_PROJECT, false), dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_3dels_LT2dups.cnv", 0, 3, 0, Integer.MAX_VALUE, 1);
//		FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new Project(Project.DEFAULT_PROJECT, false), dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_3dups_LT2dels.cnv", 0, 0, 3, 1, Integer.MAX_VALUE);
//		FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new Project(Project.DEFAULT_PROJECT, false), dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_2dels_2dups.cnv", 0, 2, 2, Integer.MAX_VALUE, Integer.MAX_VALUE);
//		FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new Project(Project.DEFAULT_PROJECT, false), dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_3anythings.cnv", 3, 0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
//		FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new Project(Project.DEFAULT_PROJECT, false), dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_5anythings.cnv", 5, 0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
//		FilterCalls.filterOnSegments(dir+root+"_0kb_5SNP_10.0_3anythings.cnv", dir+root+"_CommonInGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, true);
		
//		FilterCalls.union(dir+root+"_0kb_5SNP_10.0_3anythings.cnv", dir+root+"_100kb_20SNP_10.0_CNPstatusIgnored.cnv", dir+"unionOfConservativeAndCommon.cnv");
//		FilterCalls.filterOnSegments(dir+"unionOfConservativeAndCommon.cnv", dir+"unionOfConservativeAndCommonInGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, true);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		int delSize = DEFAULT_MIN_SIZE_KB;
		int dupSize = DEFAULT_MIN_SIZE_KB;
		int number = DEFAULT_MIN_NUM_SNPS;
		double score = DEFAULT_MIN_SCORE;
		boolean filterRegions = DEFAULT_FILTER_REGIONS_FLAG;
		int inOutIgnore = COMMON_IGNORED;
//		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\penncnv\\again_noGenderProblems\\";
//		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\quantisnp\\noGenderProblems\\";
//		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\osteo\\";
//		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\Singleton\\";
//		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\penncnv\\";
		String dir = "";
		String in = "conf.cnv";
		String out = getFilename("conf", delSize, dupSize, number, score, inOutIgnore);
		boolean standards = false;
		boolean tracks = false;
//		String pedfile = "plink.fam";
		String pedfile = null;
		String segs = "";
		boolean excludeSegs = false;
//		boolean sigOverlap = false;
		boolean breakCent = DEFAULT_BREAK_CENTROMERE;

		String usage = 
		"vis.cnv.FilterCalls requires 0-1 arguments\n"+
		"   (1) directory (i.e. dir="+dir+" (default))\n"+
		"   (2) file in (i.e. in="+in+" (default))\n"+
		"   (3) file out (i.e. out="+out+" (default))\n"+
		"   (4) minimum size of a deletion (in kb) (i.e. delSize="+delSize+" (default))\n"+
		"   (5) minimum size of a duplication (in kb) (i.e. dupSize="+dupSize+" (default))\n"+
		"   (6) minimum number of SNPs (i.e. number="+number+" (default))\n"+
		"   (7) minimum score (i.e. score="+score+" (default))\n"+
		"   (8) filter out cnvs in known problematicRegions (i.e. filterRegions="+filterRegions+" (default))\n"+
		"   (9) pedfile to use as a filter (i.e. ped="+pedfile+" (default))\n"+
		"   (10) if CNV spans centromere, break into two spanning actual markers (i.e. breakCentromere="+breakCent+" (default))\n"+
		"\n OR\n\n"+
		"   (1) perform all standard filters (i.e. -std (not the default))\n"+
		"   (2) make UCSC tracks as well (i.e. -tracks (not the default))\n"+
		"\n OR\n\n"+
		"   (1) keep only CNVs overlapping these segments (i.e. segs=gene_region.dat (not the default))\n"+
		"   (2) exlcude instead of include (i.e. -excludeSegs (not the default))\n"+
//		"   (3) requre a significant overlap to filter (i.e. -sigOverlap (not the default))\n"+
		"";

		System.out.println();
		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("in=")) {
				in = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				out = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("delSize=")) {
				delSize = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("dupSize=")) {
				dupSize = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("number=")) {
				number = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("score=")) {
				score = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("cnps=")) {
				inOutIgnore = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("filterRegions=")) {
				filterRegions = args[i].split("=")[1].toLowerCase().equals("true");
				numArgs--;
			} else if (args[i].startsWith("-std")) {
				standards = true;
				numArgs--;
			} else if (args[i].startsWith("-tracks")) {
				tracks = true;
				numArgs--;
			} else if (args[i].startsWith("segs=")) {
				segs = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("breakCentromere=")) {
				breakCent = ext.parseBooleanArg(args[i]);
				numArgs--;				
			} else if (args[i].startsWith("-excludeSegs")) {
				excludeSegs = true;
				numArgs--;				
			} else {
				System.err.println("Error - don't know what to do with argument: "+args[i]);
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		FilterCalls.filterOnSegments("D:/data/GEDI/global/homoDelsOverlappingGenesOnly/", "conf.cnv", "conf_overlappingGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, false);
		System.exit(1);
		
//		breakCent = true;
//		out = "noCentromeric.cnv";
		
		try {
//			FilterCalls.filterOnSegments(dir+"conf_100kb_20SNP_10.0_CNPstatusIgnored.cnv", dir+"ConservativeGeneCentric.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, true);
//			MakeUCSCtrack.makeTrack(dir, "ConservativeGeneCentric.cnv");
			if (standards) {
				stdFilters(dir, in, tracks, pedfile);
			} else if (!segs.equals("")) {
				filterOnSegments(dir, in, out, segs, excludeSegs);
			} else {
				filter(dir, in, out, delSize, dupSize, number, score, filterRegions, DEFAULT_COMMON_IN_OUT_OR_IGNORED, pedfile, breakCent);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
