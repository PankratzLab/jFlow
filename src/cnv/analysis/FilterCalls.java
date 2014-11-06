package cnv.analysis;

import java.io.*;
import java.util.Date;
import java.util.HashSet;
import java.util.Vector;

import filesys.*;
import common.*;
import cnv.filesys.*;
import cnv.manage.UCSCtrack;
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
	
	public static void filter(String dir, String in, String out, int delSize, int dupSize, int number, double score, String filenameOfProblematicRegions, int commonInOutOrIgnore, String individualsToKeepFile, boolean breakupCentromeres, String markerSetFilenameToBreakUpCentromeres, boolean makeUCSCtrack, int build, Logger log) {
		String[] individualsToKeepList;
		if (individualsToKeepFile != null && !new File(individualsToKeepFile).exists()) {
			System.err.println("Error - could not find \""+individualsToKeepFile+"\" in directory; will not be able to filter by indiviudals");
			individualsToKeepFile = null;
		}
		individualsToKeepList = individualsToKeepFile==null?null:HashVec.loadFileToStringArray(individualsToKeepFile, false, false, new int[] {0,1}, true, false, "\t"); 

		filter(dir, in, out, delSize, dupSize, number, score, filenameOfProblematicRegions, commonInOutOrIgnore, individualsToKeepList, breakupCentromeres, markerSetFilenameToBreakUpCentromeres, makeUCSCtrack, build, log);
	}
	
	public static void filter(String dir, String in, String out, int delSize, int dupSize, int number, double score, String filenameOfProblematicRegions, int commonInOutOrIgnore, String[] individualsToKeepList, boolean breakupCentromeres, String markerSetFilenameToBreakUpCentromeres, boolean makeUCSCtrack, int build, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		CNVariant cnv;
		Segment[] problemRegions, centromereMidpoints, commonReference;
		HashSet<String> indHash;
		int countGiant, countCentromeric, countGiantCentromeric;
		int[][] centromereBoundaries;

		problemRegions = filenameOfProblematicRegions==null?new Segment[0]:Segment.loadUCSCregions(filenameOfProblematicRegions, 0, false, log);
		centromereBoundaries = Positions.determineCentromereBoundariesFromMarkerSet(markerSetFilenameToBreakUpCentromeres, build, log);
		centromereMidpoints = Positions.computeCentromereMidpoints(centromereBoundaries);
		commonReference = commonInOutOrIgnore!=COMMON_IGNORED?Segment.loadUCSCregions(Files.firstDirectoryThatExists(DEFAULT_REGION_DIRECTORIES, true, true)+DEFAULT_COMMON_CNP_REFERENCE, false):new Segment[0];
		indHash = individualsToKeepList==null?null:HashVec.loadToHashSet(individualsToKeepList);

		try {
			reader = new BufferedReader(new FileReader(dir+in));
			writer = new PrintWriter(new FileWriter(dir+out));
			System.out.println("Writing to '"+dir+out+"'");
			writer.println(reader.readLine());
			countGiant = 0;
			countCentromeric = 0;
			countGiantCentromeric = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				cnv = new CNVariant(line);
				if (((cnv.getCN() < 2 && cnv.getSize()>=delSize*1000) || (cnv.getCN() > 2 && cnv.getSize()>=dupSize*1000))
						&&cnv.getNumMarkers()>=number&&cnv.getScore()>score&&!inOneOfTheseRegions(cnv, problemRegions)) {

					if ( (commonInOutOrIgnore==COMMON_IGNORED
							||(commonInOutOrIgnore==COMMON_IN&&inOneOfTheseRegions(cnv, commonReference))
							||(commonInOutOrIgnore==COMMON_OUT&&!inOneOfTheseRegions(cnv, commonReference)) )
						&& (indHash == null
							|| indHash.contains(line[0]+"\t"+line[1])) )
					{
						if (cnv.overlaps(centromereMidpoints[cnv.getChr()])) {
							if (breakupCentromeres) {
//								System.out.println("Splitting "+cnv.getUCSClocation()+" due to overlap with "+centromereMidpoints[cnv.getChr()].getUCSClocation()+" using boundaries "+Array.toStr(centromereBoundaries[cnv.getChr()], ", "));
								line[3] = cnv.getStart()+"";
								line[4] = centromereBoundaries[cnv.getChr()][0]+"";
								writer.println(Array.toStr(line));
								line[3] = centromereBoundaries[cnv.getChr()][1]+"";
								line[4] = cnv.getStop()+"";
								writer.println(Array.toStr(line));
//								return;
							}
							countCentromeric++;
							
							if (cnv.getSize()>10000000) {
								countGiantCentromeric++;
							}
							
//							System.err.println("Warning - a CNV for "+cnv.getFamilyID()+","+cnv.getIndividualID()+" spans a centromere ("+cnv.getUCSClocation()+") with "+cnv.getNumMarkers()+" markers");
						} else {
							writer.println(Array.toStr(line));
						}
					}						
					if (cnv.getSize()>10000000||cnv.getNumMarkers()>500) {
//						System.err.println("Warning - "+cnv.getFamilyID()+","+cnv.getIndividualID()+" has a gigantic CNV spanning "+ext.prettyUpDistance(cnv.getSize(), 0)+" and "+cnv.getNumMarkers()+" markers ("+cnv.getUCSClocation()+")");
						countGiant++;
					}
				}
			}
			System.err.println("Identified "+countCentromeric+" CNVs that spanned centromeres; these were "+(breakupCentromeres?"broken up into two CNVs, one on each side of the centromere":"retained as is"));
			System.err.println("Identified "+countGiant+" gigantic CNVs ( 10+ Mb or 500+ probes ), of which "+countGiantCentromeric+" spanned a centromere");
			reader.close();
			writer.close();
			if (makeUCSCtrack) {
				UCSCtrack.makeTrack(dir+out, dir+ext.rootOf(out)+".bed.gz", log);
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+in+"\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+in+"\"");
			return;
		}
	}
	
	public static void filterOnSegments(String dir, String filein, String fileout, String segmentFile, boolean excludeInsteadOfInclude) {
		BufferedReader reader;
        PrintWriter writer;
        Segment[][] genesByChr;
        CNVariant cnv;
        SegmentLists segList;

        if (segmentFile.endsWith(".segs")) {
            segList = SegmentLists.load(segmentFile, false);
        } else if (new File(segmentFile+".segs").exists()) {
        	segList = SegmentLists.load(segmentFile+".segs", false);
        } else {
        	segList = SegmentLists.parseUCSCSegmentList(segmentFile, false);
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
			return;
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+dir+filein+"\"");
			return;
        }	
	}	

	public static void filterBasedOnNumberOfCNVsAtLocus(Project proj, String filein, String fileout, int totalRequired, int delRequired, int dupRequired, int totalLimitedTo, int delLimitedTo, int dupLimitedTo, double proportionOfProbesThatNeedToPassForFinalInclusion) {
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
		int countAcceptable;
		long time;
		
		time = new Date().getTime();
		
		markerSet = proj.getMarkerSet();
		positions = markerSet.getPositionsByChr();
		counts = new int[positions.length][][];
		acceptableSNPs = new boolean[positions.length][];
		for (int i = 0; i<positions.length; i++) {
			counts[i] = new int[positions[i].length][2];
			acceptableSNPs[i] = new boolean[positions[i].length];
        }
		
		System.out.println(ext.getTime()+"\tLoading plink file...");
		cnvs = CNVariant.loadPlinkFile(filein, false);

		System.out.println(ext.getTime()+"\tDetermining acceptability...");
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
				acceptableSNPs[i][j] = dels + dups >= totalRequired && dels >= delRequired && dups >= dupRequired && dels + dups <= totalLimitedTo && dels <= delLimitedTo && dups <= dupLimitedTo;  
            }
        }
		
		System.out.println(ext.getTime()+"\tFiltering CNVs...");
		try {
	        writer = new PrintWriter(new FileWriter(fileout));
	        writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			for (int i = 0; i<cnvs.length; i++) {
				firstSNP = Array.binarySearch(positions[cnvs[i].getChr()], cnvs[i].getStart(), true);
				lastSNP = Array.binarySearch(positions[cnvs[i].getChr()], cnvs[i].getStop(), true);
				indel = cnvs[i].getCN()<2?0:1;

				if (proportionOfProbesThatNeedToPassForFinalInclusion < 1.0) {
					countAcceptable = 0;
					for (int j = firstSNP; j <= lastSNP; j++) {
						if (acceptableSNPs[cnvs[i].getChr()][j]) {
							countAcceptable++;
						}
					}
					accepted = (double)countAcceptable / (double)(lastSNP - firstSNP + 1) > proportionOfProbesThatNeedToPassForFinalInclusion;
				} else {
					index = firstSNP;
					accepted = false;
					while (!accepted && index <= lastSNP) {
						if (acceptableSNPs[cnvs[i].getChr()][index]) {
							accepted = true;
						}
						index++;
					}
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
		
		System.out.println("Finished in " + ext.getTimeElapsed(time));
	}
	
	public static boolean inOneOfTheseRegions(CNVariant cnv, Segment[] regions) {
		for (int i = 0; i<regions.length; i++) {
			if (cnv.significantOverlap(regions[i])) {
				return true;
			}
		}
		return false;
	}

//	public static boolean spansCentromereMidpoint(CNVariant cnv, Segment[] midpoints) {
//		for (int i = 0; i<midpoints.length; i++) {
//			if (cnv.overlaps(midpoints[i])) {
//				return true;
//			}
//		}
//		return false;
//	}
	
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
	
	public static void stdFilters(String dir, String filename, boolean makeUCSCtracks, String pedfile, int build) {
		String root;
		Logger log;

		log = new Logger();
		root = ext.rootOf(filename);
		FilterCalls.filter(dir, filename, root+"_allAbove10.0_unfiltered.cnv", 1, 1, 1, 10, null, COMMON_IGNORED, pedfile, true, null, makeUCSCtracks,  build, log);
		FilterCalls.filter(dir, filename, root+"_allAbove10.0_filtered.cnv", 1, 1, 1, 10, DEFAULT_PROBLEMATIC_REGIONS, COMMON_IGNORED, pedfile, true, null, makeUCSCtracks, build, log);
		FilterCalls.filter(dir, filename, root+"_ConservativeCalls.cnv", 100, 100, 20, 10, DEFAULT_PROBLEMATIC_REGIONS, COMMON_IGNORED, pedfile, true, null, makeUCSCtracks, build, log);

		FilterCalls.filterOnSegments(dir, root+"_allAbove10.0_filtered.cnv", root+"_allAbove10.0_filtered_inGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, false);
		FilterCalls.filterOnSegments(dir, root+"_allAbove10.0_filtered.cnv", root+"_allAbove10.0_filtered_inExons.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_EXONS, false);

//		FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new Project(cnv.Launch.getDefaultDebugProjectFile(), false), dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_3dels_LT2dups.cnv", 0, 3, 0, Integer.MAX_VALUE, 1);
//		FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new Project(cnv.Launch.getDefaultDebugProjectFile(), false), dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_3dups_LT2dels.cnv", 0, 0, 3, 1, Integer.MAX_VALUE);
//		FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new Project(cnv.Launch.getDefaultDebugProjectFile(), false), dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_2dels_2dups.cnv", 0, 2, 2, Integer.MAX_VALUE, Integer.MAX_VALUE);
//		FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new Project(cnv.Launch.getDefaultDebugProjectFile(), false), dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_3anythings.cnv", 3, 0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
//		FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new Project(cnv.Launch.getDefaultDebugProjectFile(), false), dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_5anythings.cnv", 5, 0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
//		FilterCalls.filterOnSegments(dir+root+"_0kb_5SNP_10.0_3anythings.cnv", dir+root+"_CommonInGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, true);
		
//		FilterCalls.union(dir+root+"_0kb_5SNP_10.0_3anythings.cnv", dir+root+"_100kb_20SNP_10.0_CNPstatusIgnored.cnv", dir+"unionOfConservativeAndCommon.cnv");
//		FilterCalls.filterOnSegments(dir+"unionOfConservativeAndCommon.cnv", dir+"unionOfConservativeAndCommonInGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, true);
	}
	
	public static void fromParameters(String filename, Logger log) {
		Vector<String> params;
		String dir = "";

		if (Files.exists("N:/statgen/NCBI/")) {
			dir = "N:/statgen/NCBI/";
		}

		params = Files.parseControlFile(filename, "filterCNVs", new String[] {
				"dir=",
				"in=penncnv.cnv",
				"out=conf15used.cnv",
				"# minimum size of a deletions / duplications (in kb):",
				"delSize=0",
				"dupSize=0",
				"# minimum number of SNPs:",
				"number=15",
				"minScore=10.0",
				"filterFile="+dir+"problematicRegions_hg19.dat",
				"# pedfile to be used as a filter:",
				"ped=plink.fam",
				"# if CNV spans centromere, break into two spanning actual markers",
				"breakCentromere=true",
				"# make a UCSC track (.bed file) as well",
				"ucsc=true",
				"",
				"# ALTERNATIVELY, in addition to the dir/in/out and ignoring all other filters you can",
				"# keep only CNVs overlapping these segments (simply uncomment the following argument):",
				"#segs=gene_region.dat",
				"# exclude instead of include:",
				"#excludeSegsInstead=true"
		}, log);

		if (params != null) {
			params.add("log=" + log.getFilename());
			main(Array.toStringArray(params));
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		int delSize = DEFAULT_MIN_SIZE_KB;
		int dupSize = DEFAULT_MIN_SIZE_KB;
		int number = DEFAULT_MIN_NUM_SNPS;
		double score = DEFAULT_MIN_SCORE;
		String filenameOfProblematicRegions = null;
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
		String logfile = null;
		int build = 37;
		String markerSetFilenameToBreakUpCentromeres = null;

		String usage = 
		"vis.cnv.FilterCalls requires 0-1 arguments\n"+
		"   (1) directory (i.e. dir="+dir+" (default))\n"+
		"   (2) file in (i.e. in="+in+" (default))\n"+
		"   (3) file out (i.e. out="+out+" (default))\n"+
		"   (4) minimum size of a deletion (in kb) (i.e. delSize="+delSize+" (default))\n"+
		"   (5) minimum size of a duplication (in kb) (i.e. dupSize="+dupSize+" (default))\n"+
		"   (6) minimum number of SNPs (i.e. number="+number+" (default))\n"+
		"   (7) minimum score (i.e. minScore="+score+" (default))\n"+
		"   (8) filter out cnvs in known problematicRegions (i.e. filterFile="+filenameOfProblematicRegions+" (default))\n"+
		"   (9) pedfile to use as a filter (i.e. ped="+pedfile+" (default))\n"+
		"   (10) if CNV spans centromere, break into two spanning actual markers (i.e. breakCentromere="+breakCent+" (default))\n"+
		"   (11) build of the genome to use for centromeres (i.e. build="+build+" (default))\n"+
		"   (12) custom marker set to determine the last and first marker of the centromeres (i.e. markerFile=plink.bim (not the default))\n"+
		"   (13) make UCSC track as well (i.e. ucsc=true (default))\n"+
		"  OR\n"+
		"   (1) keep only CNVs overlapping these segments (i.e. segs=gene_region.dat (not the default))\n"+
		"   (2) exlcude instead of include (i.e. excludeSegsInstead=false (default))\n"+
//		"   (3) require a significant overlap to filter (i.e. -sigOverlap (not the default))\n"+
		"  OR\n"+
		"   (1) perform all standard filters (i.e. -std (not the default))\n"+
		"   (2) make UCSC tracks as well (i.e. ucsc=false (default))\n"+
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
				delSize = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("dupSize=")) {
				dupSize = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("number=")) {
				number = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("minScore=")) {
				score = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("cnps=")) {
				inOutIgnore = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("filterFile=")) {
				filenameOfProblematicRegions = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-std")) {
				standards = true;
				numArgs--;
			} else if (args[i].startsWith("ucsc=")) {
				tracks = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("segs=")) {
				segs = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("breakCentromere=")) {
				breakCent = ext.parseBooleanArg(args[i]);
				numArgs--;				
			} else if (args[i].startsWith("markerFile=")) {
				markerSetFilenameToBreakUpCentromeres = ext.parseStringArg(args[i], null);
				numArgs--;				
			} else if (args[i].startsWith("excludeSegsInstead=")) {
				excludeSegs = ext.parseBooleanArg(args[i]);
				numArgs--;				
			} else if (args[i].startsWith("ped=")) {
				pedfile = ext.parseStringArg(args[i], null);
				numArgs--;				
			} else if (args[i].startsWith("build=")) {
				build = ext.parseIntArg(args[i]);
				numArgs--;				
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], null);
				numArgs--;				
			} else {
				System.err.println("Error - don't know what to do with argument: "+args[i]);
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			return;
		}

//		dir = "D:/data/GEDI/penn_results/custom_gediBoth/";
//		filterBasedOnNumberOfCNVsAtLocus(new Project(cnv.Launch.getDefaultDebugProjectFile(), false), dir+"conf15_usedFiltered.cnv", dir+"conf15_usedFilteredRare.cnv", 0, 0, 0, 275, 275, 275, 0.50);
//		System.exit(1);
		
//		filter("D:/data/GEDI/penn_results/custom_gediBoth/", "penncnv.cnv", "conf1checkers.cnv", 0, 0, 1, 0, null, -1, "plink.fam", true, null, true, 37, new Logger());
//		System.exit(1);
		
//		FilterCalls.filterOnSegments("D:/data/GEDI/global/homoDelsOverlappingGenesOnly/", "conf.cnv", "conf_overlappingGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, false);
//		FilterCalls.filterOnSegments("D:/data/GEDI/penn_results/custom_gediBoth/conf15_usedFilteredRare/homoDels/", "conf.cnv", "conf_overlappingGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, false);
//		System.exit(1);
		
//		breakCent = true;
//		out = "noCentromeric.cnv";
		
		try {
//			FilterCalls.filterOnSegments(dir+"conf_100kb_20SNP_10.0_CNPstatusIgnored.cnv", dir+"ConservativeGeneCentric.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, true);
//			MakeUCSCtrack.makeTrack(dir, "ConservativeGeneCentric.cnv");
			if (standards) {
				stdFilters(dir, in, tracks, pedfile, build);
			} else if (!segs.equals("")) {
				filterOnSegments(dir, in, out, segs, excludeSegs);
			} else {
				filter(dir, in, out, delSize, dupSize, number, score, filenameOfProblematicRegions, DEFAULT_COMMON_IN_OUT_OR_IGNORED, pedfile, breakCent, markerSetFilenameToBreakUpCentromeres, tracks, build, new Logger(logfile));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
