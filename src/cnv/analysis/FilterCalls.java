package cnv.analysis;

import java.io.*;
import java.util.ArrayList;
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
	
	
	private static class CNVFilterNode {
		public static final int HET_DEL = 0;
		public static final int HOM_DEL = 1;
		public static final int HET_DUP = 2;
		public static final int HOM_DUP = 3;
		
		final CNVariant cnv;
		final int popCnt;
		final ArrayList<CNVariant> major = new ArrayList<CNVariant>();
		final ArrayList<CNVariant> minor = new ArrayList<CNVariant>();
		
		public CNVFilterNode(CNVariant myCNV, int pop) {
			this.cnv = myCNV;
			this.popCnt = pop;
		}
		
		public void addMajor(CNVariant cnv) {
			this.major.add(cnv);
		}
		
		public void addMinor(CNVariant cnv) {
			this.minor.add(cnv);
		}
		
		public int countType(int type, boolean lookMajor) {
			int lookingFor = -1;
			switch(type) {
				case HET_DEL:
					lookingFor = 1;
					break;
				case HOM_DEL:
					lookingFor = 0;
					break;
				case HET_DUP:
					lookingFor = 3;
					break;
				case HOM_DUP:
					lookingFor = 4;
					break;
				default: 
					// TODO show error message
					return 0;
			}
			
			int cnt = 0;
			ArrayList<CNVariant> pop = lookMajor ? major : minor;
			for (CNVariant cnv : pop) {
				// TODO output warning if CN > 4
				if (cnv.getCN() == lookingFor || (lookingFor == 4 && cnv.getCN() > lookingFor)) {
					cnt++;
				}
			}
			
			return cnt;
		}
		
		public String toString() {
			System.out.println("Major: " + major.size() + " of " + popCnt + " - " + ext.formDeci((((double)(major.size()) / ((double)popCnt)) * 100.0), 3, true) + "%");
			System.out.println("Minor: " + minor.size() + " of " + popCnt + " - " + ext.formDeci((((double)(minor.size()) / ((double)popCnt)) * 100.0), 3, true) + "%");

			return cnv.toPlinkFormat() + "\t"
					+ ext.formDeci((((double)(major.size()) / ((double)popCnt)) * 100.0), 3, true) + "\t"
					+ ext.formDeci((((double)(minor.size()) / ((double)popCnt)) * 100.0), 3, true) + "\t"
					+ "(" + this.countType(HOM_DEL, true) + ","
						  + this.countType(HET_DEL, true) + ","
						  + this.countType(HET_DUP, true) + ","
						  + this.countType(HOM_DUP, true) + ")\t"
					+ "(" + this.countType(HOM_DEL, false) + ","
					  + this.countType(HET_DEL, false) + ","
					  + this.countType(HET_DUP, false) + ","
					  + this.countType(HOM_DUP, false) + ")";
		}
		
	}
	
	public static void CNVStats(Project proj, String dir, String filenameNoExt) {
		// find .cnv and .fam file from fileroot
		
		
		
	}
	
	public static void filterLists(String cnvList, String[] cnvFiles, String outputFile) {
		CNVariant[] srcCNVs = CNVariant.loadPlinkFile(cnvList, false);
		ArrayList<CNVariant> compCNVs = new ArrayList<CNVariant>();
		HashSet<String> ids = new HashSet<String>();
		for (String file : cnvFiles) {
			CNVariant[] cnvs = CNVariant.loadPlinkFile(file, false);
			for (CNVariant cnv : cnvs) {
				compCNVs.add(cnv);
				ids.add(cnv.getFamilyID() + "\t" + cnv.getIndividualID());
			}
		}
		
		
		ArrayList<CNVFilterNode> outputNodes = new ArrayList<FilterCalls.CNVFilterNode>();
		
		for (CNVariant cnv : srcCNVs) {
			CNVFilterNode cnvNode = new CNVFilterNode(cnv, ids.size());
			outputNodes.add(cnvNode);
			
			for (CNVariant comp : compCNVs) {
				if (cnv.equalsIncludingIndividual(comp)) continue;
				int overlap = cnv.amountOfOverlapInBasepairs(comp);
				if (overlap == -1) continue;
				if (overlap >= (cnv.getSize() / 2)) {
					cnvNode.addMajor(comp);
				} else {
					cnvNode.addMinor(comp);
				}
			}
			
		}
		
		String header = Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t") + "\t%major\t%minor\tStats(M)\tStats(m)";
		
		PrintWriter writer;
		try {
			writer = new PrintWriter(new FileWriter(outputFile));
			writer.println(header);
			for (CNVFilterNode node : outputNodes) {
				writer.println(node.toString());
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void filterList(String cnvList, String cnvFile, String outputFile) {	
		/*
		 * INPUT (PLINK formatted):
		 * |  FID  |  IID  |  CHR  |  BP1  |  BP2  |  TYPE  |  SCORE  |  SITES  | 
		 * 
		 * OUTPUT:
		 * |  FID  |  IID  |  CHR  |  BP1  |  BP2  |  TYPE  |  SCORE  |  SITES  |  % cnvs overlapping >50% of source cnv  |  % cnvs overlapping <50% of source cnv  |   of major overlapping, # of hDel, hmDel, hDup, hmDup  |  of minor overlapping, # of hDel, hmDel, hDup, hmDup  |
		 */
		
		CNVariant[] srcCNVs = CNVariant.loadPlinkFile(cnvList, false);
		CNVariant[] compCNVs = CNVariant.loadPlinkFile(cnvFile, false);
		
		HashSet<String> ids = new HashSet<String>();
		for (CNVariant cnv : compCNVs) {
			ids.add(cnv.getFamilyID() + "\t" + cnv.getIndividualID());
		}
		
		ArrayList<CNVFilterNode> outputNodes = new ArrayList<FilterCalls.CNVFilterNode>();
		
		for (CNVariant cnv : srcCNVs) {
			CNVFilterNode cnvNode = new CNVFilterNode(cnv, ids.size());
			outputNodes.add(cnvNode);
			
			for (CNVariant comp : compCNVs) {
				if (cnv.equalsIncludingIndividual(comp)) continue;
				int overlap = cnv.amountOfOverlapInBasepairs(comp);
				if (overlap == -1) continue;
				if (overlap >= (cnv.getSize() / 2)) {
					cnvNode.addMajor(comp);
				} else {
					cnvNode.addMinor(comp);
				}
			}
			
		}
		
		String header = Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t") + "\t%major\t%minor\tStats(M)\tStats(m)";
		
		PrintWriter writer;
		try {
			writer = new PrintWriter(new FileWriter(outputFile));
			writer.println(header);
			for (CNVFilterNode node : outputNodes) {
				writer.println(node.toString());
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void filterExclusions(Project proj, String cnvFile) {
		PrintWriter writer;
		BufferedReader reader;
		String path = cnvFile.substring(0, cnvFile.lastIndexOf('/'));
		String fileext = cnvFile.substring(cnvFile.lastIndexOf('/'));
		String filenm = fileext.substring(0, fileext.lastIndexOf('.'));
		String newFile = path + filenm + ".excluded.cnv";
		
		String[] excludes = Array.subArray(proj.getSamples(), proj.getSamplesToExclude());
		HashSet<String> excludeSet = new HashSet<String>();
		for (String exclude : excludes) excludeSet.add(exclude);
		
		try {
			reader = new BufferedReader(new FileReader(cnvFile));
			writer = new PrintWriter(new FileWriter(newFile));
			
			while(reader.ready()) {
				String line = reader.readLine();
				if (excludeSet.contains(line.split("\t")[0])) continue;
				writer.println(line);
			}
			
			writer.flush();
			
			reader.close();
			writer.close();
		} catch (FileNotFoundException e) {
			proj.getLog().reportException(e);
		} catch (IOException e) {
			proj.getLog().reportException(e);
		}
	}
	

	public static void filterExclusions(String dir, String in, String out, String excludeFile) {
		PrintWriter writer;
		Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(dir + in, null, false);
		boolean[] remove = new boolean[cnvs.size()];
		HashSet<String> excludeList = HashVec.loadFileToHashSet(excludeFile, false);
		
		for (int i = 0; i < remove.length; i++) {
			CNVariant examine = cnvs.get(i);
			if (excludeList.contains(examine.getIndividualID())) { // remove all CNVs belonging to a listed sample ID
				remove[i] = true;
				continue;
			}
			
		}

		try {
			writer = new PrintWriter(new FileWriter(dir + out));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
			for (int i = 0; i < remove.length; i++) {
				if (!remove[i]) writer.println(cnvs.get(i).toPlinkFormat());
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	/**
	 * Given a list of CNVs, and a list of sample IDs, find all CNVs that overlap at least one CNV belonging to a sample ID on the given list.
	 * 
	 * Useful for data without controls, but with known affected samples.
	 * 
	 * @param dir location of in / out files
	 * @param in Plink-formatted CNV file for input
	 * @param out desired name of filtered CNV file
	 * @param listFile file name (full path) of the list of sample IDs
	 */
	public static void filterForAllCNVsSharedInGroup(String dir, String in, String out, String listFile) {
		PrintWriter writer;
		Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(dir + in, null, false);
		boolean[] remove = new boolean[cnvs.size()];
		HashSet<String> indivList = HashVec.loadFileToHashSet(listFile, false);
		
		for (int i = 0; i < remove.length; i++) {
			CNVariant examine = cnvs.get(i);
			
			boolean mark = true;
			for (CNVariant groupCNV : cnvs) {
				if (!indivList.contains(groupCNV.getIndividualID())) { // loop through only those belonging to listed sample IDs
					continue;
				}
				
				if (examine.overlaps(groupCNV)) {
					mark = false;
					break;
				}
			}
			remove[i] = mark;
		}
		
		try {
			writer = new PrintWriter(new FileWriter(dir + out));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
			for (int i = 0; i < remove.length; i++) {
				if (!remove[i]) writer.println(cnvs.get(i).toPlinkFormat());
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public static void filterForGroupCNVs(String dir, String in, String out, String listFile, String excludeFile, boolean excludeCommon) {
		PrintWriter writer;
		Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(dir + in, null, false);
		boolean[] remove = new boolean[cnvs.size()];
		HashSet<String> indivList = HashVec.loadFileToHashSet(listFile, false);
		HashSet<String> excludeList = excludeFile == null ? new HashSet<String>() : HashVec.loadFileToHashSet(excludeFile, false);
		
		for (int i = 0; i < remove.length; i++) {
			CNVariant examine = cnvs.get(i);
			if (!indivList.contains(examine.getIndividualID()) || excludeList.contains(examine.getIndividualID())) { // remove all CNVs not belonging to a listed sample ID
				remove[i] = true;
				continue;
			}
			
			if (excludeCommon) {
				for (int j = 0; j < remove.length; j++) {
					if (i == j) continue;
					CNVariant check = cnvs.get(j);
					if (examine.significantOverlap(check, true) && !indivList.contains(check.getIndividualID())) {
						// if the sample ID of the second CNV isn't on the list, and the current cnv is on the list, remove both
						remove[i] = true;
						remove[j] = true;
					}
				}
			}
			
		}

		try {
			writer = new PrintWriter(new FileWriter(dir + out));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
			for (int i = 0; i < remove.length; i++) {
				if (!remove[i]) writer.println(cnvs.get(i).toPlinkFormat());
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	
	public static void filterOutCommonCNVs(String dir, String in, String out, double pct, boolean checkLarger) {
		PrintWriter writer;
		Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(dir + in, null, false);
		boolean[] remove = new boolean[cnvs.size()];
		
		for (int i = 0; i < remove.length; i++) {
			CNVariant examine = cnvs.get(i);
			int overlapCnt = 0;
			
			for (int j = 0; j < remove.length; j++) {
				if (i == j) continue;
				
				if (examine.significantOverlap(cnvs.get(j), checkLarger)) {
					overlapCnt++;
				}
			}
			
			if (((double)overlapCnt) / ((double)cnvs.size()) > pct) {
				remove[i] = true;
			} else {
				remove[i] = false;
			}
		}
		
		try {
			writer = new PrintWriter(new FileWriter(dir + out));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
			for (int i = 0; i < remove.length; i++) {
				if (!remove[i]) writer.println(cnvs.get(i).toPlinkFormat());
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void filterOutCommonCNVs(String dir, String in, String out, double pct) {
		PrintWriter writer;
		Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(dir + in, null, false);
		boolean[] remove = new boolean[cnvs.size()];
		
		for (int i = 0; i < remove.length; i++) {
			CNVariant examine = cnvs.get(i);
			int overlapCnt = 0;
			
			for (int j = 0; j < remove.length; j++) {
				if (i == j) continue;
				
				if (examine.significantOverlap(cnvs.get(j))) {
					overlapCnt++;
				}
			}
			
			if (((double)overlapCnt) / ((double)cnvs.size()) > pct) {
				remove[i] = true;
			} else {
				remove[i] = false;
			}
		}
		
		try {
			writer = new PrintWriter(new FileWriter(dir + out));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
			for (int i = 0; i < remove.length; i++) {
				if (!remove[i]) writer.println(cnvs.get(i).toPlinkFormat());
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static void filter(String dir, String in, String out, int[] delSize, int[] dupSize, int[] number, double score, String filenameOfProblematicRegions, int commonInOutOrIgnore, String individualsToKeepFile, boolean breakupCentromeres, String markerSetFilenameToBreakUpCentromeres, boolean makeUCSCtrack, int build, Logger log) {
		String[] individualsToKeepList;
		if (individualsToKeepFile != null && !new File(individualsToKeepFile).exists()) {
			System.err.println("Error - could not find \""+individualsToKeepFile+"\" in directory; will not be able to filter by indiviudals");
			individualsToKeepFile = null;
		}
		individualsToKeepList = individualsToKeepFile==null?null:HashVec.loadFileToStringArray(individualsToKeepFile, false, false, new int[] {0,1}, true, false, "\t"); 

		filter(dir, in, out, delSize, dupSize, number, score, filenameOfProblematicRegions, commonInOutOrIgnore, individualsToKeepList, breakupCentromeres, markerSetFilenameToBreakUpCentromeres, makeUCSCtrack, build, log);
	}
	
	public static void filter(String dir, String in, String out, int[] delSize, int[] dupSize, int[] number, double score, String filenameOfProblematicRegions, int commonInOutOrIgnore, String[] individualsToKeepList, boolean breakupCentromeres, String markerSetFilenameToBreakUpCentromeres, boolean makeUCSCtrack, int build, Logger log) {
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
				if (	(	(cnv.getCN() == 1 && cnv.getSize() >= delSize[0] * 1000) || // heterozygous deletion
							(cnv.getCN() == 0 && cnv.getSize() >= delSize[1] * 1000) || // homozygous deletion
							(cnv.getCN() > 2 && cnv.getSize() >= dupSize[0] * 1000) // duplications
							// ignoring homozygotic duplications
						) 
						&& 
							(((cnv.getCN() == 1 || cnv.getCN() == 3 || cnv.getCN() == 4) && cnv.getNumMarkers() >= number[0]) ||
							(cnv.getCN() == 2 && cnv.getNumMarkers() >= number[1]))
						&& cnv.getScore() > score 
						&& !inOneOfTheseRegions(cnv, problemRegions)) {

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
		FilterCalls.filter(dir, filename, root+"_allAbove10.0_unfiltered.cnv", new int[]{1, 1}, new int[]{1, 1}, new int[]{1, 1}, 10, null, COMMON_IGNORED, pedfile, true, null, makeUCSCtracks,  build, log);
		FilterCalls.filter(dir, filename, root+"_allAbove10.0_filtered.cnv", new int[]{1, 1}, new int[]{1, 1}, new int[]{1, 1}, 10, DEFAULT_PROBLEMATIC_REGIONS, COMMON_IGNORED, pedfile, true, null, makeUCSCtracks, build, log);
		FilterCalls.filter(dir, filename, root+"_ConservativeCalls.cnv", new int[]{100, 100}, new int[]{100, 100}, new int[]{20, 20}, 10, DEFAULT_PROBLEMATIC_REGIONS, COMMON_IGNORED, pedfile, true, null, makeUCSCtracks, build, log);

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
				"# minimum size of heterozygous deletions / duplications (in kb):",
				"delSize=0",
				"dupSize=0",
				"# minimum size of homozygous deletions / duplications (in kb):",
				"hDelSize = 0",
				"hDupSize = 0",
				"# minimum number of heterozygous SNPs:",
				"number=15",
				"# minimum number of homozygous SNPs:",
				"hNumber=15",
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
		int hDelSize = DEFAULT_MIN_SIZE_KB;
		int hDupSize = DEFAULT_MIN_SIZE_KB;
		int number = DEFAULT_MIN_NUM_SNPS;
		int hNumber = DEFAULT_MIN_NUM_SNPS;
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
		boolean common = false;
		double pct = 0.05;
		String listFile = null;
		String[] listFiles = null;
		String excludeFile = null;
		boolean exclude = false;
		boolean group = false;
		boolean stats = false;

		String usage = 
		"vis.cnv.FilterCalls requires 2+ arguments\n"+
		"   (1) directory (i.e. dir="+dir+" (default))\n"+
		"   (2) file in (i.e. in="+in+" (default))\n"+
		"   (3) file out (i.e. out="+out+" (default))\n"+
		"   (4) minimum size of a deletion (in kb) (i.e. delSize="+delSize+" (default))\n"+
		"   (5) minimum size of a duplication (in kb) (i.e. dupSize="+dupSize+" (default))\n"+
		"   (6) (Optional) minimum size of a homozygous deletion (in kb) (i.e. hDelSize="+hDelSize+" (default))\n"+
		"   (7) (Optional) minimum size of a homozygous duplication (in kb) (i.e. hDupSize="+hDelSize+" (default))\n"+
		"   (8) minimum number of heterozygous SNPs (i.e. number="+number+" (default))\n"+
		"   (9) (Optional) minimum number of homozygous SNPs (i.e. hNumber=" + number + " (default))\n" + 
		"   (9) minimum score (i.e. minScore="+score+" (default))\n"+
		"   (10) filter out cnvs in known problematicRegions (i.e. filterFile="+filenameOfProblematicRegions+" (default))\n"+
		"   (11) pedfile to use as a filter (i.e. ped="+pedfile+" (default))\n"+
		"   (12) if CNV spans centromere, break into two spanning actual markers (i.e. breakCentromere="+breakCent+" (default))\n"+
		"   (13) build of the genome to use for centromeres (i.e. build="+build+" (default))\n"+
		"   (14) custom marker set to determine the last and first marker of the centromeres (i.e. markerFile=plink.bim (not the default))\n"+
		"   (15) make UCSC track as well (i.e. ucsc=true (default))\n"+
		"  OR\n"+
		"   (1) keep only CNVs overlapping these segments (i.e. segs=gene_region.dat (not the default))\n"+
		"   (2) exclude instead of include (i.e. excludeSegsInstead=false (default))\n"+
//		"   (3) require a significant overlap to filter (i.e. -sigOverlap (not the default))\n"+
		"  OR\n"+
		"   (1) perform all standard filters (i.e. -std (not the default))\n"+
		"   (2) make UCSC tracks as well (i.e. ucsc=false (default))\n"+
		"  OR\n" + 
		"   (1) directory (i.e. dir="+dir+" (default))\n"+
		"   (2) file in (i.e. in="+in+" (default))\n"+
		"   (3) file out (i.e. out="+out+" (default))\n"+
		"   (4) filter out CNVs that overlap with a percentage of other CNVs (i.e. -common (not the default))\n" +
		"   (5) percent threshold for filtering common CNVs (i.e. pct=.05 (default))\n" +
		"  OR\n" + 
		"   (1) directory (i.e. dir="+dir+" (default))\n"+
		"   (2) file in (i.e. in="+in+" (default))\n"+
		"   (3) file out (i.e. out="+out+" (default))\n"+
		"   (4) remove all CNVs not belonging to a group of sample IDs (i.e. list=/path/to/list.txt (not the default))\n" +
		"   (5) (Optional) remove all CNVs that overlap other CNVs not belonging to the given list of sample IDs (i.e. -exclude (not the default))" + 
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
			} else if (args[i].startsWith("hDelSize=")) {
				hDelSize = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("hDupSize=")) {
				hDupSize = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("number=")) {
				number = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("hNumber=")) {
				hNumber = Integer.parseInt(args[i].split("=")[1]);
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
			} else if (args[i].startsWith("pct=")) {
				pct = ext.parseDoubleArg(args[i]);
				numArgs--;				
			} else if (args[i].startsWith("-common")) {
				common = true;
				numArgs--;				
			} else if (args[i].startsWith("-exclude")) {
				exclude = true;
				numArgs--;				
			} else if (args[i].startsWith("-group")) { 
				group = true;
				numArgs--;
			} else if (args[i].startsWith("-stats")) { 
				stats = true;
				numArgs--;
			} else if (args[i].startsWith("list=")) {
				listFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("lists=")) {
				listFiles = args[i].split("=")[1].split(",");
				numArgs--;
			}  else if (args[i].startsWith("excludeFile=")) { 
				excludeFile = args[i].split("=")[1];
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
			} else if (common) {
				filterOutCommonCNVs(dir, in, out, pct);
			} else if (listFile != null) { 
				if (group) {
					filterForAllCNVsSharedInGroup(dir, in, out, listFile);
				} else if (stats) {
					filterList(in, listFile, out);
				} else {
					filterForGroupCNVs(dir, in, out, listFile, excludeFile, exclude);
				}
			} else if (listFiles != null) {
				filterLists(in, listFiles, out);
			} else if (excludeFile != null) {
				filterExclusions(dir, in, out, excludeFile);
			} else {
				filter(dir, in, out, new int[]{delSize, hDelSize}, new int[]{dupSize, hDupSize}, new int[]{number, hNumber}, score, filenameOfProblematicRegions, DEFAULT_COMMON_IN_OUT_OR_IGNORED, pedfile, breakCent, markerSetFilenameToBreakUpCentromeres, tracks, build, new Logger(logfile));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
