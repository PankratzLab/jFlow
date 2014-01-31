// -Xms1024M -Xmx1024M
package cnv.manage;

import java.io.*;
import java.util.*;

import cnv.filesys.ABLookup;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.Sample;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import common.*;

public class PlinkFormat {
	public static boolean createPlink(Project proj, String filenameRoot, String clusterFilterFilename, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		Hashtable<String,String> hash;
		String[] line;
		int[] indices;
		String[] markerNames, targets;
		byte[] chrs;
		int[] positions;
		boolean prob;
		MarkerSet markerSet;
		Sample fsamp;
		byte[] genotypes;
		byte genIndex;
		String genotype;
		int count;
		float gcThreshold;
		String targetMarkers;
		ClusterFilterCollection clusterFilterCollection;
		char[][] abLookup;
		String temp;

		System.out.println(ext.getTime());
		hash = new Hashtable<String,String>();
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		
		for (int i = 0; i<markerNames.length; i++) {
			if (hash.containsKey(markerNames[i])) {
				System.err.println("Warning - duplicate marker name: "+markerNames[i]);
			}
			hash.put(markerNames[i], i+"");
		}
		
		targetMarkers = proj.getFilename(Project.TARGET_MARKERS_FILENAME, false, false);
		if (new File(targetMarkers).exists()) {
			targets = HashVec.loadFileToStringArray(targetMarkers, false, false, new int[] {0}, false);
			indices = new int[targets.length];
			prob = false;
			for (int i = 0; i<targets.length; i++) {
				if (hash.containsKey(targets[i])) {
					indices[i] = Integer.parseInt(hash.get(targets[i]));
				} else {
					System.err.println("Error - target marker '"+targets[i]+"' was not found in the MarkerSet");
					prob = true;
				}
			}
			if (prob) {
				System.exit(1);
			}
		} else {
			if (!targetMarkers.equals("")) {
				System.out.println("FYI, since target markers file '"+targetMarkers+"' was not found, all markers will be exported to PLINK");
			}

//			targets = HashVec.getKeys(hash);
//			indices = new int[targets.length];
//			for (int i = 0; i<targets.length; i++) {
//				indices[i] = Integer.parseInt(hash.get(targets[i]));
//			}
			indices = Array.intArray(markerNames.length);
		}

		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+filenameRoot+".map"));
			for (int i = 0; i<indices.length; i++) {
				writer.println(chrs[indices[i]]+" "+markerNames[indices[i]]+" 0 "+positions[indices[i]]);
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: failed to write to "+filenameRoot+".map (in use?)");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error writing to "+filenameRoot+".map");
			System.exit(2);
		}

		gcThreshold = proj.getFloat(Project.GC_THRESHOLD);
		if (gcThreshold < 0) {
			System.err.println("Error - GC_THRESHOLD must be greater than zero (not "+gcThreshold+")");
		}
		if (gcThreshold > 1) {
			System.err.println("Error - GC_THRESHOLD must be less than one (not "+gcThreshold+")");
		}
		System.out.println("Using a GC threshold of "+gcThreshold+" (less than or equal to will be set to missing, greater than is kept)");
		
		clusterFilterCollection = null;
		if (clusterFilterFilename != null) {
			clusterFilterFilename = proj.getProperty(Project.PROJECT_DIRECTORY)+proj.getProperty(Project.DATA_DIRECTORY)+clusterFilterFilename;
			if (Files.exists(clusterFilterFilename, proj.getJarStatus())) {
				clusterFilterCollection = ClusterFilterCollection.load(clusterFilterFilename, proj.getJarStatus());
			} else {
				System.err.println("Error - cluster filter collection is not found at '"+clusterFilterFilename+"'");
				return false;
			}
			abLookup = new ABLookup(markerNames, proj.getFilename(Project.AB_LOOKUP_FILENAME), true, true).getLookup();
			log.report("Using "+clusterFilterFilename+" and "+proj.getProperty(Project.AB_LOOKUP_FILENAME)+" to call genotypes");
		} else {
			abLookup = null;
		}
		
		try {
			reader = new BufferedReader(new FileReader(proj.getFilename(Project.PEDIGREE_FILENAME)));
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+filenameRoot+".ped"));
			count = 1;
			while (reader.ready()) {
				count++;
				if (count % 100 == 0) {
					System.out.println(count);
				}
				
				temp = reader.readLine();
				line = temp.split(ext.determineDelimiter(temp));
				if (line.length < 7) {
					log.reportError("Error - starting at line "+(count-1)+(line.length<3?"":" (individual "+line[0]+"-"+line[1]+")")+" there are only "+line.length+" columns in pedigree file '"+proj.getFilename(Project.PEDIGREE_FILENAME)+"'.");
					log.reportError("  Pedigree files require 7 columns with no header: FID IID FA MO SEX PHENO DNA");
					log.reportError("  where DNA is the sample name associated with the genotypic data (see the "+proj.getDir(Project.SAMPLE_DIRECTORY)+" directory for examples)");
					reader.close();
					writer.close();
					return false;
				}
				writer.print(line[0]+" "+line[1]+" "+line[2]+" "+line[3]+" "+line[4]+" "+line[5]);
				if (line[6].equals(".")) {
					for (int i = 0; i<indices.length; i++) {
						writer.print(" 0 0");
					}
				} else {
					fsamp = proj.getFullSampleFromRandomAccessFile(line[6]);
					if (fsamp==null) {
						System.err.println("Error - the DNA# "+line[6]+" was listed in the pedigree file but "+line[6]+".fsamp was not found in directory: "+proj.getDir(Project.SAMPLE_DIRECTORY));
						for (int i = 0; i<indices.length; i++) {
							writer.print(" 0 0");
						}
					} else {
						if (clusterFilterFilename == null) {
							genotypes = fsamp.getForwardGenotypes(gcThreshold);
						} else {
							genotypes = MarkerSet.translateABtoForwardGenotypes(fsamp.getAB_GenotypesAfterFilters(markerNames, clusterFilterCollection, gcThreshold), abLookup);
						}
						for (int i = 0; i<indices.length; i++) {
							genIndex = genotypes[indices[i]];
							if (genIndex==0) {
								writer.print(" 0 0");
							} else {
								genotype = Sample.ALLELE_PAIRS[genIndex];
								writer.print(" "+genotype.charAt(0)+" "+genotype.charAt(1));
							}
						}
					}
				}
				writer.println();
				writer.flush();
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.getFilename(Project.PEDIGREE_FILENAME)+"\" not found");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.getFilename(Project.PEDIGREE_FILENAME)+"\"");
			System.exit(2);
		}

		System.out.println(ext.getTime());
		return true;
	}

	public static void pickTargets(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int[] indices;

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter("picksForTargetMarkers.txt"));
			line = reader.readLine().split(",", -1);
			indices = ext.indexFactors(new String[] {"Name", "Intensity Only", "Minor Freq"}, line, false, true);
			while (reader.ready()) {
				line = reader.readLine().split(",", -1);
				if (line[indices[1]].equals("0")&&!line[indices[2]].equals("")&&Double.parseDouble(line[indices[2]])>0) {
					writer.println(line[indices[0]]);
				}
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}
	}
	
	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;
		String pick = "";
		String filters = null;
		String plinkPrefix = "plinkPed";

		String usage = "\\n"+
		"cnv.manage.PlinkFormat requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"     requires pedigree file and targets file to be delineated in the Project properties file\n"+
		"   (2) prefix for plink ped files (i.e. prefix="+plinkPrefix+" (default))\n"+
		"   (3) filename of cluster filters to use during processing (i.e. filters=data/clusterFilters.ser (not the default))\n"+
		"  OR\n"+
		"   (3) pick targets to include in the PLINK file from an Illumina map file (i.e. pick=filename.csv (not the default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("prefix=")) {
				plinkPrefix = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("filters=")) {
				filters = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pick=")) {
				pick = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		try {
			if (!pick.equals("")) {
				pickTargets(pick);
			} else {
				createPlink(new Project(filename, false), plinkPrefix, filters, new Logger());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
