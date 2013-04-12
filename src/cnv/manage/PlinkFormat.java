// -Xms1024M -Xmx1024M
package cnv.manage;

import java.io.*;
import java.util.*;

import cnv.filesys.FullSample;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import common.*;

public class PlinkFormat {
	public static void createPlink(Project proj) {
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
		FullSample fsamp;
		byte[] genotypes;
		byte genIndex;
		String genotype;
		int count;
		double gcThreshold;
		float[] gcs;
		String targetMarkers;
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
		
		targetMarkers = proj.getFilename(Project.TARGET_MARKERS_FILENAME);
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
				System.err.println("Error - target markers in file '"+targetMarkers+"' not found, using all markers");
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
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"gwas.map"));
			for (int i = 0; i<indices.length; i++) {
				writer.println(chrs[indices[i]]+"\t"+markerNames[indices[i]]+"\t0\t"+positions[indices[i]]);
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: failed to write to gwas.map (in use?)");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error writing to gwas.map");
			System.exit(2);
		}

		gcThreshold = proj.getDouble(Project.GC_THRESHOLD);
		if (gcThreshold < 0) {
			System.err.println("Error - GC_THRESHOLD must be greater than zero (not "+gcThreshold+")");
		}
		if (gcThreshold > 1) {
			System.err.println("Error - GC_THRESHOLD must be less than one (not "+gcThreshold+")");
		}
		System.out.println("Using a GC threshold of "+gcThreshold+" (less than or equal to will be set to missing, greater than is kept)");
		try {
			reader = new BufferedReader(new FileReader(proj.getFilename(Project.PEDIGREE_FILENAME)));
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"gwas.ped"));
			count = 1;
			while (reader.ready()) {
				System.out.println(count++);
				temp = reader.readLine();
				line = temp.split(ext.determineDelimiter(temp));
				writer.print(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]);
				if (line[6].equals(".")) {
					for (int i = 0; i<indices.length; i++) {
						writer.print("\t0\t0");
					}
				} else {
					fsamp = proj.getFullSample(line[6]);
					if (fsamp==null) {
						System.err.println("Error - the DNA# "+line[6]+" was listed in the pedigree file but "+line[6]+".fsamp was not found in directory: "+proj.getDir(Project.SAMPLE_DIRECTORY));
						for (int i = 0; i<indices.length; i++) {
							writer.print("\t0\t0");
						}
					} else {
						genotypes = fsamp.getForwardGenotypes();
						gcs = fsamp.getGCs();
						for (int i = 0; i<indices.length; i++) {
							genIndex = genotypes[indices[i]];
							if (genIndex==0 || gcs[indices[i]] <= gcThreshold) {
								writer.print("\t0\t0");
							} else {
								genotype = FullSample.ALLELE_PAIRS[genIndex];
								writer.print("\t"+genotype.charAt(0)+"\t"+genotype.charAt(1));
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

		String usage = "\\n"+
		"cnv.manage.PlinkFormat requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"     requires pedigree file and targets file to be delineated in the Project properties file\n"+
		"  OR\n"+
		"   (1) pick targets to include in the PLINK file from an Illumina map file (i.e. pick=filename.csv (not the default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
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
				createPlink(new Project(filename, false));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
