package gwas;

import java.io.*;
import java.util.*;
import common.*;

public class ForestPlot {
	public static final String[] MARKER_ID_POSSIBILITIES = {"Marker", "SNP", "Variant"};
	
	public static final String[][] SUFFIXES = {
		{"A1", "Allele1"},
		{"A2", "Allele2"},
		{"Beta"},
		{"SE", "StdErr"},
		{"Freq", "Freq1"},
	};
	
	public static void fromParameters(String filename, Logger log) {
        String trav;
        String datafile = "hits.txt";
        String markerListFile = null;
        String studyListFile = null;
        Vector<String> params;
        boolean studyFirst = true;

		params = Files.parseControlFile(filename, "forest", new String[] {"hits_parsed.xln", "# optional file lists which SNPs to parse (default is all)", "targets=hits.txt", "# optional file listing order of studies (default is the order they are in the data file) as well as an optional second tab delimited column with what to rename the study as", "order=studies.txt", "# if the study comes after the underscore (i.e. Beta_StudyA instead of StudyA_Beta) then use the following flag:", "studyFirst=false"}, log);
		if (params != null) {
			datafile = params.elementAt(0).trim();
			for (int i = 1; i < params.size(); i++) {
	        	trav = params.elementAt(i).trim();
	        	if (trav.startsWith("targets=")) {
	        		markerListFile = ext.parseStringArg(trav, null);
    		    } else if (trav.startsWith("order=")){
	        		studyListFile = ext.parseStringArg(trav, null);
    		    } else if (trav.startsWith("studyFirst=")){
    		    	studyFirst = ext.parseBooleanArg(trav, log);
    		    } else if (!trav.startsWith("#")){
    		    	log.reportError("Error - don't know what to do with argument: "+trav);
    		    }
	        }

			generateForests(ext.parseDirectoryOfFile(filename), datafile, studyFirst, markerListFile, studyListFile, log);
		}
	}

	public static void generateForests(String dir, String datafile, boolean studyFirst, String markerListFile, String studyListFile, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, data;
		Hashtable<String, String> markers;
		Vector<String> v = new Vector<String>();
		String filename;
        String[] header, roots;
        int[][] indices;
        int markerIndex, index;
		boolean problem;
		String[][] studies;
		
		if (markerListFile != null) {
			markers = HashVec.loadFileToHashNull(dir+markerListFile, false);
		} else {
			markers = null;
		}
        
		if (studyListFile != null) {
			studies = HashVec.loadFileToStringMatrix(dir+studyListFile, false, new int[] {0,1}, "\t", false, 20, true);
		} else {
			studies = null;
		}
		
        try {
			reader = new BufferedReader(new FileReader(dir+datafile));
	        header = reader.readLine().trim().split("[\\s]+");
	        v = new Vector<String>();
	        markerIndex = -1;
	        for (int i = 0; i<header.length; i++) {
	        	if (ext.indexOfStr(header[i], MARKER_ID_POSSIBILITIES, false, true) != -1) {
	        		if (markerIndex == -1) {
	        			markerIndex = i;
	        		} else {
	        			log.reportError("Warning - more than one column is labeled "+Array.toStr(MARKER_ID_POSSIBILITIES, "/")+"; using the first ("+header[markerIndex]+"), not '"+header[i]+"'");
	        		}
	        	}	        		
	        	if (header[i].contains("_")) {
	        		if (studyFirst) {
		        		HashVec.addIfAbsent(header[i].substring(0, header[i].lastIndexOf("_")), v);
	        		} else {
		        		HashVec.addIfAbsent(header[i].substring(header[i].lastIndexOf("_")+1), v);
	        		}
	        	}
            }
	        problem = false;
	        roots = Array.toStringArray(v);
	        if (studies == null) {
	        	studies = Array.toMatrix(roots);
	        } else {
	        	problem = false;
	        	for (int i = 0; i < studies.length; i++) {
	        		if (ext.indexOfStr(studies[i][0], roots) == -1) {
	        			log.reportError("Error - there is no root called '"+studies[i][0]+"' in the datafile");
	        			problem = true;
	        		}
				}
	        	if (problem) {
	        		return;
	        	}
	        }
	        
	        indices = Matrix.intMatrix(studies.length, SUFFIXES.length, -1);
	        for (int i = 0; i<studies.length; i++) {
	        	for (int j = 0; j<SUFFIXES.length; j++) {
	        		for (int k = 0; k < SUFFIXES[j].length; k++) {
	        			index = ext.indexOfStr(studyFirst?studies[i][0]+"_"+SUFFIXES[j][k]:SUFFIXES[j][k]+"_"+studies[i][0], header);
	        			if (index != -1) {
				        	if (indices[i][j] == -1) {
				        		indices[i][j] = ext.indexOfStr(studyFirst?studies[i][0]+"_"+SUFFIXES[j][k]:SUFFIXES[j][k]+"_"+studies[i][0], header);
				        	} else {
			        			log.reportError("Warning - more than one column is labeled "+(studyFirst?studies[i][0]+"_"+Array.toStr(SUFFIXES[j], "/"):Array.toStr(SUFFIXES[j], "/")+"_"+studies[i][0])+"; using the first ("+header[indices[i][j]]+")");
				        	}
	        			}
					}
		        	if (indices[i][j] == -1) {
		        		log.reportError("Error - expecting, but did not find, a column header called "+(studyFirst?studies[i][0]+"_"+Array.toStr(SUFFIXES[j], "/"):Array.toStr(SUFFIXES[j], "/")+"_"+studies[i][0]));
		        		problem = true;
		        	}
                }
            }		
	        if (problem) {
	        	return;
	        }
	        new File(dir+"forests/").mkdirs();
	        
//	        for (int i = 0; i < studies.length; i++) {
//	        	System.out.println(Array.toStr(studies[i]));
//			}
//	        System.out.println();
//	        System.exit(1);

			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (markers == null || markers.containsKey(line[markerIndex])) {
					filename = ext.replaceWithLinuxSafeCharacters(line[markerIndex], true)+".txt";
					try {
						writer = new PrintWriter(new FileWriter(dir+"forests/"+filename));
				        writer.println("Cohort\tSNP\tA1\tA2\tBeta\tSE\tFreq");
				        for (int i = 0; i < studies.length; i++) {
//				        	writer.print((studies[i].length==1||studies[i][1]==null?studies[i][0]:studies[i][1])+"\t"+line[markerIndex]);
////				        	+"\t"+Array.toStr(Array.subArray(line, indices[i])));(studies[i].length==1||studies[i][1]==null?studies[i][0]:studies[i][1])+"\t"+line[markerIndex]+"\t"+Array.toStr(Array.subArray(line, indices[i])));
//				        	writer.println((studies[i].length==1||studies[i][1]==null?studies[i][0]:studies[i][1])+"\t"+line[markerIndex]+"\t"+Array.toStr(Array.subArray(line, indices[i])));
//				        	writer.println((studies[i].length==1||studies[i][1]==null?studies[i][0]:studies[i][1])+"\t"+line[markerIndex]+"\t"+Array.toStr(Array.subArray(line, indices[i])));
				        	data = Array.subArray(line, indices[i]);
				        	if (ext.indexOfStr(".", data) == -1) {
				        		writer.println((studies[i].length==1||studies[i][1]==null?studies[i][0]:studies[i][1])+"\t"+line[markerIndex]+"\t"+Array.toStr(data));
				        	}
						}
						writer.close();
					} catch (Exception e) {
						System.err.println("Error writing to " + filename);
						e.printStackTrace();
					}
					
					filename = ext.replaceWithLinuxSafeCharacters(line[markerIndex], true)+".R";
					try {
						writer = new PrintWriter(new FileWriter(dir+"forests/"+filename));
						writer.println("##Data directory");
						writer.println("data.dir <- \""+dir+"forests/\"");
						writer.println("### Input data");
						writer.println("input.data <- read.table(paste(data.dir, \""+ext.replaceWithLinuxSafeCharacters(line[markerIndex], true)+".txt\", sep=\"\"), header=T)");
						writer.println("## library");
						writer.println("library(rmeta)");
						writer.println("### forest plot");
						writer.println("metaplot(input.data$Beta, input.data$SE, nn=(input.data$SE)^-2, input.data$Cohort,");
						writer.println("logeffect=F, boxsize=0.4,logticks=T, col=meta.colors(\"black\"),");
						writer.println("xlab=\""+line[markerIndex]+" OR\", ylab=\"Cohort\")");
						writer.close();						
					} catch (Exception e) {
						System.err.println("Error writing to " + filename);
						e.printStackTrace();
					}					
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir+datafile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir+datafile + "\"");
			System.exit(2);
		}		
	}
	
	public static void main(String[] args) {
		fromParameters("D:/mega/forests.crf", new Logger());
	}
}
