// "AB_lookup.dat" is necessary if the files do not contain {"Allele1 - AB"}/{"Allele2 - AB}
package cnv.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JOptionPane;

import cnv.filesys.ABLookup;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.filesys.Sample;
import cnv.filesys.SampleList;
import cnv.filesys.SourceFileHeaderData;
import common.Aliases;
import common.Array;
import common.CountHash;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Sort;
import common.ext;

public class SourceFileParser implements Runnable {

//    public static final String[][] SNP_HEADER_OPTIONS = {{ "SNP Name", "rsID", "Probe Set ID", "ProbeSet", "SNP", "ProbeSetName", "Name", "Marker"}}; // TODO any problems with mixing all idents together?
    public static final String[][] SNP_HEADER_OPTIONS = { Aliases.MARKER_NAMES };
//    public static final String[][] SNP_TABLE_FIELDS = { SNP_HEADER_OPTIONS[0], { "Chr", "Chromosome" }, { "Position" } };
    public static final String[][] SNP_TABLE_FIELDS = { SNP_HEADER_OPTIONS[0], Aliases.CHRS, Aliases.POSITIONS };
    public static final String[] DELIMITERS = {",", "\t", " "};
    public static final String[] DELIMITER_DESCRIPTIONS = {"COMMA", "TAB", "SPACE"};
    public static final String OVERWRITE_OPTION_FILE = ".overwrite_option";
    public static final String CANCEL_OPTION_FILE = ".cancel_option";
    public static final String HOLD_OPTION_FILE = ".hold_option";
    public static final String FILENAME_AS_ID_OPTION = "[FILENAME_ROOT]";
	
    private Project proj;
	private String[] files;
	private String[] markerNames;
	private int[] keysKeys;
	private long fingerprint;
	private char[][] abLookup;
	private Hashtable<String,String> fixes;
	private long timeBegan;
	private int threadId;
	private String delimiter;
	private boolean splitAB;

	public SourceFileParser(Project proj, String[] files, String[] markerNames, int[] keysKeys, char[][] abLookup, String delimiter, long fingerprint, Hashtable<String,String> fixes, long timeBegan) {
		this(proj, files, markerNames, keysKeys, abLookup, delimiter, fingerprint, fixes, timeBegan, -1);
	}

	public SourceFileParser(Project proj, String[] files, String[] markerNames, int[] keysKeys, char[][] abLookup, String delimiter, long fingerprint, Hashtable<String,String> fixes, long timeBegan, int threadId) {
		this.proj = proj;
		this.files = files;
		this.markerNames = markerNames;
		this.keysKeys = keysKeys;
		this.abLookup = abLookup;
		this.delimiter = delimiter; 
		this.fingerprint = fingerprint;
		this.fixes = fixes;
		this.timeBegan = timeBegan;
		this.threadId = threadId;
		this.splitAB = false;
		
	}

	public void run() {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int count, snpIndex, sampIndex, key;
//		int[] dataIndices, genotypeIndices;
		boolean parseAtAt;

		Sample samp;
		String sampleName, temp, filename, trav, idHeader;
		float[][] data;
		byte[][] genotypes;
		boolean ignoreAB;
		Hashtable<String, Float> allOutliers;
		HashMap<String, SourceFileHeaderData> headers;
		Logger log;
		
		log = proj.getLog();
		idHeader = proj.ID_HEADER.getValue();
//		delimiter = proj.getSourceFileDelimiter();
		allOutliers = new Hashtable<String, Float>();
		headers = proj.getSourceFileHeaders();
        try {
			for (int i = 0; i < files.length; i++) {
				if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+SourceFileParser.CANCEL_OPTION_FILE).exists()) {
					return;
				}
				try {
					log.report(ext.getTime()+"\t"+(i+1)+" of "+files.length + " -- " + files[i]);
					
					SourceFileHeaderData headerData = headers.get(files[i]);
					if (headerData == null) {
					    // TODO error, missing header for source files
					}
					
					reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[i]);
					temp = null;
					for (int k = 0; k < headerData.columnHeaderLineIndex + 1; k++) {
					    // iterate
					    temp = reader.readLine();
					}
//					do {
//						line = reader.readLine().trim().split(delimiter, -1);
//					} while (reader.ready() && (ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || (!idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line) == -1)));
					
//					DATA_0  {"GC Score", "GCscore", "confidence","Confidence"}
//					DATA_1  {"X Raw"}
//					DATA_2  {"Y Raw"}
//					DATA_3  {"X", "Xvalue", "Log Ratio", "intensity_1","Signal A"}
//					DATA_4  {"Y", "Yvalue", "Strength", "intensity_2","Signal B"}
//					DATA_5  {"Theta"}
//					DATA_6  {"R"}
//					DATA_7  {"B Allele Freq"},
//					DATA_8  {"Log R Ratio"}
//					
//					GENO_0  {"Allele1 - Forward", "Allele1", "genotype1", "Allele1 - Top","Forward Strand Base Calls","Forced Call","Forced Call Codes"}, 
//					GENO_1  {"Allele2 - Forward", "Forward Strand Base Calls", "genotype2", "Allele2 - Top","Allele B","Forced Call","Forced Call Codes"}, 
//					GENO_2  {"Allele1 - AB","Call Codes","Call"}, 
//					GENO_3  {"Allele2 - AB","Call Codes","Call"}}; 
					
					if (idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)) {
						sampIndex = -7;
					} else {
						sampIndex = headerData.colSampleIdent;
					}
					snpIndex = headerData.colSnpIdent;
					
					if (headerData.colX == -1 || headerData.colY == -1) {
						log.reportError("Error - File format not consistent! At the very least the files need to contain "+Array.toStr(Sample.DATA_FIELDS[3], "/")+" and "+Array.toStr(Sample.DATA_FIELDS[4], "/"));
						return;
					}
					if (headerData.colGeno1 == -1 || headerData.colGeno2 == -1) {
						log.reportError("Error - File format not consistent! The files need to contain "+Array.toStr(Sample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(Sample.GENOTYPE_FIELDS[1], "/"));
						return;
					}
					if ((headerData.colGenoAB1 == -1 || headerData.colGenoAB2 == -1) && abLookup == null) {
						ignoreAB = true;
					} else {
						ignoreAB = false;
					}

					sampleName = "";
					data = new float[][] {
					        headerData.colGC == -1 ? null : new float[markerNames.length],
			                headerData.colXRaw == -1 ? null : new float[markerNames.length],
	                        headerData.colYRaw == -1 ? null : new float[markerNames.length],
                            headerData.colX == -1 ? null : new float[markerNames.length],
                            headerData.colY == -1 ? null : new float[markerNames.length],
                            headerData.colTheta == -1 ? null : new float[markerNames.length],
                            headerData.colR == -1 ? null : new float[markerNames.length],
                            headerData.colBAF == -1 ? null : new float[markerNames.length],
                            headerData.colLRR == -1 ? null : new float[markerNames.length],
					};
					genotypes = new byte[][] {
					        Array.byteArray(markerNames.length, (byte) 0),
					        ignoreAB ? null : Array.byteArray(markerNames.length, (byte) -1),
					}; // two sets of genotypes, the "forward" alleles and the AB alleles
					
					count = 0;
					parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
					while (reader.ready()) {
						line = reader.readLine().split(delimiter, -1);
						if (idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)) {
							trav = files[i].substring(0, files[i].indexOf(proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)));
						} else {
							if (parseAtAt && line[sampIndex].indexOf("@") == -1) {
								log.reportError("Error - " + idHeader + " '" + line[sampIndex] + "' did not contain an @ sample");
								parseAtAt = false;
							}
							trav = parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex];
						}
						if (count == 0) {
							sampleName = trav;
						} else if (!trav.equals(sampleName)) {
							log.reportError("Found more than one ID in file "+files[i]+"(found "+trav+", expecting "+sampleName+")");
							return;
						} else if (count>markerNames.length) {
							log.reportError("Error - expecting only "+markerNames.length+" markers and found more than that in file "+files[i]);
							return;
						} else if (!markerNames[count].equals(line[snpIndex]) && !markerNames[count].startsWith("Blank")) {
							log.reportError("Found "+line[snpIndex]+" at marker #"+(count+1)+" in file "+files[i]+"; expecting "+markerNames[count]);
							return;
						}
						key = keysKeys[count];

						int ind = 0;
						int col = headerData.colGC;
						try {
						    data[ind][key] = col == -1 ? null : ext.isMissingValue(line[col]) ? Float.NaN : Float.parseFloat(line[col]);
                            ind = 1;
                            col = headerData.colXRaw;
						    data[ind][key] = col == -1 ? null : ext.isMissingValue(line[col]) ? Float.NaN : Float.parseFloat(line[col]);
						    ind = 2;
						    col = headerData.colYRaw;
                            data[ind][key] = col == -1 ? null : ext.isMissingValue(line[col]) ? Float.NaN : Float.parseFloat(line[col]);
                            ind = 3;
                            col = headerData.colX;
                            data[ind][key] = col == -1 ? null : ext.isMissingValue(line[col]) ? Float.NaN : Float.parseFloat(line[col]);
                            ind = 4;
                            col = headerData.colY;
                            data[ind][key] = col == -1 ? null : ext.isMissingValue(line[col]) ? Float.NaN : Float.parseFloat(line[col]);
                            ind = 5;
                            col = headerData.colTheta;
                            data[ind][key] = col == -1 ? null : ext.isMissingValue(line[col]) ? Float.NaN : Float.parseFloat(line[col]);
                            ind = 6;
                            col = headerData.colR;
                            data[ind][key] = col == -1 ? null : ext.isMissingValue(line[col]) ? Float.NaN : Float.parseFloat(line[col]);
                            ind = 7;
                            col = headerData.colBAF;
                            data[ind][key] = col == -1 ? null : ext.isMissingValue(line[col]) ? Float.NaN : Float.parseFloat(line[col]);
                            ind = 8;
                            col = headerData.colLRR;
                            data[ind][key] = col == -1 ? null : ext.isMissingValue(line[col]) ? Float.NaN : Float.parseFloat(line[col]);
                        } catch (NumberFormatException nfe) {
                            log.reportError("Error - failed at line " + key + " to parse '" + line[col] + "' into a valid " + Array.toStr(Sample.DATA_FIELDS[ind], "/"));
                            return;
                        } catch (Exception e) {
                            log.reportError("Some other exception");
                            log.reportException(e);
                            return;
                        }
                        
                        if (proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6 || proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6_CN/* || proj.ARRAY_TYPE.getValue() == ARRAY.DBGAP*/) {
                            data[0][key] = 1 - data[0][key];
                        }
                        if (proj.XY_SCALE_FACTOR.getValue() != 1) {
                            data[3][key] = data[3][key] / proj.XY_SCALE_FACTOR.getValue().floatValue(); 
                            data[4][key] = data[4][key] / proj.XY_SCALE_FACTOR.getValue().floatValue(); 
                        }
                        
						if (line[headerData.colGeno1].length() > 1) {
							if (i == 0 && !splitAB) {
								log.reportTimeInfo("Detected genotype calls for each allele are in the same column (Such as in Affymetrix .chp format) ");
								splitAB = true;
							}
							if (splitAB) {
								String tmp0 = line[headerData.colGeno1].substring(0, 1);
								String tmp1 = line[headerData.colGeno1].substring(1, 2);
								line[headerData.colGeno1] = tmp0;
								line[headerData.colGeno2] = tmp1;
								if (!ignoreAB) {
									try{
    									String tmp2 = line[headerData.colGenoAB1].substring(0, 1);
    									String tmp3 = line[headerData.colGenoAB1].substring(1, 2);
    									line[headerData.colGenoAB1] = tmp2;
    									line[headerData.colGenoAB2] = tmp3;
									} catch (Exception e) {
										log.reportTimeError("Could not parse genotypes on line "+Array.toStr(line));
										log.reportException(e);
										return;
									}
								}
							} else {
								log.reportTimeError("Inconsistant genotype call lengths");
							}
						}else if(splitAB){
							log.reportTimeError("Detected previously that genotype calls should be split, but the calls on line "+Array.toStr(line) +" did not");
							return;
						}
						genotypes[0][key] = (byte)ext.indexOfStr(line[headerData.colGeno1]+line[headerData.colGeno2], Sample.ALLELE_PAIRS);
						
						if (genotypes[0][key] == -1) {
							if (proj.getArrayType() == ARRAY.ILLUMINA) {// Affy matrix format does not use ALLELE_PAIRS

								if (ext.indexOfStr(line[headerData.colGeno1] + line[headerData.colGeno2], Sample.ALT_NULLS) == -1) {
									log.reportError("Error - failed to lookup " + line[headerData.colGeno1] + line[headerData.colGeno2] + " for marker " + markerNames[count] + " of sample " + files[i] + "; setting to missing");
								}
								genotypes[0][key] = 0;
							} else {
								genotypes[0][key] = (byte) ext.indexOfStr(line[headerData.colGenoAB1] + line[headerData.colGenoAB2], Sample.AB_PAIRS);
								if (genotypes[0][key] == -1) {
									log.reportError("Error - failed to lookup " + line[headerData.colGeno1] + line[headerData.colGeno2] + " for marker " + markerNames[count] + " of sample " + files[i] + "; setting to missing");
									genotypes[0][key] = 0;
								}
							}
						}
						if (ignoreAB) {
							// do nothing, will need to use these files to determine AB lookup table
						} else if (abLookup == null) {
							genotypes[1][key] = (byte)ext.indexOfStr(line[headerData.colGenoAB1]+line[headerData.colGenoAB2], Sample.AB_PAIRS);
						} else {
							if (genotypes[0][key] == 0) {
								genotypes[1][key] = -1;
							} else {
								genotypes[1][key] = 0;
                                if (line[headerData.colGeno1].charAt(0) == abLookup[count][1]) {
                                    genotypes[1][key]++;
                                } else if (line[headerData.colGeno1].charAt(0) != abLookup[count][0]) {
                                    log.reportError("Error - alleles for individual '" + (sampIndex < 0 ? trav : line[sampIndex]) + "' (" + line[headerData.colGeno1] + "/" + line[headerData.colGeno2] + ") do not match up with the defined AB lookup alleles (" + abLookup[count][0] + "/" + abLookup[count][1] + ") for marker " + markerNames[count]);
                                }
                                if (line[headerData.colGeno2].charAt(0) == abLookup[count][1]) {
                                    genotypes[1][key]++;
                                } else if (line[headerData.colGeno2].charAt(0) != abLookup[count][0]) {
                                    log.reportError("Error - alleles for individual '" + (sampIndex < 0 ? trav : line[sampIndex]) + "' (" + line[headerData.colGeno1] + "/" + line[headerData.colGeno2] + ") do not match up with the defined AB lookup alleles (" + abLookup[count][0] + "/" + abLookup[count][1] + ") for marker " + markerNames[count]);
                                }
							}
						}

						count++;
					}
					reader.close();
					if (count != markerNames.length) {
						log.reportError("Error - expecting "+markerNames.length+" markers and only found "+count+" in file "+files[i]);
						return;
					}

					if (fixes.containsKey(sampleName)) {
						sampleName = fixes.get(sampleName);
					}
					
					writer = null;
					trav = ext.replaceWithLinuxSafeCharacters(sampleName, true);
					if (!trav.equals(sampleName)) {
						if (writer == null) {
							try {
								writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+"FYI_IDS_WERE_CHANGED"+threadId+".txt", true));
								if (new File(proj.PROJECT_DIRECTORY.getValue()+"FYI_IDS_WERE_CHANGED.txt").length() == 0) {
									writer.println("The following IDs were changed so that spaces are removed and so that they could be used as valid filenames:");
								}
							} catch (Exception e) {
								log.reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue()+"FYI_IDS_WERE_CHANGED"+threadId+".txt");
								log.reportException(e);
							}
						}
						writer.println(sampleName+"\t"+trav);
						writer.close();
						sampleName = trav;
					}

					
					filename = SourceFileParser.determineFilename(proj.SAMPLE_DIRECTORY.getValue(true, true), sampleName, timeBegan, log);
					if (filename == null) {
						return;
					}

					samp = new Sample(sampleName, fingerprint, data, genotypes, false);
					samp.saveToRandomAccessFile(filename, allOutliers, sampleName); //TODO sampleIndex
				} catch (FileNotFoundException fnfe) {
					log.reportError("Error: file \""+files[i]+"\" not found in current directory");
					return;
				} catch (IOException ioe) {
					log.reportError("Error reading file \""+files[i]+"\"");
					return;
				}
			}

			if (allOutliers.size()>0) {
				if (threadId >= 0) {
					if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + threadId + ".ser").exists()) {
						log.reportError("Error - the following file already exists: " + proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + threadId + ".ser");
						return;
					} else {
						Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + threadId + ".ser");
					}
				} else {
					if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers0.ser").exists()) {
						log.reportError("Error - the following file already exists: " + proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers0.ser");
						return;
					} else {
						Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers0.ser");
					}
				}
			}

			log.report(ext.getTime()+"\tfinished");
		} catch (Exception e) {
			log.reportException(e);
		}
        
        // request release of resources
        log.report(ext.getTime()+"\tCleaning up...");
        System.gc();
	}

	public static void writeToLookupFile(Project proj, Hashtable<String,String[]> hash, int fileNumber) {
        PrintWriter writer;
        String filename;
        String[] keys;
        Logger log;
    
        log = proj.getLog();
        filename = proj.PROJECT_DIRECTORY.getValue()+"alleleLookup"+(fileNumber>0?"_atFile"+fileNumber:"")+".xln";
        log.report("Writing to file {" + filename + "}...", false, true);
        try {
            writer = new PrintWriter(new FileWriter(filename));
            keys = HashVec.getKeys(hash, false, false);
            writer.println("SNP\t"+Array.toStr(Sample.ALL_STANDARD_GENOTYPE_FIELDS));
            for (int k = 0; k < keys.length; k++) {
                writer.println(keys[k]+"\t"+Array.toStr(hash.get(keys[k])));
            }
            writer.close();
        } catch (Exception e) {
            log.reportError("Error writing to " + filename);
            log.reportException(e);
        }
        log.report("done.");
        
    }

    public static char[][] getABLookup(boolean abLookupRequired, String[] markerNames, Project proj) {
        ABLookup abLookup; 
        char[][] lookup;
        Logger log;
        
        log = proj.getLog();
        if (abLookupRequired && Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
            abLookup = new ABLookup(markerNames, proj.AB_LOOKUP_FILENAME.getValue(), true, false, proj.getLog());
            lookup = abLookup.getLookup();
            if (lookup == null) {
                log.report("Warning - failed to provide columns \"" + 
                                Sample.GENOTYPE_FIELDS[2][0] + "\" / \"" +
                                Sample.GENOTYPE_FIELDS[3][0] + "\" and the specificed AB_lookup file '" + 
                                proj.getProperty(proj.AB_LOOKUP_FILENAME) + 
                                "' does not exist; you'll need reconstruct the B allele for analysis");
            } else {
                abLookup.writeToFile(proj.PROJECT_DIRECTORY.getValue()+"checkAB.xln", proj.getLog());
            }
        } else {
            lookup = null;
        }
        
        return lookup;
    }

    public static String determineFilename(String dir, String sampleName, long timeBegan, Logger log) {
        String filename, trav;
        int version, versionToOverwrite;
        String[] overwriteOptions;
        int response;
        
        version = 0;
        versionToOverwrite = -1;
        do {
            trav = sampleName + (version == 0 ? "" : "." + version);
            version++;
            filename = dir + trav + Sample.SAMPLE_DATA_FILE_EXTENSION;
            if (new File(filename).exists() && new File(filename).lastModified() < timeBegan && versionToOverwrite == -1) {
                versionToOverwrite = version - 1;
            }
        } while (new File(filename).exists());
    
    
        filename = sampleName + (versionToOverwrite == 0 ? "" : "." + versionToOverwrite) + Sample.SAMPLE_DATA_FILE_EXTENSION;
        overwriteOptions = new String[] {
                "Rename new file " + trav + Sample.SAMPLE_DATA_FILE_EXTENSION, 
                "Overwrite existing file " + filename, 
                "Overwrite this and all future files", 
                "Cancel parser"
        };
        
        if (versionToOverwrite != -1) {
            while (new File(dir + SourceFileParser.HOLD_OPTION_FILE).exists()) {
                try {
                    Thread.sleep(500);
                } catch (InterruptedException ie) {}
            }
            if (new File(dir + SourceFileParser.CANCEL_OPTION_FILE).exists()) {
                return null;
            }
            
            Files.write("", dir + SourceFileParser.HOLD_OPTION_FILE);
            if (new File(dir + SourceFileParser.OVERWRITE_OPTION_FILE).exists()) {
                response = 1;
            } else {
                
                if (!System.getProperty("os.name").startsWith("Windows")) {
                    log.reportError("Error - the same sample name '"+sampleName+"' is being parsed again and the previous file existed before the current command began. "
                            + "This happens if you inadvertently restarted the parser or if the parser was interrupted and manually restarted. "
                            + "If you would like to start from scratch, the safest thing would be to cancel now and delete all files in the sample directory. "
                            + "To automatically overwrite when there is conflict like this, then create an empty file with the following name: "
                            + dir + SourceFileParser.OVERWRITE_OPTION_FILE
                    );
                    return null;
                }
                
                do {
                    response = JOptionPane.showOptionDialog(null, 
                            "Error - the same sample name '"+sampleName+"' is being parsed again and the previous file existed before the current command began.\n"+
                            "This happens if you inadvertently restarted the parser or if the parser was interrupted and manually restarted.\n"+
                            "If you would like to start from scratch, the safest thing would be to cancel now and delete all files in the sample directory.\n"+
                            "What would you like to do?"
                        , "What to do?", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null, overwriteOptions, overwriteOptions[0]);
                } while (response == -1);
            }
            new File(dir + SourceFileParser.HOLD_OPTION_FILE).delete();
            switch (response) {
            case 0:
                return dir + trav + Sample.SAMPLE_DATA_FILE_EXTENSION;
            case 2:
                Files.write("", dir + SourceFileParser.OVERWRITE_OPTION_FILE);
            case 1:
                return dir + filename;
            case 3:
                Files.write("", dir + SourceFileParser.CANCEL_OPTION_FILE);
                return null;
            default:
                JOptionPane.showMessageDialog(null, "Should be impossible to obtain this message ("+response+")", "Error", JOptionPane.ERROR_MESSAGE);
                break;
            }
        }
    
        return dir+trav+Sample.SAMPLE_DATA_FILE_EXTENSION;
    }

    public static void parseAlleleLookupFromFinalReports(Project proj) {
        BufferedReader reader;
        String[] line;
        Hashtable<String, String[]> hash;
        String[] files;
        int[] indices;
        int snpIndex;
        String delimiter, idHeader;
        String[] alleles;
        int expIndex;
        Logger log;
    
        log = proj.getLog();
        files = Files.list(proj.SOURCE_DIRECTORY.getValue(false, true), proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);
        if (files.length == 0) {
            log.reportError("Error - no files to parse");
            return;
        }
        log.report("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" to parse");
    
        idHeader = proj.getProperty(proj.ID_HEADER);
        delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
        hash = new Hashtable<String, String[]>();
        for (int i = 0; i<files.length; i++) {
            if (new File("report").exists()) {
                SourceFileParser.writeToLookupFile(proj, hash, i+1);
                new File("report").delete();
            }
            try {
                log.report(ext.getTime()+"\t"+(i+1)+" of "+files.length);
                reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[i]);
                do {
                    line = reader.readLine().trim().split(delimiter, -1);
                } while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || (!idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line)==-1)));
    
                snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
                indices = ext.indexFactors(Sample.ALL_STANDARD_GENOTYPE_FIELDS, line, false, proj.getLog(), false, false);
                
                while (reader.ready()) {
                    line = reader.readLine().split(delimiter);
                    if (hash.containsKey(line[snpIndex])) {
                        alleles = hash.get(line[snpIndex]);
                    } else {
                        if (i!=0) {
                            log.reportError("Error - snp '"+line[snpIndex]+"' first seen in file #"+i+" ("+files[i]+") and not earlier");
                        }
                        hash.put(line[snpIndex], alleles = new String[Sample.ALL_STANDARD_GENOTYPE_FIELDS.length]);
                    }
                    for (int j = 0; j < 2; j++) {
                        if (ext.indexOfStr(line[indices[j]], Sample.ALT_NULL) == -1) {
                            if (alleles[0] == null) {
                                alleles[0] = line[indices[j]];
                                expIndex = -1;
                            } else if (line[indices[j]].equals(alleles[0])) {
                                expIndex = 0;
                            } else if (alleles[1] == null) {
                                alleles[1] = line[indices[j]];
                                expIndex = -2;
                            } else if (line[indices[j]].equals(alleles[1])) {
                                expIndex = 1;
                            } else {
                                log.reportError("Error - snp '"+line[snpIndex]+"' has a new allele in file #"+i+" ("+files[i]+"): "+line[indices[j]]+" (previously "+alleles[0]+"/"+alleles[1]+")");
                                expIndex = -9;
                            }
                            
                            for (int k = 1; k < alleles.length/2; k++) {
                                switch (expIndex) {
                                case -1:
                                    if (alleles[k*2+0] != null) {
                                        log.reportError("Error - how can this be? -1");
                                    }
                                    alleles[k*2+0] = line[indices[k*2+j]];
                                    break;
                                case 0:
                                    if (!line[indices[k*2+j]].equals(alleles[k*2+0])) {
                                        log.reportError("Error - how can this be? 0");
                                    }                                   
                                    break;
                                case -2:
                                    if (alleles[k*2+1] != null) {
                                        log.reportError("Error - how can this be? -2");
                                    }
                                    alleles[k*2+1] = line[indices[k*2+j]];
                                    break;
                                case 1:
                                    if (!line[indices[k*2+j]].equals(alleles[k*2+1])) {
                                        log.reportError("Error - how can this be? 1");
                                    }                                   
                                    break;
                                case -9:
                                    break;
                                }
                            }                           
                        }
                    }
                }
            } catch (FileNotFoundException fnfe) {
                log.reportError("Error: file \""+files[i]+"\" not found in current directory");
                return;
            } catch (IOException ioe) {
                log.reportError("Error reading file \""+files[i]+"\"");
                return;
            }
        }
        SourceFileParser.writeToLookupFile(proj, hash, 0);
    }

    public static byte[] parseGenotypes(String[] line, int[] genotypeIndices, boolean ignoreAB, char[][] abLookup, int count, String sampleName, String markerName, String filename) {
        String genotype;
        byte genoForward, genoAB;
    
        if (genotypeIndices[4] >= 0) {
            genotype = line[genotypeIndices[4]];
        } else {
            genotype = line[genotypeIndices[0]] + line[genotypeIndices[1]];
        }
        genoForward = (byte) ext.indexOfStr(genotype, Sample.ALLELE_PAIRS);
        if (genoForward == -1) {
            if (ext.indexOfStr(genotype, Sample.ALT_NULLS) == -1) {
                // System.err.println("Error - failed to lookup "+genotype+" for marker "+markerName+" of sample "+filename);
                genoForward = -9;
            } else {
                genoForward = 0;
            }
        }
        if (ignoreAB) {
            genoAB = -1;
            // do nothing, will need to use these files to determine AB lookup table
        } else if (abLookup == null) {
            if (genotypeIndices[5] >= 0) {
                genotype = line[genotypeIndices[5]];
            } else {
                genotype = line[genotypeIndices[2]] + line[genotypeIndices[3]];
            }
            genoAB = (byte) ext.indexOfStr(genotype, Sample.AB_PAIRS);
            if (genoAB == -1) {
                if (ext.indexOfStr(genotype, Sample.ALT_NULLS) == -1) {
                    System.err.println("Error - failed to lookup " + genotype + " for marker " + markerName + " of sample " + filename);
                }
            }
        } else {
            if (genoForward == 0) {
                genoAB = -1;
            } else {
                genoAB = 0;
                if (genotypeIndices[4] >= 0) {
                    genotype = line[genotypeIndices[4]];
                } else {
                    genotype = line[genotypeIndices[0]] + line[genotypeIndices[1]];
                }
                for (int j = 0; j < 2; j++) {
                    if (genotype.charAt(j) == abLookup[count][1]) {
                        genoAB++;
                    } else if (genotype.charAt(j) != abLookup[count][0]) {
                        System.err.println("Error - alleles for individual '" + sampleName + "' (" + genotype + ") do not match up with the defined AB lookup alleles (" + abLookup[count][0] + "/" + abLookup[count][1] + ") for marker " + markerName);
                    }
                }
            }
        }
    
        return new byte[] { genoForward == -9 ? 1 : genoForward, genoAB };
    }

    public static void mapFilenamesToSamples(Project proj, String filename) {
                BufferedReader reader;
                PrintWriter writer;
                String[] line;
                int sampIndex;
                String[] files;
                String idHeader, delimiter;
                boolean longFormat, done;
                String prev;
                Logger log;
    
                log = proj.getLog();
                if (!Files.exists(proj.SOURCE_DIRECTORY.getValue(false, false))) {
                    proj.message("Source directory does not exist; change SOURCE_DIRECTORY= to point to the proper files");
                    return;
                }
                
                delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
                longFormat = proj.LONG_FORMAT.getValue();
                idHeader = proj.getProperty(proj.ID_HEADER);
                log.report(ext.getTime()+"\tSearching for "+proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)+" files in: "+proj.SOURCE_DIRECTORY.getValue(false, true));
                files = Files.list(proj.SOURCE_DIRECTORY.getValue(false, true), proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);
                log.report("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" to parse");
    
                try {
                    writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+filename));
                    for (int i = 0; i<files.length; i++) {
                        try {
                            reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[i]);
                            do {
                                line = reader.readLine().trim().split(delimiter);
                            } while (reader.ready()&&(line.length<3 || ext.indexOfStr(idHeader, line)==-1));
    
                            if (!reader.ready()) {
                                log.reportError("Error - went through enitre file without finding a line containing the user-defined ID header: "+idHeader);
                                return;
                            }
                            sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];
    
    //                  ParseAffymetrix, ParseAffySNP6, and ParseDbgap ::>
    //                        line = reader.readLine().split(delimiter);
    //                        writer.println(files[i]+"\t"+line[sampIndex]+"\t"+(line[sampIndex].indexOf("@") >= 0?line[sampIndex].split("@")[0]:line[sampIndex]));
    //                        reader.close();
    //                        
    //                  ParseIllumina ::>
                            done = false;
                            prev = null;
                            while (!done) {
                                line = reader.readLine().split(delimiter);
                                if (!line[sampIndex].equals(prev)) {
                                    writer.println(files[i]+"\t"+line[sampIndex]+(line[sampIndex].indexOf("@") >= 0?"\t"+line[sampIndex].split("@")[0]:""));
                                }
                                if (!longFormat || !reader.ready()) {
                                    done = true;
                                }
                                prev = line[sampIndex];
                            }
                            reader.close();
                        } catch (FileNotFoundException fnfe) {
                            log.reportError("Error: file \""+files[i]+"\" not found in "+proj.SOURCE_DIRECTORY.getValue(false, true));
                            writer.close();
                            return;
                        } catch (IOException ioe) {
                            log.reportError("Error reading file \""+files[i]+"\"");
                            writer.close();
                            return;
                        }
                    }
                    log.report(ext.getTime());
                    writer.close();
                } catch (Exception e) {
                    log.reportException(e);
                }
            }

    @SuppressWarnings("unchecked")
	public static int createFiles(Project proj, int numThreads) {
		BufferedReader reader;
		String[] line, markerNames, files;
		Hashtable<String,Integer> markerNameHash;
		int snpIndex;
		int[] keys, keysKeys, indices;
		long fingerprint;
		Thread[] threads;
		Vector<Vector<String>> fileCabinet;
		String trav;
		int count;
		boolean abLookupRequired;
		char[][] lookup;
        Hashtable<String,String> fixes;
        String idHeader, delimiter;
        String temp;
        int[][] delimiterCounts;
        boolean complete;
        long timeBegan;
        boolean parseAtAt;
        int sampIndex;
        String sampleName;
		String[] overwriteOptions;
		int response;
		Hashtable<String, Float> allOutliers;
		Vector<String> v;
		Logger log;
		int foundIDon, foundSNPon;
		AffyProcess affyProcess = null;
		
		log = proj.getLog();
        timeBegan = new Date().getTime();
        new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+SourceFileParser.OVERWRITE_OPTION_FILE).delete();
        new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+SourceFileParser.HOLD_OPTION_FILE).delete();
        new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+SourceFileParser.CANCEL_OPTION_FILE).delete();

		if (!proj.SOURCE_DIRECTORY.getValue(false, true).equals("") && !new File(proj.SOURCE_DIRECTORY.getValue(false, true)).exists()) {
			log.reportError("Error - the Project source location is invalid: "+proj.SOURCE_DIRECTORY.getValue(false, true));
			return 0;
		}
        
		if (!new File(proj.MARKER_POSITION_FILENAME.getValue(false, false)).exists()) {
			log.reportError("Error - missing markerPositions: "+proj.MARKER_POSITION_FILENAME.getValue(false, false));
			return 0;
		}

		delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
		idHeader = proj.getProperty(proj.ID_HEADER);
		log.report(ext.getTime()+"\tSearching for "+proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)+" files in: "+proj.SOURCE_DIRECTORY.getValue(false, true));
		files = Files.list(proj.SOURCE_DIRECTORY.getValue(false, true), proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);
		
		log.report("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" with a "+proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)+" extension");
		for (int i = 0; i < files.length; i++) {
			if (files[i].startsWith("Sample_Map.csv") || files[i].startsWith("SNP_Map.csv")) {
				files = Array.removeFromArray(files, i);
				i--;
			}
		}
		
		if (files.length == 0) {
			log.reportError("Error - no files to parse; are you sure you have the right extension specified for your FinalReport files? It is currently set to \""+proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)+"\"");
			return 0;
		}
		
		abLookupRequired = false;
		log.report("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" to parse");
		fixes = new Hashtable<String,String>();
		if (new File(proj.PROJECT_DIRECTORY.getValue()+"fixes.dat").exists()) {
			log.report("Also found a 'fixes.dat' file in the project directory, which will be used to rename samples");
			fixes = HashVec.loadFileToHashString(proj.PROJECT_DIRECTORY.getValue()+"fixes.dat", false);
		} else {
			log.report("Did not find a 'fixes.dat' file; assuming there are no swapped samples to rename");
		}
		ARRAY array = proj.getArrayType();
		switch (array) {
    		case AFFY_GW6:
    			log.reportTimeWarning("Affymetrix confidence scores will be imported as (1-conf)");
    			break;
    		case AFFY_GW6_CN:
    			log.reportTimeWarning("Affymetrix confidence scores will be imported as (1-conf)");
    			log.reportTimeInfo("Initializing parser for array type " + array);
    			affyProcess = new AffyProcess(proj, Files.toFullPaths(files, proj.SOURCE_DIRECTORY.getValue(false, true)), delimiter, log);
    			affyProcess.matchCn();
    			affyProcess.combineFirst();
    			for (int i = 0; i < files.length; i++) {
    				files[i] = ext.removeDirectoryInfo(affyProcess.getCombinedOutputFiles()[i]);
    			}
    			proj.setSourceFileHeaders(SourceFileHeaderData.validate(ext.parseDirectoryOfFile(affyProcess.getCombinedOutputFiles()[0]), "." + proj.getArrayType() + ".tmp.gz", true, proj.getLog()));
    			break;
    		case ILLUMINA:
    			break;
//    		case NGS:
//    		    break;
    		default:
    			log.reportTimeError("Invalid array type " + array);
    			break;
		}

		try {
			reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]);
			log.report(ext.getTime() + "]\tFound appropriate reader for: "+proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]);
			count = 0;
			foundSNPon = -1;
			foundIDon = -2;
			log.report("Searching for header line...");
			do {
				line = reader.readLine().trim().replace("#Column header: ", "").split(delimiter, -1);
				count++;
				if (count < 20) {
					log.report(Array.toStr(line), true, true, 11);
				}
				if (ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true, false, false)[0] != -1) {
					foundSNPon = count;
				}
				if (idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION) ||ext.indexOfStr(idHeader, line) != -1) {
					foundIDon = count;
				} 
				
			} while (reader.ready() && count < 1000 && foundIDon != foundSNPon);
			
			// If we reached the end of the file, it means that we didn't find the header we are looking for
			// The most common cause of this is that the delimiter was misspecified
			// The following code checks all of the common delimiters (tab, comma, space) and determines which one to use when it tries for a second time
			if (!reader.ready() || count == 1000) {
				if (foundSNPon == -1) {
					log.reportError("Could not find a header with the following tokens: "+Array.toStr(SourceFileParser.SNP_HEADER_OPTIONS[0], " / "));
				}
				if (foundIDon == -2) {
					log.reportError("Could not find a header with the selected id type: "+idHeader);
				}
				
				log.reportError("   Perhaps the delimiter, which is currently set to \""+proj.getProperty(proj.SOURCE_FILE_DELIMITER).getDelimiter() +"\", is incorrect? This can be corrected in the file "+proj.getPropertyFilename()+". In the meantime, the most stable delimiter will be determined for you...");
				log.reportError("   OR perhaps the ID_HEADER property is incorrect; the text '" + proj.getProperty(proj.ID_HEADER) + "' should be present in the header line.");
				log.reportError("   OR perhaps the ARRAY_TYPE property is incorrect;  options are " + Array.toStr(Project.ARRAY.class, ","));

				reader.close();
				reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]);
				delimiterCounts = new int[SourceFileParser.DELIMITERS.length][count];
				for (int i = 0; i < count; i++) {
					temp = reader.readLine();
					for (int j = 0; j < SourceFileParser.DELIMITERS.length; j++) {
						delimiterCounts[j][i] = ext.countInstancesOf(temp, SourceFileParser.DELIMITERS[j]);
					}
				}
				delimiter = null;
				for (int j = 0; j < SourceFileParser.DELIMITERS.length; j++) {
					if (Array.quantWithExtremeForTie(delimiterCounts[j], 0.5) > 4 && Array.quantWithExtremeForTie(delimiterCounts[j], 0.9) - Array.quantWithExtremeForTie(delimiterCounts[j], 0.1) == 0) {
						if (delimiter == null) {
							delimiter = SourceFileParser.DELIMITERS[j];
						} else {
							proj.message("Could not auto-detect the delimiter used in the Final Reports file: could be '"+delimiter+"' or '"+SourceFileParser.DELIMITERS[j]+"'");
							return 0;
						}
					}
				}
				reader.close();
				
				if (delimiter == null) {
					proj.message("Failed to auto-detect the delimiter used in the Final Reports file; exiting");
					return 0;
				}
				
				log.reportError("      - determined delimiter to be "+SourceFileParser.DELIMITER_DESCRIPTIONS[ext.indexOfStr(delimiter, SourceFileParser.DELIMITERS)]);

				// Tries again to determine the header fields and column names
				reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]);
				do {
					line = reader.readLine().trim().split(delimiter, -1);
				} while (reader.ready() && (ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || (!idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line) == -1)));
			}

			log.report(ext.getTime() + "]\tSearching for data fields...");
			// check immediately to make sure these fields are valid
			indices = ext.indexFactors(Sample.DATA_FIELDS, line, false, true, false, false); // dataIndices // TODO turn verbose back on after removing Raw X, Raw Y, Theta and R
			if (indices[3] == -1 || indices[4] == -1) {
				log.reportError("Error - at the very least the files need to contain "+Array.toStr(Sample.DATA_FIELDS[3], "/")+" and "+Array.toStr(Sample.DATA_FIELDS[4], "/"));
				log.reportError("      - failed to see that in "+files[0]);
				log.reportError(Array.toStr(line));
				return 0;
			}
			
			// TODO check different fields depending upon Affy/Illumina flag
			log.report(ext.getTime() + "]\tSearching for other fields...");
			indices = ext.indexFactors(Sample.GENOTYPE_FIELDS, line, false, true, true, false); // genotypeIndices
			if (indices[0] == -1 || indices[1] == -1) {
				log.reportError("Error - the files need to contain "+Array.toStr(Sample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(Sample.GENOTYPE_FIELDS[1], "/"));
				return 0;
			}
			if (indices[2] == -1 || indices[3] == -1) {
				abLookupRequired = true;
			}

			if (!idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)) {
				ext.indexFactors(new String[] { idHeader }, line, false, true); // sampIndex
			}
			snpIndex = ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true, true, true)[0];

			idHeader = proj.ID_HEADER.getValue();
			reader.mark(1000);

			if (idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)) {
				sampleName = files[0].substring(0, files[0].indexOf(proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)));
				sampIndex = -7;
				parseAtAt = false;
			} else {
				sampIndex = ext.indexFactors(new String[] { idHeader }, line, false, true)[0];
				line = reader.readLine().split(delimiter);
				parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
				if (parseAtAt && line[sampIndex].indexOf("@") == -1) {
					log.reportError("Error - " + idHeader + " '" + line[sampIndex] + "' did not contain an @ sample; if your ID's do not naturally contain at symbols, then set " + proj.PARSE_AT_AT_SYMBOL + " in the project properties file to FALSE");
					parseAtAt = false;
				}
				sampleName = parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex];
			}
			if (new File(proj.MARKER_DATA_DIRECTORY.getValue(false, false)+"markers.0.mdRAF").exists()) {
				
				if (!System.getProperty("os.name").startsWith("Windows")) {
					log.reportError("Error - Marker data (at least the first file 'markers.0.mdRAF') have already been parsed. This happens if you had previously transposed the data or if the parser was interrupted and manually restarted. Please delete the marker directory to start over from scratch.");
					return 0;
				}

				overwriteOptions = new String[] {
						"Delete All", 
						"Cancel parser"
				};
				
				response = JOptionPane.showOptionDialog(null, 
						"Marker data files (at least the first file 'markers.0.mdRAF') have already been parsed.\n"+
						"This happens if you had previously transposed the data or if the parser was interrupted and manually restarted.\n"+
						"If you would like to start from scratch, select \"Delete All\" earlier files.\n"+
						"Otherwise, press Cancel.\n"+
						"What would you like to do?"
					, "Marker data exists", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null, overwriteOptions, overwriteOptions[1]);

				switch (response) {
				case -1:
					return 0;
				case 0:
					deleteAllFilesInMarkerDataDirectory(proj);
					break;
				case 1:
					return 0;
				default:
					proj.message("Should be impossible to obtain this message ("+response+")");
					return 0;
				}
			}

			if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION).exists() || new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + ext.replaceWithLinuxSafeCharacters(sampleName, true) + Sample.SAMPLE_DATA_FILE_EXTENSION).exists()) {

				overwriteOptions = new String[] {
						"Delete All", 
						"Customize", 
						"Cancel parser"
				};
				
				response = JOptionPane.showOptionDialog(null, 
						"Sample data files (at least the first sample '"+sampleName+"') have already been parsed.\n"+
						"This happens if you inadvertently restarted the parser or if the parser was interrupted and manually restarted.\n"+
						"If you would like to start from scratch, select \"Delete All\" earlier files.\n"+
						"Otherwise, Cancel or you can \"Customize\" and determine what to do on a sample-by-sample basis.\n"+
						"What would you like to do?"
					, "Samples already exist", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null, overwriteOptions, overwriteOptions[2]);

				switch (response) {
				case -1:
					return 0;
				case 0:
					deleteAllFilesInSampleDirectory(proj);
					break;
				case 1:
				    // TODO complete 'customize' option for sample parsing
					// keep "outlier.ser"
					break;
				case 2:
					return 0;
				default:
					proj.message("Should be impossible to obtain this message ("+response+")");
					return 0;
				}
			}
			
			log.report(ext.getTime() + "]\tCleaning up before continuing...");
			TransposeData.deleteOlderRafs(proj.SAMPLE_DIRECTORY.getValue(true, true), new String[] {"outliers"}, new String[] {".ser"}, true, new String[] {"outliers.ser"});

			reader.reset();	
			
			int lines = Files.countLines(proj.SOURCE_DIRECTORY.getValue(false, true)+files[0], count);
			markerNameHash = new Hashtable<String, Integer>((int) (lines + (lines / 3)) + 100); // calc to never re-balance
			log.report(ext.getTime() + "]\tReading " + lines + " marker names...");
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter);
				trav = line[snpIndex];
				if (trav.equals("") || trav.equals("0")) {
					trav = "Blank"+count;
					count++;
				}
				if (markerNameHash.containsKey(trav)) {
					
					if (!proj.getProperty(proj.LONG_FORMAT)) {
						log.reportError("The same marker ('" + trav + "') was seen twice...");
					}
					
					if (idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)) {
						if (!proj.getProperty(proj.LONG_FORMAT)){
							log.report("... this could mean that the file contains multiple samples. This should not happen when " + proj.ID_HEADER.getName() + " is set to " + SourceFileParser.FILENAME_AS_ID_OPTION);
						}else{
							log.reportError("Error - " + proj.ID_HEADER.getName() + " was set to " + SourceFileParser.FILENAME_AS_ID_OPTION + ", which is invalid for Long Format files");
						}
						return 0;
					} else if (!(parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex]).equals(sampleName)) {
						if (!proj.getProperty(proj.LONG_FORMAT)) {
							log.reportError("... and the sample name changed from " + sampleName + " to " + (parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex]) + ", so this must be a Long Format file. The property file will be updated to reflect this, and an attempt will be made to launch the Long Format file processor now.");
							proj.setProperty(proj.LONG_FORMAT, Boolean.TRUE);
							proj.saveProperties();
						}
						reader.close();
						//we should have all markers now...
						markerNames = Array.toStringArray(markerNameHash);
						keys = Markers.orderMarkers(markerNames, proj.MARKER_POSITION_FILENAME.getValue(), proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());
						if (keys == null) {
							return checkForSNP_Map(proj, log);
						} else {
							return createFilesFromLongFormat(proj, files, idHeader, fixes, delimiter, abLookupRequired, timeBegan);
						}
					} else {
						log.reportError("Error - detected the same marker ('" + trav + "') twice, but the sample name did not change. Could not determine the file format to parse");
						return 0;
					}
				} else {
					markerNameHash.put(trav, markerNameHash.size());
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \""+proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]+"\" not found in current directory");
			return 0;
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]+"\"");
			log.reportException(ioe);
			return 0;
		}
		
		markerNames = Array.toStringArray(markerNameHash);
		keys = Markers.orderMarkers(markerNames, proj.MARKER_POSITION_FILENAME.getValue(), proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());
		if (keys == null) {
			return checkForSNP_Map(proj, log);
		}
		keysKeys = Sort.quicksort(keys); // very important
		fingerprint = proj.getMarkerSet().getFingerprint();
		log.report("There are "+markerNames.length+" markers being processed (fingerprint: "+fingerprint+")");
		if (proj.XY_SCALE_FACTOR.getValue() != 1) {
			log.report("The XY_SCALE_FACTOR was set to " + proj.XY_SCALE_FACTOR.getValue() + ", scaling X and Y values...");
		}
		lookup = SourceFileParser.getABLookup(abLookupRequired, markerNames, proj);

		if (affyProcess != null) {
			affyProcess.combineAll(numThreads);
		}
		fileCabinet = new Vector<Vector<String>>();
		for (int i = 0; i<numThreads; i++) {
			fileCabinet.add(new Vector<String>());
		}
		for (int i = 0; i<files.length; i++) {
			fileCabinet.elementAt(i%numThreads).add(files[i]);
		}
		threads = new Thread[numThreads];
		for (int i = 0; i<numThreads; i++) {
			threads[i] = new Thread(new SourceFileParser(proj, fileCabinet.elementAt(i).toArray(new String[fileCabinet.elementAt(i).size()]), markerNames, keysKeys, lookup, delimiter, fingerprint, fixes, timeBegan, i));
			threads[i].start();
			try {
				Thread.sleep(100L);
			} catch (InterruptedException ex) {}
		}
		
		complete = false;
		while (!complete) {
			complete = true;
			for (int i = 0; i<numThreads; i++) {
				if (threads[i].isAlive()) {
					complete = false;
				}
			}
			if (!complete) {
				try {
					Thread.sleep(1000L);
				} catch (InterruptedException ex) {}
			}
		}
		
		log.report(ext.getTime() + "]\tWriting sample list...");
		SampleList.generateSampleList(proj).writeToTextFile(proj.PROJECT_DIRECTORY.getValue()+"ListOfSamples.txt");

		allOutliers = new Hashtable<String, Float>();


		v = new Vector<String>();

        log.report(ext.getTime() + "]\tWriting outliers list(s)...");
		for (int i = 0; i<numThreads; i++) {
			if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser").exists()) {
				allOutliers.putAll((Hashtable<String, Float>) Files.readSerial(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser"));
				new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser").delete();
			}
			if (new File(proj.PROJECT_DIRECTORY.getValue()+"FYI_IDS_WERE_CHANGED"+i+".txt").exists()) {
				v.add(proj.PROJECT_DIRECTORY.getValue()+"FYI_IDS_WERE_CHANGED"+i+".txt");
			}
		}
		if (v.size() > 0) {
			Files.cat(Array.toStringArray(v), proj.PROJECT_DIRECTORY.getValue()+"FYI_IDS_WERE_CHANGED.txt", Array.intArray(v.size(), 1), log);
			for (int i = 0; i < v.size(); i++) {
				new File(v.elementAt(i)).delete();
			}
		}
		
		// changed 6-26-15 to always write an outliers.ser, even if no outliers exist
		Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
		log.report(ext.getTime() + "]\tSource File Parsing Complete.");
		if (abLookupRequired && !Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false))) {
			return 6;
		} else {
			return 1;
		}
	}

	private static int checkForSNP_Map(Project proj, Logger log) {
		String filename;
		filename = proj.getLocationOfSNP_Map(true);
		if (filename != null) {
			log.reportError("\nSuch a file was found: "+filename);
			log.reportError("\nIn order to process it use the command:");
			log.reportError("   java -cp [package_name].jar cnv.manage.Markers proj="+proj.getPropertyFilename()+" snps="+filename);
			return 7;
		}
		return 0;
	}

	public static void deleteAllFilesInSampleDirectory(Project proj) {
		String[] filesToDelete;
		filesToDelete = Files.list(proj.SAMPLE_DIRECTORY.getValue(false, true), Sample.SAMPLE_DATA_FILE_EXTENSION, false);
		for (int i = 0; i < filesToDelete.length; i++) {
			new File(proj.SAMPLE_DIRECTORY.getValue(false, true) + filesToDelete[i]).delete();
		}
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser").delete();
	}

	public static void deleteAllFilesInMarkerDataDirectory(Project proj) {
		String[] filesToDelete;
		
		filesToDelete = Files.list(proj.MARKER_DATA_DIRECTORY.getValue(false, true), MarkerData.MARKER_DATA_FILE_EXTENSION, false);
		for (int i = 0; i < filesToDelete.length; i++) {
			new File(proj.MARKER_DATA_DIRECTORY.getValue(false, true) + filesToDelete[i]).delete();
		}
		new File(proj.MARKER_DATA_DIRECTORY.getValue(true, true) + "outliers.ser").delete();
	}

	public static int createFilesFromLongFormat(Project proj, String[] files, String idHeader, Hashtable<String, String> fixes, String delimiter, boolean abLookupRequired, long timeBegan) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, markerNames, list;
		long fingerprint;
		MarkerSet markerSet;
		Hashtable<String, Integer> markerIndices;
		int count, markerCount;
		char[][] abLookup;
		String filename;
		Hashtable<String, Float> allOutliers;
		Hashtable<String, String> renamedIDsHash;
		Logger log;

		log = proj.getLog();
        log.report("Parsing files using the Long Format algorithm");
        
        // creates and serializes the markers.bim, the returned keys are not used here as the markerIndices fill that purpose
        markerSet = proj.getMarkerSet();
        markerNames = markerSet.getMarkerNames();
		fingerprint = proj.getMarkerSet().getFingerprint();
		
		abLookup = SourceFileParser.getABLookup(abLookupRequired, markerNames, proj);
		
		markerIndices = new Hashtable<String, Integer>();
		for (int i = 0; i < markerNames.length; i++) {
			markerIndices.put(markerNames[i], new Integer(i));
		}
		
		log.report("There were "+markerNames.length+" markers present in '"+proj.MARKERSET_FILENAME.getValue(true, true)+"' that will be processed from the source files (fingerprint: "+fingerprint+")");
		
		int snpIndex, sampIndex, key;
		String trav, temp;
		int[] dataIndices, genotypeIndices;
		boolean parseAtAt;

		Sample samp;
		String sampleName;
		float[][] data;
		byte[][] genotypes;
		boolean ignoreAB;
		CountHash countHash, dupHash;
		boolean done;
		int numCols;
		
		if (Files.exists(proj.PROJECT_DIRECTORY.getValue()+"FYI_IDS_WERE_CHANGED.txt")) {
			Files.backup("FYI_IDS_WERE_CHANGED.txt", proj.PROJECT_DIRECTORY.getValue(), proj.PROJECT_DIRECTORY.getValue(), true);
		}
		
		count = markerCount = 0;
		countHash = new CountHash();
		dupHash = new CountHash();
		allOutliers = new Hashtable<String, Float>();
		renamedIDsHash = new Hashtable<String, String>();
        try {
			for (int i = 0; i<files.length; i++) {
				try {
					log.report(ext.getTime()+"\t"+(i+1)+" of "+files.length+" ("+files[i]+")");
					reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[i]);
					do {
						line = reader.readLine().trim().replace("#Column header: ", "").split(delimiter, -1);
					} while (reader.ready()&&(ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || ext.indexOfStr(idHeader, line)==-1));

					log.report("Searching: "+Array.toStr(line));
					dataIndices = ext.indexFactors(Sample.DATA_FIELDS, line, false, true, false, false);
					genotypeIndices = ext.indexFactors(Sample.GENOTYPE_FIELDS, line, false, true, false, false);
					sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];
					snpIndex = ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
					
					if (dataIndices[3] == -1 || dataIndices[4] == -1) {
						log.reportError("Error - File format not consistent! At the very least the files need to contain "+Array.toStr(Sample.DATA_FIELDS[3], "/")+" and "+Array.toStr(Sample.DATA_FIELDS[4], "/"));
						return 0;
					}
					if (genotypeIndices[0] == -1 || genotypeIndices[1] == -1) {
						log.reportError("Error - File format not consistent! The files need to contain "+Array.toStr(Sample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(Sample.GENOTYPE_FIELDS[1], "/"));
						return 0;
					}
					if ((genotypeIndices[2] == -1 || genotypeIndices[3] == -1) && abLookup == null) {
						ignoreAB = true;
					} else {
						ignoreAB = false;
					}

					done = false;
					sampleName = "just starting";
					data = null;
					genotypes = null;
					parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
					numCols = -1;
					while (!done) {
						if (reader.ready()) {
							line = reader.readLine().split(delimiter);
							if (numCols == -1) {
								numCols = line.length;
							} else if (line.length != numCols) {
								log.reportError("Error - mismatched number of columns at marker index "+markerCount);
							}
							if (parseAtAt&&line[sampIndex].indexOf("@")==-1) {
								log.reportError("Error - "+idHeader+" '"+line[sampIndex]+"' did not contain an @ sample");
								parseAtAt = false;
							}
							temp = parseAtAt?line[sampIndex].substring(0, line[sampIndex].indexOf("@")):line[sampIndex];
							trav = ext.replaceWithLinuxSafeCharacters(temp, true);
							
							if (!trav.equals(temp) && !renamedIDsHash.containsKey(temp)) {
								try {
									writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+"FYI_IDS_WERE_CHANGED.txt", true));
									if (renamedIDsHash.size() == 0) {
										writer.println("The following IDs were changed so that spaces are removed and so that they could be used as valid filenames:");
									}
									writer.println(temp+"\t"+trav);
									writer.close();
								} catch (Exception e) {
									log.reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue()+"FYI_IDS_WERE_CHANGED.txt");
									log.reportException(e);
								}
								renamedIDsHash.put(temp, trav);
							}
							
						} else {
							done = true;
							trav = null;
						}

						if (done || !trav.equals(sampleName)) {
							if (!sampleName.equals("just starting")) {
								if (markerCount!=markerNames.length) {
									log.reportError("Error - expecting "+markerNames.length+" markers and only found "+markerCount+" for sample "+sampleName+"; this usually happens when the input file is truncated");
									return 0;
								}

								if (fixes.containsKey(sampleName)) {
									sampleName = fixes.get(sampleName);
								}
								
								filename = SourceFileParser.determineFilename(proj.SAMPLE_DIRECTORY.getValue(true, true), sampleName, timeBegan, log);
								if (filename == null) {
									return 0;
								}

//								log.report("Saving file sample "+sampleName+" which has data for "+markerCount+" markers"); 
								samp = new Sample(sampleName, fingerprint, data, genotypes, false);
								samp.saveToRandomAccessFile(filename, allOutliers, sampleName);
								
								count++;
								markerCount = 0;
								
								if (count % 100 == 0) {
									log.report(count+" of ???"); // TODO make this intelligent and guess the number remaining in this file based on average byte usage and the size of the file (estimate when the file is gzip compressed, the factor is probably somewhere around ~15%)
								}
							}
							if (!done) {
								if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + (fixes.containsKey(trav)?fixes.get(trav):trav) + Sample.SAMPLE_DATA_FILE_EXTENSION).exists()) {
									log.reportError("Warning - marker data must be out of order, becaue we're seeing "+trav+(fixes.containsKey(trav)?"-->"+fixes.get(trav):"")+" again at line "+count);
									samp = Sample.loadFromRandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(true, true) + (fixes.containsKey(trav)?fixes.get(trav):trav) + Sample.SAMPLE_DATA_FILE_EXTENSION, proj.JAR_STATUS.getValue());
									data = samp.getAllData();
									genotypes = samp.getAllGenotypes();
								} else {
									data = new float[Sample.DATA_FIELDS.length][];
									for (int j = 0; j<data.length; j++) {
										if (dataIndices[j] != -1) {
											data[j] = Array.floatArray(markerNames.length, Float.POSITIVE_INFINITY);
										}
				                    }
									genotypes = new byte[2][];
									genotypes[0] = Array.byteArray(markerNames.length, (byte)0);	// used to be initialized to Byte.MIN_VALUE when AB genotypes && abLookup were both absent
									if (!ignoreAB) {
										genotypes[1] = Array.byteArray(markerNames.length, (byte)-1);	// used to be initialized to Byte.MIN_VALUE when AB genotypes && abLookup were both absent
									}
								}
								sampleName = trav;
							}
						}

						if (!done) {
							countHash.add(line[snpIndex]);							
							if (markerIndices.containsKey(line[snpIndex])) {
								key = markerIndices.get(line[snpIndex]);

								for (int j = 0; j<Sample.DATA_FIELDS.length; j++) {
									try {
										if (dataIndices[j] != -1) {
											if (!(data[j][key]+"").equals("Infinity")) {
												dupHash.add(line[snpIndex]);
												log.reportError("Sample "+trav+" already has data for marker "+line[snpIndex]+" (Was the parsing restarted? Delete the old directories first)");
											}
											if (ext.isMissingValue(line[dataIndices[j]])) {
												data[j][key] = Float.NaN;
											} else {
												data[j][key] = Float.parseFloat(line[dataIndices[j]]);
											}
										}
									} catch (NumberFormatException nfe) {
										log.reportError("Error - failed to parse '"+line[dataIndices[j]]+"' into a valid "+Array.toStr(Sample.DATA_FIELDS[j], "/"));
										return 0;
									}
								}

								genotypes[0][key] = (byte)ext.indexOfStr(line[genotypeIndices[0]]+line[genotypeIndices[1]], Sample.ALLELE_PAIRS);
								if (genotypes[0][key] == -1) {
									if (ext.indexOfStr(line[genotypeIndices[0]]+line[genotypeIndices[1]], Sample.ALT_NULLS) == -1) {
										log.reportError("Error - failed to lookup "+line[genotypeIndices[0]]+line[genotypeIndices[1]]+" for marker "+line[snpIndex]+" of sample "+files[i]+"; setting to missing");
									}
									genotypes[0][key] = 0;
								}
								
//								if (!ignoreAB) {
//									genotypes[1][key] = (byte)ext.indexOfStr(line[genotypeIndices[2]]+line[genotypeIndices[3]], Sample.AB_PAIRS);
//								}
								if (ignoreAB) {
									// do nothing, will need to use these files to determine AB lookup table
								} else if (abLookup == null) {
									genotypes[1][key] = (byte)ext.indexOfStr(line[genotypeIndices[2]]+line[genotypeIndices[3]], Sample.AB_PAIRS);
								} else {
									if (genotypes[0][key] == 0) {
										genotypes[1][key] = -1;
									} else {
										genotypes[1][key] = 0;
										for (int j = 0; j<2; j++) {
											if (line[genotypeIndices[j]].charAt(0) == abLookup[key][1]) {
												genotypes[1][key]++;
											} else if (line[genotypeIndices[j]].charAt(0) != abLookup[key][0]) {
												log.reportError("Error - alleles for individual '"+line[sampIndex]+"' ("+line[genotypeIndices[0]]+"/"+line[genotypeIndices[1]]+") do not match up with the defined AB lookup alleles ("+abLookup[key][0]+"/"+abLookup[key][1]+") for marker "+markerNames[key]);
											}
		                                }
									}
								}
							}
							markerCount++;
						}
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					log.reportError("Error: file \""+files[i]+"\" not found in current directory");
					return 0;
				} catch (IOException ioe) {
					log.reportError("Error reading file \""+files[i]+"\"");
					return 0;
				}
			}

			if (allOutliers.size()>0) {
				if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser").exists()) {
					log.reportError("Error - the following file already exists: " + proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
					return 0;
				} else {
					Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
				}
			}

			log.report(ext.getTime()+"\tfinished");
		} catch (Exception e) {
			log.reportException(e);
		}

        
		log.report(ext.getTime()+"\t"+"Parsed "+count+" sample(s)");
		SampleList.generateSampleList(proj).writeToTextFile(proj.PROJECT_DIRECTORY.getValue()+"ListOfSamples.txt");

		try {
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+"ListOfMarkers.txt"));
			writer.println("Marker\tExpected\tTimesSeen\tTimesDuplicated");
			for (int j = 0; j < markerNames.length; j++) {
				writer.println(markerNames[j]+"\t1\t"+countHash.getCount(markerNames[j])+"\t"+dupHash.getCount(markerNames[j]));
				countHash.remove(markerNames[j], false);
			}
			list = countHash.getValues();
			for (int j = 0; j < list.length; j++) {
				writer.println(list[j]+"\t0\t"+countHash.getCount(markerNames[j])+"\t.");
			}			
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + "ListOfMarkers.txt");
			log.reportException(e);
		}

        new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+SourceFileParser.OVERWRITE_OPTION_FILE).delete();
        new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+SourceFileParser.HOLD_OPTION_FILE).delete();
        new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+SourceFileParser.CANCEL_OPTION_FILE).delete();
       
		if (abLookupRequired && !Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false))) {
			return 6;
		} else {
			return 1;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj;
		String filename = null;
		boolean map = false;
		int numThreads = 1;
		boolean parseAlleleLookupFromFinalReports = false;
		String mapOutput = "filenamesMappedToSamples.txt";

		String usage = "\n"+
		"cnv.manage.ParseIllumina requires 0-1 arguments\n"+
		"   (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"   (2) number of threads to use (i.e. threads="+numThreads+" (default))\n"+
		" OPTIONAL:\n"+
		"   (3) map filenames to sample IDs (i.e. -mapFiles ("+(map?"":"not the ")+"default))\n"+
		"   (4) output file for mappings (i.e. out="+mapOutput+" (default))\n"+
		" OR:\n"+
		"   (1) parse Forward/TOP/AB/etc lookup (i.e. --parseAlleleLookup (not the default))\n"+
		"";
		
		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				return;
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-mapFiles")) {
				map = true;
				numArgs--;
			} else if (args[i].startsWith("-parseAlleleLookup")) {
				parseAlleleLookupFromFinalReports = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			return;
		}

		try {
			proj = new Project(filename, false);
			if (map) {
			    SourceFileParser.mapFilenamesToSamples(proj, mapOutput);
			} else if (parseAlleleLookupFromFinalReports) {
			    SourceFileParser.parseAlleleLookupFromFinalReports(proj);
			} else {
				createFiles(proj, numThreads);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
