// "AB_lookup.dat" is necessary if the files do not contain {"Allele1 - AB"}/{"Allele2 - AB}
package cnv.manage;

import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;

import cnv.filesys.ABLookup;
import cnv.filesys.MarkerData;
import cnv.filesys.Sample;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.SampleList;
import common.*;

public class ParseIllumina implements Runnable {
	public static final String[][] SNP_HEADER_OPTIONS = {{"SNP Name", "rsID"}};
	public static final int EXP_NUM_MARKERS = 650000;
	public static final String[] FIELDS = {"Sample ID", "Sample Name"};
	public static final String[][] SNP_TABLE_FIELDS = {{"Name"}, {"Chr", "Chromosome"}, {"Position"}};
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

	public ParseIllumina(Project proj, String[] files, String[] markerNames, int[] keysKeys, char[][] abLookup, String delimiter, long fingerprint, Hashtable<String,String> fixes, long timeBegan) {
		this(proj, files, markerNames, keysKeys, abLookup, delimiter, fingerprint, fixes, timeBegan, -1);
	}

	public ParseIllumina(Project proj, String[] files, String[] markerNames, int[] keysKeys, char[][] abLookup, String delimiter, long fingerprint, Hashtable<String,String> fixes, long timeBegan, int threadId) {
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
	}

	public void run() {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int count, snpIndex, sampIndex, key;
		String trav;
		int[] dataIndices, genotypeIndices;
		boolean parseAtAt;

		Sample samp;
		String sampleName;
		float[][] data;
		byte[][] genotypes;
		boolean ignoreAB;
		String idHeader;
//		String delimiter;
		String filename;
		Hashtable<String, Float> allOutliers;
		Logger log;
		
		log = proj.getLog();
		idHeader = proj.getProperty(proj.ID_HEADER);
//		delimiter = proj.getSourceFileDelimiter();
		allOutliers = new Hashtable<String, Float>();
        try {
			for (int i = 0; i<files.length; i++) {
				if (new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)+CANCEL_OPTION_FILE).exists()) {
					return;
				}
				try {
					log.report(ext.getTime()+"\t"+(i+1)+" of "+files.length);
					reader = Files.getAppropriateReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[i]);
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready() && (ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || (!idHeader.equals(FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line) == -1)));
					
					dataIndices = ext.indexFactors(Sample.DATA_FIELDS, line, false, true, false, false);
					genotypeIndices = ext.indexFactors(Sample.GENOTYPE_FIELDS, line, false, true, false, false);
					if (idHeader.equals(FILENAME_AS_ID_OPTION)) {
						sampIndex = -7;
					} else {
						sampIndex = ext.indexFactors(new String[] { idHeader }, line, false, true)[0];
					}
					snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
					
					if (dataIndices[3] == -1 || dataIndices[4] == -1) {
						log.reportError("Error - File format not consistent! At the very least the files need to contain "+Array.toStr(Sample.DATA_FIELDS[3], "/")+" and "+Array.toStr(Sample.DATA_FIELDS[4], "/"));
						return;
					}
					if (genotypeIndices[0] == -1 || genotypeIndices[1] == -1) {
						log.reportError("Error - File format not consistent! The files need to contain "+Array.toStr(Sample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(Sample.GENOTYPE_FIELDS[1], "/"));
						return;
					}
					if ((genotypeIndices[2] == -1 || genotypeIndices[3] == -1) && abLookup == null) {
						ignoreAB = true;
					} else {
						ignoreAB = false;
					}

					sampleName = "";
					data = new float[Sample.DATA_FIELDS.length][];
					for (int j = 0; j<data.length; j++) {
						if (dataIndices[j] != -1) {
							data[j] = new float[markerNames.length];
						}
                    }
					genotypes = new byte[2][]; // two sets of genotypes, the "forward" alleles and the AB alleles
					
					// this is the for the forward alleles
					genotypes[0] = Array.byteArray(markerNames.length, (byte)0);	// used to be initialized to Byte.MIN_VALUE when AB genotypes && abLookup were both absent
					
					// this is for the AB alleles
					if (!ignoreAB) {
						genotypes[1] = Array.byteArray(markerNames.length, (byte)-1);	// used to be initialized to Byte.MIN_VALUE when AB genotypes && abLookup were both absent
					}
					
					count = 0;
//					parseAtAt = Boolean.parseBoolean(proj.getProperty(proj.PARSE_AT_AT_SYMBOL));
					parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
					while (reader.ready()) {
						line = reader.readLine().split(delimiter);
						if (idHeader.equals(FILENAME_AS_ID_OPTION)) {
							trav = files[i].substring(0, files[i].indexOf(proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)));
						} else {
							if (parseAtAt && line[sampIndex].indexOf("@") == -1) {
								log.reportError("Error - " + idHeader + " '" + line[sampIndex] + "' did not contain an @ sample");
								parseAtAt = false;
							}
							trav = parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex];
						}
						if (count==0) {
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
						for (int j = 0; j<Sample.DATA_FIELDS.length; j++) {
							try {
								if (dataIndices[j] != -1) {
									if (ext.isMissingValue(line[dataIndices[j]])) {
										data[j][key] = Float.NaN;
									} else {
										data[j][key] = Float.parseFloat(line[dataIndices[j]]);
									}
								}
							} catch (NumberFormatException nfe) {
								log.reportError("Error - failed to parse '"+line[dataIndices[j]]+"' into a valid "+Array.toStr(Sample.DATA_FIELDS[j], "/"));
								return;
							} catch (Exception e) {
								log.reportError("Some other exception");
								log.reportException(e);
								return;
							}
						}

						genotypes[0][key] = (byte)ext.indexOfStr(line[genotypeIndices[0]]+line[genotypeIndices[1]], Sample.ALLELE_PAIRS);
						if (genotypes[0][key] == -1) {
							if (ext.indexOfStr(line[genotypeIndices[0]]+line[genotypeIndices[1]], Sample.ALT_NULLS) == -1) {
								log.reportError("Error - failed to lookup "+line[genotypeIndices[0]]+line[genotypeIndices[1]]+" for marker "+markerNames[count]+" of sample "+files[i]+"; setting to missing");
							}								
							genotypes[0][key] = 0;
						}
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
									if (line[genotypeIndices[j]].charAt(0) == abLookup[count][1]) {
										genotypes[1][key]++;
									} else if (line[genotypeIndices[j]].charAt(0) != abLookup[count][0]) {
										log.reportError("Error - alleles for individual '" + (sampIndex < 0 ? trav : line[sampIndex]) + "' (" + line[genotypeIndices[0]] + "/" + line[genotypeIndices[1]] + ") do not match up with the defined AB lookup alleles (" + abLookup[count][0] + "/" + abLookup[count][1] + ") for marker " + markerNames[count]);
									}
                                }
							}
						}

						count++;
					}
					reader.close();
					if (count!=markerNames.length) {
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
								writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"FYI_IDS_WERE_CHANGED"+threadId+".txt", true));
								if (new File(proj.getProjectDir()+"FYI_IDS_WERE_CHANGED.txt").length() == 0) {
									writer.println("The following IDs were changed so that spaces are removed and so that they could be used as valid filenames:");
								}
							} catch (Exception e) {
								log.reportError("Error writing to " + proj.getProjectDir()+"FYI_IDS_WERE_CHANGED"+threadId+".txt");
								log.reportException(e);
							}
						}
						writer.println(sampleName+"\t"+trav);
						writer.close();
						sampleName = trav;
					}

					
					filename = determineFilename(proj.getDir(proj.SAMPLE_DIRECTORY, true), sampleName, timeBegan, log);
					if (filename == null) {
						return;
					}

					samp = new Sample(sampleName, fingerprint, data, genotypes, false);
//					samp.serialize(proj.getDir(proj.SAMPLE_DIRECTORY, true) + trav + Sample.SAMPLE_DATA_FILE_EXTENSION);
//					samp.saveToRandomAccessFile(filename);
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
					if (new File(proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers" + threadId + ".ser").exists()) {
						log.reportError("Error - the following file already exists: " + proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers" + threadId + ".ser");
						return;
					} else {
						Files.writeSerial(allOutliers, proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers" + threadId + ".ser");
					}
				} else {
					if (new File(proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers0.ser").exists()) {
						log.reportError("Error - the following file already exists: " + proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers0.ser");
						return;
					} else {
						Files.writeSerial(allOutliers, proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers0.ser");
					}
				}
			}

			log.report(ext.getTime()+"\tfinished");
		} catch (Exception e) {
			log.reportException(e);
		}
	}

	private static String determineFilename(String dir, String sampleName, long timeBegan, Logger log) {
		String filename, trav;
		int version, versionToOverwrite;
		String[] overwriteOptions;
		int response;
		
		version = 0;
		versionToOverwrite = -1;
		do {
			trav = sampleName+(version==0?"":"."+version);
			version++;
			filename = dir+trav+Sample.SAMPLE_DATA_FILE_EXTENSION;
			if (new File(filename).exists() && new File(filename).lastModified() < timeBegan && versionToOverwrite == -1) {
				versionToOverwrite = version-1;
			}
		} while (new File(filename).exists());


		overwriteOptions = new String[] {
				"Rename new file "+trav+Sample.SAMPLE_DATA_FILE_EXTENSION, 
				"Overwrite existing file "+sampleName+(versionToOverwrite==0?"":"."+versionToOverwrite)+Sample.SAMPLE_DATA_FILE_EXTENSION, 
				"Overwrite this and all future files", 
				"Cancel parser"
		};
		
		if (versionToOverwrite != -1) {
			while (new File(dir+HOLD_OPTION_FILE).exists()) {
				try {
					Thread.sleep(500);
				} catch (InterruptedException ie) {}
			}
			if (new File(dir+CANCEL_OPTION_FILE).exists()) {
				return null;
			}
			
			Files.write("", dir+HOLD_OPTION_FILE);
			if (new File(dir+OVERWRITE_OPTION_FILE).exists()) {
				response = 1;
			} else {
				
				if (!System.getProperty("os.name").startsWith("Windows")) {
					log.reportError("Error - the same sample name '"+sampleName+"' is being parsed again and the previous file existed before the current command began. "
							+ "This happens if you inadvertently restarted the parser or if the parser was interrupted and manually restarted. "
							+ "If you would like to start from scratch, the safest thing would be to cancel now and delete all files in the sample directory. "
							+ "To automatically overwrite when there is conflict like this, then create an empty file with the following name: "
							+ dir+OVERWRITE_OPTION_FILE
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
			new File(dir+HOLD_OPTION_FILE).delete();
			switch (response) {
			case 0:
				return dir+trav+Sample.SAMPLE_DATA_FILE_EXTENSION;
			case 2:
				Files.write("", dir+OVERWRITE_OPTION_FILE);
			case 1:
				return dir+sampleName+(versionToOverwrite==0?"":"."+versionToOverwrite)+Sample.SAMPLE_DATA_FILE_EXTENSION;
			case 3:
				Files.write("", dir+CANCEL_OPTION_FILE);
				return null;
			default:
				JOptionPane.showMessageDialog(null, "Should be impossible to obtain this message ("+response+")", "Error", JOptionPane.ERROR_MESSAGE);
				break;
			}
		}

		return dir+trav+Sample.SAMPLE_DATA_FILE_EXTENSION;
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
        long timeBegan;
        boolean parseAtAt;
        int sampIndex;
        String sampleName;
		String[] overwriteOptions;
		int response;
		boolean complete;
		Hashtable<String, Float> allOutliers;
		Vector<String> v;
		Logger log;
		
		log = proj.getLog();
        timeBegan = new Date().getTime();
        new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)+OVERWRITE_OPTION_FILE).delete();
        new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)+HOLD_OPTION_FILE).delete();
        new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)+CANCEL_OPTION_FILE).delete();

		if (!proj.getDir(proj.SOURCE_DIRECTORY).equals("")&&!new File(proj.getDir(proj.SOURCE_DIRECTORY)).exists()) {
			log.reportError("Error - the Project source location is invalid: "+proj.getDir(proj.SOURCE_DIRECTORY));
			return 0;
		}
        
		if (!new File(proj.getFilename(proj.MARKER_POSITION_FILENAME, false, false)).exists()) {
			log.reportError("Error - missing markerPositions: "+proj.getFilename(proj.MARKER_POSITION_FILENAME, false, false));
			return 0;
		}

		delimiter = proj.getSourceFileDelimiter();
		idHeader = proj.getProperty(proj.ID_HEADER);
		log.report(ext.getTime()+"\tSearching for "+proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)+" files in: "+proj.getDir(proj.SOURCE_DIRECTORY));
		files = Files.list(proj.getDir(proj.SOURCE_DIRECTORY), proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);
		
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
		if (new File(proj.getProjectDir()+"fixes.dat").exists()) {
			log.report("Also found a 'fixes.dat' file in the project directory, which will be used to rename samples");
			fixes = HashVec.loadFileToHashString(proj.getProjectDir()+"fixes.dat", false);
		} else {
			log.report("Did not find a 'fixes.dat' file; assuming there are no swapped samples to rename");
		}

		try {
			reader = Files.getAppropriateReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[0]);
			log.report("Found appropriate reader for: "+proj.getDir(proj.SOURCE_DIRECTORY)+files[0]);
			count = 0;
			do {
				line = reader.readLine().trim().split(delimiter, -1);
				count++;
//				if (count < 20) {
//					log.report(Array.toStr(line));
//				}
			} while (reader.ready() && count < 1000 && (ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || (!idHeader.equals(FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line) == -1)));
			
			// If we reached the end of the file, it means that we didn't find the header we are looking for
			// The most common cause of this is that the delimiter was misspecified
			// The following code checks all of the common delimiters (tab, comma, space) and determines which one to use when it tries for a second time
			if (!reader.ready() || count == 1000) {
				log.reportError("Could not find a header with the following tokens: "+Array.toStr(SNP_HEADER_OPTIONS[0]));
				log.reportError("   Perhaps the delimiter, which is currently set to "+proj.getProperty(proj.SOURCE_FILE_DELIMITER)+", is incorrect? This can be corrected in the file "+proj.getPropertyFilename()+". In the meantime, the most stable delimiter will be determined for you...");

				reader.close();
				reader = Files.getAppropriateReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[0]);
				delimiterCounts = new int[DELIMITERS.length][count];
				for (int i = 0; i < count; i++) {
					temp = reader.readLine();
					for (int j = 0; j < DELIMITERS.length; j++) {
						delimiterCounts[j][i] = ext.countInstancesOf(temp, DELIMITERS[j]);
					}
				}
				delimiter = null;
				for (int j = 0; j < DELIMITERS.length; j++) {
					if (Array.quantWithExtremeForTie(delimiterCounts[j], 0.5) > 4 && Array.quantWithExtremeForTie(delimiterCounts[j], 0.9) - Array.quantWithExtremeForTie(delimiterCounts[j], 0.1) == 0) {
						if (delimiter == null) {
							delimiter = DELIMITERS[j];
						} else {
							proj.message("Could not auto-detect the delimiter used in the Final Reports file: could be '"+delimiter+"' or '"+DELIMITERS[j]+"'");
							return 0;
						}
					}
				}
				reader.close();
				
				if (delimiter == null) {
					proj.message("Failed to auto-detect the delimiter used in the Final Reports file; exiting");
					return 0;
				}
				
				log.reportError("      - determined delimiter to be "+DELIMITER_DESCRIPTIONS[ext.indexOfStr(delimiter, DELIMITERS)]);

				// Tries again to determine the header fields and column names
				reader = Files.getAppropriateReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[0]);
				do {
					line = reader.readLine().trim().split(delimiter, -1);
				} while (reader.ready() && (ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1 || (!idHeader.equals(FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line) == -1)));
			}

			// check immediately to make sure these fields are valid
			indices = ext.indexFactors(Sample.DATA_FIELDS, line, false, true, false, false); // dataIndices // TODO turn verbose back on after removing Raw X, Raw Y, Theta and R
			if (indices[3] == -1 || indices[4] == -1) {
				log.reportError("Error - at the very least the files need to contain "+Array.toStr(Sample.DATA_FIELDS[3], "/")+" and "+Array.toStr(Sample.DATA_FIELDS[4], "/"));
				log.reportError("      - failed to see that in "+files[0]);
				log.reportError(Array.toStr(line));
				
				return 0;
			}
			// TODO check different fields depending upon Affy/Illumina flag
			indices = ext.indexFactors(Sample.GENOTYPE_FIELDS, line, false, true, true, false); // genotypeIndices
			if (indices[0] == -1 || indices[1] == -1) {
				log.reportError("Error - the files need to contain "+Array.toStr(Sample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(Sample.GENOTYPE_FIELDS[1], "/"));
				return 0;
			}
			if (indices[2] == -1 || indices[3] == -1) {
				abLookupRequired = true;
			}

			if (!idHeader.equals(FILENAME_AS_ID_OPTION)) {
				ext.indexFactors(new String[] { idHeader }, line, false, true); // sampIndex
			}
			snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, true, true)[0];

			idHeader = proj.getProperty(proj.ID_HEADER);
			reader.mark(1000);


			if (idHeader.equals(FILENAME_AS_ID_OPTION)) {
				sampleName = files[0].substring(0, files[0].indexOf(proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)));
				sampIndex = -7;
				parseAtAt = false;
			} else {
				sampIndex = ext.indexFactors(new String[] { idHeader }, line, false, true)[0];
				line = reader.readLine().split(delimiter);
//				parseAtAt = Boolean.parseBoolean(proj.getProperty(proj.PARSE_AT_AT_SYMBOL));
				parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
				if (parseAtAt && line[sampIndex].indexOf("@") == -1) {
					log.reportError("Error - " + idHeader + " '" + line[sampIndex] + "' did not contain an @ sample; if your ID's do not naturally contain at symbols, then set " + proj.PARSE_AT_AT_SYMBOL + " in the project properties file to FALSE");
					parseAtAt = false;
				}
				sampleName = parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex];
			}
			if (new File(proj.getDir(proj.MARKER_DATA_DIRECTORY, false, false)+"markers.0.mdRAF").exists()) {
				
				if (!System.getProperty("os.name").startsWith("Windows")) {
					log.reportError("Error - Marker data (at least the first file 'markers.0.mdRAF') have already been parsed. This happens if you had previously transposed the data or if the parser was interrupted and manually restarted. Please delete the marker directory to start over from scratch.");
					return 0;
				}

				overwriteOptions = new String[] {
						"Delete All", 
						"Cancel parser"
				};
				
				response = JOptionPane.showOptionDialog(null, 
						"Marker data (at least the first file 'markers.0.mdRAF') have already been parsed.\n"+
						"This happens if you had previously transposed the data or if the parser was interrupted and manually restarted.\n"+
						"If you would like to start from scratch, select \"Delete All\" earlier files.\n"+
						"Otherwise, press cancel.\n"+
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

			if (new File(proj.getDir(proj.SAMPLE_DIRECTORY, true) + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION).exists() || new File(proj.getDir(proj.SAMPLE_DIRECTORY, true) + ext.replaceWithLinuxSafeCharacters(sampleName, true) + Sample.SAMPLE_DATA_FILE_EXTENSION).exists()) {

				overwriteOptions = new String[] {
						"Delete All", 
						"Customize", 
						"Cancel parser"
				};
				
				response = JOptionPane.showOptionDialog(null, 
						"These data (at least the first sample '"+sampleName+"') have already been parsed.\n"+
						"This happens if you inadvertently restarted the parser or if the parser was interrupted and manually restarted.\n"+
						"If you would like to start from scratch, select \"Delete All\" earlier files.\n"+
						"Otherwise, cancel or you can \"Customize\" and determine what to do on a sample-by-sample basis.\n"+
						"What would you like to do?"
					, "Samples already exist", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null, overwriteOptions, overwriteOptions[2]);

				switch (response) {
				case -1:
					return 0;
				case 0:
					deleteAllFilesInSampleDirectory(proj);
					break;
				case 1:
					// keep "outlier.ser"
					break;
				case 2:
					return 0;
				default:
					proj.message("Should be impossible to obtain this message ("+response+")");
					return 0;
				}
			}
			
			TransposeData.deleteOlderRafs(proj.getDir(proj.SAMPLE_DIRECTORY, true), new String[] {"outliers"}, new String[] {".ser"}, true, new String[] {"outliers.ser"});

			reader.reset();	
			
//			if (Boolean.parseBoolean(proj.getProperty(proj.LONG_FORMAT))) {
//				reader.close();
//				createFilesFromLongFormat(proj, files, idHeader, fixes, delimiter, abLookupRequired, timeBegan);
//				return 0;
//			}

			markerNameHash = new Hashtable<String, Integer>(500000);
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter);
				trav = line[snpIndex];
				if (trav.equals("") || trav.equals("0")) {
					trav = "Blank"+count;
					count++;
				}
				if (markerNameHash.containsKey(trav)) {
					
//					if (!Boolean.parseBoolean(proj.getProperty(proj.LONG_FORMAT))) {
					if (!proj.getProperty(proj.LONG_FORMAT)) {
						log.reportError("The same marker ('" + trav + "') was seen twice...");
					}
					
					if (idHeader.equals(FILENAME_AS_ID_OPTION)) {
//						if (!Boolean.parseBoolean(proj.getProperty(proj.LONG_FORMAT))){
						if (!proj.getProperty(proj.LONG_FORMAT)){
							log.report("... this could mean that the file contains multiple samples. This should not happen when " + proj.ID_HEADER + " is set to " + FILENAME_AS_ID_OPTION);
						}else{
							log.reportError("Error - "+proj.ID_HEADER + " was set to " + FILENAME_AS_ID_OPTION + ", which is invalid for Long Format files");
						}
						return 0;
					} else if (!(parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex]).equals(sampleName)) {
//						if (!Boolean.parseBoolean(proj.getProperty(proj.LONG_FORMAT))) {
						if (!proj.getProperty(proj.LONG_FORMAT)) {
							log.reportError("... and the sample name changed from " + sampleName + " to " + (parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@")) : line[sampIndex]) + ", so this must be a Long Format file. The property file will be updated to reflect this, and an attempt will be made to launch the Long Format file processor now.");
//							proj.setProperty(proj.LONG_FORMAT, "TRUE");
							proj.setProperty(proj.LONG_FORMAT, Boolean.TRUE);
							proj.saveProperties();
						}
						reader.close();
						//we should have all markers now...
						markerNames = Array.toStringArray(markerNameHash);
						keys = Markers.orderMarkers(markerNames, proj.getFilename(proj.MARKER_POSITION_FILENAME), proj.getFilename(proj.MARKERSET_FILENAME, true, true), proj.getLog());
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
			log.reportError("Error: file \""+proj.getDir(proj.SOURCE_DIRECTORY)+files[0]+"\" not found in current directory");
			return 0;
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+proj.getDir(proj.SOURCE_DIRECTORY)+files[0]+"\"");
			log.reportException(ioe);
			return 0;
		}

		markerNames = Array.toStringArray(markerNameHash);
		keys = Markers.orderMarkers(markerNames, proj.getFilename(proj.MARKER_POSITION_FILENAME), proj.getFilename(proj.MARKERSET_FILENAME, true, true), proj.getLog());
		if (keys == null) {
			return checkForSNP_Map(proj, log);
		}
		keysKeys = Sort.quicksort(keys); // very important
		fingerprint = proj.getMarkerSet().getFingerprint();
		log.report("There are "+markerNames.length+" markers being processed (fingerprint: "+fingerprint+")");

		lookup = getABLookup(abLookupRequired, markerNames, proj);

		fileCabinet = new Vector<Vector<String>>();
		for (int i = 0; i<numThreads; i++) {
			fileCabinet.add(new Vector<String>());
		}
		for (int i = 0; i<files.length; i++) {
			fileCabinet.elementAt(i%numThreads).add(files[i]);
		}
		threads = new Thread[numThreads];
		for (int i = 0; i<numThreads; i++) {
			threads[i] = new Thread(new ParseIllumina(proj, fileCabinet.elementAt(i).toArray(new String[fileCabinet.elementAt(i).size()]), markerNames, keysKeys, lookup, delimiter, fingerprint, fixes, timeBegan, i));
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
		
		SampleList.generateSampleList(proj).writeToTextFile(proj.getProjectDir()+"ListOfSamples.txt");

		allOutliers = new Hashtable<String, Float>();


		v = new Vector<String>();
		
		for (int i = 0; i<numThreads; i++) {
			if (new File(proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers" + i + ".ser").exists()) {
				allOutliers.putAll((Hashtable<String, Float>) Files.readSerial(proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers" + i + ".ser"));
				new File(proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers" + i + ".ser").delete();
			}
			if (new File(proj.getProjectDir()+"FYI_IDS_WERE_CHANGED"+i+".txt").exists()) {
				v.add(proj.getProjectDir()+"FYI_IDS_WERE_CHANGED"+i+".txt");
			}
		}
		if (v.size() > 0) {
			Files.cat(Array.toStringArray(v), proj.getProjectDir()+"FYI_IDS_WERE_CHANGED.txt", Array.intArray(v.size(), 1), log);
			for (int i = 0; i < v.size(); i++) {
				new File(v.elementAt(i)).delete();
			}
		}
		
		if (allOutliers.size()>0) {
			Files.writeSerial(allOutliers, proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers.ser");
		}
		if (abLookupRequired && !Files.exists(proj.getFilename(proj.AB_LOOKUP_FILENAME, false, false))) {
			return 6;
		} else {
			return 1;
		}
	}

	private static int checkForSNP_Map(Project proj, Logger log) {
		String filename;
		filename = proj.getLocationOfSNP_Map();
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
		filesToDelete = Files.list(proj.getDir(proj.SAMPLE_DIRECTORY), Sample.SAMPLE_DATA_FILE_EXTENSION, false);
		for (int i = 0; i < filesToDelete.length; i++) {
			new File(proj.getDir(proj.SAMPLE_DIRECTORY) + filesToDelete[i]).delete();
		}
		new File(proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers.ser").delete();
	}

	public static void deleteAllFilesInMarkerDataDirectory(Project proj) {
		String[] filesToDelete;
		
		filesToDelete = Files.list(proj.getDir(proj.MARKER_DATA_DIRECTORY), MarkerData.MARKER_DATA_FILE_EXTENSION, false);
		for (int i = 0; i < filesToDelete.length; i++) {
			new File(proj.getDir(proj.MARKER_DATA_DIRECTORY) + filesToDelete[i]).delete();
		}
		new File(proj.getDir(proj.MARKER_DATA_DIRECTORY, true) + "outliers.ser").delete();
	}

	public static char[][] getABLookup(boolean abLookupRequired, String[] markerNames, Project proj) {
		ABLookup abLookup; 
		char[][] lookup;
		
		if (abLookupRequired && Files.exists(proj.getFilename(proj.AB_LOOKUP_FILENAME))) {
			abLookup = new ABLookup(markerNames, proj.getFilename(proj.AB_LOOKUP_FILENAME), true, false, proj.getLog());
			lookup = abLookup.getLookup();
            if (lookup == null) {
            	proj.getLog().reportError("Warning - filed to provide columns \""+Sample.GENOTYPE_FIELDS[2][0]+"\" / \""+Sample.GENOTYPE_FIELDS[3][0]+"\" and the specificed AB_lookup file '"+proj.getProperty(proj.AB_LOOKUP_FILENAME)+"' does not exist; you'll need reconstruct the B allele for analysis");
            } else {
            	abLookup.writeToFile(proj.getProjectDir()+"checkAB.xln", proj.getLog());
            }
		} else {
			lookup = null;
		}
		
		return lookup;
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
        
		//Markers.orderMarkers(null, proj.getFilename(proj.MARKER_POSITION_FILENAME), proj.getFilename(proj.MARKERSET_FILENAME, true, true), proj.getLog());
        // creates and serializes the markers.bim, the returned keys are not used here as the markerIndices fill that purpose
        markerSet = proj.getMarkerSet();
        markerNames = markerSet.getMarkerNames();
		fingerprint = proj.getMarkerSet().getFingerprint();
		
		abLookup = getABLookup(abLookupRequired, markerNames, proj);
		
		markerIndices = new Hashtable<String, Integer>();
		for (int i = 0; i < markerNames.length; i++) {
			markerIndices.put(markerNames[i], new Integer(i));
		}
		
		log.report("There were "+markerNames.length+" markers present in '"+proj.getFilename(proj.MARKERSET_FILENAME, true, true)+"' that will be processed from the source files (fingerprint: "+fingerprint+")");
		
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
		
		if (Files.exists(proj.getProjectDir()+"FYI_IDS_WERE_CHANGED.txt")) {
			Files.backup("FYI_IDS_WERE_CHANGED.txt", proj.getProjectDir(), proj.getProjectDir(), true);
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
					reader = Files.getAppropriateReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[i]);
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || ext.indexOfStr(idHeader, line)==-1));

					log.report("Searching: "+Array.toStr(line));
					dataIndices = ext.indexFactors(Sample.DATA_FIELDS, line, false, true, false, false);
					genotypeIndices = ext.indexFactors(Sample.GENOTYPE_FIELDS, line, false, true, false, false);
					sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];
					snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
					
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
//					parseAtAt = Boolean.parseBoolean(proj.getProperty(proj.PARSE_AT_AT_SYMBOL));
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
									writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"FYI_IDS_WERE_CHANGED.txt", true));
									if (renamedIDsHash.size() == 0) {
										writer.println("The following IDs were changed so that spaces are removed and so that they could be used as valid filenames:");
									}
									writer.println(temp+"\t"+trav);
									writer.close();
								} catch (Exception e) {
									log.reportError("Error writing to " + proj.getProjectDir()+"FYI_IDS_WERE_CHANGED.txt");
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
								
								filename = determineFilename(proj.getDir(proj.SAMPLE_DIRECTORY, true), sampleName, timeBegan, log);
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
								if (new File(proj.getDir(proj.SAMPLE_DIRECTORY, true) + (fixes.containsKey(trav)?fixes.get(trav):trav) + Sample.SAMPLE_DATA_FILE_EXTENSION).exists()) {
									log.reportError("Warning - marker data must be out of order, becaue we're seeing "+trav+(fixes.containsKey(trav)?"-->"+fixes.get(trav):"")+" again at line "+count);
									samp = Sample.loadFromRandomAccessFile(proj.getDir(proj.SAMPLE_DIRECTORY, true) + (fixes.containsKey(trav)?fixes.get(trav):trav) + Sample.SAMPLE_DATA_FILE_EXTENSION, proj.getJarStatus());
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
				if (new File(proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers.ser").exists()) {
					log.reportError("Error - the following file already exists: " + proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers.ser");
					return 0;
				} else {
					Files.writeSerial(allOutliers, proj.getDir(proj.SAMPLE_DIRECTORY, true) + "outliers.ser");
				}
			}

			log.report(ext.getTime()+"\tfinished");
		} catch (Exception e) {
			log.reportException(e);
		}

        
		log.report(ext.getTime()+"\t"+"Parsed "+count+" sample(s)");
		SampleList.generateSampleList(proj).writeToTextFile(proj.getProjectDir()+"ListOfSamples.txt");

		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"ListOfMarkers.txt"));
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

        new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)+OVERWRITE_OPTION_FILE).delete();
        new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)+HOLD_OPTION_FILE).delete();
        new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)+CANCEL_OPTION_FILE).delete();
       
		if (abLookupRequired && !Files.exists(proj.getFilename(proj.AB_LOOKUP_FILENAME, false, false))) {
			return 6;
		} else {
			return 1;
		}
	}

	public static void mapFilenamesToSamples(Project proj, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int sampIndex;
		String[] files;
		String idHeader, delimiter;
		Logger log;
		boolean longFormat, done;
		String prev;

		log = proj.getLog();
		if (!Files.exists(proj.getDir(proj.SOURCE_DIRECTORY, false, false))) {
			proj.message("Source directory does not exist; change SOURCE_DIRECTORY= to point to the proper files");
			return;
		}
		
		delimiter = proj.getSourceFileDelimiter();
		longFormat = proj.getBoolean(proj.LONG_FORMAT);
		idHeader = proj.getProperty(proj.ID_HEADER);
		log.report(ext.getTime()+"\tSearching for "+proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)+" files in: "+proj.getDir(proj.SOURCE_DIRECTORY));
		files = Files.list(proj.getDir(proj.SOURCE_DIRECTORY), proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);
		log.report("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" to parse");

		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+filename));
			for (int i = 0; i<files.length; i++) {
				try {
					reader = Files.getAppropriateReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[i]);
					do {
						line = reader.readLine().trim().split(delimiter);
					} while (reader.ready()&&(line.length<3 || ext.indexOfStr(idHeader, line)==-1));

					if (!reader.ready()) {
						log.reportError("Error - went through enitre file without finding a line containing the user-defined ID header: "+idHeader);
						return;
					}
					sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];

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
					log.reportError("Error: file \""+files[i]+"\" not found in "+proj.getDir(proj.SOURCE_DIRECTORY));
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
		files = Files.list(proj.getDir(proj.SOURCE_DIRECTORY), proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);
		if (files.length == 0) {
			log.reportError("Error - no files to parse");
			return;
		}
		log.report("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" to parse");

		idHeader = proj.getProperty(proj.ID_HEADER);
		delimiter = proj.getSourceFileDelimiter();
		hash = new Hashtable<String, String[]>();
		for (int i = 0; i<files.length; i++) {
			if (new File("report").exists()) {
				writeToLookupFile(proj, hash, i+1);
				new File("report").delete();
			}
			try {
				log.report(ext.getTime()+"\t"+(i+1)+" of "+files.length);
				reader = Files.getAppropriateReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[i]);
				do {
					line = reader.readLine().trim().split(delimiter, -1);
				} while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || (!idHeader.equals(FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line)==-1)));

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
		writeToLookupFile(proj, hash, 0);
	}
	
	public static void writeToLookupFile(Project proj, Hashtable<String,String[]> hash, int fileNumber) {
		PrintWriter writer;
		String[] keys;
		Logger log;

		log = proj.getLog();
		log.report("Writing to file...", false, true);
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"alleleLookup"+(fileNumber>0?"_atFile"+fileNumber:"")+".xln"));
			keys = HashVec.getKeys(hash, false, false);
			writer.println("SNP\t"+Array.toStr(Sample.ALL_STANDARD_GENOTYPE_FIELDS));
			for (int k = 0; k < keys.length; k++) {
				writer.println(keys[k]+"\t"+Array.toStr(hash.get(keys[k])));
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + proj.getProjectDir()+"alleleLookup.xln");
			log.reportException(e);
		}
		log.report("done");
		
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
				mapFilenamesToSamples(proj, mapOutput);
			} else if (parseAlleleLookupFromFinalReports) {
				parseAlleleLookupFromFinalReports(proj);
			} else {
				createFiles(proj, numThreads);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
