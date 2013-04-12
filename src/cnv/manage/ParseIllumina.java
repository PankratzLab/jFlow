// might want to migrate the allele lookup/conversion to a new class since it would apply to Affy data too, same potentially for generateMarkerPositions if markerPositions is used in Affy algorithm
// "AB_lookup.dat" is necessary if the files do not contain {"Allele1 - AB"}/{"Allele2 - AB}
package cnv.manage;

import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;

import cnv.filesys.ABLookup;
import cnv.filesys.FullSample;
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

	private Project proj;
	private String[] files;
	private String[] markerNames;
	private int[] keysKeys;
	private long fingerprint;
	private char[][] abLookup;
	private Hashtable<String,String> fixes;

	public ParseIllumina(Project proj, String[] files, String[] markerNames, int[] keysKeys, char[][] abLookup, long fingerprint, Hashtable<String,String> fixes) {
		this.proj = proj;
		this.files = files;
		this.markerNames = markerNames;
		this.keysKeys = keysKeys;
		this.abLookup = abLookup;
		this.fingerprint = fingerprint;
		this.fixes = fixes;
	}

	public void run() {
		BufferedReader reader;
		String[] line;
		int count, version, snpIndex, sampIndex, key;
		String trav;
		int[] dataIndices, genotypeIndices;
		boolean parseAtAt;

		FullSample samp;
		String sampleName;
		float[][] data;
		byte[][] genotypes;
		boolean ignoreAB;
		String idHeader;
		String delimiter;
		
//		try {
//			PrintWriter writer = new PrintWriter(new FileWriter(files[0]+"_list.xln"));
//			for (int j = 0; j<files.length; j++) {
//				writer.println(files[j]);
//            }
//	        writer.close();
//	        return;
//        } catch (Exception e) {
//	        System.err.println("Error writing to "+files[0]+"_list.xln");
//	        e.printStackTrace();
//        }

		idHeader = proj.getProperty(Project.ID_HEADER);
		delimiter = proj.getSourceFileDelimiter();
        try {
			for (int i = 0; i<files.length; i++) {
				try {
					System.out.println(ext.getTime()+"\t"+(i+1)+" of "+files.length);
//					reader = new BufferedReader(new FileReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[i]));
					reader = Files.getAppropriateReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[i]);
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || ext.indexOfStr(idHeader, line)==-1));

					dataIndices = ext.indexFactors(FullSample.DATA_FIELDS, line, false, true, false, false);
					genotypeIndices = ext.indexFactors(FullSample.GENOTYPE_FIELDS, line, false, true, false, false);
					sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];
					snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
					
					if (dataIndices[3] == -1 || dataIndices[4] == -1) {
						System.err.println("Error - File format not consistent! At the very least the files need to contain "+Array.toStr(FullSample.DATA_FIELDS[3], "/")+" and "+Array.toStr(FullSample.DATA_FIELDS[4], "/"));
						return;
					}
					if (genotypeIndices[0] == -1 || genotypeIndices[1] == -1) {
						System.err.println("Error - File format not consistent! The files need to contain "+Array.toStr(FullSample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(FullSample.GENOTYPE_FIELDS[1], "/"));
						return;
					}
					if ((genotypeIndices[2] == -1 || genotypeIndices[3] == -1) && abLookup == null) {
						ignoreAB = true;
					} else {
						ignoreAB = false;
					}

					sampleName = "";
					data = new float[FullSample.DATA_FIELDS.length][];
					for (int j = 0; j<data.length; j++) {
						if (dataIndices[j] != -1) {
							data[j] = new float[markerNames.length];
						}
                    }
					genotypes = new byte[2][];
					genotypes[0] = Array.byteArray(markerNames.length, Byte.MIN_VALUE);
					if (!ignoreAB) {
						genotypes[1] = Array.byteArray(markerNames.length, Byte.MIN_VALUE);
					}
					
					count = 0;
					parseAtAt = Boolean.parseBoolean(proj.getProperty(Project.PARSE_AT_AT_SYMBOL));
					while (reader.ready()) {
						line = reader.readLine().split(delimiter);
						if (parseAtAt&&line[sampIndex].indexOf("@")==-1) {
							System.err.println("Error - "+idHeader+" '"+line[sampIndex]+"' did not contain an @ sample");
							parseAtAt = false;
						}
						trav = parseAtAt?line[sampIndex].substring(0, line[sampIndex].indexOf("@")):line[sampIndex];

						if (count==0) {
							sampleName = trav;
						} else if (!trav.equals(sampleName)) {
							System.err.println("Found more than one ID in file "+files[i]+"(found "+trav+", expecting "+sampleName+")");
							return;
						} else if (count>markerNames.length) {
							System.err.println("Error - expecting only "+markerNames.length+" markers and found more than that in file "+files[i]);
							return;
						} else if (!markerNames[count].equals(line[snpIndex]) && !markerNames[count].startsWith("Blank")) {
							System.err.println("Found "+line[snpIndex]+" at marker #"+(count+1)+" in file "+files[i]+"; expecting "+markerNames[count]);
							return;
						}
						key = keysKeys[count];
						// System.out.println(count+"\t"+markerNames[count]+"\t"+key);
						for (int j = 0; j<FullSample.DATA_FIELDS.length; j++) {
							try {
								if (dataIndices[j] != -1) {
									if (line[dataIndices[j]].equals("")) {
										data[j][key] = Float.NaN;
									} else {
										data[j][key] = Float.parseFloat(line[dataIndices[j]]);
									}
								}
							} catch (NumberFormatException nfe) {
								System.err.println("Error - failed to parse '"+line[dataIndices[j]]+"' into a valid "+Array.toStr(FullSample.DATA_FIELDS[j], "/"));
								return;
							}
						}

						genotypes[0][key] = (byte)ext.indexOfStr(line[genotypeIndices[0]]+line[genotypeIndices[1]], FullSample.ALLELE_PAIRS);
						if (genotypes[0][key] == -1) {
							if (ext.indexOfStr(line[genotypeIndices[0]]+line[genotypeIndices[1]], FullSample.ALT_NULLS) == -1) {
								System.err.println("Error - failed to lookup "+line[genotypeIndices[0]]+line[genotypeIndices[1]]+" for marker "+markerNames[count]+" of sample "+files[i]);
							} else {
								genotypes[0][key] = 0;
							}								
						}
						if (ignoreAB) {
							// do nothing, will need to use these files to determine AB lookup table
						} else if (abLookup == null) {
							genotypes[1][key] = (byte)ext.indexOfStr(line[genotypeIndices[2]]+line[genotypeIndices[3]], FullSample.AB_PAIRS);
						} else {
							if (genotypes[0][key] == 0) {
								genotypes[1][key] = -1;
							} else {
								genotypes[1][key] = 0;
								for (int j = 0; j<2; j++) {
									if (line[genotypeIndices[j]].charAt(0) == abLookup[count][1]) {
										genotypes[1][key]++;
									} else if (line[genotypeIndices[j]].charAt(0) != abLookup[count][0]) {
										System.err.println("Error - alleles for individual '"+line[sampIndex]+"' ("+line[genotypeIndices[0]]+"/"+line[genotypeIndices[1]]+") do not match up with the defined AB lookup alleles ("+abLookup[count][0]+"/"+abLookup[count][1]+") for marker "+markerNames[count]);
									}
                                }
							}
						}

						count++;
					}
					reader.close();
					if (count!=markerNames.length) {
						System.err.println("Error - expecting "+markerNames.length+" markers and only found "+count+" in file "+files[i]);
						return;
					}

					if (fixes.containsKey(sampleName)) {
						sampleName = fixes.get(sampleName);
					}
					
					version = 0;
					do {
						trav = sampleName+(version==0?"":"."+version);
						version++;
					} while (new File(proj.getDir(Project.SAMPLE_DIRECTORY, true)+trav+".fsamp").exists());

					samp = new FullSample(sampleName, fingerprint, data, genotypes);
					samp.serialize(proj.getDir(Project.SAMPLE_DIRECTORY, true)+trav+".fsamp");
					samp.convertToSample().serialize(proj.getDir(Project.IND_DIRECTORY, true)+trav+".samp");
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[i]+"\" not found in current directory");
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[i]+"\"");
					return;
				}
			}

			SampleList.generateSampleList(proj).writeToTextFile(proj.getProjectDir()+"ListOfSamples.txt");

			System.out.println(ext.getTime()+"\tfinished");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void createFiles(Project proj, int numThreads) {
		BufferedReader reader;
		String[] line, markerNames, files;
		ArrayList<String> alNames;
		int snpIndex;
		int[] keys, keysKeys, indices;
		long fingerprint;
		Thread[] threads;
		Vector<Vector<String>> fileCabinet;
		String trav;
		int count;
//		Hashtable<String, char[]> hash;
		boolean abLookupRequired; // , problem;
		char[][] lookup;
        Hashtable<String,String> fixes;
        String idHeader, delimiter;
        ABLookup abLookup;
        String temp;
        int[][] delimiterCounts;

		if (!proj.getDir(Project.SOURCE_DIRECTORY).equals("")&&!new File(proj.getDir(Project.SOURCE_DIRECTORY)).exists()) {
			System.err.println("Error - the Project source location is invalid: "+proj.getDir(Project.SOURCE_DIRECTORY));
			return;
		}
        
		if (!new File(proj.getFilename(Project.MARKER_POSITION_FILENAME, false, false)).exists()) {
			System.err.println("Error - missing markerPositions: "+proj.getFilename(Project.MARKER_POSITION_FILENAME, false, false));
			return;
		}

		delimiter = proj.getSourceFileDelimiter();
		idHeader = proj.getProperty(Project.ID_HEADER);
		System.out.println(ext.getTime()+"\tSearching for "+proj.getProperty(Project.SOURCE_FILENAME_EXTENSION)+" files in: "+proj.getDir(Project.SOURCE_DIRECTORY));
		files = Files.list(proj.getDir(Project.SOURCE_DIRECTORY), proj.getProperty(Project.SOURCE_FILENAME_EXTENSION), false);
		
		if (files.length == 0) {
			System.err.println("Error - no files to parse");
			return;
		}

		abLookupRequired = false;
		System.out.println("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" to parse");
		fixes = new Hashtable<String,String>();
		if (new File(proj.getProjectDir()+"fixes.dat").exists()) {
			System.out.println("Also found a 'fixes.dat' file in the project directory, which will be used to rename samples");
			fixes = HashVec.loadFileToHashString(proj.getProjectDir()+"fixes.dat", false);
		} else {
			System.out.println("Did not find a 'fixes.dat' file; assuming you don't want to rename any IDs");
		}
		alNames = new ArrayList<String>();
		try {
//			reader = new BufferedReader(new FileReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[0]));
			reader = Files.getAppropriateReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[0]);
			System.out.println("Found appropriate reader for: "+proj.getDir(Project.SOURCE_DIRECTORY)+files[0]);
			count = 0;
			do {
				line = reader.readLine().trim().split(delimiter, -1);
				count++;
				if (count < 20) {
					System.out.println(Array.toStr(line));
				}
			} while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || ext.indexOfStr(idHeader, line)==-1));
			
			if (!reader.ready()) {
				System.err.println("Error - reached the end of the file without finding a line with the following tokens: "+Array.toStr(SNP_HEADER_OPTIONS[0]));
				System.err.println("      - perhaps the delimiter is set incorrectly? Determing most stable delimiter...");

				reader.close();
				reader = Files.getAppropriateReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[0]);
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
							JOptionPane.showMessageDialog(null, "Could not auto-detect the delimiter used in the Final Reports file: could be '"+delimiter+"' or '"+DELIMITERS[j]+"'", "Error", JOptionPane.ERROR_MESSAGE);
							return;
						}
					}
				}
				reader.close();
				
				if (delimiter == null) {
					JOptionPane.showMessageDialog(null, "Failed to auto-detect the delimiter used in the Final Reports file; exitting", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				}
				
				System.err.println("      - determined delimiter to be '"+delimiter+"'");

				reader = Files.getAppropriateReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[0]);
				do {
					line = reader.readLine().trim().split(delimiter, -1);
				} while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || ext.indexOfStr(idHeader, line)==-1));
				System.out.println(1);
			}

			// check immediately to make sure these fields are valid
			indices = ext.indexFactors(FullSample.DATA_FIELDS, line, false, true, true, false); // dataIndices
			if (indices[3] == -1 || indices[4] == -1) {
				System.err.println("Error - at the very least the files need to contain "+Array.toStr(FullSample.DATA_FIELDS[3], "/")+" and "+Array.toStr(FullSample.DATA_FIELDS[4], "/"));
				System.err.println("      - failed to see that in "+files[0]);
				System.err.println(Array.toStr(line));
				
				return;
			}
			indices = ext.indexFactors(FullSample.GENOTYPE_FIELDS, line, false, true, true, false); // genotypeIndices
			if (indices[0] == -1 || indices[1] == -1) {
				System.err.println("Error - the files need to contain "+Array.toStr(FullSample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(FullSample.GENOTYPE_FIELDS[1], "/"));
				return;
			}
			if (indices[2] == -1 || indices[3] == -1) {
				abLookupRequired = true;
			}
			
			ext.indexFactors(new String[] {idHeader}, line, false, true); // sampIndex
			snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, true, true)[0];

			if (Boolean.parseBoolean(proj.getProperty(Project.LONG_FORMAT))) {
				createFilesFromLongFormat(proj, files, idHeader, fixes, delimiter, abLookupRequired);
				return;
			}

			count = 1;
			while (reader.ready()) {
				trav = reader.readLine().trim().split(delimiter)[snpIndex];
				if (trav.equals("") || trav.equals("0")) {
					trav = "Blank"+count;
					count++;
				}
				alNames.add(trav);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.getDir(Project.SOURCE_DIRECTORY)+files[0]+"\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.getDir(Project.SOURCE_DIRECTORY)+files[0]+"\"");
			return;
		}

//		new File(proj.getDir(Project.SAMPLE_DIRECTORY, true)).mkdirs();
//		new File(proj.getDir(Project.IND_DIRECTORY)).mkdirs();
//		new File(proj.getDir(Project.DATA_DIRECTORY)).mkdirs();

		markerNames = Array.toStringArray(alNames);
		keys = Markers.orderMarkers(markerNames, proj.getFilename(Project.MARKER_POSITION_FILENAME), proj.getFilename(Project.MARKERSET_FILENAME, true, true));
		if (keys == null) {
			return;
		}
		keysKeys = Sort.quicksort(keys); // very important
		fingerprint = proj.getMarkerSet().getFingerprint();
		System.out.println("There are "+markerNames.length+" markers being processed (fingerprint: "+fingerprint+")");
		
		if (abLookupRequired) {
			abLookup = new ABLookup(markerNames, proj.getProjectDir()+ABLookup.DEFAULT_AB_FILE, true, true);
			lookup = abLookup.getLookup();
            if (lookup == null) {
    			System.err.println("Warning - filed to provide columns \""+FullSample.GENOTYPE_FIELDS[2][0]+"\" / \""+FullSample.GENOTYPE_FIELDS[3][0]+"\" and there is no file \""+ABLookup.DEFAULT_AB_FILE+"\" to compensate; you'll need reconstruct the B allele for analysis");
            } else {
            	abLookup.writeToFile(proj.getProjectDir()+"checkAB.xln");
            }
		} else {
			lookup = null;
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
			threads[i] = new Thread(new ParseIllumina(proj, fileCabinet.elementAt(i).toArray(new String[fileCabinet.elementAt(i).size()]), markerNames, keysKeys, lookup, fingerprint, fixes));
//			threads[i] = new Thread(new ParseIllumina(proj, fileCabinet.elementAt(i).toArray(new String[fileCabinet.elementAt(i).size()]), null, null, null, 0));
			threads[i].start();
			try {
				Thread.sleep(100L);
			} catch (InterruptedException ex) {}
		}
	}

	public static void createFilesFromLongFormat(Project proj, String[] files, String idHeader, Hashtable<String, String> fixes, String delimiter, boolean abLookupRequired) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, markerNames, list;
		long fingerprint;
		MarkerSet markerSet;
		Hashtable<String, Integer> markerIndices;

        System.out.println("Parsing files using the Long Format algorithm");

        Markers.orderMarkers(null, proj.getFilename(Project.MARKER_POSITION_FILENAME), proj.getFilename(Project.MARKERSET_FILENAME, true, true));
        markerSet = proj.getMarkerSet();
        markerNames = markerSet.getMarkerNames();
		fingerprint = proj.getMarkerSet().getFingerprint();
		
		markerIndices = new Hashtable<String, Integer>();
		for (int i = 0; i < markerNames.length; i++) {
			markerIndices.put(markerNames[i], new Integer(i));
		}
		
		System.out.println("There were "+markerNames.length+" markers present in '"+proj.getFilename(Project.MARKERSET_FILENAME, true, true)+"' that will be processed from the source files (fingerprint: "+fingerprint+")");
		
		if (abLookupRequired) {
//			JOptionPane.showMessageDialog(null, "Error - This project requires an AB lookup, but this has not yet been implemented for the Long Format", "Error", JOptionPane.ERROR_MESSAGE);
			System.err.println("Error - This project requires an AB lookup, but this has not yet been implemented for the Long Format");
		}
		

		int snpIndex, sampIndex, key;
		String trav;
		int[] dataIndices, genotypeIndices;
		boolean parseAtAt;

		FullSample samp;
		String sampleName;
		float[][] data;
		byte[][] genotypes;
		boolean ignoreAB;
		CountHash countHash, dupHash;
		boolean done;
		
		countHash = new CountHash();
		dupHash = new CountHash();
        try {
			for (int i = 0; i<files.length; i++) {
				try {
					System.out.println(ext.getTime()+"\t"+(i+1)+" of "+files.length+" ("+files[i]+")");
//					reader = new BufferedReader(new FileReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[i]));
					reader = Files.getAppropriateReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[i]);
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || ext.indexOfStr(idHeader, line)==-1));

					System.err.println("Searching: "+Array.toStr(line));
					dataIndices = ext.indexFactors(FullSample.DATA_FIELDS, line, false, true, false, false);
					genotypeIndices = ext.indexFactors(FullSample.GENOTYPE_FIELDS, line, false, true, false, false);
					sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];
					snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
					
					if (dataIndices[3] == -1 || dataIndices[4] == -1) {
						System.err.println("Error - File format not consistent! At the very least the files need to contain "+Array.toStr(FullSample.DATA_FIELDS[3], "/")+" and "+Array.toStr(FullSample.DATA_FIELDS[4], "/"));
						return;
					}
					if (genotypeIndices[0] == -1 || genotypeIndices[1] == -1) {
						System.err.println("Error - File format not consistent! The files need to contain "+Array.toStr(FullSample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(FullSample.GENOTYPE_FIELDS[1], "/"));
						return;
					}
//					if ((genotypeIndices[2] == -1 || genotypeIndices[3] == -1) && abLookup == null) {
					if ((genotypeIndices[2] == -1 || genotypeIndices[3] == -1)) {
						ignoreAB = true;
					} else {
						ignoreAB = false;
					}

					done = false;
					sampleName = "just starting";
					data = null;
					genotypes = null;
					parseAtAt = Boolean.parseBoolean(proj.getProperty(Project.PARSE_AT_AT_SYMBOL));
					while (!done) {
						if (reader.ready()) {
							line = reader.readLine().split(delimiter);
							if (parseAtAt&&line[sampIndex].indexOf("@")==-1) {
								System.err.println("Error - "+idHeader+" '"+line[sampIndex]+"' did not contain an @ sample");
								parseAtAt = false;
							}
							trav = parseAtAt?line[sampIndex].substring(0, line[sampIndex].indexOf("@")):line[sampIndex];
						} else {
							done = true;
							trav = null;
						}

						if (done || !trav.equals(sampleName)) {
							if (!sampleName.equals("just starting")) {
								if (fixes.containsKey(sampleName)) {
									sampleName = fixes.get(sampleName);
								}
								
								samp = new FullSample(sampleName, fingerprint, data, genotypes);
								samp.serialize(proj.getDir(Project.SAMPLE_DIRECTORY, true)+sampleName+".fsamp");
								samp.convertToSample().serialize(proj.getDir(Project.IND_DIRECTORY, true)+sampleName+".samp");
							}
							if (new File(proj.getDir(Project.SAMPLE_DIRECTORY, true)+trav+".fsamp").exists()) {
								samp = FullSample.load(proj.getDir(Project.SAMPLE_DIRECTORY, true)+(fixes.containsKey(trav)?fixes.get(trav):trav)+".fsamp", proj.getJarStatus());
								data = samp.getAllData();
								genotypes = samp.getAllGenotypes();
							} else {
								data = new float[FullSample.DATA_FIELDS.length][];
								for (int j = 0; j<data.length; j++) {
									if (dataIndices[j] != -1) {
										data[j] = Array.floatArray(markerNames.length, Float.POSITIVE_INFINITY);
									}
			                    }
								genotypes = new byte[2][];
								genotypes[0] = Array.byteArray(markerNames.length, Byte.MIN_VALUE);
								if (!ignoreAB) {
									genotypes[1] = Array.byteArray(markerNames.length, Byte.MIN_VALUE);
								}
							}							
							sampleName = trav;
						}

						if (!done) {
							countHash.add(line[snpIndex]);							
							if (markerIndices.containsKey(line[snpIndex])) {
								key = markerIndices.get(line[snpIndex]);

								for (int j = 0; j<FullSample.DATA_FIELDS.length; j++) {
									try {
										if (dataIndices[j] != -1) {
											if (!(data[j][key]+"").equals("Infinity")) {
												dupHash.add(line[snpIndex]);
												System.err.println("Sample "+trav+" already has data for marker "+line[snpIndex]+" (Was the parsing restarted? Delete the old directories first)");
											}
											data[j][key] = Float.parseFloat(line[dataIndices[j]]);
										}
									} catch (NumberFormatException nfe) {
										System.err.println("Error - failed to parse"+line[dataIndices[j]]+" into a valid "+FullSample.DATA_FIELDS[j]);
										return;
									}
								}

								genotypes[0][key] = (byte)ext.indexOfStr(line[genotypeIndices[0]]+line[genotypeIndices[1]], FullSample.ALLELE_PAIRS);
								if (genotypes[0][key] == -1) {
									if (ext.indexOfStr(line[genotypeIndices[0]]+line[genotypeIndices[1]], FullSample.ALT_NULLS) == -1) {
										System.err.println("Error - failed to lookup "+line[genotypeIndices[0]]+line[genotypeIndices[1]]+" for marker "+line[snpIndex]+" of sample "+files[i]);
									} else {
										genotypes[0][key] = 0;
									}								
								}
								if (!ignoreAB) {
									genotypes[1][key] = (byte)ext.indexOfStr(line[genotypeIndices[2]]+line[genotypeIndices[3]], FullSample.AB_PAIRS);
								}
//								if (ignoreAB) {
//									// do nothing, will need to use these files to determine AB lookup table
//								} else if (abLookup == null) {
//									genotypes[1][key] = (byte)ext.indexOfStr(line[genotypeIndices[2]]+line[genotypeIndices[3]], FullSample.AB_PAIRS);
//								} else {
//									if (genotypes[0][key] == 0) {
//										genotypes[1][key] = -1;
//									} else {
//										genotypes[1][key] = 0;
//										for (int j = 0; j<2; j++) {
//											if (line[genotypeIndices[j]].charAt(0) == abLookup[count][1]) {
//												genotypes[1][key]++;
//											} else if (line[genotypeIndices[j]].charAt(0) != abLookup[count][0]) {
//												System.err.println("Error - alleles for individual '"+line[sampIndex]+"' ("+line[genotypeIndices[0]]+"/"+line[genotypeIndices[1]]+") do not match up with the defined AB lookup alleles ("+abLookup[count][0]+"/"+abLookup[count][1]+") for marker "+markerNames[count]);
//											}
//		                                }
//									}
//								}
							}
						}
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[i]+"\" not found in current directory");
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[i]+"\"");
					return;
				}
			}

			System.out.println(ext.getTime()+"\tfinished");
		} catch (Exception e) {
			e.printStackTrace();
		}

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
			System.err.println("Error writing to " + "ListOfMarkers.txt");
			e.printStackTrace();
		}

	}

	public static void mapFilenamesToSamples(Project proj, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int sampIndex;
		String[] files;
		String idHeader, delimiter;

		delimiter = proj.getSourceFileDelimiter();
		idHeader = proj.getProperty(Project.ID_HEADER);
		System.out.println(ext.getTime()+"\tSearching for "+proj.getProperty(Project.SOURCE_FILENAME_EXTENSION)+" files in: "+proj.getDir(Project.SOURCE_DIRECTORY));
		files = Files.list(proj.getDir(Project.SOURCE_DIRECTORY), proj.getProperty(Project.SOURCE_FILENAME_EXTENSION), false);
		System.out.println("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" to parse");

		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+filename));
			for (int i = 0; i<files.length; i++) {
				try {
//					reader = new BufferedReader(new FileReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[i]));
					reader = Files.getAppropriateReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[i]);
					do {
						line = reader.readLine().trim().split(delimiter);
					} while (reader.ready()&&(line.length<3 || ext.indexOfStr(idHeader, line)==-1));

					if (!reader.ready()) {
						System.err.println("Error - went through enitre file without finding a line containing the user-defined ID header: "+idHeader);
						return;
					}
					sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];

					line = reader.readLine().split(delimiter);
					writer.println(files[i]+"\t"+line[sampIndex]+"\t"+(line[sampIndex].indexOf("@") >= 0?line[sampIndex].split("@")[0]:line[sampIndex]));
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[i]+"\" not found in "+proj.getDir(Project.SOURCE_DIRECTORY));
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[i]+"\"");
					return;
				}
			}
			System.out.println(ext.getTime());
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	
	public static void parseAlleleLookup(Project proj) {
		BufferedReader reader;
		String[] line;
		Hashtable<String, String[]> hash;
		String[] files;
		int[] indices;
		int snpIndex;
		String delimiter, idHeader;
		String[] alleles;
		int expIndex;
		
		files = Files.list(proj.getDir(Project.SOURCE_DIRECTORY), proj.getProperty(Project.SOURCE_FILENAME_EXTENSION), false);
		if (files.length == 0) {
			System.err.println("Error - no files to parse");
			return;
		}
		System.out.println("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" to parse");

		idHeader = proj.getProperty(Project.ID_HEADER);
		delimiter = proj.getSourceFileDelimiter();
		hash = new Hashtable<String, String[]>();
		for (int i = 0; i<files.length; i++) {
			if (new File("report").exists()) {
				writeToLookupFile(proj, hash, i+1);
				new File("report").delete();
			}
			try {
				System.out.println(ext.getTime()+"\t"+(i+1)+" of "+files.length);
//				reader = new BufferedReader(new FileReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[i]));
				reader = Files.getAppropriateReader(proj.getDir(Project.SOURCE_DIRECTORY)+files[i]);
				do {
					line = reader.readLine().trim().split(delimiter, -1);
				} while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || ext.indexOfStr(idHeader, line)==-1));

				snpIndex = ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
				indices = ext.indexFactors(FullSample.ALL_STANDARD_GENOTYPE_FIELDS, line, false, new Logger(null), false, false);
				
				while (reader.ready()) {
					line = reader.readLine().split(delimiter);
					if (hash.containsKey(line[snpIndex])) {
						alleles = hash.get(line[snpIndex]);
					} else {
						if (i!=0) {
							System.err.println("Error - snp '"+line[snpIndex]+"' first seen in file #"+i+" ("+files[i]+") and not earlier");
						}
						hash.put(line[snpIndex], alleles = new String[FullSample.ALL_STANDARD_GENOTYPE_FIELDS.length]);
					}
					for (int j = 0; j < 2; j++) {
						if (ext.indexOfStr(line[indices[j]], FullSample.ALT_NULL) == -1) {
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
								System.err.println("Error - snp '"+line[snpIndex]+"' has a new allele in file #"+i+" ("+files[i]+"): "+line[indices[j]]+" (previously "+alleles[0]+"/"+alleles[1]+")");
								expIndex = -9;
							}
							
							for (int k = 1; k < alleles.length/2; k++) {
								switch (expIndex) {
								case -1:
									if (alleles[k*2+0] != null) {
										System.err.println("Error - how can this be? -1");
									}
									alleles[k*2+0] = line[indices[k*2+j]];
									break;
								case 0:
									if (!line[indices[k*2+j]].equals(alleles[k*2+0])) {
										System.err.println("Error - how can this be? 0");
									}									
									break;
								case -2:
									if (alleles[k*2+1] != null) {
										System.err.println("Error - how can this be? -2");
									}
									alleles[k*2+1] = line[indices[k*2+j]];
									break;
								case 1:
									if (!line[indices[k*2+j]].equals(alleles[k*2+1])) {
										System.err.println("Error - how can this be? 1");
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
				System.err.println("Error: file \""+files[i]+"\" not found in current directory");
				return;
			} catch (IOException ioe) {
				System.err.println("Error reading file \""+files[i]+"\"");
				return;
			}
		}
		writeToLookupFile(proj, hash, 0);
	}
	
	public static void writeToLookupFile(Project proj, Hashtable<String,String[]> hash, int fileNumber) {
		PrintWriter writer;
		String[] keys;
		
		System.out.print("Writing to file...");
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"alleleLookup"+(fileNumber>0?"_atFile"+fileNumber:"")+".xln"));
			keys = HashVec.getKeys(hash, false, false);
			writer.println("SNP\t"+Array.toStr(FullSample.ALL_STANDARD_GENOTYPE_FIELDS));
			for (int k = 0; k < keys.length; k++) {
				writer.println(keys[k]+"\t"+Array.toStr(hash.get(keys[k])));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + proj.getProjectDir()+"alleleLookup.xln");
			e.printStackTrace();
		}
		System.out.println("done");
		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj;
		String filename = Project.DEFAULT_PROJECT;
		boolean map = false;
		int numThreads = 1;
//		boolean parseABlookup = false;
		boolean parseAlleleLookup = false;
		String mapOutput = "filenamesMappedToSamples.txt";

		String usage = "\n"+
		"cnv.manage.ParseIllumina requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) number of threads to use (i.e. threads="+numThreads+" (default))\n"+
		" OPTIONAL:\n"+
		"   (3) map filenames to sample IDs (i.e. -mapFiles ("+(map?"":"not the ")+"default))\n"+
		"   (4) output file for mappings (i.e. out="+mapOutput+" (default))\n"+
//		" OR:\n"+
//		"   (1) parse AB lookup (i.e. --parseAB (not the default))\n"+
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
//			} else if (args[i].startsWith("-parseAB")) {
//				parseABlookup = true;
//				numArgs--;
			} else if (args[i].startsWith("-parseAlleleLookup")) {
				parseAlleleLookup = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			return;
		}

//		proj = null;
		proj = new Project(filename, false);
//
//		if (!proj.getDir(Project.SOURCE_DIRECTORY).equals("")&&!new File(proj.getDir(Project.SOURCE_DIRECTORY)).exists()) {
//			System.err.println("Error - the project source location is invalid: "+proj.getDir(Project.SOURCE_DIRECTORY));
//			return;
//		}

//		fileToConvert = "D:\\LOAD\\paper\\second submission\\genotype_counts\\top_hits.logistic.4pc_genotypes.txt";
//		fileToConvert = "D:\\LOAD\\paper\\second submission\\genotype_counts\\top_hits.logistic.4pc_full_apoe_genotypes.txt";
//		fileToConvert = "D:\\LOAD\\paper\\second submission\\genotype_counts\\top_hits.logistic.full_apoe_genotypes.txt";
//		fileToConvert = "D:\\LOAD\\paper\\second submission\\genotype_counts\\top_hits.logistic.nocovars_genotypes.txt";
		
//		lookupFile = "D:\\LOAD\\alleleLookup.txt";
		
		try {
//			if (parseABlookup) {
//				ABLookup.parseABlookup(proj);
//			} else 
			if (map) {
				mapFilenamesToSamples(proj, mapOutput);
			} else if (parseAlleleLookup) {
				parseAlleleLookup(proj);
			} else {
				createFiles(proj, numThreads);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
