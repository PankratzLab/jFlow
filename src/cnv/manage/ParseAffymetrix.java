package cnv.manage;


import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;

import cnv.filesys.Sample;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.SampleList;
import cnv.manage.NewParseIllumina.ParseConstants;
import common.*;

public class ParseAffymetrix implements Runnable {
	private Project proj;
	private String[] files;
	private String[] markerNames;
	private int[] keysKeys;
	private long fingerprint;
	private char[][] abLookup;
	private Hashtable<String,String> fixes;
	private long timeBegan;
	private int threadId;

	public ParseAffymetrix(Project proj, String[] files, String[] markerNames, int[] keysKeys, char[][] abLookup, long fingerprint, Hashtable<String,String> fixes, long timeBegan) {
		this(proj, files, markerNames, keysKeys, abLookup, fingerprint, fixes, timeBegan, -1);
	}

	public ParseAffymetrix(Project proj, String[] files, String[] markerNames, int[] keysKeys, char[][] abLookup, long fingerprint, Hashtable<String,String> fixes, long timeBegan, int threadId) {
		this.proj = proj;
		this.files = files;
		this.markerNames = markerNames;
		this.keysKeys = keysKeys;
		this.abLookup = abLookup;
		this.fingerprint = fingerprint;
		this.fixes = fixes;
		this.timeBegan = timeBegan;
		this.threadId = threadId;

		//TODO This seems not to be necessarily, because there is a check of "samples/.sampRAF".
//		if (Files.list(proj.getDir(proj.MARKER_DATA_DIRECTORY, true), MarkerData.MARKER_DATA_FILE_EXTENSION, proj.getJarStatus()).length>0) {
//			System.err.println("Error - Refusing to create new SampleList until the plots directory is either deleted or emptied; altering the SampleList will invalidate those files");
//			System.exit(1);
//		}

		//TODO There is a more complex scenario: what happens if you just want to add some new samples?
//		TransposeData.backupOlderFiles(proj.getDir(proj.SAMPLE_DIRECTORY, true), new String [] {"outliers.ser", ".sampRAF"}, true);
//		TransposeData.deleteOlderRafs(proj.getDir(proj.SAMPLE_DIRECTORY, true), new String[] {"outliers"}, new String[] {".ser"}, true);
	}

	public void run() {
		BufferedReader reader;
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
		String delimiter;
		String filename;
		Hashtable<String, Float> allOutliers;
		String sourceExtension;
		byte[] genos;
		
//		try {
//			PrintWriter writer = new PrintWriter(new FileWriter(files[0]+"_list.xln"));
//			for (int j = 0; j<files.length; j++) {
//				writer.println(files[j]);
//         }
//	        writer.close();
//	        return;
//     } catch (Exception e) {
//	        System.err.println("Error writing to "+files[0]+"_list.xln");
//	        e.printStackTrace();
//     }

		idHeader = proj.getProperty(proj.ID_HEADER);
		sourceExtension = proj.getProperty(proj.SOURCE_FILENAME_EXTENSION);
		delimiter = proj.getSourceFileDelimiter();
		allOutliers = new Hashtable<String, Float>();
		try {
			for (int i = 0; i<files.length; i++) {
				if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+ParseConstants.CANCEL_OPTION_FILE).exists()) {
					return;
				}
				try {
					System.out.println(ext.getTime()+"\t"+(i+1)+" of "+files.length + " -- " + files[i]);
//					reader = new BufferedReader(new FileReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[i]));
					reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[i]);
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready()&&(ext.indexFactors(ParseConstants.SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || (!idHeader.equals(ParseConstants.FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line)==-1)));

					dataIndices = ext.indexFactors(Sample.DATA_FIELDS, line, false, true, false, false);
					genotypeIndices = ext.indexFactors(Sample.GENOTYPE_FIELDS, line, false, true, false, false);
					if (idHeader.equals(ParseConstants.FILENAME_AS_ID_OPTION)) {
						sampIndex = -7;
					} else {
						sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];
					}
					snpIndex = ext.indexFactors(ParseConstants.SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
					
					if (dataIndices[3] == -1 || dataIndices[4] == -1) {
						System.err.println("Error - File format not consistent! At the very least the files need to contain "+Array.toStr(Sample.DATA_FIELDS[3], "/")+" and "+Array.toStr(Sample.DATA_FIELDS[4], "/"));
						return;
					}
					if (genotypeIndices[4] == -1 && (genotypeIndices[0] == -1 || genotypeIndices[1] == -1)) {
						System.err.println("Error - File format not consistent! The files need to contain "+Array.toStr(Sample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(Sample.GENOTYPE_FIELDS[1], "/") +" or "+Sample.GENOTYPE_FIELDS[4] +" (for single token calls)");
						return;
					}
					if (genotypeIndices[5] == -1 && (genotypeIndices[2] == -1 || genotypeIndices[3] == -1) && abLookup == null) {
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
					genotypes = new byte[2][];
					genotypes[0] = Array.byteArray(markerNames.length, (byte)0);	// used to be initialized to Byte.MIN_VALUE when AB genotypes && abLookup were both absent
					if (!ignoreAB) {
						genotypes[1] = Array.byteArray(markerNames.length, (byte)-1);	// used to be initialized to Byte.MIN_VALUE when AB genotypes && abLookup were both absent
					}
					
					count = 0;
//					parseAtAt = Boolean.parseBoolean(proj.getProperty(proj.PARSE_AT_AT_SYMBOL));
					parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
					while (reader.ready()) {
						line = reader.readLine().split(delimiter);
						if (idHeader.equals(ParseConstants.FILENAME_AS_ID_OPTION)) {
							trav = files[i].substring(0, files[i].indexOf(sourceExtension));
						} else {
							if (parseAtAt&&line[sampIndex].indexOf("@")==-1) {
								System.err.println("Error - "+idHeader+" '"+line[sampIndex]+"' did not contain an @ sample");
								parseAtAt = false;
							}
							trav = parseAtAt?line[sampIndex].substring(0, line[sampIndex].indexOf("@")):line[sampIndex];
						}

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
						for (int j = 0; j<Sample.DATA_FIELDS.length; j++) {
							try {
								if (dataIndices[j] != -1) {
									if (line[dataIndices[j]].equals("")) {
										data[j][key] = Float.NaN;
									} else {
										data[j][key] = Float.parseFloat(line[dataIndices[j]]);
									}
								}
							} catch (NumberFormatException nfe) {
								System.err.println("Error - failed to parse '"+line[dataIndices[j]]+"' into a valid "+Array.toStr(Sample.DATA_FIELDS[j], "/"));
								return;
							}
						}
						
						genos = ParseConstants.parseGenotypes(line, genotypeIndices, ignoreAB, abLookup, count, sampleName, markerNames[count], files[i]);
						genotypes[0][key] = genos[0];
						if (!ignoreAB) {
							genotypes[1][key] = genos[1];
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

					filename = ParseConstants.determineFilename(proj.SAMPLE_DIRECTORY.getValue(true, true), sampleName, timeBegan, proj.getLog());
					if (filename == null) {
						return;
					}

					samp = new Sample(sampleName, fingerprint, data, genotypes, true);
//					samp.serialize(proj.getDir(proj.SAMPLE_DIRECTORY, true) + trav + Sample.SAMPLE_DATA_FILE_EXTENSION);
//					samp.saveToRandomAccessFile(filename);
					samp.saveToRandomAccessFile(filename, allOutliers, sampleName); //TODO sampleIndex
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[i]+"\" not found in current directory");
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[i]+"\"");
					return;
				}
			}

			if (allOutliers.size()>0) {
				if (threadId >= 0) {
					if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + threadId + ".ser").exists()) {
						System.err.println("Error - the following file already exists: " + proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + threadId + ".ser");
						System.exit(1);
					} else {
						Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + threadId + ".ser");
					}
				} else {
					if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers0.ser").exists()) {
						System.err.println("Error - the following file already exists: " + proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers0.ser");
						System.exit(1);
					} else {
						Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers0.ser");
					}
				}
			}

			SampleList.generateSampleList(proj).writeToTextFile(proj.PROJECT_DIRECTORY.getValue()+"ListOfSamples.txt");

			System.out.println(ext.getTime()+"\tfinished");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	@SuppressWarnings("unchecked")
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
        int lineCount;
		boolean abLookupRequired;
		char[][] lookup;
		Hashtable<String, String> fixes;
		String idHeader, delimiter;
		String temp;
        String testline;
		int[][] delimiterCounts;
		boolean done;
		long timeBegan;
		boolean parseAtAt;
		int sampIndex;
		String sampleName;
		String[] overwriteOptions;
		int response;
		String[] filesToDelete;
		boolean complete;
		Hashtable<String, Float> allOutliers;

		timeBegan = new Date().getTime();
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+ParseConstants.OVERWRITE_OPTION_FILE).delete();
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+ParseConstants.HOLD_OPTION_FILE).delete();
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+ParseConstants.CANCEL_OPTION_FILE).delete();

		if (!proj.SOURCE_DIRECTORY.getValue(false, true).equals("")&&!new File(proj.SOURCE_DIRECTORY.getValue(false, true)).exists()) {
			System.err.println("Error - the Project source location is invalid: "+proj.SOURCE_DIRECTORY.getValue(false, true));
			return;
		}
     
		if (!new File(proj.MARKER_POSITION_FILENAME.getValue(false, false)).exists()) {
			System.err.println("Error - missing markerPositions: "+proj.MARKER_POSITION_FILENAME.getValue(false, false));
			return;
		}

		delimiter = proj.getSourceFileDelimiter();
		idHeader = proj.getProperty(proj.ID_HEADER);
		System.out.println(ext.getTime()+"\tSearching for "+proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)+" files in: "+proj.SOURCE_DIRECTORY.getValue(false, true));
		files = Files.list(proj.SOURCE_DIRECTORY.getValue(false, true), proj.getProperty(proj.SOURCE_FILENAME_EXTENSION), false);
		
		System.out.println("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" with a "+proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)+" extension");
		for (int i = 0; i < files.length; i++) {
			if (files[i].equals("Sample_Map.csv") || files[i].equals("SNP_Map.csv")) {
				files = Array.removeFromArray(files, i);
				i--;
			}
		}
		
		if (files.length == 0) {
			System.err.println("Error - no files to parse");
			return;
		}

		abLookupRequired = false;
		System.out.println("\t\tFound "+files.length+" file"+(files.length==1?"":"s")+" to parse");
		fixes = new Hashtable<String,String>();
		if (new File(proj.PROJECT_DIRECTORY.getValue()+"fixes.dat").exists()) {
			System.out.println("Also found a 'fixes.dat' file in the project directory, which will be used to rename samples");
			fixes = HashVec.loadFileToHashString(proj.PROJECT_DIRECTORY.getValue()+"fixes.dat", false);
		} else {
			System.out.println("Did not find a 'fixes.dat' file; assuming you don't want to rename any IDs");
		}

		try {
			reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]);
			System.out.println("Found appropriate reader for: "+proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]);
			count = 0;
			do {
				line = reader.readLine().trim().split(delimiter, -1);
				count++;
			} while (reader.ready()&&(ext.indexFactors(ParseConstants.SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || (!idHeader.equals(ParseConstants.FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line)==-1)));
			
			if (!reader.ready()) {
				System.err.println("Error - reached the end of the file without finding a line with the following tokens: "+Array.toStr(ParseConstants.SNP_HEADER_OPTIONS[0]));
				System.err.println("      - perhaps the delimiter is set incorrectly? Determing most stable delimiter...");

				reader.close();
				reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]);
				delimiterCounts = new int[ParseConstants.DELIMITERS.length][count];
				for (int i = 0; i < count; i++) {
					temp = reader.readLine();
					for (int j = 0; j < ParseConstants.DELIMITERS.length; j++) {
						delimiterCounts[j][i] = ext.countInstancesOf(temp, ParseConstants.DELIMITERS[j]);
					}
				}
				delimiter = null;
				for (int j = 0; j < ParseConstants.DELIMITERS.length; j++) {
					if (Array.quantWithExtremeForTie(delimiterCounts[j], 0.5) > 4 && Array.quantWithExtremeForTie(delimiterCounts[j], 0.9) - Array.quantWithExtremeForTie(delimiterCounts[j], 0.1) == 0) {
						if (delimiter == null) {
							delimiter = ParseConstants.DELIMITERS[j];
						} else {
							proj.message("Could not auto-detect the delimiter used in the Final Reports file: could be '"+delimiter+"' or '"+ParseConstants.DELIMITERS[j]+"'");
							return;
						}
					}
				}
				reader.close();
				
				if (delimiter == null) {
					proj.message("Failed to auto-detect the delimiter used in the Final Reports file; exitting");
					return;
				}
				
				System.err.println("      - determined delimiter to be '"+delimiter+"'");

				reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]);
				do {
					line = reader.readLine().trim().split(delimiter, -1);
				} while (reader.ready()&&(ext.indexFactors(ParseConstants.SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || (!idHeader.equals(ParseConstants.FILENAME_AS_ID_OPTION) && ext.indexOfStr(idHeader, line)==-1)));
				System.out.println(1);
			}

			// check immediately to make sure these fields are valid
			indices = ext.indexFactors(Sample.DATA_FIELDS, line, false, true, true, false); // dataIndices
			if (indices[3] == -1 || indices[4] == -1) {
				System.err.println("Error - at the very least the files need to contain "+Array.toStr(Sample.DATA_FIELDS[3], "/")+" and "+Array.toStr(Sample.DATA_FIELDS[4], "/"));
				System.err.println("      - failed to see that in "+files[0]);
				System.err.println(Array.toStr(line));
				
				return;
			}
			indices = ext.indexFactors(Sample.GENOTYPE_FIELDS, line, false, true, true, false); // genotypeIndices
			if (indices[4] == -1 && (indices[0] == -1 || indices[1] == -1)) {
				System.err.println("Error - the files need to contain "+Array.toStr(Sample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(Sample.GENOTYPE_FIELDS[1], "/") +" or "+Sample.GENOTYPE_FIELDS[4] +" (for single token calls)");
				return;
			}
			if (indices[5] == -1 && (indices[2] == -1 || indices[3] == -1)) {
				abLookupRequired = true;
			}

			snpIndex = ext.indexFactors(ParseConstants.SNP_HEADER_OPTIONS, line, false, true, true, true)[0];

//			parseAtAt = Boolean.parseBoolean(proj.getProperty(proj.PARSE_AT_AT_SYMBOL));
			parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
			idHeader = proj.getProperty(proj.ID_HEADER);
			if (idHeader.equals(ParseConstants.FILENAME_AS_ID_OPTION)) {
				sampleName = files[0].substring(0, files[0].indexOf(proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)));
			} else {
				sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];
				reader.mark(1000);
				line = reader.readLine().split(delimiter);
				if (parseAtAt&&line[sampIndex].indexOf("@")==-1) {
					System.err.println("Error - "+idHeader+" '"+line[sampIndex]+"' did not contain an @ sample");
					parseAtAt = false;
				}
				sampleName = parseAtAt?line[sampIndex].substring(0, line[sampIndex].indexOf("@")):line[sampIndex];
				reader.reset();
			}
			
			if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+sampleName+Sample.SAMPLE_DATA_FILE_EXTENSION).exists()) {

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
					, "What to do?", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null, overwriteOptions, overwriteOptions[2]);

				switch (response) {
				case -1:
					break;
				case 0:
					filesToDelete = Files.list(proj.SAMPLE_DIRECTORY.getValue(false, true), Sample.SAMPLE_DATA_FILE_EXTENSION, false);
					for (int i = 0; i < filesToDelete.length; i++) {
						new File(proj.SAMPLE_DIRECTORY.getValue(false, true) + filesToDelete[i]).delete();
					}
					new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser").delete();
					break;
				case 1:
					// keep "outlier.ser"
					break;
				case 2:
					return;
				default:
					proj.message("Should be impossible to obtain this message ("+response+")");
					break;
				}
			}
			TransposeData.deleteOlderRafs(proj.SAMPLE_DIRECTORY.getValue(true, true), new String[] {"outliers"}, new String[] {".ser"}, true, new String[] {"outliers.ser"});

//			if (Boolean.parseBoolean(proj.getProperty(proj.LONG_FORMAT))) {
			if (proj.getProperty(proj.LONG_FORMAT)) {
				reader.close();
				createFilesFromLongFormat(proj, files, idHeader, fixes, delimiter, abLookupRequired, timeBegan);
				return;
			}

			done = false;
//			count = 1;
			alNames = new ArrayList<String>(500000);
			while (reader.ready() && !done) {
				trav = reader.readLine().trim().split(delimiter)[snpIndex];
				if (trav.equals("") || trav.equals("0")) {
					trav = "Blank"+count;
					count++;
				}
//				if (markerNameHash.containsKey(trav)) {
//					done = true;
//				} else {
//					markerNameHash.put(trav, new Integer(markerNameHash.size()));
//				}
				alNames.add(trav);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]+"\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.SOURCE_DIRECTORY.getValue(false, true)+files[0]+"\"");
			return;
		}

//		new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)).mkdirs();
//		new File(proj.getDir(proj.IND_DIRECTORY)).mkdirs();
//		new File(proj.getDir(proj.DATA_DIRECTORY)).mkdirs();

//		markerNames = Array.toStringArray(markerNameHash);
		markerNames = Array.toStringArray(alNames);
//		keys = Markers.orderMarkers(markerNames, proj.getFilename(proj.MARKER_POSITION_FILENAME), proj.getFilename(proj.MARKERSET_FILENAME, true, true), proj.getLog());
		keys = Markers.orderMarkers(markerNames, proj.MARKER_POSITION_FILENAME.getValue(), proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());
		if (keys == null) {
			return;
		}
		keysKeys = Sort.quicksort(keys); // very important
		fingerprint = proj.getMarkerSet().getFingerprint();
		System.out.println("There are "+markerNames.length+" markers being processed (fingerprint: "+fingerprint+")");

		lookup = ParseConstants.getABLookup(abLookupRequired, markerNames, proj);

		fileCabinet = new Vector<Vector<String>>();
		for (int i = 0; i<numThreads; i++) {
			fileCabinet.add(new Vector<String>());
		}
		for (int i = 0; i<files.length; i++) {
			fileCabinet.elementAt(i%numThreads).add(files[i]);
		}
		threads = new Thread[numThreads];
		for (int i = 0; i<numThreads; i++) {
			threads[i] = new Thread(new ParseAffymetrix(proj, fileCabinet.elementAt(i).toArray(new String[fileCabinet.elementAt(i).size()]), markerNames, keysKeys, lookup, fingerprint, fixes, timeBegan, i));
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
		
		allOutliers = new Hashtable<String, Float>();
		for (int i = 0; i<numThreads; i++) {
			if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser").exists()) {
				allOutliers.putAll((Hashtable<String, Float>) Files.readSerial(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser"));
				new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser").delete();
			}
		}
		if (allOutliers.size()>0) {
			Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
		}
		
	}

	public static void createFilesFromLongFormat(Project proj, String[] files, String idHeader, Hashtable<String, String> fixes, String delimiter, boolean abLookupRequired, long timeBegan) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, markerNames, list;
		long fingerprint;
		MarkerSet markerSet;
		Hashtable<String, Integer> markerIndices;
		int count;
		char[][] abLookup;
		String filename;
		Hashtable<String, Float> allOutliers;
        Hashtable<String, String> renamedIDsHash;
        Logger log;

        log = proj.getLog();
		System.out.println("Parsing files using the Long Format algorithm");
     
		Markers.orderMarkers(null, proj.MARKER_POSITION_FILENAME.getValue(), proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		fingerprint = proj.getMarkerSet().getFingerprint();
		
		abLookup = ParseConstants.getABLookup(abLookupRequired, markerNames, proj);
		
		markerIndices = new Hashtable<String, Integer>();
		for (int i = 0; i < markerNames.length; i++) {
			markerIndices.put(markerNames[i], new Integer(i));
		}
		
		System.out.println("There were "+markerNames.length+" markers present in '"+proj.MARKERSET_FILENAME.getValue(true, true)+"' that will be processed from the source files (fingerprint: "+fingerprint+")");
		
		int snpIndex, sampIndex, key;
		String trav;
		int[] dataIndices, genotypeIndices;
		boolean parseAtAt;

		Sample samp;
		String sampleName;
		float[][] data;
		byte[][] genotypes;
		boolean ignoreAB;
		CountHash countHash, dupHash;
		boolean done;
		byte[] genos;
		
		count = 0;
		countHash = new CountHash();
		dupHash = new CountHash();
		allOutliers = new Hashtable<String, Float>();
		try {
			for (int i = 0; i<files.length; i++) {
				try {
					System.out.println(ext.getTime()+"\t"+(i+1)+" of "+files.length+" ("+files[i]+")");
//					reader = new BufferedReader(new FileReader(proj.getDir(proj.SOURCE_DIRECTORY)+files[i]));
					reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)+files[i]);
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready()&&(ext.indexFactors(ParseConstants.SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 || ext.indexOfStr(idHeader, line)==-1));

					System.err.println("Searching: "+Array.toStr(line));
					dataIndices = ext.indexFactors(Sample.DATA_FIELDS, line, false, true, false, false);
					genotypeIndices = ext.indexFactors(Sample.GENOTYPE_FIELDS, line, false, true, false, false);
					sampIndex = ext.indexFactors(new String[] {idHeader}, line, false, true)[0];
					snpIndex = ext.indexFactors(ParseConstants.SNP_HEADER_OPTIONS, line, false, true, false, true)[0];
					
					if (dataIndices[3] == -1 || dataIndices[4] == -1) {
						System.err.println("Error - File format not consistent! At the very least the files need to contain "+Array.toStr(Sample.DATA_FIELDS[3], "/")+" and "+Array.toStr(Sample.DATA_FIELDS[4], "/"));
						return;
					}
					if (genotypeIndices[4] == -1 && (genotypeIndices[0] == -1 || genotypeIndices[1] == -1)) {
						System.err.println("Error - File format not consistent! The files need to contain "+Array.toStr(Sample.GENOTYPE_FIELDS[0], "/")+" and "+Array.toStr(Sample.GENOTYPE_FIELDS[1], "/") +" or "+Sample.GENOTYPE_FIELDS[4] +" (for single token calls)");
						return;
					}
					if (genotypeIndices[5] == -1 && (genotypeIndices[2] == -1 || genotypeIndices[3] == -1) && abLookup == null) {
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
								
								filename = ParseConstants.determineFilename(proj.SAMPLE_DIRECTORY.getValue(true, true), sampleName, timeBegan, proj.getLog());
								if (filename == null) {
									return;
								}
								
								samp = new Sample(sampleName, fingerprint, data, genotypes, true);
								samp.saveToRandomAccessFile(filename, allOutliers, sampleName);
							}
							if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + trav + Sample.SAMPLE_DATA_FILE_EXTENSION).exists()) {
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

						if (!done) {
							countHash.add(line[snpIndex]);							
							if (markerIndices.containsKey(line[snpIndex])) {
								key = markerIndices.get(line[snpIndex]);

								for (int j = 0; j<Sample.DATA_FIELDS.length; j++) {
									try {
										if (dataIndices[j] != -1) {
											if (!(data[j][key]+"").equals("Infinity")) {
												dupHash.add(line[snpIndex]);
												System.err.println("Sample "+trav+" already has data for marker "+line[snpIndex]+" (Was the parsing restarted? Delete the old directories first)");
											}
											data[j][key] = Float.parseFloat(line[dataIndices[j]]);
										}
									} catch (NumberFormatException nfe) {
										System.err.println("Error - failed to parse"+line[dataIndices[j]]+" into a valid "+Sample.DATA_FIELDS[j]);
										return;
									}
								}
								
								genos = ParseConstants.parseGenotypes(line, genotypeIndices, ignoreAB, abLookup, count, sampleName, markerNames[count], files[i]);
								genotypes[0][key] = genos[0];
								if (!ignoreAB) {
									genotypes[1][key] = genos[1];
								}
							}
						}
						count++;
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

			if (allOutliers.size()>0) {
				if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser").exists()) {
					System.err.println("Error - the following file already exists: " + proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
					System.exit(1);
				} else {
					Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
				}
			}

			System.out.println(ext.getTime()+"\tfinished");
		} catch (Exception e) {
			e.printStackTrace();
		}

     
		System.out.println("Parsed "+count+" sample(s)");
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
			System.err.println("Error writing to " + "ListOfMarkers.txt");
			e.printStackTrace();
		}

		new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+ParseConstants.OVERWRITE_OPTION_FILE).delete();
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+ParseConstants.HOLD_OPTION_FILE).delete();
		new File(proj.SAMPLE_DIRECTORY.getValue(true, true)+ParseConstants.CANCEL_OPTION_FILE).delete();
	}

	//	public static Hashtable filenamesUsedByExistingSampleFiless(Project proj) {
//		Hashtable <String, String> result = null;
//		String[] sampNames;
//		File[] files;
//		String[] fileNames;
//
//		files = new File(proj.getDir(proj.SAMPLE_DIRECTORY, true)).listFiles(new FilenameFilter() {
//			public boolean accept(File file, String filename) {
//				return filename.endsWith(Sample.SAMPLE_DATA_FILE_EXTENSION);
//			}
//		});
//
//		if (files != null) {
//			sampNames = proj.getSamples();
//			fileNames = new String[files.length];
//			result = new Hashtable<String, String>();
//			for (int i = 0; i < files.length; i++) {
//				fileNames[i] = files[i].getName();
//			}
//			for (int i = 0; i < sampNames.length; i++) {
//				for (int j = 0; j < fileNames.length; j++) {
//					if (sampNames[i].equals(fileNames[j])) {
//						result.put(sampNames[i], files[j].getName());
//					}
//				}
//			}
//		}
//
//		return result;
//	}
//
	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj;
		String filename = null;
		boolean map = false;
		int numThreads = 1;
//		boolean parseABlookup = false;
		boolean parseAlleleLookupFromFinalReports = false;
		String mapOutput = "filenamesMappedToSamples.txt";

		String usage = "\n"+
		"cnv.manage.ParseAffymetrix requires 0-1 arguments\n"+
		"   (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
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
				parseAlleleLookupFromFinalReports = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			return;
		}

		filename = "C:/workspace/Genvisis/projects/COGA_exome.properties";
//		proj = null;
		proj = new Project(filename, false);
//
//		if (!proj.getDir(proj.SOURCE_DIRECTORY).equals("")&&!new File(proj.getDir(proj.SOURCE_DIRECTORY)).exists()) {
//			System.err.println("Error - the project source location is invalid: "+proj.getDir(proj.SOURCE_DIRECTORY));
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
			    ParseConstants.mapFilenamesToSamples(proj, mapOutput);
			} else if (parseAlleleLookupFromFinalReports) {
			    ParseConstants.parseAlleleLookupFromFinalReports(proj);
			} else {
				createFiles(proj, numThreads);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
