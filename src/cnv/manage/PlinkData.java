package cnv.manage;

import java.io.*;
import java.util.*;
import java.text.SimpleDateFormat;

import common.*;
import cnv.filesys.*;
import cnv.var.SampleData;

public class PlinkData {
	public static final String FAM_DELIMITER = " ";
	public static final String BIM_DELIMITER = "\t";
	public static final String MAP_DELIMITER = "\t";

	/**
	 * Convert a Plink .ped data set to Plink .bed data set.
	 * 
	 * Note: This is the same feature as the Plink command "plink --file mydata --make-bed"
	 * 
	 * @param plinkDirAndFilenameRoot
	 * @param pedFileFilenameRoot
	 * @param isSnpMajor
	 */
	public static void convertPedSetToBedSet(String plinkDirAndFilenameRoot, boolean isSnpMajor, Logger log) {
//		RandomAccessFile outBed;
//		RandomAccessFile outFam;
		BufferedOutputStream outBed;
		PrintWriter outFamOrBim;
		Scanner inFile;
		int numMarkers;
		String[] line;
		byte genotype;
		byte[] genotypeByteStream = null;
		String[][] genotypeLetters = null;
		int index;
		int index2;
		byte[] alleles;
		String temp;

		if (log == null) {
			log = new Logger();
		}

		if (new File(plinkDirAndFilenameRoot + ".bed").exists() || new File(plinkDirAndFilenameRoot + ".bim").exists() || new File(plinkDirAndFilenameRoot + ".fam").exists()) {
			log.reportError("System abort. PLINK binary file set \"" + plinkDirAndFilenameRoot + "\" .bed/.bim/.map already exist. Please remove the file(s).");
			return;
//			log.report("Found existing PLINK .bed file set in out file directory. Deleting these files.");
//			new File(plinkDirAndFilenameRoot + ".bed").delete();
//			new File(plinkDirAndFilenameRoot + ".fam").delete();
//			new File(plinkDirAndFilenameRoot + ".bim").delete();
		}

		try {
			genotype = (byte) 0;
			numMarkers = -1;
			alleles = new byte[2];
			inFile = new Scanner(new FileInputStream(plinkDirAndFilenameRoot + ".ped"));
//			outFam = new RandomAccessFile(plinkDataDir + plinkPedFileFilenameRoot + ".fam", "rw");
			outFamOrBim = new PrintWriter(new FileOutputStream(plinkDirAndFilenameRoot + ".fam"));
//			outBed = new RandomAccessFile(plinkDataDir + plinkPedFileFilenameRoot + ".bed", "rw");
//			outBed.write(new byte[] {(byte) 108, (byte) 27, (byte) 1});
			outBed = new BufferedOutputStream(new FileOutputStream(plinkDirAndFilenameRoot + ".bed"));
			outBed.write(new byte[] {(byte) 108, (byte) 27, (byte) 1});

			while (inFile.hasNext()) {
				line = inFile.nextLine().split("\\s+");
				if (numMarkers == -1) {
					numMarkers = (line.length - 6) / 2;
					genotypeLetters = new String[numMarkers][2];
					for (int i = 0; i < genotypeLetters.length; i++) {
						genotypeLetters[i][0] = "0";
						genotypeLetters[i][1] = "0";
					}
					genotypeByteStream = new byte[(int) Math.ceil((double) numMarkers / 4)];
				}
				outFamOrBim.println(line[0] + FAM_DELIMITER + line[1] + FAM_DELIMITER + line[2] + FAM_DELIMITER + line[3] + FAM_DELIMITER + line[4] + FAM_DELIMITER + line[5]);

				index = 0;
				index2 = 0;
				genotypeByteStream[genotypeByteStream.length - 1] = (byte) 0;
				for (int i = 0; i < numMarkers; i++) {
					if (i == 30 && line[0].equals("4")) {
						System.out.println("");
					}
					for (int j = 0; j < 2; j++) {
						temp = line[7 + i + i - j];
						if (temp.equals("0")) {
							alleles[j] = (byte) -1;
						} else if (genotypeLetters[i][1].equals("0")) {
							genotypeLetters[i][1] = temp;
							alleles[j] = (byte) 0;
						} else if (genotypeLetters[i][1].equals(temp)) {
							alleles[j] = (byte) 0;
						} else {
							alleles[j] = (byte) 1;
							if (genotypeLetters[i][0].equals("0")) {
								genotypeLetters[i][0] = temp;
							}
						}
					}

					if (alleles[0] == 0 && alleles[1] == 0) {
						genotype = (byte) (genotype | 0x00 << (index * 2));
					} else if (alleles[0] == 0 && alleles[1] == 1) {
						genotype = (byte) (genotype | 0x01 << (index * 2));
					} else if (alleles[0] == 1 && alleles[1] == 0) {
						genotype = (byte) (genotype | 0x01 << (index * 2));
					} else if (alleles[0] == 1 && alleles[1] == 1) {
						genotype = (byte) (genotype | 0x03 << (index * 2));
					} else {
						genotype = (byte) (genotype | 0x02 << (index * 2));
					}
					index ++;
					
					if (index == 4) {
						index = 0;
						genotypeByteStream[index2] = genotype;
						index2 ++;
						genotype = 0;
					}
				}
				outBed.write(genotypeByteStream);
			}
			outBed.close();
			outFamOrBim.close();

			index = 0;
			inFile = new Scanner(new FileInputStream(plinkDirAndFilenameRoot + ".map"));
			outFamOrBim = new PrintWriter(new FileOutputStream(plinkDirAndFilenameRoot + ".bim"));
			while (inFile.hasNext()) {
				temp = inFile.nextLine();
				outFamOrBim.println(temp + BIM_DELIMITER + genotypeLetters[index][0] + BIM_DELIMITER + genotypeLetters[index][1]);
				index ++;
			}
			inFile.close();
			outFamOrBim.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Convert a PLINK .bed data set to PLINK .ped data set.
	 * Note: This is the same feature as the PLINK command "plink --bfile mydata --recode --out mynewdata"
	 * 
	 * @param plinkDirAndFilenameRoot
	 * @param log
	 */
	public static void convertBedSetToPedSet(String plinkDirAndFilenameRoot, Logger log) {
		if (log == null) {
			log = new Logger();
		}

		if (new File(plinkDirAndFilenameRoot + ".ped").exists() || new File(plinkDirAndFilenameRoot + ".map").exists()) {
			log.reportError("System abort. PLINK binary file set \"" + plinkDirAndFilenameRoot + "\" .ped/.map already exist. Please remove the file(s).");
			return;
//			log.report("Found existing PLINK .ped file set in out file directory. Deleting these files.");
//			new File(plinkDirAndFilenameRoot + ".ped").delete();
//			new File(plinkDirAndFilenameRoot + ".map").delete();
		}

		convertFamBedToPed(plinkDirAndFilenameRoot, convertBimToMap(plinkDirAndFilenameRoot), log);
	}
	
	/**
	 * Load a PLINK .bim file and return the list of alleles.
	 * @param plinkDirAndFilenameRoot
	 * @return
	 */
	public static char[][] loadBimAlleles(String plinkDirAndFilenameRoot) {
		Vector<char[]> allelesTmp;
		Scanner reader;
		String[] line;
		char[][] alleles;
		
		allelesTmp = new Vector<char[]>();
		try {
			reader = new Scanner(new FileInputStream(plinkDirAndFilenameRoot + ".bim"));
			reader.nextLine();
			while (reader.hasNext()) {
				line = reader.nextLine().split("\t");
				allelesTmp.add(new char[] {line[4].charAt(0), line[5].charAt(0)});
			}
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		alleles = new char[allelesTmp.size()][2];
		for (int i = 0; i < alleles.length; i++) {
			alleles[i] = allelesTmp.elementAt(i);
		}

		return alleles;
	}
	
	/**
	 * Convert a PLINK .bim file to PLINK .map file and return the list of alleles.
	 * This is normally used as part of PlinkData.convertBedToPed(...)
	 * @param plinkDirAndFilenameRoot
	 * @return
	 */
	public static char[][] convertBimToMap(String plinkDirAndFilenameRoot) {
		Vector<char[]> allelesTmp;
		Scanner reader;
		String[] line;
		PrintWriter writer;
		char[][] alleles;
		
		allelesTmp = new Vector<char[]>();
		try {
			reader = new Scanner(new FileInputStream(plinkDirAndFilenameRoot + ".bim"));
			reader.nextLine();
			
			writer = new PrintWriter(new FileOutputStream(plinkDirAndFilenameRoot + ".map"));
			while (reader.hasNext()) {
				line = reader.nextLine().split("\t");
				writer.println(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3]);
				allelesTmp.add(new char[] {line[4].charAt(0), line[5].charAt(0)});
			}
			reader.close();
			writer.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		alleles = new char[allelesTmp.size()][2];
		for (int i = 0; i < alleles.length; i++) {
			alleles[i] = allelesTmp.elementAt(i);
		}

		return alleles;
	}

	/**
	 * Convert PLINK .fam and .bed file to PLINK .ped file.
	 * This is normally used as part of PlinkData.convertBedSetToPedSet().
	 * @param plinkDirAndFilenameRoot
	 * @param alleles
	 */
	public static void convertFamBedToPed(String plinkDirAndFilenameRoot, char[][] alleles, Logger log) {
		Scanner famReader;
		String famLine;
		BufferedInputStream bedReader;
		byte[] isSnpMajor;
		byte[] bedByteStream;
		byte bedByte;
		int bedByteStreamIndex;
		int bytesPerSampleOrMarkerInBed;
		int bytesPerSampleInPed;
		int index;
		int iteration;
		PrintWriter printWriter;
		RandomAccessFile rafWriter;
		long seekLocation;
		Vector<String> fam;
		int nFams;
		byte[] famIDsLength;
		String newLine;

		if (log == null) {
			log = new Logger();
		}
		newLine = "\n";
		try {
			famReader = new Scanner(new FileInputStream(plinkDirAndFilenameRoot + ".fam"));
			famReader.nextLine();

			isSnpMajor = new byte[3];
			bedByteStreamIndex = 3;
			bedReader = new BufferedInputStream(new FileInputStream(plinkDirAndFilenameRoot + ".bed"));
			bedReader.read(isSnpMajor);

			if (isSnpMajor[2] == 0) {	// Individual Major
				bytesPerSampleOrMarkerInBed = (int) Math.ceil((float) alleles.length / 4);
				bedByteStream = new byte[bytesPerSampleOrMarkerInBed];
				printWriter = new PrintWriter(new FileOutputStream(plinkDirAndFilenameRoot + ".ped"));
				while (famReader.hasNext()) {
					famLine = famReader.nextLine();
					printWriter.print(famLine);
					bedReader.read(bedByteStream, bedByteStreamIndex, bytesPerSampleOrMarkerInBed);
					index = 0;
					for (int i = 0; i < bedByteStream.length; i++) {
						iteration = Math.min(alleles.length - index, 4);
						for (int j = 0; j < iteration; j++) {
							bedByte = (byte) ((bedByteStream[i] >> (j * 2)) & 0x03);
							if (bedByte == 0) {
								printWriter.print("\t" + alleles[index][0] + "\t" + alleles[index][0]);
	
							} else if (bedByte == 1) {
								printWriter.print("\t" + alleles[index][0] + "\t" + alleles[index][1]);
	
							} else if (bedByte == 2) {
								printWriter.print("\t\t");
	
							} else if (bedByte == 3) {
								printWriter.print("\t" + alleles[index][1] + "\t" + alleles[index][1]);
							}
							index ++;
						}
					}
					printWriter.println("");
					bedByteStreamIndex += bytesPerSampleOrMarkerInBed;
				}
				printWriter.close();

			} else {	// SNP Major
				bytesPerSampleInPed = alleles.length * 2;
				fam = loadFamOrBim(plinkDirAndFilenameRoot + ".fam", 0, -1, log);
				nFams = fam.size();
				famIDsLength = new byte[nFams];

				rafWriter = new RandomAccessFile(plinkDirAndFilenameRoot + ".ped", "w");
				rafWriter.write(fam.elementAt(0).getBytes());
				for (int i = 1; i < nFams; i++) {
					rafWriter.skipBytes(bytesPerSampleInPed);
					bedByteStream = (newLine + fam.elementAt(i)).getBytes();
					rafWriter.write(bedByteStream);
					famIDsLength[i] = (byte) bedByteStream.length;
				}
				rafWriter.skipBytes(bytesPerSampleInPed);
				rafWriter.write(newLine.getBytes());

				bytesPerSampleOrMarkerInBed = (int) Math.ceil((float) nFams / 4);
				bedByteStream = new byte[bytesPerSampleOrMarkerInBed];
				for (int i = 0; i < alleles.length; i++) {
					bedReader.read(bedByteStream, bedByteStreamIndex, bytesPerSampleOrMarkerInBed);
					index = 0;
					for (int j = 0; j < bedByteStream.length; j++) {
						iteration = Math.min(nFams - index, 4);
						for (int l = 0; l < iteration; l++) {
							seekLocation = famIDsLength[index] + i * 4;
							seekLocation += (bytesPerSampleInPed + famIDsLength[index]);
							rafWriter.seek(seekLocation);
							bedByte = (byte) ((bedByteStream[j] >> (l * 2)) & 0x03);
							if (bedByte == 0) {
								rafWriter.write(("\t" + alleles[index][0] + "\t" + alleles[index][0]).getBytes());
							} else if (bedByte == 1) {
								rafWriter.write(("\t" + alleles[index][0] + "\t" + alleles[index][1]).getBytes());
							} else if (bedByte == 2) {
								rafWriter.write(("\t \t ").getBytes());
							} else if (bedByte == 3) {
								rafWriter.write(("\t" + alleles[index][1] + "\t" + alleles[index][1]).getBytes());
							}
							index ++;
						}
					}
					bedByteStreamIndex += bytesPerSampleOrMarkerInBed;
				}
				rafWriter.close();
			}

			famReader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Load a PLINK .fam or .bim file into a String vector
	 * @param famOrBimFileFullPath full path of the .fam or .bim file
	 * @param indexOfStartSampleOrMarker index of the start sample, if this is a .fam file, or index of the start marker to load, if this is a .bim file
	 * @param numberOfSamplesOrMarkersToLoad number of samples to load, if this is a .fam file, or number of markers to load, if this is a .bim file.
	 * @param log
	 * @return
	 */
	public static Vector<String> loadFamOrBim(String famOrBimFileFullPath, int indexOfStartSampleOrMarker, int numberOfSamplesOrMarkersToLoad, Logger log) {
		Scanner reader;
		Vector<String> data;

		data = new Vector<String>();
		try {
			reader = new Scanner(new FileInputStream(famOrBimFileFullPath));

			while (indexOfStartSampleOrMarker > 0 && reader.hasNext()) {
				reader.nextLine();
				indexOfStartSampleOrMarker --;
			}

			if (numberOfSamplesOrMarkersToLoad <= 0) {
				while (reader.hasNext()) {
					data.add(reader.nextLine());
					numberOfSamplesOrMarkersToLoad --;
				}

			} else {
				while (numberOfSamplesOrMarkersToLoad > 0 && reader.hasNext()) {
					data.add(reader.nextLine());
					numberOfSamplesOrMarkersToLoad --;
				}
			}

			reader.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return null;
		}
		
		return data;
	}

	/**
	 * Load a PLINK .fam file into a String array.
	 * @param famFileFullPath
	 * @param indexOfStartSample
	 * @param numberOfSamplesToLoad
	 * @param log
	 * @return
	 */
	public static String[][] loadFamToArray(String famFileFullPath, int indexOfStartSample, int numberOfSamplesToLoad, Logger log) {
		Vector<String> fam;
		String[][] out;

		fam = loadFamOrBim(famFileFullPath, indexOfStartSample, numberOfSamplesToLoad, log);
		out = new String[fam.size()][];
		for (int i = 0; i < out.length; i++) {
			out[i] = fam.elementAt(i).split("\\s+");
		}
		
		return out;
	}


	/**
	 * Load a PLINK .bim file into a String vector
	 * @param famFileFullPath
	 * @param indexOfStartSample
	 * @param numberOfSamplesToLoad
	 * @param log
	 * @return
	 */
	public static String[][] loadBimToArray(String famFileFullPath, int indexOfStartSample, int numberOfSamplesToLoad, Logger log) {
		Vector<String> bim;
		String[][] out;

		bim = loadFamOrBim(famFileFullPath, indexOfStartSample, numberOfSamplesToLoad, log);
		out = new String[bim.size()][];
		for (int i = 0; i < out.length; i++) {
			out[i] = bim.elementAt(i).split("\\s+");
		}
		
		return out;
	}

//	public static void loadPedData(String plinkDataDir, String bedFileFilenameRoot, int[] indicesOfMarkersToLoad, int[] indicesOfSamplesToLoad) {
//		Scanner inFamOrBim;
//		String[] line;
//		Vector<Byte> chr;
//		Vector<String> markerName;
//		Vector<Double> morgans;
//		Vector<Integer> pos;
//		int nMarkers;
//		Vector<String> famId;
//		Vector<String> sampId;
//		Vector<String> paternalId;
//		Vector<String> maternalId;
//		Vector<Byte> sex;
//		Vector<String> phenotype;
//		byte[] genotypes;
//
//		chr = new Vector<Byte>();
//		markerName = new Vector<String>();
//		morgans = new Vector<Double>();
//		pos = new Vector<Integer>();
//		try {
//			while (inFamOrBim.hasNext()) {
//				line = inFamOrBim.nextLine().split("\t");
//				chr.add(Byte.parseByte(line[0]));	//TODO to verify chr 23, 24, 25, and 26
//				markerName.add(line[1]);
//				morgans.add(Double.parseDouble(line[2]));
//				pos.add(Integer.parseInt(line[3]));
//			}
//			inFamOrBim.close();
//			nMarkers = markerName.size();
//
//			famId = new Vector<String>();
//			sampId = new Vector<String>();
//			paternalId = new Vector<String>();
//			maternalId = new Vector<String>();
//			sex = new Vector<Byte>();
//			phenotype = new Vector<String>();
//			genotypes = new byte[nMarkers];
//			inFamOrBim = new Scanner(new FileInputStream(plinkDataDir + "plink.ped"));
//			while (inFamOrBim.hasNext()) {
//				line = inFamOrBim.nextLine().split("\t");
//				famId.add(line[0]);
//				sampId.add(line[1]);
//				paternalId.add(line[2]);
//				maternalId.add(line[3]);
//				sex.add(line[4]);
//				phenotype.add(line[5]);
//				for (int i = 0; i < genotypes.length; i++) {
//					if (line[i + 6].equalsIgnoreCase("0")) {
//						genotypes[i] = -1;
//					} else if (line[i + 6].equalsIgnoreCase("")) {
//						genotypes[i] = 0;
//					} else if (line[i + 6].equalsIgnoreCase("")) {
//						genotypes[i] = 1;
//					} else if (line[i + 6].equalsIgnoreCase("")) {
//						genotypes[i] = 2;
//					}
//				}
//			}
//			inFamOrBim.close();
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
//	}

//	public static void loadBedData(String plinkDataDir, String bedFileFilenameRoot, int[] indicesOfMarkersToLoad, int[] indicesOfSamplesToLoad) {
//		Scanner inFamOrBim;
//		String[] line;
//		Vector<String> sampId;
//		Vector<Byte> sex;
//		Vector<String> cluster;
//		BufferedInputStream inBed;
//
//		int nSamples;
//		int nMarkers;
//		int iterations;
//		int index;
//		byte[] genotypes;
//		byte buffer;
//		ClusterFilterCollection clusterFilterCollection;
//		MarkerData[] markerData;
//
//		sampId = new Vector<String>();
//		sex = new Vector<Byte>();
//		cluster = new Vector<String>();
//		inFamOrBim = new Scanner(new FileInputStream(plinkDataDir + bedFileFilenameRoot + ".bed"));
//		while (inFamOrBim.hasNext()) {
//			line = inFamOrBim.nextLine().split("\t");
//			sampId.add(line[0] + "\t" + line[1]);
//			sex.add(line[4].equals("1")? (byte) 0 : (byte) 1);
//			cluster.add(line[5]);
//		}
//
//		inBed = new BufferedInputStream(new FileInputStream(plinkDataDir + bedFileFilenameRoot + ".bed"));
//
//		if (markerList == null) {
//			nMarkers = proj.getSamples().length;
//			markerList = proj.getMarkerNames();
//		} else {
//			nMarkers = markerList.length;
//		}
//		nSamples = proj.getSamples().length;
//		clusterFilterCollection = proj.getClusterFilterCollection();
//		
//		try {
//			in = new RandomAccessFile(PlinkBinaryDataDir + "plink.bed", "rw");
//			outStream = new byte[1];
//			outStream[0] = (byte) 108;	// 0b01101100
//			outStream[1] = (byte) 27;	// 0b00011011
//			outStream[2] = (byte) 1;	//0b00000001
//	
//			iterations = (int) Math.ceil((double) nSamples / 4);
//			markerData = new MarkerDataLoader(proj, markerList, 0).markerData;
//			for (int i = 0; i < markerList.length; i++) {
//				genotypes = markerData[i].getAbGenotypesAfterFilters(clusterFilterCollection, markerList[i], (float) gcThreshold / 100);
//				index = 0;
//				for (int j = 0; j < iterations; j++) {
//					buffer = (byte) 0;
//					for (int k = 0; k < 4 && index < nSamples; k++) {
//						if (genotypes[index] == (byte) 0) {
//							buffer = (byte) (buffer | (byte) 2 << (2 * k));
//						} else if (genotypes[index] == (byte) 1) {
//							buffer = (byte) (buffer | (byte) 0 << (2 * k));
//						} else if (genotypes[index] == (byte) 2) {
//							buffer = (byte) (buffer | (byte) 1 << (2 * k));
//						} else if (genotypes[index] == (byte) 3) {
//							buffer = (byte) (buffer | (byte) 3 << (2 * k));
//						} else {
//							System.err.println("Unrecognized genotype: " + genotypes[index] + " for marker " + markerList[i] + ", sample " + index);
//							System.exit(1);
//						}
//						index ++;
//					}
//					outStream[j] = buffer;
//				}
//				in.write(outStream);
//			}
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}
	
	public static boolean saveGenvisisToPlinkPedSet(Project proj, String filenameRoot, String clusterFiltersFilename, String targetMarkersFilename) {
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
		Hashtable<Integer,Integer> invalidAbLookups;
		Logger log;

		log = proj.getLog();
		log.report(ext.getTime());
		hash = new Hashtable<String,String>();
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		
		String PROG_KEY = "PLINKEXPORT";
		proj.progressMonitor.beginIndeterminateTask(PROG_KEY, "Exporting marker data for PLINK analysis", ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		
		for (int i = 0; i<markerNames.length; i++) {
			if (hash.containsKey(markerNames[i])) {
				proj.message("Warning - duplicate marker name: "+markerNames[i]);
			}
			hash.put(markerNames[i], i+"");
		}
		
		targetMarkers = targetMarkersFilename;
		if (targetMarkers != null && new File(targetMarkers).exists()) {
			targets = HashVec.loadFileToStringArray(targetMarkers, false, false, new int[] {0}, false);
			indices = new int[targets.length];
			prob = false;
			for (int i = 0; i<targets.length; i++) {
				if (hash.containsKey(targets[i])) {
					indices[i] = Integer.parseInt(hash.get(targets[i]));
				} else {
					proj.message("Error - target marker '"+targets[i]+"' was not found in the MarkerSet");
					prob = true;
				}
			}
			if (prob) {
				return false;
			}
		} else {
			if (!targetMarkers.equals("")) {
				proj.message("FYI, since target markers file '"+targetMarkers+"' was not found, all markers will be exported to PLINK");
			}

			indices = Array.intArray(markerNames.length);
		}
		
		proj.progressMonitor.updateTask(PROG_KEY);
		proj.progressMonitor.beginDeterminateTask(PROG_KEY + "_MAPEXPORT", "Exporting marker data to .map file", indices.length, ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		try {
		    log.report(ext.getTime() + "]\tWriting " + filenameRoot + ".map");
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+filenameRoot+".map"));
			for (int i = 0; i<indices.length; i++) {
				writer.println(chrs[indices[i]]+" "+markerNames[indices[i]]+" 0 "+positions[indices[i]]);
				proj.progressMonitor.updateTask(PROG_KEY + "_MAPEXPORT");
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			proj.message("Error: failed to write to "+filenameRoot+".map (in use?)");
			proj.progressMonitor.endTask(PROG_KEY + "_MAPEXPORT");
			proj.progressMonitor.endTask(PROG_KEY);
			return false;
		} catch (IOException ioe) {
			proj.message("Error writing to "+filenameRoot+".map");
			proj.progressMonitor.endTask(PROG_KEY + "_MAPEXPORT");
			proj.progressMonitor.endTask(PROG_KEY);
			return false;
		}
		proj.progressMonitor.endTask(PROG_KEY + "_MAPEXPORT");
		
		proj.progressMonitor.updateTask(PROG_KEY);

		gcThreshold = proj.GC_THRESHOLD.getValue().floatValue();
		if (gcThreshold < 0) {
			proj.message("Error - GC_THRESHOLD must be greater than zero (not "+gcThreshold+")");
		}
		if (gcThreshold > 1) {
			proj.message("Error - GC_THRESHOLD must be less than one (not "+gcThreshold+")");
		}
		log.report("Using a GC threshold of "+gcThreshold+" (less than or equal to will be set to missing, greater than is kept)");
		
		clusterFilterCollection = null;
		if (clusterFiltersFilename != null) {
			clusterFiltersFilename = proj.getProperty(proj.PROJECT_DIRECTORY)+proj.getProperty(proj.DATA_DIRECTORY)+clusterFiltersFilename;
			if (Files.exists(clusterFiltersFilename, proj.JAR_STATUS.getValue())) {
				clusterFilterCollection = ClusterFilterCollection.load(clusterFiltersFilename, proj.JAR_STATUS.getValue());
			} else {
				proj.message("Error - cluster filter collection is not found at '"+clusterFiltersFilename+"'");
	            proj.progressMonitor.endTask(PROG_KEY);
				return false;
			}
			abLookup = new ABLookup(markerNames, proj.AB_LOOKUP_FILENAME.getValue(), true, true, proj.getLog()).getLookup();
			log.report("Using "+clusterFiltersFilename+" and "+proj.getProperty(proj.AB_LOOKUP_FILENAME)+" to call genotypes");
		} else {
			abLookup = null;
		}
		
		int exp = Files.countLines(proj.PEDIGREE_FILENAME.getValue(), 0);
		proj.progressMonitor.beginDeterminateTask(PROG_KEY + "_PEDEXPORT", "Exporting sample data to .ped file", exp, ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		
		try {
			reader = new BufferedReader(new FileReader(proj.PEDIGREE_FILENAME.getValue()));
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+filenameRoot+".ped"));
			log.report(ext.getTime() + "]\tWriting " + filenameRoot + ".ped");
			count = 1;
			invalidAbLookups = new Hashtable<Integer, Integer>();
			while (reader.ready()) {
				count++;
//				if (count % 100 == 0) {
//					System.out.println(count);
//				}
				
				temp = reader.readLine();
				line = temp.split(ext.determineDelimiter(temp));
				if (line.length < 7) {
					proj.message("Error - starting at line "+(count-1)+(line.length<3?"":" (individual "+line[0]+"-"+line[1]+")")+" there are only "+line.length+" columns in pedigree file '"+proj.PEDIGREE_FILENAME.getValue()+"'.\n"+
								"  Pedigree files require 7 columns with no header: FID IID FA MO SEX PHENO DNA\n"+
								"  where DNA is the sample name associated with the genotypic data (see the "+proj.SAMPLE_DIRECTORY.getValue(false, true)+" directory for examples)");
					reader.close();
					writer.close();
		            proj.progressMonitor.endTask(PROG_KEY + "_PEDEXPORT");
		            proj.progressMonitor.endTask(PROG_KEY);
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
						log.reportError("Error - the DNA# "+line[6]+" was listed in the pedigree file but "+line[6]+Sample.SAMPLE_DATA_FILE_EXTENSION+" was not found in directory: "+proj.SAMPLE_DIRECTORY.getValue(false, true));
						for (int i = 0; i<indices.length; i++) {
							writer.print(" 0 0");
						}
					} else {
						if (clusterFiltersFilename == null) {
							genotypes = fsamp.getForwardGenotypes(gcThreshold);
						} else {
							genotypes = MarkerSet.translateABtoForwardGenotypes(fsamp.getAB_GenotypesAfterFilters(markerNames, clusterFilterCollection, gcThreshold), abLookup);
						}
						for (int i = 0; i<indices.length; i++) {
							genIndex = genotypes[indices[i]];
							if (genIndex==0) {
								writer.print(" 0 0");
							} else {
								if (genIndex < 0) {
									if (!invalidAbLookups.containsKey(indices[i])) {
										log.reportError("Error - marker '"+markerNames[indices[i]]+"' was manually reclustered and requires a previously undefined AB lookup code ("+abLookup[indices[i]][0]+"/"+abLookup[indices[i]][1]+"); alleles will be set to missing for anyone with an invalid allele");
										invalidAbLookups.put(indices[i], invalidAbLookups.size());
									}
									writer.print(" 0 0");
								} else {
									genotype = Sample.ALLELE_PAIRS[genIndex];
									writer.print(" "+genotype.charAt(0)+" "+genotype.charAt(1));
								}
							}
						}
					}
				}
				writer.println();
				writer.flush();
				
	            proj.progressMonitor.updateTask(PROG_KEY + "_PEDEXPORT");
			}
			reader.close();
			writer.close();

			proj.progressMonitor.endTask(PROG_KEY + "_PEDEXPORT");
			
			if (invalidAbLookups.size() > 0) {
				proj.message("There "+(invalidAbLookups.size()==1?" was one marker ":"were "+invalidAbLookups.size()+" markers")+" with an invalid set of AB lookup codes that had been manually reclustered and now needs a full complement. Run \"java -cp Genvisis.jar cnv.filesys.ABLookup -h\" for options on how to fill these in, and check "+proj.getProperty(proj.DATA_DIRECTORY)+"invalid_AB_codes.out for a list of variants that this affects.");
				try {
					writer = new PrintWriter(new FileWriter(proj.DATA_DIRECTORY.getValue(false, true)+"invalid_AB_codes.out"));
					writer.println("MarkerNames\tA\tB");
					indices = Array.toIntArray(invalidAbLookups);
					for (int i = 0; i < indices.length; i++) {
						writer.println(markerNames[indices[i]]+"\t"+abLookup[indices[i]][0]+"\t"+abLookup[indices[i]][1]);
					}
					writer.close();
				} catch (Exception e) {
					proj.message("Error writing to " + proj.DATA_DIRECTORY.getValue(false, true)+"invalid_AB_codes.out");
					log.reportException(e);
				}
			}

            
		} catch (FileNotFoundException fnfe) {
			proj.message("Error: file \""+proj.PEDIGREE_FILENAME.getValue()+"\" not found");
            proj.progressMonitor.endTask(PROG_KEY + "_PEDEXPORT");
            proj.progressMonitor.endTask(PROG_KEY);
            return false;
		} catch (IOException ioe) {
			proj.message("Error reading file \""+proj.PEDIGREE_FILENAME.getValue()+"\"");
			proj.progressMonitor.endTask(PROG_KEY + "_PEDEXPORT");
            proj.progressMonitor.endTask(PROG_KEY);
			return false;
		}

		log.report(ext.getTime());
        proj.progressMonitor.endTask(PROG_KEY);
        
		return true;
	}
	

	/**
	 * Convert Genvisis data into a PLINK .bed data set.
	 * 
	 * @param proj
	 * @param markerList
	 * @param clusterFilterFileName
	 * @param gcThreshold
	 */
	public static boolean saveGenvisisToPlinkBedSet(Project proj, String plinkPrefix, String clusterFilterFileName, String targetMarkersFileName, float gcThreshold, boolean isSnpMajor) {
		String[] targetMarkers;
		int[] indicesOfTargetSamplesInProj;
		int[] indicesOfTargetMarkersInProj;
		byte[] chrsOfTargetMarkers;
		int[] posOfTargetMarkers;
		String[] allSamplesInProj;
		char[][] abLookup;
		String[] targetSamples;
		String outFileDirAndFilenameRoot;
		Logger log;
		
		log = proj.getLog();
		outFileDirAndFilenameRoot = proj.PROJECT_DIRECTORY.getValue() + plinkPrefix;
//		if (new File(outFileDirAndFilenameRoot + ".bed").exists() || new File(outFileDirAndFilenameRoot + ".bim").exists() || new File(outFileDirAndFilenameRoot + ".fam").exists()) {
//			log.reportError("System abort. PLINK binary file set \"" + outFileDirAndFilenameRoot + "\" .bed/.bim/.fam already exist. Please remove the file(s).");
//			return false;
//		}
//		if (new File(outFileDirAndFilenameRoot + ".bed").exists() || new File(outFileDirAndFilenameRoot + ".bim").exists() || new File(outFileDirAndFilenameRoot + ".fam").exists()) {
//			log.report("Found existing PLINK .ped file set in out file directory. Deleting these files.");
//			new File(outFileDirAndFilenameRoot + ".bed").delete();
//			new File(outFileDirAndFilenameRoot + ".bim").delete();
//			new File(outFileDirAndFilenameRoot + ".fam").delete();
//		}
		
		String PROG_KEY = "PLINKBINARYEXPORT";
		proj.progressMonitor.beginIndeterminateTask(PROG_KEY, "Creating .fam file", ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		targetSamples = createFamFile(proj, outFileDirAndFilenameRoot);
		proj.progressMonitor.endTask(PROG_KEY);
		
		allSamplesInProj = proj.getSamples();
		if (targetSamples != null) {
		    proj.progressMonitor.beginIndeterminateTask(PROG_KEY, "Loading sample data", ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
			indicesOfTargetSamplesInProj = getSortedIndicesOfTargetSamplesInProj(allSamplesInProj, targetSamples, log);
			targetSamples = new String[indicesOfTargetSamplesInProj.length];
			for (int i = 0; i < indicesOfTargetSamplesInProj.length; i++) {
				targetSamples[i] = allSamplesInProj[ indicesOfTargetSamplesInProj[i] ];
			}
	        proj.progressMonitor.endTask(PROG_KEY);
		} else {
			targetSamples = allSamplesInProj;
//			indicesOfTargetSamplesInProj = null;
			indicesOfTargetSamplesInProj = new int[allSamplesInProj.length];
			for (int i = 0; i < allSamplesInProj.length; i++) {
			    indicesOfTargetSamplesInProj[i] = i;
			}
		}

//		allMarkersInProj = proj.getMarkerNames();
		targetMarkers = proj.getTargetMarkers(targetMarkersFileName);
		if (targetMarkers != null) {
            proj.progressMonitor.beginDeterminateTask(PROG_KEY, "Loading marker data", targetMarkers.length, ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
			indicesOfTargetMarkersInProj = new int[targetMarkers.length];
			chrsOfTargetMarkers = new byte[targetMarkers.length];
			posOfTargetMarkers = new int[targetMarkers.length];
			getIndicesOfTargetMarkers(proj, targetMarkers, indicesOfTargetMarkersInProj, chrsOfTargetMarkers, posOfTargetMarkers);
			
            proj.progressMonitor.endTask(PROG_KEY);
		} else {
			indicesOfTargetMarkersInProj = null;
			targetMarkers = proj.getMarkerNames();
			chrsOfTargetMarkers = proj.getMarkerSet().getChrs();
			posOfTargetMarkers = proj.getMarkerSet().getPositions();
		}

		if (gcThreshold < 0) {
//			gcThreshold = proj.getFloat(proj.GC_THRESHOLD);
			gcThreshold = proj.GC_THRESHOLD.getValue().floatValue();
		}
		
		proj.progressMonitor.beginIndeterminateTask(PROG_KEY, "Creating .bed file", ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		if (isSnpMajor) {
			abLookup = createBedFileSnpMajor10KperCycle(proj, targetMarkers, indicesOfTargetMarkersInProj, targetSamples, indicesOfTargetSamplesInProj, clusterFilterFileName, gcThreshold, outFileDirAndFilenameRoot, log);
//			abLookup = createBedFileSnpMajorAllInMemory(proj, targetMarkers, indicesOfTargetMarkersInProj, targetSamples, indicesOfTargetSamplesInProj, bedFilenameRoot, gcThreshold, bedFilenameRoot);
		} else {
			abLookup = createBedFileIndividualMajor(proj, targetSamples, targetMarkers, indicesOfTargetMarkersInProj, clusterFilterFileName, gcThreshold, outFileDirAndFilenameRoot);
		}
		proj.progressMonitor.endTask(PROG_KEY);
		
		if (abLookup == null) {
			log.reportError("Error - failed to create PLINK files due to lack of an AB lookup file");
			new File(outFileDirAndFilenameRoot+".fam").delete();
			return false;
		}

		proj.progressMonitor.beginIndeterminateTask(PROG_KEY, "Creating .bim file", ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		createBimFile(targetMarkers, chrsOfTargetMarkers, posOfTargetMarkers, abLookup, outFileDirAndFilenameRoot, log);
		proj.progressMonitor.endTask(PROG_KEY);

		return true;
	}

	public static int[] getSortedIndicesOfTargetSamplesInProj(String[] allSampInProj, String[] targetSamples, Logger log) {
//		String[] sampleNames;
		int[] indicesOfTargetSampInProj;
		Hashtable<String, Integer> hash;
		Enumeration<String> keys;
		hash = new Hashtable<String, Integer>();
		boolean found;
//		sampleNames = proj.getSamples();

//		indices = new int[sampList.length];
		for (int i = 0; i < targetSamples.length; i++) {
			if (hash.containsKey(targetSamples[i])) {
				log.reportError("Warning - duplicate sample id in the list of samples to include: " + targetSamples[i]);
			} else {
				found = false;
				for (int j = 0; j < allSampInProj.length; j++) {
					if (targetSamples[i].equals(allSampInProj[j])) {
						hash.put(targetSamples[i], j);
						found = true;
						break;
					}
				}
				if (! found) {
					log.reportError("Warning - The following from the target sample list was not found in the list of all samples in the project: " + targetSamples[i]);
				}
			}
		}

		indicesOfTargetSampInProj = new int[hash.size()];
		keys = hash.keys();
		for (int i = 0; i < indicesOfTargetSampInProj.length; i++) {
			indicesOfTargetSampInProj[i] = hash.get(keys.nextElement());
		}

		Arrays.sort(indicesOfTargetSampInProj);
		return indicesOfTargetSampInProj;
	}

//	public static int[] getIndicesOfMarkersToInclude(Project proj, log) {
	public static void getIndicesOfTargetMarkers(Project proj, String[] inputTargetMarkers, int[] outputIndicesOfTargetMarkers, byte[] outputChrOfTagetMarkers, int[] outputPosOfTargetMarkers) {
		String[] allMarkersInProj;
//		Hashtable<String,String> hash;
		MarkerSet markerSet;
		byte[] chrs;
		int[] positions;
		boolean[] found;
		Logger log = proj.getLog();
	
//		hash = new Hashtable<String,String>();
		if (inputTargetMarkers == null) {
			log.reportError("inputTargetMarkers is null");

		} else {
			markerSet = proj.getMarkerSet();
			allMarkersInProj = markerSet.getMarkerNames();
			chrs = markerSet.getChrs();
			positions = markerSet.getPositions();
			found = new boolean[inputTargetMarkers.length];
			log.report("Loading " + inputTargetMarkers.length + " target markers");
			HashMap<String, Integer> markerPositions = new HashMap<String, Integer>();
			for (int j = 0; j < allMarkersInProj.length; j++) {
			    markerPositions.put(allMarkersInProj[j], j);
			}
			for (int i = 0; i < inputTargetMarkers.length; i++) {
//				if (hash.containsKey(inputTargetMarkers[i])) {
//					System.err.println("Warning - duplicate marker name: " + inputTargetMarkers[i] + " in targetMarkers.txt");
//				}
//				for (int j = 0; j < allMarkersInProj.length; j++) {
				    if (markerPositions.containsKey(inputTargetMarkers[i])) {
				        outputIndicesOfTargetMarkers[i] = markerPositions.get(inputTargetMarkers[i]);
				        found[i] = true;
				    }
//					if (inputTargetMarkers[i].equals(allMarkersInProj[j])) {
//						outputIndicesOfTargetMarkers[i] = j;
//						found[i] = true;
//						break;
//					}
//				}
//				hash.put(allMarkersInProj[i], i+"");
				proj.progressMonitor.updateTask("PLINKBINARYEXPORT");
			}
			Arrays.sort(outputIndicesOfTargetMarkers); 
			for (int i = 0; i < found.length; i++) {
				if (! found[i]) {
					log.reportError("Warning - the following marker from target marker list is not found in whole project's marker list: " + inputTargetMarkers[i]);
					break;
				} else {
					inputTargetMarkers[i] = allMarkersInProj[outputIndicesOfTargetMarkers[i]];
					outputChrOfTagetMarkers[i] = chrs[outputIndicesOfTargetMarkers[i]];
					outputPosOfTargetMarkers[i] = positions[outputIndicesOfTargetMarkers[i]];
				}
			}
		}

//		chrs = markerSet.getChrs();
//		positions = markerSet.getPositions();
////		targetMarkers = proj.getFilename(proj.TARGET_MARKERS_FILENAME, false, false);
////		if (new File(targetMarkers).exists()) {
////			targets = HashVec.loadFileToStringArray(targetMarkers, false, false, new int[] {0}, false);
//		if (inputTargetMarkers != null) {
//			indicesOfMarkersSelected = new int[inputTargetMarkers.length];
//			prob = false;
//			for (int i = 0; i < inputTargetMarkers.length; i++) {
//				if (hash.containsKey(inputTargetMarkers[i])) {
//					indicesOfMarkersSelected[i] = Integer.parseInt(hash.get(inputTargetMarkers[i]));
//				} else {
//					System.err.println("Error - target marker '" + inputTargetMarkers[i] + "' was not found in the MarkerSet");
//					prob = true;
//				}
//			}
//			if (prob) {
//				System.exit(1);
//			}
//
////			chrOfMarkersSelected = new byte[indicesOfMarkersSelected.length];
////			posOfMarkersSelected = new int[indicesOfMarkersSelected.length];
////			for (int i = 0; i < indicesOfMarkersSelected.length; i++) {
////				chrOfMarkersSelected[i] = chrs[indicesOfMarkersSelected[i]];
////				posOfMarkersSelected[i] = positions[indicesOfMarkersSelected[i]];
////			}
//
//		} else {
//			if (! inputTargetMarkers.equals("")) {
//				System.out.println("FYI, since target markers file '" + inputTargetMarkers + "' was not found, all markers will be exported to PLINK");
//			}
//
////			targets = HashVec.getKeys(hash);
////			indices = new int[targets.length];
////			for (int i = 0; i<targets.length; i++) {
////				indices[i] = Integer.parseInt(hash.get(targets[i]));
////			}
//			indicesOfMarkersSelected = Array.intArray(allMarkersInProj.length);
////			chrOfMarkersSelected = chrs;
////			posOfMarkersSelected = positions;
//		}
//
//		return indicesOfMarkersSelected;
	}

	/**
	 * Convert Genvisis data into binary PLINK .bed format (Individual Major, or in our term - organized by samples).
	 * This is normally used as part of PlinkData.createBinaryFileSetFromGenvisisData()
	 * @param proj
	 * @param targetSamples samples selected to convert 
	 * @param targetMarkers markers selected to convert
	 * @param indicesOfTargetMarkers
	 * @param clusterFilterFileName
	 * @param gcThreshold
	 * @param bedDirAndFilenameRoot
	 * @param log
	 * @return
	 */
	public static char[][] createBedFileIndividualMajor(Project proj, String[] targetSamples, String[] targetMarkers, int[] indicesOfTargetMarkers, String clusterFilterFileName, float gcThreshold, String bedDirAndFilenameRoot) {
//		String[] markList;
		RandomAccessFile out;
		byte[] outStream;
		byte[] genotypes;
//		byte[] genotypesOfTargetMarkers;
		ClusterFilterCollection clusterFilterCollection;
		Sample fsamp;
		char[][] abLookup = null;
//		MarkerSet markerSet = proj.getMarkerSet();
//		String[] markerNames;

//		allMarkersInProj = proj.getMarkerNames();
//		markList = new String[indicesOfSelectedMarkers.length];
//		for (int i = 0; i < indicesOfSelectedMarkers.length; i++) {
//			markList[i] = allMarkersInProj[indicesOfSelectedMarkers[i]];
//		}

        clusterFilterCollection = proj.getClusterFilterCollection();

        String PROG_KEY = "EXPORTBINARYBEDBATCH";
        proj.progressMonitor.beginDeterminateTask(PROG_KEY, "Exporting data to .bed file", targetSamples.length, ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
        
		try {
			if (clusterFilterCollection == null) {
				abLookup = null;
			} else {
				abLookup = new ABLookup(targetMarkers, proj.AB_LOOKUP_FILENAME.getValue(), true, true, proj.getLog()).getLookup();
			}
			
			out = new RandomAccessFile(bedDirAndFilenameRoot + ".bed", "rw");
			outStream = new byte[3];
			outStream[0] = (byte) 108;	// 0b01101100
			outStream[1] = (byte) 27;	// 0b00011011
			outStream[2] = (byte) 0;	// 0b00000000 <-- Be careful here
			out.write(outStream);

			for (int i = 0; i < targetSamples.length; i++) {
				fsamp = proj.getFullSampleFromRandomAccessFile(targetSamples[i]);

				if (fsamp == null) {
					System.err.println("Error - the DNA# " + targetSamples[i] + " was listed in the pedigree file but " + targetSamples[i] + Sample.SAMPLE_DATA_FILE_EXTENSION+ " was not found in directory: " + proj.SAMPLE_DIRECTORY.getValue(false, true));
					genotypes = new byte[1];

				} else {
					if (clusterFilterCollection == null) {
						genotypes = fsamp.getAB_Genotypes(indicesOfTargetMarkers);
					} else {
						genotypes = fsamp.getAB_GenotypesAfterFilters(targetMarkers, indicesOfTargetMarkers, clusterFilterCollection, gcThreshold);
//						genotypes = proj.getMarkerSet().translateABtoForwardGenotypes(genotypes, abLookup);
					}
					out.write(encodePlinkBedBytesForASingleMarkerOrSample(genotypes));

//					genotypesOfTargetMarkers = new byte[indicesOfTargetMarkers.length];
//					for (int j = 0; j < indicesOfTargetMarkers.length; j++) {
//						genotypesOfTargetMarkers[j] = genotypes[indicesOfTargetMarkers[i]];
//					}
//					out1.write(encodePlinkBedBytesForASingleMarkOrSamp(genotypesOfTargetMarkers));

//					if (clusterFilterCollection == null) {
//						genotypes = fsamp.getForwardGenotypes(gcThreshold);
//					} else {
//						genotypes = markerSet.translateABtoForwardGenotypes(fsamp.getAB_GenotypesAfterFilters(targetMarkers, clusterFilterCollection, gcThreshold), abLookup);
//					}

//					alleles = new String[genotypesSelected.length][2];
//					for (int j = 0; j<indicesOfTargetMarkers.length; j++) {
//						genotypesSelected[j] = genotypes[indicesOfTargetMarkers[i]];
//
//						//TODO Something does not make sense here. alleles[j] seems to be wrong here.
//						if (genotypesSelected[j] == 0) {
//							alleles[j] = new String[] {"0", "0"};
//						} else {
//							tmp = Sample.ALLELE_PAIRS[genotypesSelected[j]];
//							alleles[j][0] = tmp.charAt(0) + "";
//							alleles[j][1] = tmp.charAt(1) + "";
//						}
//					}
				}
				
				proj.progressMonitor.updateTask(PROG_KEY);
			}

			out.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		proj.progressMonitor.endTask(PROG_KEY);
		
		return abLookup;
	}

	/**
	 * Convert Genvisis data to PLINK .bed format (SNP Major, or in our term - organized by markers)
	 * This is normally used as part of PlinkData.createBinaryFileSetFromGenvisisData()
	 * @param proj
	 * @param targetMarkers markers selected to convert
	 * @param indicesOfTargetMarkersInProj
	 * @param targetSamples samples selected to convert
	 * @param indicesOfTargetSamplesInProj
	 * @param clusterFilterFileName
	 * @param gcThreshold
	 * @param bedDirAndFilenameRoot
	 * @param log
	 * @return
	 */
	public static char[][] createBedFileSnpMajorAllInMemory(Project proj, String[] targetMarkers, int[] indicesOfTargetMarkersInProj, String[] targetSamples, int[] indicesOfTargetSamplesInProj, String clusterFilterFileName, float gcThreshold, String bedDirAndFilenameRoot) {
		RandomAccessFile out1;
		byte[] outStream;
		byte[] genotypes;
		byte[] genotypesOfTargetSamples;
		ClusterFilterCollection clusterFilterCollection;
		MarkerDataLoader markerDataLoader;
		char[][] abLookup = null;
		Hashtable<String, Vector<String>> batches;
		MarkerData[] markerData;
		String[] filenames;
		Vector<String> v;
		int[] indicesOfMarkersInFileForCurrentFile;
		String[] markersOfThisFile;
		String[] temp;
		Hashtable<String, Float> outliersHash;
		long sampleFingerPrint;
		String[] allSamplesInProj;
		int[] indicesOfMarkersInProjForCurrentFile;
		String markerDataDir;
		long startTime;
		
		startTime = new Date().getTime();

		try {
	        clusterFilterCollection = proj.getClusterFilterCollection();
			if (clusterFilterCollection == null) {
				abLookup = null;
			} else {
				abLookup = new ABLookup(targetMarkers, proj.AB_LOOKUP_FILENAME.getValue(), true, true, proj.getLog()).getLookup();
			}
			
			sampleFingerPrint = proj.getSampleList().getFingerprint();
			allSamplesInProj = proj.getSamples();
			markerDataDir = proj.MARKER_DATA_DIRECTORY.getValue(false, true);
			markerDataLoader = new MarkerDataLoader(proj, targetMarkers, 0);
			outliersHash = MarkerDataLoader.loadOutliers(proj);
			batches = markerDataLoader.getBatches();
			filenames = HashVec.getKeys(batches);
			genotypesOfTargetSamples = new byte[indicesOfTargetSamplesInProj.length];

			out1 = new RandomAccessFile(bedDirAndFilenameRoot + ".bed", "rw");
			outStream = new byte[3];
			outStream[0] = (byte) 108;	// 0b01101100
			outStream[1] = (byte) 27;	// 0b00011011
			outStream[2] = (byte) 1;	// 0b00000001 <-- be careful here
			out1.write(outStream);

//			for (int i = 0; i < targetMarkers.length; i++) {
			for (int i = 0; i < filenames.length; i++) {
				v = batches.get(filenames[i]);
				markersOfThisFile = new String[v.size()];
				indicesOfMarkersInFileForCurrentFile = new int[markersOfThisFile.length];
				indicesOfMarkersInProjForCurrentFile = new int[markersOfThisFile.length];
				for (int j = 0; j < v.size(); j++) {
					temp = v.elementAt(j).split("\t");
					markersOfThisFile[j] = temp[0];
					indicesOfMarkersInFileForCurrentFile[j] = Integer.parseInt(temp[1]);
					if(indicesOfTargetMarkersInProj != null) {
						for (int k = 0; k < indicesOfTargetMarkersInProj.length; k++) {
							if (markersOfThisFile[j].equals(targetMarkers[k])) {
								indicesOfMarkersInProjForCurrentFile[j] = indicesOfTargetMarkersInProj[k];
								break;
							}
						}
					} else {
						for (int k = 0; k < targetMarkers.length; k++) {
							if (markersOfThisFile[j].equals(targetMarkers[k])) {
								indicesOfMarkersInProjForCurrentFile[j] = k;
								break;
							}
						}
					}
				}
				markerData = MarkerDataLoader.loadFromRAF(null, null, null, allSamplesInProj, markerDataDir + filenames[i], indicesOfMarkersInProjForCurrentFile, indicesOfMarkersInFileForCurrentFile, false, true, false, false, true, sampleFingerPrint, outliersHash, proj.getLog());

				for (int j = 0; j < markerData.length; j++) {
					genotypes = markerData[j].getAbGenotypesAfterFilters(clusterFilterCollection, markersOfThisFile[j], 0);
					for (int k = 0; k < indicesOfTargetSamplesInProj.length; k++) {
						genotypesOfTargetSamples[k] = genotypes[indicesOfTargetSamplesInProj[k]];
					}
					out1.write(encodePlinkBedBytesForASingleMarkerOrSample(genotypesOfTargetSamples));
				}

				//TODO Something does not make sense here. alleles[j] seems to be wrong here.
//				genotypesSelected = new byte[indicesOfSelectedSamples.length];
//				for (int j = 0; j < indicesOfSelectedSamples.length; j++) {
//					genotypesSelected[j] = genotypes[indicesOfSelectedSamples[i]];

//					//TODO Something does not make sense here. alleles[j] seems to be wrong here.
//					if (genotypesSelected[j] == 0) {
//						abLookup[j] = new char[] {0, 0};
//
//					} else {
//						tmp = Sample.ALLELE_PAIRS[genotypesSelected[j]];
//						abLookup[j][0] = tmp.charAt(0);
//						abLookup[j][0] = tmp.charAt(1);
//					}

//				}

//				out1.write(encodePlinkBedBytesForASingleMarkOrSamp(genotypes));
			}
			out1.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		proj.getLog().report("Finished creating binary PLINK files in "+ext.getTimeElapsed(startTime));

		return abLookup;
	}

	/**
	 * Convert Genvisis data to PLINK .bed format (SNP Major, or in our term - organized by markers)
	 * This is normally used as part of PlinkData.createBinaryFileSetFromGenvisisData()
	 * @param proj
	 * @param targetMarkers markers selected to convert
	 * @param indicesOfTargetMarkersInProj
	 * @param targetSamples samples selected to convert
	 * @param indicesOfTargetSamplesInProj
	 * @param clusterFilterFileName
	 * @param gcThreshold
	 * @param bedDirAndFilenameRoot
	 * @param log
	 * @return
	 */
	public static char[][] createBedFileSnpMajor10KperCycle(Project proj, String[] targetMarkers, int[] indicesOfTargetMarkersInProj, String[] targetSamples, int[] indicesOfTargetSamplesInProj, String clusterFilterFileName, float gcThreshold, String bedDirAndFilenameRoot, Logger log) {
		RandomAccessFile out;
		byte[] outStream;
		byte[] genotypes;
		byte[] genotypesOfTargetSamples;
		ClusterFilterCollection clusterFilterCollection;
		MarkerDataLoader markerDataLoader;
		char[][] abLookup;
		Hashtable<String, Vector<String>> batches;
		MarkerData[] markerData;
		String[] filenames;
		Vector<String> v;
		int[] indicesOfMarkersInFileForCurrentFile;
		String[] markersOfThisFile;
		String[] temp;
		Hashtable<String, Float> outliersHash;
		long sampleFingerPrint;
		String[] allSamplesInProj;
		int[] indicesOfMarkersInProjForCurrentFile;
		String dir;
		long startTime, subTime;
		Hashtable<String,Integer> hash;
		HashMap<String, Integer> projHash;
		int targetIndex;

		startTime = new Date().getTime();

		if (clusterFilterFileName == null) {
	        clusterFilterCollection = null;
		} else {
			clusterFilterCollection = ClusterFilterCollection.load(clusterFilterFileName, proj.JAR_STATUS.getValue());
		}
		
        if (clusterFilterCollection == null) {
			abLookup = new char[targetMarkers.length][];
		} else if (Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false))) {
			abLookup = new ABLookup(targetMarkers, proj.AB_LOOKUP_FILENAME.getValue(), true, true, proj.getLog()).getLookup();
		} else {
			proj.message("Error - could not find AB lookup file '"+proj.AB_LOOKUP_FILENAME.getValue()+"'; this file needs to be created, as it is not otherwise possible to export to PLINK when there are cluster filters.");
			return null;
		}
		
		hash = new Hashtable<String,Integer>();
		for (int i = 0; i<targetMarkers.length; i++) {
			if (hash.containsKey(targetMarkers[i])) {
				System.err.println("Warning - duplicate marker name: "+targetMarkers[i]);
			}
			hash.put(targetMarkers[i], i);
		}
		
		projHash = new HashMap<String, Integer>();
		if (indicesOfTargetMarkersInProj != null) {
    		for (int i = 0; i < indicesOfTargetMarkersInProj.length; i++) {
    		    if (projHash.containsKey(targetMarkers[i])) {
                    System.err.println("Warning - duplicate marker name: "+targetMarkers[i]);
                }
    		    projHash.put(targetMarkers[i], indicesOfTargetMarkersInProj[i]);
    		}
		} else {
		    for (int i = 0; i < targetMarkers.length; i++) {
		        if (projHash.containsKey(targetMarkers[i])) {
		            System.err.println("Warning - duplicate marker name: "+targetMarkers[i]);
		        }
		        projHash.put(targetMarkers[i], i);
		    }
		}

		sampleFingerPrint = proj.getSampleList().getFingerprint();
		allSamplesInProj = proj.getSamples();
		dir = proj.MARKER_DATA_DIRECTORY.getValue(false, true);
		subTime = new Date().getTime();
		markerDataLoader = new MarkerDataLoader(proj, targetMarkers, 0);
		log.report("MarkerDataLoader initialized in "+ext.getTimeElapsed(subTime));
		outliersHash = MarkerDataLoader.loadOutliers(proj);
		batches = markerDataLoader.getBatches();
		filenames = HashVec.getKeys(batches);
		genotypesOfTargetSamples = new byte[indicesOfTargetSamplesInProj.length];
		
		String PROG_KEY = "EXPORTBINARYBEDBATCH";
		
		int expUpdateCount = 0;
		for (int i = 0; i < filenames.length; i++) {
		    expUpdateCount += batches.get(filenames[i]).size();
		}
		proj.progressMonitor.beginDeterminateTask(PROG_KEY, "Exporting data to .bed file", expUpdateCount, ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		
		try {
			out = new RandomAccessFile(bedDirAndFilenameRoot + ".bed", "rw");
			outStream = new byte[3];
			outStream[0] = (byte) 108;	// 0b01101100
			outStream[1] = (byte) 27;	// 0b00011011
			outStream[2] = (byte) 1;	// 0b00000001 <-- be careful here
			out.write(outStream);

			subTime = new Date().getTime();

			for (int i = 0; i < filenames.length; i++) {
			    proj.progressMonitor.changeTaskLabelWithUpdate(PROG_KEY, "Exporting data to .bed file from marker file ... " + filenames[i]);
				
				v = batches.get(filenames[i]);
				markersOfThisFile = new String[v.size()];
				indicesOfMarkersInFileForCurrentFile = new int[markersOfThisFile.length];
				indicesOfMarkersInProjForCurrentFile = new int[markersOfThisFile.length];
				for (int j = 0; j < v.size(); j++) {
					temp = v.elementAt(j).split("\t");
					markersOfThisFile[j] = temp[0];
					indicesOfMarkersInFileForCurrentFile[j] = Integer.parseInt(temp[1]);
					indicesOfMarkersInProjForCurrentFile[j] = projHash.get(temp[0]);
				}
				
				subTime = new Date().getTime();
				proj.progressMonitor.beginIndeterminateTask(PROG_KEY + filenames[i] + "_load", "Loading marker data from file ... " + filenames[i], ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
				markerData = MarkerDataLoader.loadFromRAF(null, null, null, allSamplesInProj, dir + filenames[i], /*indicesOfTargetMarkersInProj*/indicesOfMarkersInProjForCurrentFile, indicesOfMarkersInFileForCurrentFile, false, true, false, false, true, sampleFingerPrint, outliersHash, proj.getLog());
				proj.progressMonitor.endTask(PROG_KEY + filenames[i] + "_load");
				
				subTime = new Date().getTime();
				for (int j = 0; j < markerData.length; j++) {
					targetIndex = hash.get(markersOfThisFile[j]);
					genotypes = markerData[j].getAbGenotypesAfterFilters(clusterFilterCollection, markersOfThisFile[j], 0);
					for (int k = 0; k < indicesOfTargetSamplesInProj.length; k++) {
						genotypesOfTargetSamples[k] = genotypes[indicesOfTargetSamplesInProj[k]];
					}
					if (abLookup[targetIndex] == null) {
						abLookup[targetIndex] = markerData[j].getAB_AlleleMappings();
					}
					
					out.write(encodePlinkBedBytesForASingleMarkerOrSample(genotypesOfTargetSamples));
//					if (j > 0 && j % mod == 0) {
//					    proj.progressMonitor.changeTaskLabelWithUpdate(PROG_KEY + filenames[i] + "_export", "Exported " + j + " of " + markerData.length + " markers from file " + filenames[i]);
//					}
					proj.progressMonitor.updateTask(PROG_KEY);
				}
			}
			out.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		proj.progressMonitor.endTask(PROG_KEY);
		proj.getLog().report("Finished creating binary PLINK files in "+ext.getTimeElapsed(startTime));

		return abLookup;
	}

	/**
	 * Convert Genvisis data to PLINK .bim format (SNP Major, or in our term - organized by markers)
	 * This is normally used as part of PlinkData.createBinaryFileSetFromGenvisisData()
	 * @param targetMarkers
	 * @param chrsOfTargetMarkers
	 * @param posOfTargetMarkers
	 * @param abLookup
	 * @param bimDirAndFilenameRoot
	 * @param log
	 * @return
	 */
	public static boolean createBimFile(String[] targetMarkers, byte[] chrsOfTargetMarkers, int[] posOfTargetMarkers, char[][] abLookup, String bimDirAndFilenameRoot, Logger log) {
		PrintWriter writer;

		if (abLookup == null) {
			System.err.println("Error - abLookup cannot be null; failed to create .bim file");
			return false;
		}
		
		try {
			writer = new PrintWriter(new FileWriter(bimDirAndFilenameRoot + ".bim"));
			for (int i = 0; i < targetMarkers.length; i++) {
				writer.println(chrsOfTargetMarkers[i] + "\t" + targetMarkers[i] + "\t0\t" + posOfTargetMarkers[i] + "\t" + abLookup[i][0] + "\t" + abLookup[i][1]); //TODO alleles[][] matching chrs[]
			}
			writer.close();

		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: failed to write to gwas.map (in use?)");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error writing to gwas.map");
			System.exit(2);
		}

		return true;
	}

	/**
	 * Convert Genvisis data to PLINK .fam format (SNP Major, or in our term - organized by markers)
	 * This is normally used as part of PlinkData.createBinaryFileSetFromGenvisisData()
	 * @param proj
	 * @param famDirAndFilenameRoot
	 * @param log
	 * @return
	 */
	public static String[] createFamFile(Project proj, String famDirAndFilenameRoot) {
		BufferedReader reader;
		PrintWriter writer;
		int count;
		String temp;
		String[] line;
		Vector<String> dna;
		String[] allSamples;
		String filename;
		Logger log;
		
		log = proj.getLog();
		allSamples = proj.getSamples();
		dna = new Vector<String>();

		try {
			filename = proj.PEDIGREE_FILENAME.getValue();
			if (!new File(filename).exists()) {
				log.reportError("Error - pedigree file ('"+filename+"') is not found.  Cannot create .fam file.");
				return null;
			}
			reader = new BufferedReader(new FileReader(proj.PEDIGREE_FILENAME.getValue()));
			writer = new PrintWriter(new FileWriter(famDirAndFilenameRoot+".fam"));
			count = 1;
			while (reader.ready()) {
				count++;
				temp = reader.readLine().trim();
				line = temp.split(ext.determineDelimiter(temp));
				if (temp.equals("")) {
					// then do nothing
				} else if (line.length < 7) {
					log.reportError("Error - starting at line "+(count-1)+(line.length<3?"":" (individual "+line[0]+"-"+line[1]+")")+" there are only "+line.length+" columns in pedigree file '"+proj.PEDIGREE_FILENAME.getValue()+"'.");
					log.reportError("  Pedigree files require 7 columns with no header: FID IID FA MO SEX PHENO DNA");
					log.reportError("  where DNA is the sample name associated with the genotypic data (see the "+proj.SAMPLE_DIRECTORY.getValue(false, true)+" directory for examples)");
					reader.close();
					writer.close();
					return null;
				} else if (ext.isMissingValue(line[6])) {
//					dna.add(null);
				} else if (ext.indexOfStr(line[6], allSamples) == -1) {
					log.reportError("Warning - sample '" + line[6] + "' from '" + proj.PEDIGREE_FILENAME.getValue() + "' is not found in the project's list of samples, and is ignored.");
					if (line.length != 7) {
						log.reportError("      check to make sure that there are no spaces in your IDs; as this will be parsed as a new column; for example there are "+line.length+" columns here, and we only want 7");
					}
					//dna.add(null);
				} else {
					dna.add(line[6]);
					writer.println(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]);
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + famDirAndFilenameRoot+".fam" + "\" not found in current directory");
			return null;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + famDirAndFilenameRoot+".fam" + "\"");
			return null;
		}

		return Array.toStringArray(dna);
	}

	/**
	 * Transpose an array of byte steam read in from a PLINK .bed file from individual-major format to SNP-major format.
	 * 
	 * Note: This method does not address the memory management issue. The caller of this method is supposed to implement
	 * the optimization of memory usage.
	 * 
	 * This method is part of the process to convert a PLINK .bed format data from individual-major (Sample lead, or
	 * organized by sample) to SNP-major (Marker lead, or organized by marker). This method is normally used together
	 * with another method reading in and writing to PLINK .bed files.
	 * 
	 * This method can also reversely convert PLINK .bed format data from SNP-major (Marker lead, or organized by marker)
	 * to individual-major (Sample lead, or organized by sample). Just treat all the variables labeled with Mark as if
	 * they were labeled with Samp, and also Samp with Mark.
	 *  
	 * The index parameters are all inclusive, meaning samples identified by indexOfStartSamp and indesOfEndSamp are all
	 * going to be included in the output, and so do the markers.
	 * 
	 * @param inputSampBytes array of a PLINK .bed stream, the input of this method
	 * @param indexOfStartSamp index basing on all the samples in the input array, starting from 0. Please not to confuse
	 * with index of the input array, because 1 single element of the array contains up to 4 samples.
	 * @param indexOfEndSamp index basing on all the samples in the .bed file, starting from. 0 Please not to confuse with
	 * index of the input array, because 1 single element of the array contains up to 4 samples.
	 * @param indexOfStartMark index basing on all the markers in the .bed file, starting from 0. Please not to confuse
	 * with index of the output array, because 1 single element of the array contains up to 4 markers.
	 * @param indexOfEndMark index basing on all the markers in the .bed file, starting from 0. Please not to confuse
	 * with index of the output array, because 1 single element of the array contains up to 4 markers.
	 * @param outputMarkBytes array of a PLINK .bed stream, the output of this method
	 * @param numberOfTotalMarkersInProj totoal number of markers in the PLINK .bed file
	 * @param numberOfTotalSamplesInProj totoal number of samples in the PLINK .bed file
	 */
	public static void transposeBedBytes(byte[] inputSampBytes, int indexOfStartSamp, int indexOfEndSamp, int indexOfStartMark, int indexOfEndMark, byte[] outputMarkBytes, int numberOfTotalMarkersInProj, int numberOfTotalSamplesInProj) {
//		int iSamp;
//		int iStartByteSamp;
//		byte iStartBitSamp;
//		int iEndByteSamp;
//		byte iEndBitSamp;

		int nBytesPerMark;
		int nBytesPerSamp;
		int offsetMark;
		int offsetSamp;
		int iByteSamp;
		byte iBitSamp;
		int iByteMark;
		byte iBitMark;

//		iStartByteSamp = (int) Math.ceil((double) iStartMark / 4);
//		iStartBitSamp = (byte) (iStartMark % 4);
//		iEndByteSamp = (int) Math.ceil((double) iEndMark / 4);
//		iEndBitSamp = (byte) (iEndMark % 4);

		nBytesPerMark = (int) Math.ceil((double) numberOfTotalSamplesInProj / 4);
		nBytesPerSamp = (int) Math.ceil((double) numberOfTotalMarkersInProj / 4);
//		iByteMark = iStartMark * bytesPerMark + iStartSamp / 4;
		offsetMark = indexOfStartMark * nBytesPerMark;
		offsetSamp = indexOfStartSamp * nBytesPerSamp;

		for (int i = indexOfStartSamp; i < indexOfEndSamp; i++) {
			iByteMark = offsetMark + i / 4;
			iBitMark = (byte) (i % 4);
			
			for (int j = indexOfStartMark; j < indexOfEndMark; j++) {
				iByteSamp = offsetSamp + j / 4;
				iBitSamp = (byte) (j % 4);

				outputMarkBytes[iByteMark] = (byte) ((outputMarkBytes[iByteMark] & (0xff - (0x03 << (iBitMark * 2)))) | (inputSampBytes[iByteSamp] & (0x03 << (iBitSamp * 2))));
			}
		}

//		for (int i = iStartByteSamp; i < iEndByteSamp - 1; i++) {
//			iBitMark = (byte) (iStartSamp % 4);
//			for (int j = iStartBitSamp; j < 4; j++) {
//				outputMarkBytes[iByteMark] = (byte) ((outputMarkBytes[iByteMark] & (0xff - (0x03 << (iBitMark * 2)))) | (inputASingleSamplesBytes[i] & (0x03 << (j * 2))));
//				iByteMark += bytesPerMark;
//			}
//			iStartBitSamp = 0;
//			iBitMark = 0;
//		}

//		for (int j = iStartBitSamp; j <= 4; j++) {
//			iByteMark = i * bytesPerMark + j / 4;
//			outputMarkBytes[iByteMark] = (byte) ((outputMarkBytes[iByteMark] & 0x03) | (inputASingleSamplesBytes[i] & 0x3));
//		}
	}

	/**
	 * Transpose one byte of data read in from a PLINK .bed file from individual-major format to SNP-major format.
	 *  
	 * This method is part of the process to convert a PLINK .bed format data from individual-major (Sample lead, or organized
	 * by sample) to SNP-major (Marker lead, or organized by marker). This method is normally used together with another
	 * method to read in a PLINK .bed file and write to another one, and looping over the byte stream loaded from the
	 * PLINK .bed file.
	 * 
	 * This method can also reversely convert PLINK .bed format data from SNP-major (Marker lead, or organized by marker) to
	 * individual-major (Sample lead, or organized by sample). Just treat all the variables labeled with Mark as if they were
	 * labeled with Samp, and also Samp with Mark.
	 *  
	 * The index parameters are all inclusive, meaning samples identified by indexOfStartSamp and indexOfEndSamp are all going
	 * to be included in the output, and so do the markers.
	 * 
	 * @param inputSampByte
	 * @param indexOfCurrentSampInProj
	 * @param indexOfEndSamp
	 * @param indexOfStartMark
	 * @param indexOfEndMark
	 * @param outputMarkBytes
	 * @param numberOfBytesPerSamp
	 * @param numberOfBytesPerMark
	 */
	public static void transposeOneBedByte(byte inputSampByte, int indexOfCurrentSampInProj, int indexOfEndSamp, int indexOfStartMark, int indexOfEndMark, byte[] outputMarkBytes, int numberOfBytesPerSamp, int numberOfBytesPerMark) {
		int offsetMark;
		int iByteMark;
		byte iBitMark;
		byte iBitSamp;

		offsetMark = indexOfStartMark * numberOfBytesPerMark;
		iByteMark = offsetMark + indexOfCurrentSampInProj / 4;
		iBitMark = (byte) (indexOfCurrentSampInProj % 4);

		for (int j = indexOfStartMark; j < indexOfEndMark; j++) {
			iBitSamp = (byte) (j % 4);

			outputMarkBytes[iByteMark] = (byte) ((outputMarkBytes[iByteMark] & (0xff - (0x03 << (iBitMark * 2)))) | (inputSampByte & (0x03 << (iBitSamp * 2))));
		}
	}

	/**
	 * Translate Genvisis genotype genotype to PLINK genotype. This is specifically for markers on chromosome X of males.
	 *  
	 * @param genvisisGenotype
	 * @param chromosome
	 * @param sex
	 * @return
	 */
	public static byte convertGenvisisGenotypeToPlink (byte genvisisGenotype, byte chromosome, byte sex) {
		byte plinkGenotype;
		if (chromosome == 23 && genvisisGenotype == 2) {
			plinkGenotype = 1;
		} else if (chromosome == 24 && genvisisGenotype == 2) {
			plinkGenotype = 1;
		} else {
			plinkGenotype = genvisisGenotype;
		}
		return plinkGenotype;
	}

	/**
	 * Translate PLINK genotype genotype to Genvisis genotype. This is specifically for markers on chromosome X of males.
	 *  
	 * @param genvisisGenotype
	 * @param chromosome
	 * @return
	 */
	public static byte convertPlinkGenotypeToGenvisis (byte plinkGenotype, byte chromosome) {
		byte genvisisGenotype;
		if (chromosome == 23 && plinkGenotype == 1) {
			genvisisGenotype = 2;
		} else if (chromosome == 24 && plinkGenotype == 1) {
			genvisisGenotype = 2;
		} else {
			genvisisGenotype = plinkGenotype;
		}
		return genvisisGenotype;
	}

	/**
	 * Translate Genvisis genotype genotype to PLINK genotype. This is specifically for markers on chromosome X of males.
	 *  
	 * @param genvisisGenotype
	 * @param chromosome
	 * @return
	 */
	public static byte convertGenvisisGenotypeToPlink (byte genvisisGenotype, byte chromosome) {
		byte plinkGenotype;
		if (chromosome == 23 && genvisisGenotype == 2) {
			plinkGenotype = 1;
		} else if (chromosome == 24 && genvisisGenotype == 2) {
			plinkGenotype = 1;
		} else {
			plinkGenotype = genvisisGenotype;
		}
		return plinkGenotype;
	}

	/**
	 * Convert the array of genotypes into an array of PLINK .bed byte stream.
	 * 
	 * @param genotype
	 * @return
	 */
	public static byte[] encodePlinkBedBytesForASingleMarkerOrSample (byte[] genotype) {
		int iBytes;
		byte[] result;
		int nBytes;
		byte shift;

		nBytes = (int) Math.ceil((double) genotype.length / 4);
		iBytes = -1;
		result = new byte[nBytes];

//		for (int i = 0; i < result.length; i++) {
//			result[i] = (byte) 0xAA;	//initilize the array to be 0b10101010, the null genotype defined by PLINK bed data.
//		}

		try {
			for (int i = 0; i < genotype.length; i++) {
				shift = (byte) ((i % 4) * 2);
				if (shift == 0) {
					iBytes ++;
				}
				result[iBytes] = (byte) ((result[iBytes] & (~(0x03 << shift))) | (encodeLastTwoBitsOfABedByte(genotype[i]) << shift));
//				displayBits(result[iBytes]);
			}
		} catch (Elision e) {
			e.printStackTrace();
		}
		
		return result;
	}

	/**
	 * This is part of the method encodePlinkBedBytesForASingleMarkOrSamp(byte[] )
	 * 
	 * @param genotype
	 * @return
	 * @throws Elision
	 */
	public static byte encodeLastTwoBitsOfABedByte (byte genotype) throws Elision {
		byte bedByte;

		if (genotype == (byte) 0) {
			bedByte = (byte) 0x00;

		} else if (genotype == (byte) 1) {
			bedByte = (byte) 0x02;

		} else if (genotype == (byte) 2) {
			bedByte = (byte) 0x03;

		} else if (genotype == (byte) -1) {
			bedByte = (byte) 0x01;

		} else {
			throw new Elision("Unrecognized genotype: " + genotype + ". Please use 0 for A/A, 1 for A/B, 2 for B/B, and -1 for null.");
		}

		return bedByte;
	}

	/**
	 * Convert an array of byte stream from a PLINK .bed file into an array of genotypes with each element corresponding to one single sample
	 * @param bedBytes
	 * @param startIndex
	 * @param indicesOfSamplesOrMarkers
	 * @return
	 */
	public static byte[] decodeBedBytesOfASingleMarkOrSamp (byte[] bedBytes, int startIndex, int[] indicesOfSamplesOrMarkers) {
		byte[] genotypes;
		int indexBedBytes;
		int indexBedByte;

		genotypes = new byte[indicesOfSamplesOrMarkers.length];
		try {
			for (int i = 0; i <= genotypes.length; i++) {
				indexBedBytes = indicesOfSamplesOrMarkers[i] / 4;
				indexBedByte = indicesOfSamplesOrMarkers[i] % 4;
				genotypes[i] = decodeLastTwoBitsOfABedByte((byte) (bedBytes[indexBedBytes] >> (indexBedByte * 2)));
			}
		} catch (Elision e) {
			e.printStackTrace();
		}

		return genotypes;
	}

	/**
	 * Convert an array of byte stream from a PLINK .bed file to an array of genotypes with each element corresponding to one single sample
	 * @param bedBytes
	 * @param startIndex
	 * @param nSamplesOrMarkers
	 * @return
	 */
	public static byte[] decodeBedBytesOfASingleMarkOrSamp (byte[] bedBytes, int startIndex, int nSamplesOrMarkers) {
		byte[] genotypes;
		int indexSampOrMark;
		int endIndex;

		indexSampOrMark = 0;
		endIndex = (int) Math.ceil((double) nSamplesOrMarkers / 4);
		genotypes = new byte[nSamplesOrMarkers];
		try {
			for (int i = startIndex; i <= endIndex; i++) {
				decodeBedByte(bedBytes[i], genotypes, indexSampOrMark);
				indexSampOrMark += 4;
			}
		} catch (Elision e) {
			e.printStackTrace();
		}

		return genotypes;
	}

	/**
	 * Convert a single byte from a PLINK .bed file to the genotype of up to 4 samples and put them into the output array. 
	 * @param inputOneByteFromBed
	 * @param outputGenotypes
	 * @param startIndexOfOutput
	 * @throws Elision
	 */
	public static void decodeBedByte (byte inputOneByteFromBed, byte[] outputGenotypes, int startIndexOfOutput) throws Elision {
		for (int i = 0; i < 4 && startIndexOfOutput < outputGenotypes.length; i++) {
			outputGenotypes[startIndexOfOutput + i] = decodeLastTwoBitsOfABedByte((byte) (inputOneByteFromBed >> (2 * i)));
			startIndexOfOutput ++;
		}
	}

	/**
	 * Convert a type from a PLINK .bed file to an array of genotypes with each element corresponding to a single sample.
	 * @param bedByte
	 * @return
	 * @throws Elision
	 */
	public static byte[] decodeBedByte(byte bedByte) throws Elision {
		byte[] genotypes;

		genotypes = new byte[4];
		for (int k = 0; k < 4; k++) {
			genotypes[k] = decodeLastTwoBitsOfABedByte((byte) (bedByte >> (2 * k)));
		}

		return genotypes;
	}

	/**
	 * This is part of the method decodeBedByte(byte)
	 * @param bedByte
	 * @return
	 * @throws Elision
	 */
	public static byte decodeLastTwoBitsOfABedByte (byte bedByte) throws Elision {
		byte genotype;

		bedByte = (byte) (bedByte &  0x03);

		if (bedByte == (byte) 0) {
			genotype = (byte) 0;

		} else if (bedByte == (byte) 2) {
			genotype = (byte) 1;

		} else if (bedByte == (byte) 3) {
			genotype = (byte) 2;

		} else if (bedByte == (byte) 1) {
			genotype = (byte) -1;

		} else {
			throw new Elision("Unrecognized genotype: " + bedByte);
		}

		return genotype;
	}

	/**
	 * A utility to show the bit maps of several bytes read in from a PLINK .bed file.
	 * @param bedSetDirAndFilenameRoot
	 * @param nSamplesToLoad
	 */
	public static void showBitMapOfBedFileWithPedLayout(String bedSetDirAndFilenameRoot, int indexOfStartMarker, int nMarkersToLoad, int indexOfStartSample, int nSamplesToLoad, Logger log) {
		RandomAccessFile reader;
		Vector<String> fam; 
		byte[] readBuffer;
		byte[][] readBufferSnpMajor;
		int totalMarkers;
		int bytesPerSample;
		int bytesPerMarker;
		long location;
		byte index;
		String out;
		
		log.report("Bitmap of " + bedSetDirAndFilenameRoot);
		try {
			nSamplesToLoad = Math.max(nSamplesToLoad, 0);
			fam = loadFamOrBim(bedSetDirAndFilenameRoot + ".fam", indexOfStartSample, nSamplesToLoad, log);
			totalMarkers = loadFamOrBim(bedSetDirAndFilenameRoot + ".bim", indexOfStartSample, nSamplesToLoad, log).size();

			reader = new RandomAccessFile(bedSetDirAndFilenameRoot + ".bed", "r");
			readBuffer = new byte[3];
			reader.read(readBuffer);
			indexOfStartSample = Math.max(indexOfStartSample, 0);
			if (readBuffer[2] == 1) {	// SNP-major
				bytesPerMarker = (int) Math.ceil((float)fam.size() / 4);
				if (nMarkersToLoad <= 0) {
					nMarkersToLoad = totalMarkers - indexOfStartMarker;
				}
				readBufferSnpMajor = new byte[nMarkersToLoad][bytesPerMarker];
				location = 3 + indexOfStartMarker * bytesPerMarker;
				reader.seek(location);
				for (int i = 0; i < nMarkersToLoad; i++) {
					reader.read(readBufferSnpMajor[i]);
//					location += bytesPerMarker;
				}
				nSamplesToLoad = fam.size();
				for (int i = 0; i < nSamplesToLoad; i++) {
					out = fam.elementAt(indexOfStartSample);
					location = indexOfStartSample / 4;
					index = (byte) ((indexOfStartSample % 4) * 2);
					for (int j = 0; j < nMarkersToLoad; j++) {
						out += showBitsLittleEndian(readBufferSnpMajor[j][(int) location], index);
					}
					log.report(out);
					indexOfStartSample ++;
				}

			} else if (readBuffer[2] == 0) {	// individual-major
				bytesPerSample = (int) Math.ceil((float)totalMarkers / 4);
				readBuffer = new byte[bytesPerSample];
				location = 3 + Math.max(indexOfStartSample, 0) * bytesPerSample;
				while (nSamplesToLoad > 0) {
					reader.seek(location);
					reader.read(readBuffer);
					showBitMapOfByteArrayInGroupsOfTwo(readBuffer, log);
					location += bytesPerSample;
					nSamplesToLoad --;
				}

			} else {
				log.reportError("Error - unrecognized flag at the 3rd byte of the .bed file: " + readBuffer[2] + " (should be either 1 for SNP-major or 0 for individual-major)");

			}

//			if (nMarkersToLoad < 0) {
//				nMarkersToLoad = 10;
//			}


		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * A utility to show the bit maps of several bytes read in from a PLINK .bed file.
	 * 
	 * @param bedFileName
	 * @param startByte with 0 meaning the 1st byte, 1 the 2nd, and etc. 
	 * @param nBytes
	 * @param log
	 */
	public static void showBitMapOfBedFile(String bedFileName, int startByte, int nBytes, Logger log) {
		RandomAccessFile in;
		byte[] readBuffer;
		
		log.report("Bitmap of " + bedFileName);
		try {
			in = new RandomAccessFile(bedFileName, "r");
			if (nBytes < 3) {
				nBytes = 10;
			}
			if (startByte >= 3) {
				readBuffer = new byte[3];
				in.read(readBuffer);
				showBitMapOfByteArray(readBuffer, 0, log);
				readBuffer = new byte[(int) Math.min(nBytes, in.length())];
				in.seek(startByte);
			} else {
				readBuffer = new byte[(int) Math.min(nBytes, in.length())];
			}
			in.read(readBuffer);
			showBitMapOfByteArray(readBuffer, startByte, log);


		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * A utility to show the bit map of a byte array;
	 * @param data
	 * @param startByte
	 * @param nBytes
	 * @param log
	 */
	public static void showBitMapOfByteArrayInGroupsOfTwo(byte[] data, Logger log) {
		for (int i = 0; i < data.length; i++) {
			log.report(showBitsLittleEndian(data[i]));
		}
	}

	/**
	 * A utility to show the bit map of a byte.
	 * @param data
	 */
	public static String showBitsLittleEndian(byte data) {
		String out;

		out = "";
		for (int j = 0; j < 4; j++) {
			out += ("\t" + (((data & 0x80) >> (j*2)) & 0x03));
		}

		return out;
	}

	/**
	 * A utility to show the bit map of a byte.
	 * @param data
	 */
	public static String showBitsLittleEndian(byte data, byte index) {
		return "\t" + ((data >> (index + 1)) & 0x01) + ((data >> index) & 0x01);
	}

	/**
	 * A utility to show the bit map of a byte array;
	 * @param data
	 * @param startByte
	 * @param nBytes
	 * @param log
	 */
	public static void showBitMapOfByteArray(byte[] data, int indexOfTheFirstElement, Logger log) {
		for (int i = 0; i < data.length; i++) {
			log.report((indexOfTheFirstElement + i) + ":\t" + showBits(data[i]));
		}
	}

	/**
	 * A utility to show the bit map of a byte.
	 * @param data
	 */
	public static String showBits(byte data) {
		String out;

		out = "";
		for (int j = 0; j < 8; j++) {
			out += ((j==4? " " : "") + ((data & 0x80) >> 7));
			data = (byte) (data << 1);
		}

		return out;
	}

	/**
	 * Generate MarkerLookup file from PLINK data sets. This is part of the process to convert PLINK data sets to Genvisis data.
	 * @param plinkFileRoot
	 * @return
	 */
	public static MarkerLookup parseMarkerLookup(String plinkFileRoot) {
		BufferedReader reader;
		String[] line;
		Hashtable<String, String> hash;
		int count;
		
		hash = new Hashtable<String, String>();
		try {
			reader = new BufferedReader(new FileReader(plinkFileRoot+".bim"));
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("\\s+");
				hash.put(line[1], ":\t"+count+"\t"+line[0] + "\t" + line[3] + "\t" + line[4]+"\t"+line[5]);
				count++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + plinkFileRoot+".bim" + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + plinkFileRoot+".bim" + "\"");
			System.exit(2);
		}

		return new MarkerLookup(hash);
	}

	/**
	 * Match samples in a PLINK data set with samples in a Genvisis project. This is part of the process to convert PLINK data sets to Genvisis data.
	 * @param proj
	 * @param plinkFileRoot
	 * @param log
	 * @return
	 */
	public static int[] parseSampleIndicesForProject(Project proj, String plinkFileRoot) {
		BufferedReader reader;
		String[] line;
		int count;
		String[] finalSampleIDs, allIDs;
		String famIndID, sampleID;
		SampleData sampleData;
		int[] sampleIndices;
		
		finalSampleIDs = proj.getSamples();
		sampleData = proj.getSampleData(SampleData.BASIC_CLASSES.length, false);
		
		sampleIndices = Array.intArray(finalSampleIDs.length, -1);
		try {
			reader = new BufferedReader(new FileReader(plinkFileRoot+".fam"));
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("\\s+");
				famIndID = line[0]+"\t"+line[1];
				allIDs = sampleData.lookup(famIndID);
				if (allIDs == null) {
					proj.getLog().report("Warning - sample in PLINK file "+plinkFileRoot+".fam that is not in the project's sampleData file: "+famIndID);
				} else {
					sampleID = allIDs[0];
					sampleIndices[ext.indexOfStr(sampleID, finalSampleIDs)] = count;
				}
				count++;				
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + plinkFileRoot+".fam" + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + plinkFileRoot+".fam" + "\"");
			System.exit(2);
		}
		
		return sampleIndices;
	}

	/**
	 * Generate the list of sample indices for a PLINK data set.
	 * This is part of the process to convert PLINK data to Genvisis data.
	 * 
	 * @param plinkFileRoot
	 * @param log
	 * @return
	 */
	public static int[] parseSampleIndicesAll(String plinkFileRoot, Logger log) {
		return Array.intArray(Files.countLines(plinkFileRoot+"fam", 0));
	}

	/**
	 * (This method is dropped, because of a conceptual error - PLINK data do not have intensity x and y)
	 * 
	 * This is a method to load data from a PLINK data set into MarkerData[]. This is part of the process to conver PLINK data into Genvisis data.
	 * 
	 * @param allMarkersInProj
	 * @param allChrsInProj
	 * @param allPositionsInProj
	 * @param allSamplesInProj
	 * @param bedFileName
	 * @param markersIndicesInProj
	 * @param markersIndicesInFile
	 * @param sampleFingerprint
	 * @param sampIndices
	 * @return
	 */
	@SuppressWarnings("resource")
	public static MarkerData[] loadBedUsingRAF(String[] allMarkersInProj, byte[] allChrsInProj, int[] allPositionsInProj, String[] allSamplesInProj, String bedFileName, int[] markersIndicesInProj, int[] markersIndicesInFile, long sampleFingerprint, int[] sampIndices) {
		MarkerData[] result;
		RandomAccessFile in;
		int nSampTotallyInBed;
		int nBytesPerMarkInBed;
		int indexBedBytes;
		int indexBedBits;
		byte[] bytesOfOneMarkerInBed;
		byte[] genotypes;

		result = new MarkerData[markersIndicesInProj.length];
		try {
			in = new RandomAccessFile(bedFileName, "r");

			//Test the PLINK file header && Marker Dominant or Sample Dominant?
			byte[] isSnpMajor = new byte[3];
			in.read(isSnpMajor);
			if (isSnpMajor[2] == 0) {
				// sample dominant
				System.err.println("Error - .bed file is sample-dominant.");
			} else {
				
			}

			// Sort the markerList for sequential loading

			// Load from bim and fam
			nSampTotallyInBed = 0;
			nBytesPerMarkInBed = (int) Math.ceil((double) nSampTotallyInBed / 4);

			// Load the data
			bytesOfOneMarkerInBed = new byte[nBytesPerMarkInBed];
			genotypes = new byte[sampIndices.length];
			for (int i = 0; i < markersIndicesInProj.length; i++) {
				// Should we do sequential reading?
				in.seek(3 + i * nBytesPerMarkInBed);
				in.read(bytesOfOneMarkerInBed);
				for (int j = 0; j < sampIndices.length; j++) {
					indexBedBytes = sampIndices[j] / 4;
					indexBedBits = sampIndices[i] % 4;
					genotypes[j] = decodeLastTwoBitsOfABedByte((byte) (bytesOfOneMarkerInBed[indexBedBytes] >> (indexBedBits * 2)));
				}

		        result[i] = new MarkerData(allMarkersInProj[markersIndicesInProj[i]], allChrsInProj[markersIndicesInProj[i]], allPositionsInProj[markersIndicesInProj[i]], 0l, null, null, null, null, null, null, null, null, null, genotypes, null);
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (Elision e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return result;
	}

//	public static String loadBedFile(String bedFileName, int[] markersIndicesInFile, int[] sampIndices) {
//		String result;
//		RandomAccessFile in;
//		int nSampTotallyInBed;
//		int nBytesPerMarkInBed;
//		int indexBedBytes;
//		int indexBedBits;
//		byte[] bytesOfOneMarkerInBed;
//		byte[] genotypes;
//
//		try {
//			in = new RandomAccessFile(bedFileName, "r");
//			//Test the PLINK file header && Marker Dominant or Sample Dominant?
//
//			// Sort the markerList for sequential loading
//
//			// Load from bim and fam
//			nSampTotallyInBed = 0;
//			nBytesPerMarkInBed = (int) Math.ceil((double) nSampTotallyInBed / 4);
//
//			// Load the data
//			bytesOfOneMarkerInBed = new byte[nBytesPerMarkInBed];
//			genotypes = new byte[sampIndices.length];
//			for (int i = 0; i < markersIndicesInProj.length; i++) {
//				// Should we do sequential reading?
//				in.seek(3 + i * nBytesPerMarkInBed);
//				in.read(bytesOfOneMarkerInBed);
//				for (int j = 0; j < sampIndices.length; j++) {
//					indexBedBytes = sampIndices[j] / 4;
//					indexBedBits = sampIndices[i] % 4;
//					genotypes[j] = decodeLastTwoBitsOfABedByte((byte) (bytesOfOneMarkerInBed[indexBedBytes] >> (indexBedBits * 2)));
//				}
//
//		        result[i] = new MarkerData(allMarkersInProj[markersIndicesInProj[i]], allChrsInProj[markersIndicesInProj[i]], allPositionsInProj[markersIndicesInProj[i]], 0l, null, null, null, null, null, null, null, null, null, genotypes, null);
//			}
//
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		} catch (Elision e) {
//			e.printStackTrace();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//
//		return result;
//	}

	public static void main (String[] args) {
		String conversionToRun = null;
		String projPropertyFileFullPath = null;
		Project proj;
		String plinkDataDirAndFilenameRoot = "plink";
		boolean isSnpMajor = true;
		int startByte = 100;
		int nBytes = -1;
		float gcThreshold = 0.15f;
		int indexOfStartMarker = 0;
		int nMarkers = -1;
		int indexOfStartSample = 0;
		int nSamples = -1;
		Logger log;

		String usage = "\n"+
				"cnv.manage.PlinkData provides methods to do the following:\n" +
				"   (1) Conversion from Genvisis data to PLINK text (.ped) data set;\n" +
				"   (2) Conversion from Genvisis data to PLINK binary (.bed) data set with SNP major;\n" +
				"   (3) Conversion from Genvisis data to PLINK binary (.bed) data set with individual major;\n" +
				"   (4) Conversion from PLINK text (.ped) data set to PLINK binary (.bed) data set with SNP major;\n" +
				"   (5) Conversion from PLINK text (.ped) data set to PLINK binary (.bed) data set with individual major;\n" +
				"   (6) Conversion from PLINK binary (.bed) data set with SNP major to PLINK text (.ped) data set;\n" +
				"   (7) Show the bit map of a PLINK binary (.bed) data set;\n" +
				"   (8) Read in genotype from PLINK text (.ped) data set (not available in command line mode);\n" +
				"   (9) Read in genotype from PLINK binary (.bed) data set (not available in command line mode);\n" +
				"\n" +
				"To export from Genvisis to PLINK text (.ped) data set, the following arguments are required:\n" +
				"   (1) the command (i.e. -genvisisToPed (not the default));\n" +
				"   (2) the Genvisis project's property file location (i.e. proj=~/projects/default.properties (not the default));\n" +
				"   (3) the GC threshold (i.e. gcthreshold=" + gcThreshold + " (default));\n" +
				"   (4) the directory and filename root for the final PLINK data set (i.e. plinkdata=" + plinkDataDirAndFilenameRoot + " (default));\n" +
				"\n" +
				"To export from Genvisis to PLINK binary (.bed) data set, the following arguments are required:\n" +
				"   (1) the command (i.e. -genvisisToBed (not the default));\n" +
				"   (2) the Genvisis project's property file location (i.e. proj=" + projPropertyFileFullPath + " (default));\n" +
				"   (3) the GC threshold (i.e. gcthreshold=" + gcThreshold + " (default));\n" +
				"   (4) is the .bed file going to be SNP Major (i.e. issnpmajor=" + isSnpMajor + " (default));\n" +
				"   (5) the directory and filename root for the final PLINK data set (i.e. plinkdata=" + plinkDataDirAndFilenameRoot + " (default));\n" +
				"Note: the following specified by the project's property file are also required:" +
				"   (6) the text file \"TARGET_MARKERS_FILENAME\";\n" +
				"   (7) the text file \"PEDIGREE_FILENAME\";\n" +
				"\n" +
				"To convert PLINK text (.ped) data to PLINK binary (.bed) data, the following arguments are required:\n" +
				"   (1) the command (i.e. -pedToBed (not the default));\n" +
				"   (2) the directory and filename root for the PLINK data, which will be the same for both input/output (i.e. plinkdata=" + plinkDataDirAndFilenameRoot + " (default));\n" +
				"   (3) is the .bed file going to be SNP Major (i.e. issnpmajor=" + isSnpMajor + " (default));\n" +
				"\n" +
				"To convert PLINK binary (.bed) data to PLINK text (.ped) data, the following arguments are required:\n" +
				"   (1) the command (i.e. -bedToPed (not the default));\n" +
				"   (2) the directory and filename root for the PLINK data, which will be the same for both input/output (i.e. plinkdata=" + plinkDataDirAndFilenameRoot + " (default));\n" +
				"\n" +
//				"To show bitmap of PLINK binary (.bed) data, the following arguments are required:\n" +
//				"   (1) the command (i.e. -bitmap (not the default));\n" +
//				"   (2) the PLINK binary (.bed) data set's directory and filename root (i.e. plinkdata=" + plinkDataDirAndFilenameRoot + " (default));\n" +
//				"   (3) the starting byte in .bed file to show (i.e. startbyte=" + startByte + " (default));\n" +
//				"   (4) the number of bytes in .bed file to show (i.e. nbytes=" + nBytes + " (default, with a negative number meaning all));\n" +
//				"\n" +
//				"To show bitmap of PLINK binary (.bed) data in a layout similiar to .ped file, the following arguments are required:\n" +
//				"   (1) the command (i.e. -bitmapInPedLayout (not the default));\n" +
//				"   (2) the PLINK data's directory and filename root for both the text set and binary set (i.e. plinkdata=" + plinkDataDirAndFilenameRoot + " (default));\n" +
//				"   (3) the index of the starting marker in the PLINK binary (.bed) data set to show (i.e. startmarker=" + indexOfStartMarker + " (default));\n" +
//				"   (4) the number of markers in the PLINK binary (.bed) data set to show (i.e. nmarkers=" + nMarkers + " (default, with a negative number meaning all));\n" +
//				"   (5) the index of the starting sample in the PLINK binary (.bed) data set to show (i.e. startsample=" + indexOfStartSample + " (default));\n" +
//				"   (6) the number of samples in the PLINK binary (.bed) data set to show (i.e. nsamples=" + nSamples + " (default, with a negative number meaning all));"+
				"\n";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].equalsIgnoreCase("-genvisisToPed")) {
				conversionToRun = args[i];
			} else if (args[i].equalsIgnoreCase("-genvisisToBed")) {
				conversionToRun = args[i];
			} else if (args[i].equalsIgnoreCase("-pedToBed")) {
				conversionToRun = args[i];
			} else if (args[i].equalsIgnoreCase("-bedToPed")) {
				conversionToRun = args[i];
			} else if (args[i].startsWith("-bitmap")) {
				conversionToRun = args[i];
			} else if (args[i].startsWith("-bitmapInPedLayout")) {
				conversionToRun = args[i];
			} else if (args[i].startsWith("proj=")) {
				projPropertyFileFullPath = args[i].split("=")[1];
			} else if (args[i].startsWith("plinkdata=")) {
				plinkDataDirAndFilenameRoot = args[i].split("=")[1];
			} else if (args[i].startsWith("issnpmajor=")) {
				isSnpMajor = Boolean.parseBoolean(args[i].split("=")[1]);
			} else if (args[i].startsWith("gcthreshold=")) {
				gcThreshold = Float.parseFloat(args[i].split("=")[1]);
			} else if (args[i].startsWith("startbyte=")) {
				startByte = Integer.parseInt(args[i].split("=")[1]);
			} else if (args[i].startsWith("nbytes=")) {
				nBytes = Integer.parseInt(args[i].split("=")[1]);
			} else if (args[i].startsWith("startmarker=")) {
				indexOfStartMarker = Integer.parseInt(args[i].split("=")[1]);
			} else if (args[i].startsWith("nmarkers=")) {
				nMarkers = Integer.parseInt(args[i].split("=")[1]);
			} else if (args[i].startsWith("startsample=")) {
				indexOfStartSample = Integer.parseInt(args[i].split("=")[1]);
			} else if (args[i].startsWith("nsamples=")) {
				nSamples = Integer.parseInt(args[i].split("=")[1]);
			} else {
				System.err.println("Error - invalid argument: "+args[i]);
			}
		}

		if  (conversionToRun == null) {
			log = new Logger();
			log.report(usage);

		} else if (conversionToRun.equals("-genvisisToPed")) {
			proj = new Project(projPropertyFileFullPath, false);
			log = proj.getLog();
			log.report(ext.getTime()+"\tConverting Genvisis to PLINK text (.ped) data set.");
			PlinkData.saveGenvisisToPlinkPedSet(proj, plinkDataDirAndFilenameRoot, proj.DATA_DIRECTORY.getValue(false, true) + proj.CLUSTER_FILTER_COLLECTION_FILENAME, null);

		} else if (conversionToRun.equals("-genvisisToBed")) {
			proj = new Project(projPropertyFileFullPath, false);
			log = proj.getLog();
			log.report(ext.getTime()+"\tConverting from Genvisis to PLINK binary (.bed) data set.");
			saveGenvisisToPlinkBedSet(proj, plinkDataDirAndFilenameRoot, proj.DATA_DIRECTORY.getValue(false, true) + proj.CLUSTER_FILTER_COLLECTION_FILENAME, proj.TARGET_MARKERS_FILENAMES.getValue()[0], gcThreshold, true);
			
		} else if (conversionToRun.equals("-pedToBed")) {
			log = new Logger(ext.parseDirectoryOfFile(plinkDataDirAndFilenameRoot) + "PlinkData_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
			log.report(ext.getTime()+"\tConverting from PLINK text (.ped) data set to PLINK binary (.bed) data set.");
			convertPedSetToBedSet(plinkDataDirAndFilenameRoot, isSnpMajor, log);

		} else if (conversionToRun.equals("-bedToPed")) {
			log = new Logger(ext.parseDirectoryOfFile(plinkDataDirAndFilenameRoot) + "PlinkData_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
			log.report(ext.getTime()+"\tConverting from PLINK binary (.bed) data set to PLINK text (.ped) data set.");
			convertBedSetToPedSet(plinkDataDirAndFilenameRoot, log);

		} else if (conversionToRun.equals("-bitmap")) {
			log = new Logger(ext.parseDirectoryOfFile(plinkDataDirAndFilenameRoot) + "PlinkData_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
			log.report(ext.getTime()+"\tShowing bitmap.");
			showBitMapOfBedFile(plinkDataDirAndFilenameRoot + ".bed", startByte, nBytes, log);

		} else if (conversionToRun.equals("-bitmapInPedLayout")) {
			log = new Logger(ext.parseDirectoryOfFile(plinkDataDirAndFilenameRoot) + "PlinkData_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
			log.report(ext.getTime()+"\tShowing bitmap.");
			showBitMapOfBedFileWithPedLayout(plinkDataDirAndFilenameRoot, indexOfStartMarker, nMarkers, indexOfStartSample, nSamples, log);

		} else {
			log = new Logger();
			log.report("Unrecognized command conversion=" + conversionToRun);
		}

		log.report(ext.getTime()+"\tFinished PlinkData.java.");
		
	}
}
