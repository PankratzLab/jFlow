package cnv.manage;

import java.io.*;
import java.util.*;
import common.*;

import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerLookup;
import cnv.filesys.Project;
import cnv.filesys.SampleList;
import cnv.var.SampleData;

public class PlinkData {
	public static final String FAM_DELIMITER = " ";
	public static final String BIM_DELIMITER = "\t";
	public static final String MAP_DELIMITER = "\t";

	public static void convertPedToBed(String plinkDataDir, String pedFileNameStem, boolean isMarkDominant) {
//		RandomAccessFile outBed;
//		RandomAccessFile outFam;
		BufferedOutputStream outBed;
		PrintWriter outFamOrBim;
		Scanner inFile;
		int nMarks;
		String[] line;
		byte genotype;
		byte[] genotypeByteStream = null;
		String[][] genotypeLetters = null;
		int index;
		int index2;
		byte[] alleles;
		String temp;

		try {
			genotype = (byte) 0;
			nMarks = -1;
			alleles = new byte[2];
			inFile = new Scanner(new FileInputStream(plinkDataDir + pedFileNameStem + ".ped"));
//			outBed = new RandomAccessFile(plinkDataDir + plinkPedFileNameStem + ".bed", "rw");
//			outBed.write(new byte[] {(byte) 108, (byte) 27, (byte) 1});
			outBed = new BufferedOutputStream(new FileOutputStream(plinkDataDir + pedFileNameStem + ".bed"));
			outBed.write(new byte[] {(byte) 108, (byte) 27, (byte) 1});
//			outFam = new RandomAccessFile(plinkDataDir + plinkPedFileNameStem + ".fam", "rw");
			outFamOrBim = new PrintWriter(new FileOutputStream(plinkDataDir + pedFileNameStem + ".fam"));

			while (inFile.hasNext()) {
				line = inFile.nextLine().split("[\t\\s]");
				if (nMarks == -1) {
					nMarks = (line.length - 6) / 2;
					genotypeLetters = new String[nMarks][2];
					for (int i = 0; i < genotypeLetters.length; i++) {
						genotypeLetters[i][0] = "0";
						genotypeLetters[i][1] = "0";
					}
					genotypeByteStream = new byte[(int) Math.ceil((double) nMarks / 4)];
				}
				outFamOrBim.println(line[0] + FAM_DELIMITER + line[1] + FAM_DELIMITER + line[2] + FAM_DELIMITER + line[3] + FAM_DELIMITER + line[4] + FAM_DELIMITER + line[5]);

				index = 0;
				index2 = 0;
				genotypeByteStream[genotypeByteStream.length - 1] = (byte) 0;
				for (int i = 0; i < nMarks; i++) {
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
			inFile = new Scanner(new FileInputStream(plinkDataDir + pedFileNameStem + ".map"));
			outFamOrBim = new PrintWriter(new FileOutputStream(plinkDataDir + pedFileNameStem + ".bim"));
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

//	public static void convertBedToPed(String plinkDataDir, String bedFileNameStem) {
////		BufferedInputStream inBed;
//		RandomAccessFile inBed;
//		PrintWriter outFile;
////		Scanner inFamOrBim;
//		Scanner inBimOrFam;
//		int nMarks;
//		int nSamps;
//		int nBedBytesPerMarker;
//		byte[] genotypeByteStream;
//		String[] line;
//		byte tmp;
//		byte[] genotypes;
//		Vector<String[]> genotypeLetters;
//		Vector<String> famColumns;
//		int index;
//		int index2;
//		byte[] alleles;
//		String temp;
//
//		try {
//			genotypeLetters = new Vector<String[]>();
//			inBimOrFam = new Scanner(new FileInputStream(plinkDataDir + bedFileNameStem + ".bim"));
//			outFile = new PrintWriter(new FileOutputStream(plinkDataDir + bedFileNameStem + ".map"));
//			while(inBimOrFam.hasNext()) {
//				line = inBimOrFam.nextLine().split("[\t]");
//				outFile.println(line[0] + BIM_DELIMITER + line[1] + BIM_DELIMITER + line[2] + BIM_DELIMITER + line[3]);
//				genotypeLetters.add(new String[] {line[4], line[5]});
//			}
//			inBimOrFam.close();
//			outFile.close();
//
//			famColumns = new Vector<String>();
//			inBimOrFam = new Scanner(new FileInputStream(plinkDataDir + bedFileNameStem + ".fam"));
//			while (inBimOrFam.hasNext()) {
//				famColumns.add(inBimOrFam.nextLine());
//			}
//			inBimOrFam.close();
//
//			nMarks = genotypeLetters.size();
//			nSamps = famColumns.size();
//			nBedBytesPerMarker = (int) Math.ceil((double) nSamps / 4);
//
//			inBed = new RandomAccessFile(plinkDataDir + bedFileNameStem + ".bed", "r");
//			outFile = new PrintWriter(new FileOutputStream(plinkDataDir + bedFileNameStem + ".ped"));
//			for (int i = 0; i < nSamps; i++) {
//				temp = famColumns.elementAt(i);
//				for (int j = 0; j < nMarks; j++) {
//					inBed.seek(3 + (long) i * nSamps * nBedBytesPerMarker);
//					tmp = inBed.readByte();
//					temp += ;
//				}
//				outFile.println(temp);
//			}
//
//
//			inBed.close();
//			outFile.close();
//
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}
//
//	public static void loadPedData(String plinkDataDir, String bedFileNameStem, int[] indeciesOfMarkersToLoad, int[] indeciesOfSamplesToLoad) {
//		Scanner inFamOrBim;
//		String[] line;
//		Vector<Byte> chr;
//		Vector<String> markerName;
//		Vector<Double> morgans;
//		Vector<Integer> pos;
//		int nMarks;
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
//			nMarks = markerName.size();
//
//			famId = new Vector<String>();
//			sampId = new Vector<String>();
//			paternalId = new Vector<String>();
//			maternalId = new Vector<String>();
//			sex = new Vector<Byte>();
//			phenotype = new Vector<String>();
//			genotypes = new byte[nMarks];
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
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//	}
//
//	public static void dumpPlinkBinaryData(String PlinkDataDir) {
//	}
//
//	public static void loadBedData(String plinkDataDir, String bedFileNameStem, int[] indeciesOfMarkersToLoad, int[] indeciesOfSamplesToLoad) {
//		Scanner inFamOrBim;
//		String[] line;
//		Vector<String> sampId;
//		Vector<Byte> sex;
//		Vector<String> cluster;
//		BufferedInputStream inBed;
//
//		int nSamps;
//		int nMarks;
//		int iterations;
//		int index;
//		byte[] genotypes;
//		byte buffer;
//		ClusterFilterCollection clusterFilterCollection;
//		MarkerData[] markData;
//
//		sampId = new Vector<String>();
//		sex = new Vector<Byte>();
//		cluster = new Vector<String>();
//		inFamOrBim = new Scanner(new FileInputStream(plinkDataDir + bedFileNameStem + ".bed"));
//		while (inFamOrBim.hasNext()) {
//			line = inFamOrBim.nextLine().split("\t");
//			sampId.add(line[0] + "\t" + line[1]);
//			sex.add(line[4].equals("1")? (byte) 0 : (byte) 1);
//			cluster.add(line[5]);
//		}
//
//		inBed = new BufferedInputStream(new FileInputStream(plinkDataDir + bedFileNameStem + ".bed"));
//
//		if (markerList == null) {
//			nMarks = proj.getSamples().length;
//			markerList = proj.getMarkerNames();
//		} else {
//			nMarks = markerList.length;
//		}
//		nSamps = proj.getSamples().length;
//		clusterFilterCollection = ClusterFilterCollection.load(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, true), proj.getJarStatus()); 
//		
//		try {
//			in = new RandomAccessFile(PlinkBinaryDataDir + "plink.bed", "rw");
//			outStream = new byte[1];
//			outStream[0] = (byte) 108;	// 0b01101100
//			outStream[1] = (byte) 27;	// 0b00011011
//			outStream[2] = (byte) 1;	//0b00000001
//	
//			iterations = (int) Math.ceil((double) nSamps / 4);
//			markData = new MarkerDataLoader(proj, markerList, 0).markerData;
//			for (int i = 0; i < markerList.length; i++) {
//				genotypes = markData[i].getAbGenotypesAfterFilters(clusterFilterCollection, markerList[i], (float) gcThreshold / 100);
//				index = 0;
//				for (int j = 0; j < iterations; j++) {
//					buffer = (byte) 0;
//					for (int k = 0; k < 4 && index < nSamps; k++) {
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
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//	}
//
//	public static void generateBedFile(Project proj, String[] markerList, String clusterFilterFileName, int gcThreshold) {
//		RandomAccessFile out;
//		byte[] outStream, genotypes;
//		ClusterFilterCollection clusterFilterCollection;
//		MarkerDataLoader markDataLoader;
//		
//		if (markerList == null) {
//			proj.getSamples().length;
//			markerList = proj.getMarkerNames();
//		} else {
//		}
//
//		clusterFilterCollection = ClusterFilterCollection.load(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, true), proj.getJarStatus()); 
//		
//		try {
//			out = new RandomAccessFile("plink.bed", "rw");
//			outStream = new byte[1];
//			outStream[0] = (byte) 108;	// 0b01101100
//			outStream[1] = (byte) 27;	// 0b00011011
//			outStream[2] = (byte) 1;	// 0b00000001
//			out.write(outStream);
//	
//			markDataLoader = new MarkerDataLoader(proj, markerList, 0);
//			for (int i = 0; i < markerList.length; i++) {
//				genotypes = markDataLoader.getMarkerData(i).getAbGenotypesAfterFilters(clusterFilterCollection, markerList[i], 0);
//				out.write(encodePlinkBedBytesForASingleMarkOrSamp(genotypes));
//			}
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}


	/**
	 * Transpose individual-major byte stream to SNP-major one
	 */
	public static void transposeBedByte(byte[] inputSampBytes, int iStartSamp, int iEndSamp, byte[] outputMarkBytes, int nAllMarksInProj, int nAllSampsInProj) {
//		for (int i = 0; i < originalBytes.length; i++) {
//			for (int j = 0; j < transposedBytes.length; j++) {
//				transposedBytes[j] = originalBytes[i];
//			}
//		}

		int iSampStream1;
		int iBytes2;
		int iSampStream3;
		int iBytes4;
		int iMarkStream;
		int bytesPerMark;

		bytesPerMark = (int) Math.ceil((double) nAllSampsInProj / 4);

		iSampStream1 = (int) Math.ceil((double) iStartSamp / 4);
		iBytes2 = iStartSamp % 4;
		iSampStream3 = (int) Math.ceil((double) iEndSamp / 4);
		iBytes4 = iEndSamp % 4;
		
		for (int i = iSampStream1; i < iSampStream3; i++) {
			for (int j = iBytes2; j <= 4; j++) {
				iMarkStream = i * bytesPerMark + j / 4;
				outputMarkBytes[iMarkStream] = (byte) ((outputMarkBytes[iMarkStream] & 0x03) | (inputSampBytes[i] & 0x3));
			}
		}
	}

	public static void transposeOneBedByte(byte originalBytes, int nSampOrMarkInThisByte, byte[] transposedBytes, int iMarksOrSampsInProj, int nAllMarksInProj, int nAllSampsInProj) {
		if(nSampOrMarkInThisByte >= 4) {
			nSampOrMarkInThisByte = 4;
		}

		for (int i = 0; i < nSampOrMarkInThisByte; i++) {
			transposedBytes[iMarksOrSampsInProj / 4 + i] = (byte) (((transposedBytes[i]) & 0x00) | ((originalBytes >> (2 * i)) & 0x00));
		}
	}

	public static byte[] encodePlinkBedBytesForASingleMarkOrSamp (byte[] genotype) {
		int iBytes;
		byte[] result;
		int nBytes;

		nBytes = (int) Math.ceil((double) genotype.length / 4);
		iBytes = -1;
		result = new byte[nBytes];

		for (int i = 0; i < result.length; i++) {
			result[i] = (byte) 0xAA;	//initilize the array to be 0b10101010, the null genotype defined by plink bed data.
		}

		try {
			for (int i = 0; i < genotype.length; i++) {
				if (i % 4 == 0) {
					iBytes ++;
				}
				result[iBytes] = (byte) (result[iBytes] << 2 | encodeLastTwoBitsOfABedByte(genotype[i]));
			}
		} catch (Elision e) {
			e.printStackTrace();
		}
		
		return result;
	}

	public static byte encodeLastTwoBitsOfABedByte (byte genotype) throws Elision {
		byte bedByte;

		if (genotype == (byte) 0) {
			bedByte = (byte) 0x00;

		} else if (genotype == (byte) 1) {
			bedByte = (byte) 0x01;

		} else if (genotype == (byte) 2) {
			bedByte = (byte) 0x03;

		} else if (genotype == (byte) -1) {
			bedByte = (byte) 0x02;

		} else {
			throw new Elision("Unrecognized genotype: " + genotype + ". Please use 0 for A/A, 1 for A/B, 2 for B/B, and -1 for null.");
		}

		return bedByte;
	}

	public static byte[] decodeBedBytesOfASingleMarkOrSamp (byte[] bedBytes, int startIndex, int[] indeciesOfSampsOrMarks) {
		byte[] genotypes;
		int indexBedBytes;
		int indexBedByte;

		genotypes = new byte[indeciesOfSampsOrMarks.length];
		try {
			for (int i = 0; i <= genotypes.length; i++) {
				indexBedBytes = indeciesOfSampsOrMarks[i] / 4;
				indexBedByte = indeciesOfSampsOrMarks[i] % 4;
				genotypes[i] = decodeLastTwoBitsOfABedByte((byte) (bedBytes[indexBedBytes] >> (indexBedByte * 2)));
			}
		} catch (Elision e) {
			e.printStackTrace();
		}

		return genotypes;
	}

	public static byte[] decodeBedBytesOfASingleMarkOrSamp (byte[] bedBytes, int startIndex, int nSampsOrMarks) {
		byte[] genotypes;
		int indexSampOrMark;
		int endIndex;

		indexSampOrMark = 0;
		endIndex = (int) Math.ceil((double) nSampsOrMarks / 4);
		genotypes = new byte[nSampsOrMarks];
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

	public static void decodeBedByte (byte inputOneByteFromBed, byte[] outputGenotypes, int startIndexOfOutput) throws Elision {
		for (int i = 0; i < 4 && startIndexOfOutput < outputGenotypes.length; i++) {
			outputGenotypes[startIndexOfOutput + i] = decodeLastTwoBitsOfABedByte((byte) (inputOneByteFromBed >> (2 * i)));
			startIndexOfOutput ++;
		}
	}

	public static byte[] decodeBedByte(byte bedByte) throws Elision {
		byte[] genotypes;

		genotypes = new byte[4];
		for (int k = 0; k < 4; k++) {
			genotypes[k] = decodeLastTwoBitsOfABedByte((byte) (bedByte >> (2 * k)));
		}

		return genotypes;
	}

	public static byte decodeLastTwoBitsOfABedByte (byte bedByte) throws Elision {
		byte genotype;

		bedByte = (byte) (bedByte &  0x03);

		if (bedByte == (byte) 0) {
			genotype = (byte) 0;

		} else if (bedByte == (byte) 1) {
			genotype = (byte) 1;

		} else if (bedByte == (byte) 3) {
			genotype = (byte) 2;

		} else if (bedByte == (byte) 2) {
			genotype = (byte) -1;

		} else {
			throw new Elision("Unrecognized genotype: " + bedByte);
		}

		return genotype;
	}

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
				line = reader.readLine().trim().split("[\\s]+");
				hash.put(line[1], ":\t"+count+"\t"+line[4]+"\t"+line[5]);
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

	public static int[] parseSampleIndicesForProject(Project proj, String plinkFileRoot, Logger log) {
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
				line = reader.readLine().trim().split("[\\s]+");
				famIndID = line[0]+"\t"+line[1];
				allIDs = sampleData.lookup(famIndID);
				if (allIDs == null) {
					log.report("Warning - sample in plink file "+plinkFileRoot+".fam that is not in the project's sampleData file: "+famIndID);
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

	public static int[] parseSampleIndicesAll(String plinkFileRoot, Logger log) {
		return Array.intArray(Files.countLines(plinkFileRoot+"fam", false));
	}

	@SuppressWarnings("resource")
	public static MarkerData[] loadBedUsingRAF(String[] allMarkersInProj, byte[] allChrsInProj, int[] allPositionsInProj, String[] allSamplesInProj, String bedFileName, int[] markersIndicesInProj, int[] markersIndicesInFile, long sampleFingerprint, int[] sampsIndecies) {
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
			//Test the plink file header && Marker Dominant or Sample Dominant?

			// Sort the markerList for sequential loading

			// Load from bim and fam
			nSampTotallyInBed = 0;
			nBytesPerMarkInBed = (int) Math.ceil((double) nSampTotallyInBed / 4);

			// Load the data
			bytesOfOneMarkerInBed = new byte[nBytesPerMarkInBed];
			genotypes = new byte[sampsIndecies.length];
			for (int i = 0; i < markersIndicesInProj.length; i++) {
				// Should we do sequential reading?
				in.seek(3 + i * nBytesPerMarkInBed);
				in.read(bytesOfOneMarkerInBed);
				for (int j = 0; j < sampsIndecies.length; j++) {
					indexBedBytes = sampsIndecies[j] / 4;
					indexBedBits = sampsIndecies[i] % 4;
					genotypes[j] = decodeLastTwoBitsOfABedByte((byte) (bytesOfOneMarkerInBed[indexBedBytes] >> (indexBedBits * 2)));
				}

//		        result[i] = new MarkerData(allMarkersInProj[markersIndicesInProj[i]], allChrsInProj[markersIndicesInProj[i]], allPositionsInProj[markersIndicesInProj[i]], fingerprint, null, null, null, null, null, null, null, null, null, genotypes, null);
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

	public static void main (String[] args) {
		String fileDir = "N:/statgen/Genvisis/plinkTesting/";
		String fileNameRoot = "plink";

		convertPedToBed(fileDir, fileNameRoot, true);
	}
}
