package cnv.manage;

import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;

import common.*;

import cnv.filesys.ABLookup;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerLookup;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
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
		int numMarkers;
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
			numMarkers = -1;
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

	/**
	 * Unfinished
	 * @param proj
	 * @param markerList
	 * @param clusterFilterFileName
	 * @param gcThreshold
	 */
	public static void saveToPedFile(Project proj, String[] markerList, String clusterFilterFileName, int gcThreshold, String pedFilenameRoot) {
		PrintWriter out;
		ClusterFilterCollection clusterFilterCollection;
		MarkerDataLoader markDataLoader;
		String[] sampleList;
		
		if (markerList == null) {
			markerList = proj.getMarkerNames();
		}

		try {
			markDataLoader = new MarkerDataLoader(proj, markerList, 0);
			clusterFilterCollection = ClusterFilterCollection.load(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, true), proj.getJarStatus());
			sampleList = proj.getSamples();
			out = new PrintWriter(new FileOutputStream(pedFilenameRoot + ".map"));
			for (int i = 0; i < sampleList.length; i++) {
				markDataLoader.getMarkerData(i).getAbGenotypesAfterFilters(clusterFilterCollection, markerList[i], 0);
			}
			out.close();

			out = new PrintWriter(new FileOutputStream(pedFilenameRoot + ".ped"));
			for (int i = 0; i < markerList.length; i++) {
				markDataLoader.getMarkerData(i).getAbGenotypesAfterFilters(clusterFilterCollection, markerList[i], 0);
			}
			out.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Unfinished
	 * @param proj
	 * @param markerList
	 * @param clusterFilterFileName
	 * @param gcThreshold
	 */
	public static boolean createBinaryFileSetFromGenvisisData(Project proj, String clusterFilterFileName, float gcThreshold, String bedFilenameRoot, boolean isSnpMajor, Logger log) {
		String[] targetMarkers;
		int[] indicesOfTargetSampsInProj;
		int[] indicesOfTargetMarksInProj;
		byte[] chrsOfTargetMarks;
		int[] posOfTargetMarks;
		String[] allSampsInProj;
		char[][] abLookup;
		String[] targetSamps;
		

		bedFilenameRoot = proj.getProjectDir() + bedFilenameRoot;
		if (new File(bedFilenameRoot + ".bed").exists() || new File(bedFilenameRoot + ".bim").exists() || new File(bedFilenameRoot + ".fam").exists()) {
			log.reportError("System abort. Plink binary file set \"" + bedFilenameRoot + "\" .bed/.bim/.fam already exist. Please remove the file or use another file name to try again.");
			return false;
		}

		targetSamps = createFamFile(proj, bedFilenameRoot, log);
		allSampsInProj = proj.getSamples();
		if (targetSamps != null) {
			indicesOfTargetSampsInProj = getSortedIndicesOfTargetSampsInProj(allSampsInProj, targetSamps, log);
			targetSamps = new String[indicesOfTargetSampsInProj.length];
			for (int i = 0; i < indicesOfTargetSampsInProj.length; i++) {
				targetSamps[i] = allSampsInProj[ indicesOfTargetSampsInProj[i] ];
			}
		} else {
			indicesOfTargetSampsInProj = null;
			targetSamps = allSampsInProj;
		}

//		allMarksInProj = proj.getMarkerNames();
		targetMarkers = proj.getTargetMarkers(log);
		if (targetMarkers != null) {
			indicesOfTargetMarksInProj = new int[targetMarkers.length];
			chrsOfTargetMarks = new byte[targetMarkers.length];
			posOfTargetMarks = new int[targetMarkers.length];
			getIndicesOfTargetMarks(proj, targetMarkers, indicesOfTargetMarksInProj, chrsOfTargetMarks, posOfTargetMarks, log);

		} else {
			indicesOfTargetMarksInProj = null;
			targetMarkers = proj.getMarkerNames();
			chrsOfTargetMarks = proj.getMarkerSet().getChrs();
			posOfTargetMarks = proj.getMarkerSet().getPositions();
		}

		if (isSnpMajor) {
			abLookup = createBedFileSnpMajor10KperCycle(proj, targetMarkers, indicesOfTargetMarksInProj, targetSamps, indicesOfTargetSampsInProj, bedFilenameRoot, gcThreshold, bedFilenameRoot, log);
//			abLookup = createBedFileSnpMajorAllInMemory(proj, targetMarkers, indicesOfTargetMarksInProj, targetSamps, indicesOfTargetSampsInProj, bedFilenameRoot, gcThreshold, bedFilenameRoot, log);
		} else {
			abLookup = createBedFileIndividualMajor(proj, targetSamps, targetMarkers, indicesOfTargetMarksInProj, bedFilenameRoot, gcThreshold, bedFilenameRoot, log);
		}
		
		if (abLookup == null) {
			log.reportError("Error - failed to create Plink files due to lack of an AB lookup file");
			new File(bedFilenameRoot+".fam").delete();
			return false;
		}

		createBimFile(targetMarkers, chrsOfTargetMarks, posOfTargetMarks, abLookup, bedFilenameRoot, log);

		return true;
	}

	public static int[] getSortedIndicesOfTargetSampsInProj(String[] allSampInProj, String[] targetSamps, Logger log) {
//		String[] sampleNames;
		int[] indicesOfTargetSampInProj;
		Hashtable<String, Integer> hash;
		Enumeration<String> keys;
		hash = new Hashtable<String, Integer>();
		boolean found;
//		sampleNames = proj.getSamples();

//		indices = new int[sampList.length];
		for (int i = 0; i < targetSamps.length; i++) {
			if (hash.containsKey(targetSamps[i])) {
				log.reportError("Warning - duplicate sample id in the list of samples to include: " + targetSamps[i]);
			} else {
				found = false;
				for (int j = 0; j < allSampInProj.length; j++) {
					if (targetSamps[i].equals(allSampInProj[j])) {
						hash.put(targetSamps[i], j);
						found = true;
						break;
					}
				}
				if (! found) {
					log.reportError("Warning - The following from the target sample list was not found in the list of all samples in the project: " + targetSamps[i]);
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

//	public static int[] getIndicesOfMarksToInclude(Project proj, log) {
	public static void getIndicesOfTargetMarks(Project proj, String[] inputTargetMarkers, int[] outputIndicesOfTargetMarks, byte[] outputChrOfTagetMarks, int[] outputPosOfTargetMarks, Logger log) {
		String[] allMarksInProj;
//		Hashtable<String,String> hash;
		MarkerSet markerSet;
		byte[] chrs;
		int[] positions;
		boolean[] found;
	
//		hash = new Hashtable<String,String>();
		if (inputTargetMarkers == null) {
			log.reportError("inputTargetMarkers is null");

		} else {
			markerSet = proj.getMarkerSet();
			allMarksInProj = markerSet.getMarkerNames();
			chrs = markerSet.getChrs();
			positions = markerSet.getPositions();
			found = new boolean[inputTargetMarkers.length];
			for (int i = 0; i < inputTargetMarkers.length; i++) {
//				if (hash.containsKey(inputTargetMarkers[i])) {
//					System.err.println("Warning - duplicate marker name: " + inputTargetMarkers[i] + " in targetMarks.txt");
//				}
				for (int j = 0; j < allMarksInProj.length; j++) {
					if (inputTargetMarkers[i].equals(allMarksInProj[j])) {
						outputIndicesOfTargetMarks[i] = j;
						found[i] = true;
						break;
					}
				}
//				hash.put(allMarksInProj[i], i+"");
			}
			Arrays.sort(outputIndicesOfTargetMarks);
			for (int i = 0; i < found.length; i++) {
				if (! found[i]) {
					log.reportError("Warning - the following marker from target marker list is not found in whole project's marker list: " + inputTargetMarkers[i]);
					break;
				} else {
					inputTargetMarkers[i] = allMarksInProj[outputIndicesOfTargetMarks[i]];
					outputChrOfTagetMarks[i] = chrs[outputIndicesOfTargetMarks[i]];
					outputPosOfTargetMarks[i] = positions[outputIndicesOfTargetMarks[i]];
				}
			}
		}

//		chrs = markerSet.getChrs();
//		positions = markerSet.getPositions();
////		targetMarkers = proj.getFilename(Project.TARGET_MARKERS_FILENAME, false, false);
////		if (new File(targetMarkers).exists()) {
////			targets = HashVec.loadFileToStringArray(targetMarkers, false, false, new int[] {0}, false);
//		if (inputTargetMarkers != null) {
//			indicesOfMarksSelected = new int[inputTargetMarkers.length];
//			prob = false;
//			for (int i = 0; i < inputTargetMarkers.length; i++) {
//				if (hash.containsKey(inputTargetMarkers[i])) {
//					indicesOfMarksSelected[i] = Integer.parseInt(hash.get(inputTargetMarkers[i]));
//				} else {
//					System.err.println("Error - target marker '" + inputTargetMarkers[i] + "' was not found in the MarkerSet");
//					prob = true;
//				}
//			}
//			if (prob) {
//				System.exit(1);
//			}
//
////			chrOfMarksSelected = new byte[indicesOfMarksSelected.length];
////			posOfMarksSelected = new int[indicesOfMarksSelected.length];
////			for (int i = 0; i < indicesOfMarksSelected.length; i++) {
////				chrOfMarksSelected[i] = chrs[indicesOfMarksSelected[i]];
////				posOfMarksSelected[i] = positions[indicesOfMarksSelected[i]];
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
//			indicesOfMarksSelected = Array.intArray(allMarksInProj.length);
////			chrOfMarksSelected = chrs;
////			posOfMarksSelected = positions;
//		}
//
//		return indicesOfMarksSelected;
	}

	public static char[][] createBedFileIndividualMajor(Project proj, String[] targetSamps, String[] targetMarks, int[] indicesOfTargetMarks, String clusterFilterFileName, float gcThreshold, String bedDirAndFilenameRoot, Logger log) {
//		String[] markList;
		RandomAccessFile out1;
		byte[] outStream;
		byte[] genotypes;
//		byte[] genotypesOfTargetMarks;
		ClusterFilterCollection clusterFilterCollection;
		Sample fsamp;
		char[][] abLookup = null;
//		MarkerSet markerSet = proj.getMarkerSet();
//		String[] markerNames;

//		allMarksInProj = proj.getMarkerNames();
//		markList = new String[indicesOfSelectedMarks.length];
//		for (int i = 0; i < indicesOfSelectedMarks.length; i++) {
//			markList[i] = allMarksInProj[indicesOfSelectedMarks[i]];
//		}

		if (Files.exists(clusterFilterFileName, proj.getJarStatus())) {
			clusterFilterCollection = ClusterFilterCollection.load(clusterFilterFileName, proj.getJarStatus());
		} else {
			clusterFilterCollection = null;
		}

		try {
			clusterFilterCollection = ClusterFilterCollection.load(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, true), proj.getJarStatus()); 
			if (clusterFilterCollection == null) {
				abLookup = null;
			} else {
				abLookup = new ABLookup(targetMarks, proj.getFilename(Project.AB_LOOKUP_FILENAME), true, true).getLookup();
			}
			
			out1 = new RandomAccessFile(bedDirAndFilenameRoot + ".bed", "rw");
			outStream = new byte[3];
			outStream[0] = (byte) 108;	// 0b01101100
			outStream[1] = (byte) 27;	// 0b00011011
			outStream[2] = (byte) 0;	// 0b00000000 <-- Be careful here
			out1.write(outStream);

			for (int i = 0; i < targetSamps.length; i++) {
				fsamp = proj.getFullSampleFromRandomAccessFile(targetSamps[i]);

				if (fsamp == null) {
					System.err.println("Error - the DNA# " + targetSamps[i] + " was listed in the pedigree file but " + targetSamps[i] + ".fsamp was not found in directory: " + proj.getDir(Project.SAMPLE_DIRECTORY));
					genotypes = new byte[1];

				} else {
					if (clusterFilterCollection == null) {
						genotypes = fsamp.getAB_Genotypes(indicesOfTargetMarks);
					} else {
						genotypes = fsamp.getAB_GenotypesAfterFilters(targetMarks, indicesOfTargetMarks, clusterFilterCollection, gcThreshold);
//						genotypes = proj.getMarkerSet().translateABtoForwardGenotypes(genotypes, abLookup);
					}
					out1.write(encodePlinkBedBytesForASingleMarkOrSamp(genotypes));

//					genotypesOfTargetMarks = new byte[indicesOfTargetMarks.length];
//					for (int j = 0; j < indicesOfTargetMarks.length; j++) {
//						genotypesOfTargetMarks[j] = genotypes[indicesOfTargetMarks[i]];
//					}
//					out1.write(encodePlinkBedBytesForASingleMarkOrSamp(genotypesOfTargetMarks));

//					if (clusterFilterCollection == null) {
//						genotypes = fsamp.getForwardGenotypes(gcThreshold);
//					} else {
//						genotypes = markerSet.translateABtoForwardGenotypes(fsamp.getAB_GenotypesAfterFilters(targetMarks, clusterFilterCollection, gcThreshold), abLookup);
//					}

//					alleles = new String[genotypesSelected.length][2];
//					for (int j = 0; j<indicesOfTargetMarks.length; j++) {
//						genotypesSelected[j] = genotypes[indicesOfTargetMarks[i]];
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
			}

			out1.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return abLookup;
	}

	public static char[][] createBedFileSnpMajorAllInMemory(Project proj, String[] targetMarks, int[] indicesOfTargetMarksInProj, String[] targetSamps, int[] indicesOfTargetSampsInProj, String clusterFilterFileName, float gcThreshold, String bedDirAndFilenameRoot, Logger log) {
		RandomAccessFile out1;
		byte[] outStream;
		byte[] genotypes;
		byte[] genotypesOfTargetSamps;
		ClusterFilterCollection clusterFilterCollection;
		MarkerDataLoader markDataLoader;
		char[][] abLookup = null;
		Hashtable<String, Vector<String>> batches;
		MarkerData[] markerData;
		String[] filenames;
		Vector<String> v;
		int[] indicesOfMarksInFileForCurrentFile;
		String[] marksOfThisFile;
		String[] temp;
		Hashtable<String, Float> outliersHash;
		long sampleFingerPrint;
		String[] allSampsInProj;
		int[] indicesOfMarksInProjForCurrentFile;
		String dir;
		long startTime;
		
		startTime = new Date().getTime();

		try {
			if (Files.exists(clusterFilterFileName, proj.getJarStatus())) {
				clusterFilterCollection = ClusterFilterCollection.load(clusterFilterFileName, proj.getJarStatus());
			} else {
				clusterFilterCollection = null;
			}
			if (clusterFilterCollection == null) {
				abLookup = null;
			} else {
				abLookup = new ABLookup(targetMarks, proj.getFilename(Project.AB_LOOKUP_FILENAME), true, true).getLookup();
			}
			
			sampleFingerPrint = proj.getSampleList().getFingerprint();
			allSampsInProj = proj.getSamples();
			dir = proj.getDir(Project.MARKER_DATA_DIRECTORY);
			markDataLoader = new MarkerDataLoader(proj, targetMarks, 0);
			outliersHash = MarkerDataLoader.loadOutliers(proj);
			batches = markDataLoader.getBatches();
			filenames = HashVec.getKeys(batches);
			genotypesOfTargetSamps = new byte[indicesOfTargetSampsInProj.length];

			out1 = new RandomAccessFile(bedDirAndFilenameRoot + ".bed", "rw");
			outStream = new byte[3];
			outStream[0] = (byte) 108;	// 0b01101100
			outStream[1] = (byte) 27;	// 0b00011011
			outStream[2] = (byte) 1;	// 0b00000001 <-- be careful here
			out1.write(outStream);

//			for (int i = 0; i < targetMarks.length; i++) {
			for (int i = 0; i < filenames.length; i++) {
				v = batches.get(filenames[i]);
				marksOfThisFile = new String[v.size()];
				indicesOfMarksInFileForCurrentFile = new int[marksOfThisFile.length];
				indicesOfMarksInProjForCurrentFile = new int[marksOfThisFile.length];
				for (int j = 0; j < v.size(); j++) {
					temp = v.elementAt(j).split("\t");
					marksOfThisFile[j] = temp[0];
					indicesOfMarksInFileForCurrentFile[j] = Integer.parseInt(temp[1]);
					if(indicesOfTargetMarksInProj != null) {
						for (int k = 0; k < indicesOfTargetMarksInProj.length; k++) {
							if (marksOfThisFile[j].equals(targetMarks[k])) {
								indicesOfMarksInProjForCurrentFile[j] = indicesOfTargetMarksInProj[k];
								break;
							}
						}
					} else {
						for (int k = 0; k < targetMarks.length; k++) {
							if (marksOfThisFile[j].equals(targetMarks[k])) {
								indicesOfMarksInProjForCurrentFile[j] = k;
								break;
							}
						}
					}
				}
				markerData = MarkerDataLoader.loadFromRAF(null, null, null, allSampsInProj, dir + filenames[i], indicesOfMarksInProjForCurrentFile, indicesOfMarksInFileForCurrentFile, false, true, false, false, true, sampleFingerPrint, outliersHash);

				for (int j = 0; j < markerData.length; j++) {
					genotypes = markerData[j].getAbGenotypesAfterFilters(clusterFilterCollection, marksOfThisFile[j], 0);
					for (int k = 0; k < indicesOfTargetSampsInProj.length; k++) {
						genotypesOfTargetSamps[k] = genotypes[indicesOfTargetSampsInProj[k]];
					}
					out1.write(encodePlinkBedBytesForASingleMarkOrSamp(genotypesOfTargetSamps));
				}

				//TODO Something does not make sense here. alleles[j] seems to be wrong here.
//				genotypesSelected = new byte[indicesOfSelectedSamps.length];
//				for (int j = 0; j < indicesOfSelectedSamps.length; j++) {
//					genotypesSelected[j] = genotypes[indicesOfSelectedSamps[i]];

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
		
		log.report("Finished creating binary PLINK files in "+ext.getTimeElapsed(startTime));

		return abLookup;
	}

	public static char[][] createBedFileSnpMajor10KperCycle(Project proj, String[] targetMarkers, int[] indicesOfTargetMarkersInProj, String[] targetSamps, int[] indicesOfTargetSampsInProj, String clusterFilterFileName, float gcThreshold, String bedDirAndFilenameRoot, Logger log) {
		RandomAccessFile out;
		byte[] outStream;
		byte[] genotypes;
		byte[] genotypesOfTargetSamps;
		ClusterFilterCollection clusterFilterCollection;
		MarkerDataLoader markDataLoader;
		char[][] abLookup;
		Hashtable<String, Vector<String>> batches;
		MarkerData[] markerData;
		String[] filenames;
		Vector<String> v;
		int[] indicesOfMarksInFileForCurrentFile;
		String[] marksOfThisFile;
		String[] temp;
		Hashtable<String, Float> outliersHash;
		long sampleFingerPrint;
		String[] allSampsInProj;
		int[] indicesOfMarksInProjForCurrentFile;
		String dir;
		long startTime, subTime;
		Hashtable<String,Integer> hash;
		int targetIndex;

		startTime = new Date().getTime();

		if (Files.exists(clusterFilterFileName, proj.getJarStatus())) {
			clusterFilterCollection = ClusterFilterCollection.load(clusterFilterFileName, proj.getJarStatus());
		} else {
			clusterFilterCollection = null;
		}
		if (Files.exists(proj.getFilename(Project.AB_LOOKUP_FILENAME))) {
			abLookup = new ABLookup(targetMarkers, proj.getFilename(Project.AB_LOOKUP_FILENAME), true, true).getLookup();
		} else if (clusterFilterCollection == null) {
			abLookup = new char[targetMarkers.length][];
		} else {
			JOptionPane.showMessageDialog(null, "Error - could not find AB lookup file '"+proj.getFilename(Project.AB_LOOKUP_FILENAME)+"'; this file needs to be created, as it is not otherwise possible to export to Plink when there are cluster filters, ", "Error", JOptionPane.ERROR_MESSAGE);
			return null;
		}
		
		hash = new Hashtable<String,Integer>();
		for (int i = 0; i<targetMarkers.length; i++) {
			if (hash.containsKey(targetMarkers[i])) {
				System.err.println("Warning - duplicate marker name: "+targetMarkers[i]);
			}
			hash.put(targetMarkers[i], i);
		}

		sampleFingerPrint = proj.getSampleList().getFingerprint();
		allSampsInProj = proj.getSamples();
		dir = proj.getDir(Project.MARKER_DATA_DIRECTORY);
		markDataLoader = new MarkerDataLoader(proj, targetMarkers, 0);
		outliersHash = MarkerDataLoader.loadOutliers(proj);
		batches = markDataLoader.getBatches();
		filenames = HashVec.getKeys(batches);
		genotypesOfTargetSamps = new byte[indicesOfTargetSampsInProj.length];

		try {
			out = new RandomAccessFile(bedDirAndFilenameRoot + ".bed", "rw");
			outStream = new byte[3];
			outStream[0] = (byte) 108;	// 0b01101100
			outStream[1] = (byte) 27;	// 0b00011011
			outStream[2] = (byte) 1;	// 0b00000001 <-- be careful here
			out.write(outStream);

			subTime = new Date().getTime();
//			for (int i = 0; i < targetMarks.length; i++) {
			for (int i = 0; i < filenames.length; i++) {
				
				v = batches.get(filenames[i]);
				marksOfThisFile = new String[v.size()];
				indicesOfMarksInFileForCurrentFile = new int[marksOfThisFile.length];
				indicesOfMarksInProjForCurrentFile = new int[marksOfThisFile.length];
				System.out.println("Prepping batch "+(i+1)+" of "+filenames.length+" which has "+v.size()+" elements");
				for (int j = 0; j < v.size(); j++) {
					temp = v.elementAt(j).split("\t");
					marksOfThisFile[j] = temp[0];
					indicesOfMarksInFileForCurrentFile[j] = Integer.parseInt(temp[1]);
					indicesOfMarksInProjForCurrentFile[j] = hash.get(temp[0]);
				}
				System.out.println("done prepping in "+ext.getTimeElapsed(subTime));
				
				subTime = new Date().getTime();
				System.out.println("Starting batch "+(i+1)+" of "+filenames.length);
				markerData = MarkerDataLoader.loadFromRAF(null, null, null, allSampsInProj, dir + filenames[i], indicesOfMarksInProjForCurrentFile, indicesOfMarksInFileForCurrentFile, false, true, false, false, true, sampleFingerPrint, outliersHash);
				System.out.println("Done loading in "+ext.getTimeElapsed(subTime));

				subTime = new Date().getTime();
				for (int j = 0; j < markerData.length; j++) {
					targetIndex = hash.get(marksOfThisFile[j]);
					genotypes = markerData[j].getAbGenotypesAfterFilters(clusterFilterCollection, marksOfThisFile[j], 0);
					if (ext.indexOfStr(marksOfThisFile[j], new String[] {"rs6650104", "rs3115860", "rs17160939"}) >= 0) {
						markerData[j].dump(marksOfThisFile[j]+".xln", allSampsInProj);
					}
					for (int k = 0; k < indicesOfTargetSampsInProj.length; k++) {
						genotypesOfTargetSamps[k] = genotypes[indicesOfTargetSampsInProj[k]];
					}
					if (abLookup[targetIndex] == null) {
						abLookup[targetIndex] = markerData[j].getAB_AlleleMappings();
						if (ext.indexOfStr(marksOfThisFile[j], new String[] {"rs6650104", "rs3115860", "rs17160939"}) >= 0) {
							abLookup[targetIndex] = markerData[j].getAB_AlleleMappings();
							System.out.println(marksOfThisFile[j]+"\t"+abLookup[targetIndex][0]+"\t"+abLookup[targetIndex][1]);
						}
					}
					
					out.write(encodePlinkBedBytesForASingleMarkOrSamp(genotypesOfTargetSamps));
				}
				System.out.println("Done writing in "+ext.getTimeElapsed(subTime));
			}
			out.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		log.report("Finished creating binary PLINK files in "+ext.getTimeElapsed(startTime));

		return abLookup;
	}

	public static boolean createBimFile(String[] targetMarks, byte[] chrsOfTargetMarks, int[] posOfTargetMarks, char[][] abLookup, String bimDirAndFilenameRoot, Logger log) {
		PrintWriter writer;

		if (abLookup == null) {
			System.err.println("Error - abLookup cannot be null; failed to create .bim file");
			return false;
		}
		
		try {
			writer = new PrintWriter(new FileWriter(bimDirAndFilenameRoot + ".bim"));
			for (int i = 0; i < targetMarks.length; i++) {
				writer.println(chrsOfTargetMarks[i] + "\t" + targetMarks[i] + "\t0\t" + posOfTargetMarks[i] + "\t" + abLookup[i][0] + "\t" + abLookup[i][1]); //TODO alleles[][] matching chrs[]
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

	public static String[] createFamFile(Project proj, String famDirAndFilenameRoot, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		int count;
		String temp;
		String[] line;
		Vector<String> dna;
		String[] allSamples;
		String filename;
		
		allSamples = proj.getSamples();
		dna = new Vector<String>();

		try {
			filename = proj.getFilename(Project.PEDIGREE_FILENAME);
			if (!new File(filename).exists()) {
				log.reportError("Error - pedigree file ('"+filename+"') is not found. Program aborted.");
				return null;
			}
			reader = new BufferedReader(new FileReader(proj.getFilename(Project.PEDIGREE_FILENAME)));
			writer = new PrintWriter(new FileWriter(famDirAndFilenameRoot+".fam"));
			count = 1;
			while (reader.ready()) {
				count++;
				if (count % 100 == 0) {
					System.out.println(count);
				}
				temp = reader.readLine().trim();
				line = temp.split(ext.determineDelimiter(temp));
				writer.println(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]);
				line = temp.split("[\\s]+");
				if (temp.equals("")) {
					// then do nothing
				} else if (line.length < 7) {
					log.reportError("Error - starting at line "+(count-1)+(line.length<3?"":" (individual "+line[0]+"-"+line[1]+")")+" there are only "+line.length+" columns in pedigree file '"+proj.getFilename(Project.PEDIGREE_FILENAME)+"'.");
					log.reportError("  Pedigree files require 7 columns with no header: FID IID FA MO SEX PHENO DNA");
					log.reportError("  where DNA is the sample name associated with the genotypic data (see the "+proj.getDir(Project.SAMPLE_DIRECTORY)+" directory for examples)");
					reader.close();
					writer.close();
					return null;
				} else if (ext.isMissingValue(line[6])) {
					dna.add(null);
				} else if (ext.indexOfStr(line[6], allSamples) == -1) {
					log.reportError("Error - sample '"+line[6]+"' was not found in the list of the projects samples; all marker data will be set to missing for individual "+line[0]+"-"+line[1]);
					dna.add(null);
				} else {
					dna.add(line[6]);
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + famDirAndFilenameRoot+".fam" + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + famDirAndFilenameRoot+".fam" + "\"");
			System.exit(2);
		}

		return Array.toStringArray(dna);
	}

	/**
	 * Note: This method does not address the memory management issue. The caller of this method is supposed to implement the optimization of memory usage.
	 * 
	 * Transpose individual-major data (Sample lead data) to SNP-major one (Marker lead data), or vice versa. If you want to apply this to the latter, just
	 * treat all the variables labeled with Mark as if they were labeled with Samp, and also Samp with Mark.
	 *  
	 * The index parameters are all inclusive, meaning iStartSamp and iEndSamp are all going to be included in the output.
	 * 
	 */
	public static void transposeBedBytes(byte[] inputSampBytes, int iStartSamp, int iEndSamp, int iStartMark, int iEndMark, byte[] outputMarkBytes, int nTotalMarksInProj, int nTotalSampsInProj) {
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

		nBytesPerMark = (int) Math.ceil((double) nTotalSampsInProj / 4);
		nBytesPerSamp = (int) Math.ceil((double) nTotalMarksInProj / 4);
//		iByteMark = iStartMark * bytesPerMark + iStartSamp / 4;
		offsetMark = iStartMark * nBytesPerMark;
		offsetSamp = iStartSamp * nBytesPerSamp;

		for (int i = iStartSamp; i < iEndSamp; i++) {
			iByteMark = offsetMark + i / 4;
			iBitMark = (byte) (i % 4);
			
			for (int j = iStartMark; j < iEndMark; j++) {
				iByteSamp = offsetSamp + j / 4;
				iBitSamp = (byte) (j % 4);

				outputMarkBytes[iByteMark] = (byte) ((outputMarkBytes[iByteMark] & (0xff - (0x03 << (iBitMark * 2)))) | (inputSampBytes[iByteSamp] & (0x03 << (iBitSamp * 2))));
			}
		}

//		for (int i = iStartByteSamp; i < iEndByteSamp - 1; i++) {
//			iBitMark = (byte) (iStartSamp % 4);
//			for (int j = iStartBitSamp; j < 4; j++) {
//				outputMarkBytes[iByteMark] = (byte) ((outputMarkBytes[iByteMark] & (0xff - (0x03 << (iBitMark * 2)))) | (inputASingleSampsBytes[i] & (0x03 << (j * 2))));
//				iByteMark += bytesPerMark;
//			}
//			iStartBitSamp = 0;
//			iBitMark = 0;
//		}

//		for (int j = iStartBitSamp; j <= 4; j++) {
//			iByteMark = i * bytesPerMark + j / 4;
//			outputMarkBytes[iByteMark] = (byte) ((outputMarkBytes[iByteMark] & 0x03) | (inputASingleSampsBytes[i] & 0x3));
//		}
	}

	public static void transposeOneBedByte(byte inputSampByte, int iCurrentSampInProj, int iEndSamp, int iStartMark, int iEndMark, byte[] outputMarkBytes, int nBytesPerSamp, int nBytesPerMark) {
		int offsetMark;
		int iByteMark;
		byte iBitMark;
		byte iBitSamp;

		offsetMark = iStartMark * nBytesPerMark;
		iByteMark = offsetMark + iCurrentSampInProj / 4;
		iBitMark = (byte) (iCurrentSampInProj % 4);

		for (int j = iStartMark; j < iEndMark; j++) {
			iBitSamp = (byte) (j % 4);

			outputMarkBytes[iByteMark] = (byte) ((outputMarkBytes[iByteMark] & (0xff - (0x03 << (iBitMark * 2)))) | (inputSampByte & (0x03 << (iBitSamp * 2))));
		}
	}


	public static byte[] encodePlinkBedBytesForASingleMarkOrSamp (byte[] genotype) {
		int iBytes;
		byte[] result;
		int nBytes;
		byte shift;

		nBytes = (int) Math.ceil((double) genotype.length / 4);
		iBytes = -1;
		result = new byte[nBytes];

//		for (int i = 0; i < result.length; i++) {
//			result[i] = (byte) 0xAA;	//initilize the array to be 0b10101010, the null genotype defined by plink bed data.
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

	public static byte[] decodeBedBytesOfASingleMarkOrSamp (byte[] bedBytes, int startIndex, int[] indicesOfSampsOrMarks) {
		byte[] genotypes;
		int indexBedBytes;
		int indexBedByte;

		genotypes = new byte[indicesOfSampsOrMarks.length];
		try {
			for (int i = 0; i <= genotypes.length; i++) {
				indexBedBytes = indicesOfSampsOrMarks[i] / 4;
				indexBedByte = indicesOfSampsOrMarks[i] % 4;
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

	public static void displayBits(byte data) {
		for (int j = 0; j < 8; j++) {
			System.out.print((j==4? " " : "") + ((data & 0x80) >> 7));
			data = (byte) (data << 1);
		}
		System.out.println("");
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
			//Test the plink file header && Marker Dominant or Sample Dominant?

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


	public static void displayBitMapOfBedFile(String bedFileName, int nBytes) {
		RandomAccessFile in;
		byte[] readBuffer;
		
		try {
			in = new RandomAccessFile(bedFileName, "r");
			readBuffer = new byte[(int) Math.min(nBytes, in.length())];
			in.read(readBuffer);
			for (int i = 0; i < readBuffer.length; i++) {
				displayBits(readBuffer[i]);
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void main (String[] args) {
		new File("D:/data/practice/plinkers.bed").delete();
		new File("D:/data/practice/plinkers.fam").delete();
		new File("D:/data/practice/plinkers.bim").delete();
		createBinaryFileSetFromGenvisisData(new Project("D:/home/npankrat/projects/practice.properties", false), null, 0.15f, "plinkers", true, new Logger());
		
//		displayBitMapOfBedFile("D:/GEDI_exome/gwas_plink.bed", 30);
//
//		String fileDir = "N:/statgen/Genvisis/plinkTesting/";
//		String fileNameRoot = "plink";
//
//		convertPedToBed(fileDir, fileNameRoot, true);
	}
}
