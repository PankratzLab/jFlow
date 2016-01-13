package cnv.manage;

import java.io.*;
import java.util.*;

import javax.swing.JFileChooser;
import filesys.Segment;
import cnv.var.CNVariant;
import common.*;

public class ExportCNVsToPedFormat {
	private static final String MATRIX_FORMAT = "MATRIX";
	private static final String PLINK_TEXT_FORMAT = "PLINK_TEXT";
	private static final String PLINK_BINARY_FORMAT = "PLINK_BINARY_FORMAT";
	private static final String RFGLS_FORMAT = "RFGLS";

	/*
	 * Convert a cnv data file into
	 */
	// Step 1a
	public static void convertFromRawCNV(String cnvFilename, String pedFilename, String outputRoot, String endOfLine, boolean rfglsOutput, boolean includeDele, boolean includeDupl, boolean ordered, boolean collapsed, boolean homozygous, boolean excludeMonomorphicLoci, int markersPerFile, int win) {
		
	}

	// Step 1b
	public static void convertFromTrioOrJointCNVs(String cnvFilename, String pedFilename, String outputRoot, String endOfLine, boolean rfglsOutput, boolean includeDele, boolean includeDupl, boolean ordered, boolean collapsed, boolean homozygous, boolean excludeMonomorphicLoci, int markersPerFile, int win) {
		
	}

	// Step 2
	public static void export(CNVariant[] cnvs, String[][] pedstructure, String outputRoot, String endOfLine, boolean rfglsOutput, boolean includeDele, boolean includeDupl, boolean ordered, boolean collapsed, boolean homozygous, boolean excludeMonomorphicLoci, int markersPerFile, int win) {
		
	}
	
	public static void export(String cnvFilename, String pedFilename, String outputRoot, String endOfLine, String fileFormat, boolean includeDele, boolean includeDupl, boolean ordered, boolean collapsed, boolean homozygousOnly, boolean excludeMonomorphicLoci, int markersPerFile, int windowInBasepairs, Logger log) {
		PrintWriter writer;
		CNVariant[] cnvs;
		Hashtable<String,String> allChrPosHash;
		Hashtable<String,String> sampleListHashFromCnvOrPedData;
		String[][] dnaMapping;
		Hashtable<String,Vector<String>> mzTwins;
		Vector<String> mzPairs;
		String[] allChrPosKeys;
		byte[][] currentCNs;
		byte[] previousCNs;
		Segment[] allChrPosSegs, currentChrPosSegs;
		String[] line;
		byte currentCN;
		long time;
		String[] tempSampleList;
		String[] finalSampleList;
		int cnvIndex;
		boolean done;
		String outputFilename;
		int numMarkers;
		int currentSampleIndex;
		int indexOfCurrentSeg;
		int countValidLoci;
		int fileNumber;
//		Logger log;
		Vector<CNVariant> cnVector;
		CNVariant cnv;
		int numBaseCNVs;

		log.report("Generating files for "+outputRoot);
		if (!Files.exists(ext.parseDirectoryOfFile(outputRoot), false)) {
			log.report("Created directory: '"+ext.parseDirectoryOfFile(outputRoot)+"'");
			new File(ext.parseDirectoryOfFile(outputRoot)).mkdirs();
		}

		time = new Date().getTime();
		mzTwins = new Hashtable<String, Vector<String>>();
		if (pedFilename != null) {
			line = Files.getHeaderOfFile(pedFilename, new Logger());
			ext.checkHeader(line, new String[] {"FAMID", "ID", "DNA"}, true);
			sampleListHashFromCnvOrPedData = HashVec.loadFileToHashString(pedFilename, new int[] {0,1}, new int[] {-7}, pedFilename.endsWith(".csv"), "\t", true, false, false);
			dnaMapping = HashVec.loadFileToStringMatrix(pedFilename, true, new int[] {0,1,2}, "[\\s]+", false, 10000, false);
			for (int i = 0; i < dnaMapping.length; i++) {
				if (!dnaMapping[i][1].equals(dnaMapping[i][2])) {
					HashVec.addToHashVec(mzTwins, dnaMapping[i][0]+"\t"+dnaMapping[i][2], dnaMapping[i][0]+"\t"+dnaMapping[i][1], false);
				}
			}
		} else {
			sampleListHashFromCnvOrPedData = null;
		}
		
		log.report("Matched "+mzTwins.size()+" twins");
		
		cnVector = CNVariant.loadPlinkFile(cnvFilename, sampleListHashFromCnvOrPedData, false);
		if (mzTwins.size() > 0) {
			numBaseCNVs = cnVector.size();
			for (int i = 0; i < numBaseCNVs; i++) {
				cnv = cnVector.elementAt(i);
				if (mzTwins.containsKey(cnv.getFamilyID()+"\t"+cnv.getIndividualID())) {
					mzPairs = mzTwins.get(cnv.getFamilyID()+"\t"+cnv.getIndividualID());
					for (int j = 0; j < mzPairs.size(); j++) {
						cnv = cnVector.elementAt(i).clone();
						line = mzPairs.elementAt(j).split("[\\s]+");
						cnv.setFamilyID(line[0]);
						cnv.setIndividualID(line[1]);
						cnVector.add(cnv);
					}					
				}
			}
			
		}
		cnvs = CNVariant.toCNVariantArray(cnVector);
		log.report("Loaded file in " + ext.getTimeElapsed(time));
		
		//Generate chrPosition and sampleHash, to be used for the rows and columns of the final result.
		time = new Date().getTime();
		allChrPosHash = new Hashtable<String,String>();
		if (pedFilename==null) {
			sampleListHashFromCnvOrPedData = new Hashtable<String,String>();
		}
		for (int i=0; i<cnvs.length; i++) {
			allChrPosHash.put(cnvs[i].getChr()+"\t"+cnvs[i].getStart(), "");
			allChrPosHash.put(cnvs[i].getChr()+"\t"+cnvs[i].getStop(), "");
			allChrPosHash.put(cnvs[i].getChr()+"\t"+(cnvs[i].getStop()+1), "");
			if (pedFilename==null && !sampleListHashFromCnvOrPedData.containsKey(cnvs[i].getFamilyID()+"\t"+cnvs[i].getIndividualID())) {
				sampleListHashFromCnvOrPedData.put(cnvs[i].getFamilyID()+"\t"+cnvs[i].getIndividualID(), sampleListHashFromCnvOrPedData.size()+"");
			}
		}
		allChrPosKeys = HashVec.getKeys(allChrPosHash);
		allChrPosSegs = new Segment[allChrPosKeys.length];
		for (int i = 0; i < allChrPosSegs.length; i++) {
			line = allChrPosKeys[i].split("\t");
			allChrPosSegs[i] = new Segment(Byte.parseByte(line[0]), Integer.parseInt(line[1])-windowInBasepairs, Integer.parseInt(line[1])+windowInBasepairs);
		}
		log.report("Generated hashtable of positions in " + ext.getTimeElapsed(time));

		time = new Date().getTime();
		allChrPosSegs = Segment.sortSegments(allChrPosSegs);
		log.report("Sorted positions in " + ext.getTimeElapsed(time));

		time = new Date().getTime();
		cnvs = CNVariant.sort(cnvs);
		log.report("Sorted CNVariants in " + ext.getTimeElapsed(time));

		
		cnvIndex = 0;
		fileNumber = 0;
		countValidLoci = 0;
		writer = null;
		outputFilename = "very first file";
		for (int startPosition = 0; startPosition < allChrPosSegs.length; startPosition+=markersPerFile) {
			numMarkers = Math.min(markersPerFile, allChrPosSegs.length-startPosition);
			log.report("Calculating CN for positions "+(startPosition+1)+" through "+(startPosition+numMarkers));

			currentChrPosSegs = new Segment[numMarkers];
			for (int i = 0; i < currentChrPosSegs.length; i++) {
				currentChrPosSegs[i] = allChrPosSegs[startPosition+i];
			}
			
			currentCNs = new byte[currentChrPosSegs.length][sampleListHashFromCnvOrPedData.size()];
			previousCNs = Array.byteArray(sampleListHashFromCnvOrPedData.size(), Byte.MIN_VALUE);

			while (cnvs[cnvIndex].getChr() < currentChrPosSegs[0].getChr() || (cnvs[cnvIndex].getChr() == currentChrPosSegs[0].getChr() && cnvs[cnvIndex].getStop() < currentChrPosSegs[0].getStart())) {
				cnvIndex++;
			}

			done = false;				
			for (int i = cnvIndex; !done && i < cnvs.length; i++) {
				// determine column from sample hash
				currentSampleIndex = sampleListHashFromCnvOrPedData.containsKey(cnvs[i].getFamilyID()+"\t"+cnvs[i].getIndividualID())?Integer.parseInt(sampleListHashFromCnvOrPedData.get(cnvs[i].getFamilyID()+"\t"+cnvs[i].getIndividualID())):-1;
				if (cnvs[i].getChr() > currentChrPosSegs[currentChrPosSegs.length-1].getChr() || (cnvs[i].getChr() == currentChrPosSegs[currentChrPosSegs.length-1].getChr() && cnvs[i].getStart() > currentChrPosSegs[currentChrPosSegs.length-1].getStop())) {
					done = true;
				} else if (currentSampleIndex >= 0) {
					indexOfCurrentSeg = Math.max(Segment.binarySearchForStartPositions(cnvs[i], currentChrPosSegs), 0);
					while (indexOfCurrentSeg < currentChrPosSegs.length && currentChrPosSegs[indexOfCurrentSeg].overlaps(cnvs[i])) {
						currentCNs[indexOfCurrentSeg][currentSampleIndex] = (byte) (cnvs[i].getCN()-2);
						indexOfCurrentSeg++;
					}

				}
			}

			tempSampleList = HashVec.getKeys(sampleListHashFromCnvOrPedData, false, false);
			finalSampleList = new String[sampleListHashFromCnvOrPedData.size()];
			for (int i = 0; i < tempSampleList.length; i++) {
				finalSampleList[Integer.parseInt(sampleListHashFromCnvOrPedData.get(tempSampleList[i]))] = ext.replaceAllWith(tempSampleList[i], "\t", "-");
			}

			for (int i = 0; i < currentCNs.length; i++) {
				for (int j = 0; j < currentCNs[i].length; j++) {
					currentCN = currentCNs[i][j];
					if (!includeDele && currentCN < 0) {
						currentCN = 0;
					}
					if (!includeDupl && currentCN > 0) {
						currentCN = 0;
					}
					if (!ordered) {
						currentCN = (byte)Math.abs(currentCN);
					}
					if (homozygousOnly && Math.abs(currentCN) == 1) {
						currentCN = 0;
					}
					if ((collapsed || homozygousOnly) && Math.abs(currentCN) == 2) {
						currentCN /= 2;
					}
					currentCNs[i][j] = currentCN;
				}
			}

			time = new Date().getTime();
			if (fileFormat.equals(PLINK_TEXT_FORMAT) || fileFormat.equals(PLINK_BINARY_FORMAT)) {
				// dig into whether we need to convert this into MarkerData + an AB_Lookup object, or if this can be modularized
			} else {
				try {
					for (int i = 0; i<currentChrPosSegs.length; i++) {
						if (!excludeMonomorphicLoci || Array.min(currentCNs[i]) < 0 || Array.max(currentCNs[i]) > 0 && !Array.equals(currentCNs[i], previousCNs)) {
							if (countValidLoci % markersPerFile == 0) {
								if (writer != null) {
									writer.close();
									if (fileFormat.equals(RFGLS_FORMAT)) {
										convertToRfglsFormat(outputRoot, fileNumber, endOfLine, log);
									}
									fileNumber++;
								}
								outputFilename = outputRoot+"_"+fileNumber;
								writer = new PrintWriter(new FileWriter(outputFilename));
								writer.print("markerName\t"+Array.toStr(finalSampleList));
								writer.print(endOfLine);
							}
							
							writer.print(currentChrPosSegs[i].getChr()+":"+currentChrPosSegs[i].getStart());
							for (int j=0; j<finalSampleList.length; j++) {
								writer.print("\t"+currentCNs[i][j]);
							}
							writer.print(endOfLine);
							countValidLoci++;
							previousCNs = currentCNs[i];
						}					
					}
				} catch (Exception e) {
					log.reportError("Error writing to '" + outputFilename + "'");
					log.reportException(e);
				}
				if (writer != null) {
					writer.close();
					if (fileFormat.equals(RFGLS_FORMAT)) {
						convertToRfglsFormat(outputRoot, fileNumber, endOfLine, log);
					}
				}
			}
			log.report("    ...finished in " + ext.getTimeElapsed(time));
		}
	}
	
	private static void convertToRfglsFormat(String root, int fileNumber, String endOfLine, Logger log) {
		BufferedReader reader;
		PrintWriter writer, mapWriter;
		String[] line, ids;
		String dir;
		long time;
		
		dir = ext.parseDirectoryOfFile(root);
		
		try {
			reader = new BufferedReader(new FileReader(root+"_"+fileNumber));
			ids = reader.readLine().trim().split("[\\s]+");
			try {
				new File(dir+"rfgls/").mkdirs();
				writer = new PrintWriter(new FileWriter(dir+"rfgls/ids"));
				for (int i = 1; i < ids.length; i++) {
					line = ids[i].split("-");
					if (line.length != 2) {
						log.reportError("Error recreating pedfile since there are hyphens in the FID/IID: "+ids[i]);
						reader.close();
						writer.close();
						return;
					}
					writer.print(line[0]+"\t"+line[1]+endOfLine);
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + dir+"rfgls/ids");
				e.printStackTrace();
			}
				
			writer = new PrintWriter(new FileWriter(dir+"rfgls/"+ext.removeDirectoryInfo(root+"_byMarker_"+fileNumber)));
			mapWriter = new PrintWriter(new FileWriter(dir+"rfgls/"+ext.removeDirectoryInfo(root+"_map_"+fileNumber)));

			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				mapWriter.print(line[0]+endOfLine);
				for (int i = 1; i < line.length; i++) {
					writer.print((i==1?"":"\t")+line[i]);
				}
				writer.print(endOfLine);
			}			
			writer.close();
			mapWriter.close();
			reader.close();
			
			time = new Date().getTime();
			log.report("Transposing "+ext.removeDirectoryInfo(root+"_byMarker_"+fileNumber));
			Files.transpose(dir+"rfgls/"+ext.removeDirectoryInfo(root+"_byMarker_"+fileNumber), "\t", dir+"rfgls/"+ext.removeDirectoryInfo(root+"_bySample_"+fileNumber), "\t", log);
			log.report("\t...finished in "+ext.getTimeElapsed(time));
			new File(dir+"rfgls/"+ext.removeDirectoryInfo(root+"_byMarker_"+fileNumber)).delete();
			
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + root+"_"+fileNumber + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + root+"_"+fileNumber + "\"");
			System.exit(2);
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String cnvFilename = "penncnv.cnv";
		String pedFilename = "pedigree.dat";
		String outputFilename = "cnv_matrix";
		JFileChooser fileChooser = new JFileChooser(); // TODO Should this be open up in any particular directory?
		File input;
		String output;
		boolean includeDele = true;
		boolean includeDupl = true;
		boolean ordered = true;
		boolean collapsed = false; 
		boolean homozygous = false;
		boolean excludeMonomorphicLoci = false;
		int lociPerFile = 5000;
		int window = 0;
		String endOfLine = "\r\n";
//		boolean rfglsOutput = false;
		String fileFormat = MATRIX_FORMAT;

		String usage = "\n" +
				"cnv.analysis.CnvBySample requires the following arguments\n" +
				"   (1) cnv filename (i.e. cnv=" + cnvFilename + " (default))\n" +
				"   (2) pedigree filename (i.e. ped=" + pedFilename + " (default))\n" +
				"   (3) output format (i.e. format=" + fileFormat + " (default; options are "+MATRIX_FORMAT+", "+PLINK_TEXT_FORMAT+", "+PLINK_BINARY_FORMAT+", and "+RFGLS_FORMAT+"))\n" +
				"   (4) output filename (i.e. out=" + outputFilename + " (default))\n" +
				"   (5) to include Deletion or not (i.e. del=" + includeDele + " (default))\n" +
				"   (6) to include Duplication or not (i.e. dup=" + includeDupl + " (default))\n" +
				"   (7) to use ordered value or not (i.e. ord=" + ordered + " (default))\n" +
				"   (8) to use collapsed value or not (i.e. coll=" + collapsed + " (default))\n" +
				"   (9) to only use homozygous variants (i.e. homozygousOnly=" + homozygous + " (default))\n" +
				"   (10) exclude monomorphic loci (i.e. excludeMonomorphicLoci=" + excludeMonomorphicLoci + " (default))\n" +
				"   (11) number of loci per file in the output file (i.e. lociPerFile=" + lociPerFile + " (default))\n" +
				"   (12) width of the window (i.e. win=" + window + " (default))\n" +
				"   (13) use linux line endings (i.e. -linux (not the default))\n" +
				" Note: RFGLS format is sample dominant matrix w/ separate map and ped file\n" +
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("cnv=")) {
				cnvFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ped=")) {
				pedFilename = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outputFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("del=")) {
				includeDele = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("dup=")) {
				includeDupl = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("ord=")) {
				ordered = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("coll=")) {
				collapsed = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("homozygousOnly=")) {
				homozygous = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("excludeMonomorphicLoci=")) {
				homozygous = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("lociPerFile=")) {
				lociPerFile = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("win=")) {
				window = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-linux")) {
				endOfLine = "\n";
				numArgs--;
			} else if (args[i].startsWith("format=")) {
				fileFormat = ext.parseStringArg(args[i], MATRIX_FORMAT);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}

		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		//TODO To remove this block of code, which is for testing only.
//		cnvFilename = "C:/projects/Gedi/data/filtered.cnv";
//		pedFilename = null;
//		pedFilename = "C:/projects/Geti/filtered.fam";
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "ordered", true, true, true, false, markersPerFile, window);
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "ordered_collapsed", true, true, true, true, markersPerFile, window);
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "del_only", true, false, true, true, markersPerFile, window);
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "dup_only", false, true, true, false, markersPerFile, window);
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "unordered", true, true, false, true, markersPerFile, window);

//		cnvFilename = "C:/projects/Gedi/data/filtered.cnv";
//		cnvFilename = "D:/data/GEDI/penn_results/rfgls/data/conf15_usedFilteredRare.cnv";
//		lociPerFile = 1000;
//		pedFilename = ext.parseDirectoryOfFile(cnvFilename)+"IQ_IDs4Nathan121114.txt";
//		endOfLine = "\n";
//		rfglsOutput = true;
//		export(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "../results/deletionsAndDuplications", endOfLine, rfglsOutput, true, true, false, true, false, true, lociPerFile, window);
//		export(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "../results/deletionsOnly", endOfLine, rfglsOutput, true, false, false, true, false, true, lociPerFile, window);
//		export(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "../results/homozygousDeletionsOnly", endOfLine, rfglsOutput, true, false, false, true, true, true, lociPerFile, window);
//		System.exit(1);

		try {
			int fileAction = fileChooser.showOpenDialog(null);
			if (fileAction == JFileChooser.APPROVE_OPTION) {
				input = fileChooser.getSelectedFile();
				output = input.getPath() + "cnvBySample.txt";
				//TODO this line should be further elaborated.
				export(input.toString(), null, output, endOfLine, fileFormat, true, true, true, false, false, false, lociPerFile, window, new Logger());
			} else {
				export(cnvFilename, pedFilename, outputFilename, endOfLine, fileFormat, includeDele, includeDupl, ordered, collapsed, homozygous, excludeMonomorphicLoci, lociPerFile, window, new Logger());
//				cnvBySample("C:/projects/Geti/filtered.cnv", null, "C:/projects/Geti/cnvBySample.txt", endOfLine, rfglsOutput, true, false, true, false, saveIntermediateFiles, markersPerFile);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
