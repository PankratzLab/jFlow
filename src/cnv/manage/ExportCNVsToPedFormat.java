package cnv.manage;

import java.io.*;
import java.util.*;

import javax.swing.JFileChooser;
import filesys.Segment;
import cnv.var.CNVariant;
import common.*;

public class ExportCNVsToPedFormat {

	/*
	 * Convert a cnv data file into a ped file
	 */
	public static void export(String cnvFilename, String pedFilename, String outputRoot, String endOfLine, boolean includeDele, boolean includeDupl, boolean ordered, boolean collapsed, boolean homozygous, boolean excludeMonomorphicLoci, int markersPerFile, int win) {
		PrintWriter writer;
		CNVariant[] cnvs;
		Hashtable<String,String> allChrPosHash;
		Hashtable<String,String> sampleListHashFromCnvOrPedData;
		String[] allChrPosKeys;
		byte[][] currentCNs;
		Segment[] allChrPosSegs, currentChrPosSegs;
		String[] line;
		byte currentCN;
//		CNVariantHash cnvHash;
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

		System.out.println("Generating files for "+outputRoot);

		time = new Date().getTime();
		if (pedFilename != null) {
			sampleListHashFromCnvOrPedData = HashVec.loadFileToHashString(pedFilename, new int[] {0,1}, new int[] {-7}, pedFilename.endsWith(".csv"), "\t", false, false, false);
		} else {
			sampleListHashFromCnvOrPedData = null;
		}
		
		cnvs = CNVariant.loadPlinkFile(cnvFilename, sampleListHashFromCnvOrPedData, false);
		System.out.println("Loaded file in " + ext.getTimeElapsed(time));
		
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
			allChrPosSegs[i] = new Segment(Byte.parseByte(line[0]), Integer.parseInt(line[1])-win, Integer.parseInt(line[1])+win);
		}
		System.out.println("Generated hashtable of positions in " + ext.getTimeElapsed(time));

		time = new Date().getTime();
		allChrPosSegs = Segment.sortSegments(allChrPosSegs);
		System.out.println("Sorted positions in " + ext.getTimeElapsed(time));

//		cnvHash = new CNVariantHash(filename, structureType, jar)
		
		time = new Date().getTime();
		cnvs = CNVariant.sort(cnvs);
		System.out.println("Sorted CNVariants in " + ext.getTimeElapsed(time));

		
		cnvIndex = 0;
		countValidLoci = 0;
		writer = null;
		outputFilename = "very first file";
		for (int startPosition = 0; startPosition < allChrPosSegs.length; startPosition+=markersPerFile) {
			numMarkers = Math.min(markersPerFile, allChrPosSegs.length-startPosition);
			System.out.print("Calculating CN for positions "+(startPosition+1)+" through "+(startPosition+numMarkers));

			currentChrPosSegs = new Segment[numMarkers];
			for (int i = 0; i < currentChrPosSegs.length; i++) {
				currentChrPosSegs[i] = allChrPosSegs[startPosition+i];
			}
			
//			System.out.println("positions.length: "+positions.length+"\tsampleHash.size(): "+sampleHash.size()+"\tmatrix size: "+positions.length*sampleHash.size()+"\theap Size: "+Runtime.getRuntime().totalMemory());
			currentCNs = new byte[currentChrPosSegs.length][sampleListHashFromCnvOrPedData.size()];

			while (cnvs[cnvIndex].getChr() < currentChrPosSegs[0].getChr() || (cnvs[cnvIndex].getChr() == currentChrPosSegs[0].getChr() && cnvs[cnvIndex].getStop() < currentChrPosSegs[0].getStart())) {
				cnvIndex++;
			}
//			System.out.println("Starting at "+cnvs[startCNV].getUCSClocation());

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


//			if (saveIntermediateFiles) {
//				time = new Date().getTime();
//				new SerialByteMatrix(result).serialize(cnvFilename+"."+startPosition+"_bmatrix.ser");
//				new SerialStringArray(ids).serialize(cnvFilename+".ids.ser");
//				new SerialStringArray(positionNames).serialize(cnvFilename+".positions.ser");
//				System.out.println("Serialized intermediate files in " + ext.getTimeElapsed(time));
//			}

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
					if (homozygous && Math.abs(currentCN) == 1) {
						currentCN = 0;
					}
					if (collapsed && Math.abs(currentCN) == 2) {
						currentCN /= 2;
					}
					currentCNs[i][j] = currentCN;
				}
			}

			time = new Date().getTime();
			try {
				for (int i = 0; i<currentChrPosSegs.length; i++) {
					if (!excludeMonomorphicLoci || Array.min(currentCNs[i]) < 0 || Array.max(currentCNs[i]) > 0) {
						if (countValidLoci % markersPerFile == 0) {
							if (writer != null) {
								writer.close();
							}
							outputFilename = outputRoot+(allChrPosSegs.length<markersPerFile?"":"_"+countValidLoci)+".dat";
//							System.out.println("Beginning export to "+outputFilename);
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
					}
				}
			} catch (Exception e) {
				System.err.println("Error writing to '" + outputFilename + "'");
				e.printStackTrace();
			}
			System.out.println("    ...finished in " + ext.getTimeElapsed(time));
		}
		if (writer != null) {
			writer.close();
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String cnvFilename = "penncnv.cnv";
		String pedFilename = "pedigree.dat";
		String outputFilename = "cnv_matrix";
		JFileChooser fileChooser = new JFileChooser();
		File input;
		String output;
		boolean includeDele = true;
		boolean includeDupl = true;
		boolean ordered = true;
		boolean collapsed = false; 
		boolean homozygous = false;
		boolean excludeMonomorphicLoci = false;
		int markersPerFile = 5000;
		int window = 0;
		String endOfLine = "\r\n";

		String usage = "\n" +
				"cnv.analysis.CnvBySample requires the following arguments\n" +
				"   (1) cnv filename (i.e. cnv=" + cnvFilename + " (default))\n" +
				"   (2) pedegree filename (i.e. ped=" + pedFilename + " (default))\n" +
				"   (3) output filename (i.e. out=" + outputFilename + " (default))\n" +
				"   (4) to include Deletion or not (i.e. del=" + includeDele + " (default))\n" +
				"   (5) to include Duplication or not (i.e. dup=" + includeDupl + " (default))\n" +
				"   (6) to use ordered value or not (i.e. ord=" + ordered + " (default))\n" +
				"   (7) to use collapsed value or not (i.e. coll=" + collapsed + " (default))\n" +
				"   (8) to only use homozygous variants (i.e. homozygousOnly=" + homozygous + " (default))\n" +
				"   (9) to only use homozygous variants (i.e. excludeMonomorphicLoci=" + excludeMonomorphicLoci + " (default))\n" +
				"   (10) number of markers per file in the output file (i.e. markersPerFile=" + markersPerFile + " (default))\n" +
				"   (11) width of the window (i.e. win=" + window + " (default))\n" +
				"   (12) use linux line endings (i.e. -linux (not the default))\n" +
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("cnv=")) {
				cnvFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ped=")) {
				pedFilename = args[i].split("=")[1];
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
			} else if (args[i].startsWith("markersPerFile=")) {
				markersPerFile = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("win=")) {
				window = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-linux")) {
				endOfLine = "\n";
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
		cnvFilename = "C:/projects/Gedi/data/filtered.cnv";
		pedFilename = null;
//		pedFilename = "C:/projects/Geti/filtered.fam";
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "ordered", true, true, true, false, markersPerFile, window);
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "ordered_collapsed", true, true, true, true, markersPerFile, window);
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "del_only", true, false, true, true, markersPerFile, window);
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "dup_only", false, true, true, false, markersPerFile, window);
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "unordered", true, true, false, true, markersPerFile, window);

		endOfLine = "\n";
//		cnvBySample(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "../results/deletionsAndDuplications", endOfLine, true, true, false, true, false, true, markersPerFile, window);
		export(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "../results/deletionsOnly", endOfLine, true, false, false, true, false, true, markersPerFile, window);
		export(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "../results/homozygousDeletionsOnly", endOfLine, true, false, false, true, true, true, markersPerFile, window);
		System.exit(1);

		try {
			int fileAction = fileChooser.showOpenDialog(null);
			if (fileAction == JFileChooser.APPROVE_OPTION) {
				input = fileChooser.getSelectedFile();
				output = input.getPath() + "cnvBySample.txt";
				//TODO this line should be further elaborated.
				export(input.toString(), null, output, endOfLine, true, true, true, false, false, false, markersPerFile, window);
			} else {
				export(cnvFilename, pedFilename, outputFilename, endOfLine,includeDele, includeDupl, ordered, collapsed, homozygous, excludeMonomorphicLoci, markersPerFile, window);
//				cnvBySample("C:/projects/Geti/filtered.cnv", null, "C:/projects/Geti/cnvBySample.txt", endOfLine, true, false, true, false, saveIntermediateFiles, markersPerFile);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
