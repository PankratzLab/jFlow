package cnv.analysis;

import java.io.*;
import java.util.*;

import javax.sound.sampled.Line;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

//import cnv.analysis.FilterCalls;
import cnv.filesys.FullSample;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.var.CNVariant;
import cnv.var.CNVariantHash;
import cnv.var.SampleData;
import common.*;
import filesys.Segment;
import filesys.SerialByteMatrix;
import filesys.SerialHash;
import filesys.SerialStringArray;

public class CnvBySampleID {
	
	public static void cnvBySampleId(String cnvFilename, String pedFilename, String outputRoot, boolean includeDele, boolean includeDupl, boolean ordered, boolean collapsed, int markersPerFile, int win) {
		CNVariant[] cnvs;
		Hashtable<String,String> chrPositionHash;
		Hashtable<String,String> sampleHash;
		String[] chrPosition;
		byte[][] result;
		Segment[] allPositions, positions;
		String[] line;
		byte tmp;
//		CNVariantHash cnvHash;
		long time;
		String[] keys;
		String[] ids;
		int startCNV;
		boolean done;
		String outputFilename;
		int numMarkers;
		String sample;
		int sampleIndex;
		int iter;

		System.out.println("Generating files for "+outputRoot);

		time = new Date().getTime();
		if (pedFilename != null) {
			sampleHash = HashVec.loadFileToHashString(pedFilename, new int[] {0,1}, new int[] {-7}, pedFilename.endsWith(".csv"), "\t", false, false, false);
		} else {
			sampleHash = null;
		}
		
		cnvs = CNVariant.loadPlinkFile(cnvFilename, sampleHash, false);
		System.out.println("Loaded file in " + ext.getTimeElapsed(time));
		
		//Generate chrPosition and sampleHash, to be used for the rows and columns of the final result.
		time = new Date().getTime();
		chrPositionHash = new Hashtable<String,String>();
		if (pedFilename==null) {
			sampleHash = new Hashtable<String,String>();
		}
		for (int i=0; i<cnvs.length; i++) {
			chrPositionHash.put(cnvs[i].getChr()+"\t"+cnvs[i].getStart(), "");
			chrPositionHash.put(cnvs[i].getChr()+"\t"+cnvs[i].getStop(), "");
			chrPositionHash.put(cnvs[i].getChr()+"\t"+(cnvs[i].getStop()+1), "");
			if (pedFilename==null && !sampleHash.containsKey(cnvs[i].getFamilyID()+"\t"+cnvs[i].getIndividualID())) {
				sampleHash.put(cnvs[i].getFamilyID()+"\t"+cnvs[i].getIndividualID(), sampleHash.size()+"");
			}
		}
		chrPosition = HashVec.getKeys(chrPositionHash);
		allPositions = new Segment[chrPosition.length];
		for (int i = 0; i < allPositions.length; i++) {
			line = chrPosition[i].split("\t");
			allPositions[i] = new Segment(Byte.parseByte(line[0]), Integer.parseInt(line[1])-win, Integer.parseInt(line[1])+win);
		}
		System.out.println("Generated hashtable of positions in " + ext.getTimeElapsed(time));

		time = new Date().getTime();
		allPositions = Segment.sortSegments(allPositions);
		System.out.println("Sorted positions in " + ext.getTimeElapsed(time));

//		cnvHash = new CNVariantHash(filename, structureType, jar)
		
		time = new Date().getTime();
		cnvs = CNVariant.sort(cnvs);
		System.out.println("Sorted CNVariants in " + ext.getTimeElapsed(time));

		
		startCNV = 0;
		for (int startPosition = 0; startPosition < allPositions.length; startPosition+=markersPerFile) {
			System.out.print("Exporting "+(startPosition+1)+" of "+allPositions.length+" marker positions");
			numMarkers = Math.min(markersPerFile, allPositions.length-startPosition);
			outputFilename = outputRoot+(allPositions.length<markersPerFile?"":"_"+startPosition)+".dat";

			positions = new Segment[numMarkers];
			for (int i = 0; i < positions.length; i++) {
				positions[i] = allPositions[startPosition+i];
			}
			
//			System.out.println("positions.length: "+positions.length+"\tsampleHash.size(): "+sampleHash.size()+"\tmatrix size: "+positions.length*sampleHash.size()+"\theap Size: "+Runtime.getRuntime().totalMemory());
			result = new byte[positions.length][sampleHash.size()];

			while (cnvs[startCNV].getChr() < positions[0].getChr() || (cnvs[startCNV].getChr() == positions[0].getChr() && cnvs[startCNV].getStop() < positions[0].getStart())) {
				startCNV++;
			}
//			System.out.println("Starting at "+cnvs[startCNV].getUCSClocation());

			done = false;				
			for (int i = startCNV; !done && i < cnvs.length; i++) {
				// determine column from sample hash
				sample = sampleHash.get(cnvs[i].getFamilyID()+"\t"+cnvs[i].getIndividualID());
				if (cnvs[i].getChr() > positions[positions.length-1].getChr() || (cnvs[i].getChr() == positions[positions.length-1].getChr() && cnvs[i].getStart() > positions[positions.length-1].getStop())) {
					done = true;
				} else if (sample != null) {
					sampleIndex = Integer.parseInt(sample);
					iter = Math.max(Segment.binarySearchForStartPositions(cnvs[i], positions), 0);
					while (iter < positions.length && positions[iter].overlaps(cnvs[i])) {
						result[iter][sampleIndex] = (byte) (cnvs[i].getCN()-2);
						iter++;
					}

				}
			}

			ids = new String[sampleHash.size()];
			keys = HashVec.getKeys(sampleHash, false, false);
			for (int i = 0; i < keys.length; i++) {
				ids[Integer.parseInt(sampleHash.get(keys[i]))] = ext.replaceAllWith(keys[i], "\t", "-");
			}


//			if (saveIntermediateFiles) {
//				time = new Date().getTime();
//				new SerialByteMatrix(result).serialize(cnvFilename+"."+startPosition+"_bmatrix.ser");
//				new SerialStringArray(ids).serialize(cnvFilename+".ids.ser");
//				new SerialStringArray(positionNames).serialize(cnvFilename+".positions.ser");
//				System.out.println("Serialized intermediate files in " + ext.getTimeElapsed(time));
//			}
			
			for (int j = 0; j < result.length; j++) {
				for (int i = 0; i < result[j].length; i++) {
					tmp = result[j][i];
					if (!includeDele && tmp < 0) {
						tmp = 0;
					}
					if (!includeDupl && tmp > 0) {
						tmp = 0;
					}
					if (!ordered) {
						tmp = (byte)Math.abs(tmp);
					}
					if (collapsed && Math.abs(tmp) == 2) {
						tmp /= 2;
					}
					result[j][i] = tmp;
				}
			}

			time = new Date().getTime();
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(outputFilename));
				writer.println("markerName\t"+Array.toStr(ids));
				for (int i = 0; i<positions.length; i++) {
					writer.print(positions[i].getChr()+":"+positions[i].getStart());
					for (int j=0; j<ids.length; j++) {
						writer.print("\t"+result[i][j]);
					}
					writer.println();
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to '" + outputFilename + "'");
				e.printStackTrace();
			}
			System.out.println("    ...finished in " + ext.getTimeElapsed(time));
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
		int markersPerFile = 5000;
		int window = 0;

		String usage = "\n" +
				"cnv.analysis.CnvBySampleID requires the following arguments\n" +
				"   (1) cnv filename (i.e. cnv=" + cnvFilename + " (default))\n" +
				"   (2) pedegree filename (i.e. ped=" + pedFilename + " (default))\n" +
				"   (3) output filename (i.e. out=" + outputFilename + " (default))\n" +
				"   (4) to include Deletion or not (i.e. del=" + includeDele + " (default))\n" +
				"   (5) to include Duplication or not (i.e. dup=" + includeDupl + " (default))\n" +
				"   (6) to use ordered value or not (i.e. ord=" + ordered + " (default))\n" +
				"   (7) to use collapsed value or not (i.e. coll=" + collapsed + " (default))\n" +
				"   (8) number of markers per file in the output file (i.e. markersPerFile=" + markersPerFile + " (default))\n" +
				"   (9) width of the window (i.e. win=" + window + " (default))\n" +
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
			} else if (args[i].startsWith("markersPerFile=")) {
				markersPerFile = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("win=")) {
				window = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}

		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		cnvFilename = "C:/projects/Geti/filtered.cnv";
		pedFilename = null;
//		pedFilename = "C:/projects/Geti/filtered.fam";
		cnvBySampleId(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "ordered", true, true, true, false, markersPerFile, window);
		cnvBySampleId(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "ordered_collapsed", true, true, true, true, markersPerFile, window);
		cnvBySampleId(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "del_only", true, false, true, true, markersPerFile, window);
		cnvBySampleId(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "dup_only", false, true, true, false, markersPerFile, window);
		cnvBySampleId(cnvFilename, pedFilename, ext.parseDirectoryOfFile(cnvFilename) + "unordered", true, true, false, true, markersPerFile, window);
		System.exit(1);

		try {
			int fileAction = fileChooser.showOpenDialog(null);
			if (fileAction == JFileChooser.APPROVE_OPTION) {
				input = fileChooser.getSelectedFile();
				output = input.getPath() + "cnvBySampleId.txt";
				cnvBySampleId(input.toString(), null, output, true, true, true, false, markersPerFile, window);
			} else {
				cnvBySampleId(cnvFilename, pedFilename, outputFilename, includeDele, includeDupl, ordered, collapsed, markersPerFile, window);
//				cnvBySampleId("C:/projects/Geti/filtered.cnv", null, "C:/projects/Geti/cnvBySampleId.txt", true, false, true, false, saveIntermediateFiles, markersPerFile);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
