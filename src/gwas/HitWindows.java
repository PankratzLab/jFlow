package gwas;

import java.io.*;
//import java.util.*;
import common.*;

public class HitWindows {
	public static void determine(String filename, String outfile, float indexThreshold, int windowMinSizePerSide, float windowExtensionThreshold) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, header;
		int count;
		String[] markerNames;
		byte[] chrs;
		int[] positions;
		float[] pvals;
		int startIndex, stopIndex, offset, minIndex;
		float minPval;
		int region;
		
		count = Files.countLines(filename, true);
		System.out.println("Parsing "+count+" lines");
		markerNames = new String[count];
		chrs = new byte[count];
		positions = new int[count];
		pvals = new float[count];
		try {
			reader = new BufferedReader(new FileReader(filename));
			header = reader.readLine().trim().split("[\\s]+");
			ext.checkHeader(header, new String[] {"MarkerName", "Chr", "Position", "Pvalue"}, false);
			count = 0;
			System.out.println("Parsing...");
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				markerNames[count] = line[0];
				chrs[count] = Byte.parseByte(line[1]);
				positions[count] = Integer.parseInt(line[2]);
				pvals[count] = line[3].equals(".")?999:Float.parseFloat(line[3]);
				count++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
		
		try {
			writer = new PrintWriter(new FileWriter(outfile));
			startIndex = -1;
			stopIndex = -1;
			region = 1;
			System.out.println("Starting search...");
			for (int i = 0; i < markerNames.length; i++) {
				if (pvals[i] < indexThreshold) {
					startIndex = i;
					minIndex = i;
					minPval = pvals[i];					
					offset = 0;
					while (chrs[startIndex] == chrs[startIndex-offset-1] && positions[startIndex] - windowMinSizePerSide <= positions[startIndex-offset-1]) {
						offset++;
						if (pvals[startIndex-offset] < windowExtensionThreshold) {
							startIndex -= offset;
							offset = 0;
						}
					}
					System.out.println(markerNames[i]+"\t"+region);
					
					stopIndex = i;
					offset = 0;
					while (chrs[stopIndex] == chrs[stopIndex+offset+1] && positions[stopIndex] + windowMinSizePerSide >= positions[stopIndex+offset+1]) {
						offset++;
						if (pvals[stopIndex+offset] < indexThreshold) {
							System.out.println(markerNames[stopIndex+offset]+"\t"+region);
						}
						if (pvals[stopIndex+offset] < windowExtensionThreshold) {
							stopIndex += offset;
							offset = 0;
						}
						if (pvals[stopIndex] < minPval) {
							minIndex = stopIndex;
							minPval = pvals[stopIndex];
						}
					}
//					System.out.println(markerNames[minIndex]+"\t"+chrs[minIndex]+"\t"+positions[minIndex]+"\tchr"+chrs[startIndex]+":"+(positions[startIndex]-windowMinSizePerSide)+":"+(positions[stopIndex]+windowMinSizePerSide));
					writer.println(markerNames[minIndex]+"\t"+chrs[minIndex]+"\t"+positions[minIndex]+"\tchr"+chrs[startIndex]+":"+(positions[startIndex]-windowMinSizePerSide)+":"+(positions[stopIndex]+windowMinSizePerSide));
					i = stopIndex+offset;
					region++;
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + outfile);
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "HitWindows.dat";
		String outfile = "hits.out";
		float indexThreshold = (float)0.00000005;
		int windowMinSizePerSide = 500000;
		float windowExtensionThreshold = (float)0.00001;

		String usage = "\n" + 
		"gwas.HitWindows requires 0-1 arguments\n" + 
		"   (1) input filename (i.e. file=" + filename + " (default))\n" + 
		"   (1) filename of output (i.e. out=" + outfile + " (default))\n" + 
		"   (2) p-value threshold for index SNPs (i.e. indexThresh=" + indexThreshold + " (default))\n" + 
		"   (3) minimum num pb per side of window (i.e. minWinSize=" + windowMinSizePerSide + " (default))\n" + 
		"   (4) p-value threshold to extend the window (i.e. winThresh=" + windowExtensionThreshold + " (default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("indexThresh=")) {
				indexThreshold = ext.parseFloatArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("minWinSize=")) {
				windowMinSizePerSide = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("winThresh=")) {
				windowExtensionThreshold = ext.parseFloatArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		try {
//			determine(filename, outfile, indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("D:/mega/FromMike.11032011/pvals1.txt", "D:/mega/FromMike.11032011/pvals1.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("D:/mega/FromMike.11032011/pvals2.txt", "D:/mega/FromMike.11032011/pvals2.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("D:/mega/FromMike.11032011/fixed_together.txt", "D:/mega/FromMike.11032011/fixed_together.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("D:/mega/FromMike.11032011/all_together.txt", "D:/mega/FromMike.11032011/all_together.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("C:/CARe_data/conditionalMeta/uni_input.txt", "C:/CARe_data/conditionalMeta/uniqueRegions.out", (float)0.00001, windowMinSizePerSide, (float)0.0001);
			determine("D:/tWork/Consortium/Megas/input.txt", "D:/tWork/Consortium/Megas/uniqueRegions.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
