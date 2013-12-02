package gwas;

import java.io.*;
import java.util.Hashtable;
import java.util.Vector;

import common.*;
import filesys.Segment;
import filesys.SnpMarkerSet;

// currently crashes if first or last marker passes significance threshold or is in window of a SNP that does
public class HitWindows {
	public static void determine(String filename, String outfile, float indexThreshold, int windowMinSizePerSide, float windowExtensionThreshold) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, header;
		int count;
		String[] markerNames;
		byte[] chrs;
		int[] positions;
		double[] pvals;
		int startIndex, stopIndex, offset, minIndex;
		double minPval;
		int region;
		int numSig, numSuggestive;
		
		count = Files.countLines(filename, true);
		System.out.println("Parsing "+count+" lines");
		markerNames = new String[count];
		chrs = new byte[count];
		positions = new int[count];
		pvals = new double[count];
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
				pvals[count] = ext.isMissingValue(line[3])?999:Double.parseDouble(line[3]);
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
			writer.println("Region\tMarkerName\tChr\tPosition\tp-value\tRegion+Window\tRegionStart\tRegionStop\tNumSigMarkers\tNumSuggestiveMarkers\tNumTotalMarkers\tSizeOfRegion");
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
					numSig = numSuggestive = 1;
					while (startIndex-offset-1 >= 0 && chrs[startIndex] == chrs[startIndex-offset-1] && positions[startIndex] - windowMinSizePerSide*2 <= positions[startIndex-offset-1]) { // *2 required to ensure that there are no overlapping SNPs 500kb after last hit and 500kb before next hit is technically a 1M region that should be merged 
						offset++;
						if (pvals[startIndex-offset] < windowExtensionThreshold) {
							startIndex -= offset;
							offset = 0;
							numSuggestive++;
						}
					}
					System.out.println(markerNames[i]+"\t"+region+"\t"+numSig+"\t"+numSuggestive);
									
					stopIndex = i;
					offset = 0;
					while (stopIndex+offset+1 < markerNames.length && chrs[stopIndex] == chrs[stopIndex+offset+1] && positions[stopIndex] + windowMinSizePerSide*2 >= positions[stopIndex+offset+1]) {
						offset++;
						if (pvals[stopIndex+offset] < indexThreshold) {
//							System.out.println(markerNames[stopIndex+offset]+"\t"+region);
							numSig++;
						}
						if (pvals[stopIndex+offset] < windowExtensionThreshold) {
							stopIndex += offset;
							offset = 0;
							numSuggestive++;
						}
						if (pvals[stopIndex] < minPval) {
							minIndex = stopIndex;
							minPval = pvals[stopIndex];
						}
					}
					System.out.println(markerNames[minIndex]+"\t"+region+"\t"+numSig+"\t"+numSuggestive);
//					System.out.println(markerNames[minIndex]+"\t"+chrs[minIndex]+"\t"+positions[minIndex]+"\tchr"+chrs[startIndex]+":"+(positions[startIndex]-windowMinSizePerSide)+":"+(positions[stopIndex]+windowMinSizePerSide));
					writer.println(region+"\t"+markerNames[minIndex]+"\t"+chrs[minIndex]+"\t"+positions[minIndex]+"\t"+pvals[minIndex]+"\tchr"+chrs[startIndex]+":"+Math.max(positions[startIndex]-windowMinSizePerSide, 1)+"-"+(positions[stopIndex]+windowMinSizePerSide)+"\t"+positions[startIndex]+"\t"+positions[stopIndex]+"\t"+numSig+"\t"+numSuggestive+"\t"+(stopIndex-startIndex+1)+"\t"+(positions[stopIndex]-positions[startIndex]+1));
					
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
	
	public static void generateHitsLookup(String inputHits, int window, String outputFile, String mapfile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Hashtable<String, Vector<String>> hash;
		Vector<String> v;
		int count;
		long time;
		Segment[][] segs;
		SnpMarkerSet markerSet;
		Logger log;
		int[] indices;
		String[] header, traits, markerNames, chrPositions;
		Segment variant;
		
		log = new Logger();

		v = new Vector<String>();
		hash = new Hashtable<String, Vector<String>>();
		try {
			reader = new BufferedReader(new FileReader(inputHits));
			header = reader.readLine().trim().split("[\\s]+");
			indices = ext.indexFactors(new String[][] {{"Trait"}, Aliases.CHRS, Aliases.POSITIONS}, header, false, true, true, true);
			if (!Array.equals(indices, new int[] {0,1,2})) {
				log.reportError("Error - currently expecting format: Trait\tChr\tPosition");
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				HashVec.addIfAbsent(line[0], v);
				HashVec.addToHashVec(hash, line[0], line[1]+"\t"+line[1], false);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + inputHits + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + inputHits + "\"");
			System.exit(2);
		}
		
		traits = Array.toStringArray(v);
		segs = new Segment[traits.length][];
		for (int i = 0; i < traits.length; i++) {
			v = hash.get(traits[i]);
			segs[i] = new Segment[v.size()];
			for (int j = 0; j < segs[i].length; j++) {
				line = v.elementAt(j).split("[\\s]+");
				segs[i][j] = new Segment(Positions.chromosomeNumber(line[0]), Integer.parseInt(line[1])-window, Integer.parseInt(line[1])+window);
			}
		}
		
		
		markerSet = new SnpMarkerSet(mapfile, false, log);
		markerNames = markerSet.getMarkerNames();
		chrPositions = markerSet.getChrAndPositions();
		for (int m = 0; m < 10; m++) {
			line = chrPositions[m].split("[\\s]+");
			variant = new Segment(Byte.parseByte(line[0]), Integer.parseInt(line[1]), Integer.parseInt(line[1])+1);
			v = new Vector<String>();
			for (int i = 0; i < segs.length; i++) {
				for (int j = 0; j < segs[i].length; j++) {
					if (variant.overlaps(segs[i][j])) {
						v.add(traits[i]);
					}
				}
			}
			System.out.println(markerSet.getMarkerNames()[m]+"\t"+(v.size()==0?".":Array.toStr(Array.toStringArray(v), "/")));
		}
		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "input.dat";
		String outfile = "hits.out";
		float indexThreshold = (float)0.00000005;
		int windowMinSizePerSide = 500000; // 500kb each side is technically a 1M window until the next hit region, but we now take this into consideration in the main algorithm
		float windowExtensionThreshold = (float)0.00000005; // (float)0.00001;
		String knownHits = null;
		String map = "markers.dat";

		String usage = "\n" + 
		"gwas.HitWindows requires 0-1 arguments\n" + 
		"   (1) input filename (i.e. file=" + filename + " (default))\n" + 
		"   (2) filename of output (i.e. out=" + outfile + " (default))\n" + 
		"   (3) p-value threshold for index SNPs (i.e. indexThresh=" + indexThreshold + " (default))\n" + 
		"   (4) minimum num bp per side of window (i.e. minWinSize=" + windowMinSizePerSide + " (default))\n" + 
		"   (5) p-value threshold to extend the window (i.e. winThresh=" + windowExtensionThreshold + " (default))\n" + 
		" OR\n" + 
		"   (1) list of known hits, 3 columns=trait+chr+position (i.e. knownHits=filenameOfKnownHits.dat (not the default))\n" + 
		"   (2) window around hit to extend (i.e. minWinSize=" + windowMinSizePerSide + " (default))\n" + 
		"   (3) map file for lookup (i.e. map=" + map + " (default))\n" + 
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
			} else if (args[i].startsWith("knownHits=")) {
				knownHits = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("map=")) {
				map = ext.parseStringArg(args[i], null);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		generateHitsLookup("D:/ExomeChip/Hematology/00src/CHARGE-RBC/knownHits.dat", 200000, "D:/ExomeChip/Hematology/00src/CHARGE-RBC/hitLookup.dat", "D:/ExomeChip/Hematology/00src/CHARGE-RBC/ExomeChipV5_wMAF.csv");
		System.exit(1);
		
		try {
			if (knownHits != null) {
				generateHitsLookup(knownHits, windowMinSizePerSide, outfile, map);
			} else {
				determine(filename, outfile, indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
			}
//			determine("D:/mega/FromMike.11032011/pvals1.txt", "D:/mega/FromMike.11032011/pvals1.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("D:/mega/FromMike.11032011/pvals2.txt", "D:/mega/FromMike.11032011/pvals2.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("D:/mega/FromMike.11032011/fixed_together.txt", "D:/mega/FromMike.11032011/fixed_together.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("D:/mega/FromMike.11032011/all_together.txt", "D:/mega/FromMike.11032011/all_together.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("C:/CARe_data/conditionalMeta/uni_input.txt", "C:/CARe_data/conditionalMeta/uniqueRegions.out", (float)0.00001, windowMinSizePerSide, (float)0.0001);
//			determine("D:/tWork/Consortium/Megas/input.txt", "D:/tWork/Consortium/Megas/uniqueRegions.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("D:/mega/filtered/OnlyThoseInRef/input.dat", "D:/mega/filtered/OnlyThoseInRef/uniqueRegions.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("D:/mega/filtered/OnlyThoseInRef/originalWrongRE2_input.dat", "D:/mega/filtered/OnlyThoseInRef/originalWrongRE2_uniqueRegions.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
//			determine("D:/home/npankrat/jProjects/Ruhi/bioinf/input2.txt", "D:/home/npankrat/jProjects/Ruhi/bioinf/my.out.xln", 0.00001f, 500000, 0.0001f);
			
			
//			determine("/home/input.dat", "/home/uniqueRegions.out", indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
