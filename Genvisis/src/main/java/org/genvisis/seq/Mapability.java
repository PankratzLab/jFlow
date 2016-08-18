package org.genvisis.seq;

import java.io.*;

import org.genvisis.common.*;
import org.genvisis.filesys.Segment;

public class Mapability {
	public static final String[] REQS = {"Chr", "Position", "totalreads", "notduplicated", "brokenmatepairs", "mapqualitygtzero", "notproperlypaired"};
	
	public static void parseBeds(String filename) {
		PrintWriter writer;
		String[] header;
		String[][] data;
		int[] indices;
		int pos;
		double value, sum, mean, max;
		
		header = Files.getHeaderOfFile(filename, "\t", new Logger());
		data = HashVec.loadFileToStringMatrix(filename, true, Array.arrayOfIndices(header.length), false);
		indices = ext.indexFactors(REQS, header, false, true);
	
		max = -1;
		sum = 0;
		for (int i = 0; i < data.length; i++) {
			value = Double.parseDouble(data[i][indices[3]]);
			sum += value;
			if (value > max) {
				max = value;
			}
		}
		mean = sum / data.length;
		
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_coverage.bed"));
			writer.println("track name=totalUniqueReads description=\"Coverage relative to mean (mean=500)\" useScore=1");
			for (int i = 0; i < data.length; i++) {
				pos = Integer.parseInt(data[i][indices[1]]);
				writer.println("chr"+data[i][indices[0]]+" "+(pos-50)+" "+(pos+50)+" "+pos+" "+(int)Math.min(Double.parseDouble(data[i][indices[3]])/mean*500, 1000));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename)+"_coverage.bed");
			e.printStackTrace();
		}
		
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_brokenMatePairs.bed"));
			writer.println("track name=percentBrokenMatePairs description=\"Percent of unique reads with broken mate pairs (scaled from 0-1000) \" useScore=1");
			for (int i = 0; i < data.length; i++) {
				pos = Integer.parseInt(data[i][indices[1]]);
				writer.println("chr"+data[i][indices[0]]+" "+(pos-50)+" "+(pos+50)+" "+pos+" "+(int)(Double.parseDouble(data[i][indices[4]])/Double.parseDouble(data[i][indices[3]])*1000));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename)+"_coverage.bed");
			e.printStackTrace();
		}
		
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_mapQualityZero.bed"));
			writer.println("track name=percentMapQualityZero description=\"Percent of unique reads with map quality of zero (scaled from 0-1000) \" useScore=1");
			for (int i = 0; i < data.length; i++) {
				pos = Integer.parseInt(data[i][indices[1]]);
				writer.println("chr"+data[i][indices[0]]+" "+(pos-50)+" "+(pos+50)+" "+pos+" "+(int)((Double.parseDouble(data[i][indices[3]])-Double.parseDouble(data[i][indices[5]]))/Double.parseDouble(data[i][indices[3]])*1000));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename)+"_coverage.bed");
			e.printStackTrace();
		}

		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_notProperlyPaired.bed"));
			writer.println("track name=percentNotProperlyPaired description=\"Percent of unique reads that are not properly paired (scaled from 0-1000) \" useScore=1");
			for (int i = 0; i < data.length; i++) {
				pos = Integer.parseInt(data[i][indices[1]]);
				writer.println("chr"+data[i][indices[0]]+" "+(pos-50)+" "+(pos+50)+" "+pos+" "+(int)(Double.parseDouble(data[i][indices[6]])/Double.parseDouble(data[i][indices[3]])*1000));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename)+"_coverage.bed");
			e.printStackTrace();
		}
		
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_compositePoorQualty.bed"));
			writer.println("track name=CompositePoorQuality description=\"Composite Poor Quality = MapQualityZero + NotProperlyPaired\" useScore=1");
			for (int i = 0; i < data.length; i++) {
				pos = Integer.parseInt(data[i][indices[1]]);
				writer.println("chr"+data[i][indices[0]]+" "+(pos-50)+" "+(pos+50)+" "+pos+" "+(int)((Double.parseDouble(data[i][indices[3]])-Double.parseDouble(data[i][indices[5]])+Double.parseDouble(data[i][indices[6]]))/Double.parseDouble(data[i][indices[3]])*1000));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename)+"_coverage.bed");
			e.printStackTrace();
		}
	}
	
	public static void assignPass(String filename, String[] passes) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String[][] data;
		Segment[][] segs;
		Segment seg;
		boolean overlaps;
		
		segs = new Segment[passes.length][];
		for (int i = 0; i < passes.length; i++) {
			data = HashVec.loadFileToStringMatrix(passes[i], false, new int[] {0,1,2}, false);
			segs[i] = new Segment[data.length];
			for (int j = 0; j < data.length; j++) {
				segs[i][j] = new Segment("chr"+data[j][0]+":"+data[j][1]+"-"+data[j][2]);
			}
		}
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename+"_positions.xln"));
			line = reader.readLine().trim().split("[\\s]+");
			ext.checkHeader(line, new String[] {"Chr", "Position"}, new int[] {0,1}, false, new Logger(), true);
			for (int i = 0; i < passes.length; i++) {
				line = Array.insertStringAt(ext.rootOf(passes[i]), line, 2+i);
			}
			writer.println(Array.toStr(line));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				seg = new Segment("chr"+line[0]+":"+line[1]+"-"+(Integer.parseInt(line[1])+1));
				for (int i = 0; i < passes.length; i++) {
					overlaps = false;
					for (int j = 0; j < segs[i].length; j++) {
						if (seg.overlaps(segs[i][j])) {
							overlaps = true;
						}
					}
					line = Array.insertStringAt(overlaps?"1":"0", line, 2+i);
				}
				writer.println(Array.toStr(line));
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
//		String filename = "D:\\tWork\\SequencingProjectWithCIDR\\Coverage\\all.out";
		String filename = "positions.txt";
		boolean parse = false;
		String[] passes = new String[] {"1stPass.bed", "2ndPass.bed"};

		String usage = "\n" + 
		"seq.Mapability requires 0-1 arguments\n" + 
		"   (1) filename (i.e. file=" + filename + " (default))\n" + 
		"   (2) parse beds (i.e. -parse (not the default))\n" + 
		" OR\n" + 
		"   (2) assign passes to file (i.e. passes=passOne.bed,passTwo.bed (not the default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-parse")) {
				parse = true;
				numArgs--;
			} else if (args[i].startsWith("passes=")) {
				passes = args[i].split("=")[1].split(",");
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (parse) {
				parseBeds(filename);
			} else if (passes != null) {
				assignPass(filename, passes);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
