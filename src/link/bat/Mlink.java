package link.bat;

import java.io.*;
import java.util.*;

import common.*;
import link.LinkageMap;

public class Mlink {
	public static final double DEFAULT_FREQ = 0.005;

	public static final double[] DEFAULT_MODEL = {0.03, 0.80, 0.80};

	// public static final double DEFAULT_FREQ = 0.01;
	// public static final double[] DEFAULT_MODEL = {0.0, 1.0, 1.0};
	public static final int MAX_MARKERS = 49;

	public static final String[] THETAS = {" 0.000", " 0.010", " 0.050", " 0.100", " 0.200", " 0.300", " 0.400"};

	public static final boolean DEFAULT_PED = true;

	public static void doChrome(int start, int stop, String filename, boolean pedformat) {
		for (int i = start; i<=stop; i++) {
			doChrome(i, filename, pedformat);
		}
	}

	public static void doChrome(int chr, String filename, boolean pedformat) {
		double[][] models;

		if (filename==null) {
			doChrome(chr, DEFAULT_FREQ, DEFAULT_MODEL, pedformat, -1);
		} else {
			models = LinkageMap.parseModels(filename);
			for (int i = 0; i<models.length; i++) {
				doChrome(chr, models[i][0], new double[] {models[i][1], models[i][2], models[i][3]}, pedformat, i+1);
			}
		}
	}

	public static void doChrome(int chr, double freq, double[] model, boolean pedformat, int iteration) {
		PrintWriter writer = null;
		int numFiles;
		String chrome = ext.chrome(chr);
		String dir;
		LinkageMap map;
		String[] markerNames;

		dir = (iteration==-1?"":iteration+"/")+"chrom"+chrome+"/";
		map = new LinkageMap("map"+chrome+"mlink.dat");
		markerNames = map.getMarkerNames();
		map.alterPenetrance("map"+chrome+"mlink.dat", freq, model[0], model[1], model[2], false);

		numFiles = (int)(Math.ceil(markerNames.length/(double)MAX_MARKERS));
		breakApart(chr, dir, pedformat);

		if (pedformat) {
			try {
				writer = new PrintWriter(new FileWriter("batch"));
				writer.println("cd chrom"+chrome);
				for (int i = 0; i<numFiles; i++) {
					writer.println("be pedin"+chrome+"-"+(i+1));
					writer.println("./pedin"+chrome+"-"+(i+1)+" > /dev/null");
				}
				writer.println("cd ..");
				writer.close();
				Runtime.getRuntime().exec("be batch").waitFor();
				System.out.println("Running chromosome "+chr);
				System.out.println("(if this runs for what seems like forever, check that there are no Mendelian errors: pedcheck  -p "+"re_chrom"+chrome+".pre"+" -d "+"map"+chrome+".dat"+" -o error.out  -4)");
				Runtime.getRuntime().exec("./batch >> batch.log").waitFor();
				parseStreams(markerNames, map.getCumulativePositions(false), chr, dir);
			} catch (InterruptedException ie) {
				System.err.println("Error running file \""+"batch"+"\"");
				System.exit(2);
			} catch (IOException ioe) {
				if (ioe.getMessage().startsWith("Cannot run program")) {
					// System.err.println("Not in linux; analysis skipped (hope
					// you have stored results!)");
					System.err.println("Not in linux; analysis skipped and summary");
				} else {
					System.err.println("Error writing file \""+"batch"+"\"");
					ioe.printStackTrace();
					System.exit(2);
				}
			}
		}

	}

	public static void parseStreams(String[] markerNames, double[] cumulative_cM_positions, int chr, String dir) {
		BufferedReader reader = null;
		PrintWriter writer;
		String temp;
		String chrome = ext.chrome(chr);
		int numFiles = (int)(Math.ceil(markerNames.length/(double)MAX_MARKERS));
		double[][] lods = new double[markerNames.length][THETAS.length];

		try {
			for (int i = 0; i<numFiles; i++) {
				reader = new BufferedReader(new FileReader(dir+"stream"+chrome+"-"+(i+1)+".out"));
				for (int j = 0; j<(i==numFiles-1?(markerNames.length-(numFiles-1)*MAX_MARKERS):MAX_MARKERS); j++) {
					for (int k = 0; k<53; k++) {
						temp = reader.readLine();
						if ((k==2||k==39)&&!temp.equals("MLINK")) {
							System.err.println("Error - change in file format, expecting MLINK after a @");
						}
					}
					for (int k = 0; k<THETAS.length; k++) {
						temp = reader.readLine();
						if (!temp.equals(THETAS[k])) {
							System.err.println("Error - change in file format, expecting theta of '"+THETAS[k]+"', found '"+temp+"'");
						}
						reader.readLine();
						reader.readLine();
						lods[i*MAX_MARKERS+j][k] = Double.parseDouble(reader.readLine());
					}
					reader.readLine();
				}
				reader.close();
			}

			writer = new PrintWriter(new FileWriter("chrom"+chrome+".xls"));
			writer.println("Chr.\tpos\tMarker\t"+Array.toStr(THETAS)+"\tThetaMax");
			for (int i = 0; i<markerNames.length; i++) {
				writer.print(chr+"\t"+cumulative_cM_positions[i]+"\t"+markerNames[i]);
				for (int j = 0; j<THETAS.length; j++) {
					writer.print("\t"+ext.formDeci(lods[i][j], 2));
				}
				writer.println("\t"+Array.max(lods[i]));

			}

			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(2);
		}

	}

	public static int breakApart(int chr, String destination, boolean pedformat) {
		BufferedReader reader = null;
		PrintWriter[] writers;
		String[] line, pedinfo;
		String temp, trav;
		int numMarkers;
		String chrome = ext.chrome(chr);
		int numFiles;
		Hashtable<String,String> pedStruct2;

		if (!destination.endsWith("/")) {
			destination += "/";
		}
		new File(destination).mkdirs();

		try {
			reader = new BufferedReader(new FileReader("map"+chrome+"mlink.dat"));

			temp = reader.readLine();
			line = temp.trim().split("[\\s]+");
			numMarkers = Integer.valueOf(line[0]).intValue()-1;
			numFiles = (int)(Math.ceil(numMarkers/(double)MAX_MARKERS));
			writers = new PrintWriter[numFiles];
			for (int i = 0; i<numFiles; i++) {
				writers[i] = new PrintWriter(new FileWriter(destination+"map"+chrome+"-"+(i+1)+".dat"));
				line[0] = i==numFiles-1?""+(numMarkers-(numFiles-1)*MAX_MARKERS+1):""+(MAX_MARKERS+1);
				writers[i].println(Array.toStr(line, " "));
			}
			ext.writeToAll(reader.readLine(), writers);
			reader.readLine();
			for (int i = 0; i<numFiles-1; i++) {
				writers[i].println(Array.toStr(Array.stringArraySequence(MAX_MARKERS+1, ""), " "));
			}
			writers[numFiles-1].println(Array.toStr(Array.stringArraySequence(numMarkers-(numFiles-1)*MAX_MARKERS+1, ""), " "));
			for (int i = 0; i<(chr==23?5:4); i++) {
				ext.writeToAll(reader.readLine(), writers);
			}

			for (int i = 0; i<numMarkers; i++) {
				writers[(int)Math.floor(i/MAX_MARKERS)].println(reader.readLine());
				writers[(int)Math.floor(i/MAX_MARKERS)].println(reader.readLine());
			}
			ext.writeToAll(reader.readLine(), writers);
			line = reader.readLine().split("[\\s]+");
			for (int i = 0; i<line.length-3; i++) {
				if (i%MAX_MARKERS==0) {
					writers[(int)Math.floor(i/MAX_MARKERS)].print(line[0]);
				} else {
					writers[(int)Math.floor(i/MAX_MARKERS)].print(" "+line[i]);
				}
			}
			ext.writeToAll(" "+line[line.length-3]+" "+line[line.length-2]+" "+line[line.length-1], writers);
			ext.writeToAll(reader.readLine(), writers);

			for (int i = 0; i<numFiles; i++) {
				writers[i].close();
			}
			reader.close();

			if (!new File("re_chrom"+chrome+".pre").exists()) {
				System.err.println("Error - failed to find "+"re_chrom"+chrome+".pre"+" for processing");
				System.exit(1);
			}
			reader = new BufferedReader(new FileReader("re_chrom"+chrome+".pre"));
			for (int i = 0; i<numFiles; i++) {
				writers[i] = new PrintWriter(new FileWriter(destination+"chrom"+chrome+"-"+(i+1)+(pedformat?".ped":".pre")));
			}
			if (pedformat) {
				pedStruct2 = procPedStruct(new File("re_chrom"+chrome+".ped").exists()?"re_chrom"+chrome+".ped":searchForPedStruct(), true);
			} else {
				pedStruct2 = procPedStruct("re_chrom"+chrome+".pre", false);
			}
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				temp = pedStruct2.get(line[0]+"\t"+line[1]);
				if (temp==null) {
					System.err.println("Error - .ped file is out of sync with "+"re_chrom"+chrome+".pre (could not find individual "+line[0]+"-"+line[1]+")");
					System.exit(1);
				}
				pedinfo = temp.trim().split("[\\s]+");
				trav = pedinfo[0]+"\t"+pedinfo[pedformat?7:4];
				if (!trav.equals(line[0]+"\t"+line[4])) {
					System.err.println("Error - .ped file does not match "+"re_chrom"+chrome+".pre"+" (found '"+line[0]+"/"+line[4]+"', expecting '"+trav.replace('\t', '/')+"'); aborting");
					System.err.println("        This has happened before when makeped reorders the individuals, but we made an automatic lookup table, so they must not matchup");
					System.exit(1);
				}
				for (int j = 0; j<numFiles; j++) {
					writers[j].print(temp);
				}
				for (int j = 0; j<numMarkers; j++) {
					writers[(int)Math.floor(j/MAX_MARKERS)].print("\t"+line[6+j*2+0]+"\t"+line[6+j*2+1]);
				}
				ext.writeToAll(pedformat?"‘  Ped: "+line[0]+"  Per: "+line[1]:"", writers);
			}
			for (int i = 0; i<numFiles; i++) {
				writers[i].close();
			}
			reader.close();

			if (pedformat) {
				for (int i = 0; i<numFiles; i++) {
					pedinXX(destination+"pedin"+chrome+"-"+(i+1), i==numFiles-1?(numMarkers-(numFiles-1)*MAX_MARKERS):MAX_MARKERS, "chrom"+chrome+"-"+(i+1)+".ped", "map"+chrome+"-"+(i+1)+".dat", "final"+chrome+"-"+(i+1)+".out", "stream"+chrome+"-"+(i+1)+".out");
				}
			}

			return numFiles;
		} catch (IOException ioe) {
			System.err.println("Error processing chromosome "+chr);
			ioe.printStackTrace();
			return -1;
		}
	}

	public static String searchForPedStruct() {
		String[] files = new File(".").list(new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(".ped");
			}
		});

		if (files.length==0) {
			System.err.println("Error - could not find a single .ped file to use as a template");
			System.exit(1);
			return null;
		} else {
			System.out.println("FYI, using "+files[0]+" as the template for the .ped structure");
			return files[0];
		}
	}

	public static Hashtable<String,String> procPedStruct(String filename, boolean ped) {
		BufferedReader reader;
		String temp;
		String[] line;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		int count, spaces;
		boolean on;

		try {
			reader = new BufferedReader(new FileReader(filename));
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				on = false;
				count = spaces = 0;
				while (spaces<(ped?10:6)) {
					if ((temp.charAt(count)==' '||temp.charAt(count)=='\t')&&on==false) {
						spaces++;
						on = true;
					} else {
						on = false;
					}
					count++;
				}
				hash.put(line[line.length-3]+"\t"+line[line.length-1], temp.substring(0, count-1));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		return hash;
	}

	public static void pedinXX(String batchfile, int numMarkers, String pedfile, String mapfile, String finalfile, String streamfile) {
		PrintWriter writer = null;

		try {
			writer = new PrintWriter(new FileWriter(batchfile));
			writer.println(":");
			writer.println("#");
			writer.println("#\t**");
			writer.println("#\t**  LCF  Generated  Command  File");
			writer.println("#\t**");
			writer.println("#\t**  Version : V4.0-00 ( 7-Jan-1991) - UNIX/Bourne");
			writer.println("#\t**  Created : 14-Feb-101  15:57:11");
			writer.println("#\t**");
			writer.println("#");
			writer.println("");
			writer.println("trap \"");
			writer.println("      if [ $1x != 'nodelete'x -a $1x != 'NODELETE'x -a $2x != 'nodelete'x -a $2x != 'NODELETE'x ]");
			writer.println("      then");
			writer.println("        rm -f datafile.dat final.dat ipedfile.dat lsp.log lsp.stm lsp.tmp\\");
			writer.println("              outfile.dat pedfile.dat recfile.dat speedfile.dat stream.dat tempdat.dat\\");
			writer.println("              tempped.dat *.tpd *.tdt *.mpd *.mdt");
			writer.println("      fi;");
			writer.println("      exit");
			writer.println("     \" 0 1 2 3 15");
			writer.println("");
			writer.println("rm -f lsp.log");
			writer.println("rm -f lsp.stm");
			writer.println("rm -f lsp.tmp");
			writer.println("rm -f final.dat");
			writer.println("rm -f stream.dat");
			writer.println("rm -f outfile.dat");
			writer.println("rm -f recfile.dat");
			writer.println("");
			writer.println("if [ -f "+finalfile+" ]");
			writer.println("then");
			writer.println("  mv "+finalfile+" "+finalfile+".bak");
			writer.println("fi");
			writer.println("cp /dev/null "+finalfile+"");
			writer.println("if [ $? != '0' ]");
			writer.println("then");
			writer.println("  exit 1");
			writer.println("fi");
			writer.println("");
			writer.println("if [ -f "+streamfile+" ]");
			writer.println("then");
			writer.println("  mv "+streamfile+" stream.bak");
			writer.println("fi");
			writer.println("cp /dev/null "+streamfile+"");
			writer.println("if [ $? != '0' ]");
			writer.println("then");
			writer.println("  exit 1");
			writer.println("fi");
			writer.println("");
			for (int j = 2; j<=numMarkers+1; j++) {
				writer.println("echo");
				writer.println("echo");
				writer.println("echo '******************************************************************************'");
				writer.println("echo 'Run   "+(j-1)+" - MLINK : p1 p"+j+"'");
				writer.println("echo '******************************************************************************'");
				writer.println("echo");
				writer.println("cp "+pedfile+" "+pedfile+".tpd");
				writer.println("cp "+mapfile+" "+mapfile+".tdt");
				writer.println("lsp \\");
				writer.println("\t  mlink \\");
				writer.println("\t  "+pedfile+".tpd "+mapfile+".tdt \\");
				writer.println("\t  2 \\");
				writer.println("\t   1  "+j+" \\");
				writer.println("\t  0 \\");
				writer.println("\t  0 \\");
				writer.println("\t  1 2 1 6 \\");
				writer.println("\t  0.01 2 1 \\");
				writer.println("\t  0.05 2 1 \\");
				writer.println("\t  0.1 2 1 \\");
				writer.println("\t  0.2 2 1 \\");
				writer.println("\t  0.3 2 1 \\");
				writer.println("\t  0.4 2 1");
				writer.println("if [ $? = '0' -o $? = '1' ]");
				writer.println("then");
				writer.println("  cat lsp.log >> "+finalfile+"");
				writer.println("  cat lsp.stm >> "+streamfile+"");
				writer.println("  unknown");
				writer.println("  if [ $? = '0' ]");
				writer.println("  then");
				writer.println("    mlink");
				writer.println("    if [ $? = '0' ]");
				writer.println("    then");
				writer.println("      cat outfile.dat >> "+finalfile+"");
				writer.println("      cat stream.dat  >> "+streamfile+"");
				writer.println("    fi");
				writer.println("  fi");
				writer.println("fi");
				writer.println("");
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+mapfile+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error parsing batch file for "+pedfile+" and "+mapfile+"");
			System.exit(2);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		int start = 1;
		int stop = 1;
		int chr = -1;
		// int chr = 23;
		String filename = null;
		boolean pedformat = DEFAULT_PED;

		String usage = "\n"+"park.mlink requires 0-4 arguments\n"+"   (1) chromosome start number (i.e. start="+start+" (default))\n"+"   (2) chromosome stop number (i.e. stop="+stop+" (default))\n"+" OR simply a single chromosome number\n"+"   (1) chromosome number (i.e. 2)\n"+" PLUS\n"+"   (3) filename of models (i.e. file=models.dat (default is a dominant model)\n"+"   (4) use PED format (i.e. ped="+pedformat+" (default, required to run mlink)\n"+"\n"+"\n"+"   NOTE #1: you'll need both a .pre and .ped, as makped truncates after around 50 markers worth of data\n"+"   NOTE #2: makeped often reorders individuals, but thankfully you created a lookup function to address this\n"+"   NOTE #3: makeped treats the last line differently than the rest so you'll have to mannually delete the newline character until you come up with a work around\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("start=")) {
				start = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("stop=")) {
				stop = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ped=")) {
				pedformat = args[i].split("=")[1].equalsIgnoreCase("true");
				numArgs--;
			} else {
				if (chr==-1) {
					try {
						chr = Integer.parseInt(args[i]);
					} catch (NumberFormatException nfe) {
						System.err.println("Error - invalid chromosome number: "+args[i]);
						System.exit(1);
					}
					numArgs--;
				} else {
					System.err.println("Error - tried to specify more than one flagless option ("+chr+" and "+args[i]+")");
				}
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (chr!=-1) {
				doChrome(chr, filename, pedformat);
			} else {
				doChrome(start, stop, filename, pedformat);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
