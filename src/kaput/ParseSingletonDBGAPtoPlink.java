package kaput;

import java.io.*;
import java.util.*;

import common.*;

public class ParseSingletonDBGAPtoPlink {
	// public static final String WINDOWS_DIRECTORY = "C:\\Documents and
	// Settings\\npankrat\\My Documents\\jProjects\\park\\singleton_dbGaP";
	public static final String WINDOWS_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\PD-singleton\\10_dbGaP\\allSamples_results\\whatgives\\";

	public static final String LINUX_DIRECTORY = "/home/npankrat/PD-singleton/dbgap/IND";

	public static final String PED_DIRECTORY = "peds/";

	public static final String MARKER_FILE = "markers.dat";

	public static final String PED_FILE = "pedfile.dat";

	public static final String TARGET_MARKERS = "markerListWithIndices.dat";

	public static final int ALLELES_COL = 6;

	public static final double GC_SCORE_CUTOFF = 0.25;

	public static void parse(String filename) {
		BufferedReader reader;
		PrintWriter writer = null;
		String[] line;
		int count;
		String id = "", trav;

		if (new File(WINDOWS_DIRECTORY).exists()) {
			trav = WINDOWS_DIRECTORY;
		} else if (new File(LINUX_DIRECTORY).exists()) {
			trav = LINUX_DIRECTORY;
		} else {
			trav = "";
			System.err.println("Error - could not resolve directory to parse (none of the following worked)");
			System.err.println(WINDOWS_DIRECTORY);
			System.err.println(LINUX_DIRECTORY);
			System.exit(1);
		}

		File[] files = new File(trav).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.startsWith("phg000012.case")||filename.startsWith("phg000012.control");
			}
		});

		new File(PED_DIRECTORY).mkdir();

		Vector<String> v = new Vector<String>();
		System.out.println(ext.getTime());
		try {
			for (int i = 0; i<files.length; i++) {
				System.out.println(files[i].getName());
				try {
					reader = new BufferedReader(new FileReader(files[i]));
					for (int j = 0; j<3; j++) {
						reader.readLine();
					}
					id = reader.readLine().split("[\\s]+")[1];
					writer = new PrintWriter(new FileWriter(PED_DIRECTORY+id+".ped"));
					writer.print(id);

					trav = files[i].getName();
					if (!trav.endsWith(trav)) {
						System.err.println("Error - sample name is '"+trav+"' within the file, but the filename ends with "+trav.substring(trav.lastIndexOf(".")+1));

					}
					reader.readLine();
					count = 0;
					while (reader.ready()) {
						line = reader.readLine().split(",");
						if (i==0) {
							v.add(line[2]+(line[2].equals(line[1])?"":"\t"+line[1]));
						} else if (!v.elementAt(count).split("\t")[0].equals(line[2])) {
							System.err.println("Found "+line[2]+" -- expecting "+v.elementAt(count).split("\t")[0]+" in file "+files[i].getName());
							System.exit(1);
						}
						if (line[ALLELES_COL-1].equals("NA")) { // this is
							// the
							// column,
							// not the
							// index, so
							// you need
							// to
							// subtract
							// 1
							writer.print("\t0\t0");
						} else {
							writer.print("\t"+line[ALLELES_COL-1].charAt(0)+"\t"+line[ALLELES_COL-1].charAt(1)); // this
							// is
							// the
							// column,
							// not
							// the
							// index,
							// so
							// you
							// need
							// to
							// subtract
							// 1
						}
						count++;
					}
					reader.close();
					writer.println();
					writer.close();
					if (i==0) {
						System.out.println("There are "+v.size()+" markers being read in");
						writer = new PrintWriter(new FileWriter(MARKER_FILE));
						for (int j = 0; j<v.size(); j++) {
							writer.println(v.elementAt(j));
						}
						writer.close();
					}
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[i].getName()+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[i].getName()+"\"");
					System.exit(2);
				}
			}
			System.out.println(ext.getTime());
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "ParseSingletonDBGAPtoPlink.dat";

		String usage = "\\n"+"park.gwa.ParseSingletonDBGAPtoPlink requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			parse(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
