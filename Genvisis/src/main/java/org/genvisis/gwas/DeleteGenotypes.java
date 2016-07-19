package gwas;

import java.io.*;
import java.util.*;

import common.*;

public class DeleteGenotypes {
	public static final String WINDOWS_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\mosaics\\original\\";
	public static final String LINUX_DIRECTORY = "/archive/parkinsons/gwas/Genotype and intensity data files/";
	public static final String BACKUP_DIRECTORY = "backup/";
	public static final String OUTPUT_DIRECTORY = "updated/";
	public static final String LOOKUP = "lookup.txt";

	@SuppressWarnings("unchecked")
	public static void deleteGenotypes(String filename) {
		BufferedReader reader;
		PrintWriter writer = null;
		String[] line, marks;
		int count;
		String id = "", dir, trav;
		Hashtable<String,Hashtable<String,String>> hash = new Hashtable<String,Hashtable<String,String>>();
		Hashtable<String,String> markers, markerList;
		Hashtable<String,Vector<String>> lookup;
		Vector<String> sources;

		try {
			reader = new BufferedReader(new FileReader(filename));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (line[0].equals("for")) {
					markerList = new Hashtable<String,String>();
					while (reader.ready()) {
						markerList.put(reader.readLine().split("[\\s]+")[0], "");
					}
					for (int i = 1; i<line.length; i++) {
						hash.put(line[i], markerList);
					}
				} else {
					if (line[0].indexOf("__")>0) {
						line[0] = line[0].substring(0, line[0].indexOf("__"));
					}
					HashVec.addToHashHash(hash, line[0], line[1], "");
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		if (new File(WINDOWS_DIRECTORY).exists()) {
			dir = WINDOWS_DIRECTORY;
		} else if (new File(LINUX_DIRECTORY).exists()) {
			dir = LINUX_DIRECTORY;
		} else {
			dir = "";
			System.err.println("Error - could not resolve directory to parse (none of the following worked)");
			System.err.println(WINDOWS_DIRECTORY);
			System.err.println(LINUX_DIRECTORY);
			System.exit(1);
		}

		if (!new File(LOOKUP).exists()) {
			createLookupTable(1);
		}

		lookup = HashVec.loadFileToHashVec(LOOKUP, 1, new int[] {0}, "\t", false, false);

		new File(OUTPUT_DIRECTORY).mkdir();

		marks = HashVec.getKeys(hash);
		System.out.println(ext.getTime());
		try {
			for (int i = 0; i<marks.length; i++) {
				sources = lookup.get(marks[i]);
				if (sources==null) {
					System.err.println("Error - Could not find '"+marks[i]+"' in the lookup file. Either this is not a valid local_id or the lookup table is outdated (delete and I'll create a new one) or the directory I'm pointed too ('"+dir+"') is not the one you wanted.");
					System.exit(1);
				}
				for (int j = 0; j<sources.size(); j++) {
					markers = (Hashtable<String,String>)hash.get(marks[i]).clone();
					if (!new File(dir+sources.elementAt(j)).exists()) {
						System.err.println("Error - Could not find file '"+dir+sources.elementAt(j)+"'. Either this lookup table is outdated (delete and I'll create a new one) or the directory I'm pointed too ('"+dir+"') is not the one you wanted.");
						System.exit(1);
					}
					try {
						reader = new BufferedReader(new FileReader(Files.backup(sources.elementAt(j), dir, dir+BACKUP_DIRECTORY)));
						writer = new PrintWriter(new FileWriter(dir+sources.elementAt(j)));
						for (int k = 0; k<11; k++) {
							writer.println(reader.readLine());
						}
						count = 0;
						while (reader.ready()) {
							line = reader.readLine().split(",");
							trav = line[2].substring(0, line[2].indexOf("@"));
							if (count==0) {
								id = trav;
								if (!marks[i].startsWith(id)) {
									System.err.println("Error - incorrect match up. Is the lookup table incorrect? Opened '"+dir+sources.elementAt(j)+"' and found '"+id+"' instead of '"+marks[i]+"'");
									System.exit(1);
								}
							}
							if (markers.containsKey(line[0])) {
								line[4] = "0.0000";
								for (int k = 5; k<=12; k++) {
									line[k] = "-";
								}
								markers.remove(line[0]);
							}
							writer.println(Array.toStr(line, ","));
							count++;
						}
						reader.close();
						writer.close();
					} catch (IOException ioe) {
						System.err.println("Error reading file '"+dir+sources.elementAt(j)+"'");
						System.exit(2);
					}

				}
			}
			System.out.println(ext.getTime());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void createLookupTable(int serial) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String trav;
		Vector<String> v = new Vector<String>();
		int version;
		String dir, serialNumber;

		if (new File(WINDOWS_DIRECTORY).exists()) {
			dir = WINDOWS_DIRECTORY;
		} else if (new File(LINUX_DIRECTORY).exists()) {
			dir = LINUX_DIRECTORY;
		} else {
			dir = "";
			System.err.println("Error - could not resolve directory to parse (none of the following worked)");
			System.err.println(WINDOWS_DIRECTORY);
			System.err.println(LINUX_DIRECTORY);
			System.exit(1);
		}

		if (!new File(dir).exists()) {
			System.err.println("Error - could not resolve directory:");
			System.err.println(dir);
			System.exit(1);
		}

		File[] files = new File(dir).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".csv");
			}
		});

		try {
			writer = new PrintWriter(new FileWriter(LOOKUP));
			for (int i = 0; i<files.length; i++) {
				try {
					System.out.println(files[i].getName());
					reader = new BufferedReader(new FileReader(files[i]));
					for (int j = 0; j<11; j++) {
						reader.readLine();
					}
					line = reader.readLine().split(",");
					serialNumber = line[2];
					trav = line[2].substring(0, line[2].indexOf("@"));

					version = 0;
					if (v.contains(trav+(version==0?"":"__"+version))) {
						version++;
					}
					trav += (version==0?"":"__"+version);
					if (serial==1) {
						writer.println(files[i].getName()+"\t"+trav);
					} else if (serial==2) {
						writer.println(files[i].getName()+"\t"+trav+"\t"+serialNumber);
					} else {
						System.err.println("Error - don't know whether to include serial numbers or not");
						System.exit(1);
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[i].getName()+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[i].getName()+"\"");
					System.exit(2);
				}
			}
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "deleteThese.txt";
		int lookup = 0;

		String usage = "\\n"+
		"park.gwa.DeleteGenotypes requires 0-1 arguments\n"+
		"   (1) filename (i.e. file="+filename+" (default))\n"+
		"   (2) --createLookupWithoutSerialNumbers (not the default)\n"+
		"   (3) --createLookupWithSerialNumbers (not the default)\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("--createLookupWithout")) {
				lookup = 1;
				numArgs--;
			} else if (args[i].startsWith("--createLookupWith")) {
				lookup = 2;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (lookup>0) {
				createLookupTable(lookup);
			} else {
				deleteGenotypes(filename);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
