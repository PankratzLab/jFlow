package org.genvisis.dead;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class RemoveMarkersFromResults {
	public static final String DEFAULT_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\gaw\\results\\";

	public static void remove(String dir, String filename, final String suffix) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, current, future;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		int snpIndex;

		File[] files = new File(dir.equals("")?".":dir).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(suffix);
			}
		});

		try {
			reader = new BufferedReader(new FileReader(dir+filename));
			line = reader.readLine().trim().split("[\\s]+");
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				hash.put(line[0], "");
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+filename+"\"");
			System.exit(2);
		}

		for (int i = 0; i<files.length; i++) {
			current = files[i].getAbsolutePath();
			if (current.contains("\\")) {
				new File(current.substring(0, current.lastIndexOf("\\"))+"\\backup").mkdirs();
				future = current.substring(0, current.lastIndexOf("\\"))+"\\backup"+current.substring(current.lastIndexOf("\\"));
			} else {
				new File(current.substring(0, current.lastIndexOf("/"))+"/backup").mkdirs();
				future = current.substring(0, current.lastIndexOf("/"))+"/backup"+current.substring(current.lastIndexOf("/"));
			}
			System.out.println(current);
			files[i].renameTo(new File(future));
			try {
				reader = new BufferedReader(new FileReader(future));
				writer = new PrintWriter(new FileWriter(current));

				temp = reader.readLine();
				writer.println(temp);
				line = temp.trim().split("[\\s]+");
				snpIndex = ext.indexOfStr("SNP", line);
				if (snpIndex==-1) {
					System.err.println("Error - could not find a column named \"SNP\" in file "+current);
				}
				while (reader.ready()) {
					temp = reader.readLine();
					line = temp.trim().split("[\\s]+");
					if (!hash.containsKey(line[snpIndex])) {
						writer.println(temp);
					}
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""+future+"\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \""+future+"\"");
				System.exit(2);
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "dropThese.txt";
		String suffix = "logistic";
		String dir = new File(DEFAULT_DIR).exists()?DEFAULT_DIR:"";

		String usage = "\\n"+"park.gwa.RemoveMarkersFromResults requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"   (2) suffix of files to drop from (i.e. suffix="+suffix+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("suffix=")) {
				suffix = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			remove(dir, filename, suffix);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
