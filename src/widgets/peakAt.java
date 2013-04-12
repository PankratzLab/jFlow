package widgets;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import common.Files;
import common.ext;

public class peakAt {
	public static final int DEFAULT_NUM_LINES = 1000;

	public static void peak(String filename) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String temp;
		int numLines;
		boolean transpose = false;
		String[] line;
		LinkedList<String> ll;
		boolean tail, counting;
		int count;
		long time;
		String outputFilename;
		InputStreamReader isReader;
		
		outputFilename = "PeakAt_"+ext.rootOf(filename, true);

		temp = filename.substring(0, filename.length()-(".peakAt").length());
		tail = filename.indexOf(".tail.") > 0;
		counting = filename.indexOf(".count.") > 0 || filename.indexOf(".wc.") > 0;
		
		if (temp.indexOf(".")>0) {
			try {
				if (temp.substring(temp.lastIndexOf(".")+1).startsWith("trans")) {
					transpose = true;
					temp = temp.substring(0, temp.lastIndexOf("."));
				} else {
					transpose = false;
				}
			} catch (Exception e) {
				transpose = false;
			}

			try {
				numLines = Integer.parseInt(temp.substring(temp.lastIndexOf(".")+1));
				outputFilename = "PeakAt_"+(tail?"last":"first")+numLines+"_"+ext.rootOf(ext.rootOf(filename, true), true);
			} catch (NumberFormatException nfe) {
				numLines = DEFAULT_NUM_LINES;
				outputFilename = "PeakAt_"+(tail?"last":"first")+numLines+"_"+ext.rootOf(filename, true);
			}

		} else {
			numLines = DEFAULT_NUM_LINES;
		}

		try {
			writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(filename)+outputFilename));
			if (counting) {
				System.out.println("Counting the number of rows.");
	            time = new Date().getTime();
				count = Files.countLines(filename, true);
				writer.println(count);
				writer.println("Counted "+count+" rows in "+ext.getTimeElapsed(time));
			} else {
				isReader = null;
				if (outputFilename.endsWith(".gz")) {
					isReader = new InputStreamReader(new GZIPInputStream(new FileInputStream(filename)));
				} else if (outputFilename.endsWith(".zip")) {
					isReader = new InputStreamReader(new ZipInputStream(new FileInputStream(filename)));
				} else {
					isReader = new FileReader(filename);
				}

//				reader = new BufferedReader(new FileReader(filename));
				reader = new BufferedReader(isReader);
				if (transpose) {
					System.out.println("Taking the first "+numLines+" columns of all rows.");
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						for (int i = 0; i<Math.min(numLines, line.length); i++) {
							writer.print((i==0?"":"\t")+line[i]);
						}
						writer.println();
						writer.flush();
					}
				} else if (tail) {
					ll = new LinkedList<String>();
					while (reader.ready()) {
						ll.addLast(reader.readLine());
						if (ll.size() > numLines) {
							ll.removeFirst();
						}
					}
					while (!ll.isEmpty()) {
						writer.println(ll.removeFirst());
					}
				} else {
					for (int i = 0; i<numLines&&reader.ready(); i++) {
						writer.println(reader.readLine());
					}
				}
				reader.close();
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "peakAt.dat";

		String usage = "\n"+
		"widgets.peakAt requires 0-1 arguments\n"+
		"   (1) filename (i.e. file="+filename+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else {
				filename = args[i];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			peak(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
