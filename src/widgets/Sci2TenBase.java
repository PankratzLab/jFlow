package widgets;

import java.io.*;

import common.*;

public class Sci2TenBase {
	public Sci2TenBase(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename+"-better.txt"));
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				for (int i = 0; i<line.length; i++) {
					writer.print((i==0?"":"\t")+ext.prettyP(line[i]));
				}
				writer.println();
			}
			reader.close();
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
		// String filename = "Sci2TenBase.txt";
		String filename = "makeMePretty.txt";

		String usage = "\n"+"park.Sci2TenBase requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default)\n"+"";

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
			new Sci2TenBase(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
