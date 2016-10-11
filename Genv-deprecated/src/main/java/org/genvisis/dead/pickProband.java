package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.park.tools;

public class pickProband {
	public pickProband() throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		String trav;
		int pro, vpdpro;

		hash = tools.getBestPDdx();

		try {
			reader = tools.getNinfoReader(1, false);
			writer = new PrintWriter(new FileWriter("probands.csv"));
			writer.println("FamID,IndID,proband,VPDproband");
			reader.readLine();
			trav = "";
			pro = vpdpro = 1;
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (!line[0].equals(trav)) {
					pro = vpdpro = 1;
					trav = line[0];
				}

				writer.print(line[0] + "," + line[1]);
				if (pro == 1 && tools.isAffected(hash, line[0] + "\t" + line[1])) {
					writer.print(",1");
					pro = 0;
				} else {
					writer.print(",0");
				}
				if (vpdpro == 1	&& hash.containsKey(line[0] + "\t" + line[1])
						&& hash.get(line[0] + "\t" + line[1]).equals("VPD")) {
					writer.print(",1");
					vpdpro = 0;
				} else {
					writer.print(",0");
				}
				writer.println();
			}
			reader.close();
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error reading ninfo1 file");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;

		String usage = "\n" + "park.pickProband requires no arguments\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			new pickProband();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
