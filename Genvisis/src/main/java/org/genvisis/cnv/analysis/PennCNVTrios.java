package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.CLI;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

/**
 * Utility for building input list for penncnv trio analysis
 * (http://penncnv.openbioinformatics.org/en/latest/user-guide/trio/)
 *
 */
public class PennCNVTrios {

	public static void main(String... args) {
		final String trios = "trios";
		final String pennData = "pennData";
		CLI c = new CLI(PennCNVTrios.class);
		c.addArg(trios, "List of known trios", true);
		c.addArg(pennData, "Directory containing zipped, exported penncnv data", true);

		c.parseWithExit(args);

		buildInputList(c.get(trios), c.get(pennData));
	}

	private static void buildInputList(String trios, String pennDir) {
		try {
			BufferedReader r = Files.getAppropriateReader(trios);
			PrintWriter w = Files.getAppropriateWriter("trioInput.txt");
			String header = r.readLine();
			int[] dnaIdxs = ext.indexFactors(	new String[] {"DNA", "FA_DNA", "MO_DNA"}, header.split("\t"),
																				false, true);
			while (r.ready()) {
				String[] line = r.readLine().split("\t");

				boolean valid = true;
				for (int i=0; i<dnaIdxs.length; i++) {
					valid = valid && isValidDNA(line[dnaIdxs[i]]);
				}
				if (!valid) {
					continue;
				}
				StringBuilder sb = new StringBuilder();
				for (int i=0; i<dnaIdxs.length; i++) {
					sb.append("`gunzip -c ");
					sb.append(pennDir);
					sb.append(line[dnaIdxs[i]]);
					sb.append(".gz`");
					if (i+1 < dnaIdxs.length) {
						sb.append("\t");
					}
				}
				w.println(sb.toString());
			}

			r.close();
			w.flush();
			w.close();
		} catch (FileNotFoundException exc) {
			// TODO Auto-generated catch block
			exc.printStackTrace();
		} catch (IOException exc) {
			// TODO Auto-generated catch block
			exc.printStackTrace();
		}


	}

	private static boolean isValidDNA(String s) {
		return !s.isEmpty() && !s.equals("0") && !s.equals(".");
	}
}
