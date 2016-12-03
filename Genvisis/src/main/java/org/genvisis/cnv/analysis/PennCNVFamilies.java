package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.CLI;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

/**
 * Utility for building input list for penncnv trio
 * (http://penncnv.openbioinformatics.org/en/latest/user-guide/trio/) and joint
 * (http://penncnv.openbioinformatics.org/en/latest/user-guide/joint/) analysis.
 *
 */
public class PennCNVFamilies {

	public static void main(String... args) {
		final String trios = "trios";
		final String pennData = "pennData";
		final String jointSize = "joint";

		CLI c = new CLI(PennCNVFamilies.class);
		c.addArg(trios, "List of known trios", true);
		c.addArg(jointSize, "Joint chunk size. Allow ~1hr/trio. Will create [#samples]/jointSize input files.", true, CLI.Arg.NUMBER);
		c.addArg(pennData, "Directory containing zipped, exported penncnv data", true);

		c.parseWithExit(args);

		buildInputList(c.get(trios), c.get(pennData), c.getI(jointSize));
	}

	private static void buildInputList(String trios, String pennDir, int chunkSize) {
		final String out = "trioInput.txt";
		try {
			BufferedReader r = Files.getAppropriateReader(trios);
			PrintWriter writeTrios = Files.getAppropriateWriter(out);

			int jointChunk = 1;
			PrintWriter writeJoints = Files.getAppropriateWriter("jointInput" + jointChunk + ".txt");

			String header = r.readLine();
			int[] dnaIdxs = ext.indexFactors(	new String[] {"DNA", "FA_DNA", "MO_DNA"}, header.split("\t"),
																				false, true);
			int sample = 0;

			while (r.ready()) {
				String[] line = r.readLine().split("\t");

				boolean valid = true;
				for (int i = 0; i < dnaIdxs.length; i++) {
					valid = valid && isValidDNA(line[dnaIdxs[i]]);
				}
				if (!valid) {
					continue;
				}
				sample++;
				StringBuilder sb = new StringBuilder();
				for (int i = 0; i < dnaIdxs.length; i++) {
					sb.append("`gunzip -c ");
					sb.append(pennDir);
					sb.append(line[dnaIdxs[i]]);
					sb.append(".gz`");
					if (i + 1 < dnaIdxs.length) {
						sb.append("\t");
					}
				}
				writeTrios.println(sb.toString());
				writeJoints.println(sb.toString());
				if (sample >= chunkSize) {
					sample = 0;
					jointChunk++;
					writeJoints.flush();
					writeJoints.close();

					writeJoints = Files.getAppropriateWriter("jointInput" + jointChunk + ".txt");
				}
			}

			r.close();
			writeTrios.flush();
			writeTrios.close();
			writeJoints.flush();
			writeJoints.close();
		} catch (FileNotFoundException exc) {
			// TODO Auto-generated catch block
			exc.printStackTrace();
		} catch (IOException exc) {
			// TODO Auto-generated catch block
			exc.printStackTrace();
		}


		System.out.println("Finished writing penncnv trio inputs to: " + out);
	}

	private static boolean isValidDNA(String s) {
		return !s.isEmpty() && !s.equals("0") && !s.equals(".");
	}
}
