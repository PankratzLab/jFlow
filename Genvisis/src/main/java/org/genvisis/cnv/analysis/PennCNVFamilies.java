package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.gwas.PlinkMendelianChecker;

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
		final String sexSpecific = "sex";

		CLI c = new CLI(PennCNVFamilies.class);
		c.addArg(trios, "List of known trios", true);
		c.addArg(jointSize,
						 "Joint chunk size. Allow ~1hr/trio. Will create [#samples]/jointSize input files.",
						 true, CLI.Arg.NUMBER);
		c.addArg(pennData, "Directory containing zipped, exported penncnv data", true);
		c.addFlag(sexSpecific, "Use samples from sex-specific subdirectores?");

		c.parseWithExit(args);

		buildInputList(c.get(trios), c.get(pennData), c.getI(jointSize), c.has(sexSpecific));
	}

	private static void buildInputList(String trios, String pennDir, int chunkSize,
																		 boolean sexSpecific) {
		final String out = "trioInput.txt";
		int jointChunk = 1;
		try {
			BufferedReader r = Files.getAppropriateReader(trios);
			PrintWriter writeTrios = Files.getAppropriateWriter(out);

			PrintWriter writeJoints = Files.getAppropriateWriter("jointInput" + jointChunk + ".txt");

			String[] header = r.readLine().split("\t");
			int[] dnaIdxs = ext.indexFactors(new String[] {"FA_DNA", "MO_DNA", "DNA"}, header, false,
																			 true);

			int sexIndex = ext.indexOfStr("SEX", header, false, true);

			int sample = 0;
			int[] sexes = {1, 2, 0};

			while (r.ready()) {
				String[] line = r.readLine().split("\t");

				boolean valid = true;
				for (int i = 0; i < dnaIdxs.length; i++) {
					valid = valid && PlinkMendelianChecker.isValidDNA(line[dnaIdxs[i]]);
				}
				if (!valid) {
					continue;
				}
				sample++;
				StringBuilder sb = new StringBuilder();
				sexes[2] = Integer.parseInt(line[sexIndex]);
				for (int i = 0; i < dnaIdxs.length; i++) {
					sb.append("`gunzip -c ");
					String dir = pennDir;
					if (sexSpecific) {
						dir += "sexSpecific" + File.separator;
						dir += sexes[i] == 1 ? "male" : "female";
						dir += File.separator;
					}

					String pennDataFile = new StringBuilder().append(dir).append(line[dnaIdxs[i]])
																									 .append(".gz").toString();


					valid = valid && new File(pennDataFile).exists();

					sb.append(pennDataFile);
					sb.append("`");
					if (i + 1 < dnaIdxs.length) {
						sb.append("\t");
					}
				}

				if (!valid) {
					continue;
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
		System.out.println("Wrote " + jointChunk + " joint input files");
	}

}
