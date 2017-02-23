package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.cnv.analysis.pca.BetaOptimizer;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class processChargeBetas {

	public static void main(String[] args) {
		String betaDir = "C:/data/misc/ChargeBetas/";
		String[] betaFiles = Files.list(betaDir, "", ".txt", true, false, true);
		Logger log = new Logger(betaDir + "log.log");
		for (String betaFile : betaFiles) {
			String[] top = Files.getFirstNLinesOfFile(betaFile, 2, log);
			if (!top[0].contains("NCBI dbGaP analysis accession")) {
				throw new IllegalArgumentException();
			}
			if (!top[1].contains("Name:")) {
				throw new IllegalArgumentException();
			}
			String acc = top[0].split("\t")[1].replaceAll(" ", "_");
			String type = top[1].split("\t")[1].replaceAll(" ", "_");
			String file = betaDir + acc + "_" + type.replaceAll("/", "_") + ".beta";
			String[] required = new String[] {"SNP ID", "Allele1", "Allele2", "|&beta;|", "P-value"};
			try {

				BufferedReader reader = Files.getAppropriateReader(betaFile);
				PrintWriter writer = new PrintWriter(new FileWriter(file));
				int[] indices = null;
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					if (line[0].equals("ID")) {
						writer.println(ArrayUtils.toStr(BetaOptimizer.BETA_HEADER));
						indices = ext.indexFactors(required, line, true, false);
					}

					else if (indices != null) {
						writer.println(ArrayUtils.toStr(ArrayUtils.subArray(line, indices)));
					}
				}
				writer.close();

				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + betaFile + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + betaFile + "\"");
				return;
			}
		}
	}

}
