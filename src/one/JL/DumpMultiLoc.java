package one.JL;

import java.io.FileWriter;
import java.io.PrintWriter;

import common.Logger;

import filesys.GeneData;
import filesys.GeneTrack;

public class DumpMultiLoc {

	public static void dumpMultiLoc(String geneTrackFile, String outputFile, Logger log) {
		GeneTrack geneTrack = GeneTrack.load(geneTrackFile, false);
		GeneData[][] genes = geneTrack.getGenes();

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
			for (int i = 0; i < genes.length; i++) {
				for (int j = 0; j < genes[i].length; j++) {
					if (genes[i][j].getMultiLoc() > 0) {
						writer.println(genes[i][j].getChr() + "\t" + genes[i][j].getStart() + "\t" + genes[i][j].getStop());
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outputFile);
			log.reportException(e);
		}

	}

}
