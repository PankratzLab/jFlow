package one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import common.Files;
import common.Logger;
import common.ext;
import cnv.var.LocusSet;
import filesys.Segment;

public class filterByLocPval {

	public static void filter(String segFile, String locPvalFile, Logger log) {
		LocusSet<Segment> set = LocusSet.loadSegmentSetFromFile(segFile, 0, 1, 2, 0, true, true, 0, log);
		String output = ext.rootOf(locPvalFile, false) + ".union.txt";
		try {
			BufferedReader reader = Files.getAppropriateReader(locPvalFile);

			PrintWriter writer = new PrintWriter(new FileWriter(output));

			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("[\\s]+");
				Segment curSeg = new Segment(line[0].replaceAll("\\.\\.", "-"));
				
				if (set.getOverlappingIndices(curSeg) != null) {
					double maf = Double.parseDouble(line[1]);
					if (maf > 0 && maf < 0.05) {
						System.out.println(set.getLoci()[set.getOverlappingIndices(curSeg)[0]].getUCSClocation());
						System.out.println(curSeg.getUCSClocation());
//						try {
//							Thread.sleep(1000);
//						} catch (InterruptedException ie) {
//						}
						writer.println(line[2]);
					}
				}

			}
			writer.close();

			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + locPvalFile + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + locPvalFile + "\"");
			return;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String segFile = "D:/data/Project_Tsai_21_25_26_spector/meregedCaptureRegions.txt";
		String locPvalFile = "D:/data/Project_Tsai_21_25_26_spector/NORMALIZED_GC_CORRECTED/NORMALIZED_GC_CORRECTED/maf.p.loc.txt";

		Logger log = new Logger(ext.rootOf(locPvalFile, false) + ".log");
		filter(segFile, locPvalFile, log);
	}
}
