package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;

public class filterByLocPval {

	public static void filter(String segFile, String locPvalFile, boolean mafOnly, Logger log) {
		LocusSet<Segment> set = LocusSet.loadSegmentSetFromFile(segFile, 0, 1, 2, 0, true, true, 0,
																														log);
		String output = ext.rootOf(locPvalFile, false) + (mafOnly ? ".maf.txt" : ".maf.union.txt");
		try {
			BufferedReader reader = Files.getAppropriateReader(locPvalFile);

			PrintWriter writer = Files.openAppropriateWriter(output);

			while (reader.ready()) {
				String[] line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
				Segment curSeg = new Segment(line[0].replaceAll("\\.\\.", "-"));

				if (set.getOverlappingIndices(curSeg) != null || mafOnly) {
					double maf = Double.parseDouble(line[1]);
					if (maf > 0 && maf < 0.05) {
						// System.out.println(set.getLoci()[set.getOverlappingIndices(curSeg)[0]].getUCSClocation());
						// System.out.println(curSeg.getUCSClocation());
						// try {
						// Thread.sleep(1000);
						// } catch (InterruptedException ie) {
						// }
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
		// int numArgs = args.length;
		String segFile = "D:/data/Project_Tsai_21_25_26_spector/meregedCaptureRegions.txt";
		String locPvalFile = "D:/data/Project_Tsai_21_25_26_spector/NORMALIZED_GC_CORRECTED/NORMALIZED_GC_CORRECTED/p.loc.txt";

		// String locPvalFile =
		// "D:/data/Project_Tsai_21_25_26_spector/NORMALIZED_GC_CORRECTED/NORMALIZED_GC_CORRECTED/maf.p.loc.filt.txt";
		//
		// String segFile =
		// "D:/data/Project_Tsai_21_25_26_spector/MergeAric/EPP_V_ARIC_hists_ARIC_FILTER.MAF_0.segments";
		Logger log = new Logger(ext.rootOf(locPvalFile, false) + ".log");
		filter(segFile, locPvalFile, true, log);
		filter(segFile, locPvalFile, false, log);

		segFile = "D:/data/Project_Tsai_21_25_26_spector/MergeAric/EPP_V_ARIC_hists_ARIC_FILTER.MAF_0.segments";
		locPvalFile = "D:/data/Project_Tsai_21_25_26_spector/NORMALIZED_GC_CORRECTED/NORMALIZED_GC_CORRECTED/p.loc.filt.txt";
		filter(segFile, locPvalFile, true, log);
		filter(segFile, locPvalFile, false, log);

	}
}
