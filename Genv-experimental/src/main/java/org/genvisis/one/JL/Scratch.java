package org.genvisis.one.JL;

import java.util.ArrayList;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class Scratch {

	public static void main(String[] args) {
		String dir = "C:/data/ARIC/scratch/";
		Logger log = new Logger(dir + "log");
		String[] pcMarks = HashVec.loadFileToStringArray(dir
																												+ "autosomal_PC_markers.oneHitWonders_20.txt",
																											false, new int[] {0}, false);
		String[] mito = HashVec.loadFileToStringArray(dir	+ "gw6_MT_USE.oneHitWonders_20.txt", false,
																									new int[] {0}, false);

		String[] combo = Array.concatAll(pcMarks, mito);
		log.reportTimeInfo("PCMarks " + pcMarks.length);
		log.reportTimeInfo("Mito" + mito.length);

		String pos = dir + "markerPositionsHG18.txt";
		String mset = ext.rootOf(pos, false) + ".ser";
		Markers.orderMarkers(null, pos, mset, log);
		MarkerSet markerSet = MarkerSet.load(mset, false);

		int[] pcmarks =
									ext.indexLargeFactors(pcMarks, markerSet.getMarkerNames(), true, log, true, true);
		int[][] indices = markerSet.getIndicesByChr();
		ArrayList<String> outMarks = new ArrayList<String>();
		outMarks.add("SNP Name\tChr\tPosition");
		for (int pcmark : pcmarks) {
			outMarks.add(markerSet.getMarkerNames()[pcmark]	+ "\t" + markerSet.getChrs()[pcmark] + "\t"
										+ markerSet.getPositions()[pcmark]);
		}
		for (int i = 0; i < indices[26].length; i++) {
			outMarks.add(markerSet.getMarkerNames()[indices[26][i]]	+ "\t"
										+ markerSet.getChrs()[indices[26][i]] + "\t"
										+ markerSet.getPositions()[indices[26][i]]);
		}

		String out = ext.addToRoot(pos, ".PC.mito.subset");
		Files.writeIterable(outMarks, out);
		Markers.orderMarkers(combo, out, out + "check", log);

		// String[] us = HashVec.loadFileToStringArray(dir + "custom.gcmodel", false, new int[] { 0 },
		// false);
		// String[] them = HashVec.loadFileToStringArray(dir + "SNPThem", false, new int[] { 0 },
		// false);
		//
		// String[] usgc = HashVec.loadFileToStringArray(dir + "custom.gcmodel", false, new int[] { 3 },
		// false);
		// String[] themgc = HashVec.loadFileToStringArray(dir + "SNPThem", false, new int[] { 3 },
		// false);
		//
		// int[] indices = ext.indexLargeFactors(us, them, false, new Logger(), false, false);
		// ArrayList<String> go = new ArrayList<String>();
		// for (int i = 0; i < indices.length; i++) {
		// if (indices[i] >= 0) {
		// go.add(us[i] + "\t" + them[indices[i]] + "\t" + usgc[i] + "\t" + themgc[indices[i]]);
		// }
		// }
		// Files.writeArrayList(go, dir + "merge");

		// Project proj = new Project("C:/workspace/Genvisis/projects/Cushings.properties", false);
		// generateGCModel(proj, new ReferenceGenome(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(),
		// proj.getLog()), 1000000);
		//
		//
		// Segment test = new Segment("chr1:874792-874792");
		// System.out.println(test.overlaps(test));
		// System.exit(1);
		//
		// int[] numStates = new int[] { 5, 10, 15, 20, 50, 100 };
		//
		// for (int i = 0; i < numStates.length; i++) {
		// // bin(numStates[i]);
		// }
		//
		// // System.out.println(new Segment("chr1:883539-883539").getSize());
		// Project proj = new Project("C:/workspace/Genvisis/projects/Cushings.properties", false);
		// String colDir = proj.PROJECT_DIRECTORY.getValue() + "Colors/";
		// String colFile = colDir + "markerColors.txt";
		// proj.MARKER_COLOR_KEY_FILENAMES.setValue(new String[] { colFile });
		// proj.saveProperties();
		// // ColorManager<String> cm = ColorExt.getColorManager(proj,
		// proj.MARKER_COLOR_KEY_FILENAMES.getValue()[0]);
		//
		// String[] names = proj.getMarkerNames();
		// ArrayList<String> out = new ArrayList<String>();
		// GcModel gcModel = GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(), true,
		// proj.getLog());
		// System.out.println(Array.mean(gcModel.getGcs()));
		// // Color[] color = new Co
		//
		// out.add("MarkerName\tCLASS=MARKER_COLOR;OFF_TARGET=Blue;VARIANT_SITE=RED;ON_TARGET=Green");
		// // for (int i = 0; i < names.length; i++) {
		// // String markClass = "ON_TARGET";
		// // if (names[i].contains(NGS_MARKER_TYPE.OFF_TARGET.getFlag())) {
		// // markClass = NGS_MARKER_TYPE.OFF_TARGET.getFlag();
		// // } else if (names[i].contains(BamImport.VARIANT_SITE_FLAG)) {
		// // markClass = BamImport.OFF_TARGET_FLAG;
		// // }
		// // out.add(names[i] + "\t" + markClass);
		// // }
		// Files.writeList(Array.toStringArray(out), colFile);
	}

}
