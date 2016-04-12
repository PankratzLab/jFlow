package one.JL;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import seq.manage.BamImport;
import seq.manage.ReferenceGenome;
import seq.manage.BamImport.NGS_MARKER_TYPE;
import stats.Stats;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;
import cnv.filesys.MarkerSet.PreparedMarkerSet;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.hmm.PennHmm.ViterbiResult;
import cnv.manage.Markers;
import cnv.qc.GcAdjustor.GcModel;
import cnv.var.CNVariant;
import cnv.var.CNVariant.CNVBuilder;
import cnv.var.LocusSet;
import cnv.var.LocusSet.TO_STRING_TYPE;
import filesys.Segment;

public class Scratch {

	private static void bin(int numStates) {
		Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);
		GcModel gcModel = GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(), true, proj.getLog());

		double[] gcs = gcModel.getGcs();
		OpdfGaussian opdfGaussian = new OpdfGaussian(Array.mean(gcs, true), Math.pow(Array.stdev(gcs, true), 2));
		int zeroState = 5;
		int[] stateSequence = new int[gcs.length];
		for (int i = 0; i < gcs.length; i++) {
			double cdf = opdfGaussian.cdf(new ObservationReal(gcs[i]));
			// Stats.ttestOneSample(mean, stdev, n, expected)
			int state = (int) Math.round((double) cdf * numStates) - 1;
			stateSequence[i] = state + zeroState;
			// int state =
		}
		PreparedMarkerSet markerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());
		int[][] indices = markerSet.getIndicesByChr();
		ArrayList<CNVariant> allGc = new ArrayList<CNVariant>();

		for (int i = 0; i < indices.length; i++) {
			if (indices[i].length > 0) {
				int[] currentStates = Array.subArray(stateSequence, indices[i]);
				int[] positions = Array.subArray(markerSet.getPositions(), indices[i]);
				byte chr = (byte) i;
				ViterbiResult vtr = new ViterbiResult(currentStates, null);
				LocusSet<CNVariant> tmp = vtr.analyzeStateSequence(proj, "7139072086_R01C01", "7139072086_R01C01", chr, positions, null, 0, false, true);
				ArrayList<int[]> indicesStates = vtr.getIndexStateChange();
				for (int j = 0; j < tmp.getLoci().length; j++) {
					CNVBuilder builder = new CNVBuilder(tmp.getLoci()[j]);
					int[] startStop = indicesStates.get(j);
					int[] stateIndices = Array.subArray(indices[i], startStop[0], startStop[1] + 1);
					double[] clean = Array.removeNaN(Array.subArray(gcs, stateIndices));
					builder.score(Array.median(clean));
					if (j % 2 == 0) {
						builder.cn(0);
					} else {
						builder.cn(4);
					}
					allGc.add(builder.build());
				}
			}
		}

		String outDir = proj.PROJECT_DIRECTORY.getValue() + "GC_Bins/";
		new File(outDir).mkdirs();
		String outFile = outDir + "gcBins_" + numStates + ".bins";
		LocusSet<CNVariant> binset = new LocusSet<CNVariant>(allGc.toArray(new CNVariant[allGc.size()]), true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		binset.writeRegions(outFile, TO_STRING_TYPE.REGULAR, true, proj.getLog());
	}

	private static String generateGCModel(Project proj, ReferenceGenome referenceGenome, int buffer) {
		String gcFile = ext.addToRoot(proj.GC_MODEL_FILENAME.getValue(), ".buffer_" + buffer);
		if (!Files.exists(gcFile)) {
			MarkerSet markerSet = proj.getMarkerSet();
			String[] markerNames = markerSet.getMarkerNames();

			try {
				PrintWriter writer = new PrintWriter(new FileWriter(gcFile));
				String[] header = new String[] { "Name", "Chr", "Position", "GC" };
				writer.println(Array.toStr(header));
				for (int i = 0; i < markerNames.length; i++) {
					if (i % 10000 == 0) {
						proj.getLog().reportTimeInfo("Loaded gc content for " + (i + 1) + " bins");
					}
					writer.println(markerNames[i] + "\t" + markerSet.getChrs()[i] + "\t" + markerSet.getPositions()[i] + "\t" + ReferenceGenome.getPercentGC(referenceGenome.getSequenceFor(new Segment(markerNames[i].split("\\|")[0]).getBufferedSegment(buffer), buffer > 100000)));
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + gcFile);
				proj.getLog().reportException(e);
			}
		} else {
			proj.getLog().reportFileExists(gcFile);
		}
		return gcFile;
	}

	public static void main(String[] args) {
		String dir = "C:/data/ARIC/scratch/";
		Logger log = new Logger(dir + "log");
		String[] pcMarks = HashVec.loadFileToStringArray(dir +"autosomal_PC_markers.oneHitWonders_20.txt", false, new int[] { 0 }, false);
		String[] mito = HashVec.loadFileToStringArray(dir + "gw6_MT_USE.oneHitWonders_20.txt" , false, new int[] { 0 }, false);

		String[] combo = Array.concatAll(pcMarks, mito);
		log.reportTimeInfo("PCMarks " + pcMarks.length);
		log.reportTimeInfo("Mito" + mito.length);

		String pos = dir + "markerPositionsHG18.txt";
		String mset = ext.rootOf(pos, false) + ".ser";
		Markers.orderMarkers(null, pos, mset, log);
		MarkerSet markerSet = MarkerSet.load(mset, false);

		int[] pcmarks = ext.indexLargeFactors(pcMarks, markerSet.getMarkerNames(), true, log, true, true);
		int[][] indices = markerSet.getIndicesByChr();
		ArrayList<String> outMarks = new ArrayList<String>();
		outMarks.add("SNP Name\tChr\tPosition");
		for (int i = 0; i < pcmarks.length; i++) {
			outMarks.add(markerSet.getMarkerNames()[pcmarks[i]] + "\t" + markerSet.getChrs()[pcmarks[i]] + "\t" + markerSet.getPositions()[pcmarks[i]]);
		}
		for (int i = 0; i < indices[26].length; i++) {
			outMarks.add(markerSet.getMarkerNames()[indices[26][i]] + "\t" + markerSet.getChrs()[indices[26][i]] + "\t" + markerSet.getPositions()[indices[26][i]]);
		}

		String out = ext.addToRoot(pos, ".PC.mito.subset");
		Files.writeArrayList(outMarks, out);
		Markers.orderMarkers(combo, out, out + "check", log);

		// String[] us = HashVec.loadFileToStringArray(dir + "custom.gcmodel", false, new int[] { 0 }, false);
		// String[] them = HashVec.loadFileToStringArray(dir + "SNPThem", false, new int[] { 0 }, false);
		//
//		String[] usgc = HashVec.loadFileToStringArray(dir + "custom.gcmodel", false, new int[] { 3 }, false);
//		String[] themgc = HashVec.loadFileToStringArray(dir + "SNPThem", false, new int[] { 3 }, false);
//
//		int[] indices = ext.indexLargeFactors(us, them, false, new Logger(), false, false);
//		ArrayList<String> go = new ArrayList<String>();
//		for (int i = 0; i < indices.length; i++) {
//			if (indices[i] >= 0) {
//				go.add(us[i] + "\t" + them[indices[i]] + "\t" + usgc[i] + "\t" + themgc[indices[i]]);
//			}
//		}
//		Files.writeArrayList(go, dir + "merge");

		// Project proj = new Project("C:/workspace/Genvisis/projects/Cushings.properties", false);
//		generateGCModel(proj, new ReferenceGenome(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(), proj.getLog()), 1000000);
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
		// // ColorManager<String> cm = ColorExt.getColorManager(proj, proj.MARKER_COLOR_KEY_FILENAMES.getValue()[0]);
		//
		// String[] names = proj.getMarkerNames();
		// ArrayList<String> out = new ArrayList<String>();
		// GcModel gcModel = GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(), true, proj.getLog());
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
