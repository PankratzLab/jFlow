package org.genvisis.cnv.analysis;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.GcAdjustor.GCAdjustorBuilder;
import org.genvisis.cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import org.genvisis.cnv.qc.GcAdjustor.GcModel;
import org.genvisis.cnv.qc.GcAdjustorParameter;
import org.genvisis.cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;

/**
 * Generate marker-dominant LRR SD values, potentially including any number of GC correction models.
 */
public final class MarkerStats {

	private MarkerStats() {
		// Prevent instantiation of utility class
	}

	private static void generateMarkerStats(String projectPath, String outName) {
		Project proj = new Project(projectPath, false);
		String outFile = proj.PROJECT_DIRECTORY.getValue() + outName;
		Hashtable<String, Integer> markerIndices = proj.getMarkerIndices();
		Logger log = proj.getLog();
		String gcModel = proj.GC_MODEL_FILENAME.getValue(false, false);
		String gcModelName = new File(gcModel).getName();

		List<String> outHeader = new ArrayList<String>();
		outHeader.add("Marker name");
		outHeader.add("Chr");
		outHeader.add("Pos");
		outHeader.add("GC: " + gcModelName);
		outHeader.add("LRR SD");
		outHeader.add("BAF Mean");
		outHeader.add("BAF SD");

		List<GcAdjustorParameters> params = getParams(proj, outHeader, gcModel);

		Hashtable<String, Vector<String>> gcModelHash = HashVec.loadFileToHashVec(gcModel, true, "Name",
		                                                                          "GC");

		PrintWriter writer = Files.getAppropriateWriter(outFile);
		writer.println(Array.toStr(outHeader, "\t"));

		MDL mdl = new MDL(proj, proj.getMarkerSet(), proj.getMarkerNames());

		try {
			while (mdl.hasNext()) {
				MarkerData marker = mdl.next();
				final List<String> line = new ArrayList<String>();
				String markerName = marker.getMarkerName();
				int markerIndexInProject = markerIndices.get(markerName);
				float bafAvg = Array.mean(marker.getBAFs(), true);
				float bafSd = Array.stdev(marker.getBAFs(), true);
				float lrrSd = Array.stdev(marker.getLRRs(), true);
				line.add(markerName);
				line.add(String.valueOf(marker.getChr()));
				line.add(String.valueOf(marker.getPosition()));
				line.add(gcModelHash.get(markerName).get(0));
				line.add(String.valueOf(bafAvg));
				line.add(String.valueOf(bafSd));
				line.add(String.valueOf(lrrSd));
				for (GcAdjustorParameters ps : params) {
					float[][] lrrbaf = marker.getGCCorrectedLRRBAF(ps, markerIndexInProject, log);
					float gcLrrSd = Array.stdev(lrrbaf[1], true);
					line.add(String.valueOf(gcLrrSd));
					// NB: looks like GC correction doesn't currently affect BAF calculation. Not clear why both come back?
				}

				writer.println(Array.toStr(line, "\t"));
			}
		} finally {
			mdl.shutdown();
			writer.close();
		}

		log.report("Finished writing output to: " + outFile);
	}

	private static List<GcAdjustorParameters> getParams(Project proj, List<String> header,
	                                                    String gcModel) {
		List<GcAdjustorParameters> list = new ArrayList<GcAdjustorParameters>();
		final Logger log = proj.getLog();
		String[] gcParamsFile = proj.GC_CORRECTION_PARAMETERS_FILENAMES.getValue();

		if (gcParamsFile == null || gcParamsFile.length == 0) {
			log.reportTimeWarning("No GC-corrected parameter files found. Will attempt to create default for GC model file: "
			                      + gcModel);
			GCAdjustorBuilder builder = new GcAdjustor.GCAdjustorBuilder();
			// TODO return value is a list of outputs which I guess we can deserialize?
			// GCAdjustorBuilders
			GcModel model = GcModel.populateFromFile(gcModel, false, log);
			String cents = proj.CUSTOM_CENTROIDS_FILENAME.getValue();
			String[] centParams = Files.exists(cents) ? new String[] {cents} : null;
			gcParamsFile = GcAdjustorParameter.generateAdjustmentParameters(proj, new GCAdjustorBuilder[] {builder},
			                                                 centParams,
			                                                 new GC_CORRECTION_METHOD[] {GC_CORRECTION_METHOD.GENVISIS_GC},
			                                                 model, new String[] {"GC_Correction" + File.separator},
			                                                 proj.NUM_THREADS.getValue(), false)[0][0];
		}
		for (String gaf : gcParamsFile) {
			String name = new File(gaf).getName();
			header.add(name + ": LRR SD");
			list.add(GcAdjustorParameters.readSerial(gaf, log));
		}
		return list;
	}

	public static void main(String... args) {
		CLI c = new CLI(MarkerStats.class);
		final String outFile = "marker_lrr_sd.xln";

		c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ, true);
		c.addArgWithDefault(CLI.ARG_OUTFILE, CLI.DESC_OUTFILE, outFile);
		c.parseWithExit(args);

		generateMarkerStats(c.get(CLI.ARG_PROJ), c.get(CLI.ARG_OUTFILE));
	}
}
