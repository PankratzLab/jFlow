package org.genvisis.one.JL;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Map;

import org.genvisis.cnv.annotation.markers.MarkerBlastAnnotation;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ArrayUtils;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Not efficient.
 *
 */
public class BlastResults {

	private static void getEm(Project proj) {
		Map<String, MarkerBlastAnnotation> blastResults = MarkerBlastAnnotation.initForMarkers(proj.getMarkerNames());

		String outDir = proj.PROJECT_DIRECTORY.getValue() + "blastSummary/";
		new File(outDir).mkdirs();
		String outSum = outDir + "blastSummary.txt";
		try {

			PrintWriter writer = Files.openAppropriateWriter(outSum);

			writer.println("MarkerName\tPerfectMatch\tNonPerfectOnTarget\tOffTarget\tTotalAlignments");
			VCFFileReader reader = new VCFFileReader(new File(proj.BLAST_ANNOTATION_FILENAME.getValue()),
																							 true);
			int index = 0;
			for (VariantContext vc : reader) {
				// proj.getLog().reportTimeInfo(vc.getID());
				MarkerBlastAnnotation blastResult = blastResults.get(vc.getID());
				if (blastResult != null) {
					blastResult.parseAnnotation(vc, proj.getLog());
					boolean pm = blastResult.hasPerfectMatch(proj.getLog());
					int numOff = blastResult.getNumOffTarget(proj.getLog());
					int numOnNonPerf = blastResult.getNumOnTargetNonPerfect(proj.getLog());
					int numTotal = ArrayUtils.sum(blastResult.getAlignmentHistogram(proj));
					writer.println(vc.getID() + "\t" + pm + "\t" + numOnNonPerf + "\t" + numOff + "\t"
												 + numTotal);
				} else {
					proj.getLog().reportError("Cannot find marker " + vc.getID() + " in project");
				}
				index++;
			}
			reader.close();

			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + outSum);
			proj.getLog().reportException(e);
		}

	}

	public static void main(String[] args) {
		Project proj = new Project("C:/workspace/Genvisis/projects/singapore_PCs.properties", false);
		getEm(proj);
	}

}
