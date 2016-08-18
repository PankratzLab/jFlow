package org.genvisis.one.JL;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.genvisis.cnv.annotation.MarkerBlastAnnotation;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Not efficient.
 *
 */
public class BlastResults {

	private static void getEm(Project proj) {
		MarkerBlastAnnotation[] blastResults = MarkerBlastAnnotation
				.initForMarkers(proj.getMarkerNames());

		String outDir = proj.PROJECT_DIRECTORY.getValue() + "blastSummary/";
		new File(outDir).mkdirs();
		String outSum = outDir + "blastSummary.txt";
		try {

			PrintWriter writer = new PrintWriter(new FileWriter(outSum));

			writer.println("MarkerName\tPerfectMatch\tNonPerfectOnTarget\tOffTarget\tTotalAlignments");
			VCFFileReader reader = new VCFFileReader(new File(proj.BLAST_ANNOTATION_FILENAME.getValue()), true);
			int index = 0;
			for (VariantContext vc : reader) {
				//proj.getLog().reportTimeInfo(vc.getID());
				blastResults[index].parseAnnotation(vc, proj.getLog());
				boolean pm = blastResults[index].hasPerfectMatch(proj.getLog());
				int numOff = blastResults[index].getNumOffTarget(proj.getLog());
				int numOnNonPerf = blastResults[index].getNumOnTargetNonPerfect(proj.getLog());
				int numTotal = Array.sum(blastResults[index].getAlignmentHistogram(proj));
				writer.println(vc.getID() + "\t" + pm + "\t" + numOnNonPerf + "\t" + numOff+"\t"+numTotal);
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
