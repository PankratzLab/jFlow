package one.JL;

import java.util.ArrayList;

import cnv.annotation.BlastAnnotationLoader;
import cnv.annotation.BlastAnnotationWriter;
import cnv.annotation.BlastAnnotationLoader.MarkerBlastResult;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import common.Array;
import common.Files;

/**
 * @author lane0212 Tests for the annotation writing/loading
 */
public class BlastAnnotationTesting {

	public static void test() {
		Project proj = new Project("/home/pankrat2/lanej/projects/aric_exome.properties", false);
		String outfile = proj.PROJECT_DIRECTORY.getValue() + "Blasts/blast.anno.vcf.gz";
		String[] blastResultFiles = Files.list("/home/pankrat2/shared/aric_exome_chip/Blasts/", "GPL18544_humanexome-12v1_a.csv.blasted.ws.30.rep.0.tmp", null, true, false, true);
		int minAlignmentLength = proj.getArrayType().getProbeLength() - 10;
		int maxGaps = 0;
		int maxMismatches = 0;

		BlastAnnotationWriter blastAnnotation = new BlastAnnotationWriter(proj, outfile, blastResultFiles, minAlignmentLength, maxGaps, maxMismatches, 15);
		blastAnnotation.summarizeResultFiles();
		blastAnnotation.close();

		BlastAnnotationLoader blastAnnotationLoader = new BlastAnnotationLoader(proj, outfile, true);
		String[] testMarkers = Array.subArray(proj.getMarkerNames(), Array.intArray(100));
		ArrayList<String> testItOut = new ArrayList<String>();
		for (int i = 0; i < testMarkers.length; i++) {
			testItOut.add(testMarkers[i]);
		}
		MarkerSet markerSet = proj.getMarkerSet();
		for (int i = 0; i < markerSet.getPositionsByChr().length; i++) {
			if (markerSet.getPositionsByChr()[i].length - 1 >= 0) {
				testItOut.add(proj.getMarkerNames()[markerSet.getPositionsByChr()[i][markerSet.getPositionsByChr()[i].length - 1]]);
			}
		}
		MarkerBlastResult[] markerBlastResults = blastAnnotationLoader.loadBlastAnnotationsFor(Array.toStringArray(testItOut));
		for (int i = 0; i < markerBlastResults.length; i++) {
			for (int j = 0; j < markerBlastResults[i].getAnnotationLists().length; j++) {
				for (int j2 = 0; j2 < markerBlastResults[i].getAnnotationLists()[j].size(); j2++) {
					System.out.println(markerBlastResults[i].getbTypes()[j].getName());
					System.out.println(testItOut.get(i) + "\t" + markerBlastResults[i].getAnnotationLists()[j].get(j2).getCigar().toString());
				}
			}
		}
	}

	public static void main(String[] args) {
		test();
	}
}
