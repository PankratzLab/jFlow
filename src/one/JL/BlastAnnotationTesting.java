package one.JL;

import java.io.File;
import java.util.ArrayList;

import cnv.annotation.BlastAnnotationLoader;
import cnv.annotation.BlastAnnotationWriter;
import cnv.annotation.BlastAnnotationLoader.MarkerBlastResult;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import common.Array;
import common.Files;
import common.ext;

/**
 * @author lane0212 Tests for the annotation writing/loading
 */
public class BlastAnnotationTesting {

	public static void testLoad(Project proj, String annoFile) {

		// Reading
		System.out.println("getting testers");
		ArrayList<String> t = getTestMarks(proj);
		String[] markers = proj.getMarkerNames();
		for (int i = 300; i < 4000; i++) {
			t.add(markers[i]);
		}
		t.add(markers[20]);
		t.add(markers[200000]);
		t.add(markers[100000]);
		t.add(markers[50000]);
		t.add(markers[200001]);
		t.add(markers[2]);
		t.add(markers[200003]);


		System.out.println("Finished getting testers");

		long time = System.currentTimeMillis();

		proj.getLog().reportTimeInfo("Loading "+t.size()+" markers");
		BlastAnnotationLoader blastAnnotationLoader = new BlastAnnotationLoader(proj, annoFile, true);

		MarkerBlastResult[] markerBlastResults = blastAnnotationLoader.loadBlastAnnotationsFor(Array.toStringArray(t));
		proj.getLog().reportTimeElapsed(time);
		for (int i = 0; i < markerBlastResults.length; i++) {
			for (int j = 0; j < markerBlastResults[i].getAnnotationLists().length; j++) {
				for (int j2 = 0; j2 < markerBlastResults[i].getAnnotationLists()[j].size(); j2++) {
					// System.out.println();
					// System.out.println(markerBlastResults[i].getbTypes()[j].getName() + "\t" + t.get(i) + "\t" + markerBlastResults[i].getAnnotationLists()[j].get(j2).getRefLoc().getUCSClocation() + "\t" + markerBlastResults[i].getAnnotationLists()[j].get(j2).getCigar().toString());
				}
			}
		}
		blastAnnotationLoader.close();
	}

	public static void test(Project proj, String annoFile) {

		new File(ext.parseDirectoryOfFile(annoFile)).mkdirs();
		String[] blastResultFiles = Files.list("/home/pankrat2/shared/aric_exome_chip/Blasts/", "GPL18544_humanexome-12v1_a.csv.blasted.ws.30.rep.0.tmp", null, true, false, true);
		int minAlignmentLength = proj.getArrayType().getProbeLength() - 10;
		int maxGaps = 0;
		int maxMismatches = 0;

		// Writing:
		BlastAnnotationWriter blastAnnotation = new BlastAnnotationWriter(proj, annoFile, blastResultFiles, minAlignmentLength, maxGaps, maxMismatches, 15);
		blastAnnotation.summarizeResultFiles();
		blastAnnotation.close();

	}

	private static ArrayList<String> getTestMarks(Project proj) {
		ArrayList<String> t = new ArrayList<String>();
		MarkerSet markerSet = proj.getMarkerSet();

		for (int i = 0; i < markerSet.getIndicesByChr().length; i++) {
			if (markerSet.getIndicesByChr()[i].length - 1 >= 0) {
				t.add(proj.getMarkerNames()[markerSet.getIndicesByChr()[i][markerSet.getIndicesByChr()[i].length - 1]]);
			}
		}
		return t;
	}

	public static void main(String[] args) {
		Project proj = new Project("/home/pankrat2/lanej/projects/aric_exome.properties", false);
		String annoFile = proj.PROJECT_DIRECTORY.getValue() + "TestBlastLoad/blast.anno.vcf.gz";
	//	test(proj,annoFile);
		testLoad(proj, annoFile);
	}
}
