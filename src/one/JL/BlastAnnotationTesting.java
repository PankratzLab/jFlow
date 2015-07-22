package one.JL;

import htsjdk.variant.vcf.VCFHeaderLineType;

import java.io.File;
import java.util.ArrayList;

import cnv.annotation.Annotation;
import cnv.annotation.AnnotationData;
import cnv.annotation.AnnotationFileWriter;
import cnv.annotation.BlastAnnotationLoader;
import cnv.annotation.BlastAnnotationTypes;
import cnv.annotation.BlastAnnotationWriter;
import cnv.annotation.LocusAnnotation;
import cnv.annotation.BlastAnnotationLoader.MarkerBlastResult;
import cnv.annotation.LocusAnnotation.Builder;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import common.Array;
import common.Files;
import common.ext;
import filesys.Segment;

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

		proj.getLog().reportTimeInfo("Loading " + t.size() + " markers");
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
		int maxGaps = 10;
		int maxMismatches = 10;

		// Writing:
		BlastAnnotationWriter blastAnnotation = new BlastAnnotationWriter(proj, annoFile, blastResultFiles, minAlignmentLength, maxGaps, maxMismatches, 15);
		blastAnnotation.summarizeResultFiles();
		blastAnnotation.close();

		blastAnnotation = new BlastAnnotationWriter(proj, annoFile, blastResultFiles, minAlignmentLength, maxGaps, maxMismatches, 15);
		blastAnnotation.summarizeResultFiles();
		blastAnnotation.close();
		// (Project proj, Annotation[] annotations, String annotationFilename, boolean overWriteExisting)
		AnnotationFileWriter test = new AnnotationFileWriter(proj, new AnnotationData[] { new AnnotationData(VCFHeaderLineType.String, "TestAdd", "TestAddidtion", "DSF", ".") }, annoFile, false) {
		};
		LocusAnnotation[] testAdd = getTestAddition(proj);
		for (int i = 0; i < testAdd.length; i++) {
			test.write(testAdd[i]);
		}
		test.close();
	}

	/**
	 * @param proj
	 * @return initialized blast summaries for all markers
	 */
	private static LocusAnnotation[] getTestAddition(Project proj) {
		// ABLookup abLookup = null;
		// if (Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
		// proj.getLog().reportTimeInfo("Ref and alt alleles will be determined by " + proj.AB_LOOKUP_FILENAME.getValue());
		// abLookup = new ABLookup(proj.getMarkerNames(), proj.AB_LOOKUP_FILENAME.getValue(), true, true, proj.getLog());
		// } else {
		// proj.getLog().reportTimeWarning(proj.AB_LOOKUP_FILENAME.getValue() + " did not exist so ref and alt alleles will be in-accurate");
		// }
		MarkerSet markerSet = proj.getMarkerSet();
		byte[] chrs = markerSet.getChrs();
		int[] pos = markerSet.getPositions();
		String[] markerNames = proj.getMarkerNames();
		LocusAnnotation[] anDatas = new LocusAnnotation[markerNames.length];
		for (int i = 0; i < anDatas.length; i++) {
			Builder builder = new Builder();
			// if (abLookup != null) {
			// String ref = abLookup.getLookup()[i][0] + "";
			// String alt = abLookup.getLookup()[i][1] + "";
			// if (!alt.equals("B") && !alt.equals("I") && !alt.equals("D")) {
			// builder.ref(ref);
			// builder.alts(new String[] { alt });
			// }
			// }
			builder.annotations(new AnnotationData[] { new AnnotationData(VCFHeaderLineType.String, "TestAdd", "TestAddidtion", "DSF_"+markerNames[i], ".") });
			Segment markerSeg = new Segment(chrs[i], pos[i], pos[i]);
			anDatas[i] = builder.build(markerNames[i], markerSeg);
		}
		return anDatas;
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
		test(proj, annoFile);
		testLoad(proj, annoFile);
	}
}
