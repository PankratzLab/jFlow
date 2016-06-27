package cnv.qc;

import java.util.ArrayList;

import cnv.annotation.AnnotationParser;
import cnv.annotation.MarkerAnnotationLoader;
import cnv.annotation.MarkerBlastAnnotation;
import cnv.annotation.MarkerGCAnnotation;
import cnv.annotation.AnnotationFileLoader.QUERY_ORDER;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import common.Logger;

public class MarkerBlastQC {

	private static void getOneHitWonders(Project proj, String blastVCF, double crossHybePercent, Logger log) {
		MarkerSet markerSet = proj.getMarkerSet();
		String[] markerNames = markerSet.getMarkerNames();
		MarkerAnnotationLoader markerAnnotationLoader = new MarkerAnnotationLoader(proj, null, proj.BLAST_ANNOTATION_FILENAME.getValue(), proj.getMarkerSet(), true);
		markerAnnotationLoader.setReportEvery(500000);
		MarkerGCAnnotation[] gcAnnotations = MarkerGCAnnotation.initForMarkers(proj, markerNames, markerAnnotationLoader.getMarkerSet(), markerAnnotationLoader.getIndices());
		MarkerBlastAnnotation[] blastResults = MarkerBlastAnnotation.initForMarkers(markerNames);
		ArrayList<AnnotationParser[]> parsers = new ArrayList<AnnotationParser[]>();
		parsers.add(gcAnnotations);
		parsers.add(blastResults);

		markerAnnotationLoader.fillAnnotations(null, parsers, QUERY_ORDER.ONE_PER_IN_ORDER);

		for (int i = 0; i < blastResults.length; i++) {
			MarkerBlastAnnotation current = blastResults[i];
			
			int[] alignmentHistogram = current.getAlignmentHistogram(proj);
			int sub =(int) Math.round((double)crossHybePercent*alignmentHistogram.length);
			if(current.hasPerfectMatch(log)){
				
			}
		}

	}

}
