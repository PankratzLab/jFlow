package one.JL;

import cnv.annotation.MarkerBlastAnnotation;
import cnv.filesys.Project;

/**
 * Not efficient.
 *
 */
public class BlastResults {

	private static void getEm(Project proj) {
		MarkerBlastAnnotation[] blastResults = MarkerBlastAnnotation
				.initForMarkers(proj.getMarkerNames());

	}

}
