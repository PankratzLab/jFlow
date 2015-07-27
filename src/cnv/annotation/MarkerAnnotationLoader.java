package cnv.annotation;

import java.util.Hashtable;
import java.util.List;

import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import filesys.Segment;

/**
 * @author lane0212 Class that concentrates on loading annotations for specific markers
 */
public class MarkerAnnotationLoader extends AnnotationFileLoader {
	private byte[] chrs;
	private int[] pos;
	private Hashtable<String, Integer> markerIndices;

	public MarkerAnnotationLoader(Project proj, String annotationFilename, MarkerSet markerSet, boolean indexRequired) {
		super(proj, null, annotationFilename, indexRequired);
		if (markerSet == null) {
			markerSet = proj.getMarkerSet();
		}
		this.chrs = markerSet.getChrs();
		this.pos = markerSet.getPositions();
		this.markerIndices = proj.getMarkerIndices();
	}

	public void fillAnnotations(final String[] markers, List<AnnotationParser[]> parsers) {
		Segment[] markerSegments = getSegmentsForMarkers(markers);
		query(markerSegments, parsers);
	}

	private Segment[] getSegmentsForMarkers(final String[] markers) {
		Segment[] segs = new Segment[markers.length];
		for (int i = 0; i < segs.length; i++) {
			int markerIndex = markerIndices.get(markers[i]);
			Segment markerSeg = new Segment(chrs[markerIndex], pos[markerIndex], pos[markerIndex]);
			segs[i] = markerSeg;
		}
		return segs;
	}

}
