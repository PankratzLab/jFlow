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
	private MarkerSet markerSet;
	private Hashtable<String, Integer> indices;
	private byte[] chrs;
	private int[] pos;
	private Hashtable<String, Integer> markerIndices;

	/**
	 * @param proj
	 * @param annotationFilename
	 * @param markerSet
	 *            will be loaded if null, can save time if it is already available
	 * @param indexRequired
	 *            should always be true for now
	 */
	public MarkerAnnotationLoader(Project proj, String annotationFilename, MarkerSet markerSet, boolean indexRequired) {
		super(proj, null, annotationFilename, indexRequired);
		if (markerSet == null) {
			markerSet = proj.getMarkerSet();
		}
		this.indices = proj.getMarkerIndices();
		this.markerSet = markerSet;
		this.chrs = markerSet.getChrs();
		this.pos = markerSet.getPositions();
		this.markerIndices = proj.getMarkerIndices();
	}

	public MarkerSet getMarkerSet() {
		return markerSet;
	}

	public Hashtable<String, Integer> getIndices() {
		return indices;
	}

	/**
	 * @param markers
	 * @param parsersQueries
	 *            typically each entry in the {@link AnnotationParser} array represents a single marker
	 */
	public void fillAnnotations(final String[] markers, List<AnnotationParser[]> parsersQueries) {
		Segment[] markerSegments = getSegmentsForMarkers(markers);
		query(markerSegments, parsersQueries);
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
