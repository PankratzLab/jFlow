package cnv.annotation;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Hashtable;

import common.Array;
import common.ArraySpecialList.ArrayBlastAnnotationList;
import common.ext;
import cnv.annotation.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import filesys.Segment;

/**
 * @author lane0212 Loads summarized blast results in annotation format
 */
public class BlastAnnotationLoader extends AnnotationFileLoader {
	private byte[] chrs;
	private int[] pos;
	private Hashtable<String, Integer> markerIndices;

	/**
	 * @param proj
	 * @param annotationFilename
	 *            the file created by a {@link BlastAnnotationWriter} run
	 * @param indexRequired
	 *            , should always be set to true
	 */
	public BlastAnnotationLoader(Project proj, String annotationFilename, boolean indexRequired) {
		super(proj, BlastAnnotationTypes.getBaseAnnotations(), annotationFilename, indexRequired);
		MarkerSet markerSet = proj.getMarkerSet();
		this.chrs = markerSet.getChrs();
		this.pos = markerSet.getPositions();
		this.markerIndices = proj.getMarkerIndices();
	}

	public MarkerBlastResult[] loadBlastAnnotationsFor(String[] markers) {
		MarkerBlastResult[] blastAnnotations = new MarkerBlastResult[markers.length];
		for (int i = 0; i < blastAnnotations.length; i++) {
			blastAnnotations[i] = new MarkerBlastResult(BLAST_ANNOTATION_TYPES.values());
		}

		Segment[] segs = getSegmentsForMarkers(markers);
		
		AnnotationQuery annotationQuery = getAnnotationQuery(segs);
		
		boolean[] found = Array.booleanArray(markers.length, false);
		
		while (annotationQuery.hasNext()) {
			VariantContext vc = annotationQuery.next();
			String id = vc.getID();
			int annoIndex = ext.indexOfStr(id, markers);
			
			if (annoIndex < 0) {
				proj.getLog().reportTimeWarning("Query has returned un-desired marker " + id + ", ignoring");
			} else {
				found[annoIndex] = true;
				for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {//each annotation type has a separate key in the file
					
					String[] groups = vc.getCommonInfo().getAttributeAsString(BLAST_ANNOTATION_TYPES.values()[i].toString(), ".").split(BLAST_ANNOTATION_TYPES.values()[i].getSep());
					for (int j = 0; j < groups.length; j++) {
						String[] segCigar = groups[j].split("/");
						blastAnnotations[annoIndex].getAnnotationLists()[i].add(new BlastAnnotation(TextCigarCodec.getSingleton().decode(segCigar[0]), new Segment(segCigar[1])));
					}
				}
			}

		}
		if (Array.booleanArraySum(found) != markers.length) {
			String error = markers.length + " markers were expected to be loaded, but only " + Array.booleanArraySum(found) + " markers were found";
			proj.getLog().reportTimeError(error);
			throw new IllegalStateException(error);
		} else {
			proj.getLog().reportTimeInfo("Loaded " + markers.length + " marker annotations");
		}
		return blastAnnotations;
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

	public static class MarkerBlastResult {
		private BLAST_ANNOTATION_TYPES[] bTypes;
		private ArrayBlastAnnotationList[] annotationLists;

		public MarkerBlastResult(BLAST_ANNOTATION_TYPES[] bTypes) {
			super();
			this.bTypes = bTypes;
			this.annotationLists = new ArrayBlastAnnotationList[bTypes.length];
			for (int i = 0; i < annotationLists.length; i++) {
				annotationLists[i] = new ArrayBlastAnnotationList(100);
			}
		}

		public BLAST_ANNOTATION_TYPES[] getbTypes() {
			return bTypes;
		}

		public ArrayBlastAnnotationList[] getAnnotationLists() {
			return annotationLists;
		}

	}

}
