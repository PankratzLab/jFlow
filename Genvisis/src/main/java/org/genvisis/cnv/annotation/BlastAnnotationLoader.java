package org.genvisis.cnv.annotation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;

import org.genvisis.cnv.annotation.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import org.genvisis.cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.ArraySpecialList.ArrayBlastAnnotationList;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * @author lane0212 Loads summarized blast results in annotation format
 */
public class BlastAnnotationLoader extends AnnotationFileLoader {
	private final byte[] chrs;
	private final int[] pos;
	private final Hashtable<String, Integer> markerIndices;

	/**
	 * @param proj
	 * @param annotationFilename the file created by a {@link BlastAnnotationWriter} run
	 * @param indexRequired , should always be set to true
	 */
	public BlastAnnotationLoader(Project proj, String annotationFilename, boolean indexRequired) {
		super(proj, null, BlastAnnotationTypes.getBaseAnnotations(), annotationFilename, indexRequired);
		MarkerSet markerSet = proj.getMarkerSet();
		chrs = markerSet.getChrs();
		pos = markerSet.getPositions();
		markerIndices = proj.getMarkerIndices();
	}

	/**
	 * @param markers
	 * @param otherQueries these queries can be null, but if not the appropriate annotations will be
	 *        parsed
	 * @return
	 */
	public MarkerBlastResult[] loadBlastAnnotationsFor(	String[] markers,
																											AnnotationParser[]... otherQueries) {

		if (Array.unique(markers).length != markers.length) {
			String error = "Internal error, markers for blast annotation retrieval must be unique";
			proj.getLog().reportTimeError(error);
			throw new IllegalArgumentException(error);
		}

		MarkerBlastResult[] blastAnnotations = initResults(markers);

		Segment[] segs = getSegmentsForMarkers(markers);

		AnnotationQuery annotationQuery = getAnnotationQuery(segs);

		boolean[] found = Array.booleanArray(markers.length, false);

		int count = 0;
		while (annotationQuery.hasNext()) {
			count++;
			if (count % proj.MAX_MARKERS_LOADED_PER_CYCLE.getValue() == 0) {
				proj.getLog().reportTimeInfo("Loaded " + count + " annotations");
			}
			VariantContext vc = annotationQuery.next();
			if (otherQueries != null) {
				for (AnnotationParser[] annotationParsers : otherQueries) {
					for (AnnotationParser annotationParser : annotationParsers) {
						if (annotationParser.shouldAnnotateWith(vc, proj.getLog())) {
							annotationParser.parseAnnotation(vc, proj.getLog());
						}
					}
				}
			}
			String id = vc.getID();
			int annoIndex = ext.indexOfStr(id, markers);
			if (annoIndex < 0) {
				proj.getLog()
						.reportTimeWarning("Query has returned un-desired marker " + id + ", ignoring");
			} else {
				found[annoIndex] = true;
				blastAnnotations[annoIndex].parseAnnotation(vc, proj.getLog());
			}
		}
		if (Array.booleanArraySum(found) != markers.length) {
			String error = markers.length	+ " markers were expected to be loaded, but only "
											+ Array.booleanArraySum(found) + " markers were found";
			for (int i = 0; i < found.length; i++) {
				if (!found[i]) {
					error += "\nMissing " + markers[i];
				}
			}
			proj.getLog().reportTimeError(error);
			proj.getLog().reportTimeError(Array.toStr(markers));
			// throw new IllegalStateException(error);
		} else {
			proj.getLog().reportTimeInfo("Loaded " + markers.length + " marker annotations");
		}
		return blastAnnotations;
	}

	private MarkerBlastResult[] initResults(String[] markers) {
		MarkerBlastResult[] blastAnnotations = new MarkerBlastResult[markers.length];
		for (int i = 0; i < blastAnnotations.length; i++) {
			blastAnnotations[i] = new MarkerBlastResult(markers[i], BLAST_ANNOTATION_TYPES.values(), 100);
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

	/**
	 * @author lane0212 Probably could be its own class
	 */
	public static class MarkerBlastResult implements AnnotationParser {
		private final BLAST_ANNOTATION_TYPES[] bTypes;
		private final ArrayBlastAnnotationList[] annotationLists;
		private final String markerName;
		private boolean found;

		public MarkerBlastResult(	String markerName, BLAST_ANNOTATION_TYPES[] bTypes,
															int initialCapacity) {
			super();
			this.markerName = markerName;
			this.bTypes = bTypes;
			annotationLists = new ArrayBlastAnnotationList[bTypes.length];
			for (int i = 0; i < annotationLists.length; i++) {
				annotationLists[i] = new ArrayBlastAnnotationList(initialCapacity);
			}
		}

		public boolean hasPerfectMatch(Logger log) {
			return getAnnotationsFor(BLAST_ANNOTATION_TYPES.PERFECT_MATCH, log).size() > 0;
		}

		public int getNumOffTarget(Logger log) {
			return getAnnotationsFor(BLAST_ANNOTATION_TYPES.OFF_T_ALIGNMENTS, log).size();
		}

		public ArrayList<BlastAnnotation> getAnnotationsFor(BLAST_ANNOTATION_TYPES btype, Logger log) {
			int testIndex = getAnnotationIndexFor(btype, log);
			return annotationLists[testIndex];
		}

		public BLAST_ANNOTATION_TYPES[] getbTypes() {
			return bTypes;
		}

		public ArrayBlastAnnotationList[] getAnnotationLists() {
			return annotationLists;
		}

		private int getAnnotationIndexFor(BLAST_ANNOTATION_TYPES bType, Logger log) {
			int index = -1;
			for (int i = 0; i < bTypes.length; i++) {
				if (bTypes[i] == bType) {
					index = i;
					break;
				}
			}
			if (index < 0) {
				String error = "Internal error: Annotation does not contain " + bType;
				log.reportTimeError(error);
				throw new IllegalStateException(error);
			}
			return index;
		}

		@Override
		public void parseAnnotation(VariantContext vc, Logger log) {
			for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {// each annotation type has
																																				// a separate key in the
																																				// file
				String info = vc.getCommonInfo()
												.getAttributeAsString(BLAST_ANNOTATION_TYPES.values()[i].toString(),
																							BLAST_ANNOTATION_TYPES.values()[i].getDefaultValue());
				if (!info.equals(BLAST_ANNOTATION_TYPES.values()[i].getDefaultValue())) {
					List<String> groups = Arrays.asList(info.replaceAll("\\[", "").replaceAll("\\]", "")
																									.split("\\s*,\\s*"));
					for (String group : groups) {
						String[] segCigarStrand = group.split("/");
						annotationLists[i].add(new BlastAnnotation(	TextCigarCodec.decode(segCigarStrand[0]),
																												new Segment(segCigarStrand[1]),
																												Strand.valueOf(segCigarStrand[2]),
																												PROBE_TAG.valueOf(segCigarStrand[3]),
																												Double.parseDouble(segCigarStrand[4])));
					}
				}
			}
		}

		@Override
		public boolean shouldAnnotateWith(VariantContext vc, Logger log) {
			return markerName.equals(vc.getID());
		}

		@Override
		public void setFound(boolean found) {
			this.found = found;
		}

		@Override
		public boolean isFound() {
			return found;
		}

	}
}
