/**
 * 
 */
package org.genvisis.cnv.annotation.segments;

import java.util.ArrayList;
import java.util.List;

import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BedOps;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.bed.BEDFeature;

/**
 * Note this class should only be used for annotating CNV calls from WES data (exon calls) still
 * pretty ugly
 */
public class WESMappabilityAnnotator extends BEDFileAnnotator {

	private final GeneAnnotator geneAnnotator;
	private final Logger log;



	/**
	 * @param geneAnnotator used for getting exon bounderies
	 * @param mapBed mappability bed file
	 * @param log
	 */
	public WESMappabilityAnnotator(GeneAnnotator geneAnnotator, String mapBed, Logger log) {
		super(mapBed);
		this.geneAnnotator = geneAnnotator;
		this.log = log;
	}


	/**
	 * 
	 *
	 * 
	 * @param geneAnnotator this is required for accessing exon start stops
	 * @param build the genome build for the mappability track
	 * @param log
	 * @return
	 */
	public static WESMappabilityAnnotator getDefaultAnnotator(GeneAnnotator geneAnnotator,
																														GENOME_BUILD build, Logger log) {
		String mapBed = Resources.genome(build, log).get100MerMappabilityTrack().get();
		BedOps.verifyBedIndex(mapBed, log);
		return new WESMappabilityAnnotator(geneAnnotator, mapBed, log);
	}



	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.genvisis.cnv.annotation.segments.SegmentAnnotator#annotate(org.genvisis.filesys.Segment,
	 * org.genvisis.cnv.annotation.segments.SegmentAnotation)
	 */
	@Override
	public void annotate(Segment segment, SegmentAnotation segmentAnotation) {
		if (segmentAnotation.getAttributes()
												.containsKey(SegmentAnnotationKeys.MAPPABILITY.toString())) {
			String error = "This annotation has already been used for Gene lookup";
			log.reportError(error);
			throw new IllegalArgumentException(error);
		}
		// get all mapability segments overlapping segment of interest
		CloseableIterator<BEDFeature> iterator = query(	Positions.getChromosomeUCSC(segment.getChr(),
																																								true),
																										segment.getStart(), segment.getStop());
		List<String> values = new ArrayList<String>();
		double numBases = 0;
		double cumulativeMapScore = 0;
		boolean found = false;

		while (iterator.hasNext()) {
			MappabilityFeature feature = new MappabilityFeature(iterator.next(), log);

			double mapScore = feature.getMappability(log);

			// Since this is WES, we limit only to those overlapping exons, for array, we will limit to
			// markers
			GeneData[] geneDatas = geneAnnotator.getGeneSet().getOverLappingLoci(feature);// should
																																										// typically be
																																										// one gene
			if (geneDatas != null && geneDatas.length > 0) {
				for (GeneData geneData : geneDatas) {
					int[][] exons = geneData.getExonBoundaries();
					for (int[] exon : exons) {
						Segment exSeq = new Segment(feature.getChr(), exon[0], exon[1]);
						int overlap = exSeq.amountOfOverlapInBasepairs(segment);// has to overlap original
																																		// segment
						if (overlap > 0) {
							numBases += overlap;
							cumulativeMapScore += mapScore * overlap; // per nucleotide overlap
							found = true;
						}
					}
				}
			}
		}
		if (found) {
			double averageMapScore = cumulativeMapScore / numBases;
			values.add(Double.toString(averageMapScore));
		} else {
			values.add(SegmentAnnotationKeys.MAPPABILITY.getMissingValue());
		}
		segmentAnotation.getAttributes().put(SegmentAnnotationKeys.MAPPABILITY.toString(), values);
	}
}
