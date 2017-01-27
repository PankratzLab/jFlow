/**
 * 
 */
package org.genvisis.cnv.annotation.segments;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.genvisis.cnv.manage.Resources;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;

/**
 *
 *
 * Relies on annotations being present from {@link GeneAnnotator} and annotates each gene with GDI
 * score,
 *
 * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4640721/
 *
 */
public class GDIAnnotator implements SegmentAnnotator {


	private static final String[] REQUIRED_HEADER = new String[] {"GENE", "GDI-Percentile", "GDI",
																																"GDI-Phred"};


	/**
	 * Key is gene name, value is GDI Percentile
	 */
	private final Map<String, String> gdiLookup;
	private final Logger log;

	private GDIAnnotator(Map<String, String> gdiLookup, Logger log) {
		super();
		this.gdiLookup = gdiLookup;
		this.log = log;
	}



	public Map<String, String> getGdiLookup() {
		return gdiLookup;
	}



	/**
	 * @param log
	 * @return GDI annotator stored in resources
	 */
	public static GDIAnnotator getDefaultGDIAnnotator(Logger log) {
		String gdiFile = Resources.annotation(log).getGDI().get();
		int[] indices = ext.indexFactors(	REQUIRED_HEADER, Files.getHeaderOfFile(gdiFile, log), true,
																			false);

		Map<String, String> attributes = new HashMap<String, String>();

		if (Array.countIf(indices, -1) > 0) {
			throw new IllegalArgumentException("Invalid header in "+ gdiFile + " , require "
																					+ Array.toStr(REQUIRED_HEADER));
		}

		try {
			BufferedReader reader = Files.getAppropriateReader(gdiFile);
			reader.readLine();// skip header
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				if (attributes.containsKey(line[indices[0]])) {
					log.reportTimeWarning("Multiple entries for " + line[indices[0]] + " in " + gdiFile);
				} else {
					attributes.put(	line[indices[0]],
													line[indices[1]] + ";" + line[indices[2]] + ";" + line[indices[3]]);
				}
			}
			reader.close();
		} catch (FileNotFoundException e) {
			log.reportException(e);

			return null;
		} catch (IOException e) {
			log.reportException(e);
			return null;

		}

		return new GDIAnnotator(attributes, log);
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
		if (!segmentAnotation.getAttributes().containsKey(SegmentAnnotationKeys.GENE.toString())) {
			String error = "This annotation requires previous gene lookup";
			log.reportError(error);
			throw new IllegalArgumentException(error);

		}


		if (segmentAnotation.getAttributes().containsKey(SegmentAnnotationKeys.GDI.toString())) {
			String error = "This annotation has already been used for GDI lookup";
			log.reportError(error);
			throw new IllegalArgumentException(error);
		}



		List<String> genes = segmentAnotation	.getAttributes()
																					.get(SegmentAnnotationKeys.GENE.toString());
		List<String> gdis = new ArrayList<String>();
		for (String gene : genes) {
			if (gdiLookup.containsKey(gene)) {
				gdis.add(gdiLookup.get(gene));
			} else {
				gdis.add(SegmentAnnotationKeys.GDI.getMissingValue());
			}
		}
		segmentAnotation.getAttributes().put(SegmentAnnotationKeys.GDI.toString(), gdis);
	}
}


