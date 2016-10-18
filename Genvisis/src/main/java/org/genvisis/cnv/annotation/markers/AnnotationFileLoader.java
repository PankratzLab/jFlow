package org.genvisis.cnv.annotation.markers;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.ReadingFilePrep;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * @author lane0212 Does tabix-based loading for annotations
 */
public abstract class AnnotationFileLoader extends AnnotationFile implements ReadingFilePrep {
	private final boolean indexRequired;
	private final boolean valid;
	private int reportEvery;

	/**
	 * @param proj
	 * @param annotations these annotations will be checked for in the annotation file, can be null...
	 * @param annotationFilename the file must be in proper vcf format and contain all annotations
	 *        requested
	 * @param indexRequired this should always be true, but may be optional at a later time
	 */
	public AnnotationFileLoader(Project proj, AnalysisParams[] params, Annotation[] annotations,
															String annotationFilename, boolean indexRequired) {
		super(proj, annotationFilename);
		setAnnotations(annotations);
		setParams(params);
		this.indexRequired = indexRequired;
		valid = validate();
		reportEvery = -1;

	}

	public void setReportEvery(int reportEvery) {
		this.reportEvery = reportEvery;
	}

	public enum QUERY_TYPE {
		/**
		 * It is assumed that the {@link AnnotationParser} will be found exactly once
		 */
		ONE_TO_ONE,
		/**
		 * No assumptions are made on whether the {@link AnnotationParser} will be found
		 */
		DISCRETE_LIST;
	}

	/**
	 * @param segs
	 * @param parsersQueries checks and annotates these queries
	 */
	public void query(Segment[] segs, List<Map<String, ? extends AnnotationParser>> parsersQueries,
										QUERY_TYPE queryType) {

		AnnotationQuery annotationQuery = getAnnotationQuery(segs);
		int index = 0;
		while (annotationQuery.hasNext()) {
			VariantContext vc = annotationQuery.next();
			index++;
			if (reportEvery > 0 && index % reportEvery == 0) {
				proj.getLog().reportTimeInfo(index + " -loaded");
			}
			for (Map<String, ? extends AnnotationParser> parsers : parsersQueries) {
				AnnotationParser parser = parsers.get(vc.getID());
				if (parser != null) {
					parser.parseAnnotation(vc, proj.getLog());
					parser.setFound(true);
				}
			}
		}
		validateSearch(parsersQueries, queryType);
	}

	/**
	 * Make sure that every {@link AnnotationParser } was found
	 */
	private void validateSearch(List<Map<String, ? extends AnnotationParser>> parsersQueries,
															QUERY_TYPE queryType) {
		switch (queryType) {
			case DISCRETE_LIST:
				break;
			case ONE_TO_ONE:
				boolean allFound = true;
				for (Map<String, ? extends AnnotationParser> parsers : parsersQueries) {
					for (AnnotationParser parser : parsers.values()) {
						if (!parser.isFound()) {
							allFound = false;
						}
					}

				}
				if (!allFound) {
					String error = "Did not find all queries for type " + queryType
												 + " , missing annotations or internal bug";
					proj.getLog().reportError(error);
					throw new IllegalStateException(error);
				}
				break;
			default:
				break;
		}
	}

	public AnnotationQuery getAnnotationQuery() {
		return getAnnotationQuery(null);
	}

	public AnnotationQuery getAnnotationQuery(Segment[] segs) {
		if (valid) {
			AnnotationQuery annotationIterator = new AnnotationQuery(annotationFilename, segs,
																															 indexRequired, proj.getLog());
			return annotationIterator;
		} else {
			proj.getLog().reportError("Invalid loader...");
			return null;
		}
	}

	@Override
	public void init() {

	}

	@Override
	public boolean validate() {
		if (!Files.exists(annotationFilename)) {
			proj.getLog().reportError("Could not find annotation file " + annotationFilename);
			return false;
		} else {
			try {
				VCFFileReader reader = new VCFFileReader(new File(annotationFilename), indexRequired);
				VCFHeader vcfHeader = reader.getFileHeader();// doing this will trigger the htsjdk file
																										 // format checks
				if (params != null) {
					for (AnalysisParams param : params) {
						if (vcfHeader.getMetaDataLine(param.getKey()) != null) {
							param.parseHeaderLine(vcfHeader.getMetaDataLine(param.getKey()));
						} else {
							proj.getLog().reportTimeWarning("Could not find parameters for " + param.getKey());
						}
					}
				}
				boolean hasAllAnno = true;
				if (annotations != null) {
					for (int i = 0; i < annotations.length; i++) {
						if (!vcfHeader.hasInfoLine(annotations[i].getName())) {
							proj.getLog().reportError("Could not find annotation " + annotations[i].getName()
																				+ " in " + annotationFilename);
							hasAllAnno = false;
						}
					}
					reader.close();
				}
				return hasAllAnno;

			} catch (TribbleException trib) {
				proj.getLog()
						.reportError("Index was required and failed to load it for " + annotationFilename);

				return false;
			} catch (Exception e) {
				proj.getLog()
						.reportError("Could not properly initialize reader for  " + annotationFilename);
				proj.getLog().reportException(e);
				return false;

			}
		}
	}

	@Override
	public void close() {

	}

	/**
	 * @author lane0212 Class that returns each {@link VariantContext} by the {@link Segment} query
	 */
	public static class AnnotationQuery implements Iterator<VariantContext> {
		private final VCFFileReader vcfFileReader;
		private int currentIndex;
		private final QueryInterval[] queryIntervals;
		private final VCFHeader vcfHeader;
		private CloseableIterator<VariantContext> currentIterator;

		/**
		 * @param annotationFile the file to load from
		 * @param segs can be null, if not null, only these regions will be returned. Otherwise the
		 *        iterator will traverse the entire file
		 * @param requireIndex should always be true
		 * @param log
		 */
		public AnnotationQuery(String annotationFile, Segment[] segs, boolean requireIndex,
													 Logger log) {
			super();
			vcfFileReader = new VCFFileReader(new File(annotationFile), requireIndex);
			vcfHeader = vcfFileReader.getFileHeader();
			currentIndex = 0;
			queryIntervals = segs == null ? null : VCFOps.convertSegsToQI(segs, vcfHeader, 0, true, log);
			currentIterator = queryIntervals == null ? vcfFileReader.iterator()
																							 : vcfFileReader.query(vcfHeader.getSequenceDictionary()
																																							.getSequence(queryIntervals[currentIndex].referenceIndex)
																																							.getSequenceName(),
																																		 queryIntervals[currentIndex].start,
																																		 queryIntervals[currentIndex].end);
		}

		@Override
		public boolean hasNext() {
			boolean hasNext = false;
			if (currentIterator.hasNext()) {
				hasNext = true;
			} else if (queryIntervals != null) {// try the next interval if they exist
				while (!currentIterator.hasNext()) {
					currentIterator.close();
					currentIndex++;
					if (currentIndex >= queryIntervals.length) {
						break;
					}
					currentIterator = vcfFileReader.query(
																								vcfHeader.getSequenceDictionary()
																												 .getSequence(queryIntervals[currentIndex].referenceIndex)
																												 .getSequenceName(),
																								queryIntervals[currentIndex].start,
																								queryIntervals[currentIndex].end);
				}
				hasNext = currentIterator.hasNext();
			}
			if (!hasNext) {
				vcfFileReader.close();
			}
			return hasNext;

		}

		@Override
		public VariantContext next() {
			return currentIterator.next();
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}

}
