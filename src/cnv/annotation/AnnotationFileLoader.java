package cnv.annotation;

import java.util.Iterator;
import java.util.List;

import seq.manage.VCFOps;
import filesys.Segment;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import common.Files;
import common.Logger;
import cnv.filesys.Project;
import cnv.filesys.ReadingFilePrep;

/**
 * @author lane0212 Does tabix-based loading for annotations
 */
public abstract class AnnotationFileLoader extends AnnotationFile implements ReadingFilePrep {
	private boolean indexRequired;
	private boolean valid;
	private int reportEvery;

	/**
	 * @param proj
	 * @param annotations
	 *            these annotations will be checked for in the annotation file, can be null...
	 * @param annotationFilename
	 *            the file must be in proper vcf format and contain all annotations requested
	 * @param indexRequired
	 *            this should always be true, but may be optional at a later time
	 */
	public AnnotationFileLoader(Project proj, Annotation[] annotations, String annotationFilename, boolean indexRequired) {
		super(proj, annotationFilename);
		setAnnotations(annotations);
		this.indexRequired = indexRequired;
		this.valid = validate();
		this.reportEvery = -1;

	}

	public void setReportEvery(int reportEvery) {
		this.reportEvery = reportEvery;
	}

	public enum QUERY_ORDER {
		/**
		 * It is assumed that the {@link AnnotationParser} will be found exactly once and in order of the file
		 */
		ONE_PER_IN_ORDER, /**
		 * No assumptions are made for order or whether the {@link AnnotationParser} will be found
		 */
		NO_ORDER;
	}

	/**
	 * @param segs
	 * @param parsersQueries
	 *            checks and annotates these queries
	 */
	public void query(Segment[] segs, List<AnnotationParser[]> parsersQueries, QUERY_ORDER queryType) {
		int[][] search = getSearchSpace(parsersQueries, queryType);

		AnnotationQuery annotationQuery = getAnnotationQuery(segs);
		int index = 0;
		while (annotationQuery.hasNext()) {
			VariantContext vc = annotationQuery.next();
			index++;
			if (reportEvery > 0 && index % reportEvery == 0) {
				proj.getLog().reportTimeInfo(index + " -loaded");
			}
			int searchIndex = 0;
			for (AnnotationParser[] parsers : parsersQueries) {
				int stopSearch = parsers.length;
				int startSearch = 0;
				switch (queryType) {
				case NO_ORDER:
					for (int i = startSearch; i < stopSearch; i++) {

						if (parsers[i].shouldAnnotateWith(vc, proj.getLog())) {
							parsers[i].parseAnnotation(vc, proj.getLog());
							parsers[i].setFound(true);
							search[searchIndex][0] = search[searchIndex][0] + 1;
							search[searchIndex][1] = search[searchIndex][1] + 1;
						} else {

						}
					}
					break;
				case ONE_PER_IN_ORDER:
					if (parsers[index - 1].shouldAnnotateWith(vc, proj.getLog())) {
						String error = "Query type was set to " + queryType + " but the annotation at index " + (index - 1) + " did not match";
						proj.getLog().reportTimeError(error);
						throw new IllegalStateException(error);
					} else {
						parsers[index - 1].parseAnnotation(vc, proj.getLog());
						parsers[index - 1].setFound(true);
					}
					break;
				default:
					String error = "Invalid query " + queryType;
					proj.getLog().reportTimeInfo(error);
					throw new IllegalArgumentException(error);
				}

				switch (queryType) {
				case NO_ORDER:
					break;
				case ONE_PER_IN_ORDER:
					//searchIndex++;
					break;
				default:
					String error = "Invalid query " + queryType;
					proj.getLog().reportTimeInfo(error);
					throw new IllegalArgumentException(error);
				}
			}
		}
		validateSearch(parsersQueries, queryType);
	}

	/**
	 * Make sure that every {@link AnnotationParser } was found
	 */
	private void validateSearch(List<AnnotationParser[]> parsersQueries, QUERY_ORDER queryType) {
		switch (queryType) {
		case NO_ORDER:
			break;
		case ONE_PER_IN_ORDER:
			boolean allFound = true;
			for (int i = 0; i < parsersQueries.size(); i++) {
				for (int j = 0; j < parsersQueries.get(i).length; j++) {
					if (!parsersQueries.get(i)[j].isFound()) {
						allFound = false;
					}
				}

			}
			if (!allFound) {
				String error = "Did not find all queries for type " + queryType + " , missing annotations or internal bug";
				proj.getLog().reportTimeError(error);
				throw new IllegalStateException(error);
			}
			break;
		default:
			break;
		}
	}

	private int[][] getSearchSpace(List<AnnotationParser[]> parsersQueries, QUERY_ORDER queryType) {
		int[][] search = new int[parsersQueries.size()][];
		for (int i = 0; i < parsersQueries.size(); i++) {
			search[i] = new int[2];
			switch (queryType) {
			case NO_ORDER:
				search[i][0] = 0;
				search[i][1] = parsersQueries.get(i).length;

				break;
			case ONE_PER_IN_ORDER:
				search[i][0] = 0;
				search[i][1] = 2;

				break;
			default:
				break;

			}
		}
		return search;
	}

	public AnnotationQuery getAnnotationQuery() {
		return getAnnotationQuery(null);
	}

	public AnnotationQuery getAnnotationQuery(Segment[] segs) {
		if (valid) {
			AnnotationQuery annotationIterator = new AnnotationQuery(annotationFilename, segs, indexRequired, proj.getLog());
			return annotationIterator;
		} else {
			proj.getLog().reportTimeError("Invalid loader...");
			return null;
		}
	}

	@Override
	public void init() {

	}

	@Override
	public boolean validate() {
		if (!Files.exists(annotationFilename)) {
			proj.getLog().reportTimeError("Could not find annotation file " + annotationFilename);
			return false;
		} else {
			try {
				VCFFileReader reader = new VCFFileReader(annotationFilename, indexRequired);
				VCFHeader vcfHeader = reader.getFileHeader();// doing this will trigger the htsjdk file format checks

				boolean hasAllAnno = true;
				if (annotations != null) {
					for (int i = 0; i < annotations.length; i++) {
						if (!vcfHeader.hasInfoLine(annotations[i].getName())) {
							proj.getLog().reportTimeError("Could not find annotation " + annotations[i].getName() + " in " + annotationFilename);
							hasAllAnno = false;
						}
					}
					reader.close();
				}
				return hasAllAnno;

			} catch (TribbleException trib) {
				proj.getLog().reportTimeError("Index was required and failed to load it for " + annotationFilename);

				return false;
			} catch (Exception e) {
				proj.getLog().reportTimeError("Could not properly initialize reader for  " + annotationFilename);
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
		private VCFFileReader vcfFileReader;
		private int currentIndex;
		private QueryInterval[] queryIntervals;
		private VCFHeader vcfHeader;
		private CloseableIterator<VariantContext> currentIterator;

		/**
		 * @param annotationFile
		 *            the file to load from
		 * @param segs
		 *            can be null, if not null, only these regions will be returned. Otherwise the iterator will traverse the entire file
		 * @param requireIndex
		 *            should always be true
		 * @param log
		 */
		public AnnotationQuery(String annotationFile, Segment[] segs, boolean requireIndex, Logger log) {
			super();
			this.vcfFileReader = new VCFFileReader(annotationFile, requireIndex);
			this.vcfHeader = vcfFileReader.getFileHeader();
			this.currentIndex = 0;
			this.queryIntervals = segs == null ? null : VCFOps.convertSegsToQI(segs, vcfHeader, 0, true, log);
			this.currentIterator = queryIntervals == null ? vcfFileReader.iterator() : vcfFileReader.query(vcfHeader.getSequenceDictionary().getSequence(queryIntervals[currentIndex].referenceIndex).getSequenceName(), queryIntervals[currentIndex].start, queryIntervals[currentIndex].end);
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
					currentIterator = vcfFileReader.query(vcfHeader.getSequenceDictionary().getSequence(queryIntervals[currentIndex].referenceIndex).getSequenceName(), queryIntervals[currentIndex].start, queryIntervals[currentIndex].end);
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
	}

}
