package seq.manage;

import java.io.File;
import java.io.IOException;
import java.util.List;

import common.Logger;
import common.Positions;
import filesys.Segment;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReader.Indexing;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * Class for common bamFile manips
 *
 */
public class BamOps {
	public static final String BAM_EXT = ".bam";
	public static final String BAI_EXT = ".bai";

	/**
	 * @param bamOrSam
	 *            .bam or .sam file
	 * @param stringency
	 *            Stringency validation for the records
	 * @return new reader
	 */
	public static SamReader getDefaultReader(String bamOrSam, ValidationStringency stringency) {
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
		samReaderFactory.validationStringency(stringency);
		SamReader reader = samReaderFactory.open(new File(bamOrSam));
		return reader;
	}

	/**
	 * @param segs
	 *            Genvisis type segments to search
	 * @param sFileHeader
	 *            an {@link SAMFileHeader}
	 * @param bpBuffer
	 *            bp buffer to be added to the segments
	 * @param optimize
	 *            perform an extra step to make a more efficient query
	 * @param log
	 * @return array of {@link QueryInterval} that can be queried by a bamfile reader
	 */
	public static QueryInterval[] convertSegsToQI(Segment[] segs, SAMFileHeader sFileHeader, int bpBuffer, boolean optimize, Logger log) {
		QueryInterval[] qIntervals = new QueryInterval[segs.length];
		segs = Segment.sortSegments(segs);
		for (int i = 0; i < qIntervals.length; i++) {
			String sequenceName = Positions.getChromosomeUCSC(segs[i].getChr(), true);
			int referenceIndex = sFileHeader.getSequenceIndex(sequenceName);
			if (referenceIndex < 0) {
				log.reportError("Error - could not find " + sequenceName + " in the sequence dictionary, halting");
				return null;
			}
			qIntervals[i] = new QueryInterval(referenceIndex, segs[i].getStart() - bpBuffer, segs[i].getStop() + bpBuffer);
		}
		if (optimize) {
			qIntervals = QueryInterval.optimizeIntervals(qIntervals);
		}
		return qIntervals;
	}

	public static Segment[] converQItoSegs(QueryInterval[] qIntervals, SAMFileHeader sFileHeader, Logger log) {
		Segment[] segs = new Segment[qIntervals.length];
		for (int i = 0; i < segs.length; i++) {
			segs[i] = new Segment(Positions.chromosomeNumber(sFileHeader.getSequence(qIntervals[i].referenceIndex).getSequenceName()), qIntervals[i].start, qIntervals[i].end);
		}
		return segs;
	}

	public static SAMFileHeader getHeader(String bamfile) {
		SamReader reader = getDefaultReader(bamfile, ValidationStringency.STRICT);
		SAMFileHeader samFileHeader = reader.getFileHeader();
		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return samFileHeader;

	}

	public static String[] getSampleNames(String[] bamFiles) {
		String[] sampleNames = new String[bamFiles.length];
		for (int i = 0; i < sampleNames.length; i++) {
			sampleNames[i] = getSampleName(bamFiles[i]);
		}
		return sampleNames;
	}

	/**
	 * @author lane0212 Stores some simple counts that can be quickly retrieved from the index file
	 */
	public static class BamIndexStats {
		private int alignedRecordCount;
		private int unalignedRecordCount;

		public BamIndexStats(int alignedRecordCount, int unalignedRecordCount) {
			super();
			this.alignedRecordCount = alignedRecordCount;
			this.unalignedRecordCount = unalignedRecordCount;
		}

		public int getAlignedRecordCount() {
			return alignedRecordCount;
		}

		public int getUnalignedRecordCount() {
			return unalignedRecordCount;
		}
	}

	public static BamIndexStats getBamIndexStats(String bamFile) {
		SamReader reader = BamOps.getDefaultReader(bamFile, ValidationStringency.STRICT);
		return getBamIndexStats(reader);
	}

	public static BamIndexStats getBamIndexStats(SamReader reader) {
		BAMIndexMetaData[] result = getIndexMetaData(reader);
		int alignedRecordCount = 0;
		int unalignedRecordCount = 0;
		for (int i = 0; i < result.length; i++) {
			alignedRecordCount += result[i].getAlignedRecordCount();
			unalignedRecordCount += result[i].getUnalignedRecordCount();
		}

		return new BamIndexStats(alignedRecordCount, unalignedRecordCount);

	}

	public static BAMIndexMetaData[] getIndexMetaData(SamReader reader) {
		Indexing index = reader.indexing();
		BAMIndex bamIndex = index.getIndex();
		List<SAMSequenceRecord> records = reader.getFileHeader().getSequenceDictionary().getSequences();
		BAMIndexMetaData[] result = new BAMIndexMetaData[records.size()];
		for (int i = 0; i < result.length; i++) {
			result[i] = bamIndex.getMetaData(i);
		}
		return result;

		// bamIndex.g
		//
		// bamIndex.getMetaData(0).

		// int nRefs = bamIndex.getNumberOfReferences();
		//
		//
		// AbstractBAMFileIndex index = (AbstractBAMFileIndex) bam.getIndex();
		// // read through all the bins of every reference.
		// BAMIndexMetaData[] result = new BAMIndexMetaData[nRefs == 0 ? 1 : nRefs];
		// for (int i = 0; i < nRefs; i++) {
		// result[i] = index.getMetaData(i);
		// }
		//
		// if (result[0] == null){
		// result[0] = new BAMIndexMetaData();
		// }
		// final Long noCoordCount = index.getNoCoordinateCount();
		// if (noCoordCount != null) // null in old index files without metadata
		// result[0].setNoCoordinateRecordCount(noCoordCount);
		//
		// return result;
		// }
		//
	}

	public static String getSampleName(String bamFile) {
		SamReader reader = getDefaultReader(bamFile, ValidationStringency.STRICT);
		String sample = reader.getFileHeader().getReadGroups().get(0).getSample();
		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return sample;
	}
}
