package seq.manage;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;

import common.Logger;
import common.Positions;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import filesys.Segment;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
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
	public static QueryInterval[] convertSegsToQI(Segment[] segs, SAMFileHeader sFileHeader, int bpBuffer, boolean optimize, boolean chr, Logger log) {
		QueryInterval[] qIntervals = new QueryInterval[segs.length];
		segs = Segment.sortSegments(segs);
		for (int i = 0; i < qIntervals.length; i++) {
			String sequenceName = Positions.getChromosomeUCSC(segs[i].getChr(), chr);
			int referenceIndex = sFileHeader.getSequenceIndex(sequenceName);
			if (referenceIndex < 0) {
				referenceIndex = sFileHeader.getSequenceIndex(sequenceName +"T");// MT
				if (referenceIndex < 0) {
					log.reportError("Error - could not find " + sequenceName + " in the sequence dictionary, halting");
					return null;
				}
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

	public static String[] getAllBarCodes(String[] bams, Logger log) {
		HashSet<String> unique = new HashSet<String>();
		for (int i = 0; i < bams.length; i++) {
			unique.addAll(getBarcodesFor(bams[i], log));
		}
		return unique.toArray(new String[unique.size()]);
	}

	public static ArrayList<String> getBarcodesFor(String bam, Logger log) {
		ArrayList<String> barcodes = new ArrayList<String>();
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
		samReaderFactory.validationStringency(ValidationStringency.LENIENT);
		SamReader reader = samReaderFactory.open(new File(bam));
		List<SAMReadGroupRecord> rgs = reader.getFileHeader().getReadGroups();
		HashSet<String> barcodesUnique = new HashSet<String>();
		for (SAMReadGroupRecord samReadGroupRecord : rgs) {
			String[] id = samReadGroupRecord.getId().split("_");
			String[] tmpCodes = id[id.length - 3].split("-");
			if (tmpCodes.length != 2) {
				throw new IllegalArgumentException("Could not parse barcodes for RG " + samReadGroupRecord + " in bam file " + bam);

			} else {
				for (int i = 0; i < tmpCodes.length; i++) {
					if (tmpCodes[i].replaceAll("A", "").replaceAll("C", "").replaceAll("T", "").replaceAll("G", "").length() != 0) {
						throw new IllegalArgumentException("Invalid barcode " + tmpCodes[i]);
					} else {
						barcodesUnique.add(tmpCodes[i]);
					}
				}
			}

		}
		barcodes.addAll(barcodesUnique);
		return barcodes;
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

	public static int estimateReadSize(String bamFile, int numReads, Logger log) {
		SamReader reader = getDefaultReader(bamFile, ValidationStringency.STRICT);
		int readsize = 0;
		SAMRecordIterator iterator = reader.iterator();
		
		int readsScanned = 0;
		while (iterator.hasNext()) {
			SAMRecord samRecord = iterator.next();
			if (!samRecord.getReadUnmappedFlag() && samRecord.getCigar().getCigarElements().size() == 1) {

				readsize += samRecord.getReadLength();
				readsScanned++;

			}
			if (readsScanned > numReads) {
				break;
			}
		}
		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if (readsScanned > 0) {
			double avg = (double) readsize / readsScanned;
			return (int) avg;
		}
		return -1;
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

	private static class SampleNameExtractor implements Callable<SampleNameExtractor> {
		private String bamFile;
		private String sampleName;

		public SampleNameExtractor(String bamFile) {
			super();
			this.bamFile = bamFile;
		}

		@Override
		public SampleNameExtractor call() throws Exception {
			this.sampleName = getSampleName(bamFile);
			return this;
		}

		public String getBamFile() {
			return bamFile;
		}

		public String getName() {
			return sampleName;
		}

	}

	private static class SampleNameProducer implements Producer<SampleNameExtractor> {

		private String[] bamFiles;
		private int index;

		public SampleNameProducer(String[] bamFiles) {
			super();
			this.bamFiles = bamFiles;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < bamFiles.length;
		}

		@Override
		public Callable<SampleNameExtractor> next() {
			SampleNameExtractor ex = new SampleNameExtractor(bamFiles[index]);
			index++;
			return ex;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	/**
	 * Designed to take the sample names from {@link VCFOps#getSamplesInFile(htsjdk.variant.vcf.VCFFileReader)} and match to an array of bams
	 * 
	 * @param samples
	 *            samples to match
	 * @param variantSets
	 *            variant sets that may be appended to the vcf sample names
	 * @param bams
	 *            the bam files
	 * @param log
	 * @return Hashtable of the sample -> bam file mapping
	 */
	public static Hashtable<String, String> matchToVcfSamplesToBamFiles(String[] samples, Set<String> variantSets, String[] bams, int numThreads, Logger log) {
		Hashtable<String, String> matched = new Hashtable<String, String>();
		Hashtable<String, String> bamSamples = new Hashtable<String, String>();
		SampleNameProducer producer = new SampleNameProducer(bams);
		WorkerTrain<SampleNameExtractor> train = new WorkerTrain<SampleNameExtractor>(producer, numThreads, 10, log);

		while (train.hasNext()) {
			SampleNameExtractor ex = train.next();
			String bamSamp = ex.getName();
			if (bamSamples.containsKey(bamSamp)) {
				throw new IllegalArgumentException("Bams must be sample unique");
			} else {
				bamSamples.put(bamSamp, ex.getBamFile());
				if (variantSets != null) {
					for (String set : variantSets) {
						bamSamples.put(bamSamp + set, ex.getBamFile());
					}
				}
			}
		}

		for (int i = 0; i < samples.length; i++) {
			if (!bamSamples.containsKey(samples[i])) {
				log.reportTimeWarning("Did not find matching bam file for " + samples[i]);
			} else {
				if (matched.contains(samples[i])) {
					throw new IllegalArgumentException("Multiple bam files matched sample " + samples[i] + ", perhaps because of variant sets?");
				}
				matched.put(samples[i], bamSamples.get(samples[i]));
			}
		}

		log.reportTimeInfo("Found matching bam files for" + matched.size() + " of " + samples.length + " samples");

		return matched;

	}

}
