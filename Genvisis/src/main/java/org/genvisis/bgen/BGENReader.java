package org.genvisis.bgen;

import java.io.Closeable;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import org.genvisis.bgen.BGENReader.BGENRecord;

// ZSTD library:
// https://github.com/prestodb/presto/tree/master/presto-orc/src/main/java/com/facebook/presto/orc/zstd
// or
// https://github.com/airlift/aircompressor

public class BGENReader implements Closeable, Iterable<BGENRecord> {

	private BGENReader() {/**/}

	/**
	 * Open a BGEN file.<br />
	 * <br />
	 * If {@code buildMap} is set to {@code TRUE}, the reader will scan through the file and parse the
	 * variant record header blocks, storing the data into a genomically-ordered set of record headers
	 * along with their associated locations in file. This makes later reading of the file extremely
	 * fast, but carries a significant initial time cost. Intended usage should therefore be
	 * thoughtfully considered and if a streaming model is the best-fit, the {@code buildMap}
	 * parameter should be set to {@code FALSE} and the various {@code query()} methods used instead.<br />
	 * <br />
	 * If a serialized file of the map data exists (built with
	 * {@link BGENTools#serializeMapInfo(BGENReader, String)}), {@code buildMap} should be set to
	 * {@code FALSE} and after this method, {@link BGENTools#loadMapInfo(BGENReader, String)} should
	 * be called.
	 * 
	 * 
	 * @param file BGEN file
	 * @param buildMap Scan through file and build map information
	 * @throws IOException
	 */
	public static BGENReader open(String file, boolean buildMap) throws IOException {
		BGENReader reader = new BGENReader();
		reader.openFile(file, buildMap);
		return reader;
	}

	/**
	 * Indicates {@code close()} has been called.
	 */
	private boolean closed = false;
	/**
	 * Filename of BGEN file
	 */
	private String filename;
	/**
	 * Underlying file object
	 */
	private RandomAccessFile raf;
	/**
	 * Start of records in bytes
	 */
	private long recordStart;
	/**
	 * Length of header in bytes
	 */
	private long headerLength;
	/**
	 * Number of records (from header)
	 */
	private long recordCount;
	/**
	 * Number of samples (from header)
	 */
	private long sampleCount;
	/**
	 * Magic bytes!
	 */
	private byte[] magicBytes;
	/**
	 * Compression, layout version, and sample inclusion flags.
	 */
	private long[] flags;
	/**
	 * Compression used in this file
	 */
	private COMPRESSION compress;
	/**
	 * Layout version of this file
	 */
	private LAYOUT layout;
	/**
	 * Flag for sample names in the file
	 */
	private boolean hasSampleNames;
	/**
	 * List of sample names
	 */
	private String[] samples;
	/**
	 * Internal map of byte locations to the metadata located at those locations
	 */
	private Map<Long, BGENRecordMetaData> ptrMap;
	/**
	 * Genomically-ordered chromosome-mapped variant metadata
	 */
	Map<Integer, TreeSet<BGENRecordMetaData>> chrSets;

	/**
	 * 
	 * @param file
	 * @param buildMap
	 * @throws IOException
	 */
	private void openFile(String file, boolean buildMap) throws IOException {
		filename = file;
		raf = new RandomAccessFile(file, "r");
		ptrMap = new HashMap<>();
		chrSets = new HashMap<Integer, TreeSet<BGENRecordMetaData>>();
		readHeader();
		if (buildMap) {
			readMap();
		}
	}

	private void readHeader() throws IOException {
		byte[] read;
		read = new byte[4];
		raf.read(read);
		setRecordStart(BGENBitMath.bytesToLong(read, true) + 4);

		raf.read(read);
		headerLength = BGENBitMath.bytesToLong(read, true);

		raf.read(read);
		recordCount = BGENBitMath.bytesToLong(read, true);

		raf.read(read);
		sampleCount = BGENBitMath.bytesToLong(read, true);

		raf.read(magicBytes = new byte[4]);

		if (!testMagicBytes(magicBytes)) {
			raf.close();
			throw new IOException(
														"Magic bytes were incorrect, expected {'b','g','e','n'} or {'0','0','0','0'}, found {"
																+
																((char) magicBytes[0]) + "," +
																((char) magicBytes[1]) + "," +
																((char) magicBytes[2]) + "," +
																((char) magicBytes[3]) + "}.");
		}

		if (headerLength - 20 > 0) {
			for (int i = 0; i < headerLength - 20; i++) {
				raf.read();
			}
		}

		raf.read(read);
		flags = readFlags(read);
		checkFlagsOrException(flags);
		compress = COMPRESSION.values()[(int) flags[0]];
		layout = LAYOUT.values()[(int) flags[1] - 1];
		hasSampleNames = flags[2] == 1;

		if (hasSampleNames) {
			samples = new String[(int) sampleCount];
			raf.read(read);
			long numBytesSampIDBlock = BGENBitMath.bytesToLong(read, true);
			long readB = 8;
			raf.read(read);
			long sampleCount2 = BGENBitMath.bytesToLong(read, true);
			if (sampleCount != sampleCount2) {
				raf.close();
				throw new IOException(
															"Sample count in header block must equal sample count in sample block, found "
																	+ sampleCount + " in header and " + sampleCount2
																	+ " in sample block.");
			}

			byte[] lenByt = new byte[2];
			byte[] identByt = new byte[0];
			for (int i = 0; i < sampleCount; i++) {
				raf.read(lenByt);
				readB += lenByt.length;
				int len = (int) BGENBitMath.bytesToLong(lenByt, true);
				if (len != identByt.length) {
					identByt = new byte[len];
				}
				raf.read(identByt);
				readB += len;
				samples[i] = new String(identByt, StandardCharsets.US_ASCII);
			}
			if (readB != numBytesSampIDBlock) {
				throw new IOException(
															"Number of bytes recorded for sample ID header does not match the number of bytes read.  Expected "
																	+ numBytesSampIDBlock + " | Read " + readB);
			}
		} else {
			initSamples();
		}
		raf.seek(getRecordStartByte());
	}

	private void initSamples() {
		samples = new String[(int) sampleCount];
		for (int i = 0; i < samples.length; i++) {
			samples[i] = "(unknown_sample_" + (1 + i) + ")";
		}
	}

	private final boolean testMagicBytes(byte[] magicBytes) {
		boolean len = magicBytes.length == 4;
		boolean b1 = magicBytes[0] == 'b' || magicBytes[0] == 0;
		boolean b2 = magicBytes[1] == 'g' || magicBytes[1] == 0;
		boolean b3 = magicBytes[2] == 'e' || magicBytes[2] == 0;
		boolean b4 = magicBytes[3] == 'n' || magicBytes[3] == 0;
		return len && b1 && b2 && b3 && b4;
	}

	private long[] readFlags(byte[] read) {
		long compressionFlag;
		long layoutFlag;
		long hasSampFlag;

		boolean b0, b1, b2, b3, b4, b5, b6;

		long ord = BGENBitMath.bytesToLong(read, true);
		b0 = BGENBitMath.getBit(ord, 0);
		b1 = BGENBitMath.getBit(ord, 1);

		compressionFlag = BGENBitMath.bitsToInt(false, b0, b1);

		b2 = BGENBitMath.getBit(ord, 2);
		b3 = BGENBitMath.getBit(ord, 3);
		b4 = BGENBitMath.getBit(ord, 4);
		b5 = BGENBitMath.getBit(ord, 5);
		layoutFlag = BGENBitMath.bitsToInt(false, b2, b3, b4, b5);

		b6 = BGENBitMath.getBit(ord, 31);
		hasSampFlag = b6 ? 1 : 0;

		return new long[] {compressionFlag, layoutFlag, hasSampFlag};
	}

	/**
	 * Read variant map info, skipping gen-/hapl-otype info
	 * 
	 * @throws IOException
	 */
	private void readMap() throws IOException {
		Iterator<BGENRecord> mapIter = new BGENIterators.BGENIterator(this, this.raf, false);
		while (mapIter.hasNext()) {
			mapIter.next();
		}
		reset();
		setupLookup();
	}

	/**
	 * Catalog variant chr/pos
	 */
	private void setupLookup() {
		for (BGENRecordMetaData rec : ptrMap.values()) {
			TreeSet<BGENRecordMetaData> chr = chrSets.get(rec.getChr());
			if (chr == null) {
				chr = new TreeSet<>();
				chrSets.put(rec.getChr(), chr);
			}
			chr.add(rec);
		}
	}

	private BGENRecordMetaData initRecord(long filePointer, RandomAccessFile in, LAYOUT layout)
																																														 throws IOException {
		if (in.getFilePointer() != filePointer) {
			in.seek(filePointer);
		}

		BGENRecordMetaData r = new BGENRecordMetaData();
		r.ptrByt = filePointer;
		int numAll;

		byte[] read = new byte[2]; // 2 bytes to determine length of field
		byte[] readV = new byte[0]; // init to 0 length, will vary based on 2byte values

		if (layout == LAYOUT.V1) {
			byte[] nB = new byte[4];
			in.read(nB);
			r.N = BGENBitMath.bytesToLong(nB, true);
		} else {
			r.N = sampleCount;
		}

		// id
		in.read(read);
		readV = new byte[(int) BGENBitMath.bytesToLong(read, true)];
		in.read(readV);
		r.setId(new String(readV, StandardCharsets.US_ASCII));

		// rs_id
		in.read(read);
		readV = new byte[(int) BGENBitMath.bytesToLong(read, true)];
		in.read(readV);
		r.setRsId(new String(readV, StandardCharsets.US_ASCII));

		// chr
		in.read(read);
		readV = new byte[(int) BGENBitMath.bytesToLong(read, true)];
		in.read(readV);
		r.setChr(Integer.parseInt(new String(readV, StandardCharsets.US_ASCII)));

		// pos
		read = new byte[4];
		in.read(read);
		r.setPos(BGENBitMath.bytesToLong(read, true));

		// number of alleles
		if (layout == LAYOUT.V1) {
			numAll = 2;
		} else if (layout == LAYOUT.V2) {
			read = new byte[2];
			in.read(read);
			numAll = (int) BGENBitMath.bytesToLong(read, true);
		} else {
			throw new UnsupportedOperationException("Only V1 and V2 layouts are currently supported.");
		}

		// alleles
		r.setAlleles(new String[numAll]);
		if (read.length != 4) {
			read = new byte[4];
		}
		for (int i = 0; i < numAll; i++) {
			in.read(read);
			readV = new byte[(int) BGENBitMath.bytesToLong(read, true)];
			in.read(readV);
			r.getAlleles()[i] = new String(readV, StandardCharsets.US_ASCII).intern();
		}

		if (layout == LAYOUT.V2) {
			in.read(read);
			r.blockLength = BGENBitMath.bytesToLong(read, true);
		} else {
			if (compress == COMPRESSION.NONE) {
				r.blockLength = 6 * r.N;
			} else {
				in.read(read);
				r.blockLength = BGENBitMath.bytesToLong(read, true);
			}
		}

		r.lenByt = in.getFilePointer() - filePointer;
		return r;
	}

	BGENRecord readRecord(RandomAccessFile raf, BGENRecordMetaData mapData) throws IOException {
		BGENRecord rec = new BGENRecord(mapData);
		if (raf.getFilePointer() != mapData.ptrByt + mapData.lenByt) {
			raf.seek(mapData.ptrByt + mapData.lenByt);
		}
		switch (layout) {
			case V2:
				rec.data = BGENTools.readLayout2Record(false, raf, mapData, compress);
				break;
			case V1:
				rec.data = BGENTools.readLayout1Record(false, raf, mapData, compress);
				break;
			default:
				throw new UnsupportedOperationException("Layout version " + layout + " is not supported.");
		}
		return rec;
	}

	BGENRecord readNextRecord(RandomAccessFile raf, long location, boolean readFull)
																																									throws IOException {
		BGENRecordMetaData r = getMetaData(raf, location, true);
		return readRecord(raf, r);
	}

	/**
	 * Synchronized on {@code this} if the record metadata is not already cached.
	 * 
	 * @param filePointer Location of variant header to read
	 * @param save Cache metadata if true
	 * @return {@link BGENRecordMetaData}
	 * @throws IOException
	 */
	BGENRecordMetaData getMetaData(RandomAccessFile raf, long filePointer, boolean save)
																																											throws IOException {
		BGENRecordMetaData brmd = ptrMap.get(filePointer);
		if (brmd == null) {
			if (raf.getFilePointer() != filePointer) {
				// allow non-linear reads
				raf.seek(filePointer);
			}
			brmd = initRecord(filePointer, raf, layout);
			if (save) {
				ptrMap.put(filePointer, brmd);
			}
		} else {
			raf.skipBytes((int) brmd.lenByt);
		}
		return brmd;
	}

	private static void checkFlagsOrException(long[] flags) throws IOException {
		if (flags[0] < 0 || flags[0] > 2) {
			throw new IOException("Compression flag must be one of 0, 1, or 2 - found " + flags[0]);
		}
		if (flags[1] < 1) {
			throw new IOException("Only BGEN layout version 1 or greater is supported.");
		}
		if (flags[1] == 1 && flags[0] == 2) {
			throw new IOException("BGEN layout version 1 does not support ZSTD compression.");
		}
	}

	/**
	 * Initialize a new RandomAccessFile object and seek to the start of the first record.
	 * 
	 * @return RandomAccessFile
	 * @throws IOException
	 */
	RandomAccessFile initFilePointer() throws IOException {
		RandomAccessFile raf1 = new RandomAccessFile(filename, "r");
		raf1.seek(recordStart);
		return raf1;
	}

	@Override
	public Iterator<BGENRecord> iterator() {
		return new BGENIterators.BGENIterator(this, true);
	}

	private boolean checkAndReportIfEmpty() {
		if (chrSets.isEmpty()) {
			System.err.println("Map info lookups not built - no records available for query.");
			return false;
		}
		return true;
	}

	/**
	 * Request an Iterable of all records for the specified variant IDs. Will return a
	 * {@link BGENIterators.BGENIterable} wrapped around either a
	 * {@link BGENIterators.BGENRecordIterator} if map info is present, or a
	 * {@link BGENIterators.BGENRecordQueryIterator} if not.
	 * 
	 * @param variantIDs Collections of variant IDs
	 * @return Iterable
	 */
	public Iterable<BGENRecord> query(Collection<String> variantIDs) {
		boolean hasMap = checkAndReportIfEmpty();
		Iterator<BGENRecord> iter = new BGENIterators.EmptyIterator();
		if (!hasMap) {
			iter = new BGENIterators.BGENRecordQueryIterator(this, variantIDs);
		} else {
			List<BGENRecordMetaData> metas = new ArrayList<>();
			for (Entry<Integer, TreeSet<BGENRecordMetaData>> ent : chrSets.entrySet()) {
				for (BGENRecordMetaData rec : ent.getValue()) {
					if (variantIDs.contains(rec.id) || variantIDs.contains(rec.rsId)) {
						metas.add(rec);
					}
				}
			}
			if (!metas.isEmpty()) {
				iter = new BGENIterators.BGENRecordIterator(this, metas);
			}
		}
		return new BGENIterators.BGENIterable(iter);
	}

	/**
	 * @see {@link BGENReader#query(int)}
	 * 
	 * @param chr String chromosome
	 * @return Iterable
	 * @throws NumberFormatException if {@code chr} is not a valid number
	 */
	public Iterable<BGENRecord> query(String chr) {
		return query(Integer.parseInt(chr));
	}

	/**
	 * Request an iterator of all records on the given chromosome. Will return a
	 * {@link BGENIterators.BGENIterable} wrapped around either a
	 * {@link BGENIterators.BGENRecordIterator} if map info is present, or a
	 * {@link BGENIterators.BGENRegionQueryIterator} if not.
	 * 
	 * @param chr Chromosome to query
	 * @return Iterable
	 */
	public Iterable<BGENRecord> query(int chr) {
		boolean hasMap = checkAndReportIfEmpty();
		Iterator<BGENRecord> iter = new BGENIterators.EmptyIterator();
		if (!hasMap) {
			HashMap<Integer, List<int[]>> rgns = new HashMap<>();
			rgns.put(chr, null);
			iter = new BGENIterators.BGENRegionQueryIterator(this, rgns);
		} else {
			Set<BGENRecordMetaData> chrData;
			if ((chrData = chrSets.get(chr)) == null) {
				iter = new BGENIterators.EmptyIterator();
			} else {
				iter = new BGENIterators.BGENRecordIterator(this, chrData);
			}
		}
		return new BGENIterators.BGENIterable(iter);
	}

	/**
	 * @see {@link BGENReader#query(int, int, int)}
	 * @param chr Chromosome to query
	 * @param start Region start (inclusive)
	 * @param stop Region stop (inclusive)
	 * @return Iterable
	 * @throws NumberFormatException if {@code chr} is not a valid number
	 */
	public Iterable<BGENRecord> query(String chr, int start, int stop) {
		return query(Integer.parseInt(chr), start, stop);
	}

	/**
	 * Request an iterator of all records within the specified regions. Will return a
	 * {@link BGENIterators.BGENIterable} wrapped around either a
	 * {@link BGENIterators.BGENRecordIterator} if map info is present, or a
	 * {@link BGENIterators.BGENRegionQueryIterator} if not.
	 * 
	 * @param chr Chromosome to query
	 * @param start Region start (inclusive)
	 * @param stop Region stop (inclusive)
	 * @return Iterable
	 */
	public Iterable<BGENRecord> query(Map<Integer, List<int[]>> regions) {
		boolean hasMap = checkAndReportIfEmpty();
		Iterator<BGENRecord> iter = new BGENIterators.EmptyIterator();
		if (!hasMap) {
			iter = new BGENIterators.BGENRegionQueryIterator(this, regions);
		} else {
			Set<BGENRecordMetaData> chrData = new HashSet<>();
			List<BGENRecordMetaData> chrList = new ArrayList<>();
			Set<Integer> chrs = regions.keySet();
			for (Integer chr : chrs) {
				if (chrSets.containsKey(chr)) {
					chrData.addAll(chrSets.get(chr));
				}
			}
			if (chrData.isEmpty()) {
				iter = new BGENIterators.EmptyIterator();
			} else {
				for (BGENRecordMetaData rec : chrData) {
					if (!chrs.contains(rec.chr))
						continue;
					for (int[] rgn : regions.get(rec.chr)) {
						if (rec.getPos() >= rgn[0] || rec.getPos() <= rgn[1]) {
							chrList.add(rec);
						}
					}
				}
				iter = new BGENIterators.BGENRecordIterator(this, chrList);
			}
		}
		return new BGENIterators.BGENIterable(iter);
	}

	/**
	 * Request an iterator of all records on the given chromosome within a specified region. Will
	 * return a {@link BGENIterators.BGENIterable} wrapped around either a
	 * {@link BGENIterators.BGENRecordIterator} if map info is present, or a
	 * {@link BGENIterators.BGENRegionQueryIterator} if not.
	 * 
	 * @param chr Chromosome to query
	 * @param start Region start (inclusive)
	 * @param stop Region stop (inclusive)
	 * @return Iterable
	 */
	public Iterable<BGENRecord> query(int chr, int start, int stop) {
		boolean hasMap = checkAndReportIfEmpty();
		Iterator<BGENRecord> iter = new BGENIterators.EmptyIterator();
		if (!hasMap) {
			HashMap<Integer, List<int[]>> rgns = new HashMap<>();
			ArrayList<int[]> rgn = new ArrayList<>();
			rgn.add(new int[] {start, stop});
			rgns.put(chr, rgn);
			iter = new BGENIterators.BGENRegionQueryIterator(this, rgns);
		} else {
			Set<BGENRecordMetaData> chrData;
			List<BGENRecordMetaData> chrList = new ArrayList<>();
			if ((chrData = chrSets.get(chr)) == null) {
				iter = new BGENIterators.EmptyIterator();
			} else {
				for (BGENRecordMetaData rec : chrData) {
					if (rec.getPos() >= start || rec.getPos() <= stop) {
						chrList.add(rec);
					}
				}
				iter = new BGENIterators.BGENRecordIterator(this, chrList);
			}
		}
		return new BGENIterators.BGENIterable(iter);
	}

	/**
	 * @see {@link BGENReader#query(int, List<int[]>)}
	 * @param chr Chromosome to query
	 * @param regions List of int arrays of {@code start1, stop1}, {start2, stop2}, ... , {startN,
	 *        stopN}}} (start and stop both inclusive)
	 * @return Iterable
	 * @throws NumberFormatException if {@code chr} is not a valid number
	 */
	public Iterable<BGENRecord> query(String chr, List<int[]> regions) {
		return query(Integer.parseInt(chr), regions);
	}


	/**
	 * Request an iterator of all records on the given chromosome within the specified regions. Will
	 * return a {@link BGENIterators.BGENIterable} wrapped around either a
	 * {@link BGENIterators.BGENRecordIterator} if map info is present, or a
	 * {@link BGENIterators.BGENRegionQueryIterator} if not.
	 * 
	 * @param chr Chromosome to query
	 * @param regions int arrays of {@code start1, stop1}, {start2, stop2}, ... , {startN, stopN}}}
	 *        (start and stop both inclusive)
	 * @return Iterable
	 * @throws NumberFormatException if {@code chr} is not a valid number
	 */
	public Iterable<BGENRecord> query(int chr, List<int[]> regions) {
		boolean hasMap = checkAndReportIfEmpty();
		Iterator<BGENRecord> iter = new BGENIterators.EmptyIterator();
		if (!hasMap) {
			HashMap<Integer, List<int[]>> rgns = new HashMap<>();
			rgns.put(chr, regions);
			iter = new BGENIterators.BGENRegionQueryIterator(this, rgns);
		} else {
			Set<BGENRecordMetaData> chrData;
			List<BGENRecordMetaData> chrList = new ArrayList<>();
			if ((chrData = chrSets.get(chr)) == null) {
				iter = new BGENIterators.EmptyIterator();
			} else {
				for (BGENRecordMetaData rec : chrData) {
					for (int[] rgn : regions) {
						if (rec.getPos() >= rgn[0] || rec.getPos() <= rgn[1]) {
							chrList.add(rec);
						}
					}
				}
				iter = new BGENIterators.BGENRecordIterator(this, chrList);
			}
		}
		return new BGENIterators.BGENIterable(iter);
	}

	@Override
	public void close() throws IOException {
		closed = true;
		raf.close();
	}

	/**
	 * Reset the reader to the beginning of the first BGEN record.<br />
	 * No-op if closed.
	 * 
	 * @throws IOException
	 */
	public void reset() throws IOException {
		if (closed) {
			return;
		}
		raf.seek(getRecordStartByte());
	}

	public long getRecordCount() {
		return recordCount;
	}

	public long getSampleCount() {
		return sampleCount;
	}

	public COMPRESSION getCompression() {
		return compress;
	}

	public LAYOUT getLayout() {
		return layout;
	}

	public boolean hasSampleNames() {
		return hasSampleNames;
	}

	/**
	 * @return Array of String sample names, or {@code null} if {@link #hasSampleNames()} is
	 *         {@code false}.
	 */
	public String[] getSamples() {
		return samples;
	}

	public long getRecordStartByte() {
		return recordStart;
	}

	public void setRecordStart(long recordStart) {
		this.recordStart = recordStart;
	}

	public enum COMPRESSION {
		NONE,
		ZLIB,
		ZSTD;
	}

	public enum LAYOUT {
		V1,
		V2;
	}

	public final static class BGENRecordMetaData implements Comparable<BGENRecordMetaData> {
		String id;
		String rsId;
		private int chr;
		private long pos;
		private String[] alleles;

		long N;
		long ptrByt;
		long lenByt;
		long blockLength;

		@Override
		public int compareTo(BGENRecordMetaData o) {
			int comp = Integer.compare(getChr(), o.getChr());
			if (comp != 0)
				return comp;
			comp = Long.compare(getPos(), o.getPos());
			if (comp != 0)
				return comp;
			comp = getRsID().compareTo(o.getRsID());
			if (comp != 0)
				return comp;
			comp = getID().compareTo(o.getID());
			return comp;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (int) (N ^ (N >>> 32));
			result = prime * result + Arrays.hashCode(getAlleles());
			result = prime * result + getChr();
			result = prime * result + ((getID() == null) ? 0 : getID().hashCode());
			result = prime * result + (int) (lenByt ^ (lenByt >>> 32));
			result = prime * result + (int) (getPos() ^ (getPos() >>> 32));
			result = prime * result + (int) (ptrByt ^ (ptrByt >>> 32));
			result = prime * result + ((getRsID() == null) ? 0 : getRsID().hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			BGENRecordMetaData other = (BGENRecordMetaData) obj;
			if (N != other.N)
				return false;
			if (!Arrays.equals(getAlleles(), other.getAlleles()))
				return false;
			if (getChr() != other.getChr())
				return false;
			if (getID() == null) {
				if (other.getID() != null)
					return false;
			} else if (!getID().equals(other.getID()))
				return false;
			if (lenByt != other.lenByt)
				return false;
			if (getPos() != other.getPos())
				return false;
			if (ptrByt != other.ptrByt)
				return false;
			if (getRsID() == null) {
				if (other.getRsID() != null)
					return false;
			} else if (!getRsID().equals(other.getRsID()))
				return false;
			return true;
		}

		public String getRsID() {
			return rsId;
		}

		public String getID() {
			return id;
		}

		public String[] getAlleles() {
			return alleles;
		}

		public long getPos() {
			return pos;
		}

		public int getChr() {
			return chr;
		}

		private void setRsId(String rsId) {
			this.rsId = rsId;
		}

		private void setId(String id) {
			this.id = id;
		}

		private void setAlleles(String[] alleles) {
			this.alleles = alleles;
		}

		private void setPos(long pos) {
			this.pos = pos;
		}

		private void setChr(int chr) {
			this.chr = chr;
		}

	}

	public final static class BGENRecord {
		public BGENRecord(BGENRecordMetaData metaData) {
			this.metaData = metaData;
		}

		private final BGENRecordMetaData metaData;
		private double[][] data;

		public BGENRecordMetaData getMetaData() {
			return metaData;
		}

		public double[][] getData() {
			return data;
		}
	}

}
