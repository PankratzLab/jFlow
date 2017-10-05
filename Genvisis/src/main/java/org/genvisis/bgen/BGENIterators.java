package org.genvisis.bgen;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public final class BGENIterators {

	private BGENIterators() {}

	static final class EmptyIterator implements Iterator<BGENReader.BGENRecord> {
		@Override
		public boolean hasNext() {
			return false;
		}

		@Override
		public BGENReader.BGENRecord next() {
			return null;
		}

		@Override
		public void remove() {}

	}

	/**
	 * Iterator over a specific set of BGENRecords
	 * 
	 */
	static final class BGENRecordIterator implements Iterator<BGENReader.BGENRecord> {
		final BGENReader reader;
		final List<BGENReader.BGENRecordMetaData> records;
		int read = 0;

		BGENRecordIterator(BGENReader reader, Collection<BGENReader.BGENRecordMetaData> records) {
			this.reader = reader;
			this.records = records instanceof List ? (List<BGENReader.BGENRecordMetaData>) records
																						: new ArrayList<>(records);
		}

		@Override
		public boolean hasNext() {
			return read < records.size();
		}

		@Override
		public BGENReader.BGENRecord next() {
			try {
				return reader.readRecord(records.get(read++));
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

	}

	/**
	 * BGENQueryIterator searches for the next record based on the subclass implementation of
	 * <code>findNextRecord()</code>.<br />
	 * The <code>hasNext()</code> method will scan forward through metadata blocks searching for the
	 * next variant record. The <code>next()</code> method then reads the full variant information.
	 * This implies a few things: <br/>
	 * <ul>
	 * <li>(a) as BGEN files are not ordered, <code>hasNext()</code> will only return false when the
	 * end of the file is reached,</li>
	 * <li>
	 * (b) <code>hasNext()</code> may block for a indeterminate amount of time as it scans for the
	 * next valid variant record,</li>
	 * <li>and (c) the <code>next()</code> method itself is fast, as it already has the next variant
	 * record ready.</li>
	 * </ul>
	 * 
	 */
	static abstract class BGENQueryIterator implements Iterator<BGENReader.BGENRecord> {
		protected final BGENReader reader;
		private BGENReader.BGENRecordMetaData next = null;
		private boolean end = false;

		/**
		 * @see {@link BGENQueryIterator}
		 *
		 * @param reader
		 * @param regions
		 */
		BGENQueryIterator(BGENReader reader) {
			this.reader = reader;
			try {
				this.reader.reset();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

		@Override
		public boolean hasNext() {
			if (end)
				return false;
			try {
				next = findNextRecord();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			if (next == null) {
				end = true;
			}
			return next != null;
		}

		/**
		 * Search for the next record to load, or return null if the end of the file is reached
		 * 
		 * @return
		 * @throws IOException
		 */
		protected abstract BGENReader.BGENRecordMetaData findNextRecord() throws IOException;

		@Override
		public BGENReader.BGENRecord next() {
			try {
				return next == null ? null : reader.readRecord(next);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

	}

	/**
	 * BGENQueryIterator subclass for extracting specific records based on their IDs
	 * 
	 * @see {@link BGENQueryIterator}
	 *
	 */
	static final class BGENRecordQueryIterator extends BGENQueryIterator {
		private final Collection<String> variantIDs;

		BGENRecordQueryIterator(BGENReader reader, Collection<String> variantIDs) {
			super(reader);
			this.variantIDs = variantIDs;
		}

		@Override
		protected BGENReader.BGENRecordMetaData findNextRecord() throws IOException {
			BGENReader.BGENRecordMetaData r = null;
			boolean inVars = false;
			do {
				r = this.reader.getMetaData(this.reader.getRAF().getFilePointer(), false);
				inVars = (variantIDs.contains(r.id) || variantIDs.contains(r.rsId));
			} while (!inVars);
			if (!inVars)
				return null;
			return r;
		}

	}

	/**
	 * BGENQueryIterator subclass for extracting records within regions.
	 * 
	 * @see {@link BGENQueryIterator}
	 *
	 */
	static final class BGENRegionQueryIterator extends BGENQueryIterator {
		private final Map<Integer, int[][]> regions;

		BGENRegionQueryIterator(BGENReader reader, Map<Integer, int[][]> regions) {
			super(reader);
			this.regions = regions;
		}

		protected BGENReader.BGENRecordMetaData findNextRecord() throws IOException {
			BGENReader.BGENRecordMetaData r = null;
			boolean inRegion = false;
			do {
				r = this.reader.getMetaData(this.reader.getRAF().getFilePointer(), false);
				inRegion = regions.keySet().contains(r.getChr());
				for (int[] rgn : regions.get(r.getChr())) {
					if (r.getPos() >= rgn[0] && r.getPos() <= rgn[1]) {
						inRegion = true;
						break;
					}
				}
			} while (!inRegion);
			if (!inRegion)
				return null;
			return r;
		}
	}

	/**
	 * Iterate over all records in the given BGENReader, reading either just map info or the full
	 * record.
	 * 
	 */
	static final class BGENIterator implements Iterator<BGENReader.BGENRecord> {
		private final BGENReader reader;
		private final long N;
		private long read = 0;
		/**
		 * Read data if true or metadata only if false
		 */
		private boolean readFully = true;

		/**
		 * @see {@link BGENIterator}
		 * 
		 * @param reader
		 * @param readInFull Read data if true or metadata only if false
		 */
		BGENIterator(BGENReader reader, boolean readInFull) {
			this.reader = reader;
			this.N = reader.getRecordCount();
			this.readFully = readInFull;
		}

		@Override
		public boolean hasNext() {
			return read < N;
		}

		@Override
		public BGENReader.BGENRecord next() {
			try {
				read++;
				return reader.readNextRecord(readFully);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

	}

	/**
	 * Wrapper
	 */
	static class BGENIterable implements Iterable<BGENReader.BGENRecord> {
		final Iterator<BGENReader.BGENRecord> iter;

		public BGENIterable(Iterator<BGENReader.BGENRecord> iter) {
			this.iter = iter;
		}

		@Override
		public Iterator<BGENReader.BGENRecord> iterator() {
			return this.iter;
		}
	}

}
