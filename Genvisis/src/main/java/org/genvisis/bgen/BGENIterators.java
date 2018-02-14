package org.genvisis.bgen;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.genvisis.bgen.BGENReader.BGENRecord;
import org.genvisis.bgen.BGENReader.BGENRecordMetaData;

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

  static abstract class BGENRAFIterator implements Iterator<BGENReader.BGENRecord> {

    final BGENReader reader;
    final RandomAccessFile raf;

    public BGENRAFIterator(BGENReader reader, RandomAccessFile raf) {
      this.reader = reader;
      this.raf = raf;
    }

    public BGENRAFIterator(BGENReader reader) {
      this.reader = reader;
      try {
        this.raf = reader.initFilePointer();
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }

  }

  /**
   * Iterator over a specific set of BGENRecords
   */
  static final class BGENRecordIterator extends BGENRAFIterator {

    private Iterator<BGENRecordMetaData> iterator;

    BGENRecordIterator(BGENReader reader, Collection<BGENReader.BGENRecordMetaData> records) {
      super(reader);
      this.iterator = records.iterator();
    }

    @Override
    public boolean hasNext() {
      boolean has = iterator.hasNext();
      if (!has) {
        try {
          this.raf.close();
        } catch (IOException e) {
          throw new RuntimeException(e);
        }
      }
      return has;
    }

    @Override
    public BGENReader.BGENRecord next() {
      try {
        return reader.readRecord(this.raf, iterator.next(), false);
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }

  }

  /**
   * BGENQueryIterator searches for the next record based on the subclass implementation of
   * {@code findNextRecord()}.<br />
   * The {@code hasNext()} method will scan forward through metadata blocks searching for the next
   * variant record. The {@code next()} method then reads the full variant information. <br/>
   * This implies a few things:
   * <ul>
   * <li>(a) as BGEN files are not ordered, {@code hasNext()} will only return false when the end of
   * the file is reached,</li>
   * <li>(b) {@code hasNext()} may block for a indeterminate amount of time as it scans for the next
   * valid variant record,</li>
   * <li>and (c) the {@code next()} method itself is fast, as it already has the next variant record
   * ready.</li>
   * </ul>
   */
  static abstract class BGENQueryIterator extends BGENRAFIterator {

    private BGENReader.BGENRecordMetaData next = null;
    private boolean end = false;

    /**
     * @see {@link BGENQueryIterator}
     * @param reader
     * @param regions
     */
    BGENQueryIterator(BGENReader reader) {
      super(reader);
    }

    @Override
    public boolean hasNext() {
      if (end) return false;
      if (next == null) {
        try {
          next = findNextRecord();
        } catch (IOException e) {
          throw new RuntimeException(e);
        }
      }
      if (next == null) {
        end = true;
        if (end) {
          try {
            this.raf.close();
          } catch (IOException e) {
            throw new RuntimeException(e);
          }
        }
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
        BGENRecord rec = hasNext() ? reader.readRecord(this.raf, next, false) : null;
        next = null;
        return rec;
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }

  }

  /**
   * BGENQueryIterator subclass for extracting specific records based on their IDs
   * 
   * @see {@link BGENQueryIterator}
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
        if (this.raf.getFilePointer() == this.raf.length()) return null;
        r = this.reader.getMetaData(this.raf, this.raf.getFilePointer(), false);
        inVars = (variantIDs.contains(r.id) || variantIDs.contains(r.rsId));
        if (!inVars) {
          this.raf.skipBytes((int) r.blockLength);
        }
      } while (!inVars);
      if (!inVars) return null;
      return r;
    }

  }

  /**
   * BGENQueryIterator subclass for extracting records within regions.
   * 
   * @see {@link BGENQueryIterator}
   */
  static final class BGENRegionQueryIterator extends BGENQueryIterator {

    private final Map<Integer, List<int[]>> regions;

    BGENRegionQueryIterator(BGENReader reader, Map<Integer, List<int[]>> regions) {
      super(reader);
      this.regions = regions;
    }

    protected BGENReader.BGENRecordMetaData findNextRecord() throws IOException {
      BGENReader.BGENRecordMetaData r = null;
      boolean inRegion = false;
      do {
        if (this.raf.getFilePointer() == this.raf.length()) return null;
        r = this.reader.getMetaData(this.raf, this.raf.getFilePointer(), false);
        inRegion = regions.containsKey(r.getChr());
        if (inRegion && regions.get(r.getChr()) != null) {
          inRegion = false;
          for (int[] rgn : regions.get(r.getChr())) {
            if (r.getPos() >= rgn[0] && r.getPos() <= rgn[1]) {
              inRegion = true;
              break;
            }
          }
          if (!inRegion) {
            break;
          }
        }
        if (!inRegion) {
          this.raf.skipBytes((int) r.blockLength);
        }
      } while (!inRegion);
      if (!inRegion) return null;
      return r;
    }
  }

  /**
   * Iterate over all records in the given BGENReader, reading either just map info or the full
   * record.
   */
  static final class BGENIterator extends BGENRAFIterator {

    private final long N;
    private long start = 0;
    private long read = 0;

    /**
     * Read data if true or metadata only if false
     */
    private boolean readFully = true;
    private boolean saveMeta = true;
    private boolean closeWhenDone = true;

    /**
     * @see {@link BGENIterator}
     * @param reader
     * @param readInFull Read data if true or metadata only if false
     */
    BGENIterator(BGENReader reader, boolean readInFull) {
      super(reader);
      this.N = reader.getRecordCount();
      this.readFully = readInFull;
    }

    public BGENIterator(BGENReader reader, RandomAccessFile raf, boolean readInFull,
                        boolean saveMetaData, boolean closeWhenDone) {
      this(reader, raf, readInFull, saveMetaData, closeWhenDone, 0);
    }

    public BGENIterator(BGENReader reader, RandomAccessFile raf, boolean readInFull,
                        boolean saveMetaData, boolean closeWhenDone, long start) {
      super(reader, raf);
      this.N = reader.getRecordCount();
      this.readFully = readInFull;
      this.saveMeta = saveMetaData;
      this.closeWhenDone = closeWhenDone;
      this.start = start;
      skipAhead();
    }

    private void skipAhead() {
      while (read < start) {
        read++;
        try {
          reader.readNextRecord(this.raf, this.raf.getFilePointer(), false, false);
        } catch (IOException e) {
          throw new RuntimeException(e);
        }
      }
    }

    @Override
    public boolean hasNext() {
      return read < N;
    }

    @Override
    public BGENReader.BGENRecord next() {
      try {
        read++;
        BGENRecord rec = reader.readNextRecord(this.raf, this.raf.getFilePointer(), readFully,
                                               saveMeta);
        return rec;
      } catch (IOException e) {
        throw new RuntimeException(e);
      } finally {
        if (read == N && closeWhenDone) {
          try {
            this.raf.close();
          } catch (IOException e) {
            throw new RuntimeException(e);
          }
        }
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
