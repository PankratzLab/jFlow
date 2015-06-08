package filesys.rao;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.zip.GZIPOutputStream;

import common.Logger;

/**
 * Class to write java objecs in a random manner
 *
 */
public class RAOWriter {
	public static String RAO_EXT = ".rao";
	public static String RAO_INDEX_EXT = ".idx";
	private String fullPathToFile;
	private long offset;
	private RandomAccessProducer rProducer;
	private RandomAccessFile file;
	private FileOutputStream fos;
	private RAOIndex index;
	private OutputStream os;
	private Logger log;

	/**
	 * @param fullPathToFile
	 *            where the {@link RAObject } objects will be stored, must end with {@link RAOWriter#RAO_EXT}
	 * @param rProducer
	 *            dishes up {@link RAObject } objects to be serialized, compressed, and written
	 * @param log
	 */
	public RAOWriter(String fullPathToFile, RandomAccessProducer rProducer, Logger log) {
		super();
		this.fullPathToFile = fullPathToFile;
		this.index = getIndex();
		this.rProducer = rProducer;
		this.log = log;
		this.offset = 0;
	}

	public RAOWriter(String fullPathToFile, Logger log) throws IOException {
		this.fullPathToFile = fullPathToFile;
		this.file = new RandomAccessFile(fullPathToFile, "rw");
		this.os = new FileOutputStream(file.getFD());
		// this.os = new GZIPOutputStream(fos);
		this.index = getIndex();
		this.log = log;
		this.offset = 0;
	}

	public void close() throws IOException {
		os.close();
		// fos.close();
		file.close();
		index.serialize();
	}

	public void add(RAObject raObject) throws IOException {
		offset += addObject(raObject);
	}

	private RAOIndex getIndex() {
		String indexFile = fullPathToFile + RAO_INDEX_EXT;
		return new RAOIndex(indexFile, new Hashtable<String, ArrayList<Long>>());
	}

	public WriteComplete writeToFile() {
		WriteComplete writeComplete = new WriteComplete(fullPathToFile, fullPathToFile + RAO_INDEX_EXT, false);
		if (!fullPathToFile.endsWith(RAO_EXT)) {
			String error = "Invalid file extension, " + fullPathToFile + " must end with " + RAO_EXT;
			log.reportTimeError(error);
		} else {
			try {
				RandomAccessFile file = new RandomAccessFile(fullPathToFile, "rw");
				FileOutputStream fos = new FileOutputStream(file.getFD());
				GZIPOutputStream gos = new GZIPOutputStream(fos);
				index = writeToFile(gos, index);
				index.serialize();
				fos.close();
				file.close();
				writeComplete.setWritten(true);
			} catch (FileNotFoundException e) {
				log.reportFileNotFound(fullPathToFile);
				log.reportException(e);
				e.printStackTrace();
			} catch (IOException e) {
				log.reportException(e);
				e.printStackTrace();
			}

		}
		return writeComplete;
	}

	private RAOIndex writeToFile(OutputStream fos, RAOIndex index) throws IOException {
		while (rProducer.hasNext()) {
			offset += addObject(rProducer.next());
		}
		return index;
	}

	private long addObject(RAObject rObject) throws IOException {
		String[] indexKeys = rObject.getIndexKeys();
		for (int i = 0; i < indexKeys.length; i++) {
			if (!index.getIndex().containsKey(indexKeys[i])) {
				index.getIndex().put(indexKeys[i], new ArrayList<Long>());
			}
		}
		ByteArrayOutputStream bos = RAOExt.convertAndCompress(rObject);
		os.write(bos.toByteArray());
		bos.close();
		long size = (long) bos.size();
		if (size > index.getMaxSize()) {
			index.setMaxSize(index.getMaxSize() + size + 1);
		}
		if (offset % 10000 == 0) {
			System.out.println(offset + "\t" + size + "\t" + bos.size());
		}

		for (int i = 0; i < indexKeys.length; i++) {
			index.getIndex().get(indexKeys[i]).add(offset);
		}
		return size;
	}

	public interface RandomAccessProducer extends Iterator<RAObject> {

	}

	/**
	 * Stores the write status from the {@link RAOWriter#writeToFile()} method
	 *
	 */
	public static class WriteComplete {
		private String dataFile;
		private String indexFile;
		private boolean written;

		public WriteComplete(String dataFile, String indexFile, boolean written) {
			super();
			this.dataFile = dataFile;
			this.indexFile = indexFile;
			this.written = written;
		}

		public String getDataFile() {
			return dataFile;
		}

		public void setDataFile(String dataFile) {
			this.dataFile = dataFile;
		}

		public String getIndexFile() {
			return indexFile;
		}

		public void setIndexFile(String indexFile) {
			this.indexFile = indexFile;
		}

		public boolean getWritten() {
			return written;
		}

		public void setWritten(boolean written) {
			this.written = written;
		}

	}

}
