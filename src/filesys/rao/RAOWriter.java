package filesys.rao;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
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
	private RandomAccessProducer rProducer;
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
		this.rProducer = rProducer;
		this.log = log;
	}

	public WriteComplete writeToFile() {
		WriteComplete writeComplete = new WriteComplete(fullPathToFile, fullPathToFile + RAO_INDEX_EXT, false);
		if (!fullPathToFile.endsWith(RAO_EXT)) {
			String error = "Invalid file extension, " + fullPathToFile + " must end with " + RAO_EXT;
			log.reportTimeError(error);
		} else {
			try {
				String indexFile = fullPathToFile + RAO_INDEX_EXT;
				RandomAccessFile file = new RandomAccessFile(fullPathToFile, "rw");
				FileOutputStream fos = new FileOutputStream(file.getFD());
				GZIPOutputStream gos = new GZIPOutputStream(fos);
				BlockCompressedOutputStream bcos = new BlockCompressedOutputStream(fullPathToFile, 9);
				RAOIndex index = new RAOIndex(indexFile, new Hashtable<String, ArrayList<Long>>());
				index = writeToFile(bcos, index);
				index.serialize();
				// gos.close();
				bcos.close();
				System.exit(1);
				fos.close();
				file.close();
				BlockCompressedInputStream bcis = new BlockCompressedInputStream(new File(fullPathToFile));
				//bcis.getFileBlock(bgzfOffset);
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

	private RAOIndex writeToFile(BlockCompressedOutputStream fos, RAOIndex index) throws IOException {
		long offset = 0;
		long maxSize = 0;

		while (rProducer.hasNext()) {
			RAObject rObject = rProducer.next();
			String[] indexKeys = rObject.getIndexKeys();
			for (int i = 0; i < indexKeys.length; i++) {
				if (!index.getIndex().containsKey(indexKeys[i])) {
					index.getIndex().put(indexKeys[i], new ArrayList<Long>());

				}
			}
			// ObjectOutputStream oos = new ObjectOutputStream(bos);
			// oos.writeObject(rObject);
			// oos.flush();
			// ByteArrayOutputStream byteOutput = new ByteArrayOutputStream();
			// GZIPOutputStream gzipOutput = new GZIPOutputStream(byteOutput);
			// gzipOutput.write(bos.toByteArray());
			ByteArrayOutputStream bos = RAOExt.convertAndCompress(rObject);
			fos.write(bos.toByteArray());
			//bos.writeTo(fos);
			bos.close();
			// for (int i = 0; i < bos.toByteArray().length; i++) {
			// System.out.println(bos.toByteArray()[i]);
			// }
			long size = (long) bos.size();
			if (size > maxSize) {
				maxSize = size + 1;
			}
			if (offset % 10000 == 0) {
				System.out.println(offset + "\t" + size + "\t" + bos.size());
			}

			for (int i = 0; i < indexKeys.length; i++) {
				index.getIndex().get(indexKeys[i]).add(offset);
			}
			offset += size;
		}
		index.setMaxSize(maxSize);
		return index;
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
