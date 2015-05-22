package filesys.rao;

import java.util.ArrayList;
import java.util.Hashtable;

import common.Files;
import common.Logger;

/**
 * Serialized index for accessing java objects
 *
 */
public class RAOIndex implements RAObject {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String indexFileName;
	private Hashtable<String, ArrayList<Long>> index;
	private long maxSize;

	public String getFileName() {
		return indexFileName;
	}

	public RAOIndex(String fileName, Hashtable<String, ArrayList<Long>> index) {
		super();
		this.indexFileName = fileName;
		this.index = index;
		this.maxSize = maxSize;
	}

	public void setFileName(String fileName) {
		this.indexFileName = fileName;
	}

	public long getMaxSize() {
		return maxSize;
	}

	public void setMaxSize(long maxSize) {
		this.maxSize = maxSize;
	}

	public Hashtable<String, ArrayList<Long>> getIndex() {
		return index;
	}

	public void setIndex(Hashtable<String, ArrayList<Long>> index) {
		this.index = index;
	}

	@Override
	public String[] getIndexKeys() {
		return new String[] { indexFileName };
	}

	public void serialize() {
		Files.writeSerial(this, indexFileName);
		//TODO change this
	}

	public static RAOIndex load(String indexFileName, Logger log) {
		return (RAOIndex) Files.readSerial(indexFileName, false, log, false);
	}

}
