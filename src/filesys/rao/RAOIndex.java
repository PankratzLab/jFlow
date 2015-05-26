package filesys.rao;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;
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
		Files.writeSerial(this, indexFileName, true);
	}

	public static RAOIndex load(String indexFileName, Logger log) {
		return (RAOIndex) Files.readSerial(indexFileName, false, log, false, true);
	}

	public RAOIndex() {
	}

	@Override
	public void writeExternal(ObjectOutput out) throws IOException {
		out.writeObject(indexFileName);
		out.writeObject(index);
		out.writeLong(maxSize);
		// TODO Auto-generated method stub

	}

	@SuppressWarnings("unchecked")
	@Override
	public void readExternal(ObjectInput in) throws IOException, ClassNotFoundException {
		this.indexFileName = (String) in.readObject();
		this.index = (Hashtable<String, ArrayList<Long>>) in.readObject();
		this.maxSize = in.readLong();
	}

}
