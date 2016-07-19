package filesys.rao;

import java.io.Serializable;

public interface RAObject extends Serializable {
	/**
	 * @return the keys that will be used for the indexed retrieval of data
	 */
	public String[] getIndexKeys();
}
