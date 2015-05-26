package filesys.rao;

import java.io.Externalizable;

public interface RAObject extends Externalizable {
	/**
	 * @return the keys that will be used for the indexed retrieval of data
	 */
	public String[] getIndexKeys();
}
