package common;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Hashtable;

/**
 * Author: Rohit Sinha
 * Class to maintain a Boolean Map of Chr Positions.
 * The outer Hashtable houses the data for each of the chromosomes (Byte) separately. The inner Hashtable flags each position (Integer) as being inside of a DHS region (Boolean, true/false)
 */
public class ChrPositionBooleanMap implements Serializable {

	public static final long serialVersionUID = 1L;

	Hashtable<Byte, HashSet<Integer>> chrPositionMap = new Hashtable<Byte, HashSet<Integer>>();

	public Hashtable<Byte, Hashtable<Integer, Boolean>> getChrPositionMap() {
		return chrPositionMap;
	}

	public ChrPositionBooleanMap(ChrPositionBooleanMap chrPositionMapObject) {
		this.chrPositionMap = chrPositionMapObject.chrPositionMap;
	}

	public void setChrPositionMap(Hashtable<Byte, Hashtable<Integer, Boolean>> chrPositionMap) {
		this.chrPositionMap = chrPositionMap;
	}

	public void writeToFile(String filepath) throws IOException {

		ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(filepath));
		out.writeObject(this);
		out.close();
	}

	public void readFromFile(String filepath) throws IOException, ClassNotFoundException {
		FileInputStream in = new FileInputStream(filepath);
		ObjectInputStream reader = new ObjectInputStream(in);
		this.setChrPositionMap(((ChrPositionBooleanMap) reader.readObject()).chrPositionMap);
		reader.close();
	}
}
