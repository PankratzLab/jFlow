package cnv.filesys;

import java.io.Serializable;
import java.util.Hashtable;

import common.Files;

public class MarkerLookup implements Serializable {
	public static final long serialVersionUID = 1L;

	private Hashtable<String,String> hash;

	public MarkerLookup(Hashtable<String,String> hash) {
		this.hash = hash;
	}

	public boolean contains(String markerName) {
		return hash.containsKey(markerName);
	}

	public String get(String markerName) {
		return hash.get(markerName);
	}

	public int getSize() {
		return hash.size();
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static MarkerLookup load(String filename, boolean jar) {
		return (MarkerLookup)Files.readSerial(filename, jar, true);
	}
}
