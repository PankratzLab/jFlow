package org.genvisis.cnv.filesys;

import java.io.Serializable;
import java.util.Hashtable;

import org.genvisis.common.HashVec;
import org.genvisis.common.SerializedFiles;

public class MarkerLookup implements Serializable {
	public static final long serialVersionUID = 1L;

	private final Hashtable<String, String> hash;

	public MarkerLookup(Hashtable<String, String> hash) {
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
		SerializedFiles.writeSerial(this, filename);
	}

	public static MarkerLookup load(String filename, boolean jar) {
		return (MarkerLookup) SerializedFiles.readSerial(filename, jar, true);
	}

	public String[] getMarkerList() {
		return HashVec.getKeys(hash, false, false);
	}

	public String getFirstMarkerDataRafFilename() {
		return hash.elements().nextElement().split("\t")[0];
	}

	public String[] getMarkerDataRafFilenames() {
		Hashtable<String, String> filenames;
		String[] line;
		String[] listOfMarkersInMarkerLookup;

		filenames = new Hashtable<String, String>();

		listOfMarkersInMarkerLookup = getMarkerList();
		for (String element : listOfMarkersInMarkerLookup) {
			line = get(element).split("[\\s]+");
			if (!filenames.containsKey(line[0])) {
				filenames.put(line[0], "");
			}
		}

		return HashVec.getKeys(filenames, false, false);
	}

}
