package cnv.filesys;

import java.io.Serializable;
import common.Files;

// TODO phase out or retain for demos?
public class MarkerDataCollection implements Serializable {
	public static final long serialVersionUID = 1L;

	private String filename;
	private long fingerprint;
	private String[] markerNames;
	private MarkerData[] collection;

	public MarkerDataCollection(String filename, long fingerprint, String[] markerNames, MarkerData[] collection) {
		this.filename = filename;
		this.fingerprint = fingerprint;
		this.markerNames = markerNames;
		this.collection = collection;
	}

	public String getFilename() {
		return filename;
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public String[] getMarkerNames() {
		return markerNames;
	}

	public MarkerData[] getCollection() {
		return collection;
	}

	public void serialize() {
		Files.writeSerial(this, filename);
	}

	public static MarkerDataCollection load(String filename, boolean jar) {
		return (MarkerDataCollection)Files.readSerial(filename, jar, true);
	}
}
