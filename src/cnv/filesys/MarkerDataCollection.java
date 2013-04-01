package cnv.filesys;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.util.Date;
import java.util.Hashtable;

import common.Files;

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

	/**
	 * So far only TransposeDate.extractForDemo(...) needs to save MarkerDataCollection, and serialize() works better for it.
	 * So, we are not going to write code for saveToRandomAccessFile() and loadFromRandomAccessFile();
	 */
	public static void saveToRandomAccessFile() {
	}

	/**
	 * So far only TransposeDate.extractForDemo(...) needs to save MarkerDataCollection, and serialize() works better for it.
	 * So, we are not going to write code for saveToRandomAccessFile() and loadFromRandomAccessFile();
	 */
	public static void loadFromRandomAccessFile(String filename, boolean jar) {
	}
}
