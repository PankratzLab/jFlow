package cnv.filesys;

import java.io.*;
import java.util.Hashtable;

import common.Elision;
import common.Files;

public class OutlierHashtable implements Serializable {
	private static final long serialVersionUID = 1L;
	public static final String[] DATA_ITEMS = new String[] {};
	public static final int MAX_NUM_OF_SAMPLES_BIG_HASH = 1000000;
	public static final int MAX_NUM_OF_SAMPLES_SMALL_HASH = 10000;	//Need to make this dynamic

	private long sampFingerPrintFromProj;
	private long markFingerPrintFromProj;
	Hashtable <Integer, Float> outlierHashtableSmall;
	Hashtable <Long, Float> outlierHashtableBig;

	public OutlierHashtable(long sampleFingerPrintOfTheProj, long markerFingerPrintOfTheProj) {
		this.sampFingerPrintFromProj = sampleFingerPrintOfTheProj;
		this.markFingerPrintFromProj = markerFingerPrintOfTheProj;
		this.outlierHashtableBig = new Hashtable<Long, Float>();
	}

    public void add(int markerIndexInProj, int sampleIndexInProj, byte dataItem, float value) {
    	if (outlierHashtableSmall == null) {
    		outlierHashtableBig.put((long) (markerIndexInProj * MAX_NUM_OF_SAMPLES_BIG_HASH + sampleIndexInProj * 10 + dataItem), value);
    	} else {
    		outlierHashtableSmall.put(markerIndexInProj * MAX_NUM_OF_SAMPLES_SMALL_HASH + sampleIndexInProj * 10 + dataItem, value);
    	}
    }

	public float getValue(int markerIndexInProj, int sampleIndexInProj, byte dataItem) {
    	if (outlierHashtableSmall == null) {
    		return outlierHashtableBig.get((long) markerIndexInProj * MAX_NUM_OF_SAMPLES_BIG_HASH + sampleIndexInProj * 10 + dataItem);
    	} else {
    		return outlierHashtableSmall.get(markerIndexInProj * MAX_NUM_OF_SAMPLES_SMALL_HASH + sampleIndexInProj * 10 + dataItem);
    	}
    }

	public long getSampFingerPrint() {
		return sampFingerPrintFromProj;
    }

	public long getMarkFingerPrint() {
		return markFingerPrintFromProj;
    }

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
    }

	public static OutlierHashtable load(String filename) {
		return (OutlierHashtable) Files.readSerial(filename);
    }
}
