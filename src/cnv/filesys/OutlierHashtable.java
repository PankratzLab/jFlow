package cnv.filesys;

import java.io.*;
import java.util.Hashtable;

import common.Elision;

public class OutlierHashtable {

	public static final byte REDUCED_PRECISION_XY_NUM_BYTES = 2;

	Hashtable <Long, Float> outlierHashtable;

	public OutlierHashtable() {
		this.outlierHashtable = new Hashtable<Long, Float>();
	}

    public Hashtable <Long, Float> getOutlierHashtable() {
        return this.outlierHashtable;
    }

}
