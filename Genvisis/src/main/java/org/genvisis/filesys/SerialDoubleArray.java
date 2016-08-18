package org.genvisis.filesys;

import java.io.Serializable;

import org.genvisis.common.*;

public class SerialDoubleArray implements Serializable {
	public static final long serialVersionUID = 1L;
	private double[] array;

	public SerialDoubleArray(double[] array) {
		this.array = array;
	}
	
	public double[] getArray() {
		return array;
	}
	
	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}

	public static SerialDoubleArray load(String filename, boolean jar) {
		return (SerialDoubleArray)SerializedFiles.readSerial(filename, jar, true);
	}
}
