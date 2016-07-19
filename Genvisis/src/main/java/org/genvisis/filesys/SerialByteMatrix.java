package org.genvisis.filesys;

import java.io.Serializable;

import org.genvisis.common.*;

public class SerialByteMatrix implements Serializable {
	public static final long serialVersionUID = 1L;
	private byte[][] matrix;

	public byte[][] getMatrix() {
		return matrix;
	}
	
	public SerialByteMatrix(byte[][] matrix) {
		this.matrix = matrix;
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static SerialByteMatrix load(String filename, boolean jar) {
		return (SerialByteMatrix)Files.readSerial(filename, jar, true);
	}
}
