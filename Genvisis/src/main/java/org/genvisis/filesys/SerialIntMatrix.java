package org.genvisis.filesys;

import java.io.Serializable;

import org.genvisis.common.*;

public class SerialIntMatrix implements Serializable {
	public static final long serialVersionUID = 1L;
	private int[][] matrix;

	public int[][] getMatrix() {
		return matrix;
	}
	
	public SerialIntMatrix(int[][] matrix) {
		this.matrix = matrix;
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static SerialIntMatrix load(String filename, boolean jar) {
		return (SerialIntMatrix)Files.readSerial(filename, jar, true);
	}
}
