package filesys;

import java.io.Serializable;
import common.*;

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
		Files.writeSerial(this, filename);
	}

	public static SerialDoubleArray load(String filename, boolean jar) {
		return (SerialDoubleArray)Files.readSerial(filename, jar, true);
	}
}
