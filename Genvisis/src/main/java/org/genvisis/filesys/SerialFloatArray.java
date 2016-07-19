package org.genvisis.filesys;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;

import org.genvisis.common.*;

public class SerialFloatArray implements Serializable {
	public static final long serialVersionUID = 1L;
	private float[] array;

	public SerialFloatArray(float[] array) {
		this.array = array;
	}
	
	public float[] getArray() {
		return array;
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static SerialFloatArray load(String filename, boolean jar) {
		return (SerialFloatArray)Files.readSerial(filename, jar, true);
	}
	
	public static void dump(String filename) {
		float[] all = load(filename, false).getArray();
		
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(filename+".xln"));
			for (int i = 0; i < all.length; i++) {
				writer.println(all[i]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename+".xln");
			e.printStackTrace();
		}
	}
}
