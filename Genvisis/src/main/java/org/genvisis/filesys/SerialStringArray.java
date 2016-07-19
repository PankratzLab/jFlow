package org.genvisis.filesys;

import java.io.Serializable;

import org.genvisis.common.*;

public class SerialStringArray implements Serializable {
	public static final long serialVersionUID = 1L;
	private String[] array;

	public SerialStringArray(String[] array) {
		this.array = array;
	}
	
	public String[] getArray() {
		return array;
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static SerialStringArray load(String filename, boolean jar) {
		return (SerialStringArray)Files.readSerial(filename, jar, true);
	}

	public static void dump(String filename) {
		Files.writeList(load(filename, false).getArray(), filename+".xln");
	}
}
