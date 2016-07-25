package org.genvisis.common;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.ObjectStreamClass;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class SerializedFiles {

	public static void writeSerial(Object o, String filename) {
		writeSerial(o, filename, false);
	}

	public static void writeSerial(Object o, String filename, boolean gzip) {
		ObjectOutputStream oos;
	
		try {
			if (gzip) {
				oos = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(filename)));
			} else {
				oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(filename)));
			}
			oos.writeObject(o);
			oos.flush();
			oos.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static Object readSerial(String filename) {
		return readSerial(filename, false, new Logger(), true);
	}

	public static Object readSerial(String filename, boolean jar, boolean kill) {
		return readSerial(filename, jar, new Logger(), kill);
	}

	/**
	 * Checks if a serialized version of the non serial filename exists
	 * 
	 * @param altDir
	 *            if the file should be in a different directory (optional)
	 * @param nonSerialFilename
	 * @return
	 */
	public static boolean serializedVersionExists(String altDir, String nonSerialFilename) {
		return Files.exists(getSerializedFileName(altDir, nonSerialFilename));
	}

	/**
	 * returns a serialized version of the non serial filename
	 * 
	 * @param altDir
	 *            if the file should be in a different directory (optional)
	 * @param nonSerialFilename
	 * @return
	 */
	public static String getSerializedFileName(String altDir, String nonSerialFilename) {
		if (altDir == null) {
			return ext.rootOf(nonSerialFilename, false) + Files.SERIALIZED_FILE_EXTENSION;
		} else {
			return altDir + ext.rootOf(nonSerialFilename, true) + Files.SERIALIZED_FILE_EXTENSION;
		}
	}

	public static Object readSerial(String filename, boolean jar, Logger log, boolean kill) {
		return readSerial(filename, jar, log, kill, filename.endsWith(".gz"));
	}

	public static Object readSerial(String filename, boolean jar, Logger log, boolean kill, boolean gzipped){
			InputStream in;
			ObjectInputStream ois = null;
			Object o = null;
	
			try {
				if (jar) {
					in = new BufferedInputStream(ClassLoader.getSystemResourceAsStream(filename));
				} else if (gzipped) {
					in = new GZIPInputStream(new FileInputStream(filename));
				} else {
					in = new BufferedInputStream(new FileInputStream(filename));
				}
				ois = new ObjectInputStream(in);
				o = ois.readObject();
				ois.close();
			} catch (Exception e) {
				try {
					if (ois != null) {
						ois.close();
					}
					o = readSerialFixClassname(filename, jar, log, kill, gzipped);
				} catch (Exception e2) {
					log.reportError("Error - failed to load " + filename);
					log.reportException(e2);
					if (kill) {
						System.exit(1);
					}
				}

			} 
	
			return o;
		}
	
	private static synchronized Object readSerialFixClassname(String filename, boolean jar, Logger log, boolean kill, boolean gzipped) throws FileNotFoundException, IOException, ClassNotFoundException {
		InputStream in;
		ObjectInputStream ois;
		Object o;
		
		if (jar) {
			in = new BufferedInputStream(ClassLoader.getSystemResourceAsStream(filename));
		} else if (gzipped) {
			in = new GZIPInputStream(new FileInputStream(filename));
		} else {
			in = new BufferedInputStream(new FileInputStream(filename));
		}
		ois = new ObjectInputStream(in) {
			@Override
			protected ObjectStreamClass readClassDescriptor() throws IOException, ClassNotFoundException {
				ObjectStreamClass resultClassDescriptor = null;
				try {
					resultClassDescriptor = super.readClassDescriptor();
					Class.forName(resultClassDescriptor.getName());
				} catch (ClassNotFoundException cnfe) {
					String fileClassDesc = resultClassDescriptor.getName();
					int startInsert = 0;
					for (int i = 0; i < fileClassDesc.length(); i++) {
						if (Character.isLowerCase(fileClassDesc.charAt(i))){
							startInsert = i;
							break;
						}
					}
					String convertedClassDesc = fileClassDesc.substring(0, startInsert) + "org.genvisis." + fileClassDesc.substring(startInsert);
					Class<?> convertedClass = Class.forName(convertedClassDesc);
					log.reportTimeWarning("The Class (" + fileClassDesc +  ") for the Serialized Object " + filename + " cannot be resolved, attempting to use " + convertedClassDesc);
					resultClassDescriptor = ObjectStreamClass.lookup(convertedClass);
				}
				return resultClassDescriptor;
			}
		};
		o = ois.readObject();
		ois.close();
		log.report("Succesfully deserialized " + filename + " to " + o.getClass().getName());
		writeSerial(o, filename, gzipped);
		log.report("Succesfully rewrote " + filename + " as a serialized " + o.getClass().getName());
		return o;
	}
}
