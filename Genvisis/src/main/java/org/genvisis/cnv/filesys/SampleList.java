package org.genvisis.cnv.filesys;

import java.io.*;
import java.util.Hashtable;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class SampleList implements Serializable {
	public static final long serialVersionUID = 1L;

	private long fingerprint;
	private String[] samples;

	public SampleList(String[] samples) {
		this.samples = samples;
		this.fingerprint = MarkerSet.fingerprint(samples);
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public String[] getSamples() {
		return samples;
	}

	public void writeToTextFile(String filename) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (int i = 0; i<samples.length; i++) {
				writer.println(samples[i]);
			}
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error writing to "+filename);
			ioe.printStackTrace();
		}
	}

	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}

	public static SampleList load(String filename, boolean jar) {
		return (SampleList)SerializedFiles.readSerial(filename, jar, true);
	}

	public static SampleList generateSampleList(Project proj) {
		String[] files, samples;
		SampleList list;
		int[] keys;
		int countAt;
		Logger log;
		
		log = proj.getLog();

//		if (Files.list(proj.getDir(Project.MARKER_DATA_DIRECTORY, true), MarkerData.MARKER_DATA_FILE_EXTENSION, proj.getJarStatus()).length>0) {
//			System.err.println("Error - Refusing to create new SampleList until the plots directory is either deleted or emptied; altering the SampleList will invalidate those files");
//			System.exit(1);
//		}

		if (Files.exists(proj.SAMPLE_DIRECTORY.getValue(false, true), false)) {
			files = Files.list(proj.SAMPLE_DIRECTORY.getValue(false, true), Sample.SAMPLE_FILE_EXTENSION, false);
		} else {
			log.reportError("Error - failed to find the SAMPLE_DIRECTORY ("+proj.SAMPLE_DIRECTORY.getValue(false, true)+"); no SampleList could be generated");
			return null;
		}

		countAt = 0;
		keys = Sort.quicksort(files);
		samples = new String[files.length];
		for (int i = 0; i<samples.length; i++) {
			samples[i] = files[keys[i]].substring(0, files[keys[i]].lastIndexOf("."));
			if (samples[i].contains("@")) {
				countAt++;
			}
		}
		list = new SampleList(samples);
		if (samples.length > 0) {
			list.serialize(proj.SAMPLELIST_FILENAME.getValue(true, true));
		} else {
			log.reportError("Error - there are no samples in the samples directory; parsing must have failed, so cannot create a SampleList");
		}
		if (countAt > 0) {
			proj.getLog().report("Note - "+countAt+" ("+(Double.parseDouble(ext.prettyP((double)countAt/(double)samples.length))*100)+"%) of your Sample IDs contain the @ symbol, which is often used when concatenating the sample's bar code. If you would like these to be stripped, then set "+proj.PARSE_AT_AT_SYMBOL.getName()+"=TRUE in the properties file, delete the samples/ directory and reparse the data");
		}

		return list;
	}

	public static void serializeOutliers(Project proj) {
		String[] samples = proj.getSamples();
		String sampleDir = proj.SAMPLE_DIRECTORY.getValue(false, true);
		
		Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();
		
		proj.getLog().report("Extracting outliers from " + samples.length + " Sample files");
		for (int i = 0; i < samples.length; i++) {
			try {
				String sampleFile = sampleDir + samples[i] + Sample.SAMPLE_FILE_EXTENSION;
				Hashtable<String, Float> outliers = Sample.loadOutOfRangeValuesFromRandomAccessFile(sampleFile);
				
				// TODO duplicates?
				if (outliers != null) {
					for (java.util.Map.Entry<String, Float> entry : outliers.entrySet()) {
						String[] keyParts = entry.getKey().split("\t");
						allOutliers.put(keyParts[0] + "\t" + samples[i] + "\t" + keyParts[1], entry.getValue());
					}
				}
				// for duplicates, loop through outliers and use allOutliers.putIfAbsent and compare returned value
			} catch (Exception e) {
				proj.getLog().reportException(e);
			}
			
		}
		SerializedFiles.writeSerial(allOutliers, sampleDir + "outliers.ser");
		
	}
	
	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = null;
		boolean outliers = false;
		
		String usage = "\n"+
		"filesys.SampleList requires 1+ argument\n"+
		"   (1) project properties filename (i.e. proj="+org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"   (2) (Optional) use -outliers argument to generate an 'outliers.ser' file (not the default)\n" +
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-outliers")) {
				outliers = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		if (outliers) {
			serializeOutliers(new Project(filename, false));
		} else {
			try {
				generateSampleList(new Project(filename, false));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

}
