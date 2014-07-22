package cnv.filesys;

import java.io.*;

import common.Files;
import common.Sort;

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
		Files.writeSerial(this, filename);
	}

	public static SampleList load(String filename, boolean jar) {
		return (SampleList)Files.readSerial(filename, jar, true);
	}

	public static SampleList generateSampleList(Project proj) {
		String[] files, samples;
		SampleList list;
		int[] keys;

//		if (Files.list(proj.getDir(Project.MARKER_DATA_DIRECTORY, true), MarkerData.MARKER_DATA_FILE_EXTENSION, proj.getJarStatus()).length>0) {
//			System.err.println("Error - Refusing to create new SampleList until the plots directory is either deleted or emptied; altering the SampleList will invalidate those files");
//			System.exit(1);
//		}

		if (Files.exists(proj.getDir(Project.SAMPLE_DIRECTORY), false)) {
			files = new File(proj.getDir(Project.SAMPLE_DIRECTORY)).list(new FilenameFilter() {
				public boolean accept(File file, String filename) {
					return filename.endsWith(Sample.SAMPLE_DATA_FILE_EXTENSION);
				}
			});
		} else {
			System.err.println("Error - failed to find the SAMPLE_DIRECTORY ("+proj.getDir(Project.SAMPLE_DIRECTORY)+"); no SampleList could be generated");
			return null;
		}
			
		keys = Sort.quicksort(files);
		samples = new String[files.length];
		for (int i = 0; i<samples.length; i++) {
			samples[i] = files[keys[i]].substring(0, files[keys[i]].lastIndexOf("."));
		}
		list = new SampleList(samples);
		list.serialize(proj.getFilename(Project.SAMPLELIST_FILENAME, true, true));

		return list;
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = cnv.Launch.getDefaultDebugProjectFile();

		String usage = "\n"+
		"filesys.SampleList requires 1 argument\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		try {
			generateSampleList(new Project(filename, false));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
