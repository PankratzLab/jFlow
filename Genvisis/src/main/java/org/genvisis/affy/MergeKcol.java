package org.genvisis.affy;

import java.io.BufferedReader;
// import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class MergeKcol implements Runnable {
	public static final String[][] SNP_HEADER_OPTIONS = {{"SNP", "rsID", "ProbeSetName", "Name"}};
	public static final String FILENAME_AS_ID_OPTION = "[FILENAME_ROOT]";

	private final String[] files;
	// private long timeBegan;
	// private int threadId;
	private final String kColDir;

	public MergeKcol(String[] files, long timeBegan) {
		this(files, timeBegan, -1, "");
	}

	public MergeKcol(String[] files, long timeBegan, int threadId, String kColDir) {
		this.files = files;
		// this.timeBegan = timeBegan;
		// this.threadId = threadId;
		this.kColDir = kColDir;
	}

	@Override
	public void run() {

		BufferedReader reader;
		PrintWriter writer;
		String delimiter;
		delimiter = "\t";
		String[] line;
		// common folder output by apt-genotype, present in each subdirectory of Source
		// check source directory
		String[] dirList = Files.listDirectories(kColDir, false);

		int counts = 0;

		for (int j = 0; j < files.length; j++) {
			writer = Files.getAppropriateWriter(kColDir + files[j]);
			System.out.println("merging files " + (j + 1) + " of " + files.length);

			for (int i = 0; i < dirList.length; i++) {
				try {
					reader = Files.getAppropriateReader(kColDir + dirList[i] + "/" + files[j]);
					// filter comments
					do {
						line = reader.readLine().trim().split(delimiter, -1);
					} while (reader.ready() && (ext.indexFactors(	SNP_HEADER_OPTIONS, line, false, true, false,
																												false)[0] == -1));
					// if its the first directory, print the header

					if (i == 0) {
						writer.println(Array.toStr(line));
					}
					while (reader.ready()) {
						line = reader.readLine().trim().split(delimiter, -1);
						writer.println(Array.toStr(line));
						counts++;
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""	+ files[j] + "\" not found in " + kColDir
															+ dirList[i]);
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + kColDir + dirList[i] + files[j] + "\"");
					return;
				}
			}
			writer.close();
			System.out.println("Merged file contains " + counts + " IDs");
			counts = 0;
		}
	}

	public static void combineKcolFiles(int numThreads, String kColDir, String prefix,
																			String suffix) {
		Vector<Vector<String>> fileCabinet;
		Thread[] threads;
		boolean complete;
		long timeBegan;
		// common folder output by apt-genotype, present in each subdirectory of Source
		// check source directory
		timeBegan = new Date().getTime();

		String[] dirList = Files.listDirectories(kColDir, false);
		String[] files = Files.list(kColDir + dirList[0], prefix, suffix, false, false);
		System.out.println(files.length);
		System.out.println(dirList.length);
		fileCabinet = new Vector<Vector<String>>();
		for (int i = 0; i < numThreads; i++) {
			fileCabinet.add(new Vector<String>());
		}
		for (int i = 0; i < files.length; i++) {
			fileCabinet.elementAt(i % numThreads).add(files[i]);
		}
		threads = new Thread[numThreads];
		for (int i = 0; i < numThreads; i++) {
			threads[i] = new Thread(new MergeKcol(fileCabinet	.elementAt(i)
																												.toArray(new String[fileCabinet	.elementAt(i)
																																												.size()]),
																						timeBegan, i, kColDir));
			threads[i].start();
			try {
				Thread.sleep(100L);
			} catch (InterruptedException ex) {
			}
		}

		complete = false;
		while (!complete) {
			complete = true;
			for (int i = 0; i < numThreads; i++) {
				if (threads[i].isAlive()) {
					complete = false;
				}
			}
			if (!complete) {
				try {
					Thread.sleep(1000L);
				} catch (InterruptedException ex) {
				}
			}
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		int numThreads = 8;
		String kColDir = "";
		String suffix = null;
		String prefix = "gw6_split";
		String usage = "\n"	+ "affy.MergeChprequires 0-1 arguments\n"
										+ "   (2) number of threads to use (i.e. threads=" + numThreads
										+ " (default))\n";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				return;
			} else if (arg.startsWith("threads=")) {
				numThreads = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("inputDir=")) {
				kColDir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("prefix=")) {
				prefix = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		// filename = "C:/workspace/Genvisis/projects/dbGaP_test_CEL.properties";
		try {
			combineKcolFiles(numThreads, kColDir, prefix, suffix);
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Done Merging Files");
	}
}
