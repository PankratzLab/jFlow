// -Xms1024M -Xmx1024M
package org.genvisis.common;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import org.genvisis.filesys.SerialHash;
import org.genvisis.parse.GenParser;

import com.google.common.collect.Lists;
import com.google.common.io.Closeables;
import com.google.common.primitives.Ints;

// class DeleteLater implements Runnable {
// private String filename;
// private int timeToWait;
//
// public DeleteLater(String filename, int timeToWait) {
// this.filename = filename;
// this.timeToWait = timeToWait;
// }
// public void run() {
// try {
// Thread.sleep(timeToWait);
// } catch (InterruptedException e) {}
// new File(filename).delete();
// }
// }

public class Files {
	public static final String ROOT_DIRECTORY = "/home/npankrat/"; // galileo
	// public static final String ROOT_DIRECTORY = "/export/home/npankrat/"; // alcatraz
	// public static final String ROOT_DIRECTORY = "/state/partition1/npankrat/"; // indiviudal nodes
	public static final String JAVA = "/usr/java/latest/bin/java";
	public static final String JCP = JAVA	+ " -jar /home/npankrat/"
																		+ org.genvisis.common.PSF.Java.GENVISIS;
	public static final String SERIALIZED_FILE_EXTENSION = ".ser";

	public static String getJarLocation() {
		try {
			return new File(Files.class	.getProtectionDomain().getCodeSource().getLocation().toURI()
																	.getPath()).getAbsolutePath();
		} catch (URISyntaxException e) {
			return "~/" + org.genvisis.common.PSF.Java.GENVISIS;
		}
	}

	public static String getRunString() {
		return getRunString(-1);
	}

	public static String getRunString(int memInMeg) {
		return getRunString(memInMeg, false);
	}

	public static String getRunString(int memInMeg, boolean interpretAsGig) {
		int mem = memInMeg == -1 ? -1 : interpretAsGig ? memInMeg / 1024 : memInMeg;
		return "java "	+ (memInMeg == -1 ? "" : "-Xmx" + mem + (interpretAsGig ? "G" : "M")) + " -jar "
						+ getJarLocation();
	}

	public static void batchIt(	String root_batch_name, String init, int numBatches, String commands,
															String[] iterations) {
		String[][] iters = new String[iterations.length][1];

		for (int i = 0; i < iters.length; i++) {
			iters[i][0] = iterations[i];
		}

		batchIt(root_batch_name, init, numBatches, commands, iters);
	}

	public static void batchIt(	String root_batch_name, int sleep, int start, int stop, int numBatches,
															String commands) {
		PrintWriter[] writers = new PrintWriter[numBatches];
		PrintWriter writer;
		String batchName;

		try {
			for (int i = 0; i < numBatches; i++) {
				writers[i] = new PrintWriter(new FileWriter(getBatchName(root_batch_name, i, numBatches)));
				if (!Files.isWindows()) {
					writers[i].println("#/bin/sh\n");
				}
				if (sleep > 0) {
					writers[i].println("sleep " + sleep + "\n");
				}
			}
			for (int i = start; i <= stop; i++) {
				writer = writers[i % numBatches];
				writer.println(ext.insertNumbers(commands, i));
			}
			if (numBatches > 1) {
				writer = new PrintWriter(new FileWriter("master"));
			} else {
				writer = null;
			}
			for (int i = 0; i < numBatches; i++) {
				writers[i].close();
				batchName = getBatchName(root_batch_name, i, numBatches);
				chmod(batchName, i == 0);
				if (numBatches > 1) {
					writer.println("nohup ./" + batchName + " > " + (i + 1) + ".out &");
				}

			}
			if (numBatches > 1) {
				writer.close();
			}

		} catch (IOException ioe) {
			throw new RuntimeException("Problem creating batch files named " + root_batch_name);
		}
	}

	private static String getBatchName(String root, int i, int numBatches) {
		String batchName;

		batchName = numBatches == 1 ? root : root + "." + (i + 1);
		if (Files.isWindows()) {
			batchName += ".bat";
		}

		return batchName;
	}

	public static void qsub(String root, int start, int stop, String commands) {
		qsub(root, start, stop, commands, -1, -1);
	}

	public static void qsub(String root, int start, int stop, String commands, int memRequiredInMb,
													double walltimeRequestedInHours) {
		qsub("", root, start, stop, commands, memRequiredInMb, walltimeRequestedInHours, null);
	}

	public static void qsub(String dir, String root, int start, int stop, String commands,
													int memRequiredInMb, double walltimeRequestedInHours, String queue) {
		String[] lines;

		lines = qsub(	dir, "chr#_" + root, start, stop, commands, null, memRequiredInMb,
									walltimeRequestedInHours, null);

		if (lines.length > 1) {
			writeArray(lines, dir + "master." + (root == null ? "qsub" : root));
			chmod(dir + "master." + (root == null ? "qsub" : root));
		}
	}

	private static void writeQsubHeader(PrintWriter writer, String filename,
																			int totalMemoryRequestedInMb, double walltimeRequestedInHours,
																			int numProcs, String nodeToUse) {
		writeQsubHeader(writer, filename, totalMemoryRequestedInMb, walltimeRequestedInHours, numProcs,
										nodeToUse, true);
	}


	private static void writeQsubHeader(PrintWriter writer, String filename,
																			int totalMemoryRequestedInMb, double walltimeRequestedInHours,
																			int numProcs, String nodeToUse, boolean defualtMods) {
		Vector<String> params;
		int hours, minutes;

		writer.println("#!/bin/bash");
		writer.println("#$ -cwd");
		writer.println("#$ -S /bin/bash");
		writer.println("#PBS -e $PBS_JOBNAME.$PBS_JOBID.e");
		writer.println("#PBS -o $PBS_JOBNAME.$PBS_JOBID.o");
		params = new Vector<String>();
		if (totalMemoryRequestedInMb > 0) {
			params.add("mem=" + totalMemoryRequestedInMb + "mb");
		}
		if (walltimeRequestedInHours > 0) {
			hours = (int) Math.floor(walltimeRequestedInHours);
			minutes = (int) Math.ceil((walltimeRequestedInHours - hours) * 60.0);
			params.add("walltime=" + ext.formNum(hours, 2) + ":" + ext.formNum(minutes, 2) + ":00");
		}
		params.add("nodes=1:ppn=" + (numProcs <= 0 ? "1" : numProcs));
		if (params.size() > 0) {
			writer.println("#PBS -m ae"); // send mail when aborts or ends (add b, as in #PBS -m abe, for
																		// begins as well)
			writer.println("#PBS -l " + Array.toStr(Array.toStringArray(params), ","));
		}


		if (nodeToUse != null) {
			writer.println("#$ -q *@" + nodeToUse);
		}
		writer.println();

		if (defualtMods) {
			writer.println(Array.toStr(PSF.Load.getAllModules(), "\n"));
			writer.println(PSF.Cmd.getProfilePLWithOutput(ext.rootOf(filename, false) + ".profile"));
		}

		writer.println();
		writer.println("echo \"start " + ext.rootOf(filename) + " at: \" `date`");
		writer.println("/bin/hostname");

	}

	public static String[] qsub(String dir, String filenameFormat, int start, int stop,
															String commands, String[] patterns, int totalMemoryRequestedInMb,
															double walltimeRequestedInHours, String nodeToUse) {
		return qsub(dir, filenameFormat, start, stop, commands, patterns, totalMemoryRequestedInMb,
								walltimeRequestedInHours, nodeToUse, null);
	}

	@Deprecated
	public static String[] qsub(String dir, String filenameFormat, int start, int stop,
															String commands, String[] patterns, int totalMemoryRequestedInMb,
															double walltimeRequestedInHours, String nodeToUse, String queueName) {
		PrintWriter writer;
		String filename;
		String[] lines;
		Vector<String> v;

		if (dir == null) {
			dir = "";
		} else if (!dir.equals("") && !new File(dir).exists()) {
			System.err.println("Error - directory '"	+ dir
													+ "' does not exist, cannot create batches there");
		}

		v = new Vector<String>();
		for (int i = start; i <= stop; i++) {
			// filename = ext.parseDirectoryOfFile(prefix,
			// true)+(prefix==null?"":ext.removeDirectoryInfo(prefix))+(start!=stop||(start>0&&start<25)?(prefix==null||prefix.equals("")?i:(prefix.endsWith("chr")?"":".")+i):"")+(root==null?"":"_"+root)+".qsub";
			filename = ext.insertNumbers(filenameFormat, i)
									+ (filenameFormat.endsWith(".qsub") ? "" : ".qsub");
			try {
				writer = new PrintWriter(new FileWriter(dir + filename));
				writeQsubHeader(writer, filename, totalMemoryRequestedInMb, walltimeRequestedInHours, 1,
												nodeToUse);
				if (patterns == null) {
					writer.println(ext.insertNumbers(commands, i));
				} else {
					writer.println("/bin/hostname");
					writer.println(">plug");
					writer.println("rep=0");
					writer.println("total_reps=0");
					writer.println("while [ -e \"plug\" ]; do ");
					writer.println("    rep=$(java -jar "	+ ROOT_DIRECTORY
													+ org.genvisis.common.PSF.Java.GENVISIS
													+ " common.Files -nextRep patterns=" + Array.toStr(patterns, ",")
													+ " lastRep=$rep wait=1000)");
					writer.println("    echo \"Beginning replicate $rep\"");
					lines = commands.split("\n");
					for (String line : lines) {
						writer.println("    " + ext.replaceAllWith(line, "#", "$rep"));
					}
					writer.println("    rm $rep.taken");
					writer.println("    total_reps=$(expr $total_reps + 1)");
					writer.println("done");
					writer.println("echo \"Performed $total_reps replicate(s) within this thread\"");
				}
				writer.println("echo \"end " + ext.rootOf(filename) + " at: \" `date`");
				writer.close();
				chmod(dir + filename, false);
				v.add("qsub " + (queueName == null ? "" : "-q " + queueName + " ") + filename);
			} catch (IOException ioe) {
				throw new RuntimeException("Problem creating " + dir + filename);
			}
		}
		return Array.toStringArray(v);
	}

	// TODO get rid of memoryPerProcRequestedInMb at some point
	public static void qsubMultiple(String chunkFilename, String[] jobs, int numJobsToForce,
																	int memoryPerProcRequestedInMb, int totalMemoryRequestedInMb,
																	double walltimeRequestedInHours) {
		PrintWriter writer;
		int numProcs;

		if (jobs.length > 16) {
			System.err.println("Error - not optimized past 16 processes; due some testing and implement");
			return;
		}

		try {
			writer = new PrintWriter(new FileWriter(chunkFilename));

			if (numJobsToForce <= 0) {
				numProcs = jobs.length;
			} else if (numJobsToForce < jobs.length) {
				System.err.println("Error - cannot force fewer jobs than are provided ("	+ numJobsToForce
														+ "<" + jobs.length + ")");
				writer.close();
				return;
			} else {
				numProcs = numJobsToForce;
			}
			writeQsubHeader(writer, chunkFilename, totalMemoryRequestedInMb, walltimeRequestedInHours,
											numProcs, null);


			for (int j = 0; j < jobs.length; j++) {
				writer.println("pbsdsh -n " + j + " " + jobs[j] + " &");

			}
			writer.println("wait");
			writer.println("echo \"end at: \" `date`");
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + chunkFilename);
			e.printStackTrace();
		}
	}

	public static void qsubMultiple(Vector<String> jobNamesWithAbsolutePaths, IntVector jobSizes,
																	String batchDir, String batchRoot, int maxJobsPerBatch,
																	boolean forceMaxJobs, String queueName,
																	int memoryPerProcRequestedInMb, int totalMemoryRequestedInMb,
																	double walltimeRequestedInHours) {
		String[] files;
		int count;
		Vector<String> v;

		if (jobSizes == null) {
			files = Array.toStringArray(jobNamesWithAbsolutePaths);
		} else {
			int[] jobOrder = Sort.getReverseIndices(Ints.toArray(jobSizes));
			files = Sort.getOrdered(jobNamesWithAbsolutePaths, jobOrder);
		}
		count = 0;
		v = new Vector<String>();
		new File(batchDir).mkdirs();
		v.add("cd " + batchDir);
		while (count * maxJobsPerBatch < files.length) {
			Files.qsubMultiple(batchDir	+ ext.removeDirectoryInfo(batchRoot) + "." + count,
													Array.subArray(	files, count * maxJobsPerBatch,
																					Math.min((count + 1) * maxJobsPerBatch, files.length)),
													forceMaxJobs ? maxJobsPerBatch : -1, memoryPerProcRequestedInMb,
													totalMemoryRequestedInMb, walltimeRequestedInHours);
			v.add("qsub "	+ (queueName == null ? "" : "-q " + queueName + " ")
						+ ext.removeDirectoryInfo(batchRoot) + "." + count);
			count++;
		}
		Files.writeArray(Array.toStringArray(v), ext.parseDirectoryOfFile(batchRoot)	+ "master."
																							+ ext.removeDirectoryInfo(batchRoot));
		Files.chmod(ext.parseDirectoryOfFile(batchRoot)	+ "master."
								+ ext.removeDirectoryInfo(batchRoot));
	}

	public static void qsubExecutor(String dirToSwitchToBeforeRunning,
																	Vector<String> commandsWithAbsolutePaths, IntVector jobSizes,
																	String batchRoot, int numProcs, int totalMemoryRequestedInMb,
																	double walltimeRequestedInHours) {
		String[] commands;

		if (jobSizes == null) {
			commands = Array.toStringArray(commandsWithAbsolutePaths);
		} else {
			int[] jobOrder = Sort.getReverseIndices(Ints.toArray(jobSizes));
			commands = Sort.getOrdered(commandsWithAbsolutePaths, jobOrder);
		}

		Files.writeArray(commands, batchRoot + ".chain");
		Files.qsub(batchRoot	+ ".pbs",
								"cd "			+ dirToSwitchToBeforeRunning + "\njava -jar ~/"
													+ org.genvisis.common.PSF.Java.GENVISIS + " one.ScriptExecutor file="
													+ batchRoot + ".chain threads=" + numProcs,
								totalMemoryRequestedInMb, walltimeRequestedInHours, numProcs);
	}

	public static void batchIt(	String root_batch_name, String init, int numBatches, String commands,
															String[][] iterations) {
		PrintWriter[] writers = new PrintWriter[numBatches];
		PrintWriter writer;
		String trav;
		String[] lines;

		try {
			for (int i = 0; i < numBatches; i++) {
				writers[i] = new PrintWriter(new FileWriter(i == 0 && numBatches == 1	? root_batch_name
																																							: root_batch_name
																																								+ "." + (i + 1)));
				if (init != null) {
					writers[i].println(init);
					writers[i].println();
				}
			}
			for (int i = 0; i < iterations.length; i++) {
				trav = commands;
				for (int j = 0; j < iterations[i].length; j++) {
					trav = ext.replaceAllWith(trav, "[%" + j + "]", iterations[i][j]);
				}
				writer = writers[i % numBatches];
				lines = trav.split("\\n");
				for (String line : lines) {
					writer.println(line);
				}
			}
			if (numBatches > 1) {
				writer = new PrintWriter(new FileWriter("master"));
			} else {
				writer = null;
			}
			for (int i = 0; i < numBatches; i++) {
				writers[i].close();
				chmod(i == 0	&& numBatches == 1 ? root_batch_name : root_batch_name + "." + (i + 1),
							i == 0);
				if (numBatches > 1) {
					writer.println("nohup ./" + root_batch_name + "." + (i + 1) + " > " + (i + 1) + ".out &");
				}

			}
			if (numBatches > 1) {
				writer.close();
			}
		} catch (IOException ioe) {
			throw new RuntimeException("Problem creating batch files named " + root_batch_name);
		}
	}

	public static void qsub(String root_batch_name, String dirToSwitchToBeforeRunning, int numBatches,
													String commands, String[][] iterations, int totalMemoryRequestedInMb,
													double walltimeRequestedInHours) {
		qsub(	root_batch_name, dirToSwitchToBeforeRunning, numBatches, commands, iterations,
					totalMemoryRequestedInMb, walltimeRequestedInHours, null);
	}

	public static void qsub(String root_batch_name, String dirToSwitchToBeforeRunning, int numBatches,
													String commands, String[][] iterations, int totalMemoryRequestedInMb,
													double walltimeRequestedInHours, String queueName) {
		PrintWriter[] writers;
		PrintWriter writer;
		String trav;
		String[] lines;
		int index;
		String dir;

		if (iterations.length == 0) {
			System.out.println("No iterations specified for root " + root_batch_name);
			return;
		}

		if (numBatches <= 0) {
			numBatches = iterations.length;
		}

		dir = ext.parseDirectoryOfFile(root_batch_name, true);
		new File(dir).mkdirs();
		writers = new PrintWriter[numBatches];
		try {
			if (numBatches > 1) {
				writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(root_batch_name)
																									+ "master."
																								+ ext.removeDirectoryInfo(root_batch_name)));
				if (!dir.equals("./") && !dir.equals("")) {
					writer.println("cd " + dir);
				}
				for (int i = 0; i < numBatches; i++) {
					writers[i] = new PrintWriter(new FileWriter(root_batch_name + "_" + (i + 1) + ".qsub"));
					writeQsubHeader(writers[i], root_batch_name + "_" + (i + 1) + ".qsub",
													totalMemoryRequestedInMb, walltimeRequestedInHours, 1, null);
					if (dirToSwitchToBeforeRunning != null) {
						writers[i].println("cd " + dirToSwitchToBeforeRunning);
					}
					writer.println("qsub "	+ (queueName == null ? "" : "-q " + queueName + " ")
													+ ext.removeDirectoryInfo(root_batch_name) + "_" + (i + 1) + ".qsub");
				}
				writer.close();
				chmod(ext.parseDirectoryOfFile(root_batch_name)	+ "master."
							+ ext.removeDirectoryInfo(root_batch_name));
			} else {
				// writer = null;
				writers[0] = new PrintWriter(new FileWriter(root_batch_name + ".qsub"));
				writeQsubHeader(writers[0], root_batch_name + ".qsub", totalMemoryRequestedInMb,
												walltimeRequestedInHours, 1, null);
			}


			for (int i = 0; i < iterations.length; i++) {
				trav = commands;
				for (int j = 0; j < iterations[i].length; j++) {
					trav = ext.replaceAllWith(trav, "[%" + j + "]", iterations[i][j]);
				}
				lines = trav.split("\\n");
				index = i % numBatches;
				for (String line : lines) {
					writers[index].println(line);
				}
			}

			for (int i = 0; i < numBatches; i++) {
				writers[i].println("echo \"end "	+ ext.removeDirectoryInfo(root_batch_name) + "_" + (i + 1)
														+ " at: \" `date`");
				writers[i].close();
				chmod(root_batch_name + "_" + (i + 1) + ".qsub");
			}
			if (numBatches == 1) {
				chmod(root_batch_name + ".qsub");
			}

			// if (numBatches > 1) {
			// writer.close();
			// chmod(ext.parseDirectoryOfFile(root_batch_name)+"master." +
			// ext.removeDirectoryInfo(root_batch_name));
			// }
		} catch (IOException ioe) {
			ioe.printStackTrace();
			throw new RuntimeException("Problem creating batch files named " + root_batch_name);
		}
	}



	public static void qsub(String root, String commands, String[][] iterations) {
		qsub(root, commands, iterations, -1, 10, -1);
	}

	public static void qsub(String root, String commands, String[][] iterations,
													int totalMemoryRequestedInMb, double walltimeRequestedInHours,
													int numProcs) {
		qsub(	root, commands, iterations, totalMemoryRequestedInMb, walltimeRequestedInHours, numProcs,
					null);
	}

	public static void qsub(String root, String commands, String[][] iterations,
													int totalMemoryRequestedInMb, double walltimeRequestedInHours,
													int numProcs, String queueName) {
		PrintWriter writer;
		String filename, trav;

		try {
			if (root == null) {
				filename = "master.qsub";
			} else {
				trav = root;
				for (int j = 0; j < iterations[0].length; j++) {
					trav = ext.replaceAllWith(trav, "[%" + j + "]", "");
				}
				if (ext.rootOf(trav).equals("")) {
					filename = ext.parseDirectoryOfFile(trav) + "master.qsub";
				} else {
					filename = ext.parseDirectoryOfFile(trav) + "master." + ext.rootOf(trav) + ".qsub";
				}

			}
			writer = new PrintWriter(new FileWriter(filename));
			for (int i = 0; i < iterations.length; i++) {
				if (root == null) {
					filename = i + ".qsub";
				} else if (ext.rootOf(root).equals("")) {
					filename = root + i + ".qsub";
				} else {
					trav = root;
					for (int j = 0; j < iterations[i].length; j++) {
						trav = ext.replaceAllWith(trav, "[%" + j + "]", iterations[i][j]);
					}
					if (trav.equals(root)) {
						filename = root + "_" + i + ".qsub";
					} else {
						filename = trav + ".qsub";
					}
				}

				trav = commands;
				for (int j = 0; j < iterations[i].length; j++) {
					trav = ext.replaceAllWith(trav, "[%" + j + "]", iterations[i][j]);
				}
				qsub(filename, trav, totalMemoryRequestedInMb, walltimeRequestedInHours, numProcs);
				writer.println("qsub " + (queueName == null ? "" : "-q " + queueName + " ") + filename);
			}
			writer.close();
			Files.chmod("master." + (root == null ? "qsub" : root));
		} catch (IOException ioe) {
			throw new RuntimeException("Problem creating master batch file");
		}
	}

	public static void qsub(String filename, String command, int totalMemoryRequestedInMb,
													double walltimeRequestedInHours, int numProcs) {
		PrintWriter writer;
		String[] lines;

		lines = command.split("\\n");
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writeQsubHeader(writer, filename, totalMemoryRequestedInMb, walltimeRequestedInHours,
											numProcs, null);

			for (String line : lines) {
				writer.println(line);
			}
			writer.println("echo \"end " + ext.rootOf(filename) + " at: \" `date`");
			writer.close();
			chmod(filename, false);
		} catch (IOException ioe) {
			throw new RuntimeException("Problem creating " + filename);
		}
	}

	public static void execListAdd(	Vector<String> execList, String commands, String[] iterations,
																	Logger log) {
		execListAdd(execList, commands, Matrix.toMatrix(iterations), log);
	}

	public static void execListAdd(	Vector<String> execList, String commands, String[][] iterations,
																	Logger log) {
		String trav;

		if (iterations == null || iterations.length == 0) {
			log.reportError("No iterations specified to add to the execList");
			return;
		}

		if (execList == null) {
			log.reportError("execList is null");
			return;
		}

		if (commands == null) {
			log.reportError("commands is null, cannot be added to the execList");
			return;
		}

		if (commands.split("\\n").length > 1) {
			log.reportError("Addition to execList has multiple lines, not appropriate for a serial executor");
			return;
		}

		for (String[] iteration : iterations) {
			trav = commands;
			for (int j = 0; j < iteration.length; j++) {
				trav = ext.replaceAllWith(trav, "[%" + j + "]", iteration[j]);
			}
			execList.addElement(trav);
		}
	}

	public static boolean chmod(String filename, boolean verbose) {
		if (Files.isWindows()) {
			if (verbose) {
				System.err.println("chmod not attempted on windows platform");
			}
			return false;
		} else {
			try {
				Runtime.getRuntime().exec("chmod +x " + filename).waitFor();
				return true;
			} catch (Exception e) {
				if (verbose) {
					System.err.println("Warning - chmod failed; "	+ filename
															+ " was not given execute permissions");
				}
				return false;
			}
		}
	}

	public static boolean chmod(String filename) {
		return chmod(filename, true);
	}

	public static boolean copyFileUsingFileChannels(String source, String dest, Logger log) {
		return copyFileUsingFileChannels(new File(source), new File(dest), log);
	}

	/**
	 * Method should handle gzip and other binary /compressed formats like .ser
	 *
	 * @param source source file
	 * @param dest destination file
	 * @param log
	 * @return if the copy was a success
	 */
	public static boolean copyFileUsingFileChannels(File source, File dest, Logger log) {
		FileChannel inputChannel = null;
		FileChannel outputChannel = null;
		boolean copy = false;
		try {
			inputChannel = new FileInputStream(source).getChannel();
			outputChannel = new FileOutputStream(dest).getChannel();
			outputChannel.transferFrom(inputChannel, 0, inputChannel.size());
			inputChannel.close();
			outputChannel.close();
			copy = true;
		} catch (FileNotFoundException e) {
			log.reportFileNotFound(source.getPath());
			log.reportException(e);
			e.printStackTrace();
		} catch (IOException e) {
			log.reportException(e);
			e.printStackTrace();
		}
		return copy;
	}


	// causes trouble with Serialized data
	public static boolean copyFile(String from, String to) {
		FileReader in;
		FileWriter out;
		int c;

		try {
			in = new FileReader(from);
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error - Cannot find " + from + " in current directory");
			return false;
		}

		try {
			out = new FileWriter(to);

			while ((c = in.read()) != -1) {
				out.write(c);
			}

			in.close();
			out.close();

			new File(to).setLastModified(new File(from).lastModified());

			return true;
		} catch (Exception e) {
			return false;
		}
	}

	public static boolean copyFileFromJar(String from, String to) {
		InputStream is;
		OutputStream os;

		try {
			is = ClassLoader.getSystemResourceAsStream(from);
			os = new FileOutputStream(to);

			byte[] buffer = new byte[4096];
			int length;
			while ((length = is.read(buffer)) > 0) {
				os.write(buffer, 0, length);
			}
			os.close();
			is.close();

			return true;
		} catch (Exception e) {
			return false;
		}
	}

	// super slow, but copies exactly
	public static boolean copyFileExactlyByteByByte(String from, String to) {
		FileInputStream in;
		FileOutputStream out;
		int c;

		in = null;
		out = null;
		try {
			in = new FileInputStream(from);
			out = new FileOutputStream(to);

			while ((c = in.read()) != -1) {
				out.write(c);
			}

			in.close();
			out.close();

			new File(to).setLastModified(new File(from).lastModified());

			return true;

		} catch (Exception e) {
			return false;
		}
	}

	public static boolean identicalFiles(String file1, String file2, boolean checkTime) {
		FileReader in1;
		FileReader in2;

		if (!new File(file1).exists() || !new File(file2).exists()) {
			System.err.println("Error - Cannot find one or both of the files ("	+ file1 + " or " + file2
													+ ")");
			return false;
		}

		if (checkTime && new File(file1).lastModified() != new File(file2).lastModified()) {
			return false;
		}

		try {
			in1 = new FileReader(file1);
			in2 = new FileReader(file2);

			while (in1.ready() && in2.ready()) {
				if (in1.read() != in2.read()) {
					in1.close();
					in2.close();
					return false;
				}
			}

			in1.close();
			in2.close();

			return true;
		} catch (Exception e) {
			return false;
		}
	}

	public static String backup(String filename, String sourceDir, String backupDir) {
		return backup(filename, sourceDir, backupDir, true);
	}

	public static String backup(String filename, String sourceDir, String backupDir,
															boolean deleteOriginal) {
		int index = filename.indexOf(".") > 0 ? filename.lastIndexOf(".") : filename.length();
		String root = filename.substring(0, index);
		String exten = filename.substring(index);
		String lastBackup, newBackup, finalBackup;
		int count = 1;

		new File(backupDir).mkdirs();

		while (new File(backupDir + root + "(" + count + ")" + exten).exists()) {
			count++;
		}

		lastBackup = root + "(" + (count - 1) + ")" + exten;
		newBackup = root + "(" + count + ")" + exten;
		if (count > 1 && identicalFiles(sourceDir + filename, backupDir + lastBackup, false)) {
			finalBackup = lastBackup;
			if (deleteOriginal) {
				new File(sourceDir + filename).delete();
			}
		} else {
			if (deleteOriginal) {
				new File(sourceDir + filename).renameTo(new File(backupDir + newBackup));
			} else {
				copyFile(sourceDir + filename, backupDir + newBackup);
			}
			finalBackup = newBackup;
		}

		return finalBackup;
	}

	public static void backupAndMove(String filename, String source, String target, String backup) {
		int index = filename.indexOf(".") > 0 ? filename.lastIndexOf(".") : filename.length();
		String root = filename.substring(0, index);
		String exten = filename.substring(index);
		int count = 1;

		if (!new File(target + filename).exists()
				|| !identicalFiles(source + filename, target + filename, false)) {
			while (new File(backup + root + "(" + count + ")" + exten).exists()) {
				count++;
			}

			new File(target + filename).renameTo(new File(backup + root + "(" + count + ")" + exten));
			new File(source + filename).renameTo(new File(target + filename));
		} else {
			new File(source + filename).delete();
		}
	}

	/**
	 * Searches the current directory and then the alternate location for a file and then opens it
	 * with the appropriate reader
	 *
	 * @param filename the name of the filename to search for
	 * @param alt_locations the alternate directories in which to search for the file
	 * @return a BufferedReader for this file
	 * @throws FileNotFoundException if the file does not exist in any of the locations
	 */
	public static BufferedReader getReader(	String filename,
																					String alt_location) throws FileNotFoundException {
		return getReader(filename, new String[] {alt_location});
	}

	/**
	 * Searches the current directory and then all of the alternate locations until it finds the file
	 * and then opens it with the appropriate reader
	 *
	 * @param filename the name of the filename to search for
	 * @param alt_locations the alternate directories in which to search for the file
	 * @return a BufferedReader for this file
	 * @throws FileNotFoundException if the file does not exist in any of the locations
	 */
	public static BufferedReader getReader(	String filename,
																					String[] alt_locations) throws FileNotFoundException {
		BufferedReader reader;

		reader = null;
		for (int i = -1; reader == null && i < alt_locations.length; i++) {
			reader = getAppropriateReader(i == -1 ? filename : alt_locations[i] + "/" + filename);
		}

		if (reader == null) {
			throw new FileNotFoundException("Error: file \""	+ filename
																			+ "\" not found in current directory, or any of the alternate directories");
		}

		return reader;
	}

	public static BufferedReader getAppropriateReader(String filename) throws FileNotFoundException {
		BufferedReader reader = new BufferedReader(getAppropriateInputStreamReader(filename));
		try {
			if (!reader.ready()) {
				System.err.println("Error - " + filename + " was empty");
			}
		} catch (IOException e) {
			System.err.println("Error accessing '" + filename + "'");
			e.printStackTrace();
		}
		return reader;
	}

	public static InputStreamReader getAppropriateInputStreamReader(String filename) throws FileNotFoundException {
		InputStream is = null;
		InputStreamReader isReader = null;

		if (!exists(filename, false) && !exists(filename, true)) {
			throw new FileNotFoundException("File '" + filename + "' was no where to be found");
		}

		if (!exists(filename, false)) {
			return null;
		}

		if (filename.endsWith(".gz")) {
			try {
				if (!Files.checkJVMUpToDateApprox()) {
					System.err.println("\nYOUR VERSION OF JAVA IS OUT OF DATE; reading gzipped files may fail.");
				}
			} catch (Exception e) {
			}

			try {
				is = new GZIPInputStream(new FileInputStream(filename));
			} catch (IOException e) {
				System.err.println("Error accessing '" + filename + "'");
				e.printStackTrace();
			}
		} else if (filename.endsWith(".zip")) {
			try {
				is = new ZipInputStream(new FileInputStream(filename));
			} catch (IOException e) {
				System.err.println("Error accessing '" + filename + "'");
				e.printStackTrace();
			}
		} else {
			is = new FileInputStream(filename);
		}
		if (is == null) {
			return null;
		}

		if (filename.contains(".utf8.")) {
			try {
				isReader = new InputStreamReader(is, "UTF-8");
			} catch (UnsupportedEncodingException e) {
				e.printStackTrace();
			}
		} else {
			isReader = new InputStreamReader(is);
		}

		return isReader;
	}

	public static PrintWriter getWriter(String filename) {
		return getAppropriateWriter(filename);
	}

	/**
	 * @param filenames initialize writers for these filenames using
	 *        {@link Files#getAppropriateWriter(String)}
	 * @return
	 */
	public static PrintWriter[] getAppropriateWriters(String[] filenames) {
		PrintWriter[] writers = new PrintWriter[filenames.length];
		for (int i = 0; i < writers.length; i++) {
			writers[i] = getAppropriateWriter(filenames[i]);
		}
		return writers;
	}

	/**
	 * @param writers close all of em
	 */
	public static void closeAllWriters(PrintWriter[] writers) {
		for (PrintWriter writer : writers) {
			writer.close();
		}
	}

	public static PrintWriter getAppropriateWriter(String filename) {
		PrintWriter writer = null;

		try {
			if (filename.endsWith(".gz")) {
				writer = new PrintWriter(new GZIPOutputStream(new FileOutputStream(filename)));
			} else if (filename.endsWith(".zip")) {
				writer = new PrintWriter(new ZipOutputStream(new FileOutputStream(filename)));
			} else {
				writer = new PrintWriter(new FileWriter(filename));
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""	+ filename
													+ "\" could not be written to (it's probably open)");
			return null;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			return null;
		}

		return writer;
	}

	public static void mergeFromParameters(String filename, Logger log) {
		List<String> params;
		String file1, file2, mergedFile;
		int lookup1, lookup2;
		int[] indices1, indices2;
		boolean keepRowsUniqueToFile1, keepRowsUniqueToFile2;

		params = Files.parseControlFile(filename, "merge",
																		new String[] {"sourcefile.txt numFiles=6 sizeOfHeader=1 blockSize=1 root=list ext=.dat"},
																		log);
		if (params != null) {
			file1 = params.get(0);
			lookup1 = Integer.parseInt(params.get(1));
			indices1 = Array.toIntArray(params.get(2).split("[\\s]+"));
			keepRowsUniqueToFile1 = params.get(3).trim().toLowerCase().equals("true");
			file2 = params.get(4);
			lookup2 = Integer.parseInt(params.get(5));
			indices2 = Array.toIntArray(params.get(6).split("[\\s]+"));
			keepRowsUniqueToFile2 = params.get(7).trim().toLowerCase().equals("true");
			mergedFile = params.get(8);

			try {
				Files.merge(file1, lookup1, indices1, keepRowsUniqueToFile1, file2, lookup2, indices2,
										keepRowsUniqueToFile2, mergedFile);
			} catch (Elision e) {
				log.reportError("Error merging files '" + file1 + "' and '" + file2 + "'");
				log.reportException(e);
			} catch (Exception e) {
				log.reportError("Error merging files '" + file1 + "' and '" + file2 + "'");
				log.reportException(e);
			}
		}
	}

	public static void merge(	String file1, int lookup1, int[] indices1, boolean keepRowsUniqueToFile1,
														String file2, int lookup2, int[] indices2,
														boolean keepRowsUniqueToFile2, String mergedFile) throws Elision {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
		Hashtable<String, String[]> used = new Hashtable<String, String[]>();
		String[] point;
		String[] keys;
		String trav;

		try {
			reader = new BufferedReader(new FileReader(file2));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				point = new String[indices2.length];
				for (int i = 0; i < point.length; i++) {
					point[i] = line[indices2[i]];
				}
				hash.put(line[lookup2], point);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			throw new Elision("Error: file \"" + file2 + "\" not found in current directory");
		} catch (IOException ioe) {
			ioe.printStackTrace();
			throw new Elision("Error reading file \"" + file2 + "\"");
		}

		writer = Files.getWriter(mergedFile);
		try {
			reader = new BufferedReader(new FileReader(file1));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				trav = "";
				for (int i = 0; i < indices1.length; i++) {
					trav += (i == 0 ? "" : "\t") + line[indices1[i]];
				}
				if (used.containsKey(line[lookup1])) {
					writer.println(trav + "\t" + Array.toStr(used.get(line[lookup1])));
				} else if (hash.containsKey(line[lookup1])) {
					writer.println(trav + "\t" + Array.toStr(hash.get(line[lookup1])));
					used.put(line[lookup1], hash.remove(line[lookup1]));
				} else if (keepRowsUniqueToFile1) {
					point = Array.stringArray(indices2.length, ".");
					if (Ints.indexOf(indices2, lookup2) != -1) {
						point[Ints.indexOf(indices2, lookup2)] = line[lookup1];
					}
					writer.println(trav + "\t" + Array.toStr(point));
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			throw new Elision("Error: file \"" + file1 + "\" not found in current directory");
		} catch (IOException ioe) {
			ioe.printStackTrace();
			throw new Elision("Error reading file \"" + file1 + "\"");
		}
		if (keepRowsUniqueToFile2) {
			keys = HashVec.getKeys(hash);
			for (String key : keys) {
				line = Array.stringArray(indices1.length, ".");
				if (Ints.indexOf(indices1, lookup1) != -1) {
					line[Ints.indexOf(indices1, lookup1)] = key;
				}
				writer.println(Array.toStr(line) + "\t" + Array.toStr(hash.get(key)));
			}
		}

		writer.close();
	}

	public static void mergeSNPlistsFromParameters(String filename, Logger log) {
		List<String> params;
		String file1, file2, mergedFile;
		int lookup1, lookup2;
		int[] indices1;
		boolean keepRowsUniqueToFile1, keepRowsUniqueToFile2;

		params = Files.parseControlFile(filename, "merge",
																		new String[] {"sourcefile.txt numFiles=6 sizeOfHeader=1 blockSize=1 root=list ext=.dat"},
																		log);
		if (params != null) {
			file1 = params.get(0);
			lookup1 = Integer.parseInt(params.get(1));
			indices1 = Array.toIntArray(params.get(2).split("[\\s]+"));
			keepRowsUniqueToFile1 = params.get(3).trim().toLowerCase().equals("true");
			file2 = params.get(4);
			lookup2 = Integer.parseInt(params.get(5));
			if (!params.get(6).equals(lookup2 + "")) {
				log.reportError("FYI - ignoring other indices for file2; with mergeSNPs, only presence is kept");
			}
			keepRowsUniqueToFile2 = params.get(7).trim().toLowerCase().equals("true");
			mergedFile = params.get(8);

			try {
				Files.mergeSNPLists(file1, lookup1, indices1, keepRowsUniqueToFile1, file2, lookup2,
														keepRowsUniqueToFile2, mergedFile);
			} catch (Elision e) {
				log.reportError("Error merging files '" + file1 + "' and '" + file2 + "'");
				log.reportException(e);
			} catch (Exception e) {
				log.reportError("Error merging files '" + file1 + "' and '" + file2 + "'");
				log.reportException(e);
			}
		}
	}

	public static void mergeSNPLists(	String file1, int lookup1, int[] indices1,
																		boolean keepRowsUniqueToFile1, String file2, int lookup2,
																		boolean keepRowsUniqueToFile2,
																		String mergedFile) throws Elision {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String trav;
		int count, num, index;
		int[] rsNums;
		boolean[] used;

		System.out.println("Loading file2 into memory");
		try {
			reader = new BufferedReader(new FileReader(file2));
			count = 0;
			while (reader.ready()) {
				reader.readLine();
				count++;
			}
			reader.close();

			reader = new BufferedReader(new FileReader(file2));
			rsNums = new int[count];
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line[lookup2].startsWith("rs")) {
					try {
						rsNums[count] = Integer.parseInt(line[lookup2].substring(2));
					} catch (NumberFormatException nfe) {
						System.err.println("Error - while it starts with rs, this is an invalid rs number: "
																+ line[lookup2]);
					}
				} else {
					System.err.println("Error - invalid rs number: " + line[lookup2]);
				}
				count++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			throw new Elision("Error: file \"" + file2 + "\" not found in current directory");
		} catch (IOException ioe) {
			ioe.printStackTrace();
			throw new Elision("Error reading file \"" + file2 + "\"");
		}

		System.out.println("Sorting rs numbers");
		Arrays.sort(rsNums);
		used = new boolean[rsNums.length];
		System.out.println("Writing merged filed");
		writer = Files.getWriter(mergedFile);
		try {
			reader = new BufferedReader(new FileReader(file1));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				trav = "";
				for (int i = 0; i < indices1.length; i++) {
					trav += (i == 0 ? "" : "\t") + line[indices1[i]];
				}

				if (line[lookup1].startsWith("rs")) {
					try {
						num = Integer.parseInt(line[lookup1].substring(2));
					} catch (NumberFormatException nfe) {
						System.err.println("Error - while it starts with rs, this is an invalid rs number: "
																+ line[lookup1]);
						num = -1;
					}
				} else {
					System.err.println("Error - invalid rs number: " + line[lookup1]);
					num = -1;
				}

				index = Array.binarySearch(rsNums, num, true);
				if (index >= 0) {
					writer.println(trav + "\t1");
					used[index] = true;
				} else if (keepRowsUniqueToFile1) {
					writer.println(trav + "\t0");
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			throw new Elision("Error: file \"" + file1 + "\" not found in current directory");
		} catch (IOException ioe) {
			ioe.printStackTrace();
			throw new Elision("Error reading file \"" + file1 + "\"");
		}
		if (keepRowsUniqueToFile2) {
			for (int i = 0; i < rsNums.length; i++) {
				if (!used[i]) {
					line = Array.stringArray(indices1.length, ".");
					if (Ints.indexOf(indices1, lookup1) != -1) {
						line[Ints.indexOf(indices1, lookup1)] = "rs" + rsNums[i];
					}
					writer.println(Array.toStr(line) + "\t0");
				}
			}
		}

		writer.close();
	}

	public static void combine(	String[] keys, String[] fileParameters, String unit,
															String outputFilename, Logger log, boolean ignoreCase) {
		combine(keys, fileParameters, new String[fileParameters.length][], unit, ".", outputFilename,
						log, ignoreCase, true, false, null);
	}

	public static void combine(	String[] keys, String[] fileParameters, String[][] headers,
															String unit, String missingValue, String outputFilename, Logger log,
															boolean ignoreCase, boolean finalHeader, boolean hideIndex) {
		combine(keys, fileParameters, headers, unit, missingValue, outputFilename, log, ignoreCase,
						finalHeader, hideIndex, null);
	}

	/**
	 * @param altHeaderMap Map containing keys as filename_headerInFile and values of header in output
	 */
	public static void combine(	String[] keys, String[] fileParameters, String[][] headers,
															String unit, String missingValue, String outputFilename, Logger log,
															boolean ignoreCase, boolean finalHeader, boolean hideIndex,
															Hashtable<String, Hashtable<String, String>> altHeaderMap) {
		PrintWriter writer;
		String[] line, colNames;
		Hashtable<String, String> hash;
		Hashtable<String, String[]> sHash;
		int trav;
		String[][][] data;
		GenParser parser;
		String serializedFilename;
		boolean serializing;
		String delimiter;
		ArrayList<String> fileNames = new ArrayList<String>(fileParameters.length);
		delimiter = Files.suggestDelimiter(outputFilename, log);

		hash = new Hashtable<String, String>();
		for (int i = 0; i < keys.length; i++) {
			hash.put(ignoreCase ? keys[i].toLowerCase() : keys[i], i + "");
		}
		data = new String[fileParameters.length][keys.length + 1][];

		if (log.getLevel() > 8) {
			log.report(ext.addCommas(Runtime.getRuntime().maxMemory()) + " memory available");
		}

		try {
			for (int i = 0; i < fileParameters.length; i++) {
				line = fileParameters[i].trim().split("[\\s]+");
				serializedFilename = GenParser.parseSerialFilename(line);
				fileNames.add(line[0]);
				if (Files.exists(serializedFilename, false)) {
					if (log.getLevel() > 8) {
						log.report("Loading pre-serialized data from '" + line[0] + "'");
					}
					sHash = SerialHash.loadSerializedStringArrayHash(serializedFilename);

					data[i][keys.length] = sHash.get("!colNames");

					for (int j = 0; j < keys.length; j++) {
						data[i][j] = Array.stringArray(data[i][keys.length].length, missingValue);
					}

					for (int k = 0; k < keys.length; k++) {
						if (sHash.containsKey(keys[k])) {
							data[i][k] = sHash.get(keys[k]);
						}
					}
					sHash = null;
				} else {
					if (!Files.exists(line[0])) {
						log.reportError("File " + line[0] + " not found");
						return;
					}
					if (log.getLevel() > 8) {
						log.report("Loading data from '" + line[0] + "'");
					}
					parser = new GenParser(line, log);
					serializing = parser.toBeSerialized();
					if (serializing) {
						sHash = new Hashtable<String, String[]>();
					} else {
						sHash = null;
					}
					if (headers != null	&& headers[i] != null
							&& !ext.checkHeader(parser.getOriginalColumnNames(), headers[i], true, log, false)) {
						log.reportError("Error - unexpected header for file " + line[0]);
						System.exit(1);
					}

					for (int j = 0; j < keys.length; j++) {
						if (parser.forcingFailCodes()) {
							data[i][j] = Array.subArray(parser.getFailCodes(), 1);
						} else {
							data[i][j] = Array.stringArray(parser.getNumColumns() - 1, missingValue);
						}
					}

					colNames = parser.getColumnNames();
					data[i][keys.length] = new String[colNames.length - 1];
					for (int j = 0; j < colNames.length - 1; j++) {
						data[i][keys.length][j] = colNames[j + 1];
					}
					if (serializing) {
						sHash.put("!colNames", data[i][keys.length]);
					}
					while (parser.ready()) {
						line = parser.nextLine();
						if (line != null && hash.containsKey(ignoreCase ? line[0].toLowerCase() : line[0])) {
							trav = Integer.parseInt(hash.get(ignoreCase ? line[0].toLowerCase() : line[0]));
							for (int j = 1; j < line.length; j++) {
								data[i][trav][j - 1] = line[j];
							}
						}
						if (serializing && line != null) {
							sHash.put(line[0], Array.subArray(line, 1));
						}
					}
					parser.reportTruncatedLines();
					parser.close();
					if (serializing) {
						SerialHash.createSerializedStringArrayHash(serializedFilename, sHash);
						sHash = null;
					}
				}
			}
		} catch (OutOfMemoryError oome) {
			log.reportError("Uh oh! Ran out of memory while combining files!");
		} catch (Exception e) {
			log.reportException(e);
		}

		try {
			writer = new PrintWriter(new FileWriter(outputFilename));
			if (finalHeader) {
				if (!hideIndex) {
					writer.print(unit);
				}
				for (int i = 0; i < data.length; i++) {
					try {
						if (altHeaderMap != null) {
							String[] tmp = new String[data[i][keys.length].length];
							for (int j = 0; j < data[i][keys.length].length; j++) {
								String f = ext.removeDirectoryInfo(fileNames.get(i));
								if (altHeaderMap.containsKey(f)
										&& altHeaderMap.get(f).contains(data[i][keys.length][j])) {
									tmp[j] = altHeaderMap.get(f).get(data[i][keys.length][j]);
								} else {
									tmp[j] = data[i][keys.length][j];
								}
							}
							writer.print((hideIndex && i == 0 ? "" : delimiter) + Array.toStr(tmp, delimiter));
						} else {
							writer.print((hideIndex && i == 0 ? "" : delimiter)
														+ Array.toStr(data[i][keys.length], delimiter));
						}
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
				writer.println();
			}
			for (int j = 0; j < keys.length; j++) {
				if (!hideIndex) {
					writer.print(keys[j]);
				}
				for (int i = 0; i < data.length; i++) {
					writer.print((hideIndex && i == 0 ? "" : delimiter) + Array.toStr(data[i][j], delimiter));
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outputFilename);
			log.reportException(e);
		}
	}

	public static void combineWithLessMemory(	String[] keys, String[] fileParameters,
																						String[][] headers, String unit, String missingValue,
																						String outputFilename, Logger log, boolean ignoreCase,
																						boolean finalHeader, boolean hideIndex,
																						boolean keepIntermediateFiles) {
		BufferedReader[] readers;
		PrintWriter writer;
		String[] line, colNames;
		Hashtable<String, String> hash;
		int trav;
		String[][] data;
		GenParser parser;
		String delimiter;
		int counter;
		Hashtable<String, String[]> sHash;
		boolean serializing;
		String serializedFilename;

		delimiter = Files.suggestDelimiter(outputFilename, log);

		log.report("Data will be parsed for each file separately and then merged");
		hash = new Hashtable<String, String>();
		for (int i = 0; i < keys.length; i++) {
			hash.put(ignoreCase ? keys[i].toLowerCase() : keys[i], i + "");
		}

		line = null;
		for (int i = 0; i < fileParameters.length; i++) {
			try {
				data = new String[keys.length + 1][];

				line = fileParameters[i].trim().split(ext.REGEX_TO_SPLIT_SPACES_NOT_IN_QUOTES);
				line[0] = line[0].replaceAll("\"", "");

				serializedFilename = GenParser.parseSerialFilename(line);
				if (Files.exists(serializedFilename, false)) {
					if (log.getLevel() > 8) {
						log.report("Loading pre-serialized data from '" + line[0] + "'");
					}
					sHash = SerialHash.loadSerializedStringArrayHash(serializedFilename);

					data[0] = sHash.get("!colNames");

					for (int j = 0; j < keys.length; j++) {
						data[j + 1] = Array.stringArray(data[0].length, missingValue);
					}

					for (int k = 0; k < keys.length; k++) {
						if (sHash.containsKey(keys[k])) {
							data[k + 1] = sHash.get(keys[k]);
						}
					}
					sHash = null;

					Files.writeMatrix(data, "file." + (i + 1) + ".temp", delimiter);
				} else {
					if (!Files.exists(line[0])) {
						log.reportError("File " + line[0] + " not found");
						return;
					}
					if (log.getLevel() > 8) {
						log.report("Parsing data from '" + line[0] + "'");
					}

					parser = new GenParser(line, log);

					serializing = parser.toBeSerialized();
					if (serializing) {
						sHash = new Hashtable<String, String[]>();
					} else {
						sHash = null;
					}
					if (headers != null	&& headers[i] != null
							&& !ext.checkHeader(parser.getOriginalColumnNames(), headers[i], true, log, false)) {
						log.reportError("Error - unexpected header for file " + line[0]);
						System.exit(1);
					}

					for (int j = 0; j < keys.length; j++) {
						if (parser.forcingFailCodes()) {
							data[j + 1] = Array.subArray(parser.getFailCodes(), 1);
						} else {
							data[j + 1] = Array.stringArray(parser.getNumColumns() - 1, missingValue);
						}
					}

					colNames = parser.getColumnNames();
					data[0] = new String[colNames.length - 1];
					for (int j = 0; j < colNames.length - 1; j++) {
						data[0][j] = colNames[j + 1];
					}

					if (serializing) {
						sHash.put("!colNames", data[keys.length]);
					}
					while (parser.ready()) {
						line = parser.nextLine();
						if (line != null && hash.containsKey(ignoreCase ? line[0].toLowerCase() : line[0])) {
							trav = Integer.parseInt(hash.get(ignoreCase ? line[0].toLowerCase() : line[0]));
							for (int j = 1; j < line.length; j++) {
								data[trav + 1][j - 1] = line[j];
							}
						}
						if (serializing && line != null) {
							sHash.put(line[0], Array.subArray(line, 1));
						}
					}
					parser.reportTruncatedLines();
					parser.close();
					if (serializing) {
						SerialHash.createSerializedStringArrayHash(serializedFilename, sHash);
						sHash = null;
					}
					Files.writeMatrix(data, "file." + (i + 1) + ".temp", delimiter);
				}
			} catch (OutOfMemoryError oome) {
				log.reportError("Uh oh! Ran out of memory parsing file '"
												+ fileParameters[i].trim().split("[\\s]+")[0] + "' at line:");
				log.reportError(Array.toStr(line, "/"));
			} catch (Exception e) {
				log.reportException(e);
			}
		}

		try {
			writer = new PrintWriter(new FileWriter(outputFilename));
			readers = new BufferedReader[fileParameters.length];
			for (int i = 0; i < readers.length; i++) {
				readers[i] = new BufferedReader(new FileReader("file." + (i + 1) + ".temp"));
				if (finalHeader) {
					if (!hideIndex && i == 0) {
						writer.print(unit);
					}
					writer.print((hideIndex && i == 0 ? "" : delimiter) + readers[i].readLine());
				} else {
					readers[i].readLine();
				}
			}
			if (finalHeader) {
				writer.println();
			}

			counter = 0;
			while (readers[0].ready()) {
				for (int i = 0; i < readers.length; i++) {
					if (!hideIndex && i == 0) {
						writer.print(keys[counter]);
					}
					writer.print((hideIndex && i == 0 ? "" : delimiter) + readers[i].readLine());
				}
				writer.println();
				counter++;
			}
			writer.close();
			for (int i = 0; i < readers.length; i++) {
				readers[i].close();
				if (!keepIntermediateFiles) {
					new File("file." + (i + 1) + ".temp").delete();
				}
			}
		} catch (Exception e) {
			log.reportError("Error writing to " + outputFilename);
			log.reportException(e);
		}
	}


	/**
	 * Looks up the values for a set of keys in memory
	 *
	 * @param keys - the keys to lookup (if they contain tabs, then they will be matched on the first
	 *        N set of columns in the lookup values)
	 * @param lookupValues - Matrix of String with the lookup values
	 * @param missingValue - the String to use if the key is not present in the lookup values
	 * @param ignoreCase - ignore case when matching keys to values
	 * @param hideIndex - do not return the keys as the first column of the result
	 * @param log - the Logger to report any errors
	 * @return
	 * @throws Elision
	 */
	public static String[][] combineInMemory(	String[] keys, String[][] lookupSource,
																						String missingValue, boolean ignoreCase,
																						boolean hideIndex, Logger log) throws Elision {
		Hashtable<String, Integer> hash;
		int numColsInKey, numColsInValues;
		String[][] data;
		String key;
		int count;

		numColsInKey = keys[0].split("\t", -1).length;
		numColsInValues = lookupSource[0].length;

		if (numColsInKey > 1) {
			log.report("The keys contain tabs, so the lookup will be performed on the first "
									+ numColsInKey + " columns in the lookup source data");
		}

		count = 0;
		hash = new Hashtable<String, Integer>(keys.length);
		for (int i = 0; i < lookupSource.length; i++) {
			if (lookupSource[i].length > 0) {
				if (lookupSource[i].length < numColsInKey + 1) {
					throw new Elision("The lookup values require at least "	+ (numColsInKey + 1)
														+ " columns (" + numColsInKey + " for the key and at least one value"
														+ "), but row " + (i + 1) + " has only " + lookupSource[i].length);
				}
				if (lookupSource[i].length != numColsInValues) {
					throw new Elision("The first row of the lookup values had "	+ numColsInValues
														+ " columns while row " + (i + 1) + " has "
														+ keys[i].split("\t", -1).length);
				}

				key = lookupSource[i][0];
				for (int j = 1; j < numColsInKey; j++) {
					key += "\t" + lookupSource[i][j];
				}
				if (ignoreCase) {
					key = key.toLowerCase();
				}

				if (hash.containsKey(key)) {
					if (count < 5 && log.getLevel() > 8) {
						log.reportError("Warning duplicate key in lookup values ('"	+ key
														+ "'); only the first value will be used:");
						log.reportError(Array.toStr(lookupSource[hash.get(key)]));
						log.reportError("...and not this one:");
						log.reportError(Array.toStr(lookupSource[i]));
					} else if (count == 5) {
						log.reportError("...");
					}
					count++;
				}
				hash.put(key, i);
			}
		}
		if (count > 5) {
			log.reportError("There were " + count + " instances of duplicate keys in the lookup values");
		}
		data = new String[keys.length][];

		if (log.getLevel() > 8) {
			log.report(ext.addCommas(Runtime.getRuntime().maxMemory()) + " memory available");
		}

		count = 0;
		try {
			for (int i = 0; i < keys.length; i++) {
				key = ignoreCase ? keys[i].toLowerCase() : keys[i];
				if (key.split("\t", -1).length != numColsInKey) {
					throw new Elision("The first row of the keys had "	+ numColsInKey + " column"
														+ (numColsInKey == 1 ? "" : "s") + " while row " + (i + 1) + " has "
														+ keys[i].split("\t", -1).length);
				}

				if (hash.containsKey(key)) {
					data[i] = lookupSource[hash.get(key)];
				} else {
					data[i] = Array.stringArray(numColsInValues, missingValue);
					count++;
				}
			}
		} catch (OutOfMemoryError oome) {
			log.reportError("Uh oh! Ran out of memory while combining files!");
			return null;
		} catch (Exception e) {
			log.reportException(e);
			return null;
		}

		if (count > 0 && log.getLevel() > 8) {
			log.report("\nOf the "	+ data.length + " keys, " + count + " " + (count == 1 ? "was" : "were")
									+ " not found among the lookup values\n");

			if (data.length == count) {
				log.report("The first few keys were: ");
				for (int i = 0; i < Math.min(keys.length, 10); i++) {
					log.report(keys[i]);
				}
				log.report("");
				log.report("The first few lookup values were: ");
				for (int i = 0; i < Math.min(lookupSource.length, 10); i++) {
					key = lookupSource[i][0];
					for (int j = 1; j < numColsInKey; j++) {
						key += "\t" + lookupSource[i][j];
					}
					log.report(key);
				}
			}
		}


		return data;
	}

	public static void filterByKeys(String keysFile, String fileIn, String fileOut, int keyColumn,
																	boolean preserveKeyOrder) throws IOException {
		if (!Files.exists(keysFile)) {
			System.err.println("Error - key file " + keysFile + " doesn't exist!");
			return;
		}
		if (!Files.exists(fileIn)) {
			System.err.println("Error - input data file " + fileIn + " doesn't exist!");
			return;
		}
		if (Files.exists(fileOut)) {
			System.err.println("Error - output file " + fileOut + " already exists!");
			return;
		}
		if (keyColumn < 0) {
			System.err.println("Error - invalid key column index: " + keyColumn);
			return;
		}

		BufferedReader reader;
		PrintWriter writer;
		String line, delim;

		String[] keys = HashVec.loadFileToStringArray(keysFile, false, new int[] {0}, false);
		HashMap<String, Integer> keyMap = new HashMap<String, Integer>();
		for (int i = 0; i < keys.length; i++) {
			keyMap.put(keys[i], i);
		}
		String[] lines = new String[keys.length];

		reader = Files.getAppropriateReader(fileIn);
		writer = preserveKeyOrder ? null : Files.getAppropriateWriter(fileOut);
		line = delim = null;
		while ((line = reader.readLine()) != null) {
			if (delim == null) {
				delim = ext.determineDelimiter(line);
			}
			String key = line.split(delim, -1)[keyColumn];
			if (keyMap.containsKey(key)) {
				if (!preserveKeyOrder) {
					writer.println(line);
				} else {
					lines[keyMap.get(key)] = line;
				}
			}
		}
		if (!preserveKeyOrder) {
			writer.flush();
			writer.close();
		}
		reader.close();
		if (preserveKeyOrder) {
			lines = Array.removeMissingValues(lines);
			Files.writeArray(lines, fileOut);
		}
	}

	public static void generateTables(String outputFile, String[] files, String[] fileDescriptions,
																		String[][] parameters, Logger log) {
		PrintWriter writer;
		String[] values;
		String[][] data;
		Vector<String> filters;
		Vector<Vector<String>> output;
		double[] array, means;
		int[] counts;
		boolean stdev, blank, percent;
		int sf, percentile;

		// delimiters = new String[files.length];
		// headers = new String[files.length][];
		output = new Vector<Vector<String>>();
		for (String file : files) {
			// headers[i] = getHeaderOfFile(files[i], delimiters[i], log);
			output.add(new Vector<String>());
		}

		try {
			writer = new PrintWriter(new FileWriter(outputFile));
			writer.println("\t" + Array.toStr(fileDescriptions) + "\tAll");
			for (String[] parameter : parameters) {
				writer.print(parameter[0]);
				stdev = false;
				blank = false;
				percent = false;
				sf = 4;
				values = new String[] {parameter[2]};
				filters = new Vector<String>();
				if (parameter[1].equals("mean") || parameter[1].equals("stdev")) {
					// stdev = false;
					// blank = false;
					// percent = false;
					// sf = 4;
					// values = new String[] {parameters[i][2]};
					// filters = new Vector<String>();
					for (int j = 3; j < parameter.length; j++) {
						if (parameter[j].equals("-stdev")) {
							stdev = true;
						} else if (parameter[j].startsWith("-sf=")) {
							sf = ext.parseIntArg(parameter[j]);
						} else if (parameter[j].equals("-blank")) {
							blank = true;
						} else if (parameter[j].equals("-percent")) {
							percent = true;
						} else {
							filters.add(parameter[j]);
						}
					}
					means = new double[files.length + 1];
					counts = new int[files.length];
					for (int j = 0; j < files.length; j++) {
						data = generateDataset(	files[j], determineDelimiter(files[j], log), values,
																		Array.toStringArray(filters), log);
						array = Array.toDoubleArray(Matrix.extractColumn(data, 0));
						means[j] = Array.mean(array);
						counts[j] = data.length;
						if (counts[j] > 0) {
							means[files.length] += means[j] * counts[j];
						}
						if (parameter[1].equals("mean")) {
							writer.print("\t" + (counts[j] > 0
																										? (percent	? ext.formDeci(means[j] * 100, sf)
																																+ "%"
																															: ext.formDeci(	means[j],
																																							sf)
																																+ (stdev	? " (+/- "
																																							+ ext.formDeci(	Array.stdev(array),
																																														sf)
																																						+ ")"
																																					: ""))
																									: (blank ? "" : ".")));
						} else {
							writer.print("\t" + (counts[j] > 0	? ext.formDeci(Array.stdev(array), sf)
																									: (blank ? "" : ".")));
						}
					}
					if (parameter[1].equals("mean")) {
						writer.print("\t" + (Array.sum(counts) > 0
																													? (percent	? ext.formDeci(means[files.length]
																																											/ Array.sum(counts)
																																										* 100, sf)
																																			+ "%"
																																		: ext.formDeci(means[files.length]
																																											/ Array.sum(counts),
																																										sf))
																												: (blank ? "" : ".")));
					} else {
						// TODO calculate the overall stdev of crossing files.
					}

				} else if (parameter[1].equals("count")) {
					// blank = false;
					// values = new String[] {parameters[i][2]};
					// filters = new Vector<String>();
					for (int j = 3; j < parameter.length; j++) {
						if (parameter[j].equals("-blank")) {
							blank = true;
						} else {
							filters.add(parameter[j]);
						}
					}
					counts = new int[files.length];
					for (int j = 0; j < files.length; j++) {
						data = generateDataset(	files[j], determineDelimiter(files[j], log), values,
																		Array.toStringArray(filters), log);
						counts[j] = data.length;
						writer.print("\t" + (counts[j] == 0 ? (blank ? "" : "0") : counts[j]));
					}
					writer.print("\t" + (Array.sum(counts) == 0 ? (blank ? "" : "0") : Array.sum(counts)));

				} else if (parameter[1].contains("percentile")) {
					percentile = Integer.parseInt(parameter[1].split(" ")[0].trim());
					for (int j = 3; j < parameter.length; j++) {
						if (parameter[j].equals("-blank")) {
							blank = true;
						} else if (parameter[j].startsWith("-sf=")) {
							sf = ext.parseIntArg(parameter[j]);
						} else if (parameter[j].equals("-percent")) {
							percent = true;
						} else {
							filters.add(parameter[j]);
						}
					}
					counts = new int[files.length];
					for (int j = 0; j < files.length; j++) {
						data = generateDataset(	files[j], determineDelimiter(files[j], log), values,
																		Array.toStringArray(filters), log);
						array = Array.toDoubleArray(Matrix.extractColumn(data, 0));
						counts[j] = data.length;
						if (counts[j] > 0) {
							Arrays.sort(array);
						}
						if (percentile == 100) {
							writer.print("\t" + (counts[j] == 0 ? (blank ? "" : "0") : array[array.length - 1]));
						} else {
							writer.print("\t" + (counts[j] == 0	? (blank ? "" : "0")
																									: array[(int) (((double) percentile / 100)
																																	* array.length + .5)]));
						}
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outputFile);
			log.reportException(e);
		}

	}

	// consider using this function in db.comp (i.e. .ctl runs)
	public static String[][] generateDataset(	String filename, String delimiter, String[] values,
																						String[] filters, Logger log) {
		BufferedReader reader;
		String[] line, header, filterTargets, data;
		int[] valueIndices, filterIndices;
		Vector<String[]> v;
		boolean include;

		v = new Vector<String[]>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			header = reader.readLine().trim().split(delimiter);
			valueIndices = ext.indexFactors(values, header, true, log, true, true);
			filterIndices = new int[filters.length];
			filterTargets = new String[filters.length];
			for (int i = 0; i < filters.length; i++) {
				filterIndices[i] = ext.indexOfStr(filters[i].split("=")[0], header, true, true, log, true);
				filterTargets[i] = filters[i].split("=")[1];
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter);
				include = true;
				for (int i = 0; i < filterIndices.length; i++) {
					if (filterTargets[i].equalsIgnoreCase("validDouble")) {
						if (!ext.isValidDouble(line[filterIndices[i]])) {
							include = false;
						}
					} else if (filterTargets[i].equalsIgnoreCase("validInteger")) {
						if (!ext.isValidDouble(line[filterIndices[i]])) {
							include = false;
						}
					} else if (!line[filterIndices[i]].equals(filterTargets[i])) {
						include = false;
					}
				}
				if (include) {
					data = new String[values.length];
					for (int i = 0; i < values.length; i++) {
						data[i] = line[valueIndices[i]];
					}
					v.add(data);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		return Matrix.toStringArrays(v);
	}

	public static void transpose(	String filein, String delimiterIn, String fileout,
																String delimiterOut, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[][] matrix;
		int size = -1;
		int lineCount;
		long time;

		time = new Date().getTime();
		try {
			log.report("Counting lines in " + filein);
			lineCount = Files.countLines(filein, 0);
			log.report("Line count = " + lineCount);
			log.report("Loading " + filein + " into memory");
			reader = new BufferedReader(new FileReader(filein));
			matrix = new String[lineCount][];
			for (int i = 0; i < matrix.length; i++) {
				if (lineCount > 20 && i % (matrix.length / 20) == 0) {
					log.report(ext.getTime()	+ "\t" + Math.round(100 * (float) i / matrix.length)
											+ "% loaded; ", false, true);
					log.memoryPercentFree();
				}
				matrix[i] = reader.readLine().trim().split(delimiterIn);
				if (size == -1) {
					size = matrix[i].length;
				} else if (matrix[i].length != size) {
					log.reportError("Error transposing file: different number of columns (previously: "	+ size
													+ "; line " + (i + 1) + ": " + matrix[i].length + ")");
					log.reportError("   Make sure that the delimiter is set correctly");
				}
			}
			reader.close();

			if (fileout == null) {
				fileout = filein + "-transposed.xln";
			}

			log.report("Writing out to " + fileout);
			writer = new PrintWriter(new FileWriter(fileout));
			for (int j = 0; j < size; j++) {
				if (size > 20 && j % (size / 20) == 0) {
					log.report(ext.getTime()	+ "\t" + Math.round(100 * (float) j / size) + "% loaded; ",
											false, true);
					log.memoryPercentFree();
				}
				for (int i = 0; i < matrix.length; i++) {
					writer.print((i == 0 ? "" : delimiterOut) + matrix[i][j]);
				}
				writer.println();
			}

			writer.close();

		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filein + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filein + "\"");
			return;
		} catch (OutOfMemoryError oome) {
			log.reportError("Ran out of memory");
			return;
		}
		log.report("Finished in " + ext.getTimeElapsed(time));
	}

	public static void transposeHuge(String filename, int step, Logger log) {
		transposeHuge(filename, Files.determineDelimiter(filename, log),
									filename + "-transposeHuge.xln", "\t", step, log);
	}

	public static void transposeHuge(	String filename, String delimiterIn, String fileout,
																		String delimiterOut, int step, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String[][] matrix;
		int size = 1, inc = 0;
		int lineCount;

		log.report("Counting lines in " + filename);
		lineCount = Files.countLines(filename, 0);
		matrix = new String[lineCount][];
		try {
			writer = new PrintWriter(new FileWriter(fileout));
			while (inc < size) {
				log.report(ext.getTime()	+ "\t" + Math.round(100 * (float) inc / size) + "% processed; ",
										false, true);
				log.memoryPercentFree();

				reader = new BufferedReader(new FileReader(filename));
				for (int i = 0; i < matrix.length; i++) {
					line = reader.readLine().trim().split(delimiterIn);
					if (size == 1) {
						size = line.length;
					} else if (line.length != size) {
						log.reportError("Error transposing file: different number of columns (previously: "
														+ size + "; line " + (i + 1) + ": " + line.length + ")");
						log.reportError("   Make sure that the delimiter is set correctly");
					}
					matrix[i] = Array.subArray(line, inc, inc + step > size ? size : inc + step);
				}
				reader.close();

				log.report(ext.getTime() + "\tWriting");
				for (int j = 0; j < (inc + step > size ? size - inc : step); j++) {
					for (int i = 0; i < matrix.length; i++) {
						writer.print((i == 0 ? "" : delimiterOut) + matrix[i][j]);
					}
					writer.println();
					writer.flush();
				}
				inc += step;
			}

			writer.close();

		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			return;
		} catch (OutOfMemoryError oome) {
			log.reportError("Ran out of memory");
			return;
		}
	}

	public static void transposeFromParameters(String filename, Logger log) {
		List<String> params;
		String infile, outfile;
		String[] line;
		boolean huge;
		int step = 1000;
		long memoryAvailable, filesize;
		int colCount, suggestedStep;

		params = Files.parseControlFile(filename, "transpose",
																		new String[] {"sourcefile.txt out=newfile.txt forceHuge=false forceStep=-1",
																									"# (step is only used for huge"},
																		log);
		if (params != null) {
			line = params.remove(0).trim().split("[\\s]+");

			huge = false;
			infile = line[0];
			outfile = ext.rootOf(filename) + "_transposed.xln";
			for (int i = 1; i < line.length; i++) {
				if (line[i].startsWith("out=")) {
					outfile = line[i].split("=")[1];
				} else if (line[i].startsWith("forceHuge=")) {
					huge = ext.parseBooleanArg(line[i]);
				} else if (line[i].startsWith("forceStep=")) {
					step = ext.parseIntArg(line[i]);
					System.out.println("Set step to " + step);
				} else {
					log.reportError("Error - don't know what to do with argument: " + line[i]);
				}
			}

			memoryAvailable = log.memoryMax();
			filesize = new File(infile).length();
			log.report("The size of the file is " + ext.prettyUpSize(filesize, 1));
			if (huge || memoryAvailable < filesize * 20) {
				log.report("Using the (slower) huge file transpose given these constraints");

				colCount = Files.getHeaderOfFile(infile, log).length;
				log.report("There are " + colCount + " columns in the file");
				suggestedStep = (int) ((double) memoryAvailable	/ 40
																/ ((double) filesize / (double) colCount));
				log.report("( "	+ memoryAvailable + " / 40 ) / ( " + filesize + " / " + colCount + " ) = "
										+ suggestedStep);
				if (step == -1) {
					step = suggestedStep;
				}

				log.report("Using the huge transpose option with a step of "	+ step + " which will take "
										+ (int) (Math.ceil((double) colCount / (double) step)) + " cycles");
				transposeHuge(infile, Files.determineDelimiter(infile, log), outfile,
											Files.suggestDelimiter(outfile, log), step, log);
			} else {
				transpose(infile, Files.determineDelimiter(infile, log), outfile,
									Files.suggestDelimiter(outfile, log), log);
			}
		}
	}

	public static void extractColumns(String filename, int[] columns) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String delimiter;

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename + "-extracted.xln"));
			delimiter = determineDelimiter(filename, new Logger());
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter);
				for (int i = 0; i < columns.length; i++) {
					writer.print((i == 0 ? "" : "\t") + line[columns[i]]);
				}
				writer.println();
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	/**
	 * Searches all of the directories in the array to see if any of them exist
	 *
	 * @param dirs the array of directories to search
	 * @param verbose whether to report an error if none of the locations exists
	 * @param kill whether to System.exit() if no valid directory is found
	 * @return String first directory that exists, otherwise null
	 */
	public static String firstDirectoryThatExists(String[] dirs, boolean verbose, boolean kill,
																								Logger log) {
		for (int i = 0; i < dirs.length; i++) {
			dirs[i] = ext.verifyDirFormat(dirs[i]);
			if (new File(dirs[i]).exists()) {
				return dirs[i];
			}
		}

		if (verbose) {
			log.reportError("Error - none of the following directories exists..."
											+ (kill ? " terminating" : ""));
			for (String dir : dirs) {
				log.reportError("   " + dir);
			}
		}
		if (kill) {
			System.exit(1);
		}

		return null;
	}

	/**
	 * Searches all of the directories in the array to see if it contains the specified file
	 *
	 * @param filename the filename to search for
	 * @param dirs the array of directories to search
	 * 
	 * @return String the full path to the file of interest if it exists in one of the directories,
	 *         otherwise null
	 */
	public static String firstPathToFileThatExists(String filename, String... dirs) {
		return firstPathToFileThatExists(dirs, filename, false, false, new Logger());
	}

	/**
	 * Searches all of the directories in the array to see if it contains the specified file
	 *
	 * @param dirs the array of directories to search
	 * @param filename the filename to search for
	 * @param verbose whether to report an error if none of the locations exists or if none of the
	 *        locations contains the specified file
	 * @param kill whether to System.exit() if no valid directory is found
	 * @return String the full path to the file of interest if it exists in one of the directories,
	 *         otherwise null
	 */
	public static String firstPathToFileThatExists(	String[] dirs, String filename, boolean verbose,
																									boolean kill, Logger log) {
		for (int i = 0; i < dirs.length; i++) {
			dirs[i] = ext.verifyDirFormat(dirs[i]);
			if (new File(dirs[i] + filename).exists()) {
				return dirs[i] + filename;
			}
		}

		if (verbose) {
			log.reportError("Error - none of the following directories contains the file '"	+ filename
											+ "' " + (kill ? "... terminating" : ""));
			for (String dir : dirs) {
				log.reportError("   " + dir);
			}
		}
		if (kill) {
			System.exit(1);
		}

		return null;
	}

	public static BufferedReader getReader(	String filename, boolean jar, boolean verbose,
																					boolean kill) {
		return getReader(filename, jar, verbose, new Logger(), kill);
	}

	public static BufferedReader getReader(	String filename, boolean jar, boolean verbose, Logger log,
																					boolean kill) {
		try {
			if (Files.exists(filename, jar)) {
				if (jar) {
					return new BufferedReader(new InputStreamReader(ClassLoader.getSystemResourceAsStream(filename)));
				} else {
					return getAppropriateReader(filename);
				}
			} else {
				if (verbose) {
					log.reportError("Error - file not found: " + filename);
				}
				if (kill) {
					System.exit(1);
				}
			}
		} catch (Exception e) {
			if (verbose) {
				log.reportException(e);
			}
		}
		if (kill) {
			System.exit(1);
		}
		return null;
	}

	public static int getSize(String filename, boolean jar) {
		int size = -1;

		try {
			if (Files.exists(filename, jar)) {
				if (jar) {
					size = ClassLoader.getSystemResourceAsStream(filename).available();
				} else {
					size = (int) new File(filename).length();
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		return size;
	}

	public static double[] getSizeScaled(String filename, boolean jar) {
		int sz = getSize(filename, jar);
		if (sz == -1) {
			return new double[] {sz, -1};
		}

		double temp = sz, dec = 0;
		while (temp > 1024) {
			temp = temp / 1024;
			dec++;
		}
		return new double[] {temp, dec};
	}

	public static String getSizeScaledString(String filename, boolean jar) {
		String[] bases = {"b", "kb", "Mb", "Gb", "Tb", "Eb"};

		double[] sz = getSizeScaled(filename, jar);
		if (sz[0] == -1) {
			return "-1";
		} else {
			return ext.formDeci(sz[0], 3) + bases[(int) sz[1]];
		}
	}

	public static boolean exists(String dir, String[] filenames) {
		return exists(dir, filenames, false);
	}

	public static boolean exists(String dir, String[] filenames, boolean treatEmptyAsMissing) {
		boolean result;

		result = true;
		for (String filename : filenames) {
			if (!exists(dir + filename, false, treatEmptyAsMissing)) {
				result = false;
				break;
			}
		}

		return result;
	}

	public static boolean exists(String dir, Iterable<String> filenames) {
		return exists(dir, filenames, false);
	}

	public static boolean exists(	String dir, Iterable<String> filenames,
																boolean treatEmptyAsMissing) {
		boolean result;

		result = true;
		for (String filename : filenames) {
			if (!exists(dir + filename, false, treatEmptyAsMissing)) {
				result = false;
				break;
			}
		}

		return result;
	}

	public static boolean exists(String filename) {
		return exists(filename, false);
	}

	public static boolean exists(String filename, boolean jar) {
		return exists(filename, jar, false);
	}

	public static boolean exists(String filename, boolean jar, boolean treatEmptyAsMissing) {
		if (filename == null) {
			return false;
		}
		if (jar) {
			try {
				ClassLoader.getSystemResourceAsStream(filename).close();
				return !treatEmptyAsMissing || getSize(filename, jar) != 0;
			} catch (Exception e) {
				return false;
			}
		} else {
			return new File(filename).exists() && (!treatEmptyAsMissing || getSize(filename, jar) != 0);
		}
	}

	/**
	 * @param dir directory, set to "" if full paths
	 * @param filenames array of filenames to check...all will be checked
	 * @param verbose report an error to the log for each filename missing
	 * @param log
	 * @return true if all files exist
	 */
	public static boolean checkAllFiles(String dir, String[] filenames, boolean verbose, Logger log) {
		boolean result = true;
		for (int i = 0; i < filenames.length; i++) {
			if (!Files.exists(dir + filenames[i])) {
				if (verbose & log != null) {
					log.reportError("Error - could not find file " + dir + filenames[i]);
				}
				result = false;
			}
		}
		return result;
	}

	public static boolean isDirectory(String handle) {
		return exists(handle) && new File(handle).isDirectory();
	}

	public static boolean checkJVMUpToDateApprox() {
		int currVer = 8;
		String[] version = System.getProperty("java.version").split("\\.");
		if (version[0].equals("1") && Integer.parseInt(version[1]) < currVer) {
			return false;
		}
		return true;
	}

	public static boolean isWindows() {
		return System.getProperty("os.name").startsWith("Windows");
	}

	public static String[] listFullPaths(String directory, final String suffix, boolean jar) {
		return list(directory, null, suffix, false, jar, true);
	}

	public static String[] list(String directory, final String suffix, boolean jar) {
		return list(directory, null, suffix, false, jar);
	}

	public static String[] list(String directory, final String prefix, final String suffix,
															final boolean caseSensitive, boolean jar) {
		return list(directory, prefix, suffix, caseSensitive, jar, false);
	}

	// These variables need to be final in order to work in the FilenameFilter
	// The only illegal character in all operating systems is the colon :
	// so this was chosen to signify NOT
	public static String[] list(String directory, final String prefix, final String suffix,
															final boolean caseSensitive, boolean jar, final boolean fullPath) {
		if (directory == null || directory.length() == 0) {
			directory = "./";
		}

		if (jar) {
			try {
				// System.err.println("I haven't been able to get listFiles() to work inside a jar file
				// yet");

				URL repositoryURL = ClassLoader.getSystemResource("common/Files.class");
				String repositoryPath = repositoryURL.getPath();
				URI jarURI = new URI(repositoryPath.substring(0, repositoryPath.indexOf('!')));
				JarFile jarFile = new JarFile(new File(jarURI));
				Vector<String> v = new Vector<String>();

				Enumeration<JarEntry> entries = jarFile.entries();
				while (entries.hasMoreElements()) {
					String entryName = entries.nextElement().getName();
					boolean passes = true;

					if (!entryName.startsWith(directory)) {
						passes = false;
					}

					if (prefix != null && !prefix.equals("")) {
						if (prefix.startsWith(":")) {
							if (entryName.toLowerCase().startsWith(prefix.toLowerCase())) {
								passes = false;
							}
						} else if (!entryName.toLowerCase().startsWith(prefix.toLowerCase())) {
							passes = false;
						}
					}

					if (suffix != null && !suffix.equals("")) {
						if (suffix.startsWith(":")) {
							if (entryName.toLowerCase().endsWith(suffix.toLowerCase())) {
								passes = false;
							}
						} else if (!entryName.toLowerCase().endsWith(suffix.toLowerCase())) {
							passes = false;
						}
					}

					if (passes) {
						String trav = entryName.substring(directory.length());
						if (trav.startsWith("/")) {
							trav = trav.substring(1);
						}
						if (!trav.contains("/")) {
							if (fullPath) {
								v.add(directory + trav);
							} else {
								v.add(trav);
							}
						}
					}
				}

				jarFile.close();

				return Array.toStringArray(v);
			} catch (Exception e) {
				System.err.println("Error reading files in jar file");
				e.printStackTrace();
				return null;
			}
		} else {
			String[] files;

			files = new File(directory).list(new FilenameFilter() {
				@Override
				public boolean accept(File file, String filename) {
					boolean passes = true;
					String pre, suf;

					if (new File(file, filename).isDirectory()) {
						return false;
					}

					pre = prefix == null ? null : (caseSensitive ? prefix : prefix.toLowerCase());
					suf = suffix == null ? null : (caseSensitive ? suffix : suffix.toLowerCase());
					filename = caseSensitive ? filename : filename.toLowerCase();

					if (pre != null && !pre.equals("")) {
						if (pre.startsWith(":")) {
							if (filename.toLowerCase().startsWith(pre.substring(1).toLowerCase())) {
								passes = false;
							}
						} else if (!filename.toLowerCase().startsWith(pre.toLowerCase())) {
							passes = false;
						}
					}

					if (suf != null && !suf.equals("")) {
						if (suf.startsWith(":")) {
							if (filename.toLowerCase().endsWith(suf.substring(1).toLowerCase())) {
								passes = false;
							}
						} else if (!filename.toLowerCase().endsWith(suf.toLowerCase())) {
							passes = false;
						}
					}

					return passes;
				}
			});

			if (fullPath) {
				for (int i = 0; i < files.length; i++) {
					files[i] = directory + files[i];
				}
			}

			if (files == null) {
				return new String[0];
			} else {
				return files;
			}
		}
	}

	public static String[] listDirectories(String directory, boolean jar) {
		if (jar) {
			try {
				URL repositoryURL = ClassLoader.getSystemResource("common/Files.class");
				String repositoryPath = repositoryURL.getPath();
				URI jarURI = new URI(repositoryPath.substring(0, repositoryPath.indexOf('!')));
				JarFile jarFile = new JarFile(new File(jarURI));
				Vector<String> v = new Vector<String>();

				Enumeration<JarEntry> entries = jarFile.entries();
				while (entries.hasMoreElements()) {
					String entryName = entries.nextElement().getName();

					if (entryName.startsWith(directory)) {
						String trav = entryName.substring(directory.length());
						if (trav.startsWith("/")) {
							trav = trav.substring(1);
						}
						if (trav.contains("/")) {
							HashVec.addIfAbsent(trav.substring(0, trav.indexOf("/")), v, true);
						}
					}
				}

				jarFile.close();

				return Array.toStringArray(v);
			} catch (Exception e) {
				System.err.println("Error reading files in jar file");
				e.printStackTrace();
				return null;
			}
		} else {
			return new File(directory).list(new FilenameFilter() {
				@Override
				public boolean accept(File file, String filename) {
					return new File(file, filename).isDirectory();
				}
			});
		}
	}

	public static void summarizeAllFilesInDirectory(String dir) {
		PrintWriter writer;
		String[] data;

		data = listAllFilesInTree(dir, false);
		try {
			writer = new PrintWriter(new FileWriter(dir + "summaryOfFiles.xln"));
			writer.println("Full filename\tdirectory\tfilename\troot\textension");
			for (String element : data) {
				writer.println(element	+ "\t" + ext.parseDirectoryOfFile(element) + "\t"
												+ ext.removeDirectoryInfo(element) + "\t" + ext.rootOf(element) + "\t"
												+ element.substring(element.lastIndexOf(".") + 1));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir + "summaryOfFiles.xln");
			e.printStackTrace();
		}
	}

	public static void summarizeDirectoryFromParameters(String filename, Logger log) {
		String[][] params;

		params = parseControlFile(filename, false, "dir", new String[] {"directory"}, log);
		if (params != null) {
			summarizeAllFilesInDirectory(params[0][0]);
		}
	}

	public static String[] listAllFilesInTree(String dir, boolean jar) {
		Vector<String> allFiles;

		allFiles = new Vector<String>();
		traverseTree(dir, "", allFiles, jar);
		return Array.toStringArray(allFiles);
	}

	private static void traverseTree(String root, String dir, Vector<String> allFiles, boolean jar) {
		String[] dirs;

		dirs = listDirectories(root + dir, jar);
		for (String dir2 : dirs) {
			traverseTree(root, dir + dir2 + "/", allFiles, jar);
		}

		HashVec.addAllInArrayToVector(Array.addPrefixSuffixToArray(	list(root + dir, null, jar), dir,
																																null),
																	allFiles);
	}

	public static void copySpecificFiles(String filename, Logger log) {
		String[][] params;
		String sourceDir, targetDir;
		boolean lower;

		lower = false;
		params = parseControlFile(filename, false, "copy",
															new String[] {"sourceDirectory/", "targetDirectory/ -toLowerCase",
																						"file1.txt", "file2.txt"},
															log);
		if (params != null) {
			sourceDir = params[0][0];
			targetDir = params[1][0];
			for (int i = 1; i < params[1].length; i++) {
				if (params[1][i].equals("-toLowerCase")) {
					lower = true;
				} else {
					log.reportError("Error - don't know what to do with argument: " + params[1][i]);
				}
			}
			new File(targetDir).mkdirs();
			for (int i = 2; i < params.length; i++) {
				Files.copyFile(	new File(sourceDir + params[i][0]).toString(),
												new File(targetDir + (lower	? params[i][0].toLowerCase()
																										: params[i][0])).toString());
			}
		}
	}

	public static void more(String filename) {
		BufferedReader reader;

		try {
			reader = new BufferedReader(new FileReader(filename));
			while (reader.ready()) {
				System.out.println(reader.readLine());
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public static void catFilesFromParameters(String filename, Logger log) {
		List<String> params;
		String[] line, files;
		int[] skips;
		String out;
		boolean addFilename = false;

		// get all files in the directory, excluding the crf itself and its corresponding log
		files = Files.list("./", ":" + ext.rootOf(filename), ":.crf", false, false);
		files = Array.addStrToArray("outfile.xln", files, 0);
		files =
					Array.addStrToArray("# include add_filename_as_first_column after the output filename if you want it",
															files, 1);
		for (int i = 3; i < files.length; i++) {
			files[i] += " skip=1";
		}

		params = parseControlFile(filename, "cat", files, log);
		if (params != null) {
			out = params.remove(0);
			if (out.contains("add_filename_as_first_column")) {
				out = out.substring(0, out.indexOf("add_filename_as_first_column")).trim();
				addFilename = true;
			}
			files = new String[params.size()];
			skips = Array.intArray(params.size(), 0);
			for (int i = 0; i < files.length; i++) {
				line = params.get(i).trim().split("[\\s]+");
				files[i] = line[0];
				for (int j = 1; j < line.length; j++) {
					if (line[j].startsWith("skip=")) {
						skips[i] = Integer.parseInt(line[j].split("=")[1]);
					}
				}
			}
			cat(files, out, skips, addFilename, log);
		}
	}

	public static void renameFilesFromParameters(String filename, Logger log) {
		List<String> params;
		String[] line, files;
		String[][] matchingFilenames;
		boolean problem;

		// get all files in the directory, excluding the crf itself and its corresponding log
		files = Files.list("./", ":" + ext.rootOf(filename), ":.crf", false, false);
		params = parseControlFile(filename, "rename",
															Array.toStringArray(Matrix.transpose(new String[][] {files, files}),
																									"\t"),
															log);
		if (params != null) {
			matchingFilenames = new String[params.size()][2];
			problem = false;
			for (int i = 0; i < matchingFilenames.length; i++) {
				line = params.get(i).trim().split("\t", -1);
				if (line.length != 2) {
					log.reportError("Error - skipping this line, invalid number of arguments (needs to be tab delimited): "
													+ Array.toStr(line, "/"));
					problem = true;
				} else {
					matchingFilenames[i] = line;
					if (!Files.exists(line[0])) {
						log.reportError("Error - file '"	+ line[0] + "' does not exist, cannot be renamed to "
														+ line[1]);
						problem = true;
					}
					if (Files.exists(line[1])) {
						log.reportError("Error - file '"	+ line[1]
														+ "' already exists, and will not be overwritten by the contents of "
														+ line[0]);
						problem = true;
					}
				}
			}
			if (!problem) {
				for (String[] matchingFilename : matchingFilenames) {
					log.report("Renaming " + matchingFilename[0] + " to " + matchingFilename[1]);
					new File(matchingFilename[0]).renameTo(new File(matchingFilename[1]));
				}
			}

		}
	}

	/**
	 * Sort of like the linux paste function, files must be the same length, and in the same order
	 */
	public static String[][] paste(	String[] orginalFiles, String finalFile, int[] columns,
																	int keyColumn, String[] newColumnNameTags, String[] skipKeys,
																	Logger log) {
		if (newColumnNameTags.length != orginalFiles.length) {
			log.reportError("Name tags must be the same length as files");
			return null;
		}
		try {
			PrintWriter writer = Files.getAppropriateWriter(finalFile);
			BufferedReader[] readers = new BufferedReader[orginalFiles.length];
			String[][] newColumns = new String[orginalFiles.length][columns.length];
			for (int i = 0; i < orginalFiles.length; i++) {
				readers[i] = Files.getAppropriateReader(orginalFiles[i]);
				readers[i].readLine();
				String[] curHeader = Files.getHeaderOfFile(orginalFiles[i], log);
				for (int j = 0; j < columns.length; j++) {
					String newColumn = curHeader[columns[j]] + "_" + newColumnNameTags[i];
					writer.print((i == 0 && j == 0 ? "" : "\t") + newColumn);
					newColumns[i][j] = newColumn;
				}
			}
			writer.println();

			while (readers[0].ready()) {
				String[] line = readers[0].readLine().trim().split("\t");
				String key = line[keyColumn];
				if (skipKeys == null || ext.indexOfStr(key, skipKeys) < 0) {
					writer.print(Array.toStr(Array.subArray(line, columns)));
					for (int i = 1; i < readers.length; i++) {
						String[] oLine = readers[i].readLine().trim().split("\t");
						if (!key.equals(oLine[keyColumn])) {
							String error = "Files must be in the same order as defined by the key column";
							log.reportError(error);
							throw new IllegalStateException(error);
						} else {
							writer.print("\t" + Array.toStr(Array.subArray(oLine, columns)));
						}
					}
					writer.println();
				} else {
					for (int i = 1; i < readers.length; i++) {
						readers[i].readLine();
					}
				}
			}
			for (int i = 0; i < readers.length; i++) {
				if (readers[i].ready()) {
					String error = orginalFiles[i] + " still has more lines";

					log.reportError(error);
					throw new IllegalStateException(error);

				} else {
					readers[i].close();
				}
			}
			writer.close();
			return newColumns;
		} catch (Exception e) {
			log.reportError("Error writing to " + finalFile);
			log.reportException(e);
		}
		return null;
	}

	public static void renameFilesUsingSubstitutionsFromParameters(String filename, Logger log) {
		List<String> params;
		String[] line, dirs, filenames, newFilenames;
		String[][] substitutions;
		boolean echo, extract;
		int count;
		long summedSizes = 0;
		long[] filesizes;
		int[] order, indices;

		echo = extract = false;
		params = parseControlFile(filename, "subs",
															new String[] {"# echo reports what it's going to change without actually doing it",
																						"echo=true",
																						"# Extract the largest file from each subdirectory and script the deletion of the subdirectory",
																						"extract=false",
																						"# Examples: two tabbed delimited columns represents a replacement, anything with a single column will just be deleted",
																						"studyname.\tStudyName.", ".GREEN.IDAT\t_Grn.idat",
																						".RED.IDAT\t_Red.idat", "justDeleteThisPhrase",
																						"replaceThis\twithThis"},
															log);
		if (params != null) {
			substitutions = new String[params.size()][];
			count = 0;
			while (count < params.size()) {
				line = params.get(count).trim().split("\t", -1);
				if (line[0].startsWith("echo=")) {
					echo = ext.parseBooleanArg(line[0]);
					params.remove(count);
				} else if (line[0].startsWith("extract=")) {
					extract = ext.parseBooleanArg(line[0]);
					params.remove(count);
				} else {
					if (line.length == 1) {
						line = Array.addStrToArray("", line);
					}
					if (line.length > 2) {
						System.out.println("Error more than two tokens in: " + Array.toStr(line, " / "));
					}
					substitutions[count] = line;
					count++;
				}
			}
			if (substitutions.length > count) {
				boolean[] keeps = Array.booleanArray(substitutions.length, true);
				for (int i = count; i < keeps.length; i++) {
					keeps[i] = false;
				}
				substitutions = Matrix.subset(substitutions, keeps);
			}

			if (extract) {
				// get all sub-directories
				dirs = Files.listDirectories("./", false);
				for (String dir : dirs) {
					filenames = Files.list(dir, null, null, false, false);
					filesizes = new long[filenames.length];
					for (int j = 0; j < filenames.length; j++) {
						filesizes[j] = new File(dir + "/" + filenames[j]).length();
					}
					order = Sort.getReverseIndices(filesizes);
					if (echo) {
						System.out.println("Searching " + dir);
					}

					for (int j = 1; j < order.length; j++) {
						summedSizes += filesizes[order[j]];
					}

					if (filesizes[order[0]] < summedSizes) {
						System.out.println("  a dominant file could not be found; extraction skipped for this directory");
					} else {
						for (int j = 0; j < Math.min(order.length, 1); j++) {
							if (echo) {
								System.out.println("  Moving " + filenames[order[j]] + " one level up");
							} else {
								new File(dir + "/" + filenames[order[j]]).renameTo(new File(filenames[order[j]]));
							}
						}
						for (int j = 1; j < order.length; j++) {
							if (echo) {
								System.out.println("  deleting " + filenames[order[j]]);
							} else {
								new File(dir + "/" + filenames[order[j]]).delete();
							}
						}
						if (echo) {
							System.out.println(" deleting " + dir);
						} else {
							new File(dir).delete();
						}
					}
				}
			}

			// get all files in the directory, excluding the crf itself and its corresponding log
			filenames = Files.list("./", ":" + ext.rootOf(filename), ":.crf", false, false);
			newFilenames = new String[filenames.length];

			for (int i = 0; i < filenames.length; i++) {
				newFilenames[i] = ext.replaceAllWith(filenames[i], substitutions);
				indices = ext.indicesWithinString(".", newFilenames[i]);
				for (int j = 0; j < indices.length - 1; j++) {
					newFilenames[i] = newFilenames[i].substring(0, indices[j] + 1)
															+ newFilenames[i]	.substring(indices[j] + 1, indices[j] + 2)
																							.toUpperCase()
														+ newFilenames[i].substring(indices[j] + 2);
				}
				if (!filenames[i].equals(newFilenames[i])) {
					if (echo) {
						System.out.println("Renaming '" + filenames[i] + "' to '" + newFilenames[i] + "'");
					} else {
						new File(filenames[i]).renameTo(new File(newFilenames[i]));
					}
				}
			}
		}
	}

	public static void cat(String[] originalFiles, String finalFile, int[] skips, Logger log) {
		cat(originalFiles, finalFile, skips, false, log);
	}

	// can pass skips as null if there are no skips to be made
	// can pass skips as an empty array (new int[0]) if the first file should pass a header but the
	// rest should be skipped
	public static void cat(	String[] originalFiles, String finalFile, int[] skips, boolean addFilename,
													Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String trav;
		boolean problem;
		String delimiter;

		delimiter = finalFile.endsWith(".csv") || finalFile.endsWith(".csv.gz") ? "," : "\t";

		if (skips != null) {
			if (skips.length == 0) {
				skips = Array.addIntToArray(0, Array.intArray(originalFiles.length - 1, 1), 0);
			}
			if (skips.length != originalFiles.length) {
				log.reportError("Error - mismatched length of arrays for the files and number of lines to skip for file "
												+ finalFile + "; aborting...");
				return;
			}
		}

		problem = false;
		for (int i = 0; i < originalFiles.length; i++) {
			if (originalFiles[i] == null) {
				log.reportError("Error - can't cat if the filename is null");
				problem = true;
			} else if (!new File(originalFiles[i]).exists()) {
				log.reportError("Error - missing file '" + originalFiles[i] + "'");
				problem = true;
			}
		}
		if (problem) {
			return;
		}

		try {
			writer = Files.getAppropriateWriter(finalFile);
			for (int i = 0; i < originalFiles.length; i++) {
				try {

					reader = new BufferedReader(new FileReader(originalFiles[i]));
					for (int j = 0; skips != null && j < skips[i]; j++) {
						reader.readLine();
					}
					while (reader.ready()) {
						trav = reader.readLine();
						writer.println((addFilename ? originalFiles[i] + delimiter : "") + trav);
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					log.reportError("Error: file \""	+ originalFiles[i]
													+ "\" not found in current directory");
				} catch (IOException ioe) {
					log.reportError("Error reading file \"" + originalFiles[i] + "\"");
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + finalFile);
			e.printStackTrace();
		}
	}

	// this is much faster than any c=in.read() method would be (~3.7 times faster for a large .mldose
	// file)
	// however it is almost twice as slow as the LineNumberReader below
	public static int countLinesSlow(String filename, boolean dontCountFirstLine) {
		BufferedReader reader;
		int count;

		try {
			reader = getAppropriateReader(filename);
			for (count = 0; reader.ready() && reader.readLine() != null; count++) {
				;
			}
			reader.close();
			return count - (dontCountFirstLine ? 1 : 0);
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
		}

		return -1;
	}

	/**
	 * Put a tester at {@link NumLineTest}. Let JL know if this is messed up
	 *
	 */
	public static int countLines(String filename, int numberOfLinesNotToCount) {
		try {
			LineNumberReader lnr = new LineNumberReader(getAppropriateInputStreamReader(filename));
			long num = -1;
			while (num != 0) {
				num = lnr.skip(8192);
			}
			int lines = lnr.getLineNumber();// I think we do not need the +1, seems similar to a
																			// array.length
			lnr.close();
			return lines - numberOfLinesNotToCount;
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
		}
		return -1;
	}

	public static boolean countIfMoreThanNLines(String filename, int n) {
		BufferedReader reader;
		int count;

		try {
			reader = getAppropriateReader(filename);
			count = 0;
			while (reader.ready()) {
				reader.readLine();
				count++;
				if (count > n) {
					return true;
				}
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
		}

		return false;
	}

	public static void write(String str, String filename) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println(str);
			writer.flush();
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}
	}

	public static void writeArray(String[] entries, String filename) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (String element : entries) {
				writer.println(element);
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}
	}

	public static void writeIterable(Iterable<String> entries, String filename) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (String element : entries) {
				writer.println(element);
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}
	}

	public static void writeMatrix(String[][] matrix, String filename, String delimiterToUse) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (String[] element : matrix) {
				writer.println(Array.toStr(element, delimiterToUse));
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}
	}

	public static boolean isFileReady(String filename, long timeOfFirstCheck,
																		long timeIntervalBetweenChecking) {
		File file;
		long fileSize;

		file = new File(filename);
		try {
			Thread.sleep(timeOfFirstCheck);
			if (!file.exists()) {
				System.out.println("Cannot find the file:" + timeOfFirstCheck + "_InvVar");
				return false;
			}

			fileSize = -1;
			while (fileSize != file.length()) {
				fileSize = file.length();
				Thread.sleep(timeIntervalBetweenChecking);
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		return true;
	}

	public static void splitFile(String[] line) {
		String filename = line[0];
		int numSplits = 2;
		int numLinesCopiedToAll = 1;
		int blockSize = 1;
		String rootForNewFiles = "list";
		String extForNewFiles = ".dat";
		boolean allowUnevenBlocks = false;

		for (int i = 1; i < line.length; i++) {
			if (line[i].startsWith("numFiles=")) {
				numSplits = ext.parseIntArg(line[i]);
			} else if (line[i].startsWith("sizeOfHeader=")) {
				numLinesCopiedToAll = ext.parseIntArg(line[i]);
			} else if (line[i].startsWith("blockSize=")) {
				blockSize = ext.parseIntArg(line[i]);
			} else if (line[i].startsWith("root=")) {
				rootForNewFiles = ext.parseStringArg(line[i], "");
			} else if (line[i].startsWith("ext=")) {
				extForNewFiles = ext.parseStringArg(line[i], "");
			} else if (line[i].startsWith("allowUnevenBlocks=")) {
				allowUnevenBlocks = ext.parseBooleanArg(line[i]);
			} else {
				System.err.println("Error - unknown argument");
				return;
			}
		}

		splitFile(filename, numSplits, numLinesCopiedToAll, blockSize, rootForNewFiles, extForNewFiles,
							allowUnevenBlocks);
	}

	public static void splitFileFromParameters(String filename, Logger log) {
		String[][] params;

		params = parseControlFile(filename, false, "split",
															new String[] {"sourcefile.txt numFiles=6 sizeOfHeader=1 blockSize=1 root=list ext=.dat"},
															log);
		if (params != null) {
			splitFile(params[0]);
		}
	}

	public static void splitFile(	String filename, int numSplits, int numLinesCopiedToAll,
																int blockSize, String rootForNewFiles, String extForNewFiles,
																boolean allowUnevenBlocks) {
		BufferedReader reader;
		PrintWriter writer;
		int count;
		String header;

		for (int i = 1; i <= numSplits; i++) {
			new File(rootForNewFiles + i + extForNewFiles).delete();
		}

		try {
			reader = new BufferedReader(new FileReader(filename));
			header = "";
			for (int i = 0; i < numLinesCopiedToAll; i++) {
				header += reader.readLine() + "\n";
			}
			count = 0;
			writer = null;
			while (reader.ready()) {
				if (writer != null) {
					writer.close();
				}
				writer = new PrintWriter(new FileWriter(rootForNewFiles	+ (count % numSplits + 1)
																								+ extForNewFiles, true));
				if (new File(rootForNewFiles + (count % numSplits + 1) + extForNewFiles).length() == 0) {
					writer.print(header);
				}
				for (int i = 0; i < blockSize; i++) {
					if (reader.ready()) {
						writer.println(reader.readLine());
					} else {
						if (!allowUnevenBlocks) {
							System.out.println("Error - invalid block size or trailing whitespace; last block size was only "
																	+ i);
							writer.close();
							try {
								writer = new PrintWriter(new FileWriter(filename + "_BLOCK_ERROR.log"));
								writer.println("Error - invalid block size or trailing whitespace; last block size was only "
																+ i);
								writer.close();
							} catch (Exception e) {
								System.err.println("Error writing to " + filename + ".log");
								e.printStackTrace();
							}
							i = blockSize;
						}
					}
				}
				count++;
			}
			reader.close();
			if (writer != null) {
				writer.close();
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public static void splitFilesByLine(String dir, String suffix, int start, int stop) {
		BufferedReader[] readers;
		PrintWriter writer;
		String[] files;


		files = list(dir, suffix, false);
		readers = new BufferedReader[files.length];
		try {
			for (int i = 0; i < files.length; i++) {
				readers[i] = new BufferedReader(new FileReader(dir + files[i]));
			}

			for (int i = 0; i < stop; i++) {
				writer = new PrintWriter(new FileWriter(dir + (i + 1) + ".xln"));
				for (int j = 0; j < files.length; j++) {
					writer.println(readers[j].readLine());
				}
				writer.close();
			}

			for (int i = 0; i < files.length; i++) {
				readers[i].close();
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + "\"");
			System.exit(2);
		}



	}

	// className keeps track of which method is altering the file
	public static String getBakFilename(String base, String className) {
		return getBakFilename(base, className, false);
	}

	public static String getBakFilename(String filename, String className, boolean renameOldToNew) {
		String bakFilename;
		int count = 0;

		if (className.indexOf(".") > 0) {
			className = className.substring(className.indexOf(".") + 1);
		}
		do {
			count++;
			bakFilename = filename + "-" + className + "(" + count + ").bak";
		} while (new File(bakFilename).exists());

		if (renameOldToNew) {
			if (new File(filename).exists()) {
				new File(filename).renameTo(new File(bakFilename));
			} else {
				System.err.println("Error - " + filename + " was not found in current directory");
			}
		}

		return bakFilename;
	}

	public static int findNextRep(String pattern) {
		return findNextRep(new String[] {pattern}, -1, 0);
	}

	public static int findNextRep(String[] patterns, int numDigits, int lastKnownRep) {
		boolean done;
		int count;

		count = Math.max(lastKnownRep - 1, 0);
		do {
			count++;
			done = true;
			for (String pattern : patterns) {
				if (new File(ext.insertNumbers(pattern, count, numDigits)).exists()) {
					done = false;
				}
			}
		} while (!done);

		return count;
	}

	public static int findNextRepSafely(String[] patterns, int numDigits, int lastKnownRep,
																			int patienceInMilliseconds) {
		int trav, rand;
		boolean done;
		String[] files;
		int[] rands;

		trav = lastKnownRep;
		rand = (int) (Math.random() * 100000000);
		do {
			trav++;
			done = true;
			trav = findNextRep(Array.addStrToArray("#.taken", patterns), numDigits, trav);
			Files.writeArray(new String[0], trav + "." + rand + ".reserved");
			if (new File(trav + ".taken").exists()) {
				done = false;
			} else {
				for (int repeat = 0; done && repeat < 2; repeat++) {
					try {
						Thread.sleep(patienceInMilliseconds / 2);
					} catch (InterruptedException ie) {
					}

					if (new File(trav + ".taken").exists()) {
						done = false;
					} else {
						files = list("./", trav + ".", ".reserved", false, false);
						if (files.length > 1) {
							rands = new int[files.length];
							for (int i = 0; i < rands.length; i++) {
								rands[i] = Integer.parseInt(files[i].substring(files[i].indexOf(".")	+ 1,
																																files[i].lastIndexOf(".")));
							}
							Arrays.sort(rands);
							if (rand != rands[0]) {
								new File(trav + "." + rand + ".reserved").delete();
								done = false;
							}
						}
					}
				}
			}
		} while (!done);

		Files.writeArray(new String[0], trav + ".taken");
		new File(trav + "." + rand + ".reserved").delete();

		return trav;
	}

	// public static void deleteLater(String filename, int timeToWait) {
	// Thread thread;
	//
	// thread = new Thread(new DeleteLater(filename, timeToWait));
	// thread.setDaemon(true);
	// thread.start();
	// }
	//
	public static long getMostRecentUpdate(String dir) {
		String[] files;
		long max;

		if (!new File(dir).exists()) {
			return -1;
		}

		dir = ext.verifyDirFormat(dir);

		max = -1;
		files = list(dir, null, false);
		for (String file : files) {
			max = Math.max(max, new File(dir + file).lastModified());
		}

		return max;
	}



	public static boolean headerOfFileContainsAll(String filename, String[] toSearch, Logger log) {
		return headerOfFileContainsAll(filename, toSearch, true, log);
	}

	/**
	 * @param filename search the header of this file
	 * @param toSearch strings to search the header for
	 * @param verbose report missing strings
	 * @param log
	 * @return true if all strings are present
	 */
	public static boolean headerOfFileContainsAll(String filename, String[] toSearch, boolean verbose,
																								Logger log) {
		String[] header = getHeaderOfFile(filename, log);
		boolean has = true;
		if (toSearch.length > header.length) {
			has = false;
		} else if (Array.countIf(	ext.indexFactors(toSearch, header, true, log, verbose, false),
															-1) > 0) {
			if (verbose) {
				log.reportError("Searched header " + Array.toStr(header));
			}
			has = false;
		}
		return has;
	}

	public static String[] getHeaderOfFile(String filename, Logger log) {
		return getHeaderOfFile(filename, null, log);
	}

	public static String[] getHeaderOfFile(String filename, String delimiter, Logger log) {
		return getHeaderOfFile(filename, delimiter, null, log);
	}



	// if you want to force the removal of quotes in a comma-delimited file, then set delimiter to
	// ",!"
	/**
	 * @param filename
	 * @param delimiter
	 * @param patternToSkip lines starting with this pattern will not be included, set to null to skip
	 *        this filter
	 * @param log
	 * @return
	 */
	public static String[] getHeaderOfFile(	String filename, String delimiter,
																					String[] patternsThatIfPresentAtTheStartOfTheLineWillSkipLine,
																					Logger log) {
		String[] lines;

		lines = getFirstNLinesOfFile(	filename, 1, patternsThatIfPresentAtTheStartOfTheLineWillSkipLine,
																	log);
		if (lines == null) {
			return null;
		}

		if (delimiter == null) {
			delimiter = ext.determineDelimiter(lines[0]);
		}

		if (delimiter.startsWith(",")) {
			return ext.splitCommasIntelligently(lines[0], delimiter.endsWith("!"), log);
		}

		return lines[0].trim().split(delimiter);
	}

	public static String[] getFirstNLinesOfFile(String filename, int nLines, Logger log) {
		return getFirstNLinesOfFile(filename, nLines, null, log);
	}

	/**
	 * @param filename
	 * @param nLines
	 * @param patternsToSkip lines starting with this pattern will not be included, set to null to
	 *        skip this filter
	 * @param log
	 * @return
	 */
	public static String[] getFirstNLinesOfFile(String filename, int nLines, String[] patternsToSkip,
																							Logger log) {
		BufferedReader reader;
		String[] lines;
		Vector<String> v;
		int count;

		lines = null;
		try {
			reader = getAppropriateReader(filename);
			if (reader == null && log != null) {
				log.reportError("Error: file \"" + filename + "\" not found in current directory");
			}
			v = new Vector<String>();
			count = 0;
			while (reader.ready() && count < nLines) {
				String line = reader.readLine();
				boolean add = true;
				if (patternsToSkip != null) {
					for (String element : patternsToSkip) {
						if (element != null && line.startsWith(element)) {
							add = false;
							break;
						}
					}
				}
				if (add) {
					v.add(line);
					count++;
				}

			}
			reader.close();

			return Array.toStringArray(v);

		} catch (FileNotFoundException fnfe) {
			if (log != null) {
				log.reportError("Error: file \"" + filename + "\" not found in current directory");
			}
		} catch (IOException ioe) {
			if (log != null) {
				log.reportError("Error reading file \"" + filename + "\"");
			}
		}

		return lines;
	}

	/**
	 * @param filename
	 * @param delimiter
	 * @param containing
	 * @param log
	 * @return the first line containing the containing array, null if all are not found
	 */
	public static String[] getLineContaining(	String filename, String delimiter, String[] containing,
																						Logger log) {

		try {
			BufferedReader reader = Files.getAppropriateReader(filename);

			String line;
			while ((line = reader.readLine()) != null) {
				String[] tmp = line.trim().split(delimiter);
				if (Array.countIf(ext.indexFactors(containing, tmp, false, log, false, false), -1) == 0) {
					reader.close();
					return tmp;
				}
			}
			reader.close();
			log.reportError("Could not find the header containing "	+ Array.toStr(containing)
											+ " in file " + filename);
			return null;
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			return null;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			return null;
		}
	}

	public static String[][] parseControlFile(String filename, boolean tab, String command,
																						String[] sampleCode, Logger log) {
		List<String> list;

		list = parseControlFile(filename, command, sampleCode, log);
		if (list == null) {
			return null;
		} else {
			return Array.splitStrings(Array.toStringArray(list), tab);
		}
	}

	public static List<String> parseControlFile(String filename, String command,
																								String[] sampleCode, Logger log) {
		return parseControlFile(filename, command, Arrays.asList(sampleCode), log);
	}

	public static List<String> parseControlFile(String filename, String command,
																								List<String> sampleCode, Logger log) {
		List<String> outputCode = Lists.newArrayList(sampleCode);
		outputCode.add(0, command);

		if (new File(filename).length() == 0) {
			Files.writeIterable(outputCode, filename);
			return null;
		}
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(filename));
			String[] line = reader.readLine().trim().split("[\\s]+");
			if (!line[0].equalsIgnoreCase(command)) {
				log.reportError("Error - file must start with the line '" + command + "'");
				reader.close();
				return null;
			}
			if (!reader.ready()) {
				reader.close();
				Files.writeIterable(outputCode, filename);
				return null;
			} else {
				List<String> ret = Lists.newArrayList();
				while (reader.ready()) {
					String temp = reader.readLine();
					if (!temp.startsWith("#") && !temp.trim().split("[\\s]+")[0].equals("")) {
						ret.add(temp);
					}
				}
				reader.close();
				return ret;
			}
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			log.reportException(fnfe);
			return null;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			log.reportException(ioe);
			return null;
		} finally {
			Closeables.closeQuietly(reader);
		}
	}

	public static void makeQsub(String filename, boolean multiple, int start, int stop,
															boolean separate, String[] patterns,
															boolean changeToCurrentWorkingDirectoryFirst) {
		String[] lines, qsubs;
		Vector<String> v;
		int numThreads;

		lines = HashVec.loadFileToStringArray(filename, false, false, null, false);

		if (multiple) { // could be more elegant and incorporated with those below if desired
			System.out.println("Creating a ScriptExecutor");
			numThreads = Math.min(24, Files.countLines(filename, 0));
			qsub("scriptExecutorFor_"	+ ext.removeDirectoryInfo(filename),
						"cd "																	+ ext.parseDirectoryOfFile(filename) + "\n"
																									+ Files.getRunString()
																									+ " one.ScriptExecutor file="
																									+ ext.removeDirectoryInfo(filename) + " threads="
																									+ numThreads,
						63000, 12, numThreads);
		} else if (separate) {
			v = new Vector<String>();
			for (int i = 0; i < lines.length; i++) {
				qsubs = qsub(	"", ext.rootOf(filename) + (i + 1) + ".#", start, stop,
											(changeToCurrentWorkingDirectoryFirst ? "cd " + ext.pwd() + "\n" : "")
																																							+ lines[i],
											patterns, 63488, 24, null);
				v.add(qsubs[0]);
			}
			writeArray(Array.toStringArray(v), "master." + ext.rootOf(filename));
			Files.chmod("master." + ext.rootOf(filename));
		} else {
			if (changeToCurrentWorkingDirectoryFirst) {
				lines = Array.addStrToArray("cd " + ext.pwd(), lines);
			}
			if (start == stop) {
				qsub(ext.rootOf(filename) + ".qsub", Array.toStr(lines, "\n"), 63488, 24, 1);
			} else {
				qsubs = qsub(	"", ext.rootOf(filename) + "#", start, stop, Array.toStr(lines, "\n"),
											patterns, 63488, 24, null);
				if (qsubs.length > 1) {
					writeArray(qsubs, "master." + ext.rootOf(filename));
					Files.chmod("master." + ext.rootOf(filename));
				}
			}

		}
	}

	/**
	 * Function to get path of serialized file from a filepath
	 *
	 * @param filepath the fileopath
	 * @return the filepath of serialized file obtained by adding ".ser" at end
	 */
	public static String getSerializedFilepath(String filepath) {
		return (filepath + SERIALIZED_FILE_EXTENSION);
	}

	// public static String determineDelimiterFromFileExt(String filename, Logger log) {
	// if (filename.endsWith(".csv") || filename.endsWith(".csv.gz")) {
	// return ",";
	// }
	//
	// if (filename.endsWith(".xln") || filename.endsWith(".xln.gz")) {
	// return "\t";
	// }
	//
	// if (filename.endsWith(".dat") || filename.endsWith(".dat.gz")) {
	// return "\t";
	// }
	//
	// if (filename.endsWith(".txt") || filename.endsWith(".txt.gz")) {
	// return "\t";
	// }
	//
	// return "\t";
	// }
	//
	public static String determineDelimiter(String filename, Logger log) {
		if (filename.endsWith(".csv") || filename.endsWith(".csv.gz")) {
			return ",";
		}

		if (filename.endsWith(".xln") || filename.endsWith(".xln.gz")) {
			return "\t";
		}

		return determineDelimiterFromFirstLine(filename, log);
	}

	public static String suggestDelimiter(String filename, Logger log) {
		if (filename.endsWith(".csv") || filename.endsWith(".csv.gz")) {
			return ",";
		}

		return "\t";
	}

	private static String determineDelimiterFromFirstLine(String filename, Logger log) {
		String[] lines = Files.getFirstNLinesOfFile(filename, 1, null, log);
		if (lines == null || lines.length == 0) {
			return "";
		}
		return ext.determineDelimiter(lines[0]);
	}

	public static void moveFilesMoved(String filesMoved, String dir) {
		String[] files;

		files = HashVec.loadFileToStringArray(filesMoved, false, new int[] {0}, false);
		new File(dir + "moved/").mkdirs();
		for (String file : files) {
			new File(dir + file).renameTo(new File(dir + "moved/" + file));
		}
	}

	public static void makeVLookupReadyFile(String in, String out, int[] keyIndices,
																					int[] valueIndices) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String key, value;

		try {
			reader = new BufferedReader(new FileReader(in));
			writer = new PrintWriter(new FileWriter(out));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				key = "";
				for (int i = 0; i < keyIndices.length; i++) {
					key += (i == 0 ? "" : "^") + line[keyIndices[i]];
				}
				value = "";
				for (int i = 0; i < valueIndices.length; i++) {
					value += (i == 0 ? "" : "^") + line[valueIndices[i]];
				}
				writer.println(key + "\t" + value);
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + in + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + in + "\"");
			System.exit(2);
		}
	}

	/**
	 * Function to remove exntention from a the filenaem of a give file path
	 *
	 * @param filePath the path of the file
	 * @return the name of the file without exntention
	 */
	public static String removeExtention(String filePath) {
		File f = new File(filePath);
		// if it's a directory, don't remove the extention
		if (f.isDirectory()) {
			return f.getName();
		}
		String name = f.getName();
		String path = f.getParent();
		// if it is a hidden file
		if (name.startsWith(".")) {
			// if there is no extn, do not rmove one...
			if (name.lastIndexOf('.') == name.indexOf('.')) {
				return name;
			}
		}
		// if there is no extention, don't do anything
		if (!name.contains(".")) {
			return name;
		}
		// Otherwise, remove the last 'extension type thing'
		if (path != null) {
			return path + File.separator + name.substring(0, name.lastIndexOf('.'));
		} else {
			return name.substring(0, name.lastIndexOf('.'));
		}
	}

	// doesn't work, better to just search to see if the file you want exists, or check the home
	// directory for clarification
	// public static String getHostname() {
	// System.out.println("Tring to get hostname");
	// CmdLine.run("/bin/hostname > hostname.out", "./");
	// System.out.println("...done");
	// if (Files.exists("hostname.out")) {
	// System.out.println("Success!");
	// return Files.getFirstNLinesOfFile("hostname.out", 1, new Logger())[0];
	// }
	// System.out.println("Failed");
	//
	// return "failed to get hostname";
	// }

	public static String getNextAvailableFilename(String pattern) {
		String filename;
		int count;

		count = 0;
		do {
			count++;
			filename = ext.replaceAllWith(pattern, "#", count + "");
		} while (exists(filename));

		return filename;
	}

	public static void closeAll(BufferedReader[] readers) {
		for (BufferedReader reader : readers) {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	public static void closeAll(PrintWriter[] writers) {
		for (PrintWriter writer : writers) {
			if (writer != null) {
				writer.close();
			}
		}
	}

	/**
	 * Warning - does not check existence,null, etc.. just appends and moves on
	 *
	 * @param fileNamesNoDirectory these are fileNames with directory info removed (i.e from
	 *        Files.list)
	 * @param directoryToAppend append this directory info to each file in fileNamesNoDirectory
	 * @return the appended fileNames
	 */
	public static String[] toFullPaths(String[] fileNamesNoDirectory, String directoryToAppend) {
		String[] fullPaths = new String[fileNamesNoDirectory.length];
		for (int i = 0; i < fileNamesNoDirectory.length; i++) {
			fullPaths[i] = directoryToAppend + fileNamesNoDirectory[i];
		}
		return fullPaths;
	}

	/**
	 * trys to determine if the path is relative or absolute
	 *
	 * @param path
	 */
	public static boolean isRelativePath(String path) {
		return (!path.startsWith("/") && path.indexOf(":") < 0);
	}

	/**
	 * Attempts to grab all file paths of a type from an ftp site
	 *
	 * @param ftpdirAddress a remote directory of an ftp site
	 * @param type the type of file to collect
	 * @param log
	 * @return
	 */
	public static String[] parseRemoteFTPFiles(String ftpdirAddress, String type, Logger log) {
		ArrayList<String> remoteVcfs = new ArrayList<String>();
		URL url;
		if (!ftpdirAddress.startsWith("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/")) {
			log.reportTimeWarning("Did not detect that "	+ ftpdirAddress
														+ " was an ftp address starting with ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/");
			log.reportTimeWarning("\t this parsing method is therefore un-tested");
		}
		try {
			url = new URL(ftpdirAddress);
			BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
			String inputLine;

			while ((inputLine = in.readLine()) != null) {
				String[] line = inputLine.trim().split("[\\s]+");
				int possibleVcfIndex = ext.indexOfEndsWith(type, line, false);
				if (possibleVcfIndex > 0) {
					remoteVcfs.add(ftpdirAddress + line[possibleVcfIndex]);
				}
			}
			in.close();

		} catch (MalformedURLException e) {
			log.reportError("malformed URL " + ftpdirAddress);
			e.printStackTrace();
		} catch (IOException e) {
			log.reportIOException(ftpdirAddress);
			e.printStackTrace();
		}
		return remoteVcfs.toArray(new String[remoteVcfs.size()]);
	}

	public static boolean programExists(String programName) {
		String cmd = "";
		if (Files.isWindows()) {
			cmd = "where " + programName;
		} else {
			cmd = "which " + programName;
		}
		boolean exists = false;
		try {
			if (Runtime.getRuntime().exec(cmd).waitFor() == 0) {
				exists = true;
			}
		} catch (InterruptedException e1) {
		} catch (IOException e1) {
		}
		return exists;
	}

	public static void appendStringToFile(String filename, String str) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(filename, true));
			writer.println(str);
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			System.err.println(e);
		}
	}

	public static void replaceAll(String filename, String outfile, String replacements, Logger log) {
		BufferedReader reader;
		String[] line;
		String temp;
		Vector<String[]> list = new Vector<String[]>();
		boolean problem = false;

		try {
			reader = Files.getAppropriateReader(replacements);
			while (reader.ready()) {
				temp = reader.readLine();
				if (!temp.startsWith("#") && !temp.trim().equals("")) {
					line = temp.split("\t", -1);
					if (line.length == 1) {
						log.report("Will delete all instances of \"" + line[0] + "\"");
					} else if (line[1].indexOf(line[0]) != -1) {
						log.reportError("Error - the replacement cannot be an extension of itself otherwise it will result in an infinite loop");
						log.reportError("\"" + line + "\" to \"" + line[1] + "\"");
						problem = true;
					}
					if (line.length > 2) {
						log.reportError("Warning - the following line has more than two TABs in it, only the first and second will be used");
						log.reportError(temp);
					}

					if (!line[0].equals("")) {
						list.add(line.length > 1	? new String[] {line[0], line[1]}
																			: new String[] {line[0], ""});
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + replacements + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + replacements + "\"");
			return;
		}

		if (!problem) {
			replaceAll(filename, outfile, Matrix.toStringArrays(list), log);
		}
	}

	public static void replaceAll(String filename, String outfile, String[][] replacements,
																Logger log) {
		BufferedReader reader;
		PrintWriter writer;

		try {
			reader = Files.getAppropriateReader(filename);
			writer = Files.getAppropriateWriter(outfile);
			while (reader.ready()) {
				writer.println(ext.replaceAllWith(reader.readLine(), replacements));
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			log.reportException(ioe);
			return;
		}
	}

	public static void replaceAllFromParameters(String filename, Logger log) {
		List<String> params;

		params = Files.parseControlFile(filename, "replaceAll",
																		new String[] {"file=input.txt.gz", "out=outputFile.txt.gz",
																									"# the swap/replacement file is two tab delimited columns of what to search for (first col) and what to replace it with (second col)",
																									"swap=replacements.txt"},
																		log);

		if (params != null) {
			params.add("log=" + log.getFilename());
			main(Array.toStringArray(params));
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		int start = 1;
		int stop = 1;
		boolean separate = false;
		String[] patterns = null;
		int numDigits = -1;
		boolean findNextRep = false;
		int lastKnownRep = 0;
		int patienceInMilliseconds = 1000;
		String transpose = null;
		boolean commaDelimitedIn = false;
		boolean commaDelimitedOut = false;
		String dir = null;
		String outfile = null;
		String keyfile = null;
		boolean multiple = false;
		boolean cwd = false;
		boolean wc = false;
		String replacements = null;
		Logger log;
		String logfile = null;
		boolean filter = false;
		boolean ord = true;
		int keyIndex = -1;

		String usage = "\n"	+ "common.Files requires 0-1 arguments\n"
										+ "   (1) filename to convert to a .qsub (i.e. file=batchFile (not the default))\n"
										+ "   (2) make it multithreaded (i.e. -multiple (not the default))\n"
										+ "   (3) (optional) chr/rep to start with (use # or ## within file to designate where to use) (i.e. start="
										+ start + " (default))\n" + "   (4) (optional) chr/rep to end with (i.e. stop="
										+ stop + " (default))\n"
										+ "   (5) separate each line into a separate file (i.e. separate=" + separate
										+ " (default))\n"
										+ "   (6) (optional) don't stop until plug is pulled, counting based on patterns (i.e. patterns=perm#.log;perm#.assoc;perm#.assoc.mperm (not the default))\n"
										+ "   (7) (optional) change to current working directory at the top of [each] script (i.e. -cwd (not the default))\n"
										+ "  OR\n" + "   (1) find next rep safely (i.e. -nextRep (not the default))\n"
										+ "   (2) (required) patterns to match when incrementing rep (i.e. patterns=perm#.log;perm#.assoc;perm#.assoc.mperm (not the default))\n"
										+ "   (3) passing the last known rep, speeds things up (i.e. lastRep="
										+ lastKnownRep + " (default))\n"
										+ "   (4) time to wait in milliseconds, in order to ensure no ties (i.e. wait="
										+ patienceInMilliseconds + " (default))\n" + "  OR\n"
										+ "   (1) transpose file (i.e. transpose=file.txt (not the default))\n"
										+ "   (2) input file is comma-delimited (i.e. commaDelimitedIn="
										+ commaDelimitedIn + " (default))\n"
										+ "   (3) name of output file (i.e. out=[input]-transposed.xln (default))\n"
										+ "   (3) output file is comma-delimited (i.e. commaDelimitedOut="
										+ commaDelimitedOut + " (default))\n" + "  OR\n"
										+ "   (1) move files already successfully (i.e. currently hard coded (not the default))\n"
										+ "  OR\n"
										+ "   (1) list all files in directory and all its subdirectories (i.e. dir=./ (not the default))\n"
										+ "  OR\n"
										+ "   (1) count the number of lines in the file (i.e. wc (not the default))\n"
										+ "   (2) filename to count (i.e. file=large_file.txt (not the default))\n"
										+ "  OR\n"
										+ "   (1) Replace all instances of a set of Strings in a file (i.e. swap=replacements.txt (not the default))\n"
										+ "       [requires two tab delimited columns of what to search for (first col) and what to replace it with (second col)]\n"
										+ "   (2) name of input file (i.e. file=input.txt (not the default))\n"
										+ "   (3) name of output file (i.e. out=output.txt (not the default))\n"
										+ "  OR\n"
										+ "   (1) filter a file by a set of keys (i.e. -filter (not the default))\n"
										+ "   (2) input data file (i.e. file=plink.tdt (not the default))\n"
										+ "   (3) key file (i.e. keys=gwasHits.txt (not the default))\n"
										+ "   (4) output file name (i.e. out=gwasHits.tdt (not the default))\n"
										+ "   (5) index of key column in data file (i.e. keyIndex=7 (not the default))\n"
										+ "   (6) order results by key order (i.e. ordered=TRUE (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-multiple")) {
				multiple = true;
				numArgs--;
			} else if (arg.startsWith("start=")) {
				start = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("stop=")) {
				stop = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("separate=")) {
				separate = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("patterns=")) {
				patterns = arg.substring(arg.indexOf("=") + 1).split(",");
				numArgs--;
			} else if (arg.startsWith("-cwd")) {
				cwd = true;
				numArgs--;
			} else if (arg.startsWith("-nextRep")) {
				findNextRep = true;
				numArgs--;
			} else if (arg.startsWith("-filter")) {
				filter = true;
				numArgs--;
			} else if (arg.startsWith("wc")) {
				wc = true;
				numArgs--;
			} else if (arg.startsWith("lastRep=")) {
				lastKnownRep = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("wait=")) {
				patienceInMilliseconds = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("transpose=")) {
				transpose = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("commaDelimitedIn=")) {
				commaDelimitedIn = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("commaDelimitedOut=")) {
				commaDelimitedOut = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("out=")) {
				outfile = ext.parseStringArg(arg, null);
				numArgs--;
			} else if (arg.startsWith("keys=")) {
				keyfile = ext.parseStringArg(arg, null);
				numArgs--;
			} else if (arg.startsWith("keyIndex=")) {
				keyIndex = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("ordered=")) {
				ord = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("swap=")) {
				replacements = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			return;
		}
		try {
			// String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\imputation\\MACH
			// comparison\\Mach_chr21\\";
			// String filename = "trans.txt";

			// String dir = "C:\\Documents and Settings\\npankrat\\My
			// Documents\\gwas\\imputation\\01_proxy_assoc\\";
			// String filename = "MACH_step2_chr21.mldose";

			// String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\imputation\\Boston
			// Method\\comparison\\minicomp\\";
			// String filename = "genotypes_chr21_CEU_r22_nr.b36_fwd.phase.txt";

			// String filename = "transBack.txt";
			// String filename = "trans.txt";
			// transpose(dir+filename);
			// transposeHuge(dir+filename, 1000);

			// findNextRep = true;
			// patterns = new String[] {"perms#.qsub"};

			// String filesMoved = "D:/data/GEDI/all.out";
			// String directory = "D:/data/GEDI/penn_data/";
			// moveFilesMoved(filesMoved, directory);
			// System.exit(1);

			log = new Logger(logfile);

			if (filter) {
				filterByKeys(keyfile, filename, outfile, keyIndex, ord);
			} else if (wc) {
				long time = new java.util.Date().getTime();
				log.report("Counted " + countLines(args[0], 0) + " lines in " + ext.getTimeElapsed(time));
			} else if (findNextRep && patterns != null) {
				System.out.println(findNextRepSafely(	patterns, numDigits, lastKnownRep,
																							patienceInMilliseconds));
			} else if (transpose != null) {
				transpose(transpose, commaDelimitedIn ? "," : "[\\s]+", outfile,
									commaDelimitedOut ? "," : "[\\s]+", log);
			} else if (replacements != null) {
				replaceAll(filename, outfile, replacements, log);
			} else if (filename != null) {
				makeQsub(	new File(filename).getAbsolutePath(), multiple, start, stop, separate, patterns,
									cwd);
			} else if (dir != null) {
				summarizeAllFilesInDirectory(dir);
			} else {
				log.reportError(usage);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
