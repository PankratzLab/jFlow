package org.genvisis.qsub;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import org.genvisis.cnv.gui.QueuePicker;
import org.genvisis.cnv.gui.UITools;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.PSF;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

import com.google.common.primitives.Ints;

public class Qsub {

	/**
	 * Display a GUI with fields for memory, wall-time, and # of processors.  <br /><br />
	 * GUI is populated with info from the default properties file ({@link QueueProperties.PROPERTIES_FILE}).
	 * 
	 * @param filename Suggested filename of qsub file to be created
	 * @param command Single string of command(s) to be written to a qsub file
	 * @return Actual filename of created file
	 */
	public static String qsubGUI(String filename, String command) {
		return qsubGUI(QueueProperties.PROPERTIES_FILE, filename, command);
	}
	
	/**
	 * Display a GUI with fields for memory, wall-time, and # of processors.  <br /><br />
	 * GUI is populated with info from the given properties file.  
	 * 
	 * @param queuePropFile Queue properties file - it is intended that this point to a 
	 * single queue.properties file, ideally located wherever the Genvisis executable is located.
	 * @param filename Filename of qsub file to be created
	 * @param command Single string of command(s) to be written to a qsub file
	 */
	public static String qsubGUI(String queuePropFile, String filename, String command) {
		QueuePicker qp = new QueuePicker(filename);
		UITools.centerComponent(qp);
		List<JobQueue> queues = QueueProperties.getJobQueues(queuePropFile);
		qp.populate(queues);
		qp.setModal(true);
		qp.setVisible(true);
		
		if (qp.wasCancelled()) {
			qp.dispose();
			return null;
		}
		
		String file = qp.getFilename();
		int mem = qp.getMemoryInMb();
		int wall = qp.getWalltimeHours();
		int proc = qp.getProcessors();
		
		qp.dispose();
		
		qsub(file, command, mem, wall, proc);
		return file;
	}
	
	/**
	 * Create a qsub file. <br /><br />
	 * Uses values for mem/proc/walltime based on the default JobQueue as defined
	 *  in the default queue properties file ({@link QueueProperties.PROPERTIES_FILE}).  
	 * If this file is missing and default values cannot be parsed from available 
	 * queue information, Genvisis defaults will be used.
	 * 
	 * @param filename Filename of qsub file to be created
	 * @param command Single string of command(s) to be written to a qsub file
	 */
	public static void qsubDefaults(String filename, String command) {
		qsubDefaults(QueueProperties.PROPERTIES_FILE, filename, command);
	}
	
	/**
	 * Create a qsub file. <br /><br />
	 * Uses values for mem/proc/walltime based on the default JobQueue as defined
	 *  in the specified queue properties file.  
	 * If this file is missing and default values cannot be parsed from available 
	 * queue information, Genvisis defaults will be used.
	 * 
	 * @param queuePropFile Queue properties file - it is intended that this point to a 
	 * single queue.properties file, ideally located wherever the Genvisis executable is located.
	 * @param filename Filename of qsub file to be created
	 * @param command Single string of command(s) to be written to a qsub file
	 */
	public static void qsubDefaults(String queuePropFile, String filename, String command) {
		List<JobQueue> queues = QueueProperties.getJobQueues(queuePropFile);
		JobQueue def = QueuesParser.findSensibleDefault(queues);
		qsubDefaults(def, filename, command);
	}
	
	/**
	 * Create a qsub file using the default mem/wall/proc values in the given JobQueue.  
	 * If these values are missing or not set, use Genvisis defaults.
	 * 
	 * @param queue A JobQueue
	 * @param filename Filename of qsub file to be created
	 * @param command Single string of command(s) to be written to a qsub file
	 */
	public static void qsubDefaults(JobQueue queue, String filename, String command) {
		int mem = Files.PBS_MEM;
		int proc = Files.PBS_PROC;
		int wall = Files.PBS_WALL;
		if (queue != null) {
			if (queue.getDefaultMem() != -1) {
				mem = (int) queue.getDefaultMem();
			}
			if (queue.getDefaultProc() != -1) {
				proc = queue.getDefaultProc();
			}
			if (queue.getDefaultWalltime() != -1) {
				wall = queue.getDefaultWalltime();
			}
		}
		qsub(filename, command, mem, wall, proc);
	}
	
	public static void makeQsub(String filename, boolean multiple, int start, int stop,
															boolean separate, String[] patterns,
															boolean changeToCurrentWorkingDirectoryFirst) {
		String[] lines, qsubs;
		List<String> v;
		int numThreads;
	
		lines = HashVec.loadFileToStringArray(filename, false, false, null, false);
	
		if (multiple) { // could be more elegant and incorporated with those below if desired
			System.out.println("Creating a ScriptExecutor");
			numThreads = Math.min(24, Files.countLines(filename, 0));
			Qsub.qsub("scriptExecutorFor_" + ext.removeDirectoryInfo(filename),
					 "cd " + ext.parseDirectoryOfFile(filename) + "\n" + Files.getRunString()
																																		 + " one.ScriptExecutor file="
																																		 + ext.removeDirectoryInfo(filename)
																																		 + " threads=" + numThreads,
					 63000, 12, numThreads);
		} else if (separate) {
			v = new ArrayList<String>();
			for (int i = 0; i < lines.length; i++) {
				qsubs = Qsub.qsub("", ext.rootOf(filename) + (i + 1) + ".#", start, stop,
										 (changeToCurrentWorkingDirectoryFirst ? "cd " + ext.pwd() + "\n" : "")
																																						 + lines[i],
										 patterns, Files.PBS_MEM, Files.PBS_PROC, null);
				v.add(qsubs[0]);
			}
			Files.writeArray(ArrayUtils.toStringArray(v), "master." + ext.rootOf(filename));
			Files.chmod("master." + ext.rootOf(filename));
		} else {
			if (changeToCurrentWorkingDirectoryFirst) {
				lines = ArrayUtils.addStrToArray("cd " + ext.pwd(), lines);
			}
			if (start == stop) {
				Qsub.qsub(ext.rootOf(filename) + ".qsub", ArrayUtils.toStr(lines, "\n"), Files.PBS_MEM, Files.PBS_PROC, 1);
			} else {
				qsubs = Qsub.qsub("", ext.rootOf(filename) + "#", start, stop, ArrayUtils.toStr(lines, "\n"),
										 patterns, Files.PBS_MEM, Files.PBS_PROC, null);
				if (qsubs.length > 1) {
					Files.writeArray(qsubs, "master." + ext.rootOf(filename));
					Files.chmod("master." + ext.rootOf(filename));
				}
			}
		}
	}

	public static void qsub(String filename, String command, int totalMemoryRequestedInMb,
													double walltimeRequestedInHours, int numProcs) {
		PrintWriter writer;
		String[] lines;
	
		lines = command.split("\\n");
		writer = Files.getAppropriateWriter(filename);
		if (writer == null) {
			return;
		}
		
		Qsub.writeQsubHeader(writer, filename, totalMemoryRequestedInMb, walltimeRequestedInHours,	numProcs, null);
		
		boolean rewriteJavaCmd = (totalMemoryRequestedInMb / 1024) > 1; // default Java heap size is min(1/4 mem avail, 1GB)
		
		for (String line : lines) {
			if (line.startsWith("java ")) {
				if (!line.contains("-Xmx") && rewriteJavaCmd) {
					int memG = (totalMemoryRequestedInMb / 1024);
					line = line.replace("java ", "java -Xmx" + memG + "G ");
				}
			}
			writer.println(line);
		}
		writer.println("echo \"end " + ext.rootOf(filename) + " at: \" `date`");
		writer.flush();
		writer.close();
		Files.chmod(filename, false);
	}

	public static void writeQsubHeader(PrintWriter writer, String filename,
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
			writer.println("#PBS -l " + ArrayUtils.toStr(ArrayUtils.toStringArray(params), ","));
		}
	
	
		if (nodeToUse != null) {
			writer.println("#$ -q *@" + nodeToUse);
		}
		writer.println();
	
		if (defualtMods) {
			writer.println(ArrayUtils.toStr(PSF.Load.getAllModules(), "\n"));
			writer.println(PSF.Cmd.getProfilePLWithOutput(ext.rootOf(filename, false) + ".profile"));
		}
	
		writer.println();
		writer.println("echo \"start " + ext.rootOf(filename) + " at: \" `date`");
		writer.println("/bin/hostname");
	
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
	
		lines = qsub(dir, "chr#_" + root, start, stop, commands, null, memRequiredInMb,
								 walltimeRequestedInHours, null);
	
		if (lines.length > 1) {
			Files.writeArray(lines, dir + "master." + (root == null ? "qsub" : root));
			Files.chmod(dir + "master." + (root == null ? "qsub" : root));
		}
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
			System.err.println("Error - directory '" + dir
												 + "' does not exist, cannot create batches there");
		}
	
		v = new Vector<String>();
		for (int i = start; i <= stop; i++) {
			// filename = ext.parseDirectoryOfFile(prefix,
			// true)+(prefix==null?"":ext.removeDirectoryInfo(prefix))+(start!=stop||(start>0&&start<25)?(prefix==null||prefix.equals("")?i:(prefix.endsWith("chr")?"":".")+i):"")+(root==null?"":"_"+root)+".qsub";
			filename = ext.insertNumbers(filenameFormat, i)
								 + (filenameFormat.endsWith(".qsub") ? "" : ".qsub");
			try {
				writer = Files.openAppropriateWriter(dir + filename);
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
					writer.println("    rep=$(java -jar " + Files.ROOT_DIRECTORY
												 + org.genvisis.common.PSF.Java.GENVISIS
												 + " common.Files -nextRep patterns=" + ArrayUtils.toStr(patterns, ",")
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
				Files.chmod(dir + filename, false);
				v.add("qsub " + (queueName == null ? "" : "-q " + queueName + " ") + filename);
			} catch (IOException ioe) {
				throw new RuntimeException("Problem creating " + dir + filename);
			}
		}
		return ArrayUtils.toStringArray(v);
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
			writer = Files.openAppropriateWriter(chunkFilename);
	
			if (numJobsToForce <= 0) {
				numProcs = jobs.length;
			} else if (numJobsToForce < jobs.length) {
				System.err.println("Error - cannot force fewer jobs than are provided (" + numJobsToForce
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

	public static void qsubMultiple(List<String> jobNamesWithAbsolutePaths, IntVector jobSizes,
																	String batchDir, String batchRoot, int maxJobsPerBatch,
																	boolean forceMaxJobs, String queueName,
																	int memoryPerProcRequestedInMb, int totalMemoryRequestedInMb,
																	double walltimeRequestedInHours) {
		String[] files;
		int count;
		List<String> v;
	
		if (jobSizes == null) {
			files = ArrayUtils.toStringArray(jobNamesWithAbsolutePaths);
		} else {
			int[] jobOrder = Sort.getReverseIndices(Ints.toArray(jobSizes));
			files = Sort.getOrdered(jobNamesWithAbsolutePaths, jobOrder);
		}
		count = 0;
		v = new ArrayList<String>();
		new File(batchDir).mkdirs();
		v.add("cd " + batchDir);
		while (count * maxJobsPerBatch < files.length) {
			Qsub.qsubMultiple(batchDir + ext.removeDirectoryInfo(batchRoot) + "." + count,
												 ArrayUtils.subArray(files, count * maxJobsPerBatch,
																						 Math.min((count + 1) * maxJobsPerBatch, files.length)),
												 forceMaxJobs ? maxJobsPerBatch : -1, memoryPerProcRequestedInMb,
												 totalMemoryRequestedInMb, walltimeRequestedInHours);
			v.add("qsub " + (queueName == null ? "" : "-q " + queueName + " ")
						+ ext.removeDirectoryInfo(batchRoot) + "." + count);
			count++;
		}
		Files.writeArray(ArrayUtils.toStringArray(v), ext.parseDirectoryOfFile(batchRoot) + "master."
																									+ ext.removeDirectoryInfo(batchRoot));
		Files.chmod(ext.parseDirectoryOfFile(batchRoot) + "master."
								+ ext.removeDirectoryInfo(batchRoot));
	}

	public static void qsubExecutor(String dirToSwitchToBeforeRunning,
																	List<String> commandsWithAbsolutePaths, IntVector jobSizes,
																	String batchRoot, int numProcs, int totalMemoryRequestedInMb,
																	double walltimeRequestedInHours) {
		String[] commands;
	
		if (jobSizes == null) {
			commands = ArrayUtils.toStringArray(commandsWithAbsolutePaths);
		} else {
			int[] jobOrder = Sort.getReverseIndices(Ints.toArray(jobSizes));
			commands = Sort.getOrdered(commandsWithAbsolutePaths, jobOrder);
		}
	
		Files.writeArray(commands, batchRoot + ".chain");
		qsub(batchRoot + ".pbs",
							 "cd " + dirToSwitchToBeforeRunning + "\njava -jar ~/"
																	 + org.genvisis.common.PSF.Java.GENVISIS
																	 + " one.ScriptExecutor file=" + batchRoot + ".chain threads="
																	 + numProcs,
							 totalMemoryRequestedInMb, walltimeRequestedInHours, numProcs);
	}

	public static void qsub(String root_batch_name, String dirToSwitchToBeforeRunning, int numBatches,
													String commands, String[][] iterations, int totalMemoryRequestedInMb,
													double walltimeRequestedInHours) {
		qsub(root_batch_name, dirToSwitchToBeforeRunning, numBatches, commands, iterations,
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
				writer = Files.openAppropriateWriter(ext.parseDirectoryOfFile(root_batch_name)
																								+ "master."
																								+ ext.removeDirectoryInfo(root_batch_name));
				if (!dir.equals("./") && !dir.equals("")) {
					writer.println("cd " + dir);
				}
				for (int i = 0; i < numBatches; i++) {
					writers[i] = Files.openAppropriateWriter(root_batch_name + "_" + (i + 1) + ".qsub");
					writeQsubHeader(writers[i], root_batch_name + "_" + (i + 1) + ".qsub",
													totalMemoryRequestedInMb, walltimeRequestedInHours, 1, null);
					if (dirToSwitchToBeforeRunning != null) {
						writers[i].println("cd " + dirToSwitchToBeforeRunning);
					}
					writer.println("qsub " + (queueName == null ? "" : "-q " + queueName + " ")
												 + ext.removeDirectoryInfo(root_batch_name) + "_" + (i + 1) + ".qsub");
				}
				writer.close();
				Files.chmod(ext.parseDirectoryOfFile(root_batch_name) + "master."
							+ ext.removeDirectoryInfo(root_batch_name));
			} else {
				// writer = null;
				writers[0] = Files.openAppropriateWriter(root_batch_name + ".qsub");
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
				writers[i].println("echo \"end " + ext.removeDirectoryInfo(root_batch_name) + "_" + (i + 1)
													 + " at: \" `date`");
				writers[i].close();
				Files.chmod(root_batch_name + "_" + (i + 1) + ".qsub");
			}
			if (numBatches == 1) {
				Files.chmod(root_batch_name + ".qsub");
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
		qsub(root, commands, iterations, totalMemoryRequestedInMb, walltimeRequestedInHours, numProcs,
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
			writer = Files.openAppropriateWriter(filename);
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
	
}
