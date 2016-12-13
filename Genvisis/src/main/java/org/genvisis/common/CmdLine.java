package org.genvisis.common;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class CmdLine {
	public static class Command {
		List<String> commandList;
		Collection<String> necessaryInputFiles;
		Collection<String> expectedOutputFiles;
		String dir;

		/**
		 * @param commandList the commands as a List of Strings, spaces will be inserted between each
		 *        element
		 * @param necessaryInputFiles check these files for existence, will fail if they do not all
		 *        exist
		 * @param expectedOutputFiles check these files for existence, will skip the command if all
		 *        exist
		 * @param dir directory to run the command in, and also directory to check for existing files
		 */
		public Command(	List<String> commandList, Collection<String> necessaryInputFiles,
										Collection<String> expectedOutputFiles, String dir) {
			super();
			this.commandList = commandList;
			this.necessaryInputFiles = necessaryInputFiles;
			this.expectedOutputFiles = expectedOutputFiles;
			this.dir = dir;
		}

		/**
		 * @param commandArray the String array of commands with associated commands grouped at a given
		 *        index (String[] commands = new
		 *        String[]{"myFavoriteCommand","input=one.txt","output=2.txt"}
		 * @param necessaryInputFiles will check these files for existence, will fail if they do not all
		 *        exist
		 * @param expectedOutputFiles will check these files for existence, will skip the command if all
		 *        exist
		 * @param dir directory to run the command in, and also directory to check for existing files
		 */
		public Command(	String[] commandArray, String[] necessaryInputFiles,
										String[] expectedOutputFiles, String dir) {
			this(	Arrays.asList(commandArray),
						necessaryInputFiles == null ? null : Arrays.asList(necessaryInputFiles),
						expectedOutputFiles == null ? null : Arrays.asList(expectedOutputFiles), dir);
		}

		/**
		 * @param commandArray the String array of commands with associated commands grouped at a given
		 *        index (String[] commands = new
		 *        String[]{"myFavoriteCommand","input=one.txt","output=2.txt"}
		 * @param necessaryInputFiles will check these files for existence, will fail if they do not all
		 *        exist
		 * @param expectedOutputFiles will check these files for existence, will skip the command if all
		 *        exist
		 */
		public Command(	String[] commandArray, String[] necessaryInputFiles,
										String[] expectedOutputFiles) {
			this(commandArray, necessaryInputFiles, expectedOutputFiles, "");
		}

		public boolean runCommand(Logger log) {
			return runCommand(true, false, false, true, log);
		}

		/**
		 * @param verbose
		 * @param overWriteExistingOutput
		 * @param skipReporting if set to false, the cmdline reporting will be sent to the log,
		 *        otherwise it will be ignored. Nice to report when error checking, but there is a high
		 *        probability of reduced performance and strange output if always used
		 * @param log
		 * @return true on success, false on failure
		 */
		public boolean runCommand(boolean verbose, boolean overWriteExistingOutput,
															boolean skipReporting, boolean treatEmptyAsMissing, Logger log) {
			boolean success = false;
			if (expectedOutputFiles == null
						|| !Files.exists(dir, expectedOutputFiles, treatEmptyAsMissing)
					|| overWriteExistingOutput) {
				if (necessaryInputFiles == null
						|| Files.exists(dir, necessaryInputFiles, treatEmptyAsMissing)) {
					if (verbose) {
						log.report(ext.getTime()	+ " Info - running command "
												+ IterableUtils.toStr(commandList, " "));
					}
					if (run(commandList, dir, null, null, (skipReporting ? null : log), false)) {
						if (expectedOutputFiles != null
								&& !Files.exists(dir, expectedOutputFiles, treatEmptyAsMissing)) {
							log.reportError("Error - the command "	+ IterableUtils.toStr(commandList, " ")
															+ " appeared to run, but could not find all necessary output files:"
															+ IterableUtils.toStr(expectedOutputFiles, "\n"));
						} else {
							if (verbose) {
								log.report(ext.getTime()	+ " Info - finished running command "
														+ IterableUtils.toStr(commandList, " "));
							}
							success = true;
						}
					} else {
						log.reportError("Error - the command "	+ IterableUtils.toStr(commandList, " ")
														+ " has failed");
					}
				} else {
					log.reportError("Error - could not find all necessary input files:\n"
													+ IterableUtils.toStr(necessaryInputFiles, "\n"));
				}
			} else {
				if (verbose) {
					log.report(ext.getTime()
											+ " Info - all of the expected output files exist and the overwrite option was not flagged, skipping:");
					log.report("COMMAND SKIPPED: " + IterableUtils.toStr(commandList, " "));
				}
				success = true;
			}
			return success;
		}



	}

	public static String getCmdLocation(String commmand) {
		Map<String, String> env = System.getenv();
		for (String envName : env.keySet()) {
			System.out.format("%s=%s%n", envName, env.get(envName));
		}

		if (env.containsKey(commmand)) {
			return env.get(commmand);
		} else {
			return null;
		}
	}

	/**
	 * @param commandArray an array representing the command to run
	 * @param batFile where the the command will be written
	 * @param verbose report the command written to the batFile
	 * @param log
	 *
	 *
	 * @return String[] of the batFile
	 */

	public static String[] prepareBatchForCommandLine(String[] commandArray, String batFile,
																										boolean verbose, Logger log) {
		if (verbose) {
			log.report(ext.getTime()	+ " Info - running command " + Array.toStr(commandArray, " ")
									+ "\nUsing file " + batFile);
		}
		Files.write(Array.toStr(commandArray, " "), batFile);
		Files.chmod(batFile);
		return new String[] {batFile};
	}

	public static boolean run(String command, String dir) {
		return run(command, dir, null);
	}

	public static boolean run(String command, String dir, PrintStream os) {
		return run(command, dir, os, false);
	}

	public static boolean run(String command, String dir, PrintStream os,
														boolean ignoreIllegalStateExceptions) {
		return run(command, dir, os, os, new Logger(), ignoreIllegalStateExceptions);
	}

	public static boolean run(Collection<String> commands, String dir, PrintStream inOs,
														PrintStream errOS, Logger log, boolean ignoreIllegalStateExceptions) {
		return run(	commands.toArray(new String[commands.size()]), dir, inOs, errOS, log,
								ignoreIllegalStateExceptions);
	}

	public static boolean run(String command, String dir, PrintStream inOs, PrintStream errOS,
														Logger log, boolean ignoreIllegalStateExceptions) {
//		StringTokenizer st = new StringTokenizer(command, " \t\n\r\f");
//		String[] cmdarray = new String[st.countTokens()];
//		for (int i = 0; st.hasMoreTokens(); i++) {
//			cmdarray[i] = st.nextToken();
//		}
		String regex = "[\"\']([^\"\']*)[\"\']|(\\S+)";
		Matcher m = Pattern.compile(regex).matcher(command);
		ArrayList<String> cmdList = new ArrayList<String>();
		while (m.find()) {
			if (m.group(1) != null) {
				cmdList.add(m.group(1));
			} else {
				cmdList.add(m.group(2));
			}
		}
		String[] cmdarray = cmdList.toArray(new String[cmdList.size()]);
		return run(cmdarray, dir, inOs, errOS, log, ignoreIllegalStateExceptions);
	}

	public static boolean run(String[] commandArray, String dir, PrintStream inOs, PrintStream errOS,
														Logger log, boolean ignoreIllegalStateExceptions) {
		Process proc;
		InputStream in, err;
		// PrintWriter writer;
		boolean finish;
		boolean noError;
		String charSet = "UTF-8";

		noError = true;

		if (dir.equals("")) {
			dir = "./";
		}

		for (String command : commandArray) {
			if (command.contains(">") || command.contains("|")) {
				if (Files.isWindows()) {
					log.reportError("FYI - the Runtime.exec command will likely not work, since it contains a pipe or redirect, write command to a file and exec that instead");
					break;
				} else if (!commandArray[0].startsWith("/bin/")) { 
					log.reportTimeWarning("FYI - the command may not work as it contains a pipe or redirect and isn't prefaced by a shell invocation (such as \"/bin/bash\").  Try prefacing the command with \"/bin/bash -c\" and wrapping the command in quotes.");
					break;
				}
			}
		}
		
		try {
			proc = Runtime.getRuntime().exec(commandArray, null, new File(dir));
			// if (logfile != null) {
			// writer = new PrintWriter(new FileWriter(dir+logfile));
			// } else {
			// writer = null;
			// }
			in = proc.getInputStream();
			err = proc.getErrorStream();
			finish = false;
			byte[] b;
			while (!finish) {
				try {
					while (in.available() > 0 || err.available() > 0) {
						while (in.available() > 0) {
							b = new byte[in.available()];
							in.read(b);
							if (inOs != null) {
								inOs.print(new String(b, charSet));
							}
							if (log != null) {
								log.report(new String(b, charSet), false, true);
							} /*
								 * else { }
								 */
							b = null;
						}
						while (err.available() > 0) {
							b = new byte[err.available()];
							err.read(b);
							if (errOS != null) {
								errOS.print(new String(b, charSet));
							}
							if (log != null) {
								log.report(new String(b, charSet), false, true);
							}/* else  
								 * else { }
								 */
							b = null;
						}
					}
					proc.exitValue();
					finish = true;
				} catch (IllegalThreadStateException e) {
					try {
						Thread.sleep(10);
					} catch (InterruptedException e2) {
					}
					if (ignoreIllegalStateExceptions) {
						// System.err.println("Ignoring illegal state exception");
						finish = true;
						noError = false;
					}
				}
			}
			// if (writer != null) {
			// writer.close();
			// }
		} catch (IOException ioe) {
			String message;

			message = ioe.getMessage();
			if (message.startsWith("Cannot run program ")) {
				message = message.substring(20);
				log.reportError("Error - The program \""	+ message.substring(0, message.indexOf("\""))
												+ "\" is not installed or is not accessible "
												+ message.substring(message.indexOf("("), message.indexOf(")") + 1));
				noError = false;
			} else {
				log.reportException(ioe);
				noError = false;
			}
		}

		return noError;
	}


	/**
	 * @param commandArray the String array of commands with associated commands grouped at a given
	 *        index (String[] commands = new
	 *        String[]{"myFavoriteCommand","input=one.txt","output=2.txt"}
	 * @param dir directory to run the command in, and also directory to check for existing files
	 * @param necessaryInputFiles will check these files for existence, will fail if they do not all
	 *        exist
	 * @param expectedOutputFiles will check these files for existence, will skip the command if all
	 *        exist and overWriteExistingOutput is not flagged
	 * @param verbose
	 * @param overWriteExistingOutput
	 * @param skipReporting if set to true, the cmdline reporting will be sent to the log, otherwise
	 *        it will be ignored. Nice to report when error checking, but there is a high probability
	 *        of reduced performance and strange output if always used
	 *
	 * @param log
	 * @return
	 */
	public static boolean runCommandWithFileChecks(	String[] commandArray, String dir,
																									String[] necessaryInputFiles,
																									String[] expectedOutputFiles, boolean verbose,
																									boolean overWriteExistingOutput,
																									boolean skipReporting, Logger log) {
		return runCommandWithFileChecks(commandArray, dir, necessaryInputFiles, expectedOutputFiles,
																		verbose, overWriteExistingOutput, skipReporting, true, log);
	}

	/**
	 * @param commandArray the String array of commands with associated commands grouped at a given
	 *        index (String[] commands = new
	 *        String[]{"myFavoriteCommand","input=one.txt","output=2.txt"}
	 * @param dir directory to run the command in, and also directory to check for existing files
	 * @param necessaryInputFiles will check these files for existence, will fail if they do not all
	 *        exist
	 * @param expectedOutputFiles will check these files for existence, will skip the command if all
	 *        exist and overWriteExistingOutput is not flagged
	 * @param verbose
	 * @param overWriteExistingOutput
	 * @param skipReporting if set to true, the cmdline reporting will be sent to the log, otherwise
	 *        it will be ignored. Nice to report when error checking, but there is a high probability
	 *        of reduced performance and strange output if always used
	 * @param treatEmptyAsMissing if true, 0 sized files will be set to missing. If 0 sized files are
	 *        used as placeholders, set to false
	 * @param log
	 * @return
	 */
	public static boolean runCommandWithFileChecks(	String[] commandArray, String dir,
																									String[] necessaryInputFiles,
																									String[] expectedOutputFiles, boolean verbose,
																									boolean overWriteExistingOutput,
																									boolean skipReporting,
																									boolean treatEmptyAsMissing, Logger log) {
		return new Command(	commandArray, necessaryInputFiles, expectedOutputFiles,
												dir).runCommand(verbose, overWriteExistingOutput, skipReporting,
																				treatEmptyAsMissing, log);
	}

	/**
	 * @param commandList the commands as a List of Strings, spaces will be inserted between each
	 *        element
	 * @param dir directory to run the command in, and also directory to check for existing files
	 * @param necessaryInputFiles will check these files for existence, will fail if they do not all
	 *        exist
	 * @param expectedOutputFiles will check these files for existence, will skip the command if all
	 *        exist and overWriteExistingOutput is not flagged
	 * @param verbose
	 * @param overWriteExistingOutput
	 * @param treatEmptyAsMissing if true, 0 sized files will be set to missing. If 0 sized files are
	 *        used as placeholders, set to false
	 * @param log
	 * @return
	 */
	public static boolean runCommandWithFileChecks(	List<String> commandList, String dir,
																									Collection<String> necessaryInputFiles,
																									Collection<String> expectedOutputFiles,
																									boolean verbose, boolean overWriteExistingOutput,
																									boolean treatEmptyAsMissing, Logger log) {
		return new Command(	commandList, necessaryInputFiles, expectedOutputFiles,
												dir).runCommand(verbose, overWriteExistingOutput, true, treatEmptyAsMissing,
																				log);
	}


	/**
	 * @param commandList the commands as a List of Strings, spaces will be inserted between each
	 *        element
	 * @param dir directory to run the command in, and also directory to check for existing files
	 * @param necessaryInputFiles will check these files for existence, will fail if they do not all
	 *        exist
	 * @param expectedOutputFiles will check these files for existence, will skip the command if all
	 *        exist and overWriteExistingOutput is not flagged
	 * @param verbose
	 * @param overWriteExistingOutput
	 * @param skipReporting if set to true, the cmdline reporting will be sent to the log, otherwise
	 *        it will be ignored. Nice to report when error checking, but there is a high probability
	 *        of reduced performance and strange output if always used
	 * @param treatEmptyAsMissing if true, 0 sized files will be set to missing. If 0 sized files are
	 *        used as placeholders, set to false
	 * @param log
	 * @return
	 */
	public static boolean runCommandWithFileChecks(	List<String> commandList, String dir,
																									Collection<String> necessaryInputFiles,
																									Collection<String> expectedOutputFiles,
																									boolean verbose, boolean overWriteExistingOutput,
																									boolean skipReporting,
																									boolean treatEmptyAsMissing, Logger log) {
		return new Command(	commandList, necessaryInputFiles, expectedOutputFiles,
												dir).runCommand(verbose, overWriteExistingOutput, skipReporting,
																				treatEmptyAsMissing, log);
	}

	public static boolean runDefaults(String command, String dir) {
		return run(command, dir, System.out, System.err, null, false);
	}

	public static boolean runDefaults(String command, String dir, Logger log) {
		return run(command, dir, System.out, System.err, log, false);
	}
	
	public static boolean runDefaults(Collection<String> command, String dir) {
		return run(command, dir, System.out, System.err, null, false);
	}

	public static boolean runDefaults(Collection<String> command, String dir, Logger log) {
		return run(command, dir, System.out, System.err, log, false);
	}



}
