package org.genvisis.common;

import java.io.*;
import java.util.Map;

public class CmdLine {
    public static boolean runDefaults(String command, String dir) {
        return run(command, null, dir, System.out, System.err, null, false);
    }
    public static boolean runDefaults(String command, String dir, Logger log) {
        return run(command, null, dir, System.out, System.err, log, false);
    }
    public static boolean run(String command, String dir) {
		return run(command, dir, null);
	}
	public static boolean run(String command, String dir, PrintStream os) {
		return run(command, dir, os, false);
	}
	public static boolean run(String command, String dir, PrintStream os, boolean ignoreIllegalStateExceptions) {
		return run(command, null, dir, os, os, null, ignoreIllegalStateExceptions);
	}
	
	public static boolean run(String command, String[] commandArray, String dir, PrintStream inOs, PrintStream errOS, Logger log, boolean ignoreIllegalStateExceptions) {
		Process proc;
		InputStream in, err;
//        PrintWriter writer;
        boolean finish;
        boolean noError;
        
        noError = true;
		
        if (dir.equals("")) {
        	dir = "./";
        }
        
        try {
			if (Files.isWindows() && (command.contains(">") || command.contains("|"))) {
				log.reportError("FYI - the Runtime.exec command will likely not work, since it contains a pipe, write command to a .bat file and exec that instead");
			}
			if (commandArray != null) {
				proc = Runtime.getRuntime().exec(commandArray, null, new File(dir));
			} else {
				proc = Runtime.getRuntime().exec(command, null, new File(dir));
			}
//			if (logfile != null) {
//				writer = new PrintWriter(new FileWriter(dir+logfile));
//			} else {
//				writer = null;
//			}
			in = proc.getInputStream();
			err = proc.getErrorStream();
			finish = false;
			byte[] b;
			while (!finish) {
				try {
					while (in.available()>0||err.available()>0) {
						while (in.available()>0) {
						    b = new byte[in.available()];
						    in.read(b);
							if (log != null) {
								log.report(new String(b, "UTF-8"),false,true);
							} else if (inOs != null) {
								inOs.print(new String(b, "UTF-8"));
							}/* else {
							}*/
							b = null;
						}
						while (err.available() > 0) {
						    b = new byte[err.available()];
						    err.read(b);
							if (log != null) {
								log.report(new String(b, "UTF-8"),false,true);
							} else if (errOS != null) {
								errOS.print(new String(b, "UTF-8"));
							}/* else {
							}*/
							b = null;
						}
					}
					proc.exitValue();
					finish = true;
				} catch (IllegalThreadStateException e) {
					try {
						Thread.sleep(10);
					} catch (InterruptedException e2) {}
					if (ignoreIllegalStateExceptions) {
//						System.err.println("Ignoring illegal state exception");
						finish = true;
						noError = false;
					}
				}
			}
//			if (writer != null) {
//				writer.close();
//			}
		} catch (IOException ioe) {
			String message;
			
			message = ioe.getMessage();
			if (message.startsWith("Cannot run program ")) {
				message = message.substring(20);
				log.reportError("Error - The program \""+message.substring(0, message.indexOf("\""))+"\" is not installed or is not accessible "+message.substring(message.indexOf("("), message.indexOf(")")+1));
				noError = false;
			} else {
				log.reportException(ioe);
				noError = false;
			}
		}
		
		return noError;
	}

	/**
	 * @param commandArray
	 *            the String array of commands with associated commands grouped at a given index (String[] commands = new String[]{"myFavoriteCommand","input=one.txt","output=2.txt"}
	 * @param dir
	 *            directory to run the command in, and also directory to check for existing files
	 * @param neccesaryInputFiles
	 *            will check these files for existence, will fail if they do not all exist
	 * @param expectedOutputFiles
	 *            will check these files for existence, will skip the command if all exist and overWriteExistingOutput is not flagged
	 * @param verbose
	 * @param overWriteExistingOutput
	 * @param skipReporting
	 *            if set to true, the cmdline reporting will be sent to the log, otherwise it will be ignored. Nice to report when error checking, but there is a high probability of reduced performance and strange output if always used
	 * @param log
	 * @return
	 */
	public static boolean runCommandWithFileChecks(String[] commandArray, String dir, String[] neccesaryInputFiles, String[] expectedOutputFiles, boolean verbose, boolean overWriteExistingOutput, boolean skipReporting, Logger log) {
		boolean success = false;
		if (expectedOutputFiles == null || !Files.exists(dir, expectedOutputFiles, true) || (Files.exists(dir, expectedOutputFiles) && overWriteExistingOutput)) {
			if (neccesaryInputFiles == null || Files.exists(dir, neccesaryInputFiles)) {
				if (verbose) {
					log.report(ext.getTime() + " Info - running command " + Array.toStr(commandArray, " "));
				}
				if (run(Array.toStr(commandArray, " "), commandArray, dir, null, null, (skipReporting ? null : log), false)) {
					if (expectedOutputFiles != null && !Files.exists(dir, expectedOutputFiles)) {
						log.reportError("Error - the command " + Array.toStr(commandArray, " ") + " appeared to run, but could not find all neccesary output files:" + Array.toStr(expectedOutputFiles, "\n"));
					} else {
						if (verbose) {
							log.report(ext.getTime() + " Info - finished running command " + Array.toStr(commandArray, " "));
						}
						success = true;
					}
				} else {
					log.reportError("Error - the command " + Array.toStr(commandArray, " ") + " has failed");
				}
			} else {
				log.reportError("Error - could not find all necessary input files:\n" + Array.toStr(neccesaryInputFiles, "\n"));
			}
		} else {
			if (verbose) {
				log.report(ext.getTime() + " Info - all of the expected output files exist and the overwrite option was not flagged, skipping:");
				log.report("COMMAND SKIPPED: " + Array.toStr(commandArray, " "));
			}
			success = true;
		}
		return success;
	}
	
	/**
	 * @param commandArray
	 *            an array representing the command to run
	 * @param batFile
	 *            where the the command will be written
	 * @param verbose
	 *            report the command written to the batFile
	 * @param log
	 * 
	 * 
	 * @return String[] of the batFile
	 */

	public static String[] prepareBatchForCommandLine(String[] commandArray, String batFile, boolean verbose, Logger log) {
		if (verbose) {
			log.report(ext.getTime() + " Info - running command " + Array.toStr(commandArray, " ") + "\nUsing file " + batFile);
		}
		Files.write(Array.toStr(commandArray, " "), batFile);
		Files.chmod(batFile);
		return new String[] { batFile };
	}
	
	public static String getCmdLocation(String commmand) {
		Map<String, String> env = System.getenv();
		 for (String envName : env.keySet()) {
	            System.out.format("%s=%s%n",
	                              envName,
	                              env.get(envName));
	        }
		
		if (env.containsKey(commmand)) {
			return env.get(commmand);
		} else {
			return null;
		}
	}
	
	

	
}
