package common;

import java.io.*;

public class CmdLine {
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
			if (System.getProperty("os.name").startsWith("Windows") && (command.contains(">") || command.contains("|"))) {
				System.err.println("FYI - the Runtime.exec command will likely not work, since it contains a pipe, write command to a .bat file and exec that instead");
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

			while (!finish) {
				try {
					while (in.available()>0||err.available()>0) {
						while (in.available()>0) {
							if (log != null) {
								log.report((char) in.read()+"",false,true);
							} else if (inOs != null) {
								inOs.print((char) in.read());
							} else {
								in.read();
							}
						}
						while (err.available() > 0) {
							if (log != null) {
								log.report((char) err.read() + "",false,true);
							} else if (errOS != null) {
								errOS.print((char) err.read());
							} else {
								err.read();
							}
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
				System.err.println("Error - The program \""+message.substring(0, message.indexOf("\""))+"\" is not installed or is not accessible "+message.substring(message.indexOf("("), message.indexOf(")")+1));
				noError = false;
			} else {
				ioe.printStackTrace();
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
		if (expectedOutputFiles == null || !Files.exists(dir, expectedOutputFiles) || (Files.exists(dir, expectedOutputFiles) && overWriteExistingOutput)) {
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
}
