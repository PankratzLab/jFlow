package org.genvisis.one;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

public class ScriptExecutor {

	ExecutorService executor;

	public ScriptExecutor(int numThreads) {
		int threads = numThreads;
		if (threads == -1) {
			threads = Runtime.getRuntime().availableProcessors();
		}
		executor = Executors.newFixedThreadPool(numThreads);
	}

	public static boolean outLogExistsComplete(String outLog, String outLogEndToken) {
		if (!Files.exists(outLog)) {
			return false;
		}
		String[] lines = HashVec.loadFileToStringArray(outLog, false, null, false);
		return !(outLogEndToken == null) && lines.length > 0
					 && (lines[lines.length - 1].contains(outLogEndToken));
	}

	public void run(String file, String outLogEndToken) throws IOException {
		String dirTemp = ext.parseDirectoryOfFile(file);
		if (!dirTemp.endsWith("/")) {
			dirTemp = dirTemp + "/";
		}
		final String dir = dirTemp;
		final String[] lines = HashVec.loadFileToStringArray(file, false, null, false);
		final String fileRoot = ext.rootOf(file, true);
		final PrintWriter outWriter = Files.openAppropriateWriter(dir + fileRoot + ".outlog");

		new File(dir + "output/").mkdirs();
		for (int i = 0; i < lines.length; i++) {
			final int index = i;
			final String cmdLine = lines[index];

			String outLogTemp = "output/" + fileRoot + ".log_" + index + ".out";
			boolean foundCompleteOutLog = outLogExistsComplete(dir + outLogTemp, outLogEndToken);
			if (foundCompleteOutLog) {
				outWriter.println("Ignoring runnable " + index + "/" + lines.length + "; output logfile ["
													+ outLogTemp + "] contains given token [" + outLogEndToken + "]");
				continue;
			}
			final String outLog = ensureWriteableOutLog(dir + "output/" + fileRoot + ".log_" + index);
			Runnable executable = new Runnable() {
				@Override
				public void run() {
					try {
						(new File(outLog)).createNewFile();
						PrintStream ps = new PrintStream(outLog);
						CmdLine.run(cmdLine, dir, ps);
					} catch (Exception e) {
						e.printStackTrace(outWriter);
					}
				}
			};
			outWriter.println("Queuing runnable " + index + "/" + lines.length);
			try {
				executor.execute(executable);
			} catch (Exception e) {
				e.printStackTrace(outWriter);
			}
		}
		executor.shutdown();
		try {
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		outWriter.flush();
		outWriter.close();
	}


	private String ensureWriteableOutLog(String outLog) {
		String tempLogFile = outLog + ".out";
		File outLogFile = new File(tempLogFile);
		if (outLogFile.exists() && !outLogFile.canWrite()) {
			int cnt = 1;
			do {
				tempLogFile = outLog + "_" + cnt + ".out";
				outLogFile = new File(tempLogFile);
				cnt++;
			} while (outLogFile.exists() && !outLogFile.canWrite());
		}
		return tempLogFile;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String inputFile = "input.txt";
		String outlogToken = null;
		int numThreads = 1;

		String usage = "one.ScriptExecutor requires 3 arguments\n"
									 + "   (1) filename with one command per line (i.e. file="
									 + inputFile
									 + " (default))\n"
									 + "   (2) token indicating successful completion in last line of log files (i.e. token="
									 + outlogToken + " (default))\n" + "   (3) number of threads (i.e. threads="
									 + numThreads + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				inputFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("token=")) {
				outlogToken = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("threads=")) {
				numThreads = Integer.valueOf(arg.split("=")[1]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			ScriptExecutor executor = new ScriptExecutor(numThreads);
			executor.run(inputFile, outlogToken);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


}
