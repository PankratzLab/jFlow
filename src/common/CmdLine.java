package common;

import java.io.*;

public class CmdLine {
	public static void stdin() {
		try {
			new BufferedReader(new InputStreamReader(System.in)).readLine();
		} catch (IOException ioe) {
			System.err.println("Error trying to get input");
		}		
	}
	
	public static void run(String command, String dir) {
		run(command, dir, null);
	}
	
	public static void run(String command, String dir, PrintStream os) {
		run(command, dir, os, false);
	}
	
	public static boolean run(String command, String dir, PrintStream os, boolean ignoreIllegalStateExceptions) {
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
			proc = Runtime.getRuntime().exec(command, null, new File(dir));
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
							if (os != null) {
								os.print((char)in.read());
							} else {
								in.read();
							}
						}
						while (err.available()>0) {
							if (os != null) {
								os.print((char)err.read());
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
			} else {
				ioe.printStackTrace();
			}
		}
		
		return noError;

	}
}
