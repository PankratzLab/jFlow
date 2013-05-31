package common;

import java.io.*;

public class Compress {
	public static void compressAllTextFilesInDirectory(String dirin, String dirout, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] files;
		
		if (dirin == null) {
			log.reportError("Error - dir cannot be null");
			return;
		}

		if (!Files.exists(dirin)) {
			log.reportError("Error - invalid directory: "+dirin);
			return;
		}
		
		dirin = ext.verifyDirFormat(dirin);
		if (dirout == null) {
			dirout = dirin+"compressed/";
		}
		if (!Files.exists(dirout)) {
			new File(dirout).mkdirs();
		}
		
		files = Files.list(dirin, null, false);
		log.report("Found "+files.length+" files to compress");
		for (int i = 0; i < files.length; i++) {
			log.report((i+1)+" of "+files.length+"\t"+files[i]);
			try {
				reader = new BufferedReader(new FileReader(dirin+files[i]));
				writer = Files.getAppropriateWriter(dirout+files[i]+".gz");
				while (reader.ready()) {
					writer.println(reader.readLine());
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + dirin+files[i] + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + dirin+files[i] + "\"");
				System.exit(2);
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dirin = "./";
		String dirout = null;
		String logfile = null;
		Logger log;

		String usage = "\n" + 
		"common.Compress requires 0-1 arguments\n" + 
		"   (1) directory containing files to be gzipped (i.e. dirin=" + dirin + " (default))\n" + 
		"   (2) directory where the compressed files should be written to (i.e. dirout=[dirin]/compressed/ (default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dirin=")) {
				dirin = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dirout=")) {
				dirin = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			
			dirin = "I:/GEDI/02_After_reclustering/00src/";
			dirout = "C:/GEDI_compressed/";
			log = new Logger(logfile);
			compressAllTextFilesInDirectory(dirin, dirout, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
