package stats;

//import java.io.*;
import java.util.*;
import common.*;

public class Rscript {
	public static final String[] RSCRIPT_EXECS = {
//		"/soft/R/3.0.1/bin/Rscript", // MSI
		"/soft/R/2.15.1/bin/Rscript", // Itasca nodes can only see this
		"/share/apps/R-3.0.1/bin/Rscript", // psych
		"/share/apps/src/R-3.0.1/bin/Rscript", // alcatraz
	};
	
	public static final String[] R_EXECS = {
//		"/soft/R/3.0.1/bin/R", // MSI
		"/soft/R/2.15.1/bin/R", // Itasca nodes can only see this
		"/share/apps/R-3.0.1/bin/R", // psych
		"/share/apps/src/R-3.0.1/bin/R", // alcatraz
	};
	
	public static String getRscriptExecutable(Logger log) {
		if (System.getProperty("os.name").startsWith("Windows")) {
			return "Rscript";
		} else {
			for (int i = 0; i < RSCRIPT_EXECS.length; i++) {
				if (Files.exists(RSCRIPT_EXECS[i])) {
					return RSCRIPT_EXECS[i];
				}
			}
		}

		log.report("Warning - could not determine hostname, assuming R executbale is simply 'Rscript'");
		
		return "Rscript";
	}

	public static String getRExecutable(Logger log) {
		if (System.getProperty("os.name").startsWith("Windows")) {
			return "R";
		} else {
			for (int i = 0; i < R_EXECS.length; i++) {
				if (Files.exists(R_EXECS[i])) {
					return R_EXECS[i];
				}
			}
		}

		log.report("Warning - could not determine hostname, assuming R executbale is simply 'Rscript'");
		
		return "R";
	}
	
	private static void batchUp(String dir, Logger log) {
		Vector<String> v;
		String[] files;
		String root;
		
		dir = ext.verifyDirFormat(dir);
		
		if (dir.equals("") || dir.equals("./")) {
			dir = ext.pwd();
		}
		
		files = Files.list(dir, ".R", false);
		
		v = new Vector<String>();
		for (int i = 0; i < files.length; i++) {
			root = ext.rootOf(dir+files[i], false);
			Files.qsub(root+".qsub", getRscriptExecutable(log)+" --no-save "+root+".R", 4000, 2, 1);
			v.add("qsub "+root+".qsub");
		}
		
		Files.writeList(Array.toStringArray(v), dir+"master");
		Files.chmod(dir+"master");
	}
	
	/**
	 * @param toVector
	 *            String[] that will be converted to c("toVector[0]",toVector[1]...);
	 * @return
	 */
	public static String generateRVector(String[] toVector) {
		String rV = "c(";
		for (int i = 0; i < toVector.length; i++) {
			rV += (i == 0 ? "" : ",") +"\""+ toVector[i] +"\"";
		}
		rV += ")\n";
		return rV;
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		boolean batchUp = false;
		String dir = "./";
		String logfile = null;
		Logger log;

		String usage = "\n" + 
				"stats.Rscript requires 0-1 arguments\n" + 
				"   (1) create qsub files for all .R files in the directory (i.e. -batchUp (not the default))\n" + 
				"   (2) directory to search for the .R files (i.e. dir="+dir+" (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("-batchUp")) {
				batchUp = true;
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
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
			log = new Logger(logfile);
			if (batchUp) {
				batchUp(dir, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
