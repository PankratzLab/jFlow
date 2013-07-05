package gwas;

import java.io.*;
import java.util.*;
import common.*;

public class SummarizePhenotype {
	public static void parse(String filename, String idList, Logger log) {
		String dbFile = null;
		String[][] parameters;
		
		if (idList == null) {
			dbFile = filename;
		} else {
			// create dbFile by filtering filename using the idList
		}
		
		parameters = new String[][] {
				{"# with genotypes and ICAM levels", "count", "icam", "icam=ValidDouble", "BMI=ValidDouble", "-blank"},
				{"mean ICAM levels", "mean", "icam", "icam=ValidDouble", "BMI=ValidDouble", "-blank", "sf=2", "-stdev"},
				{"mean age", "mean", "Age", "icam=ValidDouble", "BMI=ValidDouble", "-blank", "sf=2", "-stdev"},
				{"% male", "mean", "Male", "icam=ValidDouble", "BMI=ValidDouble", "-blank", "sf=2", "-percent"},
				{"mean BMI", "mean", "BMI", "icam=ValidDouble", "BMI=ValidDouble", "-blank", "sf=2", "-stdev"},
		};
		
		
		Files.generateTables(dbFile, new String[] {dbFile}, new String[] {ext.rootOf(filename)}, parameters, log);
		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "pheno.csv";
		String idList = null;
		String logfile = null;
		Logger log;

		String usage = "\n" + 
		"gwas.SummarizePhenotype requires 0-1 arguments\n" + 
		"   (1) name of phenotype file (i.e. file=" + filename + " (default))\n" + 
		"   (2) (optional) name of file with list of IDs actually used (i.e. idList=idsUsed.txt (not the default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("pheno=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("idList=")) {
				filename = args[i].split("=")[1];
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
			parse(filename, idList, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
