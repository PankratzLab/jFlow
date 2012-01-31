import java.io.*;

import link.TrimFam;
import mining.Transformations;
import bioinformatics.*;
import common.*;
import parse.*;
import gwas.*;
import db.*;


public class Launch {
	public static final String[] LAUNCH_TYPES = {"hits", "dummy", "counts", "miss", "indep", "genes", "plink", "simpleM", "score", "parse", "ucsc", "split", "cat", "db", "merge", "mergeSNPs", "trimFam", "freq - computes weighted allele frequency", "uniform - creates a hits control file where each file listed has the same column names, only with a different prefix", "metal", "transpose", "forest"};

	public static void run(String filename, Logger log) throws Elision { // dir is necessary for UpdateCRFdb
		BufferedReader reader;
		String temp;

		try {
			reader = Files.getReader(filename, false, true, log, true);
			temp = reader.readLine();
			
			log.report("Launching option '"+temp+"'");

			if (temp.equals("hits")) {
				LookupTable.fromParameters(filename, log);
			} else if (temp.equals("dummy")) {
				DummyDataset.createFromParameters(filename, log);
			} else if (temp.equals("counts")) {
				DummyDataset.createReverseFromParameters(filename, log);
			} else if (temp.equals("miss")) {
				MarkerQC.parseParameters(filename, log);
			} else if (temp.equals("indep")) {
				IndependentSNPs.selectFromParameters(filename, log);
			} else if (temp.equals("genes")) {
				MapGenesToSNPs.filter(filename, log);
			} else if (temp.equals("plink")) {
				CreateDatabaseFromPlink.createDatabase(filename, log);
			} else if (temp.equals("simpleM")) {
				CreateDatabaseFromPlink.createCountsMatrix(filename, log);
			} else if (temp.equals("score")) {
				CreateDatabaseFromPlink.createScore(filename, log);
			} else if (temp.equals("parse")) {
				GenParser.parseFile(filename, log);
			} else if (temp.equals("ucsc")) {
				UCSCtrack.describeFromParameters(filename, log);
			} else if (temp.equals("split")) {
				Files.splitFileFromParamters(filename, log);
			} else if (temp.equals("cat")) {
				Files.catFilesFromParamters(filename, log);
			} else if (temp.equals("db")) {
				CreateMarkerDatabase.createFromParameters(filename, log);
			} else if (temp.equals("merge")) {
				Files.mergeFromParameters(filename, log);
			} else if (temp.equals("mergeSNPs")) {
				Files.mergeSNPlistsFromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("trimFam")) {
				TrimFam.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("freq")) {
				Metal.calculateWeightedAlleleFrequencyFromParamters(filename, log);
			} else if (temp.equalsIgnoreCase("uniform")) {
				Metal.generateUniformsFromParamters(filename, log);
			} else if (temp.equalsIgnoreCase("metal")) {
				Metal.generateInputFile(filename, log);
			} else if (temp.equalsIgnoreCase("transpose")) {
				Transformations.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("forest")) {
				ForestPlot.fromParameters(filename, log);
			} else {
				log.reportError("Error - '"+temp+"' is an invalid launch type, options include:");
				log.reportError(Array.toStr(LAUNCH_TYPES, "\n"));
			}
		} catch (Exception e) {
			log.reportError(e.getMessage());
			log.reportException(e);
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
//		String filename = "trim_withScores.crf";
		String filename = "singleton.crf";
		boolean suppress = false;

		String usage = "\n"+
		"park.crfDB requires 0-1 arguments\n"+
		"   (1) filename (i.e. "+filename+" (default)\n"+
		"   (2) suppress Windows stdin (i.e. -suppress (not the default)\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("-suppress")) {
				suppress = true;
				numArgs--;
			} else {
				filename = args[i];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			run(filename, new Logger(ext.rootOf(filename)+".log"));
		} catch (Exception e) {
			e.printStackTrace();
		}
		if (!suppress && System.getProperty("os.name").startsWith("Windows")) {
			System.err.println();
			System.err.println("Press any key to close this window");
			new BufferedReader(new InputStreamReader(System.in)).readLine();
		}
	}
}
