package gwas;

import java.io.*;

import common.*;

public class Emim {
	private static void setTo(String runType) {
		String filenameOriginal;
		Logger log;
		
		if (!Files.exists("emimparams.dat")) {
			System.err.println("Error - emimparams.dat is not in the current directory; aborting rewrite");
			return;
		}

		log = new Logger();
		filenameOriginal = Files.backup("emimparams.dat", "./", "./", true);
		
		if (runType.equals("C")) {
			replaceLines(filenameOriginal, "emimparams.dat", new String[][] {
					{"0   << estimate R1 (0=no, 1=yes)", "1   << estimate R1 (0=no, 1=yes)"}, 
					{"0   << estimate R2 (0=no, 1=yes)", "1   << estimate R2 (0=no, 1=yes)"},
					{"1   << estimate S1 (0=no, 1=yes)", "0   << estimate S1 (0=no, 1=yes)"}, 
					{"1   << estimate S2 (0=no, 1=yes)", "0   << estimate S2 (0=no, 1=yes)"}
				}, log);
		}
		
		if (runType.equals("CM")) {
			replaceLines(filenameOriginal, "emimparams.dat", new String[][] {
					{"0   << estimate R1 (0=no, 1=yes)", "1   << estimate R1 (0=no, 1=yes)"}, 
					{"0   << estimate R2 (0=no, 1=yes)", "1   << estimate R2 (0=no, 1=yes)"},
					{"0   << estimate S1 (0=no, 1=yes)", "1   << estimate S1 (0=no, 1=yes)"}, 
					{"0   << estimate S2 (0=no, 1=yes)", "1   << estimate S2 (0=no, 1=yes)"}
				}, log);
		}

		if (runType.equals("M")) {
			replaceLines(filenameOriginal, "emimparams.dat", new String[][] {
					{"1   << estimate R1 (0=no, 1=yes)", "0   << estimate R1 (0=no, 1=yes)"}, 
					{"1   << estimate R2 (0=no, 1=yes)", "0   << estimate R2 (0=no, 1=yes)"},
					{"0   << estimate S1 (0=no, 1=yes)", "1   << estimate S1 (0=no, 1=yes)"}, 
					{"0   << estimate S2 (0=no, 1=yes)", "1   << estimate S2 (0=no, 1=yes)"}
				}, log);
		}
	}

	private static void scriptAll(String plinkPrefix) {
		String commands;
		
		// TODO generate sexChrMarkers.txt		
		
		commands = "plink --noweb --bfile "+plinkPrefix+" --exclude sexChrMarkers.txt --make-bed --out emimPrep\n"+
				"plink --noweb --bfile emimPrep --mendel\n"+
				"plink --noweb --bfile emimPrep --tdt\n"+
				"premim -cg -a -rout risksnplist.txt emimPrep.bed\n"+
				"\n"+
				"jcp gwas.Emim run=C\n"+
				"emim\n"+
				"mv emimsummary.out emimsummary_C.out\n"+
				"mv emimresults.out emimresults_C.out\n"+
				"cp emimparams.dat emimparams_C.dat\n"+
				"\n"+
				"jcp gwas.Emim run=CM\n"+
				"emim\n"+
				"mv emimsummary.out emimsummary_CM.out\n"+
				"mv emimresults.out emimresults_CM.out\n"+
				"cp emimparams.dat emimparams_CM.dat\n"+
				"\n"+
				"jcp gwas.Emim run=M\n"+
				"emim\n"+
				"mv emimsummary.out emimsummary_M.out\n"+
				"mv emimresults.out emimresults_M.out\n"+
				"cp emimparams.dat emimparams_M.dat\n"+
				"# run jcp gwas.ResultsPackager"+"\n";
		
		Files.qsub(plinkPrefix+"_runEmim.pbs", commands, 5000, 24, 1);
	}
	
	public static void parse(String dir, String hweFile, double pValueThreshold) {
		String resultsFileChild, resultsFileChildMom, resultsFileMom, resultsFileTdt, mapFile, mendelErrorFile, outfile;

		resultsFileChild = dir+"emimsummary_C.out";
		resultsFileMom = dir+"emimsummary_M.out";
		resultsFileChildMom = dir+"emimsummary_CM.out";
		resultsFileTdt = dir+"plink.tdt";
		mapFile = dir+"plink.bim";
		mendelErrorFile = dir+"plink.lmendel";
		outfile = dir+"results_pVals.xln";

		ResultsPackager.parseEmimFormat(resultsFileChild, resultsFileMom, resultsFileChildMom, resultsFileTdt, mapFile, mendelErrorFile, hweFile, pValueThreshold, outfile, new Logger("EMIMparser.log"));
	}
	
	public static void replaceLines(String filenameOriginal, String filenameWithReplacements, String[][] relacements, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		
		try {
			reader = Files.getAppropriateReader(filenameOriginal);
			writer = new PrintWriter(new FileWriter(filenameWithReplacements));
			while (reader.ready()) {
				writer.println(ext.replaceAllWith(reader.readLine(), relacements));
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filenameOriginal + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filenameOriginal + "\"");
			return;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String runType = "C";
		String dir = "./";
		double pValueThreshold = 0.0001;
		String hweFile = null;
		String plinkPrefix = null;

		String usage = "\n" +
		"gwas.Emim requires 0-1 arguments\n" +
		"   (1) run type (either C, CM, or M) (i.e. run=" + runType + " (default))\n" +
		"  OR\n" +
		"   (1) generate script that runs the full process (i.e. script=plinkPrefix (not the default))\n" + 
		"  OR\n" +
		"   (1) directory to parse (i.e. parse=" + dir + " (default))\n" + 
		"   (2) p-value threshold to filter on (i.e. pThreshold=" + pValueThreshold + " (default))\n" + 
		"   (3) (optional) plink.hwe file to merge with results (i.e. hwe=" + hweFile + " (default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("run=")) {
				runType = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("plinkPrefix=")) {
				plinkPrefix = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("parse=")) {
				dir = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("pThreshold=")) {
				pValueThreshold = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("hwe=")) {
				hweFile = ext.parseStringArg(args[i], null);
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
			if (dir != null) {
				parse(dir, hweFile, pValueThreshold);
			} else if (plinkPrefix != null) {
				scriptAll(plinkPrefix);
			} else {
				setTo(runType);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
