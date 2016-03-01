package gwas;

import java.io.*;

import common.*;

public class Emim {
	private static void setTo(String runType, boolean allelic) {
		String filenameOriginal;
		Logger log;
		
		if (!Files.exists("emimparams.dat")) {
			System.err.println("Error - emimparams.dat is not in the current directory; aborting rewrite");
			return;
		}

		log = new Logger();
		filenameOriginal = Files.backup("emimparams.dat", "./", "./", true);
		
		if (runType.equals("C")) {
			if (allelic) {
				replaceLines(filenameOriginal, "emimparams.dat", new String[][] {
						{"0   << estimate R1 (0=no, 1=yes)", "1   << estimate R1 (0=no, 1=yes)"},
						{"0   << estimate R2 (0=no, 1=yes)", "1   << estimate R2 (0=no, 1=yes)"},
						{"0   << R2=R1 (0=no, 1=yes)", "1   << R2=R1 (0=no, 1=yes)"},
						{"1   << estimate S1 (0=no, 1=yes)", "0   << estimate S1 (0=no, 1=yes)"},
						{"1   << estimate S2 (0=no, 1=yes)", "0   << estimate S2 (0=no, 1=yes)"},
						{"0   << S2=S1 (0=no, 1=yes)", "1   << S2=S1 (0=no, 1=yes)"}
	
					}, log);
			} else {
				replaceLines(filenameOriginal, "emimparams.dat", new String[][] {
						{"0   << estimate R1 (0=no, 1=yes)", "1   << estimate R1 (0=no, 1=yes)"},
						{"0   << estimate R2 (0=no, 1=yes)", "1   << estimate R2 (0=no, 1=yes)"},
						{"1   << R2=R1 (0=no, 1=yes)", "0   << R2=R1 (0=no, 1=yes)"},
						{"1   << estimate S1 (0=no, 1=yes)", "0   << estimate S1 (0=no, 1=yes)"},
						{"1   << estimate S2 (0=no, 1=yes)", "0   << estimate S2 (0=no, 1=yes)"},
						{"1   << S2=S1 (0=no, 1=yes)", "0   << S2=S1 (0=no, 1=yes)"}
	
					}, log);
			}
			
		}
		
		if (runType.equals("CM")) {
			if (allelic) {
				replaceLines(filenameOriginal, "emimparams.dat", new String[][] {
						{"0   << estimate R1 (0=no, 1=yes)", "1   << estimate R1 (0=no, 1=yes)"}, 
						{"0   << estimate R2 (0=no, 1=yes)", "1   << estimate R2 (0=no, 1=yes)"},
						{"0   << R2=R1 (0=no, 1=yes)", "1   << R2=R1 (0=no, 1=yes)"},
						{"0   << estimate S1 (0=no, 1=yes)", "1   << estimate S1 (0=no, 1=yes)"}, 
						{"0   << estimate S2 (0=no, 1=yes)", "1   << estimate S2 (0=no, 1=yes)"},
						{"0   << S2=S1 (0=no, 1=yes)", "1   << S2=S1 (0=no, 1=yes)"}
					}, log);
			} else {
				replaceLines(filenameOriginal, "emimparams.dat", new String[][] {
						{"0   << estimate R1 (0=no, 1=yes)", "1   << estimate R1 (0=no, 1=yes)"}, 
						{"0   << estimate R2 (0=no, 1=yes)", "1   << estimate R2 (0=no, 1=yes)"},
						{"1   << R2=R1 (0=no, 1=yes)", "0   << R2=R1 (0=no, 1=yes)"},
						{"0   << estimate S1 (0=no, 1=yes)", "1   << estimate S1 (0=no, 1=yes)"}, 
						{"0   << estimate S2 (0=no, 1=yes)", "1   << estimate S2 (0=no, 1=yes)"},
						{"1   << S2=S1 (0=no, 1=yes)", "0   << S2=S1 (0=no, 1=yes)"}
					}, log);
			}
		}

		if (runType.equals("M")) {
			if (allelic) {
				replaceLines(filenameOriginal, "emimparams.dat", new String[][] {
						{"1   << estimate R1 (0=no, 1=yes)", "0   << estimate R1 (0=no, 1=yes)"}, 
						{"1   << estimate R2 (0=no, 1=yes)", "0   << estimate R2 (0=no, 1=yes)"},
						{"0   << R2=R1 (0=no, 1=yes)", "1   << R2=R1 (0=no, 1=yes)"},
						{"0   << estimate S1 (0=no, 1=yes)", "1   << estimate S1 (0=no, 1=yes)"}, 
						{"0   << estimate S2 (0=no, 1=yes)", "1   << estimate S2 (0=no, 1=yes)"},
						{"0   << S2=S1 (0=no, 1=yes)", "1   << S2=S1 (0=no, 1=yes)"}
					}, log);
			} else {
				replaceLines(filenameOriginal, "emimparams.dat", new String[][] {
						{"1   << estimate R1 (0=no, 1=yes)", "0   << estimate R1 (0=no, 1=yes)"}, 
						{"1   << estimate R2 (0=no, 1=yes)", "0   << estimate R2 (0=no, 1=yes)"},
						{"1   << R2=R1 (0=no, 1=yes)", "0   << R2=R1 (0=no, 1=yes)"},
						{"0   << estimate S1 (0=no, 1=yes)", "1   << estimate S1 (0=no, 1=yes)"}, 
						{"0   << estimate S2 (0=no, 1=yes)", "1   << estimate S2 (0=no, 1=yes)"},
						{"1   << S2=S1 (0=no, 1=yes)", "0   << S2=S1 (0=no, 1=yes)"}
					}, log);
			}
		}
	}
	
	private static void listSexMarkers(String bimFile, String sexFile) throws NumberFormatException, IOException {
	    PrintWriter writer;
	    BufferedReader reader;
	    String line;
	    String[] parts;
	    int chr;
	    
	    writer = Files.getAppropriateWriter(sexFile);
	    reader = Files.getAppropriateReader(bimFile);
	    line = null;
	    while ((line = reader.readLine()) != null) {
	        parts = line.trim().split("\t", -1);
	        chr = Integer.parseInt(parts[0]);
	        if (chr > 22) {
	            writer.println(parts[1]);
	        }
	    }
	    writer.flush();
	    writer.close();
	    reader.close();
	}

	protected static void scriptAllInDir(String runDir, String plinkDirAndRoot, String relativePlinkRoot, String excludeFile, String keepFile, double pThreshold, String resultPrefix) {
	    String commands;
        String currDir = ext.verifyDirFormat(runDir);
        
        if (currDir.charAt(0) != '/' && !currDir.contains(":")) {
            currDir = (new File("./" + currDir)).getAbsolutePath() + "/";
        }
        
        if (excludeFile.equals("GEN")) {
            try {
                listSexMarkers(plinkDirAndRoot + ".bim", currDir + "sexChrMarkers.txt");
                excludeFile = "sexChrMarkers.txt";
            } catch (Exception e) {
                excludeFile = null;
                e.printStackTrace();
            }
        }
        
        commands = "cd " + currDir + "\n" + 
                "plink2 --noweb --bfile " 
                        + relativePlinkRoot
                        + (excludeFile != null ? " --exclude " + excludeFile : "")  
                        + (keepFile != null ? " --keep " + keepFile : "")  
                        + " --make-bed --out emimPrep\n"+
                "mv plink.log plink_prep.log\n"+
                "\n"+
                "plink2 --noweb --bfile emimPrep --mendel\n"+
                "mv plink.log plink_mendel.log\n"+
                "\n"+
                "plink2 --noweb --bfile emimPrep --tdt\n"+
                "mv plink.log plink_tdt.log\n"+
                "\n"+
                "plink2 --noweb --bfile emimPrep --hardy\n"+
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
                "\n"+
                "jcp gwas.Emim parse=./ hwe=plink.hwe pThreshold=" + pThreshold + 
                "\n";
            if (resultPrefix != null) {
                commands += "if [ -f results_pVals.xln ] ; then \n" +
                            "    mv results_pVals.xln " + resultPrefix + "_results_pVals.xln\n" + 
                            "fi;\n" + 
                            "";
            }
        
        
        Files.qsub(currDir + ext.rootOf(plinkDirAndRoot, true)+"_runEmim.pbs", commands, 5000, 24, 1);
	}
	
	public static void scriptAll(String plinkPrefix, String excludeFile, String keepFile, double pThreshold) {
	    scriptAllInDir("./", plinkPrefix, plinkPrefix, excludeFile, keepFile, pThreshold, null);
//		String commands;
//		String currDir;
//		
//		if (excludeFile.equals("GEN")) {
//		    try {
//                listSexMarkers(plinkPrefix + ".bim", "sexChrMarkers.txt");
//                excludeFile = "sexChrMarkers.txt";
//            } catch (NumberFormatException | IOException e) {
//                excludeFile = null;
//                e.printStackTrace();
//            }
//		}
//		
//		currDir = (new File("./")).getAbsolutePath();
//		commands = "cd " + currDir + "\n" + 
//		        "plink2 --noweb --bfile " 
//		                + plinkPrefix 
//		                + (excludeFile != null ? " --exclude " + excludeFile : "")  
//		                + (keepFile != null ? " --keep " + keepFile : "")  
//		                + " --make-bed --out emimPrep\n"+
//				"plink2 --noweb --bfile emimPrep --mendel\n"+
//				"plink2 --noweb --bfile emimPrep --tdt\n"+
//				"plink2 --noweb --bfile emimPrep --hardy\n"+
//				"premim -cg -a -rout risksnplist.txt emimPrep.bed\n"+
//				"\n"+
//				"jcp gwas.Emim run=C\n"+
//				"emim\n"+
//				"mv emimsummary.out emimsummary_C.out\n"+
//				"mv emimresults.out emimresults_C.out\n"+
//				"cp emimparams.dat emimparams_C.dat\n"+
//				"\n"+
//				"jcp gwas.Emim run=CM\n"+
//				"emim\n"+
//				"mv emimsummary.out emimsummary_CM.out\n"+
//				"mv emimresults.out emimresults_CM.out\n"+
//				"cp emimparams.dat emimparams_CM.dat\n"+
//				"\n"+
//				"jcp gwas.Emim run=M\n"+
//				"emim\n"+
//				"mv emimsummary.out emimsummary_M.out\n"+
//				"mv emimresults.out emimresults_M.out\n"+
//				"cp emimparams.dat emimparams_M.dat\n"+
//				"\n"+
//				"jcp gwas.Emim parse=./ hwe=plink.hwe pThreshold=" + pThreshold;
//		
//		Files.qsub(plinkPrefix+"_runEmim.pbs", commands, 5000, 24, 1);
	}
	
	public static void parse(String dir, String hweFile, double pValueThreshold) {
		String resultsFileChild, resultsFileChildMom, resultsFileMom, resultsFileTdt, mapFile, mendelErrorFile, outfile;

		resultsFileChild = dir+"emimsummary_C.out";
		resultsFileMom = dir+"emimsummary_M.out";
		resultsFileChildMom = dir+"emimsummary_CM.out";
		resultsFileTdt = dir+"plink.tdt";
		mapFile = dir+"emimPrep.bim";
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
		String dir = null;
		double pValueThreshold = 0.0001;
		String hweFile = null;
		String plinkPrefix = null;
		boolean allelic = true;
		String excludeFile = "GEN";
		String keepFile = null;

		String usage = "\n" +
		"gwas.Emim requires 0-1 arguments\n" +
		"   (1) run type (either C, CM, or M) (i.e. run=" + runType + " (default))\n" +
		"   (2) allelic instead of 2df genotypic test (i.e. allelic=" + allelic + " (default))\n" +
		"  OR\n" +
		"   (1) generate script that runs the full process (i.e. script=plinkPrefix (not the default))\n" + 
		"   (2) p-value threshold to filter on (piped to parse method) (i.e. pThreshold=" + pValueThreshold + " (default))\n" + 
		"   (2) (optional) a file of markers to exclude - the default will generate a file of sex markers, as Emim won't parse these. (i.e. exclude=drops.dat (not the default))\n" + 
		"   (3) (optional) a file of tab-delimited FID/IID pairs to keep (i.e. keep=completeTrios.dat (not the default))\n" + 
		"  OR\n" +
		"   (1) directory to parse (i.e. parse=./ (not the default))\n" + 
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
			} else if (args[i].startsWith("exclude=")) {
			    excludeFile = ext.parseStringArg(args[i], null);
			    numArgs--;
			} else if (args[i].startsWith("keep=")) {
			    keepFile = ext.parseStringArg(args[i], null);
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
			} else if (args[i].startsWith("allelic=")) {
				allelic = ext.parseBooleanArg(args[i]);
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
				scriptAll(plinkPrefix, excludeFile, keepFile, pValueThreshold);
			} else {
				setTo(runType, allelic);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
