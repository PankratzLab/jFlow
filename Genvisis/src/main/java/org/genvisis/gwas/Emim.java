package org.genvisis.gwas;

import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.genvisis.common.*;

public class Emim {
	
	public enum EMIM_MODEL {
		GENOTYPIC("Genotypic", 2),
		DOMINANT("Dominant", 1),
		ADDITIVE("Additive", 1);
		
		private String name;
		private int degreesFreedom;
		
		EMIM_MODEL(String name, int degreesFreedom){
			this.name = name;
			this.degreesFreedom = degreesFreedom;
		}
		
		public String toString() {
			return name;
		}
		
		public int getDegreesOfFreedom() {
			return degreesFreedom;
		}
		
		public static Set<EMIM_MODEL> valueSet() {
			HashSet<EMIM_MODEL> models = new HashSet<EMIM_MODEL>();
			for (EMIM_MODEL model : EMIM_MODEL.values()) models.add(model);
			return models;
		}
	}
	
	private enum EMIM_PARAM {
		ESTIMATE_R1("   << estimate R1 (0=no, 1=yes)"),
		ESTIMATE_R2("   << estimate R2 (0=no, 1=yes)"),
		R2_EQUALS_R1("   << R2=R1 (0=no, 1=yes)"),
		R2_EQUALS_R1_SQUARED("   << R2=R1squared\t(0=no, 1=yes)"),
		ESTIMATE_S1("   << estimate S1 (0=no, 1=yes)"),
		ESTIMATE_S2("   << estimate S2 (0=no, 1=yes)"),
		S2_EQUALS_S1("   << S2=S1 (0=no, 1=yes)"),
		S2_EQUALS_S1_SQUARED("   << S2=S1squared\t(0=no, 1=yes)");
		
		private String lineSuffix;
		
		EMIM_PARAM(String lineSuffix){
			this.lineSuffix = lineSuffix;
		}
		
		String[] getReplacement(boolean setTo) {
			return setTo ? new String[] {"0" + lineSuffix, "1" + lineSuffix} : new String[] {"1" + lineSuffix, "0" + lineSuffix};
		}
	}
	
	private static void generateEmimParams(String filenameOriginal, HashSet<EMIM_PARAM> setParams, Logger log){
		String[][] replacements = new String[EMIM_PARAM.values().length][];
		for (int i = 0; i < EMIM_PARAM.values().length; i++) {
			EMIM_PARAM param = EMIM_PARAM.values()[i];
			replacements[i] = param.getReplacement(setParams.contains(param));
		}
		
		replaceLines(filenameOriginal, "emimparams.dat", replacements, log);
	}
	
	private static void setTo(String runType, EMIM_MODEL model) {
		String filenameOriginal;
		Logger log;
		
		if (!Files.exists("emimparams.dat")) {
			System.err.println("Error - emimparams.dat is not in the current directory; aborting rewrite");
			return;
		}

		log = new Logger();
		filenameOriginal = Files.backup("emimparams.dat", "./", "./", true);
		
		HashSet<EMIM_PARAM> setParams = new HashSet<EMIM_PARAM>();
		
		if (runType.equals("C")) {
			switch (model) {
			case GENOTYPIC:
				setParams.add(EMIM_PARAM.ESTIMATE_R1);
				setParams.add(EMIM_PARAM.ESTIMATE_R2);
				break;
			case DOMINANT: 
				setParams.add(EMIM_PARAM.R2_EQUALS_R1);
				break;
			case ADDITIVE:
				setParams.add(EMIM_PARAM.R2_EQUALS_R1_SQUARED);
				break;
			default:
				log.reportError("Invalid EMIM model: " + model.toString());
				return;
			}
			
		} else if (runType.equals("CM")) {
			switch (model) {
			case GENOTYPIC:
				setParams.add(EMIM_PARAM.ESTIMATE_R1);
				setParams.add(EMIM_PARAM.ESTIMATE_R2);
				setParams.add(EMIM_PARAM.ESTIMATE_S1);
				setParams.add(EMIM_PARAM.ESTIMATE_S2);
				break;
			case DOMINANT: 
				setParams.add(EMIM_PARAM.R2_EQUALS_R1);
				setParams.add(EMIM_PARAM.S2_EQUALS_S1);
				break;
			case ADDITIVE:
				setParams.add(EMIM_PARAM.R2_EQUALS_R1_SQUARED);
				setParams.add(EMIM_PARAM.S2_EQUALS_S1_SQUARED);
				break;
			default:
				log.reportError("Invalid EMIM model: " + model.toString());
				return;
			}
		} else if (runType.equals("M")) {
			switch (model) {
			case GENOTYPIC:
				setParams.add(EMIM_PARAM.ESTIMATE_S1);
				setParams.add(EMIM_PARAM.ESTIMATE_S2);
				break;
			case DOMINANT: 
				setParams.add(EMIM_PARAM.S2_EQUALS_S1);
				break;
			case ADDITIVE:
				setParams.add(EMIM_PARAM.S2_EQUALS_S1_SQUARED);
				break;
			default:
				log.reportError("Invalid EMIM model: " + model.toString());
				return;
			}
		} else {
			log.reportError("Invalid Run Type: " + runType);
			return;
		}
		
		generateEmimParams(filenameOriginal, setParams, log);
		
		System.out.println("emimparams.dat set to " + model.toString() + " model of " + runType + " effect");
		
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

	protected static boolean scriptAllInDir(String runDir, String plinkDirAndRoot, String relativePlinkRoot, String excludeFile, String keepFile, double pThreshold, Set<EMIM_MODEL> models, String resultPrefix, Logger log) {
	    String commands;
        String currDir = ext.verifyDirFormat(runDir);
        boolean forceRun = false, forceParse = false;
        
        if (currDir.charAt(0) != '/' && !currDir.contains(":")) {
            currDir = (new File("./" + currDir)).getAbsolutePath() + "/";
        }
        
        if (excludeFile.equals("GEN")) {
            excludeFile = "sexChrMarkers.txt";
            if (Files.exists(currDir + excludeFile)) {
            	log.report(currDir + excludeFile + " already exists, skipping exclude file generation");
            }
            else {
            	try {
            		listSexMarkers(plinkDirAndRoot + ".bim", currDir + excludeFile);
                	forceRun = true;
                	forceParse = true;
            	} catch (Exception e) {
            		excludeFile = null;
            		e.printStackTrace();
            	}
            }
        }
        
        commands = "";
        
        if (!forceRun && Files.exists(currDir + "plink_prep.log") && Files.exists(currDir + "emimPrep.bed") && Files.exists(currDir + "emimPrep.bim") && Files.exists(currDir + "emimPrep.fam")) {
        	log.report(currDir + "emimPrep PLINK files already exist, skipping PLINK file generation");
        } else {
        	forceRun = true;
        	forceParse = true;
        	commands += "plink2 --noweb --bfile " + relativePlinkRoot
                        	    + (excludeFile != null ? " --exclude " + excludeFile : "")  
                        	    + (keepFile != null ? " --keep " + keepFile : "")  
                        	    + " --make-bed --out emimPrep\n"+
                        "mv emimPrep.log plink_prep.log\n"+
                        "\n";
        }
        
        if (!forceRun && Files.exists(currDir + "plink_mendel.log") && Files.exists(currDir + "plink.mendel") && Files.exists(currDir + "plink.lmendel") && Files.exists(currDir + "plink.fmendel") && Files.exists(currDir + "plink.imendel")) {
        	log.report(currDir + "plink.*mendel files already exist, skipping mendelian error calculation");
        } else {
        	forceParse = true;
        	commands += "plink2 --noweb --bfile emimPrep --mendel\n"+
                    	"mv plink.log plink_mendel.log\n"+
                    	"\n";
        }
        
        if (!forceRun && Files.exists(currDir + "plink_tdt.log") && Files.exists(currDir + "plink.tdt") ){
        	log.report(currDir + "plink.tdt already exists, skipping TDT");
        } else {
        	forceParse = true;
        	commands += "plink2 --noweb --bfile emimPrep --tdt --ci 0.95\n"+
                    	"mv plink.log plink_tdt.log\n"+
                    	"\n";
        }
        
        if (!forceRun && Files.exists(currDir + "plink_hwe.log") && Files.exists(currDir + "plink.hwe")) {
        	log.report(currDir + "plink.hwe already exists, skipping Hardy-Weinberg Equilibrium calculation");
        } else {
        	forceParse = true;
        	commands += "plink2 --noweb --bfile emimPrep --hardy\n"+
        				"mv plink.log plink_hwe.log\n"+
        				"\n";
        }
        
        if (!forceRun && Files.exists(currDir + "plink_freq.log") && Files.exists(currDir + "plink.frq")) {
        	log.report(currDir + "plink.frq already exists, skipping Minor Allele Frequency calculation");
        } else {
        	forceParse = true;
        	commands += "plink2 --noweb --bfile emimPrep --freq\n"+
        				"mv plink.log plink_freq.log\n"+
        				"\n";
        }

        if (!forceRun &&
        		Files.exists(currDir + "premim.log") && 
        		Files.exists(currDir + "risksnplist.txt") &&
        		Files.exists(currDir + "emimparams.dat") &&
        		Files.exists(currDir + "emimmarkers.dat") &&
        		Files.exists(currDir + "caseparenttrios.dat") &&
        		Files.exists(currDir + "caseparents.dat") &&
        		Files.exists(currDir + "casemotherduos.dat") &&
        		Files.exists(currDir + "casefatherduos.dat") &&
        		Files.exists(currDir + "casemothers.dat") &&
        		Files.exists(currDir + "casefathers.dat") &&
        		Files.exists(currDir + "cases.dat") &&
        		Files.exists(currDir + "conparents.dat") &&
        		Files.exists(currDir + "conmotherduos.dat") &&
        		Files.exists(currDir + "confatherduos.dat") &&
        		Files.exists(currDir + "cons.dat")) {
        	log.report("Outputs of PREMIM in " + currDir + " already exist, skipping PREMIM");
        } else {
        	forceRun = true;
        	forceParse = true;
        	commands += "premim -cg -a -rout risksnplist.txt emimPrep.bed\n";
        }
                
        
        for (EMIM_MODEL model : models){
        	boolean skipModel = !forceRun &&
        			Files.exists(currDir + (resultPrefix == null ? "" : resultPrefix + "_") + "results_pVals_" + model.toString() + ".xln") &&
        			Files.exists(currDir + "emimsummary_C_" + model.toString() + ".out") &&
        			Files.exists(currDir + "emimparams_C_" + model.toString() + ".dat") &&
        			Files.exists(currDir + "emimsummary_CM_" + model.toString() + ".out") &&
        			Files.exists(currDir + "emimparams_CM_" + model.toString() + ".dat") &&
        			Files.exists(currDir + "emimsummary_M_" + model.toString() + ".out") &&
        			Files.exists(currDir + "emimparams_M_" + model.toString() + ".dat");
        	if (skipModel) {
        		log.report("Results already exist in " + currDir + " for " + model.toString() + " model, skipping " + model.toString() + " EMIM");
        	} else {
	        	commands += "\n"+
	        	            Files.getRunString() + " gwas.Emim run=C model=" + model.toString() + "\n"+
	                    	"emim\n"+
	                    	"mv emimsummary.out emimsummary_C_" + model.toString() + ".out\n"+
	                    	"rm emimresults.out\n"+
	                    	"cp emimparams.dat emimparams_C_" + model.toString() + ".dat\n"+
	                    	"\n"+
	                    	Files.getRunString() + " gwas.Emim run=CM model=" + model.toString() + "\n"+
	                    	"emim\n"+
	                    	"mv emimsummary.out emimsummary_CM_" + model.toString() + ".out\n"+
	                    	"rm emimresults.out\n"+
	                    	"cp emimparams.dat emimparams_CM_" + model.toString() + ".dat\n"+
	                    	"\n"+
	                    	Files.getRunString() + " gwas.Emim run=M model=" + model.toString() + "\n"+
	                    	"emim\n"+
	                    	"mv emimsummary.out emimsummary_M_" + model.toString() + ".out\n"+
	                    	"rm emimresults.out\n"+
	                    	"cp emimparams.dat emimparams_M_" + model.toString() + ".dat\n"+
	                    	"\n";
        	}
        	
        	if (!skipModel || forceParse) {
        		commands += Files.getRunString() + " gwas.Emim parse=./" 
	                      + " hwe=plink.hwe" 
	                      + " frq=plink.frq"
	                      + " pThreshold=" + pThreshold 
	                      + " model=" + model.toString() 
	                      + (resultPrefix == null ? "" : " resultPrefix=" + resultPrefix) +
	                      "\n\n";
        	}
        }
        
        if (commands.equals("")) {
        	log.report("Warning - No commands were generated for " + currDir + ", remove directory to re-run EMIM Pipeline from scratch");
        	return false;
        }
        commands = "cd " + currDir + "\n" + commands;
        
        Files.qsub(currDir + ext.rootOf(plinkDirAndRoot, true)+"_runEmim.pbs", commands, 8000, 6, 1);
        
        return true;
	}
	
	public static boolean scriptAll(String plinkPrefix, String excludeFile, String keepFile, double pThreshold) {
	    return scriptAllInDir("./", plinkPrefix, plinkPrefix, excludeFile, keepFile, pThreshold, EMIM_MODEL.valueSet(), null, new Logger());
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
	
	public static void parse(String dir, String resultPrefix, String hweFile, String frqFile, double pValueThreshold, EMIM_MODEL model) {
		String resultsFileChild, resultsFileChildMom, resultsFileMom, resultsFileTdt, mapFile, mendelErrorFile, outfile;

		resultsFileChild = dir+"emimsummary_C_" + model.toString() + ".out";
		resultsFileMom = dir+"emimsummary_M_" + model.toString() + ".out";
		resultsFileChildMom = dir+"emimsummary_CM_" + model.toString() + ".out";
		resultsFileTdt = dir+"plink.tdt";
		mapFile = dir+"emimPrep.bim";
		mendelErrorFile = dir+"plink.lmendel";
		outfile = dir + (resultPrefix == null ? "" : resultPrefix + "_") + "results_pVals_" + model.toString() + ".xln";
		
		ResultsPackager.parseEmimFormat(resultsFileChild,
										resultsFileMom,
										resultsFileChildMom, 
										resultsFileTdt, 
										mapFile, 
										mendelErrorFile, 
										hweFile,
										frqFile,
										pValueThreshold, 
										outfile,
										new Logger("EMIMparser.log"));
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
		double pValueThreshold = 1.1;
		String hweFile = null;
		String frqFile = null;
		String plinkPrefix = null;
		EMIM_MODEL model = EMIM_MODEL.DOMINANT;
		String excludeFile = "GEN";
		String keepFile = null;
		String resultPrefix = null;

		String usage = "\n" +
		"gwas.Emim requires 0-1 arguments\n" +
		"   (1) run type (either C, CM, or M) (i.e. run=" + runType + " (default))\n" +
		"   (2) model (" + Array.toStr(EMIM_MODEL.values(), ",") + ") (i.e. model=" + model.toString() + " (default))\n" +
		"  OR\n" +
		"   (1) generate script that runs the full process (i.e. script=plinkPrefix (not the default))\n" + 
		"   (2) p-value threshold to filter on (piped to parse method) (i.e. pThreshold=" + pValueThreshold + " (default))\n" + 
		"   (3) (optional) a file of markers to exclude - the default will generate a file of sex markers, as Emim won't parse these. (i.e. exclude=drops.dat (not the default))\n" + 
		"   (4) (optional) a file of tab-delimited FID/IID pairs to keep (i.e. keep=completeTrios.dat (not the default))\n" + 
		"  OR\n" +
		"   (1) directory to parse (i.e. parse=./ (not the default))\n" + 
		"   (2) file prefix for results (i.e. resultPrefix=" + resultPrefix + " (default))\n" + 
		"   (3) p-value threshold to filter on (i.e. pThreshold=" + pValueThreshold + " (default))\n" + 
		"   (4) model " + Arrays.toString(EMIM_MODEL.values()) + " (i.e. model=" + model.toString() + " (default))\n" +
		"   (5) (optional) plink.hwe file to merge with results (i.e. hwe=" + hweFile + " (default))\n" + 
		"   (6) (optional) plink.frq file to merge with results (i.e. frq=" + frqFile + " (default))\n" + 
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
			} else if (args[i].startsWith("frq=")) {
				frqFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("resultPrefix=")) {
				resultPrefix = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("model=")) {
				String modelString = ext.parseStringArg(args[i], null);
				model = null;
				for (EMIM_MODEL m : EMIM_MODEL.values()){
					if (modelString.equals(m.toString())) {
						model = m;
						break;
					}
				}
				if (model == null) System.err.println("Invalid model: " + modelString);
				else numArgs--;
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
				parse(dir, resultPrefix, hweFile, frqFile, pValueThreshold, model);
			} else if (plinkPrefix != null) {
				scriptAll(plinkPrefix, excludeFile, keepFile, pValueThreshold);
			} else {
				setTo(runType, model);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
