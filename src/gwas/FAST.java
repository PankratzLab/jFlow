package gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.math.BigInteger;
import java.net.URLDecoder;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import parse.GenParser;
import common.Aliases;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

public class FAST {
	
	private static final String COUNT_SYMB = "<>";
	private static final String CHARGE_FORMAT = " 'SNP.id'=Markername 'Chr'=Chr 'Pos'=Pos $#" + COUNT_SYMB + "=N 'Coded.Allele'=Effect_allele 'NonCoded.Allele'=Other_allele 'Coded.Af'=EAF 'Qual'=Imp_info 'Beta' 'Se'=SE 'Pval'=Pvalue";
	
	public static final String[] FORMATS = new String[]{CHARGE_FORMAT}; 
	public static final String DEFAULT_FORMAT = CHARGE_FORMAT;
	
	public static final String RUN_SCRIPT_NAME = "master_run.qsub";
	public static final String PROCESS_SCRIPT_NAME = "master_process.qsub";
	
	public static final String DATA_BUILD_1000G = "APR12";
	public static final String PROCESSED_RESULT_FILE_EXT = ".csv";
	
	public static final int QSUB_RAM_MB = 10000;
	public static final int QSUB_TIME_HRS = 8;
	public static final int QSUB_THREADS = 24;

	String FAST_LOC = "FAST";
	String dir = "/home/pankarne/chandap/ARIC.whites.impute2/";
	String indivFile = "~/ordered9489.indiv";
	String traitFile = "~/ordered9489.trait";
	String filePattern = ".impute2.gz";
	String runDir = "/home/pankarne/saonlib2/1000genomes/LnFXI/FAST/FAST_pC/";
	int covarCount = 4;
	String study = null;
	String pop = null;
	String factor = null;
	int sex = -2;
	boolean isLinear;
	
	private static FilenameFilter dirFilter = new FilenameFilter() {
        @Override
        public boolean accept(File dir, String name) {
            return (new File(dir, name)).isDirectory();
        }
    };
    private static FilenameFilter dataFileFilter = new FilenameFilter() {
        @Override
        public boolean accept(File dir, String name) {
            return name.endsWith(PROCESSED_RESULT_FILE_EXT) && name.contains(DATA_BUILD_1000G);
        }
    };
    private static FilenameFilter resultsFileFilter = new FilenameFilter() {
        @Override
        public boolean accept(File dir, String name) {
            return name.endsWith(".Linear.txt") || name.endsWith(".Logistic.txt") ;
        }
    };
    
    public static void processAndPrepareMETAL(final String studyDir, final String dataDefFile) {
        HashMap<String, HashMap<String, DataDefinitions>> dataDefs;
        try {
            dataDefs = parseDataDefinitionsFile(dataDefFile);
        } catch (IOException e) {
            throw new RuntimeException(e); // pass along
        }
        processAndPrepareMETAL(studyDir, dataDefs);
    }
    
	public static void processAndPrepareMETAL(final String studyDir, final HashMap<String, HashMap<String, DataDefinitions>> dataDefs/*, final boolean gcMetal*/) {
	    String tempName = ext.verifyDirFormat(studyDir);
	    final String studyName = ext.rootOf(tempName.substring(0, tempName.length() - 1), true);
	    final HashMap<String, DataDefinitions> popDefs = dataDefs.get(studyName);
	    File[] factorDirs = (new File(studyDir)).listFiles(dirFilter);
	    System.out.println(ext.getTime() + "]\tProcessing study " + studyName + " with " + factorDirs.length + " threads");
	    
	    ExecutorService executor = Executors.newFixedThreadPool(factorDirs.length); 

	    final double pvalThresh = 0.001;
	    
	    for (final File factorDir : factorDirs) {
            Runnable parseMetalRunnable = new Runnable() {
                @Override
                public void run() {
                    ArrayList<String> metalAnalyses = new ArrayList<String>();
                    String factorName = ext.rootOf(factorDir.getName(), true);
                    Logger factorLog = new Logger(studyDir + factorName +".log");
                    factorLog.report(ext.getTime() + "]\tBegin processing for factor " + factorName);
                    File[] popDirs = factorDir.listFiles(dirFilter);
                    factorLog.report(ext.getTime() + "]\tFound " + popDirs.length + " populations");
                    StringBuilder metalFileContents = new StringBuilder(writeMetalCRF(factorName, pvalThresh/*, gcMetal*/));
                    int foundCount = 0;
                    for (File popDir : popDirs) {
                        String popName = ext.rootOf(popDir.getName(), true);
                        
                        auto : {
                            String[] names = popDir.list(dataFileFilter);                
                            boolean parsed = false;
                            // names.length should only ever be 1 or 0
                            if (names.length == 0) {
                                factorLog.report(ext.getTime() + "]\tNo parsed datafile found for population " + popName + " - will attempt to parse results files");
                                File resultsDir = new File(popDir, "output/"); 
                                if (!resultsDir.exists()) {
                                    factorLog.report(ext.getTime() + "]\tNo output directory found for population " + popName + "!");
                                    break auto;
                                } else {
                                    try {
                                        String resultsDirPath = ext.verifyDirFormat(resultsDir.getAbsolutePath());
                                        String midOut = "concatenated.result";
                                        String finalOut = buildFinalFilename(studyName, popName, factorName, -1);
                                        String traitFile = ext.verifyDirFormat(popDir.getAbsolutePath()) + studyName + "_" + popName + "_" + factorName + ".trait";

                                        File f = new File(resultsDirPath + midOut);
                                        if (!f.exists() || f.length() <= 0) {
                                            concatResults(resultsDirPath, midOut, pvalThresh, true, true);
                                        }
                                        File f2 = new File(resultsDirPath + midOut);
                                        System.out.println("FILEPATH: " + f2.getAbsolutePath());
                                        System.out.println("FILESIZE: " + f2.length());
                                        if (f2.exists() && f2.length() > 0) {
                                            runParser(DEFAULT_FORMAT, resultsDirPath + midOut, resultsDirPath + "../" + finalOut, countValid(traitFile));
                                            factorLog.report(ext.getTime() + "]\tParsing complete.");
                                            parsed = true;
                                            names = popDir.list(dataFileFilter);
                                        } else {
                                            factorLog.reportError(ext.getTime() + "]\tError - concatenated result file is either missing or empty; cannot create a final results file without results!");
                                            break auto;
                                        }
                                    } catch (IOException e) {
                                        factorLog.report(ext.getTime() + "]\tError occurred while counting valid individuals in .trait file:");
                                        factorLog.report("\t" + e.getMessage());
                                    }
                                }
                            }
                            // Check again, in case we parsed new files
                            if (names.length == 0 && parsed) {
                                factorLog.report(ext.getTime() + "]\tError - Parsing failed; do final result files exist?");
                                break auto;
                                // uh-oh; parsing failed!
                            } else {
                                foundCount += names.length;
                                for (String name : names) {
                                    metalFileContents.append(popName).append("/").append(name).append("\t").append(popDefs.get(popName).gc).append("\n");
                                }
                            }
                        }
        
                        File femaleDir = new File(popDir, "female/");
                        File maleDir = new File(popDir, "male/");
        
                        StringBuilder metaSex = new StringBuilder(writeMetalCRF(factorName, pvalThresh/*, gcMetal*/));
                        if (femaleDir.exists() && femaleDir.isDirectory() && maleDir.exists() && maleDir.isDirectory()) {
                            String[] dataFilesF = femaleDir.list(dataFileFilter);
                            String[] dataFilesM = maleDir.list(dataFileFilter);
                            boolean parsedF = false;
                            boolean parsedM = false;
                            if (dataFilesF.length == 0) {
                                File femaleResultsDir = new File(popDir, "female/output/");
                                if (!femaleResultsDir.exists()) {
                                    factorLog.report(ext.getTime() + "]\tNo female-specific output directory found for population " + popName + "!");
                                } else {
                                    try {
                                        String resultsDirPathFemale = femaleResultsDir.getAbsolutePath();
                                        String midOutF = "concatenated.result";
                                        String finalOutF = buildFinalFilename(studyName, popName, factorName, 0);
                                        String traitFileF = ext.verifyDirFormat(femaleDir.getAbsolutePath()) + studyName + "_" + popName + "_" + factorName + "_female.trait";
                                        concatResults(resultsDirPathFemale, midOutF, pvalThresh, true, true);
                                        if (Files.exists(resultsDirPathFemale + midOutF) && Files.getSize(resultsDirPathFemale + midOutF, false) > 0) { 
                                            runParser(DEFAULT_FORMAT, ext.verifyDirFormat(resultsDirPathFemale) + midOutF, ext.verifyDirFormat(resultsDirPathFemale) + "../" + finalOutF, countValid(traitFileF));
                                            factorLog.report(ext.getTime() + "]\tParsing complete.");
                                            parsedF = true;
                                            dataFilesF = femaleDir.list(dataFileFilter);
                                        } else {
                                            factorLog.reportError(ext.getTime() + "]\tError - concatenated result file is either missing or empty; cannot create a final results file without results!");
                                        }
                                    } catch (IOException e) {
                                        factorLog.report(ext.getTime() + "]\tError occurred while counting valid individuals in .trait file:");
                                        factorLog.report("\t" + e.getMessage());
                                    }
                                }
                            }
                            if (dataFilesF.length == 0 && parsedF) {
                                factorLog.report(ext.getTime() + "]\tParsing failed - do final result files exist?");
                                // uh-oh; parsing failed!
                            } else {
                                for (String dataF : dataFilesF) {
                                    metaSex.append("female/").append(dataF).append("\t").append(popDefs.get(popName).gc).append("\n");
                                }
                            }
                            if (dataFilesM.length == 0) {
                                File maleResultsDir = new File(popDir, "male/output/");
                                if (!maleResultsDir.exists()) {
                                    factorLog.report(ext.getTime() + "]\tNo male-specific output directory found for population " + popName + "!");
                                } else {
                                    try {
                                        String resultsDirPathMale = ext.verifyDirFormat(maleResultsDir.getAbsolutePath());
                                        String midOutM = "concatenated.result";
                                        String finalOutM = buildFinalFilename(studyName, popName, factorName, 1);
                                        String traitFileM = ext.verifyDirFormat(maleDir.getAbsolutePath()) + studyName + "_" + popName + "_" + factorName + "_male.trait";
                                        concatResults(resultsDirPathMale, midOutM, pvalThresh, true, true);
                                        if (Files.exists(resultsDirPathMale + midOutM) && Files.getSize(resultsDirPathMale + midOutM, false) > 0) {
                                            runParser(DEFAULT_FORMAT, resultsDirPathMale + midOutM, resultsDirPathMale + "../" + finalOutM, countValid(traitFileM));
                                            factorLog.report(ext.getTime() + "]\tParsing complete.");
                                            parsedM = true;
                                            dataFilesM = maleDir.list(dataFileFilter);
                                        } else {
                                            factorLog.reportError(ext.getTime() + "]\tError - concatenated result file is either missing or empty; cannot create a final results file without results!");
                                        }
                                    } catch (IOException e) {
                                        factorLog.report(ext.getTime() + "]\tError occurred while counting valid individuals in .trait file:");
                                        factorLog.report("\t" + e.getMessage());
                                    }
                                }
                            }
                            if (dataFilesM.length == 0 && parsedM) {
                                factorLog.report(ext.getTime() + "]\tParsing failed - do final result files exist?");
                                // uh-oh; parsing failed!
                            } else {
                                for (String dataM : dataFilesM) {
                                    metaSex.append("male/").append(dataM).append("\t").append(popDefs.get(popName).gc).append("\n");
                                }
                            }
                            if (dataFilesF.length >= 1 && dataFilesM.length >= 1) {
                                String metalName = "metal_" + factorName + "_" + popName + "_sex.crf";
                                Files.write(metaSex.toString(), ext.verifyDirFormat(popDir.getAbsolutePath()) + metalName);
                                metalAnalyses.add(ext.verifyDirFormat(popDir.getAbsolutePath()) + metalName);
                            }
                        }
                    }
                    if (popDirs.length > 1 && foundCount > 1) {
                        String metalName = "metal_" + factorName + ".crf";
                        Files.write(metalFileContents.toString(), ext.verifyDirFormat(factorDir.getAbsolutePath()) + metalName);
                        metalAnalyses.add(ext.verifyDirFormat(factorDir.getAbsolutePath()) + metalName);
                    }
                    factorLog.report(ext.getTime() + "]\tProcessing complete - will now run " + metalAnalyses.size() + " METAL analyses.");
                    for (String metalCRF : metalAnalyses) { 
                        try {
                            factorLog.report("Running METAL analysis: " + metalCRF);
                            // Runtime.exec doesn't play well with '~', so we have to find the location of the park.jar file
                            String path = FAST.class.getProtectionDomain().getCodeSource().getLocation().getPath();
                            String decodedPath = path;
                            try {
                                decodedPath = URLDecoder.decode(path, "UTF-8");
                            } catch (UnsupportedEncodingException e) {
                                // TODO Auto-generated catch block
                                e.printStackTrace();
                            }
                            String metalDir = ext.parseDirectoryOfFile(metalCRF);
                            CmdLine.run(null, new String[]{"java", "-cp", decodedPath, "Launch", metalCRF}, metalDir, System.out, System.err, factorLog, false);
                            factorLog.report(ext.getTime() + "]\tRunning HitWindows analysis on METAL results...");
                            try {
                                String[][] results1 = HitWindows.determine(metalDir + "topHits.xln", 0.00000005f, 500000, 0.000005f, new String[0]);
                                Files.writeMatrix(results1, metalDir + factorName + "_topHitWindows.out", "\t");
                            } catch (Exception e) {
                                factorLog.report("ERROR - " + e.getMessage());
                            }
                            // Can't run HitWindows on InvVar1.out; columns Chr and Pos are missing 
//                            try {
//                                String[][] results2 = HitWindows.determine(metalDir + factorName + "_InvVar1.out", 0.00000005f, 500000, 0.000005f, new String[0]);
//                                Files.writeMatrix(results2, metalDir + factorName + "_InvVar1_hitWindows.out", "\t");
//                            } catch (Exception e) {
//                                System.out.println("ERROR - " + e.getMessage());
//                            }
                            factorLog.report(ext.getTime() + "]\tHitWindows analysis of METAL results complete!");
                        } catch (Exception e) {
                            factorLog.report("ERROR - " + e.getMessage());
                            continue;
                        }
                    }
                }
            };
            try {
                executor.execute(parseMetalRunnable);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        executor.shutdown();
        try {
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
	}
	
	
	private static void prepareMETAL(String studyDir, String dataFile/*, boolean gcMetal*/) {
	    HashMap<String, HashMap<String, DataDefinitions>> dataDefs = null;
	    try {
            dataDefs = parseDataDefinitionsFile(dataFile);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        final String studyName = ext.rootOf(ext.verifyDirFormat(studyDir).substring(0, ext.verifyDirFormat(studyDir).length() - 1), true);
        final HashMap<String, DataDefinitions> popDefs = dataDefs.get(studyName);
        
	    File[] factorDirs = (new File(studyDir)).listFiles(dirFilter);
	    
	    StringBuilder runMetal = new StringBuilder();
	    
	    double pvalThresh = 0.001;
	    
	    for (File factorDir : factorDirs) {
	        String factorName = ext.rootOf(factorDir.getName(), true);
	        StringBuilder metaFileContents = new StringBuilder(writeMetalCRF(factorName, pvalThresh/*, gcMetal*/));
	        int foundCount = 0;
	        File[] popDirs = factorDir.listFiles(dirFilter);
	        for (File popDir : popDirs) {
	            String popName = ext.rootOf(popDir.getName(), true);
	            String[] names = popDir.list(dataFileFilter);
	            for (String name : names) {
	                metaFileContents.append(popName).append("/").append(name).append("\t").append(popDefs.get(popName).gc).append("\n");
	                foundCount++;
	            }
	            File femDir = new File(popDir, "female");
	            File malDir = new File(popDir, "male");
	            if (femDir.exists() && femDir.isDirectory() && malDir.exists() && malDir.isDirectory()) {
	                StringBuilder metaSex = new StringBuilder(writeMetalCRF(factorName, pvalThresh/*, gcMetal*/));
	                String[] dataFilesF = femDir.list(dataFileFilter);
	                String[] dataFilesM = malDir.list(dataFileFilter);
	                for (String dataF : dataFilesF) {
	                    metaSex.append("female/").append(dataF).append("\n");
	                }
	                for (String dataM : dataFilesM) {
	                    metaSex.append("male/").append(dataM).append("\n");
	                }
	                if (dataFilesF.length >= 1 && dataFilesM.length >= 1) {
	                    String metalName = "metal_" + factorName + "_" + popName + "_sex.crf";
	                    Files.write(metaSex.toString(), ext.verifyDirFormat(popDir.getAbsolutePath()) + metalName);
	                    runMetal.append("cd ").append(ext.verifyDirFormat(popDir.getAbsolutePath())).append("\n");
	                    runMetal.append("java -cp ~/park.jar Launch ").append(metalName).append("\n");
	                }
	            }
	        }
	        if (popDirs.length > 1 && foundCount > 1) {
	            String metalName = "metal_" + factorName + ".crf";
	            Files.write(metaFileContents.toString(), ext.verifyDirFormat(factorDir.getAbsolutePath()) + metalName);
                runMetal.append("cd ").append(ext.verifyDirFormat(factorDir.getAbsolutePath())).append("\n");
                runMetal.append("java -cp ~/park.jar Launch ").append(metalName).append("\n");
	        }
	    }

//        Files.write(runMetal.toString(), ext.verifyDirFormat(studyDir) + "runMETALAnalyses.sh");
        Files.qsub(ext.verifyDirFormat(studyDir) + "master_runMETAL.qsub", runMetal.toString(), QSUB_RAM_MB, QSUB_TIME_HRS, QSUB_THREADS);
	}
	
	private static String writeMetalCRF(String factor, double pvalThresh/*, boolean gc*/) {
	    StringBuilder metalCRF = new StringBuilder("metal\n")
                                	    .append(factor).append("\n")
//                                	    .append("build=37\ngenomic_control=" + (gc ? "TRUE" : "FALSE") + "\nhits_p<="+pvalThresh+"\n");
                                	    .append("build=37\nhits_p<="+pvalThresh+"\n");
	    return metalCRF.toString();
	}

    public static String[] prepareFAST(String traitDir, String dataFile, String runDir, boolean isLinear, boolean run, String qsubQueue) throws IOException {
		HashMap<String, HashMap<String, HashMap<String, String>>> traits = loadTraitFiles(traitDir);
		HashMap<String, HashMap<String, DataDefinitions>> data = parseDataDefinitionsFile(dataFile);
		ArrayList<String> dirs = new ArrayList<String>();
		// TODO ensure 1-1 keymapping between traits and data maps
		// TODO ensure 1-1 keymapping between study.pop in both maps
		
		traitDir = ext.verifyDirFormat(traitDir);
		runDir = ext.verifyDirFormat(runDir);
		
		/*
		 mkdir STUDY
		 for each FACTOR:
		    mkdir FACTOR
		    for each POP:
		    	mkdir POP
				cp TRAITFILE > POP_DIR
		*/
		for (java.util.Map.Entry<String, HashMap<String, HashMap<String, String>>> entry : traits.entrySet()) {
			String study = entry.getKey();
			HashMap<String, HashMap<String, String>> factorToPopToTraitMap = entry.getValue();
			
			(new File(runDir + study)).mkdir();
			
			for (java.util.Map.Entry<String, HashMap<String, String>> factors : factorToPopToTraitMap.entrySet()) {
				String factor = factors.getKey();
				HashMap<String, String> popToTraitMap = factors.getValue();
				File factorDir = new File(runDir + study+"/"+factor+"/");
				(factorDir).mkdir();
				for (java.util.Map.Entry<String, String> popEntry : popToTraitMap.entrySet()) {
					String pop = popEntry.getKey();
					String traitFile = popEntry.getValue();
	                (new File(runDir + study+"/"+factor+"/"+pop)).mkdir();
	                Files.copyFile(traitDir + traitFile, runDir + study+"/"+factor+"/"+pop + "/" + traitFile);
				}
			}
			
		}
		
		StringBuilder masterRunScript = new StringBuilder();
		StringBuilder masterProcessScript = new StringBuilder();
		
		for (java.util.Map.Entry<String, HashMap<String, HashMap<String, String>>> entry : traits.entrySet()) {
			String study = entry.getKey();
			HashMap<String, HashMap<String, String>> factorToPopToTraitMap = entry.getValue();
			
			HashMap<String, DataDefinitions> popToDataDef = data.get(study);
			
			for (java.util.Map.Entry<String, HashMap<String, String>> factors : factorToPopToTraitMap.entrySet()) {
				String factor = factors.getKey();
				HashMap<String, String> popToTraitMap = factors.getValue();
				
				for (java.util.Map.Entry<String, String> popEntry : popToTraitMap.entrySet()) {
					String pop = popEntry.getKey();
					String traitFile = popEntry.getValue();
					String dir = new StringBuilder(runDir).append(study).append("/").append(factor).append("/").append(pop).append("/").toString();
					
					DataDefinitions dataDef = popToDataDef.get(pop);
					int covars = countCovars(traitDir + traitFile);
					
					FAST fastRun = new FAST("FAST", dataDef.dataDir, dataDef.indivFile, traitDir + traitFile, dataDef.dataSuffix, dir, covars, isLinear);
					fastRun.study = study;
					fastRun.factor = factor;
					fastRun.pop = pop;
					fastRun.sex = -1;
					fastRun.run();
					
					dirs.add(dir);
					
					masterRunScript.append("cd ").append(dir).append("\n");
					masterRunScript.append("qsub ");
					if (qsubQueue != null) {
					    masterRunScript.append("-q ").append(qsubQueue).append(" ");
					}
					masterRunScript.append(RUN_SCRIPT_NAME).append("\n");
					masterProcessScript.append("cd ").append(dir).append("\n");
					masterProcessScript.append("qsub ");
                    if (qsubQueue != null) {
                        masterProcessScript.append("-q ").append(qsubQueue).append(" ");
                    }
					masterProcessScript.append(PROCESS_SCRIPT_NAME).append("\n");
					
					if (dataDef.sexDir != null) {
					    String maleTraitFile = sexCopyTraitFile(dir + "male/", traitDir + traitFile, true);
					    String femaleTraitFile = sexCopyTraitFile(dir + "female/", traitDir + traitFile, false);
	                    FAST fastRunMale = new FAST("FAST", dataDef.sexDir, dataDef.indivFile, maleTraitFile, dataDef.sexSuffix, dir + "male/", covars, isLinear);
	                    fastRun.study = study;
	                    fastRun.factor = factor;
	                    fastRun.pop = pop;
	                    fastRun.sex = 1;
	                    fastRunMale.run();
	                    FAST fastRunFemale = new FAST("FAST", dataDef.sexDir, dataDef.indivFile, femaleTraitFile, dataDef.sexSuffix, dir + "female/", covars, isLinear);
	                    fastRun.study = study;
	                    fastRun.factor = factor;
	                    fastRun.pop = pop;
	                    fastRun.sex = 0;
	                    fastRunFemale.run();
	                    masterRunScript.append("cd ").append(dir).append("male/\n");
	                    masterRunScript.append("qsub ");
	                    if (qsubQueue != null) {
	                        masterRunScript.append("-q ").append(qsubQueue).append(" ");
	                    }
	                    masterRunScript.append(RUN_SCRIPT_NAME + "\n");
	                    masterRunScript.append("cd ").append(dir).append("female/\n");
	                    masterRunScript.append("qsub ");
	                    if (qsubQueue != null) {
	                        masterRunScript.append("-q ").append(qsubQueue).append(" ");
	                    }
	                    masterRunScript.append(RUN_SCRIPT_NAME + "\n");
	                    masterProcessScript.append("cd ").append(dir).append("male/\n");
	                    masterProcessScript.append("qsub ");
	                    if (qsubQueue != null) {
	                        masterProcessScript.append("-q ").append(qsubQueue).append(" ");
	                    }
	                    masterProcessScript.append(PROCESS_SCRIPT_NAME + "\n");
	                    masterProcessScript.append("cd ").append(dir).append("female/\n");
	                    masterProcessScript.append("qsub ");
	                    if (qsubQueue != null) {
	                        masterProcessScript.append("-q ").append(qsubQueue).append(" ");
	                    }
	                    masterProcessScript.append(PROCESS_SCRIPT_NAME + "\n");
					}
				}
			}
			String metalCmd = "java -cp ~/park.jar gwas.FAST rundir=" + runDir + study + " data=" + dataFile + (qsubQueue == null ? "" : " qsub=" + qsubQueue) + " -process";
			Files.qsub(runDir + "step4_" + study + "_metaAnalyzeFAST.qsub", metalCmd, QSUB_RAM_MB, QSUB_TIME_HRS, QSUB_THREADS);
			Files.write("qsub" + (qsubQueue == null ? "" : " -q " + qsubQueue) + " step4_" + study + "_metaAnalyzeFAST.qsub", runDir + "step4_" + study + "_metaAnalyzeFAST.sh");
			Files.chmod(runDir + "step4_" + study + "_metaAnalyzeFAST.sh");
		}
		
		Files.write(masterRunScript.toString(), runDir+"step2_runFAST.sh");
        Files.write(masterProcessScript.toString(), runDir+"step3_processFAST.sh");
		Files.chmod(runDir+"step2_runFAST.sh");
		Files.chmod(runDir+"step3_processFAST.sh");
		
		
		if (run) {
		    CmdLine.run("./step2_runFAST.sh", runDir);
		}
		
		return dirs.toArray(new String[dirs.size()]);
	}

    public FAST(String FASTloc, String dataDir, String indivFile, String traitFile, String dataFileSuffix, String runDir, int covarCount, boolean isLinear) {
		this.FAST_LOC = FASTloc;
		this.dir = ext.verifyDirFormat(dataDir);
		this.indivFile = indivFile;
		this.traitFile = traitFile;
		this.filePattern = dataFileSuffix;
		this.runDir = ext.verifyDirFormat(runDir);
		this.isLinear = isLinear;
		this.covarCount = covarCount;
	}
	
	public void run() throws IOException {
		String[] dataFiles = (new File(dir)).list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(filePattern);
			}
		});
		
		PrintWriter scriptInputWriter = new PrintWriter(new FileWriter(runDir + "input.txt"));
		
		for (int i = 0; i < dataFiles.length; i++) {
			String chr = dataFiles[i].substring(3, 5);
			if (chr.charAt(1) == '.') {
				chr = "" + chr.charAt(0);
			}
			if (chr.charAt(0) == 'X' || chr.charAt(0) == 'x') {
			    chr = "23";
			}
		
			StringBuilder fastString = new StringBuilder(FAST_LOC)
				.append(" --mode genotype --impute2-geno-file ")
				.append(dir)
				.append(dataFiles[i])
				.append(" --impute2-info-file ")
				.append(dir)
				.append(dataFiles[i].substring(0, dataFiles[i].length() - 3)) // TODO assuming files are gzipped
				.append("_info --indiv-file ") // assuming info files end with _info
				.append(indivFile)
				.append(" --trait-file ")
				.append(traitFile)
				.append(" --num-covariates ")
				.append(covarCount)
				.append(isLinear ? " --linear-snp " : " --logistic-snp")
				.append(" --chr ")
				.append(chr)
				.append(" --out-file ")
				.append(runDir)
				.append("output/")
				.append(dataFiles[i].substring(0, dataFiles[i].length() - 3))
				.append(".out");
			
			scriptInputWriter.println(fastString.toString());
			
		}
		
		scriptInputWriter.flush();
		scriptInputWriter.close();
		
		String command = "java -cp ~/park.jar one.ScriptExecutor file=\""+runDir+"input.txt\" token=took threads="+QSUB_THREADS;
		String procFileOut = buildFinalFilename();
		String processCommand = "cd \"" + runDir + "\"\njava -cp ~/park.jar gwas.FAST -convert -concat -writePVals -hitWindows out=\"" + procFileOut + "\" results=\""+runDir+"output/\" trait=\""+traitFile + "\"";
		Files.qsub(runDir + RUN_SCRIPT_NAME, command, QSUB_RAM_MB, QSUB_TIME_HRS, QSUB_THREADS);
		Files.qsub(runDir + PROCESS_SCRIPT_NAME, processCommand, QSUB_RAM_MB, QSUB_TIME_HRS, QSUB_THREADS);
		(new File(runDir + "output/")).mkdirs();
	}
	
	private String buildFinalFilename() {
	    return FAST.buildFinalFilename(this.study, this.pop, this.factor, this.sex);
	}
	
	private static String buildFinalFilename(String study, String pop, String factor, int sex) {
	    // TODO could cache this result
	    StringBuilder procFileOut = new StringBuilder();
	    if (study != null) {
            procFileOut.append(study).append("_");
        }
        if (pop != null) {
            procFileOut.append(pop).append("_");
        }
        if (factor != null) {
            procFileOut.append(factor).append("_");
        }
        procFileOut.append(DATA_BUILD_1000G).append("_");
        if (sex == 0) {
            procFileOut.append("CHRXF_");
        } else if (sex == 1) {
            procFileOut.append("CHRXM_");
        } else {
            procFileOut.append("AUTO_");
        }
        procFileOut.append((new SimpleDateFormat("ddMMMyyyy")).format(new Date()).toUpperCase());
        procFileOut.append(PROCESSED_RESULT_FILE_EXT);
        return procFileOut.toString();
    }

    public static HashMap<String, HashMap<String, HashMap<String, String>>> loadTraitFiles(String traitDir) {
    	String[] files = (new File(traitDir)).list(new FilenameFilter() {
    		@Override
    		public boolean accept(File dir, String name) {
    			return name.split("_").length == 3 && name.endsWith(".trait");
    		}
    	});
    	
    	HashMap<String, HashMap<String, HashMap<String, String>>> studyToFactorToPopToFile = new HashMap<String, HashMap<String,HashMap<String,String>>>();
    	
    	for (String file : files) {
    		String[] pts = file.substring(0, file.lastIndexOf(".")).split("_");
    		String study = pts[0];
    		String pop = pts[1];
    		String factor = pts[2];
    		HashMap<String, HashMap<String, String>> factorMap = studyToFactorToPopToFile.get(study);
    		if (factorMap == null) {
    			factorMap = new HashMap<String, HashMap<String,String>>();
    			studyToFactorToPopToFile.put(study, factorMap);
    		}
    		HashMap<String, String> popMap = factorMap.get(factor);
    		if (popMap == null) {
    			popMap = new HashMap<String, String>();
    			factorMap.put(factor, popMap);
    		}
    		popMap.put(pop, file);
    	}
    	return studyToFactorToPopToFile;
    }

    public static HashMap<String, HashMap<String, DataDefinitions>> parseDataDefinitionsFile(String file) throws IOException {
    	HashMap<String, HashMap<String, DataDefinitions>> defs = new HashMap<String, HashMap<String, DataDefinitions>>();
    	
    	BufferedReader reader = Files.getAppropriateReader(file);
    	String line = null;
    	while ((line = reader.readLine()) != null) {
    		String[] parts = line.split("\t");
    		// currently index based
    		DataDefinitions dd = new DataDefinitions();
    		dd.study = parts[0];
    		dd.popcode = parts[1];
    		dd.dataDir = parts[2];
    		dd.dataSuffix = parts[3];
    		if (parts.length == 5) {
    		    dd.indivFile = parts[4];
    		} else if (parts.length == 6) {
    		    dd.indivFile = parts[4];
    		    dd.gc = Float.parseFloat(parts[5]);
    		} else if (parts.length == 7) {
    		    dd.sexDir = parts[4];
    		    dd.sexSuffix = parts[5];
    		    dd.indivFile = parts[6];
    		} else if (parts.length == 8) {
                dd.sexDir = parts[4];
                dd.sexSuffix = parts[5];
                dd.indivFile = parts[6];
    		    dd.gc = Float.parseFloat(parts[7]);
    		} else {
    		    throw new RuntimeException("ERROR - malformed data.txt file!  Valid columns are: | STUDY | POP-CODE | DATA DIRECTORY | DATAFILE SUFFIX | [SEX DATAFILE DIRECTORY | SEX DATAFILE SUFFIX] | INDIV FILE | [GC VALUE] |");
    		}
    		
    		HashMap<String, DataDefinitions> defsMap = defs.get(dd.study);
    		if (defsMap == null) {
    			defsMap = new HashMap<String, DataDefinitions>();
    			defs.put(dd.study, defsMap);
    		}
    		defsMap.put(dd.popcode, dd);
    	}
    	return defs;
    }

    private static int countCovars(String traitFile) {
//		#Fam_ID	Ind_ID	Dad_ID	Mom_ID	Sex	Phenotype	Age	PC1	PC2	Sex
		String[] hdr = Files.getHeaderOfFile(traitFile, null);
		return hdr.length - 6;
	}

    private static String sexCopyTraitFile(String destDir, String traitFile, boolean male) throws IOException {
        (new File(destDir)).mkdirs();
        BufferedReader reader = Files.getAppropriateReader(traitFile);
        String newFile = destDir + ext.rootOf(traitFile, true) + "_" + (male ? "male" : "female") + ".trait";
        PrintWriter writer = Files.getAppropriateWriter(newFile);
        
        String line = null;
        writer.println(reader.readLine());
        while((line = reader.readLine()) != null) {
            String sexStr = line.split("\t")[4];
            if ((Integer.parseInt(sexStr) == 1 && male) ||
                    ((Integer.parseInt(sexStr) != 1 && !male))) {
                writer.println(line);
            } else {
                String[] parts = line.split("\t");
                String newLine = parts[0] + "\t" + parts[1] + "\t" + parts[2] + "\t" + parts[3] + "\t" + parts[4] + "\t" + parts[5]; 
                for (int i = 0; i < parts.length - 6; i++) {
                    newLine += "\tNA";  
                }
                writer.println(newLine);
                
            }
        }
        writer.flush();
        writer.close();
        reader.close();
        
        return newFile;
    }

    public static class DataDefinitions {
    	public String study;
    	public String popcode;
    	public String dataDir;
    	public String dataSuffix;
    	public String sexDir;
    	public String sexSuffix;
    	public String indivFile;
    	public double gc = -1; // -1=OFF, -9=ON, other values used as given
    }

    private static void runParser(String FORMAT, String concattedResultsFile, String outFileName, int count) {
    	System.out.println(ext.getTime() + "]\tParsing results file according to given FORMAT...");
    	String finalFormat = concattedResultsFile + " tab out=" + outFileName + FORMAT.replace(COUNT_SYMB, count + "");
    	String[] args = ext.removeQuotes(finalFormat).trim().split("[\\s]+");
    	GenParser.parse(args, new Logger());
    	System.out.println(ext.getTime() + "]\tParsing complete!");
    }

    private static void concatResults(String resultsDirectory, String resultsFile, double pvalThreshold, boolean writePValThresh, boolean runHitWindows) {
		String resultsDir = ext.verifyDirFormat(resultsDirectory);
		String[] filenames = (new File(resultsDir)).list(resultsFileFilter);
		Arrays.sort(filenames, new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				String[] pts1 = o1.split("\\.");
				String[] pts2 = o2.split("\\.");
				
				Integer chr1 = pts1[0].substring(3).charAt(0) == 'X' ? 23 : pts1[0].substring(3).charAt(0) == 'Y' ? 24 : Integer.valueOf(pts1[0].substring(3));
				Integer chr2 = pts2[0].substring(3).charAt(0) == 'X' ? 23 : pts2[0].substring(3).charAt(0) == 'Y' ? 24 : Integer.valueOf(pts2[0].substring(3));
				
				int chrComp = chr1.compareTo(chr2);
				if (chrComp != 0) return chrComp;
				
				BigInteger pos1 = new BigInteger(pts1[1]);
				BigInteger pos2 = new BigInteger(pts2[1]);

				int posComp1 = pos1.compareTo(pos2);
				if (posComp1 != 0) return posComp1;
				
				BigInteger pos12 = new BigInteger(pts1[2]);
				BigInteger pos22 = new BigInteger(pts2[2]);
				
				return pos12.compareTo(pos22);
			}
		});
		
		PrintWriter writer = Files.getAppropriateWriter(resultsDir + resultsFile);
		int[] indices = null;
		
		PrintWriter pvalWriter = writePValThresh ? Files.getAppropriateWriter(resultsDir + "meetsPVal_" + ".out") : null;
		boolean first = true;
		System.out.print(ext.getTime() + "]\tConcatenating results files: <");
		for (String str : filenames) {
			BufferedReader reader;
			try {
				reader = Files.getAppropriateReader(resultsDir + str);
				String line = reader.readLine();
				if (first) {
					if (writePValThresh) {
						indices = ext.indexFactors(new String[][] {Aliases.PVALUES}, line.split("[\\s]+"), false, true, true, false);
					}
					writer.println(line);
					first = false;
				}
				while((line = reader.readLine()) != null) {
					if (writePValThresh && indices[0] != -1) {
						double pval = Double.parseDouble(line.split("[\\s]+")[indices[0]]);
						if (pval <= pvalThreshold) {
							pvalWriter.println(line);
						}
					}
					writer.println(line);
				}
				reader.close();
				System.out.print("-");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		System.out.println(">");
		writer.flush();
		writer.close();
		if (writePValThresh) {
			pvalWriter.flush();
			pvalWriter.close();
		}
		System.out.println(ext.getTime() + "]\tConcatenation complete!");
		if (runHitWindows) {
		    File f = new File(resultsDir + resultsFile);
            System.out.println("FILEPATH: " + f.getAbsolutePath());
            System.out.println("FILESIZE: " + f.length());
            // TODO PROBLEMS WITH FILES.EXISTS and FILES.SIZE!!!!!!!
		    if (f.exists() && f.length() > 0) {
    			System.out.println(ext.getTime() + "]\tRunning HitWindows analysis...");
    			String[][] results = HitWindows.determine(resultsDir + resultsFile, 0.00000005f, 500000, 0.000005f, new String[0]);
    			Files.writeMatrix(results, resultsDir + "hits.out", "\t");
    			System.out.println(ext.getTime() + "]\tHitWindows analysis complete!");
		    } else {
		        System.err.println(ext.getTime() + "]\tError - Can't run HitWindows; input file [" + resultsDir + resultsFile + "] either doesn't exist or is empty!");
		    }
		}
	}
	
	private static int countValid(String traitFile) throws IOException {
	    BufferedReader reader = Files.getAppropriateReader(traitFile);
	    reader.readLine();
	    String line = null;
	    int count = 0;
	    read: while((line = reader.readLine()) != null) {
	        String[] pts = line.split("\t");
	        for (String str : pts) {
	            if (ext.isMissingValue(str)) {
	                continue read;
	            }
	        }
	        count++;
	    }
	    return count;
	}
	
	private static void extract(String data, String snpList, String outfileBase) {
	    HashMap<String, HashMap<String, DataDefinitions>> defs;
	    try {
            defs = parseDataDefinitionsFile(data);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
            return;
        }
	    
	    HashMap<String, String[]> studyPopIDs = new HashMap<String, String[]>();
	    HashMap<String, HashMap<String, String[]>> studyPopDataPerSNP = new HashMap<String, HashMap<String, String[]>>();
	    HashMap<String, HashMap<String, String>> studyPopInfoPerSNP = new HashMap<String, HashMap<String, String>>();
	    HashSet<String> snpSet = new HashSet<String>();
	    
	    String[] snps = null;
	    if (new File(snpList).exists()) {
	        snps = HashVec.loadFileToStringArray(snpList, false, new int[]{0}, false);
	    } else {
	        snps = snpList.split(",");
	    }
	    
	    for (String snp : snps) {
	        snpSet.add(snp);
	    }
	    
	    for (Entry<String, HashMap<String, DataDefinitions>> studyEntry : defs.entrySet()) {
	        String study = studyEntry.getKey();
	        
	        for (Entry<String, DataDefinitions> popEntry : studyEntry.getValue().entrySet()) {
	            String pop = popEntry.getKey();
	            final DataDefinitions dataDefs = popEntry.getValue();
	            
	            String[] iids = HashVec.loadFileToStringArray(dataDefs.indivFile, false, new int[]{0}, false);
	            studyPopIDs.put(study + "\t" + pop, iids);
	            
    	        String[] chrDataFiles = (new File(dataDefs.dataDir)).list(new FilenameFilter() {
                    @Override
                    public boolean accept(File dir, String name) {
                        // TODO assuming datafile names start with "chr#."
                        return name.startsWith("chr") && name.endsWith(dataDefs.dataSuffix);
                    }
                });
                
    	        HashMap<String, String> dataFilesPerSNP = new HashMap<String, String>();
    	        HashMap<String, Integer> dataLinePerSNP = new HashMap<String, Integer>();
                HashSet<String> tempSnps = new HashSet<String>(snpSet);
                try {
                
                    search : for (String dFile : chrDataFiles) {
                        String infoFile = ext.verifyDirFormat(dataDefs.dataDir) + dFile.substring(0, dFile.length() - 3) + "_info";
                        
                        BufferedReader reader = Files.getAppropriateReader(infoFile);
                        int cnt = 0;
                        String line = null;
                        reader.readLine(); // header
                        while ((line = reader.readLine()) != null) {    
                            String mkr = line.split("[\\s]+")[1];
                            if (tempSnps.contains(mkr)) {
                                dataFilesPerSNP.put(mkr, ext.verifyDirFormat(dataDefs.dataDir) + dFile);
                                dataLinePerSNP.put(mkr, cnt);
                                HashMap<String, String> snpInfoMap = studyPopInfoPerSNP.get(study + "\t" + pop);
                                if (snpInfoMap == null) {
                                    snpInfoMap = new HashMap<String, String>();
                                    studyPopInfoPerSNP.put(study + "\t" + pop, snpInfoMap);
                                }
                                snpInfoMap.put(mkr, line);
                                tempSnps.remove(mkr);
                                if (tempSnps.isEmpty()) {
                                    reader.close();
                                    break search;
                                }
                            }
                            cnt++;
                        }
                        reader.close();
                    }
                    
                    if (!tempSnps.isEmpty()) {
                        System.err.println("Error - couldn't find the following markers: {" + tempSnps.toString() + "} for " + study + "_" + pop);
                    }
                    
                    for (String snp : snpSet) {
                        BufferedReader reader = Files.getAppropriateReader(dataFilesPerSNP.get(snp));
                        String line = reader.readLine();
                        int lineCnt = dataLinePerSNP.get(snp).intValue();
                        for (int i = 0; line != null && i < lineCnt; i++) {
                            line = reader.readLine();
                        }
                        reader.close();
                        HashMap<String, String[]> snpDataMap = studyPopDataPerSNP.get(study + "\t" + pop);
                        if (snpDataMap== null) {
                            snpDataMap = new HashMap<String, String[]>();
                            studyPopDataPerSNP.put(study + "\t" + pop, snpDataMap);
                        }
                        snpDataMap.put(snp, line.split("[\\s]+"));
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
	        }
	    }
	    
	    int offset = 5;
	    for (Entry<String, HashMap<String, String[]>> dataEntry : studyPopDataPerSNP.entrySet()) {
	        String studyPop = dataEntry.getKey();
	        HashMap<String, String[]> line = dataEntry.getValue();
	        String[] ids = studyPopIDs.get(studyPop);
	        
	        ArrayList<String> snpOrder = new ArrayList<String>(line.keySet());
	        PrintWriter writer = Files.getAppropriateWriter(outfileBase + "_" + studyPop.replaceAll("\t", "_") + ".data");
	        
	        StringBuilder sb = new StringBuilder();
            sb.append("IID");
            for (String snp : snpOrder) {
                sb.append("\t").append(snp);
            }
            writer.println(sb.toString());
            sb = new StringBuilder();
            for (String snp : snpOrder) {
                for (int i = 0; i < offset; i++) {
                    sb.append(line.get(snp)[i]);
                }
            }
	        for (int i = 0; i < ids.length; i++) {
	            sb = new StringBuilder();
	            sb.append(ids[i]);
	            for (String snp : snpOrder) {
	                double geno2 = Double.parseDouble(line.get(snp)[offset + (3 * i) + 1]);
	                double geno3 = Double.parseDouble(line.get(snp)[offset + (3 * i) + 2]);
	                double geno = (geno2 + (2 * geno3));
	                
	                sb.append("\t").append(geno);
	            }
	            
	            writer.println(sb.toString());
	        }
	        writer.flush();
	        writer.close();
	    }
	    
        String infoHeader = "snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0";
	    for (Entry<String, HashMap<String, String>> infoEntry : studyPopInfoPerSNP.entrySet()) {
	        String studyPop = infoEntry.getKey();
	        HashMap<String, String> line = infoEntry.getValue();
	        
	        PrintWriter writer = Files.getAppropriateWriter(outfileBase + "_" + studyPop.replaceAll("\t", "_") + ".info");
	        writer.println(infoHeader);
	        for (String info : line.values()) {
	            writer.println(info);
	        }
	        writer.flush();
	        writer.close();
	    }
	    
	    
	}
	
	public static void main(String[] args) {
//	    results="F:/FAST analysis/FVIII/output/" out=finalResults.txt -concat
//	    trait="F:/FAST analysis/construct test/" data="F:/FAST analysis/construct test/data.txt" -prep
//	    data="F:/FAST analysis/construct test/ARIC/" -metal
		int numArgs = args.length;
		String fast = "~/FAST";
		String data = "~/data/";
		String indiv = "~/indiv.txt";
		String trait = "~/trait.txt";
		String suffix = ".impute2.gz";
		String run = ext.verifyDirFormat(System.getProperty("user.dir"));
		String snps = null;
		int covars = 0;
		String results = "~/FAST/output/";
		String out = "finalResults.txt";
		boolean concat = false;
//		boolean gc = true;
		String qsub = null;
		
		int format = 0;
		boolean convert = false;
		int count = -1;
		
		double pval = 0.0001;
		boolean printPVals = false;
		boolean runHitWindows = false;
		
		boolean prep = false;
		boolean metal = false; 
		boolean runFAST = false; 
		boolean process = false;
		boolean linear = true;
		
		String usage = "gwas.FAST requires 2-8 arguments:\n" + 
		                "   (1) Data file defining input files, in tab-delimited format (i.e. data=data.txt (not the default))\n" +
		                "         Data file must be tab-delimited in the following order:\n" +
		                "             POPULATION\n" + 
		                "             FACTOR\n" + 
		                "             DATA FOLDER\n" + 
		                "             DATA FILE SUFFIX\n" +
		                "             (Optional:)\n" + 
		                "                 SEX DATA FOLDER\n" +
		                "                 SEX DATA FILE SUFFIX\n" +
		                "             PATH TO .indiv FILE \n" + 
		                "             (Optional:)\n" + 
		                "                 GC CORRECTION VALUE (-1=OFF, -9=ON, other values used as given)\n" + 
		                "   (2) Path to folder containing .trait files (i.e. trait=~/traits/ (not the default))\n" + 
		                "   (3) Full-path to the directory in which you want to run these scripts (must include a folder named 'output') (i.e. rundir=" + run + " (default))\n" +
		                "   (4) -prep flag\n" +
	                    "   (5) OPTIONAL: -run flag to run FAST analyses after preparing FAST scripts\n" + 
		                "   (6) OPTIONAL: specify the batch queue through which to run qsub files (i.e. qsub=" + qsub + " (default))" +
		                " OR: \n" +
		                "   (1) Path to population folder containing sub-folders for FAST analyses (i.e. rundir=~/FAST/ARIC/ (not the default))\n" +
		                "   (2) Data file defining input files, in tab-delimited format (i.e. data=data.txt (not the default))\n" +
		                "   (3) -process flag\n" + 
		                "\n" +
		                "  These two options (-prep and -process) are, given no errors, the only commands needed to run multiple FAST analyses from start to finish.\n" + 
		                "  However, FAST includes other options for partial processing:\n" + 
		                "\n" + 
					   "   (1) Full-path to FAST script (including /FAST) (i.e. fast=" + fast + " (default))\n" + 
					   "   (2) Full-path to data directory (i.e. data=" + data + " (default))\n" +
					   "   (3) Full-path to .indiv file (i.e. indiv=" + indiv + " (default))\n" +
					   "   (4) Full-path to .trait file (i.e. trait=" + trait + " (default))\n" +
					   "   (5) Suffix by which to identify data files in the data directory (i.e. suffix=" + suffix + " (default))\n" +
					   "   (6) Full-path to the directory in which you want to run these scripts (must include a folder named 'output') (i.e. rundir=" + run + " (default))\n" +
					   "   (7) Number of covariates in .trait file (i.e. covars=" + covars + " (default))\n" +
					   " OR \n" +
					   "   (1) Flag to indicate results processing is desired (i.e. -concat (not the default))\n" +
					   "   (2) Path to directory with results files (i.e. results=" + results + " (default))\n" +
					   "   (3) Desired name of concatenated result file (i.e. out=" + out + " (default))\n" +
					   "   (4) -writePVals \n" +
					   "   (5) P-Value threshold (i.e. pval=" + pval + "\n" + 
					   "   (6) -hitWindows \n" + 
					   " OR \n" +
					   "   (1) Flag to indicate format conversion processing is desired (i.e. -convert (not the default))\n" +
					   "   (2) Path to concatenated result files (i.e. results=" + results + " (default))\n" +
					   "   (3) Desired name of processed result file (i.e. out=" + out + " (default))\n" +
					   "   (4) Format flag: (i.e. format=" + format + " (default))\n" + 
					   "              0: CHARGE format \n" +
					   "   (5) Number of individuals in analysis (i.e. count=" + count + " (not the default))\n" +
					   "     (5a) OPTIONAL: specify a .trait file instead of a count value, and the non-NaN and non-NA will be summed as the count value (i.e. trait= (not the default))\n" + 
					   " OR \n" +
					   "   -concat and -convert can be combined:\n" +
					   "   (1) Both -concat and -convert flags\n" +
					   "   (2) Path to directory with results files (i.e. results=" + results + " (default))\n" +
					   "   (3) Desired name of processed result file (i.e. out=" + out + " (default))\n" +
					   "   (4) Format flag: (i.e. format=" + format + " (default))\n" + 
					   "           FORMATS:\n" + 
					   "               0: CHARGE format \n" +
					   "   (5) Number of individuals in analysis (i.e. count=" + count + " (not the default))\n" +
					   "     (5a) OPTIONAL: specify a .trait file instead of a count value, and the non-NaN and non-NA will be summed as the count value (i.e. trait= (not the default))\n" + 
					   "   (6) -writePVals \n" +
					   "   (7) P-Value threshold (i.e. pval=" + pval + "\n" + 
					   "   (8) -hitWindows \n" + 
					   " OR \n" +
					   "   (1) Path to study directory with fully-parsed results files (i.e. rundir=" + data + " (default))\n" +
                       "   (2) Data file defining input files, in tab-delimited format (i.e. data=data.txt (not the default))\n" +
					   "   (3) -metal flag to create meta-analysis scripts to run METAL program\n " +
					   "\n" +
					   " FAST also provides a function to extract SNP-specific genotype probabilities and info data from data files:\n" + 
                       "   (1) Data file defining input files, in tab-delimited format (i.e. data=data.txt (not the default))\n" +
                       "   (2) List of SNPs (i.e. snps=rs1000001,rs10000002,rs10000004 (not the default))\n" +
                       "   (3) Desired name of processed result file (i.e. out=" + out + " (default))\n";
		
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("fast=")) {
				fast = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("data=")) {
				data = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("snps=")) {
			    snps = args[i].split("=")[1];
			    numArgs--;
			} else if (args[i].startsWith("indiv=")) {
				indiv = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("trait=")) {
				trait = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("suffix=")) {
				suffix = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("rundir=")) {
				run = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("results=")) {
				results = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("format=")) {
				format = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("count=")) {
				count = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				out = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("covars=")) {
				covars = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("pval=")) {
				pval = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("qsub=")) {
			    qsub = args[i].split("=")[1];
			    numArgs--;
//			} else if (args[i].startsWith("gcMetal=")) {
//			    gc = ext.parseBooleanArg(args[i]);
//			    numArgs--;
			} else if (args[i].startsWith("-concat")) {
				concat = true;
				numArgs--;
			} else if (args[i].startsWith("-convert")) {
				convert = true;
				numArgs--;
			} else if (args[i].startsWith("-writePVals")) {
				printPVals = true;
				numArgs--;
			} else if (args[i].startsWith("-hitWindows")) {
				runHitWindows = true;
				numArgs--;
			} else if (args[i].startsWith("-prep")) {
			    prep = true;
			    numArgs--;
			} else if (args[i].startsWith("-run")) {
			    runFAST = true;
			    numArgs--;
			} else if (args[i].startsWith("-metal")) {
			    metal = true;
			    numArgs--;
			} else if (args[i].startsWith("-process")) {
			    process = true;
			    numArgs--;
			} else if (args[i].startsWith("-logistic")) { 
			    linear = false;
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
		    if (metal) {
		        prepareMETAL(run, data);
		    } else if (prep) {
		        prepareFAST(trait, data, run, linear, runFAST, qsub);
		    } else if (process) {
		        processAndPrepareMETAL(run, data);
		    } else if (concat && convert) {
				String midOut = "concatenated.result";
				concatResults(results, midOut, pval, printPVals, runHitWindows);
				runParser(FORMATS[format], ext.verifyDirFormat(results) + midOut, ext.verifyDirFormat(results) + "../" + out, count == -1 ? countValid(trait) : count);
			} else if (concat) {
				concatResults(results, out, pval, printPVals, runHitWindows);
			} else if (convert) {
				runParser(FORMATS[format], results, out, count == -1 ? countValid(trait) : count);
			} else if (snps != null) {
			    extract(data, snps, out);
			} else {
				new FAST(fast, data, indiv, trait, suffix, run, covars, linear).run();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
	

