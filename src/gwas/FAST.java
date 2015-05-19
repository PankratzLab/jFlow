package gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;

import parse.GenParser;
import common.Aliases;
import common.Files;
import common.Logger;
import common.ext;

public class FAST {
	
	private static final String COUNT_SYMB = "<>";
	private static final String CHARGE_FORMAT = " 'SNP.id'=Markername 'Chr'=Chr 'Pos'=Pos $#" + COUNT_SYMB + "=N 'Coded.Allele'=Effect_allele 'NonCoded.Allele'=Other_allele 'Coded.Af'=EAF 'Qual'=Imp_info 'Beta' 'Se'=SE 'Pval'=Pvalue";
	
	public static final String[] FORMATS = new String[]{CHARGE_FORMAT}; 
	
	public static final String RUN_SCRIPT_NAME = "master_run.qsub";
	public static final String PROCESS_SCRIPT_NAME = "master_process.qsub";
	
	public static final String DATA_BUILD_1000G = "APR12";
	public static final String PROCESSED_RESULT_FILE_EXT = ".csv.gz";
	
	public static final int METAL_QSUB_RAM_MB = 10000;
	public static final int METAL_QSUB_TIME_HRS = 8;
	public static final int METAL_QSUB_THREADS = 24;
	
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
	
	private static void prepareMETAL(String studyDir) {
	    FilenameFilter dirFilter = new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return (new File(dir, name)).isDirectory();
            }
        };
        FilenameFilter dataFileFilter = new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(PROCESSED_RESULT_FILE_EXT) && name.contains(DATA_BUILD_1000G);
            }
        };
	    File[] factorDirs = (new File(studyDir)).listFiles(dirFilter);
	    
	    StringBuilder runMetal = new StringBuilder();
	    
	    for (File factorDir : factorDirs) {
	        String factorName = ext.rootOf(factorDir.getName(), true);
	        StringBuilder metaFileContents = new StringBuilder(writeMetalCRF(factorName));
	        int foundCount = 0;
	        File[] popDirs = factorDir.listFiles(dirFilter);
	        for (File popDir : popDirs) {
	            String popName = ext.rootOf(popDir.getName(), true);
	            String[] names = popDir.list(dataFileFilter);
	            for (String name : names) {
	                metaFileContents.append(popName).append("/").append(name).append("\n");
	                foundCount++;
	            }
	            File femDir = new File(popDir, "female");
	            File malDir = new File(popDir, "male");
	            if (femDir.exists() && femDir.isDirectory() && malDir.exists() && malDir.isDirectory()) {
	                StringBuilder metaSex = new StringBuilder(writeMetalCRF(factorName));
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
        Files.qsub(ext.verifyDirFormat(studyDir) + "master_runMETAL.qsub", runMetal.toString(), METAL_QSUB_RAM_MB, METAL_QSUB_TIME_HRS, METAL_QSUB_THREADS);
	}
	
	private static String writeMetalCRF(String factor) {
	    StringBuilder metalCRF = new StringBuilder("metal\n")
                                	    .append(factor).append("\n")
                                	    .append("build=37\ngenomic_control=TRUE\nhits_p<=0.0001\n");
	    return metalCRF.toString();
	}

    private static void prepareFAST(String traitDir, String dataFile) throws IOException {
    		HashMap<String, HashMap<String, HashMap<String, String>>> traits = loadTraitFiles(traitDir);
    		HashMap<String, HashMap<String, DataDefinitions>> data = parseFile(dataFile);
    
    		String runDir = ext.verifyDirFormat(System.getProperty("user.dir"));
    		
    		// TODO ensure 1-1 keymapping between traits and data maps
    		// TODO ensure 1-1 keymapping between study.pop in both maps
    		
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
    			
    			(new File(study)).mkdir();
    			
    			for (java.util.Map.Entry<String, HashMap<String, String>> factors : factorToPopToTraitMap.entrySet()) {
    				String factor = factors.getKey();
    				HashMap<String, String> popToTraitMap = factors.getValue();
    				File factorDir = new File(study+"/"+factor+"/");
    				(factorDir).mkdir();
    //				writeMetalCRF(factorDir.getAbsolutePath(), factor, false);
    				for (java.util.Map.Entry<String, String> popEntry : popToTraitMap.entrySet()) {
    					String pop = popEntry.getKey();
    					String traitFile = popEntry.getValue();
    	                (new File(study+"/"+factor+"/"+pop)).mkdir();
    	                Files.copyFile(traitDir + traitFile, study+"/"+factor+"/"+pop + "/" + traitFile);
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
    					
    					DataDefinitions dataDef = popToDataDef.get(pop);
    					int covars = countCovars(traitDir + traitFile);
    					FAST fastRun = new FAST("FAST", dataDef.dataDir, dataDef.indivFile, runDir+study+"/"+factor+"/"+pop+"/"+traitFile, dataDef.dataSuffix, runDir+study+"/"+factor+"/"+pop, covars);
    					fastRun.study = study;
    					fastRun.factor = factor;
    					fastRun.pop = pop;
    					fastRun.sex = -1;
    					fastRun.run();
    					
    					masterRunScript.append("cd ").append(runDir).append(study).append("/").append(factor).append("/").append(pop).append("/\n");
    					masterRunScript.append("qsub " + RUN_SCRIPT_NAME + "\n");
    					masterProcessScript.append("cd ").append(runDir).append(study).append("/").append(factor).append("/").append(pop).append("/\n");
    //					java -cp ~/park.jar gwas.FAST -convert -concat -writePVals -hitWindows out=ARIC_EA_F7_APR12_AUTO_18MAY2015.csv.gz results=/home/pankarne/shared/1000G/FAST_withSex/ARIC/F7/EA/output/ trait=/home/pankarne/shared/1000G/FAST_withSex/ARIC/F7/EA/ARIC_EA_F7.trait
    
    					
    					
    					if (dataDef.sexDir != null) {
    					    String maleTraitFile = sexCopyTraitFile(study+"/"+factor+"/"+pop+"/male/", traitDir + traitFile, true);
    					    String femaleTraitFile = sexCopyTraitFile(study+"/"+factor+"/"+pop+"/female/", traitDir + traitFile, false);
    	                    FAST fastRunMale = new FAST("FAST", dataDef.sexDir, dataDef.indivFile, runDir+maleTraitFile, dataDef.sexSuffix, runDir+study+"/"+factor+"/"+pop+"/male/", covars);
    	                    fastRun.study = study;
    	                    fastRun.factor = factor;
    	                    fastRun.pop = pop;
    	                    fastRun.sex = 1;
    	                    fastRunMale.run();
    	                    FAST fastRunFemale = new FAST("FAST", dataDef.sexDir, dataDef.indivFile, runDir+femaleTraitFile, dataDef.sexSuffix, runDir+study+"/"+factor+"/"+pop+"/female/", covars);
    	                    fastRun.study = study;
    	                    fastRun.factor = factor;
    	                    fastRun.pop = pop;
    	                    fastRun.sex = 0;
    	                    fastRunFemale.run();
    	                    masterRunScript.append("cd ").append(runDir).append(study).append("/").append(factor).append("/").append(pop).append("/male/\n");
    	                    masterRunScript.append("qsub " + RUN_SCRIPT_NAME + "\n");
    	                    masterRunScript.append("cd ").append(runDir).append(study).append("/").append(factor).append("/").append(pop).append("/female/\n");
    	                    masterRunScript.append("qsub " + RUN_SCRIPT_NAME + "\n");
    //	                    writeMetalCRF(runDir+study+"/"+factor+"/"+pop+"/", factor, true);
    					}
    				}
    			}
    		}
    		
    		Files.write(masterRunScript.toString(), runDir+"runFAST.sh");
    		Files.chmod(runDir+"runFAST.sh");
    		
    	}

    public FAST(String FASTloc, String dataDir, String indivFile, String traitFile, String dataFileSuffix, String runDir, int covarCount) {
		this.FAST_LOC = FASTloc;
		this.dir = ext.verifyDirFormat(dataDir);
		this.indivFile = indivFile;
		this.traitFile = traitFile;
		this.filePattern = dataFileSuffix;
		this.runDir = ext.verifyDirFormat(runDir);
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
				.append(dataFiles[i].substring(0, dataFiles[i].length() - 3))
				.append("_info --indiv-file ")
				.append(indivFile)
				.append(" --trait-file ")
				.append(traitFile)
				.append(" --num-covariates ")
				.append(covarCount)
				.append(" --linear-snp ")
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
		
		int threads = 24;
		String command = "java -cp ~/park.jar one.ScriptExecutor file=\""+runDir+"input.txt\" token=took threads="+threads;
		String procFileOut = buildFinalFilename();
		String processCommand = "cd \"" + runDir + "\"\njava -cp ~/park.jar gwas.FAST -convert -concat -writePVals -hitWindows out=\"" + procFileOut + "\" results=\""+runDir+"output/\" trait=\""+traitFile + "\"";
		Files.qsub(runDir + RUN_SCRIPT_NAME, command, 10000, 8, threads);
		Files.qsub(runDir + PROCESS_SCRIPT_NAME, processCommand, 10000, 8, threads);
		(new File(runDir + "output/")).mkdirs();
	}
	
	private String buildFinalFilename() {
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

    private static HashMap<String, HashMap<String, HashMap<String, String>>> loadTraitFiles(String traitDir) {
    	String[] files = (new File(traitDir)).list(new FilenameFilter() {
    		@Override
    		public boolean accept(File dir, String name) {
    			return name.split("_").length == 3 && name.endsWith(".trait");
    		}
    	});
    	System.out.println("Found " + files.length + " trait files");
    	
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

    private static HashMap<String, HashMap<String, DataDefinitions>> parseFile(String file) throws IOException {
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
    		} else {
    		    dd.sexDir = parts[4];
    		    dd.sexSuffix = parts[5];
    		    dd.indivFile = parts[6];
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

    static class DataDefinitions {
    	String study;
    	String popcode;
    	String dataDir;
    	String dataSuffix;
    	String sexDir;
    	String sexSuffix;
    	String indivFile;
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
		String[] filenames = (new File(resultsDir)).list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(".Linear.txt");
			}
		});
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
			System.out.println(ext.getTime() + "]\tRunning HitWindows analysis...");
			String[][] results = HitWindows.determine(resultsDir + resultsFile, 0.00000005f, 500000, 0.000005f, new String[0]);
			Files.writeMatrix(results, resultsDir + "hits.out", "\t");
			System.out.println(ext.getTime() + "]\tHitWindows analysis complete!");
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
		String run = "~/runFAST/";
		int covars = 0;
		String results = "~/FAST/output/";
		String out = "finalResults.txt";
		boolean concat = false;
		
		int format = 0;
		boolean convert = false;
		int count = -1;
		
		double pval = 0.0001;
		boolean printPVals = false;
		boolean runHitWindows = false;
		
		boolean prep = false;
		boolean metal = false; 
		
		String usage = "gwas.FAST requires 3-8 arguments\n" +
					   "   (1) Path to directory containing .trait files (i.e. trait=C:/traitFiles/ (not the default))\n" +
					   "   (2) Location of data definitions file (i.e. data=C:/dataFile.txt (not the default))\n" + 
					   "   (3) -prep flag \n" + 
					   " OR \n" +
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
					   "   (1) Path to study directory with fully-parsed results files (i.e. data=" + data + " (default))\n" +
					   "   (2) -metal flag to create meta-analysis scripts to run METAL program\n " +
					   "";

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
			} else if (args[i].startsWith("-metal")) {
			    metal = true;
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
		        prepareMETAL(data);
		    } else if (prep) {
		        prepareFAST(trait, data);
		    } else if (concat && convert) {
				String midOut = "concatenated.result";
				concatResults(results, midOut, pval, printPVals, runHitWindows);
				runParser(FORMATS[format], ext.verifyDirFormat(results) + midOut, ext.verifyDirFormat(results) + "../" + out, count == -1 ? countValid(trait) : count);
			} else if (concat) {
				concatResults(results, out, pval, printPVals, runHitWindows);
			} else if (convert) {
				runParser(FORMATS[format], results, out, count == -1 ? countValid(trait) : count);
			} else {
				new FAST(fast, data, indiv, trait, suffix, run, covars).run();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
	

