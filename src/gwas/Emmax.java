package gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import stats.ProbDist;

import common.Array;
import common.Files;
import common.Logger;
import common.Matrix;
import common.ext;

public class Emmax {
	public static final String[] KINSHIP_LABELS = new String[] {"hBN", "hIBS"};
	public static final String[] KINSHIP_EXTENSIONS = new String[] {".hBN.kinf", ".hIBS.kinf"};
	public static final String[] KINSHIP_PARAMETERS = new String[] {"emmax-kin -v -h -d 10", "emmax-kin -v -h -s -d 10"};

	public static void generateScripts(String commandFullPath, String kinshipGenoDirAndNameRoot, String analysisGenoDirAndNameRoot, String phenoCovDir, String phenoNameExtOrFullPath, String covNameExtOrFullPath, String outDir, String batchDir, int qsubMemInMBs, double qsubWalltimeInHours, Logger log) {
		String[] phenoFiles = null;
		String[] covFiles = null;
		String phenoLabel;
		Vector<String> scripts;
		
		if (log==null) {
			log = new Logger();
		}

		if (!Files.exists(kinshipGenoDirAndNameRoot+".bed")) {
			log.reportError("Kinship genotype files not found at " + kinshipGenoDirAndNameRoot);
			System.exit(1);
		}

		if (!Files.exists(analysisGenoDirAndNameRoot+".bed")) {
			log.reportError("Analysis genotype files not found at " + analysisGenoDirAndNameRoot);
			System.exit(1);
		}
		
		phenoFiles = Files.list(phenoCovDir, phenoNameExtOrFullPath, false);
		if (phenoFiles.length == 0) {
			log.reportError("No pheno file is found at " + phenoCovDir);
			return;
		}
		
		covFiles = new String[phenoFiles.length];
		for (int i = 0; i < covFiles.length; i++) {
			covFiles[i] = phenoFiles[i].substring(0, phenoFiles[i].lastIndexOf(phenoNameExtOrFullPath)) + covNameExtOrFullPath;
			if (! new File(phenoCovDir + covFiles[i]).exists()) {
				covFiles[i] = null;
			}
		}

		scripts = new Vector<String>();
		for (int l = 0; l < phenoFiles.length; l++) {
			log.report("Prepping scripts for "+phenoFiles[l]+" and "+(covFiles[l] == null || Files.exists(phenoCovDir + covFiles[l])?"zero":(Files.getHeaderOfFile(phenoCovDir + covFiles[l], log).length-2))+" covariates");
			for (int i = 0; i < KINSHIP_LABELS.length; i++) {
				phenoLabel = phenoFiles[l].substring(0, phenoFiles[l].lastIndexOf(phenoNameExtOrFullPath)) + "_" + ext.rootOf(kinshipGenoDirAndNameRoot) + "_" + ext.rootOf(analysisGenoDirAndNameRoot) + "_" + KINSHIP_LABELS[i];
				new File(phenoCovDir+phenoLabel).mkdirs();
				scripts.add("cd "+phenoCovDir+phenoLabel+"/\n"
						   + "plink --bfile "+kinshipGenoDirAndNameRoot+" --keep ../"+phenoFiles[l]+" --transpose --recode12 --output-missing-genotype 0 --out forKinship --noweb\n"
						   + "plink --bfile "+analysisGenoDirAndNameRoot+" --keep ../"+phenoFiles[l]+" --freq --out mafCheck --noweb\n"
						   + "awk '$5*$6>5 {print $2}' mafCheck.frq > polymorphicSNPs.txt\n"
						   + "plink --bfile "+analysisGenoDirAndNameRoot+" --keep ../"+phenoFiles[l]+" --extract polymorphicSNPs.txt --transpose --recode12 --output-missing-genotype 0 --out forAnalysis --noweb\n"
						   + ext.parseDirectoryOfFile(commandFullPath) + KINSHIP_PARAMETERS[i] + " forKinship\n"
						   + commandFullPath
						   		+ " -t forAnalysis" 
						   		+ " -k forKinship"+KINSHIP_EXTENSIONS[i]
				   				+ " -p " + phenoCovDir + phenoFiles[l]
		   						+ (covFiles[l] == null?"":" -c " + phenoCovDir + covFiles[l])
		   						+ " -o " + outDir + phenoLabel);
			}
		}

		batchDir = ext.verifyDirFormat(batchDir);
		if (!Files.exists(batchDir)) {
			new File(batchDir).mkdirs();
		}
		// TODO consolidate with new method
		Files.qsub(batchDir+"runEmmax", batchDir, -1, "[%0]", Matrix.toMatrix(Array.toStringArray(scripts)), qsubMemInMBs, qsubWalltimeInHours);
//		System.out.println("scripts.size(): " + scripts.size() + "\nbatchDir: " + batchDir);
		Files.qsubMultiple(batchDir+"chunkSB_emmax", Array.stringArraySequence(scripts.size(), batchDir+"runEmmax_", ".qsub"), 16, -1, qsubMemInMBs, qsubWalltimeInHours);
	}
	
	public static void parseResults(String resultDir, double pValThreshold) {
		String[] fileNames;
		String fileName;
		String label;
		String[] fileNameRootsTemp;
		int index;
		boolean found;
		Hashtable<String, String[]> modelList;
		BufferedReader reader;
		String line;
		int numSamples = 0;
		int[] numMarkers;
		Vector<Double> pVals;
		double[] pvals;
		double[] lambda;
		int[] numMarkersSig;
		Enumeration<String> keys;
		String key;
		String[] trav;
		PrintWriter writer;
		
		fileNames = Files.list(resultDir, ".log", false);
		modelList = new Hashtable<String, String[]>();
		for (int i = 0; i < fileNames.length; i++) {
			found = false;
			fileName = ext.rootOf(fileNames[i]);
			for (int j = 0; j < KINSHIP_LABELS.length; j++) {
				if(fileNames[i].contains("_" + KINSHIP_LABELS[j])) {
					found = true;
					index = fileName.indexOf("_" + KINSHIP_LABELS[j]);
					label = fileName.substring(0, index) + fileName.substring(index +  KINSHIP_LABELS[j].length() + 1);
					if (! modelList.containsKey(label)) {
						fileNameRootsTemp = new String[KINSHIP_LABELS.length];
						fileNameRootsTemp[j] = fileName;
						modelList.put(label, fileNameRootsTemp);
					} else {
						modelList.get(label)[j] = fileName;
					}
					break;
				}
			}
			if (! found) {
				modelList.put(fileName, new String[] {fileName});
			}
		}

		keys = modelList.keys();
		try {
			writer = new PrintWriter(resultDir + "resultsSummary.xln");
			line = "Pheno\t#Samples";
			for (int i = 0; i < KINSHIP_LABELS.length; i++) {
				line += "\t#Markers_" + KINSHIP_LABELS[i] + "\tLambda_" + KINSHIP_LABELS[i] + "\t#MarkersSig_" + KINSHIP_LABELS[i];
			}
			writer.println(line);
			while (keys.hasMoreElements()) {
				key = keys.nextElement();
				fileNameRootsTemp = modelList.get(key);
				numMarkersSig = new int[fileNameRootsTemp.length];
				lambda = new double[fileNameRootsTemp.length];
				numMarkers = new int[fileNameRootsTemp.length];
				for (int i = 0; i < fileNameRootsTemp.length; i++) {
					pVals = new Vector<Double>();
					if (fileNameRootsTemp[i] != null) {
							reader = new BufferedReader(new FileReader(resultDir + fileNameRootsTemp[i] + ".log"));
							while(reader.ready()) {
								line = reader.readLine();
								if (line.contains("Reading the phenotype file")) {
									while (reader.ready()) {
										line = reader.readLine();
										//TODO should consider and tolerate the scenario when there is no " rows and" here
										if (line.contains(" rows and")) {
											numSamples = Integer.parseInt(line.trim().split(" rows and")[0]);
											break;
										}
									}
								}
							}
							reader.close();
							if (new File(resultDir + fileNameRootsTemp[i] + ".ps").exists()) {
								reader = new BufferedReader(new FileReader(resultDir + fileNameRootsTemp[i] + ".ps"));
								while (reader.ready()) {
									trav = reader.readLine().split("\t");
									if(ext.isValidDouble(trav[1])) {
										numMarkers[i] ++;
										pVals.add(Double.parseDouble(trav[2]));
										if(Double.parseDouble(trav[2]) <= pValThreshold) {
											numMarkersSig[i] ++;
										}
									}
								}
								reader.close();
								pvals = new double[pVals.size()];
								for (int j = 0; j < pvals.length; j++) {
									pvals[j] = pVals.elementAt(j);
								}
								lambda[i] = ProbDist.ChiDistReverse(Array.median(pvals), 1)/ProbDist.ChiDistReverse(0.50, 1);
							}
					}
				}
				line = key + "\t" + numSamples;
				for (int i = 0; i < fileNameRootsTemp.length; i++) {
					line += "\t" + numMarkers[i] + "\t" + ext.formDeci(lambda[i],4) + "\t" + numMarkersSig[i];
				}
				writer.println(line);
			}
			writer.close();
			System.out.println("Results summary is ready at: " + resultDir + "resultsSummary.xln");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String commandFullPath;
		String phenoAndCovDir;
		String phenoNameExtOrFullPath;
		String kinshipGenoDirAndNameRoot, analysisGenoDirAndNameRoot;
//		String genoNameExtOrFullPath;
		String covNameExtOrFullPath;
		String outDir;
		String batchDir = "";
		int qsubMemInMBs;
		double qsubWalltimeInHours;
		Logger log;
		boolean isExisted;
		boolean isSucceeded;
		boolean isParseResults;
		double psig;

//		parseResults("C:/projects/BOSS/toDelete/dichotomic/", .000001);
//		parseResults("C:/projects/BOSS/toDelete/ver4_1391samples_intercept/", .000001);
//		System.exit(0);
		
		isParseResults = false;
		commandFullPath = "/home/pankrat2/shared/bin/emmax-beta-07Mar2010/emmax -v -d 10";
		phenoAndCovDir = "/home/pankrat2/shared/boss/emmax/phenoCovs4/";
//		phenoAndCovDir = "D:/boss/test/";
		phenoNameExtOrFullPath = "_pheno.txt";
		covNameExtOrFullPath = "_covars.txt";
		kinshipGenoDirAndNameRoot = "/home/pankrat2/shared/boss/BOSS_plinkData/ldPruned/plink";
		analysisGenoDirAndNameRoot = "/home/pankrat2/shared/boss/BOSS_plinkData/analyzeThis/analyzeThis";
//		genoDirAndNameRoot = "D:/boss/BOSS_plinkData/analyzeThis/analyzeThis";
//		genoNameExtOrFullPath = "plink";
		outDir = "/home/pankrat2/shared/boss/emmax/results/phenoCovs4/";
//		outDir = "D:/boss/temp/";
		batchDir = "/home/pankrat2/shared/boss/emmax/scripts/phenoCovs4/";
//		batchDir = "D:/boss/temp/";
		qsubMemInMBs = 15000;
		qsubWalltimeInHours = 12;

		psig = 0.000001;

		String usage = "\n" +
		"gwas.Emmax requires 6 - 8 arguments\n" +
		"   (1) full path of the command (i.e. command=" + commandFullPath + " (default))\n" +
		"   (2) directory of the pheno and covariate files (i.e. phenodir=" + phenoAndCovDir + " (default))\n" +
		"   (3) name extension of the pheno files, or full path of a specific pheno file (i.e. phenoext=" + phenoNameExtOrFullPath + " (default))\n" +
		"   (4) name extension of the covariate files, or full path of a specific covariate file (i.e. covext=" + covNameExtOrFullPath + " (default))\n" +
		"   (5) directory and root of the genotype files for kinship (i.e. kinshipGenoDir=" + kinshipGenoDirAndNameRoot + " (default))\n" +
		"   (6) directory and root of the genotype files for analysis (i.e. analysisGenoDir=" + analysisGenoDirAndNameRoot + " (default))\n" +
		"   (7) directory of output files (i.e. outdir=" + outDir + " (default))\n" +
		"   (8) directory of batch files (i.e. batchdir=" + batchDir + " (default))\n" +
		"   (9) (optional) qsub memory size (i.e. qsubmem=" + qsubMemInMBs + " (default; in megabytes))\n" +
		"   (10) (optional) qsub walltime (i.e. qsubwalltime=" + qsubWalltimeInHours + " (default; in hours))\n" +
		"Or\n" +
		"   (1) to parse results (i.e. parsereslt=" + isParseResults + " (default))\n" +
		"   (2) directory of the results to parse (i.e. outdir=" + outDir + " (default))\n" +
		"   (3) threshold for significant p-value (i.e. psig=" + psig + " (default))\n" +
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("command=")) {
				commandFullPath = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("phenodir=")) {
				phenoAndCovDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("phenoext=")) {
				phenoNameExtOrFullPath = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("covext=")) {
				covNameExtOrFullPath = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("kinshipGenoDir=")) {
				kinshipGenoDirAndNameRoot = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("analysisGenoDir=")) {
				analysisGenoDirAndNameRoot = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outdir=")) {
				outDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("batchdir=")) {
				batchDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("qsubmem=")) {
				qsubMemInMBs = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("qsubwalltime=")) {
				qsubWalltimeInHours = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("parseresult=")) {
				isParseResults = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("psig=")) {
				psig = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		isExisted = false;
		isSucceeded = true;
		try {
			if(!new File(outDir).exists()) {
				if (! new File(outDir).mkdir()) {
					log = new Logger();
					isSucceeded = false;
				} else {
					log = new Logger(outDir + "Emmax_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
					isExisted = true;
				}
			} else {
				log = new Logger(outDir + "Emmax_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
			}

			if (isParseResults) {
				parseResults(outDir, psig);
			} else {
				batchDir = ext.verifyDirFormat(batchDir);
				log.report("Genvisis (R) 2014. \nEmmax analysis "
							+ (new SimpleDateFormat("MM/dd/yyyy HH:mm:ss").format(new Date()))
							+ "\n-Kinship genotype files directory and name root: " +  kinshipGenoDirAndNameRoot
							+ "\n-Analysis genotype files directory and name root: " +  analysisGenoDirAndNameRoot
//							+ "\n-Geno file name extension or full path: " + genoNameExtOrFullPath
							+ "\n-Pheno and covariate files directory: " + phenoAndCovDir
							+ "\n-Pheno file name extension or full path: " + phenoNameExtOrFullPath
							+ "\n-Covariate file name extension or full path: " + covNameExtOrFullPath
							+ "\n-Output files directory: " + outDir
							+ "\n-Batch files directory: " + batchDir);
	
				if (isExisted) {
					log.reportError("Warning --- Directory " + outDir + " already exists. Existing files might be reused or overwritten.");
				} else if (! isSucceeded) {
					log.reportError("Warning --- Cannot create the directory " + outDir);
				} else {
					log.report("Creating result directory " + outDir);
				}
	
//				jobNamesWithAbsolutePaths = generateScripts(commandFullPath, genoDir, genoNameExtOrFullPath, phenoAndCovDir, phenoNameExtOrFullPath, covNameExtOrFullPath, outDir, batchDir, log);
				generateScripts(commandFullPath, kinshipGenoDirAndNameRoot, analysisGenoDirAndNameRoot, phenoAndCovDir, phenoNameExtOrFullPath, covNameExtOrFullPath, outDir, batchDir, qsubMemInMBs, qsubWalltimeInHours, log);
				log.report("\nEmmax scripts are finished.");
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
