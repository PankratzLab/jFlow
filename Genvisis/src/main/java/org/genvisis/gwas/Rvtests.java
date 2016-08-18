package org.genvisis.gwas;

import java.io.File;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;

public class Rvtests {
	public static void generateRvtestsScript(String dirOfRvtestsCommands, String dirOfPedFiles, String fullpathToVcf, String scriptsDir, String resultsDir, Logger log) {
		String[] files;
		String command = null;
		String[][] iterations;

		if (! scriptsDir.endsWith("/")) {
			scriptsDir += "/";
			log.report("Warning - a '/' symbol has been added to the following entry " + scriptsDir);
		}
		if (! resultsDir.endsWith("/")) {
			resultsDir += "/";
			log.report("Warning - a '/' symbol has been added to the following entry " + resultsDir);
		}

		files = Files.list(dirOfPedFiles, ".ped", false);
		iterations = new String[files.length][2];
		for (int i = 0; i < files.length; i++) {
			iterations[i][0] = files[i];
			iterations[i][1] = files[i].split("_pheno")[0];
		}

		command = dirOfRvtestsCommands + "vcf2kinship"
										+ " --ped " + dirOfPedFiles + "[%0]"
										+ " --inVcf " + fullpathToVcf
										+ " --bn"
										+ " --out " + resultsDir + "[%1]"
										+ " --xLabel X --xHemi"
										+ "\n"
										+ dirOfRvtestsCommands + "rvtest"
										+ " --pheno " + dirOfPedFiles + "[%0]"
										+ " --inVcf " + fullpathToVcf
										+ " --out " + resultsDir + "[%1]"
										+ " --kinship " + resultsDir + "[%1].kinship"
										+ " --xLabel X --meta score,cov,dominant,recessive";

//		Files.qsub("rvtests", command, Matrix.toMatrix(files), -1, 6, -1);
		Files.qsub(scriptsDir + "[%1]", command, iterations, -1, 6, -1);
	}

	public static void renameToNamingScheme(String resultsDir, String[] filenameConversionDictionary, String[] segmentIndices, Logger log) {
		String[] filenames, tmp1;
		String newFilename, tmp2;
		Vector<String> filenameAdditionalSegments;
		Vector<String[]> conversionTable;
		int[] indices;
		int index;

		if (! resultsDir.endsWith("/")) {
			resultsDir += "/";
			log.report("Warning - a '/' symbol has been added to the following entry " + resultsDir);
		}

		index = -1;
		filenameAdditionalSegments = new Vector<String>();
		indices = new int[segmentIndices.length];
		for (int i = 0; i < indices.length; i++) {
			try {
				indices[i] = Integer.parseInt(segmentIndices[i]);
			} catch (NumberFormatException e) {
				indices[i] = index;
				filenameAdditionalSegments.add(segmentIndices[i]);
				index --;
			}
		}

		conversionTable = new Vector<String[]> ();
		for (int i = 0; i < filenameConversionDictionary.length; i++) {
			tmp1 = filenameConversionDictionary[i].split("=");
			if (tmp1.length == 2) {
				conversionTable.add(tmp1);
			}
		}
		for (int i = 0; i < conversionTable.size(); i++) {
			for (int j = 0; j < conversionTable.size(); j++) {
				if (i != j && conversionTable.elementAt(i)[1].contains(conversionTable.elementAt(j)[0])) {
					log.reportError("Error - chained conversion: " + conversionTable.elementAt(i)[0] + "=" + conversionTable.elementAt(i)[1] + ";\t" + conversionTable.elementAt(j)[0] + "=" + conversionTable.elementAt(j)[1]);
					System.exit(0);
				}
			}
		}

		filenames = Files.list(resultsDir, null, false);
		for (int i = 0; i < filenames.length; i++) {
			newFilename = filenames[i];
			for (int j = 0; j < conversionTable.size(); j++) {
				tmp1 = conversionTable.elementAt(j);
				newFilename = newFilename.replaceAll(tmp1[0], tmp1[1]);
			}
			tmp1 = newFilename.substring(0, newFilename.indexOf(".")).split("_");
			tmp2 = newFilename.substring(newFilename.indexOf("."));
			newFilename = "";
			for (int j = 0; j < indices.length; j++) {
				if (indices[j] < 0) {
					newFilename += ("_" + filenameAdditionalSegments.elementAt(-indices[j] - 1).replaceAll("'", ""));
				} else {
					newFilename += ("_" + tmp1[indices[j]]);
				}
			}
			newFilename = newFilename.substring(1) + tmp2;
			new File(resultsDir, filenames[i]).renameTo(new File(resultsDir, newFilename));
		}
	}

	/**
	 * This is a program specifically for 
	 * @param args
	 */
	public static void main(String[] args) {
		int numArgs = args.length;
		boolean isScripts, isRename;
		String dirOfRvtestsCommands, dirOfPedFiles, fullpathToVcf, rScriptDir, resultsDir, filenameConversionTable, filenameSegmentIndices;
		String commandIsScript, commandRvExecutable, commandPedDir, commandVcf, commandResultsDir, commandScriptDir, commandIsRename;
		Logger log;

		commandIsScript = "-script";
		commandRvExecutable = "exec=";
		commandPedDir = "peddir=";
		commandVcf = "vcf=";
		commandResultsDir = "outdir=";
		commandScriptDir = "scriptdir=";
		commandIsRename = "-rename";
		isScripts = false;
		dirOfRvtestsCommands = "/home/pankrat2/shared/skatMeta/giant_hematology/rvtests/executable/";
		dirOfPedFiles = "/home/pankrat2/shared/skatMeta/giant_hematology/pedigrees/";
		fullpathToVcf = "/home/pankrat2/shared/skatMeta/giant_hematology/vcf_files/CARDIA_AA_CHR1_23_STRANDALIGNED_FINAL_nochr.vcf.gz";
		rScriptDir = "/home/pankrat2/shared/skatMeta/giant_hematology/scripts/";
		resultsDir = "/home/pankrat2/shared/skatMeta/giant_hematology/results/";
		
		isRename = false;
		filenameConversionTable = "RBC,HCT,HGB,MCV,MCH,MCHC,RDW,WBC_TOTAL=WBC,WBC_BASO=BAS,WBC_EOS=EOS,WBC_LYMPH=LYM,WBC_MONO=MON,WBC_NEUTRO=NEU,PLT,MPV";
		filenameSegmentIndices = "0,2,COHORT,1,'121114',np";

		String usage = "\nTo generate scripts for all the .ped files in a single directory:"
					+ "\n   (1) command for generating scripts (i.e. " + commandIsScript + " (default))"
					+ "\n   (2) directory of the RvTests executable files (i.e. " + commandRvExecutable + dirOfRvtestsCommands + " (default))"
					+ "\n   (3) directory of the .ped files (i.e. " + commandPedDir + dirOfPedFiles + " (default))"
					+ "\n   (4) directory of the condition files (i.e. " + commandVcf + fullpathToVcf + " (default))"
					+ "\n   (5) directory to output the RvTests results files (i.e. " + commandResultsDir + resultsDir + " (default))"
					+ "\n   (6) directory to output the scripts (i.e. " + commandScriptDir + rScriptDir + " (default))"
					+ "\n"
					+ "\nTo rename the output files per naming scheme:"
					+ "\n   (1) command for renaming file (i.e. " + commandIsScript + " (default))"
					+ "\n   (2) directory to output the RvTests results files (i.e. " + commandResultsDir + resultsDir + " (default))"
					+ "\n   (3) indices of file name segements (i.e. " + resultsDir + " (default))"
					+ "\n"
					+ "";

//		resultsDir = "D:/RvTests/";
		resultsDir = "/home/pankrat2/shared/skatMeta/giant_hematology/results/";
		filenameConversionTable = "RBC,HCT,HGB,MCV,MCH,MCHC,RDW,WBC_TOTAL=WBC,WBC_BASO=BAS,WBC_EOS=EOS,WBC_LYMPH=LYM,WBC_MONO=MON,WBC_NEUTRO=NEU,PLT,MPV";
		filenameSegmentIndices = "0,2,COHORT,1,'121114',np";
		isRename = true;

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith(commandIsScript)) {
				isScripts = true;
				numArgs--;
			} else if (args[i].startsWith(commandRvExecutable)) {
				dirOfRvtestsCommands = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandPedDir)) {
				dirOfPedFiles = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandVcf)) {
				fullpathToVcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandResultsDir)) {
				resultsDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandScriptDir)) {
				rScriptDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commandIsRename)) {
				isRename = true;
				numArgs--;
			} else if (args[i].startsWith(filenameSegmentIndices)) {
				filenameSegmentIndices = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}

		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		if (isScripts) {
			log = new Logger();
			generateRvtestsScript(dirOfRvtestsCommands, dirOfPedFiles, fullpathToVcf, rScriptDir, resultsDir, log);
		} else if (isRename) {
			log = new Logger();
			renameToNamingScheme(resultsDir, filenameConversionTable.split(","), filenameSegmentIndices.split(","), log);
		} else {
			log = new Logger();
			log.reportError("No command executed, due to none of the following is specified: " + commandIsScript);
		}

		log.report("Program completed.");
	}

}
