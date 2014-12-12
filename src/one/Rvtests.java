package one;

import java.text.SimpleDateFormat;
import java.util.Date;
import common.Files;
import common.Logger;

public class Rvtests {
	public static void generateRvtestsScript(String dirOfRvtestsCommands, String dirOfPedFiles, String fullpathToVcf, String scriptDir, String resultsDir) {
		String[] files;
		String command = null;
		String[][] iterations;

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
		Files.qsub(scriptDir + "[%1]", command, iterations, -1, 6, -1);
	}


	/**
	 * This is a program specifically for 
	 * @param args
	 */
	public static void main(String[] args) {
		int numArgs = args.length;
		boolean isScripts;
		String dirOfRvtestsCommands, dirOfPedFiles, fullpathToVcf, rScriptDir, resultsDir;
		String commandIsScript, commandRvExecutable, commandPedDir, commandVcf, commandResultsDir, commandScriptDir;
		Logger log;

		commandIsScript = "--script";
		commandRvExecutable = "exec=";
		commandPedDir = "peddir=";
		commandVcf = "vcf=";
		commandResultsDir = "outdir=";
		commandScriptDir = "scriptdir=";

		isScripts = false;
		dirOfRvtestsCommands = "/home/pankrat2/shared/skatMeta/giant_hematology/rvtests/executable/";
		dirOfPedFiles = "/home/pankrat2/shared/skatMeta/giant_hematology/pedigrees/";
		fullpathToVcf = "/home/pankrat2/shared/skatMeta/giant_hematology/vcf_files/CARDIA_AA_CHR1_23_STRANDALIGNED_FINAL_nochr.vcf.gz";
		rScriptDir = "/home/pankrat2/shared/skatMeta/giant_hematology/scripts/";
		resultsDir = "/home/pankrat2/shared/skatMeta/giant_hematology/results";
		isScripts = true;

		String usage = "\nTo generate scripts for all the .ped files in a single directory:"
					+ "\n   (1) command for generating scripts (i.e. " + commandIsScript + " (default))"
					+ "\n   (2) directory of the RvTests executable files (i.e. " + commandRvExecutable + dirOfRvtestsCommands + " (default))"
					+ "\n   (3) directory of the .ped files (i.e. " + commandPedDir + dirOfPedFiles + " (default))"
					+ "\n   (4) directory of the condition files (i.e. " + commandVcf + fullpathToVcf + " (default))"
					+ "\n   (5) directory to output the RvTests results files (i.e. " + commandResultsDir + resultsDir + " (default))"
					+ "\n   (6) directory to output the scripts (i.e. " + commandScriptDir + rScriptDir + " (default))"
					+ "\n"
					+ "";

		dirOfRvtestsCommands = "/home/pankrat2/shared/skatMeta/giant_hematology/rvtests/executable/";
		dirOfPedFiles = "/home/pankrat2/shared/skatMeta/giant_hematology/pedigrees/";
		fullpathToVcf = "/home/pankrat2/shared/skatMeta/giant_hematology/vcf_files/CARDIA_AA_CHR1_23_STRANDALIGNED_FINAL_nochr.vcf.gz";
		rScriptDir = "/home/pankrat2/shared/skatMeta/giant_hematology/scripts/";
		resultsDir = "/home/pankrat2/shared/skatMeta/giant_hematology/results";
		isScripts = true;

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
			generateRvtestsScript(dirOfRvtestsCommands, dirOfPedFiles, fullpathToVcf, rScriptDir, resultsDir);
		} else {
			log = new Logger();
			log.reportError("No command executed, due to none of the following is specified: " + commandIsScript);
		}

		log.report("Program completed.");
	}

}
