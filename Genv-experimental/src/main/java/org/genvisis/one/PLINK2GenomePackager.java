package org.genvisis.one;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.util.HashMap;
import java.util.Hashtable;

import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class PLINK2GenomePackager {

	// plink2 --threads 8 --dosage fullList.txt list format=1 Zout --fam gedi_exome_plink.fam --covar
	// GEDI_covars.dat --pheno GEDI_pheno_mtPC0_exome.dat

	int QSUB_RAM = 10000;
	int QSUB_HRS = 4;
	int QSUB_PRC = 8;

	void setup(String dir, String fileList, String pmAllFile, String qsubQueue) {

		File dirFile = new File(dir);

		if (!dirFile.exists()) {
			System.err.println("Error - run directory { " + dir + " } does not exist");
			return;
		}

		String[] covars = dirFile.list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith("_covars.dat");
			}
		});
		if (covars.length == 0) {
			System.err.println(ext.getTime()
												 + "]\tWARNING - no covariate files found.  If covar files exist, ensure they are named *_covars.dat and try again.  Otherwise, skipping covariate analysis.");
			// return;
		}


		HashMap<String, String> covarPrefices = new HashMap<String, String>();

		for (String covarFile : covars) {
			String[] parts = covarFile.split("_");
			if (parts.length == 1) {
				System.err.println(ext.getTime() + "]\tError - covariate file {" + covarFile
													 + "} is named incorrectly.  Correct naming is *_covars.dat for covariate files.");
				continue;
			} else if (parts.length == 2) {
				covarPrefices.put(parts[0], covarFile);
			} else if (parts.length > 2) {
				covarPrefices.put(covarFile.substring(0, covarFile.length() - "_covars.dat".length()),
													covarFile);
			}
		}

		String[] phenos = dirFile.list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith("_pheno.dat");
			}
		});
		if (phenos.length == 0) {
			System.err.println(ext.getTime()
												 + "]\tError - no phenotype files found.  Ensure that files are named *_pheno.dat and try again.");
			return;
		}

		String[] famList = dirFile.list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(".fam");
			}
		});
		if (famList.length == 0) {
			System.err.println(ext.getTime() + "]\tError - no .fam file found.");
			return;
		} else if (famList.length > 1) {
			System.err.println(ext.getTime() + "]\tError - multiple .fam files found.");
			return;
		}

		String famFile = famList[0];

		HashMap<String, String> plinkRuns = new HashMap<String, String>();

		for (String pheno : phenos) {
			String[] parts = pheno.split("_");
			if (parts.length == 1) {
				System.err.println(ext.getTime() + "]\tError - phenotype file {" + pheno
													 + "} is named incorrectly.  Correct naming is *_pheno.dat for phenotype files.");
				continue;
			} else if (parts.length == 2) {
				// if (covarPrefices.containsKey(parts[0])) {
				if (setupForScript(dir, fileList, famFile, covarPrefices.get(parts[0]), pheno)) {
					String phenoPrefix = pheno.substring(0, pheno.length() - "_pheno.dat".length());
					String phenoDir = dir + phenoPrefix + "/";
					String script = createScript(fileList, famFile, covarPrefices.get(parts[0]), pheno,
																			 phenoDir);
					plinkRuns.put(phenoDir, script);
				}
				// } else {
				// System.err.println(ext.getTime() + "]\tError - no covariate file found for phenotype file
				// {" + pheno + "}.");
				// continue;
				// }
			} else if (parts.length > 2) {
				String phenoCovarCheck = pheno.substring(0, pheno.length() - "_pheno.dat".length());
				String covar = covarPrefices.get(phenoCovarCheck);
				if (covar == null) {
					covar = covarPrefices.get(parts[0]);
				}
				// if (covarPrefices.containsKey(phenoCovarCheck)) {
				if (setupForScript(dir, fileList, famFile, covar/* Prefices.get(phenoCovarCheck) */,
													 pheno)) {
					String phenoPrefix = pheno.substring(0, pheno.length() - "_pheno.dat".length());
					String phenoDir = dir + phenoPrefix + "/";
					String script = createScript(fileList, famFile, covar/* Prefices.get(phenoCovarCheck) */,
																			 pheno, phenoDir);
					plinkRuns.put(phenoDir, script);
				}
				// } else if (covarPrefices.containsKey(parts[0])) {
				// if (setupForScript(dir, fileList, famFile, covarPrefices.get(parts[0]), pheno)) {
				// String phenoPrefix = pheno.substring(0, pheno.length() - "_pheno.dat".length());
				// String phenoDir = dir + phenoPrefix + "/";
				// String script = createScript(fileList, famFile, covarPrefices.get(parts[0]), pheno,
				// phenoDir);
				// plinkRuns.put(phenoDir, script);
				// }
				// } else {
				// System.err.println(ext.getTime() + "]\tError - no covariate file found for phenotype file
				// {" + pheno + "}.");
				// continue;
				// }
			}
		}

		if (plinkRuns.size() == 0) {
			System.err.println(ext.getTime() + "]\tError - no valid PLINK2 scripts found.");
			return;
		}

		for (java.util.Map.Entry<String, String> plinkRun : plinkRuns.entrySet()) {
			Files.qsub(plinkRun.getKey() + "runPlink2.qsub", plinkRun.getValue(), QSUB_RAM, QSUB_HRS,
								 QSUB_PRC);
		}

		StringBuilder sb = new StringBuilder("cd ").append(dir).append("\n");
		for (String plinkDir : plinkRuns.keySet()) {
			sb.append("cd ").append(plinkDir).append("\n").append("qsub");
			if (qsubQueue != null) {
				sb.append(" -q ").append(qsubQueue);
			}
			sb.append(" runPlink2.qsub").append("\n").append("cd ..").append("\n");
		}
		Files.write(sb.toString(), dir + "masterRun.sh");
		Files.chmod(dir + "masterRun.sh");

		setupProcess(dir, plinkRuns, pmAllFile, qsubQueue);
	}

	boolean setupForScript(String dir, String fileList, String famFile, String covarFile,
												 String phenoFile) {

		String phenoPrefix = phenoFile.substring(0, phenoFile.length() - "_pheno.dat".length());

		String phenoDir = dir + phenoPrefix + "/";
		File phenoDirFile = new File(phenoDir);
		if (phenoDirFile.exists()) {
			return false;
		}
		boolean createdPhenoDir;
		createdPhenoDir = phenoDirFile.mkdir();
		if (!createdPhenoDir) {
			System.err.println(ext.getTime() + "]\tError - couldn't create directory {" + phenoDir
												 + "}.");
			return false;
		}
		boolean fileCopy;
		if (covarFile != null) {
			fileCopy = Files.copyFile(dir + covarFile, phenoDir + covarFile);
			if (!fileCopy) {
				System.err.println(ext.getTime()
													 + "]\tError - couldn't copy covariate file into phenotype subdirectory {"
													 + phenoDir + "}.");
				return false;
			}
		}
		fileCopy = Files.copyFile(dir + phenoFile, phenoDir + phenoFile);
		if (!fileCopy) {
			System.err.println(ext.getTime()
												 + "]\tError - couldn't copy phenotype file into phenotype subdirectory {"
												 + phenoDir + "}.");
			return false;
		}
		fileCopy = Files.copyFile(dir + famFile, phenoDir + famFile);
		if (!fileCopy) {
			System.err.println(ext.getTime()
												 + "]\tError - couldn't copy .fam file into phenotype subdirectory {"
												 + phenoDir + "}.");
			return false;
		}

		return true;
	}

	String createScript(String fileList, String famFile, String covarFile, String phenoFile,
											String phenoDir) {
		String script = "cd " + phenoDir + "\nplink2 --threads " + QSUB_PRC + " --dosage ../" + fileList
										+ " list format=1 Zout --fam " + famFile
										+ (covarFile != null ? " --covar " + covarFile : "") + " --pheno " + phenoFile;
		return script;
	}

	void setupProcess(String dir, HashMap<String, String> plinkRuns, String pmAllFile,
										String qsubQueue) {
		StringBuilder masterProc = new StringBuilder();
		for (String plinkDir : plinkRuns.keySet()) {
			StringBuilder procCmd = new StringBuilder();
			procCmd.append("cd " + plinkDir + "\n").append(Files.getRunString())
						 .append(" one.PLINK2GenomePackager -process dir=").append(plinkDir).append(" pmFile=")
						 .append(pmAllFile)
						 .append(" N=`grep -o -E '[0-9]+ people pass filters and QC' *.o | sed 's/.*://g' | sed 's/[^0-9]*//g'`")
						 .append("\n");
			Files.qsub(plinkDir + "process.qsub", procCmd.toString(), 5000, 3, 1);
			masterProc.append("cd ").append(plinkDir).append("\n");
			masterProc.append("qsub");
			if (qsubQueue != null) {
				masterProc.append(" -q ").append(qsubQueue);
			}
			masterProc.append(" process.qsub").append("\n");
		}
		Files.write(masterProc.toString(), dir + "masterProcess.sh");
		Files.chmod(dir + "masterProcess.sh");
	}

	void process(String phenoDir, String pmAllFile, String N) {
		if (N == null || "".equals(N) || !ext.isValidInteger(N)) {
			System.err.println(ext.getTime() + "]\tError - specified value of N {" + N + "} is invalid.");
			return; // continue without N?
		}

		if (!(new File(pmAllFile)).exists()) {
			System.err.println(ext.getTime() + "]\tError - specified PM_ALL file {" + pmAllFile
												 + "} doesn't exist.");
			return;
		}

		String file = phenoDir + "plink.assoc.dosage.gz";
		String newHeader = "SNP\tCHR\tPOS\tA1\tA2\tN\tFRQ\tINFO\tBETA\tSE\tP";


		System.out.println(ext.getTime() + "]\tLoading PM_ALL file...");
		Hashtable<String, String> pmAll = HashVec.loadFileToHashString(pmAllFile, 1, new int[] {0, 3},
																																	 "\t", false, false);
		System.out.println(ext.getTime() + "]\tPM_ALL file successfully loaded.");

		System.out.println(ext.getTime() + "]\tProcessing results file...");
		String temp = phenoDir.substring(0, phenoDir.length() - 1);
		temp = temp.substring(temp.lastIndexOf("/") + 1, temp.length());
		String outFile = phenoDir + temp + ".plink.assoc.dosage.gz";
		try {
			BufferedReader reader = Files.getAppropriateReader(file);
			System.out.println(ext.getTime() + "]\tWriting results to {" + outFile + "}...");
			PrintWriter writer = Files.getAppropriateWriter(outFile);
			writer.println(newHeader);
			String line = null;
			reader.readLine();
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				String[] parts = line.split("[\\s]+");
				String mkrChrPos = pmAll.get(parts[0]);

				StringBuilder sb = new StringBuilder();
				sb.append(parts[0]).append("\t");
				if (mkrChrPos != null) {
					sb.append(mkrChrPos).append("\t");
				} else {
					sb.append(".\t.\t");
				}
				sb.append(parts[1]).append("\t").append(parts[2]).append("\t");
				sb.append(N);
				for (int i = 3; i < parts.length; i++) {
					sb.append("\t").append(parts[i]);
				}
				writer.println(sb.toString());
			}
			writer.flush();
			writer.close();
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		System.out.println(ext.getTime() + "]\tProcessing Complete!");

		double pvalThresh = 0.001;
		boolean gc = true;

		System.out.println(ext.getTime() + "]\tRunning METAL...");
		StringBuilder metalCRF = new StringBuilder("metal\n").append(temp).append("\n")
																												 .append("build=37\n")
																												 .append("hits_p<=" + pvalThresh + "\n")
																												 .append("genomic_control="
																																 + (gc ? "TRUE" : "FALSE") + "\n")
																												 .append(temp + ".plink.assoc.dosage.gz");
		String metalName = "metal_" + temp + ".crf";
		String metalDir = (new File(phenoDir)).getAbsolutePath();
		Files.write(metalCRF.toString(), ext.verifyDirFormat(metalDir) + metalName);
		String path = PLINK2GenomePackager.class.getProtectionDomain().getCodeSource().getLocation()
																						.getPath();
		String decodedPath = path;
		try {
			decodedPath = URLDecoder.decode(path, "UTF-8");
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		CmdLine.run(new String[] {"java", "-jar", decodedPath, "Launch", metalName}, metalDir,
								System.out, System.err, new Logger(), false);
		System.out.println(ext.getTime() + "]\tMETAL Complete!");
	}



	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = ext.pwd();
		String fileList = "fullList.txt";
		boolean process = false;
		String N = null;
		String pmAll = "D:/test/PM_all";
		String qsub = null;

		String usage = "\n" + "one.PLINK2GenomePackager requires 3-4 arguments\n"
									 + "   (1) directory (i.e. dir=" + dir + " (default))\n"
									 + "   (2) file list of source data files, indexed (i.e. file=" + fileList
									 + " (default))\n"
									 + "   (3) location of PM_ALL file (i.e. pmFile=PM_ALL (not the default))\n"
									 + " OR:\n" + "   (1) directory of PLINK2 run (i.e. dir=" + dir
									 + " (not the default))\n"
									 + "   (2) location of PM_ALL file (i.e. pmFile=PM_ALL (not the default))\n"
									 + "   (3) Number of individuals in analysis (i.e. N=1000 (not the default))\n"
									 + "   (4) -process flag" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("file=")) {
				fileList = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("pmFile=")) {
				pmAll = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("qsub=")) {
				qsub = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("N=")) {
				N = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-process")) {
				process = true;
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (process) {
				new PLINK2GenomePackager().process(dir, pmAll, N);
			} else {
				new PLINK2GenomePackager().setup(dir, fileList, pmAll, qsub);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
