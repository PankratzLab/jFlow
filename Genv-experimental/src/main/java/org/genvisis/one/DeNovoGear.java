package org.genvisis.one;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Scanner;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class DeNovoGear {

	public static void generateScripts(	String pedigreeFileFullPath, String bamFilesDir,
																			String bcfFilesDir, String denovogearResultDir,
																			String scriptFileDir, String qsubLogsDir,
																			String bamFilePrefix, String bamFileSuffix, Logger log) {
		Scanner reader;
		PrintWriter writer;
		String[] line;
		// String[] qsubFiles;
		Vector<String> qsubFilesVec;
		// Vector <String[]> iterationsVec;
		// String[][] iterations;
		String command;
		// Vector<String> jobNamesWithAbsolutePaths;

		// command = "samtools mpileup -gDf /home/pankrat2/shared/bin/ref/hg19_canonical.fa " +
		// bamFilesDir + bamFilePrefix + "[%0].bam " + bamFilesDir + bamFilePrefix + "[%1].bam " +
		// bamFilesDir + bamFilePrefix + "[%2].bam > " + bcfFilesDir + "[%0].bcf\n"
		// + "~/bin/denovogear/build/src/denovogear dnm auto --ped " + bcfFilesDir + "[%0].ped --bcf " +
		// bcfFilesDir + "[%0].bcf > " + denovogearResultDir + "[%0].txt";
		// iterationsVec = new Vector<String[]>();
		qsubFilesVec = new Vector<String>();
		try {
			reader = new Scanner(new File(pedigreeFileFullPath));
			while (reader.hasNext()) {
				line = reader.nextLine().split("\t");
				if (!(new File(bamFilesDir + bamFilePrefix + line[1] + ".bam").exists()
								|| new File(bamFilesDir + bamFilePrefix + line[2] + ".bam").exists()
							|| new File(bamFilesDir + bamFilePrefix + line[3] + ".bam").exists())) {
					log.reportError("Not all bam files are found for the trio "	+ line[1] + "; " + line[2]
													+ "; " + line[3] + "\nSystem quit with error.");
					System.exit(1);
				} else if (!new File(denovogearResultDir + line[1] + ".txt").exists()) {
					writer = new PrintWriter(new FileOutputStream(bcfFilesDir + line[1] + ".ped"));
					if (bamFileSuffix != null) {
						writer.println(line[0]	+ "\t" + line[0] + "C\t" + line[0] + "D\t" + line[0]
														+ "M\t-1\t-1");
					} else {
						writer.println(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t-1\t-1");
					}
					writer.close();
					command = "cd " + qsubLogsDir;
					if (!new File(bcfFilesDir + line[1] + ".bcf").exists()) {
						// command += "\nsamtools mpileup -gDf /home/pankrat2/shared/bin/ref/hg19_canonical.fa "
						// + bamFilesDir + bamFilePrefix + line[1] + ".bam " + bamFilesDir + bamFilePrefix +
						// line[2] + ".bam " + bamFilesDir + bamFilePrefix + line[3] + ".bam > " + bcfFilesDir +
						// line[1] + ".bcf";
						command += "\nsamtools mpileup -gDf /home/spectorl/xuz2/hg19_canonical.fa "
													+ bamFilesDir + bamFilePrefix + line[1] + ".bam " + bamFilesDir
												+ bamFilePrefix + line[2] + ".bam " + bamFilesDir + bamFilePrefix + line[3]
												+ ".bam > " + bcfFilesDir + line[1] + ".bcf";
					}
					command += "\n~/bin/denovogear/build/src/denovogear dnm auto --ped "	+ bcfFilesDir
											+ line[1] + ".ped --bcf " + bcfFilesDir + line[1] + ".bcf > "
											+ denovogearResultDir + line[1] + ".txt";
					Files.qsub(scriptFileDir + line[1] + ".qsub", command, 2000, 24, 1);
					qsubFilesVec.add(scriptFileDir + line[1] + ".qsub");
				}
			}
			reader.close();

			writer = new PrintWriter(new FileOutputStream(scriptFileDir + "master_run"));
			writer.println("cd " + qsubLogsDir);
			for (int i = 0; i < qsubFilesVec.size(); i++) {
				writer.println("qsub " + qsubFilesVec.elementAt(i));
			}
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		// iterations = new String[iterationsVec.size()][3];
		// for (int i = 0; i < iterations.length; i++) {
		// iterations[i] = iterationsVec.elementAt(i);
		// }

		// Files.qsub(scriptFileDir + "trio_", qsubLogsDir, iterations.length, command, iterations,
		// 2000, 24);
		// qsubFiles = Files.list(scriptFileDir, "trio", ".qsub", false, false);
		// jobNamesWithAbsolutePaths = new Vector<String>(qsubFiles.length);
		// for (int i = 0; i < qsubFiles.length; i++) {
		// jobNamesWithAbsolutePaths.add(scriptFileDir+qsubFiles[i]);
		// }
		// Files.qsubMultiple(jobNamesWithAbsolutePaths, null, scriptFileDir, "qsubmultiple_",
		// qsubFiles.length<16?qsubFiles.length:16, true, null, -1, 22000, 24);
		Files.qsubMultiple(	qsubFilesVec, null, scriptFileDir, "qsubmultiple_",
												qsubFilesVec.size() < 4 ? qsubFilesVec.size() : 4, true, null, -1, 22000,
												24);
	}

	public static void main(String[] args) {
		String pedigreeFileFullPath, bamFilesDir, bcfFilesDir, denovogearResultDir, scriptFileDir,
				qsubLogsDir, bamFilePrefix, bamFileSuffix;
		String[] commands;
		Logger log;

		// pedigreeFileFullPath = "/home/pankrat2/shared/logan/denovogear/all.ped";
		// bamFilesDir = "/home/pankrat2/shared/logan/denovogear/data_source_bam/";
		// bcfFilesDir = "/home/pankrat2/shared/logan/denovogear/data_converted_bcf/";
		// denovogearResultDir = "/home/pankrat2/shared/logan/denovogear/results/";
		// scriptFileDir = "/home/pankrat2/shared/logan/denovogear/scripts/";
		// qsubLogsDir = "/home/pankrat2/shared/logan/denovogear/logs/";
		// bamFilePrefix = "dedup_";
		pedigreeFileFullPath = "/home/spectorl/xuz2/denovo/allTrios.ped";
		bamFilesDir = "/home/spectorl/shared/exome_processing/bam/";
		bcfFilesDir = "/scratch/bcfs/";
		denovogearResultDir = "/home/spectorl/xuz2/denovo/results/";
		scriptFileDir = "/home/spectorl/xuz2/denovo/scripts/";
		qsubLogsDir = "/home/spectorl/xuz2/denovo/logs/";
		bamFilePrefix = "rrd_";
		bamFileSuffix = "_L00?";

		commands = new String[] {	"ped=", "bamdir=", "bcfdir=", "out=", "scriptdir=", "qsublogsdir=",
															"bamprefix=", "bamsuffiex="};
		// commandVariables = new String[] {pedigreeFileFullPath, bamFilesDir, bcfFilesDir,
		// denovogearResultDir};
		String usage = "\n"
											+ "cnv.analysis.DeNovoCNV analysis normally contains the following steps, split into several rounds, within each of which all the steps will be finished by just calling the program once.\n"
										+ "   (1) full path of the pedigree file (i.e. " + commands[0]
										+ pedigreeFileFullPath + " (the default));\n"
										+ "   (2) directory of the source data bam files (i.e. " + commands[1]
										+ bamFilesDir + " (the default));\n"
										+ "   (3) directory to place the interim data bcf files (i.e. " + commands[2]
										+ bcfFilesDir + " (the default));\n"
										+ "   (4) directory to place the results of DeNovoGear (i.e. " + commands[3]
										+ denovogearResultDir + " (the default));\n"
										+ "   (5) directory to place the script files generated by this program (i.e. "
										+ commands[4] + scriptFileDir + " (the default));\n"
										+ "   (6) directory to place the log files for the scripts to run (i.e. "
										+ commands[5] + qsubLogsDir + " (the default));\n"
										+ "   (7) prefix of bam files (i.e. " + commands[6] + bamFilePrefix
										+ " (the default));\n" + "   (8) suffix of bam files (i.e. " + commands[7]
										+ bamFileSuffix + " (the default));\n" + "";

		for (String arg : args) {
			if (arg.startsWith(commands[0])) {
				pedigreeFileFullPath = arg.split("=")[1];
			} else if (arg.startsWith(commands[1])) {
				bamFilesDir = arg.split("=")[1];
			} else if (arg.startsWith(commands[2])) {
				bcfFilesDir = arg.split("=")[1];
			} else if (arg.startsWith(commands[3])) {
				denovogearResultDir = arg.split("=")[1];
			} else if (arg.startsWith(commands[4])) {
				scriptFileDir = arg.split("=")[1];
			} else if (arg.startsWith(commands[5])) {
				qsubLogsDir = arg.split("=")[1];
			} else if (arg.startsWith(commands[6])) {
				bamFilePrefix = arg.split("=")[1];
			} else if (arg.startsWith(commands[7])) {
				bamFileSuffix = arg.split("=")[1];
			} else {
				System.err.println(usage);
				System.err.println("Error - invalid argument: " + arg);
				System.exit(1);
			}
		}

		log = new Logger(ext.parseDirectoryOfFile(pedigreeFileFullPath)	+ "Genvisis_"
											+ (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
		log.report(new SimpleDateFormat("MMM dd HH:mm:ss").format(new Date())
								+ ": starting generating qsub files for DeNovoGear.java");

		generateScripts(pedigreeFileFullPath, bamFilesDir, bcfFilesDir, denovogearResultDir,
										scriptFileDir, qsubLogsDir, bamFilePrefix, bamFileSuffix, log);

		log.report(new SimpleDateFormat("MMM dd HH:mm:ss").format(new Date())
								+ ": completed generating qsub files for DeNovoGear.java.");
	}

}
