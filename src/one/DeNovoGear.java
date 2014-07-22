package one;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Scanner;
import java.util.Vector;

import common.Files;
import common.Logger;
import common.ext;

public class DeNovoGear {

	public static void generateScripts(String pedegreeFileFullPath, String bamFilesDir, String bcfFilesDir, String denovogearResultDir, String scriptFileDir, String qsubLogsDir, String bamFilePrefix, Logger log) {
		Scanner reader;
		PrintWriter writer;
		String[] line;
//		String[] qsubFiles;
		Vector<String> qsubFilesVec;
//		Vector <String[]> iterationsVec;
//		String[][] iterations;
		String command;
//		Vector<String> jobNamesWithAbsolutePaths;

//		command = "samtools mpileup -gDf /home/pankrat2/shared/bin/ref/hg19_canonical.fa " + bamFilesDir + bamFilePrefix + "[%0].bam " + bamFilesDir + bamFilePrefix + "[%1].bam " + bamFilesDir + bamFilePrefix + "[%2].bam > " + bcfFilesDir + "[%0].bcf\n"
//				+ "~/bin/denovogear/build/src/denovogear dnm auto --ped " + bcfFilesDir + "[%0].ped --bcf " + bcfFilesDir + "[%0].bcf > " + denovogearResultDir + "[%0].txt";
//		iterationsVec = new Vector<String[]>();
		qsubFilesVec = new Vector<String>();
		try {
			reader = new Scanner(new File(pedegreeFileFullPath));
			while (reader.hasNext()) {
				line = reader.nextLine().split("\t");
				if (!	(  new File(bamFilesDir + bamFilePrefix + line[1] + ".bam").exists()
						|| new File(bamFilesDir + bamFilePrefix + line[2] + ".bam").exists()
						|| new File(bamFilesDir + bamFilePrefix + line[3] + ".bam").exists()
						)) {
					log.reportError("Not all bam files are found for the trio " + line[1] + "; " + line[2] + "; " + line[3] + "\nSystem quit with error.");
					System.exit(1);
				} else if (! new File(denovogearResultDir + line[1] + ".txt").exists()) {
					writer = new PrintWriter(new FileOutputStream(bcfFilesDir + line[1] + ".ped"));
					writer.println(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t-1\t-1");
					writer.close();
					command = "cd " + qsubLogsDir;
					if (! new File(bcfFilesDir + line[1] + ".bcf").exists()){
						command += "\nsamtools mpileup -gDf /home/pankrat2/shared/bin/ref/hg19_canonical.fa " + bamFilesDir + bamFilePrefix + line[1] + ".bam " + bamFilesDir + bamFilePrefix + line[2] + ".bam " + bamFilesDir + bamFilePrefix + line[3] + ".bam > " + bcfFilesDir + line[1] + ".bcf";
					}
					command += "\n~/bin/denovogear/build/src/denovogear dnm auto --ped " + bcfFilesDir + line[1] + ".ped --bcf " + bcfFilesDir + line[1] + ".bcf > " + denovogearResultDir + line[1] + ".txt";
					Files.qsub(scriptFileDir + line[1] + ".qsub", command, 2000, 24, 1);
					qsubFilesVec.add(scriptFileDir +  line[1] + ".qsub");
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

//		iterations = new String[iterationsVec.size()][3];
//		for (int i = 0; i < iterations.length; i++) {
//			iterations[i] = iterationsVec.elementAt(i);
//		}

//		Files.qsub(scriptFileDir + "trio_", qsubLogsDir, iterations.length, command, iterations, 2000, 24);
//		qsubFiles = Files.list(scriptFileDir, "trio", ".qsub", false, false);
//		jobNamesWithAbsolutePaths = new Vector<String>(qsubFiles.length);
//		for (int i = 0; i < qsubFiles.length; i++) {
//			jobNamesWithAbsolutePaths.add(scriptFileDir+qsubFiles[i]);
//		}
//		Files.qsubMultiple(jobNamesWithAbsolutePaths, null, scriptFileDir, "qsubmultiple_", qsubFiles.length<16?qsubFiles.length:16, true, null, -1, 22000, 24);
		Files.qsubMultiple(qsubFilesVec, null, scriptFileDir, "qsubmultiple_", qsubFilesVec.size() < 16 ? qsubFilesVec.size() : 16, true, null, -1, 22000, 24);
	}

	public static void main(String[] args) {
		String pedegreeFileFullPath, bamFilesDir, bcfFilesDir, denovogearResultDir, scriptFileDir, qsubLogsDir, bamFilePrefix;
		String[] commands;
		Logger log;

//		pedegreeFileFullPath = "/home/pankrat2/shared/logan/denovogear/all.ped";
//		bamFilesDir = "/home/pankrat2/shared/logan/denovogear/data_source_bam/";
//		bcfFilesDir = "/home/pankrat2/shared/logan/denovogear/data_converted_bcf/";
//		denovogearResultDir = "/home/pankrat2/shared/logan/denovogear/results/";
//		scriptFileDir = "/home/pankrat2/shared/logan/denovogear/scripts/";
//		qsubLogsDir = "/home/pankrat2/shared/logan/denovogear/logs/";
//		bamFilePrefix = "dedup_";
		pedegreeFileFullPath = "/home/pankrat2/shared/logan/denovogear/alltrios.ped";
		bamFilesDir = "/scratch/trio_bams/";
		bcfFilesDir = "/scratch/trio_bcfs/";
		denovogearResultDir = "/home/pankrat2/shared/logan/denovogear/results/";
		scriptFileDir = "/home/pankrat2/shared/logan/denovogear/scripts/";
		qsubLogsDir = "/home/pankrat2/shared/logan/denovogear/logs/";
		bamFilePrefix = "dedup_";

		commands = new String[] {"ped=", "bamdir=", "bcfdir=", "out=", "scriptdir=", "qsublogsdir=", "bamprefix="};
//		commandVariables = new String[] {pedegreeFileFullPath, bamFilesDir, bcfFilesDir, denovogearResultDir};
		String usage = "\n"+
				"cnv.analysis.DeNovoCNV analysis normally contains the following steps, split into several rounds, within each of which all the steps will be finished by just calling the program once." +
				"1st round:" +
				"   (1) full path of the pedegree file (i.e. " + commands[0] + pedegreeFileFullPath + " (the default));" +
				"   (2) directory of the source data bam files (i.e. " + commands[1] + bamFilesDir + " (the default));" +
				"   (3) directory to place the interim data bcf files (i.e. " + commands[2] + bcfFilesDir + " (the default));" +
				"   (4) directory to place the results of DeNovoGear (i.e. " + commands[3] + denovogearResultDir + " (the default));" +
				"   (5) directory to place the script files generated by this program (i.e. " + commands[4] + scriptFileDir + " (the default));" +
				"   (6) directory to place the log files for the scripts to run (i.e. " + commands[5] + qsubLogsDir + " (the default));" +
				"   (7) prefix of bam files (i.e. " + commands[6] + bamFilePrefix + " (the default));" +
				"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].startsWith(commands[0])) {
				pedegreeFileFullPath = args[i].split("=")[1];
			} else if (args[i].startsWith(commands[1])) {
				bamFilesDir = args[i].split("=")[1];
			} else if (args[i].startsWith(commands[2])) {
				bcfFilesDir = args[i].split("=")[1];
			} else if (args[i].startsWith(commands[3])) {
				denovogearResultDir = args[i].split("=")[1];
			} else if (args[i].startsWith(commands[4])) {
				scriptFileDir = args[i].split("=")[1];
			} else if (args[i].startsWith(commands[5])) {
				qsubLogsDir = args[i].split("=")[1];
			} else if (args[i].startsWith(commands[6])) {
				bamFilePrefix = args[i].split("=")[1];
			} else {
				System.err.println(usage);
				System.err.println("Error - invalid argument: "+args[i]);
				System.exit(1);
			}
		}

		log = new Logger(ext.parseDirectoryOfFile(pedegreeFileFullPath) + "Genvisis_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
		log.report(new SimpleDateFormat("MMM dd HH:mm:ss").format(new Date()) + ": starting generating qsub files for DeNovoGear.java");

		generateScripts(pedegreeFileFullPath, bamFilesDir, bcfFilesDir, denovogearResultDir, scriptFileDir, qsubLogsDir, bamFilePrefix, log);

		log.report(new SimpleDateFormat("MMM dd HH:mm:ss").format(new Date()) + ": completed generating qsub files for DeNovoGear.java.");
	}

}
