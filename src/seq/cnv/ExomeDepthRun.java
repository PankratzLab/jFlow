package seq.cnv;

import java.io.File;

import seq.cnv.ExomeDepth.ExomeDepthAnalysis;
import seq.manage.BamOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.ext;

/**
 * @author lane0212 Handles sample and global exclusions for the reference set, need in particular when related samples are run together (i.e trios)
 */
public class ExomeDepthRun {

	public static void runExomeDepth(String bams, String vpopFile, String outputDir, String outputRoot, String rLoc, boolean somaticMode, int numthreads, Logger log) {
		VcfPopulation vpop = null;
		String[] allReferenceBamFiles = Files.isDirectory(bams) ? Files.listFullPaths(bams, BamOps.BAM_EXT, false) : HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, true);
		outputDir = outputDir == null ? ext.parseDirectoryOfFile(bams) : outputDir;
		new File(outputDir).mkdirs();
		ExomeDepth exomeDepth = new ExomeDepth(allReferenceBamFiles, allReferenceBamFiles, outputDir, outputRoot, rLoc, log);
		if (!Files.exists(exomeDepth.getCountFile())) {
			log.reportTimeWarning("Did not find " + exomeDepth.getCountFile() + ", generating it now (takes a long time)");
			exomeDepth.generateCountFile();
		} else {
			log.reportTimeWarning("Using existing count file " + exomeDepth.getCountFile());
		}
		if (vpopFile == null) {
			log.reportTimeWarning("A vpopulation file was not provided, sample specific and global exclusions will not be applied");
		} else {
			if (somaticMode) {
				vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.TUMOR_NORMAL, log);
			}else{
				vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.EXOME_DEPTH, log);
			}
			vpop.report();
			System.exit(1);

			exomeDepth.parseVpop(vpop);
		}
		ExomeDepthAnalysis[] eDepthAnalysis = ExomeDepth.callCNVs(exomeDepth, outputDir, outputRoot, numthreads, log);
		log.reportTimeInfo("Finished running exome depth for " + eDepthAnalysis.length + " .bam files");
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String bams = "bams/";
		String outputDir = null;
		String outputRoot = "ExomeDepth";
		int numthreads = 1;
		String vpopFile = null;
		String logfile = null;
		String Rloc = null;
		boolean somaticMode =false;
		Logger log;

		String usage = "\n" + "seq.analysis.ExomeDepth requires 0-1 arguments\n";
		usage += "   (1) full path to a directory of or file of bams (i.e. bams=" + bams + " (default))\n" + "";
		usage += PSF.Ext.getOutputDirCommand(2, outputDir);
		usage += "   (3) output root command (i.e. root=" + outputRoot + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(4, numthreads);
		usage += "   (5) full path to a v population file, individuals with the same population will not be used as ref(i.e. vpop= (no default))\n" + "";
		usage += "   (6) alternative R location (i.e. rDir= (no default))\n" + "";
		usage += "   (7) somatic mode (i.e. somaticMode="+somaticMode+" (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("bams=")) {
				bams = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("vpop=")) {
				vpopFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
				outputDir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("rDir=")) {
				Rloc = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("root=")) {
				outputRoot = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("somaticMode=")) {
				somaticMode = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
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
			log = new Logger(logfile);
			runExomeDepth(bams, vpopFile, outputDir, outputRoot, Rloc, somaticMode, numthreads, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
