package gwas;

import java.io.*;
import java.util.*;

import common.*;

public class Qc {
    
    /** A rough listing of the Folders created by fullGamut */
    public static String[] FOLDERS_CREATED = {"markerQC/", "sampleQC/", "ldPruning/", "genome/", "ancestry/"};
    /** A rough listing of the files created, by folder, by fullGamut */
    public static String[][] FILES_CREATED = {
        {"plink.bed", "freq.frq", "missing.imiss", "test.missing.missing", "hardy.hwe", "mishap.missing.hap", "gender.assoc", "gender.missing", "miss_drops.dat"},
        {"plink.bed", "missing.imiss"},
        {"plink.bed", "plink.prune.in"},
        {"plink.bed", "plink.genome", "plink.genome_keep.dat"},
        {"plink.bed", "unrelateds.txt"}
    };
    
	public static void fullGamut(String dir, String plinkPrefix, boolean keepGenomeInfoForRelatedsOnly, Logger log) {
		long time;
		
		time = new Date().getTime();
		
		dir = ext.verifyDirFormat(dir);

		String plink = plinkPrefix == null ? "plink" : plinkPrefix;
		
		new File(dir+"markerQC/").mkdirs();
		log.report(ext.getTime() + "]\tRunning --mind 0.1");
		CmdLine.runDefaults("plink --bfile ../" + plink + " --mind 0.1 --make-bed --noweb", dir+"markerQC/", log);
		
		if (!Files.exists(dir+"markerQC/freq.frq")) {
			log.report(ext.getTime() + "]\tRunning --freq");
			CmdLine.runDefaults("plink --bfile " + plink + " --maf 0 --geno 1 --mind 1 --freq --out freq --noweb", dir+"markerQC/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"markerQC/missing.imiss")) {
			log.report(ext.getTime() + "]\tRunning --missing");
			CmdLine.runDefaults("plink --bfile " + plink + " --maf 0 --geno 1 --mind 1 --missing --out missing --noweb", dir+"markerQC/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"markerQC/test.missing.missing")) {
			log.report(ext.getTime() + "]\tRunning --test-missing");
			CmdLine.runDefaults("plink --bfile " + plink + " --maf 0 --geno 1 --mind 1 --test-missing --out test.missing --noweb", dir+"markerQC/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"markerQC/hardy.hwe")) {
			log.report(ext.getTime() + "]\tRunning --hardy");
			CmdLine.runDefaults("plink --bfile " + plink + " --maf 0 --geno 1 --mind 1 --hardy --out hardy --noweb", dir+"markerQC/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"markerQC/mishap.missing.hap")) {
			log.report(ext.getTime() + "]\tRunning --test-mishap");
			CmdLine.runDefaults("plink --bfile " + plink + " --maf 0 --geno 1 --mind 1 --test-mishap --out mishap --noweb", dir+"markerQC/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"markerQC/gender.assoc")) {
			log.report(ext.getTime() + "]\tRunning --assoc gender");
			CmdLine.runDefaults("plink --bfile " + plink + " --maf 0 --geno 1 --mind 1 --pheno plink.fam --mpheno 3 --assoc --out gender --noweb", dir+"markerQC/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"markerQC/gender.missing")) {
			log.report(ext.getTime() + "]\tRunning --test-missing gender");
			CmdLine.runDefaults("plink --bfile " + plink + " --maf 0 --geno 1 --mind 1 --pheno plink.fam --mpheno 3 --test-missing --out gender --noweb", dir+"markerQC/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		
		if (!Files.exists(dir+"markerQC/miss_drops.dat")) {
			new File(dir+"markerQC/miss.crf").delete();
			int runCode1 = MarkerQC.parseParameters(dir+"markerQC/miss.crf", log, false); // first generates a file with the defaults
			if (runCode1 != 0) {
			    // ERROR!  TODO not sure if we should quit here; for now, continue;
			}
			Files.writeList(Array.addStrToArray("dir="+dir+"markerQC/", HashVec.loadFileToStringArray(dir+"markerQC/miss.crf", false, new int[] {0}, false), 1), dir+"markerQC/miss.crf");
			int runCode2 = MarkerQC.parseParameters(dir+"markerQC/miss.crf", log, false); // second runs the filtering
			if (runCode2 != 0) {
			    // ERROR  TODO not sure if we should quit here; for now, continue;
			}
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }

		new File(dir+"sampleQC/").mkdirs();
		if (!Files.exists(dir+"sampleQC/" + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --exclude miss_drops.dat");
			CmdLine.runDefaults("plink --bfile ../" + plink + " --exclude ../markerQC/miss_drops.dat --make-bed --noweb", dir+"sampleQC/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"sampleQC/missing.imiss")) {
			log.report(ext.getTime() + "]\tRunning --missing");
			CmdLine.runDefaults("plink --bfile " + plink + " --maf 0 --geno 1 --mind 1 --missing --out missing --noweb", dir+"sampleQC/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		
		new File(dir+"ldPruning/").mkdirs();
		if (!Files.exists(dir+"ldPruning/" + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --mind 0.05 (removes samples with callrate <95% for the markers that did pass QC)");
			CmdLine.runDefaults("plink --bfile ../sampleQC/" + plink + " --mind 0.05 --make-bed --noweb", dir+"ldPruning/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"ldPruning/" + plink + ".prune.in")) {
			log.report(ext.getTime() + "]\tRunning --indep-pairwise 50 5 0.3");
			CmdLine.runDefaults("plink --noweb --bfile " + plink + " --indep-pairwise 50 5 0.3", dir+"ldPruning/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }

		new File(dir+"genome/").mkdirs();
		if (!Files.exists(dir+"genome/" + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --extract " + plink + ".prune.in");
			CmdLine.runDefaults("plink --bfile ../ldPruning/" + plink + " --extract ../ldPruning/" + plink + ".prune.in --make-bed --noweb", dir+"genome/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"genome/" + plink + ".genome")) {
			log.report(ext.getTime() + "]\tRunning --genome"+(keepGenomeInfoForRelatedsOnly?" --min 0.1":""));
			CmdLine.runDefaults("plink --noweb --bfile " + plink + " --genome"+(keepGenomeInfoForRelatedsOnly?" --min 0.1":""), dir+"genome/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!keepGenomeInfoForRelatedsOnly && !Files.exists(dir + "genome/mds20.mds")) {
			log.report(ext.getTime() + "]\tRunning --mds-plot 20");
			CmdLine.runDefaults("plink --bfile " + plink + " --read-genome " + plink + ".genome --cluster --mds-plot 20 --out mds20 --noweb", dir+"genome/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"genome/" + plink + ".genome_keep.dat")) {
			log.report(ext.getTime() + "]\tRunning flagRelateds");
			Plink.flagRelateds(dir+"genome/" + plink + ".genome", dir+"genome/" + plink + ".fam", dir+"markerQC/missing.imiss", dir+"genome/lrr_sd.xln", Plink.FLAGS, Plink.THRESHOLDS, 4, false);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		
		new File(dir+"ancestry/").mkdirs();
		if (!Files.exists(dir+"ancestry/unrelateds.txt")) {
			log.report(ext.getTime() + "]\tCopying genome/" + plink + ".genome_keep.dat to ancestry/unrelateds.txt");
			Files.copyFile(dir+"genome/" + plink + ".genome_keep.dat", dir+"ancestry/unrelateds.txt");
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		if (!Files.exists(dir+"ancestry/" + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --extract " + plink + ".prune.in (again, this time to ancestry/)");
			CmdLine.runDefaults("plink --bfile ../genome/" + plink + " --make-bed --noweb", dir+"ancestry/", log);
		}
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		
		ancestry(dir+"ancestry/");
		
		
		System.out.println("Finished this round in " + ext.getTimeElapsed(time));
	}
	
	public static void ancestry(String dir) {
		Logger log;
		
		dir = ext.verifyDirFormat(dir);
		log = new Logger(dir+"ancestry.log");

		if (!Files.exists(dir+"unrelateds.txt")) {
			log.reportError("Error - need a file called unrelateds.txt with FID and IID pairs before we can proceed");
			return;
		}
		
		// TODO to be continued...
		
	}
	
	public static void fromParameters(String filename, Logger log) {
		Vector<String> params;

		params = Files.parseControlFile(filename, "gwas.Qc", new String[] { "dir=./", "# Make keepGenomeInfoForRelatedsOnly=false if the sample size is small and you want to run MDS plot as well", "keepGenomeInfoForRelatedsOnly=true"}, log);

		if (params != null) {
			params.add("log=" + log.getFilename());
			main(Array.toStringArray(params));
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "./";
		boolean keepGenomeInfoForRelatedsOnly = true;
		String logfile = null;
		Logger log;
		
		String usage = "\n" +
				"gwas.Qc requires 0-1 arguments\n" +
				"   (1) directory with plink.* files (i.e. dir=" + dir + " (default))\n" + 
				"   (2) if no MDS will be run, smaller file (i.e. keepGenomeInfoForRelatedsOnly=" + keepGenomeInfoForRelatedsOnly + " (default))\n" + 
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "./");
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("keepGenomeInfoForRelatedsOnly=")) {
				keepGenomeInfoForRelatedsOnly = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (logfile == null) {
			log = new Logger(dir+"fullGamutOfMarkerAndSampleQC.log");
		} else {
			log = new Logger(logfile);
		}
		try {
			fullGamut(dir, null, keepGenomeInfoForRelatedsOnly, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
