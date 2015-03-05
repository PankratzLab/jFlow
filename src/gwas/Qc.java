package gwas;

import java.io.*;
import java.util.*;
import common.*;

public class Qc {
	public static void fullGamut(String dir, boolean keepGenomeInfoForRelatedsOnly) {
		Logger log;
		long time;
		
		time = new Date().getTime();
		
		dir = ext.verifyDirFormat(dir);
		log = new Logger(dir+"fullGamutOfMarkerAndSampleQC.log");

		
		new File(dir+"QC/").mkdirs();
		if (!Files.exists(dir+"QC/freq.frq")) {
			log.report("Running --freq");
			CmdLine.run("plink --bfile ../plink --maf 0 --geno 1 --mind 1 --freq --out freq --noweb", dir+"QC/");
		}
		if (!Files.exists(dir+"QC/missing.imiss")) {
			log.report("Running --missing");
			CmdLine.run("plink --bfile ../plink --maf 0 --geno 1 --mind 1 --missing --out missing --noweb", dir+"QC/");
		}
		if (!Files.exists(dir+"QC/test.missing.missing")) {
			log.report("Running --test-missing");
			CmdLine.run("plink --bfile ../plink --maf 0 --geno 1 --mind 1 --test-missing --out test.missing --noweb", dir+"QC/");
		}
		if (!Files.exists(dir+"QC/hardy.hwe")) {
			log.report("Running --hardy");
			CmdLine.run("plink --bfile ../plink --maf 0 --geno 1 --mind 1 --hardy --out hardy --noweb", dir+"QC/");
		}
		if (!Files.exists(dir+"QC/mishap.missing.hap")) {
			log.report("Running --test-mishap");
			CmdLine.run("plink --bfile ../plink --maf 0 --geno 1 --mind 1 --test-mishap --out mishap --noweb", dir+"QC/");
		}
		if (!Files.exists(dir+"QC/gender.assoc")) {
			log.report("Running --assoc gender");
			CmdLine.run("plink --bfile ../plink --maf 0 --geno 1 --mind 1 --pheno ../plink.fam --mpheno 3 --assoc --out gender --noweb", dir+"QC/");
		}
		if (!Files.exists(dir+"QC/gender.missing")) {
			log.report("Running --test-missing gender");
			CmdLine.run("plink --bfile ../plink --maf 0 --geno 1 --mind 1 --pheno ../plink.fam --mpheno 3 --test-missing --out gender", dir+"QC/");
		}
		
		if (!Files.exists(dir+"QC/miss_drops.dat")) {
			new File(dir+"QC/miss.crf").delete();
			MarkerQC.parseParameters(dir+"QC/miss.crf", log); // first generates a file with the defaults
			Files.writeList(Array.addStrToArray("dir="+dir+"QC/", HashVec.loadFileToStringArray(dir+"QC/miss.crf", false, new int[] {0}, false), 1), dir+"QC/miss.crf");
			MarkerQC.parseParameters(dir+"QC/miss.crf", log); // second runs the filtering
		}
		
		new File(dir+"ldPruning/").mkdirs();
		if (!Files.exists(dir+"ldPruning/plink.bed")) {
			log.report("Running --exclude miss_drops.dat");
			CmdLine.run("plink --bfile ../plink --exclude ../QC/miss_drops.dat --make-bed --noweb", dir+"ldPruning/");
		}
		if (!Files.exists(dir+"ldPruning/plink.prune.in")) {
			log.report("Running --indep-pairwise 50 5 0.3");
			CmdLine.run("plink --bfile plink --indep-pairwise 50 5 0.3", dir+"ldPruning/");
		}

		new File(dir+"genome/").mkdirs();
		if (!Files.exists(dir+"genome/plink.bed")) {
			log.report("Running --extract plink.prune.in");
			CmdLine.run("plink --bfile ../plink --extract ../ldPruning/plink.prune.in --make-bed --noweb", dir+"genome/");
		}
		if (!Files.exists(dir+"genome/plink.genome")) {
			log.report("Running --genome"+(keepGenomeInfoForRelatedsOnly?" --min 0.1":""));
			CmdLine.run("plink --bfile plink --genome"+(keepGenomeInfoForRelatedsOnly?" --min 0.1":""), dir+"genome/");
		}
		if (!keepGenomeInfoForRelatedsOnly) {
			log.report("Running --mds-plot 20");
			CmdLine.run("plink --bfile plink --read-genome plink.genome --cluster --mds-plot 20 --out mds20", dir+"genome/");
			
		}
		if (!Files.exists(dir+"genome/plink.genome_keep.dat")) {
			log.report("Running flagRelateds");
			Plink.flagRelateds(dir+"genome/plink.genome", dir+"genome/plink.fam", dir+"QC/missing.imiss", dir+"genome/lrr_sd.xln", Plink.FLAGS, Plink.THRESHOLDS, 4);
		}
		
		new File(dir+"ancestry/").mkdirs();
		if (!Files.exists(dir+"ancestry/unrelateds.txt")) {
			log.report("Copying genome/plink.genome_keep.dat to ancestry/unrelateds.txt");
			Files.copyFile(dir+"genome/plink.genome_keep.dat", dir+"ancestry/unrelateds.txt");
		}
		if (!Files.exists(dir+"ancestry/plink.bed")) {
			log.report("Running --extract plink.prune.in (again, this time to ancestry/)");
			CmdLine.run("plink --bfile ../plink --extract ../ldPruning/plink.prune.in --make-bed --noweb", dir+"ancestry/");
		}
		
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
		try {
			fullGamut(dir, keepGenomeInfoForRelatedsOnly);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
