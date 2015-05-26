package seq.analysis;

import java.io.File;

import seq.manage.VCFOps;
import seq.manage.VCFOps.VcfPopulation;
import common.Files;
import common.Logger;
import common.PSF;
import common.ext;

public class Ancestry {

	public static void runAncestry(String vcf, String hapMapRsIds, String g1000RsIds, String vpopFile, String outputDirectory) {
		new File(outputDirectory).mkdirs();

		Logger log = new Logger(outputDirectory + "ancestry.log");
		String curVCF = VCFOps.extractIDs(vcf, g1000RsIds, outputDirectory, true, true, log);
		curVCF = VCFOps.extractIDs(curVCF, hapMapRsIds, outputDirectory, true, true, log);
		String[] filenames = VcfPopulation.splitVcfByPopulation(curVCF, vpopFile, log);
		for (int i = 0; i < filenames.length; i++) {
			if (filenames[i].endsWith(".USE.vcf.gz")) {
				curVCF = VCFOps.removeFilteredVariants(vcf, true, true, log);
				Files.copyFile(curVCF, ext.parseDirectoryOfFile(curVCF) + "ancestry.vcf.gz");
				Files.copyFile(curVCF + ".tbi", ext.parseDirectoryOfFile(curVCF) + "ancestry.vcf.gz.tbi");
				curVCF = ext.parseDirectoryOfFile(curVCF) + "ancestry.vcf.gz";
				VCFOps.vcfGwasQC(curVCF, log);
				return;
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = "avcf.vcf.gz";
		String hapMapRsIds = null;
		String g1000RsIds = null;
		String vpopFile = null;
		String outputDirectory = null;

		String usage = "\n" + "seq.analysis.ancestry requires 0-1 arguments\n";
		usage += "   (1) vcf (i.e. vcf=" + vcf + " (default))\n" + "";
		usage += "   (2) file of hapMap rs ids to extract (i.e. hapMap=" + vcf + " (default))\n" + "";
		usage += "   (3) file of 1000 genome rs ids to extract (i.e. g1000=" + vcf + " (default))\n" + "";
		usage += "   (4) vpopFile  (i.e. vpop=" + vpopFile + " (default))\n" + "";
		usage += PSF.Ext.getOutputDirCommand(4, "");

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcf=")) {
				vcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("hapMap=")) {
				hapMapRsIds = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("g1000=")) {
				g1000RsIds = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("vpop=")) {
				vpopFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
				outputDirectory = args[i].split("=")[1];
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
			runAncestry(vcf, hapMapRsIds, g1000RsIds, vpopFile, outputDirectory);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
