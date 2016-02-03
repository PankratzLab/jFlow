package seq.analysis;

import java.io.File;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import seq.analysis.GATK.GenotypeRefiner;
import seq.manage.VCFOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import common.HashVec;
import common.Logger;
import common.ext;

public class GenotypeRefinement {

	private static void refineGenotypes(String vcf, String ped, String outputDir, GATK gatk) {
		new File(outputDir).mkdirs();
		Logger log = new Logger(outputDir + "gtRefine.log");
		log.reportTimeWarning("Subsetting vcf to samples found in " + ped);
		String[] sampsinPed = HashVec.loadFileToStringArray(ped, false, new int[] { 1 }, true);
		String[] vcfSamps = VCFOps.getSamplesInFile(vcf);
		if (!ext.containsAll(sampsinPed, vcfSamps)) {
			throw new IllegalArgumentException("All pedigree samples must be in ped file");
		}
		Set<String> sub = new HashSet<String>();
		Hashtable<String, Set<String>> pop = new Hashtable<String, Set<String>>();
		pop.put(ext.rootOf(ped), sub);
		for (int i = 0; i < sampsinPed.length; i++) {
			pop.get(ext.rootOf(ped)).add(sampsinPed[i]);
		}
		VcfPopulation vpop = new VcfPopulation(pop, pop, POPULATION_TYPE.ANY, log);
		String vpopTmp = outputDir + ext.rootOf(ped) + ".vpop";
		vpop.dump(vpopTmp);
		String subVcf = VcfPopulation.splitVcfByPopulation(vcf, vpopTmp, true, true, log)[0];

		GenotypeRefiner genotypeRefiner = gatk.refineGenotypes(subVcf, ped, outputDir, log);
		
		//TODO, filter low qual, high qual denovo -> pipe to tally
		
		
		System.out.println(genotypeRefiner);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String referenceGenomeFasta = "hg19_canonical.fa";
		String supportingSnps = "dbsnp_138.hg19.vcf";
		String gatkLocation = "GATK_3_5/";

		String outputDir = "genotypeRefine/";
		String inputVcf = null;
		String ped = null;
		String usage = "\n" + "seq.analysis.DeNovoMatic requires 0-1 arguments\n";
		usage += "   (1) full path to a reference genome (i.e. ref=" + referenceGenomeFasta + " (default))\n" + "";
		usage += "   (2) output root directory (i.e. outputDir=" + outputDir + " (default))\n" + "";
		usage += "   (3) gatk directory (i.e. gatk=" + gatkLocation + " (default))\n" + "";
		usage += "   (4) input vcf (i.e. vcf= (no default))\n" + "";
		usage += "   (5) full to a ped file (i.e. ped= (no default))\n" + "";
		usage += "   (6) full path to a file of supporting snps (i.e. supportSnps= (no default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("ref=")) {
				referenceGenomeFasta = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ped=")) {
				ped = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("vcf=")) {
				inputVcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("gatk=")) {
				gatkLocation = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outputDir=")) {
				outputDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("supportSnps=")) {
				supportingSnps = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		new File(outputDir).mkdirs();
		Logger log = new Logger(outputDir + "TN.log");
		GATK gatk = new GATK(gatkLocation, referenceGenomeFasta, null, null, null, true, false, true, log);
		gatk.setSupportingSnps(supportingSnps);
		refineGenotypes(inputVcf, ped, outputDir, gatk);
	}

}
