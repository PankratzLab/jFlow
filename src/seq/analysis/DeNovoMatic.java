package seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;

import javax.jms.IllegalStateException;

import seq.analysis.GATK.MutectTumorNormal;
import seq.analysis.Mutect2.MUTECT_RUN_TYPES;
import seq.manage.BamOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

/**
 * Hijack somatic calling to screen for denovo's
 *
 */
public class DeNovoMatic {
	public static void run(String vpopFile, String fileOfBams, String outputDir, String ponVcf, GATK gatk, MUTECT_RUN_TYPES type, int numThreads, int numSampleThreads, Logger log) throws IllegalStateException {
		new File(outputDir).mkdirs();
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
		String matchFile = outputDir + ext.rootOf(vpopFile) + ".matchedbams.txt";
		if (!Files.exists(matchFile)) {
			prepareMatchedBamFile(vpop, fileOfBams, matchFile, log);
		}

		String callDir = outputDir + "mutect2Calls/";
		new File(callDir).mkdirs();
		// Mutect2.run(null, matchFile, callDir, ponVcf, gatk, type, numThreads, numSampleThreads, log);
		MutectTumorNormal[] results = Mutect2.callSomatic(matchFile, outputDir, ponVcf, gatk, numThreads, numSampleThreads, false, log);
		System.out.println(results.length);
	}

	private static void prepareMatchedBamFile(VcfPopulation vpop, String fileOfBams, String output, Logger log) {
		String[] bams = HashVec.loadFileToStringArray(fileOfBams, false, new int[] { 0 }, true);
		ArrayList<String> toWrite = new ArrayList<String>();
		Hashtable<String, String> match = new Hashtable<String, String>();
		for (int i = 0; i < bams.length; i++) {
			match.put(BamOps.getSampleName(bams[i]), bams[i]);
		}
		for (String family : vpop.getSubPop().keySet()) {
			Set<String> fam = vpop.getSubPop().get(family);
			if (fam.size() != 3) {
				throw new IllegalArgumentException("trio " + family + " had " + fam.size() + " members");
			}
			String off = null;
			String p1 = null;
			String p2 = null;
			for (String ind : fam) {
				String key = ind.replaceAll(".variant.*", "");
				if (vpop.getPopulationForInd(ind, RETRIEVE_TYPE.SUPER)[0].equals("OFFSPRING")) {
					off = match.get(key);
				} else if (vpop.getPopulationForInd(ind, RETRIEVE_TYPE.SUPER)[0].equals("CONTROL")) {
					if (p1 == null) {
						p1 = match.get(key);
					} else {
						p2 = match.get(key);
					}
				} else {
					throw new IllegalArgumentException("Super pop can only be OFFSPRING or CONTROL and found "+vpop.getPopulationForInd(ind, RETRIEVE_TYPE.SUPER)[0]);
				}
			}
			if (off == null || p1 == null || p2 == null) {
				throw new IllegalArgumentException("Could not detect all bam files for trio " + family);
			}

			toWrite.add(p1 + "\t" + off);
			toWrite.add(p2 + "\t" + off);
		}
		Files.writeList(Array.toStringArray(toWrite), output);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String referenceGenomeFasta = "hg19_canonical.fa";
		String knownSnps = "dbsnp_138.hg19.vcf";
		String gatkLocation = "GATK_3_5/";
		String cosmic = "b37_cosmic_v54_120711.hg19_chr.vcf";
		String regions = "AgilentCaptureRegions.bed";
		String ponVCF = null;
		int numthreads = 1;
		int numSampleThreads = 4;
		String outputDir = "mutect/";
		String vpopFile = null;
		String fileOfBams = null;

		String usage = "\n" + "seq.analysis.DeNovoMatic requires 0-1 arguments\n";
		usage += "   (2) full path to a reference genome (i.e. ref=" + referenceGenomeFasta + " (default))\n" + "";
		usage += "   (3) known dbsnp snps (i.e. knownSnps=" + knownSnps + " (default))\n" + "";
		usage += "   (4) cosmic snps (i.e. cosmic=" + cosmic + " (default))\n" + "";
		usage += "   (5) regions for calling (i.e. regions=" + regions + " (default))\n" + "";
		usage += "   (6) output root directory (i.e. outputDir=" + outputDir + " (default))\n" + "";
		usage += "   (8) number of threads (i.e. numthreads=" + numthreads + " (default))\n" + "";
		usage += "   (9) gatk directory (i.e. gatk=" + gatkLocation + " (default))\n" + "";
		usage += "   (11) pon vcf (i.e. ponVcf= (no default))\n" + "";
		usage += "   (12) number of threads per sample (i.e. numSampleThreads=" + numSampleThreads + " (default))\n" + "";
		usage += "   (13) full to a vpop file (i.e. vpop= (no default))\n" + "";
		usage += "   (14) full path to a file of bams (i.e. bams= (no default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("ref=")) {
				referenceGenomeFasta = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("vpop=")) {
				vpopFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bams=")) {
				fileOfBams = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("gatk=")) {
				gatkLocation = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ponVCF=")) {
				ponVCF = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("knownSnps=")) {
				knownSnps = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cosmic=")) {
				cosmic = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("regions=")) {
				regions = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outputDir=")) {
				outputDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("numthreads=")) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("numSampleThreads=")) {
				numSampleThreads = ext.parseIntArg(args[i]);
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
		GATK gatk = new GATK(gatkLocation, referenceGenomeFasta, knownSnps, regions, cosmic, true, false, true, log);
		try {
			run(vpopFile, fileOfBams, outputDir, ponVCF, gatk, MUTECT_RUN_TYPES.CALL_SOMATIC, numthreads, numSampleThreads, log);
		} catch (IllegalStateException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
