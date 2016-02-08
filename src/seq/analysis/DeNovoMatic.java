package seq.analysis;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import javax.jms.IllegalStateException;

import seq.analysis.GATK.MutectTumorNormal;
import seq.analysis.Mutect2.MUTECT_RUN_TYPES;
import seq.manage.BamOps;
import seq.manage.VCFOps;
import seq.manage.VCOps;
import seq.manage.VCFOps.HEADER_COPY_TYPE;
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
	private static final String[] ACCEPTED_FILTER_BYPASS = new String[] { "str_contraction" };

	public static void run(String vpopFile, String fileOfBams, String outputDir, String ponVcf, double freqFilter, GATK gatk, MUTECT_RUN_TYPES type, int numThreads, int numSampleThreads, Logger log) throws IllegalStateException {
		new File(outputDir).mkdirs();

		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.DENOVO, log);
		vpop.report();
		String matchFile = outputDir + ext.rootOf(vpopFile) + ".matchedbams.txt";
		if (!Files.exists(matchFile) || Files.countLines(matchFile, 0) != vpop.getSuperPop().get(VcfPopulation.OFFSPRING).size() * 2) {
			prepareMatchedBamFile(vpop, fileOfBams, matchFile, log);
		}

		// Mutect2.run(null, matchFile, callDir, ponVcf, gatk, type, numThreads, numSampleThreads, log);
		MutectTumorNormal[] results = Mutect2.callSomatic(matchFile, outputDir, ponVcf, gatk, null, null, numThreads, numSampleThreads, false, log);

		MergeFamResult[] resultsMerge = combineResultsForFamilies(gatk, vpop, results, outputDir, numThreads, log);
		String[] finalFilesToMerge = new String[resultsMerge.length];
		for (int i = 0; i < resultsMerge.length; i++) {
			resultsMerge[i].scanForDenovo();
			finalFilesToMerge[i] = resultsMerge[i].getPotentialDenovoVcf();
		}

		String mergeDenovoOut = outputDir + ext.rootOf(vpopFile) + ".merge.denovo.vcf";

		gatk.mergeVCFs(finalFilesToMerge, mergeDenovoOut, numThreads, false, log);
		String annotatedVcf = GATK_Genotyper.annotateOnlyWithDefualtLocations(mergeDenovoOut,null, false, false, log);
		log.reportTimeInfo("Filtering " + annotatedVcf);
		String annoFreq = VCFOps.getAppropriateRoot(annotatedVcf, false) + ".freq_" + freqFilter + ".func.vcf.gz";
		filterByFreq(annotatedVcf, annoFreq, freqFilter, log);
		log.reportTimeInfo("FIN");

	}

	private static void filterByFreq(String invcf, String outVcf, double freq, Logger log) {

		VCFFileReader reader = new VCFFileReader(invcf, true);
		VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, outVcf, VCFOps.DEFUALT_WRITER_OPTIONS, log);
		int total = 0;
		int pass = 0;
		for (VariantContext vc : reader) {
			total++;
			if (!vc.hasAttribute("PopFreqMax")) {
				throw new IllegalArgumentException("Method expects PopFreqMax annotations for all variants");
			}
			String f = vc.getAttributeAsString("PopFreqMax", "-1");

			if (f.equals(".")) {
				f = "-1";
			}
			double freqVc = Double.parseDouble(f);
			if (freqVc < freq && VCOps.isHighModLowSNP_EFFImpact(vc)) {

				writer.add(vc);
				pass++;
			}
		}
		writer.close();
		reader.close();
		log.reportTimeInfo("Scanned " + total + " variants and " + pass + " passed freq threshold of " + freq);
	}

	private static class MergeFamResult {
		private MutectTumorNormal[] famResults;
		private String mergedVCF;
		private String potentialDenovoVcf;
		private String off;
		private HashSet<String> offCombo;
		private String p1;
		private String p2;
		private Logger log;

		public MergeFamResult(MutectTumorNormal[] famResults, String mergedVCF, String off, String p1, String p2, Logger log) {
			super();
			this.famResults = famResults;
			this.mergedVCF = mergedVCF;
			this.potentialDenovoVcf = VCFOps.getAppropriateRoot(mergedVCF, false) + ".denovo.vcf.gz";
			this.off = off;
			this.offCombo = new HashSet<String>();
			offCombo.add(off + ".variant");
			offCombo.add(off + ".variant2");
			this.p1 = p1;
			this.p2 = p2;
			this.log = log;
		}

		public String getPotentialDenovoVcf() {
			return potentialDenovoVcf;
		}

		private void scanForDenovo() {

			if (!VCFOps.existsWithIndex(potentialDenovoVcf)) {
				VCFFileReader reader = new VCFFileReader(mergedVCF, true);
				VariantContextWriter writer = VCFOps.initWriter(potentialDenovoVcf, VCFOps.DEFUALT_WRITER_OPTIONS, reader.getFileHeader().getSequenceDictionary());
				VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
				int numPossible = 0;
				int numTotal = 0;
				int numByPassFilters = 0;
				for (VariantContext vc : reader) {
					numTotal++;
					VariantContext vcSubOff = VCOps.getSubset(vc, offCombo);
					Genotype g1 = vcSubOff.getGenotype(0);
					Genotype g2 = vcSubOff.getGenotype(1);

					if (g1.sameGenotype(g2)) {// called somatic in both files
						writer.add(vc);
						numPossible++;
						// Set<String> filters = vc.getFilters();
						// if (filters.size() == 0 || filters.size() == 1) {
						// boolean ok = true;
						// for (String filt : filters) {
						// if (ext.indexOfStr(filt, ACCEPTED_FILTER_BYPASS) < 0) {
						// ok = false;
						// }
						// }
						// if (ok) {
						// if(filters.size()==1){
						// numByPassFilters++;
						// }
						//
						// }
						// }
					}
					if (numTotal % 100000 == 0) {
						log.reportTimeInfo("Detected " + numPossible + " possible denovos after scanning " + numTotal);
					}
				}
				reader.close();
				writer.close();
				log.reportTimeInfo("Detected " + numPossible + " possible denovos after scanning " + numTotal);
				if (numByPassFilters > 0) {
					log.reportTimeWarning(numByPassFilters + " possible denovos were allowed past mutects filters " + Array.toStr(ACCEPTED_FILTER_BYPASS));

				}
			}
		}
	}

	private static MergeFamResult[] combineResultsForFamilies(GATK gatk, VcfPopulation vpop, MutectTumorNormal[] results, String outputDir, int numThreads, Logger log) {
		Hashtable<String, ArrayList<MutectTumorNormal>> offSpringMatch = new Hashtable<String, ArrayList<MutectTumorNormal>>();
		if (results.length % 2 != 0) {
			throw new IllegalArgumentException("Expecting even number of results");

		}
		ArrayList<MergeFamResult> mergeResults = new ArrayList<MergeFamResult>();
		for (int i = 0; i < results.length; i++) {
			String currentVCF = results[i].getReNamedFilteredVCF();
			String[] inds = VCFOps.getSamplesInFile(results[i].getReNamedFilteredVCF());
			if (inds.length != 2) {
				throw new IllegalArgumentException("Expecting two samples per vcf, found " + inds.length + " in " + currentVCF);
			}
			boolean found = false;
			for (int j = 0; j < inds.length; j++) {
				String fam = vpop.getPopulationForInd(inds[j], RETRIEVE_TYPE.SUB)[0];
				if (vpop.getPopulationForInd(inds[j], RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.OFFSPRING)) {
					if (!offSpringMatch.containsKey(fam)) {
						offSpringMatch.put(fam, new ArrayList<MutectTumorNormal>());
					}
					found = true;
					offSpringMatch.get(fam).add(results[i]);
				}
			}
			if (!found) {
				throw new IllegalArgumentException("Could not find matching entries for " + currentVCF);
			}
		}

		for (String fam : vpop.getSubPop().keySet()) {
			String outputMerge = outputDir + fam + ".merge.vcf.gz";
			log.reportTimeInfo("Combining variants for " + fam + " to " + outputMerge);
			MutectTumorNormal[] famResults = offSpringMatch.get(fam).toArray(new MutectTumorNormal[offSpringMatch.get(fam).size()]);

			ArrayList<String> toMerge = new ArrayList<String>();
			for (int i = 0; i < famResults.length; i++) {
				toMerge.add(famResults[i].getReNamedOutputVCF());
			}
			String[] vcfsToMerge = Array.toStringArray(toMerge);
			if (vcfsToMerge.length != 2) {
				throw new IllegalArgumentException("Internal error, need to merge two vcfs");
			}
			if (!Files.exists(outputMerge))
				gatk.mergeVCFs(vcfsToMerge, outputMerge, numThreads, false, log);
			String[] famMembers = vpop.getOffP1P2ForFam(fam);
			MergeFamResult mergeFamResult = new MergeFamResult(famResults, outputMerge, famMembers[0], famMembers[1], famMembers[2], log);
			mergeResults.add(mergeFamResult);
		}
		return mergeResults.toArray(new MergeFamResult[mergeResults.size()]);
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
				// String key = ind.replaceAll(".variant.*", "");
				String key = ind;
				if (vpop.getPopulationForInd(ind, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.OFFSPRING)) {
					off = match.get(key);
				} else if (vpop.getPopulationForInd(ind, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.PARENTS)) {
					if (p1 == null) {
						p1 = match.get(key);
					} else {
						p2 = match.get(key);
					}
				} else {
					throw new IllegalArgumentException("Super pop can only be " + VcfPopulation.OFFSPRING + "  or " + VcfPopulation.OFFSPRING + " and found " + vpop.getPopulationForInd(ind, RETRIEVE_TYPE.SUPER)[0]);
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
		double freqFilter = .01;

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
		usage += "   (15) full path to a file of supporting snps (i.e. supportSnps= (no default))\n" + "";

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
			run(vpopFile, fileOfBams, outputDir, ponVCF, freqFilter, gatk, MUTECT_RUN_TYPES.CALL_SOMATIC, numthreads, numSampleThreads, log);
		} catch (IllegalStateException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
