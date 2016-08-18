package org.genvisis.seq.analysis;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import javax.jms.IllegalStateException;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.GATK.MutectTumorNormal;
import org.genvisis.seq.analysis.GATK_Genotyper.ANNOVCF;
import org.genvisis.seq.analysis.Mutect2.MUTECT_RUN_TYPES;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFTumorNormalOps;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.genvisis.seq.manage.VCOps.GENOTYPE_INFO;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.FILTER_GENERATION_TYPE;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilterPass;

/**
 * Hijack somatic calling to screen for denovo's
 *
 */
public class DeNovoMatic {

	public static void run(String vpopFile, String fileOfBams, String outputDir, String ponVcf, double freqFilter, GATK gatk, MUTECT_RUN_TYPES type, int numThreads, int numSampleThreads, ANNOVCF annoVCF, String finalVcf, String tparams, Logger log) throws IllegalStateException {
		new File(outputDir).mkdirs();

		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.DENOVO, log);
		vpop.report();
		String matchFile = outputDir + ext.rootOf(vpopFile) + ".matchedbams.txt";
		if (!Files.exists(matchFile) || Files.countLines(matchFile, 0) != vpop.getSuperPop().get(VcfPopulation.OFFSPRING).size() * 2) {
			prepareMatchedBamFile(vpop, fileOfBams, matchFile, log);
		}

		// Mutect2.run(null, matchFile, callDir, ponVcf, gatk, type, numThreads, numSampleThreads, log);
		MutectTumorNormal[] results = Mutect2.callSomatic(matchFile, outputDir, ponVcf, gatk, null, null, null, numThreads, numSampleThreads, false, log);

		MergeFamResult[] resultsMerge = prepareResultsForFamilies(gatk, vpop, results, outputDir, numThreads, log);
		String[] finalFilesToMerge = new String[resultsMerge.length];

		WorkerHive<MergeFamResult> hive = new WorkerHive<DeNovoMatic.MergeFamResult>(numThreads, 10, log);
		hive.addCallables(resultsMerge);
		hive.execute(true);
		ArrayList<MergeFamResult> resultsDenovo = hive.getResults();
		for (int i = 0; i < resultsDenovo.size(); i++) {
			finalFilesToMerge[i] = resultsDenovo.get(i).getPotentialDenovoVcf();
		}

		String mergeDenovoOut = outputDir + ext.rootOf(vpopFile) + ".merge.denovo.vcf.gz";

		gatk.mergeVCFs(finalFilesToMerge, mergeDenovoOut, numThreads, false, log);
		String renamed = VCFOps.getAppropriateRoot(mergeDenovoOut, false) + ".renamed.vcf";
		if (!Files.exists(renamed)) {
			VCFTumorNormalOps.renameMergeVCF(mergeDenovoOut, renamed);
		}
		String annotatedVcf = GATK_Genotyper.annotateOnlyWithDefualtLocations(renamed, annoVCF, false, false, log);
		String mergeFinal = VCFOps.getAppropriateRoot(annotatedVcf, false) + ".merged.vcf.gz";
		if (finalVcf != null) {
			gatk.mergeVCFs(new String[] { annotatedVcf, finalVcf }, mergeFinal, numThreads, false, log);
			if (tparams != null) {
				Mutect2.runTally(tparams, FILTER_GENERATION_TYPE.EHQ_DNM, false, log, mergeFinal);// failure taken care of
				Mutect2.runTally(tparams, FILTER_GENERATION_TYPE.HQ_DNM, false, log, mergeFinal);
			}
		}
		// log.reportTimeInfo("Filtering " + annotatedVcf);
		// String annoFreq = VCFOps.getAppropriateRoot(annotatedVcf, false) + ".freq_" + freqFilter + ".func.vcf.gz";
		// filterByFreq(annotatedVcf, annoFreq, freqFilter, log);
		// log.reportTimeInfo("FIN");

	}

	private static class MergeFamResult implements Callable<MergeFamResult> {
		// private MutectTumorNormal[] famResults;
		private GATK gatk;
		private String[] vcfsToMerge;
		private String mergedVCF;
		private String potentialDenovoVcf;
		private String off;
		private HashSet<String> offCombo;
		private Logger log;

		public MergeFamResult(GATK gatk, MutectTumorNormal[] famResults, String mergedVCF, String[] vcfsToMerge, String off, String p1, String p2, Logger log) {
			super();
			this.gatk = gatk;
			this.vcfsToMerge = vcfsToMerge;
			// this.famResults = famResults;
			this.mergedVCF = mergedVCF;
			this.potentialDenovoVcf = VCFOps.getAppropriateRoot(mergedVCF, false) + ".denovo.vcf.gz";
			this.off = off;
			this.offCombo = new HashSet<String>();
			offCombo.add(off + ".variant");
			offCombo.add(off + ".variant2");
			// this.p1 = p1;
			// this.p2 = p2;
			this.log = log;
		}

		@Override
		public MergeFamResult call() throws Exception {
			if (!Files.exists(mergedVCF)) {
				gatk.mergeVCFs(vcfsToMerge, mergedVCF, 1, false, log);
			}
			scanForDenovo();
			return this;
		}

		public String getPotentialDenovoVcf() {
			return potentialDenovoVcf;
		}

		private void scanForDenovo() {
			VariantContextFilter filter = FilterNGS.generateFilter(FILTER_GENERATION_TYPE.TN, Double.NaN, false, log);
			
			if (!VCFOps.existsWithIndex(potentialDenovoVcf)) {
				VCFFileReader reader = new VCFFileReader(new File(mergedVCF), true);
				VariantContextWriter writer = VCFOps.initWriter(potentialDenovoVcf, VCFOps.DEFUALT_WRITER_OPTIONS, reader.getFileHeader().getSequenceDictionary());
				HashSet<VCFHeaderLine> newHeader = new HashSet<VCFHeaderLine>();
				// for (VCFHeaderLine vcfHeaderLine : reader.getFileHeader().getMetaDataInInputOrder()) {
				// newHeader.add(vcfHeaderLine);
				// }
				newHeader.addAll(reader.getFileHeader().getFormatHeaderLines());
				newHeader.addAll(reader.getFileHeader().getInfoHeaderLines());
				newHeader.addAll(reader.getFileHeader().getContigLines());
				newHeader.addAll(reader.getFileHeader().getOtherHeaderLines());
				newHeader.addAll(reader.getFileHeader().getFilterLines());

				VCFFormatHeaderLine hqInBoth = new VCFFormatHeaderLine("HQ_DNM", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Variant status for custom denovo filters in both MO -> off direction and FA -> OFF direction");
				VCFFormatHeaderLine ehqInBoth = new VCFFormatHeaderLine("EHQ_DNM", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Variant status for custom and mutect denovo filters in both MO -> off direction and FA -> OFF direction");

				VCFFormatHeaderLine hqNonTransmissionP1 = new VCFFormatHeaderLine("HQ_P1_NT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Variant  status for non-transmission of alleles filters in  P1 -> off direction");
				VCFFormatHeaderLine hqNonTransmissionP2 = new VCFFormatHeaderLine("HQ_P2_NT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Variant  status for non-transmission of alleles filters in  P2 -> off direction");

				newHeader.add(hqInBoth);
				newHeader.add(ehqInBoth);
				newHeader.add(hqNonTransmissionP1);
				newHeader.add(hqNonTransmissionP2);

				ArrayList<String> originalAtts = new ArrayList<String>();
				for (VCFFormatHeaderLine vcfFormatHeaderLine : reader.getFileHeader().getFormatHeaderLines()) {
					originalAtts.add(vcfFormatHeaderLine.getID());
					VCFFormatHeaderLine newFormatP1 = new VCFFormatHeaderLine(vcfFormatHeaderLine.getID() + "_P1", vcfFormatHeaderLine.isFixedCount() ? vcfFormatHeaderLine.getCount() : 1, vcfFormatHeaderLine.getType(), vcfFormatHeaderLine.getDescription());
					VCFFormatHeaderLine newFormatP2 = new VCFFormatHeaderLine(vcfFormatHeaderLine.getID() + "_P2", vcfFormatHeaderLine.isFixedCount() ? vcfFormatHeaderLine.getCount() : 1, vcfFormatHeaderLine.getType(), vcfFormatHeaderLine.getDescription());
					// if (!vcfFormatHeaderLine.getID().equals("GT") && !vcfFormatHeaderLine.getID().equals("AD")) {
					// newHeader.remove(vcfFormatHeaderLine);
					// }
					newHeader.add(newFormatP1);
					newHeader.add(newFormatP2);

				}
				HashSet<String> offFinal = new HashSet<String>();
				offFinal.add(off);
				VCFHeader header = new VCFHeader(newHeader, offFinal);

				writer.writeHeader(header);
				int numPossible = 0;
				int numTotal = 0;
				int numHQ = 0;
				int numEHQ = 0;
				for (VariantContext vc : reader) {
					numTotal++;
					VariantContext vcSubOff = VCOps.getSubset(vc, offCombo);
					Genotype g1 = vcSubOff.getGenotype(0);
					Genotype g2 = vcSubOff.getGenotype(1);
					VariantContextBuilder vBuilder = new VariantContextBuilder(vcSubOff);

					if (g1.sameGenotype(g2)) {// called somatic in both files
						VariantContextFilterPass pass1 = filter.filter(VCOps.getSubset(vcSubOff, g1.getSampleName(), VC_SUBSET_TYPE.SUBSET_STRICT));
						VariantContextFilterPass pass2 = filter.filter(VCOps.getSubset(vcSubOff, g2.getSampleName(), VC_SUBSET_TYPE.SUBSET_STRICT));
						GenotypeBuilder g1bBuilder = new GenotypeBuilder(g1);
						g1bBuilder.name(off);
						Hashtable<String, Object> map = new Hashtable<String, Object>();
						map.put("HQ_P1_NT", pass1.getTestPerformed().replaceAll(":", "_").replaceAll(" ", "_"));
						map.put("HQ_P2_NT", pass2.getTestPerformed().replaceAll(":", "_").replaceAll(" ", "_"));

						if (pass1.passed() && pass2.passed()) {
							map.put("HQ_DNM", true + "");
							numHQ++;
							String mutF1 = g1.getAnyAttribute(GENOTYPE_INFO.MUTECT_FILTERS.getFlag()).toString();
							String mutF2 = g2.getAnyAttribute(GENOTYPE_INFO.MUTECT_FILTERS.getFlag()).toString();

							if (mutF1.equals("[PASS]") && mutF2.equals("[PASS]")) {
								map.put("EHQ_DNM", true + "");
								numEHQ++;
							} else {
								map.put("EHQ_DNM", mutF1 + "_" + mutF2);
							}
						} else {
							map.put("HQ_DNM", pass1.passed() + "_" + pass2.passed());
						}

						for (String att : originalAtts) {
							if (g1.hasAnyAttribute(att)) {
								map.put(att + "_P1", g1.getAnyAttribute(att));
							}
							if (g2.hasAnyAttribute(att)) {
								map.put(att + "_P2", g2.getAnyAttribute(att));

							}
						}
						g1bBuilder.attributes(map);
						ArrayList<Genotype> gtypes = new ArrayList<Genotype>();
						gtypes.add(g1bBuilder.make());
						vBuilder.genotypes(gtypes);
						writer.add(vBuilder.make());
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
				log.reportTimeInfo("Detected " + numPossible + " possible denovos (" + numHQ + " HQ, " + numEHQ + " EHQ)  after scanning " + numTotal);

			}
		}
	}

	private static MergeFamResult[] prepareResultsForFamilies(GATK gatk, VcfPopulation vpop, MutectTumorNormal[] results, String outputDir, int numThreads, Logger log) {
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
			// if (!Files.exists(outputMerge)) {
			// gatk.mergeVCFs(vcfsToMerge, outputMerge, numThreads, false, log);
			// }
			String[] famMembers = vpop.getOffP1P2ForFam(fam);
			MergeFamResult mergeFamResult = new MergeFamResult(gatk, famResults, outputMerge, vcfsToMerge, famMembers[0], famMembers[1], famMembers[2], log);
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
		ANNOVCF annoVCF = null;
		String finalVCF = null;
		String tparams = null;
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
			} else if (args[i].startsWith(GATK_Genotyper.EXTRA_VCF_ANNOTATIONS)) {
				annoVCF = ANNOVCF.fromArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("finalVCF")) {
				finalVCF = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("tparams=")) {
				tparams = args[i].split("=")[1];
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
		GATK gatk = new GATK.Mutect(gatkLocation, referenceGenomeFasta, knownSnps, regions, cosmic, true, false, log);
		try {
			run(vpopFile, fileOfBams, outputDir, ponVCF, freqFilter, gatk, MUTECT_RUN_TYPES.CALL_SOMATIC, numthreads, numSampleThreads, annoVCF, finalVCF, tparams, log);
		} catch (IllegalStateException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	// private static void filterByFreq(String invcf, String outVcf, double freq, Logger log) {
	//
	// VCFFileReader reader = new VCFFileReader(invcf, true);
	// VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, outVcf, VCFOps.DEFUALT_WRITER_OPTIONS, log);
	// int total = 0;
	// int pass = 0;
	// for (VariantContext vc : reader) {
	// total++;
	// if (!vc.hasAttribute("PopFreqMax")) {
	// throw new IllegalArgumentException("Method expects PopFreqMax annotations for all variants");
	// }
	// String f = vc.getAttributeAsString("PopFreqMax", "-1");
	//
	// if (f.equals(".")) {
	// f = "-1";
	// }
	// double freqVc = Double.parseDouble(f);
	// if (freqVc < freq && VCOps.isHighModLowSNP_EFFImpact(vc)) {
	//
	// writer.add(vc);
	// pass++;
	// }
	// }
	// writer.close();
	// reader.close();
	// log.reportTimeInfo("Scanned " + total + " variants and " + pass + " passed freq threshold of " + freq);
	// }

}
