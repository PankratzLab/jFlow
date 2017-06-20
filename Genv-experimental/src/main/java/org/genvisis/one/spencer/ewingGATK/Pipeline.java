package org.genvisis.one.spencer.ewingGATK;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;
import org.genvisis.qsub.Qsub;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class Pipeline {

	private static final String BAM = ".bam";
	private static final String BAI = ".bai";
	private static final String PBS = ".pbs";

	private static final String QSUB = "qsub -A $MYGROUP -W group_list=$MYGROUP ";

	private static final int PANKRATZ_CORES = 24;
	private static final int PANKRATZ_MEM_GB = 252;

	private static final int SMALL_CORES = 24;
	private static final int SMALL_MEM_GB = 62;

	private static final Joiner filenameJoiner = Joiner.on("_");

	private static final String REF_GENO = "/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/WUSTL_RefGenome/all_sequences.fa";

	private final String bamDir;
	private final Set<String> ids;


	public Pipeline(String bamDir) {
		this.bamDir = bamDir;
		String[] bams = Files.list(bamDir, "", ".bam", true, false);
		ids = Arrays.stream(bams).map(ext::rootOf).filter(id -> !id.contains("."))
								.filter(id -> !id.contains("_")).collect(Collectors.toSet());
	}

	public void generatePipeline() {
		// Map<String, String> toIndex = indexBAMs();
		// Map<String, String> toRecal = runRecal();
		// Map<String, String> toCall = runHaplotypeCaller();

		// Set<String> pbsToIndex = Sets.newTreeSet();
		//
		// for (Map.Entry<String, String> indexEntry : toIndex.entrySet()) {
		// String id = indexEntry.getKey();
		// String pbs = indexEntry.getValue();
		//
		// pbsToIndex.add(QSUB + pbs);
		// }
		// Files.writeIterable(pbsToIndex, bamDir + "../pbsScripts/pbsToIndex");

		// Map<String, String> toRecal = runRecal();
		// Map<String, String> toCall = runHaplotypeCaller();

		// runHaplotypeCallerItasca();

		runCombineGVCFs();

	}

	private Map<String, String> indexBAMs() {
		Map<String, String> toIndex = Maps.newHashMap();
		for (String id : ids) {
			if (!fileExists(bamIndex(id))) {
				StringJoiner cmd = new StringJoiner("\n");
				cmd.add(cdCmd());
				cmd.add("");
				cmd.add("java -jar /home/pankrat2/public/bin/picard/picard.jar BuildBamIndex INPUT="
								+ bam(id));
				cmd.add("");
				// cmd.add(QSUB + recalPBS(id));
				Qsub.qsubGb(indexPBS(id), cmd.toString(), SMALL_MEM_GB / SMALL_CORES, 2, 1);
				toIndex.put(id, indexPBS(id));
			}
		}
		return toIndex;
	}

	private Map<String, String> runRecal() {
		final int threads = 8;
		Map<String, String> toRecal = Maps.newHashMap();
		for (String id : ids) {
			// if (!fileExists(recal1(id), recalPBS(id))) {
			if (true) {
				StringJoiner cmd = new StringJoiner("\n");
				cmd.add(cdCmd());
				cmd.add("");
				cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
					 .add("\t-T BaseRecalibrator \\")
					 .add("\t-R " + REF_GENO + " \\")
					 .add("\t-I " + bam(id) + " \\")
					 .add("\t-knownSites /home/pankrat2/public/bin/ref/dbsnp_138.b37.vcf \\")
					 .add("\t-knownSites /home/pankrat2/public/bin/ref/Mills_and_1000G_gold_standard.indels.b37.vcf \\")
					 .add("\t-o " + recal1(id) + " \\")
					 .add("\t-nct 4");
				cmd.add("");
				cmd.add("qsub -q batch " + callerPBS(id));
				// cmd.add(QSUB + callerPBS(id));
				// cmd.add("");
				// cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
				// .add("\t-T BaseRecalibrator \\")
				// .add("\t-R " + REF_GENO + " \\")
				// .add("\t-I " + bam(id) + " \\")
				// .add("\t-knownSites /home/pankrat2/public/bin/ref/dbsnp_138.b37.vcf \\")
				// .add("\t-knownSites
				// /home/pankrat2/public/bin/ref/Mills_and_1000G_gold_standard.indels.b37.vcf \\")
				// .add("\t-o " + recal2(id) + " \\")
				// .add("\t-BQSR " + recal1(id))
				// .add("\t-nct " + threads);
				// cmd.add("");
				// cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
				// .add("\t-T AnalyzeCovariates \\")
				// .add("\t-R " + REF_GENO + " \\")
				// .add("\t-before " + recal1(id) + " \\")
				// .add("\t-after " + recal2(id) + " \\")
				// .add("\t-plots " + recalPlots(id));
				Qsub.qsubGb(recalPBS(id), cmd.toString(), 22, 20, threads);
				toRecal.put(id, recalPBS(id));
			}
		}
		return toRecal;
	}

	private Map<String, String> runHaplotypeCaller() {
		final int threads = 4;
		Map<String, String> toCall = Maps.newHashMap();
		for (String id : ids) {
			// if (!fileExists(callerGvcf(id), callerPBS(id))) {
			if (!fileExists(callerGvcf(id), callerGvcfIndex(id))) {
				StringJoiner cmd = new StringJoiner("\n");
				cmd.add(cdCmd());
				cmd.add("");
				cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
					 .add("\t-T HaplotypeCaller \\")
					 .add("\t-R " + REF_GENO + " \\")
					 .add("\t-I " + bam(id) + " \\")
					 .add("\t-BQSR " + recal1(id) + " \\")
					 .add("\t--genotyping_mode DISCOVERY \\")
					 .add("\t-o " + callerGvcf(id) + " \\")
					 .add("\t--dbsnp /home/pankrat2/public/bin/ref/dbsnp_138.b37.vcf \\")
					 .add("\t-ERC GVCF \\")
					 .add("\t-nct " + threads);

				Qsub.qsubGb(callerPBS(id), cmd.toString(), 20, 40, threads);
				toCall.put(id, callerPBS(id));
			}
		}
		return toCall;
	}

	private void runHaplotypeCallerItasca() {
		final int threads = 4;
		Iterator<String> idIter = ids.iterator();
		Set<String> itascaIDs = HashVec.loadFileToHashSet("/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/bams/itascaIDs",
																											false);
		Set<String> pbsScripts = Sets.newHashSet();
		while (idIter.hasNext()) {
			Set<String> idsToBatch = Sets.newHashSetWithExpectedSize(3);
			int i = 0;
			while (i < 3 && idIter.hasNext()) {
				String id = idIter.next();
				if (!fileExists(callerGvcf(id)) && itascaIDs.contains(id)) {
					idsToBatch.add(id);
					i++;
				}
			}

			StringJoiner cmd = new StringJoiner("\n");
			cmd.add(cdCmd());
			for (String id : idsToBatch) {
				cmd.add("");
				cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
					 .add("\t-T HaplotypeCaller \\")
					 .add("\t-R " + REF_GENO + " \\")
					 .add("\t-I " + bam(id) + " \\")
					 .add("\t-BQSR " + recal1(id) + " \\")
					 .add("\t--genotyping_mode DISCOVERY \\")
					 .add("\t-o " + callerGvcf(id) + " \\")
					 .add("\t--dbsnp /home/pankrat2/public/bin/ref/dbsnp_138.b37.vcf \\")
					 .add("\t-ERC GVCF \\")
					 .add("\t-nct " + threads + " &");
				cmd.add("");
			}

			cmd.add("wait");
			String pbs = callerPBS(Joiner.on("_").join(idsToBatch));
			Qsub.qsubGb(pbs, cmd.toString(), 62, 40, 16);
			pbsScripts.add(QSUB + " -q sb " + pbs);
		}
		Files.writeIterable(pbsScripts,
												"/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/bams/runItascaBatch");

	}

	private void runSplitGVCFs() {

		final int threads = 1;
		List<String> commands = Lists.newArrayList();
		int cmdCount = 0;
		int commandScriptCount = 0;
		for (String id : ids) {
			boolean write = false;
			StringJoiner cmd = new StringJoiner("\n");
			cmd.add(cdCmd());
			cmd.add("");
			for (int i = 2; i < 26; i++) {

				// if (!fileExists(callerGvcf(id), callerPBS(id))) {
				String chr = Integer.toString(i);
				if (i == 23)
					chr = "X";
				if (i == 24)
					chr = "Y";
				if (i == 25)
					chr = "MT";
				if (!fileExists(splitGvcf(id, chr), splitGvcfIndex(id, chr))) {
					write = true;
					cmd.add("");
					cmd.add("java -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
						 .add("\t-T SelectVariants \\")
						 .add("\t-R " + REF_GENO + " \\")
						 .add("\t-V " + callerGvcf(id) + " \\")
						 .add("\t-o " + splitGvcf(id, chr) + " \\")
						 .add("\t-L " + chr);
					cmd.add("");

				}
			}
			if (write) {
				Qsub.qsub(splitPBS(id), cmd.toString(), 2500, 24, threads);
				commands.add("qsub " + (cmdCount < 48 ? "-q pankratz " : "") + splitPBS(id));
				if (cmdCount++ > 400) {
					Files.writeIterable(commands,
															"/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/pbsScripts/splitGVCFJobs_"
																				+ commandScriptCount);
					commands.clear();
					commandScriptCount++;
					cmdCount = 0;
				}
			}
		}
		Files.writeIterable(commands,
												"/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/pbsScripts/splitGVCFJobs_"
																	+ commandScriptCount);
	}

	private void runCombineGVCFs() {

		final int threads = 1;

		for (int i = 1; i < 26; i++) {

			String chr = intToChr(i);
			for (int batch = 1; batch <= 6; batch++) {
				String gvcfList = "/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/pbsScripts/chr"
													+ chr + "_gVCFs_batch_" + batch + ".list";
				// Files.writeIterable(ids.stream().map(id -> splitGvcf(id,
				// chr)).collect(Collectors.toList()),
				// gvcfList);

				String pbsScript = "/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/pbsScripts/CG_chr"
													 + chr + "_gVCFs_batch_" + batch + "_combineGVCFs.pbs";

				if (!fileExists(combinedGVCF(chr, batch), combinedGVCFIndex(chr, batch))) {

					StringJoiner cmd = new StringJoiner("\n");
					cmd.add(cdCmd());
					cmd.add("");
					cmd.add("java -Xmx14G -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
						 .add("\t-T CombineGVCFs \\")
						 .add("\t-R " + REF_GENO + " \\")
						 .add("\t-V " + gvcfList + " \\")
						 .add("\t-o " + combinedGVCF(chr, batch) + " \\")
						 .add("\t-L " + chr);
					cmd.add("");


					Qsub.qsubGb(pbsScript, cmd.toString(), 16, 96, threads);
				}

			}
		}

	}

	private static String combinedGVCF(String chr, int batch) {
		return "ES_combined_chr" + chr + "_batch_" + batch + "_presplit.g.vcf.gz";
	}

	private static String combinedGVCFIndex(String chr, int batch) {
		return combinedGVCF(chr, batch) + ".tbi";
	}

	private void runJointGenotyping() {


		final int threads = 1;
		for (int i = 1; i < 26; i++) {
			StringJoiner cmd = new StringJoiner("\n");
			cmd.add(cdCmd());
			cmd.add("");
			final String chr = intToChr(i);

			String chrGVCFs = "/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/chr" + chr
												+ "_gVCFs.list";

			Files.writeIterable(ids.stream()
														 .map(s -> "/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/gVCFs_byChr/"
																			 + s + "_chr" + chr + "_recalibrated.snps.indels.g.vcf.gz")
														 .collect(Collectors.toList()),
													chrGVCFs);

			if (!fileExists(jointVCF(chr), jointVCFIndex(chr))) {
				cmd.add("");
				cmd.add("java -Xmx15G -jar /home/pankrat2/public/bin/GATK_3.7/GenomeAnalysisTK.jar \\")
					 .add("\t-T GenotypeGVCFs \\")
					 .add("\t-R " + REF_GENO + " \\")
					 .add("\t-V " + chrGVCFs + " \\")
					 .add("\t-o " + jointVCF(chr) + " \\")
					 .add("\t-L " + chr);
				cmd.add("");
				Qsub.qsubGb(jointGenoPBS(chr), cmd.toString(), 19, 440, threads);
			}
		}

	}

	private String intToChr(int i) {
		String chr = Integer.toString(i);
		if (i == 23)
			chr = "X";
		else if (i == 24)
			chr = "Y";
		else if (i == 25)
			chr = "MT";
		return chr;
	}


	private String cdCmd() {
		return "cd " + bamDir;
	}

	private String bam(String id) {
		return id + BAM;
	}

	private String bamIndex(String id) {
		return id + BAI;
	}

	private String recal1(String id) {
		return id + ".recal_data.table";
	}

	private String recal2(String id) {
		return id + ".post_recal_data.table";
	}

	private String recalPlots(String id) {
		return id + ".recalibration_plots.pdf";
	}

	private String[] recalOuts(String id) {
		return new String[] {recal1(id), recal2(id), recalPlots(id)};
	}

	private String callerGvcf(String id) {
		return "../gVCFs/" + id + "_recalibrated.snps.indels.g.vcf.gz";
	}

	private String splitGvcf(String id, String chr) {
		return "../gVCFs_byChr/" + id + "_chr" + chr + "_recalibrated.snps.indels.g.vcf.gz";
	}

	private String jointVCF(String chr) {
		return "/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/vcfs/ES_JointGenotypes_chr"
					 + chr + "_splitGVCFs.vcf.gz";
	}

	private String jointVCFIndex(String chr) {
		return jointVCF(chr) + ".tbi";
	}

	private String callerGvcfIndex(String id) {
		return "../gVCFs/" + id + "_recalibrated.snps.indels.g.vcf.gz.tbi";
	}

	private String splitGvcfIndex(String id, String chr) {
		return splitGvcf(id, chr) + ".tbi";
	}

	private String indexPBS(String id) {
		return "../pbsScripts/I_" + id + "_index" + PBS;
	}

	private String recalPBS(String id) {
		return "../pbsScripts/R_" + id + "_recal" + PBS;
	}

	private String callerPBS(String id) {
		return "../pbsScripts/C_" + id + "_haplotypeCaller" + PBS;
	}

	private String splitPBS(String id) {
		return "../pbsScripts/S_" + id + "_splitGVCFs" + PBS;
	}

	private String jointGenoPBS(String chr) {
		return "../pbsScripts/J_chr" + chr + "_jointGenotyping" + PBS;
	}

	private boolean fileExists(String... filenames) {
		return Files.exists(bamDir, filenames);
	}

	public static void main(String[] args) {
		new Pipeline("/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/bams/").generatePipeline();
	}

}
