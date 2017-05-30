package org.genvisis.one.spencer.ewingGATK;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;
import org.genvisis.qsub.Qsub;

import com.google.common.base.Joiner;
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
					idsToBatch.add(idIter.next());
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
					 .add("\t-nct " + threads + "&");
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

	private String callerGvcfIndex(String id) {
		return "../gVCFs/" + id + "_recalibrated.snps.indels.g.vcf.gz.tbi";
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

	private boolean fileExists(String... filenames) {
		return Files.exists(bamDir, filenames);
	}

	public static void main(String[] args) {
		new Pipeline("/panfs/roc/pankrat2-SpectorEwingSarcoma/SpectorEwingSarcoma/bams/").generatePipeline();
	}

}
