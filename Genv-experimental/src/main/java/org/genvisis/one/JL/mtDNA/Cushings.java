package org.genvisis.one.JL.mtDNA;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.JL.ssh.SSH;
import org.genvisis.seq.analysis.GATK;
import org.genvisis.seq.analysis.SimpleTallyGene;
import org.genvisis.seq.analysis.TumorNormalSummary;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.mtdna.GenBankMtDNA;
import org.genvisis.seq.manage.mtdna.RCRS;
import org.genvisis.seq.manage.mtdna.VCFOpsMT;
import org.genvisis.seq.manage.mtdna.VCFOpsMT.MT_GENOME;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class Cushings {
	// VCF validation
	// java -jar /Users/Kitty/bin/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T
	// ValidateVariants -R /Volumes/Beta/ref/GRCh37_canon.fa -V
	// /Volumes/Beta/data/Cushings/mito/testRcrs/joint_genotypes_tsai_21_25_26_28_spector.chrM.rcrs.poly.disease.vcf
	// --warnOnErrors
	public static void main(String[] args) {

		boolean summarize = false;
		if (summarize) {
			String cVcf = "/Volumes/Beta/data/Cushings/mito/testRcrs/joint_genotypes_tsai_21_25_26_28_spector.chrM.vcf";
			String polyVCF = "/Volumes/Beta/data/Cushings/mito/testRcrs/polymorphismsMTVariants.vcf";
			String diseaseVCF = "/Volumes/Beta/data/Cushings/mito/testRcrs/diseaseMTVariants.vcf";
			String diseaseTrimVCF = ext.addToRoot(diseaseVCF, ".trim");
			String rcrsC = ext.addToRoot(cVcf, ".rcrs");
			String outDir = ext.parseDirectoryOfFile(cVcf);
			Logger log = new Logger(outDir + "log.log");
			GenBankMtDNA.fixGenBankVCF(diseaseVCF, diseaseTrimVCF, log);

			RCRS.writeRef(outDir, log);

			VCFOpsMT.convertHg19ToRCRS(cVcf, rcrsC, new Logger());

			GATK gatk = new GATK("/Users/Kitty/bin/GenomeAnalysisTK-3.6/", "/Volumes/Beta/ref/GRCh37_canon.fa", true,
					true, log);

			String outAnno = ext.addToRoot(rcrsC, ".poly");
			gatk.annotateWithAnotherVCF(rcrsC, polyVCF, outAnno, new String[] { "AF", "AC" }, "mtPoly", "MT", 1);

			String outAnnoD = ext.addToRoot(outAnno, ".disease");
			gatk.annotateWithAnotherVCF(outAnno, diseaseTrimVCF, outAnnoD,
					new String[] { "AF", "AC", "DiseaseStatus", "Disease" }, "mtPoly", "MT", 1);

			String finalOut = ext.addToRoot(outAnnoD, ".conv");
			VCFOpsMT.convertContigs(outAnnoD, finalOut, MT_GENOME.HG19, log);

			String remoteVcf = "/home/pankrat2/lanej/tmp/cushings/" + ext.removeDirectoryInfo(finalOut);
			// SSH.copyLocalToRemote(finalOut, remoteVcf, log);
			String annotateCommand = "java -jar /home/pankrat2/lanej/genvisis.jar one.JL.quickAnno vcf=" + remoteVcf;
			System.out.println(annotateCommand);
		} else {

			String annotated = "/home/pankrat2/lanej/tmp/cushings/joint_genotypes_tsai_21_25_26_28_spector.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.vcf";
			String annoVcf = "/Volumes/Beta/data/Cushings/mito/joint_genotypes_tsai_21_25_26_28_spector.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.vcf";
			// System.exit(1);

			String outputDir = "/Volumes/Beta/data/Cushings/mito/tumorNormal/";
			String vpop = outputDir + "TN.form.vpop";
			VcfPopulation tnVpop = VcfPopulation.load(vpop, POPULATION_TYPE.TUMOR_NORMAL, new Logger());
			VCFFileReader reader = new VCFFileReader(new File(annoVcf), true);
			HashMap<String, Integer> tumorGCounts = new HashMap<String, Integer>();

			HashMap<String, Integer> nonTumorGCounts = new HashMap<String, Integer>();

			for (VariantContext vc : reader) {
				if (vc.getStart() == 302) {
					for (Genotype g : vc.getGenotypes()) {
						String key = g.getGenotypeString();
						if (tnVpop.getTumorSamples().contains(g.getSampleName())) {
							if (tumorGCounts.containsKey(key)) {
								tumorGCounts.put(key, tumorGCounts.get(key) + 1);
							} else {
								tumorGCounts.put(key, 1);
							}
						} else {
							if (nonTumorGCounts.containsKey(key)) {
								nonTumorGCounts.put(key, nonTumorGCounts.get(key) + 1);
							} else {
								nonTumorGCounts.put(key, 1);
							}
						}
					}

				}
			}
			System.out.println(tumorGCounts.toString());
			System.out.println(nonTumorGCounts.toString());

			ArrayList<String> summary = new ArrayList<String>();
			summary.add("Genotype\tTumorCount\tNormalCount");

			for (String key : tumorGCounts.keySet()) {
				if (nonTumorGCounts.containsKey(key)) {
					summary.add(key + "\t" + tumorGCounts.get(key) + "\t" + nonTumorGCounts.get(key));
				} else {
					summary.add(key + "\t" + tumorGCounts.get(key) + "\t0");
				}
			}
			Files.writeIterable(summary, ext.parseDirectoryOfFile(vpop) + "D310Summary.txt");
			// System.exit(1);
			// SSH.copyRemoteToLocal(annoVcf, annotated, log);
			// SSH.copyRemoteToLocal(annoVcf + ".idx", annotated + ".idx", log);
			if (Files.exists(annoVcf)) {
				//
				SimpleTallyGene.run(annoVcf, new double[] { 1.2, .01, 0.0 });
			}
			TumorNormalSummary.main(new String[] { annoVcf });

		}

	}

	// SSH.runRemoteCommand(makeDir, log);
	// String annotateCommand = "java -jar /home/pankrat2/lanej/genvisis.jar
	// one.JL.quickAnno vcf=" + remoteVcf;
	// SSH.runRemoteCommand(annotateCommand, log);
	// String annoVcf = ext.addToRoot(finalOut, ".anno");
	// SSH.copyRemoteToLocal(remoteVcf + ".anno.vcf", annoVcf, log);
	// // =
	// //
	// "/Volumes/Beta/data/Cushings/mito/joint_genotypes_tsai_21_25_26_28_spector.chrM.posAdjust_-1.hg19_multianno.eff.gatk.sed1000g.posAdjust_1.disease.poly.vcf";
	// //

}
