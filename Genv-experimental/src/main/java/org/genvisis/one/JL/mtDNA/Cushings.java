package org.genvisis.one.JL.mtDNA;

import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.JL.ssh.SSH;
import org.genvisis.seq.analysis.GATK;
import org.genvisis.seq.analysis.SimpleTallyGene;
import org.genvisis.seq.manage.mtdna.GenBankMtDNA;
import org.genvisis.seq.manage.mtdna.RCRS;
import org.genvisis.seq.manage.mtdna.VCFOpsMT;
import org.genvisis.seq.manage.mtdna.VCFOpsMT.MT_GENOME;

public class Cushings {
	// VCF validation
	// java -jar /Users/Kitty/bin/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T
	// ValidateVariants -R /Volumes/Beta/ref/GRCh37_canon.fa -V
	// /Volumes/Beta/data/Cushings/mito/testRcrs/joint_genotypes_tsai_21_25_26_28_spector.chrM.rcrs.poly.disease.vcf
	// --warnOnErrors
	public static void main(String[] args) {
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

		GATK gatk = new GATK("/Users/Kitty/bin/GenomeAnalysisTK-3.6/", "/Volumes/Beta/ref/GRCh37_canon.fa", true, true,
				log);

		String outAnno = ext.addToRoot(rcrsC, ".poly");
		gatk.annotateWithAnotherVCF(rcrsC, polyVCF, outAnno, new String[] { "AF", "AC" }, "mtPoly", "MT", 1);

		String outAnnoD = ext.addToRoot(outAnno, ".disease");
		gatk.annotateWithAnotherVCF(outAnno, diseaseTrimVCF, outAnnoD,
				new String[] { "AF", "AC", "DiseaseStatus", "Disease" }, "mtPoly", "MT", 1);

		String finalOut = ext.addToRoot(outAnnoD, ".conv");
		VCFOpsMT.convertContigs(outAnnoD, finalOut, MT_GENOME.HG19, log);

		String remoteVcf = "/home/pankrat2/lanej/tmp/cushings/" + ext.removeDirectoryInfo(finalOut);
		String makeDir = "mkdir -p " + ext.parseDirectoryOfFile(remoteVcf);

		SSH.runRemoteCommand(makeDir, log);
		SSH.copyLocalToRemote(finalOut, remoteVcf, log);
		String annotateCommand = "java -jar /home/pankrat2/lanej/genvisis.jar one.JL.quickAnno vcf=" + remoteVcf;
		SSH.runRemoteCommand(annotateCommand, log);
		String annoVcf = ext.addToRoot(finalOut, ".anno");
		SSH.copyRemoteToLocal(remoteVcf + ".anno.vcf", annoVcf, log);
		// =
		// "/Volumes/Beta/data/Cushings/mito/joint_genotypes_tsai_21_25_26_28_spector.chrM.posAdjust_-1.hg19_multianno.eff.gatk.sed1000g.posAdjust_1.disease.poly.vcf";
		//
		SimpleTallyGene.run(annoVcf);

	}

}
