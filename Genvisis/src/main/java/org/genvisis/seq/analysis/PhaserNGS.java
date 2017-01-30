/**
 * 
 */
package org.genvisis.seq.analysis;

import java.io.File;
import java.util.ArrayList;

import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.seq.manage.BamOps;

/**
 * @author Kitty
 * 
 *         Wrapper for https://github.com/secastel/phaser -seems to be one of few that will phase
 *         indels: Paper at http://www.nature.com/articles/ncomms12817
 *
 *         Probably want pip/brew for install
 */
public class PhaserNGS {

	private static final String VCF = "vcf";
	private static final String BAM_DIR = "bamDir";
	private static final String PHASER_PY = "phaser";

	private PhaserNGS() {

	}

	private static void run(String phaserPy, String vcf, String bamDir, String outDir) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "phaser.log");
		String[] bams = Files.listFullPaths(bamDir, ".bam", false);
		// VCFOps.extractSegments( vcf, "/Volumes/Beta/ref/mica.targets.bed", 10, bamDir,
		// ext.parseDirectoryOfFile(vcf) + "extracted/", false, true, false, 1,
		// log);
		// System.exit(1);
		for (String bam : bams) {
			phase(phaserPy, bam, vcf, outDir, log);
		}
	}

	private static String phase(String phaserPy, String bam, String vcf, String outDir, Logger log) {
		ArrayList<String> command = new ArrayList<String>();

		// --bam ../bams/C018307-02_S129_L001_001.sorted.not.dedup.realigned.recal.bam --vcf
		// ../gatk/gatk.mica.nodDup.chr6.vcf --sample C018307-02 --paired_end 1 --mapq 255 --baseq 10
		// --o phaser_test_case
		String sampleName = BamOps.getSampleName(bam, log);
		String rootOut = outDir + sampleName + ".phase";

		command.add("python");
		command.add(phaserPy);

		command.add("--bam");
		command.add(bam);
		command.add("--vcf");
		command.add(vcf);
		command.add("--sample");
		command.add(sampleName);
		command.add("--paired_end");
		command.add(Integer.toString(0));
		command.add("--mapq");
		command.add(Integer.toString(20));
		command.add("--baseq");
		command.add(Integer.toString(0));
		command.add("--o");
		command.add(rootOut);
		command.add("--pass_only");
		command.add(Integer.toString(0));
		command.add("--include_indels");
		command.add(Integer.toString(1));
		command.add("--cc_threshold");
		command.add(Double.toString(0));
		command.add("--remove_dups");
		command.add(Integer.toString(0));

		String vcfOut = rootOut + ".vcf";
		CmdLine.runCommandWithFileChecks(	ArrayUtils.toStringArray(command), "", new String[] {vcf},
																			new String[] {vcfOut}, true, true, false, log);
		return vcfOut;

	}

	public static void main(String[] args) {
		CLI c = new CLI(PhaserNGS.class);
		c.addArgWithDefault(VCF, "vcf file for phasing", null);
		c.addArgWithDefault(BAM_DIR, "directory of bams", null);
		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, null);
		c.addArgWithDefault(PHASER_PY, "full path to phaser.py", null);


		c.parseWithExit(args);
		run(c.get(PHASER_PY), c.get(VCF), c.get(BAM_DIR), c.get(CLI.ARG_OUTDIR));
	}

}
