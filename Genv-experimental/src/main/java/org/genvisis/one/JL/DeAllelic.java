package org.genvisis.one.JL;

import java.io.File;
import org.genvisis.CLI;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.GATK_Genotyper;
import org.genvisis.seq.manage.VCFOps;

/**
 * Test ground for finally handling triallelic sites
 */
public class DeAllelic {

  //
  // bcftools norm -Ou -m -any input.vcf.gz |
  // bcftools norm -Ou -f human_g1k_v37.fasta |
  // bcftools annotate -Ob -x ID \
  // -I +'%CHROM:%POS:%REF:%ALT' |

  //

  private static void run(String vcf, String ref, String outDir) {
    Logger log = new Logger(ext.parseDirectoryOfFile(vcf) + "log.log");
    String outputVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".norm.vcf";
    new File(outDir).mkdirs();
    StringBuilder command = new StringBuilder();
    command.append("bcftools norm -Ov -m -any " + vcf + "| bcftools norm -Ov -f " + ref);
    command.append("| bcftools annotate -Ov -x ID -I +'%CHROM:%POS:%REF:%ALT'" + " > " + outputVCF);
    String bat = vcf + ".run";
    Files.write(command.toString(), bat);
    Files.chmod(bat);
    CmdLine.run(bat, ext.parseDirectoryOfFile(vcf));

    GATK_Genotyper.annotateOnlyWithDefualtLocations(outputVCF, null, PSF.Ext.DEFAULT_MEMORY_MB,
                                                    true, false, new Logger());
  }

  public static void main(String[] args) {
    CLI c = new CLI(DeAllelic.class);

    c.addArgWithDefault(CLI.ARG_VCF, CLI.DESC_VCF, "a.vcf");
    c.addArgWithDefault(CLI.ARG_REFERENCE_GENOME, CLI.DESC_REFERENCE_GENOME, "ref.fa");
    c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "na");

    c.parseWithExit(args);
    run(c.get(CLI.ARG_VCF), c.get(CLI.ARG_REFERENCE_GENOME), c.get(CLI.ARG_OUTDIR));

  }

}
