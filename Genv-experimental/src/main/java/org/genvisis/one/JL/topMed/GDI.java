package org.genvisis.one.JL.topMed;

import java.io.File;
import java.util.List;
import org.genvisis.CLI;
import org.genvisis.common.Logger;
import org.genvisis.one.JL.topMed.TOPMedUtils.GeneImpact;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Prototype to compute GDI using TOPMed
 */
public class GDI {

  private static void run(CLI c) {
    String outDir = c.get(CLI.ARG_OUTDIR);
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "gdi.log");
    String vcf = c.get(CLI.ARG_VCF);

    log.reportTimeInfo("computing GDI from " + vcf);
    VCFFileReader reader = new VCFFileReader(new File(vcf));

    CloseableIterator<VariantContext> iter = reader.iterator();
    while (iter.hasNext()) {
      VariantContext vc = iter.next();
      List<GeneImpact> geneImpacts = TOPMedUtils.getAllGeneImpacts(vc);
    }
    reader.close();

  }

  public static void main(String[] args) {
    CLI c = new CLI(GDI.class);

    c.addArg(CLI.ARG_VCF, CLI.DESC_VCF);
    //    c.addArg("bed","bed file containing ");

    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR);

    c.parseWithExit(args);

    run(c);
  }

}
