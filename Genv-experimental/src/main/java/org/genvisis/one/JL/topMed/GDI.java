package org.genvisis.one.JL.topMed;

import java.io.File;
import java.util.List;
import org.genvisis.CLI;
import org.genvisis.common.Logger;
import org.genvisis.one.JL.topMed.TOPMedUtils.GeneImpact;
import org.genvisis.one.JL.topMed.TOPMedUtils.IMPACT;
import org.genvisis.seq.manage.VCFOps;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;

// : (i) “raw” GDI, calculated for each human gene by first multiplying each variant’s CADD score by
// the corresponding variant’s number of alleles in the 1,000 Genomes Project (a total of 610,160
// missense/nonsense/frameshift/in-frame indels/splice variants, with a MAF < 0.5, from a total of
// 20,243,313 alleles), then summing up all (CADD × allele count) products for one gene; (ii) the
// “CADD-normalized” gene-level model of accumulated mutational damage, calculated as in i, with
// each CADD score divided by the expected (median) CADD score of a variant with a similar allele
// frequency (Fig. S2A)
// http://www.pnas.org/content/112/44/13615.long
// http://www.pnas.org/content/112/44/13615.long#F2

/**
 * Prototype to compute GDI using TOPMed
 */
public class GDI {

  private static final String TRIM_TO_GENES = "trimToGenes";

  private static void trimToImpactVariants(CLI c) {
    String outDir = c.get(CLI.ARG_OUTDIR);
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "gdi.log");
    String vcf = c.get(CLI.ARG_VCF);

    log.reportTimeInfo("selecting variants from " + vcf);
    VCFFileReader reader = new VCFFileReader(new File(vcf));

    String outputVcf = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".func.vcf.gz";
    VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputVcf);
    builder.setReferenceDictionary(reader.getFileHeader().getSequenceDictionary());

    VariantContextWriter writer = builder.build();
    writer.writeHeader(reader.getFileHeader());

    CloseableIterator<VariantContext> iter = reader.iterator();
    int numTotal = 0;
    int numUsed = 0;

    while (iter.hasNext()) {
      numTotal++;
      VariantContext vc = iter.next();
      if (vc.getAlternateAlleles().size() != 1) {
        reader.close();
        throw new IllegalArgumentException("Must not have multiple alternate alleles");
      }
      List<GeneImpact> geneImpacts = TOPMedUtils.getAllGeneImpacts(vc);
      boolean useImpact = false;
      for (GeneImpact g : geneImpacts) {
        if (g.impact.ordinal() > IMPACT.LOW.ordinal()) {
          useImpact = true;
          break;
        }
      }
      if (useImpact) {

        double af = Double.parseDouble(vc.getAttributeAsString("AF", "0"));
        if (af < 0.5) {//with a MAF < 0.5, alt is annotated so we make sure it is minor
          writer.add(vc);
          numUsed++;
        }
      }
      if (numTotal % 1000000 == 0) {
        log.reportTimeInfo("processed " + numTotal + " variants, retained " + numUsed);
      }
    }
    log.reportTimeInfo("processed " + numTotal + " variants, retained " + numUsed);
    reader.close();
  }

  public static void main(String[] args) {
    CLI c = new CLI(GDI.class);

    c.addFlag(TRIM_TO_GENES, "trim the .vcf to variants with a gene annotation");
    c.addArg(CLI.ARG_VCF, CLI.DESC_VCF);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR);

    c.parseWithExit(args);

    if (c.has(TRIM_TO_GENES)) {
      trimToImpactVariants(c);
    }
  }

}
