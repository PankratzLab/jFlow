package org.genvisis.one.JL.topMed;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;
import java.util.StringJoiner;
import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.one.JL.topMed.TOPMedUtils.GeneImpact;
import org.genvisis.one.JL.topMed.TOPMedUtils.IMPACT;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;
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

// For each human gene g with n minor alleles (missense/nonsense/frameshift/in-frame indels/splice,
// with MAF < 0.5), we calculated GDIg, the cumulative predicted damage to the exonic regions of the
// gene, by multiplying the CADD raw C score by the count f of each allele A and dividing the
// expected CE score for a variant with a corresponding allele frequency (calculated by the 1,000
// Genomes Project Phase 3 median C score for each allele frequency slot; Fig. S2A) and then summing
// the results: GDIg=∑nA=1(CACE)fA
// . We calculated a homogenized Phred I-score for each metric, indicating the ranking of the gene
// of interest i relative to all other human genes (T = 19,558 genes in our analyses), with values
// ranging from 0 (lowest Phred score: human gene with lowest GDI) to 42.91 (highest Phred score:
// most damaged human gene): Ii=−10[log10(i/T)]
// . The GDI model described had the highest performance, in three alternative forms: (i) raw GDI,
// as above but without expected CADD score normalization: GDIg=∑nA=1CAfA
// ; (ii) observed/expected GDI-normalized model, with raw GDI calculated by dividing the observed
// CADD score values CO by the expected CADD score CE, determined by calculating the 1,000 Genomes
// Project Phase 3 median C score for each allele frequency slot: GDIg=(∑nA=1COfA)/(∑nA=1CEfA)
// ; and (iii) gene size-normalized GDI, where the raw GDI score was normalized by dividing by the
// CDS length L of the gene’s canonical transcript: GDIg=(∑nA=1COfA)/Lg
// .
//

// http://www.pnas.org/content/112/44/13615.long
// http://www.pnas.org/content/112/44/13615.long#F2

/**
 * Prototype to compute GDI using TOPMed
 */
public class GDI {

  private static final String TRIM_TO_GENES = "trimToGenes";
  private static final String COMPUTE_RAW_GDI = "computeRawGDI";
  private static final String EXTRACT_GDI_INFO = "extractGDIinfo";

  private static void trimToImpactVariants(CLI c) {
    String outDir = c.get(CLI.ARG_OUTDIR);
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "gdi.log");
    String vcf = c.get(CLI.ARG_VCF);

    String outputVcf = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".func.vcf.gz";
    if (!Files.exists(outputVcf)) {
      log.reportTimeInfo("selecting variants from " + vcf);

      VCFFileReader reader = new VCFFileReader(new File(vcf));

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
  }

  private static final String[] BASE = new String[] {"CHROM", "START", "END", "ID", "REF", "ALT"};
  private static final String[] ANNOS = new String[] {"CADD_raw", "CADD_phred",
                                                      "CADD_raw_rankscore", "AF", "AC",
                                                      "SNPEFF_GENE_NAME", "SNPEFF_IMPACT"};

  private static void extractGDIComponents(CLI c) {

    String outDir = c.get(CLI.ARG_OUTDIR);
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "gdi.log");
    String vcf = c.get(CLI.ARG_VCF);

    log.reportTimeInfo("computing GDI using variants from " + vcf);

    try (VCFFileReader reader = new VCFFileReader(new File(vcf))) {
      String outputTmp = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".cadd.af.txt.gz";
      CloseableIterator<VariantContext> iter = reader.iterator();
      int numTotal = 0;
      int numUsed = 0;
      int numSkipped = 0;

      PrintWriter writer = Files.getAppropriateWriter(outputTmp);
      writer.println(ArrayUtils.toStr(BASE) + "\t" + ArrayUtils.toStr(ANNOS));

      while (iter.hasNext()) {
        numTotal++;
        VariantContext vc = iter.next();
        if (!vc.isFiltered()) {

          try {
            Double.parseDouble(vc.getAttributeAsString("CADD_raw", "."));
            StringJoiner out = new StringJoiner("\t");
            out.add(vc.getContig());
            out.add(Integer.toString(vc.getStart()));
            out.add(Integer.toString(vc.getEnd()));
            out.add(vc.getID());
            out.add(vc.getReference().getBaseString());
            out.add(vc.getAlternateAlleles().get(0).getBaseString());

            writer.println(out.toString() + "\t"
                           + ArrayUtils.toStr(VCOps.getAnnotationsFor(ANNOS, vc, ".")));
            numUsed++;
          } catch (NumberFormatException nfe) {
            numSkipped++;
          }

          if (numTotal % 10000 == 0) {
            log.reportTimeInfo("processed " + numTotal + " variants, retained " + numUsed
                               + " skipped " + numSkipped + " for missing CADD");
          }
        }
      }
    }
  }

  public static void main(String[] args) {
    CLI c = new CLI(GDI.class);

    c.addFlag(TRIM_TO_GENES, "trim the .vcf to variants to be used in GDI scores");
    c.addFlag(EXTRACT_GDI_INFO, "extract components of the GDI score ");

    c.addArg(CLI.ARG_VCF, CLI.DESC_VCF);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR);

    c.parseWithExit(args);

    if (c.has(TRIM_TO_GENES)) {
      trimToImpactVariants(c);
    } else if (c.has(EXTRACT_GDI_INFO)) {
      extractGDIComponents(c);
    }

  }

}
