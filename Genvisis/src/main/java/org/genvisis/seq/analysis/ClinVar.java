/**
 * 
 */
package org.genvisis.seq.analysis;

import java.io.File;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.seq.analysis.mutMap.VCMut;
import org.genvisis.seq.manage.CreateNonSiteOnlyVcf;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.core.CLI;
import org.pankratzlab.shared.filesys.Segment;
import org.genvisis.seq.manage.VCOps;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import scala.actors.threadpool.Arrays;

/**
 * Class for processing ClinVar data wget
 * ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar//vcf_GRCh37/clinvar_20170530.vcf.gz wget
 * ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar//vcf_GRCh37/clinvar_20170530.vcf.gz.tbi
 * https://www.ncbi.nlm.nih.gov/clinvar/docs/faq/ Why doesn't the VCF file contain all the data in
 * the XML file? ClinVar's VCF files are currently limited to records that have been assigned rs# in
 * dbSNP. Thus there may be two types of gaps: The variant is not in scope for dbSNP (i.e. the
 * length is greater than 50 bp or the exact location is not known). The variant has not yet been
 * processed by dbSNP. The ClinVar VCF and related VCF set can be retrieved from ClinVar's ftp site:
 */
public class ClinVar {

  /// types of CLINSIGss
  // ##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Variant Clinical Significance, 0 -
  // Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic,
  // 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other">

  private enum CLNSIG {

    UNCERTAIN_SIGNIFICANCE(0),
    NOT_PROVIDED(1),
    BENIGN(2),
    LIKELY_BENIGN(3),
    LIKELY_PATHOGENIC(4),
    PATHOGENIC(5),
    DRUG_RESPONSE(6),
    HISTOCOMPATIBILITY(7),
    OTHER(255);

    private int flag;

    private CLNSIG(int flag) {
      this.flag = flag;
    }

    private static CLNSIG getCLNSIG(int flag) {
      for (CLNSIG c : CLNSIG.values()) {
        if (flag == c.flag) {
          return c;
        }
      }

      throw new IllegalArgumentException("Invalid flag " + flag);
    }

  }
  private enum CLNSIG_SETS {
    @SuppressWarnings("unchecked")
    PATH_SET(new HashSet<CLNSIG>(Arrays.asList(new CLNSIG[] {CLNSIG.PATHOGENIC}))),
    @SuppressWarnings("unchecked")
    PATH_LIKELY_SET(new HashSet<CLNSIG>(Arrays.asList(new CLNSIG[] {CLNSIG.PATHOGENIC,
                                                                    CLNSIG.LIKELY_PATHOGENIC})));

    private HashSet<CLNSIG> set;

    private CLNSIG_SETS(HashSet<CLNSIG> set) {
      this.set = set;
    }
  }

  private static Set<CLNSIG> getCLNSIGs(VariantContext vc) {
    String[] anno = VCOps.getAnnotationsFor(new String[] {"CLNSIG"}, vc, ".")[0].split("\\|");
    HashSet<CLNSIG> clnsigs = new HashSet<>();
    for (String a : anno) {
      try {
        clnsigs.add(CLNSIG.getCLNSIG(Integer.parseInt(a)));

      } catch (NumberFormatException nfe) {
        String[] tmp = a.replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(" ", "").split(",");
        for (String t : tmp) {
          clnsigs.add(CLNSIG.getCLNSIG(Integer.parseInt(t)));
        }
      }
    }
    return clnsigs;
  }

  private static void run(String vcf, Segment seg, String outputDir, double[] mafs) {
    new File(outputDir).mkdirs();
    Logger log = new Logger(outputDir + "clinvar.log");
    String subsetHG19Vcf = outputDir + VCFOps.getAppropriateRoot(vcf, true) + ".clinvar."
                           + seg.getUCSClocation().replaceAll(":", "_").replaceAll("-", "_")
                           + "vcf.gz";
    log.reportTimeInfo("writing to " + subsetHG19Vcf);
    subsetVcf(vcf, subsetHG19Vcf, seg, log);

    String nonSiteOnlyVcf = CreateNonSiteOnlyVcf.createNonSiteOnlyVcf(subsetHG19Vcf, false, true);

    log.reportTimeInfo("non-site only vcf created at " + nonSiteOnlyVcf);
    String annotatedVCF = GATK_Genotyper.annotateOnlyWithDefualtLocations(nonSiteOnlyVcf, null,
                                                                          PSF.Ext.DEFAULT_MEMORY_MB,
                                                                          true, false, log);
    String preppedVcf = outputDir + "clinvar.init."
                        + seg.getUCSClocation().replaceAll(":", "_").replaceAll("-", "_")
                        + ".vcf.gz";
    Files.copyFileUsingFileChannels(new File(annotatedVCF), new File(preppedVcf), log);
    Files.copyFileUsingFileChannels(new File(annotatedVCF + ".tbi"), new File(preppedVcf + ".tbi"),
                                    log);
    annotatedVCF = preppedVcf;
    log.reportTimeInfo("annotated vcf created at " + annotatedVCF);

    Hashtable<String, Set<String>> fake = new Hashtable<>();
    fake.put("FAKE", new HashSet<>());
    fake.get("FAKE").add(VCFOps.getSamplesInFile(annotatedVCF)[0]);

    VcfPopulation fakeVpop = new VcfPopulation(fake, fake, POPULATION_TYPE.ANY, new Logger());

    fakeVpop.dump(annotatedVCF + ".vpop");
    Files.write(seg.getChromosomeUCSC() + "\t" + seg.getStart() + "\t" + seg.getStop(),
                annotatedVCF + ".seq");
    VCFOps.extractSegments(annotatedVCF, annotatedVCF + ".seq", 100, null, outputDir, false, true,
                           true, 1, log);

    for (CLNSIG_SETS set : CLNSIG_SETS.values()) {
      String filt = filterByClnsigSet(annotatedVCF, set, log);
      VCFOps.extractSegments(filt, annotatedVCF + ".seq", 100, null, outputDir, false, true, true,
                             1, log);
      for (double maf : mafs) {
        VCMut.run(filt, outputDir, annotatedVCF + ".vpop", maf, false, seg);
      }
    }

  }

  private static void subsetVcf(String inputVCF, String outputVCF, Segment seg, Logger log) {
    VCFFileReader reader = new VCFFileReader(new File(inputVCF), true);

    VariantContextWriter writer = VCFOps.initWriter(outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
                                                    VCFOps.getSequenceDictionary(reader));

    VCFHeader newVCFHeader = reader.getFileHeader();
    SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome(GENOME_BUILD.HG19,
                                                                      log).getIndexedFastaSequenceFile()
                                                                          .getSequenceDictionary();
    newVCFHeader.setSequenceDictionary(samSequenceDictionary);
    writer.writeHeader(newVCFHeader);
    CloseableIterator<VariantContext> iterator = reader.query(seg.getChromosomeUCSC()
                                                                 .replaceAll("chr", ""),
                                                              seg.getStart(), seg.getStop());
    int num = 0;
    while (iterator.hasNext()) {
      num++;
      VariantContext vc = iterator.next();
      VariantContextBuilder builder = new VariantContextBuilder(vc);
      builder.chr(seg.getChromosomeUCSC());
      writer.add(builder.make());
    }
    reader.close();
    writer.close();
    log.reportTimeInfo(num + " variants written to " + outputVCF);
  }

  private static String filterByClnsigSet(String vcf, CLNSIG_SETS set, Logger log) {
    VCFFileReader reader = new VCFFileReader(new File(vcf), true);
    String output = VCFOps.getAppropriateRoot(vcf, false) + set.toString() + ".vcf.gz";
    VariantContextWriter writer = VCFOps.initWriter(output, VCFOps.DEFUALT_WRITER_OPTIONS,
                                                    VCFOps.getSequenceDictionary(reader));

    VCFHeader newVCFHeader = reader.getFileHeader();

    writer.writeHeader(newVCFHeader);
    int num = 0;
    for (VariantContext vc : reader) {
      Set<CLNSIG> cur = getCLNSIGs(vc);
      for (CLNSIG c : cur) {
        if (set.set.contains(c)) {
          writer.add(vc);
          num++;
          break;
        }
      }
    }

    reader.close();
    writer.close();
    log.reportTimeInfo(num + " variants written to " + output);
    return output;
  }

  /**
   * @param args
   */
  public static void main(String[] args) {

    CLI c = new CLI(ClinVar.class);
    c.addArgWithDefault("clinvarVCF", "ClinVarVCF", "clinvar_20170530.vcf.gz");
    c.addArgWithDefault("segment", "segment of interest", "chr11:108056992-108276393");
    c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "/clinVar");
    c.addArgWithDefault("mafs", "comma-delimited maf values to use", "0, .001, .01, .05, 1.2");
    c.parseWithExit(args);
    run(c.get("clinvarVCF"), new Segment(c.get("segment")), c.get(CLI.ARG_OUTDIR),
        ArrayUtils.toDoubleArray(c.get("mafs").split(",")));
  }

}
