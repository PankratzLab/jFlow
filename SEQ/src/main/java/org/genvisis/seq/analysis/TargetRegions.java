package org.genvisis.seq.analysis;

import java.io.File;
import java.io.PrintWriter;
import java.util.Set;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCOps.ALT_ALLELE_CONTEXT_TYPE;
import org.genvisis.seq.manage.VCOps.GENOTYPE_INFO;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.shared.filesys.LocusSet;
import org.pankratzlab.shared.filesys.Positions;
import org.pankratzlab.shared.filesys.Segment;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class TargetRegions<T extends Segment> {

  private static final String[] TO_REPORT = new String[] {"SNPEFF_GENE_NAME", "SNPEFF_EFFECT",
                                                          "SNPEFF_IMPACT", "AAChange.refGene",
                                                          "SNPEFF_EXON_ID", "culprit", "snp138",
                                                          "esp6500si_all", "g10002014oct_all"};

  private final String vcfFile;
  private final LocusSet<T> targetRegions;
  private final VcfPopulation vpop;
  private final Logger log;

  public TargetRegions(String vcfFile, LocusSet<T> targetRegions, VcfPopulation vpop, Logger log) {
    super();
    this.vcfFile = vcfFile;
    this.targetRegions = targetRegions;
    this.vpop = vpop;
    this.log = log;
  }

  public void summarizeRegions(String fullPathToOutput, String toMatchVCF,
                               String[] toMatchAnnotations) {
    VCFFileReader reader = new VCFFileReader(new File(vcfFile), true);
    T[] regions = targetRegions.getLoci();

    String[] subpop = vpop.getSubPop().keySet()
                          .toArray(new String[vpop.getSubPop().keySet().size()]);
    try {
      PrintWriter writer = Files.openAppropriateWriter(fullPathToOutput);
      writer.print("CHR\tStart\tStop\tRef\tAlt\tFILTER");
      for (String element : subpop) {
        writer.print("\t" + element + "_AAC");
        writer.print("\t" + element + "_SamplesWithAlt");
        writer.print("\t" + element + "_NumExcluded");

      }
      writer.print("\t" + ArrayUtils.toStr(TO_REPORT));
      if (toMatchVCF != null && toMatchAnnotations != null) {
        for (String toMatchAnnotation : toMatchAnnotations) {
          writer.print("\t" + ext.rootOf(toMatchVCF) + "_" + toMatchAnnotation);
        }
      }
      for (String element : subpop) {
        writer.print("\t" + element + "_AVG_DP");
        writer.print("\t" + element + "_AVG_GQ");
      }
      writer.println();
      for (T region : regions) {
        CloseableIterator<VariantContext> cIterator = reader.query(Positions.getChromosomeUCSC(region.getChr(),
                                                                                               true),
                                                                   region.getStart(),
                                                                   region.getStop());

        while (cIterator.hasNext()) {
          VariantContext vc = cIterator.next();
          writer.print(Positions.getChromosomeUCSC(VCOps.getSegment(vc).getChr(), true));
          writer.print("\t" + VCOps.getSegment(vc).getStart());
          writer.print("\t" + VCOps.getSegment(vc).getStop());
          writer.print("\t" + vc.getReference().getBaseString());
          writer.print("\t" + vc.getAlternateAlleles().toString());
          writer.print("\t" + vc.getFilters().toString());

          for (String element : subpop) {
            writer.print("\t" + VCOps.getAAC(vc, vpop.getSubPop().get(element)));
            VariantContext vcAlt = VCOps.getAltAlleleContext(VCOps.getSubset(vc,
                                                                             vpop.getSubPop()
                                                                                 .get(element)),
                                                             null, null,
                                                             ALT_ALLELE_CONTEXT_TYPE.ALL, log);
            Set<String> tmpSamps = vcAlt.getSampleNames();
            writer.print("\t");
            int index = 0;
            for (String altSamp : tmpSamps) {
              VariantContext curContext = VCOps.getSubset(vcAlt, altSamp,
                                                          VC_SUBSET_TYPE.SUBSET_STRICT);
              writer.print((index > 0 ? "|" : "") + altSamp + ":GQ="
                           + curContext.getGenotype(0).getGQ() + ":AD="
                           + ArrayUtils.toStr(curContext.getGenotype(0).getAD(), ","));
              index++;
            }
            int numExcluded = 0;
            for (String altSamp : vcAlt.getSampleNames()) {
              if (vpop.getPopulationForInd(altSamp,
                                           RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.EXCLUDE)) {
                numExcluded++;
              }
            }
            writer.print("\t" + numExcluded);
          }
          writer.print("\t" + ArrayUtils.toStr(VCOps.getAnnotationsFor(TO_REPORT, vc, ".")));
          if (toMatchVCF != null && toMatchAnnotations != null) {
            VariantContext vcMatch = VCFOps.lookupExactVariant(toMatchVCF, vc, log);
            if (vcMatch == null) {
              for (String toMatchAnnotation : toMatchAnnotations) {
                writer.print("\tNA");
              }
            } else {
              writer.print("\t" + ArrayUtils.toStr(VCOps.getAnnotationsFor(toMatchAnnotations,
                                                                           vcMatch, ".")));
            }
          }
          for (String element : subpop) {
            writer.print("\t" + VCOps.getAvgGenotypeInfo(vc, vpop.getSubPop().get(element),
                                                         GENOTYPE_INFO.DP, log));
            writer.print("\t" + VCOps.getAvgGenotypeInfo(vc, vpop.getSubPop().get(element),
                                                         GENOTYPE_INFO.GQ, log));
          }

          writer.println();
        }
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + fullPathToOutput);
      log.reportException(e);
    }
    reader.close();
  }

  public static void test() {
    String vcf = "D:/data/Project_Tsai_21_25_26_spector/joint_genotypes_tsai_21_25_26_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf.gz";
    String vpopFile = "D:/data/Project_Tsai_21_25_26_spector/candidateGenes/USP8/usp8.vpop";
    String toMatchVCF = "D:/data/CHARGE/CHARGE_MAFS/charge_fibrinogen_mafs_and_counts.xln.hg19_multianno.eff.gatk.sed.vcf";
    String[] toMatchAnnotations = new String[] {"MAF_blacks", "MAF_whites"};
    String output = ext.parseDirectoryOfFile(vpopFile) + "targetRegions.txt";
    Logger log = new Logger(ext.parseDirectoryOfFile(vpopFile) + "target.log");
    VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
    vpop.report();
    Segment[] segs = new Segment[] {new Segment("chr15:50714579-50795277")};
    LocusSet<Segment> set = new LocusSet<Segment>(segs, true, log) {

      /**
       *
       */
      private static final long serialVersionUID = 1L;
    };
    TargetRegions<Segment> targetRegions = new TargetRegions<>(vcf, set, vpop, log);
    targetRegions.summarizeRegions(output, toMatchVCF, toMatchAnnotations);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "TargetRegions.dat";
    // String logfile = null;
    // Logger log;

    String usage = "\n" + "seq.analysis.TargetRegions requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      }
      // else if (args[i].startsWith("log=")) {
      // logfile = args[i].split("=")[1];
      // numArgs--;
      // }
      else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      // log = new Logger(logfile);
      test();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
