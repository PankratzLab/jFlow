package org.genvisis.one.JL;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;

import org.genvisis.seq.analysis.VCFSimpleTally;
import org.genvisis.seq.analysis.mutMap.VCMut;
import org.genvisis.seq.analysis.mutMap.VCMut.MutInd;
import org.genvisis.seq.analysis.mutMap.VCMut.PARSE_METHOD;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.filesys.Segment;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class TestMutInd {

  public static void main(String[] args) {
    String vcf = "/Volumes/Beta/data/Cushings/joint_genotypes_tsai_21_25_26_28_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.merge_ARIC.hg19_multianno.eff.gatk.anno_charge.sed1000g.vcf.gz";
    String outputDir = "/Volumes/Beta/data/Cushings/ATM/mutMap/";
    new File(outputDir).mkdirs();
    Logger log = new Logger(outputDir + "log.log");
    VcfPopulation vpop = VcfPopulation.load("/Volumes/Beta/data/Cushings/genesets/CUSHING_FREQ_V3.vpop",
                                            POPULATION_TYPE.ANY, log);
    double maf = 0.01;
    run(vcf, outputDir, log, vpop, maf);

  }

  private static void run(String vcf, String outputDir, Logger log, VcfPopulation vpop,
                          double maf) {
    VariantContextFilter filter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {},
                                                           new VARIANT_FILTER_BOOLEAN[] {VARIANT_FILTER_BOOLEAN.FAILURE_FILTER},
                                                           new String[] {"Freq"},
                                                           new String[] {FilterNGS.getPopFreqFilterString(maf)},
                                                           log);

    Map<String, String> effmap = VCMut.getSNPEFFOntologyMap(new Logger());

    Segment seg = new Segment("chr11:108056992-108276393");

    for (String pop : vpop.getSuperPop().keySet()) {
      if ("CUSHING_FREQ_V3".equals(pop) || "ARIC".equals(pop)) {
        System.out.println(pop);
        String output = outputDir + VCFOps.getAppropriateRoot(vcf, true) + "_" + pop + ".txt";
        PrintWriter writer = Files.getAppropriateWriter(output);
        writer.println(ArrayUtils.toStr(MutInd.getMutIndMapHeader()));
        VCFFileReader reader = new VCFFileReader(new File(vcf));
        CloseableIterator<VariantContext> iter = reader.query(seg.getChromosomeUCSC(),
                                                              seg.getStart(), seg.getStop());
        while (iter.hasNext()) {
          VariantContext vc = iter.next();
          VCMut vcMut = VCMut.fromVariantContext(vc, PARSE_METHOD.SNP_EFF, effmap);
          if (vcMut.isHasFullInfo() && VCFSimpleTally.filterCHARGEAndTOPMed(vc, maf)) {
            List<MutInd> inds = vcMut.parseToInds(vpop.getSuperPop().get(pop), filter, -4);
            for (MutInd ind : inds) {
              writer.println(ArrayUtils.toStr(ind.getMutIndMapFormat()));

            }
          }

        }
        iter.close();
        reader.close();
        writer.close();
      }
    }
  }

}
