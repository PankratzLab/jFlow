package org.genvisis.seq.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import javax.jms.IllegalStateException;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.HEADER_COPY_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author lane0212 Quick summary of a region geared towards tumor normal
 */
public class TumorNormalSummary {

  private static final String[] BASE_OUT = new String[] {"CHROM", "POS", "ID", "REF", "FULL_ALT",
                                                         "ALT", "HIGH||MODERATE||LOW", "TN_PAIR",
                                                         "NORMAL_SAMPLE", "TUMOR_SAMPLE",
                                                         "NORMAL_GENOTYPE", "TUMOR_GENOTYPE",
                                                         "NORMAL_GQ", "TUMOR_GQ", "NORMAL_HAS_ALT",
                                                         "TUMOR_HAS_ALT", "MIN_GQ", "TN_MATCH"};

  private static void run(String vcf, String vpopFile, String outputDir, Segment seg, String name,
                          int buffer, Logger log) {
    VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.TUMOR_NORMAL, log);
    Set<String> all = new HashSet<String>();
    all.addAll(vpop.getSuperPop().get(VcfPopulation.TUMOR));
    all.addAll(vpop.getSuperPop().get(VcfPopulation.NORMAL));

    log.reportTimeInfo("Detected " + vpop.getSubPop().keySet().size() + " Tumor normal pairs");
    String outRoot = outputDir + ext.rootOf(vpop.getFileName()) + "_" + name;
    String outVCF = outRoot + ".vcf.gz";
    String outSummary = outRoot + ".summary.txt";
    VCFFileReader reader = new VCFFileReader(new File(vcf), true);
    VariantContextWriter writer = VCFOps.initWriter(outVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
                                                    reader.getFileHeader().getSequenceDictionary());
    VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
    Segment segBuffer = seg.getBufferedSegment(buffer);
    CloseableIterator<VariantContext> iter = reader.query(
                                                          Positions.getChromosomeUCSC(segBuffer.getChr(),
                                                                                      true),
                                                          segBuffer.getStart(),
                                                          segBuffer.getStop());
    try {
      PrintWriter writerSummary = new PrintWriter(new FileWriter(outSummary));
      String[][] annos = VCFOps.getAnnotationKeys(vcf, log);
      writerSummary.println(Array.toStr(BASE_OUT) + "\t" + Array.toStr(annos[1]));
      writerSummary.println(Array.toStr(BASE_OUT) + "\t" + Array.toStr(annos[0]));

      while (iter.hasNext()) {
        VariantContext vcFull = iter.next();
        VariantContext vc = VCOps.getSubset(vcFull, all, VC_SUBSET_TYPE.SUBSET_STRICT);
        if (!vc.isMonomorphicInSamples()) {
          writer.add(vc);
          for (String tnPair : vpop.getSubPop().keySet()) {
            Set<String> samps = vpop.getSubPop().get(tnPair);
            VariantContext vcTNPair = VCOps.getSubset(vc, samps);
            if (!vcTNPair.isMonomorphicInSamples()) {
              String tumor = null;
              String normal = null;
              for (String samp : samps) {
                if (vpop.getPopulationForInd(samp,
                                             RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.TUMOR)) {
                  tumor = samp;
                } else if (vpop.getPopulationForInd(samp,
                                                    RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.NORMAL)) {
                  normal = samp;
                } else {
                  writerSummary.close();
                  throw new IllegalStateException("Unknown types");
                }
              }
              // private static final String[] BASE_OUT = new String[] { "CHROM", "POS", "ID",
              // "REF", "ALT", "TN_PAIR", "NORMAL_GENOTYPE", "TUMOR_GENOTYPE", "TN_MATCH" };
              Genotype gTumor = vc.getGenotype(tumor);
              Genotype gNormal = vc.getGenotype(normal);

              StringBuilder builder = new StringBuilder();
              builder.append(vc.getContig());
              builder.append("\t" + vc.getStart());
              builder.append("\t" + vc.getID());
              builder.append("\t" + vc.getReference().getDisplayString());
              builder.append("\t" + vcFull.getAlternateAlleles());
              builder.append("\t" + vc.getAlternateAlleles());
              String impact = VCOps.getSNP_EFFImpact(vc);
              builder.append("\t" + (impact.equals("HIGH") || impact.equals("MODERATE")
                                     || impact.equals("LOW")));
              builder.append("\t" + tnPair);
              builder.append("\t" + normal);
              builder.append("\t" + tumor);
              builder.append("\t" + gNormal.toString());
              builder.append("\t" + gTumor.toString());
              builder.append("\t" + gNormal.getGQ());
              builder.append("\t" + gTumor.getGQ());
              builder.append("\t" + (gNormal.isCalled() && !gNormal.isHomRef()));
              builder.append("\t" + (gTumor.isCalled() && !gTumor.isHomRef()));

              builder.append("\t" + Math.min(gNormal.getGQ(), gTumor.getGQ()));
              builder.append("\t" + gTumor.sameGenotype(gNormal));
              builder.append("\t" + Array.toStr(VCOps.getAnnotationsFor(annos[0], vc, ".")));
              writerSummary.println(builder.toString());
            }
          }
        }
      }

      writerSummary.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + outSummary);
      log.reportException(e);
    }
    writer.close();
    reader.close();

  }

  public static void main(String[] args) {
    String vcf =
               "D:/data/Project_Tsai_21_25_26_28_spector/joint_genotypes_tsai_21_25_26_28_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.merge_ARIC.hg19_multianno.eff.gatk.anno_charge.sed1000g.vcf.gz";
    String outputDir = "D:/data/Project_Tsai_21_25_26_28_spector/TumorNormal/";
    String vpop = outputDir + "TN.vpop";

    Logger log = new Logger(outputDir + "TN.log");
    Segment[] segs = new Segment[] {new Segment("chr15:50714579-50795277"),
                                    new Segment("chr8:143543377-143628368")};
    String[] names = new String[] {"USP8", "BAI1"};
    int buffer = 300;
    for (int i = 0; i < names.length; i++) {
      run(vcf, vpop, outputDir, segs[i], names[i], buffer, log);
    }
  }
}

// String subsetVcf = outputDir + ext.rootOf(vpop.getFileName()) + ".vcf.gz";
// if (!Files.exists(subsetVcf)) {
// log.reportTimeInfo("Generating fast- query vcf");
// Set<String> all = new HashSet<String>();
// all.addAll(vpop.getSuperPop().get(VcfPopulation.TUMOR));
// all.addAll(vpop.getSuperPop().get(VcfPopulation.NORMAL));
// Hashtable<String, Set<String>> allHash = new Hashtable<String, Set<String>>();
// allHash.put(ext.rootOf(vpopFile), all);
// VcfPopulation allpop = new VcfPopulation(allHash, allHash, POPULATION_TYPE.ANY, log);
// String tmpOut = outputDir + "tmp.vpop";
// allpop.dump(tmpOut);
// VcfPopulation.splitVcfByPopulation(vcf, tmpOut, true, true, log);
// }
