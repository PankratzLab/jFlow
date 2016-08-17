package org.genvisis.seq.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;

import org.genvisis.common.Logger;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.genvisis.seq.manage.VCOps.GENOTYPE_INFO;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Consolidating some common operations for tumor normal calling, summarizing filtering, and etc
 *
 */
public class VCFTumorNormalOps {
  public static class TNSample {
    private final String tumorSample;
    private final String normalSample;
    private final String tumorBam;
    private final String normalBam;

    public TNSample(String tumorSample, String normalSample, String tumorBam, String normalBam) {
      super();
      this.tumorSample = tumorSample;
      this.normalSample = normalSample;
      this.tumorBam = tumorBam;
      this.normalBam = normalBam;
    }

    public String getNormalBam() {
      return normalBam;
    }

    public String getNormalSample() {
      return normalSample;
    }

    public String getTumorBam() {
      return tumorBam;
    }

    public String getTumorSample() {
      return tumorSample;
    }

  }

  /**
   * 
   * NOTE, we also add "variant level" annotations to the genotypes format to preserve this
   * information after merging<br>
   * NOTE, really designed for mutect2
   * 
   * @param vcf vcf to rename
   * @param tumorSamp the tumor sample
   * @param tumorDef tumorSamp will replace this sample
   * @param normalSamp the normal sample
   * @param normalDef normalSamp will replace this sample
   * @param output output vcf
   * @param log
   */

  private static final String NORMAL_TAG = "_NORMAL";

  public static TNSample[] matchSamples(String[] bamFiles, VcfPopulation vpop, Logger log) {
    if (vpop.getType() != POPULATION_TYPE.TUMOR_NORMAL) {
      throw new IllegalArgumentException("Vpop must be " + POPULATION_TYPE.TUMOR_NORMAL);
    }

    Hashtable<String, String> all = new Hashtable<String, String>();
    for (String bamFile : bamFiles) {
      all.put(BamOps.getSampleName(bamFile), bamFile);
    }
    ArrayList<TNSample> tnSamples = new ArrayList<TNSample>();
    for (String tnPair : vpop.getSubPop().keySet()) {
      Set<String> samps = vpop.getSubPop().get(tnPair);
      String tumor = null;
      String normal = null;
      for (String samp : samps) {
        if (vpop.getPopulationForInd(samp, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.TUMOR)) {
          tumor = samp;
        } else if (vpop.getPopulationForInd(samp, RETRIEVE_TYPE.SUPER)[0]
            .equals(VcfPopulation.NORMAL)) {
          normal = samp;
        } else {
          throw new IllegalArgumentException("Unknown types");
        }
      }

      if (!all.containsKey(tumor) || !all.containsKey(normal)) {
        throw new IllegalArgumentException(
            "Could not find bam file for Tumor " + tumor + " or for Normal " + normal);
      } else {
        TNSample tSample = new TNSample(tumor, normal, all.get(tumor), all.get(normal));
        tnSamples.add(tSample);
      }
    }
    // if (analysisBams.size() < bamFiles.length) {
    // Files.writeList(Array.toStringArray(analysisBams), outputDir + ext.rootOf(vpop.getFileName()
    // + ".analysis.bams.txt"));
    // }

    log.reportTimeInfo("Matched " + tnSamples.size() + " tumor normals with bams ");
    return tnSamples.toArray(new TNSample[tnSamples.size()]);
  }

  private static Genotype rename(Genotype g, String newName) {
    GenotypeBuilder builder = new GenotypeBuilder(g);
    builder.name(newName);
    return builder.make();
  }

  public static void renameAndTransferInfo(String vcf, String tumorSamp, String tumorDef,
      String normalSamp, String normalDef, String output, String outputFiltered, Logger log) {

    if (VCFOps.getSamplesInFile(vcf).length != 2) {
      throw new IllegalArgumentException("This method is only designed for tumor normal renaming");
    }
    VCFFileReader reader = new VCFFileReader(new File(vcf), false);
    VariantContextWriter writer = VCFOps.initWriter(output, VCFOps.DEFUALT_WRITER_OPTIONS,
        reader.getFileHeader().getSequenceDictionary());
    VariantContextWriter writerFiltered = VCFOps.initWriter(outputFiltered,
        VCFOps.DEFUALT_WRITER_OPTIONS, reader.getFileHeader().getSequenceDictionary());

    Set<String> samps = new HashSet<String>();
    samps.add(normalSamp);
    samps.add(tumorSamp);
    ArrayList<String> attsRemoveVCAddGT = new ArrayList<String>();
    Set<VCFHeaderLine> newHeaderLines = new HashSet<VCFHeaderLine>();
    Collection<VCFInfoHeaderLine> infos = reader.getFileHeader().getInfoHeaderLines();
    for (VCFInfoHeaderLine vcfInfoHeaderLine : infos) {
      if (!vcfInfoHeaderLine.getID().equals("DB")
          && vcfInfoHeaderLine.getType() != VCFHeaderLineType.Flag) {// keep dbsnp
        if (!vcfInfoHeaderLine.getID().equals("PON")) {
          // oldHeaderLine.remove(vcfInfoHeaderLine);// remove "variant" level annotations
          attsRemoveVCAddGT.add(vcfInfoHeaderLine.getID());

          VCFFormatHeaderLine newFormat = new VCFFormatHeaderLine(vcfInfoHeaderLine.getID(),
              vcfInfoHeaderLine.isFixedCount() ? vcfInfoHeaderLine.getCount() : 1,
              vcfInfoHeaderLine.getType(), vcfInfoHeaderLine.getDescription());

          newHeaderLines.add(newFormat);// transfer to genotype level annotations for merging
        }
      }
      newHeaderLines.add(vcfInfoHeaderLine);
    }

    VCFFormatHeaderLine ADPreserve = new VCFFormatHeaderLine(GENOTYPE_INFO.AD_TUMOR.getFlag(),
        VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer,
        "Allelic depths for the ref and alt alleles in the order listed for the tumor sample");
    VCFFormatHeaderLine filterPreserve = new VCFFormatHeaderLine(
        GENOTYPE_INFO.MUTECT_FILTERS.getFlag(), VCFHeaderLineCount.UNBOUNDED,
        VCFHeaderLineType.String, "Filters applied by mutect to somatic calls");

    newHeaderLines.add(ADPreserve);
    newHeaderLines.add(filterPreserve);
    newHeaderLines.addAll(reader.getFileHeader().getFormatHeaderLines());
    ArrayList<String> attsToTransferFromNormal = new ArrayList<String>();
    for (VCFFormatHeaderLine vcfFormatHeaderLine : reader.getFileHeader().getFormatHeaderLines()) {
      VCFFormatHeaderLine normal = new VCFFormatHeaderLine(vcfFormatHeaderLine.getID() + NORMAL_TAG,
          vcfFormatHeaderLine.isFixedCount() ? vcfFormatHeaderLine.getCount() : 1,
          vcfFormatHeaderLine.getType(),
          vcfFormatHeaderLine.getDescription() + " in the normal sample from mutect calls");
      newHeaderLines.add(normal);
      attsToTransferFromNormal.add(vcfFormatHeaderLine.getID());
    }

    newHeaderLines.addAll(reader.getFileHeader().getOtherHeaderLines());
    newHeaderLines.addAll(reader.getFileHeader().getContigLines());
    newHeaderLines.addAll(reader.getFileHeader().getFilterLines());
    final VCFHeader outHeader = new VCFHeader(newHeaderLines, samps);
    writer.writeHeader(outHeader);
    writerFiltered.writeHeader(outHeader);
    int index = 0;
    int pass = 0;
    for (VariantContext vc : reader) {
      VariantContextBuilder builder = new VariantContextBuilder(vc);
      index++;
      if (index % 10000 == 0) {
        log.reportTimeInfo("Parsed " + index + " total variants, " + pass + " variants were PASS");
      }
      ArrayList<Genotype> renamed = new ArrayList<Genotype>();
      Genotype normal = rename(vc.getGenotype(normalDef), normalSamp);
      Genotype tumor = rename(vc.getGenotype(tumorDef), tumorSamp);
      tumor = transferFormat(tumor, normal, vc, attsRemoveVCAddGT, attsToTransferFromNormal);
      Hashtable<String, Object> map = new Hashtable<String, Object>();
      map.putAll(vc.getAttributes());

      renamed.add(tumor);
      renamed.add(normal);
      builder.genotypes(renamed);
      if (!renamed.get(0).sameGenotype(vc.getGenotype(tumorDef))) {
        reader.close();
        writer.close();
        writerFiltered.close();
        throw new IllegalStateException("Improprer rename");
      }
      builder.genotypes(renamed);
      if (!renamed.get(1).sameGenotype(vc.getGenotype(normalDef))) {
        reader.close();
        writer.close();
        throw new IllegalStateException("Improprer rename");
      }
      for (String key : attsRemoveVCAddGT) {
        map.remove(key);
      }
      builder.attributes(map);
      VariantContext vcRename = builder.make(true);
      writer.add(vcRename);
      if (!vcRename.isFiltered() || (vcRename.getFilters().size() == 1
          && vcRename.getFilters().contains("str_contraction"))) {
        writerFiltered.add(vcRename);
        pass++;
      }
    }
    log.reportTimeInfo("Re-named and indexed " + vcf + " to " + output);
    reader.close();
    writer.close();
    writerFiltered.close();

  }

  /**
   * GATK appends .variant## to each sample when merging individual tumor normal calls, this will
   * re-name the samples to the original (removes .variant.*)
   */
  public static void renameMergeVCF(String inputVCF, String outputVCF) {
    VCFFileReader reader = new VCFFileReader(new File(inputVCF), true);
    Set<String> samps = new HashSet<String>();
    String[] sampIn = VCFOps.getSamplesInFile(inputVCF);
    for (String element : sampIn) {
      String fix = element.replaceAll(".variant.*", "");
      samps.add(fix);
    }

    final VCFHeader outHeader =
        new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder(), samps);

    VariantContextWriter writer = VCFOps.initWriter(outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
        reader.getFileHeader().getSequenceDictionary());
    writer.writeHeader(outHeader);
    for (VariantContext vc : reader) {
      VariantContextBuilder builder = new VariantContextBuilder(vc);
      ArrayList<Genotype> renameGeno = new ArrayList<Genotype>();
      for (Genotype g : vc.getGenotypes()) {
        GenotypeBuilder gBuilder = new GenotypeBuilder(g);
        gBuilder.name(g.getSampleName().replaceAll(".variant.*", ""));
        renameGeno.add(gBuilder.make());
      }
      builder.genotypes(renameGeno);
      VariantContext vcFilt = builder.make();
      if (!vcFilt.isMonomorphicInSamples()) {
        writer.add(vcFilt);
      }
    }
    reader.close();
    writer.close();
  }

  public static void renameTumorNormalVCF(String vcf, String tumorSamp, String normalSamp,
      String output, String outputFiltered, Logger log) {
    renameAndTransferInfo(vcf, tumorSamp, "TUMOR", normalSamp, "NORMAL", output, outputFiltered,
        log);
  }

  private static Genotype transferFormat(Genotype gTumor, Genotype gNormal, VariantContext vc,
      List<String> attsToadd, List<String> attsToTransferFromNormal) {
    GenotypeBuilder builder = new GenotypeBuilder(gTumor);
    for (String att : attsToadd) {
      if (vc.hasAttribute(att)) {
        builder.attribute(att, vc.getAttribute(att));
      }
    }
    if (gTumor.hasAD()) {
      builder.attribute(GENOTYPE_INFO.AD_TUMOR.getFlag(), gTumor.getAnyAttribute("AD"));
    }

    for (String att : attsToTransferFromNormal) {
      if (gNormal.hasAnyAttribute(att)) {
        builder.attribute(att + NORMAL_TAG, gNormal.getAnyAttribute(att));
      }
    }
    if (vc.isFiltered()) {
      builder.attribute(GENOTYPE_INFO.MUTECT_FILTERS.getFlag(), vc.getFilters().toString());
    } else {
      HashSet<String> noFilt = new HashSet<String>();
      noFilt.add("PASS");
      builder.attribute(GENOTYPE_INFO.MUTECT_FILTERS.getFlag(), noFilt.toString());
    }

    return builder.make();
  }

  // double mapQ = 0;
  // double ssc = 0;
  // try {
  // if (g.hasAnyAttribute("MQ")) {
  // mapQ = Double.parseDouble(g.getAnyAttribute("MQ").toString());
  // }
  // if (g.hasAnyAttribute("SSC")) {
  // ssc = Double.parseDouble(g.getAnyAttribute("SSC").toString());
  // }
  // } catch (NumberFormatException nfe) {
  //
  // }
  // if (mapQ < 40 || ssc < 40) {
  // gBuilder.alleles(GenotypeOps.getNoCall());
  // }

}
