package org.genvisis.one.JL.topMed;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import org.genvisis.CLI;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.seq.manage.VCFOps;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class AddWGSA {

  //  Important notice:
  //    1. To facilitate gene-based genotype-phenotype association analysis, the most "deleterious" consequence of the variant for its corresponding gene or genes were identified according to
  //            http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences. Each vairant may contain multiple rows in freeze.6.chr1-22X.snp.general.gz
  //            and freeze.6.chr1-22X.indel.general.gz, one row for a gene. To obtain only one row for a variant, I recommend to use the row with "Y" in the unique_variant column,
  //            which corresponds to the most "deleterious" consequence of the variant across genes it impacts.
  //    2. The coordinates of hg38 were converted to those of hg19 via CrossMap (https://doi.org/10.1093/bioinformatics/btt730). There are a few cases in which the reference allele of the hg19 coordinate does not agree with
  //            that of the hg38 coordinate. Those variants will have a "N" in column ref_hg19_equals_ref_hg38 and a "." in column alt_hg19.
  //    3. When an annotation resource has both hg38 and hg19 data, the hg38 will be used in priority. When an annotation resource only has hg19 data, the annotation will be based on the CrossMap converted hg19 coordiate
  //            (and ref_hg19 and alt_hg19 if applicable).
  //    4. Because there are hundreds of annotation columns for each variant, I recommend read "List of resources v0.75.docx" first before resorting the column description files (description.txt)
  private static final String WGSA_SNP_FILE = "wgsaSnpFile";
  private static final String WGSA_INDEL_FILE = "wgsaIndelFile";

  private static final String CHR = "chr";
  private static final String POS = "pos";
  private static final String REF = "ref";
  private static final String ALT = "alt";
  private static final String DB_SNP = "rs_dbSNP150";
  private static final String TOPMED_FRZ6_AC = "TOPMed_frz6_AC";
  private static final String UNIQUE_VARIANT = "unique_variant";

  private static final List<String> VARIANT_KEYS = Arrays.asList(CHR, POS, REF, ALT);

  private static LinkedHashSet<String> getHeaderOfInterest(Set<String> header) {
    LinkedHashSet<String> interest = new LinkedHashSet<>();
    for (String h : header) {
      if (!VARIANT_KEYS.contains(h)) {
        interest.add(h);
      }
    }
    return interest;

  }

  private static void run(String vcf, String outDir, String wgsaSNP, String wgsaIndel) {
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "wgsa.log");
    log.reportTimeInfo("writing results to " + outDir);

    try (BufferedReader readerAnnotationSNP = Files.getAppropriateReader(wgsaSNP)) {
      try (BufferedReader readerAnnotationIndel = Files.getAppropriateReader(wgsaIndel)) {

        Map<String, Integer> lookupSNP = processHeader(readerAnnotationSNP);

        Set<String> interestSNP = getHeaderOfInterest(lookupSNP.keySet());

        Map<String, Integer> lookupIndel = processHeader(readerAnnotationIndel);

        Set<String> interestIndel = getHeaderOfInterest(lookupIndel.keySet());

        String outputVcf = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".wgsa.vcf.gz";
        try (VCFFileReader reader = new VCFFileReader(new File(vcf), false)) {

          VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputVcf);
          if (reader.getFileHeader().getSequenceDictionary() != null) {
            builder.setReferenceDictionary(reader.getFileHeader().getSequenceDictionary());
          }

          VariantContextWriter writer = builder.build();
          LinkedHashSet<String> unique = new LinkedHashSet<>();
          unique.addAll(interestSNP);
          unique.addAll(interestIndel);

          addHeader(unique, reader.getFileHeader(), writer, log);
          int num = 0;
          int numSNPsTransferred = 0;
          int numIndelsTransferred = 0;
          WGSALine wSNPLine = processLine(lookupSNP, readerAnnotationSNP, interestSNP, log);
          WGSALine wIndelLine = processLine(lookupIndel, readerAnnotationIndel, interestIndel, log);
          for (VariantContext vc : reader) {
            num++;
            if (num % 1 == 0) {
              log.reportTimeInfo("Processed " + num + " variants\nNum SNP annotations transferred="
                                 + numSNPsTransferred + "\nNum INDEL annotations transferred="
                                 + numIndelsTransferred);

              //              if (vc.getStart() == 5031125) {
              //
              log.reportTimeInfo(vc.toStringWithoutGenotypes() + "\t" + wIndelLine.line[1] + "\t"
                                 + wIndelLine.line[2] + "\t" + wIndelLine.line[3] + "\t"
                                 + vc.getReference().getBaseString() + "\t"
                                 + vc.getAlternateAllele(0).getBaseString() + "\t"
                                 + vc.getReference().getBaseString().equals(wIndelLine.line[2])
                                 + vc.getAlternateAllele(0).getBaseString()
                                     .equals(wIndelLine.line[3]));
              log.reportTimeInfo(vc.toStringWithoutGenotypes() + "\t" + wSNPLine.line[1]);
            }
            //          }
            //            WANT
            //            CTTTTTT
            //            CTTTTTTT
            //            CTTTTTTT
            //            
            //            
            //            
            //            CTTTTTTT
            //            CTTTTTT
            //            CTTTTTT
            //            CTTTTTTT
            //            CTTTTTTT
            VariantContextBuilder vcBuilder = new VariantContextBuilder(vc);
            if (validate(lookupIndel, vc, wIndelLine.line)) {
              numSNPsTransferred++;
              transferAnnotations(vc, vcBuilder, wIndelLine.toAdd);
              wIndelLine = processLine(lookupIndel, readerAnnotationIndel, interestIndel, log);
            } else if (validate(lookupSNP, vc, wSNPLine.line)) {
              transferAnnotations(vc, vcBuilder, wSNPLine.toAdd);
              numIndelsTransferred++;
              wSNPLine = processLine(lookupSNP, readerAnnotationSNP, interestSNP, log);
            }
            writer.add(vcBuilder.make());
          }
          writer.close();
        }
      }
    } catch (

    IOException e) {
      log.reportException(e);

    }
  }

  static void transferAnnotations(VariantContext vc, VariantContextBuilder vcBuilder,
                                  Map<String, String> toAdd) {
    for (Entry<String, String> entry : toAdd.entrySet()) {
      vcBuilder.attribute(entry.getKey(), entry.getValue());
    }
    if (".".equals(vc.getID())) {
      vcBuilder.id(toAdd.get(DB_SNP));
    }
  }

  private static Map<String, Integer> processHeader(BufferedReader readerAnnotation) throws IOException {
    Map<String, Integer> lookup = new HashMap<>();
    String[] header = readerAnnotation.readLine().trim().split("\t");
    for (int i = 0; i < header.length; i++) {
      header[i] = header[i].replaceAll("=", ".eq.");
      lookup.put(header[i], i);
    }
    return lookup;
  }

  private static class WGSALine {

    private final Map<String, String> toAdd;
    private final String[] line;

    /**
     * @param toAdd
     * @param line
     */
    private WGSALine(Map<String, String> toAdd, String[] line) {
      super();
      this.toAdd = toAdd;
      this.line = line;
    }

  }

  private static WGSALine processLine(Map<String, Integer> lookup, BufferedReader reader,
                                      Set<String> interest, Logger log) throws IOException {
    String[] line;

    while ((line = reader.readLine().trim().split("\t")) != null) {
      if (!line[lookup.get(TOPMED_FRZ6_AC)].equals(".")
          && line[lookup.get(UNIQUE_VARIANT)].equals("Y")) break;
    }

    if (line == null) {
      log.reportTimeWarning("Returning null WGSALine");
      return null;
    }

    Map<String, String> toAdd = new HashMap<>();
    for (String key : interest) {
      toAdd.put(key, line[lookup.get(key)]);
    }
    return new WGSALine(toAdd, line);
  }

  private static boolean validate(Map<String, Integer> lookupIndel, VariantContext vc,
                                  String[] line) {
    if (!vc.getContig().equals("chr" + line[lookupIndel.get(CHR)])) {
      throw new IllegalArgumentException("invalid contig , " + vc.getContig() + "\t"
                                         + line[lookupIndel.get(CHR)] + "\t"
                                         + vc.toStringWithoutGenotypes());
    }
    if (vc.getStart() != Integer.parseInt(line[lookupIndel.get(POS)])) {
      //      throw new IllegalArgumentException("invalid position , " + vc.getStart() + "\t"
      //                                         + line[lookupIndel.get(POS)] + "\t"
      //                                         + vc.toStringWithoutGenotypes());
      return false;
    }
    if (!vc.getReference().getBaseString().equals(line[lookupIndel.get(REF)])) {
      //      throw new IllegalArgumentException("invalid ref , " + vc.getReference() + "\t"
      //                                         + line[lookupIndel.get(REF)] + "\t"
      //                                         + vc.toStringWithoutGenotypes());
      return false;
    }
    if (!vc.getAlternateAllele(0).getBaseString().equals(line[lookupIndel.get(ALT)])) {
      //      throw new IllegalArgumentException("invalid alt , " + vc.getReference() + "\t"
      //                                         + line[lookupIndel.get(ALT)] + "\t"
      //                                         + vc.toStringWithoutGenotypes());
      return false;
    }
    return true;
  }

  private static void addHeader(Set<String> newHeaderLines, VCFHeader header,
                                VariantContextWriter writer, Logger log) {

    List<VCFInfoHeaderLine> vcfInfos = new ArrayList<>();

    for (String key : newHeaderLines) {
      log.reportTimeInfo("adding header " + key);
      VCFInfoHeaderLine infoN = new VCFInfoHeaderLine(key, 1, VCFHeaderLineType.String,
                                                      key + " annotation transferred from WGSA");
      vcfInfos.add(infoN);
    }
    for (VCFInfoHeaderLine vcfInfo : vcfInfos) {
      header.addMetaDataLine(vcfInfo);
    }
    writer.writeHeader(new VCFHeader(header.getMetaDataInInputOrder(), new HashSet<String>()));
  }

  public static void main(String[] args) {
    CLI c = new CLI(AddWGSA.class);
    c.addArg(CLI.ARG_VCF, CLI.DESC_VCF);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR);
    c.addArg(WGSA_SNP_FILE, "corresponding wgsa snp annotation file");
    c.addArg(WGSA_INDEL_FILE, "corresponding wgsa indel annotation file");
    c.parseWithExit(args);
    run(c.get(CLI.ARG_VCF), c.get(CLI.ARG_OUTDIR), c.get(WGSA_SNP_FILE), c.get(WGSA_INDEL_FILE));
  }

}
