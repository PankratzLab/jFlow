package org.genvisis.one.JL.topMed;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.StringJoiner;

import org.genvisis.seq.manage.VCFOps;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Adds WGSA annotations to a vcf <br>
 * Note: We enforce that every variant in the .vcf has an annotation<br>
 * Note: We allow annotations that are not present in the .vcf (which are written at the end)
 */
public class AddWGSA {

  // Important notice:
  // 1. To facilitate gene-based genotype-phenotype association analysis, the most "deleterious"
  // consequence of the variant for its corresponding gene or genes were identified according to
  // http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences. Each vairant may
  // contain multiple rows in freeze.6.chr1-22X.snp.general.gz
  // and freeze.6.chr1-22X.indel.general.gz, one row for a gene. To obtain only one row for a
  // variant, I recommend to use the row with "Y" in the unique_variant column,
  // which corresponds to the most "deleterious" consequence of the variant across genes it impacts.
  // 2. The coordinates of hg38 were converted to those of hg19 via CrossMap
  // (https://doi.org/10.1093/bioinformatics/btt730). There are a few cases in which the reference
  // allele of the hg19 coordinate does not agree with
  // that of the hg38 coordinate. Those variants will have a "N" in column ref_hg19_equals_ref_hg38
  // and a "." in column alt_hg19.
  // 3. When an annotation resource has both hg38 and hg19 data, the hg38 will be used in priority.
  // When an annotation resource only has hg19 data, the annotation will be based on the CrossMap
  // converted hg19 coordiate
  // (and ref_hg19 and alt_hg19 if applicable).
  // 4. Because there are hundreds of annotation columns for each variant, I recommend read "List of
  // resources v0.75.docx" first before resorting the column description files (description.txt)
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

  private enum VARIANT_TYPE {
    INDEL, SNP;
  }

  private static LinkedHashSet<String> getHeaderOfInterest(Set<String> header) {
    LinkedHashSet<String> interest = new LinkedHashSet<>();
    for (String h : header) {
      if (!VARIANT_KEYS.contains(h) && (h.equals(DB_SNP) || h.equals(TOPMED_FRZ6_AC))) {
        interest.add(h);
      }
    }
    return interest;

  }

  private static class WGASMap {

    private final Map<String, WGSALine> map;
    private final Set<String> interest;

    /**
     * @param map
     * @param lookup
     * @param interest
     */
    private WGASMap(Map<String, WGSALine> map, Set<String> interest) {
      super();
      this.map = map;
      this.interest = interest;
    }

  }

  private static WGASMap loadMap(String wgsaFile, VARIANT_TYPE type, Logger log) {
    try (BufferedReader readerAnnotation = Files.getAppropriateReader(wgsaFile)) {
      Map<String, Integer> lookup = processHeader(readerAnnotation);
      Set<String> interest = getHeaderOfInterest(lookup.keySet());
      Map<String, WGSALine> map = new HashMap<>();
      String tmp;

      int count = 0;

      while ((tmp = readerAnnotation.readLine()) != null) {
        String[] line = tmp.trim().split("\t");
        count++;
        if (count % 100000 == 0) {
          log.reportTimeInfo(count + " annotations loaded from " + wgsaFile);
        }

        if (line[lookup.get(UNIQUE_VARIANT)].equals("Y")) {
          Map<String, String> toAdd = new HashMap<>();
          for (String key : interest) {
            toAdd.put(key, line[lookup.get(key)]);
          }
          WGSALine wLine = new WGSALine(toAdd);
          map.put(getWGSAKey(lookup, line, type), wLine);
        }
      }
      return new WGASMap(map, interest);

    } catch (IOException e) {
      log.reportException(e);
    }
    throw new IllegalStateException("cannot load annotations");
  }

  private static void run(String vcf, String outDir, String wgsaSNP, String wgsaIndel) {
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "wgsa.log");
    log.reportTimeInfo("writing results to " + outDir);
    WGASMap snpMap = loadMap(wgsaSNP, VARIANT_TYPE.SNP, log);
    WGASMap indelMap = loadMap(wgsaIndel, VARIANT_TYPE.INDEL, log);

    String outputVcf = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".wgsa.vcf.gz";
    try (VCFFileReader reader = new VCFFileReader(new File(vcf), false)) {

      VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputVcf);
      if (reader.getFileHeader().getSequenceDictionary() != null) {
        builder.setReferenceDictionary(reader.getFileHeader().getSequenceDictionary());
      }

      builder.setOption(Options.USE_ASYNC_IO);
      VariantContextWriter writer = builder.build();
      LinkedHashSet<String> unique = new LinkedHashSet<>();
      unique.addAll(snpMap.interest);
      unique.addAll(indelMap.interest);

      addHeader(unique, reader.getFileHeader(), writer, log);
      int num = 0;
      int numSNPsTransferred = 0;
      int numIndelsTransferred = 0;

      for (VariantContext vc : reader) {
        num++;
        if (num % 100000 == 0) {
          log.reportTimeInfo("Processed " + num + " variants\nNum SNP annotations transferred="
                             + numSNPsTransferred + "\nNum INDEL annotations transferred="
                             + numIndelsTransferred);

        }

        VariantContextBuilder vcBuilder = new VariantContextBuilder(vc);

        String vcKey = getVCKey(vc);
        if (vc.isIndel()) {
          numIndelsTransferred++;

          if (!indelMap.map.containsKey(vcKey)) {
            log.reportTimeInfo(vc.toStringWithoutGenotypes());
            throw new IllegalArgumentException("Key not found");
          }
          transferAnnotations(vc, vcBuilder, indelMap.map.get(vcKey).toAdd);
          indelMap.map.remove(vcKey);
        } else {
          numSNPsTransferred++;

          if (!snpMap.map.containsKey(vcKey)) {
            log.reportTimeInfo(vc.toStringWithoutGenotypes());
            throw new IllegalArgumentException("Key not found");
          }
          transferAnnotations(vc, vcBuilder, snpMap.map.get(vcKey).toAdd);
          snpMap.map.remove(vcKey);
        }
        writer.add(vcBuilder.make());

      }
      writer.close();
    }
    String outputSnpLeftover = outDir + VCFOps.getAppropriateRoot(vcf, true)
                               + ".wgsa.snp.not.in.vcf.gz";
    dumpLeftovers(outputSnpLeftover, snpMap);
    String outputIndelLeftover = outDir + VCFOps.getAppropriateRoot(vcf, true)
                                 + ".wgsa.indel.not.in.vcf.gz";
    dumpLeftovers(outputIndelLeftover, indelMap);

  }

  private static void dumpLeftovers(String outFile, WGASMap map) {
    try (PrintWriter writer = Files.getAppropriateWriter(outFile)) {
      writer.println("KEY\t" + ArrayUtils.toStr(map.interest));

      for (String key : map.map.keySet()) {
        StringJoiner out = new StringJoiner("\t");
        out.add(key);
        for (String ann : map.interest) {
          out.add(map.map.get(key).toAdd.get(ann));
        }
        writer.println(out.toString());
      }
    }
  }

  private static void transferAnnotations(VariantContext vc, VariantContextBuilder vcBuilder,
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

    /**
     * @param toAdd
     * @param line
     */
    private WGSALine(Map<String, String> toAdd) {
      super();
      this.toAdd = toAdd;
    }

  }

  private static String getVCKey(VariantContext vc) {
    String base = vc.getContig() + "_" + vc.getStart() + "_" + vc.getReference().getBaseString()
                  + "_" + vc.getAlternateAllele(0).getBaseString();

    if (vc.isIndel()) {
      return VARIANT_TYPE.INDEL.toString() + "_" + base;
    } else {
      return VARIANT_TYPE.SNP.toString() + "_" + base;
    }
  }

  private static String getWGSAKey(Map<String, Integer> lookup, String[] line, VARIANT_TYPE type) {
    String base = "chr" + line[lookup.get(CHR)] + "_" + line[lookup.get(POS)] + "_"
                  + line[lookup.get(REF)] + "_" + line[lookup.get(ALT)];

    return type.toString() + "_" + base;

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
    writer.writeHeader(header);
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
