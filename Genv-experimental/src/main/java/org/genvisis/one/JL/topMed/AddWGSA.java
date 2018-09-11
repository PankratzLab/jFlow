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

  private static final String WGSA_SNP_FILE = "wgsaSnpFile";
  private static final String WGSA_INDEL_FILE = "wgsaIndelFile";

  private static final String CHR = "chr";
  private static final String POS = "pos";
  private static final String REF = "ref";
  private static final String ALT = "alt";
  private static final String DB_SNP = "rs_dbSNP150";

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

          addHeader(unique, reader.getFileHeader(), writer);

          for (VariantContext vc : reader) {

            VariantContextBuilder vcBuilder = new VariantContextBuilder(vc);
            Map<String, String> toAdd = new HashMap<>();

            if (vc.isIndel()) {
              processLine(lookupIndel, readerAnnotationIndel, interestIndel, vc, toAdd);

            } else {
              processLine(lookupSNP, readerAnnotationSNP, interestSNP, vc, toAdd);
            }
            for (Entry<String, String> entry : toAdd.entrySet()) {
              vcBuilder.attribute(entry.getKey(), entry.getValue());
            }
            if (".".equals(vc.getID())) {
              vcBuilder.id(toAdd.get(DB_SNP));
            }
            writer.add(vcBuilder.make());
          }
          writer.close();
        }
      }
    } catch (IOException e) {
      log.reportException(e);

    }
  }

  private static Map<String, Integer> processHeader(BufferedReader readerAnnotationIndel) throws IOException {
    Map<String, Integer> lookup = new HashMap<>();
    String[] headerIndel = readerAnnotationIndel.readLine().trim().split("\t");
    for (int i = 0; i < headerIndel.length; i++) {
      lookup.put(headerIndel[i], i);
    }
    return lookup;
  }

  private static void processLine(Map<String, Integer> lookup, BufferedReader reader,
                                  Set<String> interest, VariantContext vc,
                                  Map<String, String> toAdd) throws IOException {
    String[] line = reader.readLine().trim().split("\t");
    validate(lookup, vc, line);

    for (String key : interest) {
      toAdd.put(key, line[lookup.get(key)]);
    }
  }

  private static void validate(Map<String, Integer> lookupIndel, VariantContext vc, String[] line) {
    if (!vc.getContig().equals(line[lookupIndel.get(CHR)])) {
      throw new IllegalArgumentException("invalid contig , " + vc.getContig() + "\t"
                                         + line[lookupIndel.get(CHR)] + "\t"
                                         + vc.toStringWithoutGenotypes());
    }
    if (vc.getStart() != Integer.parseInt(line[lookupIndel.get(POS)])) {
      throw new IllegalArgumentException("invalid position , " + vc.getStart() + "\t"
                                         + line[lookupIndel.get(POS)] + "\t"
                                         + vc.toStringWithoutGenotypes());
    }
    if (!vc.getReference().getBaseString().equals(line[lookupIndel.get(REF)])) {
      throw new IllegalArgumentException("invalid ref , " + vc.getReference() + "\t"
                                         + line[lookupIndel.get(REF)] + "\t"
                                         + vc.toStringWithoutGenotypes());
    }
    if (!vc.getAlternateAllele(0).getBaseString().equals(line[lookupIndel.get(ALT)])) {
      throw new IllegalArgumentException("invalid alt , " + vc.getReference() + "\t"
                                         + line[lookupIndel.get(ALT)] + "\t"
                                         + vc.toStringWithoutGenotypes());
    }
  }

  private static void addHeader(Set<String> newHeaderLines, VCFHeader header,
                                VariantContextWriter writer) {

    List<VCFInfoHeaderLine> vcfInfos = new ArrayList<>();

    for (String key : newHeaderLines) {
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

    run(c.get(CLI.ARG_VCF), c.get(CLI.ARG_OUTDIR), c.get(WGSA_SNP_FILE), c.get(WGSA_INDEL_FILE));
  }

}
