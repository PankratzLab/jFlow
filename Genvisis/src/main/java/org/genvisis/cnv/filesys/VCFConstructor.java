package org.genvisis.cnv.filesys;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import org.genvisis.cnv.Resources;
import org.genvisis.cnv.Resources.GENOME_BUILD;
import org.genvisis.seq.ReferenceGenome;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.CLI.Arg;
import org.pankratzlab.utils.filesys.SnpMarkerSet;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.Sort;
import org.pankratzlab.common.ext;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class VCFConstructor {

  Logger log;
  String markersFile;
  String dosageFile;
  String outputFile;

  String[] markers;
  int[][] locations;
  Map<String, List<Allele>> alleles;
  Map<String, Collection<Genotype>> genotypes;
  String[] ids;
  double windowForHet = HET_WINDOW_STRICT;

  private static final double HET_WINDOW_STRICT = 0.1;

  protected void readInputFile() {
    if (markersFile == null || "".equals(markersFile) || !Files.exists(markersFile)) {
      String error = markersFile == null ? "Marker list file hasn't been set yet. Please specify marker list file and try again."
                                         : "".equals(markersFile) ? "Marker list file is blank - was it set properly? Please fix marker list file name and try again."
                                                                  : !Files.exists(markersFile) ? "Marker list file {"
                                                                                                 + markersFile
                                                                                                 + "} doesn't exist. Please fix marker list file name and try again."
                                                                                               : "Problem loading marker list file {"
                                                                                                 + markersFile
                                                                                                 + "}";
      throw new UnsupportedOperationException(error);
    }
    markers = HashVec.loadFileToStringArray(markersFile, false, new int[] {0}, true, false, "\t");
    SnpMarkerSet markerSet = new SnpMarkerSet(markers);
    markerSet.parseSNPlocations(log);
    locations = markerSet.getChrAndPositionsAsInts();
    parseAlleles(markerSet.getMarkerNames());
  }

  protected void readDosageFile() {
    if (dosageFile == null || "".equals(dosageFile) || !Files.exists(dosageFile)) {
      String error = dosageFile == null ? "Dosage file hasn't been set yet. Please specify dosage file and try again."
                                        : "".equals(dosageFile) ? "Dosage file is blank - was it set properly? Please fix dosage file name and try again."
                                                                : !Files.exists(dosageFile) ? "Dosage file {"
                                                                                              + dosageFile
                                                                                              + "} doesn't exist. Please fix dosage file name and try again."
                                                                                            : "Problem loading dosage file {"
                                                                                              + dosageFile
                                                                                              + "}";
      throw new UnsupportedOperationException(error);
    }
    // TODO refactor to use DosageData - adds complexity and requirements (map file, id file, etc)
    // for FRZ5 data, all info is included in dosage file, so using DosageData is more than what's
    // needed
    // DosageData dd = new DosageData(dosageFile, null, null, false, log);
    String[][] matr = HashVec.loadFileToStringMatrix(dosageFile, false, null);
    ids = ArrayUtils.subArray(Matrix.extractColumn(matr, 1), 1);
    String[] snps = ArrayUtils.subArray(matr[0], 2);

    genotypes = new HashMap<>();
    for (int i = 0; i < snps.length; i++) {
      String mkr = snps[i];
      List<Allele> all = alleles.get(mkr);
      int m = i + 2;
      Collection<Genotype> genos = new ArrayList<>();
      for (int id = 0; id < ids.length; id++) {
        int ind = id + 1;
        String indDose = matr[ind][m];
        Genotype g;
        if (ext.isMissingValue(indDose)) {
          g = createGeno(GenotypeType.NO_CALL, id, all);
        } else {
          try {
            double dose = Double.parseDouble(indDose);
            if (dose >= 0 + windowForHet) {
              if (dose <= 2 - windowForHet) {
                g = createGeno(GenotypeType.HOM_VAR, id, all);
              } else {
                g = createGeno(GenotypeType.HET, id, all);
              }
            } else {
              g = createGeno(GenotypeType.HOM_REF, id, all);
            }
          } catch (NumberFormatException e) {
            log.reportError("Couldn't parse dosage {" + indDose + "} for id {" + ids[id]
                            + "}. Setting to missing/no-call.");
            g = createGeno(GenotypeType.NO_CALL, id, all);
          }
        }
        genos.add(g);
      }
      genotypes.put(mkr, genos);
    }
  }

  private Genotype createGeno(GenotypeType call, int id, List<Allele> all) {
    switch (call) {
      case HOM_REF:
        return GenotypeBuilder.create(ids[id],
                                      ArrayUtils.toList(new Allele[] {all.get(0), all.get(0)}));
      case HET:
        return GenotypeBuilder.create(ids[id],
                                      ArrayUtils.toList(new Allele[] {all.get(0), all.get(1)}));
      case HOM_VAR:
        return GenotypeBuilder.create(ids[id],
                                      ArrayUtils.toList(new Allele[] {all.get(1), all.get(1)}));
      case NO_CALL:
      default:
        return GenotypeBuilder.create(ids[id], ArrayUtils.toList(new Allele[] {Allele.NO_CALL,
                                                                               Allele.NO_CALL}));
    }
  }

  private void parseAlleles(String[] markerNames) {
    alleles = new HashMap<>();
    for (int i = 0; i < markerNames.length; i++) {
      if (markerNames[i].startsWith("rs")) {
        log.reportError("Cannot parse alleles from RS marker");
        alleles.put(markerNames[i], new ArrayList<Allele>());
      } else {
        String[] pts = markerNames[i].split(":");
        ArrayList<Allele> a = new ArrayList<>();
        a.add(Allele.create(pts[2], true));
        a.add(Allele.create(pts[3], false));
        alleles.put(markerNames[i], a);
      }
    }
  }

  public boolean prep() {
    try {
      readInputFile();
    } catch (UnsupportedOperationException e) {
      log.reportError(e.getMessage());
      return false;
    }
    try {
      readDosageFile();
    } catch (UnsupportedOperationException e) {
      log.reportError(e.getMessage());
      return false;
    }
    return true;
  }

  public void build() {
    VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputFile);
    builder.clearOptions();
    builder.setOption(Options.INDEX_ON_THE_FLY);
    HashSet<VCFHeaderLine> lines = new HashSet<>();
    VCFFormatHeaderLine format = new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "GT");
    lines.add(format);

    List<String> idNames = new ArrayList<>();
    for (String s : ids) {
      idNames.add(s);
    }
    VCFHeader vcfHeader = new VCFHeader(lines, idNames);

    SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome(Resources.genome(GENOME_BUILD.HG19,
                                                                                       log)
                                                                               .getFASTA()
                                                                               .getAbsolute(),
                                                                      log).getIndexedFastaSequenceFile()
                                                                          .getSequenceDictionary();

    builder.setReferenceDictionary(samSequenceDictionary);
    vcfHeader.setSequenceDictionary(samSequenceDictionary);
    VariantContextWriter writer = builder.build();
    vcfHeader.hasGenotypingData();
    writer.writeHeader(vcfHeader);

    int[] sortedIndices = Sort.getSort2DIndices(Matrix.extractColumn(locations, 0),
                                                Matrix.extractColumn(locations, 1));
    for (int m = 0; m < markers.length; m++) {
      int i = sortedIndices[m];
      VariantContextBuilder builderVc = new VariantContextBuilder();
      builderVc.chr("chr" + locations[i][0]);
      int len = alleles.get(markers[i]).get(0).length() - 1;
      builderVc.alleles(alleles.get(markers[i]));
      builderVc.start(locations[i][1]);
      builderVc.stop(locations[i][1] + len);
      builderVc.id(markers[i]);
      builderVc.genotypes(genotypes.get(markers[i]));

      writer.add(builderVc.make());
    }
    writer.close();
  }

  public static void main(String[] args) {
    String markersFile = "markers.txt";
    String dosageFile = "data.db.xln.gz";
    String outputFile = "data.vcf";
    double window = HET_WINDOW_STRICT;

    Object[][] argSet = {{"markers", "Marker list file", markersFile, Arg.STRING},
                         {"data", "Dosage data file", dosageFile, Arg.STRING},
                         {CLI.ARG_OUTFILE, CLI.DESC_OUTFILE, outputFile, Arg.STRING},
                         {"window", "Window inset for heterozygous calls (0 + window, 2 - window)",
                          window, Arg.NUMBER},};

    CLI cli = new CLI(VCFConstructor.class);

    for (Object[] arg : argSet) {
      cli.addArgWithDefault((String) arg[0], (String) arg[1], arg[2] + "", (Arg) arg[3]);
    }

    cli.parseWithExit(args);

    VCFConstructor c = new VCFConstructor();
    c.markersFile = cli.get((String) argSet[0][0]);
    c.dosageFile = cli.get((String) argSet[1][0]);
    c.outputFile = cli.get((String) argSet[2][0]);
    c.windowForHet = cli.getD((String) argSet[3][0]);
    c.log = new Logger();
    c.readInputFile();
    c.readDosageFile();
    c.build();

    //// to create an index file from an existing vcf, use:
    // VCFOps.verifyIndex("F:/temp/variantviewer/output.vcf", new Logger());

  }

}
