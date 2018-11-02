package org.genvisis.seq.manage;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.apache.commons.math3.util.Pair;
import org.genvisis.CLI;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.Maths.COMPARISON;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class NGSBinSNPSelector {

  public static final double MAPPING_QUALITY_FILTER = .9;
  private static final boolean DEBUG = true;

  String outputDir;
  String inputDir;
  String format;
  String token;
  String mappingQualityFile;
  String afTag;
  String mqTag;
  Set<Integer> chrs;
  LocusSet<Segment> bins;
  Logger log = new Logger();
  MappabilityCompute mc;

  private abstract static class Filter {

    public boolean use = true;

    public Filter(String tag) {
      this.tag = tag;
    }

    String tag;

    public void setUse(boolean use) {
      this.use = use;
    }

    public abstract boolean passes(VariantContext vc);

  }

  private static class ComparisonFilter extends Filter {

    public ComparisonFilter(String tag) {
      super(tag);
      this.comparisons = new ArrayList<>();
    }

    public ComparisonFilter addComparison(COMPARISON comp, double value) {
      this.comparisons.add(new Pair<>(comp, value));
      return this;
    }

    List<Pair<COMPARISON, Double>> comparisons;

    @Override
    public boolean passes(VariantContext vc) {
      String v = vc.getAttributes().get(tag).toString();
      if (v.contains(",")) {
        // multi-allelic
        return false;
      }
      double val = Double.parseDouble(v);
      for (Pair<COMPARISON, Double> comp : comparisons) {
        if (!comp.getFirst().check(val, comp.getSecond())) {
          return false;
        }
      }
      return true;
    }

  }

  private static class TagFilter extends Filter {

    boolean required;

    public TagFilter(String tag, boolean required) {
      super(tag);
      this.required = required;
    }

    @Override
    public boolean passes(VariantContext vc) {
      return vc.getFilters().contains(tag) && required;
    }

  }

  private static class MappingQualityFilter extends ComparisonFilter {

    MappabilityCompute mc;

    public MappingQualityFilter(String tag, MappabilityCompute mc) {
      super(tag);
      this.mc = mc;
    }

    @Override
    public boolean passes(VariantContext vc) {
      double v = mc.getAverageMap(vc);
      for (Pair<COMPARISON, Double> comp : comparisons) {
        if (!comp.getFirst().check(v, comp.getSecond())) {
          return false;
        }
      }
      return true;
    }

  }

  private List<FilterSet> parseFilters(String filterFile, MappabilityCompute mc) {
    List<FilterSet> listFS = new ArrayList<>();
    String[] filterDefs = HashVec.loadFileToStringArray(filterFile, false, null, false);
    for (String fsd : filterDefs) {
      listFS.add(parseFilterSet(fsd));
    }
    return listFS;
  }

  private FilterSet parseFilterSet(String filterLine) {
    String[] filterParts = filterLine.split(";");
    List<Filter> filters = new ArrayList<>();
    for (String filter : filterParts) {
      Filter f = parseFilter(filter);
      filters.add(f);
    }
    return new FilterSet(filters);
  }

  private Filter parseFilter(String filter) {
    String filterDef = filter;
    // test for usage tag
    boolean use = !filterDef.endsWith("!");
    if (!use) {
      filterDef = filter.substring(0, filter.length() - 1);
    }

    String[] pts = filterDef.split("\\|");
    if (pts.length != 2) {
      throw new RuntimeException("Invalid filter definition: " + filterDef);
    }

    Filter newFilter;
    // test for tag
    if (pts[1].charAt(0) == '+' || pts[1].charAt(0) == '-') {
      newFilter = new TagFilter(pts[0], pts[1].charAt(0) == '+');
    } else {
      // check for multiple thresholds
      String[] compSets = pts[1].split(",");
      List<Pair<COMPARISON, Double>> comps = parseComparisons(filterDef, compSets);
      if (pts[0].equalsIgnoreCase(mqTag)) {
        newFilter = new MappingQualityFilter(pts[0], mc);
      } else {
        newFilter = new ComparisonFilter(pts[0]);
      }
      for (Pair<COMPARISON, Double> pr : comps) {
        ((ComparisonFilter) newFilter).addComparison(pr.getFirst(), pr.getSecond());
      }
    }
    newFilter.setUse(use);

    return newFilter;
  }

  private List<Pair<COMPARISON, Double>> parseComparisons(String fullFilter, String[] compSets) {
    List<Pair<COMPARISON, Double>> comps = new ArrayList<>();
    for (String comp : compSets) {
      // greedy test for comparison
      COMPARISON compar = COMPARISON.forSymbol(comp.substring(0, 2));
      int sz = 2;
      if (compar == null) {
        // couldn't find a two-char comp, check for single char
        compar = COMPARISON.forSymbol(comp.substring(0, 1));
        sz = 1;
      }
      if (compar == null) {
        throw new RuntimeException("Invalid filter definition: " + fullFilter);
      }
      double compVal = Double.parseDouble(comp.substring(sz));
      comps.add(new Pair<>(compar, compVal));
    }
    return comps;
  }

  private List<FilterSet> parseFilters(MappabilityCompute mc) {
    String[] filters = {"SVM|-!;AF|>.4,<.6;MQ|>=.99;AVGDP|>20",
                        "SVM|-!;AF|>.4,<.6;MQ|>=.95;AVGDP|>20",
                        "SVM|-!;AF|>.35,<.65;MQ|>=.9;AVGDP|>20",
                        "SVM|-!;AF|>.3,<.7;MQ|>=.85;AVGDP|>20",
                        "SVM|-!;AF|>.2,<.8;MQ|>=.8;AVGDP|>15",
                        "SVM|-!;AF|>.05,<.95;MQ|>=.8;AVGDP|>10",};
    List<FilterSet> listFS = new ArrayList<>();
    for (String fsd : filters) {
      FilterSet f = parseFilterSet(fsd);
      listFS.add(f);
    }
    return listFS;
  }

  private static class FilterSet {

    public FilterSet(List<Filter> filters) {
      this.filters = new ArrayList<>(filters);
    }

    List<Filter> filters;

    public List<Filter> check(VariantContext vc) {
      List<Filter> failures = new ArrayList<>();
      for (Filter f : filters) {
        if (!f.passes(vc)) {
          failures.add(f);
        }
      }
      return failures;
    }

    public VariantContext best(List<VariantContext> passing, String afTag) {
      double best = 1;
      VariantContext bestVC = null;
      for (VariantContext vc : passing) {
        double afD = Math.abs(Double.parseDouble(vc.getAttributes().get(afTag).toString()) - .5);
        if (afD < best) {
          bestVC = vc;
          best = afD;
        }
      }
      return bestVC;
    }

  }

  String[] hdr = {"BIN", "ID", "CHR", "POS", "REF", "ALT", "MAPQ", "MAF", "PASSED", "LEVEL",
                  "TOTAL"};

  public void run() {
    mc = null;
    if (mappingQualityFile != null) {
      mc = MappabilityCompute.open(mappingQualityFile);
      log.reportTime("Mapping Quality Data Loaded");
    }
    List<FilterSet> filters = parseFilters(mc);

    Map<Integer, String> allVCFs = new HashMap<>();

    for (int c : chrs) {
      String vcf = inputDir + format.replace(token, Integer.toString(c));
      if (!Files.exists(vcf)) {
        if (Files.exists(inputDir + format.replace(token, "0" + Integer.toString(c)))) {
          vcf = inputDir + format.replace(token, "0" + Integer.toString(c));
        }
      }
      if (!Files.exists(vcf)) continue;
      allVCFs.put(c, vcf);
    }

    File exFl = new File(allVCFs.values().iterator().next());
    VariantContextWriter vcfWriter = new htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder().setOutputFile(outputDir
                                                                                                                          + "selected.vcf")
                                                                                                           .setReferenceDictionary(VCFFileReader.getSequenceDictionary(exFl))
                                                                                                           .setOutputFileType(OutputType.VCF)
                                                                                                           .build();
    VCFFileReader temp = new VCFFileReader(exFl);
    VCFHeader vcfHeader = temp.getFileHeader();
    vcfHeader.addMetaDataLine(new VCFInfoHeaderLine("BINSTART", 1, VCFHeaderLineType.Integer,
                                                    "NGS Bin Start"));
    vcfHeader.addMetaDataLine(new VCFInfoHeaderLine("BINSTOP", 1, VCFHeaderLineType.Integer,
                                                    "NGS Bin Stop"));
    vcfWriter.writeHeader(vcfHeader);
    temp.close();

    for (Entry<Integer, String> vcf : allVCFs.entrySet()) {
      long startN = System.nanoTime();
      PrintWriter writerDEBUG = DEBUG ? Files.getAppropriateWriter(outputDir + "all_chr"
                                                                   + vcf.getKey() + ".xln")
                                      : null;
      PrintWriter writer = Files.getAppropriateWriter(outputDir + "selected_chr" + vcf.getKey()
                                                      + ".xln");
      writer.println(ArrayUtils.toStr(hdr, "\t"));
      if (writerDEBUG != null) {
        writerDEBUG.println(ArrayUtils.toStr(hdr, "\t"));
      }
      VCFHeader header;

      try (VCFFileReader reader = new VCFFileReader(new File(vcf.getValue()),
                                                    Files.exists(vcf + ".tbi"))) {
        header = reader.getFileHeader();

        // discover which chrs are in this vcf file
        BidiMap<String, Integer> contigMap = new DualHashBidiMap<>();
        boolean useChrPrepend = header.getContigLines().get(0).getID().startsWith("chr");
        header.getContigLines().forEach(vch -> {
          // ensure parsability
          contigMap.put(vch.getID(), (int) Positions.chromosomeNumber(vch.getID()));
        });

        for (Segment bin : bins.getLoci()) {
          // no variants in file for this bin's chromosome  
          if (bin.getChr() != vcf.getKey() || !contigMap.containsValue((int) bin.getChr())
              || !chrs.contains((int) bin.getChr())) {
            //            System.err.println("Error - bin contig " + bin.getChr()
            //                               + " was not present in the VCF!");
            continue;
          }

          List<VariantContext> iter = reader.query((useChrPrepend ? "chr" : "") + bin.getChr(),
                                                   bin.getStart(), bin.getStop())
                                            .toList();

          int count = iter.size();
          VariantContext selected = null;
          int passedCount = 0;
          int filterLevel = -1;

          for (int i = 0; i < filters.size(); i++) {
            FilterSet filterSet = filters.get(i);
            List<VariantContext> passed = new ArrayList<>();
            for (VariantContext vc : iter) {
              List<Allele> alts = vc.getAlternateAlleles();
              if (alts.size() > 1 || alts.get(0).getBaseString().length() > 1
                  || vc.getReference().getBaseString().length() > 1) {
                // skip multi-allelic snps and indels 
                continue;
              }
              List<Filter> failedFilters = filterSet.check(vc);
              if (failedFilters.size() == 0) {
                passed.add(vc);
                continue;
              }
              boolean trueFail = false;
              for (Filter f : failedFilters) {
                if (f.use) {
                  trueFail = true;
                  break;
                }
              }
              if (!trueFail) {
                passed.add(vc);
              }
            }
            passedCount = passed.size();
            filterLevel = i;
            if (passed.size() == 0) {
              continue;
            }
            if (passed.size() == 1) {
              selected = passed.get(0);
            } else {
              selected = filterSet.best(passed, afTag);
            }
            if (DEBUG) {
              for (VariantContext vc : passed) {
                String id = vc.getID();
                if (ext.isMissingValue(id)) {
                  id = vc.getContig() + ":" + vc.getStart();
                }
                StringBuilder lineOut = new StringBuilder();
                lineOut.append(bin.getUCSClocation()).append("\t");
                lineOut.append(id).append("\t");
                lineOut.append(vc.getContig()).append("\t");
                lineOut.append(vc.getStart()).append("\t");
                lineOut.append(vc.getReference().getBaseString()).append("\t");
                lineOut.append(ArrayUtils.toStr(vc.getAlternateAlleles(), ",")).append("\t");
                lineOut.append(mc == null ? "." : mc.getAverageMap(vc)).append("\t");
                lineOut.append(vc.getAttributes().get(afTag).toString()).append("\t");
                lineOut.append(passedCount).append("\t");
                lineOut.append((i + 1)).append("\t");
                lineOut.append(count);
                writerDEBUG.println(lineOut.toString());
              }
            }
            break;
          }

          if (selected != null) {
            String id = selected.getID();
            if (ext.isMissingValue(id)) {
              id = selected.getContig() + ":" + selected.getStart();
            }
            StringBuilder lineOut = new StringBuilder();
            lineOut.append(bin.getUCSClocation()).append("\t");
            lineOut.append(id).append("\t");
            lineOut.append(selected.getContig()).append("\t");
            lineOut.append(selected.getStart()).append("\t");
            lineOut.append(selected.getReference().getBaseString()).append("\t");
            lineOut.append(ArrayUtils.toStr(selected.getAlternateAlleles(), ",")).append("\t");
            lineOut.append(mc == null ? "." : mc.getAverageMap(selected)).append("\t");
            lineOut.append(selected.getAttributes().get(afTag).toString()).append("\t");
            lineOut.append(passedCount).append("\t");
            lineOut.append((filterLevel + 1)).append("\t");
            lineOut.append(count);
            writer.println(lineOut.toString());

            VariantContext vc = new VariantContextBuilder().alleles(selected.getAlleles())
                                                           .chr(selected.getContig()).id(id)
                                                           .start(selected.getStart())
                                                           .stop(selected.getEnd()).noGenotypes()
                                                           .attribute("BINSTART", bin.getStart())
                                                           .attribute("BINSTOP", bin.getStop())
                                                           .make();

            vcfWriter.add(vc);
          } else {
            writer.println(bin.getUCSClocation() + "\t.\t.\t.\t.\t.\t.\t.\t0\t.\t" + count);
          }
          writer.flush();
          if (DEBUG) {
            writerDEBUG.flush();
          }
        }

      }

      writer.close();
      if (DEBUG) {
        writerDEBUG.close();
      }
      vcfWriter.close();
      log.reportTime("Finished processing chr " + vcf.getKey() + " in "
                     + ext.getTimeElapsedNanos(startN) + "!");
    }
  }

  public static void main(String[] args) {
    CLI cli = new CLI(NGSBinSNPSelector.class);

    cli.addArg("out", "Output directory");
    cli.addArg("in", "Input directory");
    cli.addArg("format", "File format with special token where the chromosome number goes");
    cli.addArg("token", "Special filename token to replace with chromosome number", "##", false);
    cli.addArg("mapFile", "Mapping quality file", false);
    cli.addArg("mq", "Mapping quality tag in the filters definitions file.", "MQ", false);
    cli.addArg("af", "Allele frequency tag in the VCF and filters definitions file", "AF", false);
    cli.addArg("chrs", "Comma-delimited list of chromosomes to process.", false);
    cli.addArg("bed", "BED file with predefined regions to use.", false);
    cli.addArg("bin", "Bin size to use when splitting up reference genome, default value 1000.",
               false);
    cli.addGroup("bed", "bin");

    cli.parseWithExit(args);

    Logger log = new Logger();
    NGSBinSNPSelector selector = new NGSBinSNPSelector();

    selector.outputDir = cli.get("out");
    selector.inputDir = cli.get("in");
    selector.format = cli.get("format");
    selector.token = cli.get("token");
    selector.mappingQualityFile = cli.has("mapFile") ? cli.get("mapFile") : null;
    selector.afTag = cli.has("af") ? cli.get("af") : "AF";
    selector.mqTag = cli.has("mq") ? cli.get("mq") : "MQ";
    int bin = 1000;
    String bedFile = null;
    if (cli.has("bin")) {
      bin = cli.getI("bin");
    } else if (cli.has("bed")) {
      bedFile = cli.get("bed");
    } else {
      log.reportTime("No BED file specified, nor was a bin size specified, so the default bin size of "
                     + bin + " will be used.");
    }

    Set<Integer> chrs = null;
    if (cli.has("chrs")) {
      String[] v = cli.get("chrs").split(",");
      chrs = new HashSet<>();
      for (String vi : v) {
        chrs.add(Integer.parseInt(vi));
      }
    } else {
      chrs = Sets.newHashSet(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                             21, 22, 23, 24, 25);
    }
    selector.chrs = chrs;
    if (bedFile != null) {
      selector.bins = new BEDFileReader(bedFile, false).loadAll(log).getStrictSegmentSet();
    } else {
      selector.bins = new ReferenceGenome(GENOME_BUILD.HG19, new Logger()).getBins(bin, chrs);
    }
    selector.run();
  }

}
