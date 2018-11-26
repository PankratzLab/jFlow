package org.genvisis.seq.manage.mosdepth;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import org.apache.commons.collections4.BidiMap;
import org.apache.commons.math3.util.Pair;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.Maths.COMPARISON;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public abstract class NGSBinSNPSelector {

  public static final double MAPPING_QUALITY_FILTER = .9;

  String outputFile;

  String filterDefsFile;
  String mappingQualityFile;
  String afTag;
  String mqTag;
  Set<Integer> chrs;
  LocusSet<Segment> bins;
  Logger log = new Logger();
  MappabilityCompute mc;
  List<FilterSet> filters;

  abstract static class Filter {

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
    String[] filterDefs = HashVec.loadFileToStringArray(filterFile, false, null, false);
    return parseFilters(filterDefs, mc);
  }

  private List<FilterSet> parseFilters(String[] filterDefs, MappabilityCompute mc) {
    List<FilterSet> listFS = new ArrayList<>();
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

    String[] filters = {"SVM|-;AF|>.45,<.55;MQ|>=.99;AVGDP|>20",
                        "SVM|-;AF|>.44,<.56;MQ|>=.99;AVGDP|>20",
                        "SVM|-;AF|>.43,<.57;MQ|>=.99;AVGDP|>20",
                        "SVM|-;AF|>.42,<.58;MQ|>=.99;AVGDP|>20",
                        "SVM|-;AF|>.41,<.59;MQ|>=.99;AVGDP|>20",
                        "SVM|-;AF|>.4,<.6;MQ|>=.99;AVGDP|>20",
                        "SVM|-;AF|>.4,<.6;MQ|>=.98;AVGDP|>20",
                        "SVM|-;AF|>.4,<.6;MQ|>=.97;AVGDP|>20",
                        "SVM|-;AF|>.4,<.6;MQ|>=.96;AVGDP|>20",
                        "SVM|-;AF|>.4,<.6;MQ|>=.95;AVGDP|>20",
                        "SVM|-;AF|>.39,<.61;MQ|>=.95;AVGDP|>20",
                        "SVM|-;AF|>.38,<.62;MQ|>=.95;AVGDP|>20",
                        "SVM|-;AF|>.37,<.63;MQ|>=.95;AVGDP|>20",
                        "SVM|-;AF|>.36,<.64;MQ|>=.95;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.95;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.94;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.93;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.92;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.91;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.9;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.89;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.88;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.87;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.86;AVGDP|>20",
                        "SVM|-;AF|>.35,<.65;MQ|>=.85;AVGDP|>20",
                        "SVM|-;AF|>.34,<.66;MQ|>=.85;AVGDP|>20",
                        "SVM|-;AF|>.33,<.67;MQ|>=.85;AVGDP|>20",
                        "SVM|-;AF|>.32,<.68;MQ|>=.85;AVGDP|>20",
                        "SVM|-;AF|>.31,<.69;MQ|>=.85;AVGDP|>20",
                        "SVM|-;AF|>.3,<.7;MQ|>=.85;AVGDP|>20",
                        "SVM|-;AF|>.3,<.7;MQ|>=.84;AVGDP|>20",
                        "SVM|-;AF|>.3,<.7;MQ|>=.83;AVGDP|>20",
                        "SVM|-;AF|>.3,<.7;MQ|>=.82;AVGDP|>20",
                        "SVM|-;AF|>.3,<.7;MQ|>=.81;AVGDP|>20", "SVM|-;AF|>.3,<.7;MQ|>=.8;AVGDP|>20",
                        "SVM|-;AF|>.28,<.72;MQ|>=.8;AVGDP|>20",
                        "SVM|-;AF|>.26,<.74;MQ|>=.8;AVGDP|>20",
                        "SVM|-;AF|>.24,<.76;MQ|>=.8;AVGDP|>20",
                        "SVM|-;AF|>.22,<.78;MQ|>=.8;AVGDP|>20",
                        "SVM|-;AF|>.2,<.8;MQ|>=.8;AVGDP|>20", "SVM|-;AF|>.2,<.8;MQ|>=.8;AVGDP|>19",
                        "SVM|-;AF|>.2,<.8;MQ|>=.8;AVGDP|>18", "SVM|-;AF|>.2,<.8;MQ|>=.8;AVGDP|>17",
                        "SVM|-;AF|>.2,<.8;MQ|>=.8;AVGDP|>16", "SVM|-;AF|>.2,<.8;MQ|>=.8;AVGDP|>15",
                        "SVM|-;AF|>.05,<.95;MQ|>=.8;AVGDP|>10",};
    return parseFilters(filters, mc);
  }

  static class FilterSet {

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

  void setup() {
    mc = null;
    if (mappingQualityFile != null) {
      mc = MappabilityCompute.open(mappingQualityFile);
      log.reportTime("Mapping Quality Data Loaded");
    }
    if (filterDefsFile != null) {
      filters = parseFilters(filterDefsFile, mc);
    } else {
      filters = parseFilters(mc);
    }
  }

  VariantContextWriter openWriter(String exampleVCF) {
    File exFl = new File(exampleVCF);
    VariantContextWriter vcfWriter = new htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder().setOutputFile(outputFile)
                                                                                                           .setReferenceDictionary(VCFFileReader.getSequenceDictionary(exFl))
                                                                                                           .setOutputFileType(OutputType.VCF)
                                                                                                           .build();
    VCFFileReader temp = new VCFFileReader(exFl);
    VCFHeader vcfHeader = temp.getFileHeader();
    vcfHeader.addMetaDataLine(new VCFInfoHeaderLine("BINSTART", 1, VCFHeaderLineType.Integer,
                                                    "NGS Bin Start"));
    vcfHeader.addMetaDataLine(new VCFInfoHeaderLine("BINSTOP", 1, VCFHeaderLineType.Integer,
                                                    "NGS Bin Stop"));
    vcfHeader.addMetaDataLine(new VCFInfoHeaderLine("LEVEL", 1, VCFHeaderLineType.Integer,
                                                    "Level of filter used"));
    vcfHeader.addMetaDataLine(new VCFInfoHeaderLine("PASSED", 1, VCFHeaderLineType.Integer,
                                                    "Number of snps at the filter level specified"));
    vcfHeader.addMetaDataLine(new VCFInfoHeaderLine("COUNT", 1, VCFHeaderLineType.Integer,
                                                    "Total number of snps in bin"));
    vcfWriter.writeHeader(vcfHeader);
    temp.close();
    return vcfWriter;
  }

  void writeSelectedForBin(Segment bin, boolean useChrPrepend, BidiMap<String, Integer> contigMap,
                           VCFFileReader reader, VariantContextWriter vcfWriter) {
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
        if (!vc.isBiallelic()) {
          // skip multi-allelic snps
          continue;
        }
        if (vc.isFiltered()) {
          // skip VQSRTranche snps
          continue;
        }
        //        if (vc.isIndel()) {
        //          // skip indels
        //          continue;
        //        }
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
      break;
    }

    if (selected != null) {
      String id = selected.getID();
      if (ext.isMissingValue(id)) {
        id = selected.getContig() + ":" + selected.getStart();
      }

      VariantContext vc = new VariantContextBuilder(selected).alleles(selected.getAlleles())
                                                             .chr(selected.getContig()).id(id)
                                                             .start(selected.getStart())
                                                             .stop(selected.getEnd()).noGenotypes()
                                                             .attribute("BINSTART", bin.getStart())
                                                             .attribute("BINSTOP", bin.getStop())
                                                             .attribute("LEVEL", filterLevel)
                                                             .attribute("PASSED", passedCount)
                                                             .attribute("COUNT", count).make();

      vcfWriter.add(vc);
    }
  }

  public abstract void run();

}
