package org.genvisis.seq.manage;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.apache.commons.math3.util.Pair;
import org.genvisis.CLI;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.Maths.COMPARISON;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class NGSBinSNPSelector {

  public static final double MAPPING_QUALITY_FILTER = .9;
  private static final boolean DEBUG = true;

  String outputDir;
  String inputDir;
  String format;
  String token;
  String mappingQualityFile;
  String afTag;
  Set<Integer> chrs;
  LocusSet<Segment> bins;
  Logger log = new Logger();

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

    public ComparisonFilter(String tag, COMPARISON comp, double compValue) {
      super(tag);
      this.comparisons = new ArrayList<>();
      this.comparisons.add(new Pair<COMPARISON, Double>(comp, compValue));
    }

    public ComparisonFilter(String tag, COMPARISON comp1, double compValue1, COMPARISON comp2,
                            double compValue2) {
      super(tag);
      this.comparisons = new ArrayList<>();
      this.comparisons.add(new Pair<COMPARISON, Double>(comp1, compValue1));
      this.comparisons.add(new Pair<COMPARISON, Double>(comp2, compValue2));
    }

    List<Pair<COMPARISON, Double>> comparisons;

    @Override
    public boolean passes(VariantContext vc) {
      String af = vc.getAttributes().get(tag).toString();
      if (af.contains(",")) {
        // multi-allelic
        return false;
      }
      double val = Double.parseDouble(af);
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

  private static class MappingQualityFilter extends Filter {

    MappabilityCompute mc;
    double threshold;
    COMPARISON comp;

    public MappingQualityFilter(String tag, MappabilityCompute mc, COMPARISON comp, double thresh) {
      super(tag);
      this.mc = mc;
      this.comp = comp;
      this.threshold = thresh;
    }

    @Override
    public boolean passes(VariantContext vc) {
      return comp.check(mc.getAverageMap(vc), threshold);
    }

  }

  private List<FilterSet> parseFilters(MappabilityCompute mc) {
    List<FilterSet> listFS = new ArrayList<>();
    {
      Filter svmFilter = new TagFilter("SVM", false);
      svmFilter.setUse(false);
      Filter mafFilter = new ComparisonFilter("AF", COMPARISON.GT, .4, COMPARISON.LT, .6);
      Filter depFilter = new ComparisonFilter("AVGDP", COMPARISON.GT, 20);
      if (mc != null) {
        Filter mqFilter = new MappingQualityFilter("MQ", mc, COMPARISON.GTE, .95);
        mqFilter.setUse(false);
        listFS.add(new FilterSet(svmFilter, mafFilter, mqFilter, depFilter));
      } else {
        listFS.add(new FilterSet(svmFilter, mafFilter, depFilter));
      }
    }
    {
      Filter svmFilter = new TagFilter("SVM", false);
      svmFilter.setUse(false);
      Filter mafFilter = new ComparisonFilter("AF", COMPARISON.GT, .35, COMPARISON.LT, .65);
      Filter depFilter = new ComparisonFilter("AVGDP", COMPARISON.GT, 20);
      if (mc != null) {
        Filter mqFilter = new MappingQualityFilter("MQ", mc, COMPARISON.GTE, .9);
        mqFilter.setUse(false);
        listFS.add(new FilterSet(svmFilter, mafFilter, mqFilter, depFilter));
      } else {
        listFS.add(new FilterSet(svmFilter, mafFilter, depFilter));
      }
    }
    {
      Filter svmFilter = new TagFilter("SVM", false);
      svmFilter.setUse(false);
      Filter mafFilter = new ComparisonFilter("AF", COMPARISON.GT, .3, COMPARISON.LT, .7);
      Filter depFilter = new ComparisonFilter("AVGDP", COMPARISON.GT, 20);
      if (mc != null) {
        Filter mqFilter = new MappingQualityFilter("MQ", mc, COMPARISON.GTE, .85);
        mqFilter.setUse(false);
        listFS.add(new FilterSet(svmFilter, mafFilter, mqFilter, depFilter));
      } else {
        listFS.add(new FilterSet(svmFilter, mafFilter, depFilter));
      }
    }
    {
      Filter svmFilter = new TagFilter("SVM", false);
      svmFilter.setUse(false);
      Filter mafFilter = new ComparisonFilter("AF", COMPARISON.GT, .2, COMPARISON.LT, .8);
      Filter depFilter = new ComparisonFilter("AVGDP", COMPARISON.GT, 15);
      if (mc != null) {
        Filter mqFilter = new MappingQualityFilter("MQ", mc, COMPARISON.GTE, .8);
        mqFilter.setUse(false);
        listFS.add(new FilterSet(svmFilter, mafFilter, mqFilter, depFilter));
      } else {
        listFS.add(new FilterSet(svmFilter, mafFilter, depFilter));
      }
    }
    {
      Filter svmFilter = new TagFilter("SVM", false);
      svmFilter.setUse(false);
      Filter mafFilter = new ComparisonFilter("AF", COMPARISON.GT, .05, COMPARISON.LT, .95);
      Filter depFilter = new ComparisonFilter("AVGDP", COMPARISON.GT, 10);
      if (mc != null) {
        Filter mqFilter = new MappingQualityFilter("MQ", mc, COMPARISON.GTE, .8);
        mqFilter.setUse(false);
        listFS.add(new FilterSet(svmFilter, mafFilter, mqFilter, depFilter));
      } else {
        listFS.add(new FilterSet(svmFilter, mafFilter, depFilter));
      }
    }

    return listFS;
  }

  private static class FilterSet {

    public FilterSet(Filter... filters) {
      this.filters = new ArrayList<>();
      for (Filter f : filters) {
        this.filters.add(f);
      }
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

  public void run() {
    MappabilityCompute mc = MappabilityCompute.open(mappingQualityFile);
    log.reportTime("Mapping Quality Data Loaded");
    List<FilterSet> filters = parseFilters(mc);

    for (int c : chrs) {
      long startN = System.nanoTime();
      PrintWriter writerDEBUG = DEBUG ? Files.getAppropriateWriter(outputDir + "all_chr" + c
                                                                   + ".xln")
                                      : null;
      PrintWriter writer = Files.getAppropriateWriter(outputDir + "selected_chr" + c + ".xln");
      writer.println("BIN\tID\tMAP_Q\tMAF\tPASSED\tLEVEL\tTOTAL");
      if (writerDEBUG != null) {
        writerDEBUG.println("BIN\tID\tMAP_Q\tMAF\tPASSED\tLEVEL\tTOTAL");
      }

      VCFHeader header;
      String vcf = format.replace(token, Integer.toString(c));
      if (!Files.exists(vcf)) {
        if (Files.exists(format.replace(token, "0" + Integer.toString(c)))) {
          vcf = format.replace(token, "0" + Integer.toString(c));
        }
      }
      try (VCFFileReader reader = new VCFFileReader(new File(vcf), Files.exists(vcf + ".tbi"))) {
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
          if (!contigMap.containsValue((int) bin.getChr())) {
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
            if (passed.size() == 0) {
              continue;
            }
            if (passed.size() == 1) {
              selected = passed.get(0);
            } else {
              selected = filterSet.best(passed, afTag);
            }
            filterLevel = i;
            if (DEBUG) {
              for (VariantContext vc : passed) {
                String id = vc.getID();
                if (ext.isMissingValue(id)) {
                  id = vc.getContig() + ":" + vc.getStart();
                }
                writerDEBUG.println(bin.getUCSClocation() + "\t" + id + "\t"
                                    + (mc == null ? "." : mc.getAverageMap(vc)) + "\t"
                                    + vc.getAttributes().get(afTag).toString() + "\t" + passedCount
                                    + "\t" + (i + 1) + "\t" + count);
              }
            }
            break;
          }

          if (selected != null) {
            String id = selected.getID();
            if (ext.isMissingValue(id)) {
              id = selected.getContig() + ":" + selected.getStart();
            }

            writer.println(bin.getUCSClocation() + "\t" + id + "\t"
                           + (mc == null ? "." : mc.getAverageMap(selected)) + "\t"
                           + selected.getAttributes().get(afTag).toString() + "\t" + passedCount
                           + "\t" + (filterLevel + 1) + "\t" + count);
          } else {
            writer.println(bin.getUCSClocation() + "\t.\t.\t.\t0\t.\t" + count);
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
      log.reportTime("Finished processing chr " + c + " in " + ext.getTimeElapsedNanos(startN)
                     + "!");
    }
  }

  public static void main(String[] args) {
    CLI cli = new CLI(NGSBinSNPSelector.class);

    cli.addArg("out", "Output directory");
    cli.addArg("in", "Input directory");
    cli.addArg("format", "File format with special token where the chromosome number goes");
    cli.addArg("token", "Special filename token to replace with chromosome number", "##", false);
    cli.addArg("mapQ", "Mapping quality file", false);
    cli.addArg("af", "Allele frequency tag in the VCF and filters definitions file", "AF", false);
    cli.addArg("chrs", "Comma-delimited list of chromosomes to process.", false);

    cli.parseWithExit(args);

    NGSBinSNPSelector selector = new NGSBinSNPSelector();

    selector.outputDir = cli.get("out");
    selector.inputDir = cli.get("in");
    selector.format = cli.get("format");
    selector.token = cli.get("token");
    selector.mappingQualityFile = cli.has("mapQ") ? cli.get("mapQ") : null;
    selector.afTag = cli.has("af") ? cli.get("af") : "AF";

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

    selector.bins = new ReferenceGenome(GENOME_BUILD.HG19, new Logger()).getBins(1000, chrs);
    selector.run();
  }

}
