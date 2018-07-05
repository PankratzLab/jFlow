package org.genvisis.one.JL.topMed;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.seq.manage.VCFOps;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * 
 *
 */
public class SummarizePopAF {

  private static final String MAFS = "mafs";
  private static final String ANN = "ANN";
  private static final String ANN_DELIM = ",";
  private static final String ANN_DELIM_SUB = "\\|";
  private static final String ANN_BLANK = ".";

  private static final int IMPACT_INDEX = 2;
  private static final int GENE_INDEX = 3;

  /**
   * If this is used elsewhere, can make this customizable
   */
  private static final List<String> FREQ_ANNOS = Arrays.asList("AF_EM_POP_DEF_African_Americans",
                                                               "AF_EM_POP_DEF_EAS",
                                                               "AF_EM_POP_DEF_Hispanics",
                                                               "AF_EM_POP_DEF_SAS",
                                                               "AF_EM_POP_DEF_Whites", "AF");

  private enum IMPACT {
    MODIFIER, LOW, MODERATE, HIGH;
  }

  private static class MAFGeneCounter {

    private final Map<String, Map<IMPACT, Integer>> countMap;
    private final double maf;

    /**
     * @param maf
     */
    private MAFGeneCounter(String population, double maf) {
      super();
      this.countMap = new HashMap<>();
      this.maf = maf;
    }

    private void init(String gene) {
      if (!countMap.containsKey(gene)) {
        EnumMap<IMPACT, Integer> subMap = new EnumMap<>(IMPACT.class);
        for (IMPACT im : IMPACT.values()) {
          subMap.put(im, 0);
        }
        countMap.put(gene, subMap);
      }
    }

    private void count(IMPACT impact, String gene, double af) {

      if (af < maf) {
        int current = countMap.get(gene).get(impact);
        countMap.get(gene).put(impact, current + 1);
      }
    }

  }

  private static void run(CLI c) {

    String outDir = c.get(CLI.ARG_OUTDIR);
    String vcf = c.get(CLI.ARG_VCF);
    double[] mafs = ArrayUtils.toDoubleArray(c.get(MAFS).split(","));
    Map<String, List<MAFGeneCounter>> mCounters = new HashMap<>();
    for (String pop : FREQ_ANNOS) {
      mCounters.put(pop, new ArrayList<>());
      for (double d : mafs) {
        mCounters.get(pop).add(new MAFGeneCounter(pop, d));
      }
    }
    new File(outDir).mkdirs();

    Logger log = new Logger(outDir + "log.log");

    VCFFileReader reader = new VCFFileReader(new File(vcf), true);
    int numScanned = 0;
    HashSet<String> allGenes = new HashSet<>();
    for (VariantContext vc : reader) {
      numScanned++;
      if (numScanned % 100000 == 0) {
        log.reportTimeInfo(Integer.toString(numScanned) + " variants scanned");
      }
      String ann = vc.getAttributeAsString(ANN, ANN_BLANK);
      if (!ANN_BLANK.equals(ann)) {
        String[] anns = ann.split(ANN_DELIM);
        IMPACT maxImpact = null;
        String gene = null;
        for (String a : anns) {
          String[] subAnn = a.split(ANN_DELIM_SUB);
          IMPACT impact = IMPACT.valueOf(subAnn[IMPACT_INDEX]);
          if (maxImpact == null || impact.ordinal() > maxImpact.ordinal()) {

            maxImpact = impact;
            gene = subAnn[GENE_INDEX];
          }
        }

        allGenes.add(gene);
        double maxAF = -1;
        for (String afAnno : FREQ_ANNOS) {
          double af = vc.getAttributeAsDouble(afAnno, -1);
          if (af > maxAF) {
            maxAF = af;
          }
        }
        for (String afAnno : FREQ_ANNOS) {
          double af = vc.getAttributeAsDouble(afAnno, -1);
          for (MAFGeneCounter mCounter : mCounters.get(afAnno)) {
            mCounter.init(gene);
            if (af > 0) {
              mCounter.count(maxImpact, gene, maxAF);
            }
          }
        }
      }
    }
    reader.close();

    String outFile = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".geneCounts.txt";

    StringJoiner full = new StringJoiner("\n");
    StringJoiner header = new StringJoiner("\t");
    header.add("GENE");
    for (String pop : mCounters.keySet()) {
      for (MAFGeneCounter mCounter : mCounters.get(pop)) {
        for (IMPACT im : IMPACT.values()) {
          String reportPop = pop.equals("AF") ? "AF_ALL" : pop;
          header.add("POP_" + reportPop + "_MAF_" + mCounter.maf + "_IMPACT_" + im + "_COUNT");

        }
      }
    }
    full.add(header.toString());
    for (String gene : allGenes) {
      StringJoiner geneJoiner = new StringJoiner("\t");
      geneJoiner.add(gene);
      for (String pop : mCounters.keySet()) {
        for (MAFGeneCounter mCounter : mCounters.get(pop)) {
          for (IMPACT im : IMPACT.values()) {
            geneJoiner.add(Integer.toString(mCounter.countMap.get(gene).get(im)));
          }
        }
      }
      full.add(geneJoiner.toString());
    }
    log.reportTimeInfo("Writing to " + outFile);
    Files.write(full.toString(), outFile);

  }

  public static void main(String[] args) {
    CLI c = new CLI(SummarizePopAF.class);

    c.addArg(CLI.ARG_VCF, CLI.DESC_VCF);

    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR);
    c.addArg(MAFS, "comma delimited list of mafs to summarize, note that maf 0 will be singletons");

    c.parseWithExit(args);

    run(c);

  }

}
