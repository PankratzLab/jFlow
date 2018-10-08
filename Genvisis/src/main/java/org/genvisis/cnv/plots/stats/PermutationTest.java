package org.genvisis.cnv.plots.stats;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.apache.commons.cli.ParseException;
import org.genvisis.cnv.filesys.CNVariant;
import org.genvisis.cnv.plots.PlotUtilities;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomicPosition;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.CLI;
import org.pankratzlab.shared.collect.RangeMultimap;
import org.pankratzlab.shared.collect.TreeRangeSetMultimap;
import org.pankratzlab.shared.stats.LeastSquares;
import org.pankratzlab.shared.stats.MpermResults;
import org.pankratzlab.shared.stats.MpermResults.TestType;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Range;

public class PermutationTest {

  public enum Filter {
    DELETIONS {

      public boolean include(int cn) {
        return cn < 0;
      }
    },
    DUPLICATIONS {

      public boolean include(int cn) {
        return cn > 0;
      }
    },
    ALL {

    },
    HOMOZYGOUS_DELETION {

      public boolean include(int cn) {
        return cn == -2;
      }
    };

    public boolean include(int cn) {
      return cn != 0;
    }
  }

  /**
   * Compute and write out the results of a Permutation Test over the given CNVs and phenotype
   * 
   * @param cnvfile File containing CNV information
   * @param mapfile .map file containing regions to evaluate (form: chr/name/start/stop)
   * @param pheno Map from Sample ID to phenotype value
   * @param windowSize Size of window to check around regions
   * @param mperm Number of permutations to run
   * @param outfile Name of the output file
   * @param twoSided Whether to run a two-sided test (applicable for plink only)
   * @param filter CNV filter to apply (eg deletions, duplications, etc)
   * @param type Type of test to run (eg plink, Mann-Whitney)
   */
  public static void mperm(String cnvfile, String mapfile, Map<String, Double> pheno,
                           int windowSize, int mperm, String outfile, boolean twoSided,
                           Filter filter, MpermResults.TestType type, int numThreads, Logger log) {
    long time;
    Map<String, Integer> sampleSet = new HashMap<String, Integer>();

    // save keys and indices so we have a consistent index into our table
    List<Double> valueList = new ArrayList<Double>();
    int count = 0;
    for (String k : pheno.keySet()) {
      sampleSet.put(k, count++);
      valueList.add(pheno.get(k));
    }

    time = new Date().getTime();
    // read in cnv file
    CNVariant[] cnvs = CNVariant.loadPlinkFile(cnvfile);
    log.report("Read CNVariants in " + ext.getTimeElapsed(time));

    time = new Date().getTime();
    String[][] mat;
    // create list of regions to check using the given map file and window size
    if (mapfile == null) {
      List<String[]> tempMat = new ArrayList<String[]>();
      Set<String> seen = new HashSet<String>();
      for (int i = 0; i < cnvs.length; i++) {
        CNVariant c = cnvs[i];
        String[] start = {c.getChr() + "", c.getChr() + ":" + c.getStart(), "0", c.getStart() + ""};
        String[] stop = {c.getChr() + "", c.getChr() + ":" + c.getStop(), "0", c.getStop() + ""};
        String[] onePlus = {c.getChr() + "", c.getChr() + ":" + (c.getStop() + 1), "0",
                            c.getStop() + 1 + ""};

        // check if we've seen this named region before, if not, add it to our matrix
        if (!seen.contains(start[1])) {
          seen.add(start[1]);
          tempMat.add(start);
        }
        if (!seen.contains(stop[1])) {
          seen.add(stop[1]);
          tempMat.add(stop);
        }
        if (!seen.contains(onePlus[1])) {
          seen.add(onePlus[1]);
          tempMat.add(onePlus);
        }
      }
      mat = Matrix.toStringArrays(tempMat);
    } else {
      mat = HashVec.loadFileToStringMatrix(mapfile, false, null);
    }
    TreeRangeSetMultimap<GenomicPosition, Integer> regions = TreeRangeSetMultimap.create();
    // If the first position column is 0, we're looking at a single position. Otherwise, we're looking at a gene/region and need to include both the start and stop index (2 and 3)
    int startIndex = mat[0][2].equals("0") ? 3 : 2;
    for (int i = 0; i < mat.length; i++) {
      regions.put(Range.closed(new GenomicPosition(Byte.parseByte(mat[i][0]),
                                                   Integer.parseInt(mat[i][startIndex])
                                                                              - windowSize),
                               new GenomicPosition(Byte.parseByte(mat[i][0]),
                                                   Integer.parseInt(mat[i][3]) + (windowSize == 0 ? 1
                                                                                                  : windowSize))),
                  i);
    }
    log.report("Read map in " + ext.getTimeElapsed(time));

    time = new Date().getTime();
    int currentSampleIndex;
    int[][] cnMatrix = new int[mat.length][sampleSet.size()];
    for (int i = 0; i < cnvs.length; i++) {
      // determine index of sample to use as column index in the resulting matrix 
      String id = cnvs[i].getFamilyID() + "\t" + cnvs[i].getIndividualID();
      if (sampleSet.containsKey(id)) {
        currentSampleIndex = sampleSet.get(id);

        RangeMultimap<GenomicPosition, Integer, ImmutableSet<Integer>> overlap = regions.subRangeMap(Range.closed(new GenomicPosition(cnvs[i].getChr(),
                                                                                                                                      cnvs[i].getStart()),
                                                                                                                  new GenomicPosition(cnvs[i].getChr(),
                                                                                                                                      cnvs[i].getStop())));
        for (Range<GenomicPosition> k : overlap.keySet()) {
          for (Integer j : overlap.get(k)) {
            cnMatrix[j][currentSampleIndex] = filter.include(cnvs[i].getCN() - 2) ? 1 : 0;
          }
        }
      }
    }

    log.report("Built PED matrix in " + ext.getTimeElapsed(time));

    MpermResults[] out = new MpermResults[cnMatrix.length];
    String[] outHeader = new String[] {"CHR", "BP", "SNP", "NCNV", "M1", "M0", "EMP1", "EMP2"};
    time = new Date().getTime();
    // make a results object for each region we're testing
    for (int i = 0; i < cnMatrix.length; i++) {
      out[i] = new MpermResults(Byte.parseByte(mat[i][0]), Integer.parseInt(mat[i][3]), mat[i][1],
                                valueList, cnMatrix[i], mperm, twoSided, type);
    }
    log.report("Computed initial results in " + ext.getTimeElapsed(time));

    time = new Date().getTime();

    ExecutorService pool = Executors.newFixedThreadPool(numThreads);
    for (int i = 0; i < mperm; i++) {
      Collections.shuffle(valueList);
      List<Double> vals = new ArrayList<Double>(valueList);
      pool.execute(new Runnable() {

        public void run() {
          double bestStat = type.worst();

          // run permutation for this shuffled phenotype value list
          for (int j = 0; j < out.length; j++) {
            MpermResults res = out[j];
            synchronized (res) {
              double stat_j = res.permute(vals);

              // track the maximum value we've seen this permutation 
              bestStat = type.betterStat(stat_j, bestStat);
            }
          }

          // increment r2, the number of greater differences across all permutations
          for (int j = 0; j < out.length; j++) {
            MpermResults res = out[j];
            synchronized (res) {
              res.updateR2(bestStat);
            }
          }
        }
      });
    }

    pool.shutdown();
    try {
      pool.awaitTermination(10, TimeUnit.MINUTES);
    } catch (InterruptedException e) {

      e.printStackTrace();
    }

    log.report("Finished permutations in " + ext.getTimeElapsed(time));

    PrintWriter writer = Files.getAppropriateWriter(outfile);
    writer.println(ArrayUtils.toStr(outHeader));

    // print results
    for (int i = 0; i < out.length; i++) {
      writer.println(out[i]);
    }

    writer.close();
  }

  /**
   * Create the map from sample ID to phenotype, regressing out covariates as necessary. Assumes
   * phenoFile and covarFile are in the same ID order (eg the output of PhenoPrep)
   */
  public static Map<String, Double> prepPheno(String phenoFile, String covarFile, Logger log) {
    Map<String, Double> residuals = new HashMap<String, Double>();
    double[] deps;
    double[][] indeps;
    List<String> ids = new ArrayList<String>();
    List<Double> pheno = new ArrayList<Double>();
    List<double[]> covars = new ArrayList<double[]>();
    double[] res;

    try {
      BufferedReader phenoReader = Files.getAppropriateReader(phenoFile);

      String phenoHeader = phenoReader.readLine();

      String delim = ext.determineDelimiter(phenoHeader);
      String[] line = phenoHeader.split(delim);
      boolean fidiid = line[0].equals("FID") && line[1].equals("IID");

      while (phenoReader.ready()) {
        line = phenoReader.readLine().split(delim);
        String id = line[0] + (fidiid ? "\t" + line[1] : "");
        ids.add(id);
        pheno.add(Double.parseDouble(line[line.length - 1]));
      }

      deps = new double[pheno.size()];
      int i = 0;
      for (Double d : pheno) {
        deps[i++] = d;
      }

      if (covarFile != null) {
        BufferedReader covarReader = Files.getAppropriateReader(covarFile);
        String covarHeader = covarReader.readLine();

        delim = ext.determineDelimiter(covarHeader);
        line = covarHeader.split(delim);

        // do we have the same ID column headers? If not, abort
        if (fidiid ^ (line[0].equals("FID") && line[1].equals("IID"))) {
          log.reportError("ID columns do not match. Aborting.");
          return null;
        } else {
          int rangeStart = fidiid ? 2 : 1;
          while (covarReader.ready()) {
            line = covarReader.readLine().split(delim);
            line = ArrayUtils.subArray(line, rangeStart);
            covars.add(ArrayUtils.toDoubleArray(line));
          }

          indeps = Matrix.toDoubleArrays(covars);

          LeastSquares reg = new LeastSquares(deps, indeps, null, false, true);
          res = reg.getResiduals();
        }
      } else {
        res = deps;
      }

      if (ids.size() != res.length) {
        log.reportError("Error - lost a few rows in regression model. Aborting.");
        return null;
      }

      // map ids to residualized phenotype value
      for (int j = 0; j < res.length; j++) {
        residuals.put(ids.get(j), res[j]);
      }

    } catch (IOException e) {
      e.printStackTrace();
    }

    return residuals;
  }

  public static void main(String[] args) throws ParseException {
    CLI cli = new CLI("PermutationTest");
    cli.addArg("phenoFile", "phenotype file");
    cli.addArg("covarFile", "file containing covariates");
    cli.addArg("outFile", "output file name");
    cli.addArg("cnvFile", "cnv file name");
    cli.addArg("mapFile", "map file name");
    cli.addArgWithDefault("window", "window size in kb", 0);
    cli.addArg("mperm", "number of permutations");
    cli.addFlag("twoSided", "Run a two-sided test (plink only)");
    cli.addArgWithDefault("filter",
                          "CNV filter to be applied (eg deletions, duplications, homozygousDeletions, all [default])",
                          "all");
    cli.addArgWithDefault("type", "Test type (eg \"plink\", Mann-Whitney (\"mw\"))", "plink");
    cli.addArgWithDefault("threads", "number of threads to use", 1);
    cli.addFlag("parseResults",
                "run hitWindows, manhattan plot, and qq plot on the resulting file");

    cli.parse(args);
    Filter filter;
    String f = cli.get("filter");
    if (f.equals("deletions")) {
      filter = Filter.DELETIONS;
    } else if (f.equals("duplications")) {
      filter = Filter.DUPLICATIONS;
    } else if (f.equals("homozygousDeletions")) {
      filter = Filter.HOMOZYGOUS_DELETION;
    } else {
      filter = Filter.ALL;
    }

    Logger log = new Logger("mperm.log");
    Map<String, Double> pheno = prepPheno(cli.get("phenoFile"), cli.get("covarFile"), log);
    if (pheno == null) {
      log.reportError("Failed to create phenotye map.");
      return;
    }

    MpermResults.TestType type;
    String t = cli.get("type");
    if (t.equals("logistic")) {
      type = MpermResults.TestType.LOGISTIC;
      // Logistic regression is currently inefficient, and does not handle covars correctly as we're regressing them out by default
      // It is disabled for now.
      log.reportError("Logistic regression is not currently implemented. Consider running as plink or mw.");
      return;
    } else if (t.equals("mw")) {
      type = MpermResults.TestType.MW;
    } else {
      type = MpermResults.TestType.PLINK;
    }

    mperm(cli.get("cnvFile"), cli.get("mapFile"), pheno, cli.getI("window"), cli.getI("mperm"),
          cli.get("outFile"), cli.has("twoSided"), filter, type, cli.getI("threads"), log);
    if (cli.has("parseResults")) {
      String outFile = cli.get("outFile");
      PlotUtilities.runHitWindows(outFile);

      PlotUtilities.createManPlotScreenshot(outFile);

      PlotUtilities.createQQPlotScreenshot(outFile);

    }
  }
}
