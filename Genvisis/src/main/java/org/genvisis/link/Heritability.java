package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.CountHash;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.FamilyStructure;
import org.genvisis.parse.GenParser;
import org.genvisis.stats.Correlation;
import org.genvisis.stats.ICC;
import org.genvisis.stats.RegressionModel;

import com.google.common.primitives.Doubles;

public class Heritability {
  public static final String DEFAULT_MERLIN_EXEC = "merlin";
  public static final String DEFAULT_SOLAR_EXEC = "/home/pankrat2/public/bin/solarEclipse/solar";

  private static String computeWithMerlin(String dir, String pedfile, String pheno, String covars,
                                          String prefix, String merlinExec, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String> phenoHash, covarHash, seen, referenced;
    String[] keys, covarHeader;
    String temp;
    String estimate;

    estimate = ".";
    dir = ext.verifyDirFormat(dir);

    line = Files.getHeaderOfFile(dir + pheno, log);
    phenoHash = HashVec.loadFileToHashString(dir + pheno, new int[] {0, 1}, new int[] {2},
                                             pheno.endsWith(".csv"), "\t", true, false, false);
    covarHeader = null;
    if (covars != null) {
      // System.err.println("Error - covars not implemented for Merlin yet, will eventually
      // construct a residual and run that instead; though it looks like it may work in the .dat
      // file");
      covarHeader = Files.getHeaderOfFile(dir + covars, log);
      if (!covarHeader[0].equals("FID") || !covarHeader[1].equals("IID")) {
        System.err.println("Error - warning covariates file expected to start with FID and IID; may not merge properly without");
      }
      covarHash = HashVec.loadFileToHashString(dir + covars, new int[] {0, 1},
                                               Array.subArray(Array.arrayOfIndices(covarHeader.length),
                                                              2),
                                               covars.endsWith(".csv"), "\t", true, false, false);
      covarHeader = Array.subArray(covarHeader, 2);
    } else {
      covarHash = new Hashtable<String, String>();
      covarHeader = new String[0];
    }

    seen = new Hashtable<String, String>();
    referenced = new Hashtable<String, String>();
    try {
      reader = new BufferedReader(new FileReader(pedfile));
      reader.mark(10000);
      line = reader.readLine().trim().split("[\\s]+");
      if (FamilyStructure.likelyPedHeader(line)) {
        if (log.getLevel() > 8) {
          log.report("Header row detected in the pedigree file");
        }
      } else {
        if (log.getLevel() > 8) {
          log.report("No header row detected in the pedigree file");
        }
        reader.reset();
      }
      writer = new PrintWriter(new FileWriter(dir + "re_chrom01.pre"));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        writer.println(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4]
                       + "\t"
                       + (phenoHash.containsKey(line[0] + "\t"
                                                + line[1]) ? phenoHash.get(line[0] + "\t" + line[1])
                                                           : "x")
                       + "\t"
                       + (covarHash.containsKey(line[0] + "\t"
                                                + line[1]) ? covarHash.get(line[0] + "\t" + line[1])
                                                           : Array.toStr(Array.stringArray(covarHeader.length,
                                                                                           "x")))
                       + "\t1\t1");
        seen.put(line[0] + "\t" + line[1], "");
        if (!line[2].equals("0")) {
          referenced.put(line[0] + "\t" + line[2], "1");
        }
        if (!line[3].equals("0")) {
          referenced.put(line[0] + "\t" + line[3], "2");
        }
      }
      reader.close();

      keys = HashVec.getKeys(referenced);
      for (int i = 0; i < keys.length; i++) {
        if (!seen.containsKey(keys[i])) {
          writer.println(keys[i] + "\t0\t0\t" + referenced.get(keys[i]) + "\t"
                         + (phenoHash.containsKey(keys[i]) ? phenoHash.get(keys[i]) : "x")
                         + "\t0\t0");
        }
      }

      writer.close();
      new LinkageMap(1, new String[] {"rs007"}, 2, new double[] {1}, false,
                     false).createFileInDir(dir);
      Merlin.createMerlinFiles(dir, 1, dir + "chr01", Merlin.QUANTITATIVE_TRAIT, covarHeader);
      // System.exit(1);
      // if (Files.exists("map01.dat")) {
      // new File("map01.dat").delete();
      // }

      try {
        CmdLine.run(merlinExec + " -d chr01.dat -p re_chrom01.pre -m chr01.map --vc"
                    + (covarHeader.length > 0 ? " --useCovariates" : "")
                    + " --bits 1000 --megabytes 8000 --tabulate --step 5 --markerNames --information --prefix "
                    + prefix, dir, new PrintStream(new File(dir + prefix + "_merlin.log")));
      } catch (IOException ioe) {
        log.reportError("Error - merlin failed to run/not available");
      }

      try {
        reader = new BufferedReader(new FileReader(dir + prefix + "_merlin.log"));
        while (reader.ready()) {
          temp = reader.readLine();
          if (temp.contains("h2 = ")) {
            estimate = temp.substring(temp.indexOf("h2 = ") + 5, temp.length() - 1);
            log.report("Merlin estimate: " + estimate);
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        log.reportError("Error - merlin failed to run");
        fnfe.printStackTrace();
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + dir + prefix + ".log"
                        + "\"; merlin probably ran into trouble");
      }
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + pedfile + "\" not found in current directory");
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + pedfile + "\"");
    }

    return estimate;
  }

  private static String[] computeWithSolar(String dir, String pedfile, String pheno, String covars,
                                           String prefix, String solarExec, Logger log) {
    BufferedReader reader;
    String temp;
    String[] ids, params, line, covarsHeader;
    String trait;

    dir = ext.verifyDirFormat(dir);

    GenParser.parse(new String[] {pedfile, "out=" + dir + "/" + prefix + "_fam.csv", "0=FAMID",
                                  "1=ID", "2=FA", "3=MO", "4=SEX"},
                    log);
    ids = HashVec.loadFileToStringArray(pedfile, false, new int[] {1}, false);
    if (Array.unique(ids).length != ids.length) {
      System.err.println("Error - lookup for solar heritabilty screen currently requires unique IDs (found "
                         + (ids.length - Array.unique(ids).length)
                         + " that were not unique); aborting solar run");
      return null;
    }

    params = new String[covars == null ? 2 : 3];
    params[0] = pedfile + " 1 0=FAMID 1=ID skip=0";
    params[1] = dir + pheno + " 1 2";
    trait = Files.getHeaderOfFile(dir + pheno, log)[2];
    if (covars != null) {
      covarsHeader = Files.getHeaderOfFile(dir + covars, log);
      line = Array.toStringArray(Array.arrayOfIndices(covarsHeader.length));
      line = Array.removeFromArray(line, 0);
      params[2] = dir + covars + " " + Array.toStr(line, " ");
    } else {
      covarsHeader = null;
    }
    Files.combine(ids, params, null, null, "", dir + "/" + prefix + "_ptypes.csv", log, true, true,
                  true);
    Solar.cleanPhenotypeFile(dir + "/" + prefix + "_ptypes.csv");

    // Files.write("echo -e \"load ped "+prefix+"_fam.csv\\nautomodel "+prefix+"_ptypes.csv
    // "+trait+"\\npolygenic -screen\\nquit\\n\" | "+solarExec+" > "+prefix+"_solar.log",
    // dir+"/batch");
    Files.write("echo -e \"load ped " + prefix + "_fam.csv\\nload phenotype "
                + prefix + "_ptypes.csv\\ntrait " + trait + "\\n" + (covarsHeader != null
                                                                                          ? "covariates "
                                                                                            + Array.toStr(Array.subArray(covarsHeader,
                                                                                                                         2),
                                                                                                          " ")
                                                                                            + "\\n"
                                                                                          : "")
                + "polygenic -screen\\nquit\\n\" | " + solarExec + " > " + prefix + "_solar.log",
                dir + "/batch");
    System.out.println(dir);
    if (!Files.isWindows()) {
      CmdLine.run("chmod +x batch", dir);
      CmdLine.run("./batch", dir);
    }
    String[] solarEstimate = new String[] {"NaN", "NaN", "NaN", "NaN", "NaN"};
    try {
      reader = new BufferedReader(new FileReader(dir + trait + "/polygenic.out"));
      while (reader.ready()) {
        temp = reader.readLine();
        if (temp.contains("H2r is ")) {
          String estString = temp.substring(temp.indexOf("H2r is ") + 7);
          log.report("Solar  estimate: " + estString);
          String[] toParse = estString.trim().split("[\\s]+");
          try {
            double h2 = Double.parseDouble(toParse[0]);
            solarEstimate[0] = h2 + "";
            double pval = Double.parseDouble(toParse[3]);
            solarEstimate[1] = pval + "";
          } catch (NumberFormatException nfe) {
            log.reportTimeError("Could not parse " + estString);
            reader.close();
            return null;
          }

          temp = reader.readLine();
          while (!temp.contains("Loglikelihoods")) {
            log.report(temp);
            if (temp.contains("H2r Std. Error:")) {
              String[] tmpString = temp.trim().split("[\\s]+");
              try {
                double stError = Double.parseDouble(tmpString[3]);
                solarEstimate[2] = stError + "";
              } catch (NumberFormatException nfe) {
                log.reportTimeError("Could not parse " + temp);
                reader.close();
                return null;
              }
            }
            temp = reader.readLine();
          }

          while (reader.ready()) {
            temp = reader.readLine();
            if (temp.contains("Kurtosis")) {
              if (temp.contains("Residual Kurtosis is")) {
                String kurtString = temp.replaceAll(".*Residual Kurtosis is ", "");
                kurtString = kurtString.replaceAll(" .*", "");
                kurtString = kurtString.replaceAll(",.*", "");
                try {
                  double kurt = Double.parseDouble(kurtString);
                  solarEstimate[3] = kurt + "";

                } catch (NumberFormatException nfe) {
                  log.reportTimeError("Could not parse " + temp);
                  reader.close();
                  return null;
                }
                if (temp.contains("high")) {
                  solarEstimate[4] = "HIGH";
                } else if (temp.contains("normal")) {
                  solarEstimate[4] = "NORMAL";
                } else {
                  log.reportTimeError("Could not kurtosis type for " + temp);
                }
              }
              log.report(temp);
            }
          }
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error - solar failed to run");
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + trait
                      + "/polygenic.out\"; solar probably ran into trouble");
    }
    return solarEstimate;
  }

  public static void getHeritabilitiesOfAllPhenosInADir(String dir, String pedigreeFile,
                                                        String covars, Logger log) {
    String[] phenos, results;
    String merlinExec;

    merlinExec = DEFAULT_MERLIN_EXEC;
    phenos = Files.list(dir, ".txt", false);
    results = new String[phenos.length];
    for (int i = 0; i < phenos.length; i++) {
      results[i] = computeWithMerlin(dir, pedigreeFile, phenos[i], covars, ext.rootOf(phenos[i]),
                                     merlinExec, log);
    }
    Files.writeArray(results, dir + "heritabilities.txt");
  }

  public static void fromParameters(String filename, boolean skipExtra, Logger log) {
    PrintWriter writer, summary;
    String[] line;
    Hashtable<String, String> famIdHash;
    Vector<String> params;
    Vector<String[]> models;
    String pedigreeFile, dbFile;
    String merlinExec, solarExec;
    String temp;
    String[][] data;
    String root;
    String dbDelimiter;
    String[] dbHeader;
    boolean[] use;
    int numNotInPed;
    int[] indices;
    CountHash counter;
    String dir;
    String merlinEstimate;
    String[] solarEstimate;
    int numOfAllSamples;
    int numOfFamiliesSizedTwoOrAbove;
    int numOfFamiliesSizedOne;

    dir = ext.parseDirectoryOfFile(filename);

    params = Files.parseControlFile(filename, "heritability",
                                    new String[] {"ped=plink.fam", "db=input.xln",
                                                  "# name of model (will be used to make a subdirectory) and dependent variable, followed by all covariates to be included",
                                                  "model1 pta age sex PC1 PC2",
                                                  "model2 icam1 age sex center", "",
                                                  "# optional specifications of program locations; simply correct and un-comment",
                                                  "# merlin_exec=/bin/merlin",
                                                  "# solar_exec=/bin/solar"},
                                    log);

    if (params != null) {
      pedigreeFile = null;
      dbFile = null;
      models = new Vector<String[]>();
      merlinExec = DEFAULT_MERLIN_EXEC;
      solarExec = DEFAULT_SOLAR_EXEC;
      for (int i = 0; i < params.size(); i++) {
        temp = params.elementAt(i);

        if (temp.startsWith("ped=")) {
          pedigreeFile = ext.parseStringArg(temp, "failed_to_define_pedigree_file");
        } else if (temp.startsWith("db=")) {
          dbFile = ext.parseStringArg(temp, "failed_to_define_db_file");
        } else if (temp.startsWith("merlin_exec=")) {
          merlinExec = ext.parseStringArg(temp, DEFAULT_MERLIN_EXEC);
        } else if (temp.startsWith("solar_exec=")) {
          solarExec = ext.parseStringArg(temp, DEFAULT_SOLAR_EXEC);
        } else {
          models.add(temp.trim().split("[\\s]+"));
        }
      }

      if (!Files.exists(pedigreeFile)) {
        log.report("Error - could not find pedigree file '" + pedigreeFile + "'; aborting");
        return;
      }

      if (!Files.exists(dbFile)) {
        log.report("Error - could not find database file '" + dbFile + "'; aborting");
        return;
      }

      famIdHash = HashVec.loadFileToHashString(pedigreeFile, new int[] {1}, new int[] {0},
                                               Files.determineDelimiter(pedigreeFile, log)
                                                    .equals(","),
                                               null, false, false, false);

      log.setLevel(8);

      try {
        summary = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + "_summary.xln"));
        summary.println("Model\tMerlin_est.\tSolar_est.\tSolar_p\tSolar_StdError\tSolar_Kurt\tSolar_KurtWarning\tn_Samples\tn_Families\tn_Families_size>1\tAverage_size_families_siez>1\tn_Families_size=1\tPearson_Correl_Z_Siblings\tPearson_Correl_Z_ParentOffspringPairs\tPearson_Correl_Z_Trios\tPearson_Correl_Pval_Siblings\tPearson_Correl_Pval_ParentOffspringPairs\tPearson_Correl_Pval_Trios\tICC_Siblings\tICC_ParentOffspringPairs\tICC_Trios\tN_Siblings\tN_ParentOffspringPairs\tN_Trios");
        for (int i = 0; i < models.size(); i++) {
          line = models.elementAt(i);
          if (line.length < 2) {
            log.reportError("Error - need at least two arguments (model->directory and the model, even if there are no more covariates), found:");
            log.reportError(Array.toStr(line));
          } else {
            root = line[0];
            root = ext.replaceWithLinuxSafeCharacters(root, false);

            new File(dir + root).mkdirs();

            dbDelimiter = Files.determineDelimiter(dbFile, log);
            dbHeader = Files.getHeaderOfFile(dbFile, dbDelimiter, log);
            line[0] = dbHeader[0];
            indices = ext.indexFactors(line, dbHeader, false, log, true, false);
            if (Array.min(indices) == -1) {
              summary.close();
              return;
            }
            data = HashVec.loadFileToStringMatrix(dbFile, true, indices, dbDelimiter, false, 1000,
                                                  false);
            //FIXME this is guaranteed to cause a NPE. What is supposed to be computed here?
            use = RegressionModel.getRowsWithCompleteData(null,
                                                          Matrix.prune(data, null,
                                                                       Array.subArray(Array.arrayOfIndices(line.length),
                                                                                      1),
                                                                       log),
                                                          log);

            numNotInPed = 0;
            counter = new CountHash();
            HashSet<String> validIDs = new HashSet<String>();
            Vector<String> deps = new Vector<String>();
            Vector<double[]> indeps = new Vector<double[]>();
            int count = 0;
            HashMap<String, Integer> subIndexMap = new HashMap<String, Integer>();
            try {
              System.out.println(dir + root + "/pheno.dat");
              writer = new PrintWriter(new FileWriter(dir + root + "/pheno.dat"));
              writer.println("FID\tIID\t" + line[1]);
              for (int j = 0; j < use.length; j++) {
                if (use[j]) {
                  if (famIdHash.containsKey(data[j][0])) {
                    String fidiid = famIdHash.get(data[j][0]) + "\t" + data[j][0];
                    writer.println(fidiid + "\t" + data[j][1]);
                    counter.add(famIdHash.get(data[j][0]));
                    validIDs.add(fidiid);
                    deps.add(data[j][1]);
                    subIndexMap.put(fidiid, count++);
                  } else {
                    if (numNotInPed == 0) {
                      log.reportError("Warning - the following samples were not found in the pedigree file:");
                    }
                    if (numNotInPed == 10) {
                      log.reportError("...");
                    } else if (numNotInPed < 10) {
                      log.reportError(data[j][0]);
                    }
                    numNotInPed++;
                    use[j] = false;
                  }
                }
              }
              writer.close();
            } catch (Exception e) {
              System.err.println("Error writing to " + root + "/pheno.dat");
              e.printStackTrace();
            }
            // counter.getCounts();

            if (line.length > 2) {
              try {
                writer = new PrintWriter(new FileWriter(root + "/covars.dat"));
                writer.println("FID\tIID\t" + Array.toStr(Array.subArray(line, 2)));
                for (int j = 0; j < use.length; j++) {
                  if (use[j]) {
                    writer.println(famIdHash.get(data[j][0]) + "\t" + data[j][0] + "\t"
                                   + Array.toStr(Array.subArray(data[j], 2)));
                    indeps.add(Array.toDoubleArray(Array.subArray(data[j], 2)));
                  }
                }
                writer.close();
              } catch (Exception e) {
                System.err.println("Error writing to " + root + "/pheno.dat");
                e.printStackTrace();
              }

              if (numNotInPed > 0) {
                log.reportError("There were " + numNotInPed
                                + " sample(s) with valid and complete phenotype data that were not in the pedigree file");
              }
            }
            double sibICC = Double.NaN;
            double poICC = Double.NaN;
            double trioICC = Double.NaN;
            double[] sibCorrel = Array.doubleArray(2, Double.NaN);
            double[] poCorrel = Array.doubleArray(2, Double.NaN);
            double[] trioCorrel = Array.doubleArray(2, Double.NaN);
            int cntSib = 0;
            int cntPO = 0;
            int cntTrio = 0;
            if (!skipExtra) {
              FamilyStructure ped = new FamilyStructure(pedigreeFile, false);
              ArrayList<String[]> sibList = Pedigree.PedigreeUtils.loadSibs(ped, true, null,
                                                                            validIDs, true); // double
                                                                                             // actual
                                                                                             // due
                                                                                             // to
                                                                                             // bidirectionality
              ArrayList<String[]> poPairs = Pedigree.PedigreeUtils.loadPOPairs(ped, true, null,
                                                                               validIDs, true);
              ArrayList<int[]> trios = Pedigree.PedigreeUtils.loadCompleteTrios(ped, null, validIDs,
                                                                                true);

              double[] resids = RegressionModel.processDeps(deps);
              if (indeps.size() > 0) {
                RegressionModel model = RegressionModel.determineAppropriate(resids,
                                                                             RegressionModel.processIndeps(indeps),
                                                                             false, true);
                resids = model.getResiduals();
              }


              /* sibs : */ {
                double[][] correlationData;
                ArrayList<Double> sibICCData = new ArrayList<Double>();
                ArrayList<String> sibICCResponseIDs = new ArrayList<String>();
                ArrayList<Double> correl1 = new ArrayList<Double>();
                ArrayList<Double> correl2 = new ArrayList<Double>();
                HashSet<String> used = new HashSet<String>();
                for (int k = 0; k < sibList.size(); k++) {
                  if (used.contains(sibList.get(k)[0] + "\t" + sibList.get(k)[1])
                      || used.contains(sibList.get(k)[0] + "\t" + sibList.get(k)[1])) {
                    // already found
                    continue;
                  }
                  Integer resid1 = subIndexMap.get(sibList.get(k)[0]);
                  Integer resid2 = subIndexMap.get(sibList.get(k)[1]);
                  if (resid1 == null || resid2 == null) {
                    continue; // skip
                  }
                  cntSib++;
                  used.add(sibList.get(k)[0] + "\t" + sibList.get(k)[1]);
                  correl1.add(resids[resid1.intValue()]);
                  correl2.add(resids[resid2.intValue()]);
                  sibICCData.add(resids[resid1]);
                  sibICCData.add(resids[resid2]);
                  sibICCResponseIDs.add(sibList.get(k)[0] + "\t" + sibList.get(k)[1]);
                  sibICCResponseIDs.add(sibList.get(k)[0] + "\t" + sibList.get(k)[1]);
                }
                correlationData =
                                new double[][] {Doubles.toArray(correl1), Doubles.toArray(correl2)};
                sibCorrel = Correlation.Pearson(correlationData);
                ICC sibICCAnalysis = new ICC(Doubles.toArray(sibICCData),
                                             Array.toStringArray(sibICCResponseIDs), null, null,
                                             true, log);
                sibICCAnalysis.computeICC();
                sibICC = sibICCAnalysis.getICC();
              }
              /* po : */ {
                double[][] correlationData;
                ArrayList<Double> poICCData = new ArrayList<Double>();
                ArrayList<String> poICCResponseIDs = new ArrayList<String>();
                ArrayList<Double> correl1 = new ArrayList<Double>();
                ArrayList<Double> correl2 = new ArrayList<Double>();
                for (int k = 0; k < poPairs.size(); k++) {
                  Integer resid1 = subIndexMap.get(poPairs.get(k)[0]);
                  Integer resid2 = subIndexMap.get(poPairs.get(k)[1]);
                  if (resid1 == null || resid2 == null) {
                    continue; // skip
                  }
                  cntPO++;
                  correl1.add(resids[resid1.intValue()]);
                  correl2.add(resids[resid2.intValue()]);
                  poICCData.add(resids[resid1.intValue()]);
                  poICCData.add(resids[resid2.intValue()]);
                  poICCResponseIDs.add(poPairs.get(k)[0] + "\t" + poPairs.get(k)[1]);
                  poICCResponseIDs.add(poPairs.get(k)[0] + "\t" + poPairs.get(k)[1]);
                }
                correlationData =
                                new double[][] {Doubles.toArray(correl1), Doubles.toArray(correl2)};
                poCorrel = Correlation.Pearson(correlationData);
                ICC poICCAnalysis = new ICC(Doubles.toArray(poICCData),
                                            Array.toStringArray(poICCResponseIDs), null, null, true,
                                            log);
                poICCAnalysis.computeICC();
                poICC = poICCAnalysis.getICC();
              }
              /* trios : */ {
                double[][] correlationData;
                ArrayList<Double> trioICCData = new ArrayList<Double>();
                ArrayList<String> trioICCResponseIDs = new ArrayList<String>();
                ArrayList<Double> correl1 = new ArrayList<Double>();
                ArrayList<Double> correl2 = new ArrayList<Double>();
                for (int k = 0; k < trios.size(); k++) {
                  int iidInd = trios.get(k)[0];
                  int faInd = trios.get(k)[1];
                  int moInd = trios.get(k)[2];
                  String fidiid = famIdHash.get(ped.getIID(iidInd)) + "\t" + ped.getIID(iidInd);
                  String fidiid2 = famIdHash.get(ped.getIID(faInd)) + "\t" + ped.getIID(faInd);
                  String fidiid3 = famIdHash.get(ped.getIID(moInd)) + "\t" + ped.getIID(moInd);
                  Integer resid1 = subIndexMap.get(fidiid);
                  Integer resid2 = subIndexMap.get(fidiid2);
                  Integer resid3 = subIndexMap.get(fidiid3);
                  if (resid1 == null || resid2 == null || resid3 == null) {
                    continue; // skip
                  }
                  cntTrio++;
                  correl1.add(resids[resid1]);
                  double p1 = resids[resid2];
                  double p2 = resids[resid3];
                  correl2.add((p1 + p2) / 2);
                  trioICCData.add(resids[resid1]);
                  trioICCData.add((p1 + p2) / 2);
                  trioICCResponseIDs.add(trios.get(k)[0] + "\t" + trios.get(k)[1] + "\t"
                                         + trios.get(k)[2]);
                  trioICCResponseIDs.add(trios.get(k)[0] + "\t" + trios.get(k)[1] + "\t"
                                         + trios.get(k)[2]);
                }
                correlationData =
                                new double[][] {Doubles.toArray(correl1), Doubles.toArray(correl2)};
                trioCorrel = Correlation.Pearson(correlationData);
                ICC trioICCAnalysis = new ICC(Doubles.toArray(trioICCData),
                                              Array.toStringArray(trioICCResponseIDs), null, null,
                                              true, log);
                trioICCAnalysis.computeICC();
                trioICC = trioICCAnalysis.getICC();
              }
            }
            log.report("Heritability for " + root);
            log.report(line[1] + " ~ " + Array.toStr(Array.subArray(line, 2), " "));

            merlinEstimate = computeWithMerlin(dir + root, pedigreeFile, "pheno.dat",
                                               line.length > 2 ? "covars.dat" : null, root,
                                               merlinExec, log);
            solarEstimate = computeWithSolar(dir + root, pedigreeFile, "pheno.dat",
                                             line.length > 2 ? "covars.dat" : null, root, solarExec,
                                             log);
            if (solarEstimate == null) {
              solarEstimate = new String[] {"NaN", "NaN", "NaN", "NaN", "NaN"};
            }
            numOfAllSamples = counter.getTotalCount();
            numOfFamiliesSizedOne = counter.getSizeOfCountEquals(1);
            numOfFamiliesSizedTwoOrAbove = counter.getSizeOfCountGreaterThan(2);
            log.report("Number of samples: " + numOfAllSamples + "\nNumber of families: "
                       + counter.getSize() + "\nNumber of families of size>=2: "
                       + numOfFamiliesSizedTwoOrAbove + "\nAverage size of families of size>=2: "
                       + ext.formDeci((numOfAllSamples - numOfFamiliesSizedOne)
                                      / (float) numOfFamiliesSizedTwoOrAbove, 3)
                       + "\nNumber of families of size=1: " + numOfFamiliesSizedOne);
            // summary.println(root + "\t" + merlinEstimate + "\t" + solarEstimate[0] + "\t" +
            // solarEstimate[1] + "\t" + numOfAllSamples + "\t" + counter.getSize() + "\t" +
            // numOfFamiliesSizedTwoOrAbove + "\t" + String.format("%.3", ((float) (numOfAllSamples
            // - numOfFamiliesSizedOne)) / numOfFamiliesSizedTwoOrAbove) + "\t" +
            // numOfFamiliesSizedOne);
            summary.println(root + "\t" + merlinEstimate + "\t" + Array.toStr(solarEstimate) + "\t"
                            + numOfAllSamples + "\t" + counter.getSize() + "\t"
                            + numOfFamiliesSizedTwoOrAbove + "\t"
                            + ext.formDeci((float) (numOfAllSamples - numOfFamiliesSizedOne)
                                           / numOfFamiliesSizedTwoOrAbove, 3)
                            + "\t" + numOfFamiliesSizedOne + "\t" + sibCorrel[0] + "\t"
                            + poCorrel[0] + "\t" + trioCorrel[0] + "\t" + sibCorrel[1] + "\t"
                            + poCorrel[1] + "\t" + trioCorrel[1] + "\t" + sibICC + "\t" + poICC
                            + "\t" + trioICC + "\t" + cntSib + "\t" + cntPO + "\t" + cntTrio);
            summary.flush();
            log.report("");
          }
        }
        summary.close();
      } catch (Exception e) {
        System.err.println("Error writing to " + ext.rootOf(filename, false) + "_summary.xln");
        e.printStackTrace();
      }


    }
  }



  public static class HeritabilityProducer extends AbstractProducer<String> {

    @Override
    public boolean hasNext() {
      // TODO Auto-generated method stub
      return false;
    }

    @Override
    public Callable<String> next() {
      // TODO Auto-generated method stub
      return null;
    }
  }

  /**
   * For internal generation of heritability .crfs
   *
   * @param pedfile ped file
   * @param db db file
   * @param crf full path to the desired crf
   * @param modelRoot will be concatenated with models[i]
   * @param models
   * @param log
   */
  public static void developCrf(String pedfile, String db, String crf, String modelRoot,
                                String[] models, Logger log) {
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(crf));
      writer.println("heritability");
      writer.println("ped=" + pedfile);
      writer.println("db=" + db);
      for (String model : models) {
        writer.println(modelRoot + "_" + model + "\t" + model);
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + crf);
      log.reportException(e);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String phenoDir = "N:/statgen/Mitochondrial_CN_heritability/phenos/";
    String pheno = "pheno.dat";
    String covars = "covars.dat";
    String pedfile = "plink.fam";
    // String method = "merlin";
    // String method = "solar";
    String method = "both";
    Logger log;
    String logfile = null;
    String controlFile = null;
    String prefix = "vc-chr01";

    // fromParameters("N:/statgen/BOSS/phenotypes/PhenoPrep/word_recogition/heritability_noCovars.crf",
    // new Logger());
    // System.exit(0);

    String usage = "\n" + "link.Heritability requires 0-1 arguments\n"
                   + "   (1) pedigree filename (i.e. ped=" + pedfile + " (default))\n"
                   + "   (2) phenotype filename (i.e. pheno=" + pheno + " (default))\n"
                   + "   (3) covariates filename (i.e. covars=" + covars
                   + " (default; set to null if not needed))\n"
                   + "   (4) method: merlin/solar/both (i.e. method=" + method + " (default))\n"
                   + "   (5) (optional) prefix for the files (i.e. prefix=" + prefix
                   + " (default))\n" + " OR:\n"
                   + "   (1) control file (i.e. crf=pheno.crf (not the default))\n" + "";

    controlFile = null;
    method = "";
    phenoDir = null;

    phenoDir = "N:/statgen/Mitochondrial_CN_heritability/phenos/";
    pedfile = "N:/statgen/Mitochondrial_CN_heritability/pedigrees/pedigree.dat";
    covars = null;

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("ped=")) {
        pedfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("pheno=")) {
        pheno = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("covars=")) {
        covars = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("method=")) {
        method = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("prefix=")) {
        method = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("crf=")) {
        controlFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("phenodir=")) {
        phenoDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    // if (args.length == 0) {
    // controlFile = "heritabilityTest.crf";
    // controlFile = "N:/statgen/BOSS/phenotypes/PhenoPrep/MMSE/heritability.crf";
    // }
    try {
      log = new Logger(logfile);
      if (controlFile != null) {
        fromParameters(controlFile, false, new Logger(ext.rootOf(controlFile, false) + ".log"));
      } else {
        if (phenoDir != null) {
          getHeritabilitiesOfAllPhenosInADir(phenoDir, pedfile, covars, log);
        } else {
          if (method.equalsIgnoreCase("merlin") || method.equalsIgnoreCase("both")) {
            computeWithMerlin(pedfile, pheno, covars, prefix, "", DEFAULT_MERLIN_EXEC, log);
          }
          if (method.equalsIgnoreCase("solar") || method.equalsIgnoreCase("both")) {
            computeWithSolar(pedfile, pheno, covars, prefix, "", DEFAULT_SOLAR_EXEC, log);
          }
        }
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
