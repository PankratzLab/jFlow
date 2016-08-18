package org.genvisis.gwas;

import com.google.common.primitives.Ints;

import java.io.File;
import java.util.Vector;

import org.genvisis.common.Aliases;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.filesys.MetaAnalysisParams;
import org.genvisis.stats.Rscript;

public class SkatMeta {
  public static final String[] ALGORITHMS =
      {"singlesnpMeta", "burdenMeta", "skatMeta", "skatOMeta"};

  public static final String[][] UNIT_OF_ANALYSIS =
      {Aliases.MARKER_NAMES, Aliases.GENE_UNITS, Aliases.GENE_UNITS, Aliases.GENE_UNITS};

  public static final boolean[] SINGLE_VARIANTS = {true, false, false, false};

  public static final String[][] HEADER_TYPES =
      {{"gene", "Name", "p", "maf", "nmiss", "ntotal", "beta", "se"}, // Single SNP
       {"gene", "p", "beta", "se", "cmafTotal", "cmafUsed", "nsnpsTotal", "nsnpsUsed", "nmiss"}, // Burden
                                                                                                 // Test
       {"gene", "p", "Qmeta", "cmaf", "nmiss", "nsnps"}, // SKAT test
       {"gene", "p", "Qmeta", "cmaf", "nmiss", "nsnps"} // SKAT-O test (not verified)
      };

  public static String getRscriptExecutable(MetaAnalysisParams maps, Logger log) {
    if (maps != null && maps.getRExec() != null) {
      return maps.getRExec();
    } else {
      return Rscript.getRscriptExecutable(log);
    }
  }

  private static void determineObjectNames(String dir, MetaAnalysisParams maps, Logger log) {
    String[] lines, files;
    String[][] iterations;
    String root, commands, filename;
    Vector<String> v, remaining;

    files = Files.list("./", null, ".rdata", false, false);
    log.report("There are " + files.length + " total .Rdata files");

    dir = new File(dir).getAbsolutePath() + "/";

    v = new Vector<String>();
    remaining = new Vector<String>();
    new File("batchChecks/").mkdir();
    for (String file : files) {
      root = ext.rootOf(file);
      if (!Files.exists("batchChecks/" + root + ".object")) {
        lines = new String[] {"load(\"" + file + "\")", "name <- ls()",
                              // "write.table(name, \"name.txt\", sep=\"\t\")",
                              "fileConn<-file(\"batchChecks/" + root + ".object\")",
                              "writeLines(c(name), fileConn)", "close(fileConn)",};
        filename = "batchChecks/" + root + ".R";
        Files.writeList(lines, filename);
        v.add(filename);
        remaining.add(file);
      }
    }
    log.report("There are " + v.size() + " .Rdata files remaining to interrogate:\n"
               + Array.toStr(Array.toStringArray(v), "\n"));

    if (v.size() > 0) {
      commands = getRscriptExecutable(maps, log) + " --no-save [%0]";
      iterations = Matrix.toMatrix(Array.toStringArray(v));
      Files.qsub("batchChecks/checkObject", dir, -1, commands, iterations, 4000, 1);
      Files.batchIt("master.checkObjectAll", null, 1, commands, iterations);
    }
  }

  public static String getObjectName(String dir, String filename) {
    return HashVec.loadFileToStringArray(dir + "batchChecks/" + ext.rootOf(filename) + ".object",
                                         false, new int[] {0}, false)[0];
  }

  public static String[][][] identifySet(MetaAnalysisParams maps, String[] files, Logger log) {
    String[][][] finalSets;
    boolean[] picks, used;
    int numMatches, index;
    String[][] phenotypes, maskedPhenos;
    String[] studies;
    String[][] races;

    used = Array.booleanArray(files.length, false);

    index = ext.indexOfStr(maps.getSnpInfoFilename(), files);
    if (index >= 0) {
      used[index] = true;
    }

    phenotypes = maps.getPhenotypesWithFilenameAliases(true);
    maskedPhenos = maps.getPhenotypesWithFilenameAliases(false);
    studies = maps.getStudies();
    races = maps.getRacesWithFilenameAliases();

    finalSets = new String[phenotypes.length][][]; // [pheno][study][race] <- all files meeting
                                                   // criteria
    for (int i = 0; i < phenotypes.length; i++) {
      finalSets[i] = Matrix.stringMatrix(studies.length, races.length, "<missing>");
      log.report("For " + phenotypes[i][0] + " identified:", true, true);
      log.report("\tStudy\t" + Array.toStr(Matrix.extractColumn(races, 0)));
      for (int j = 0; j < studies.length; j++) {
        log.report("\t" + studies[j], false, true);
        for (int k = 0; k < races.length; k++) {
          picks = Array.booleanArray(files.length, false);
          for (int f = 0; f < files.length; f++) {
            if (files[f].contains(studies[j]) && ext.containsAny(files[f], maskedPhenos[i])
                && ext.containsAny(files[f], races[k])) {
              picks[f] = true;
              if (finalSets[i][j][k].equals("<missing>")) {
                finalSets[i][j][k] = files[f];
              } else {
                finalSets[i][j][k] += ";" + files[f];
              }
              if (used[f]) {
                log.reportError("Error - file '" + files[f] + "' matches to " + studies[j] + "/"
                                + phenotypes[i][0] + " but was already picked for another purpose");
              }
              used[f] = true;
            }
          }
          numMatches = Array.booleanArraySum(picks);
          if (numMatches == 0) {
            // log.reportError("Warning - could not find a match for
            // "+studies[j]+"/"+phenotypes[i][0]+"/"+races[k][0]);
          } else if (numMatches > 1) {
            // log.reportError("Error - found multiple matched for
            // "+studies[j]+"/"+phenotypes[i][0]);
            // log.reportError(Array.toStr(Array.subArray(files, picks), "\n"));
          }
          // log.report(" "+finalSets[i][j], true, false);
          log.report("\t" + numMatches, false, true);
        }
        log.report("");
      }
      log.report("");
    }

    numMatches = Array.booleanArraySum(used);
    if (numMatches != files.length) {
      log.reportError("Warning - did not find a match for the following file(s):");
      for (int i = 0; i < files.length; i++) {
        if (!used[i]) {
          log.reportError("  " + files[i]);
        }
      }
    }

    return finalSets;
  }

  public static int getMaxChr() {
    String chrom;
    int maxChr;

    maxChr = 0;
    for (int chr = 1; chr <= 24; chr++) {
      chrom = chr == 23 ? "X" : (chr == 24 ? "Y" : chr + "");
      if (Files.exists("snpInfos/snpInfo_chr" + chrom + ".RData")) {
        maxChr = chr;
      }
    }

    return maxChr;
  }

  public static void splitAll(String dir, MetaAnalysisParams maps) {
    String[] files;
    String[][][] finalSets;
    Logger log;
    String localDir;
    String filename;
    Vector<String> toBeSplit, commands;
    String objectName, snpInfoName, chrom, subsetObject;
    String[][] phenotypes;
    String[] studies;
    String[][] races;
    String snpInfoFile;
    String chromName;
    String geneName;
    boolean problem;
    int maxChr;
    IntVector chrsToDo;
    IntVector jobSizes;
    Vector<String> jobNames;

    problem = false;
    if (dir == null || dir.equals("")) {
      dir = new File("").getAbsolutePath() + "/";
    }

    log = new Logger(dir + "splitAll.log");
    new File(dir + "batchSplits/").mkdir();
    files = Files.list(dir, null, ".Rdata", false, false);

    phenotypes = maps.getPhenotypesWithFilenameAliases();
    studies = maps.getStudies();
    races = maps.getRacesWithFilenameAliases();
    snpInfoFile = maps.getSnpInfoFilename();
    chromName = maps.getChromName();
    geneName = maps.getGeneName();

    if (ext.indexOfStr(snpInfoFile, files) == -1) {
      log.reportError("Error - could not find SNP Info file '" + snpInfoFile + "'; aborting");
      return;
    }

    if (Files.exists(dir + "batchChecks/" + ext.rootOf(snpInfoFile) + ".object")) {
      snpInfoName = getObjectName(dir, snpInfoFile);
    } else {
      log.reportError("Error - could not find file '" + dir + "batchChecks/"
                      + ext.rootOf(snpInfoFile) + ".object" + "'");
      snpInfoName = "UNKNOWN_SNP_INFO_OBJECT_NAME";
      problem = true;
    }


    commands = new Vector<String>();
    commands.add("load(\"" + dir + snpInfoFile + "\")");
    commands.add("ls()");

    commands.add("chroms <- unique(" + snpInfoName + "$" + chromName + ")");
    commands.add("write.table( chroms, \"chroms.csv\", sep=\",\")");

    commands.add("for (chr in chroms) {");
    commands.add("  snps_on_chr <- " + snpInfoName + "[" + snpInfoName + "$CHROM == chr,]");
    commands.add("  filename <- paste(\"snpInfos/snpInfo_chr\", chr, \".RData\", sep='')");
    commands.add("  save(snps_on_chr, file=filename, compress=\"bzip2\")");
    commands.add("}");

    filename = dir + "batchSplits/splitChrs.R";
    Files.writeList(Array.toStringArray(commands), filename);

    new File(dir + "snpInfos/").mkdirs();
    Files.qsub(dir + "batchSplits/" + ext.rootOf(filename) + ".qsub",
               "cd " + dir + "\n" + getRscriptExecutable(maps, log) + " --no-save " + filename,
               5000, 0.25, 1);

    toBeSplit = new Vector<String>();
    toBeSplit.add("# make sure to run \"qsub " + ext.rootOf(filename) + ".qsub\" first!!!");
    toBeSplit.add("cd batchSplits/");

    jobNames = new Vector<String>();
    jobSizes = new IntVector();

    dir = ext.verifyDirFormat(dir);
    finalSets = identifySet(maps, files, log);
    for (int i = 0; i < phenotypes.length; i++) {
      for (int j = 0; j < studies.length; j++) {
        for (int k = 0; k < races.length; k++) {
          if (!finalSets[i][j][k].equals("<missing>")) {
            localDir =
                dir + "objects/" + studies[j] + "/" + races[k][0] + "/" + phenotypes[i][0] + "/";
            new File(localDir).mkdirs();

            files = finalSets[i][j][k].split(";");
            for (int f = 0; f < files.length; f++) {
              commands = new Vector<String>();
              commands.add("load(\"" + dir + snpInfoFile + "\")");
              commands.add("load(\"" + dir + files[f] + "\")");
              if (Files.exists(dir + "batchChecks/" + ext.rootOf(files[f]) + ".object")) {
                objectName = getObjectName(dir, files[f]);
              } else {
                log.reportError("Error - could not find file '" + dir + "batchChecks/"
                                + ext.rootOf(files[f]) + ".object" + "'");
                objectName = "UNKNOWN_SKAT_COHORT_OBJECT_NAME";
                problem = true;
              }
              commands.add("ls()");

              chrsToDo = new IntVector();
              maxChr = getMaxChr();
              for (int chr = 1; chr <= maxChr; chr++) {
                chrom = chr == 23 ? "X" : (chr == 24 ? "Y" : chr + "");
                subsetObject =
                    studies[j] + "_" + races[k][0] + "_" + phenotypes[i][0] + "_chr" + chrom;
                if (!Files.exists(localDir + subsetObject + "_f" + f + ".RData")
                    && !Files.exists(localDir + subsetObject + ".RData")) {
                  chrsToDo.add(chr);
                }
              }
              if (chrsToDo.size() != 0 && chrsToDo.size() != maxChr) {
                log.reportError("Warning - for " + studies[j] + ";" + races[k][0] + "/"
                                + phenotypes[i][0] + ", missing chr(s) "
                                + ext.listWithCommas(Array.toStringArray(Ints.toArray(chrsToDo))));
                log.reportError("        - if batch job was killed in the middle, suggest deleting the last attempted chromosome, in case it was incomplete");
              }

              for (int c = 0; c < chrsToDo.size(); c++) {
                chrom = chrsToDo.elementAt(c) == 23 ? "X"
                                                    : (chrsToDo.elementAt(c) == 24 ? "Y"
                                                                                   : chrsToDo.elementAt(c)
                                                                                     + "");
                subsetObject =
                    studies[j] + "_" + races[k][0] + "_" + phenotypes[i][0] + "_chr" + chrom;

                // filter for the gene names present on the chromosome
                commands.add("genes <- unique(" + snpInfoName + "[" + snpInfoName + "$" + chromName
                             + " == \"" + chrom + "\", \"" + geneName + "\"])");

                // take the intersect of those actually present in the skatCohort object
                commands.add("idx <- intersect(genes, names(" + objectName + "))");

                // create the skatCohort subset
                commands.add(subsetObject + " <- " + objectName + "[idx]");

                // make sure the dataset has the skatCohort class
                commands.add("class(" + subsetObject + ") <- \"skatCohort\"");

                // save the new file
                commands.add("save(" + subsetObject + ", file=\"" + localDir + subsetObject + "_f"
                             + f + ".RData\", compress=\"bzip2\")");

                // free up memory
                commands.add("rm(" + subsetObject + ")");
                commands.add("");
              }

              if (chrsToDo.size() > 0) {
                filename = dir + "batchSplits/" + studies[j] + "_" + races[k][0] + "_"
                           + phenotypes[i][0] + "_f" + f + ".R";
                Files.writeList(Array.toStringArray(commands), filename);

                Files.qsub(dir + "batchSplits/" + ext.rootOf(filename) + ".qsub",
                           "cd " + dir + "\n" + getRscriptExecutable(maps, log) + " --no-save "
                                                                                  + filename,
                           10000, 0.5, 1);
                toBeSplit.add("qsub " + ext.rootOf(filename) + ".qsub");
                jobNames.add(dir + "batchSplits/" + ext.rootOf(filename) + ".qsub");
                jobSizes.add((int) (new File(dir + files[f]).length()
                                    + chrsToDo.size() * 2 / maxChr
                                      * new File(dir + files[f]).length()));
              }
            }
          }
        }
      }
    }

    if (problem) {
      log.reportError("   need to first run using the -determineObjectNames option");
      return;
    }

    toBeSplit.add("# make sure to run \"SkatMeta -consolidate\" after everything else is run!!!");
    Files.writeList(Array.toStringArray(toBeSplit), dir + "master.toBeSplit");
    Files.chmod(dir + "master.toBeSplit");

    log.report("");
    log.report("Make sure to run \"qsub splitChrs.qsub\" first!!!");

    Files.qsubMultiple(jobNames, jobSizes, "chunks/", "chunkSplit", 8, true, null, -1, 22000, 2);
  }

  public static void consolidate(String dir, MetaAnalysisParams maps) {
    String[] files;
    String[][][] finalSets;
    Logger log;
    String localDir;
    String filename;
    String chrom, subsetObject;
    String[][] phenotypes;
    String[] studies;
    String[][] races;
    boolean problem;
    long fileSize, largestFileSize;
    int[][][][] finalSelections;
    int count;
    int maxChr;

    problem = false;
    if (dir == null || dir.equals("")) {
      dir = new File("").getAbsolutePath() + "/";
    }

    log = new Logger(dir + "consolidateAll.log");
    files = Files.list(dir, null, ".Rdata", false, false);

    phenotypes = maps.getPhenotypesWithFilenameAliases();
    studies = maps.getStudies();
    races = maps.getRacesWithFilenameAliases();
    maxChr = getMaxChr();

    count = 0;
    dir = ext.verifyDirFormat(dir);
    finalSets = identifySet(maps, files, log);
    finalSelections = new int[finalSets.length][finalSets[0].length][finalSets[0][0].length][];
    for (int iter = 0; iter < 2; iter++) {
      for (int i = 0; i < phenotypes.length; i++) {
        System.out.println(phenotypes[i][0]);
        for (int j = 0; j < studies.length; j++) {
          System.out.println("  " + studies[j]);
          for (int k = 0; k < races.length; k++) {
            if (!finalSets[i][j][k].equals("<missing>")) {
              System.out.println("    " + races[k][0]);
              localDir =
                  dir + "objects/" + studies[j] + "/" + races[k][0] + "/" + phenotypes[i][0] + "/";
              files = finalSets[i][j][k].split(";");
              if (iter == 0) {
                finalSelections[i][j][k] = Array.intArray(maxChr, -1);
              }

              for (int chr = 1; chr <= maxChr; chr++) {

                chrom = chr == 23 ? "X" : (chr == 24 ? "Y" : chr + "");
                subsetObject =
                    studies[j] + "_" + races[k][0] + "_" + phenotypes[i][0] + "_chr" + chrom;

                largestFileSize = 0;
                for (int f = 0; f < files.length; f++) {
                  filename = subsetObject + "_f" + f + ".RData";

                  if (iter == 0) {
                    if (Files.exists(localDir + filename)) {
                      fileSize = new File(localDir + filename).length();
                      if (fileSize > largestFileSize) {
                        largestFileSize = fileSize;
                        finalSelections[i][j][k][chr - 1] = f;
                      }
                      count++;
                      // } else {
                    } else if (!Files.exists(localDir + subsetObject + ".RData")) {
                      System.err.println("Error - could not find '" + subsetObject + "_f" + f
                                         + ".RData' in " + localDir);
                      problem = true;
                    }
                  } else {
                    if (f == finalSelections[i][j][k][chr - 1]) {
                      new File(localDir + filename).renameTo(new File(localDir + subsetObject
                                                                      + ".RData"));
                    } else {
                      new File(localDir + filename).delete();
                    }
                  }
                }


              }
            }
          }
        }
      }
      if (iter == 0) {
        if (problem) {
          if (count == 0) {
            log.reportError("\n   discrepancies found; possible explanations are that the R parser was not run or did not complete, or the original .RData files were either removed or added to; if the latter, then suggest rerunning the parsers for those study/pheno pairs; no consolidation will occur\n");
          } else {
            log.reportError("\n   did not find a single file with a _f# extension; either nothing has been run or everything has already been processed by this algorithm; check the objects/ directory\n");
          }
          return;
        } else {
          log.reportError("\nEverything seems to be in order, proceeding with consolidation\n");
        }
      }
    }
  }

  public static void runAll(String dir, MetaAnalysisParams maps, boolean forceMeta) {
    String[] files;
    String[][][] finalSets;
    Logger log;
    String localDir;
    String root, filename, objectFilename, outputFilename;
    Vector<String> toBeRunIndividually, toBeRunMetad, commands, objects;
    int count;
    String originalObjectName, objectName, snpInfoFile, snpInfoName, chrom;
    String[][] phenotypes, races;
    String[] studies;
    String snpName;
    String[][] methods;
    String functionFlagName, geneName;
    boolean runningByChr;
    IntVector jobSizes;
    Vector<String> jobNames;
    int[] infoSizes;
    int maxChr;

    if (dir == null || dir.equals("")) {
      dir = new File("").getAbsolutePath() + "/";
    }

    log = new Logger(dir + "runAll.log");
    phenotypes = maps.getPhenotypesWithFilenameAliases();
    studies = maps.getStudies();
    races = maps.getRacesWithFilenameAliases();
    snpName = maps.getVariantName();
    methods = maps.getMethods();
    functionFlagName = maps.getFunctionFlagName();
    geneName = maps.getGeneName();
    runningByChr = maps.runningByChr();
    snpInfoFile = maps.getSnpInfoFilename();

    files = Files.list(dir, null, ".Rdata", false, false);
    finalSets = identifySet(maps, files, log);

    maxChr = getMaxChr();
    jobSizes = new IntVector();
    jobNames = new Vector<String>();
    infoSizes = new int[maxChr + 2];

    if (runningByChr) {
      for (int chr = 1; chr <= maxChr; chr++) {
        chrom = chr == 23 ? "X" : (chr == 24 ? "Y" : chr + "");
        filename = "snpInfos/snpInfo_chr" + chrom + ".RData";
        if (!Files.exists(filename)) {
          log.reportError("Error - could not find SNP Info file '" + filename + "'; aborting");
          return;
        } else {
          infoSizes[chr] = (int) new File(filename).length();
        }

      }
      snpInfoName = "snps_on_chr";
    } else {
      if (!Files.exists(snpInfoFile)) {
        log.reportError("Error - could not find SNP Info file '" + snpInfoFile + "'; aborting");
        return;
      }
      if (Files.exists(dir + "batchChecks/" + ext.rootOf(snpInfoFile) + ".object")) {
        snpInfoName = getObjectName(dir, snpInfoFile);
      } else {
        log.reportError("Error - could not find file '" + dir + "batchChecks/"
                        + ext.rootOf(snpInfoFile) + ".object" + "'");
        return;
      }
    }

    toBeRunIndividually = new Vector<String>();
    toBeRunIndividually.add("cd batchRuns/");
    toBeRunMetad = new Vector<String>();
    toBeRunMetad.add("cd batchRuns/");
    new File(dir + "batchRuns/").mkdir();
    dir = ext.verifyDirFormat(dir);
    for (String[] phenotype : phenotypes) {
      for (String[] method : methods) {
        localDir = dir + phenotype[0] + "/" + method[0] + "/";
        new File(localDir).mkdirs();
        for (String[] race : races) {
          localDir = dir + phenotype[0] + "/" + race[0] + "/" + method[0] + "/";
          new File(localDir).mkdirs();
        }
      }
    }

    // Primary analysis by study/race
    for (int i = 0; i < phenotypes.length; i++) {
      for (int j = 0; j < studies.length; j++) {
        for (int k = 0; k < races.length; k++) {
          if (!finalSets[i][j][k].equals("<missing>")) {
            localDir =
                dir + "objects/" + studies[j] + "/" + races[k][0] + "/" + phenotypes[i][0] + "/";

            for (int chr = 1; chr <= (runningByChr ? maxChr : 1); chr++) {
              chrom = chr == 23 ? "X" : (chr == 24 ? "Y" : chr + "");

              if (runningByChr) {
                objectName =
                    studies[j] + "_" + races[k][0] + "_" + phenotypes[i][0] + "_chr" + chrom;
                objectFilename = localDir + objectName + ".RData";
                snpInfoFile = "snpInfos/snpInfo_chr" + chrom + ".RData";
              } else {
                objectFilename = dir + finalSets[i][j][k];
                if (objectFilename.contains(";")) {
                  log.reportError("Error - more than one file is mapped to " + studies[j] + "_"
                                  + races[k][0] + "_" + phenotypes[i][0] + ": " + objectFilename);
                  return;
                } else if (Files.exists(dir + "batchChecks/" + ext.rootOf(objectFilename)
                                        + ".object")) {
                  objectName = getObjectName(dir, objectFilename);
                } else {
                  log.reportError("Error - could not find file '" + dir + "batchChecks/"
                                  + ext.rootOf(objectFilename) + ".object" + "'");
                  objectName = "UNKNOWN_SKAT_COHORT_OBJECT_NAME";
                }
              }

              if (!Files.exists(objectFilename)) {
                log.reportError("Error - missing object file: '" + objectFilename + "'");
                if (Files.exists(ext.addToRoot(objectFilename, "_f0"))) {
                  log.reportError("     - however did find '"
                                  + ext.removeDirectoryInfo(ext.addToRoot(objectFilename, "_f0"))
                                  + "'; so try running -consolidate");
                  return;
                }
              }

              commands = new Vector<String>();
              commands.add("library(skatMeta)");
              commands.add("load(\"" + dir + snpInfoFile + "\")");
              commands.add("load(\"" + objectFilename + "\")");
              commands.add("ls()");
              count = 0;
              for (String[] method : methods) {
                root = studies[j] + "_" + races[k][0] + "_" + phenotypes[i][0] + "_" + method[0];
                outputFilename = dir + phenotypes[i][0] + "/" + races[k][0] + "/" + method[0] + "/"
                                 + root + (runningByChr ? "_chr" + chrom : "") + ".csv";
                if (!Files.exists(outputFilename) || new File(outputFilename).length() == 0) {
                  if (new File(objectFilename).length() > 1024) {
                    commands.add("results <- " + method[2] + "(" + objectName + ", SNPInfo="
                                 + (SINGLE_VARIANTS[ext.indexOfStr(method[2], ALGORITHMS)]
                                    || functionFlagName == null ? snpInfoName
                                                                : "subset(" + snpInfoName + ", "
                                                                  + functionFlagName + "==TRUE)")
                                 + ", snpNames = \"" + snpName + "\"" + ", aggregateBy=\""
                                 + geneName + "\""
                                 + (method.length > 3
                                    && ext.isValidDouble(method[3]) ? ", mafRange = c(0,"
                                                                      + method[3] + ")"
                                                                      + (method.length > 4 ? ", "
                                                                                             + Array.toStr(Array.subArray(method,
                                                                                                                          4),
                                                                                                           ", ")
                                                                                           : "")
                                                                    : (method.length > 3 ? ", "
                                                                                           + Array.toStr(Array.subArray(method,
                                                                                                                        3),
                                                                                                         ", ")
                                                                                         : ""))
                                 + ")");
                    commands.add("write.table( results, \"" + outputFilename
                                 + "\", sep=\",\", row.names = F)");
                    count++;
                  } else {
                    Files.write(Array.toStr(getHeaderForMethod(method), ","), outputFilename);
                  }
                }
              }
              if (count > 0) {
                count = 0;
                do {
                  filename = dir + "batchRuns/" + studies[j] + "_" + races[k][0] + "_"
                             + phenotypes[i][0] + (runningByChr ? "_chr" + chrom : "")
                             + (count == 0 ? "" : "_" + count) + ".R";
                  count++;
                } while (Files.exists(filename));
                Files.writeList(Array.toStringArray(commands), filename);

                Files.qsub(dir + "batchRuns/" + ext.rootOf(filename) + ".qsub",
                           "cd " + dir + "\n" + getRscriptExecutable(maps, log) + " --no-save "
                                                                                + filename,
                           5000, 1, 1);
                toBeRunIndividually.add("qsub " + ext.rootOf(filename) + ".qsub");
                jobNames.add(dir + "batchRuns/" + ext.rootOf(filename) + ".qsub");
                jobSizes.add(infoSizes[chr]);
              }
            }
          }
        }
      }
    }

    Files.writeList(Array.toStringArray(toBeRunIndividually), dir + "master.toBeRunIndividually");
    Files.chmod(dir + "master.toBeRunIndividually");
    System.err.println("qsubing multiple individual runs");
    Files.qsubMultiple(jobNames, jobSizes, "chunks/", "chunkRun", 16, true, "sb", -1, 62000, 2);
    System.err.println("multiple individual runs done");


    jobNames = new Vector<String>();
    jobSizes = new IntVector();
    // Meta-analysis stratified by race
    for (int i = 0; i < phenotypes.length; i++) {
      for (int k = 0; k < races.length; k++) {
        for (int chr = 1; chr <= (runningByChr ? maxChr : 1); chr++) {
          chrom = chr == 23 ? "X" : (chr == 24 ? "Y" : chr + "");
          commands = new Vector<String>();
          commands.add("library(skatMeta)");
          if (runningByChr) {
            snpInfoFile = "snpInfos/snpInfo_chr" + chrom + ".RData";
          }
          commands.add("load(\"" + dir + snpInfoFile + "\")");

          objects = new Vector<String>();
          for (int j = 0; j < studies.length; j++) {
            if (!finalSets[i][j][k].equals("<missing>")) {
              if (runningByChr) {
                localDir = dir + "objects/" + studies[j] + "/" + races[k][0] + "/"
                           + phenotypes[i][0] + "/";
                objectName =
                    studies[j] + "_" + races[k][0] + "_" + phenotypes[i][0] + "_chr" + chrom;
                if (new File(localDir + objectName + ".RData").length() > 1024) {
                  commands.add("load(\"" + localDir + objectName + ".RData" + "\")");
                  objects.add(objectName);
                }
              } else {
                objectFilename = finalSets[i][j][k];
                originalObjectName = getObjectName(dir, objectFilename);
                objectName = studies[j] + "_" + races[k][0] + "_" + phenotypes[i][0];
                if (new File(dir + objectFilename).length() > 1024) {
                  commands.add("load(\"" + dir + objectFilename + "\")");
                  commands.add(objectName + " <- " + originalObjectName);
                  commands.add("rm(\"" + originalObjectName + "\")");
                  objects.add(objectName);
                }
              }
            }
          }
          commands.add("ls()");
          commands.add("");
          count = 0;
          for (String[] method : methods) {
            root = races[k][0] + "_" + phenotypes[i][0] + "_" + method[0];
            outputFilename = dir + phenotypes[i][0] + "/" + races[k][0] + "/" + method[0] + "/"
                             + root + (runningByChr ? "_chr" + chrom : "") + ".csv";
            if (forceMeta || !Files.exists(outputFilename)
                || new File(outputFilename).length() == 0) {
              if (objects.size() > 0) {
                commands.add("results <- " + method[2] + "("
                             + Array.toStr(Array.toStringArray(objects), ", ") + ", SNPInfo="
                             + (SINGLE_VARIANTS[ext.indexOfStr(method[2], ALGORITHMS)]
                                || functionFlagName == null ? snpInfoName
                                                            : "subset(" + snpInfoName + ", "
                                                              + functionFlagName + "==TRUE)")
                             + ", snpNames = \"" + snpName + "\"" + ", aggregateBy=\"" + geneName
                             + "\""
                             + (method.length > 3
                                && ext.isValidDouble(method[3]) ? ", mafRange = c(0," + method[3]
                                                                  + ")"
                                                                  + (method.length > 4 ? ", "
                                                                                         + Array.toStr(Array.subArray(method,
                                                                                                                      4),
                                                                                                       ", ")
                                                                                       : "")
                                                                : (method.length > 3 ? ", "
                                                                                       + Array.toStr(Array.subArray(method,
                                                                                                                    3),
                                                                                                     ", ")
                                                                                     : ""))
                             + ")");
                commands.add("write.table( results, \"" + outputFilename
                             + "\", sep=\",\", row.names = F)");
                commands.add("");
                count++;
              } else {
                Files.write(Array.toStr(getHeaderForMethod(method), ","), outputFilename);
              }

            }
          }
          if (count > 0) {
            count = 0;
            do {
              filename = dir + "batchRuns/" + races[k][0] + "_" + phenotypes[i][0]
                         + (runningByChr ? "_chr" + chrom : "") + (count == 0 ? "" : "_" + count)
                         + ".R";
              count++;
            } while (Files.exists(filename));
            Files.writeList(Array.toStringArray(commands), filename);

            Files.qsub(dir + "batchRuns/" + ext.rootOf(filename) + ".qsub",
                       "cd " + dir + "\n" + getRscriptExecutable(maps, log) + " --no-save "
                                                                            + filename,
                       25000, 2, 1);
            toBeRunMetad.add("qsub " + ext.rootOf(filename) + ".qsub");
            jobNames.add(dir + "batchRuns/" + ext.rootOf(filename) + ".qsub");
            jobSizes.add(infoSizes[chr]);
          }
        }
      }

      // Meta-analysis of all races
      for (int chr = 1; chr <= (runningByChr ? maxChr : 1); chr++) {
        chrom = chr == 23 ? "X" : (chr == 24 ? "Y" : chr + "");
        commands = new Vector<String>();
        commands.add("library(skatMeta)");
        if (runningByChr) {
          snpInfoFile = "snpInfos/snpInfo_chr" + chrom + ".RData";
        }
        commands.add("load(\"" + dir + snpInfoFile + "\")");

        objects = new Vector<String>();
        for (int j = 0; j < studies.length; j++) {
          for (int k = 0; k < races.length; k++) {
            if (!finalSets[i][j][k].equals("<missing>")) {
              if (runningByChr) {
                localDir = dir + "objects/" + studies[j] + "/" + races[k][0] + "/"
                           + phenotypes[i][0] + "/";
                objectName =
                    studies[j] + "_" + races[k][0] + "_" + phenotypes[i][0] + "_chr" + chrom;
                if (new File(localDir + objectName + ".RData").length() > 1024) {
                  commands.add("load(\"" + localDir + objectName + ".RData" + "\")");
                  objects.add(objectName);
                }
              } else {
                objectFilename = finalSets[i][j][k];
                originalObjectName = getObjectName(dir, objectFilename);
                objectName = studies[j] + "_" + races[k][0] + "_" + phenotypes[i][0];
                if (new File(dir + objectFilename).length() > 1024) {
                  commands.add("load(\"" + dir + objectFilename + "\")");
                  commands.add(objectName + " <- " + originalObjectName);
                  commands.add("rm(\"" + originalObjectName + "\")");
                  objects.add(objectName);
                }
              }
            }
          }
        }
        commands.add("ls()");
        commands.add("");
        count = 0;
        for (String[] method : methods) {
          root = phenotypes[i][0] + "_" + method[0];
          outputFilename = dir + phenotypes[i][0] + "/" + method[0] + "/" + root
                           + (runningByChr ? "_chr" + chrom : "") + ".csv";
          if (forceMeta || !Files.exists(outputFilename)
              || new File(outputFilename).length() == 0) {
            if (objects.size() > 0) {
              commands.add("results <- " + method[2] + "("
                           + Array.toStr(Array.toStringArray(objects), ", ") + ", SNPInfo="
                           + (SINGLE_VARIANTS[ext.indexOfStr(method[2], ALGORITHMS)]
                              || functionFlagName == null ? snpInfoName
                                                          : "subset(" + snpInfoName + ", "
                                                            + functionFlagName + "==TRUE)")
                           + ", snpNames = \"" + snpName + "\"" + ", aggregateBy=\"" + geneName
                           + "\""
                           + (method.length > 3
                              && ext.isValidDouble(method[3]) ? ", mafRange = c(0," + method[3]
                                                                + ")"
                                                                + (method.length > 4 ? ", "
                                                                                       + Array.toStr(Array.subArray(method,
                                                                                                                    4),
                                                                                                     ", ")
                                                                                     : "")
                                                              : (method.length > 3 ? ", "
                                                                                     + Array.toStr(Array.subArray(method,
                                                                                                                  3),
                                                                                                   ", ")
                                                                                   : ""))
                           + ")");
              commands.add("write.table( results, \"" + outputFilename
                           + "\", sep=\",\", row.names = F)");
              commands.add("");
              count++;
            } else {
              Files.write(Array.toStr(getHeaderForMethod(method), ","), outputFilename);
            }
          }
        }
        if (count > 0) {
          count = 0;
          do {
            filename = dir + "batchRuns/" + phenotypes[i][0] + (runningByChr ? "_chr" + chrom : "")
                       + (count == 0 ? "" : "_" + count) + ".R";
            count++;
          } while (Files.exists(filename));
          Files.writeList(Array.toStringArray(commands), filename);

          Files.qsub(dir + "batchRuns/" + ext.rootOf(filename) + ".qsub",
                     "cd " + dir + "\n" + getRscriptExecutable(maps, log) + " --no-save "
                                                                          + filename,
                     30000, 2, 1);
          toBeRunMetad.add("qsub " + ext.rootOf(filename) + ".qsub");
          jobNames.add(dir + "batchRuns/" + ext.rootOf(filename) + ".qsub");
          jobSizes.add(infoSizes[chr]);
        }
      }
    }
    Files.writeList(Array.toStringArray(toBeRunMetad), dir + "master.toBeMetaAnalyzed");
    Files.chmod(dir + "master.toBeMetaAnalyzed");
    System.err.println("qsubing multiple meta runs");
    Files.qsubMultiple(jobNames, jobSizes, "chunks/", "chunkMeta", 16, true, "sb", -1, 62000, 2);
    System.err.println("multiple meta runs done");
  }

  public static String[] getHeaderForMethod(String[] method) {
    return HEADER_TYPES[ext.indexOfStr(method[2], ALGORITHMS)];
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String logfile = null;
    Logger log;
    String dir = "";
    boolean determineObjectNames = false;
    boolean splitAll = false;
    boolean runAll = false;
    String mapsFile = "metaAnalysis.params";
    MetaAnalysisParams maps;
    boolean consolidate = false;
    boolean forceMeta = false;

    String usage = "\n" + "gwas.SkatMeta requires 0-1 arguments\n" + "   (0) directory (i.e. dir="
                   + dir + " (default))\n" + "   (1) filename of MetaAnalysisParameters (i.e. maps="
                   + mapsFile
                   + " (default; create an empty file of this name to populate with examples))\n"
                   + " AND\n"
                   + "   (2) determine object names (i.e. -determineObjectNames (not the default))\n"
                   + " OR\n" + "   (3) split all (i.e. -splitAll (not the default))\n" + " OR\n"
                   + "   (3) consolidate split files (i.e. -consolidate (not the default))\n"
                   + " OR\n" + "   (3) run all (i.e. -runAll (not the default))\n"
                   + "   (4) force the meta-analysis to be redone even if meta-analysis output files exist (i.e. -forceMeta (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("maps=")) {
        mapsFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-determineObjectNames")) {
        determineObjectNames = true;
        numArgs--;
      } else if (arg.startsWith("-splitAll")) {
        splitAll = true;
        numArgs--;
      } else if (arg.startsWith("-consolidate")) {
        consolidate = true;
        numArgs--;
      } else if (arg.startsWith("-runAll")) {
        runAll = true;
        numArgs--;
      } else if (arg.startsWith("-forceMeta")) {
        forceMeta = true;
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

    try {
      log = new Logger(logfile);
      maps = new MetaAnalysisParams(mapsFile, log);
      if (determineObjectNames) {
        determineObjectNames(dir, maps, log);
      } else if (splitAll) {
        splitAll(dir, maps);
      } else if (consolidate) {
        consolidate(dir, maps);
      } else if (runAll) {
        runAll(dir, maps, forceMeta);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
