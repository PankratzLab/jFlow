package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import org.genvisis.cnv.manage.ExportCNVsToPedFormat;
import org.genvisis.cnv.plots.ForestPlot;
import org.genvisis.cnv.plots.QQPlot;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.qsub.Qsub;
import com.google.common.primitives.Booleans;

public class EmimPipeline {

  private static final String[] QQ_PLOT_EFFECTS = {"C", "CM-C", "M"};
  private static final String FOREST_PLOT_DISPLAY_ORDER_NAME = "forestPlotDisplayOrder.txt";

  static class PopFileData {

    String[] pops;
    LinkedHashMap<String, Integer> idsIndexMap;
    boolean[][] inclArr;
  }

  private static PopFileData loadPopFile(String popFile) throws IOException {
    BufferedReader reader = Files.getAppropriateReader(popFile);
    String line = reader.readLine();
    String[] hdrs = line.split("\t", -1);
    ArrayList<String> pops = new ArrayList<>();
    ArrayList<Integer> popIndices = new ArrayList<>();
    for (int i = 0; i < hdrs.length; i++) {
      if (!"".equals(hdrs[i])) {
        popIndices.add(i);
        pops.add(hdrs[i]);
      }
    }
    ArrayList<String> ids = new ArrayList<>();
    @SuppressWarnings("unchecked")
    ArrayList<Boolean>[] incls = new ArrayList[pops.size()];
    for (int i = 0; i < incls.length; i++) {
      incls[i] = new ArrayList<>();
    }
    while ((line = reader.readLine()) != null) {
      String[] incl = line.split("\t", -1);
      for (int i = 0; i < incls.length; i++) {
        incls[i].add(Integer.parseInt(incl[popIndices.get(i)]) == 1);
      }
      ids.add(incl[0] + "\t" + incl[1]);
    }
    reader.close();
    reader = null;

    PopFileData popFileObj = new PopFileData();
    popFileObj.pops = pops.toArray(new String[pops.size()]);
    popFileObj.idsIndexMap = new LinkedHashMap<>();
    for (int i = 0; i < ids.size(); i++) {
      popFileObj.idsIndexMap.put(ids.get(i), i);
    }
    popFileObj.inclArr = new boolean[incls.length][pops.size()];
    for (int i = 0; i < incls.length; i++) {
      popFileObj.inclArr[i] = Booleans.toArray(incls[i]);
    }

    hdrs = null;
    ids = null;
    incls = null;
    pops = null;
    popIndices = null;

    return popFileObj;
  }

  private static void generateKeepsFile(String keepsFile, PopFileData popData, int pop,
                                        PopFileData subPopData, int subPop, Logger log) {
    if (Files.exists(keepsFile)) {
      log.report(keepsFile + " already exists, skipping Keeps File generation");
      return;
    }
    PrintWriter writer = Files.getAppropriateWriter(keepsFile);
    for (String s : popData.idsIndexMap.keySet()) {
      if (popData.inclArr[pop][popData.idsIndexMap.get(s)]) {
        if (subPopData == null || (subPopData.idsIndexMap.containsKey(s)
                                   && subPopData.inclArr[subPop][subPopData.idsIndexMap.get(s)])) {
          writer.println(s);
        }
      }
    }
    writer.flush();
    writer.close();
  }

  private static boolean generateFolderStructureAndKeepsFiles(String runDir, String[] cnvFiles,
                                                              String[] plinkRoots,
                                                              PopFileData popFile,
                                                              PopFileData subPopFile, Logger log) {
    String[] fileroots = combineRoots(cnvFiles, plinkRoots);
    for (String fileroot : fileroots) {
      String dir = runDir + fileroot + "/";
      if (Files.exists(dir)) {
        log.report(dir + " already exists, will not be created");
      } else if (!new File(dir).mkdir()) {
        log.reportError("Could not create " + dir);
        return false;
      }

      for (int p = 0; p < popFile.pops.length; p++) {
        String popDir = dir + ext.replaceWithLinuxSafeCharacters(popFile.pops[p], true) + "/";
        if (Files.exists(popDir)) {
          log.report(popDir + " already exists, will not be created");
        } else if (!new File(popDir).mkdir()) {
          log.reportError("Could not create " + popDir);
          return false;
        }
        generateKeepsFile(popDir + "keeps.txt", popFile, p, null, 0, log);

        for (int sP = 0; sP < subPopFile.pops.length; sP++) {
          String subPopDir = popDir + ext.replaceWithLinuxSafeCharacters(subPopFile.pops[sP], true)
                             + "/";
          if (Files.exists(subPopDir)) {
            log.report(subPopDir + " already exists, will not be created");
          } else if (!new File(subPopDir).mkdir()) {
            log.reportError("Could not create " + subPopDir);
            return false;
          }
          generateKeepsFile(subPopDir + "keeps.txt", popFile, p, subPopFile, sP, log);
        }
      }
    }

    return true;
  }

  static void setup(String runDir, String[] cnvFiles, String[] plinkRoots, String pedFile,
                    String popFile, String subPopFile, String riskAlleleFile, double pThreshold,
                    Set<Emim.EMIM_MODEL> models, boolean phaseWithShapeit, String qsubQueue,
                    Logger log1) {
    ArrayList<String> pbsFiles = new ArrayList<>();
    Logger log = log1 == null ? new Logger() : log1;
    PopFileData popData, subPopData;
    popFile = new File(popFile).getAbsolutePath();
    subPopFile = new File(subPopFile).getAbsolutePath();
    try {
      popData = loadPopFile(popFile);
      subPopData = loadPopFile(subPopFile);
    } catch (IOException e) {
      e.printStackTrace();
      log.reportException(e);
      return;
    }
    if (!checkRoots(cnvFiles, plinkRoots)) {
      log.reportError("Cannot have cnv roots and plink roots that have the same name");
      return;
    }

    if (Files.exists(runDir + FOREST_PLOT_DISPLAY_ORDER_NAME)) {
      log.report(FOREST_PLOT_DISPLAY_ORDER_NAME + " already exists, will not be created");
    } else {
      log.reportTimeInfo("Generating naive Forest Plot Display Order file...");
      ForestPlot.generateNaiveOrderFile(subPopData.pops, runDir + FOREST_PLOT_DISPLAY_ORDER_NAME);
    }
    log.reportTimeInfo("Generating folder structure in run directory...");
    if (!generateFolderStructureAndKeepsFiles(runDir, cnvFiles, plinkRoots, popData, subPopData,
                                              log)) {
      log.reportError("Failed to generate folder structure");
      return;
    }

    if (cnvFiles != null) {
      for (String cnvFile : cnvFiles) {
        String cnvRoot = ext.rootOf(cnvFile, true);
        String cnvDir = cnvRoot + "/";
        String plinkRoot = runDir + cnvRoot + "_0";
        if (!Files.exists(plinkRoot + ".bim") || !Files.exists(plinkRoot + ".bed")
            || !Files.exists(plinkRoot + ".fam")) {
          log.reportTimeInfo("Exporting cnv files for root " + plinkRoot + " to PLINK files...");
          ExportCNVsToPedFormat.export(cnvFile, pedFile, runDir + cnvRoot, "\r\n",
                                       ExportCNVsToPedFormat.PLINK_BINARY_FORMAT, true, true, false,
                                       false, false, false, Integer.MAX_VALUE, 0, log);
        }
        if (!Files.exists(plinkRoot + ".bim") || !Files.exists(plinkRoot + ".bed")
            || !Files.exists(plinkRoot + ".fam")) {
          log.reportError("ERROR - couldn't find exported PLINK files for CNV root " + cnvRoot
                          + " in directory " + runDir);
          continue;
        }
        String relativePlinkRoot = "../" + cnvRoot + "_0";
        String resultFile = cnvRoot;

        log.reportTimeInfo("Generating EMIM files for cnv root " + plinkRoot + "...");
        String pbsFile = Emim.scriptAllInDir(runDir + cnvDir, plinkRoot, relativePlinkRoot, "GEN",
                                             null, riskAlleleFile, pThreshold, models,
                                             phaseWithShapeit, resultFile, log);
        if (pbsFile != null) {
          pbsFiles.add(pbsFile);
        }

        for (String pop : popData.pops) {
          String popDir = cnvDir + ext.replaceWithLinuxSafeCharacters(pop, true) + "/";
          resultFile = cnvRoot + "_" + ext.replaceWithLinuxSafeCharacters(pop, true);
          pbsFile = Emim.scriptAllInDir(runDir + popDir, plinkRoot, "../" + relativePlinkRoot,
                                        "GEN", "keeps.txt", riskAlleleFile, pThreshold, models,
                                        phaseWithShapeit, resultFile, log);
          if (pbsFile != null) {
            pbsFiles.add(pbsFile);
          }

          for (String pop2 : subPopData.pops) {
            String subPopDir = popDir + ext.replaceWithLinuxSafeCharacters(pop2, true) + "/";
            resultFile = cnvRoot + "_" + ext.replaceWithLinuxSafeCharacters(pop, true) + "_"
                         + ext.replaceWithLinuxSafeCharacters(pop2, true);
            pbsFile = Emim.scriptAllInDir(runDir + subPopDir, plinkRoot,
                                          "../../" + relativePlinkRoot, "GEN", "keeps.txt",
                                          riskAlleleFile, pThreshold, models, phaseWithShapeit,
                                          resultFile, log);
            if (pbsFile != null) {
              pbsFiles.add(pbsFile);
            }
          }
        }
      }
    }
    if (plinkRoots != null) {
      for (String plinkRoot : plinkRoots) {
        String plinkDir = plinkRoot + "/";
        String resultFile = plinkRoot;

        log.reportTimeInfo("Generating EMIM files for PLINK root " + plinkRoot + "...");
        String pbsFile = Emim.scriptAllInDir(runDir + plinkDir, plinkRoot, "../" + plinkRoot, "GEN",
                                             null, riskAlleleFile, pThreshold, models,
                                             phaseWithShapeit, resultFile, log);
        if (pbsFile != null) {
          pbsFiles.add(pbsFile);
        }

        for (String pop : popData.pops) {
          String popDir = plinkDir + ext.replaceWithLinuxSafeCharacters(pop, true) + "/";
          resultFile = plinkRoot + "_" + ext.replaceWithLinuxSafeCharacters(pop, true);
          pbsFile = Emim.scriptAllInDir(runDir + popDir, plinkRoot, "../../" + plinkRoot, "GEN",
                                        "keeps.txt", riskAlleleFile, pThreshold, models,
                                        phaseWithShapeit, resultFile, log);
          if (pbsFile != null) {
            pbsFiles.add(pbsFile);
          }

          for (String pop2 : subPopData.pops) {
            String subPopDir = popDir + ext.replaceWithLinuxSafeCharacters(pop2, true) + "/";
            resultFile = plinkRoot + "_" + ext.replaceWithLinuxSafeCharacters(pop, true) + "_"
                         + ext.replaceWithLinuxSafeCharacters(pop2, true);
            pbsFile = Emim.scriptAllInDir(runDir + subPopDir, plinkRoot, "../../../" + plinkRoot,
                                          "GEN", "keeps.txt", riskAlleleFile, pThreshold, models,
                                          phaseWithShapeit, resultFile, log);
            if (pbsFile != null) {
              pbsFiles.add(pbsFile);
            }
          }
        }
      }
    }

    if (pbsFiles.size() > 0) {
      log.reportTimeInfo("Writing QSUB files...");
      writeQsubs(pbsFiles, runDir, qsubQueue);
      // for (String pbsFile : pbsFiles) {
      // CmdLine.run("qsub " + (qsubQueue != null ? "-q " + qsubQueue : "") + pbsFile,
      // ext.parseDirectoryOfFile(pbsFile));
      // }
      log.reportTimeInfo("Pipeline prepared, run " + runDir + "runPBSFiles.sh to submit QSUB jobs");

    } else {
      log.reportError("No pipeline pbs scripts generated, to re-run EMIM remove/rename the existing directories");
    }

    String processCommand = "cd " + runDir + "\n";

    processCommand += Files.getRunString() + " gwas.EmimPipeline -process dir=" + runDir;
    if (cnvFiles != null) {
      processCommand += " cnvs=" + ArrayUtils.toStr(cnvFiles, ",");
    }
    if (plinkRoots != null) {
      processCommand += " plink=" + ArrayUtils.toStr(plinkRoots, ",");
    }
    for (Emim.EMIM_MODEL model : Emim.EMIM_MODEL.optionalSet()) {
      processCommand += " " + model.toString() + "=" + models.contains(model);
    }
    processCommand += " pop=" + popFile + " subPop=" + subPopFile;
    if (log1 != null) {
      processCommand += " log=" + log1.getFilename();
    }
    Qsub.qsub(runDir + "/processResults.pbs", processCommand, 45000, 6, 1);

    log.report("PBS script for post-pipeline processing generated, submit " + runDir
               + "processResults.pbs after the pipeline has completed to process results");

  }

  private static boolean checkRoots(String[] cnvFiles, String[] plinkRoots) {
    HashSet<String> rootSet = new HashSet<>();
    if (plinkRoots != null) {
      for (String root : plinkRoots) {
        if (!rootSet.add(root)) {
          return false;
        }
      }
    }
    if (cnvFiles != null) {
      for (String cnv : cnvFiles) {
        if (!rootSet.add(ext.rootOf(cnv, true))) {
          return false;
        }
      }
    }
    return true;
  }

  private static String[] combineRoots(String[] cnvFiles, String[] plinkRoots) {
    int length = (cnvFiles != null ? cnvFiles.length : 0)
                 + (plinkRoots != null ? plinkRoots.length : 0);
    String[] fileroots = new String[length];
    int i = 0;
    if (cnvFiles != null) {
      for (String cnvFile : cnvFiles) {
        fileroots[i++] = ext.rootOf(cnvFile);
      }
    }
    if (plinkRoots != null) {
      for (String plinkRoot : plinkRoots) {
        fileroots[i++] = plinkRoot;
      }
    }
    return fileroots;
  }

  private static String getResultsFilename(String baseDir, String fileroot, Emim.EMIM_MODEL model) {
    return getResultsFilename(baseDir, fileroot, null, null, model);
  }

  private static String getResultsFilename(String baseDir, String fileroot, String popName,
                                           Emim.EMIM_MODEL model) {
    return getResultsFilename(baseDir, fileroot, popName, null, model);
  }

  private static String getResultsFilename(String baseDir, String fileroot, String popName,
                                           String subPopName, Emim.EMIM_MODEL model) {
    String file = fileroot;
    if (popName != null) {
      popName = ext.replaceWithLinuxSafeCharacters(popName, true);
      file += "_" + popName;
      if (subPopName != null) {
        subPopName = ext.replaceWithLinuxSafeCharacters(subPopName, true);
        file += "_" + subPopName;
      }
    }
    file += "_results_pVals_" + model.toString() + ".xln";
    return getResultsDirectory(baseDir, fileroot, popName, subPopName) + file;
  }

  private static String getResultsDirectory(String baseDir, String fileroot) {
    return getResultsDirectory(baseDir, fileroot, null, null);
  }

  private static String getResultsDirectory(String baseDir, String fileroot, String popName) {
    return getResultsDirectory(baseDir, fileroot, popName, null);
  }

  private static String getResultsDirectory(String baseDir, String fileroot, String popName,
                                            String subPopName) {
    String dir = baseDir + fileroot + "/";
    if (popName != null) {
      popName = ext.replaceWithLinuxSafeCharacters(popName, true);
      dir += popName + "/";
      if (subPopName != null) {
        subPopName = ext.replaceWithLinuxSafeCharacters(subPopName, true);
        dir += subPopName + "/";
      }
    }
    return dir;
  }

  private static void process(String dir, String[] cnvFiles, String[] plinkRoots, String popFile,
                              String subPopFile, double thresh, Set<Emim.EMIM_MODEL> models,
                              Logger log1) throws IOException {
    Logger log = log1 == null ? new Logger() : log1;
    String finalDir = ext.verifyDirFormat(dir);

    PopFileData popData, subPopData;
    try {
      popData = loadPopFile(popFile);
      subPopData = loadPopFile(subPopFile);
    } catch (IOException e) {
      e.printStackTrace();
      log.reportException(e);
      return;
    }

    String[] fileroots = combineRoots(cnvFiles, plinkRoots);

    boolean firstModel = true;
    for (Emim.EMIM_MODEL model : models) {
      if (!model.isOptional()) {
        continue;
        // POO effects are processed with all genotype models
      }
      String finalOut = finalDir + "final_results_pVals_" + model.toString() + ".xln";

      String[] qqHeaders = new String[QQ_PLOT_EFFECTS.length + (firstModel ? 3 : 0)];
      for (int i = 0; i < QQ_PLOT_EFFECTS.length; i++) {
        qqHeaders[i] = "pVal_" + QQ_PLOT_EFFECTS[i] + "_df" + model.getDegreesOfFreedom();
      }
      if (firstModel) {
        qqHeaders[qqHeaders.length - 1] = "tdt_P";
        qqHeaders[qqHeaders.length - 2] = "pVal_Ip_df1";
        qqHeaders[qqHeaders.length - 3] = "pVal_Im_df1";
      }

      PrintWriter writer = Files.getAppropriateWriter(finalOut);
      boolean first = true;
      for (String fileroot : fileroots) {

        writeResults(getResultsFilename(finalDir, fileroot, model), writer, fileroot, "#N/A",
                     "#N/A", thresh, first, log);
        first = false;

        String[] populationResultFiles = new String[popData.pops.length];

        for (int p = 0; p < popData.pops.length; p++) {
          populationResultFiles[p] = getResultsFilename(finalDir, fileroot, popData.pops[p], model);
          writeResults(populationResultFiles[p], writer, fileroot, popData.pops[p], "#N/A", thresh,
                       first, log);

          String[] subPopResultFiles = new String[subPopData.pops.length];

          for (int sP = 0; sP < subPopData.pops.length; sP++) {
            subPopResultFiles[sP] = getResultsFilename(finalDir, fileroot, popData.pops[p],
                                                       subPopData.pops[sP], model);
            writeResults(subPopResultFiles[sP], writer, fileroot, popData.pops[p],
                         subPopData.pops[sP], thresh, first, log);
          }

          generateQQPlots(subPopResultFiles, qqHeaders,
                          getResultsDirectory(finalDir, fileroot, popData.pops[p]),
                          fileroot + "_subpopulations_of_" + ext.replaceWithLinuxSafeCharacters(popData.pops[p],
                                                                                                true),
                          model, log);
        }

        generateQQPlots(populationResultFiles, qqHeaders, getResultsDirectory(finalDir, fileroot),
                        fileroot, model, log);
      }
      writer.flush();
      writer.close();
      firstModel = false;
    }
  }

  private static void generateQQPlots(String[] resultFiles, String[] headers, String baseDir,
                                      String baseName, Emim.EMIM_MODEL model, Logger log) {
    String[][] resultFilesWithCols = new String[headers.length][resultFiles.length];
    for (int i = 0; i < resultFiles.length; i++) {
      if (!Files.exists(resultFiles[i])) {
        log.reportError(resultFiles[i] + " does not exist, QQ Plot cannot be generated");
        return;
      }
      int[] resultFileIndices = ext.indexFactors(headers,
                                                 Files.getHeaderOfFile(resultFiles[i], log), true,
                                                 log, false);
      for (int j = 0; j < headers.length; j++) {
        if (resultFileIndices[j] == -1) {
          log.reportError("Could not find " + headers[j] + " in " + resultFiles[i]);
          return;
        }
        resultFilesWithCols[j][i] = resultFiles[i] + "," + resultFileIndices[j];
      }
    }
    for (int i = 0; i < headers.length; i++) {
      String modelLabel = model.toString() + "_";
      if (headers[i].equals("tdt_P") || headers[i].equals("pVal_Im_df1")
          || headers[i].equals("pVal_Ip_df1")) {
        modelLabel = "";
      }
      String label = baseName + "_" + modelLabel + headers[i];
      QQPlot qqPlot = QQPlot.loadPvals(resultFilesWithCols[i], label, false, true, false, -1, -1,
                                       false, Float.MAX_VALUE, log);
      qqPlot.screenCap(baseDir + label + "_QQPlot.png");
    }

  }

  private static void generateForestPlots(String dir, String[] cnvFiles, String[] plinkRoots,
                                          String popFile, String subPopFile, String forestMarkers,
                                          Set<Emim.EMIM_MODEL> models, Logger log1) {
    Logger log = log1 == null ? new Logger() : log1;
    String finalDir = ext.verifyDirFormat(dir);

    PopFileData popData, subPopData;
    try {
      popData = loadPopFile(popFile);
      subPopData = loadPopFile(subPopFile);
    } catch (IOException e) {
      e.printStackTrace();
      log.reportException(e);
      return;
    }

    String forestPlotDisplayOrderFile = dir + FOREST_PLOT_DISPLAY_ORDER_NAME;
    if (!Files.exists(forestPlotDisplayOrderFile)) {
      ForestPlot.generateNaiveOrderFile(subPopData.pops, forestPlotDisplayOrderFile);
    }

    String[] fileRoots = combineRoots(cnvFiles, plinkRoots);

    for (String fileRoot : fileRoots) {

      for (int p = 0; p < popData.pops.length; p++) {

        for (Emim.EMIM_MODEL model : models) {
          if (!model.isOptional()) {
            // POO Effects are included with all genotypic models
            continue;
          }
          String[][] resultFiles = new String[subPopData.pops.length + 1][];
          resultFiles[subPopData.pops.length] = new String[] {".",
                                                              getResultsFilename(finalDir, fileRoot,
                                                                                 popData.pops[p],
                                                                                 model)};
          for (int sP = 0; sP < subPopData.pops.length; sP++) {
            resultFiles[sP] = new String[] {subPopData.pops[sP],
                                            getResultsFilename(finalDir, fileRoot, popData.pops[p],
                                                               subPopData.pops[sP], model)};
          }
          String forestParameterFile = getResultsDirectory(finalDir, fileRoot, popData.pops[p])
                                       + ext.replaceWithLinuxSafeCharacters(popData.pops[p], true)
                                       + "_forestplot_" + model.toString() + ".xln";
          ResultsPackager.getForestPlotParameterFile(resultFiles, forestMarkers, "MarkerName",
                                                     new String[] {"TDT", "EMIM Child Effect",
                                                                   "EMIM Maternal Effect",
                                                                   "EMIM Maternal POO Effect",
                                                                   "EMIM Paternal POO Effect"},
                                                     new String[][] {{"tdt_OR", "tdt_U95", "tdt_P"},
                                                                     {"C_lnR1", "C_se_lnR1",
                                                                      "pVal_C_df" + model.getDegreesOfFreedom()},
                                                                     {"M_lnS1", "M_se_lnS1",
                                                                      "pVal_M_df" + model.getDegreesOfFreedom(),
                                                                      "pVal_CM-C_df" + model.getDegreesOfFreedom()},
                                                                     {"POO_lnIm", "POO_se_lnIm",
                                                                      "pVal_Im_df1"},
                                                                     {"POO_lnIp", "POO_se_lnIp",
                                                                      "pVal_Ip_df1"}},
                                                     new String[][] {{"", "", "p"}, {"", "", "p"},
                                                                     {"", "", "p",
                                                                      "removing Child Effect p"},
                                                                     {"", "", "p"}, {"", "", "p"}},
                                                     forestParameterFile, log);
          ForestPlot fp = new ForestPlot(ext.rootOf(forestParameterFile, false) + ".input", log);
          fp.loadMarkerFile();
          fp.loadOrderFile(forestPlotDisplayOrderFile, true);
          fp.screenCapAll("ForestPlots/", true, false);
        }
      }
    }

  }

  private static void writeResults(String file, PrintWriter writer, String root1, String root2,
                                   String root3, double thresh, boolean includeHeader,
                                   Logger log) throws IOException {
    if (!Files.exists(file)) {
      log.reportError(file + " does not exist and will not be included in the summarized results");
      return;
    }
    BufferedReader reader = Files.getAppropriateReader(file);

    String line = includeHeader ? null : reader.readLine();
    int pValueColumn = includeHeader ? -1 : ext.indexOfStr("tdt_P", line.split("\t"));
    boolean isHeader = true;
    while ((line = reader.readLine()) != null) {
      if (includeHeader && isHeader) {
        isHeader = false;
        writer.println("DataFileRoot\tPopulation\tSubPopulation\t" + line);
        pValueColumn = ext.indexOfStr("tdt_P", line.split("\t"));
      } else {
        String val = line.split("\t")[pValueColumn];
        if (!ext.isMissingValue(val) && Double.parseDouble(val) < thresh) {
          writer.println(root1 + "\t" + root2 + "\t" + root3 + "\t" + line);
        }
      }
    }
    reader.close();
  }

  private static void writeQsubs(List<String> qsubFiles, String runDir, String qsubQueue) {
    ArrayList<String> qsubCommands = new ArrayList<>();
    for (String qsub : qsubFiles) {

      String dir = ext.parseDirectoryOfFile(qsub);
      qsubCommands.add("cd " + dir);
      qsubCommands.add("qsub " + (qsubQueue != null ? "-q " + qsubQueue : "")
                       + ext.rootOf(qsub, true) + ".pbs");
      String ret = "cd ";
      for (int i = 0; i < dir.split("/").length; i++) {
        ret += "../";
      }
      qsubCommands.add(ret);

    }

    Files.writeIterable(qsubCommands, runDir + "runPBSFiles.sh");
    Files.chmod(runDir + "runPBSFiles.sh");
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String runDir = "./";
    String[] cnvFiles = null;
    String[] plinkRoots = null;
    String pedFile = "./pedigree.dat";
    String popFile = "./pops.xln";
    String subPopFile = "./subPops.xln";
    String riskAlleleFile = null;
    double pThreshold = 1.1;
    String logFile = null;
    Logger log = null;
    String qsub = null;
    boolean process = false;
    boolean forest = false;
    Set<Emim.EMIM_MODEL> models = new HashSet<>(Emim.EMIM_MODEL.valueSet());
    boolean phaseWithShapeit = false;
    String forestMarkers = "./gwasHits.txt";

    String usage = "\n" + "gwas.EmimPipeline requires 2-8 arguments\n"
                   + "   (1) run directory (i.e. dir=" + runDir + " (default))\n" + " AND\n"
                   + "   (2a) cnv files (i.e. cnvs=cnvFile1.cnv,cnvFile2.cnv (not the default))\n"
                   + " AND/OR \n"
                   + "   (2b) PLINK fileroots (i.e. plink=plink1,plink2 (not the default))\n"
                   + " AND\n" + "   (3) population file (i.e. pop=" + popFile + " (default))\n"
                   + " AND\n" + "   (4) subpopulation file (i.e. subPop=" + subPopFile
                   + " (default))\n" + " AND\n"
                   + "   (5) desired risk allele file (i.e. riskAlleles=forceRiskAllele.txt (not the default))\n"
                   + " AND\n" + "   (6) p-value threshold to filter on (i.e. pThreshold="
                   + pThreshold + " (default))\n"
                   + " AND, if desired (though the script to run this will be created automatically)\n"
                   + "   (7) -process flag to consolidate results after PBS files have completed (i.e. -process (not the default))\n"
                   + " AND/OR \n"
                   + "   (8) -forest flag to generate forest plot parameters for a set of markers (i.e. -forest (not the default))\n"
                   + " AND\n"
                   + "   (9) markers to use for forest plot parameter generation (i.e. forestMarkers=./gwasHits.txt (default))\n"
                   + " AND\n" + "   (10) Phase with shapeit in PREMIM run (i.e. phase="
                   + phaseWithShapeit + " (default))\n";
    int argNum = 11;
    for (Emim.EMIM_MODEL model : Emim.EMIM_MODEL.optionalSet()) {
      usage += " AND\n" + "   (" + argNum++ + ") Include " + model.toString() + " model (i.e. "
               + model.toString() + "=true (default))\n";
    }

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (args[i].startsWith("dir=")) {
        runDir = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("cnvs=")) {
        cnvFiles = args[i].split("=")[1].split(",");
        numArgs--;
      } else if (args[i].startsWith("plink=")) {
        plinkRoots = args[i].split("=")[1].split(",");
        numArgs--;
      } else if (args[i].startsWith("log=")) {
        logFile = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("pop=")) {
        popFile = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("subPop=")) {
        subPopFile = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("riskAlleles=")) {
        riskAlleleFile = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("pThreshold=")) {
        pThreshold = ext.parseDoubleArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith("-process")) {
        process = true;
        numArgs--;
      } else if (args[i].startsWith("-forest")) {
        forest = true;
        numArgs--;
      } else if (args[i].startsWith("forestMarkers=")) {
        forestMarkers = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("phase=")) {
        phaseWithShapeit = ext.parseBooleanArg(args[i]);
        numArgs--;
      } else {
        boolean foundModel = false;
        for (Emim.EMIM_MODEL model : Emim.EMIM_MODEL.optionalSet()) {
          if (args[i].startsWith(model.toString() + "=")) {
            foundModel = true;
            if (!ext.parseBooleanArg(args[i])) {
              models.remove(model);
            }
            numArgs--;
            break;
          }
        }
        if (!foundModel) {
          System.err.println("Error - invalid argument: " + args[i]);
        }
      }
    }
    if (args.length == 0 || numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (logFile != null) {
        log = new Logger(logFile);
      }
      runDir = ext.verifyDirFormat(new File(runDir).getCanonicalPath());
      if (process || forest) {
        if (process) {
          System.out.println("Processing EMIM results...");
          process(runDir, cnvFiles, plinkRoots, popFile, subPopFile, pThreshold, models, log);
        }
        if (forest) {
          System.out.println("Generating EMIM forest plots...");
          generateForestPlots(runDir, cnvFiles, plinkRoots, popFile, subPopFile, forestMarkers,
                              models, log);
        }
      } else {
        System.out.println("Preparing EMIM pipeline...");
        if (riskAlleleFile != null) {
          riskAlleleFile = Files.firstPathToFileThatExists(riskAlleleFile, runDir, "");
        }
        setup(runDir, cnvFiles, plinkRoots, pedFile, popFile, subPopFile, riskAlleleFile,
              pThreshold, models, phaseWithShapeit, qsub, log);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}