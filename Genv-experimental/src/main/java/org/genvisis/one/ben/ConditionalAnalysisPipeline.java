package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.genvisis.common.Aliases;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.gwas.FAST;
import org.genvisis.gwas.FAST.DataDefinitions;
import org.genvisis.one.ScriptExecutor;
import org.genvisis.stats.LeastSquares;

import com.google.common.primitives.Doubles;

public class ConditionalAnalysisPipeline {

  static class ConditionalAnalysisToolset_FAST implements Runnable {

    static final String DATA_FILE_TEMPLATE = "chr#.<1>.<2>";
    static final String DATA_FILE_DELIMITER = ".";
    static final String TEMP_GENO_FILE_EXTENSION = ".temp.gz";
    static final String TEMP_INFO_FILE_EXTENSION = ".temp_info";
    static final String MISSING_DATA = ".";
    static final int NUM_THREADS = 24;
    static final double PVAL_THRESHOLD = 0.0001;

    private static void createNewTraitFiles(final Region region, String traitDir,
                                            DataDefinitions dd, boolean baseline) {
      String[] iids = HashVec.loadFileToStringArray(dd.indivFile, false, new int[] {0}, false);
      HashMap<String, Integer> indexMap = new HashMap<String, Integer>();
      for (int i = 0; i < iids.length; i++) {
        indexMap.put(iids[i], i);
      }

      HashMap<String, HashMap<String, HashMap<String, String>>> studyToFactorToPopToFile =
          FAST.loadTraitFiles(traitDir);
      for (Entry<String, HashMap<String, HashMap<String, String>>> studyMap : studyToFactorToPopToFile.entrySet()) {
        String study = studyMap.getKey();
        if (!dd.study.equals(study)) {
          continue;
        }
        for (Entry<String, HashMap<String, String>> factorMap : studyMap.getValue().entrySet()) {
          for (Entry<String, String> popMap : factorMap.getValue().entrySet()) {
            String pop = popMap.getKey();
            if (!dd.popcode.equals(pop)) {
              continue;
            }
            regressNewTraitFile(region, traitDir, popMap.getValue(), indexMap, baseline);
          }
        }
      }

    }

    private static void dumpDosages(Region region, DataDefinitions dd) {
      PrintWriter writer;

      HashMap<String, String[]> pullFrom;
      HashMap<String, String[]> prevFrom;
      String[] dosageData, iids, prevData;

      StringBuilder sb;
      ArrayList<String> prevDataKeyList;
      int offset = 5; // column index offset to start of geno data
      double geno, geno2, geno3;

      pullFrom = region.genoData;
      prevFrom = region.prevSNPdata;
      dosageData = pullFrom.get(dd.study + "\t" + dd.popcode);

      iids = HashVec.loadFileToStringArray(dd.indivFile, false, new int[] {0}, false);
      writer = Files.getAppropriateWriter(region.analysisRootDir + region.label + "_" + dd.study
                                          + "_" + dd.popcode + "_snpDosages.txt");

      sb = new StringBuilder();

      prevDataKeyList = new ArrayList<String>();
      sb.append("IID");
      sb.append("\t").append(dd.study).append("_").append(dd.popcode).append("_")
        .append(region.indexSNP);
      for (String key : prevFrom.keySet()) {
        if (key.startsWith(dd.study + "\t" + dd.popcode)) {
          sb.append("\t").append(key.replaceAll("\t", "_"));
          prevDataKeyList.add(key);
        }
      }
      writer.println(sb.toString());

      for (int i = 0; i < iids.length; i++) {
        sb = new StringBuilder();

        sb.append(iids[i]).append("\t");

        geno2 = Double.parseDouble(dosageData[offset + (3 * i) + 1]);
        geno3 = Double.parseDouble(dosageData[offset + (3 * i) + 2]);
        geno = (geno2 + (2 * geno3));

        sb.append(geno);

        for (String key : prevDataKeyList) {
          prevData = prevFrom.get(key);

          geno2 = Double.parseDouble(prevData[offset + (3 * i) + 1]);
          geno3 = Double.parseDouble(prevData[offset + (3 * i) + 2]);
          geno = (geno2 + (2 * geno3));

          sb.append("\t").append(geno);
        }

        writer.println(sb.toString());
      }

      writer.flush();
      writer.close();
    }

    private static void dumpRegion(Region region, boolean dumpInfo) {
      PrintWriter writer =
          Files.getAppropriateWriter(region.analysisRootDir + region.label
                                     + (dumpInfo ? "_snpInfo.txt" : "_snpData.txt"));
      StringBuilder sb = new StringBuilder();
      ArrayList<String> studyPopOrder = new ArrayList<String>();
      int maxLength = 0;

      HashMap<String, String[]> pullFrom = dumpInfo ? region.infoData : region.genoData;
      HashMap<String, String[]> prevFrom = dumpInfo ? region.prevSNPinfo : region.prevSNPdata;

      int cnt = 0;
      for (Entry<String, String[]> studyPop : pullFrom.entrySet()) {
        sb.append(studyPop.getKey().replace("\t", "_") + "_" + region.indexSNP);
        studyPopOrder.add(studyPop.getKey());
        maxLength = Math.max(maxLength, studyPop.getValue().length);
        if (cnt < pullFrom.size() - 1) {
          sb.append("\t");
        }
        cnt++;
      }
      ArrayList<String> keyListOrder = new ArrayList<String>();
      for (String key : prevFrom.keySet()) {
        keyListOrder.add(key);
        sb.append("\t").append(key.replaceAll("\t", "_"));
      }

      writer.println(sb.toString());

      int currentCount = 0;
      while (currentCount < maxLength) {
        sb = new StringBuilder();
        int keyCnt = 0;
        for (String key : studyPopOrder) {
          String[] line = pullFrom.get(key);
          if (currentCount < line.length) {
            sb.append(line[currentCount]);
          }
          if (keyCnt < studyPopOrder.size() - 1) {
            sb.append("\t");
          }
          keyCnt++;
        }
        for (String key : keyListOrder) {
          sb.append("\t");
          String[] data = prevFrom.get(key);
          if (currentCount < data.length) {
            sb.append(data[currentCount]);
          }
        }
        currentCount++;
        writer.println(sb.toString());
      }

      writer.flush();
      writer.close();
    }

    private static String extractGenoAndInfoDataForRegion(final Region region,
                                                          final DataDefinitions dataDefs,
                                                          String[] dataFiles, String tempDir,
                                                          boolean baseline) {
      // TODO currently only partially reuses data files. This increases space used, but we have
      // problems reading from shared files in a multi-threaded environment
      String tempDataDir =
          tempDir + region.label + "/" /* + (baseline ? "baseline" : region.indexSNP) + "/" */
                           + dataDefs.study + "_" + dataDefs.popcode + "/";
      File dirFile = new File(tempDataDir);
      if (!dirFile.exists()) {
        if (!new File(tempDataDir).mkdirs()) {
          throw new RuntimeException("ERROR - failed to create temporary data directory {"
                                     + tempDataDir + "}");
        }
      }

      String infoHeader =
          "snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0";

      TreeSet<String> sortedFiles = new TreeSet<String>(new Comparator<String>() {
        @Override
        public int compare(String o1, String o2) {
          String[] parts1 = o1.split("\\.");
          String[] parts2 = o2.split("\\.");
          int start1 = Integer.parseInt(parts1[1]);
          // int stop1 = Integer.parseInt(parts1[2]);
          int start2 = Integer.parseInt(parts2[1]);
          // int stop2 = Integer.parseInt(parts2[2]);
          return Integer.valueOf(start1).compareTo(Integer.valueOf(start2));
        }
      });
      for (String file : dataFiles) {
        sortedFiles.add(file);
      }

      String newGenoFileName =
          "chr" + region.chr + "." + region.start + "." + region.stop + ".temp.gz";
      String newInfoFileName =
          "chr" + region.chr + "." + region.start + "." + region.stop + ".temp_info";

      boolean found = false;
      if (Files.exists(tempDataDir + newGenoFileName)
          && Files.exists(tempDataDir + newInfoFileName)) {
        log("WARNING - Skipping data file export, files already exist!");
        log("Loading geno/info data for index SNP [" + region.indexSNP + "]...");
        found = true;
      }

      PrintWriter genoWriter = null;
      PrintWriter infoWriter = null;
      if (!found) {
        genoWriter = Files.getAppropriateWriter(tempDataDir + newGenoFileName);
        infoWriter = Files.getAppropriateWriter(tempDataDir + newInfoFileName);

        infoWriter.println(infoHeader);
      }

      for (String file : sortedFiles) {
        try {
          BufferedReader genoReader = Files.getAppropriateReader(dataDefs.dataDir + file);
          BufferedReader infoReader =
              Files.getAppropriateReader(dataDefs.dataDir + file.substring(0, file.length() - 3)
                                         + "_info");

          String infoLine = infoReader.readLine();
          String delim = ext.determineDelimiter(infoLine);
          String genoLine = null;
          int count = 0;
          while ((infoLine = infoReader.readLine()) != null
                 && (genoLine = genoReader.readLine()) != null) {
            count++;
            // geno/info lines should be one to one
            String[] infoParts = infoLine.split(delim);
            if (!baseline && infoParts[1].trim().equals(region.indexSNP)) {
              region.genoData.put(dataDefs.study + "\t" + dataDefs.popcode,
                                  genoLine.split("[\\s]+"));
              region.infoData.put(dataDefs.study + "\t" + dataDefs.popcode,
                                  infoLine.split("[\\s]+"));
              // if (found) {
              // break; // don't break anymore - we have to load data for each index SNP we've
              // tracked
              // }
            } else if (region.prevSNPs.contains(dataDefs.study + "\t" + dataDefs.popcode + "\t"
                                                + infoParts[1])) {
              region.prevSNPdata.put(dataDefs.study + "\t" + dataDefs.popcode + "\t" + infoParts[1],
                                     genoLine.split("[\\s]+"));
              region.prevSNPinfo.put(dataDefs.study + "\t" + dataDefs.popcode + "\t" + infoParts[1],
                                     infoLine.split("[\\s]+"));
            }
            int mkrPos = Integer.parseInt(infoParts[2]);
            if (mkrPos < region.start || mkrPos > region.stop) {
              continue;
            }
            if (!found) {
              genoWriter.println(genoLine);
              infoWriter.println(infoLine);
            }
          }
          log("Read " + count + " data/info lines...");

          genoReader.close();
          infoReader.close();
        } catch (FileNotFoundException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        } catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      }
      if (!found) {
        infoWriter.flush();
        genoWriter.flush();
        infoWriter.close();
        genoWriter.close();
      }

      return tempDataDir;
    }

    private static String[] findDataFiles(final Region region, final DataDefinitions dataDefs) {
      String[] chrDataFiles = (new File(dataDefs.dataDir)).list(new FilenameFilter() {
        @Override
        public boolean accept(File dir, String name) {
          // TODO assuming datafile names start with "chr#."
          return name.startsWith("chr" + region.chr + DATA_FILE_DELIMITER)
                 && name.endsWith(dataDefs.dataSuffix);
        }
      });

      if (chrDataFiles.length == 0) {
        throw new RuntimeException("ERROR - no data files found in directory.  Looking for files with pattern 'chr#.<>.<>' or 'chr#' ending with suffix "
                                   + dataDefs.dataSuffix + " in directory {" + dataDefs.dataDir
                                   + "}");
      }

      // TODO generic parsing for file-name template: is it in chr#.<>.<> format [chunked], or chr#
      // format [whole_chr]?
      // assuming data file names include position
      ArrayList<String> inclDataFiles = new ArrayList<String>();
      for (String chrDataFile : chrDataFiles) {
        String[] pts = chrDataFile.split("\\.");
        int chunkStart = Integer.parseInt(pts[1]);
        int chunkStop = Integer.parseInt(pts[2]);
        boolean checkStartOverlap = region.start <= chunkStart && region.stop >= chunkStart;
        boolean checkOverlap1 = region.start >= chunkStart && region.stop <= chunkStop;
        boolean checkOverlap2 = region.start <= chunkStart && region.stop >= chunkStop;
        boolean checkStopOverlap = region.start < chunkStop && region.stop > chunkStop;
        if (checkStartOverlap || checkOverlap1 || checkOverlap2 || checkStopOverlap) {
          inclDataFiles.add(chrDataFile);
        }
      }
      log("Found " + inclDataFiles.size() + " out of " + chrDataFiles.length
          + " possibly-matching data files!");
      return inclDataFiles.toArray(new String[inclDataFiles.size()]);
    }

    private static void regressNewTraitFile(Region region, String traitDir, String traitFile,
                                            HashMap<String, Integer> iids, boolean baseline) {
      String[] pts = traitFile.substring(0, traitFile.lastIndexOf(".")).split("_");
      String study = pts[0];
      String pop = pts[1];
      String factor = pts[2];
      String newTraitFile = study + "_" + pop + "_" + factor + ".trait";
      String dir = region.analysisRootDir
                   + (baseline ? region.label + "_iter0_baseline/" : region.regionDirNameRoot);

      int offset = 5; // column index offset to start of geno data

      ArrayList<String> missing = new ArrayList<String>();

      try {
        int traitCount = Files.countLines(traitDir + traitFile, 1);
        BufferedReader reader = Files.getAppropriateReader(traitDir + traitFile);

        int phenoCol = 5;
        //
        double[] phenoData = new double[traitCount];
        // double[][] indepData = new double[traitCount][];
        ArrayList<double[]> indepDataLines = new ArrayList<double[]>();

        String line = reader.readLine(); // header
        String[] parts = line.split("[\\s]+");
        ArrayList<String> colNames = new ArrayList<String>();

        for (int i = 6; i < parts.length; i++) {
          colNames.add(parts[i]);
        }
        if (!baseline) {
          colNames.add(region.indexSNP);
          for (String snp : region.prevSNPs) {
            if (snp.contains(study) && snp.contains(pop)) {
              colNames.add(snp.split("\t")[2]);
            }
          }
        }

        int cnt = 0;
        ArrayList<Double> phenoDataList = new ArrayList<Double>();

        while ((line = reader.readLine()) != null) {
          parts = line.split("[\\s]+");

          ArrayList<Double> lineData = new ArrayList<Double>();

          String iid = parts[1];
          Integer iidIndex = iids.get(iid);
          double geno = Double.NaN;
          if (iidIndex == null) {
            missing.add(iid);
            continue;
          } else {
            for (int i = 6; i < parts.length; i++) {
              lineData.add(Double.parseDouble(parts[i]));
            }

            if (!baseline) {
              int iidInd = iidIndex.intValue();
              // double geno1 = Double.parseDouble(genoData[offset + (3 * iidInd)]);
              double geno2 =
                  Double.parseDouble(region.genoData.get(study + "\t" + pop)[offset + (3 * iidInd)
                                                                             + 1]);
              double geno3 =
                  Double.parseDouble(region.genoData.get(study + "\t" + pop)[offset + (3 * iidInd)
                                                                             + 2]);
              geno = (geno2 + (2 * geno3));
              lineData.add(geno);

              // TODO include each previous SNP? Comment out to remove
              for (String snp : region.prevSNPs) {
                if (snp.contains(study) && snp.contains(pop)) {
                  // log("Including data for: {" + snp + "}");
                  // log("Keys: {" + region.prevSNPdata.keySet() + "}");
                  geno2 =
                      Double.parseDouble(region.prevSNPdata.get(snp)[offset + (3 * iidInd) + 1]);
                  geno3 =
                      Double.parseDouble(region.prevSNPdata.get(snp)[offset + (3 * iidInd) + 2]);
                  geno = (geno2 + (2 * geno3));
                  lineData.add(geno);
                }
              }
            }

            phenoDataList.add(Double.parseDouble(parts[phenoCol]));

            indepDataLines.add(Doubles.toArray(lineData));
            cnt++;
          }

        }
        reader.close();

        phenoData = Doubles.toArray(phenoDataList); // Double.parseDouble(parts[phenoCol]);
        double[][] indepData = indepDataLines.toArray(new double[0][]);

        String[] cols = colNames.toArray(new String[colNames.size()]);

        LeastSquares lsReg = new LeastSquares(phenoData, indepData, cols, false, false); // TODO
                                                                                         // should
                                                                                         // bypass
                                                                                         // data
                                                                                         // check??
                                                                                         // - could
                                                                                         // easily
                                                                                         // have
                                                                                         // NaNs in
                                                                                         // there

        double[] resids = lsReg.getResiduals();
        if (resids.length != traitCount) {
          // TODO error!
        }

        String header = "#Fam_ID\tInd_ID\tDad_ID\tMom_ID\tSex\tPhenotype";

        PrintWriter writer = Files.getAppropriateWriter(dir + newTraitFile);
        writer.println(header);

        reader = Files.getAppropriateReader(traitDir + traitFile);
        reader.readLine(); // header
        cnt = 0;
        while ((line = reader.readLine()) != null) {
          parts = line.split("[\\s]+");
          String iid = parts[1];
          Integer iidIndex = iids.get(iid);

          if (iidIndex != null) {
            StringBuilder lineStr = new StringBuilder();
            lineStr.append(parts[0]).append("\t").append(parts[1]).append("\t").append(parts[2])
                   .append("\t").append(parts[3]).append("\t").append(parts[4]).append("\t");

            lineStr.append(resids[cnt++]);

            writer.println(lineStr.toString());
          } else {
            // lineStr.append(".");
          }

          // cnt++;
        }
        reader.close();
        writer.flush();
        writer.close();

      } catch (FileNotFoundException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }

    Region region;
    String dataFile;
    String traitDir;
    String tempDir;

    HashMap<String, HashMap<String, DataDefinitions>> dataDefs;

    boolean isBaseline;

    public ConditionalAnalysisToolset_FAST(Region region, String dataFile, String traitDir,
                                           String tempDir,
                                           HashMap<String, HashMap<String, DataDefinitions>> dataDefs,
                                           boolean baseline) {
      this.region = region;
      this.dataFile = dataFile;
      this.traitDir = traitDir;
      this.tempDir = tempDir;
      this.dataDefs = dataDefs;
      isBaseline = baseline;
    }

    @Override
    public void run() {
      String[] regionPathAndDataFile = setup();

      try {
        log("Preparing FAST analysis in directory [" + regionPathAndDataFile[0] + "]...");
        String[] analysisDirs = FAST.prepareFAST(regionPathAndDataFile[0], regionPathAndDataFile[1],
                                                 regionPathAndDataFile[0], true, false, true, null);

        log("Running " + analysisDirs.length + " FAST analyses...");
        boolean[] runs = Array.booleanArray(analysisDirs.length, false);
        for (int i = 0; i < analysisDirs.length; i++) {
          (new ScriptExecutor(NUM_THREADS)).run(analysisDirs[i] + "input.txt", "took");
          runs[i] = ScriptExecutor.outLogExistsComplete(analysisDirs[i] + "output/input.log_0.out",
                                                        "took");
          if (!runs[i]) {
            // TODO Error - FAST failed!
          }
        }

        log("Processing FAST results...");
        for (String study : dataDefs.keySet()) {
          String studyDir = regionPathAndDataFile[0] + study + "/";

          if (!Files.exists(studyDir)) {
            System.err.println(ext.getTime() + "ERROR - analysis directory [" + studyDir
                               + "] does not exist.");
            continue;
          }

          FAST.processAndPrepareMETAL(studyDir, dataDefs, true);

          HashMap<String, DataDefinitions> popDefs = dataDefs.get(study);

          // should only ever be one directory...
          String[] factorDirs = Files.listDirectories(studyDir, false);
          for (String factorDir : factorDirs) {
            String dir = studyDir + factorDir + "/";
            String file = null;
            if (popDefs.size() == 1) {
              // one population, no meta analysis from which to pull results
              file = popDefs.keySet().toArray(new String[1])[0] + "/output/concatenated.result";
            } else {
              file = factorDir + "_InvVar1.out";
            }

            if (Files.exists(dir + file)) {
              String newSNP = extractIndexSnp(dir + file, region, PVAL_THRESHOLD);

              if (isBaseline) {
                log("Baseline analysis for region " + region.label
                    + " complete!  Index snp for analysis determined to be [" + newSNP + "]");
              } else if (newSNP == null) {
                dumpRegion(region, false); // dump data
                dumpRegion(region, true); // dump info

                for (DataDefinitions dd : popDefs.values()) {
                  dumpDosages(region, dd);
                }

                log("Couldn't find a candidate SNP for iterative analysis; recursive analysis for region ["
                    + region.label + "] complete.");
              } else {
                log("Iterating analysis with most-significant SNP [" + newSNP + "]");
                Region r2 = new Region();
                r2.chr = region.chr;
                r2.start = region.start;
                r2.stop = region.stop;
                r2.label = region.label;
                r2.indexSNP = newSNP;

                String newDir =
                    r2.label + "_iter" + ((region.prevSNPs.size() / popDefs.size()) + 2) + "_"
                                + ext.replaceWithLinuxSafeCharacters(r2.indexSNP, false) + "/";
                new File(region.analysisRootDir + newDir).mkdirs();
                r2.analysisRootDir = region.analysisRootDir;
                r2.regionDirNameRoot = newDir;
                r2.prevSNPs.addAll(region.prevSNPs);
                for (String def : popDefs.keySet()) {
                  r2.prevSNPs.add(study + "\t" + def + "\t" + region.indexSNP);
                }
                (new ConditionalAnalysisToolset_FAST(r2, dataFile, traitDir, tempDir, dataDefs,
                                                     false)).run();
              }
            } else {
              log("Error - file [" + dir + file + "] not found!");
              // TODO error message - result file not found!
            }

          }

        }

      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }

    private String[] setup() {
      // Study -> PopCode -> Defs
      ArrayList<DataDefinitions> allDefs = new ArrayList<FAST.DataDefinitions>();
      for (HashMap<String, DataDefinitions> sub : dataDefs.values()) {
        allDefs.addAll(sub.values());
      }

      ArrayList<String> newDataDefs = new ArrayList<String>();
      for (DataDefinitions dd : allDefs) {
        log("Retrieving required data files...");
        String[] files = findDataFiles(region, dd);

        log("Extracting region-specific geno/info data...");
        String dir = extractGenoAndInfoDataForRegion(region, dd, files, tempDir, isBaseline);

        log("Creating new .trait files...");
        createNewTraitFiles(region, traitDir, dd, isBaseline);

        StringBuilder newDef = new StringBuilder();
        newDef.append(dd.study).append("\t").append(dd.popcode).append("\t").append(dir)
              .append("\t").append(".temp.gz").append("\t")
              // TODO sex-specific
              .append(dd.indivFile);
        newDataDefs.add(newDef.toString());
      }

      String regionDir = isBaseline ? region.label + "_iter0_baseline/" : region.regionDirNameRoot;

      log("Writing new data.txt file...");
      String newDataFile = region.analysisRootDir + regionDir
                           + "data_" + (isBaseline ? "baseline"
                                                   : ext.replaceWithLinuxSafeCharacters(region.indexSNP,
                                                                                        false))
                           + ".txt";
      Files.writeList(newDataDefs.toArray(new String[newDataDefs.size()]), newDataFile);

      return new String[] {region.analysisRootDir + regionDir, newDataFile};
    }

  }

  static class Region {
    String label;
    String indexSNP;
    int chr;
    int start;
    int stop;
    String analysisRootDir;
    String regionDirNameRoot;
    // set programmatically:
    HashMap<String, String[]> genoData = new HashMap<String, String[]>();
    HashMap<String, String[]> infoData = new HashMap<String, String[]>();

    HashSet<String> prevSNPs = new HashSet<String>();
    HashMap<String, String[]> prevSNPdata = new HashMap<String, String[]>();
    HashMap<String, String[]> prevSNPinfo = new HashMap<String, String[]>();

    @Override
    public String toString() {
      StringBuilder sb = new StringBuilder("region[id:" + label + ", SNP:" + indexSNP + ", UCSC:chr"
                                           + chr + ":" + start + ":" + stop + "]");
      return sb.toString();
    }

  }

  private static final Object PRINT_LOCK = new Object();

  private static String extractIndexSnp(String resultFile, Region region,
                                        double thresh) throws IOException {
    String[][] aliases = {Aliases.MARKER_NAMES, Aliases.PVALUES};

    BufferedReader reader = Files.getAppropriateReader(resultFile);

    String line = reader.readLine();

    String delim = ext.determineDelimiter(line);
    String[] hdr = line.split(delim);
    int[] indices = ext.indexFactors(aliases, hdr, false, true, false, false);

    double minPVal = thresh;
    // String minPValStr = null;
    String minPValSNP = null;
    while ((line = reader.readLine()) != null) {
      String[] parts = line.split(delim);
      // if (parts[indices[0]].equals(region.indexSNP) ||
      // region.prevSNPs.contains(parts[indices[0]])) {
      // continue;
      // }
      double pval = thresh;
      try {
        pval = Double.parseDouble(parts[indices[1]]);
      } catch (NumberFormatException e) {
      }
      if (pval < minPVal) {
        minPVal = pval;
        // minPValStr = parts[indices[1]];
        minPValSNP = parts[indices[0]];
      }
    }
    reader.close();

    return minPValSNP;
  }


  private static String[] extractSnpDetails(String resultFile, String snp) throws IOException {
    String[][] aliases = {Aliases.MARKER_NAMES, Aliases.EFFECTS, Aliases.STD_ERRS, Aliases.PVALUES};

    BufferedReader reader = Files.getAppropriateReader(resultFile);

    String line = reader.readLine();

    String delim = ext.determineDelimiter(line);
    String[] hdr = line.split(delim);
    int[] indices = ext.indexFactors(aliases, hdr, false, true, false, false);

    String[] deets = new String[3];
    while ((line = reader.readLine()) != null) {
      String[] parts = line.split(delim);
      if (parts[indices[0]].equals(snp)) {
        deets[0] = parts[indices[1]];
        deets[1] = parts[indices[2]];
        deets[2] = parts[indices[3]];
        break;
      }
    }
    reader.close();

    return deets;
  }

  private static String[] getIterDirs(final Region region) {
    String[] iterDirs = (new File(region.analysisRootDir)).list(new FilenameFilter() {
      @Override
      public boolean accept(File arg0, String arg1) {
        return arg1.startsWith(region.label + "_iter");
      }
    });

    Arrays.sort(iterDirs, new Comparator<String>() {
      @Override
      public int compare(String o1, String o2) {
        String[] parts1 = o1.split("_");
        String[] parts2 = o2.split("_");
        String iter1 = parts1[1].substring(4);
        String iter2 = parts2[1].substring(4);
        return new Integer(iter1).compareTo(new Integer(iter2));
      }
    });

    return iterDirs;
  }

  // private HashMap<String, HashMap<String, String>> performBaselineAnalyses(String dir, String
  // dataFile, Region region, HashMap<String, HashMap<String, DataDefinitions>> dataDefs) {
  // HashMap<String, HashMap<String, String>> studyFactorSnp = new HashMap<String,
  // HashMap<String,String>>();
  //
  // String newDir = region.label + "_iter0_baseline/";
  // if(!new File(dir + newDir).mkdirs()) {
  // // TODO failed!
  // }
  //
  // synchronized(ConditionalAnalysisToolset_FAST.PRINT_LOCK) { log("Preparing baseline FAST
  // analysis in directory [" + newDir + "]..."); }
  // try {
  // String[] analysisDirs = FAST.prepareFAST(dir, dataFile, dir + newDir, false);
  //
  // synchronized(ConditionalAnalysisToolset_FAST.PRINT_LOCK) { log("Running " + analysisDirs.length
  // + " baseline FAST analyses..."); }
  // boolean[] runs = Array.booleanArray(analysisDirs.length, false);
  // for (int j = 0; j < analysisDirs.length; j++) {
  // (new ScriptExecutor(ConditionalAnalysisToolset_FAST.NUM_THREADS)).run(analysisDirs[j] +
  // "input.txt", "took");
  // runs[j] = ScriptExecutor.outLogExistsComplete(analysisDirs[j] + "output/input.log_0.out",
  // "took");
  // if (!runs[j]) {
  // // TODO Error - FAST failed!
  // }
  // }
  //
  // synchronized(ConditionalAnalysisToolset_FAST.PRINT_LOCK) { log("Processing baseline FAST
  // results..."); }
  // for (String study : dataDefs.keySet()) {
  // HashMap<String, String> factorSnpMap = new HashMap<String, String>();
  // studyFactorSnp.put(study, factorSnpMap);
  //
  // String studyDir = dir + newDir + study + "/";
  //
  // if (!Files.exists(studyDir)) {
  // System.err.println(ext.getTime() + "ERROR - baseline analysis directory [" + studyDir + "] does
  // not exist.");
  // continue;
  // }
  //
  // FAST.processAndPrepareMETAL(studyDir);
  //
  // HashMap<String, DataDefinitions> popDefs = dataDefs.get(study);
  //
  // // should only ever be one directory...
  // String[] factorDirs = Files.listDirectories(studyDir, false);
  // for (String factorDir : factorDirs) {
  // String dir2 = studyDir + factorDir + "/";
  // String file = null;
  // if (popDefs.size() == 1) {
  // // one population, no meta analysis from which to pull results
  // file = popDefs.keySet().toArray(new String[1])[0] + "/output/concatenated.result";
  // } else {
  // file = factorDir + "_InvVar1.out";
  // }
  //
  // if (Files.exists(dir2 + file)) {
  // String newSNP = ConditionalAnalysisToolset_FAST.extractIndexSnp(dir2 + file, region,
  // ConditionalAnalysisToolset_FAST.PVAL_THRESHOLD);
  //
  // if (newSNP == null) {
  // ConditionalAnalysisToolset_FAST.dumpRegion(region);
  //
  // synchronized (ConditionalAnalysisToolset_FAST.PRINT_LOCK) {
  // log("Couldn't find a candidate SNPin baseline analysis for region [" + region.label + "]");
  // }
  // } else {
  // synchronized (ConditionalAnalysisToolset_FAST.PRINT_LOCK) {
  // log("Iterating analysis with most-significant SNP [" + newSNP + "]");
  // }
  //
  // factorSnpMap.put(factorDir, newSNP);
  // }
  // } else {
  // synchronized(ConditionalAnalysisToolset_FAST.PRINT_LOCK) {
  // log("Error - file [" + dir2 + file + "] not found!");
  // }
  // // TODO error message - result file not found!
  // }
  //
  // }
  //
  // }
  //
  //
  // } catch (IOException e) {
  // // TODO Auto-generated catch block
  // e.printStackTrace();
  // }
  //
  // return studyFactorSnp;
  // }

  public static void log(String msg) {
    synchronized (PRINT_LOCK) {
      System.out.println(ext.getTime() + "]\t" + msg);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;

    String analysisDir = "";
    String inputFile = "";
    String dataFile = "";
    String tempDataDir = "";
    String traitDir = "";

    String usage = "WRONG";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        analysisDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("input=")) {
        inputFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("data=")) {
        dataFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("temp=")) {
        tempDataDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("trait=")) {
        traitDir = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.err.println(Array.toStr(args, "\n"));
      System.exit(1);
    }
    (new ConditionalAnalysisPipeline()).setup(ext.verifyDirFormat(analysisDir), inputFile, dataFile,
                                              ext.verifyDirFormat(tempDataDir),
                                              ext.verifyDirFormat(traitDir));
  }

  private static Region[] parseSetupFile(String file) {
    String[][] rgnDefs = HashVec.loadFileToStringMatrix(file, false, new int[] {0, 1, 2}, false);
    Region[] rgns = new Region[rgnDefs.length];
    for (int i = 0; i < rgnDefs.length; i++) {
      Region rgn = new Region();
      rgn.label = rgnDefs[i][0];
      rgn.indexSNP = rgnDefs[i][1];
      int[] tmp = Positions.parseUCSClocation(rgnDefs[i][2]);
      rgn.chr = tmp[0];
      rgn.start = tmp[1];
      rgn.stop = tmp[2];
      rgns[i] = rgn;
    }
    return rgns;
  }

  private static void processAndWriteResults(final Region region,
                                             HashMap<String, HashMap<String, DataDefinitions>> dataDefs) {
    String[] iterDirs = getIterDirs(region);
    String[][] factors =
        new String[][] {Aliases.MARKER_NAMES, Aliases.EFFECTS, Aliases.STD_ERRS, Aliases.PVALUES};

    HashMap<String, HashMap<String, StringBuilder>> headerMap =
        new HashMap<String, HashMap<String, StringBuilder>>();
    HashMap<String, HashMap<String, HashMap<String, StringBuilder>>> resultsMap =
        new HashMap<String, HashMap<String, HashMap<String, StringBuilder>>>();

    for (String regionIter : iterDirs) {
      String iterMarker = regionIter.split("_")[2];

      String iterPath = ext.verifyDirFormat(region.analysisRootDir + regionIter);

      for (String study : dataDefs.keySet()) {
        String iterStudyDir = iterPath + study + "/";

        if (!Files.exists(iterStudyDir)) {
          // System.err.println(ext.getTime() + "ERROR - analysis directory [" + iterStudyDir + "]
          // does not exist.");
          continue;
        }

        HashMap<String, StringBuilder> factorHeaderMap = headerMap.get(study);
        if (factorHeaderMap == null) {
          factorHeaderMap = new HashMap<String, StringBuilder>();
          headerMap.put(study, factorHeaderMap);
        }
        HashMap<String, HashMap<String, StringBuilder>> factorResultsMap = resultsMap.get(study);
        if (factorResultsMap == null) {
          factorResultsMap = new HashMap<String, HashMap<String, StringBuilder>>();
          resultsMap.put(study, factorResultsMap);
        }

        HashMap<String, DataDefinitions> popDefs = dataDefs.get(study);

        // should only ever be one directory...
        String[] factorDirs = Files.listDirectories(iterStudyDir, false);
        for (String factorDir : factorDirs) {
          String dir = iterStudyDir + factorDir + "/";
          String file = null;
          if (popDefs.size() == 1) {
            // one population, no meta analysis from which to pull results
            file = popDefs.keySet().toArray(new String[1])[0] + "/output/concatenated.result";
          } else {
            file = factorDir + "_InvVar1.out";
          }

          if (Files.exists(dir + file)) {
            StringBuilder sb = factorHeaderMap.get(factorDir);
            if (sb == null) {
              sb = new StringBuilder("MarkerName");
              factorHeaderMap.put(factorDir, sb);
            }
            sb.append("\t").append(iterMarker).append("_beta\t").append(iterMarker).append("_SE\t")
              .append(iterMarker).append("_pval");
            factorHeaderMap.put(factorDir, sb); // because StringBuilders are immutable, we need to
                                                // replace the instance each time [is this true?]

            int[] indices = ext.indexFactors(factors, Files.getHeaderOfFile(dir + file, null),
                                             false, true, false, false);
            String[][] fileData = HashVec.loadFileToStringMatrix(dir + file, true, indices, false);

            HashMap<String, StringBuilder> markerResultsMap = factorResultsMap.get(factorDir);
            if (markerResultsMap == null) {
              markerResultsMap = new HashMap<String, StringBuilder>();
              factorResultsMap.put(factorDir, markerResultsMap);
            }

            for (String[] markerData : fileData) {
              String mkr = markerData[0];
              String beta = markerData[1];
              String se = markerData[2];
              String pval = markerData[3];

              StringBuilder markerResults = markerResultsMap.get(mkr);
              if (markerResults == null) {
                markerResults = new StringBuilder(mkr);
                markerResultsMap.put(mkr, markerResults);
              }
              markerResults.append("\t").append(beta).append("\t").append(se).append("\t")
                           .append(pval);
              markerResultsMap.put(mkr, markerResults);
            }
          }
        }
      }
    }

    for (Entry<String, HashMap<String, StringBuilder>> headerByStudy : headerMap.entrySet()) {
      String study = headerByStudy.getKey();

      for (Entry<String, StringBuilder> headerByFactorDir : headerByStudy.getValue().entrySet()) {
        String factor = headerByFactorDir.getKey();

        String file =
            region.analysisRootDir + region.label + "_" + study + "_" + factor + "_iterations.xln";

        PrintWriter writer = Files.getAppropriateWriter(file);
        writer.println(headerByFactorDir.getValue().append("\titerationMin").toString());

        for (StringBuilder sb : resultsMap.get(study).get(factor).values()) {
          String line = sb.toString();
          for (int i = 0; i < iterDirs.length; i++) {
            String iterMarker = iterDirs[i].split("_")[2];
            if (line.startsWith(iterMarker)) {
              line = line + "\t" + (i - 1);
            }
          }
          writer.println(line);
        }
        writer.flush();
        writer.close();
      }
    }

  }

  private static void processAndWriteResults2(final Region region,
                                              HashMap<String, HashMap<String, DataDefinitions>> dataDefs) {
    String[] iterDirs = getIterDirs(region);

    HashMap<String, HashMap<String, StringBuilder>> headerMap =
        new HashMap<String, HashMap<String, StringBuilder>>();
    HashMap<String, HashMap<String, ArrayList<StringBuilder>>> resultsMap =
        new HashMap<String, HashMap<String, ArrayList<StringBuilder>>>();
    ArrayList<String> popCodeOrder = new ArrayList<String>();

    for (String study : dataDefs.keySet()) {
      HashMap<String, StringBuilder> factorMap = new HashMap<String, StringBuilder>();
      headerMap.put(study, factorMap);
      String iterPath = ext.verifyDirFormat(region.analysisRootDir + iterDirs[0]);
      String iterStudyDir = iterPath + study + "/";
      String[] factorDirs = Files.listDirectories(iterStudyDir, false);
      for (String factorDir : factorDirs) {
        StringBuilder hdr1SB = new StringBuilder(study).append("\t").append(factorDir)
                                                       .append("\t\t\tMeta\tMeta\tMeta");
        String perHdr = "\tbeta\tSE\tpval";
        StringBuilder hdr2SB = new StringBuilder("Condition\tTopHit\tChr\tPosition").append(perHdr);
        HashMap<String, DataDefinitions> popDefs = dataDefs.get(study);
        for (DataDefinitions dd : popDefs.values()) {
          popCodeOrder.add(dd.popcode);
          hdr1SB.append("\t").append(dd.popcode).append("\t").append(dd.popcode).append("\t")
                .append(dd.popcode);
          hdr2SB.append(perHdr);
        }
        factorMap.put(factorDir, hdr1SB.append("\n").append(hdr2SB));
      }
    }

    for (String study : dataDefs.keySet()) {
      HashMap<String, ArrayList<StringBuilder>> map =
          new HashMap<String, ArrayList<StringBuilder>>();
      resultsMap.put(study, map);
    }

    for (String study : dataDefs.keySet()) {
      for (int i = 0; i < iterDirs.length; i++) {

        String iterMarker = iterDirs[i].split("_")[2];
        String iterPath = ext.verifyDirFormat(region.analysisRootDir + iterDirs[i]);
        String iterStudyDir = iterPath + study + "/";
        String[] factorDirs = Files.listDirectories(iterStudyDir, false);

        for (String factorDir : factorDirs) {
          ArrayList<StringBuilder> iterSBs = resultsMap.get(study).get(factorDir);
          if (iterSBs == null) {
            iterSBs = new ArrayList<StringBuilder>();
            resultsMap.get(study).put(factorDir, iterSBs);
          }

          StringBuilder resultLine = new StringBuilder();
          if (i == 0) { // baseline
            resultLine.append(iterMarker);
          } else {
            for (int j = 1; j < i + 1; j++) {
              resultLine.append(iterDirs[j].split("_")[2]);
              if (j < i) {
                resultLine.append(",");
              }
            }
          }

          String dir = iterStudyDir + factorDir + "/";

          HashMap<String, DataDefinitions> popDefs = dataDefs.get(study);

          String metaFile = null;
          if (popDefs.size() == 1) {
            // one population, no meta analysis from which to pull results
            metaFile = popDefs.keySet().toArray(new String[1])[0] + "/output/concatenated.result";
          } else {
            metaFile = factorDir + "_InvVar1.out";
          }

          HashMap<String, String> popFiles = new HashMap<String, String>();
          for (String popCode : popDefs.keySet()) {
            popFiles.put(popCode, popCode + "/output/concatenated.result");
          }

          String minSNP = null;
          try {
            double thresh =
                i == iterDirs.length - 1 ? 1d : ConditionalAnalysisToolset_FAST.PVAL_THRESHOLD;
            minSNP = extractIndexSnp(dir + metaFile, region, thresh);
            resultLine.append("\t").append(minSNP).append("\t\t\t");
            resultLine.append(Array.toStr(extractSnpDetails(dir + metaFile, minSNP)));
            for (String code : popCodeOrder) {
              String[] newParts = extractSnpDetails(dir + popFiles.get(code), minSNP);
              if (newParts[0] == null) {
                newParts[0] = ".";
              }
              if (newParts[1] == null) {
                newParts[1] = ".";
              }
              if (newParts[2] == null) {
                newParts[2] = ".";
              }
              resultLine.append("\t").append(Array.toStr(newParts));
            }
          } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
          }

          iterSBs.add(resultLine);
        }
      }
    }


    for (String study : dataDefs.keySet()) {
      String iterPath = ext.verifyDirFormat(region.analysisRootDir + iterDirs[0]);
      String iterStudyDir = iterPath + study + "/";
      String[] factorDirs = Files.listDirectories(iterStudyDir, false);
      for (String factorDir : factorDirs) {
        PrintWriter writer =
            Files.getAppropriateWriter(region.analysisRootDir + "/" + region.label + "_" + study
                                       + "_" + factorDir + "_topSNPs.xln");

        writer.println(headerMap.get(study).get(factorDir).toString());

        for (StringBuilder sb : resultsMap.get(study).get(factorDir)) {
          writer.println(sb.toString());
        }

        writer.flush();
        writer.close();
      }
    }
  }


  public static void processOnly(String analysisDir, String inputFile, String dataFile) {
    Region[] rgns = parseSetupFile(inputFile);
    for (Region rgn : rgns) {
      String newDir =
          rgn.label + "_iter1_" + (ext.replaceWithLinuxSafeCharacters(rgn.indexSNP, false)
                                      .replaceAll("_", ""))
                      + "/";
      rgn.analysisRootDir = analysisDir;
      rgn.regionDirNameRoot = newDir;
    }
    HashMap<String, HashMap<String, DataDefinitions>> dataDefsTemp = null;
    try {
      dataDefsTemp = FAST.parseDataDefinitionsFile(dataFile);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    for (Region region : rgns) {
      processAndWriteResults(region, dataDefsTemp);
      processAndWriteResults2(region, dataDefsTemp);
    }

  }

  private boolean[] buildAnalysisFolders(String dir, Region[] rgns) {
    boolean[] results = Array.booleanArray(rgns.length, false);
    for (int i = 0; i < rgns.length; i++) {
      String newDir = rgns[i].label + "_iter0_baseline/";
      results[i] = new File(dir + newDir).mkdirs();
      newDir = rgns[i].label + "_iter1_"
               + (ext.replaceWithLinuxSafeCharacters(rgns[i].indexSNP, false).replaceAll("_", ""))
               + "/";
      results[i] = results[i] && new File(dir + newDir).mkdirs();
      rgns[i].analysisRootDir = dir;
      rgns[i].regionDirNameRoot = newDir;
    }
    return results;
  }

  /*
   * Input file format: Tab-delimited columns: REGION_LABEL INDEX_SNP UCSC_REGION
   */
  private void setup(final String analysisDir, String inputFile, final String dataFile,
                     String tempDataDir, String traitDir) {
    log("Parsing regions from input file...");
    Region[] rgns = parseSetupFile(inputFile);
    log("Found " + rgns.length + " regions for analysis");
    boolean[] dirCreation = buildAnalysisFolders(analysisDir, rgns);
    log("Parsing data file..."); // TODO divorce from FAST?
    HashMap<String, HashMap<String, DataDefinitions>> dataDefsTemp = null;
    try {
      dataDefsTemp = FAST.parseDataDefinitionsFile(dataFile);
    } catch (IOException e1) {
      // TODO Auto-generated catch block
      e1.printStackTrace();
    }
    final HashMap<String, HashMap<String, DataDefinitions>> dataDefs = dataDefsTemp;

    int threads = rgns.length;
    int avail = Runtime.getRuntime().availableProcessors();
    threads = Math.min(threads, avail);
    ExecutorService executor = Executors.newFixedThreadPool(threads);

    // TODO error checking on folder creation
    for (Region rgn : rgns) {
      ConditionalAnalysisToolset_FAST run =
          new ConditionalAnalysisToolset_FAST(rgn, dataFile, traitDir, tempDataDir, dataDefs, true);
      executor.execute(run);
    }

    for (int i = 0; i < rgns.length; i++) {
      if (dirCreation[i]) {
        log("Processing region " + rgns[i].toString());
        ConditionalAnalysisToolset_FAST run =
            new ConditionalAnalysisToolset_FAST(rgns[i], dataFile, traitDir, tempDataDir, dataDefs,
                                                false);
        executor.execute(run);
      }
    }

    executor.shutdown();
    try {
      executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    for (Region region : rgns) {
      log("Writing final results files for region [" + region.label + "]");
      processAndWriteResults(region, dataDefs);
      processAndWriteResults2(region, dataDefs);
    }

  }

}


