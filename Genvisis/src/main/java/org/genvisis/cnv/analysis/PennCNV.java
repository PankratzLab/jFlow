package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;

import org.genvisis.cnv.Resources;
import org.genvisis.cnv.filesys.CNVariant;
import org.genvisis.cnv.filesys.Centroids.CENTROID_STRATEGY;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.hmm.PFB;
import org.genvisis.cnv.qc.GcAdjustor.GcModel;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.cnv.var.SampleData;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.SciStringComparator;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.qsub.Qsub;

import com.google.common.collect.ImmutableSet;

public class PennCNV {

  public static final String[] QC_HEADS = {"LRR_mean", "LRR_median", "LRR_SD", "BAF_mean",
                                           "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"};
  public static final String[] ERRORS = {"large SD for LRR", "drifting BAF values",
                                         "waviness factor values", "Small-sized CNV calls",
                                         "NoCall rate"};
  public static final String QC_SUMMARY_EXTENSION = "_QC.xln";
  public static final int MISSING_SCORE = -1;

  public static void batch(Project proj, int numChunks, List<String> execList, String pfbFile,
                           String gcmodelFile, String hmmFile, String scriptSubDir,
                           String dataSubDir, String resultsSubDir) {
    String commands;
    PrintWriter writer;
    String[] files;
    int step;
    String execDir, dataDir, resultsDir, projDir, pennDir;
    Logger log;

    log = proj.getLog();
    final String runLine = Files.getRunString() + " " + PennCNV.class.getCanonicalName() + " proj="
                           + new File(proj.getPropertyFilename()).getAbsolutePath();
    projDir = proj.PROJECT_DIRECTORY.getValue();
    execDir = proj.PENNCNV_EXECUTABLE_DIRECTORY.getValue(false, true);
    pennDir = proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true);
    dataDir = proj.PENNCNV_DATA_DIRECTORY.getValue(false, true) + dataSubDir;
    resultsDir = pennDir + resultsSubDir;

    String penncnvExecutable = execDir + "detect_cnv.pl";
    if (!Files.exists(penncnvExecutable)) {
      log.reportError("WARNING - couldn't find PennCNV executable 'detect_cnv.pl' in given directory: "
                      + execDir);
      if (Files.programExists("detect_cnv.pl")) {
        log.report("PennCNV executable 'detect_cnv.pl' found on the PATH; please check the PENNCNV_EXECUTABLE_DIRECTORY project property.");
      }
    }

    if (pfbFile != null) {
      pfbFile = ext.replaceTilde(pfbFile);
      if (!pfbFile.startsWith("/") && (pfbFile.charAt(1) != ':')) {
        pfbFile = ext.pwd() + pfbFile;
      }
      if (!Files.exists(pfbFile)) {
        log.reportError("Error - pfb file '" + pfbFile + "' does not exist; aborting");
        return;
      }
    }

    if (gcmodelFile != null) {
      gcmodelFile = ext.replaceTilde(gcmodelFile);
      if (!gcmodelFile.startsWith("/") && (gcmodelFile.charAt(1) != ':')) {
        gcmodelFile = ext.pwd() + gcmodelFile;
      }
      if (!Files.exists(gcmodelFile)) {
        log.reportError("Error - gcmodel file '" + gcmodelFile + "' does not exist; aborting");
        return;
      }
    }

    new File(resultsDir).mkdirs();
    new File(dataDir).mkdirs();

    files = new File(dataDir).list(new FilenameFilter() {

      @Override
      public boolean accept(File file, String filename) {
        return filename.endsWith(".gz");
        // file.length() > 1000 && !filename.equals("chrX") && !filename.equals("sexSpecific");
      }
    });

    if (files == null || files.length == 0) {
      log.reportError("Found zero files in " + dataDir);
      log.reportError("Will not proceed");
      return;
    }
    log.report("Found " + files.length + " files in " + dataDir);

    step = (int) Math.ceil((double) files.length / (double) numChunks);
    log.report("Which means the step for " + numChunks + " chunks would be " + step);

    for (int i = 0; i < numChunks; i++) {
      try {
        writer = Files.openAppropriateWriter(resultsDir + "list" + (i + 1) + ".txt");
        for (int j = i * step; j < Math.min(files.length, (i + 1) * step); j++) {
          if (files[j].endsWith(".gz")) {
            writer.println("`gunzip -c " + dataDir + files[j] + "`");
          } else {
            writer.println(dataDir + files[j]);
          }
        }
        writer.close();
      } catch (Exception e) {
        log.reportError("Error writing to " + resultsSubDir + "list" + (i + 1) + ".txt");
        log.reportException(e);
      }
    }

    commands = "/bin/bash -c \"" + execDir + "detect_cnv.pl -test -conf -hmm "
               + (hmmFile == null ? execDir + "lib/hhall.hmm" : hmmFile) + " -pfb "
               + (pfbFile == null ? execDir + "lib/hhall.hg18.pfb" : pfbFile) + " -gcmodel "
               + (gcmodelFile == null ? execDir + "lib/hhall.hg18.gcmodel" : gcmodelFile)
               + " -list " + resultsDir + "list[%0].txt -log " + resultsDir + "[%0].log -out "
               + resultsDir + "[%0].rawcnv > " + resultsDir + "[%0].out\"";

    new File(pennDir + scriptSubDir).mkdirs();

    if (execList == null) {
      Qsub.qsub(pennDir + scriptSubDir + "runPenn", dataDir, numChunks, commands,
                Matrix.toMatrix(ArrayUtils.stringArraySequence(numChunks, "")), 2200, 16);
    } else {
      Files.execListAdd(execList, commands, ArrayUtils.stringArraySequence(numChunks, ""), log);
    }

    Files.writeArray(new String[] {"cd " + projDir, "cat " + resultsDir + "*.log > " + resultsDir
                                                    + "penncnv.rawlog",
                                   "cat " + resultsDir + "*.rawcnv > " + resultsDir + "penncnv.rawcnv",
                                   runLine + " rawlog=" + resultsDir + "penncnv.rawlog",
                                   runLine + " rawcnv=" + resultsDir + "penncnv.rawcnv",},
                     pennDir + scriptSubDir + "assemblePenncnv");
    Files.chmod(pennDir + scriptSubDir + "assemblePenncnv");
  }

  // FIXME need to unify this method with batch
  public static void batchX(Project proj, int numChunks, List<String> execList, String pfbFile,
                            String gcmodelFile, String hmmFile, String scriptSubDir,
                            String dataSubDir, String resultsSubDir) {
    String commands;
    PrintWriter writer;
    String[] files;
    int step;
    String execDir, dataDir, resultsDir, projDir, pennDir;
    Logger log;

    log = proj.getLog();
    final String runLine = Files.getRunString() + " " + PennCNV.class.getCanonicalName() + " proj="
                           + new File(proj.getPropertyFilename()).getAbsolutePath();
    projDir = proj.PROJECT_DIRECTORY.getValue();
    execDir = proj.PENNCNV_EXECUTABLE_DIRECTORY.getValue(false, true);
    pennDir = proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true);
    dataDir = proj.PENNCNV_DATA_DIRECTORY.getValue(false, true) + dataSubDir;
    resultsDir = pennDir + resultsSubDir;

    // if (!Files.exists(proj.getFilename("SAMPLE_DATA_FILENAME", false, false),
    // proj.getJarStatus())) {
    if (!Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false))) {
      log.reportError("Error - sample data file " + proj.SAMPLE_DATA_FILENAME.getValue()
                      + " does not exist;");
      return;
    }
    SampleData sampleData;
    try {
      sampleData = proj.getSampleData(false);
    } catch (Exception e) {
      log.reportError("Error - without a sample data file, PennCNV will fail to analyze sex chromosomes");
      return;
    }

    if (pfbFile != null) {
      pfbFile = ext.replaceTilde(pfbFile);
      if (!pfbFile.startsWith("/") && (pfbFile.charAt(1) != ':')) {
        pfbFile = ext.pwd() + pfbFile;
      }
      if (!Files.exists(pfbFile)) {
        log.reportError("Error - pfb file '" + pfbFile + "' does not exist; aborting");
        return;
      }
    }

    if (gcmodelFile != null) {
      gcmodelFile = ext.replaceTilde(gcmodelFile);
      if (!gcmodelFile.startsWith("/") && (gcmodelFile.charAt(1) != ':')) {
        gcmodelFile = ext.pwd() + gcmodelFile;
      }
      if (!Files.exists(gcmodelFile)) {
        log.reportError("Error - gcmodel file '" + gcmodelFile + "' does not exist; aborting");
        return;
      }
    }

    new File(resultsDir).mkdirs();
    new File(dataDir).mkdirs();

    String newGCFile = dataDir + "chrX.gcModel";
    String newPFBFile = dataDir + "chrX.pfb";

    pfbFile = filterFile(proj, pfbFile, newPFBFile, new String[] {"23", "X"});
    gcmodelFile = filterFile(proj, gcmodelFile, newGCFile, new String[] {"23", "X"});

    String sexFileStatus = writeSexFile(proj, sampleData, dataDir, log);
    if (sexFileStatus.length() > 0) {
      log.reportError("Error - " + sexFileStatus);
      return;
    }

    files = new File(dataDir).list(new FilenameFilter() {

      @Override
      public boolean accept(File file, String filename) {
        return filename.endsWith(".gz");
        // return file.length() > 1000 && !filename.endsWith(".pfb") &&
        // !filename.endsWith(".gcmodel") && !filename.startsWith("sex_file");
      }
    });
    log.report("Found " + files.length + " files");

    if (files == null || files.length == 0) {
      log.reportError("Found zero files in " + dataDir);
      log.reportError("Will not proceed");
      return;
    }
    step = (int) Math.ceil((double) files.length / (double) numChunks);
    log.report("Which means the step for " + numChunks + " chunks would be " + step);

    for (int i = 0; i < numChunks; i++) {
      try {
        writer = Files.openAppropriateWriter(resultsDir + "list" + (i + 1) + ".txt");
        for (int j = i * step; j < Math.min(files.length, (i + 1) * step); j++) {
          if (files[j].endsWith(".gz")) {
            writer.println("`gunzip -c " + dataDir + files[j] + "`");
          } else {
            writer.println(files[j]);
          }
        }
        writer.close();
      } catch (Exception e) {
        log.reportError("Error writing to " + resultsSubDir + "list" + (i + 1) + ".txt");
        log.reportException(e);
      }
    }

    commands = "/bin/bash -c \"" + execDir + "detect_cnv.pl -test -conf -hmm "
               + (hmmFile == null ? execDir + "lib/hhall.hmm" : hmmFile) + " -pfb "
               + (pfbFile == null ? execDir + "lib/hhall.hg18.pfb" : pfbFile) + " -gcmodel "
               + (gcmodelFile == null ? execDir + "lib/hhall.hg18.gcmodel" : gcmodelFile)
               + " -chrx -sexfile " + dataDir + "sex_file.txt -list " + resultsDir
               + "list[%0].txt -log " + resultsDir + "[%0].log -out " + resultsDir
               + "[%0].rawcnv > " + resultsDir + "[%0].out\"";

    new File(pennDir + scriptSubDir).mkdirs();

    if (execList == null) {
      Qsub.qsub(pennDir + scriptSubDir + "runPennX", dataDir, numChunks, commands,
                Matrix.toMatrix(ArrayUtils.stringArraySequence(numChunks, "")), 2200, 16);
    } else {
      Files.execListAdd(execList, commands, ArrayUtils.stringArraySequence(numChunks, ""), log);
    }
    Files.writeArray(new String[] {"cd " + projDir,
                                   "cat " + resultsDir + "*.log > " + resultsDir
                                                    + "penncnvX.rawlog",
                                   "cat " + resultsDir + "*.rawcnv > " + resultsDir
                                                                         + "penncnvX.rawcnv",
                                   // don't parse warnings; the parseWarnings method isn't written
                                   // to
                                   // parse X-chromosome warnings
                                   runLine + " rawcnv=" + resultsDir + "penncnvX.rawcnv",},
                     pennDir + scriptSubDir + "assemblePenncnv");
    Files.chmod(pennDir + scriptSubDir + "assemblePenncnv");
  }

  private static String filterFile(Project proj, String fileToFilter, String outputFile,
                                   String[] chrs) {
    // TODO combine method with filterPFB - literally the same except different names/extensions
    BufferedReader reader = null;
    PrintWriter writer = null;

    try {
      reader = new BufferedReader(new FileReader(fileToFilter));
      writer = Files.openAppropriateWriter(outputFile);

      String header;
      String temp;
      String[] line;
      if (reader.ready()) {
        header = reader.readLine();
        writer.println(header);
      }
      while (reader.ready()) {
        temp = reader.readLine();
        line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        for (String chr : chrs) {
          if (line[1].equals(chr)) {
            writer.println(temp);
          }
        }
      }
    } catch (IOException e) {
      proj.getLog().reportError("Error - filtering failed for file: " + fileToFilter);
      proj.getLog().reportException(e);
      return fileToFilter;
    } finally {
      if (reader != null) {
        try {
          reader.close();
        } catch (IOException e) {
          proj.getLog()
              .reportError("Error - couldn't properly close file reader for " + fileToFilter);
          proj.getLog().reportException(e);
        }
        reader = null;
      }
      if (writer != null) {
        writer.close();
        writer = null;
      }
    }

    return outputFile;
  }

  private static String writeSexFile(Project proj, SampleData sampleData, String resultsDir,
                                     Logger log) {
    String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
    String[] header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
    int sexInd = -1;
    for (int i = 0; i < header.length; i++) {
      if (("CLASS=" + SexChecks.EST_SEX_HEADER).equalsIgnoreCase(header[i])) {
        sexInd = i;
        break;
      }
    }
    if (sexInd == -1) {
      return "no estimated sex found in sample data file - please run SexChecks with -check argument to generate the required data";
    }
    Hashtable<String, Vector<String>> sexData = HashVec.loadFileToHashVec(sampleDataFile, 0,
                                                                          new int[] {sexInd}, "\t",
                                                                          true, false);
    try {
      PrintWriter writer = Files.openAppropriateWriter(resultsDir + "sex_file.txt");
      for (Map.Entry<String, Vector<String>> lineData : sexData.entrySet()) {
        String estSexStr = lineData.getValue().get(0);
        if (!ext.isMissingValue(estSexStr)) {
          int estSex = Integer.parseInt(estSexStr);
          estSex = SexChecks.EstimatedSex.values()[estSex].getKaryotype().contains("XX") ? 2 : 1;
          writer.println(lineData.getKey() + "\t" + estSex);
        }
      }
      writer.close();
    } catch (IOException e) {
      log.reportException(e);
      return "unable to complete writing of sex_file for PennCNV";
    }
    return "";
  }

  public static void parseWarnings(Project proj, String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, data;
    String temp, sampleID = null;
    Hashtable<String, String[]> hash = new Hashtable<>();
    Vector<String> v = new Vector<>();
    SampleData sampleData;
    int err;
    double lrrSD_cutoff;
    String[] ids;
    long time;
    Logger log;

    time = new Date().getTime();
    log = proj.getLog();
    log.report("Parsing PennCNV warning...");

    sampleData = proj.getSampleData(false);
    // lrrSD_cutoff = proj.getDouble(proj.LRRSD_CUTOFF);
    lrrSD_cutoff = proj.LRRSD_CUTOFF.getValue();

    try {
      reader = new BufferedReader(new FileReader(filename));
      while (reader.ready()) {
        temp = reader.readLine();
        temp = translateDerivedSamples(temp, log);
        line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        try {
          if (temp.contains("quality summary")) {
            sampleID = line[4].substring(line[4].lastIndexOf("/") + 1, line[4].indexOf(":"));
            v.add(sampleID);
            data = new String[QC_HEADS.length + ERRORS.length];
            for (int i = 0; i < ERRORS.length; i++) {
              data[i] = ".";
            }
            if (line.length < QC_HEADS.length + 5) {
              log.reportError("Error - line doesn't have all the expected pieces:");
              log.reportError(temp);
            }
            for (int i = 0; i < QC_HEADS.length; i++) {
              data[ERRORS.length + i] = line[5 + i].split("=")[1];
            }
            hash.put(sampleID, data);
          } else if (temp.startsWith("WARNING")) {
            if (temp.contains("Small-sized CNV calls")) {
              // use old trav
            } else {
              sampleID = line[3].substring(line[3].lastIndexOf("/") + 1);
            }
            data = hash.get(sampleID);
            err = -1;
            for (int i = 0; i < ERRORS.length; i++) {
              if (temp.contains(ERRORS[i])) {
                err = i;
              }
            }
            if (err == -1) {
              log.reportError("Unknown WARNING: " + temp);
            } else {
              data[err] = "1";
            }
          }
        } catch (Exception e) {
          log.reportError("Error with: " + temp);
          log.reportException(e);
        }
      }
      reader.close();

      writer = Files.openAppropriateWriter(proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(filename)
                                           + QC_SUMMARY_EXTENSION);
      writer.print("Sample\tFID\tIID\tUse_" + ext.formDeci(lrrSD_cutoff, 2));
      for (String element : ERRORS) {
        writer.print("\t" + element);
      }
      for (String element : QC_HEADS) {
        writer.print("\t" + element);
      }
      writer.println();
      Collections.sort(v);
      for (int i = 0; i < v.size(); i++) {
        sampleID = v.get(i);
        data = hash.get(sampleID);
        ids = sampleData.lookup(sampleID);
        writer.print(sampleID + "\t" + (ids == null ? "NotInSampleData\t" + sampleID : ids[1]));
        writer.print("\t" + (data[1].equals("1") || data[2].equals("1")
                             || Double.parseDouble(data[6]) > lrrSD_cutoff ? "0" : "1"));
        writer.println("\t" + ArrayUtils.toStr(data));
      }
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + filename + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + filename + "\"");
      return;
    }
    log.report("Parsed PennCNV warnings in " + ext.getTimeElapsed(time));
  }

  public static String translateDerivedSamples(String str, Logger log) {
    String trav;
    int start, stop;

    start = str.indexOf("`");
    stop = str.lastIndexOf("`");
    if (start == -1) {
      return str;
    }

    trav = str.substring(start + 1, stop);
    if (trav.contains("`")) {
      log.reportError("Error - more than one set of quotes for: " + str);
    }
    if (trav.startsWith("gunzip -c") && trav.endsWith(".gz")) {
      trav = trav.substring(9, trav.length() - 3).trim();
    } else {
      log.reportError("Error - not currently set up to handle the following construction into a sample_ID: "
                      + trav);
    }

    return str.substring(0, start) + trav + str.substring(stop + 1);
  }

  public static void combineResults(Project proj, String[] cnvFiles, String outputFile,
                                    boolean recode, boolean removeChr11) {
    BufferedReader reader;
    PrintWriter writer;
    Logger log = proj.getLog();

    // TODO check input and output file names for .cnv extension( - error if not? or just
    // warning...?)

    java.util.HashMap<String, java.util.TreeMap<String, java.util.ArrayList<String[]>>> cnvSet = new HashMap<>();

    String temp;
    String[] line;
    String key, chr, currFile = null;
    boolean readAll = false;
    try {
      for (String cnvFile : cnvFiles) {
        currFile = cnvFile;
        reader = new BufferedReader(new FileReader(cnvFile));
        if (reader.ready()) {
          // skip header
          reader.readLine();
        }
        while (reader.ready()) {
          temp = reader.readLine();
          line = temp.split("\t");
          key = line[0] + "\t" + line[1];
          // get all CNVs for an individual:
          TreeMap<String, ArrayList<String[]>> chrSets = cnvSet.get(key);
          if (chrSets == null) {
            chrSets = new TreeMap<>();
            cnvSet.put(key, chrSets);
          }
          chr = line[2];
          // get all CNVs for a specific chromosome
          ArrayList<String[]> chrSet = chrSets.get(chr);
          if (chrSet == null) {
            chrSet = new ArrayList<String[]>() {

              private static final long serialVersionUID = 1L;

              @Override
              public boolean add(String[] e) {
                int index = ArrayUtils.binarySearch(this, e, 0, false);
                super.add(index, e);
                return true;
              }
            };
            chrSets.put(chr, chrSet);
          }
          // add CNV to list
          chrSet.add(ArrayUtils.subArray(line, 3));
        }
        reader.close();
      }
      readAll = true;
    } catch (FileNotFoundException e) {
      log.reportError("Error: file \"" + currFile + "\" not found in current directory");
      log.reportException(e);
    } catch (IOException e) {
      log.reportException(e);
    }

    if (readAll) {
      try {
        writer = Files.openAppropriateWriter(outputFile);
        writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER));
        String FIDIID;
        String cnvChr;
        for (Map.Entry<String, TreeMap<String, ArrayList<String[]>>> sample : cnvSet.entrySet()) {
          FIDIID = sample.getKey();
          for (Map.Entry<String, ArrayList<String[]>> chrSet : sample.getValue().entrySet()) {
            cnvChr = chrSet.getKey();
            if (removeChr11 && "11".equals(cnvChr)) {
              continue;
            }
            if (recode) {
              if ("1".equals(cnvChr)) {
                cnvChr = "23";
              } else if ("2".equals(cnvChr)) {
                cnvChr = "24";
              } else if ("3".equals(cnvChr)) {
                cnvChr = "25";
              } else if ("4".equals(cnvChr)) {
                cnvChr = "26";
              }
            }
            for (String[] cnv : chrSet.getValue()) {
              writer.println(FIDIID + "\t" + cnvChr + "\t" + ArrayUtils.toStr(cnv, "\t"));
            }
          }
        }
        writer.close();
      } catch (IOException e) {
        log.reportException(e);
      }
    }
  }

  public static void parseResults(Project proj, String filename, boolean denovoOnly) {
    BufferedReader reader;
    PrintWriter writer;
    PrintWriter denovoWriter = null;
    String[] line;
    String temp, trav;
    Vector<String> warnings;
    Hashtable<String, Vector<String>> pedinfo;
    int[] position;
    String score;
    SampleData sampleData;
    String famIndPair;
    Hashtable<String, String> hash;
    String[] ids;
    List<String> inds;
    String[] fams;
    long time;
    int sex;
    Logger log;

    log = proj.getLog();
    log.report("Parsing PennCNV rawcnvs...");
    time = new Date().getTime();

    if (!Files.exists(filename)) {
      log.reportError("Error - could not find file '" + filename + "'");
      return;
    }

    warnings = new Vector<>();
    sampleData = proj.getSampleData(false);
    pedinfo = new Hashtable<>();
    Pedigree ped = proj.loadPedigree();
    PrintWriter[] denoValWriter = new PrintWriter[1];
    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = Files.openAppropriateWriter(ext.rootOf(filename, false) + ".cnv");
      writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER));
      hash = new Hashtable<>();
      while (reader.ready()) {
        temp = reader.readLine();
        if (!temp.startsWith("NOTICE:")) {
          temp = translateDerivedSamples(temp, log);
          line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
          position = Positions.parseUCSClocation(line[0]);
          trav = line[4];
          trav = trav.substring(trav.lastIndexOf("/") + 1);
          ids = sampleData.lookup(trav);
          if (ids == null) {
            if (!hash.containsKey(trav)) {
              // log.reportError("Error - '"+trav+"' was not found in
              // "+proj.getFilename(proj.SAMPLE_DATA_FILENAME));
              log.reportError("Error - '" + trav + "' was not found in "
                              + proj.SAMPLE_DATA_FILENAME.getValue());
              hash.put(trav, "");
            }
            famIndPair = trav + "\t" + trav;
          } else {
            famIndPair = ids[1];
          }

          ids = famIndPair.split("\t");
          HashVec.addToHashVec(pedinfo, ids[0], ids[1], true);

          if (line.length < 8 || !line[7].startsWith("conf=")
              || line[7].toUpperCase().contains("NAN")) {
            score = Integer.toString(MISSING_SCORE);
            if (!warnings.contains(trav) && warnings.size() < 10) {
              log.reportError("Warning - no conf estimates for " + trav);
              warnings.add(trav);
            }
          } else {
            score = ext.formDeci(Double.parseDouble(line[7].substring(5)), 4, true);
          }
          boolean isDenovo = false;
          for (String s : line) {
            if (s.startsWith("statepath=33") || s.startsWith("triostate=33")) {
              isDenovo = true;
            }
          }

          String copynum = line[3].substring(line[3].indexOf("=") + 1);
          String sites = line[1].substring(7);
          StringBuilder lineOut = new StringBuilder(famIndPair).append("\t").append(position[0])
                                                               .append("\t").append(position[1])
                                                               .append("\t").append(position[2])
                                                               .append("\t").append(copynum)
                                                               .append("\t").append(score)
                                                               .append("\t").append(sites);
          if (!denovoOnly) {
            writer.println(lineOut.toString());
          }
          if (isDenovo) {
            if (denovoWriter == null) {
              denovoWriter = Files.openAppropriateWriter(ext.rootOf(filename, false)
                                                         + "_denovo.cnv");
              denovoWriter.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER));
            }
            denovoWriter.println(lineOut.toString());
            writeValidation(ped, ids, position, copynum, line, filename, denoValWriter, log);
          }
        }
      }
      reader.close();
      writer.close();
      if (denovoWriter != null) {
        denovoWriter.close();
      }
      if (denoValWriter[0] != null) {
        denoValWriter[0].close();
      }

      // FilterCalls.stdFilters(dir, ext.rootOf(filename)+".cnv", MAKE_UCSC_TRACKS);

      writer = Files.openAppropriateWriter(ext.rootOf(filename, false) + ".fam");
      fams = HashVec.getNumericKeys(pedinfo);
      for (String fam : fams) {
        inds = pedinfo.get(fam);
        Collections.sort(inds, new SciStringComparator());
        for (String ind : inds) {
          ids = sampleData.lookup(fam + "\t" + ind);
          if (ids != null) {
            sex = sampleData.getSexForIndividual(ids[0]);
          } else {
            sex = 0;
          }
          int pedIndex = ped == null ? -1 : ped.getIndIndex(fam, ind);
          String fa = pedIndex >= 0 && ped != null ? ped.getFA(pedIndex) : "0";
          String mo = pedIndex >= 0 && ped != null ? ped.getMO(pedIndex) : "0";
          writer.println(fam + "\t" + ind + "\t" + fa + "\t" + mo + "\t" + Math.max(0, sex)
                         + "\t-9");
        }
      }
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + ext.rootOf(filename, false)
                      + ".cnv\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + ext.rootOf(filename, false) + "\"");
      return;
    }

    log.report("...finished in " + ext.getTimeElapsed(time));
  }

  private static void writeValidation(Pedigree ped, String[] ids, int[] position, String copynum,
                                      String[] line, String filename, PrintWriter[] denoValWriter,
                                      Logger log) {
    int pedIndex = ped.getIndIndex(ids[0], ids[1]);
    if (pedIndex < 0) {
      return;
    }
    int faIndex = ped.getIndexOfFaInIDs(pedIndex);
    int moIndex = ped.getIndexOfMoInIDs(pedIndex);
    if (faIndex < 0 || moIndex < 0) {
      return;
    }

    String cDna = ped.getiDNA(pedIndex);
    String faDna = ped.getiDNA(faIndex);
    String moDna = ped.getiDNA(moIndex);

    String outDir = ext.parseDirectoryOfFile(filename);

    if (denoValWriter[0] == null) {
      try {
        denoValWriter[0] = Files.openAppropriateWriter(outDir + "denovoValidation.txt");
        denoValWriter[0].println("export HMMFILE=");
        denoValWriter[0].println("export PFBFILE=");
        denoValWriter[0].println();
      } catch (IOException e) {
        log.reportException(e);
        return;
      }
    }

    String[] bounds = getBounds(line);
    if (bounds[0] == null || bounds[1] == null) {
      return;
    }

    outDir += "denovo" + File.separator;
    File outFile = new File(outDir);
    if (!outFile.exists()) {
      outFile.mkdirs();
    }

    String childSource = "gunzip -c " + line[4] + ".gz";
    String faSource = childSource.replace(cDna, faDna);
    String moSource = childSource.replace(cDna, moDna);

    if (childSource.contains("sexSpecific")) {
      faSource = faSource.replaceAll("/female/", "/male/");
      moSource = moSource.replaceAll("/male/", "/female/");
    }

    String out = outDir + ids[0] + "_" + ids[1] + "_" + position[0] + "_" + position[1] + "_"
                 + position[2];
    String faFile = out + "_fa.txt";
    String moFile = out + "_mo.txt";
    String childFile = out + "_off.txt";

    StringBuilder extractLine = new StringBuilder(faSource).append(" > ").append(faFile)
                                                           .append(" && ").append(moSource)
                                                           .append(" > ").append(moFile)
                                                           .append(" && ").append(childSource)
                                                           .append(" > ").append(childFile);

    denoValWriter[0].println(extractLine.toString());

    StringBuilder sb = new StringBuilder("/home/pankrat2/shared/bin/infer_snp_allele.pl -pfbfile $PFBFILE -hmmfile $HMMFILE").append(" -denovocn ")
                                                                                                                             .append(copynum)
                                                                                                                             .append(" -startsnp ")
                                                                                                                             .append(bounds[0])
                                                                                                                             .append(" -endsnp ")
                                                                                                                             .append(bounds[1])
                                                                                                                             .append(" -outfile ")
                                                                                                                             .append(out)
                                                                                                                             .append(".gen  -logfile ")
                                                                                                                             .append(out)
                                                                                                                             .append(".log ")
                                                                                                                             .append(faFile)
                                                                                                                             .append(" ")
                                                                                                                             .append(moFile)
                                                                                                                             .append(" ")
                                                                                                                             .append(childFile);

    denoValWriter[0].println(sb.toString());
    StringBuilder cleanup = new StringBuilder("rm ").append(faFile).append(" && rm ").append(moFile)
                                                    .append(" && rm ").append(childFile);

    denoValWriter[0].println(cleanup.toString());
  }

  private static String[] getBounds(String[] line) {
    String[] bounds = new String[2];
    for (String s : line) {
      if (s.startsWith("startsnp")) {
        bounds[0] = s.split("=")[1];
      } else if (s.startsWith("endsnp")) {
        bounds[1] = s.split("=")[1];
      }
    }
    return bounds;
  }

  // Available in cnv.Launch
  // public static void fromParameters(String filename, Logger log) {
  // Vector<String> params;
  //
  // params = Files.parseControlFile(filename, "penncnv", new String[]
  // {"proj=/home/npankrat/projects/GEDI.properties", "rawcnv=all.rawcnv", "rawlog=all.log"}, log);
  //
  // if (params != null) {
  // params.add("log=" + log.getFilename());
  // main(Array.toStringArray(params));
  // }
  // }

  private static String[] getSamplesForTransform(Project proj, boolean excludeExcludeds) {
    if (excludeExcludeds) {
      return ArrayUtils.subArray(proj.getSamples(), proj.getSamplesToInclude(null));
    } else {
      return proj.getSamples();
    }
  }

  public static void doBatch(Project proj, boolean auto, boolean chrx, boolean sexCent,
                             boolean transformData, int numChunks, boolean separateQsubFiles,
                             String pfbFile, String gcmodelFile, String hmmFile,
                             boolean submitImmed, boolean createCombined, boolean useExcludes,
                             int threadCount) {
    boolean problem = false;
    Vector<String> execList;
    Logger log = proj.getLog();
    final String runLine = Files.getRunString() + " " + PennCNV.class.getCanonicalName() + " proj="
                           + new File(proj.getPropertyFilename()).getAbsolutePath();

    String dir = proj.PENNCNV_RESULTS_DIRECTORY.getValue();
    dir += "penn_scripts/";

    if (hmmFile == null || !Files.exists(hmmFile)) {
      hmmFile = Resources.cnv(log).getAllHmm().get();
    }

    if ((pfbFile == null || !Files.exists(pfbFile))
        && (pfbFile = proj.CUSTOM_PFB_FILENAME.getValue()) == null) {
      System.err.println("Error - could not find " + pfbFile);
      problem = true;
    }
    if ((gcmodelFile == null || !Files.exists(gcmodelFile))
        && (gcmodelFile = proj.GC_MODEL_FILENAME.getValue()) == null) {
      System.err.println("Error - could not find " + gcmodelFile);
      problem = true;
    }
    if (problem) {
      return;
    }

    if (separateQsubFiles) {
      execList = null;
    } else {
      execList = new Vector<>();
    }

    String[] samples = getSamplesForTransform(proj, !useExcludes);

    if (auto) {
      if (transformData) {
        log.report("Transforming data for autosomal CNV analysis");
        AnalysisFormats.exportPenncnvSamples(proj, samples, null, null,
                                             Runtime.getRuntime().availableProcessors());
      }
      log.report("Creating batch scripts for autosomal CNV analysis");
      batch(proj, numChunks, execList, pfbFile, gcmodelFile, hmmFile, "penn_scripts/", "", "");
    }
    if (chrx) {
      MarkerDetailSet ms = proj.getMarkerSet();
      if (ms == null) {
        log.reportError("Error - no marker set available.");
      } else {
        log.report("Transforming data for chromosomal CNV analysis");
        Set<String> xMarkers = ms.getChrMap().get((byte) 23).stream().map(Marker::getName)
                                 .collect(ImmutableSet.toImmutableSet());
        AnalysisFormats.exportPenncnvSamples(proj, samples, xMarkers, "chrX/",
                                             Runtime.getRuntime().availableProcessors());
      }
      log.report("Creating batch scripts for chromosomal CNV analysis");
      batchX(proj, numChunks, execList, pfbFile, gcmodelFile, hmmFile, "penn_scripts/chrX/",
             "chrX/", "chrX/");
    }
    if ((auto && chrx) || (auto && createCombined) || (chrx && createCombined)) {
      // write combine script
      String resultsDir = proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true);
      String outdir = resultsDir + "penn_scripts/";
      new File(outdir).mkdirs();
      String outfile = "combineAutoXCNVs";
      Files.writeArray(new String[] {"cd " + resultsDir,
                                     runLine + " combine=penncnv.cnv,chrX/penncnvX.cnv output=combinedAX.cnv",},
                       outdir + outfile);
      Files.chmod(outdir + outfile);
    }
    if (sexCent) {
      log.report("Transforming data for 'faked' chromosomal CNV analysis");
      // [males.pfb, females.pfb, sexSpecific.gcModel]

      String[] files = AnalysisFormats.pennCNVSexHackMultiThreaded(proj, gcmodelFile,
                                                                   CENTROID_STRATEGY.USE_CENT_IF_EXISTS_OTHERWISE_COMPUTE,
                                                                   useExcludes, threadCount);

      log.report("Creating batch scripts for 'faked' chromosomal CNV analysis");
      String scriptDir = "penn_scripts/sexSpecific/";
      batch(proj, numChunks, execList, files[0], files[2], hmmFile, scriptDir + "male/",
            "sexSpecific/male/", "sexSpecific/male/");
      batch(proj, numChunks, execList, files[1], files[2], hmmFile, scriptDir + "female/",
            "sexSpecific/female/", "sexSpecific/female/");
      // write combine script
      String resultsDir = proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true);
      String outdir = resultsDir + "penn_scripts/";
      String outfile = "combineMFCNVs";
      Files.writeArray(new String[] {"cd " + resultsDir,
                                     runLine + " combine=sexSpecific/male/penncnv.cnv output=sexSpecific/male/recodedM.cnv -recode",
                                     runLine + " combine=sexSpecific/female/penncnv.cnv output=sexSpecific/female/recodedF.cnv -recode",
                                     runLine + " combine=sexSpecific/male/recodedM.cnv,sexSpecific/female/recodedF.cnv output=combinedMF.cnv -recode",},
                       outdir + outfile);
      Files.chmod(outdir + outfile);

      if (auto) {
        outfile = "combineAMFCNVs";
        Files.writeArray(new String[] {"cd " + resultsDir,
                                       runLine + " combine=sexSpecific/male/penncnv.cnv output=sexSpecific/male/recodedM.cnv -recode",
                                       runLine + " combine=sexSpecific/female/penncnv.cnv output=sexSpecific/female/recodedF.cnv -recode",
                                       runLine + " combine=penncnv.cnv,sexSpecific/male/recodedM.cnv,sexSpecific/female/recodedF.cnv output=combinedMF.cnv",},
                         outdir + outfile);
        Files.chmod(outdir + outfile);
      }

    }

    if (execList != null) {
      Qsub.qsubExecutor(proj.PROJECT_DIRECTORY.getValue(), execList, null,
                        proj.PENNCNV_RESULTS_DIRECTORY.getValue() + "runAllPenncnv", 24, 5000, 8);
      log.report("All PennCNV files and scripts have been prepped. The next thing would be to qsub "
                 + proj.PENNCNV_RESULTS_DIRECTORY.getValue() + "runAllPenncnv.pbs");
    }
    List<String> toRun = new ArrayList<>();
    toRun.add(dir + "assemblePenncnv");
    toRun.add(dir + "chrX/assemblePenncnv");
    if (sexCent) {
      toRun.add(dir + "sexSpecific/female/assemblePenncnv");
      toRun.add(dir + "sexSpecific/male/assemblePenncnv");
      toRun.add(dir + "combineAMFCNVs");
    }
    Files.writeArray(toRun.toArray(new String[toRun.size()]), dir + "parseAllPenncnv");
    Files.chmod(dir + "parseAllPenncnv");

    log.report("Script generation complete. See: " + dir);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    String rawlog = null;
    String rawcnvs = null;
    int numChunks = 0;
    boolean transformData = true;
    boolean separateQsubs = false;
    boolean auto = true;
    boolean chrx = true;
    boolean sexCent = true;
    Project proj;
    String pfbFile = null;
    String gcmodelFile = null;
    boolean denovoOnly = false;
    boolean parsePFB = false;
    String gc5base = null;
    String logfile = null;
    String[] cnvFiles = null;
    String outputFile = null;
    boolean recode = false;
    boolean removeChr11 = true;
    boolean submit = false;
    boolean excludes = false;
    String hmmFile = null;
    int numThreads = 1;

    String usage = "\n" + "org.genvisis.cnv.analysis.PennCNV requires 0-1 arguments\n"
                   + "   (0) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + " AND\n" + "   (1) number of chunks to split everything in to (i.e. chunks="
                   + numChunks + " (default))\n"
                   + "   (2) generate seperate qsub files instead of a single executor chain (i.e. -sepqsub (not the default))\n"
                   + "   (3) generate PennCNV scripts to analyze autosomes (i.e. auto=TRUE (default))\n"
                   + "   (4) generate PennCNV scripts to analyze X Chromosome (i.e. chrx=TRUE (default))\n"
                   + "   (5) recompute centroids of chr23-26 (X, Y, XY, MT) and recode as chr1-4 in subdirectory (i.e. sexSpecificCentroids=TRUE (default))\n"
                   + "   (6) transform sample data into PennCNV data files (i.e. data=TRUE (default))\n"
                   + "   (7) number of threads to use (i.e. threads=" + numThreads + " (default))\n"
                   + "   (8) (optional) use custom pfb file (i.e. pfb=custom.pfb (not the default))\n"
                   + "   (9) (optional) use custom gcmodel file (i.e. gcmodel=custom.gcmodel (not the default))\n"
                   + "   (10) (optional) use an array specific hmm file (i.e. hmm= (no default))\n"
                   +

                   " OR\n"
                   + "   (1) compute file containing project based b allele frequencies for file using parameters in properties file (i.e. -pfb (not the default))\n"
                   + " OR\n"
                   + "   (1) compute a custom gcmodel file for the markers in this project using this file (i.e. gc5base=gc5base.txt (not the default))\n"
                   + " OR\n"
                   + "   (1) parse warnings from log file (i.e. rawlog=final.log (not the default))\n"
                   + " OR\n"
                   + "   (1) raw cnvs to parse (i.e. rawcnv=final.rawcnv (not the default))\n"
                   + "   (2) (optional) parse only de novo variants (i.e. -denovoOnly (not the default))\n"
                   + " OR\n"
                   + "   (1) a comma-separated list of .cnv files to combine together (i.e. combine=/full/path/to/cnv1.cnv,relative/path/to/cnv2.cnv (not the default))\n"
                   + "   (2) full path of the desired output file (i.e. output=/path/to/output/file.cnv (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        return;
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("chunks=")) {
        numChunks = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("-sepqsub")) {
        separateQsubs = true;
        numArgs--;
      } else if (arg.startsWith("rawlog=")) {
        rawlog = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("rawcnv=")) {
        rawcnvs = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-pfb")) {
        parsePFB = true;
        numArgs--;
      } else if (arg.startsWith("gc5base=")) {
        gc5base = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-denovoOnly")) {
        denovoOnly = true;
        numArgs--;
      } else if (arg.startsWith("pfb=")) {
        pfbFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("gcmodel=")) {
        gcmodelFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("hmm=")) {
        hmmFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("data=")) {
        transformData = Boolean.parseBoolean(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("auto=")) {
        auto = Boolean.parseBoolean(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("chrx=")) {
        chrx = Boolean.parseBoolean(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("sexSpecificCentroids=")) {
        sexCent = Boolean.parseBoolean(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("combine=")) {
        cnvFiles = arg.split("=")[1].split(",");
        numArgs--;
      } else if (arg.startsWith("output=")) {
        outputFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-recode")) {
        recode = true;
        numArgs--;
      } else if (arg.startsWith("-submit")) {
        submit = true;
        numArgs--;
      } else if (arg.startsWith("-useExcluded")) {
        excludes = true;
        numArgs--;
      } else if (arg.startsWith("threads=")) {
        numThreads = ext.parseIntArg(arg);
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      return;
    }
    try {

      // filename = "C:/workspace/Genvisis/projects/GEDI_exome.properties";
      // parsePFB = true;
      // gc5base = "C:/projects/gcModel/gc5Base.txt";

      // logfile = "penncnv/penncnv.log";
      // rawcnvs = "penncnv/penncnv.rawcnv";

      // filename = "/home/npankrat/projects/GEDI.properties";
      // batch = 60;
      // qsub = true;
      // pfbFile = "gedi.pfb";
      // gcmodelFile = "gedi.gcmodel";
      //
      // batch = 1;
      // filename = "C:/data/FarrarMike/default.properties";
      // qsub = true;
      // pfbFile = "C:/data/FarrarMike/custom.pfb";
      // gcmodelFile = "C:/data/FarrarMike/data/custom.gcmodel";
      // numThreads = 5;

      proj = new Project(filename, logfile);
      if (parsePFB) {
        PFB.populationBAF(proj);
      }
      if (gc5base != null) {
        GcModel.gcModel(proj, gc5base, proj.GC_MODEL_FILENAME.getValue(), 100);
      }
      if (numChunks > 0) {
        if (hmmFile == null || !new File(hmmFile).exists()) {
          hmmFile = Resources.cnv(proj.getLog()).getAllHmm().get();
        }
        if (pfbFile == null || !new File(pfbFile).exists()) {
          pfbFile = Resources.cnv(proj.getLog()).genome(proj.GENOME_BUILD_VERSION.getValue())
                             .getAllPfb().get();
        }
        if (gcmodelFile == null || !new File(pfbFile).exists()) {
          gcmodelFile = Resources.cnv(proj.getLog()).genome(proj.GENOME_BUILD_VERSION.getValue())
                                 .getAllGcmodel().get();
        }
        doBatch(proj, auto, chrx, sexCent, transformData, numChunks, separateQsubs, pfbFile,
                gcmodelFile, hmmFile, separateQsubs ? submit : false, recode, excludes, numThreads);
      }
      if (rawlog != null) {
        parseWarnings(proj, rawlog);
      }
      if (rawcnvs != null) {
        parseResults(proj, rawcnvs, denovoOnly);
      }

      if (cnvFiles != null && outputFile != null) {
        combineResults(proj, cnvFiles, outputFile, recode, removeChr11);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
