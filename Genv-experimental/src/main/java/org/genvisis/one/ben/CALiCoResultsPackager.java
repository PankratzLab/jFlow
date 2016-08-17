package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedHashSet;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

public class CALiCoResultsPackager {

  static class ModelData {

    public static ModelData loadModelData(String modelFile) throws IOException {
      ModelData md = new ModelData(modelFile);
      BufferedReader reader = Files.getAppropriateReader(modelFile);
      reader.readLine();
      String line = null;
      while ((line = reader.readLine()) != null) {
        line = line.trim();
        String[] parts = line.split("[\\s]+");
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < md.RS_INDICES.length; i++) {
          sb.append(parts[md.RS_INDICES[i]]);
          if (i < md.RS_INDICES.length - 1) {
            sb.append(",");
          }
        }
        md.dataMap.put(sb.toString(), parts);
      }
      return md;
    }

    int LOCUS_INDEX;
    int CHR_INDEX;
    int B37_START_INDEX;
    int B37_STOP_INDEX;
    int B36_START_INDEX;
    int B36_STOP_INDEX;
    int[] SNP_INDICES;
    int[] RS_INDICES;
    int MODEL_INDEX;

    int MODEL_DESC_INDEX;

    HashMap<String, String[]> dataMap = new HashMap<String, String[]>();

    private ModelData(String file) throws IOException {
      BufferedReader reader = Files.getAppropriateReader(file);
      String[] hdrParts = reader.readLine().split("[\\s]+");
      int overage = hdrParts.length - 8;

      LOCUS_INDEX = 0;
      CHR_INDEX = 1;
      B37_START_INDEX = 2;
      B37_STOP_INDEX = 3;
      B36_START_INDEX = 4;
      B36_STOP_INDEX = 5;
      SNP_INDICES = new int[overage / 2];
      RS_INDICES = new int[overage / 2];
      for (int i = 0; i < overage / 2; i++) {
        SNP_INDICES[i] = 6 + i;
      }
      for (int i = 0; i < overage / 2; i++) {
        RS_INDICES[i] = 6 + overage / 2 + i;
      }
      MODEL_INDEX = 6 + overage;
      MODEL_DESC_INDEX = 7 + overage;

    }
  }

  static class ModelSNP {
    String modelSNP;
    String modelName;
    int n;
    String locus;
    String chr;
    String bp36start;
    String bp36stop;
    String bp37start;
    String bp37stop;
    String[] snpIDs;
    String[] rsIDs;
    String modelID;
    String modelDesc;
    LinkedHashSet<String> keyList = new LinkedHashSet<String>();
    HashMap<String, String> addLines = new HashMap<String, String>();
    HashMap<String, String> selfLines = new HashMap<String, String>();
    HashMap<String, String> freqLines = new HashMap<String, String>();

    public ModelSNP() {}
  }

  // static final int MODEL_INDEX_LOCUS = 0;
  // static final int MODEL_INDEX_CHR = 1;
  // static final int MODEL_INDEX_B37_START = 2;
  // static final int MODEL_INDEX_B37_STOP = 3; // B37 stop
  // static final int MODEL_INDEX_B36_START = 4; // B36 start
  // static final int MODEL_INDEX_B36_STOP = 5; // B36 stop
  // static final int MODEL_INDEX_RSID = 6; // rsID
  // static final int MODEL_INDEX_ORIG_SNP = 7; // snp
  // static final int MODEL_INDEX_MODEL = 8; // model name
  // static final int MODEL_INDEX_MODEL_DESC = 9; // model desc
  static final int[] MODEL_INDICES = {
      // MODEL_INDEX_LOCUS,
      // MODEL_INDEX_CHR,
      // MODEL_INDEX_B37_START,
      // MODEL_INDEX_B37_STOP,
      // MODEL_INDEX_B36_START,
      // MODEL_INDEX_B36_STOP,
      // MODEL_INDEX_RSID,
      // MODEL_INDEX_ORIG_SNP,
      // MODEL_INDEX_MODEL,
      // MODEL_INDEX_MODEL_DESC,
  };

  static final int[] FREQ_INDICES = {0, // CHR
                                     1, // SNP
                                     2, // A1
                                     3, // A2
                                     4, // MAF
                                     5, // NCHROBS
  };

  static final int[] ASSOC_INDICES = {0, // CHR
                                      1, // SNP
                                      2, // BP
                                      3, // A1
                                      4, // TEST
                                      5, // NMISS
                                      6, // BETA
                                      7, // SE
                                      8, // L95
                                      9, // U95
                                      10, // STAT
                                      11 // P
  };

  static String runDir = "F:/GlucInul/conditional_round/round2/";
  static String resultsDir = "results/";
  static String scratchDir = "scratch/";
  static String modelFile = "round2model.xln";
  static String plinkFileRoot = "F:/GlucInul/plink";

  static String[] hdr1 = {"SNP.original", "rsID", "CHR.build36", "CHR.build37", "bp.build36",
                          "bp.build37", "analysis.locus", "analysis.race"};

  static String[] hdr2 = {"snp.coded", "snp.noncoded", "snp.BETA", "snp.SE", "snp.STAT", "snp.P",
                          "snp.L95", "snp.U95", "snp.NMISS", "snp.CAF", "snp.CAC",};

  static String[] hdrInd =
      {"index#.BETA", "index#.SE", "index#.STAT", "index#.P", "index#.L95", "index#.U95",
       "index#.NMISS", "index#.CAF", "index#.CAC", "index#.coded", "index#.noncoded"};

  private static String[] getFileNames(ModelSNP model) {
    String assoc = runDir + scratchDir + model.modelName + ".assoc.linear";
    String freq = runDir + scratchDir + model.modelName + "_freq.frq";
    return new String[] {assoc, freq};
  }

  private static String getHeader(ModelData md) {
    StringBuilder header = new StringBuilder(Array.toStr(hdr1, "\t"));
    int rsCount = md.RS_INDICES.length;
    for (int i = 0; i < rsCount; i++) {
      header.append("\tindex");
      if (rsCount > 1) {
        header.append(i + 1);
      }
      header.append(".rsID");
    }
    header.append("\t").append(Array.toStr(hdr2, "\t")).append("\t");
    for (int i = 0; i < rsCount; i++) {
      String repl = rsCount > 1 ? (i + 1) + "" : "";
      header.append(Array.toStr(hdrInd, "\t").replaceAll("#", repl)).append("\t");
    }
    return header.toString();
  }

  private static String getStringForSNP(ModelSNP model, String snp, String delim,
                                        PrintWriter logWriter) {
    String assocLine = model.addLines.get(snp);
    String[] lineParts = assocLine.split("[\\s]+");
    // String chr = lineParts[ASSOC_INDICES[0]];
    // String snp = lineParts[ASSOC_INDICES[1]];
    // String bp = lineParts[ASSOC_INDICES[2]];
    String a1assoc = lineParts[ASSOC_INDICES[3]];
    // String test = lineParts[ASSOC_INDICES[4]];
    String nmissStr = lineParts[ASSOC_INDICES[5]];
    String beta = lineParts[ASSOC_INDICES[6]];
    String se = lineParts[ASSOC_INDICES[7]];
    String l95 = lineParts[ASSOC_INDICES[8]];
    String u95 = lineParts[ASSOC_INDICES[9]];
    String stat = lineParts[ASSOC_INDICES[10]];
    String p = lineParts[ASSOC_INDICES[11]];

    String freqLine = model.freqLines.get(snp);
    String[] freqParts = freqLine.split("[\\s]+");
    String a1 = freqParts[FREQ_INDICES[2]];
    String a2 = freqParts[FREQ_INDICES[3]];
    String mafStr = freqParts[FREQ_INDICES[4]];

    String a1final = a1assoc;
    String a2final = (a1assoc.equals(a1) ? a2 : a1);
    if (!a1assoc.equals(a1)) {
      logWriter.println("Error - Mismatched Alleles! SNP: [" + snp + "] ::>> A1_a: [" + a1assoc
                        + "] <-> A1/A2: [" + a1 + "/" + a2 + "]");
    }

    int nmiss = -1;
    String nStr = "NA";
    try {
      nmiss = Integer.parseInt(nmissStr);
      nStr = "" + (model.n - nmiss);
    } catch (NumberFormatException e) {
    }

    double maf = -1;
    try {
      maf = Double.parseDouble(mafStr);
    } catch (NumberFormatException e) {
    }
    double mac = (nmiss == -1 || maf == -1 ? -1d : maf * nmiss * 2d);

    StringBuilder sb = new StringBuilder();
    for (String modSnp : model.rsIDs) {
      sb.append(modSnp).append(delim);
    }

    sb.append(a1final).append(delim).append(a2final).append(delim).append(beta).append(delim)
      .append(se).append(delim).append(stat).append(delim).append(p).append(delim).append(l95)
      .append(delim).append(u95).append(delim).append(nStr).append(delim).append(mafStr)
      .append(delim).append(mac).append(delim);

    for (String rsID : model.rsIDs) {

      freqLine = model.freqLines.get(rsID);
      freqParts = freqLine.split("[\\s]+");
      a1 = freqParts[FREQ_INDICES[2]];
      a2 = freqParts[FREQ_INDICES[3]];
      mafStr = freqParts[FREQ_INDICES[4]];

      String selfLine = model.selfLines.get(snp);
      lineParts = selfLine.split("[\\s]+");
      // String chr = lineParts[ASSOC_INDICES[0]];
      // String snp = lineParts[ASSOC_INDICES[1]];
      // String bp = lineParts[ASSOC_INDICES[2]];
      a1assoc = lineParts[ASSOC_INDICES[3]];
      // String test = lineParts[ASSOC_INDICES[4]];
      nmissStr = lineParts[ASSOC_INDICES[5]];
      beta = lineParts[ASSOC_INDICES[6]];
      se = lineParts[ASSOC_INDICES[7]];
      l95 = lineParts[ASSOC_INDICES[8]];
      u95 = lineParts[ASSOC_INDICES[9]];
      stat = lineParts[ASSOC_INDICES[10]];
      p = lineParts[ASSOC_INDICES[11]];

      nmiss = -1;
      nStr = "NA";
      try {
        nmiss = Integer.parseInt(nmissStr);
        nStr = "" + (model.n - nmiss);
      } catch (NumberFormatException e) {
      }

      maf = -1;
      try {
        maf = Double.parseDouble(mafStr);
      } catch (NumberFormatException e) {
      }
      mac = (nmiss == -1 || maf == -1 ? -1d : maf * nmiss * 2d);

      a1assoc = model.selfLines.get(rsID).split("[\\s]+")[ASSOC_INDICES[3]];
      a1final = a1assoc;
      a2final = (a1assoc.equals(a1) ? a2 : a1);
      if (!a1assoc.equals(a1)) {
        logWriter.println("Error - Mismatched Alleles! INDEXSNP: [" + rsID + "] ::>> A1_a: ["
                          + a1assoc + "] <-> A1/A2: [" + a1 + "/" + a2 + "]");
      }

      sb.append(beta).append(delim).append(se).append(delim).append(stat).append(delim).append(p)
        .append(delim).append(l95).append(delim).append(u95).append(delim).append(nStr)
        .append(delim).append(mafStr).append(delim).append(mac).append(delim).append(a1final)
        .append(delim).append(a2final).append(delim);
    }


    return sb.toString();
  }

  public static void loadDataFromAssocFile(ModelSNP model) throws IOException {
    // CHR SNP BP A1 TEST NMISS BETA SE L95 U95 STAT P
    String filename = getFileNames(model)[0];

    BufferedReader reader = Files.getAppropriateReader(filename);
    reader.readLine();
    String line = null;
    while ((line = reader.readLine()) != null) {
      line = line.trim();
      String[] lineParts = line.split("[\\s]+");
      // String chr = lineParts[ASSOC_INDICES[0]];
      String snp = lineParts[ASSOC_INDICES[1]];
      // String bp = lineParts[ASSOC_INDICES[2]];
      // String a1 = lineParts[ASSOC_INDICES[3]];
      String test = lineParts[ASSOC_INDICES[4]];
      // String nmiss = lineParts[ASSOC_INDICES[5]];
      // String beta = lineParts[ASSOC_INDICES[6]];
      // String se = lineParts[ASSOC_INDICES[7]];
      // String l95 = lineParts[ASSOC_INDICES[8]];
      // String u95 = lineParts[ASSOC_INDICES[9]];
      // String stat = lineParts[ASSOC_INDICES[10]];
      // String p = lineParts[ASSOC_INDICES[11]];

      if (test.equals("ADD")) {
        model.keyList.add(snp);
        model.addLines.put(snp, line);
      } else {
        for (String modSnp : model.rsIDs) {
          if (test.equals(modSnp)) {
            model.keyList.add(snp);
            model.selfLines.put(snp, line);
          }
        }
      }
    }
  }

  public static void loadDataFromFreqFile(ModelSNP model) throws IOException {
    String filename = getFileNames(model)[1];
    BufferedReader reader = Files.getAppropriateReader(filename);
    reader.readLine();
    String line = null;
    while ((line = reader.readLine()) != null) {
      line = line.trim();
      String[] parts = line.split("[\\s]+");
      // String chr = parts[FREQ_INDICES[0];
      String snp = parts[FREQ_INDICES[1]];
      // String a1 = parts[FREQ_INDICES[2]];
      // String a2 = parts[FREQ_INDICES[3]];
      // String maf = parts[FREQ_INDICES[4]];
      // String nchrobs = parts[FREQ_INDICES[5]];
      model.freqLines.put(snp, line);
    }
  }

  public static ModelSNP[] loadModels(ModelData md) throws IOException {
    ArrayList<ModelSNP> keys = new ArrayList<ModelSNP>();
    String[][] keyLines =
        HashVec.loadFileToStringMatrix(runDir + "conditionals.txt", false, new int[] {0, 1}, false);
    for (String[] keyLine : keyLines) {
      ModelSNP m = new ModelSNP();
      m.modelSNP = keyLine[1];
      m.modelName = keyLine[0].substring(0, keyLine[0].indexOf(".")) + "_"
                    + ext.replaceWithLinuxSafeCharacters(keyLine[1], false);
      String file = runDir + resultsDir + m.modelName + "/"
                    + m.modelName.substring(0, m.modelName.indexOf("_")) + "_pheno.dat";
      m.n = Files.countLines(file, 1);

      String[] dataLine = md.dataMap.get(m.modelSNP);
      m.locus = dataLine[md.LOCUS_INDEX];
      m.chr = dataLine[md.CHR_INDEX];
      m.bp37start = dataLine[md.B37_START_INDEX];
      m.bp37stop = dataLine[md.B37_STOP_INDEX];
      m.bp36start = dataLine[md.B36_START_INDEX];
      m.bp36stop = dataLine[md.B36_STOP_INDEX];
      m.snpIDs = new String[md.SNP_INDICES.length];
      for (int i = 0; i < md.SNP_INDICES.length; i++) {
        m.snpIDs[i] = dataLine[md.SNP_INDICES[i]];
      }
      m.rsIDs = new String[md.RS_INDICES.length];
      for (int i = 0; i < md.RS_INDICES.length; i++) {
        m.rsIDs[i] = dataLine[md.RS_INDICES[i]];
      }
      m.modelID = dataLine[md.MODEL_INDEX];
      m.modelDesc = dataLine[md.MODEL_DESC_INDEX];

      keys.add(m);
    }
    return keys.toArray(new ModelSNP[keys.size()]);
  }

  public static void main(String[] args) {
    ModelData md;
    ModelSNP[] models;
    try {
      md = ModelData.loadModelData(runDir + modelFile);
      models = loadModels(md);
    } catch (IOException e) {
      System.err.println("Error - couldn't read model data");
      return;
    }
    Hashtable<String, String> ml = parseMarkerLookup(plinkFileRoot);
    String outputFile = runDir + "results1.xln";
    PrintWriter writer = Files.getAppropriateWriter(outputFile);
    PrintWriter logWriter = Files.getAppropriateWriter(runDir + "resultLog.txt");
    writer.println(getHeader(md));
    for (int i = 0; i < models.length; i++) {
      ModelSNP model = models[i];
      logWriter.println("Processing model " + model.modelName);
      try {
        loadDataFromAssocFile(model);
        loadDataFromFreqFile(model);
        String[] lines = writeLinesForModel(model, ml, logWriter);
        for (String line : lines) {
          writer.println(line);
        }
      } catch (IOException e) {
        logWriter.println("Error - couldn't read from files for model " + model.modelName);
      }
      models[i] = null;
    }
    writer.flush();
    writer.close();
    logWriter.flush();
    logWriter.close();
  }

  public static Hashtable<String, String> parseMarkerLookup(String plinkFileRoot) {
    BufferedReader reader;
    String[] line;
    Hashtable<String, String> hash;
    int count;

    hash = new Hashtable<String, String>();
    try {
      reader = new BufferedReader(new FileReader(plinkFileRoot + ".bim"));
      count = 0;
      while (reader.ready()) {
        line = reader.readLine().trim().split("\\s+");
        hash.put(line[1],
                 ":\t" + count + "\t" + line[0] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5]);
        count++;
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + plinkFileRoot + ".bim"
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + plinkFileRoot + ".bim" + "\"");
      System.exit(2);
    }

    return hash;
  }

  public static String[] writeLinesForModel(ModelSNP model, Hashtable<String, String> markerHash,
                                            PrintWriter logWriter) {
    ArrayList<String> outputLines = new ArrayList<String>();
    String delim = "\t";
    for (String snp : model.keyList) {
      String[] snpBIMData = markerHash.get(snp).split("[\\s]+");
      StringBuilder sb = new StringBuilder();
      sb.append(snp).append(delim).append(delim).append(snpBIMData[2]).append(delim)// CHR.build36
        .append(delim)// CHR.build37
        .append(snpBIMData[1]).append(delim)// bp.build36
        .append(delim)// bp.build37
        .append(model.locus).append(delim)// analysis.locus
        .append("AA").append(delim)// analysis.race
        .append(getStringForSNP(model, snp, delim, logWriter));
      outputLines.add(sb.toString());
    }
    return outputLines.toArray(new String[outputLines.size()]);
  }
}
