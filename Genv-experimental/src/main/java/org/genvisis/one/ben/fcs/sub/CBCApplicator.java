package org.genvisis.one.ben.fcs.sub;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class CBCApplicator {
  
  private static final String[][] PANELS = {
    {"panel 1", "p1"},
    {"panel 2", "p2"},
  };
  
  private static final FilenameFilter FILTER_P1 = new FilenameFilter() {
    @Override
    public boolean accept(File dir, String name) {
      return (name.toLowerCase().contains(PANELS[0][0]) || name.toLowerCase().contains(PANELS[0][1]));
    }
  };
  private static final FilenameFilter FILTER_P2 = new FilenameFilter() {
    @Override
    public boolean accept(File dir, String name) {
      return (name.toLowerCase().contains(PANELS[1][0]) || name.toLowerCase().contains(PANELS[1][1]));
    }
  };
  
  private static final String PNL_1_HDR_START = "Lymphocytes (SSC-A v FSC-A) | Freq. of Parent"; 
  private static final String PNL_1_CBC_START = "ALYMP x10e9/L";
  
  private static final String PNL_2_HDR_START = "PBMCs (SSC-A v FSC-A) | Freq. of Parent";
  private static final String[] PNL_2_CBC_COLS = {"WBC x10e9/L", "ANEU x10e9/L", "AEOS x10e9/L", "ABASO x10e9/L"};

  private static final String PANEL_1_UNITS = PNL_1_CBC_START.split(" ")[1];
  
  private String cbcDir = "F:/Flow/CBC_processing/cbc/";
  private String dataDir = "F:/Flow/CBC_processing/data/";
  private String outDir = null;
  private String[] filesP1;
  private String[] filesP2;
  private int panel1DataColumn;
  
  private Logger log = new Logger();
  
  HashMap<String, String[]> idMap = new HashMap<String, String[]>();
  
  private void loadCBCs() {
    String[] files = new File(cbcDir).list();
    log.reportTime("Loading " + files.length + " CBC files from " + cbcDir);
    for (String file : files) {
      String[][] data = HashVec.loadFileToStringMatrix(cbcDir + file, true, null, false);
      for (String[] line : data) {
        if (!idMap.containsKey(line[0])) {
          idMap.put(line[0], line);
        } else {
          log.reportTimeError("Duplicate ID in CBC data: " + line[0]);
        }
      }
    }
    panel1DataColumn = ext.indexOfStr(PNL_1_CBC_START, Files.getHeaderOfFile(cbcDir + files[0], files[0].endsWith(".csv") ? "," : null, log));
    // TODO load panel 2 column indices
  }
  
  private void discoverDataFiles() {
    File dataDirFile = new File(dataDir);
    filesP1 = dataDirFile.list(FILTER_P1);
    filesP2 = dataDirFile.list(FILTER_P2);
    log.reportTime("Discovered " + filesP1.length + " Panel 1 files and " + filesP2.length + " Panel 2 files in " + dataDir);
  }
  
  private void runPanel1() {
    if (filesP1 == null || filesP1.length == 0) {
      log.reportTime("No Panel 1 files available!");
      return;
    }
    for (int i = 0; i < filesP1.length; i++) {
      processPanel1File(filesP1[i]);
    }
  }

  private void runPanel2() {
    if (filesP2 == null || filesP2.length == 0) {
      log.reportTime("No Panel 2 files available!");
      return;
    }
    for (int i = 0; i < filesP2.length; i++) {
      processPanel2File(filesP2[i]);
    }
  }
  
  private String getOutFile(String file) {
    return ((outDir == null || "".equals(outDir)) ? dataDir : outDir) + ext.rootOf(file, true) + "_COUNTS.xln";
  }
  
  private void processPanel1File(String file) {
    BufferedReader reader;
    PrintWriter writer;
    String outFile, line, id, cbcCnt;
    String[] header, parts, idParts, cbcData;
    StringBuilder outLine;
    
    try {
      reader = Files.getAppropriateReader(dataDir + file);
    } catch (FileNotFoundException e) {
      log.reportTimeError("Data file " + file + " not found!");
      return;
    }
    outFile = getOutFile(file);
    if (Files.exists(outFile)) {
      log.reportError("Error - output file " + outFile + " already exists!");
      try {
        reader.close();
      } catch (IOException e) {}
      return;
    }
    writer = Files.getAppropriateWriter(outFile);
    
    try {
      header = reader.readLine().replaceAll("\"", "").split("\t", -1);
      outLine = new StringBuilder();
      for (int i = 1; i < header.length; i++) {
        outLine.append("\t").append(header[i].split("\\|")[0]).append("| ").append(PANEL_1_UNITS);
      }
      writer.println(outLine.toString());
    } catch (IOException e) {
      log.reportTimeError("Error - unable to read file " + file);
      try {
        reader.close();
      } catch (IOException e1) {}
      writer.close();
      (new File(outFile)).delete();
      return;
    }
    
    log.reportTime("Processing file " + file);
    line = null;
    try {
      while((line = reader.readLine()) != null) {
        parts = line.split("\t", -1);
        if (parts[0].indexOf("_") == -1) continue; // skip mean/sd lines - TODO compute mean/sd of counts at end?
        idParts = parts[0].split("_");
        id = idParts[idParts.length - 2]; // id is always second to last token
        cbcData = idMap.get(id);
        if (cbcData == null) {
          log.reportTimeWarning("id not found: " + id);
          continue;
        }
        cbcCnt = cbcData[panel1DataColumn];
        if ("NULL".equals(cbcCnt) || "NA".equals(cbcCnt)) {
          writer.println(id + "\t" + Array.toStr(Array.stringArray(parts.length - 1, "NaN")));
        } else {
          double[] cnts = new double[parts.length - 1];
          cnts[0] = Double.parseDouble(cbcCnt);
          
          for (int i = 2; i < parts.length; i++) {
            String column = header[i];
            double pct = Double.parseDouble(parts[i].replace("%", "")) / 100;
            double parentCnt = cnts[getParentIndexPanel1(header, column) - 1];
            cnts[i - 1] = parentCnt * pct;
          }
          
          outLine = new StringBuilder(parts[0]);
          for (int i = 0; i < cnts.length; i++) {
            outLine.append("\t").append(cnts[i]);
          }
          writer.println(outLine.toString());
        }
      }
    } catch (IOException e) {
      log.reportTimeError("Error while reading file " + file + " - aborting.  Partial output file will be removed.");
      try {
        reader.close();
      } catch (IOException e1) {}
      writer.close();
      (new File(outFile)).delete();
      return;
    }
    
    writer.flush();
    writer.close();
    try {
      reader.close();
    } catch (IOException e) {}
    
  }
  
  private void processPanel2File(String file) {
    return;
  }
  
  private int getParentIndexPanel1(String[] header, String column) {
    String[] parts = column.split("\\|");
    String[] path = parts[0].split("/");
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < path.length - 1; i++) {
      sb.append(path[i]);
      if (i < path.length - 2) {
        sb.append("/");
      }
    }
    sb.append(" |").append(parts[1]);
    int val = ext.indexOfStr(sb.toString(), header); 
    if (val == -1 && sb.charAt(0) == ('\"')) {
      sb.delete(0, 1);
      val = ext.indexOfStr(sb.toString(), header); 
      if (val == -1 && sb.charAt(sb.length() - 1) == ('\"')) {
        sb.delete(sb.length() - 1, sb.length());
        val = ext.indexOfStr(sb.toString(), header); 
      }
    }
    return val;
  }
  
  public void run() {
    loadCBCs();
    discoverDataFiles();
    runPanel1();
    runPanel2();
  }
  
  public static void main(String[] args) {
    new CBCApplicator().run();
  }

}
