package org.genvisis.one.ben.fcs.sub;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;

import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class CBCApplicator implements Runnable {
  
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
  
  private static final String PNL_1_CBC_START = "ALYMP x10e9/L";
  private static final String[] PNL_2_CBC_COLS = {"WBC x10e9/L", "ANEU x10e9/L", "AEOS x10e9/L", "ABASO x10e9/L"};
  private static final String UNITS = "x10e9/L";
  
  private int panel1DataColumn;
  private int[] panel2DataColumns;
  private String cbcDir = null;
  private String dataDir = null;
  private String outDir = null;
  private String[] filesP1;
  private String[] filesP2;
  
  private Logger log = new Logger();
  
  HashMap<String, String[]> idMap = new HashMap<String, String[]>();

  private JProgressBar progressBar;
  
  public void setCBCDir(String cbcD) {
    this.cbcDir = cbcD;
  }
  
  public void setDataDir(String dataD) {
    this.dataDir = dataD;
  }
  
  public void setOutputDirectory(String outD) {
    this.outDir = outD;
  }
  
  public void setLog(Logger log) {
    this.log = log;
  }
  
  public void setProgressBar(JProgressBar prog) {
    this.progressBar = prog;
  }

  private void loadCBCs() {
    if (cbcDir == null || "".equals(cbcDir)) {
      log.reportTimeError("CBC directory not set!");
      return;
    }
    if (!Files.exists(cbcDir)) {
      log.reportTimeError("CBC directory not found: " + cbcDir);
      return;
    }
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
    String[] hdr = Files.getHeaderOfFile(cbcDir + files[0], files[0].endsWith(".csv") ? "," : null, log);
    if (hdr == null) {
      log.reportTimeError("couldn't load header of file: " + cbcDir + files[0]);
      return;
    }
    panel1DataColumn = ext.indexOfStr(PNL_1_CBC_START, hdr);
    panel2DataColumns = ext.indexFactors(PNL_2_CBC_COLS, hdr, false, false);
    if (panel1DataColumn == -1) {
      log.reportTimeError("missing column " + PNL_1_CBC_START + " from CBC file: " + cbcDir + files[0]);
      return;
    }
    if (Array.countIf(panel2DataColumns, -1) > 0) {
      log.reportTimeError("missing column " + PNL_2_CBC_COLS[Arrays.asList(panel2DataColumns).indexOf(-1)] + " CBC file: " + cbcDir + files[0]);
      return;
    }
  }
  
  private void discoverDataFiles() {
    if (dataDir == null || "".equals(dataDir)) {
      log.reportTimeError("data file directory not set!");
      return;
    }
    if (!Files.exists(dataDir)) {
      log.reportTimeError("data file directory not found: " + dataDir);
      return;
    }
    File dataDirFile = new File(dataDir);
    filesP1 = dataDirFile.list(FILTER_P1);
    filesP2 = dataDirFile.list(FILTER_P2);
    log.reportTime("Discovered " + filesP1.length + " Panel_1 files and " + filesP2.length + " Panel_2 files in " + dataDir);
  }
  
  private void runPanel1() {
    if (filesP1 == null || filesP1.length == 0) {
      log.reportTime("No Panel_1 files available!");
      return;
    }
    if (progressBar != null) {
      SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
          progressBar.setString("Panel 1");
          progressBar.setMinimum(0);
          progressBar.setMaximum(filesP1.length);
        }
      });
    }
    for (int i = 0; i < filesP1.length; i++) {
      processPanel1File(filesP1[i]);
      if (progressBar != null) {
        final int ind = i + 1;
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            progressBar.setValue(ind);
            progressBar.setString("Panel 1: (" + ind + " / " + filesP1.length + ")");
          }
        });
      }
    }
  }

  private void runPanel2() {
    if (filesP2 == null || filesP2.length == 0) {
      log.reportTime("No Panel_2 files available!");
      return;
    }
    if (progressBar != null) {
      SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
          progressBar.setString("Panel 2");
          progressBar.setMinimum(0);
          progressBar.setMaximum(filesP2.length);
        }
      });
    }
    for (int i = 0; i < filesP2.length; i++) {
      processPanel2File(filesP2[i]);
      if (progressBar != null) {
        final int ind = i + 1;
        SwingUtilities.invokeLater(new Runnable() {
          @Override
          public void run() {
            progressBar.setValue(ind);
            progressBar.setString("Panel 2: (" + ind + " / " + filesP2.length + ")");
          }
        });
      }
    }
  }
  
  private void processPanel1File(String file) {
    BufferedReader reader;
    PrintWriter writer;
    String outFile, line, id, cbcCnt, delim;
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
      log.reportTimeError("output file " + outFile + " already exists!");
      try {
        reader.close();
      } catch (IOException e) {}
      return;
    }
    writer = Files.getAppropriateWriter(outFile);
    
    try {
      line = reader.readLine();
      delim = file.endsWith("csv") ? "," : ext.determineDelimiter(line);
      header = delim.equals(",") ? ext.splitCommasIntelligently(line, true, log) : line.replaceAll("\"", "").split(delim, -1);
      outLine = new StringBuilder();
      for (int i = "".equals(header[0]) ? 1 : 0; i < header.length; i++) {
        if ("".equals(header[i])) continue;
        outLine.append("\t").append(header[i].split("\\|")[0]).append("| ").append(PNL_1_CBC_START);
      }
      writer.println(outLine.toString());
    } catch (IOException e) {
      log.reportTimeError("unable to read file " + file);
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
        parts = line.split(delim, -1);
        if (parts[0].indexOf("_") == -1) continue; // skip mean/sd lines - TODO compute mean/sd of counts at end?
        idParts = parts[0].split("_");
        id = idParts[idParts.length - 2]; // id is always second to last token
        cbcData = idMap.get(id);
        if (cbcData == null) {
          log.reportTimeWarning("CBC data not found for ID: " + id);
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
            if ("".equals(column)) continue;
            String temp = parts[i].replace("%", "");
            double pct = Double.parseDouble(temp) / 100;
            double parentCnt = cnts[getParentIndex(header, column) - 1];
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
      log.reportTimeError("Exception occurred while reading file " + file + " - aborting.  Partial output file will be removed.");
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
    BufferedReader reader;
    PrintWriter writer;
    String outFile, line, id, delim;
    double cbcCnt;
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
      log.reportTimeError("output file " + outFile + " already exists!");
      try {
        reader.close();
      } catch (IOException e) {}
      return;
    }
    writer = Files.getAppropriateWriter(outFile);
    
    try {
      line = reader.readLine();
      delim = file.endsWith("csv") ? "," : ext.determineDelimiter(line);
      header = delim.equals(",") ? ext.splitCommasIntelligently(line, true, log) : line.replaceAll("\"", "").split(delim, -1);
      outLine = new StringBuilder();
      for (int i = 1; i < header.length; i++) {
        if ("".equals(header[i])) continue;
        outLine.append("\t").append(header[i].split("\\|")[0]).append("| ").append(UNITS);
        if (i == 1) {
          outLine.append(" (");
          outLine.append(PNL_2_CBC_COLS[0].split(" ")[0]);
          for (int j = 1; j < PNL_2_CBC_COLS.length; j++) {
            outLine.append(" - ").append(PNL_2_CBC_COLS[j].split(" ")[0]);
          }
          outLine.append(")");
        }
      }
      writer.println(outLine.toString());
    } catch (IOException e) {
      log.reportTimeError("unable to read file " + file);
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
        parts = line.split(delim, -1);
        if (parts[0].indexOf("_") == -1) continue; // skip mean/sd lines - TODO compute mean/sd of counts at end?
        idParts = parts[0].split("_");
        id = idParts[idParts.length - 2]; // id is always second to last token
        cbcData = idMap.get(id);
        if (cbcData == null) {
          log.reportTimeWarning("id not found: " + id);
          continue;
        }
        cbcCnt = getCBCCountPanel2(cbcData);
        if (Double.isNaN(cbcCnt)) {
          writer.println(id + "\t" + Array.toStr(Array.stringArray(parts.length - 1, "NaN")));
        } else {
          double[] cnts = new double[parts.length - 1];
          cnts[0] = cbcCnt;
          
          for (int i = 2; i < parts.length; i++) {
            String column = header[i];
            if ("".equals(column)) continue;
            double pct = Double.parseDouble(parts[i].replace("%", "")) / 100;
            double parentCnt = cnts[getParentIndex(header, column) - 1];
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
      log.reportTimeError("Exception occurred while reading file " + file + " - aborting.  Partial output file will be removed.");
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
  
  private double getCBCCountPanel2(String[] cbcData) {
    double[] cbcs = new double[PNL_2_CBC_COLS.length];
    for (int i = 0; i < panel2DataColumns.length; i++) {
      if ("NULL".equals(cbcData[panel2DataColumns[i]]) || "NA".equals(cbcData[panel2DataColumns[i]])) {
        return Double.NaN;
      }
      cbcs[i] = Double.parseDouble(cbcData[panel2DataColumns[i]]);
    }
    double cbc = cbcs[0];
    for (int i = 1; i < cbcs.length; i++) {
      cbc -= cbcs[i];
    }
    return cbc;
  }

  private String getOutFile(String file) {
    return ((outDir == null || "".equals(outDir) || !Files.exists(outDir)) ? dataDir : outDir) + ext.rootOf(file, true) + "_COUNTS.tab.txt";
  }

  private int getParentIndex(String[] header, String column) {
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
    if (progressBar != null) {
      SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
          progressBar.setString("Done!");
        }
      });
    }
  }
  
}
