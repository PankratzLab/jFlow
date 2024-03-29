package org.genvisis.fcs.sub;

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

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

public class CBCApplicator implements Runnable {

  private static final String[][] PANELS = {
    {"panel 1", "p1"}, {"panel 2", "p2"},
  };

  private static final FilenameFilter FILTER_P1 =
      new FilenameFilter() {

        @Override
        public boolean accept(File dir, String name) {
          return (name.toLowerCase().contains(PANELS[0][0])
              || name.toLowerCase().contains(PANELS[0][1]));
        }
      };
  private static final FilenameFilter FILTER_P2 =
      new FilenameFilter() {

        @Override
        public boolean accept(File dir, String name) {
          return (name.toLowerCase().contains(PANELS[1][0])
              || name.toLowerCase().contains(PANELS[1][1]));
        }
      };

  private static final String PNL_1_CBC_START = "ALYMP x10e9/L";
  private static final String[] PNL_2_CBC_COLS = {
    "WBC x10e9/L", "ANEU x10e9/L", "AEOS x10e9/L", "ABASO x10e9/L"
  };
  private static final String UNITS = "x10e9/L";

  private int panel1DataColumn;
  private int[] panel2DataColumns;
  private String cbcDir = null;
  private String dataDir = null;
  private String outDir = null;
  private String[] filesP1;
  private String[] filesP2;

  private Logger log = new Logger();

  HashMap<String, String[]> idMap = new HashMap<>();

  private JProgressBar progressBar;

  public void setCBCDir(String cbcD) {
    cbcDir = cbcD;
  }

  public void setDataDir(String dataD) {
    dataDir = dataD;
  }

  public void setOutputDirectory(String outD) {
    outDir = outD;
  }

  public void setLog(Logger log) {
    this.log = log;
  }

  public void setProgressBar(JProgressBar prog) {
    progressBar = prog;
  }

  private void loadCBCs() {
    if (cbcDir == null || "".equals(cbcDir)) {
      log.reportError("CBC directory not set!");
      return;
    }
    if (!Files.exists(cbcDir)) {
      log.reportError("CBC directory not found: " + cbcDir);
      return;
    }
    String[] files = new File(cbcDir).list();
    log.report("Loading " + files.length + " CBC files from " + cbcDir);
    for (String file : files) {
      String[][] data = HashVec.loadFileToStringMatrix(cbcDir + file, true, null);
      for (String[] line : data) {
        if (!idMap.containsKey(line[0])) {
          idMap.put(line[0], line);
        } else {
          log.reportError("Duplicate ID in CBC data: " + line[0]);
        }
      }
    }
    String[] hdr =
        Files.getHeaderOfFile(cbcDir + files[0], files[0].endsWith(".csv") ? "," : null, log);
    if (hdr == null) {
      log.reportError("couldn't load header of file: " + cbcDir + files[0]);
      return;
    }
    panel1DataColumn = ext.indexOfStr(PNL_1_CBC_START, hdr);
    panel2DataColumns = ext.indexFactors(PNL_2_CBC_COLS, hdr, false);
    if (panel1DataColumn == -1) {
      log.reportError("missing column " + PNL_1_CBC_START + " from CBC file: " + cbcDir + files[0]);
      return;
    }
    if (ArrayUtils.countIf(panel2DataColumns, -1) > 0) {
      log.reportError(
          "missing column "
              + PNL_2_CBC_COLS[Arrays.asList(panel2DataColumns).indexOf(-1)]
              + " CBC file: "
              + cbcDir
              + files[0]);
      return;
    }
  }

  private void discoverDataFiles() {
    if (dataDir == null || "".equals(dataDir)) {
      log.reportError("data file directory not set!");
      return;
    }
    if (!Files.exists(dataDir)) {
      log.reportError("data file directory not found: " + dataDir);
      return;
    }
    File dataDirFile = new File(dataDir);
    filesP1 = dataDirFile.list(FILTER_P1);
    filesP2 = dataDirFile.list(FILTER_P2);
    log.report(
        "Discovered "
            + filesP1.length
            + " Panel_1 files and "
            + filesP2.length
            + " Panel_2 files in "
            + dataDir);
  }

  private void runPanel1() {
    if (filesP1 == null || filesP1.length == 0) {
      log.report("No Panel_1 files available!");
      return;
    }
    if (progressBar != null) {
      SwingUtilities.invokeLater(
          new Runnable() {

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
        SwingUtilities.invokeLater(
            new Runnable() {

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
      log.report("No Panel_2 files available!");
      return;
    }
    if (progressBar != null) {
      SwingUtilities.invokeLater(
          new Runnable() {

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
        SwingUtilities.invokeLater(
            new Runnable() {

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
    boolean skipFirstCol = false;

    try {
      reader = Files.getAppropriateReader(dataDir + file);
    } catch (FileNotFoundException e) {
      log.reportError("Data file " + file + " not found!");
      return;
    } catch (IOException e) {
      log.reportException(e);
      return;
    }
    outFile = getOutFile(file);
    if (Files.exists(outFile)) {
      log.reportError("output file " + outFile + " already exists!");
      try {
        reader.close();
      } catch (IOException e) {
      }
      return;
    }
    writer = Files.getAppropriateWriter(outFile);

    try {
      line = reader.readLine();
      delim =
          file.endsWith("csv") ? "," : file.endsWith("xln") ? "\t" : ext.determineDelimiter(line);
      header =
          delim.equals(",")
              ? ext.splitCommasIntelligently(line, true, log)
              : line.replaceAll("\"", "").split(delim, -1);
      outLine = new StringBuilder();
      for (int i = "".equals(header[0]) ? 1 : 0; i < header.length; i++) {
        if ("".equals(header[i])) {
          continue;
        }
        if (header[i].indexOf('|') == -1) {
          if (i == 1) {
            skipFirstCol = true;
          }
          continue;
        }
        outLine.append("\t").append(header[i].split("\\|")[0]).append("| ").append(PNL_1_CBC_START);
      }
      writer.println(outLine.toString());
    } catch (IOException e) {
      log.reportError("unable to read file " + file);
      try {
        reader.close();
      } catch (IOException e1) {
      }
      writer.close();
      (new File(outFile)).delete();
      return;
    }

    log.report("Processing file " + file);
    line = null;
    try {
      while ((line = reader.readLine()) != null) {
        parts = line.split(delim, -1);
        if (parts[0].indexOf("_") == -1) {
          continue; // skip mean/sd lines - TODO compute mean/sd of counts at end?
        }
        idParts = parts[0].split("_");
        id = idParts[idParts.length - 2]; // id is always second to last token
        cbcData = idMap.get(id);
        if (cbcData == null) {
          log.reportWarning("CBC data not found for ID: " + id);
          continue;
        }
        cbcCnt = cbcData[panel1DataColumn];
        if ("NULL".equals(cbcCnt) || "NA".equals(cbcCnt)) {
          writer.println(
              id
                  + "\t"
                  + ArrayUtils.toStr(
                      ArrayUtils.stringArray(parts.length - (skipFirstCol ? 2 : 1), "NaN")));
        } else {
          double[] cnts = new double[parts.length - (skipFirstCol ? 2 : 1)];
          cnts[0] = Double.parseDouble(cbcCnt);

          for (int i = (skipFirstCol ? 3 : 2); i < parts.length; i++) {
            String column = header[i];
            if ("".equals(column)) {
              continue;
            }
            if (column.indexOf('|') == -1) {
              continue;
            }
            String temp = parts[i].replace("%", "");
            double pct = Double.parseDouble(temp) / 100;
            double parentCnt = cnts[getParentIndex(header, column) - (skipFirstCol ? 2 : 1)];
            cnts[i - (skipFirstCol ? 2 : 1)] = parentCnt * pct;
          }

          outLine = new StringBuilder(parts[0]);
          for (double cnt : cnts) {
            outLine.append("\t").append(cnt);
          }
          writer.println(outLine.toString());
        }
      }
    } catch (IOException e) {
      log.reportError(
          "Exception occurred while reading file "
              + file
              + " - aborting.  Partial output file will be removed.");
      try {
        reader.close();
      } catch (IOException e1) {
      }
      writer.close();
      (new File(outFile)).delete();
      return;
    }

    writer.flush();
    writer.close();
    try {
      reader.close();
    } catch (IOException e) {
    }
  }

  private void processPanel2File(String file) {
    BufferedReader reader;
    PrintWriter writer;
    String outFile, line, id, delim;
    double cbcCnt;
    String[] header, parts, idParts, cbcData;
    StringBuilder outLine;
    boolean skipFirstCol = false;

    try {
      reader = Files.getAppropriateReader(dataDir + file);
    } catch (FileNotFoundException e) {
      log.reportError("Data file " + file + " not found!");
      return;
    } catch (IOException e) {
      log.reportException(e);
      return;
    }
    outFile = getOutFile(file);
    if (Files.exists(outFile)) {
      log.reportError("output file " + outFile + " already exists!");
      try {
        reader.close();
      } catch (IOException e) {
      }
      return;
    }
    writer = Files.getAppropriateWriter(outFile);

    try {
      line = reader.readLine();
      delim = file.endsWith("csv") ? "," : ext.determineDelimiter(line);
      header =
          delim.equals(",")
              ? ext.splitCommasIntelligently(line, true, log)
              : line.replaceAll("\"", "").split(delim, -1);
      outLine = new StringBuilder();
      for (int i = 1; i < header.length; i++) {
        if ("".equals(header[i])) {
          continue;
        }
        if (header[i].indexOf('|') == -1) {
          if (i == 1) {
            skipFirstCol = true;
          }
          continue;
        }
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
      log.reportError("unable to read file " + file);
      try {
        reader.close();
      } catch (IOException e1) {
      }
      writer.close();
      (new File(outFile)).delete();
      return;
    }

    log.report("Processing file " + file);
    line = null;
    try {
      while ((line = reader.readLine()) != null) {
        parts = line.split(delim, -1);
        if (parts[0].indexOf("_") == -1) {
          continue; // skip mean/sd lines - TODO compute mean/sd of counts at end?
        }
        idParts = parts[0].split("_");
        id = idParts[idParts.length - 2]; // id is always second to last token
        cbcData = idMap.get(id);
        if (cbcData == null) {
          log.reportWarning("id not found: " + id);
          continue;
        }
        cbcCnt = getCBCCountPanel2(cbcData);
        if (Double.isNaN(cbcCnt)) {
          writer.println(
              id
                  + "\t"
                  + ArrayUtils.toStr(
                      ArrayUtils.stringArray(parts.length - (skipFirstCol ? 2 : 1), "NaN")));
        } else {
          double[] cnts = new double[parts.length - (skipFirstCol ? 2 : 1)];
          cnts[0] = cbcCnt;

          for (int i = skipFirstCol ? 3 : 2; i < parts.length; i++) {
            String column = header[i];
            if ("".equals(column)) {
              continue;
            }
            if (column.indexOf('|') == -1) {
              continue;
            }
            String temp = parts[i].replace("%", "");
            double pct = Double.parseDouble(temp) / 100;
            double parentCnt = cnts[getParentIndex(header, column) - (skipFirstCol ? 2 : 1)];
            cnts[i - (skipFirstCol ? 2 : 1)] = parentCnt * pct;
          }

          outLine = new StringBuilder(parts[0]);
          for (double cnt : cnts) {
            outLine.append("\t").append(cnt);
          }
          writer.println(outLine.toString());
        }
      }
    } catch (IOException e) {
      log.reportError(
          "Exception occurred while reading file "
              + file
              + " - aborting.  Partial output file will be removed.");
      try {
        reader.close();
      } catch (IOException e1) {
      }
      writer.close();
      (new File(outFile)).delete();
      return;
    }

    writer.flush();
    writer.close();
    try {
      reader.close();
    } catch (IOException e) {
    }
  }

  private double getCBCCountPanel2(String[] cbcData) {
    double[] cbcs = new double[PNL_2_CBC_COLS.length];
    for (int i = 0; i < panel2DataColumns.length; i++) {
      if ("NULL".equals(cbcData[panel2DataColumns[i]])
          || "NA".equals(cbcData[panel2DataColumns[i]])) {
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
    return ((outDir == null || "".equals(outDir) || !Files.exists(outDir)) ? dataDir : outDir)
        + ext.rootOf(file, true)
        + "_COUNTS.tab.txt";
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

  @Override
  public void run() {
    loadCBCs();
    discoverDataFiles();
    runPanel1();
    runPanel2();
    if (progressBar != null) {
      SwingUtilities.invokeLater(
          new Runnable() {

            @Override
            public void run() {
              progressBar.setString("Done!");
            }
          });
    }
  }
}
