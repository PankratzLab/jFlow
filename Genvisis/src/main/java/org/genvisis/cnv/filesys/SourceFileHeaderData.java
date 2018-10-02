package org.genvisis.cnv.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Optional;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import org.genvisis.cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import org.genvisis.cnv.manage.SourceFileParser;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

public class SourceFileHeaderData implements Serializable {

  /**
   *
   */
  private static final long serialVersionUID = -6906302109843776908L;

  private String gsgtVersion;
  private String processingDate;
  private String content;
  private int numSnps;
  private int totalSnps;
  private int numSamples;
  private int totalSamples;
  private int numFiles = -1;
  private int currFile;
  private final int columnHeaderLineIndex;

  private int colSampleIdent = -1;
  private int colSnpIdent = -1;
  private int colGenoAB1 = -1;
  private int colGenoAB2 = -1;
  private int colGeno1 = -1;
  private int colGeno2 = -1;
  private int colX = -1;
  private int colY = -1;
  private int colTheta = -1;
  private int colR = -1;
  private int colXRaw = -1;
  private int colYRaw = -1;
  private int colBAF = -1;
  private int colLRR = -1;
  private int colGC = -1;
  private final String headerString;
  private final String delimiter;
  private final String[] cols;

  // Order from ParseIllumina
  // 0 GC
  // 1 XRAW
  // 2 YRAW
  // 3 X
  // 4 Y
  // 5 Theta
  // 6 R
  // 7 B Allele Freq
  // 8 Log R Ratio

  private SourceFileHeaderData(String file, Logger log) throws Elision, IOException {
    BufferedReader reader = Files.getAppropriateReader(file);
    String line = null;
    int lineCnt = 0;
    String delim = ",";
    while ((line = reader.readLine()) != null) {
      delim = ext.determineDelimiter(line, true); // TODO if file ends with .csv [or contains
                                                 // .csv?], can assume that delim is ','. Sim., if
                                                 // ends with .xln, can assume delim is '\t'
      if (",".equals(delim)) {
        delim = "[\\s]*,[\\s]*"; // ext.indexFactors doesn't call trim()
      }
      if ("[Data]".equals(line) || line.startsWith("rs") || line.toUpperCase().startsWith("SNP")
          || ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line.split(delim), false, true,
                              false)[0] != -1) {
        break;
      }
      String[] parts = line.trim().split(",");
      processLine(parts);
      lineCnt++;
    }
    if (!"[Data]".equals(line) && !(line.startsWith("rs") || line.toUpperCase().startsWith("SNP")
                                    || ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS,
                                                        line.split(delim), false, true,
                                                        false)[0] != -1)) {
      log.reportError("Error - malformed or missing header.");
      throw new Elision(file);
    }
    if ("[Data]".equals(line)) {
      line = reader.readLine();
      delim = file.contains(".csv") ? "[\\s]*,[\\s]*"
                                    : file.contains(".xln") ? "[\\s]*\t[\\s]*"
                                                            : ext.determineDelimiter(line, true);
      lineCnt++;
      if (!(line.startsWith("rs") || line.toUpperCase().startsWith("SNP")
            || ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line.split(delim), false, true,
                                false)[0] != -1)) {
        log.reportError("Error - malformed or missing header.  Header must start with 'rs' or 'SNP' or contain one of the following: "
                        + ArrayUtils.toStr(SourceFileParser.SNP_HEADER_OPTIONS[0]) + ".");
        throw new Elision(file);
      }
    }
    reader.close(); // done
    reader = null; // release
    String columnHeaders = line;
    delim = file.contains(".csv") ? SOURCE_FILE_DELIMITERS.COMMA.getDelimiter()
                                  : file.contains(".xln") ? SOURCE_FILE_DELIMITERS.TAB.getDelimiter()
                                                          : SOURCE_FILE_DELIMITERS.getDelimiter(ext.determineDelimiter(columnHeaders,
                                                                                                                       true))
                                                                                  .getDelimiter();
    // if (",".equals(delim)) {
    // delim = "[\\s]*,[\\s]*";
    // }
    cols = columnHeaders.split(delim);
    parseColumnsBestGuess();
    delimiter = delim;
    headerString = columnHeaders;
    columnHeaderLineIndex = lineCnt;

  }

  public static SourceFileHeaderData parseHeader(String file, Logger log) throws Elision,
                                                                          IOException {
    return new SourceFileHeaderData(file, log);
  }

  private static final String[] SAMPLE_FIELD_ID = {"Sample ID", "Sample Name"};

  private static final String[] DATA_FIELDS_GC = {"GC Score", "GCscore", "confidence",
                                                  "confidenceScore"};
  private static final String[] DATA_FIELDS_XRAW = {"X Raw"};
  private static final String[] DATA_FIELDS_YRAW = {"Y Raw"};
  private static final String[] DATA_FIELDS_X = {"X", "Xvalue", "Log Ratio", "intensity_1",
                                                 "Signal A"};
  private static final String[] DATA_FIELDS_Y = {"Y", "Yvalue", "Strength", "intensity_2",
                                                 "Signal B"};
  private static final String[] DATA_FIELDS_THETA = {"Theta"};
  private static final String[] DATA_FIELDS_R = {"R"};
  private static final String[] DATA_FIELDS_BAF = {"B Allele Freq", "BAF"};
  private static final String[] DATA_FIELDS_LRR = {"Log R Ratio", "LRR"};
  private static final String[] GENOTYPE_FIELDS_A1_FOR = {"Allele1 - Forward", "Allele1",
                                                          "genotype1", "Allele1 - Top",
                                                          "Forward Strand Base Calls",
                                                          "Forced Call", "Forced Call Codes"};
  private static final String[] GENOTYPE_FIELDS_A2_FOR = {"Allele2 - Forward", "Allele B",
                                                          "genotype2", "Allele2 - Top",
                                                          "Forward Strand Base Calls",
                                                          "Forced Call", "Forced Call Codes"};
  private static final String[] GENOTYPE_FIELDS_A1_AB = {"Allele1 - AB", "Call Codes", "Call"};
  private static final String[] GENOTYPE_FIELDS_A2_AB = {"Allele2 - AB", "Call Codes", "Call"};

  private static final String[][] LOOKUP = {/* 0 */ SourceFileParser.SNP_HEADER_OPTIONS[0],
                                            /* 1 */ SAMPLE_FIELD_ID, /* 2 */ DATA_FIELDS_GC,
                                            /* 3 */ DATA_FIELDS_XRAW, /* 4 */ DATA_FIELDS_YRAW,
                                            /* 5 */ DATA_FIELDS_X, /* 6 */ DATA_FIELDS_Y,
                                            /* 7 */ DATA_FIELDS_THETA, /* 8 */ DATA_FIELDS_R,
                                            /* 9 */ DATA_FIELDS_BAF, /* 10 */ DATA_FIELDS_LRR,
                                            /* 11 */ GENOTYPE_FIELDS_A1_FOR,
                                            /* 12 */ GENOTYPE_FIELDS_A2_FOR,
                                            /* 13 */ GENOTYPE_FIELDS_A1_AB,
                                            /* 14 */ GENOTYPE_FIELDS_A2_AB,};

  private void parseColumnsBestGuess() throws Elision {
    int[] indices = ext.indexFactors(LOOKUP, getCols(), false, true, false);
    if (indices[0] == -1) {
      throw new Elision("Error - missing SNP ID column");
    }
    setColSnpIdent(indices[0]);
    setColSampleIdent(indices[1]);
    // col_sampleIndex = indices[2];
    setColGC(indices[2]);
    setColXRaw(indices[3]);
    setColYRaw(indices[4]);
    setColX(indices[5]);
    setColY(indices[6]);
    setColTheta(indices[7]);
    setColR(indices[8]);
    setColBAF(indices[9]);
    setColLRR(indices[10]);
    setColGeno1(indices[11]);
    setColGeno2(indices[12]);
    setColGenoAB1(indices[13]);
    setColGenoAB2(indices[14]);
  }

  private void processLine(String[] parts) {
    if ("[Header]".equals(parts[0])) {
      return;
    }
    if ("File".equals(parts[0])) {
      String[] fileParts = parts[parts.length - 1].split(" of ");
      if (ext.isValidInteger(fileParts[0])) {
        currFile = Integer.parseInt(fileParts[0]);
      } else {
        // TODO error
      }
      if (ext.isValidInteger(fileParts[1])) {
        numFiles = Integer.parseInt(fileParts[1]);
      } else {
        // TODO error
      }
    }
    if ("Total Samples".equals(parts[0])) {
      if (ext.isValidInteger(parts[parts.length - 1])) {
        totalSamples = Integer.parseInt(parts[parts.length - 1]);
      } else {
        // TODO error
      }
    }
    if ("Num Samples".equals(parts[0])) {
      if (ext.isValidInteger(parts[parts.length - 1])) {
        numSamples = Integer.parseInt(parts[parts.length - 1]);
      } else {
        // TODO error
      }
    }
    if ("Total SNPs".equals(parts[0])) {
      if (ext.isValidInteger(parts[parts.length - 1])) {
        totalSnps = Integer.parseInt(parts[parts.length - 1]);
      } else {
        // TODO error
      }
    }
    if ("Num SNPs".equals(parts[0])) {
      if (ext.isValidInteger(parts[parts.length - 1])) {
        numSnps = Integer.parseInt(parts[parts.length - 1]);
      } else {
        // TODO error
      }
    }
    if ("Content".equals(parts[0])) {
      content = parts[parts.length - 1];
    }
    if ("Processing Date".equals(parts[0])) {
      processingDate = parts[parts.length - 1];
    }
    if ("GSGT Version".equals(parts[0])) {
      gsgtVersion = parts[parts.length - 1];
    }
  }

  public static HashMap<String, SourceFileHeaderData> validate(final String rawDir,
                                                               final String ext,
                                                               boolean fullValidation, Logger log,
                                                               Optional<JProgressBar> progressBar) {
    String dir = rawDir.endsWith("/")
                 || rawDir.endsWith("\\") ? rawDir
                                          : org.pankratzlab.common.ext.verifyDirFormat(rawDir);
    String[] possibleFiles = (new File(dir)).list(new FilenameFilter() {

      @Override
      public boolean accept(File dir, String name) {
        return name.endsWith(ext);
      }
    });
    // No files found - not a valid project
    if (possibleFiles.length == 0) {
      log.report("Project validation failed: no project files found discovered.");
      return null;
    }

    if (progressBar.isPresent()) {
      progressBar.get().setVisible(true);
      progressBar.get().setMinimum(0);
      progressBar.get().setMaximum(possibleFiles.length);
      progressBar.get().setString(null);
      progressBar.get().setStringPainted(true);
    }

    boolean valid = false;
    HashMap<String, SourceFileHeaderData> headers = null;
    int progCnt = 0;
    try {
      headers = new HashMap<>();
      SourceFileHeaderData exemplar = null;
      for (String possFile : possibleFiles) {
        SourceFileHeaderData frhd;
        if (!fullValidation) {
          if (exemplar == null) {
            frhd = SourceFileHeaderData.parseHeader(dir + possFile, log);
            exemplar = frhd;
          } else {
            frhd = exemplar;
          }
        } else {
          frhd = SourceFileHeaderData.parseHeader(dir + possFile, log);
        }
        headers.put(possFile, frhd);
        if (progressBar.isPresent()) {
          progressBar.get().setValue(++progCnt);
        }
      }
      if (fullValidation) {
        if (progressBar.isPresent()) {
          progressBar.get().setIndeterminate(true);
          progressBar.get().setString("Verifying...");
        }
        String error = doFullValidation(headers, log);
        if (error != null) {
          JOptionPane.showMessageDialog(null, error, "Error", JOptionPane.ERROR_MESSAGE);
          throw new Elision(error);
        }
      }
      if (progressBar.isPresent()) {
        progressBar.get().setVisible(false);
      }
      valid = true;
    } catch (Elision e) {
      log.reportError(e.getMessage());
    } catch (IOException e) {
      log.reportException(e);
    }

    return valid ? headers : null;
  }

  public static String doFullValidation(HashMap<String, SourceFileHeaderData> headers, Logger log) {
    int cnt = headers.size();
    HashMap<Integer, ArrayList<String>> totSnpsSet = new HashMap<>();
    HashMap<Integer, ArrayList<String>> headerLineIndex = new HashMap<>();
    HashMap<Integer, ArrayList<String>> sampleID = new HashMap<>();
    HashMap<Integer, ArrayList<String>> snpIndex = new HashMap<>();
    HashMap<Integer, ArrayList<String>> genoAB1 = new HashMap<>();
    HashMap<Integer, ArrayList<String>> genoAB2 = new HashMap<>();
    HashMap<Integer, ArrayList<String>> genoForward1 = new HashMap<>();
    HashMap<Integer, ArrayList<String>> genoForward2 = new HashMap<>();
    HashMap<Integer, ArrayList<String>> x = new HashMap<>();
    HashMap<Integer, ArrayList<String>> y = new HashMap<>();
    HashMap<Integer, ArrayList<String>> theta = new HashMap<>();
    HashMap<Integer, ArrayList<String>> r = new HashMap<>();
    HashMap<Integer, ArrayList<String>> xRaw = new HashMap<>();
    HashMap<Integer, ArrayList<String>> yRaw = new HashMap<>();
    HashMap<Integer, ArrayList<String>> baf = new HashMap<>();
    HashMap<Integer, ArrayList<String>> lrr = new HashMap<>();
    HashMap<Integer, ArrayList<String>> gc = new HashMap<>();
    for (java.util.Map.Entry<String, SourceFileHeaderData> entry : headers.entrySet()) {
      SourceFileHeaderData headerData = entry.getValue();
      if (headerData.getNumFiles() == -1) {
        headerData.setNumFiles(cnt);
      } else if (headerData.getNumFiles() != cnt) {
        return "Number of Files listed in Source File {" + entry.getKey()
               + "} does not equal the number of headers needing validation.  Please check source directory and extension and try again.";
      }
      ArrayList<String> files = totSnpsSet.get(headerData.getTotalSnps());
      if (files == null) {
        files = new ArrayList<>();
        totSnpsSet.put(headerData.getTotalSnps(), files);
      }
      files.add(entry.getKey());
      files = headerLineIndex.get(headerData.getColumnHeaderLineIndex());
      if (files == null) {
        files = new ArrayList<>();
        headerLineIndex.put(headerData.getColumnHeaderLineIndex(), files);
      }
      files.add(entry.getKey());
      files = sampleID.get(headerData.getColSampleIdent());
      if (files == null) {
        files = new ArrayList<>();
        sampleID.put(headerData.getColSampleIdent(), files);
      }
      files.add(entry.getKey());
      files = snpIndex.get(headerData.getColSnpIdent());
      if (files == null) {
        files = new ArrayList<>();
        snpIndex.put(headerData.getColSnpIdent(), files);
      }
      files.add(entry.getKey());
      files = genoAB1.get(headerData.getColGenoAB1());
      if (files == null) {
        files = new ArrayList<>();
        genoAB1.put(headerData.getColGenoAB1(), files);
      }
      files.add(entry.getKey());
      files = genoAB2.get(headerData.getColGenoAB2());
      if (files == null) {
        files = new ArrayList<>();
        genoAB2.put(headerData.getColGenoAB2(), files);
      }
      files.add(entry.getKey());
      files = genoForward1.get(headerData.getColGeno1());
      if (files == null) {
        files = new ArrayList<>();
        genoForward1.put(headerData.getColGeno1(), files);
      }
      files.add(entry.getKey());
      files = genoForward2.get(headerData.getColGeno2());
      if (files == null) {
        files = new ArrayList<>();
        genoForward2.put(headerData.getColGeno2(), files);
      }
      files.add(entry.getKey());
      files = x.get(headerData.getColX());
      if (files == null) {
        files = new ArrayList<>();
        x.put(headerData.getColX(), files);
      }
      files.add(entry.getKey());
      files = y.get(headerData.getColY());
      if (files == null) {
        files = new ArrayList<>();
        y.put(headerData.getColY(), files);
      }
      files.add(entry.getKey());
      files = theta.get(headerData.getColTheta());
      if (files == null) {
        files = new ArrayList<>();
        theta.put(headerData.getColTheta(), files);
      }
      files.add(entry.getKey());
      files = r.get(headerData.getColR());
      if (files == null) {
        files = new ArrayList<>();
        r.put(headerData.getColR(), files);
      }
      files.add(entry.getKey());
      files = xRaw.get(headerData.getColXRaw());
      if (files == null) {
        files = new ArrayList<>();
        xRaw.put(headerData.getColXRaw(), files);
      }
      files.add(entry.getKey());
      files = yRaw.get(headerData.getColYRaw());
      if (files == null) {
        files = new ArrayList<>();
        yRaw.put(headerData.getColYRaw(), files);
      }
      files.add(entry.getKey());
      files = baf.get(headerData.getColBAF());
      if (files == null) {
        files = new ArrayList<>();
        baf.put(headerData.getColBAF(), files);
      }
      files.add(entry.getKey());
      files = lrr.get(headerData.getColLRR());
      if (files == null) {
        files = new ArrayList<>();
        lrr.put(headerData.getColLRR(), files);
      }
      files.add(entry.getKey());
      files = gc.get(headerData.getColGC());
      if (files == null) {
        files = new ArrayList<>();
        gc.put(headerData.getColGC(), files);
      }
      files.add(entry.getKey());
    }
    int numErrors = 0;
    ArrayList<String> errorMsgs = new ArrayList<>();

    String error;
    error = checkErrors(totSnpsSet, "Total SNPs");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(headerLineIndex, "Index of header line");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(sampleID, "Sample ID column index");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(snpIndex, "SNP column index");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(genoAB1, "AB Genotype column 1");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(genoAB2, "AB Genotype column 2");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(genoForward1, "Forward Genotype column 1");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(genoForward2, "Forward Genotype column 2");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(x, "X column");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(y, "Y column");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(theta, "Theta column");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(r, "R column");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(xRaw, "X Raw column");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(yRaw, "Y Raw column");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(baf, "B Allele Freq column");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(lrr, "Log R Ratio column");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }
    error = checkErrors(gc, "GC column");
    if (!"".equals(error)) {
      numErrors++;
      errorMsgs.add(error);
    }

    if (numErrors > 0 || !errorMsgs.isEmpty()) {
      if (!errorMsgs.isEmpty()) {
        for (String s : errorMsgs) {
          log.reportError(s);
        }
      }
      return "Found " + numErrors
             + " data or column index mismatches among source files.  Please check log for more details";
    }

    return null;
  }

  private static String checkErrors(HashMap<Integer, ArrayList<String>> valueMapping,
                                    String string) {
    if (valueMapping.size() == 1) {
      return "";
    }
    StringBuilder sb = new StringBuilder("Mismatch in ").append(string).append(" between ")
                                                        .append(valueMapping.size())
                                                        .append(" sets of values: {");

    ArrayList<Integer> values = new ArrayList<>(valueMapping.keySet());
    Collections.sort(values);
    Collections.reverse(values);

    for (int i = 0, cnt = values.size(), cnt1 = cnt - 1; i < cnt; i++) {
      sb.append(values.get(i));
      if (i < cnt1) {
        sb.append(", ");
      }
    }

    sb.append("} with {");
    for (int i = 0, cnt = values.size(), cnt1 = cnt - 1; i < cnt; i++) {
      sb.append(valueMapping.get(values.get(i)).size());
      if (i < cnt1) {
        sb.append(", ");
      }
    }
    sb.append("} file(s) for each value, respectively.");
    if (valueMapping.size() == 2) {
      sb.append(" The file(s) in the second set are: ");
      ArrayList<String> files = valueMapping.get(values.get(values.size() - 1));
      for (int i = 0; i < files.size(); i++) {
        sb.append(files.get(i));
        if (i < files.size() - 1) {
          sb.append(", ");
        }
      }
      sb.append(".");
    }
    return sb.toString();
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = "D:/data/ny_registry/gedi_gwas/00src/";// null;
    String ext = ".csv.gz";
    String log = null;

    String usage = "\n" + "cnv.filesys.FinalReportHeaderData requires 2 argument\n"
                   + "   (1) Directory of FinalReport files (i.e. dir=" + dir + " (default))\n"
                   + "   (2) Extension of FinalReport files (i.e. ext=" + ext + " (default))\n";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("ext=")) {
        ext = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        log = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0 || dir == null) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      validate(dir, ext, true, log == null ? new Logger() : new Logger(log), Optional.empty());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public String getSourceFileDelimiter() {
    return delimiter;
  }

  /**
   * @return the columnHeaderLineIndex
   */
  public int getColumnHeaderLineIndex() {
    return columnHeaderLineIndex;
  }

  /**
   * @return the colSampleIdent
   */
  public int getColSampleIdent() {
    return colSampleIdent;
  }

  /**
   * @param colSampleIdent the colSampleIdent to set
   */
  public void setColSampleIdent(int colSampleIdent) {
    this.colSampleIdent = colSampleIdent;
  }

  /**
   * @return the colSnpIdent
   */
  public int getColSnpIdent() {
    return colSnpIdent;
  }

  /**
   * @param colSnpIdent the colSnpIdent to set
   */
  public void setColSnpIdent(int colSnpIdent) {
    this.colSnpIdent = colSnpIdent;
  }

  /**
   * @return the colGenoAB1
   */
  public int getColGenoAB1() {
    return colGenoAB1;
  }

  /**
   * @param colGenoAB1 the colGenoAB1 to set
   */
  public void setColGenoAB1(int colGenoAB1) {
    this.colGenoAB1 = colGenoAB1;
  }

  /**
   * @return the colGenoAB2
   */
  public int getColGenoAB2() {
    return colGenoAB2;
  }

  /**
   * @param colGenoAB2 the colGenoAB2 to set
   */
  public void setColGenoAB2(int colGenoAB2) {
    this.colGenoAB2 = colGenoAB2;
  }

  /**
   * @return the colGeno1
   */
  public int getColGeno1() {
    return colGeno1;
  }

  /**
   * @param colGeno1 the colGeno1 to set
   */
  public void setColGeno1(int colGeno1) {
    this.colGeno1 = colGeno1;
  }

  /**
   * @return the colGeno2
   */
  public int getColGeno2() {
    return colGeno2;
  }

  /**
   * @param colGeno2 the colGeno2 to set
   */
  public void setColGeno2(int colGeno2) {
    this.colGeno2 = colGeno2;
  }

  /**
   * @return the colX
   */
  public int getColX() {
    return colX;
  }

  /**
   * @param colX the colX to set
   */
  public void setColX(int colX) {
    this.colX = colX;
  }

  /**
   * @return the colY
   */
  public int getColY() {
    return colY;
  }

  /**
   * @param colY the colY to set
   */
  public void setColY(int colY) {
    this.colY = colY;
  }

  /**
   * @return the colTheta
   */
  public int getColTheta() {
    return colTheta;
  }

  /**
   * @param colTheta the colTheta to set
   */
  public void setColTheta(int colTheta) {
    this.colTheta = colTheta;
  }

  /**
   * @return the colR
   */
  public int getColR() {
    return colR;
  }

  /**
   * @param colR the colR to set
   */
  public void setColR(int colR) {
    this.colR = colR;
  }

  /**
   * @return the colXRaw
   */
  public int getColXRaw() {
    return colXRaw;
  }

  /**
   * @param colXRaw the colXRaw to set
   */
  public void setColXRaw(int colXRaw) {
    this.colXRaw = colXRaw;
  }

  /**
   * @return the colYRaw
   */
  public int getColYRaw() {
    return colYRaw;
  }

  /**
   * @param colYRaw the colYRaw to set
   */
  public void setColYRaw(int colYRaw) {
    this.colYRaw = colYRaw;
  }

  /**
   * @return the colBAF
   */
  public int getColBAF() {
    return colBAF;
  }

  /**
   * @param colBAF the colBAF to set
   */
  public void setColBAF(int colBAF) {
    this.colBAF = colBAF;
  }

  /**
   * @return the colLRR
   */
  public int getColLRR() {
    return colLRR;
  }

  /**
   * @param colLRR the colLRR to set
   */
  public void setColLRR(int colLRR) {
    this.colLRR = colLRR;
  }

  /**
   * @return the colGC
   */
  public int getColGC() {
    return colGC;
  }

  /**
   * @param colGC the colGC to set
   */
  public void setColGC(int colGC) {
    this.colGC = colGC;
  }

  /**
   * @return the cols
   */
  public String[] getCols() {
    return cols;
  }

  /**
   * @return the gsgtVersion
   */
  String getGsgtVersion() {
    return gsgtVersion;
  }

  /**
   * @return the processingDate
   */
  String getProcessingDate() {
    return processingDate;
  }

  /**
   * @return the content
   */
  String getContent() {
    return content;
  }

  /**
   * @return the numSnps
   */
  int getNumSnps() {
    return numSnps;
  }

  /**
   * @return the totalSnps
   */
  int getTotalSnps() {
    return totalSnps;
  }

  /**
   * @return the numSamples
   */
  int getNumSamples() {
    return numSamples;
  }

  /**
   * @return the totalSamples
   */
  int getTotalSamples() {
    return totalSamples;
  }

  /**
   * @return the numFiles
   */
  int getNumFiles() {
    return numFiles;
  }

  /**
   * @param numFiles the numFiles to set
   */
  private void setNumFiles(int numFiles) {
    this.numFiles = numFiles;
  }

  /**
   * @return the currFile
   */
  int getCurrFile() {
    return currFile;
  }
}
