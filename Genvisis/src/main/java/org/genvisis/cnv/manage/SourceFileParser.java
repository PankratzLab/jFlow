// "AB_lookup.dat" is necessary if the files do not contain {"Allele1 - AB"}/{"Allele2 - AB}
package org.genvisis.cnv.manage;

import java.awt.GraphicsEnvironment;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import javax.swing.JOptionPane;

import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.filesys.SourceFileHeaderData;
import org.pankratzlab.common.Aliases;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CountHash;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.Sort;
import org.pankratzlab.common.WorkerTrain;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;
import org.pankratzlab.common.ext;

import com.google.common.collect.Lists;

public class SourceFileParser implements Runnable {

  private static final String IDS_CHANGED_FILEROOT = "FYI_IDS_WERE_CHANGED";
  public static final String[][] SNP_HEADER_OPTIONS = {Aliases.MARKER_NAMES};
  protected static final String[][] SNP_TABLE_FIELDS = {SNP_HEADER_OPTIONS[0], Aliases.CHRS,
                                                        Aliases.POSITIONS};
  protected static final String[] DELIMITERS = {",", "\t", " "};
  protected static final String[] DELIMITER_DESCRIPTIONS = {"COMMA", "TAB", "SPACE"};
  public static final String OVERWRITE_OPTION_FILE = ".overwrite_option";
  protected static final String CANCEL_OPTION_FILE = ".cancel_option";
  protected static final String HOLD_OPTION_FILE = ".hold_option";
  public static final String FILENAME_AS_ID_OPTION = "[FILENAME_ROOT]";

  private final Project proj;
  private final String[] files;
  private final String[] markerNames;
  private final int[] keysKeys;
  private final long fingerprint;
  private final char[][] abLookup;
  private final Hashtable<String, String> fixes;
  private final long timeBegan;
  private final int threadId;
  private final String delimiter;
  private boolean splitAB;

  public SourceFileParser(Project proj, String[] files, String[] markerNames, int[] keysKeys,
                          char[][] abLookup, String delimiter, long fingerprint,
                          Hashtable<String, String> fixes, long timeBegan) {
    this(proj, files, markerNames, keysKeys, abLookup, delimiter, fingerprint, fixes, timeBegan,
         -1);
  }

  public SourceFileParser(Project proj, String[] files, String[] markerNames, int[] keysKeys,
                          char[][] abLookup, String delimiter, long fingerprint,
                          Hashtable<String, String> fixes, long timeBegan, int threadId) {
    this.proj = proj;
    this.files = files;
    this.markerNames = markerNames;
    this.keysKeys = keysKeys;
    this.abLookup = abLookup;
    this.delimiter = delimiter;
    this.fingerprint = fingerprint;
    this.fixes = fixes;
    this.timeBegan = timeBegan;
    this.threadId = threadId;
    splitAB = false;
  }

  @Override
  public void run() {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    int count, snpIndex, sampIndex, key;
    boolean parseAtAt;

    Sample samp;
    String sampleName, /* temp, */filename, trav, idHeader;
    float[][] data;
    byte[][] genotypes;
    Hashtable<String, Float> allOutliers;
    Map<String, SourceFileHeaderData> headers;
    Logger log;

    log = proj.getLog();
    idHeader = proj.ID_HEADER.getValue();
    allOutliers = new Hashtable<>();
    headers = proj.getSourceFileHeaders(true);
    try {
      for (int i = 0; i < files.length; i++) {

        if (Thread.currentThread().isInterrupted()) {
          return;
        }
        if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true)
                     + SourceFileParser.CANCEL_OPTION_FILE).exists()) {
          return;
        }
        try {
          log.reportTime("\tSource file parsing thread " + threadId + ": processing " + (i + 1)
                         + " of " + files.length + " -- " + files[i]);

          SourceFileHeaderData headerData = headers.get(files[i]);
          if (headerData == null) {
            // TODO error, missing header for source files
          }

          reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true)
                                              + files[i]);
          for (int k = 0; k < headerData.getColumnHeaderLineIndex() + 1; k++) {
            reader.readLine();
          }
          if (idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)) {
            sampIndex = -7;
          } else {
            sampIndex = headerData.getColSampleIdent();
          }
          snpIndex = headerData.getColSnpIdent();

          if (headerData.getColX() == -1 || headerData.getColY() == -1) {
            log.reportError("Error - source file " + files[i]
                            + " format was not consistent! At the very least the files need to contain "
                            + ArrayUtils.toStr(Sample.DATA_FIELDS[3], "/") + " and "
                            + ArrayUtils.toStr(Sample.DATA_FIELDS[4], "/"));
            return;
          }
          final boolean ignoreAB;
          if ((headerData.getColGenoAB1() == -1 || headerData.getColGenoAB2() == -1)
              && abLookup == null) {
            ignoreAB = true;
          } else {
            ignoreAB = false;
          }

          sampleName = "";
          data = new float[][] {headerData.getColGC() == -1 ? null : new float[markerNames.length],
                                headerData.getColXRaw() == -1 ? null
                                                              : new float[markerNames.length],
                                headerData.getColYRaw() == -1 ? null
                                                              : new float[markerNames.length],
                                headerData.getColX() == -1 ? null : new float[markerNames.length],
                                headerData.getColY() == -1 ? null : new float[markerNames.length],
                                headerData.getColTheta() == -1 ? null
                                                               : new float[markerNames.length],
                                headerData.getColR() == -1 ? null : new float[markerNames.length],
                                headerData.getColBAF() == -1 ? null : new float[markerNames.length],
                                headerData.getColLRR() == -1 ? null
                                                             : new float[markerNames.length],};
          genotypes = new byte[][] {ArrayUtils.byteArray(markerNames.length, (byte) 0),
                                    ignoreAB ? null
                                             : ArrayUtils.byteArray(markerNames.length, (byte) -1)};
          count = 0;
          parseAtAt = proj.PARSE_AT_AT_SYMBOL.getValue();
          String tmp;
          while ((tmp = reader.readLine()) != null) {
            line = tmp.split(delimiter, -1);
            if (idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)) {
              trav = files[i].substring(0,
                                        files[i].indexOf(proj.SOURCE_FILENAME_EXTENSION.getValue()));
            } else {
              if (parseAtAt && line[sampIndex].indexOf("@") == -1) {
                log.reportError("Error - source file " + files[i] + " header " + idHeader + " '"
                                + line[sampIndex] + "' did not contain an @ sample");
                parseAtAt = false;
              }
              trav = parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@"))
                               : line[sampIndex];
            }
            if (count == 0) {
              sampleName = trav;
            } else if (!trav.equals(sampleName)) {
              log.reportError("Found more than one ID in file " + files[i] + "(found " + trav
                              + ", expecting " + sampleName + ")");
              return;
            } else if (count > markerNames.length) {
              log.reportError("Error - expecting only " + markerNames.length
                              + " markers and found more than that in file " + files[i]);
              return;
            } else if (!markerNames[count].equals(line[snpIndex])
                       && !markerNames[count].startsWith("Blank")) {
              log.reportError("Found " + line[snpIndex] + " at marker #" + (count + 1) + " in file "
                              + files[i] + "; expecting " + markerNames[count]);
              return;
            }
            key = keysKeys[count];

            int[] indCols = {headerData.getColGC(), headerData.getColXRaw(),
                             headerData.getColYRaw(), headerData.getColX(), headerData.getColY(),
                             headerData.getColTheta(), headerData.getColR(), headerData.getColBAF(),
                             headerData.getColLRR(),};

            int ind = 0;
            try {
              for (int i1 = 0; i1 < indCols.length; i1++) {
                ind = i1;
                if (indCols[ind] != -1) {
                  data[ind][key] = ext.isMissingValue(line[indCols[ind]]) ? Float.NaN
                                                                          : Float.parseFloat(line[indCols[ind]]);
                }
              }
            } catch (NumberFormatException nfe) {
              log.reportError("Error - failed in file " + files[i] + " at line " + key
                              + " to parse '" + line[indCols[ind]] + "' into a valid "
                              + ArrayUtils.toStr(Sample.DATA_FIELDS[ind], "/"));
              return;
            } catch (Exception e) {
              log.reportError("Some other exception in file " + files[i] + ":");
              log.reportException(e);
              return;
            }

            if (proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6
                || proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6_CN/*
                                                                   * || proj.ARRAY_TYPE.getValue()
                                                                   * == ARRAY.DBGAP
                                                                   */) {
              data[0][key] = 1 - data[0][key];
            }
            if (proj.XY_SCALE_FACTOR.getValue() != 1) {
              data[3][key] = data[3][key] / proj.XY_SCALE_FACTOR.getValue().floatValue();
              data[4][key] = data[4][key] / proj.XY_SCALE_FACTOR.getValue().floatValue();
            }
            if (headerData.getColGeno1() != -1 && headerData.getColGeno2() != -1) {
              if (line[headerData.getColGeno1()].length() > 1) {
                if (i == 0 && !splitAB) {
                  log.reportTimeInfo("Detected genotype calls for each allele are in the same column (Such as in Affymetrix .chp format) ");
                  splitAB = true;
                }
              } else if (splitAB) {
                log.reportError("Detected previously that genotype calls should be split, but the calls on line "
                                + key + " --> {" + ArrayUtils.toStr(line)
                                + "} did not.  Parsing will fail for file " + files[i] + ".");
                return;
              }

              if (splitAB) {
                genotypes[0][key] = (byte) ext.indexOfStr(line[headerData.getColGeno1()],
                                                          Sample.ALLELE_PAIRS);
              } else {
                genotypes[0][key] = (byte) ext.indexOfStr(line[headerData.getColGeno1()]
                                                          + line[headerData.getColGeno2()],
                                                          Sample.ALLELE_PAIRS);
              }

              if (genotypes[0][key] == -1) {
                if (proj.getArrayType() == ARRAY.ILLUMINA) {

                  if (ext.indexOfStr(line[headerData.getColGeno1()]
                                     + line[headerData.getColGeno2()], Sample.ALT_NULLS) == -1) {
                    log.reportError("Error - failed to lookup genotype ("
                                    + line[headerData.getColGeno1()]
                                    + line[headerData.getColGeno2()] + ") for marker "
                                    + markerNames[count] + " of sample " + files[i]
                                    + "; setting to missing");
                  }
                  genotypes[0][key] = 0;
                } else if (proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6
                           || proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6_CN) {
                  genotypes[0][key] = (byte) ext.indexOfStr(line[headerData.getColGeno1()],
                                                            Sample.AB_PAIRS);
                  if (genotypes[0][key] == -1
                      && ext.indexOfStr(line[headerData.getColGeno1()], Sample.ALT_NULLS) == -1) {
                    log.reportError("Error - failed to lookup genotype ("
                                    + line[headerData.getColGeno1()] + ") for marker "
                                    + markerNames[count] + " of sample " + files[i]
                                    + "; setting to missing");

                  }
                  // Affy matrix format
                  // // does not use
                  // // ALLELE_PAIRS
                  genotypes[0][key] = 0;
                }
              }
            }
            if (headerData.getColGenoAB1() != -1 && headerData.getColGenoAB2() != -1) {
              if (ignoreAB) {
                // do nothing, will need to use these files to determine AB lookup table
              } else if (abLookup == null || headerData.getColGeno1() == -1
                         || headerData.getColGeno2() == -1
                         || proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6
                         || proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6_CN) {// FORCE AFFY HERE
                if (abLookup != null) {
                  log.reportTimeWarning("ABLookup data provided to source file parser but forward genotypes aren't present in the source file.  AB genotypes will be parsed from the source file instead.");
                }
                if (splitAB) {
                  genotypes[1][key] = (byte) ext.indexOfStr(line[headerData.getColGenoAB1()],
                                                            Sample.AB_PAIRS);
                } else {
                  genotypes[1][key] = (byte) ext.indexOfStr(line[headerData.getColGenoAB1()]
                                                            + line[headerData.getColGenoAB2()],
                                                            Sample.AB_PAIRS);
                }
              } else {
                if (genotypes[0][key] == 0) {
                  genotypes[1][key] = -1;
                } else {
                  genotypes[1][key] = 0;
                  if (line[headerData.getColGeno1()].charAt(0) == abLookup[count][1]) {
                    genotypes[1][key]++;
                  } else if (line[headerData.getColGeno1()].charAt(0) != abLookup[count][0]) {
                    log.reportError("Error - alleles for individual '"
                                    + (sampIndex < 0 ? trav : line[sampIndex]) + "' ("
                                    + line[headerData.getColGeno1()] + "/"
                                    + line[headerData.getColGeno2()]
                                    + ") do not match up with the defined AB lookup alleles ("
                                    + abLookup[count][0] + "/" + abLookup[count][1]
                                    + ") for marker " + markerNames[count]);
                  }
                  if (line[headerData.getColGeno2()].charAt(0) == abLookup[count][1]) {
                    genotypes[1][key]++;
                  } else if (line[headerData.getColGeno2()].charAt(0) != abLookup[count][0]) {
                    log.reportError("Error - alleles for individual '"
                                    + (sampIndex < 0 ? trav : line[sampIndex]) + "' ("
                                    + line[headerData.getColGeno1()] + "/"
                                    + line[headerData.getColGeno2()]
                                    + ") do not match up with the defined AB lookup alleles ("
                                    + abLookup[count][0] + "/" + abLookup[count][1]
                                    + ") for marker " + markerNames[count]);
                  }
                }
              }
            }
            count++;
          }
          reader.close();
          if (count != markerNames.length) {
            log.reportError("Error - expecting " + markerNames.length + " markers and only found "
                            + count + " in file " + files[i]);
            return;
          }

          if (fixes.containsKey(sampleName)) {
            sampleName = fixes.get(sampleName);
          }

          writer = null;
          trav = ext.replaceWithLinuxSafeCharacters(sampleName, true);
          if (!trav.equals(sampleName)) {
            if (writer == null) {
              String idsChangedFile = proj.PROJECT_DIRECTORY.getValue() + IDS_CHANGED_FILEROOT
                                      + threadId + ".txt";
              try {
                writer = Files.openAppropriateWriter(idsChangedFile, true);
                if (new File(idsChangedFile).length() == 0) {
                  writer.println("The following IDs were changed so that spaces are removed and so that they could be used as valid filenames:");
                }
              } catch (Exception e) {
                log.reportError("Error writing to " + idsChangedFile);
                log.reportException(e);
              }
            }
            writer.println(sampleName + "\t" + trav);
            writer.close();
            sampleName = trav;
          }

          filename = SourceFileParser.determineFilename(proj.SAMPLE_DIRECTORY.getValue(true, true),
                                                        sampleName, timeBegan, log);
          if (filename == null) {
            log.reportError("Could not determine file name for sample " + files[i]);
            return;
          }

          samp = new Sample(sampleName, fingerprint, data, genotypes,
                            proj.getArrayType().getCanXYBeNegative());
          samp.saveToRandomAccessFile(filename, allOutliers, sampleName);
        } catch (FileNotFoundException fnfe) {
          log.reportError("Error: file \"" + files[i] + "\" not found in current directory");
          return;
        } catch (IOException ioe) {
          log.reportError("Error reading file \"" + files[i] + "\"");
          return;
        }
      }

      if (allOutliers.size() > 0) {
        String outliersFile;
        if (threadId >= 0) {
          outliersFile = proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + threadId
                         + ".ser";
          if (new File(outliersFile).exists()) {
            log.reportError("Error - the following file already exists: " + outliersFile);
            return;
          } else {
            SerializedFiles.writeSerial(allOutliers, outliersFile);
          }
        } else {
          outliersFile = proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers0.ser";
          if (new File(outliersFile).exists()) {
            log.reportError("Error - the following file already exists: " + outliersFile);
            return;
          } else {
            SerializedFiles.writeSerial(allOutliers, outliersFile);
          }
        }
      }

    } catch (Exception e) {
      // TODO doesn't catch OOMs
      log.reportException(e);
    }

    // request release of resources
    log.reportTime("\tThread " + threadId + " cleaning up...");
    System.gc();
  }

  private static String buildColumnAssignmentsLogOutput(SourceFileHeaderData headerData) {
    StringBuilder logOutput = new StringBuilder();
    logOutput.append("Column name assignments for data import:\n");
    logOutput.append("GC: ")
             .append(headerData.getColGC() == -1 ? "[missing]"
                                                 : headerData.getCols()[headerData.getColGC()])
             .append("\n");
    logOutput.append("XRaw: ")
             .append(headerData.getColXRaw() == -1 ? "[missing]"
                                                   : headerData.getCols()[headerData.getColXRaw()])
             .append("\n");
    logOutput.append("YRaw: ")
             .append(headerData.getColYRaw() == -1 ? "[missing]"
                                                   : headerData.getCols()[headerData.getColYRaw()])
             .append("\n");
    logOutput.append("X: ")
             .append(headerData.getColX() == -1 ? "[missing]"
                                                : headerData.getCols()[headerData.getColX()])
             .append("\n");
    logOutput.append("Y: ")
             .append(headerData.getColY() == -1 ? "[missing]"
                                                : headerData.getCols()[headerData.getColY()])
             .append("\n");
    logOutput.append("Theta: ")
             .append(headerData.getColTheta() == -1 ? "[missing]"
                                                    : headerData.getCols()[headerData.getColTheta()])
             .append("\n");
    logOutput.append("R: ")
             .append(headerData.getColR() == -1 ? "[missing]"
                                                : headerData.getCols()[headerData.getColR()])
             .append("\n");
    logOutput.append("BAF: ")
             .append(headerData.getColBAF() == -1 ? "[missing]"
                                                  : headerData.getCols()[headerData.getColBAF()])
             .append("\n");
    logOutput.append("LRR: ")
             .append(headerData.getColLRR() == -1 ? "[missing]"
                                                  : headerData.getCols()[headerData.getColLRR()])
             .append("\n");
    logOutput.append("Geno1: ")
             .append(headerData.getColGeno1() == -1 ? "[missing]"
                                                    : headerData.getCols()[headerData.getColGeno1()])
             .append("\n");
    logOutput.append("Geno2: ")
             .append(headerData.getColGeno2() == -1 ? "[missing]"
                                                    : headerData.getCols()[headerData.getColGeno2()])
             .append("\n");
    logOutput.append("AB1: ")
             .append(headerData.getColGenoAB1() == -1 ? "[missing]"
                                                      : headerData.getCols()[headerData.getColGenoAB1()])
             .append("\n");
    logOutput.append("AB2: ")
             .append(headerData.getColGenoAB2() == -1 ? "[missing]"
                                                      : headerData.getCols()[headerData.getColGenoAB2()])
             .append("\n");
    return logOutput.toString();
  }

  public static char[][] getABLookup(boolean abLookupRequired, String[] markerNames, Project proj) {
    ABLookup abLookup;
    char[][] lookup;
    Logger log;

    log = proj.getLog();
    if (abLookupRequired && Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
      abLookup = new ABLookup(markerNames, proj.AB_LOOKUP_FILENAME.getValue(), true, false,
                              proj.getLog());
      lookup = abLookup.getLookup();
      if (lookup == null) {
        log.report("Warning - failed to provide columns \"" + Sample.GENOTYPE_FIELDS[2][0]
                   + "\" / \"" + Sample.GENOTYPE_FIELDS[3][0]
                   + "\" and the specificed AB_lookup file '"
                   + proj.getProperty(proj.AB_LOOKUP_FILENAME)
                   + "' does not exist; you'll need reconstruct the B allele for analysis");
      } else {
        abLookup.writeToFile(proj.PROJECT_DIRECTORY.getValue() + "checkAB.xln", proj.getLog());
      }
    } else {
      lookup = null;
    }

    return lookup;
  }

  public static String determineFilename(String dir, String sampleName, long timeBegan,
                                         Logger log) {
    String filename, trav;
    int version, versionToOverwrite;
    String[] overwriteOptions;
    int response;

    version = 0;
    versionToOverwrite = -1;
    do {
      trav = sampleName + (version == 0 ? "" : "." + version);
      version++;
      filename = dir + trav + Sample.SAMPLE_FILE_EXTENSION;
      if (new File(filename).exists() && new File(filename).lastModified() < timeBegan
          && versionToOverwrite == -1) {
        versionToOverwrite = version - 1;
      }
    } while (new File(filename).exists());

    filename = sampleName + (versionToOverwrite == 0 ? "" : "." + versionToOverwrite)
               + Sample.SAMPLE_FILE_EXTENSION;
    overwriteOptions = new String[] {"Rename new file " + trav + Sample.SAMPLE_FILE_EXTENSION,
                                     "Overwrite existing file " + filename,
                                     "Overwrite this and all future files", "Cancel parser"};

    if (versionToOverwrite != -1) {
      while (new File(dir + SourceFileParser.HOLD_OPTION_FILE).exists()) {
        try {
          Thread.sleep(500);
        } catch (InterruptedException ie) {}
      }
      if (new File(dir + SourceFileParser.CANCEL_OPTION_FILE).exists()) {
        return null;
      }

      Files.write("", dir + SourceFileParser.HOLD_OPTION_FILE);
      if (new File(dir + SourceFileParser.OVERWRITE_OPTION_FILE).exists()) {
        response = 1;
      } else {

        if (!Files.isWindows()) {
          log.reportError("Error - the same sample name '" + sampleName
                          + "' is being parsed again and the previous file existed before the current command began. "
                          + "This happens if you inadvertently restarted the parser or if the parser was interrupted and manually restarted. "
                          + "If you would like to start from scratch, the safest thing would be to cancel now and delete all files in the sample directory. "
                          + "To automatically overwrite when there is conflict like this, then create an empty file with the following name: "
                          + dir + SourceFileParser.OVERWRITE_OPTION_FILE);
          return null;
        }

        do {
          response = JOptionPane.showOptionDialog(null,
                                                  "Error - the same sample name '" + sampleName
                                                        + "' is being parsed again and the previous file existed before the current command began.\n"
                                                        + "This happens if you inadvertently restarted the parser or if the parser was interrupted and manually restarted.\n"
                                                        + "If you would like to start from scratch, the safest thing would be to cancel now and delete all files in the sample directory.\n"
                                                        + "What would you like to do?",
                                                  "What to do?", JOptionPane.DEFAULT_OPTION,
                                                  JOptionPane.QUESTION_MESSAGE, null,
                                                  overwriteOptions, overwriteOptions[0]);
        } while (response == -1);
      }
      new File(dir + SourceFileParser.HOLD_OPTION_FILE).delete();
      switch (response) {
        case 0:
          return dir + trav + Sample.SAMPLE_FILE_EXTENSION;
        case 2:
          Files.write("", dir + SourceFileParser.OVERWRITE_OPTION_FILE);
        case 1:
          return dir + filename;
        case 3:
          Files.write("", dir + SourceFileParser.CANCEL_OPTION_FILE);
          return null;
        default:
          JOptionPane.showMessageDialog(null, "Should be impossible to obtain this message ("
                                              + response + ")",
                                        "Error", JOptionPane.ERROR_MESSAGE);
          break;
      }
    }

    return dir + trav + Sample.SAMPLE_FILE_EXTENSION;
  }

  /**
   * Helper method to look up all files matching the source file extension in the source directory.
   */
  static String[] getSourceFiles(Project proj, Logger log) {
    log.report(ext.getTime() + "\tSearching for " + proj.SOURCE_FILENAME_EXTENSION.getValue()
               + " files in: " + proj.SOURCE_DIRECTORY.getValue(false, true));

    String[] files = Files.list(proj.SOURCE_DIRECTORY.getValue(false, true),
                                proj.SOURCE_FILENAME_EXTENSION.getValue());

    int unmatchedFileCount = new File(proj.SOURCE_DIRECTORY.getValue(false,
                                                                     true)).listFiles().length;
    unmatchedFileCount -= files.length;

    PSF.checkInterrupted();

    // remove known co-occurring samples
    if (proj.SOURCE_FILENAME_EXTENSION.getValue().equals(".csv")) {
      boolean[] use = new boolean[files.length];
      for (int i = 0; i < files.length; i++) {
        if (files[i].startsWith("Sample_Map.csv") || files[i].startsWith("SNP_Map.csv")) {
          use[i] = false;
        } else {
          use[i] = true;
        }
      }
      files = ArrayUtils.subArray(files, use);
    }

    if (files.length == 0) {
      log.reportError("Error - no files to parse; are you sure you have the right extension specified for your FinalReport files? It is currently set to \""
                      + proj.getProperty(proj.SOURCE_FILENAME_EXTENSION) + "\"");
    } else {
      if (unmatchedFileCount > 0) {
        log.reportError("Found " + unmatchedFileCount
                        + " file(s) in the source directory without extension: "
                        + proj.SOURCE_FILENAME_EXTENSION.getValue()
                        + " - please verify these were not supposed to be parsed.");
      }
      log.report("\t\tFound " + files.length + " file" + (files.length == 1 ? "" : "s")
                 + " to parse");
    }

    return files;
  }

  @SuppressWarnings("unchecked")
  public static int createFiles(Project proj, int numThreads) {
    BufferedReader reader;
    String[] line, markerNames, files;
    String idHeader, delimiter, temp, sampleName, trav;
    char[][] lookup;
    int[][] delimiterCounts;
    int[] keys, keysKeys, indices;
    int snpIndex, code, count, sampIndex, foundIDon, foundSNPon;
    long fingerprint;
    Vector<Vector<String>> fileCabinet;
    boolean abLookupRequired, complete, parseAtAt;
    long timeBegan;
    Hashtable<String, String> fixes;
    Hashtable<String, Float> allOutliers;
    Hashtable<String, Integer> markerNameHash;
    Vector<String> v;
    Logger log;
    AffyProcess affyProcess = null;

    log = proj.getLog();
    timeBegan = new Date().getTime();
    new File(proj.SAMPLE_DIRECTORY.getValue(true, true)
             + SourceFileParser.OVERWRITE_OPTION_FILE).delete();
    new File(proj.SAMPLE_DIRECTORY.getValue(true, true)
             + SourceFileParser.HOLD_OPTION_FILE).delete();
    new File(proj.SAMPLE_DIRECTORY.getValue(true, true)
             + SourceFileParser.CANCEL_OPTION_FILE).delete();

    if (!proj.SOURCE_DIRECTORY.getValue(false, true).equals("")
        && !new File(proj.SOURCE_DIRECTORY.getValue(false, true)).exists()) {
      log.reportError("Error - the Project source location is invalid: "
                      + proj.SOURCE_DIRECTORY.getValue(false, true));
      return 0;
    }

    PSF.checkInterrupted();
    if (!new File(proj.MARKER_POSITION_FILENAME.getValue(false, false)).exists()) {
      log.reportError("Error - missing markerPositions: "
                      + proj.MARKER_POSITION_FILENAME.getValue(false, false));
      return checkForSNP_Map(proj, log);
    }

    PSF.checkInterrupted();
    files = getSourceFiles(proj, log);
    if (files.length == 0) {
      return 0;
    }

    delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
    idHeader = proj.getProperty(proj.ID_HEADER);
    abLookupRequired = false;

    PSF.checkInterrupted();
    fixes = new Hashtable<>();
    if (new File(proj.PROJECT_DIRECTORY.getValue() + "fixes.dat").exists()) {
      log.report("Also found a 'fixes.dat' file in the project directory, which will be used to rename samples");
      fixes = HashVec.loadFileToHashString(proj.PROJECT_DIRECTORY.getValue() + "fixes.dat", false);
    } else {
      log.report("Did not find a 'fixes.dat' file; assuming there are no swapped samples to rename");
    }
    ARRAY array = proj.getArrayType();
    switch (array) {
      case AFFY_GW6:
        log.reportTimeWarning("Affymetrix confidence scores will be imported as (1-conf)");
        break;
      case AFFY_GW6_CN:
        PSF.checkInterrupted();
        log.reportTimeWarning("Affymetrix confidence scores will be imported as (1-conf)");
        log.reportTimeInfo("Initializing parser for array type " + array);
        affyProcess = new AffyProcess(proj,
                                      Files.toFullPaths(files,
                                                        proj.SOURCE_DIRECTORY.getValue(false,
                                                                                       true)),
                                      delimiter, log);
        affyProcess.matchCn();
        affyProcess.combineFirst();
        for (int i = 0; i < files.length; i++) {
          files[i] = ext.removeDirectoryInfo(affyProcess.getCombinedOutputFiles()[i]);
        }
        proj.setSourceFileHeaders(SourceFileHeaderData.validate(ext.parseDirectoryOfFile(affyProcess.getCombinedOutputFiles()[0]),
                                                                "." + proj.getArrayType() + ".tmp.gz",
                                                                true, proj.getLog(),
                                                                Optional.empty()));
        break;
      case ILLUMINA:
        break;
      // case NGS:
      // break;
      default:
        log.reportError("Invalid array type " + array);
        break;
    }

    proj.getSourceFileHeaders(true);

    PSF.checkInterrupted();
    try {
      reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[0]);
      log.report(ext.getTime() + "]\tFound appropriate reader for: "
                 + proj.SOURCE_DIRECTORY.getValue(false, true) + files[0]);

      count = 0;
      foundSNPon = -1;
      foundIDon = -2;
      log.report(ext.getTime() + "]\tSearching for header line...");
      do {
        line = reader.readLine().trim().replace("#Column header: ", "").split(delimiter, -1);
        count++;
        if (count < 20) {
          log.report(ArrayUtils.toStr(line), true, true, 11);
        }
        if (ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true,
                             false)[0] != -1) {
          foundSNPon = count;
        }
        if (idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)
            || ext.indexOfStr(idHeader, line) != -1) {
          foundIDon = count;
        }

        PSF.checkInterrupted();
      } while (reader.ready() && count < 1000 && foundIDon != foundSNPon);

      // If we reached the end of the file, it means that we didn't find the header we are looking
      // for
      // The most common cause of this is that the delimiter was misspecified
      // The following code checks all of the common delimiters (tab, comma, space) and determines
      // which one to use when it tries for a second time
      if (!reader.ready() || count == 1000) {

        reader.close();
        reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[0]);
        delimiterCounts = new int[SourceFileParser.DELIMITERS.length][count];
        for (int i = 0; i < count; i++) {
          temp = reader.readLine();
          for (int j = 0; j < SourceFileParser.DELIMITERS.length; j++) {
            delimiterCounts[j][i] = ext.countInstancesOf(temp, SourceFileParser.DELIMITERS[j]);
          }
        }
        delimiter = null;
        for (int j = 0; j < SourceFileParser.DELIMITERS.length; j++) {
          if (ArrayUtils.quantWithExtremeForTie(delimiterCounts[j], 0.5) > 4
              && ArrayUtils.quantWithExtremeForTie(delimiterCounts[j], 0.9)
                 - ArrayUtils.quantWithExtremeForTie(delimiterCounts[j], 0.1) == 0) {
            if (delimiter == null) {
              delimiter = SourceFileParser.DELIMITERS[j];
            } else {
              proj.message("Could not auto-detect the delimiter used in the Final Reports file: could be '"
                           + delimiter + "' or '" + SourceFileParser.DELIMITERS[j] + "'");
              return 0;
            }
          }
        }
        reader.close();

        if (delimiter == null) {
          if (foundSNPon == -1) {
            log.reportError("Could not find a header with the following tokens: "
                            + ArrayUtils.toStr(SourceFileParser.SNP_HEADER_OPTIONS[0], " / "));
          }
          if (foundIDon == -2) {
            log.reportError("Could not find a header with the selected id type: " + idHeader);
          }

          log.reportError("   Perhaps the delimiter, which is currently set to \""
                          + proj.getProperty(proj.SOURCE_FILE_DELIMITER).getDelimiter()
                          + "\", is incorrect? This can be corrected in the file "
                          + proj.getPropertyFilename()
                          + ". In the meantime, the most stable delimiter will be determined for you...");
          log.reportError("   OR perhaps the ID_HEADER property is incorrect; the text '"
                          + proj.getProperty(proj.ID_HEADER)
                          + "' should be present in the header line.");
          log.reportError("   OR perhaps the ARRAY_TYPE property is incorrect;  options are "
                          + ArrayUtils.toStr(Project.ARRAY.class, ","));

          proj.message("Failed to auto-detect the delimiter used in the Final Reports file; exiting");
          return 0;
        }

        log.reportError("      - determined delimiter to be "
                        + SourceFileParser.DELIMITER_DESCRIPTIONS[ext.indexOfStr(delimiter,
                                                                                 SourceFileParser.DELIMITERS)]);

        // Tries again to determine the header fields and column names
        reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[0]);
        do {
          line = reader.readLine().trim().split(delimiter, -1);
        } while (reader.ready() && (ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line,
                                                     false, true, false)[0] == -1
                                    || (!idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)
                                        && ext.indexOfStr(idHeader, line) == -1)));
      }

      PSF.checkInterrupted();
      log.report(ext.getTime() + "]\tSearching for data fields...");
      // check immediately to make sure these fields are valid
      indices = ext.indexFactors(Sample.DATA_FIELDS, line, false, true, false); // dataIndices
      // if (indices[3] == -1 || indices[4] == -1) {
      // log.reportError("Error - at the very least the files need to contain
      // "+Array.toStr(Sample.DATA_FIELDS[3], "/")+" and "+Array.toStr(Sample.DATA_FIELDS[4], "/"));
      // log.reportError(" - failed to see that in "+files[0]);
      // log.reportError(Array.toStr(line));
      // return 0;
      // }

      // TODO check different fields depending upon Affy/Illumina flag
      log.report(ext.getTime() + "]\tSearching for other fields...");
      indices = ext.indexFactors(Sample.GENOTYPE_FIELDS, line, false, true, false); // genotypeIndices
      // if (indices[0] == -1 || indices[1] == -1) {
      // log.reportError("Error - the files need to contain "+Array.toStr(Sample.GENOTYPE_FIELDS[0],
      // "/")+" and "+Array.toStr(Sample.GENOTYPE_FIELDS[1], "/"));
      // return 0;
      // }
      if (indices[2] == -1 || indices[3] == -1) {
        abLookupRequired = true;
      }

      if (!idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)) {
        ext.indexFactors(new String[] {idHeader}, line, false); // sampIndex
      }
      snpIndex = ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true, true)[0];

      idHeader = proj.ID_HEADER.getValue();
      reader.mark(1000);

      if (idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)) {
        sampleName = files[0].substring(0,
                                        files[0].indexOf(proj.getProperty(proj.SOURCE_FILENAME_EXTENSION)));
        sampIndex = -7;
        parseAtAt = false;
      } else {
        sampIndex = ext.indexFactors(new String[] {idHeader}, line, false)[0];
        line = reader.readLine().split(delimiter);
        parseAtAt = proj.getProperty(proj.PARSE_AT_AT_SYMBOL);
        if (parseAtAt && line[sampIndex].indexOf("@") == -1) {
          log.reportError("Error - " + idHeader + " '" + line[sampIndex]
                          + "' did not contain an @ sample; if your ID's do not naturally contain at symbols, then set "
                          + proj.PARSE_AT_AT_SYMBOL + " in the project properties file to FALSE");
          parseAtAt = false;
        }
        sampleName = parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@"))
                               : line[sampIndex];
      }

      PSF.checkInterrupted();
      code = checkForExistingFiles(proj, proj.SAMPLE_DIRECTORY.getValue(true, true),
                                   Sample.SAMPLE_FILE_EXTENSION);
      if (code == JOptionPane.NO_OPTION) {
        return code;
      }

      PSF.checkInterrupted();
      // log.report(ext.getTime() + "]\tCleaning up before continuing...");

      // this deletes any files that start with "outliers" and end with ".ser" that aren't
      // "outliers.ser"... why is this here?
      TransposeData.deleteOlderRafs(proj.SAMPLE_DIRECTORY.getValue(true, true),
                                    new String[] {"outliers"}, new String[] {".ser"}, true,
                                    new String[] {"outliers.ser"});

      reader.reset();

      int lines = Files.countLines(proj.SOURCE_DIRECTORY.getValue(false, true) + files[0], count);
      markerNameHash = new Hashtable<>(lines + (lines / 3) + 100); // calc to never
                                                                   // re-balance
      log.report(ext.getTime() + "]\tFound " + lines + " rows of data in the first file");
      while (reader.ready()) {
        PSF.checkInterrupted();
        line = reader.readLine().trim().split(delimiter);
        trav = line[snpIndex];
        if (trav.equals("") || trav.equals("0")) {
          trav = "Blank" + count;
          count++;
        }
        if (markerNameHash.containsKey(trav)) {

          if (!proj.LONG_FORMAT.getValue()) {
            log.reportError("The same marker ('" + trav + "') was seen twice...");
          }

          if (idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)) {
            if (!proj.LONG_FORMAT.getValue()) {
              log.report("... this could mean that the file contains multiple samples. This should not happen when "
                         + proj.ID_HEADER.getName() + " is set to "
                         + SourceFileParser.FILENAME_AS_ID_OPTION);
            } else {
              log.reportError("Error - " + proj.ID_HEADER.getName() + " was set to "
                              + SourceFileParser.FILENAME_AS_ID_OPTION
                              + ", which is invalid for Long Format files");
            }
            return 0;
          } else if (!(parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@"))
                                 : line[sampIndex]).equals(sampleName)) {
            PSF.checkInterrupted();
            if (!proj.LONG_FORMAT.getValue()) {
              log.reportError("... and the sample name changed from " + sampleName + " to "
                              + (parseAtAt ? line[sampIndex].substring(0,
                                                                       line[sampIndex].indexOf("@"))
                                           : line[sampIndex])
                              + ", so this must be a Long Format file. The property file will be updated to reflect this, and an attempt will be made to launch the Long Format file processor now.");
              proj.LONG_FORMAT.setValue(Boolean.TRUE);
              proj.saveProperties();
            }
            reader.close();
            // we should have all markers now...
            markerNames = ArrayUtils.mapToValueSortedArray(markerNameHash);
            keys = Markers.orderMarkers(markerNames, proj);
            PSF.checkInterrupted();
            if (keys == null) {
              return 0;// checkForSNP_Map(proj, log);
            } else {
              return createFilesFromLongFormat(proj, files, idHeader, fixes, delimiter,
                                               abLookupRequired, timeBegan, numThreads);
            }
          } else {
            log.reportError("Error - detected the same marker ('" + trav
                            + "') twice, but the sample name did not change. Could not determine the file format to parse");
            return 0;
          }
        } else {
          PSF.checkInterrupted();
          markerNameHash.put(trav, markerNameHash.size());
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + proj.SOURCE_DIRECTORY.getValue(false, true) + files[0]
                      + "\" not found in current directory");
      return 0;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + proj.SOURCE_DIRECTORY.getValue(false, true)
                      + files[0] + "\"");
      log.reportException(ioe);
      return 0;
    }

    PSF.checkInterrupted();
    markerNames = ArrayUtils.mapToValueSortedArray(markerNameHash);
    keys = Markers.orderMarkers(markerNames, proj);
    if (keys == null) {
      return 0;// checkForSNP_Map(proj, log);
    }

    PSF.checkInterrupted();
    keysKeys = Sort.getSortedIndices(keys); // very important
    fingerprint = proj.getMarkerSet().getFingerprint();
    log.report("There are " + markerNames.length + " markers being processed (fingerprint: "
               + fingerprint + ")");
    if (proj.XY_SCALE_FACTOR.getValue() != 1) {
      log.report("The XY_SCALE_FACTOR was set to " + proj.XY_SCALE_FACTOR.getValue()
                 + ", scaling X and Y values...");
    }
    lookup = SourceFileParser.getABLookup(abLookupRequired, markerNames, proj);

    PSF.checkInterrupted();

    log.report(buildColumnAssignmentsLogOutput(proj.getSourceFileHeaders(false).get(files[0])));

    if (affyProcess != null) {
      affyProcess.combineAll(numThreads);
    }
    fileCabinet = new Vector<>();
    for (int i = 0; i < numThreads; i++) {
      fileCabinet.add(new Vector<String>());
    }
    for (int i = 0; i < files.length; i++) {
      fileCabinet.elementAt(i % numThreads).add(files[i]);
    }

    ExecutorService executor = Executors.newFixedThreadPool(numThreads);
    List<Future<?>> futures = Lists.newArrayList();
    for (int i = 0; i < numThreads; i++) {
      String[] threadFiles = fileCabinet.elementAt(i)
                                        .toArray(new String[fileCabinet.elementAt(i).size()]);
      if (threadFiles.length == 0) continue;
      SourceFileParser parser = new SourceFileParser(proj, threadFiles, markerNames, keysKeys,
                                                     lookup, delimiter, fingerprint, fixes,
                                                     timeBegan, i);
      futures.add(executor.submit(parser));
    }
    executor.shutdown();

    if (files.length > 100) {
      try {
        // 5 minutes
        Thread.sleep(1000 * 60 * 5);
      } catch (InterruptedException e1) {}
    }
    boolean failed = false;
    watch: while (!executor.isTerminated()) {
      for (Future<?> f : futures) {
        try {
          f.get(1, TimeUnit.SECONDS);
        } catch (ExecutionException e) {
          Throwable cause = e.getCause();
          if (cause instanceof OutOfMemoryError) {
            log.reportError("Source file parsing ran out of memory - please restart Genvisis with more memory or fewer threads.");
          } else {
            log.reportError("Source file parsing encountered an error (" + cause.getMessage()
                            + ") and could not complete.  Please fix the error and restart Genvisis, or contact help@genvisis.org for assistance.");
          }
          failed = true;
          executor.shutdownNow();
          break watch;
        } catch (InterruptedException e) {
          // ignore
        } catch (TimeoutException e) {
          // ignore
        }
        int mins = files.length > 100 ? 5 : 2;
        try {
          // 5 minutes
          Thread.sleep(1000 * 60 * mins);
        } catch (InterruptedException e) {}
      }
    }
    if (failed) {
      String[] remove = Files.list(proj.SAMPLE_DIRECTORY.getValue(), "",
                                   Sample.SAMPLE_FILE_EXTENSION, false, true);
      for (String f : remove) {
        new File(f).delete();
      }
      remove = Files.list(proj.SAMPLE_DIRECTORY.getValue(), "", ".ser", false, true);
      for (String f : remove) {
        new File(f).delete();
      }
      return 0;
    }

    log.report(ext.getTime() + "]\tWriting sample list...");

    PSF.checkInterrupted();
    SampleList.generateSampleList(proj)
              .writeToTextFile(proj.PROJECT_DIRECTORY.getValue() + "ListOfSamples.txt");
    if (!proj.LONG_FORMAT.getValue()) {
      if (files.length != proj.getSamples().length) {
        proj.getLog()
            .reportError("The number of parsed samples (" + proj.getSamples().length
                         + ") does not equal the number of source files detected (" + files.length
                         + ")");
        proj.getLog()
            .reportError("Please verify source file integrity, or send us a test file to troubleshoot");
      }
    }
    allOutliers = new Hashtable<>();

    v = new Vector<>();

    PSF.checkInterrupted();
    for (int i = 0; i < numThreads; i++) {
      if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser").exists()) {
        allOutliers.putAll((Hashtable<String, Float>) SerializedFiles.readSerial(proj.SAMPLE_DIRECTORY.getValue(true,
                                                                                                                true)
                                                                                 + "outliers" + i
                                                                                 + ".ser"));
        new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers" + i + ".ser").delete();
      }
      if (new File(proj.PROJECT_DIRECTORY.getValue() + IDS_CHANGED_FILEROOT + i
                   + ".txt").exists()) {
        v.add(proj.PROJECT_DIRECTORY.getValue() + IDS_CHANGED_FILEROOT + i + ".txt");
      }
    }
    if (v.size() > 0) {
      Files.cat(ArrayUtils.toStringArray(v),
                proj.PROJECT_DIRECTORY.getValue() + IDS_CHANGED_FILEROOT + ".txt",
                ArrayUtils.intArray(v.size(), 1), log);
      for (int i = 0; i < v.size(); i++) {
        new File(v.elementAt(i)).delete();
      }
    }

    PSF.checkInterrupted();
    // changed 6-26-15 to always write an outliers.ser, even if no outliers exist
    SerializedFiles.writeSerial(allOutliers,
                                proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
    log.report(ext.getTime() + "]\tWrote " + allOutliers.size() + " outliers to file.");
    log.report(ext.getTime() + "]\tSource File Parsing Complete.");
    if (abLookupRequired && !Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false))) {
      return 6;
    } else {
      return 1;
    }

  }

  private static int checkForExistingFiles(Project proj, String dir, final String ext) {
    boolean foundFiles;
    int response = JOptionPane.YES_OPTION;
    String[] overwriteOptions;

    foundFiles = (new File(dir)).listFiles(new FilenameFilter() {

      @Override
      public boolean accept(File dir, String name) {
        return name.endsWith(ext);
      }
    }).length > 0;

    if (foundFiles) {

      overwriteOptions = new String[] {"Delete All", "Cancel Parser"};

      if (GraphicsEnvironment.isHeadless()) {
        proj.getLog()
            .reportTimeWarning("Parsed data files were found in " + dir
                               + ". Data parsing will not continue. If reparsing is desired, remove all "
                               + ext + " files from " + dir);
        response = JOptionPane.NO_OPTION;
      } else {
        response = JOptionPane.showOptionDialog(null,
                                                "Parsed data files were found in " + dir + ".\n"
                                                      + "This happens if the parser was interrupted previously, or restarted unintentionally.\n"
                                                      + "If you would like to reparse this data, select \""
                                                      + overwriteOptions[0] + "\" earlier files.\n"
                                                      + "Otherwise, select \"" + overwriteOptions[1]
                                                      + "\".\n" + "What would you like to do?",
                                                "Parsed files already exist",
                                                JOptionPane.DEFAULT_OPTION,
                                                JOptionPane.QUESTION_MESSAGE, null,
                                                overwriteOptions, overwriteOptions[1]);
      }

      switch (response) {
        case -1:
          response = JOptionPane.NO_OPTION;
          break;
        case 0:
          deleteAllFilesInDirectory(proj, dir, ext);
          response = JOptionPane.YES_OPTION;
          break;
        case 1:
          response = JOptionPane.NO_OPTION;
          break;
        default:
          proj.message("Should have been impossible to obtain this message (" + response + ")");
          response = JOptionPane.NO_OPTION;
      }
    }

    return response;
  }

  private static int checkForSNP_Map(Project proj, Logger log) {
    String filename;
    filename = proj.getLocationOfSNP_Map(true);
    if (filename != null) {
      log.reportError("\nSuch a file was found: " + filename);
      log.reportError("\nIn order to process it use the command:");
      log.reportError("   java -jar [package_name].jar cnv.manage.Markers proj="
                      + proj.getPropertyFilename() + " snps=" + filename);
      return 7;
    }
    return 0;
  }

  public static void deleteAllFilesInDirectory(Project proj, String dir, String ext) {
    String[] filesToDelete;
    filesToDelete = Files.list(dir, ext);
    for (String element : filesToDelete) {
      new File(dir + element).delete();
    }
    new File(dir + "outliers.ser").delete();
  }

  public static int createFilesFromLongFormat(Project proj, String[] files, String idHeader,
                                              Hashtable<String, String> fixes, String delimiter,
                                              boolean abLookupRequired, long timeBegan,
                                              int numThreads) {
    PrintWriter writer;
    String[] markerNames, list;
    long fingerprint;
    MarkerSetInfo markerSet;
    Hashtable<String, Integer> markerIndexMap;
    char[][] abLookup;
    Hashtable<String, Float> allOutliers;
    Hashtable<String, String> renamedIDsHash;
    Map<String, SourceFileHeaderData> headers;
    Logger log;

    log = proj.getLog();
    log.report("Parsing files using the Long Format algorithm");

    // creates and serializes the markers.bim, the returned keys are not used here as the
    // markerIndices fill that purpose
    markerSet = proj.getMarkerSet();
    markerNames = markerSet.getMarkerNames();
    fingerprint = proj.getMarkerSet().getFingerprint();

    PSF.checkInterrupted();
    abLookup = SourceFileParser.getABLookup(abLookupRequired, markerNames, proj);

    markerIndexMap = new Hashtable<>();
    for (int i = 0; i < markerNames.length; i++) {
      markerIndexMap.put(markerNames[i], new Integer(i));
    }

    log.report("There were " + markerNames.length
               + " markers present in the project that will be processed from the source files (fingerprint: "
               + fingerprint + ")");

    // int[] dataIndices, genotypeIndices;

    CountHash dupHash = new CountHash();
    CountHash countHash = new CountHash();// why are we using this?
    if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + "FYI_IDS_WERE_CHANGED.txt")) {
      Files.backup("FYI_IDS_WERE_CHANGED.txt", proj.PROJECT_DIRECTORY.getValue(),
                   proj.PROJECT_DIRECTORY.getValue(), true);
    }

    allOutliers = new Hashtable<>();
    renamedIDsHash = new Hashtable<>();
    headers = proj.getSourceFileHeaders(true);// setting to true fixed an issue parsing NGRC data
    log.report(buildColumnAssignmentsLogOutput(headers.get(files[0])));
    int count = 0;
    try {
      LongFileFormatProducer producer = new LongFileFormatProducer(proj, files, idHeader, fixes,
                                                                   delimiter, markerNames,
                                                                   fingerprint, markerIndexMap,
                                                                   abLookup, renamedIDsHash,
                                                                   headers, log);
      try (WorkerTrain<LongFormatParseResult> train = new WorkerTrain<>(producer, numThreads, 10,
                                                                        log)) {
        while (train.hasNext()) {
          LongFormatParseResult result = train.next();
          if (result.count == 0) {
            log.reportError("Encountered an error processing " + result.fileParsed + " , halting");
            train.close();
            return 0;
          }
          count += result.count;
          allOutliers.putAll(result.outliers);
          for (String key : result.countHash.getHash().keySet()) {
            countHash.add(key);
          }
        }
      }
      // for (int i = 0; i < files.length; i++) {
      // if (Thread.currentThread().isInterrupted()) {
      // throw new RuntimeException(new InterruptedException());
      // }
      //
      // parseLongFormatFile(proj, files[i], idHeader, fixes, delimiter, markerNames, fingerprint,
      // markerIndexMap, abLookup, renamedIDsHash, headers, log, countHash, i, files.length);
      // }

      PSF.checkInterrupted();

      if (allOutliers.size() > 0) {
        if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser").exists()) {// why
                                                                                             // do
                                                                                             // we
                                                                                             // care?
          log.reportError("Error - the following file already exists: "
                          + proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
          return 0;
        } else {
          SerializedFiles.writeSerial(allOutliers,
                                      proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
        }
      } else {
        SerializedFiles.writeSerial(allOutliers,
                                    proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
      }

      log.report(ext.getTime() + "\tfinished");
    } catch (Exception e) {
      log.reportException(e);
    }

    PSF.checkInterrupted();

    log.report(ext.getTime() + "\t" + "Parsed " + count + " sample(s)");
    SampleList.generateSampleList(proj)
              .writeToTextFile(proj.PROJECT_DIRECTORY.getValue() + "ListOfSamples.txt");

    PSF.checkInterrupted();
    try {
      writer = Files.openAppropriateWriter(proj.PROJECT_DIRECTORY.getValue() + "ListOfMarkers.txt");
      writer.println("Marker\tExpected\tTimesSeen\tTimesDuplicated");
      for (String markerName : markerNames) {
        writer.println(markerName + "\t1\t" + countHash.getCount(markerName) + "\t"
                       + dupHash.getCount(markerName));
        countHash.remove(markerName, false);
      }
      list = countHash.getValues();
      for (int j = 0; j < list.length; j++) {
        writer.println(list[j] + "\t0\t" + countHash.getCount(markerNames[j]) + "\t.");
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + "ListOfMarkers.txt");
      log.reportException(e);
    }

    PSF.checkInterrupted();

    new File(proj.SAMPLE_DIRECTORY.getValue(true, true)
             + SourceFileParser.OVERWRITE_OPTION_FILE).delete();
    new File(proj.SAMPLE_DIRECTORY.getValue(true, true)
             + SourceFileParser.HOLD_OPTION_FILE).delete();
    new File(proj.SAMPLE_DIRECTORY.getValue(true, true)
             + SourceFileParser.CANCEL_OPTION_FILE).delete();

    if (abLookupRequired && !Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false))) {
      return 6;
    } else {
      return 1;
    }
  }

  /**
   * Store a few things from the parsing of a long format file
   */
  private static class LongFormatParseResult {

    private final int count;
    private final Hashtable<String, Float> outliers;
    private final String fileParsed;
    private final CountHash countHash;

    public LongFormatParseResult(String fileParsed, int count, Hashtable<String, Float> outliers,
                                 CountHash countHash) {
      super();
      this.count = count;
      this.outliers = outliers;
      this.fileParsed = fileParsed;
      this.countHash = countHash;
    }
  }

  private static class LongFileFormatProducer extends AbstractProducer<LongFormatParseResult> {

    private final Project proj;
    private final String[] files;
    private final String idHeader;
    private final Hashtable<String, String> fixes;
    private final String delimiter;
    private final String[] markerNames;
    private final long fingerprint;
    private final Hashtable<String, Integer> markerIndexMap;
    private final char[][] abLookup;
    private final Hashtable<String, String> renamedIDsHash;
    private final Map<String, SourceFileHeaderData> headers;
    private final Logger log;
    private int fileIndex;

    public LongFileFormatProducer(Project proj, String[] files, String idHeader,
                                  Hashtable<String, String> fixes, String delimiter,
                                  String[] markerNames, long fingerprint,
                                  Hashtable<String, Integer> markerIndexMap, char[][] abLookup,
                                  Hashtable<String, String> renamedIDsHash,
                                  Map<String, SourceFileHeaderData> headers, Logger log) {
      super();
      this.proj = proj;
      this.files = files;
      this.idHeader = idHeader;
      this.fixes = fixes;
      this.delimiter = delimiter;
      this.markerNames = markerNames;
      this.fingerprint = fingerprint;
      this.markerIndexMap = markerIndexMap;
      this.abLookup = abLookup;
      this.renamedIDsHash = renamedIDsHash;
      this.headers = headers;
      this.log = log;
      fileIndex = 0;
    }

    @Override
    public boolean hasNext() {
      // TODO Auto-generated method stub
      return fileIndex < files.length;
    }

    @Override
    public Callable<LongFormatParseResult> next() {
      String file = files[fileIndex];
      fileIndex++;
      return new LongFileFormatWorker(proj, file, idHeader, fixes, delimiter, markerNames,
                                      fingerprint, markerIndexMap, abLookup, renamedIDsHash,
                                      headers, log, fileIndex, files.length);
    }
  }

  private static class LongFileFormatWorker implements Callable<LongFormatParseResult> {

    private final Project proj;
    private final String file;
    private final String idHeader;
    private final Hashtable<String, String> fixes;
    private final String delimiter;
    private final String[] markerNames;
    private final long fingerprint;
    private final Hashtable<String, Integer> markerIndexMap;
    private final char[][] abLookup;
    private final Hashtable<String, String> renamedIDsHash;
    private final Map<String, SourceFileHeaderData> headers;
    private final Logger log;
    private final int fileIndex;
    private final int numFiles;

    public LongFileFormatWorker(Project proj, String file, String idHeader,
                                Hashtable<String, String> fixes, String delimiter,
                                String[] markerNames, long fingerprint,
                                Hashtable<String, Integer> markerIndexMap, char[][] abLookup,
                                Hashtable<String, String> renamedIDsHash,
                                Map<String, SourceFileHeaderData> headers, Logger log,
                                int fileIndex, int numFiles) {
      super();
      this.proj = proj;
      this.file = file;
      this.idHeader = idHeader;
      this.fixes = fixes;
      this.delimiter = delimiter;
      this.markerNames = markerNames;
      this.fingerprint = fingerprint;
      this.markerIndexMap = markerIndexMap;
      this.abLookup = abLookup;
      this.renamedIDsHash = renamedIDsHash;
      this.headers = headers;
      this.log = log;
      this.fileIndex = fileIndex;
      this.numFiles = numFiles;
    }

    @Override
    public LongFormatParseResult call() throws Exception {
      return parseLongFormatFile(proj, file, idHeader, fixes, delimiter, markerNames, fingerprint,
                                 markerIndexMap, abLookup, renamedIDsHash, headers, log, fileIndex,
                                 numFiles);
    }

  }

  /**
   * @param proj
   * @param file
   * @param idHeader
   * @param fixes
   * @param delimiter
   * @param markerNames
   * @param fingerprint
   * @param markerIndexMap
   * @param abLookup
   * @param renamedIDsHash
   * @param headers
   * @param log
   * @param countHash
   * @param fileIndex
   * @param numFiles
   * @return
   * @throws Elision
   */
  private static LongFormatParseResult parseLongFormatFile(Project proj, String file,
                                                           String idHeader,
                                                           Hashtable<String, String> fixes,
                                                           String delimiter, String[] markerNames,
                                                           long fingerprint,
                                                           Hashtable<String, Integer> markerIndexMap,
                                                           char[][] abLookup,
                                                           Hashtable<String, String> renamedIDsHash,
                                                           Map<String, SourceFileHeaderData> headers,
                                                           Logger log, int fileIndex,
                                                           int numFiles) throws Elision {
    // might have some variable scope issues, but we will see.

    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String filename;
    int snpIndex;
    int sampIndex;
    String sampleIdent;
    String tempSampleIdent;
    boolean parseAtAt;
    Sample samp;
    String sampleName;
    float[][] data;
    byte[][] genotypes;
    boolean ignoreAB;
    boolean eofFound;
    int numCols;
    String tempLine;
    boolean splitAB = false;
    long timeBegan = System.currentTimeMillis();
    int sampleCount = 0;
    int lineCount = 0;
    Hashtable<String, Float> fileOutliers = new Hashtable<>();
    CountHash countHash = new CountHash();
    LongFormatParseResult result = new LongFormatParseResult(file, 0, null, countHash);
    try {
      int markerCount = 0;

      log.report(ext.getTime() + "]\t" + fileIndex + " of " + numFiles + " (" + file + ")");

      // DATA_0 {"GC Score", "GCscore", "confidence","Confidence"}
      // DATA_1 {"X Raw"}
      // DATA_2 {"Y Raw"}
      // DATA_3 {"X", "Xvalue", "Log Ratio", "intensity_1","Signal A"}
      // DATA_4 {"Y", "Yvalue", "Strength", "intensity_2","Signal B"}
      // DATA_5 {"Theta"}
      // DATA_6 {"R"}
      // DATA_7 {"B Allele Freq"},
      // DATA_8 {"Log R Ratio"}
      //
      // GENO_0 {"Allele1 - Forward", "Allele1", "genotype1", "Allele1 - Top","Forward Strand Base
      // Calls","Forced Call","Forced Call Codes"},
      // GENO_1 {"Allele2 - Forward", "Forward Strand Base Calls", "genotype2", "Allele2 -
      // Top","Allele B","Forced Call","Forced Call Codes"},
      // GENO_2 {"Allele1 - AB","Call Codes","Call"},
      // GENO_3 {"Allele2 - AB","Call Codes","Call"}};

      SourceFileHeaderData headerData = null;
      if (headers.containsKey(file)) {
        headerData = headers.get(file);
      }
      reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + file);

      if (headerData == null) {
        // TODO error, missing header for source files
        log.report("WARNING - No parsed header object found for file " + file
                   + ".  Parsing now.  If columns are non-standard, or an error occurs, please re-create the project from scratch and validate source files.");
        headerData = SourceFileHeaderData.parseHeader(file, log);
        headers.put(file, headerData);
      }

      PSF.checkInterrupted();

      for (int k = 0; k < headerData.getColumnHeaderLineIndex() + 1; k++) {
        // iterate
        reader.readLine();
      }

      PSF.checkInterrupted();
      if ((headerData.getColGenoAB1() == -1 || headerData.getColGenoAB2() == -1)
          && abLookup == null) {
        ignoreAB = true;
      } else {
        ignoreAB = false;
      }

      sampIndex = headerData.getColSampleIdent();
      snpIndex = headerData.getColSnpIdent();

      sampleName = "just starting";
      data = null;
      genotypes = null;
      parseAtAt = proj.PARSE_AT_AT_SYMBOL.getValue();
      numCols = -1;
      eofFound = false;
      while (!eofFound) {

        PSF.checkInterrupted();
        if ((tempLine = reader.readLine()) != null) {

          line = tempLine.split(delimiter);
          if (numCols == -1) {
            numCols = line.length;
          } else if (line.length != numCols) {
            log.reportError("Error - mismatched number of columns at marker index " + markerCount);
          }
          if (parseAtAt && line[sampIndex].indexOf("@") == -1) {
            log.reportError("Error - " + idHeader + " '" + line[sampIndex]
                            + "' did not contain an @ sample");
            parseAtAt = false;
          }
          tempSampleIdent = parseAtAt ? line[sampIndex].substring(0, line[sampIndex].indexOf("@"))
                                      : line[sampIndex];
          sampleIdent = ext.replaceWithLinuxSafeCharacters(tempSampleIdent, true);

          if (!sampleIdent.equals(tempSampleIdent)
              && !renamedIDsHash.containsKey(tempSampleIdent)) {
            try {
              writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()
                                                      + "FYI_IDS_WERE_CHANGED.txt", true));
              if (renamedIDsHash.size() == 0) {
                writer.println("The following IDs were changed so that spaces are removed and so that they could be used as valid filenames:");
              }
              writer.println(tempSampleIdent + "\t" + sampleIdent);
              writer.close();
            } catch (Exception e) {
              log.reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue()
                              + "FYI_IDS_WERE_CHANGED.txt");
              log.reportException(e);
            }
            renamedIDsHash.put(tempSampleIdent, sampleIdent);
          }
        } else {
          eofFound = true;
          sampleIdent = null;
          line = null;
        }

        if (eofFound || !sampleIdent.equals(sampleName)) {
          if (!sampleName.equals("just starting")) {
            if (markerCount != markerNames.length) {
              log.reportError("Error - expecting " + markerNames.length + " markers and only found "
                              + markerCount + " for sample " + sampleName
                              + "; this usually happens when the input file is truncated");
              return result;
            }

            if (fixes.containsKey(sampleName)) {
              sampleName = fixes.get(sampleName);
            }

            filename = SourceFileParser.determineFilename(proj.SAMPLE_DIRECTORY.getValue(true,
                                                                                         true),
                                                          sampleName, timeBegan, log);
            if (filename == null) {
              return result;
            }

            PSF.checkInterrupted();

            samp = new Sample(sampleName, fingerprint, data, genotypes,
                              proj.getArrayType().getCanXYBeNegative());
            samp.saveToRandomAccessFile(filename, fileOutliers, sampleName);

            sampleCount++;
            markerCount = 0;
          }
          if (!eofFound) {
            String sampleFilename = proj.SAMPLE_DIRECTORY.getValue(true, true)
                                    + (fixes.containsKey(sampleIdent) ? fixes.get(sampleIdent)
                                                                      : sampleIdent)
                                    + Sample.SAMPLE_FILE_EXTENSION;
            if (new File(sampleFilename).exists()) {
              log.reportError("Warning - marker data must be out of order, becaue we're seeing "
                              + sampleIdent
                              + (fixes.containsKey(sampleIdent) ? "-->" + fixes.get(sampleIdent)
                                                                : "")
                              + " again at line " + lineCount);
              samp = Sample.loadFromRandomAccessFile(sampleFilename);
              data = samp.getAllData();
              genotypes = samp.getAllGenotypes();
            } else {
              data = new float[][] {headerData.getColGC() == -1 ? null
                                                                : new float[markerNames.length],
                                    headerData.getColXRaw() == -1 ? null
                                                                  : new float[markerNames.length],
                                    headerData.getColYRaw() == -1 ? null
                                                                  : new float[markerNames.length],
                                    headerData.getColX() == -1 ? null
                                                               : new float[markerNames.length],
                                    headerData.getColY() == -1 ? null
                                                               : new float[markerNames.length],
                                    headerData.getColTheta() == -1 ? null
                                                                   : new float[markerNames.length],
                                    headerData.getColR() == -1 ? null
                                                               : new float[markerNames.length],
                                    headerData.getColBAF() == -1 ? null
                                                                 : new float[markerNames.length],
                                    headerData.getColLRR() == -1 ? null
                                                                 : new float[markerNames.length],};
              genotypes = new byte[2][];
              genotypes[0] = ArrayUtils.byteArray(markerNames.length, (byte) 0);
              if (!ignoreAB) {
                genotypes[1] = ArrayUtils.byteArray(markerNames.length, (byte) -1);
              }
            }
            sampleName = sampleIdent;
          }
        }
        if (!eofFound) {
          countHash.add(line[snpIndex]);
          if (markerIndexMap.containsKey(line[snpIndex])) {
            final int markerKey = markerIndexMap.get(line[snpIndex]);

            int[] indCols = {headerData.getColGC(), headerData.getColXRaw(),
                             headerData.getColYRaw(), headerData.getColX(), headerData.getColY(),
                             headerData.getColTheta(), headerData.getColR(), headerData.getColBAF(),
                             headerData.getColLRR(),};

            int ind = 0;
            try {
              for (int i1 = 0; i1 < indCols.length; i1++) {
                ind = i1;
                if (indCols[ind] != -1) {
                  data[ind][markerKey] = ext.isMissingValue(line[indCols[ind]]) ? Float.NaN
                                                                                : Float.parseFloat(line[indCols[ind]]);
                }
              }
            } catch (NumberFormatException nfe) {
              log.reportError("Error - failed at line " + lineCount + " to parse '"
                              + line[indCols[ind]] + "' into a valid "
                              + ArrayUtils.toStr(Sample.DATA_FIELDS[ind], "/"));
              return result;
            } catch (Exception e) {
              log.reportError("Some other exception");
              log.reportException(e);
              return result;
            }
            if (headerData.getColGeno1() != -1 && headerData.getColGeno2() != -1) {
              if (line[headerData.getColGeno1()].length() > 1) {
                if (fileIndex == 0 && !splitAB) {
                  log.reportTimeInfo("Detected genotype calls for each allele are in the same column (Such as in Affymetrix .chp format) ");
                  splitAB = true;
                }
                if (splitAB) {
                  String tmp0 = line[headerData.getColGeno1()].substring(0, 1);
                  String tmp1 = line[headerData.getColGeno1()].substring(1, 2);
                  line[headerData.getColGeno1()] = tmp0;
                  line[headerData.getColGeno2()] = tmp1;
                  if (!ignoreAB) {
                    try {
                      String tmp2 = line[headerData.getColGenoAB1()].substring(0, 1);
                      String tmp3 = line[headerData.getColGenoAB1()].substring(1, 2);
                      line[headerData.getColGenoAB1()] = tmp2;
                      line[headerData.getColGenoAB2()] = tmp3;
                    } catch (Exception e) {
                      log.reportError("Could not parse genotypes on line "
                                      + ArrayUtils.toStr(line));
                      log.reportException(e);
                      return result;
                    }
                  }
                } else {
                  log.reportError("Inconsistant genotype call lengths");
                }
              } else if (splitAB) {
                log.reportError("Detected previously that genotype calls should be split, but the calls on line "
                                + ArrayUtils.toStr(line) + " did not");
                return result;
              }
              genotypes[0][markerKey] = (byte) ext.indexOfStr(line[headerData.getColGeno1()]
                                                              + line[headerData.getColGeno2()],
                                                              Sample.ALLELE_PAIRS);
              if (genotypes[0][markerKey] == -1) {
                if (proj.getArrayType() == ARRAY.ILLUMINA) {// Affy matrix format does not use
                                                            // ALLELE_PAIRS

                  if (ext.indexOfStr(line[headerData.getColGeno1()]
                                     + line[headerData.getColGeno2()], Sample.ALT_NULLS) == -1) {
                    log.reportError("Error - failed to lookup genotype ("
                                    + line[headerData.getColGeno1()]
                                    + line[headerData.getColGeno2()] + ") for marker "
                                    + markerNames[markerKey] + " of sample " + sampleName
                                    + "; setting to missing");
                  }
                  genotypes[0][markerKey] = 0;
                } else {
                  genotypes[0][markerKey] = (byte) ext.indexOfStr(line[headerData.getColGenoAB1()]
                                                                  + line[headerData.getColGenoAB2()],
                                                                  Sample.AB_PAIRS);
                  if (genotypes[0][markerKey] == -1) {
                    log.reportError("Error - failed to lookup genotype ("
                                    + line[headerData.getColGeno1()]
                                    + line[headerData.getColGeno2()] + ") for marker "
                                    + markerNames[markerKey] + " of sample " + sampleName
                                    + "; setting to missing");
                    genotypes[0][markerKey] = 0;
                  }
                }
              }

              if (ignoreAB) {
                // do nothing, will need to use these files to determine AB lookup table
              } else if (abLookup == null) {
                genotypes[1][markerKey] = (byte) ext.indexOfStr(line[headerData.getColGenoAB1()]
                                                                + line[headerData.getColGenoAB2()],
                                                                Sample.AB_PAIRS);
              } else {
                if (genotypes[0][markerKey] == 0) {
                  genotypes[1][markerKey] = -1;
                } else {
                  genotypes[1][markerKey] = 0;
                  if (line[headerData.getColGeno1()].charAt(0) == abLookup[markerKey][1]) {
                    genotypes[1][markerKey]++;
                  } else if (line[headerData.getColGeno1()].charAt(0) != abLookup[markerKey][0]) {
                    log.reportError("Error - alleles for individual '"
                                    + (sampIndex < 0 ? sampleIdent : line[sampIndex]) + "' ("
                                    + line[headerData.getColGeno1()] + "/"
                                    + line[headerData.getColGeno2()]
                                    + ") do not match up with the defined AB lookup alleles ("
                                    + abLookup[markerKey][0] + "/" + abLookup[markerKey][1]
                                    + ") for marker " + markerNames[markerKey]);
                  }
                  if (line[headerData.getColGeno2()].charAt(0) == abLookup[markerKey][1]) {
                    genotypes[1][markerKey]++;
                  } else if (line[headerData.getColGeno2()].charAt(0) != abLookup[markerKey][0]) {
                    log.reportError("Error - alleles for individual '"
                                    + (sampIndex < 0 ? sampleIdent : line[sampIndex]) + "' ("
                                    + line[headerData.getColGeno1()] + "/"
                                    + line[headerData.getColGeno2()]
                                    + ") do not match up with the defined AB lookup alleles ("
                                    + abLookup[markerKey][0] + "/" + abLookup[markerKey][1]
                                    + ") for marker " + markerNames[markerKey]);
                  }
                }
              }
            }
          }
          markerCount++;
        }
        lineCount++;
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + file + "\" not found in current directory");
      return result;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + file + "\"");
      return result;
    }
    result = new LongFormatParseResult(file, sampleCount, fileOutliers, countHash);
    return result;

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    Project proj;
    String filename = null;
    boolean map = false;
    int numThreads = 1;
    boolean parseAlleleLookupFromFinalReports = false;
    String mapOutput = "filenamesMappedToSamples.txt";

    String usage = "\n" + "cnv.manage.ParseIllumina requires 0-1 arguments\n"
                   + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + "   (2) number of threads to use (i.e. " + PSF.Ext.NUM_THREADS_COMMAND
                   + numThreads + " (default))\n" + " OPTIONAL:\n"
                   + "   (3) map filenames to sample IDs (i.e. -mapFiles ("
                   + (map ? "" : "not the ") + "default))\n"
                   + "   (4) output file for mappings (i.e. out=" + mapOutput + " (default))\n"
                   + " OR:\n"
                   + "   (1) parse Forward/TOP/AB/etc lookup (i.e. --parseAlleleLookup (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        return;
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numThreads = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("-mapFiles")) {
        map = true;
        numArgs--;
      } else if (arg.startsWith("-parseAlleleLookup")) {
        parseAlleleLookupFromFinalReports = true;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      return;
    }

    try {
      proj = new Project(filename);
      if (map) {
        SourceParserUtils.mapFilenamesToSamples(proj, mapOutput);
      } else if (parseAlleleLookupFromFinalReports) {
        SourceParserUtils.parseAlleleLookupFromFinalReports(proj);
      } else {
        createFiles(proj, numThreads);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
