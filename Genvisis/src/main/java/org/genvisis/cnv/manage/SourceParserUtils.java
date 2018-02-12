package org.genvisis.cnv.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Map.Entry;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class SourceParserUtils {

  public static void mapFilenamesToSamples(Project proj, String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    int sampIndex;
    String[] files;
    String idHeader, delimiter;
    boolean longFormat, done;
    String prev;
    Logger log;

    log = proj.getLog();
    if (!Files.exists(proj.SOURCE_DIRECTORY.getValue(false, false))) {
      proj.message("Source directory does not exist; change SOURCE_DIRECTORY= to point to the proper files");
      return;
    }

    delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
    longFormat = proj.LONG_FORMAT.getValue();
    idHeader = proj.getProperty(proj.ID_HEADER);
    files = SourceFileParser.getSourceFiles(proj, log);

    try {
      writer = Files.openAppropriateWriter(proj.PROJECT_DIRECTORY.getValue() + filename);
      for (String file : files) {
        try {
          reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + file);
          do {
            line = reader.readLine().trim().split(delimiter);
          } while (reader.ready() && (line.length < 3 || ext.indexOfStr(idHeader, line) == -1));

          if (!reader.ready()) {
            log.reportError("Error - went through enitre file without finding a line containing the user-defined ID header: "
                            + idHeader);
            return;
          }
          sampIndex = ext.indexFactors(new String[] {idHeader}, line, false)[0];

          // ParseAffymetrix, ParseAffySNP6, and ParseDbgap ::>
          // line = reader.readLine().split(delimiter);
          // writer.println(files[i]+"\t"+line[sampIndex]+"\t"+(line[sampIndex].indexOf("@") >=
          // 0?line[sampIndex].split("@")[0]:line[sampIndex]));
          // reader.close();
          //
          // ParseIllumina ::>
          done = false;
          prev = null;
          while (!done) {
            line = reader.readLine().split(delimiter);
            if (!line[sampIndex].equals(prev)) {
              writer.println(file + "\t" + line[sampIndex]
                             + (line[sampIndex].indexOf("@") >= 0 ? "\t"
                                                                    + line[sampIndex].split("@")[0]
                                                                  : ""));
            }
            if (!longFormat || !reader.ready()) {
              done = true;
            }
            prev = line[sampIndex];
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          log.reportError("Error: file \"" + file + "\" not found in "
                          + proj.SOURCE_DIRECTORY.getValue(false, true));
          writer.close();
          return;
        } catch (IOException ioe) {
          log.reportError("Error reading file \"" + file + "\"");
          writer.close();
          return;
        }
      }
      log.report(ext.getTime());
      writer.close();
    } catch (Exception e) {
      log.reportException(e);
    }
  }

  public static void parseAlleleLookupFromFinalReports(Project proj) {
    BufferedReader reader;
    String[] line;
    Hashtable<String, String[]> hash;
    String[] files;
    int[] indices;
    int snpIndex;
    String delimiter, idHeader;
    String[] alleles;
    int expIndex;
    Logger log;

    log = proj.getLog();

    files = SourceFileParser.getSourceFiles(proj, log);
    if (files.length == 0) {
      return;
    }

    idHeader = proj.getProperty(proj.ID_HEADER);
    delimiter = proj.SOURCE_FILE_DELIMITER.getValue().getDelimiter();
    hash = new Hashtable<String, String[]>();
    for (int i = 0; i < files.length; i++) {
      if (new File("report").exists()) {
        SourceParserUtils.writeToLookupFile(proj, hash, i + 1);
        new File("report").delete();
      }
      try {
        log.report(ext.getTime() + "\t" + (i + 1) + " of " + files.length);
        reader = Files.getAppropriateReader(proj.SOURCE_DIRECTORY.getValue(false, true) + files[i]);
        do {
          line = reader.readLine().trim().split(delimiter, -1);
        } while (reader.ready() && (ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line,
                                                     false, true, false)[0] == -1
                                    || (!idHeader.equals(SourceFileParser.FILENAME_AS_ID_OPTION)
                                        && ext.indexOfStr(idHeader, line) == -1)));

        snpIndex = ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true,
                                    false)[0];
        indices = ext.indexFactors(Sample.ALL_STANDARD_GENOTYPE_FIELDS, line, false, proj.getLog(),
                                   false);

        while (reader.ready()) {
          line = reader.readLine().split(delimiter);
          if (hash.containsKey(line[snpIndex])) {
            alleles = hash.get(line[snpIndex]);
          } else {
            if (i != 0) {
              log.reportError("Error - snp '" + line[snpIndex] + "' first seen in file #" + i + " ("
                              + files[i] + ") and not earlier");
            }
            hash.put(line[snpIndex],
                     alleles = new String[Sample.ALL_STANDARD_GENOTYPE_FIELDS.length]);
          }
          for (int j = 0; j < 2; j++) {
            if (ext.indexOfStr(line[indices[j]], Sample.ALT_NULL) == -1) {
              if (alleles[0] == null) {
                alleles[0] = line[indices[j]];
                expIndex = -1;
              } else if (line[indices[j]].equals(alleles[0])) {
                expIndex = 0;
              } else if (alleles[1] == null) {
                alleles[1] = line[indices[j]];
                expIndex = -2;
              } else if (line[indices[j]].equals(alleles[1])) {
                expIndex = 1;
              } else {
                log.reportError("Error - snp '" + line[snpIndex] + "' has a new allele in file #"
                                + i + " (" + files[i] + "): " + line[indices[j]] + " (previously "
                                + alleles[0] + "/" + alleles[1] + ")");
                expIndex = -9;
              }

              for (int k = 1; k < alleles.length / 2; k++) {
                switch (expIndex) {
                  case -1:
                    if (alleles[k * 2 + 0] != null) {
                      log.reportError("Error - how can this be? -1");
                    }
                    alleles[k * 2 + 0] = line[indices[k * 2 + j]];
                    break;
                  case 0:
                    if (!line[indices[k * 2 + j]].equals(alleles[k * 2 + 0])) {
                      log.reportError("Error - how can this be? 0");
                    }
                    break;
                  case -2:
                    if (alleles[k * 2 + 1] != null) {
                      log.reportError("Error - how can this be? -2");
                    }
                    alleles[k * 2 + 1] = line[indices[k * 2 + j]];
                    break;
                  case 1:
                    if (!line[indices[k * 2 + j]].equals(alleles[k * 2 + 1])) {
                      log.reportError("Error - how can this be? 1");
                    }
                    break;
                  case -9:
                    break;
                }
              }
            }
          }
        }
      } catch (FileNotFoundException fnfe) {
        log.reportError("Error: file \"" + files[i] + "\" not found in current directory");
        return;
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + files[i] + "\"");
        return;
      }
    }
    SourceParserUtils.writeToLookupFile(proj, hash, 0);
  }

  public static void writeToLookupFile(Project proj, Hashtable<String, String[]> hash,
                                       int fileNumber) {
    PrintWriter writer;
    String filename;
    Logger log;

    log = proj.getLog();
    filename = proj.PROJECT_DIRECTORY.getValue() + "alleleLookup"
               + (fileNumber > 0 ? "_atFile" + fileNumber : "") + ".xln";
    log.report("Writing to file {" + filename + "}...", false, true);
    try {
      writer = Files.openAppropriateWriter(filename);
      writer.println("SNP\t" + ArrayUtils.toStr(Sample.ALL_STANDARD_GENOTYPE_FIELDS));
      for (Entry<String, String[]> entry : hash.entrySet()) {
        writer.println(entry.getKey() + "\t" + ArrayUtils.toStr(entry.getValue()));
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + filename);
      log.reportException(e);
    }
    log.report("done.");

  }

}
