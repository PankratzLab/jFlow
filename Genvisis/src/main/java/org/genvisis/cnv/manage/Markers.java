package org.genvisis.cnv.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import org.genvisis.cnv.Resources;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.qc.MarkerBlast;
import org.genvisis.cnv.qc.MarkerBlast.FILE_SEQUENCE_TYPE;
import org.pankratzlab.common.Aliases;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.Sort;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.io.Closeables;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class Markers {

  public static final int MAX_ERRORS_TO_REPORT = 30;
  public static final List<String> MARKER_POSITIONS_MISMATCHES_HEADER = ImmutableList.of("Marker",
                                                                                         "Supplied chr",
                                                                                         "Supplied pos",
                                                                                         "dbSNP chr",
                                                                                         "dbSNP pos");

  public static int[] orderMarkers(String[] markerNames, Project proj) {
    return orderMarkers(markerNames, proj, proj.getLog());
  }

  public static int[] orderMarkers(String[] markerNames, Project proj, Logger log) {
    return orderMarkers(markerNames, proj.MARKER_POSITION_FILENAME.getValue(),
                        proj.MARKERSET_FILENAME.getValue(true, false), log);
  }

  public static int[] orderMarkers(String[] markerNames, String markerDatabase, String output,
                                   Logger log) {
    return orderMarkers(markerNames, markerDatabase, output, log, true);
  }

  @Deprecated
  /**
   * @deprecated use #orderMarkers(String[], Project)
   */
  public static int[] orderMarkers(String[] markerNames, String markerDatabase, String output,
                                   Logger log, boolean allowExtraAsMissing) {
    Hashtable<String, String> snpPositions;
    byte[] chrs;
    int[] positions, keys;
    String[] line;
    Vector<String> v;
    long time;
    int logLevel;
    int[] chrCounts;
    HashSet<String> reportMarkers, databaseMarkers, databaseMarkersRef;

    logLevel = log.getLevel();
    log.setLevel(9);
    time = new Date().getTime();
    log.report(ext.getTime() + "]\tLoading marker data from " + markerDatabase);
    snpPositions = loadFileToHashString(markerDatabase, log);
    if (snpPositions == null) {
      return null;
    }
    if (markerNames == null) {
      for (String element : Aliases.MARKER_NAMES) {
        if (snpPositions.containsKey(element)) {
          snpPositions.remove(element);
        }
      }
      markerNames = HashVec.getKeys(snpPositions, false);
    }

    reportMarkers = new HashSet<>();
    for (String str : markerNames) {
      reportMarkers.add(str);
    }
    databaseMarkersRef = new HashSet<>();
    databaseMarkers = new HashSet<>();
    databaseMarkersRef.addAll(snpPositions.keySet());
    databaseMarkers.addAll(snpPositions.keySet());

    databaseMarkers.removeAll(reportMarkers);
    int posNotFoundInRpt = databaseMarkers.size();

    reportMarkers.removeAll(databaseMarkersRef);
    int rptNotFoundInPos = reportMarkers.size();

    if (posNotFoundInRpt > 0) {
      log.reportError("Error - There "
                      + (posNotFoundInRpt == 1 ? "was one" : "were " + posNotFoundInRpt)
                      + " markers found in the file of marker positions that were not listed in the FinalReport file.");
      Files.writeArray(databaseMarkers.toArray(new String[] {}),
                       ext.parseDirectoryOfFile(markerDatabase) + "markersNotInSourceFile.txt");
    }
    if (rptNotFoundInPos > 0 && !allowExtraAsMissing) {
      log.reportError("Error - There "
                      + (rptNotFoundInPos == 1 ? "was one" : "were " + rptNotFoundInPos)
                      + " markers found in the FinalReport file that were not listed in the file of marker positions; halting parse operation.");
      log.reportError("\nThe best source of complete marker positions is the SNP manifest (e.g., SNP_Map.csv from Illumina's GenomeStudio that should be exported along with the FinalReport files)");
      Files.writeArray(reportMarkers.toArray(new String[] {}),
                       ext.parseDirectoryOfFile(markerDatabase) + "markersNotInPositionsFile.txt");
    }
    databaseMarkersRef = null;
    databaseMarkers = null;
    reportMarkers = null;
    if (rptNotFoundInPos > 0 && !allowExtraAsMissing) {
      return null;
    }

    v = new Vector<>();
    log.report(ext.getTime() + "\tSorting markers by chromosome and position");
    chrs = new byte[markerNames.length];
    chrCounts = new int[Positions.CHR_CODES.length];
    positions = new int[markerNames.length];
    for (int i = 0; i < markerNames.length; i++) {
      if (snpPositions.containsKey(markerNames[i])) {
        line = snpPositions.get(markerNames[i]).split(PSF.Regex.GREEDY_WHITESPACE);
        chrs[i] = Positions.chromosomeNumber(line[0], log);
        chrCounts[chrs[i]]++;
        positions[i] = Integer.parseInt(line[1]);
      } else {
        v.add(markerNames[i]);
        chrs[i] = 0;
        chrCounts[chrs[i]]++;
        positions[i] = 0;
      }
    }
    if (!v.isEmpty()) {
      log.reportError("Error - There " + (v.size() == 1 ? "was one" : "were " + v.size())
                      + " markers found in the FinalReport file that were not listed in the file of marker positions; these will be set to missing (chr0:0) and parsing will continue.");
      Files.writeArray(ArrayUtils.toStringArray(v),
                       ext.parseDirectoryOfFile(markerDatabase) + "markersNotInPositionsFile.txt");
      // write markerPositionsNotInReport.txt
      if (!allowExtraAsMissing) {
        return null;
      }
    }

    keys = Sort.getSort2DIndices(chrs, positions);

    new MarkerSet(Sort.getOrdered(markerNames, keys), Sort.getOrdered(chrs, keys),
                  Sort.getOrdered(positions, keys)).serialize(output);

    log.report(ext.getTime() + "\tFinished loading and sorting in " + ext.getTimeElapsed(time));
    log.setLevel(logLevel);

    if (log.getLevel() >= 0) {
      log.report("\nBreakdown of marker counts by chromosome:");
      for (int i = 0; i < Positions.CHR_CODES.length; i++) {
        log.report("chr" + Positions.CHR_CODES[i] + "\t" + ext.addCommas(chrCounts[i]));
      }
      log.report("");
    }

    return keys;
  }

  public static Hashtable<String, String> loadFileToHashString(String filename, Logger log) {
    BufferedReader reader = null;
    String[] line;
    Hashtable<String, String> hash = new Hashtable<>();
    String markerName, chr, position, delimiter, temp;
    byte chrValue;
    int count, countBad, numBlankNames, numBlankChrs, numBlankPositions, numRepeatedNames,
        numInvalidChrs, numInvalidPositions, numIncompleteLines;

    delimiter = Files.determineDelimiter(filename, log);

    count = countBad = 0;
    numBlankNames = numBlankChrs = numBlankPositions = numRepeatedNames = numInvalidChrs = numInvalidPositions = numIncompleteLines = 0;
    try {
      reader = Files.getAppropriateReader(filename);
      while (reader.ready()) {
        temp = reader.readLine();
        if (delimiter.equals(",")) {
          line = ext.splitCommasIntelligently(temp, true, new Logger());
        } else if (temp.indexOf("\t") == -1) {
          line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        } else {
          line = temp.split("\t", -1);
        }
        if (count == 0 && ext.indexOfStr(line[0], Aliases.MARKER_NAMES) >= 0) {

        } else if (line.length < 3) {
          if (countBad < MAX_ERRORS_TO_REPORT) {
            log.report("Error - incomplete line at row " + count + " for marker \"" + line[0]
                       + "\"; line will not be added");
          }
          numIncompleteLines++;
        } else {
          markerName = line[0];
          chr = line[1];
          position = line[2];

          if (markerName.equals("")) {
            if (countBad < MAX_ERRORS_TO_REPORT) {
              log.reportError("Error - blank marker name at line " + count + " of " + filename);
            }
            numBlankNames++;
            countBad++;
          } else if (chr.equals("")) {
            if (countBad < MAX_ERRORS_TO_REPORT) {
              log.reportError("Error - blank chr for marker '" + markerName + "' at line " + count
                              + " of " + filename);
            }
            numBlankChrs++;
            countBad++;
          } else if (position.equals("")) {
            if (countBad < MAX_ERRORS_TO_REPORT) {
              log.reportError("Error - blank position for marker '" + markerName + "' at line "
                              + count + " of " + filename);
            }
            numBlankPositions++;
            countBad++;
          } else {
            if (hash.containsKey(markerName)) {
              log.reportError("Error - marker '" + markerName
                              + "' is already listed in the markerPositions file and is seen again at line "
                              + count + " of " + filename);
              numRepeatedNames++;
              countBad++;
            }
            chrValue = Positions.chromosomeNumber(chr, log);
            if (chrValue < 0 || chrValue > 26) {
              numInvalidChrs++;
              countBad++;
            }
            try {
              Integer.parseInt(position);
            } catch (NumberFormatException nfe) {
              if (countBad < MAX_ERRORS_TO_REPORT) {
                log.reportError("Error - invalid position (" + position + ") for marker '"
                                + markerName + "' at line " + count + " of " + filename);
              }
              numInvalidPositions++;
              countBad++;
            }
          }

          hash.put(markerName, chr + "\t" + position);
        }
        if (countBad == MAX_ERRORS_TO_REPORT) {
          log.reportError("...");
          countBad++;
        }
        count++;
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }

    log.report("\nRead in " + ext.addCommas(count) + " lines from the markerPositions file");
    if (countBad > MAX_ERRORS_TO_REPORT) {
      countBad--;
    }

    if (countBad > 0) {
      log.report("...with a total of " + ext.addCommas(countBad) + " problem"
                 + (countBad == 1 ? "" : "s"));
    }
    if (numIncompleteLines > 0) {
      log.report("...including " + ext.addCommas(numIncompleteLines) + " incomplete line"
                 + (numIncompleteLines == 1 ? "" : "s"));
    }
    log.report("Number of final valid marker positions: " + ext.addCommas(hash.size()));
    if (numBlankNames > 0) {
      log.report("Number of blank marker names: " + ext.addCommas(numBlankNames));
    }
    if (numBlankChrs > 0) {
      log.report("Number of blank chromosomes: " + ext.addCommas(numBlankChrs));
    }
    if (numBlankPositions > 0) {
      log.report("Number of blank marker positions: " + ext.addCommas(numBlankPositions));
    }
    if (numRepeatedNames > 0) {
      log.report("Number of repeated marker names: " + ext.addCommas(numRepeatedNames));
    }
    if (numInvalidChrs > 0) {
      log.report("Number of invalid chromosomes: " + ext.addCommas(numInvalidChrs));
    }
    if (numInvalidPositions > 0) {
      log.report("Number of invalid positions: " + ext.addCommas(numInvalidPositions));
    }
    log.report("");

    return hash;
  }

  public static void generateMarkerPositions(Project proj, String snpTable) {
    Logger log = proj.getLog();
    long time = new Date().getTime();

    Map<String, String> markerToChrPosLinkedMap = new LinkedHashMap<>();
    int markers = 0;
    BufferedReader reader = null;
    try {
      if (!Files.exists(snpTable) && Files.exists(proj.PROJECT_DIRECTORY.getValue() + snpTable)) {
        snpTable = proj.PROJECT_DIRECTORY.getValue() + snpTable;
      }
      String delimiter = Files.determineDelimiter(snpTable, log);
      reader = Files.getAppropriateReader(snpTable);
      int[] indices = ext.indexFactors(SourceFileParser.SNP_TABLE_FIELDS,
                                       reader.readLine().trim().split(delimiter), false, true,
                                       true);
      while (reader.ready()) {
        String[] line = reader.readLine().trim().split(delimiter);
        if (line.length <= indices[0] || line.length < indices[1] || line.length <= indices[2]) {
          log.reportTimeWarning("Skipping line with missing columns: " + ArrayUtils.toStr(line));
          continue;
        }
        String marker = line[indices[0]];
        String chr = Byte.toString(Positions.chromosomeNumber(line[indices[1]]));
        String pos = line[indices[2]];
        if (markerToChrPosLinkedMap.containsKey(marker)) {
          log.reportTimeWarning("Marker " + marker + " was duplicated in " + snpTable
                                + ", only using first occurence.");
        } else {
          markerToChrPosLinkedMap.put(marker, chr + "\t" + pos);
          markers++;
        }
      }
    } catch (FileNotFoundException fnfe) {
      proj.message("Error: file \"" + snpTable + "\" not found in "
                   + proj.PROJECT_DIRECTORY.getValue());
      return;
    } catch (IOException ioe) {
      proj.message("Error reading file \"" + snpTable + "\"");
      return;
    } finally {
      Closeables.closeQuietly(reader);
    }
    lookupMarkers(proj, markerToChrPosLinkedMap, markers, snpTable, log);
    writeMarkerPositions(proj, markerToChrPosLinkedMap, log, time);
  }

  public static void lookupMarkers(Project proj, Map<String, String> markerToChrPosLinkedMap,
                                   int numMarkers, String markerSource, Logger log) {
    int foundRsids = 0;
    int posMismatches = 0;
    List<String> mismatches = Lists.newArrayList();
    mismatches.add(Joiner.on("\t").join(MARKER_POSITIONS_MISMATCHES_HEADER));
    String dbSNP = Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), log).getDBSNP().get();
    VCFFileReader dbSNPReader = null;
    CloseableIterator<VariantContext> dbSNPIterator = null;
    try {
      dbSNPReader = new VCFFileReader(new File(dbSNP));
      dbSNPIterator = dbSNPReader.iterator();
      while (dbSNPIterator.hasNext()) {
        VariantContext vc = dbSNPIterator.next();
        String rsid = vc.getID();
        if (markerToChrPosLinkedMap.containsKey(rsid)) {
          foundRsids++;
          String dbSNPChrPos = Positions.chromosomeNumber(vc.getContig()) + "\t" + vc.getStart();
          if (!dbSNPChrPos.equals(markerToChrPosLinkedMap.get(rsid))) {
            posMismatches++;
            String oldChrPos = markerToChrPosLinkedMap.put(rsid, dbSNPChrPos);
            mismatches.add(rsid + "\t" + oldChrPos + "\t" + dbSNPChrPos);
          }
        }
      }
    } finally {
      if (dbSNPIterator != null) {
        dbSNPIterator.close();
        dbSNPReader.close();
      }
    }
    if (foundRsids == 0) {
      log.report("Markers in " + markerSource
                 + " not identified by rsid, using provided marker positions");
    } else {
      if (numMarkers - foundRsids > 0) {
        log.reportTimeWarning("An rsid match could not be made to dbSNP for "
                              + (numMarkers - foundRsids) + " of " + numMarkers
                              + " total markers in " + markerSource
                              + ", using provided marker positions for these markers.");
      }
      if (posMismatches > 0) {
        log.reportTimeWarning(posMismatches + " markers in " + markerSource
                              + " had positions that did not match the position in dbSNP, using dbSNP position for these markers.");
        String mismatchesReport = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(),
                                                "_mismatches");
        Files.writeIterable(mismatches, mismatchesReport);
        log.report("Mismatches report written to " + mismatchesReport);
      }
    }
  }

  public static void writeMarkerPositions(Project proj, Map<String, String> markerToChrPosLinkedMap,
                                          Logger log, long startTime) {
    PrintWriter writer = null;
    try {
      writer = Files.getAppropriateWriter(proj.MARKER_POSITION_FILENAME.getValue(false, false));
      writer.println("Marker\tChr\tPosition");
      for (Map.Entry<String, String> markerEntry : markerToChrPosLinkedMap.entrySet()) {
        String marker = markerEntry.getKey();
        String chrPos = markerEntry.getValue();
        writer.println(marker + "\t" + chrPos);
      }
    } finally {
      if (writer != null) {
        writer.close();
      }
    }
    log.report("Finished parsing " + proj.MARKER_POSITION_FILENAME.getValue(false, false) + " in "
               + ext.getTimeElapsed(startTime));
  }

  public static void useAlleleLookup(String filename, int alleleCol, String lookupFile, int setFrom,
                                     int setTo) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String> hash;
    String[] header, alleles;

    header = Files.getHeaderOfFile(lookupFile, "\t", new Logger());
    System.out.println("Converting from columns " + header[1 + setFrom * 2 + 0] + "/"
                       + header[1 + setFrom * 2 + 1] + " to columns " + header[1 + setTo * 2 + 0]
                       + "/" + header[1 + setTo * 2 + 1]);
    hash = HashVec.loadFileToHashString(lookupFile, 0, ArrayUtils.arrayOfIndices(header.length),
                                        "\t", true);

    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = Files.openAppropriateWriter(ext.rootOf(filename, false) + "_"
                                           + header[1 + setTo * 2 + 0] + "_"
                                           + header[1 + setTo * 2 + 1] + ".xln");
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      line = ArrayUtils.insertStringAt(header[1 + setTo * 2 + 0], line, alleleCol + 2);
      line = ArrayUtils.insertStringAt(header[1 + setTo * 2 + 1], line, alleleCol + 3);
      writer.println(ArrayUtils.toStr(line));
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        if (hash.containsKey(line[0])) {
          alleles = ArrayUtils.subArray(hash.get(line[0]).split(PSF.Regex.GREEDY_WHITESPACE), 1);
          if (line[alleleCol].equals(alleles[setFrom * 2 + 0])
              && line[alleleCol + 1].equals(alleles[setFrom * 2 + 1])) {
            line = ArrayUtils.insertStringAt(alleles[setTo * 2 + 0], line, alleleCol + 2);
            line = ArrayUtils.insertStringAt(alleles[setTo * 2 + 1], line, alleleCol + 3);
          } else if (line[alleleCol].equals(alleles[setFrom * 2 + 1])
                     && line[alleleCol + 1].equals(alleles[setFrom * 2 + 0])) {
            line = ArrayUtils.insertStringAt(alleles[setTo * 2 + 1], line, alleleCol + 2);
            line = ArrayUtils.insertStringAt(alleles[setTo * 2 + 0], line, alleleCol + 3);
          } else {
            System.err.println("Error - snp '" + line[0] + "' has alleles " + line[alleleCol] + "/"
                               + line[alleleCol + 1] + " in the file and "
                               + alleles[setFrom * 2 + 0] + "/" + alleles[setFrom * 2 + 1]
                               + " in allele lookup table");
            line = ArrayUtils.insertStringAt("XXXXX", line, alleleCol + 2);
            line = ArrayUtils.insertStringAt("XXXXX", line, alleleCol + 3);
          }
        } else {
          System.err.println("Error - snp '" + line[0] + "' not found in allele lookup table");
          line = ArrayUtils.insertStringAt("XXXXX", line, alleleCol + 2);
          line = ArrayUtils.insertStringAt("XXXXX", line, alleleCol + 3);
        }
        writer.println(ArrayUtils.toStr(line));
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    Project proj;
    String filename = null;
    String snpTable = "";
    String fileToConvert = "";
    String lookupFile = "alleleLookup.txt";
    int alleleCol = 6;
    int setFrom = 1;
    int setTo = 2;
    boolean manifestFlag = false;

    String usage = "\n" + "cnv.manage.Markers requires 0-1 arguments\n"
                   + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + "   (2) filename of SNP Table (i.e. snps=Table.csv (not the default))\n"
                   + " OR:\n"
                   + "   (2) use allele lookup to convert a file form forward to TOP strand (i.e. convert=file.txt (not the default))\n"
                   + "   (3) column of A1 in file (i.e. col=" + alleleCol + " (default))\n"
                   + "   (4) allele set to lookup from (i.e. from=" + setFrom
                   + " (default; if 1+8 columns fo 4 pairs, then use indices 0-3))\n"
                   + "   (5) allele set to lookup to (i.e. to=" + setTo + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("snps=")) {
        snpTable = arg.split("=")[1];
        numArgs--;
      } else if (arg.equals("-manifest")) {
        manifestFlag = true;
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    proj = new Project(filename);

    try {
      if (!snpTable.equals("")) {
        if (!Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false))) {
          if (manifestFlag) {
            MarkerBlast.extractMarkerPositionsFromManifest(snpTable, ARRAY.ILLUMINA,
                                                           FILE_SEQUENCE_TYPE.MANIFEST_FILE,
                                                           proj.MARKER_POSITION_FILENAME.getValue(false,
                                                                                                  false),
                                                           Files.determineDelimiter(snpTable,
                                                                                    proj.getLog()),
                                                           proj.getLog());
          } else {
            Markers.generateMarkerPositions(proj, snpTable);
          }
        } else {
          proj.getLog()
              .reportTime("Existing marker positions file already found at "
                          + proj.MARKER_POSITION_FILENAME.getValue(false, false)
                          + "; to reparse marker positions, please remove or rename this file.");
        }
      } else if (!fileToConvert.equals("")) {
        Markers.useAlleleLookup(fileToConvert, alleleCol, lookupFile, setFrom, setTo);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
