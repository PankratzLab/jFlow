// -Xms6G -Xmx6G
// expand to dynamically load/save a certain chunk of markers at a time
package org.pankratzlab.utils.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import javax.annotation.Nullable;
import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.apache.commons.compress.utils.Lists;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.IntVector;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.Sort;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.parsing.DataLine;
import org.pankratzlab.common.parsing.ExplicitIndexedFileColumn;
import org.pankratzlab.common.parsing.FileColumn;
import org.pankratzlab.common.parsing.FileParserFactory;
import org.pankratzlab.common.bioinformatics.Sequence;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.filesys.Segment;
import org.pankratzlab.common.stats.LeastSquares;
import org.pankratzlab.common.stats.LogisticRegression;
import org.pankratzlab.common.stats.RegressionModel;
import org.pankratzlab.utils.filesys.SnpMarkerSet;
import org.pankratzlab.utils.gwas.bgen.BGENReader;
import org.pankratzlab.utils.gwas.bgen.BGENTools;
import org.pankratzlab.utils.gwas.bgen.BGENReader.BGENRecord;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;
import com.google.common.primitives.Ints;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class DosageData implements Serializable {

  public class Trio implements Serializable {

    private static final long serialVersionUID = 1L;

    private final int childIndex;
    private final int fatherIndex;
    private final int motherIndex;

    public Trio(int childIndex, int fatherIndex, int motherIndex) {
      super();
      this.childIndex = childIndex;
      this.fatherIndex = fatherIndex;
      this.motherIndex = motherIndex;
    }

    /**
     * @return the childIndex
     */
    public int getChildIndex() {
      return childIndex;
    }

    /**
     * @return the fatherIndex
     */
    public int getFatherIndex() {
      return fatherIndex;
    }

    /**
     * @return the motherIndex
     */
    public int getMotherIndex() {
      return motherIndex;
    }

  }

  public static final long serialVersionUID = 1L;

  public static final int MACH_MLDOSE_FORMAT = 0;
  public static final int GEN_FORMAT = 1;
  public static final int GWAF_FORMAT = 2;
  public static final int PLINK_FORMAT = 3;
  public static final int MACH_MLPROB_FORMAT = 4;
  public static final int MINIMAC_DOSE_FORMAT = 5;
  public static final int IMPUTE2_DOSE_FORMAT = 6;
  public static final int DATABASE_DOSE_FORMAT = 7;
  public static final int BEAGLE_DOSE_FORMAT = 8;
  public static final int PLINK_BFILE_FORMAT = 9;
  public static final int FREEZE5_FORMAT = 10;
  public static final int VCF_FORMAT = 11;
  public static final int BGEN_FORMAT = 12;

  public static final int MACH_ID_TYPE = 0;
  public static final int SEPARATE_FILE_ID_TYPE = 1;
  public static final int IID_TYPE = 2;
  public static final int FID_IID_TYPE = 3;

  private static final int INDIVIDUAL_DOMINANT_FORMAT = 0;
  private static final int MARKER_DOMINANT_FORMAT = 1;
  private static final int PLINK_FORMAT_INTERNAL = 2;
  private static final int VCF_FORMAT_INTERNAL = 3;
  private static final int BGEN_FORMAT_INTERNAL = 4;

  public static final int CHR_INFO_IN_FILENAME = -2;

  public static final String CHR_REGEX = ".*?chr(\\d\\d?|\\w\\w?).*?";

  public static final String[][] HEADS = {null, null, {"id"}, {"SNP", "A1", "A2"}, null, null,
                                          {"FID", "IID"}};
  public static final String[][] LEADS = {{null, "MLDOSE"}, null, null, null, {null, "MLPROB"},
                                          {null, "DOSE"}, null};
  public static final String[] DELIMITERS = {"\t", ",", " "};

  /**
   * 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 id_type, column index where dosage/probability
   * idtype (unused), values begin, dominance format, number of columns summarizing the data, header
   * row, marker/IID index, A1 index, A2 index, chr index, pos index, delimiter, index for
   * head/lead, min number of digits, max number of digits
   */

  /** 0 1 2 3 4 5 6 7 8 9 10 11 12 13 */
  public static final int[][] PARAMETERS = {{MACH_ID_TYPE, 2, INDIVIDUAL_DOMINANT_FORMAT, 1, 0, 0,
                                             -1, -1, -1, -1, 0, 0, 3, 3}, // .mldose (MACH)
                                            {SEPARATE_FILE_ID_TYPE, 5, MARKER_DOMINANT_FORMAT, 3, 0,
                                             1, 3, 4, 0, 2, 2, 1, 0, 3}, // .gen
                                            {IID_TYPE, 1, INDIVIDUAL_DOMINANT_FORMAT, 1, 1, 0, -1,
                                             -1, -1, -1, 1, 2, 0, 3}, // .fhsR (GWAF)
                                            {FID_IID_TYPE, 3, MARKER_DOMINANT_FORMAT, 2, 1, 0, 1, 2,
                                             -1, -1, 0, 3, 3, 3}, // .dosage (PLINK)
                                            {MACH_ID_TYPE, 2, INDIVIDUAL_DOMINANT_FORMAT, 2, 0, 0,
                                             -1, -1, -1, -1, 0, 4, 3, 3}, // .mlprob (MACH)
                                            {MACH_ID_TYPE, 2, INDIVIDUAL_DOMINANT_FORMAT, 1, 0, 0,
                                             -1, -1, -1, -1, 0, 5, 3, 3}, // .dose (MINIMAC)
                                            {SEPARATE_FILE_ID_TYPE, 5, MARKER_DOMINANT_FORMAT, 3, 0,
                                             1, 3, 4, CHR_INFO_IN_FILENAME, 2, 0, 1, 3, 3}, // .impute2
                                            {FID_IID_TYPE, 3, INDIVIDUAL_DOMINANT_FORMAT, 1, 1, 0,
                                             3, 4, -1, 2, 0, 6, 3, 3}, // .db.xln
                                            {SEPARATE_FILE_ID_TYPE, 3, MARKER_DOMINANT_FORMAT, 1, 1,
                                             0, 1, 2, -1, -1, 2, 1, 4, 4}, // .dose (BEAGLE)
                                            {FID_IID_TYPE, 3, PLINK_FORMAT_INTERNAL, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1},
                                            {IID_TYPE, 1, INDIVIDUAL_DOMINANT_FORMAT, 1, 0, 0, -1,
                                             -1, -1, -1, 0, 0, 3, 3}, // freeze5
                                            {FID_IID_TYPE, 3, VCF_FORMAT_INTERNAL, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1}, // vcf
                                            {FID_IID_TYPE, 3, BGEN_FORMAT_INTERNAL, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1}, // bgen
  };

  private SnpMarkerSet markerSet;

  private String[][] ids;
  private float[][] dosageValues;
  private float[][][] genotypeProbabilities;
  private String[][] alleles;
  private ImmutableList<Trio> trios = null;
  private byte[] chrs;
  private int[] positions;
  private String labelPrepend;
  private boolean empty = false;

  private DosageData() {
    ids = null;
    markerSet = null;
    chrs = null;
    positions = null;
    alleles = null;
    genotypeProbabilities = null;
    dosageValues = null;
  }

  public DosageData(String dosageFile, String idFile, String mapFile, boolean verbose, Logger log) {
    this(dosageFile, idFile, mapFile, determineType(dosageFile), null, null, verbose, log);
  }

  public DosageData(String dosageFile, String idFile, String mapFile, int type, boolean verbose,
                    Logger log) {
    this(dosageFile, idFile, mapFile, type, null, null, verbose, log);
  }

  public DosageData(String dosageFile, String idFile, String mapFile, int type,
                    String regionsToUseFile, String markersToUseFile, boolean verbose, Logger log) {
    this(dosageFile, idFile, mapFile, PARAMETERS[type], regionsToUseFile, markersToUseFile, null,
         verbose, log);
  }

  public DosageData(String dosageFile, String idFile, String mapFile, int[] parameters,
                    String regionsToUseFile, String markersToUseFile, String markerNamePrepend,
                    boolean verbose, Logger log) {
    this(dosageFile, idFile, mapFile, parameters,
         regionsToUseFile == null ? null : loadRegions(regionsToUseFile, log),
         markersToUseFile == null ? null : loadMarkers(markersToUseFile, log), markerNamePrepend,
         verbose, log);
  }

  public DosageData(String dosageFile, String idFile, String mapFile, int[][] regionsToUse,
                    String[] markersToUse, boolean verbose, Logger log) {
    this(dosageFile, idFile, mapFile, PARAMETERS[determineType(dosageFile)], regionsToUse,
         markersToUse, null, verbose, log);
  }

  public DosageData(String dosageFile, String idFile, String mapFile, String regionsToUseFile,
                    String markersToUseFile, boolean verbose, Logger log) {
    this(dosageFile, idFile, mapFile, PARAMETERS[determineType(dosageFile)], regionsToUseFile,
         markersToUseFile, null, verbose, log);
  }

  public DosageData(String dosageFile, String idFile, String mapFile, String regionsToUseFile,
                    String markersToUseFile, String markerNamePrepend, boolean verbose,
                    Logger log) {
    this(dosageFile, idFile, mapFile, PARAMETERS[determineType(dosageFile)], regionsToUseFile,
         markersToUseFile, markerNamePrepend, verbose, log);
  }

  public DosageData(String dosageFile, String idFile, String mapFile, int[][] regions,
                    String[] markers, String markerNamePrepend, boolean verbose, Logger log) {
    this(dosageFile, idFile, mapFile, PARAMETERS[determineType(dosageFile)], regions, markers,
         markerNamePrepend, verbose, log);
  }

  public DosageData(String dosageFile, String idFile, String mapFile, int[][] regions,
                    String[] markers, Map<String, int[]> markerLocationMap,
                    String markerNamePrepend, boolean verbose, Logger log) {
    this(dosageFile, idFile, mapFile, PARAMETERS[determineType(dosageFile)], regions, markers,
         markerLocationMap, markerNamePrepend, verbose, log);
  }

  public DosageData(String dosageFile, String idFile, String mapFile, int[] parameters,
                    int[][] regions, String[] markers, String markerNamePrepend, boolean verbose,
                    Logger log) {
    this(dosageFile, idFile, mapFile, parameters, regions, markers, null, markerNamePrepend,
         verbose, log);
  }

  public DosageData(String dosageFile, String idFile, String mapFile, int[] parameters,
                    int[][] regions, String[] markers, Map<String, int[]> markerLocationMap,
                    String markerNamePrepend, boolean verbose, Logger log) {
    DosageData dd = null;
    if (parameters[2] == PLINK_FORMAT_INTERNAL) {
      dd = loadPlinkBinary(ext.parseDirectoryOfFile(dosageFile), regions, markers,
                           ext.rootOf(dosageFile, true), markerNamePrepend, true, false);
    } else if (parameters[2] == VCF_FORMAT_INTERNAL) {
      dd = loadVCF(dosageFile, mapFile, regions, markers, markerLocationMap, markerNamePrepend,
                   log);
    } else if (parameters[2] == BGEN_FORMAT_INTERNAL) {
      dd = loadBGEN(dosageFile, mapFile, idFile, regions, markers, markerNamePrepend, log);
    }
    if (dd != null) {
      alleles = dd.alleles;
      chrs = dd.chrs;
      positions = dd.positions;
      markerSet = dd.markerSet;
      ids = dd.ids;
      dosageValues = dd.dosageValues;
      genotypeProbabilities = dd.genotypeProbabilities;
      empty = dd.empty;
      labelPrepend = dd.labelPrepend;
      dd = null;
      return;
    }
    BufferedReader reader;
    String[] line;
    String[] markerNames;
    boolean[] markersToKeep;
    Hashtable<String, String> invalids;
    int keepTotal, index, firstMkr;
    if (log == null) {
      log = new Logger();
    }

    markerSet = new SnpMarkerSet(mapFile, false, log);
    markerNames = markerSet.getMarkerNames();

    boolean hasRgns = false;
    boolean hasMkrs = false;
    if (regions != null) {
      hasRgns = true;
    }
    if (markers != null) {
      hasMkrs = true;
    }
    if (hasRgns && hasMkrs) {
      log.reportError("Cannot specify both a regions file and a markers file!");
      return;
    }
    markersToKeep = filterMarkers(markerNames, regions, markers, verbose, log);
    keepTotal = ArrayUtils.booleanArraySum(markersToKeep);

    if (verbose) {
      log.report("Keeping " + keepTotal + " markers out of " + markerNames.length);
    }

    if (keepTotal == 0) {
      empty = true;
      return;
    }

    // int idType = parameters[0]; // unused
    int beginDataCol = parameters[1];
    int dominanceFormat = parameters[2];
    int numDataCols = parameters[3];
    int indexHeaderRow = parameters[4];
    int indexMkrIID = parameters[5];
    int indexA1 = parameters[6];
    int indexA2 = parameters[7];
    int indexChr = parameters[8];
    int indexPos = parameters[9];
    int delimiter = parameters[10];
    // int indexHeadLead = parameters[11]; // unused
    // int numDigitsMin = parameters[12]; // unused
    // int numDigitsMax = parameters[13]; // unused

    ids = HashVec.loadFileToStringMatrix(idFile, false, new int[] {0, 1});
    if (numDataCols == 1) {
      dosageValues = new float[keepTotal][ids.length];
    } else {
      genotypeProbabilities = new float[keepTotal][ids.length][numDataCols];
    }
    if (indexA1 != -1) {
      alleles = new String[keepTotal][2];
    }
    if (indexChr == CHR_INFO_IN_FILENAME) {
      Matcher m = Pattern.compile(CHR_REGEX).matcher(dosageFile);
      byte chr = -1;
      if (m.matches()) {
        chr = (byte) Integer.parseInt(m.group(1));
        if (verbose) {
          String msg = "Warning - the format given expects chromosome number to be part of the file name.  This was determined to be chr{"
                       + chr + "}.";
          log.report(msg);
        }
        chrs = ArrayUtils.byteArray(keepTotal, chr);
      } else {
        if (verbose) {
          String msg = "Error - the format given expects chromosome number to be part of the file name, but no chromosome number was found.  Chromosome information will not be included.";
          log.reportError(msg);
        }
      }
    } else if (indexChr != -1) {
      chrs = new byte[keepTotal];
    }
    if (indexPos != -1) {
      positions = new int[keepTotal];
    }

    if (keepTotal == 0) {
      return;
    }

    invalids = new Hashtable<>();
    try {
      reader = Files.getAppropriateReader(dosageFile);// new BufferedReader(new
                                                      // FileReader(dosageFile));
      if (indexHeaderRow == 1) {
        line = reader.readLine().trim().split(delimiter == 1 ? "," : PSF.Regex.GREEDY_WHITESPACE);
        if (dominanceFormat == INDIVIDUAL_DOMINANT_FORMAT) {
          for (int i = 0; i < markerNames.length; i++) {
            if (!markerNames[i].equals(line[beginDataCol + numDataCols * i])) {
              String msg = "Error - mismatched name at marker " + (i + 1) + " of " + dosageFile
                           + "; expecting " + markerNames[i] + " given map file " + mapFile
                           + ", found " + line[beginDataCol + numDataCols * i];
              log.reportError(msg);
              reader.close();
              return;
            }
          }
        } else if (dominanceFormat == MARKER_DOMINANT_FORMAT) {
          if (numDataCols != 2) {
            if (verbose) {
              String msg = "Warning - ignoring the header with IDs in file " + dosageFile
                           + " because it does not contain 2 columns for each individual";
              log.reportError(msg);
            }
          } else {
            for (int i = 0; i < ids.length; i++) {
              if (!ids[i][0].equals(line[beginDataCol + numDataCols * i + 0])
                  || !ids[i][1].equals(line[beginDataCol + numDataCols * i + 1])) {
                String msg = "Error - mismatched IDs at individual " + (i + 1) + " of " + dosageFile
                             + "; expecting " + ids[i][0] + "," + ids[i][1] + " given id file "
                             + idFile + ", found " + line[beginDataCol + numDataCols * i + 0] + ","
                             + line[beginDataCol + numDataCols * i + 1];
                log.reportError(msg);
                reader.close();
                return;
              }
            }
          }
        }
      }

      if (dominanceFormat == MARKER_DOMINANT_FORMAT) {
        firstMkr = 0;
        while (!markersToKeep[firstMkr]) {
          firstMkr++;
        }
        for (int i = 0; i < firstMkr; i++) {
          reader.readLine(); // skip to first marker
        }
        index = -1;
        for (int i = firstMkr; i < markerNames.length; i++) {
          String temp = reader.readLine();
          if (temp == null) {
            int rem = markerNames.length - i;
            log.reportError((rem > 0 ? "Error" : "Warning")
                            + " - Reached end of dosage file at marker " + i + ": " + markerNames[i]
                            + ", " + (rem) + " remaining expected in file: " + dosageFile);
            reader.close();
            return;
          }
          if (!markersToKeep[i]) {
            continue;
          }
          line = temp.trim().split(delimiter == 1 ? "," : PSF.Regex.GREEDY_WHITESPACE);
          if (!markerNames[i].equals(line[indexMkrIID])) {
            if (verbose) {
              String msg = "Error - mismatched name at marker " + (i + 1) + " of " + dosageFile
                           + "; expecting " + markerNames[i] + " given map file " + mapFile
                           + ", found " + line[indexMkrIID];
              log.reportError(msg);
            }
            reader.close();
            return;
          }
          index++;
          if (indexA1 != -1) {
            if (line[indexA1].length() > 1 || !Sequence.validBase(line[indexA1])) {
              if (verbose) {
                String msg = "Warning - invalid allele ('" + line[indexA1] + "') at marker "
                             + markerNames[i];
                log.reportError(msg);
              }
            }
            alleles[index][0] = line[indexA1];
          }
          if (indexA2 != -1) {
            if (line[indexA2].length() > 1 || !Sequence.validBase(line[indexA2])) {
              if (verbose) {
                String msg = "Warning - invalid allele ('" + line[indexA2] + "') at marker "
                             + markerNames[i];
                log.reportError(msg);
              }
            }
            alleles[index][1] = line[indexA2];
          }
          if (indexChr >= 0) {
            try {
              chrs[index] = Byte.parseByte(line[indexChr]);
            } catch (NumberFormatException nfe) {
              chrs[index] = -1;
              if (!invalids.containsKey(line[indexChr])) {
                if (verbose) {
                  String msg = "Warning - invalid chromosome number ('" + line[indexChr]
                               + "'), first seen at marker " + markerNames[i];
                  log.reportError(msg);
                }
                invalids.put(line[indexChr], "");
              }
            }
          }
          if (indexPos != -1) {
            try {
              positions[index] = Integer.parseInt(line[indexPos]);
            } catch (NumberFormatException nfe) {
              positions[index] = -1;
              if (!invalids.containsKey(line[indexPos])) {
                if (verbose) {
                  String msg = "Warning - invalid genome position ('" + line[indexPos]
                               + "') for marker " + markerNames[i];
                  log.reportError(msg);
                }
                invalids.put(line[indexPos], "");
              }
            }
          }

          if (line.length - beginDataCol != ids.length * numDataCols) {
            String msg = "Error - mismatched number of elements in line " + (i + 1 + indexHeaderRow)
                         + " of " + dosageFile + "; expecting " + ids.length + "*" + numDataCols
                         + "+" + beginDataCol + "[=" + (ids.length * numDataCols + beginDataCol)
                         + "], found " + line.length;
            log.reportError(msg);
            System.exit(1);
          }
          if (numDataCols == 1) {
            for (int j = 0; j < ids.length; j++) {
              dosageValues[index][j] = ext.isMissingValue(line[beginDataCol
                                                               + j]) ? Float.NaN
                                                                     : Float.parseFloat(line[beginDataCol
                                                                                             + j]);
            }
          } else {
            for (int j = 0; j < ids.length; j++) {
              genotypeProbabilities[index][j][0] = ext.isMissingValue(line[beginDataCol + j
                                                                                          * numDataCols
                                                                           + 0]) ? Float.NaN
                                                                                 : Float.parseFloat(line[beginDataCol
                                                                                                         + j
                                                                                                           * numDataCols
                                                                                                         + 0]);
              genotypeProbabilities[index][j][1] = ext.isMissingValue(line[beginDataCol + j
                                                                                          * numDataCols
                                                                           + 1]) ? Float.NaN
                                                                                 : Float.parseFloat(line[beginDataCol
                                                                                                         + j
                                                                                                           * numDataCols
                                                                                                         + 1]);
              if (numDataCols == 3) {
                genotypeProbabilities[index][j][2] = ext.isMissingValue(line[beginDataCol + j
                                                                                            * numDataCols
                                                                             + 2]) ? Float.NaN
                                                                                   : Float.parseFloat(line[beginDataCol
                                                                                                           + j
                                                                                                             * numDataCols
                                                                                                           + 2]);
                if (Math.abs((1 - genotypeProbabilities[index][j][1]
                              - genotypeProbabilities[index][j][0])
                             - genotypeProbabilities[index][j][2]) > 0.01) {
                  String msg = "Error: P(BB) does not equal [ 1 - P(AA) - P(AB) ] for individual "
                               + ids[j][0] + "," + ids[j][1] + " at marker " + markerNames[i]
                               + " which is line " + (i + 1 + indexHeaderRow) + " of " + dosageFile
                               + ": " + line[beginDataCol + j * numDataCols + 0] + " "
                               + line[beginDataCol + j * numDataCols + 1] + " "
                               + line[beginDataCol + j * numDataCols + 2];
                  log.reportError(msg);
                }
              }
            }
          }

          if (index == keepTotal - 1) {
            // found all markers, no need to continue reading
            break;
          }
        }
      } else if (dominanceFormat == INDIVIDUAL_DOMINANT_FORMAT) {
        for (int i = 0; i < ids.length; i++) {
          line = reader.readLine().trim().split(delimiter == 1 ? "," : PSF.Regex.GREEDY_WHITESPACE);
          if (line.length - beginDataCol != markerNames.length * numDataCols) {
            String msg = "Error - mismatched number of elements in line " + (i + 1 + indexHeaderRow)
                         + " of " + dosageFile + "; expecting " + markerNames.length + "*"
                         + numDataCols + "+" + beginDataCol + "[="
                         + (markerNames.length * numDataCols + beginDataCol) + "], found "
                         + line.length;
            log.reportError(msg);
            System.exit(1);
          }
          if (numDataCols == 1) {
            index = -1;
            for (int j = 0; j < markerNames.length; j++) {
              if (!markersToKeep[j]) {
                continue;
              }
              index++;
              dosageValues[index][i] = ext.isMissingValue(line[beginDataCol
                                                               + j]) ? Float.NaN
                                                                     : Float.parseFloat(line[beginDataCol
                                                                                             + j]);
              if (index == keepTotal - 1) {
                break;
              }
            }
          } else {
            index = -1;
            for (int j = 0; j < markerNames.length; j++) {
              if (!markersToKeep[j]) {
                continue;
              }
              index++;
              genotypeProbabilities[index][i][0] = ext.isMissingValue(line[beginDataCol + j
                                                                                          * numDataCols
                                                                           + 0]) ? Float.NaN
                                                                                 : Float.parseFloat(line[beginDataCol
                                                                                                         + j
                                                                                                           * numDataCols
                                                                                                         + 0]);
              genotypeProbabilities[index][i][1] = ext.isMissingValue(line[beginDataCol + j
                                                                                          * numDataCols
                                                                           + 1]) ? Float.NaN
                                                                                 : Float.parseFloat(line[beginDataCol
                                                                                                         + j
                                                                                                           * numDataCols
                                                                                                         + 1]);
              if (numDataCols == 3) {
                genotypeProbabilities[index][i][2] = ext.isMissingValue(line[beginDataCol + j
                                                                                            * numDataCols
                                                                             + 2]) ? Float.NaN
                                                                                   : Float.parseFloat(line[beginDataCol
                                                                                                           + j
                                                                                                             * numDataCols
                                                                                                           + 2]);
                if (Math.abs((1 - genotypeProbabilities[index][i][1]
                              - genotypeProbabilities[index][i][0])
                             - genotypeProbabilities[index][i][2]) > 0.01) {
                  String msg = "Error: P(BB) does not equal [ 1 - P(AA) - P(AB) ] for individual "
                               + ids[i][0] + "," + ids[i][1] + " at marker " + markerNames[j]
                               + " which is line " + (i + 1 + indexHeaderRow) + " of " + dosageFile
                               + ": " + line[beginDataCol + j * numDataCols + 0] + " "
                               + line[beginDataCol + j * numDataCols + 1] + " "
                               + line[beginDataCol + j * numDataCols + 2];
                  log.reportError(msg);
                }
              }
              if (index == keepTotal - 1) {
                break;
              }
            }
          }
        }
      }
      reader.close();

      markerNames = markerSet.getMarkerNames();
      if (markerNamePrepend != null && !"".equals(markerNamePrepend)) {
        this.labelPrepend = markerNamePrepend;
        for (int i = 0; i < markerNames.length; i++) {
          markerNames[i] = (markerNamePrepend.endsWith("_") ? markerNamePrepend
                                                            : markerNamePrepend + "_")
                           + markerNames[i];
        }
        markerSet.convertMarkerNamesToRSnumbers(markerNames, verbose, log); // resets/recreates the
        // SnpMarkerSet's internal
        // RSnumber/non-RS marker
        // lists
      }

      if (markerSet.getAlleles() == null && alleles != null) {
        markerSet.setAlleles(alleles);
      }

    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dosageFile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dosageFile + "\"");
      System.exit(2);
    }
  }

  public void analyze(String phenoFile, String phenoMissingValue, String snpList, boolean verbose,
                      Logger log) {
    PrintWriter writer, w2;
    String[] line;
    Hashtable<String, String> hash;
    HashSet<String> snps;
    int count;
    String[] traits, markerNames;
    boolean[] use, analyze;
    byte[] chrs;
    int[] positions;
    double[] deps;
    double[][] indeps;
    boolean logistic;
    RegressionModel model;
    String[][] alleles;
    double[] betas, stderrs, pvals, stats;

    traits = Files.getHeaderOfFile(phenoFile, "\t", log);
    hash = HashVec.loadFileToHashString(phenoFile, new int[] {0, 1},
                                        Arrays.copyOfRange(ArrayUtils.arrayOfIndices(traits.length),
                                                           2, traits.length),
                                        false, "\t", true, false);
    traits = ArrayUtils.subArray(traits, 2);

    markerNames = markerSet.getMarkerNames();
    chrs = markerSet.getChrs();
    positions = markerSet.getPositions();
    alleles = markerSet.getAlleles();
    analyze = ArrayUtils.booleanArray(markerNames.length, true);
    if (snpList != null) {
      snps = HashVec.loadFileToHashSet(snpList, false);
      for (int i = 0; i < markerNames.length; i++) {
        if (!snps.contains(markerNames[i])) {
          analyze[i] = false;
        }
      }
    }

    use = ArrayUtils.booleanArray(ids.length, true);
    for (int i = 0; i < ids.length; i++) {
      if (hash.containsKey(ids[i][0] + "\t" + ids[i][1])) {
        line = hash.get(ids[i][0] + "\t" + ids[i][1]).split(PSF.Regex.GREEDY_WHITESPACE);
        for (int j = 0; j < line.length; j++) {
          if (!ext.isValidDouble(line[j]) && !line[j].equals(phenoMissingValue)) {
            use[i] = false;
          }
        }
      } else {
        use[i] = false;
      }
    }
    deps = new double[ArrayUtils.booleanArraySum(use)];
    indeps = new double[deps.length][traits.length];
    log.report("There are " + deps.length + " rows with complete data", true, verbose);
    count = 0;
    for (int i = 0; i < ids.length; i++) {
      if (use[i]) {
        line = hash.get(ids[i][0] + "\t" + ids[i][1]).split(PSF.Regex.GREEDY_WHITESPACE);
        deps[count] = Double.parseDouble(line[0]);
        for (int j = 1; j < traits.length; j++) {
          indeps[count][j] = Double.parseDouble(line[j]);
        }
        count++;
      }
    }
    logistic = RegressionModel.isBinaryTrait(ArrayUtils.toStr(deps)
                                                       .split(PSF.Regex.GREEDY_WHITESPACE),
                                             log);
    log.report("Running a " + (logistic ? "logistic" : "linear") + " model for trait '" + traits[0]
               + "'", true, verbose);
    try {
      // writer = new PrintWriter(new FileWriter(ext.rootOf(phenoFile,
      // false)+(snpList==null?"":"_"+ext.rootOf(snpList,
      // false))+".results."+(logistic?"logistic":"linear")));
      // w2 = new PrintWriter(new FileWriter(ext.rootOf(phenoFile,
      // false)+(snpList==null?"":"_"+ext.rootOf(snpList, false))+".se.metal"));
      writer = Files.openAppropriateWriter(ext.rootOf(phenoFile, false) + ".results."
                                           + (logistic ? "logistic" : "linear"));
      w2 = Files.openAppropriateWriter(ext.rootOf(phenoFile, false) + ".se.metal");
      String[] arr = logistic ? PSF.Plink.LOGISTIC_SE_HEADER : PSF.Plink.LINEAR_SE_HEADER;
      line = Arrays.copyOf(arr, arr.length);
      line[1] = line[1] + "      ";
      line[2] = line[1] + "      ";
      writer.println(ArrayUtils.toStr(line));
      // w2.println("MARKER\tREF\tOTHER\tN\tDIR\tPVALUE\tbeta\tSE");
      w2.println("MarkerName\tAllele1\tAllele2\tWeight\tDirection\tP-value\tEffect\tStdErr");
      // public static final String[] LOGISTIC_SE_HEADER = {"CHR", "SNP", "BP", "A1", "TEST",
      // "NMISS", "OR", "SE", "L95", "U95", "STAT", "P"};
      for (int i = 0; i < markerNames.length; i++) {
        if (analyze[i]) {
          count = 0;
          for (int j = 0; j < ids.length; j++) {
            if (use[j]) {
              indeps[count][0] = dosageValues[i][j];
              count++;
            }
          }
          model = logistic ? new LogisticRegression(deps, indeps, false, false)
                           : new LeastSquares(deps, indeps, false, false);
          betas = model.getBetas();
          stderrs = model.getSEofBs();
          pvals = model.getSigs();
          stats = model.getStats();
          int sigfig = 4;
          // System.err.println(betas.length+"\t"+traits.length);
          if (betas.length != traits.length + 1) {
            writer.println(chrs[i] + "\t" + markerNames[i] + "\t" + positions[i] + "\t"
                           + alleles[i][0] + "\tADD\t.\t.\t.\t.\t.\t.\t.");
            w2.println(markerNames[i] + "\t" + alleles[i][0] + "\t" + alleles[i][1] + "\t"
                       + deps.length + "\t.\t.\t.\t.");
          } else {
            writer.println(chrs[i] + "\t" + markerNames[i] + "\t" + positions[i] + "\t"
                           + alleles[i][0] + "\tADD\t" + deps.length + "\t"
                           + ext.formDeci(betas[1], sigfig, true) + "\t"
                           + ext.formDeci(stderrs[1], sigfig, true) + "\t.\t.\t"
                           + (stats[1] + "    ").substring(0, 6).trim() + "\t"
                           + ext.prettyP(pvals[1], sigfig, 4, 3, true));
            w2.println(markerNames[i] + "\t" + alleles[i][0] + "\t" + alleles[i][1] + "\t"
                       + deps.length + "\t" + (betas[1] == 0 ? 0 : (betas[1] > 0 ? "+" : "-"))
                       + "\t" + ext.prettyP(pvals[1], sigfig, 4, 3, true) + "\t"
                       + ext.formDeci(betas[1], 6, true) + "\t"
                       + ext.formDeci(stderrs[1], 6, true));
            for (int j = 1; j < traits.length; j++) {
              writer.println(chrs[i] + "\t" + markerNames[i] + "\t" + positions[i] + "\t"
                             + alleles[i][0] + "\t" + traits[j] + "\t" + deps.length + "\t"
                             + ext.formDeci(betas[1 + j], sigfig, true) + "\t"
                             + ext.formDeci(stderrs[1 + j], sigfig, true) + "\t.\t.\t"
                             + (stats[1 + j] + "     ").substring(0, 6).trim() + "\t"
                             + ext.prettyP(pvals[1 + j], sigfig, 4, 3, true));
            }
          }
          writer.flush();
        }
      }
      writer.close();
      w2.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(phenoFile, false)
                         + (snpList == null ? "" : "_" + ext.rootOf(snpList, false)) + ".results."
                         + (logistic ? "logistic" : "linear"));
      e.printStackTrace();
    }
  }

  /**
   * Compute the "Best Guess" dosage value, which is the index of the largest genotype probability
   * (i.e. if the 0 index probability is highest, the dosage value becomes "0"; likewise for index 1
   * and, if it exists, index 2).
   * 
   * @param missThresh If all genoProbs are less-than or equal to this threshold, set to missing
   * @param log
   */
  public void computeDosageValueFromBestGenotype(float missThresh, float missingValue, Logger log) {
    if (genotypeProbabilities == null) {
      log.reportError("Error - cannot compute dosage values from genotype probabilities, if there are no genotype probabilities!");
    }
    dosageValues = new float[genotypeProbabilities.length][genotypeProbabilities[0].length];
    log.reportTime("Computing \"Best Guess\" dosage values from "
                   + genotypeProbabilities[0][0].length + " genotype probability values.");
    for (int i = 0; i < dosageValues.length; i++) {
      for (int j = 0; j < ids.length; j++) {
        float geno1, geno2, geno3;
        geno1 = genotypeProbabilities[i][j][0];
        geno2 = genotypeProbabilities[i][j][1];
        geno3 = genotypeProbabilities[i][j].length > 2 ? genotypeProbabilities[i][j][2]
                                                       : (1 - (geno1 + geno2));
        if (geno1 < missThresh && geno2 < missThresh && geno3 < missThresh) {
          dosageValues[i][j] = missingValue;
        } else {
          if (geno1 > geno2 && geno1 > geno3) {
            dosageValues[i][j] = 0;
          } else if (geno2 > geno1 && geno2 > geno3) {
            dosageValues[i][j] = 1;
          } else if (geno3 > geno1 && geno3 > geno2) {
            dosageValues[i][j] = 2;
          } else {
            // tie
            dosageValues[i][j] = missingValue;
          }
        }
      }
    }
  }

  public void computeDosageValues(Logger log) {
    if (genotypeProbabilities == null) {
      log.reportError("Error - cannot compute dosage values from genotype probabilities, if there are no genotype probabilities!");
      // this will error on the next line
    }
    dosageValues = new float[genotypeProbabilities.length][genotypeProbabilities[0].length];
    for (int i = 0; i < dosageValues.length; i++) {
      for (int j = 0; j < ids.length; j++) {
        dosageValues[i][j] = genotypeProbabilities[i][j][0] * 2
                             + genotypeProbabilities[i][j][1] * 1;
      }
    }
  }

  public float[][] getDosageValues() {
    return dosageValues;
  }

  public String[][] getIds() {
    return ids;
  }

  public List<Trio> getTrios() {
    if (trios == null) generateTrios();
    return trios;
  }

  private void generateTrios() {
    Map<String, Integer> indexMap = IntStream.range(0, ids.length).boxed()
                                             .collect(ImmutableMap.toImmutableMap(i -> formFidIid(ids[i][0],
                                                                                                  ids[i][1]),
                                                                                  Integer::intValue));
    ImmutableList.Builder<Trio> triosBuilder = ImmutableList.builder();
    for (int i = 0; i < ids.length; i++) {
      String[] id = ids[i];
      if (id.length > PSF.Plink.FAM_MO_INDEX && id.length > PSF.Plink.FAM_FA_INDEX) {
        String famId = id[PSF.Plink.FAM_FID_INDEX];
        String fatherId = id[PSF.Plink.FAM_FA_INDEX];
        String motherId = id[PSF.Plink.FAM_MO_INDEX];

        Integer fatherIndex = indexMap.get(formFidIid(famId, fatherId));
        Integer motherIndex = indexMap.get(formFidIid(famId, motherId));
        if (fatherIndex == null && motherIndex == null && famId.equals(id[1])) {
          fatherIndex = indexMap.get(formFidIid(fatherId, fatherId));
          motherIndex = indexMap.get(formFidIid(motherId, motherId));
        }
        if (fatherIndex != null && motherIndex != null) {
          triosBuilder.add(new Trio(i, fatherIndex, motherIndex));
        }
      }
    }
    trios = triosBuilder.build();
  }

  private static String formFidIid(String fid, String iid) {
    return fid + "\t" + iid;
  }

  public SnpMarkerSet getMarkerSet() {
    return markerSet;
  }

  public boolean isEmpty() {
    return empty;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  public void writeToFile(String filename, String mapOut, boolean allowIncompleteList,
                          boolean writeNaNsAsPeriods, int[] parameters, Logger log) {
    writeToFile(filename, mapOut, null, null, allowIncompleteList, writeNaNsAsPeriods, parameters,
                false, Float.NaN, log);
  }

  public void writeToFile(String filename, String mapOut, String extractMarkers, String regionsFile,
                          boolean allowIncompleteList, boolean writeNaNsAsPeriods, int format,
                          Logger log) {
    String[] markersToKeep = extractMarkers == null ? null
                                                    : HashVec.loadFileToStringArray(extractMarkers,
                                                                                    false,
                                                                                    new int[] {0},
                                                                                    true, false,
                                                                                    "\t");
    int[][] regions;
    if (regionsFile == null) {
      regions = null;
    } else {
      String[] rgns = HashVec.loadFileToStringArray(regionsFile, false, new int[] {0}, true, false,
                                                    "\t");
      regions = new int[rgns.length][];
      for (int i = 0; i < rgns.length; i++) {
        regions[i] = Positions.parseUCSClocation(rgns[i]);
      }
    }
    writeToFile(filename, mapOut, markersToKeep, regions, allowIncompleteList, writeNaNsAsPeriods,
                PARAMETERS[format], false, Float.NaN, log);
  }

  public void writeToFile(String filename, String mapOut, String extractMarkers, String regionsFile,
                          int format, Logger log) {
    writeToFile(filename, mapOut, extractMarkers, regionsFile, false, true, format, log);
  }

  public void writeToFile(String filename, String mapOut, String extractMarkers, String regionsFile,
                          Logger log) {
    writeToFile(filename, mapOut, extractMarkers, regionsFile, determineType(filename), log);
  }

  public void writeToFile(String filename, String mapOut, String[] markersToKeep,
                          int[][] regionsToKeep, boolean allowIncompleteList,
                          boolean writeNaNsAsPeriods, int[] parameters, boolean bestGuessOutput,
                          float bestGuessThreshold, Logger log) {
    PrintWriter writer;
    String[] line, markerNames;
    String delimiter;
    HashSet<String> keeps;
    String root;
    SnpMarkerSet newMarkerSet;
    if (log == null) {
      log = new Logger();
    }

    int dataOutType = determineType(filename);
    if (dataOutType == PLINK_BFILE_FORMAT) {
      writeToPlinkBinary(ext.parseDirectoryOfFile(filename), ext.rootOf(filename, true),
                         bestGuessOutput, bestGuessThreshold, log);
      return;
    } else if (dataOutType == VCF_FORMAT) {
      throw new UnsupportedOperationException("Writing to a .vcf file is not currently supported. Please choose a different output type and try again.");
    } else if (dataOutType == BGEN_FORMAT) {
      throw new UnsupportedOperationException("Writing to a .bgen file is not currently supported. Please choose a different output type and try again.");
    }

    if (markersToKeep == null && regionsToKeep == null) {
      keeps = null;
      markerSet.writeToFile(mapOut, SnpMarkerSet.determineType(mapOut), log);
    } else {
      if (markersToKeep != null && regionsToKeep != null) {
        log.reportError("Error - cannot specify both markersToKeep and regionsToKeep.");
        return;
      } else if (markersToKeep != null) {
        keeps = HashVec.loadToHashSet(markersToKeep);
        root = ext.rootOf(filename, false);
        if (mapOut == null) {
          mapOut = root + ".map";
        }
        newMarkerSet = markerSet.trim(markersToKeep, allowIncompleteList, false, log);
        if (!allowIncompleteList && newMarkerSet == null) {
          log.reportError("Error - failed to find all of the subset of markers, in the source file '"
                          + filename + "'");
          return;
        }
        newMarkerSet.writeToFile(mapOut, SnpMarkerSet.determineType(mapOut), log);
        try {
          writer = Files.openAppropriateWriter(root + ".ids.fam");
          for (String[] id : ids) {
            writer.println(id[0] + "\t" + id[1]);
          }
          writer.close();
        } catch (Exception e) {
          System.err.println("Error writing to " + root + ".ids.fam");
          e.printStackTrace();
        }
      } else/* if (regionsToKeep != null) */ {
        root = ext.rootOf(filename, false);
        if (mapOut == null) {
          mapOut = root + ".map";
        }
        newMarkerSet = markerSet.trim(regionsToKeep, false, log);
        keeps = HashVec.loadToHashSet(newMarkerSet.getMarkerNames());
        newMarkerSet.writeToFile(mapOut, SnpMarkerSet.determineType(mapOut), log);
        try {
          writer = Files.openAppropriateWriter(root + ".ids.fam");
          for (String[] id : ids) {
            writer.println(id[0] + "\t" + id[1]);
          }
          writer.close();
        } catch (Exception e) {
          System.err.println("Error writing to " + root + ".ids.fam");
          e.printStackTrace();
        }
      }
    }

    markerNames = markerSet.getMarkerNames();
    if (alleles == null) {
      alleles = markerSet.getAlleles();
    }
    if (chrs == null) {
      chrs = markerSet.getChrs();
      positions = markerSet.getPositions();
    }
    delimiter = DELIMITERS[parameters[10]];

    if (parameters[3] == 1 && dosageValues == null) {
      if (bestGuessOutput) {
        computeDosageValueFromBestGenotype(bestGuessThreshold, Float.NaN, log);
      } else {
        computeDosageValues(log);
      }
    }

    try {
      writer = Files.getAppropriateWriter(filename);
      if (parameters[4] == 1) {
        writer.print(ArrayUtils.toStr(HEADS[parameters[11]], delimiter));

        if (parameters[2] == MARKER_DOMINANT_FORMAT) {
          if (parameters[3] == 2 && parameters[0] == FID_IID_TYPE) {
            for (String[] id : ids) {
              writer.print(delimiter + id[0] + delimiter + id[1]);
            }
          } else if (parameters[3] == 1 && parameters[0] == IID_TYPE) {
            for (String[] id : ids) {
              writer.print(delimiter + id[1]);
            }
          } else {
            log.reportError("Error - don't know how to list IDs when there "
                            + (parameters[3] == 1 ? "is one column"
                                                  : "are " + parameters[3] + " columns")
                            + " for dosage inforation and the ID type is '" + parameters[2] + "'");
            System.exit(1);
          }
        } else if (parameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
          for (String markerName : markerNames) {
            if ((markersToKeep == null && regionsToKeep == null) || keeps.contains(markerName)) {
              for (int j = 0; j < parameters[3]; j++) {
                writer.print(delimiter + markerName);
              }
            }
          }
        }
        writer.println();
      }

      if (parameters[2] == MARKER_DOMINANT_FORMAT) {
        for (int i = 0; i < markerNames.length; i++) {
          if ((markersToKeep == null && regionsToKeep == null) || keeps.contains(markerNames[i])) {
            line = LEADS[parameters[11]] == null ? new String[parameters[1]]
                                                 : LEADS[parameters[11]];
            line[parameters[5]] = markerNames[i];
            if (parameters[6] >= 0) {
              if (alleles == null) {
                log.reportError("Error - file format requires alleles and none were supplied via the map file or the source dosage file");
                System.exit(1);
              }
              line[parameters[6]] = alleles[i][0] + "";
            }
            if (parameters[7] >= 0) {
              line[parameters[7]] = alleles[i][1] + "";
            }
            if (parameters[8] >= 0) {
              if (chrs == null) {
                log.reportError("Error - file format requires chrs/positions and none were supplied via the map file or the source dosage file");
                System.exit(1);
              }
              line[parameters[8]] = chrs[i] == -1 ? "-" : chrs[i] + "";
            }
            if (parameters[9] >= 0) {
              line[parameters[9]] = positions[i] + "";
            }
            writer.print(ArrayUtils.toStr(line, delimiter));

            if (parameters[3] == 1) {
              for (int j = 0; j < ids.length; j++) {
                writer.print(delimiter + (Float
                                               .isNaN(dosageValues[i][j]) ? (writeNaNsAsPeriods ? "." : Float.NaN) : ext.formDeci(dosageValues[i][j], parameters[13], parameters[12] == parameters[13])));
              }
            } else {
              if (genotypeProbabilities == null) {
                log.reportError("Error - cannot write genotype probabilities if they were never read in!");
                System.exit(1);
              }
              for (int j = 0; j < ids.length; j++) {
                writer.print(delimiter + (Float
                                               .isNaN(genotypeProbabilities[i][j][0]) ? (writeNaNsAsPeriods ? "." : Float.NaN) : ext.formDeci(genotypeProbabilities[i][j][0], parameters[13], parameters[12] == parameters[13])));
                writer.print(delimiter + (Float
                                               .isNaN(genotypeProbabilities[i][j][1]) ? (writeNaNsAsPeriods ? "." : Float.NaN) : ext.formDeci(genotypeProbabilities[i][j][1], parameters[13], parameters[12] == parameters[13])));
                if (parameters[3] == 3) {
                  if (genotypeProbabilities[i][j].length > 2) {
                    writer.print(delimiter + (Float
                                                   .isNaN(genotypeProbabilities[i][j][2]) ? (writeNaNsAsPeriods ? "." : Float.NaN) : ext.formDeci(genotypeProbabilities[i][j][2], parameters[13], parameters[12] == parameters[13])));
                  } else {
                    writer.print(delimiter
                                 + (Float.isNaN(genotypeProbabilities[i][j][0])
                                    || Float.isNaN(genotypeProbabilities[i][j][1]) ? (writeNaNsAsPeriods ? "." : Float.NaN) : ext.formDeci(1 - genotypeProbabilities[i][j][1] - genotypeProbabilities[i][j][0], parameters[13], parameters[12] == parameters[13])));
                  }
                }
              }
            }
            writer.println();
          }
        }
      } else if (parameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
        for (int i = 0; i < ids.length; i++) {
          line = LEADS[parameters[11]] == null ? new String[parameters[3]] : LEADS[parameters[11]];
          if (parameters[0] == MACH_ID_TYPE) {
            line[parameters[5]] = ids[i][0] + "->" + ids[i][1];
          } else if (parameters[0] == IID_TYPE) {
            line[parameters[5]] = ids[i][1];
          } else if (parameters[0] == FID_IID_TYPE) {
            line[parameters[5]] = ids[i][0] + delimiter + ids[i][1];
          } else {
            log.reportError("Error - ID type has not been defined");
            System.exit(1);
          }
          writer.print(ArrayUtils.toStr(line, delimiter));

          if (parameters[3] == 1) {
            for (int j = 0; j < markerNames.length; j++) {
              if ((markersToKeep == null && regionsToKeep == null)
                  || keeps.contains(markerNames[j])) {
                writer.print(delimiter + (Float
                                               .isNaN(dosageValues[j][i]) ? (writeNaNsAsPeriods ? "." : Float.NaN) : ext.formDeci(dosageValues[j][i], parameters[13], parameters[12] == parameters[13])));
              }
            }
          } else {
            for (int j = 0; j < markerNames.length; j++) {
              if ((markersToKeep == null && regionsToKeep == null)
                  || keeps.contains(markerNames[j])) {
                writer.print(delimiter + (Float
                                               .isNaN(genotypeProbabilities[j][i][0]) ? (writeNaNsAsPeriods ? "." : Float.NaN) : ext.formDeci(genotypeProbabilities[j][i][0], parameters[13], parameters[12] == parameters[13])));
                writer.print(delimiter + (Float
                                               .isNaN(genotypeProbabilities[j][i][1]) ? (writeNaNsAsPeriods ? "." : Float.NaN) : ext.formDeci(genotypeProbabilities[j][i][1], parameters[13], parameters[12] == parameters[13])));
                if (parameters[3] == 3) {
                  if (genotypeProbabilities[j][i].length > 2) {
                    writer.print(delimiter + (Float
                                                   .isNaN(genotypeProbabilities[j][i][2]) ? (writeNaNsAsPeriods ? "." : Float.NaN) : ext.formDeci(genotypeProbabilities[j][i][2], parameters[13], parameters[12] == parameters[13])));
                  } else {
                    writer.print(delimiter
                                 + (Float.isNaN(genotypeProbabilities[j][i][0])
                                    || Float.isNaN(genotypeProbabilities[j][i][1]) ? (writeNaNsAsPeriods ? "." : Float.NaN) : ext.formDeci(1 - genotypeProbabilities[j][i][1] - genotypeProbabilities[j][i][0], parameters[13], parameters[12] == parameters[13])));
                  }
                }
              }
            }
          }
          writer.println();
        }
      }

      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing " + filename);
      e.printStackTrace();
    }
  }

  public void writeToPlinkBinary(String dir, String plinkRoot, boolean bestGuessOutput,
                                 float bestGuessThreshold, Logger log) {
    PrintWriter writer;
    RandomAccessFile out;
    byte[] outStream;

    if (dosageValues == null && genotypeProbabilities != null) {
      if (bestGuessOutput) {
        computeDosageValueFromBestGenotype(bestGuessThreshold, -1, log);
      } else {
        computeDosageValues(log);
      }
    }
    if (dosageValues == null) {
      System.err.println("Error - Cannot write to Plink files without data!");
      return;
    }

    (new File(dir)).mkdirs();

    // FAM
    try {
      writer = Files.openAppropriateWriter(dir + plinkRoot + ".fam");
      for (String[] id : ids) {
        writer.print(id[0] + "\t" + id[1] + "\t");
        String[] parts = {id.length > 2 ? id[2] : "0", id.length > 3 ? id[3] : "0",
                          id.length > 4 ? id[4] : "0", id.length > 5 ? id[5] : "-9",};
        writer.println(ArrayUtils.toStr(parts, "\t"));
      }
      writer.close();
    } catch (IOException e) {
      log.reportError("Error writing '" + dir + plinkRoot + ".fam'");
      log.reportException(e);
      return;
    }

    try {
      // BIM
      String[] markerNames = markerSet.getMarkerNames();
      writer = Files.openAppropriateWriter(dir + plinkRoot + ".bim");
      for (int i = 0; i < markerNames.length; i++) {
        writer.println(chrs[i] + "\t" + markerNames[i] + "\t0\t" + positions[i] + "\t"
                       + alleles[i][0] + "\t" + alleles[i][1]);
      }
      writer.close();

      // BED
      out = new RandomAccessFile(dir + plinkRoot + ".bed", "rw");
      outStream = new byte[3];
      outStream[0] = (byte) 108; // 0b01101100
      outStream[1] = (byte) 27; // 0b00011011
      outStream[2] = (byte) 1; // 0b00000001 <-- be careful here
      out.write(outStream);

      for (int j = 0; j < markerNames.length; j++) {
        out.write(encodePlinkBedBytesForASingleMarkerOrSample(ArrayUtils.toByteArray(dosageValues[j])));
      }

      out.close();
    } catch (IOException e) {
      log.reportError("Error writing " + plinkRoot + ".bed and/or " + plinkRoot + ".bim in " + dir);
      log.reportException(e);
      return;
    }
  }

  private boolean[] filterMarkers(String[] markerNames, int[][] regionsToUse, String[] markersToUse,
                                  boolean verbose, Logger log) {
    boolean[] markersToKeep = ArrayUtils.booleanArray(markerNames.length, true);
    if (regionsToUse != null) {
      markerSet = markerSet.trim(regionsToUse, verbose, log);
      markersToKeep = ArrayUtils.booleanArray(markerNames.length, true);
      HashSet<String> newMarkerSet = HashVec.loadToHashSet(markerSet.getMarkerNames());
      for (int i = 0; i < markerNames.length; i++) {
        markersToKeep[i] = newMarkerSet.contains(markerNames[i]);
      }
      newMarkerSet = null;
    }
    if (markersToUse != null) {
      HashSet<String> mkrsToKeep = HashVec.loadToHashSet(markersToUse);
      HashSet<String> markers = HashVec.loadToHashSet(markerNames);
      HashSet<String> dupes = new HashSet<>();
      if (markerNames.length != markers.size()) {
        HashSet<String> hash = new HashSet<>((int) (1 + markerNames.length / .75));
        for (String mkr : markerNames) {
          if (hash.contains(mkr)) {
            dupes.add(mkr);
          } else {
            hash.add(mkr);
          }
        }
        hash = null;
      }
      markersToKeep = ArrayUtils.booleanArray(markerNames.length, false);
      for (int i = 0; i < markerNames.length; i++) {
        if (mkrsToKeep.contains(markerNames[i])) {
          markersToKeep[i] = true;
          if (!dupes.contains(markerNames[i])) {
            mkrsToKeep.remove(markerNames[i]);
          }
        }
      }
      if (mkrsToKeep.size() > 0 && verbose) {
        log.reportTimeWarning(mkrsToKeep.size()
                              + " markers listed in extract file not found in map file");
      }
      markerSet = markerSet.trim(markersToUse, true, verbose, log);
    }
    return markersToKeep;
  }

  /**
   * This is part of the method decodeBedByte(byte)
   *
   * @param bedByte
   * @return
   * @throws Elision
   */
  public static byte decodeLastTwoBitsOfABedByte(byte bedByte) throws Elision {
    byte genotype;

    bedByte = (byte) (bedByte & 0x03);

    if (bedByte == (byte) 0) {
      genotype = (byte) 0;

    } else if (bedByte == (byte) 2) {
      genotype = (byte) 1;

    } else if (bedByte == (byte) 3) {
      genotype = (byte) 2;

    } else if (bedByte == (byte) 1) {
      genotype = (byte) -1;

    } else {
      throw new Elision("Unrecognized genotype: " + bedByte);
    }

    return genotype;
  }

  /**
   * Convert a type from a PLINK .bed file to an array of genotypes with each element corresponding
   * to a single sample.
   *
   * @param bedByte
   * @return
   * @throws Elision
   */
  public static byte[] decodeBedByte(byte bedByte) throws Elision {
    byte[] genotypes;

    genotypes = new byte[4];
    for (int k = 0; k < 4; k++) {
      genotypes[k] = decodeLastTwoBitsOfABedByte((byte) (bedByte >> (2 * k)));
    }

    return genotypes;
  }

  /**
   * This is part of the method encodePlinkBedBytesForASingleMarkOrSamp(byte[] )
   *
   * @param genotype
   * @return
   * @throws Elision
   */
  public static byte encodeLastTwoBitsOfABedByte(byte genotype) throws Elision {
    byte bedByte;

    if (genotype == (byte) 0) {
      bedByte = (byte) 0x00;

    } else if (genotype == (byte) 1) {
      bedByte = (byte) 0x02;

    } else if (genotype == (byte) 2) {
      bedByte = (byte) 0x03;

    } else if (genotype == (byte) -1) {
      bedByte = (byte) 0x01;

    } else {
      throw new Elision("Unrecognized genotype: " + genotype
                        + ". Please use 0 for A/A, 1 for A/B, 2 for B/B, and -1 for null.");
    }

    return bedByte;
  }

  /**
   * Convert the array of genotypes into an array of PLINK .bed byte stream.
   *
   * @param genotype
   * @return
   */
  public static byte[] encodePlinkBedBytesForASingleMarkerOrSample(byte[] genotype) {
    int iBytes;
    byte[] result;
    int nBytes;
    byte shift;

    nBytes = (int) Math.ceil((double) genotype.length / 4);
    iBytes = -1;
    result = new byte[nBytes];

    // for (int i = 0; i < result.length; i++) {
    // result[i] = (byte) 0xAA; //initilize the array to be 0b10101010, the null genotype defined by
    // PLINK bed data.
    // }

    try {
      for (int i = 0; i < genotype.length; i++) {
        shift = (byte) ((i % 4) * 2);
        if (shift == 0) {
          iBytes++;
        }
        result[iBytes] = (byte) ((result[iBytes] & (~(0x03 << shift)))
                                 | (encodeLastTwoBitsOfABedByte(genotype[i]) << shift));
        // displayBits(result[iBytes]);
      }
    } catch (Elision e) {
      e.printStackTrace();
    }

    return result;
  }

  public static enum COMBINE_OP {
    // FAIL, DROP, OVERWRITE_IF_ALL_MISSING;
    FAIL, DROP, EITHER_IF_OTHER_MISSING, EITHER_IF_OTHER_MISSING_ALL;

  }

  public static DosageData combine(DosageData dd1, DosageData dd2, COMBINE_OP onDupeOp,
                                   boolean computeDosageAsBestGuess, float bestGuessMissThresh,
                                   Logger log) {
    if (!dd1.isEmpty() && !dd2.isEmpty()) {
      // let through
    } else if (dd1.isEmpty() && !dd2.isEmpty()) {
      log.reportError("Warning - DosageData {1} provided to combine() was empty.");
      return dd2;
    } else if (!dd1.isEmpty() && dd2.isEmpty()) {
      log.reportError("Warning - DosageData {2} provided to combine() was empty.");
      return dd1;
    } else if (dd1.isEmpty() && dd2.isEmpty()) {
      log.reportError("Warning - both DosageData objects provided to combine() were empty.");
      return dd1;
    }

    int dd1NumGeno = dd1.genotypeProbabilities == null ? (dd1.dosageValues == null ? 0 : 1)
                                                       : dd1.genotypeProbabilities[0][0].length;
    int dd2NumGeno = dd2.genotypeProbabilities == null ? (dd2.dosageValues == null ? 0 : 1)
                                                       : dd2.genotypeProbabilities[0][0].length;

    boolean dosageOverride = false;
    if ((dd1NumGeno > 1 || dd2NumGeno > 1) && dd1NumGeno != dd2NumGeno) {
      log.reportError("Warning - cannot combine different numbers of genotype probabilities - result data will be imputed dosages.");
      dosageOverride = true;
    }
    int ddNewNumGeno = dosageOverride ? 1 : Math.min(dd1NumGeno, dd2NumGeno);
    if (ddNewNumGeno == 0) {
      log.reportError("Error - cannot combine data sets when a dataset is missing both genotype and dosage data [dataset "
                      + (dd1NumGeno == 0 ? "1" : "2") + "]");
      System.exit(1);
    }

    byte missingChr = 0;
    int missingPos = 0;
    String[] missingAlleles = null;
    float missingDosage = Float.NaN;
    float missingGeno = Float.NaN;

    String[][] dd1Ids = dd1.ids;
    String[][] dd2Ids = dd2.ids;
    HashMap<String, Integer> dd1IdsAndIndices = new HashMap<>();
    HashMap<String, Integer> dd2IdsAndIndices = new HashMap<>();
    HashSet<String> duplicatedIDs = new HashSet<>();
    Map<String, String[]> idSet = new LinkedHashMap<>(); // use to ensure uniqueness and order
    for (int s = 0; s < dd1Ids.length; s++) {
      String id = dd1Ids[s][0] + "\t" + dd1Ids[s][1];
      idSet.put(id, dd1Ids[s]);
      dd1IdsAndIndices.put(id, s);
    }
    for (int s = 0; s < dd2Ids.length; s++) {
      String id = dd2Ids[s][0] + "\t" + dd2Ids[s][1];
      dd2IdsAndIndices.put(id, s);
      String[] otherIds = idSet.get(id);
      if (otherIds == null)
        idSet.put(id, dd2Ids[s]);
      else {
        duplicatedIDs.add(id);
        if (!Arrays.equals(dd2Ids[s], otherIds)) {
          boolean mismatchData = false;
          for (int i = 0; i < Math.min(dd2Ids[s].length, otherIds.length); i++) {
            if (!dd2Ids[s][i].equals(otherIds[i])) {
              mismatchData = true;
              break;
            }
          }
          if (mismatchData) {
            log.reportError("For sample " + id + " fam data to combine does not match");
          }
          if (dd2Ids[s].length > otherIds.length) idSet.put(id, dd2Ids[s]);
        }
      }
    }
    if (duplicatedIDs.size() > 0) {
      log.report(duplicatedIDs.size() + " duplicate sample IDs found, out of " + dd2Ids.length
                 + " samples present.");
    }

    String[] dd1Mkrs = dd1.markerSet.getMarkerNames();
    String[] dd2Mkrs = dd2.markerSet.getMarkerNames();
    LinkedHashSet<String> markers = new LinkedHashSet<>();
    HashMap<String, Integer> dd1MarkersAndIndices = new HashMap<>();
    for (int i = 0; i < dd1Mkrs.length; i++) {
      dd1MarkersAndIndices.put(dd1Mkrs[i], i);
      markers.add(dd1Mkrs[i]);
    }
    HashMap<String, Integer> duplicatedMarkersAndIndices = new HashMap<>();
    HashMap<String, Integer> dd2MarkersAndIndices = new HashMap<>();
    HashSet<String> droppedMarkers = new HashSet<>();
    for (int i = 0; i < dd2Mkrs.length; i++) {
      dd2MarkersAndIndices.put(dd2Mkrs[i], i);

      boolean alreadyPresentMkr = !markers.add(dd2Mkrs[i]);
      if (alreadyPresentMkr) {
        log.reportTime("Duplicate marker: " + dd2Mkrs[i]);

        duplicatedMarkersAndIndices.put(dd2Mkrs[i], i);
        if (duplicatedIDs.size() > 0) {

          if (onDupeOp == COMBINE_OP.FAIL) {
            log.reportError("Error - cannot combine data sets with duplicated marker AND sample names.  Yet.");
            // TODO combining values when marker and sample are the same?
            System.exit(1);
          } else if (onDupeOp == COMBINE_OP.DROP) {
            log.reportTimeWarning("cannot combine data sets with duplicated marker AND sample names.  Marker "
                                  + dd2Mkrs[i] + " will be dropped");
            droppedMarkers.add(dd2Mkrs[i]);
          }
        }
      }
    }

    HashSet<String> disallowedDuplicateIDsForALL = new HashSet<>();
    if (duplicatedIDs.size() > 0 && duplicatedMarkersAndIndices.size() > 0) {
      if (onDupeOp == COMBINE_OP.EITHER_IF_OTHER_MISSING_ALL) {
        // find any duplicated/shared IDs that have ANY valid data for duplicated/shared markers
        for (String s : duplicatedIDs) {
          int idInd1 = dd1IdsAndIndices.get(s);
          int idInd2 = dd2IdsAndIndices.get(s);

          boolean missingAll1 = true, missingAll2 = true;
          for (String m : duplicatedMarkersAndIndices.keySet()) {
            int mkrInd1 = dd1MarkersAndIndices.get(m);
            boolean valid1 = (dd1.genotypeProbabilities != null
                              && !checkMissing(dd1.genotypeProbabilities[mkrInd1][idInd1],
                                               missingGeno))
                             || (dd1.dosageValues != null
                                 && !checkMissing(dd1.dosageValues[mkrInd1][idInd1],
                                                  missingDosage));
            if (valid1) {
              missingAll1 = false;
              break;
            }
          }
          if (missingAll1) {
            continue; // don't care about dd2 if it's known to be already missing all from dd1
          }
          for (String m : duplicatedMarkersAndIndices.keySet()) {
            int mkrInd2 = dd2MarkersAndIndices.get(m);
            boolean valid2 = missingAll2
                             && (dd2.genotypeProbabilities != null
                                 && !checkMissing(dd2.genotypeProbabilities[mkrInd2][idInd2],
                                                  missingGeno))
                             || (dd2.dosageValues != null
                                 && !checkMissing(dd2.dosageValues[mkrInd2][idInd2],
                                                  missingDosage));

            if (valid2) {
              missingAll2 = false;
              break;
            }
          }
          if (!missingAll1 && !missingAll2) {
            disallowedDuplicateIDsForALL.add(s);
          }
        }
      }
    }

    DosageData ddNew = new DosageData();
    ddNew.ids = idSet.values().toArray(new String[0][]);
    idSet = null; // can now refer to ddNew.ids

    // don't use merge, as it sorts markers after merging
    // ddNew.markerSet = SnpMarkerSet.merge(dd1.markerSet, dd2.markerSet);
    byte[] chrSrc, chrSrc2;
    String[][] alleleSrc, alleleSrc2;
    int[] posSrc, posSrc2;

    for (String s : droppedMarkers) {
      markers.remove(s);
    }

    ddNew.alleles = new String[markers.size()][];
    ddNew.chrs = new byte[markers.size()];
    ddNew.positions = new int[markers.size()];
    int dd1NumAnnot = dd1.markerSet.getAnnotation() == null
                      || dd1.markerSet.getAnnotation().length == 0 ? 0
                                                                   : dd1.markerSet.getAnnotation()[0].length;
    int dd2NumAnnot = dd2.markerSet.getAnnotation() == null
                      || dd2.markerSet.getAnnotation().length == 0 ? 0
                                                                   : dd2.markerSet.getAnnotation()[0].length;
    String[][] annotations = new String[markers.size()][dd1NumAnnot + dd2NumAnnot];
    String[] annotLbls = new String[dd1NumAnnot + dd2NumAnnot];
    for (int i = 0; i < dd1NumAnnot; i++) {
      annotLbls[i] = dd1.labelPrepend;
    }
    for (int i = 0; i < dd2NumAnnot; i++) {
      annotLbls[i + dd1NumAnnot] = dd2.labelPrepend;
    }

    chrSrc = dd1.chrs == null ? dd1.markerSet.getChrs() : dd1.chrs;
    alleleSrc = dd1.alleles == null ? dd1.markerSet.getAlleles() : dd1.alleles;
    posSrc = dd1.positions == null ? dd1.markerSet.getPositions() : dd1.positions;

    chrSrc2 = dd2.chrs == null ? dd2.markerSet.getChrs() : dd2.chrs;
    alleleSrc2 = dd2.alleles == null ? dd2.markerSet.getAlleles() : dd2.alleles;
    posSrc2 = dd2.positions == null ? dd2.markerSet.getPositions() : dd2.positions;

    String[] mkrs = new String[markers.size()];
    Iterator<String> markerIter = markers.iterator();
    int m = 0;
    while (markerIter.hasNext()) {
      String mkr = markerIter.next();

      mkrs[m] = mkr;
      if (dd1MarkersAndIndices.containsKey(mkr)) {
        int ind = dd1MarkersAndIndices.get(mkr);
        ddNew.chrs[m] = chrSrc == null ? missingChr : chrSrc[ind];
        ddNew.alleles[m] = alleleSrc == null ? missingAlleles : alleleSrc[ind];
        ddNew.positions[m] = posSrc == null ? missingPos : posSrc[ind];
        String[] annot = dd1.markerSet.getAnnotation() == null ? new String[0]
                                                               : dd1.markerSet.getAnnotation()[ind];
        for (int i = 0; i < dd1NumAnnot; i++) {
          annotations[m][i] = annot[i];
        }
      } else if (dd2MarkersAndIndices.containsKey(mkr)) {
        int ind = dd2MarkersAndIndices.get(mkr);
        ddNew.chrs[m] = chrSrc2 == null ? missingChr : chrSrc2[ind];
        ddNew.alleles[m] = alleleSrc2 == null ? missingAlleles : alleleSrc2[ind];
        ddNew.positions[m] = posSrc2 == null ? missingPos : posSrc2[ind];
        String[] annot = dd2.markerSet.getAnnotation() == null ? new String[0]
                                                               : dd2.markerSet.getAnnotation()[ind];
        for (int i = 0; i < dd2NumAnnot; i++) {
          annotations[m][i + dd1NumAnnot] = annot[i];
        }
      } else {
        ddNew.chrs[m] = missingChr;
        ddNew.alleles[m] = missingAlleles;
        ddNew.positions[m] = missingPos;
      }
      m++;
    }

    ddNew.markerSet = new SnpMarkerSet(mkrs, ddNew.chrs, ddNew.positions, ddNew.alleles,
                                       annotations, false, false);

    if ((dosageOverride || dd1NumGeno == 1) && dd2NumGeno > 1 && dd2.dosageValues == null) {
      if (computeDosageAsBestGuess) {
        dd2.computeDosageValueFromBestGenotype(bestGuessMissThresh, Float.NaN, log);
      } else {
        dd2.computeDosageValues(log);
      }
    }
    if (dd1NumGeno > 1 && (dosageOverride || dd2NumGeno == 1) && dd1.dosageValues == null) {
      if (computeDosageAsBestGuess) {
        dd1.computeDosageValueFromBestGenotype(bestGuessMissThresh, Float.NaN, log);
      } else {
        dd1.computeDosageValues(log);
      }
    }

    ddNew.genotypeProbabilities = ddNewNumGeno > 1 ? new float[markers.size()][ddNew.ids.length][ddNewNumGeno]
                                                   : null;
    ddNew.dosageValues = ddNewNumGeno == 1 ? new float[markers.size()][ddNew.ids.length] : null;

    if (ddNewNumGeno > 1) {
      // combine genotypeProbs
      markerIter = markers.iterator();
      m = 0;
      while (markerIter.hasNext()) {
        String mkr = markerIter.next();

        for (int s = 0; s < ddNew.ids.length; s++) {
          String id = ddNew.ids[s][0] + "\t" + ddNew.ids[s][1];

          boolean idIn1, mkrIn1, idIn2, mkrIn2;
          mkrIn1 = dd1MarkersAndIndices.containsKey(mkr);
          idIn1 = dd1IdsAndIndices.containsKey(id);
          mkrIn2 = dd2MarkersAndIndices.containsKey(mkr);
          idIn2 = dd2IdsAndIndices.containsKey(id);

          float[] dd1Data = (mkrIn1
                             && idIn1) ? dd1.genotypeProbabilities[dd1MarkersAndIndices.get(mkr)][dd1IdsAndIndices.get(id)]
                                       : ArrayUtils.floatArray(ddNewNumGeno, missingGeno);
          float[] dd2Data = (mkrIn2
                             && idIn2) ? dd2.genotypeProbabilities[dd2MarkersAndIndices.get(mkr)][dd2IdsAndIndices.get(id)]
                                       : ArrayUtils.floatArray(ddNewNumGeno, missingGeno);

          if (dd1Data != null && dd2Data == null) {
            ddNew.genotypeProbabilities[m][s] = dd1Data;
          } else if (dd1Data == null && dd2Data != null) {
            ddNew.genotypeProbabilities[m][s] = dd2Data;
          } else if (dd1Data != null && dd2Data != null) {
            // "data" may be present, but also may be set to missing - check and respond
            // appropriately
            boolean miss1, miss2;
            miss1 = checkMissing(dd1Data, missingGeno);
            miss2 = checkMissing(dd2Data, missingGeno);

            if (miss1 && !miss2) {
              // valid data in dd1 and missing data in dd2
              if (onDupeOp == COMBINE_OP.EITHER_IF_OTHER_MISSING) {
                ddNew.genotypeProbabilities[m][s] = dd2Data;
              } else if (onDupeOp == COMBINE_OP.EITHER_IF_OTHER_MISSING_ALL) {
                // check that there is no data in dd1 for indiv for any duplicated/shared markers
                if (!disallowedDuplicateIDsForALL.contains(id)) {
                  ddNew.genotypeProbabilities[m][s] = dd2Data;
                } else {
                  ddNew.genotypeProbabilities[m][s] = dd1Data; // set to missingGeno
                }
              } else {
                // should have died before this point if onDupeOp is FAIL, and wouldn't reach this
                // check if onDupeOp is DROP
              }
            } else if (!miss1 && miss2) {
              // valid data in dd1 and missing data in dd2
              if (onDupeOp == COMBINE_OP.EITHER_IF_OTHER_MISSING) {
                ddNew.genotypeProbabilities[m][s] = dd1Data;
              } else if (onDupeOp == COMBINE_OP.EITHER_IF_OTHER_MISSING_ALL) {
                // check that there is no data in dd2 for indiv for any duplicated/shared markers
                if (!disallowedDuplicateIDsForALL.contains(id)) {
                  ddNew.genotypeProbabilities[m][s] = dd1Data;
                } else {
                  ddNew.genotypeProbabilities[m][s] = dd2Data; // set to missingGeno
                }
              } else {
                // should have died before this point if onDupeOp is FAIL, and wouldn't reach this
                // check if onDupeOp is DROP
              }
            } else if (miss1 && miss2) {
              // easy, set to missing, though why would this occur? Weird...
              ddNew.genotypeProbabilities[m][s] = dd1Data; // set to missingGeno
            } else if (!miss1 && !miss2) {
              log.reportError("Error - valid genotype data in both DosageData objects for ID/MKR: \""
                              + id + "\" / \"" + mkr + "\" - data will be set to missing.");
              ddNew.genotypeProbabilities[m][s] = ArrayUtils.floatArray(ddNewNumGeno, missingGeno);
            }
          } else {
            ddNew.genotypeProbabilities[m][s] = ArrayUtils.floatArray(ddNewNumGeno, missingGeno);
          }
        }

        m++;
      }
    } else if (ddNewNumGeno == 1) {

      markerIter = markers.iterator();
      m = 0;
      while (markerIter.hasNext()) {
        String mkr = markerIter.next();

        for (int s = 0; s < ddNew.ids.length; s++) {
          String id = ddNew.ids[s][0] + "\t" + ddNew.ids[s][1];

          boolean idIn1, mkrIn1, idIn2, mkrIn2;
          mkrIn1 = dd1MarkersAndIndices.containsKey(mkr) && dd1MarkersAndIndices.get(mkr) != null;
          idIn1 = dd1IdsAndIndices.containsKey(id) && dd1IdsAndIndices.get(id) != null;
          mkrIn2 = dd2MarkersAndIndices.containsKey(mkr) && dd2MarkersAndIndices.get(mkr) != null;
          idIn2 = dd2IdsAndIndices.containsKey(id) && dd2IdsAndIndices.get(id) != null;

          float dd1Data = (mkrIn1
                           && idIn1) ? dd1.dosageValues[dd1MarkersAndIndices.get(mkr)][dd1IdsAndIndices.get(id)]
                                     : missingDosage;
          float dd2Data = (mkrIn2
                           && idIn2) ? dd2.dosageValues[dd2MarkersAndIndices.get(mkr)][dd2IdsAndIndices.get(id)]
                                     : missingDosage;

          // "data" may be present, but also may be set to missing - check and respond appropriately
          boolean miss1, miss2;
          miss1 = checkMissing(dd1Data, missingDosage);
          miss2 = checkMissing(dd2Data, missingDosage);

          if (miss1 && !miss2) {
            // valid data in dd1 and missing data in dd2
            if (onDupeOp == COMBINE_OP.EITHER_IF_OTHER_MISSING) {
              ddNew.dosageValues[m][s] = dd2Data;
            } else if (onDupeOp == COMBINE_OP.EITHER_IF_OTHER_MISSING_ALL) {
              // check that there is no data in dd1 for indiv for any duplicated/shared markers
              if (!disallowedDuplicateIDsForALL.contains(id)) {
                ddNew.dosageValues[m][s] = dd2Data;
              } else {
                ddNew.dosageValues[m][s] = dd1Data; // set to missingGeno
              }
            } else {
              // should have died before this point if onDupeOp is FAIL, and wouldn't reach this
              // check if onDupeOp is DROP
            }
          } else if (!miss1 && miss2) {
            // valid data in dd1 and missing data in dd2
            if (onDupeOp == COMBINE_OP.EITHER_IF_OTHER_MISSING) {
              ddNew.dosageValues[m][s] = dd1Data;
            } else if (onDupeOp == COMBINE_OP.EITHER_IF_OTHER_MISSING_ALL) {
              // check that there is no data in dd2 for indiv for any duplicated/shared markers
              if (!disallowedDuplicateIDsForALL.contains(id)) {
                ddNew.dosageValues[m][s] = dd1Data;
              } else {
                ddNew.dosageValues[m][s] = dd2Data; // set to missingGeno
              }
            } else {
              // should have died before this point if onDupeOp is FAIL, and wouldn't reach this
              // check if onDupeOp is DROP
            }
          } else if (miss1 && miss2) {
            // easy, set to missing, though why would this occur? Weird...
            ddNew.dosageValues[m][s] = missingDosage; // set to missingGeno
          } else if (!miss1 && !miss2) {
            log.reportError("Error - valid dosage data in both DosageData objects for ID/MKR: \""
                            + id + "\" / \"" + mkr + "\" - data will be set to missing.");
            ddNew.dosageValues[m][s] = missingDosage;
          }

        }

        m++;
      }
    }

    int[] keys = Sort.getSort2DIndices(ddNew.chrs, ddNew.positions);
    ddNew.chrs = Sort.getOrdered(ddNew.chrs, keys);
    ddNew.positions = Sort.getOrdered(ddNew.positions, keys);
    if (ddNew.alleles != null) {
      ddNew.alleles = Sort.getOrdered(ddNew.alleles, keys);
    }
    ddNew.markerSet.sortMarkers();
    if (ddNew.genotypeProbabilities != null) {
      ddNew.genotypeProbabilities = Sort.getOrdered(ddNew.genotypeProbabilities, keys);
    }
    if (ddNew.dosageValues != null) {
      ddNew.dosageValues = Sort.getOrdered(ddNew.dosageValues, keys);
    }

    return ddNew;
  }

  private static boolean checkMissing(float f, float missingDose) {
    boolean missNaN = Float.isNaN(missingDose);
    if (missNaN && !Float.isNaN(f)) return false;
    if (Float.compare(f, missingDose) != 0) return false;
    return true;
  }

  private static boolean checkMissing(float[] fs, float missingGeno) {
    boolean missNaN = Float.isNaN(missingGeno);
    for (float f : fs) {
      if (missNaN && !Float.isNaN(f)) return false;
      if (Float.compare(f, missingGeno) != 0) return false;
    }
    return true;
  }

  public static void convert(String dosageFile, String idFile, String mapFile, int fromFormat,
                             String outfile, String mapOut, String extract, int toFormat,
                             boolean awk, boolean verbose, Logger log) {
    convert(dosageFile, idFile, mapFile, PARAMETERS[fromFormat], outfile, mapOut, extract,
            PARAMETERS[toFormat], awk, verbose, log);
  }

  public static void convert(String dosageFile, String idFile, String mapFile, int[] fromParameters,
                             String outfile, String mapOut, String extract, int[] toParameters,
                             boolean awk, boolean verbose, Logger log) {
    SnpMarkerSet markerSet;
    String[][] ids;
    float dosageValue;
    float[] genotypeProbabilities;
    byte[] chrs;
    int[] positions;
    String[][] alleles;
    byte chr;
    int position;
    String[] allelePair;

    BufferedReader reader;
    String[] line, markerNames;
    Hashtable<String, String> invalids;

    PrintWriter writer;
    String delimiter;
    String[] lead, markersToKeep;
    HashSet<String> keeps;
    String root;
    IntVector cols;
    int offset, numColsEach;
    boolean dos;
    String ranges;
    int start, stop;

    if (fromParameters[2] != toParameters[2]) {
      log.reportError("Conversion will have to take place in memory in order to transpose between marker dominant/individual dominant; use new DosageDate().writeToFile() instead");
      return;
    }

    markerSet = new SnpMarkerSet(mapFile, false, new Logger());
    markerNames = markerSet.getMarkerNames();
    chrs = markerSet.getChrs();
    positions = markerSet.getPositions();
    alleles = markerSet.getAlleles();
    ids = HashVec.loadFileToStringMatrix(idFile, false, new int[] {0, 1});

    if (extract == null) {
      keeps = null;
    } else {
      markersToKeep = HashVec.loadFileToStringArray(extract, false, new int[] {0}, true, false,
                                                    "\t");
      keeps = HashVec.loadToHashSet(markersToKeep);
      root = ext.rootOf(outfile, false);
      markerSet = markerSet.trim(markersToKeep, true, false, log); // allows missing markers, but
                                                                   // will list how many
      if (mapOut == null) {
        mapOut = root + ".pinfo";
        markerSet.writeToFile(mapOut, SnpMarkerSet.determineType(mapOut), log);
        mapOut = root + ".map";
      }
      markerSet.writeToFile(mapOut, SnpMarkerSet.determineType(mapOut), log);
      try {
        writer = Files.openAppropriateWriter(root + ".ids.fam");
        for (String[] id : ids) {
          writer.println(id[0] + "\t" + id[1]);
        }
        writer.close();
      } catch (Exception e) {
        System.err.println("Error writing to " + root + ".ids.fam");
        e.printStackTrace();
      }
    }

    // /** 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 */
    // /** id_type, column index where dosage/probability values begin, dominance format, number of
    // columns summarizing the data, header row, marker/IID index, A1 index, A2 index, chr index,
    // pos index, delimiter, index for head/lead, min number of digits, max number of digits */
    //
    // /** 0 1 2 3 4 5 6 7 8 9 10 11 12 13 */
    // public static final int[][] PARAMETERS = {{MACH_ID_TYPE, 2, INDIVIDUAL_DOMINANT_FORMAT, 1, 0,
    // 0, -1, -1, -1, -1, 0, 0, 3, 3}, // .mldose (MACH)
    // {SEPARATE_FILE_ID_TYPE, 5, MARKER_DOMINANT_FORMAT, 3, 0, 1, 3, 4, 0, 2, 2, 1, 0, 3}, // .gen
    // {IID_TYPE, 1, INDIVIDUAL_DOMINANT_FORMAT, 1, 1, 0, -1, -1, -1, -1, 1, 2, 0, 3}, // .fhsR
    // (GWAF)
    // {FID_IID_TYPE, 3, MARKER_DOMINANT_FORMAT, 2, 1, 0, 1, 2, -1, -1, 0, 3, 2, 2}, // .dosage
    // (PLINK)
    // {MACH_ID_TYPE, 2, INDIVIDUAL_DOMINANT_FORMAT, 2, 0, 0, -1, -1, -1, -1, 0, 4, 3, 3}, //
    // .mlprob (MACH)
    // {MACH_ID_TYPE, 2, INDIVIDUAL_DOMINANT_FORMAT, 1, 0, 0, -1, -1, -1, -1, 0, 5, 3, 3}, // .dose
    // (MINIMAC)
    // };

    if (awk && !ArrayUtils.equals(fromParameters, toParameters)) {
      log.reportError("Error - the awk option of convert is currently only available for identical file types");
    } else if (awk) {
      cols = new IntVector();
      offset = fromParameters[1];
      for (int i = 2; i <= offset; i++) {
        cols.add(i);
      }
      numColsEach = fromParameters[3];
      for (int i = 0; i < markerNames.length; i++) {
        if (keeps.contains(markerNames[i])) {
          for (int j = 0; j < numColsEach; j++) {
            cols.add(offset + i * numColsEach + j + 1);
          }
        }
      }
      ranges = ext.listRanges(Ints.toArray(cols));
      System.out.println("Extracting columns: 1," + ranges);

      dos = Files.isWindows();
      try {
        writer = Files.openAppropriateWriter("awk_command.bat");
        writer.print("awk " + (dos ? "\"" : "'") + "BEGIN { OFS = " + (dos ? "\\" : "") + "\"\\t"
                     + (dos ? "\\" : "") + "\" } { printf $1 ;");
        line = ranges.split(",");
        for (String element : line) {
          if (element.contains("-")) {
            start = Integer.parseInt(element.substring(0, element.indexOf("-")));
            stop = Integer.parseInt(element.substring(element.indexOf("-") + 1));
          } else {
            start = stop = Integer.parseInt(element);
          }
          writer.print(" for (i=" + start + "; i<=" + stop + "; i++) printf " + (dos ? "\\" : "")
                       + "\"\\t" + (dos ? "\\" : "") + "\"$i ;");
        }
        writer.println(" printf " + (dos ? "\\" : "") + "\"\\n" + (dos ? "\\" : "") + "\" }"
                       + (dos ? "\"" : "'") + " " + dosageFile + " > " + outfile);
        writer.close();
      } catch (Exception e) {
        System.err.println("Error writing to " + "awk_command.bat");
        e.printStackTrace();
      }

      // System.out.println("Running the following command to parse the new dosage
      // file:\n"+command);
      Files.chmod("awk_command.bat", false);
      CmdLine.run((dos ? "" : "./") + "awk_command.bat", "./", System.err);
    } else {
      allelePair = new String[2];
      if (fromParameters[3] > 1) {
        genotypeProbabilities = new float[fromParameters[3]];
      } else {
        genotypeProbabilities = null;
      }
      dosageValue = -999;

      delimiter = DELIMITERS[toParameters[10]];
      invalids = new Hashtable<>();
      try {
        reader = Files.getAppropriateReader(dosageFile); // new BufferedReader(new
                                                         // FileReader(dosageFile));
        writer = Files.openAppropriateWriter(outfile);

        if (fromParameters[4] == 1) {
          line = reader.readLine().trim()
                       .split(fromParameters[10] == 1 ? "," : PSF.Regex.GREEDY_WHITESPACE);
          if (fromParameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
            for (int i = 0; i < markerNames.length; i++) {
              if (!markerNames[i].equals(line[fromParameters[1] + fromParameters[3] * i])) {
                log.reportError("Error - mismatched name at marker " + (i + 1) + " of " + dosageFile
                                + "; expecting " + markerNames[i] + " given map file " + mapFile
                                + ", found " + line[fromParameters[1] + fromParameters[3] * i]);
                reader.close();
                return;
              }
            }
          } else if (fromParameters[2] == MARKER_DOMINANT_FORMAT) {
            if (fromParameters[3] != 2) {
              log.reportError("Warning - ignoring the header with IDs in file " + dosageFile
                              + " because it does not contain 2 columns for each individual");
            } else {
              for (int i = 0; i < ids.length; i++) {
                if (!ids[i][0].equals(line[fromParameters[1] + fromParameters[3] * i + 0])
                    || !ids[i][1].equals(line[fromParameters[1] + fromParameters[3] * i + 1])) {
                  log.reportError("Error - mismatched IDs at individual " + (i + 1) + " of "
                                  + dosageFile + "; expecting " + ids[i][0] + "," + ids[i][1]
                                  + " given id file " + idFile + ", found "
                                  + line[fromParameters[1] + fromParameters[3] * i + 0] + ","
                                  + line[fromParameters[1] + fromParameters[3] * i + 1]);
                  reader.close();
                  return;
                }
              }
            }
          }
        }

        if (toParameters[4] == 1) { // if there is a header row
          writer.print(ArrayUtils.toStr(HEADS[toParameters[11]], delimiter));

          if (toParameters[2] == MARKER_DOMINANT_FORMAT) {
            if (toParameters[3] == 2 && toParameters[0] == FID_IID_TYPE) {
              for (String[] id : ids) {
                writer.print(delimiter + id[0] + delimiter + id[1]);
              }
            } else if (toParameters[3] == 1 && toParameters[0] == IID_TYPE) {
              for (String[] id : ids) {
                writer.print(delimiter + id[1]);
              }
            } else {
              log.reportError("Error - don't know how to list IDs when there "
                              + (toParameters[3] == 1 ? "is one column"
                                                      : "are " + toParameters[3] + " columns")
                              + " for dosage inforation and the ID type is '" + toParameters[2]
                              + "'");
              System.exit(1);
            }
          } else if (toParameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
            for (String markerName : markerNames) {
              if (extract == null || keeps.contains(markerName)) {
                for (int j = 0; j < toParameters[3]; j++) {
                  writer.print(delimiter + markerName);
                }
              }
            }
          }
          writer.println();
        }

        if (fromParameters[2] == MARKER_DOMINANT_FORMAT) {
          for (int i = 0; i < markerNames.length; i++) {
            line = reader.readLine().trim()
                         .split(fromParameters[10] == 1 ? "," : PSF.Regex.GREEDY_WHITESPACE);
            if (!markerNames[i].equals(line[fromParameters[5]])) {
              log.reportError("Error - mismatched name at marker " + (i + 1) + " of " + dosageFile
                              + "; expecting " + markerNames[i] + " given map file " + mapFile
                              + ", found " + line[fromParameters[5]]);
              reader.close();
              return;
            }

            // gathering alleles if available
            if (fromParameters[6] != -1) {
              if (line[fromParameters[6]].length() > 1
                  || !Sequence.validBase(line[fromParameters[6]])) {
                log.reportError("Warning - invalid allele ('" + line[fromParameters[6]]
                                + "') at marker " + markerNames[i]);
              }
              allelePair[0] = line[fromParameters[6]];
            } else if (alleles != null) {
              allelePair[0] = alleles[i][0];
            } else {
              allelePair[0] = "~";
            }

            // gathering alleles if available
            if (fromParameters[7] != -1) {
              if (line[fromParameters[7]].length() > 1
                  || !Sequence.validBase(line[fromParameters[7]])) {
                log.reportError("Warning - invalid allele ('" + line[fromParameters[7]]
                                + "') at marker " + markerNames[i]);
              }
              allelePair[1] = line[fromParameters[7]];
            } else if (alleles != null) {
              allelePair[1] = alleles[i][1];
            } else {
              allelePair[1] = "~";
            }

            // gathering chromosome if available
            if (fromParameters[8] != -1) {
              try {
                chr = Byte.parseByte(line[fromParameters[8]]);
              } catch (NumberFormatException nfe) {
                chr = -1;
                if (!invalids.containsKey(line[fromParameters[8]])) {
                  log.reportError("Warning - invalid chromosome number ('" + line[fromParameters[8]]
                                  + "'), first seen at marker " + markerNames[i]);
                  invalids.put(line[fromParameters[8]], "");
                }
              }
            } else if (chrs != null) {
              chr = chrs[i];
            } else {
              chr = -9;
            }

            // gathering position if available
            if (fromParameters[9] != -1) {
              try {
                position = Integer.parseInt(line[fromParameters[9]]);
              } catch (NumberFormatException nfe) {
                position = -1;
                if (!invalids.containsKey(line[fromParameters[9]])) {
                  log.reportError("Warning - invalid genome position ('" + line[fromParameters[9]]
                                  + "') for marker " + markerNames[i]);
                  invalids.put(line[fromParameters[9]], "");
                }
              }
            } else if (positions != null) {
              position = positions[i];
            } else {
              position = 0;
            }

            if (line.length - fromParameters[1] != ids.length * fromParameters[3]) {
              log.reportError("Error - mismatched number of elements in line "
                              + (i + 1 + fromParameters[4]) + " of " + dosageFile + "; expecting "
                              + ids.length + "*" + fromParameters[3] + "+" + fromParameters[1]
                              + ", found " + line.length);
              System.exit(1);
            }

            if (extract == null || keeps.contains(markerNames[i])) {
              lead = LEADS[toParameters[11]] == null ? new String[toParameters[1]]
                                                     : LEADS[toParameters[11]];
              lead[toParameters[5]] = markerNames[i];

              // adding alleles if required
              if (toParameters[6] >= 0) {
                if (allelePair[0] == "~" || allelePair[1] == "~") {
                  log.reportError("Error - file format requires alleles and none were supplied via the map file or the source dosage file");
                  System.exit(1);
                }
                lead[toParameters[6]] = allelePair[0] + "";
              }
              // adding alleles if required
              if (toParameters[7] >= 0) {
                lead[toParameters[7]] = allelePair[1] + "";
              }
              // adding chromosome if required
              if (toParameters[8] >= 0) {
                if (chr == -9) {
                  log.reportError("Error - file format requires chrs/positions and none were supplied via the map file or the source dosage file");
                  System.exit(1);
                }
                lead[toParameters[8]] = chr == -1 ? "-" : chr + "";
              }
              // adding position if required
              if (toParameters[9] >= 0) {
                lead[toParameters[9]] = position + "";
              }
              writer.print(ArrayUtils.toStr(lead, delimiter));

              for (int j = 0; j < ids.length; j++) {
                if (fromParameters[3] == 1) {
                  dosageValue = Float.parseFloat(line[fromParameters[1] + j]);
                } else {
                  genotypeProbabilities[0] = Float.parseFloat(line[fromParameters[1]
                                                                   + j * fromParameters[3] + 0]);
                  genotypeProbabilities[1] = Float.parseFloat(line[fromParameters[1]
                                                                   + j * fromParameters[3] + 1]);
                  if (fromParameters[3] == 3) {
                    genotypeProbabilities[2] = Float.parseFloat(line[fromParameters[1]
                                                                     + j * fromParameters[3] + 2]);
                    if (Math.abs((1 - genotypeProbabilities[1] - genotypeProbabilities[0])
                                 - genotypeProbabilities[2]) > 0.01) {
                      log.reportError("Error: P(BB) does not equal [ 1 - P(AA) - P(AB) ] for individual "
                                      + ids[j][0] + "," + ids[j][1] + " at marker " + markerNames[i]
                                      + " which is line " + (i + 1 + fromParameters[4]) + " of "
                                      + dosageFile + ": "
                                      + line[fromParameters[1] + j * fromParameters[3] + 0] + " "
                                      + line[fromParameters[1] + j * fromParameters[3] + 1] + " "
                                      + line[fromParameters[1] + j * fromParameters[3] + 2]);
                    }
                  }
                }

                if (toParameters[3] == 1) {
                  if (fromParameters[3] > 1) {
                    dosageValue = genotypeProbabilities[0] * 2 + genotypeProbabilities[1] * 1;
                  }
                  writer.print(delimiter + ext.formDeci(dosageValue, toParameters[13],
                                                        toParameters[12] == toParameters[13]));
                } else {
                  writer.print(delimiter + ext.formDeci(genotypeProbabilities[0], toParameters[13],
                                                        toParameters[12] == toParameters[13]));
                  writer.print(delimiter + ext.formDeci(genotypeProbabilities[1], toParameters[13],
                                                        toParameters[12] == toParameters[13]));
                  if (toParameters[3] == 3) {
                    if (fromParameters[3] == 3) {
                      writer.print(delimiter
                                   + ext.formDeci(genotypeProbabilities[2], toParameters[13],
                                                  toParameters[12] == toParameters[13]));
                    } else {
                      writer.print(delimiter
                                   + ext.formDeci(1 - genotypeProbabilities[1]
                                                  - genotypeProbabilities[0], toParameters[13],
                                                  toParameters[12] == toParameters[13]));
                    }
                  }
                }
              }
              writer.println();
            }
          }
        } else if (fromParameters[2] == INDIVIDUAL_DOMINANT_FORMAT) {
          for (int i = 0; i < ids.length; i++) {
            line = reader.readLine().trim()
                         .split(fromParameters[10] == 1 ? "," : PSF.Regex.GREEDY_WHITESPACE);
            if (line.length - fromParameters[1] != markerNames.length * fromParameters[3]) {
              log.reportError("Error - mismatched number of elements in line "
                              + (i + 1 + fromParameters[4]) + " of " + dosageFile + "; expecting "
                              + markerNames.length + "*" + fromParameters[3] + "+"
                              + fromParameters[1] + ", found " + line.length);
              System.exit(1);
            }

            lead = LEADS[toParameters[11]] == null ? new String[toParameters[3]]
                                                   : LEADS[toParameters[11]];
            if (toParameters[0] == MACH_ID_TYPE) {
              lead[toParameters[5]] = ids[i][0] + "->" + ids[i][1];
            } else if (toParameters[0] == IID_TYPE) {
              lead[toParameters[5]] = ids[i][1];
            } else if (toParameters[0] == FID_IID_TYPE) {
              lead[toParameters[5]] = ids[i][0] + delimiter + ids[i][1];
            } else {
              log.reportError("Error - ID type has not been defined");
              System.exit(1);
            }
            writer.print(ArrayUtils.toStr(lead, delimiter));

            for (int j = 0; j < markerNames.length; j++) {
              if (extract == null || keeps.contains(markerNames[j])) {
                if (fromParameters[3] == 1) {
                  dosageValue = Float.parseFloat(line[fromParameters[1] + j]);
                } else {
                  genotypeProbabilities[0] = Float.parseFloat(line[fromParameters[1]
                                                                   + j * fromParameters[3] + 0]);
                  genotypeProbabilities[1] = Float.parseFloat(line[fromParameters[1]
                                                                   + j * fromParameters[3] + 1]);
                  if (fromParameters[3] == 3) {
                    genotypeProbabilities[2] = Float.parseFloat(line[fromParameters[1]
                                                                     + j * fromParameters[3] + 2]);
                    if (Math.abs((1 - genotypeProbabilities[1] - genotypeProbabilities[0])
                                 - genotypeProbabilities[2]) > 0.01) {
                      log.reportError("Error: P(BB) does not equal [ 1 - P(AA) - P(AB) ] for individual "
                                      + ids[i][0] + "," + ids[i][1] + " at marker " + markerNames[j]
                                      + " which is line " + (i + 1 + fromParameters[4]) + " of "
                                      + dosageFile + ": "
                                      + line[fromParameters[1] + j * fromParameters[3] + 0] + " "
                                      + line[fromParameters[1] + j * fromParameters[3] + 1] + " "
                                      + line[fromParameters[1] + j * fromParameters[3] + 2]);
                    }
                  }
                }

                if (toParameters[3] == 1) {
                  if (fromParameters[3] > 1) {
                    dosageValue = genotypeProbabilities[0] * 2 + genotypeProbabilities[1] * 1;
                  }
                  writer.print(delimiter + ext.formDeci(dosageValue, toParameters[13],
                                                        toParameters[12] == toParameters[13]));
                } else {
                  writer.print(delimiter + ext.formDeci(genotypeProbabilities[0], toParameters[13],
                                                        toParameters[12] == toParameters[13]));
                  writer.print(delimiter + ext.formDeci(genotypeProbabilities[1], toParameters[13],
                                                        toParameters[12] == toParameters[13]));
                  if (toParameters[3] == 3) {
                    if (genotypeProbabilities.length > 2) {
                      writer.print(delimiter
                                   + ext.formDeci(genotypeProbabilities[2], toParameters[13],
                                                  toParameters[12] == toParameters[13]));
                    } else {
                      writer.print(delimiter
                                   + ext.formDeci(1 - genotypeProbabilities[1]
                                                  - genotypeProbabilities[0], toParameters[13],
                                                  toParameters[12] == toParameters[13]));
                    }
                  }
                }
              }
            }
            writer.println();
          }
        }
        reader.close();
        writer.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + dosageFile + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dosageFile + "\"");
        System.exit(2);
      }
    }
  }

  public static int determineType(String dosageFileName) {
    String dosageFile = dosageFileName;
    if (dosageFile.endsWith(".gz")) {
      dosageFile = dosageFile.substring(0, dosageFile.length() - 3);
    } else if (dosageFile.endsWith(".zip")) {
      dosageFile = dosageFile.substring(0, dosageFile.length() - 4);
    }
    if (dosageFile.endsWith(".mldose")) {
      return MACH_MLDOSE_FORMAT;
    } else if (dosageFile.endsWith(".gen")) {
      return GEN_FORMAT;
    } else if (dosageFile.endsWith(".fhsR")) {
      return GWAF_FORMAT;
    } else if (dosageFile.endsWith(".dosage")) {
      return PLINK_FORMAT;
    } else if (dosageFile.endsWith(".mlprob")) {
      return MACH_MLPROB_FORMAT;
    } else if (dosageFile.endsWith(".dose")) {
      return MINIMAC_DOSE_FORMAT;
    } else if (dosageFile.endsWith(".beagle")) {
      return BEAGLE_DOSE_FORMAT;
    } else if (dosageFile.endsWith(".impute2") || dosageFile.endsWith(".imputed")) {
      return IMPUTE2_DOSE_FORMAT;
    } else if (dosageFile.endsWith(".db.xln")) {
      return DATABASE_DOSE_FORMAT;
    } else if (dosageFile.endsWith(".bed")) {
      return PLINK_BFILE_FORMAT;
    } else if (dosageFile.endsWith(".vcf")) {
      return VCF_FORMAT;
    } else if (dosageFile.endsWith(".bgen")) {
      return BGEN_FORMAT;
    } else if (dosageFile.endsWith(".frz.xln")) {
      return FREEZE5_FORMAT;
    } else {
      System.err.println("Error - format could not be deduced solely by the filename extension");
      return -1;
    }
  }

  public static DosageData load(String filename) {
    return (DosageData) SerializedFiles.readSerial(filename, true);
  }

  public static DosageData loadBGEN(String file, String indexFile, String idFile,
                                    int[][] regionsToKeep, String[] markersToKeep,
                                    String markerNamePrepend, Logger log) {
    BGENReader reader;
    try {
      reader = BGENReader.open(file, false);
    } catch (IOException e) {
      log.reportException(e);
      return null;
    }

    if (regionsToKeep != null && markersToKeep != null) {
      log.reportError("Extracting both regions and specific markers from BGEN files is not implemented.");
      return null;
    }

    if (indexFile != null && !"".equals(indexFile) && Files.exists(indexFile)) {
      try {
        BGENTools.loadMapInfo(reader, indexFile);
      } catch (ClassNotFoundException e) {
        log.reportError("Problem occured while loading index file {" + indexFile
                        + "}.  Loading will proceed without map info (likely much slower).");
        log.reportException(e);
      } catch (IOException e) {
        log.reportError("Problem occured while loading index file {" + indexFile
                        + "}.  Loading will proceed without map info (likely much slower).");
        log.reportException(e);
      }
    }

    Map<Integer, List<int[]>> regions = null;
    Collection<String> variants = null;

    if (regionsToKeep != null) {
      regions = new HashMap<>();
      for (int[] rgn : regionsToKeep) {
        List<int[]> chrRgn = regions.get(rgn[0]);
        if (chrRgn == null) {
          chrRgn = new ArrayList<>();
          regions.put(rgn[0], chrRgn);
        }
        chrRgn.add(new int[] {rgn[1], rgn[2]});
      }
    } else if (markersToKeep != null) {
      variants = new HashSet<>();
      for (String m : markersToKeep) {
        variants.add(m);
      }
    }

    if (Files.exists(idFile)) {
      String[] ids = HashVec.loadFileToStringArray(idFile, false, null, false);
      if (ids != null) {
        reader.setSamples(ids);
      }
    }

    DosageData dd = new DosageData();
    dd.ids = new String[(int) reader.getSampleCount()][];
    int count = 0;
    String id;
    for (int i = 0; i < reader.getSampleCount(); i++) {
      id = reader.getSamples()[i];
      String[] ids = id.split(PSF.Regex.GREEDY_WHITESPACE);
      if (ids.length == 2) {
        dd.ids[i] = ids;
      } else if (ids.length == 1) {
        dd.ids[i] = new String[] {id, id};
      } else {
        count++;
        dd.ids[i] = ids;
      }
    }
    if (count > 0) {
      log.reportTimeWarning(count + " IDs contained more than two parts!");
    }
    List<String> markersList = new ArrayList<>();
    List<Integer> chrs = new ArrayList<>();
    List<Integer> poss = new ArrayList<>();
    List<String[]> all = new ArrayList<>();
    List<float[][]> probs = new ArrayList<>();
    for (BGENRecord rec : (regions != null ? reader.query(regions)
                                           : (markersToKeep != null ? reader.query(variants)
                                                                    : reader))) {
      chrs.add(rec.getMetaData().getChr());
      poss.add((int) rec.getMetaData().getPos());
      markersList.add(rec.getMetaData().getRsID());
      all.add(new String[] {rec.getMetaData().getAlleles()[0], rec.getMetaData().getAlleles()[1]});
      float[][] probArr = new float[rec.getData().length][];
      for (int i = 0; i < rec.getData().length; i++) {
        probArr[i] = ArrayUtils.toFloatArray(rec.getData()[i]);
      }
      probs.add(probArr);
    }

    int numMarkers = markersList.size();
    dd.genotypeProbabilities = probs.toArray(new float[numMarkers][][]);
    dd.alleles = all.toArray(new String[numMarkers][]);
    dd.chrs = new byte[numMarkers];
    dd.positions = new int[numMarkers];
    for (int i = 0; i < numMarkers; i++) {
      dd.chrs[i] = (byte) chrs.get(i).intValue();
      dd.positions[i] = poss.get(i).intValue();
    }
    dd.labelPrepend = markerNamePrepend;
    String[] markerNames = markersList.toArray(new String[numMarkers]);

    dd.markerSet = new SnpMarkerSet(markerNames, dd.chrs, dd.positions, dd.alleles, null, false,
                                    false);

    dd.empty = markerNames.length == 0;

    markersList = null;
    chrs = null;
    poss = null;
    all = null;
    probs = null;
    markerNames = null;

    return dd;
  }

  /**
   * @param vcfFile VCF file to load
   * @param mapFile optional map file to load for IDs in VCF
   * @param regionsToKeep
   * @param markersToKeep
   * @param markerLocations
   * @param markerNamePrepend
   * @param log
   * @return
   */
  public static DosageData loadVCF(String vcfFile, @Nullable String mapFile, int[][] regionsToKeep,
                                   String[] markersToKeep, Map<String, int[]> markerLocations,
                                   String markerNamePrepend, Logger log) {
    String probTag = "GP";

    Map<String, String> segLookup = Maps.newHashMap();
    List<Segment> segsToQuery = new ArrayList<>();
    HashSet<String> markerSet = markersToKeep == null ? null : new HashSet<>();
    if (markersToKeep != null) {
      for (String s : markersToKeep) {
        markerSet.add(s);
      }
    }

    if (regionsToKeep != null) {
      for (int[] rgn : regionsToKeep) {
        segsToQuery.add(new Segment((byte) rgn[0], rgn[1], rgn[2]));
      }
    }

    if (markerLocations != null) {
      for (Entry<String, int[]> mkrLoc : markerLocations.entrySet()) {
        markerSet.add(mkrLoc.getKey());
        segsToQuery.add(new Segment((byte) mkrLoc.getValue()[0], mkrLoc.getValue()[1],
                                    mkrLoc.getValue()[1]));
        segLookup.put(mkrLoc.getValue()[0] + ":" + mkrLoc.getValue()[1], mkrLoc.getKey());
        segLookup.put("chr" + mkrLoc.getValue()[0] + ":" + mkrLoc.getValue()[1], mkrLoc.getKey());
      }
    }

    if (segsToQuery.size() > 0) {
      Collections.sort(segsToQuery);
    }

    VCFHeader header;
    List<VariantContext> keepList = Lists.newArrayList();

    try (VCFFileReader reader = new VCFFileReader(new File(vcfFile),
                                                  Files.exists(vcfFile + ".tbi"))) {
      header = reader.getFileHeader();

      // discover which chrs are in this vcf file
      BidiMap<String, Integer> contigMap = new DualHashBidiMap<>();
      header.getContigLines().forEach(vch -> {
        // ensure parsability
        contigMap.put(vch.getID(), (int) Positions.chromosomeNumber(vch.getID()));
      });

      // remove segs on unavailable chrs
      for (int i = segsToQuery.size() - 1; i >= 0; i--) {
        if (!contigMap.containsValue((int) segsToQuery.get(i).getChr())) {
          segsToQuery.remove(i);
        }
      }

      // remove markers if we have a location for them and that location isn't available
      if (markerLocations != null) {
        markerSet.removeIf(m -> {
          return markerLocations.containsKey(m)
                 && !contigMap.containsValue(markerLocations.get(m)[0]);
        });
      }

      for (Segment seg : segsToQuery) {
        CloseableIterator<VariantContext> iter = reader.query(contigMap.getKey(seg.getChr()),
                                                              seg.getStart() - 1,
                                                              seg.getStop() + 1);
        while (iter.hasNext()) {
          VariantContext vc = iter.next();
          if (markerSet.contains(vc.getID())
              || segLookup.containsKey(vc.getContig() + ":" + vc.getStart())
              || segsToQuery.contains(new Segment(Positions.chromosomeNumber(vc.getContig()),
                                                  vc.getStart(), vc.getEnd()))) {
            keepList.add(vc);
            String id = markerSet.contains(vc.getID()) ? vc.getID()
                                                       : segLookup.containsKey(vc.getContig() + ":"
                                                                               + vc.getStart()) ? segLookup.get((vc.getContig() + ":" + vc.getStart())) : vc.getID();
            markerSet.remove(id);
          }
        }
      }
    }

    if (markerSet.size() > 0) {
      log.reportTime("Scanning through " + vcfFile + " for " + markerSet.size()
                     + " remaining variants.");
      try (VCFFileReader reader = new VCFFileReader(new File(vcfFile),
                                                    Files.exists(vcfFile + ".tbi"))) {
        keepList.addAll(reader.iterator().stream().filter(vc -> {
          return markerSet.contains(vc.getID());
        }).collect(Collectors.toList()));
      }
    }

    DosageData dd = new DosageData();

    Map<ImmutableList<Object>, DataLine> mapLines = null;
    FileColumn<String> fid = new ExplicitIndexedFileColumn("FID", PSF.Plink.FAM_FID_INDEX);
    FileColumn<String> iid = new ExplicitIndexedFileColumn("IID", PSF.Plink.FAM_IID_INDEX);
    FileColumn<String> father = new ExplicitIndexedFileColumn("Father", PSF.Plink.FAM_FA_INDEX);
    FileColumn<String> mother = new ExplicitIndexedFileColumn("Mother", PSF.Plink.FAM_MO_INDEX);
    FileColumn<String> sex = new ExplicitIndexedFileColumn("Sex", PSF.Plink.FAM_SEX_INDEX);
    FileColumn<String> pheno = new ExplicitIndexedFileColumn("Pheno", PSF.Plink.FAM_AFF_INDEX);
    if (mapFile != null && !ext.isMissingValue(mapFile)) {
      try {
        mapLines = FileParserFactory.setup(mapFile, fid, iid, father, mother, sex, pheno).noHeader()
                                    .build().load(true, fid, iid);
      } catch (IOException e) {
        log.reportError("Failed to load " + mapFile + ", not loading fam data for " + vcfFile);
      }
    }
    dd.ids = new String[header.getSampleNamesInOrder().size()][];
    for (int i = 0; i < header.getSampleNamesInOrder().size(); i++) {
      String id = header.getSampleNamesInOrder().get(i);
      DataLine mapLine = mapLines.get(ImmutableList.of(id, id));
      final String[] idLine;
      if (mapLine != null) {
        idLine = new String[PSF.Plink.FAM_FIELD_COUNT];
        idLine[PSF.Plink.FAM_FID_INDEX] = id;
        idLine[PSF.Plink.FAM_IID_INDEX] = id;
        idLine[PSF.Plink.FAM_FA_INDEX] = mapLine.getString(father);
        idLine[PSF.Plink.FAM_MO_INDEX] = mapLine.getString(mother);
        idLine[PSF.Plink.FAM_SEX_INDEX] = mapLine.getString(sex);
        idLine[PSF.Plink.FAM_AFF_INDEX] = mapLine.getString(pheno);
      } else {
        idLine = new String[] {id, id};
      }
      dd.ids[i] = idLine;
    }

    int numMarkers = keepList.size();
    dd.genotypeProbabilities = new float[numMarkers][dd.ids.length][];
    dd.alleles = new String[numMarkers][2];
    dd.chrs = new byte[numMarkers];
    dd.positions = new int[numMarkers];
    dd.labelPrepend = markerNamePrepend;
    String[] markerNames = new String[numMarkers];

    VariantContext vc;
    for (int i = 0; i < numMarkers; i++) {
      vc = keepList.get(i);
      markerNames[i] = vc.getID();
      if (markerNamePrepend != null && !markerNames[i].startsWith(markerNamePrepend)) {
        markerNames[i] = markerNamePrepend + (markerNamePrepend.endsWith("_") ? "" : "_")
                         + markerNames[i];
      }
      dd.chrs[i] = Positions.chromosomeNumber(vc.getContig());
      dd.positions[i] = vc.getStart();
      dd.alleles[i] = new String[vc.getAlleles().size()];
      dd.alleles[i][0] = vc.getReference().getBaseString();
      for (int j = 1; j < dd.alleles[i].length; j++) {
        dd.alleles[i][j] = vc.getAlternateAlleles().get(j - 1).getBaseString();
      }
      GenotypesContext gc = vc.getGenotypes();
      for (int j = 0; j < gc.size(); j++) {
        Genotype g = gc.get(j);
        Object o = g.getExtendedAttribute(probTag);
        float[] probs = ArrayUtils.floatArray(3, Float.NaN);
        if (o != null) {
          String[] probStr = ((String) o).split(",");
          probs[0] = Float.parseFloat(probStr[0]);
          probs[1] = Float.parseFloat(probStr[1]);
          probs[2] = Float.parseFloat(probStr[2]);
        }
        dd.genotypeProbabilities[i][j] = probs;
      }
    }

    dd.markerSet = new SnpMarkerSet(markerNames, dd.chrs, dd.positions, dd.alleles, null, false,
                                    false);

    dd.empty = markerNames.length == 0;

    return dd;
  }

  private static boolean keepVariant(VariantContext e, Map<Integer, List<int[]>> rgnsToKeep,
                                     Set<String> markerSet) {
    if ((rgnsToKeep == null || rgnsToKeep.isEmpty()) && markerSet == null) return true;
    if (rgnsToKeep != null && !rgnsToKeep.isEmpty()) {
      int chr = Positions.chromosomeNumber(e.getContig());
      int pos = e.getStart();
      for (int[] rgn : rgnsToKeep.get(chr)) {
        if ((rgn.length == 2 && rgn[1] == pos)
            || (rgn.length == 3 && rgn[1] <= pos && rgn[2] >= pos)) {
          return true;
        }
      }
    }
    if (markerSet != null && markerSet.contains(e.getID())) {
      return true;
    }

    return false;
  }

  public static DosageData loadPlinkBinary(String dir, int[][] regionsToKeep,
                                           String[] markersToKeep, String plinkRoot,
                                           String markerNamePrepend, boolean loadMissingAsNaN,
                                           boolean includeAllMarkers) {
    DosageData dd = new DosageData();

    String[][] bimData = HashVec.loadFileToStringMatrix(dir + plinkRoot + ".bim", false,
                                                        new int[] {0, 1, 2, 3, 4, 5});
    HashSet<String> markerSet = markersToKeep == null ? null : new HashSet<>();
    if (markersToKeep != null) {
      for (String s : markersToKeep) {
        markerSet.add(s);
      }
    }
    HashSet<String> bimMarkers = new HashSet<>();
    for (String[] element : bimData) {
      bimMarkers.add(element[1]);
    }

    boolean[] markersToInclude = ArrayUtils.booleanArray(bimData.length, false);
    markers: for (int i = 0; i < bimData.length; i++) {
      if ((markerSet != null && markerSet.contains(bimData[i][1])) || includeAllMarkers) {
        markersToInclude[i] = true;
        continue;
      }
      if (regionsToKeep != null) {
        int chr = decodeChr(bimData[i][0]);
        int pos = Integer.parseInt(bimData[i][3]);
        for (int[] rgn : regionsToKeep) {
          if (rgn[0] == chr) {
            if ((rgn.length == 2 && rgn[1] == pos)
                || (rgn.length == 3 && rgn[1] <= pos && rgn[2] >= pos)) {
              markersToInclude[i] = true;
              continue markers;
            }
          }
        }
      }
    }

    dd.ids = HashVec.loadFileToStringMatrix(dir + plinkRoot + ".fam", false,
                                            new int[] {0, 1, 2, 3, 4, 5});
    int numMarkers = ArrayUtils.booleanArraySum(markersToInclude);// Files.countLines(dir +
                                                                  // plinkRoot +
                                                                  // ".bim", 0);
    dd.alleles = new String[numMarkers][2];
    dd.chrs = new byte[numMarkers];
    dd.positions = new int[numMarkers];

    String[] markerNames = new String[numMarkers];
    int index = 0;
    for (int i = 0; i < bimData.length; i++) {
      if (!markersToInclude[i]) {
        continue;
      }
      markerNames[index] = (null == markerNamePrepend ? ""
                                                      : (markerNamePrepend.endsWith("_") ? markerNamePrepend
                                                                                         : markerNamePrepend
                                                                                           + "_"))
                           + bimData[i][1];
      dd.chrs[index] = decodeChr(bimData[i][0]);
      dd.positions[index] = Integer.parseInt(bimData[i][3]);
      dd.alleles[index][0] = bimData[i][4];
      dd.alleles[index][1] = bimData[i][5];
      index++;
    }
    dd.markerSet = new SnpMarkerSet(markerNames, dd.chrs, dd.positions, dd.alleles, null, false,
                                    false);
    dd.dosageValues = new float[markerNames.length][dd.ids.length];

    RandomAccessFile in;
    try {
      in = new RandomAccessFile(dir + plinkRoot + ".bed", "r");
      byte[] magicBytes = new byte[3];
      in.read(magicBytes);
      if (magicBytes[2] == 0) {
        System.err.println("Error - .bed file is sample-dominant.");
      } else {
        int famCnt = dd.ids.length;
        int blockSize = (int) Math.ceil(famCnt / 4.0);

        index = -1;
        for (int i = 0; i < bimData.length; i++) {
          in.seek(3 + i * blockSize); // ALWAYS seek forward
          if (!markersToInclude[i]) {
            continue;
          }
          index++;
          if (dd.positions[index] == -1) {
            dd.dosageValues[index] = ArrayUtils.floatArray(dd.ids.length, -1); // TODO is this
                                                                               // correct
                                                                               // code? i.e. will
                                                                               // this
                                                                               // ever occur?
            continue;
          }

          byte[] markerBytes = new byte[blockSize];
          byte[] sampGeno = new byte[dd.ids.length];
          in.read(markerBytes);
          for (int bitInd = 0; bitInd < markerBytes.length; bitInd++) {
            byte bedByte = markerBytes[bitInd];
            byte[] genotypes = decodeBedByte(bedByte);

            for (int g = 0; g < genotypes.length; g++) {
              int idInd = bitInd * 4 + g;
              if (idInd == -1 || idInd >= sampGeno.length) {
                continue;
              }
              sampGeno[idInd] = genotypes[g];
            }
          }
          dd.dosageValues[index] = ArrayUtils.toFloatArray(sampGeno);
          if (loadMissingAsNaN) {
            for (int n = 0; n < dd.dosageValues[index].length; n++) {
              if (dd.dosageValues[index][n] == -1) {
                dd.dosageValues[index][n] = Float.NaN;
              }
            }
          }
        }
      }

      in.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    } catch (Elision e) {
      e.printStackTrace();
    }

    return dd;
  }

  private static byte decodeChr(String chrStr) {
    int chr = 0;
    try {
      chr = Integer.parseInt(chrStr);
      return (byte) chr;
    } catch (NumberFormatException e) {
      if ("X".equals(chrStr)) {
        return 23;
      }
      if ("Y".equals(chrStr)) {
        return 24;
      }
      if ("XY".equals(chrStr)) {
        return 25;
      }
      if ("MT".equals(chrStr)) {
        return 26;
      }
    }
    System.err.println("ERROR - unrecognized PLINK .bim chr: " + chrStr);
    return 0;
  }

  public static String[] loadMarkers(String markersToUseFile, Logger log) {
    String[] markers = null;

    if (markersToUseFile != null && !"".equals(markersToUseFile)) {
      if (Files.exists(markersToUseFile)) {
        markers = HashVec.loadFileToStringArray(markersToUseFile, false, new int[] {0}, true, false,
                                                "\t");
      }
    } else {
      log.reportError("Error - specified markers file: \"" + markersToUseFile
                      + "\" doesn't exist!");
    }

    return markers;
  }

  private static int[][] loadRegions(String regionsToUseFile, Logger log) {
    int[][] regions = null;

    if (regionsToUseFile != null && !"".equals(regionsToUseFile)) {
      if (Files.exists(regionsToUseFile)) {
        String[] rgns = HashVec.loadFileToStringArray(regionsToUseFile, false, new int[] {0}, true,
                                                      false, "\t");
        regions = new int[rgns.length][];
        for (int i = 0; i < rgns.length; i++) {
          regions[i] = Positions.parseUCSClocation(rgns[i]);
        }
      }
    } else {
      log.reportError("Error - specified regions file: \"" + regionsToUseFile
                      + "\" doesn't exist!");
    }

    return regions;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "chr21.mldose";
    String idFile = "chr21.fam";
    String mapFile = "chr21.mlinfo";
    String outfile = "output.dosage";
    String mapOut = null;
    int from = -1;
    int to = -1;
    String logfile = null;
    Logger log;
    String extractMkrs = null;
    String extractRgns = null;
    boolean awk = false;
    long date;

    String usage = "\n" + "filesys.DosageData requires 0-1 arguments\n"
                   + "   (1) name of dosage file to convert (i.e. dosageFile=" + filename
                   + " (default))\n" + "   (2) name of file with ids listed (i.e. idFile=" + idFile
                   + " (default))\n" + "   (3) name of associated map file (i.e. mapFile=" + mapFile
                   + " (default))\n" + "   (4) name of new dosage file (i.e. out=" + outfile
                   + " (default))\n" + "   (5) name of new map file (i.e. mapOut=" + mapOut
                   + " (default))\n" + "   (6) filename for log (i.e. log=" + logfile
                   + " (default) (optional))\n"
                   + "   (7a) name of file listing variants to extract (i.e. extractMarkers="
                   + extractMkrs + " (default) (optional))\n"
                   + "   (7b) name of file listing ranges to extract (i.e. extractRanges="
                   + extractRgns + " (default) (optional))\n"
                   + "   (8) file type of original file (i.e. from=" + from + " (default))\n"
                   + "   (9) file type of file to be generated (i.e. to=" + to + " (default))\n"
                   + "  (10) use sed/cut instead (i.e. -awk (not the default))\n" + " \n"
                   + "       Type Extension Description\n"
                   + "        -1   [any]    auto-detect from extension\n"
                   + "        0   .mldose   MACH style .mldose file\n"
                   + "        1   .gen      probabilites for all three genotypes, plus map information\n"
                   + "        2   .fhsR     input file for GWAF\n"
                   + "        3   .dosage   input file for PLINK\n"
                   + "        4   .mlprob   MACH style .mlprob file\n"
                   + "        5   .dose     output from MINIMAC\n"
                   + "        6   .impute2  output from IMPUTE2\n"
                   + "        7   .db.xln   database format\n"
                   + "        8   .dose     BEAGLE format\n"
                   + "        [Files may also be zipped or gzipped]" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dosageFile=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("idFile=")) {
        idFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("mapFile=")) {
        mapFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        outfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("mapOut=")) {
        mapOut = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("extractMarkers=")) {
        extractMkrs = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("extractRanges=")) {
        extractRgns = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("from=")) {
        from = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("to=")) {
        to = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("-awk")) {
        awk = true;
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    // filename =
    // "C:\\CARe_data\\CARe_imputed_all_llange_24mar2010\\CFS\\IBC\\whites\\00src\\leslie_lange.CFS.IBC.CEU.chr9.gen";
    // idFile =
    // "C:\\CARe_data\\CARe_imputed_all_llange_24mar2010\\CFS\\IBC\\whites\\leslie_lange.CFS.IBC.CEU.chr9.pfam";
    // mapFile =
    // "C:\\CARe_data\\CARe_imputed_all_llange_24mar2010\\CFS\\IBC\\whites\\leslie_lange.CFS.IBC.CEU.chr9.pmap";
    // outfile =
    // "C:\\CARe_data\\CARe_imputed_all_llange_24mar2010\\CFS\\IBC\\whites\\00src\\leslie_lange.CFS.IBC.CEU.chr9.fhsR";

    // filename =
    // "D:\\CARe\\CARe_imputed_all_llange_24mar2010\\ARIC\\IBC\\whites\\leslie_lange.ARIC.IBC.CEU.chr9.gen";
    // idFile =
    // "D:\\CARe\\CARe_imputed_all_llange_24mar2010\\ARIC\\IBC\\whites\\leslie_lange.ARIC.IBC.CEU.chr9.pfam";
    // mapFile =
    // "D:\\CARe\\CARe_imputed_all_llange_24mar2010\\ARIC\\IBC\\whites\\leslie_lange.ARIC.IBC.CEU.chr9.mlinfo";
    // outfile = "D:\\CARe\\CARe_imputed_all_llange_24mar2010\\ARIC\\IBC\\whites\\candi\\abo.gen";
    // extract = "D:\\CARe\\CARe_imputed_all_llange_24mar2010\\abo_snplist.txt";

    // filename = "chr4.dose";
    // mapFile = "chr4.minfo";
    // idFile = "plink.fam";
    // outfile = "first_half.dose";
    // extract = "first_half.txt";
    // awk = true;
    //
    // filename = "Fung/hla.dose";
    // mapFile = "Fung/hla.minfo";
    // idFile = "Fung/hla.ids.fam";
    // outfile = "hla_test.dose";
    // extract = "test_list.txt";
    // awk = true;

    // filename = "D:/LITE/1000Gimputation/blacks/chr9.all_parsed.xln";
    // idFile = "D:/LITE/1000Gimputation/blacks/B_ARIC_gwas_frz3v2_clean.b37.fam";
    // mapFile = "D:/LITE/1000Gimputation/blacks/ABO_blacks_fullRegion.bim";
    // outfile = "D:/LITE/1000Gimputation/blacks/ABO_blacks_fullRegion.db.xln";
    // from=6;
    // to=7;

    from = from == -1 ? determineType(filename) : from;
    to = to == -1 ? determineType(outfile) : to;

    try {
      log = new Logger(logfile);
      if (PARAMETERS[from][3] < PARAMETERS[to][3]) {
        log.reportError("Error - cannot convert dosage values to genotype probabilities");
      } else if (PARAMETERS[from][2] != PARAMETERS[to][2]) {
        log.reportError("Conversion will have to take place in memory in order to transpose between marker dominant/individual dominant");
        new DosageData(filename, idFile, mapFile, from, extractRgns, extractMkrs, true,
                       log).writeToFile(outfile, mapOut, extractMkrs, extractRgns, to, log);
      } else {
        date = new Date().getTime();
        convert(filename, idFile, mapFile, from, outfile, mapOut, extractMkrs, to, awk, true, log);
        System.out.println("Time elapsed: " + ext.getTimeElapsed(date));
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}