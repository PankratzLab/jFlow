// expand to dynamically load/save a certain chunk of markers at a time
package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.gwas.Plink;
import org.genvisis.stats.LeastSquares;
import org.genvisis.stats.LogisticRegression;
import org.genvisis.stats.RegressionModel;

public class GenotypeMatrix implements Serializable {
  public static final long serialVersionUID = 1L;
  public static final int CHARGE_S_HOUSTON_FORMAT = 0;
  public static final int CHARGE_S_BOSTON_FORMAT = 1;

  public static final int IID_TYPE = 0;

  public static final int INDIVIDUAL_DOMINANT_FORMAT = 0;
  public static final int MARKER_DOMINANT_FORMAT = 1;

  // public static final String[][] HEADS = {{}, {"id"}};
  public static final String[][] HEADS = {{}, {"CHROM", "POS"}};
  public static final String[] DELIMITERS = {"\t", ",", " "};
  public static final String[] MISSING_VALUES = {".", "NA"};

  /** 0, 1, 2, 3, 4, 5, 6 */
  /**
   * dominance format, delimiter, header row, index of which header head, marker/IID index, column
   * index where genotype counts begin, split chrPos instead of markerName
   */

  /** 0 1 2 3 4 5 6 */
  public static final int[][] PARAMETERS = {{MARKER_DOMINANT_FORMAT, 1, 1, 0, 0, 1, 0}, // .csv
                                                                                        // (ChargeS
                                                                                        // Houston
                                                                                        // format)
      {MARKER_DOMINANT_FORMAT, 0, 1, 1, 0, 2, 1}, // .txt (ChargeS Boston format)

  };

  public static int determineType(String dosageFile) {
    if (dosageFile.endsWith(".csv")) {
      return CHARGE_S_HOUSTON_FORMAT;
    } else if (dosageFile.endsWith(".txt.gz")) {
      return CHARGE_S_BOSTON_FORMAT;
    } else {
      System.err.println("Error - format could not be deduced solely by the filename extension");
      return -1;
    }
  }

  public static GenotypeMatrix load(String filename, boolean jar) {
    return (GenotypeMatrix) SerializedFiles.readSerial(filename, jar, true);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "SeqGenotypes.dat";
    String markerFile = null;
    String phenoFile = null;
    String logfile = null;
    GenotypeMatrix gens;
    Logger log;
    String dir;

    String usage = "\n" + "filesys.GenotypeMatrix requires 0-1 arguments\n"
        + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      log = new Logger(logfile);
      dir = "D:/LITE/CHARGE-S/ARIC_CHAGE_S_Freeze2/";
      if (!Files.exists(dir)) {
        dir = "C:/LITE/CHARGE-S/ARIC_CHAGE_S_Freeze2/";
      }
      filename = dir + "aric_genotypes_frz2_final.csv";
      filename = dir + "slim.csv";
      markerFile = dir + "SNPInfo_ExomeFreeze2_120810.csv";
      markerFile = dir + "slimMarkers2.csv";
      phenoFile = dir + "phenoCensored2.csv";
      phenoFile = dir + "lnFibrinogen_alanna.csv";
      if (Files.exists(filename + ".ser")) {
        System.out.println("Loading serialized version: " + filename + ".ser");
        gens = load(filename + ".ser", false);
      } else {
        System.out.println("Loading: " + filename);
        gens = new GenotypeMatrix(filename, null, markerFile, log);
        gens.serialize(filename + ".ser");
        gens.writeToPlinkFiles(ext.rootOf(filename));
      }
      System.out.println("Analyzing: " + phenoFile);
      gens.analyze(phenoFile, "NA", null, true, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private SnpMarkerSet markerSet;

  private String[] ids;

  private byte[][] values;

  public GenotypeMatrix(String genotypeFile, String idFile, String markerFile, int type,
      Logger log) {
    this(genotypeFile, idFile, markerFile, PARAMETERS[type], log);
  }

  public GenotypeMatrix(String genotypeFile, String idFile, String markerFile, int[] parameters,
      Logger log) {
    BufferedReader reader;
    String[] line;
    String[] markerNames;
    int dominance;
    String delimiter;
    boolean headerRow;
    int markerOrIdIndex, firstGenotypeCol;
    String[] headerHead;
    int markerFileType;
    boolean chrPosCombo, problem;

    if (!Files.exists(markerFile)) {
      log.reportError("Error - file \"" + markerFile + "\" not found");
      return;
    }

    log.report("Loading map from: " + markerFile);
    markerFileType = SnpMarkerSet.determineType(markerFile);
    if (markerFileType == -1) {
      markerSet = new SnpMarkerSet(markerFile, Files.determineDelimiter(markerFile, log), null,
          true, new String[][] {{"\"", ""}}, false, log);
      markerSet.sortMarkers();
    } else {
      markerSet = new SnpMarkerSet(markerFile, false, log);
    }
    markerNames = markerSet.getMarkerNames();

    try {
      // reader = new BufferedReader(new FileReader(genotypeFile));
      reader = Files.getAppropriateReader(genotypeFile);

      // can be much more complex if you want, just see DosageData
      if (idFile != null) {
        ids = HashVec.loadFileToStringArray(idFile, false, true, new int[] {0}, false, false,
            Files.determineDelimiter(idFile, log));
      }

      dominance = parameters[0];
      delimiter = DELIMITERS[parameters[1]];
      headerRow = parameters[2] == 1;
      headerHead = HEADS[parameters[3]];
      markerOrIdIndex = parameters[4];
      firstGenotypeCol = parameters[5];
      chrPosCombo = parameters[6] == 1;
      if (headerRow) {
        // line = ext.replaceAllWith(reader.readLine(), "\"", "") .trim().split(delimiter);
        if (delimiter.equals(",")) {
          line = ext.splitCommasIntelligently(reader.readLine(), true, log);
        } else {
          line = ext.replaceAllWith(reader.readLine(), "\"", "").trim().split(delimiter);
        }
        if (dominance == INDIVIDUAL_DOMINANT_FORMAT) {
          for (int i = 0; i < markerNames.length; i++) {
            if (!markerNames[i].equals(line[headerHead.length + i])) {
              log.reportError("Error - mismatched name at marker " + (i + 1) + " of " + genotypeFile
                  + "; expecting " + markerNames[i] + " given map file " + markerFile + ", found "
                  + line[headerHead.length + i]);
              System.exit(1);
            }
          }
        } else if (dominance == MARKER_DOMINANT_FORMAT) {
          problem = false;
          for (int i = 0; i < headerHead.length; i++) {
            if (!line[i].equals(headerHead[i])) {
              problem = true;
            }
          }
          if (problem) {
            log.reportError("Error - mismatched head of header row: expecting '"
                + Array.toStr(headerHead, "/") + "' but found '"
                + Array.toStr(Array.subArray(line, 0, headerHead.length), "/") + "'");
          }
          if (ids == null) {
            ids = new String[line.length - headerHead.length];
            for (int i = 0; i < ids.length; i++) {
              ids[i] = line[headerHead.length + i];
            }
          } else {
            for (int i = 0; i < ids.length; i++) {
              if (ids != null && !ids[i].equals(line[headerHead.length + i])) {
                log.reportError("Error - mismatched IDs at individual " + (i + 1) + " of "
                    + genotypeFile + "; expecting " + ids[i] + " given id file " + idFile
                    + ", found " + line[headerHead.length + i]);
                System.exit(1);
              }
            }
          }
        }
      } else if (ids == null) {
        System.err.println(
            "Error - not fully developed yet; need a header row in order to get the ids/markers");
        reader.close();
        return;
      }

      values = new byte[markerNames.length][ids.length];

      if (dominance == MARKER_DOMINANT_FORMAT) {
        for (int i = 0; i < markerNames.length; i++) {
          line = ext.replaceAllWith(reader.readLine(), "\"", "").trim().split(delimiter);
          if (chrPosCombo) {
            line[0] = "chr" + line[0] + ":" + line[1];
          }
          if (!markerNames[i].equals(line[markerOrIdIndex])) {
            log.reportError("Error - mismatched name at marker " + (i + 1) + " of " + genotypeFile
                + "; expecting " + markerNames[i] + " given map file " + markerFile + ", found "
                + line[markerOrIdIndex]);
            System.exit(1);
          }

          if (line.length - firstGenotypeCol != ids.length) {
            log.reportError("Error - mismatched number of elements in line "
                + (i + 1 + (headerRow ? 1 : 0)) + " of " + genotypeFile + "; expecting "
                + ids.length + "+" + firstGenotypeCol + ", found " + line.length);
            log.reportError("First few ids: '" + ids[0] + "', '" + ids[1] + "', '" + ids[2]
                + "'... '" + ids[ids.length - 2] + "', '" + ids[ids.length - 1] + "'");
            System.exit(1);
          }
          for (int j = 0; j < ids.length; j++) {
            if (ext.isMissingValue(line[firstGenotypeCol + j])) {
              values[i][j] = -1;
            } else {
              values[i][j] = Byte.parseByte(line[firstGenotypeCol + j]);
            }
          }
        }
      } else if (dominance == INDIVIDUAL_DOMINANT_FORMAT) {
        for (int i = 0; i < ids.length; i++) {
          line = ext.replaceAllWith(reader.readLine(), "\"", "").trim().split(delimiter);
          if (line.length - firstGenotypeCol != markerNames.length) {
            log.reportError("Error - mismatched number of elements in line "
                + (i + 1 + (headerRow ? 1 : 0)) + " of " + genotypeFile + "; expecting "
                + markerNames.length + "+" + firstGenotypeCol + ", found " + line.length);
            System.exit(1);
          }
          for (int j = 0; j < markerNames.length; j++) {
            values[j][i] = Byte.parseByte(line[firstGenotypeCol + j]);
          }
        }
      }

      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + genotypeFile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + genotypeFile + "\"");
      System.exit(2);
    }
  }

  public GenotypeMatrix(String genotypeFile, String idFile, String markerFile, Logger log) {
    this(genotypeFile, idFile, markerFile, determineType(genotypeFile), log);
  }

  public void analyze(String phenoFile, String phenoMissingValue, String snpList, boolean verbose,
      Logger log) {
    PrintWriter writer, w2;
    String[] line;
    Hashtable<String, String> hash;
    HashSet<String> snps;
    int count, countUsed;
    String[] traits, markerNames;
    boolean[] use, analyze;
    byte[] chrs;
    int[] positions;
    double[] deps;
    double[][] indeps;
    boolean logistic;
    RegressionModel model;
    char[][] alleles;
    double[] betas, stderrs, pvals, stats;
    String[] names, namesUsed;

    traits = Files.getHeaderOfFile(phenoFile, Files.determineDelimiter(phenoFile, log), log);
    names = Array.subArray(traits, 1);
    hash = HashVec.loadFileToHashString(phenoFile, new int[] {0},
        Arrays.copyOfRange(Array.arrayOfIndices(traits.length), 1, traits.length),
        phenoFile.endsWith(".csv"), "\t", true, false, false);
    traits = Array.subArray(traits, 1);
    log.report("Missing phenotype is set to '" + phenoMissingValue + "'");

    markerNames = markerSet.getMarkerNames();
    chrs = markerSet.getChrs();
    positions = markerSet.getPositions();
    alleles = markerSet.getAlleles();
    analyze = Array.booleanArray(markerNames.length, true);
    if (snpList != null) {
      snps = HashVec.loadFileToHashSet(snpList, false);
      for (int i = 0; i < markerNames.length; i++) {
        if (!snps.contains(markerNames[i])) {
          analyze[i] = false;
        }
      }
    }

    use = Array.booleanArray(ids.length, true);
    for (int i = 0; i < ids.length; i++) {
      if (hash.containsKey(ids[i])) {
        line = hash.get(ids[i]).split("[\\s]+");
        for (int j = 0; j < line.length; j++) {
          if (!ext.isValidDouble(line[j]) || line[j].equals(phenoMissingValue)) {
            use[i] = false;
          }
        }
      } else {
        use[i] = false;
      }
    }

    deps = new double[Array.booleanArraySum(use)];
    indeps = new double[deps.length][traits.length];
    log.report("There are " + deps.length + " samples with complete data", true, verbose);
    count = 0;
    for (int i = 0; i < ids.length; i++) {
      if (use[i]) {
        line = hash.get(ids[i]).split("[\\s]+");
        deps[count] = Double.parseDouble(line[0]);
        for (int j = 1; j < traits.length; j++) {
          indeps[count][j] = Double.parseDouble(line[j]);
        }
        count++;
      }
    }
    logistic = RegressionModel.isBinaryTrait(Array.toStr(deps).split("[\\s]+"), log);
    log.report(
        "Running a " + (logistic ? "logistic" : "linear") + " model for trait '" + traits[0] + "'",
        true, verbose);
    try {
      writer = new PrintWriter(new FileWriter(
          ext.rootOf(phenoFile, false) + ".results." + (logistic ? "logistic" : "linear")));
      w2 = new PrintWriter(new FileWriter(ext.rootOf(phenoFile, false) + ".se.metal"));
      String[] arr = logistic ? Plink.LOGISTIC_SE_HEADER : Plink.LINEAR_SE_HEADER;
      line = Arrays.copyOf(arr, arr.length);
      line[1] = line[1] + "      ";
      line[2] = line[1] + "      ";
      writer.println(Array.toStr(line));
      // w2.println("MARKER\tREF\tOTHER\tN\tDIR\tPVALUE\tbeta\tSE");
      w2.println("MarkerName\tAllele1\tAllele2\tWeight\tDirection\tP-value\tEffect\tStdErr");
      // public static final String[] LOGISTIC_SE_HEADER = {"CHR", "SNP", "BP", "A1", "TEST",
      // "NMISS", "OR", "SE", "L95", "U95", "STAT", "P"};

      for (int i = 0; i < markerNames.length; i++) {
        if (analyze[i]) {
          count = 0;
          countUsed = 0;
          for (int j = 0; j < ids.length; j++) {
            if (use[j]) {
              if (values[i][j] < 0) {
                indeps[count][0] = Double.NaN;
              } else {
                indeps[count][0] = values[i][j];
                countUsed++;
              }
              count++;
            }
          }
          names[0] = "SNP";
          model = logistic ? new LogisticRegression(deps, indeps, names, false, false)
              : new LeastSquares(deps, indeps, names, false, false);
          betas = model.getBetas();
          stderrs = model.getSEofBs();
          pvals = model.getSigs();
          stats = model.getStats();
          int sigfig = 4;
          // System.err.println(betas.length+"\t"+traits.length);
          namesUsed = model.getVarNames();
          if (model.getFinalDependentVariables().length != countUsed) {
            log.reportError("Mismatched number of missing values for marker '" + markerNames[i]
                + "' (" + model.getFinalDependentVariables().length + ", expected " + countUsed
                + "); might want to check missing value codes: "
                + model.getFinalIndependentVariables()[2][0]);
          }
          if (!namesUsed[1].equals("SNP")) {
            writer.println(chrs[i] + "\t" + markerNames[i] + "\t" + positions[i] + "\t"
                + alleles[i][0] + "\tADD\t.\t.\t.\t.\t.\t.\t.");
            w2.println(markerNames[i] + "\t" + alleles[i][0] + "\t" + alleles[i][1] + "\t"
                + deps.length + "\t.\t.\t.\t.");
          } else {
            writer.println(
                chrs[i] + "\t" + markerNames[i] + "\t" + positions[i] + "\t" + alleles[i][0]
                    + "\tADD\t" + countUsed + "\t" + ext.formDeci(betas[1], sigfig, true) + "\t"
                    + ext.formDeci(stderrs[1], sigfig, true) + "\t.\t.\t"
                    + (stats[1] + "    ").substring(0, 6).trim() + "\t"
                    + ext.prettyP(pvals[1], sigfig, 4, 3, true));
            w2.println(markerNames[i] + "\t" + alleles[i][0] + "\t" + alleles[i][1] + "\t"
                + countUsed + "\t" + (betas[1] == 0 ? 0 : (betas[1] > 0 ? "+" : "-")) + "\t"
                + ext.prettyP(pvals[1], sigfig, 4, 3, true) + "\t" + ext.formDeci(betas[1], 6, true)
                + "\t" + ext.formDeci(stderrs[1], 6, true));
            // for (int j = 1; j < namesUsed.length; j++) {
            for (int j = 2; j < namesUsed.length; j++) {
              // writer.println(chrs[i]+"\t"+markerNames[i]+"\t"+positions[i]+"\t"+alleles[i][0]+"\t"+namesUsed[j]+"\t"+deps.length+"\t"+ext.formDeci(betas[1+j],
              // sigfig, true)+"\t"+ext.formDeci(stderrs[1+j], sigfig,
              // true)+"\t.\t.\t"+(stats[1+j]+"
              // ").substring(0,6).trim()+"\t"+ext.prettyP(pvals[1+j], sigfig, 4, 3, true));
              writer.println(chrs[i] + "\t" + markerNames[i] + "\t" + positions[i] + "\t"
                  + alleles[i][0] + "\t" + namesUsed[j] + "\t" + deps.length + "\t"
                  + ext.formDeci(betas[j], sigfig, true) + "\t"
                  + ext.formDeci(stderrs[j], sigfig, true) + "\t.\t.\t"
                  + (stats[j] + "     ").substring(0, 6).trim() + "\t"
                  + ext.prettyP(pvals[j], sigfig, 4, 3, true));
            }
          }
          writer.flush();

          // double[] finaldeps = model.getFinalDependentVariables();
          // double[][] finalindeps = model.getFinalIndependentVariables();
          // System.out.print("deps <- c(");
          // for (int j = 0; j < finaldeps.length; j++) {
          // System.out.print((j==0?"":",")+finaldeps[j]);
          // }
          // System.out.println(")");
          // for (int k = 0; k < finalindeps[0].length; k++) {
          // System.out.print("indeps"+(k+1)+" <- c(");
          // for (int j = 0; j < finalindeps.length; j++) {
          // System.out.print((j==0?"":",")+finalindeps[j][k]);
          // }
          // System.out.println(")");
          // }
          //
          // System.exit(1);
        }
      }
      writer.close();
      w2.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(phenoFile, false) + ".results."
          + (logistic ? "logistic" : "linear"));
      e.printStackTrace();
    }
  }

  public void exportToText(String outputFile) {

  }

  public byte[][] getGenotypeCounts() {
    return values;
  }


  public String[] getIds() {
    return ids;
  }

  public SnpMarkerSet getMarkerSet() {
    return markerSet;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  public void writeToPlinkFiles(String root) {
    PrintWriter writer;
    // private SnpMarkerSet markerSet;
    // private String[] ids;
    // private byte[][] values;

    char[][] alleles;
    String[] markerNames;
    byte[] chrs;
    int[] positions;
    String[][] annotation;

    markerNames = markerSet.getMarkerNames();
    chrs = markerSet.getChrs();
    positions = markerSet.getPositions();
    alleles = markerSet.getAlleles();
    annotation = markerSet.getAnnotation();

    try {
      writer = new PrintWriter(new FileWriter(root + ".ped"));
      for (int i = 0; i < ids.length; i++) {
        writer.print(ids[i] + "\t" + ids[i] + "\t0\t0\t1\t1");
        for (int j = 0; j < values.length; j++) {
          switch (values[j][i]) {
            case 2:
              writer.print("\t" + alleles[j][0] + "\t" + alleles[j][0]);
              break;
            case 1:
              writer.print("\t" + alleles[j][0] + "\t" + alleles[j][1]);
              break;
            case 0:
              writer.print("\t" + alleles[j][1] + "\t" + alleles[j][1]);
              break;
            default:
              writer.print("\t0\t0");
              break;
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + root + ".ped");
      e.printStackTrace();
    }

    try {
      writer = new PrintWriter(new FileWriter(root + ".map"));
      for (int i = 0; i < markerNames.length; i++) {
        writer.println(chrs[i] + "\t" + markerNames[i] + "\t0\t" + positions[i]);
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + root + ".map");
      e.printStackTrace();
    }

    CmdLine.run("plink --file " + root + " --make-bed --out " + root,
        ext.parseDirectoryOfFile(root));

    try {
      writer = new PrintWriter(new FileWriter(root + ".SetID"));
      for (int i = 0; i < markerNames.length; i++) {
        writer.println(annotation[i][0] + " " + markerNames[i]);
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + root + ".SetID");
      e.printStackTrace();
    }

  }

  // main
  // " \n" +
  // " Type Extension Description\n" +
  // " -1 [any] auto-detect from extension\n" +
  // " 0 .csv CHARGE-S Houston style .csv file\n" +
  // "";
}
