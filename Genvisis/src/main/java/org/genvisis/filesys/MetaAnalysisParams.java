package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class MetaAnalysisParams {
  public static final String DEFAULT_PARAMETERS = "filesys/default_meta_anlaysis.params";
  public static final String[] KEY_PARAMETERS = {"STUDIES", "PHENOTYPES", "RACES", "METHODS",
                                                 "GROUP_ANNOTATION_PARAMS"};
  public static final String[] KEY_PROPERTIES = {"SNP_INFO_FILE=", "VARIANT_NAME=", "CHROM_NAME=",
                                                 "GENE_NAME=", "FUNCTIONAL=", "R_EXEC=",
                                                 "RUN_BY_CHR="};

  private String[] studies;
  private String[] studyGroupings;
  private String[][] phenotypesWithFilenameAliases;
  private String[][] racesWithFilenameAliases;
  private final int[][] sampleSizes;
  private String snpInfoFilename;
  private String variantName;
  private String chromName;
  private String geneName;
  private String functionFlagName;
  private String rExec;
  private String[][] methods;
  private String[][] groupAnnotationParams;
  private boolean runByChr;

  private BufferedReader reader;
  private String trav;

  public MetaAnalysisParams(String filename, Logger log) {
    String[] line;
    Vector<String> v = new Vector<String>();
    boolean problem;

    problem = false;
    runByChr = true;

    if (!Files.exists(filename)) {
      log.report("File '" + filename
                 + "' does not exist; if you create an empty text file with this same filename, then it will be filled with example parameters and instructions");
      System.exit(1);
    }

    if (new File(filename).length() == 0) {
      log.report("File '" + filename
                 + "' is being populated with example parameters and instructions; tailor to your datasets and then re-run");
      Files.copyFileFromJar(DEFAULT_PARAMETERS, filename);
      System.exit(1);
    }

    try {
      reader = new BufferedReader(new FileReader(filename));
      nextIsParam();
      while (reader.ready()) {
        if (ext.indexOfStr(trav, KEY_PARAMETERS) >= 0) {
          if (trav.equals("STUDIES")) {
            v = populateParam();
            studies = new String[v.size()];
            studyGroupings = new String[v.size()];
            for (int i = 0; i < studies.length; i++) {
              line = v.elementAt(i).split("\t");
              if (line.length > 2) {
                log.reportError("Error - more than one tab for study " + line[0] + ", please fix");
                log.reportError("Found: " + v.elementAt(i));
                problem = true;
              }
              studies[i] = line[0];
              if (line.length > 1) {
                studyGroupings[i] = line[1];
              } else {
                studyGroupings[i] = "Final";
              }
            }
          } else if (trav.equals("PHENOTYPES")) {
            phenotypesWithFilenameAliases = populateParams();
          } else if (trav.equals("RACES")) {
            racesWithFilenameAliases = populateParams();
          } else if (trav.equals("METHODS")) {
            methods = populateParams();
            for (int i = 0; i < methods.length; i++) {
              if (methods[i].length < 3) {
                log.reportError("Error - a method must have at least 3 parameters: name, grouping, algorithm, (optional) MAF threshold, (optional) additional arguments such as weighting");
                log.reportError("Found: " + v.elementAt(i));
                problem = true;
              } else if (methods[i].length == 3) {
                log.reportError("Warning - no 4th token for method " + methods[i][0]
                                + "; all markers will be included in the analysis");
              } else if (!ext.isValidDouble(methods[i][3])
                         && !methods[i][2].equals("singlesnpMeta")) {
                log.reportError("Warning - no discernable MAF threshold token for method "
                                + methods[i][0] + " since '" + methods[i][3]
                                + "' is referring to something else; all markers will be included in the analysis");
              } else if (ext.isValidDouble(methods[i][3])
                         && methods[i][2].equals("singlesnpMeta")) {
                log.reportError("Error - MAF threshold token cannot be used for singlesnpMeta (as was done for '"
                                + methods[i][0] + "')");
                problem = true;
              }
            }
          } else if (trav.equals("GROUP_ANNOTATION_PARAMS")) {
            groupAnnotationParams = populateParams();
            for (String[] groupAnnotationParam : groupAnnotationParams) {
              if (groupAnnotationParam.length != 2) {
                log.reportError("Error - additional group annotation params must have exactly 2 tokens: a method grouping and a space separated GenParser parameter set");
                log.reportError("Found: " + Array.toStr(groupAnnotationParam));
                problem = true;
              }
            }
          }
        } else if (ext.indexOfStartsWith(trav, KEY_PROPERTIES, true) >= 0) {
          if (trav.startsWith("SNP_INFO_FILE=")) {
            snpInfoFilename =
                            ext.parseStringArg(trav,
                                               "default_snp_info_file_that_probably_does_not_exist.RData");
          } else if (trav.startsWith("VARIANT_NAME=")) {
            variantName = ext.parseStringArg(trav, "Name");
          } else if (trav.startsWith("CHROM_NAME=")) {
            chromName = ext.parseStringArg(trav, "CHROM");
          } else if (trav.startsWith("GENE_NAME=")) {
            geneName = ext.parseStringArg(trav, "SKATgene");
          } else if (trav.startsWith("FUNCTIONAL=")) {
            functionFlagName = ext.parseStringArg(trav, null);
          } else if (trav.startsWith("R_EXEC=")) {
            rExec = ext.parseStringArg(trav, null);
          } else if (trav.startsWith("RUN_BY_CHR=")) {
            runByChr = ext.parseBooleanArg(trav);
          } else {
            log.reportError("Error - property '" + trav
                            + "' was defined in MetaAnalysisParams.KEY_PROPERTIES but is not yet mapped to a variable name");
            problem = true;
          }
          nextIsParam();
        } else {
          log.reportError("Error - '" + trav + "' is an unknown parameter or property");
          problem = true;
        }

      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }

    sampleSizes = null;

    if (problem) {
      System.exit(1);
    }
  }

  public Vector<String> populateParam() {
    Vector<String> v;

    v = new Vector<String>();
    while (nextIsParam()) {
      v.add(trav);
    }

    return v;
  }

  public String[][] populateParams() {
    Vector<String> v;
    String[][] params;

    v = populateParam();
    params = new String[v.size()][];
    for (int i = 0; i < params.length; i++) {
      params[i] = v.elementAt(i).split("\t");
    }

    return params;
  }

  public boolean nextIsParam() {
    try {
      if (reader.ready()) {
        do {
          if (!reader.ready()) {
            trav = null;
            return false;
          }
          trav = reader.readLine().trim();
        } while (trav.equals("") || trav.startsWith("#"));

        return ext.indexOfStr(trav, KEY_PARAMETERS) == -1
               && ext.indexOfStartsWith(trav, KEY_PROPERTIES, true) == -1;
      } else {
        return false;
      }
    } catch (IOException ioe) {
      ioe.printStackTrace();
      return false;
    }


  }

  public String[] getStudies() {
    return studies;
  }

  public boolean runningByChr() {
    return runByChr;
  }

  public String[] getStudyGroupings() {
    return studyGroupings;
  }

  public String[][] getPhenotypesWithFilenameAliases() {
    return getPhenotypesWithFilenameAliases(true);
  }

  public String[][] getPhenotypesWithFilenameAliases(boolean pruneExclamations) {
    String[][] phenotypes;

    phenotypes = Matrix.clone(phenotypesWithFilenameAliases);
    for (int i = 0; i < phenotypes.length; i++) {
      if (pruneExclamations && phenotypes[i][0].startsWith("!")) {
        phenotypes[i][0] = phenotypes[i][0].substring(1);
      }
    }

    return phenotypes;
  }

  public String[][] getRacesWithFilenameAliases() {
    return racesWithFilenameAliases;
  }

  public int[][] getSampleSizes() {
    System.err.println("Error - sample sizes are not being imported as of yet");
    return sampleSizes;
  }

  public String getSnpInfoFilename() {
    return snpInfoFilename;
  }

  public String getVariantName() {
    return variantName;
  }

  public String getChromName() {
    return chromName;
  }

  public String getGeneName() {
    return geneName;
  }

  public String getFunctionFlagName() {
    return functionFlagName;
  }

  public String getRExec() {
    return rExec;
  }

  public String[][] getMethods() {
    return methods;
  }

  public String[] getGroups() {
    String[][] methods;
    String[] groups;

    methods = getMethods();

    groups = new String[] {};
    for (String[] method : methods) {
      if (ext.indexOfStr(method[1], groups) == -1) {
        groups = Array.addStrToArray(method[1], groups);
      }
    }

    return groups;
  }

  public String[][] getGroupAnnotationParams() {
    return groupAnnotationParams;
  }
}
