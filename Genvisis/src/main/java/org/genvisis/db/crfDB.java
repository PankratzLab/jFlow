package org.genvisis.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.GregorianCalendar;
import java.util.Hashtable;
import java.util.Vector;
import java.util.regex.Pattern;

import org.genvisis.common.Array;
import org.genvisis.common.Elision;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.stats.LeastSquares;
import org.genvisis.stats.LogisticRegression;
import org.genvisis.stats.RegressionModel;

public class crfDB {
  public static boolean SHOW_ABSENTS = false;
  public static final String DELIMITERS = "$@#\\-!\\+\\/\\*";
  public static final String POSS_DIR = org.genvisis.park.tools.CRF_DIR;
  public static final String[] DONT_COMP =
      {"UniqueID", "FamID", "IndID", "Site", "BirthDate", "OnsetDate", "DiagnosisDate", "ExamDate",
       "Dx", "CONF_PD", "Affected", "Gender", "BirthMonth", "BirthDay", "BirthYear", "Race",
       "OtherLRRK2", "APOE", "MMSEcutoff", "InitialSymptom", "Ethnicity", "Onset>20",
       "levodopaChorea", "PDopinion", "Alzheimers", "Remission", "lesionMRI"};
  public static final String[] DONT_COMP_DEP =
      {"Depression", "AdjDepression", "Depressed", "MajorDepression", "MinorDepression",
       "DSMIV_Dx_Depression", "depressionBefore", "depressionSince"};
  public static final double DEFAULT_MISSING_THRESHOLD = 0.25;

  private final Logger log;

  public crfDB(String dir, String keyFilename) throws Elision {
    this(dir, keyFilename, true);
  }

  public crfDB(String dir, String keyFilename, boolean logging) throws Elision { // dir is necessary
                                                                                 // for UpdateCRFdb
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line, subline;
    String temp, output_filename = "dump", unique_id;
    double score;
    Hashtable<String, String[][]> hash = new Hashtable<String, String[][]>();
    Hashtable<String, String> colHash;
    Vector<Hashtable<String, String>> hashes;
    String[] hashNames;
    String[][] expecteds = null, allData = null, filenameKey;
    double[][][] bounds;
    Vector<String> fileDescriptions = new Vector<String>(), idCollection = new Vector<String>(),
        hashDescriptions = new Vector<String>(), specials = new Vector<String>();
    String[] data;
    int numFiles, count, missing;
    char[] operator = null;
    int[][] indices;
    int[][] functionIndicies;
    boolean doubleFun, problem, done, fits;
    Pattern delimiters = Pattern.compile("[" + DELIMITERS + "]+");
    Vector<String> recodeFrom;
    String[] recodeTo;
    boolean[] dropFromFinal;
    int countDivZero = 0;
    double missingThreshold = DEFAULT_MISSING_THRESHOLD;
    int index;
    String delimiter;

    delimiter = ",";

    log = new Logger(logging ? dir + keyFilename + ".log" : null);

    try {
      problem = false;
      try {
        reader = new BufferedReader(new FileReader(dir + keyFilename));
        while (reader.ready()) {
          done = false;
          while (!done) {
            temp = reader.readLine();

            if (temp.toLowerCase().equals("tab")) {
              delimiter = "\t";
              temp = reader.readLine();
            }

            line = temp.split("[\\s]+");
            if (line.length == 1 && line[0].equals("")) {
              done = true;
            } else {
              if (!new File(line[0]).exists() && !new File(POSS_DIR + line[0]).exists()) {
                log.reportError("Error - File '" + temp.split("[\\s]+")[0] + "' not found.");
                problem = true;
              }
              for (int i = 2; i < line.length; i++) {
                for (int j = 0; j < DELIMITERS.length(); j++) {
                  if (line[i].startsWith(DELIMITERS.charAt(j) + "")) {
                    log.reportError("Error - a variable name cannot start with a reserved delimiter character ('"
                                    + line[i] + "' for " + line[0] + ")");
                    problem = true;
                  }
                }
              }
              if (line[1].equals("HASH")) {
                hashDescriptions.add(temp);
              } else {
                fileDescriptions.add(temp);
              }
            }
          }
          line = reader.readLine().split("[\\s]+");
          output_filename = line[0];
          for (int i = 1; i < line.length; i++) {
            if (line[i].startsWith("missingThreshold=")) {
              missingThreshold = Double.parseDouble(line[i].split("=")[1]);
            } else {
              log.reportError("Error - unknown additional argument: " + line[i]);
              reader.close();
              throw new Elision();
            }
          }

          while (reader.ready()) {
            temp = reader.readLine();
            if (!temp.startsWith("//")) {
              specials.add(temp);
            }
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        log.reportError("Error: file \"" + dir + keyFilename + "\" not found in current directory");
        problem = true;
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + dir + keyFilename + "\"");
        problem = true;
      }

      numFiles = fileDescriptions.size();
      if (problem) {
        throw new Elision();
      }

      hashes = HashVec.newVecHashStringString(hashDescriptions.size());
      hashNames = new String[hashDescriptions.size()];
      for (int i = 0; i < hashDescriptions.size(); i++) {
        line = (hashDescriptions.elementAt(i)).split("[\\s]+");
        hashNames[i] = line[0];
        temp = (new File(line[0]).exists() ? "" : POSS_DIR) + line[0];
        indices = new int[][] {{Integer.parseInt(line[2]), Integer.parseInt(line[3])}};
        try {
          reader = new BufferedReader(new FileReader(temp));
          while (reader.ready()) {
            line = reader.readLine().split(delimiter);
            hashes.elementAt(i).put(line[indices[0][0]], line[indices[0][1]]);
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          log.reportError("Error: file \"" + temp + "\" not found in current directory");
          throw new Elision();
        } catch (IOException ioe) {
          log.reportError("Error reading file \"" + temp + "\"");
          throw new Elision();
        }
      }

      expecteds = new String[numFiles + 1][];
      bounds = new double[numFiles + 1][][];
      expecteds[numFiles] = new String[specials.size()];
      bounds[numFiles] = new double[specials.size()][2];
      dropFromFinal = new boolean[specials.size()];
      filenameKey = new String[numFiles + 1][2];
      filenameKey[numFiles] = new String[] {"final", "key"};
      colHash = new Hashtable<String, String>();

      for (int i = 0; i < specials.size(); i++) {
        subline = (specials.elementAt(i)).split("=");
        if (subline.length != 2) {
          log.reportError("Error - Malformed special column: " + specials.elementAt(i));
          log.reportError("  proper example is SUM=+Q1+Q2-Q3+Q4 where SUM is the name and the value adds Q1, Q2 and Q4 and subtracts Q3");
          throw new Elision();
        }

        subline = subline[0].split("\\|");
        bounds[numFiles][i] = parseBounds(subline);

        if (subline[0].startsWith("#")) {
          subline[0] = subline[0].substring(1);
          dropFromFinal[i] = true;
        }
        if (colHash.containsKey(subline[0])) {
          log.reportError("Error - duplicate variable name (" + subline[0]
                          + ") for final database (detailed in " + keyFilename + ")");
          throw new Elision();
        } else {
          colHash.put(subline[0], "1");
        }
        if (delimiters.split(subline[0]).length > 1) {
          log.reportError("Error - Final variable names may not contain characters used as operators (i.e. "
                          + DELIMITERS + ")");
          log.reportError("        " + subline[0] + " has "
                          + (delimiters.split(subline[0]).length - 1) + " of them");
          throw new Elision();
        }

        expecteds[numFiles][i] = subline[0];
      }

      for (int fn = 0; fn < numFiles; fn++) {
        line = (fileDescriptions.elementAt(fn)).split("[\\s]+");
        expecteds[fn] = new String[line.length - 2];
        bounds[fn] = new double[line.length - 2][2];
        subline = line[1].split("=");
        if (!subline[0].equals("KEY")) {
          log.reportError("Error - " + line[1]
                          + " is in the second position but does not use proper key syntax");
          log.reportError("  example - if the unique id joins columns named FAM_NO and SUBJ_NO, use KEY=FAM_NO:SUBJ_NO");
          throw new Elision();
        }
        filenameKey[fn][0] = line[0];
        filenameKey[fn][1] = subline[1];
        for (int i = 2; i < line.length; i++) {
          subline = line[i].split("\\|");
          bounds[fn][i - 2] = parseBounds(subline);
          expecteds[fn][i - 2] = subline[0];
        }

        colHash = new Hashtable<String, String>();
        for (int i = 0; i < expecteds[fn].length; i++) {
          if (colHash.containsKey(expecteds[fn][i])) {
            log.reportError("Error - duplicate column header (" + expecteds[fn][i] + ") in file '"
                            + filenameKey[fn][0] + "'");
            throw new Elision();
          } else {
            colHash.put(expecteds[fn][i], "1");
          }
        }

        reader = new BufferedReader(new FileReader((new File(line[0]).exists() ? "" : POSS_DIR)
                                                   + filenameKey[fn][0]));
        subline = reader.readLine().split(delimiter);

        for (int i = 0; i < expecteds[fn].length; i++) {
          if (i >= subline.length) {
            log.reportError("Error - too many column headers in " + keyFilename + " for file "
                            + filenameKey[fn][0]);
            temp = "      - extra columns not in actual file:";
            for (int j = subline.length; j < expecteds[fn].length; j++) {
              temp += "\t" + expecteds[fn][j];
            }
            log.reportError(temp);
            reader.close();
            throw new Elision();
          }

          if (!expecteds[fn][i].equals(ext.replaceAllWith(subline[i], "-", "_"))) {
            log.reportError("Error - column headers do not match up with expected values for file "
                            + filenameKey[fn][0]);
            log.reportError("      - expecting " + expecteds[fn][i] + " in column " + (i + 1)
                            + ", instead found, " + subline[i]);
            reader.close();
            throw new Elision();
          }
        }

        if (subline.length > expecteds[fn].length) {
          log.reportError("Error - too few column headers in " + keyFilename + " for file "
                          + filenameKey[fn][0]);
          temp = "      - additional columns not anticipated:";
          for (int j = expecteds[fn].length; j < subline.length; j++) {
            temp += "\t" + subline[j];
          }
          log.reportError(temp);
          reader.close();
          throw new Elision();
        }

        while (reader.ready()) {
          temp = reader.readLine();
          if (temp.indexOf("\"") >= 0) {
            temp = fixMarks(temp);
          }
          while (temp.indexOf(",,") != -1) {
            temp = temp.replaceAll(",,", ",.,");
          }
          line = temp.split(delimiter);
          if (line.length != expecteds[fn].length) {
            if (line.length == expecteds[fn].length - 1 && temp.endsWith(delimiter)) {
              // log.reportError(" Warning - last column is empty for
              // individual in "+filenameKey[fn][0]+", treating as
              // missing data");
              // log.reportError(" "+temp);
              temp += ".";
              line = temp.split(delimiter);
            } else {
              log.reportError("Error - expecting " + expecteds[fn].length
                              + " columns of data in file '" + filenameKey[fn][0] + "', got "
                              + line.length + " some of the time.");
              log.reportError(temp);
              reader.close();
              throw new Elision();
            }
          }

          data = new String[expecteds[fn].length];
          for (int i = 0; i < expecteds[fn].length; i++) {
            data[i] = checkData(line[i]);
          }

          subline = filenameKey[fn][1].split(":");
          if (data[ext.indexOfStr(subline[0], expecteds[fn])].equals(".")) {
            System.err.println("Warning - there is an extra blank line in " + filenameKey[fn][0]);
            continue;
          }
          if (subline.length == 1) {
            // unique_id = filenameKey[fn][1];
            unique_id = data[ext.indexOfStr(subline[0], expecteds[fn])];
          } else {
            unique_id = data[ext.indexOfStr(subline[0], expecteds[fn])] + ext.formNum(
                                                                                      Integer.valueOf(data[ext.indexOfStr(subline[1],
                                                                                                                          expecteds[fn])])
                                                                                             .intValue(),
                                                                                      3);
          }
          if (fn == 0) {
            allData = new String[numFiles + 1][];
            allData[numFiles] = new String[specials.size()];
            hash.put(unique_id, allData);
            idCollection.add(unique_id);
          } else {
            if (hash.containsKey(unique_id)) {
              allData = hash.get(unique_id);
            } else {
              if (!filenameKey[fn][0].equals("ninfo1.csv")
                  && !filenameKey[fn][0].equals("ninfo2.csv")
                  && !filenameKey[fn][0].equals("ninfo5.csv")
                  && !filenameKey[fn][0].equals("ninfo6.csv")) { // these
                // have
                // additional
                // unrequired
                // indviduals
                log.reportError("Error - Unique ID '" + unique_id + "' was in file '"
                                + filenameKey[fn][0] + "' but not in previous file(s)");
              }
              allData = new String[numFiles + 1][];
              allData[numFiles] = new String[specials.size()];
              hash.put(unique_id, allData);
              idCollection.add(unique_id);
            }
          }

          for (int i = 0; i < expecteds[fn].length; i++) {
            checkBounds(filenameKey[fn][0] + ":" + expecteds[fn][i],
                        subline.length == 1 ? unique_id : formatUniqueID(unique_id), data[i],
                        bounds[fn][i][0], bounds[fn][i][1]);
          }

          allData[fn] = data;
        }
        reader.close();
      }

      int[] keys = Sort.quicksort(Array.toStringArray(idCollection));

      for (int i = 0; i < idCollection.size(); i++) {
        unique_id = idCollection.elementAt(keys[i]);
        allData = hash.get(unique_id);
        for (int fn = 0; fn < numFiles; fn++) {
          if (allData[fn] == null) {
            if (allData[0] != null && !filenameKey[fn][0].equals("G2019S.csv") && SHOW_ABSENTS) {
              log.reportError("Error - " + unique_id + " was absent from file '"
                              + filenameKey[fn][0] + "'");
            }
          }
        }
        for (int fn = 0; fn < numFiles; fn++) {
          if (allData[fn] == null) {
            allData[fn] = new String[expecteds[fn].length];
            for (int j = 0; j < allData[fn].length; j++) {
              allData[fn][j] = ".";
            }
          }
        }
      }

      for (int i = 0; i < specials.size(); i++) {
        subline = (specials.elementAt(i)).split("=");
        temp = subline[1];

        count = 0;
        for (int j = 0; j < DELIMITERS.length(); j++) {
          if (temp.startsWith(DELIMITERS.charAt(j) + "") && DELIMITERS.charAt(j) != '\\') {
            count++;
          }
        }

        if (count != 1) {
          log.reportError("Error - All expressions defining the final variables must start with an operator");
          throw new Elision();
        }

        if (temp.startsWith("$RESIDUAL") || temp.startsWith("$PREDICTED")) {
          subline = temp.split("[,]+");
          operator = new char[1];
        } else {
          subline = delimiters.split(temp);
          operator = new char[subline.length - 1];
        }
        indices = new int[subline.length - 1][2];

        count = 0;
        for (int j = 0; j < temp.length(); j++) {
          switch (temp.charAt(j)) {
            case '+':
              if (temp.charAt(j + 1) == '!') { // ergo, can't
                // subtract a
                // negative
                operator[count] = '!';
              } else {
                operator[count] = '+';
              }
              count++;
              break;
            case '-':
              if (temp.charAt(j + 1) == '!') { // ergo, can't
                // subtract a
                // negative
                log.reportError("Error - can't subtract the inverse of a number. Fix:");
                log.reportError("        " + temp);
                throw new Elision();
              }
              operator[count] = '-';
              count++;
              break;
            case '*':
              operator[count] = '*';
              count++;
              break;
            case '/':
              operator[count] = '/';
              count++;
              break;
            case '#':
              operator[count] = '#';
              count++;
              break;
            case '@':
              operator[count] = '@';
              count++;
              break;
            case '$':
              operator[count] = '$';
              count++;
              break;
          }
        }
        if (count > subline.length - 1) {
          log.reportError("Error - Malformed special equation: " + temp);
          log.reportError("        too many operators for not enough units");
        }

        for (int j = 0; j < indices.length; j++) {
          line = subline[j + 1].split("[,]+"); // +1 because the
          // first should
          // always be blank,
          // if expression
          // starts with an
          // operator
          for (String element : line) {
            if (element.indexOf(":") >= 0) {
              indices[j] = getVector(element, expecteds, filenameKey);
              if (indices[j][0] == numFiles && indices[j][1] >= i) {
                log.reportError("Error - cannot reference a variable (final:"
                                + expecteds[numFiles][indices[j][1]]
                                + ") before it is rendered. Reorder the final variables.");
                throw new Elision();
              }
            }
          }
        }

        if (temp.startsWith("$RESIDUAL") || temp.startsWith("$PREDICTED")) {
          boolean[] used = new boolean[idCollection.size()];
          Vector<String> deps = new Vector<String>();
          Vector<double[]> indeps = new Vector<double[]>();
          Vector<String> v = new Vector<String>();
          double[] dataline;
          RegressionModel model;
          boolean logistic;

          for (int j = 0; j < idCollection.size(); j++) {
            allData = hash.get(idCollection.elementAt(j));

            used[j] = true;

            for (int[] indice : indices) {
              if (allData[indice[0]][indice[1]].equals(".")) {
                used[j] = false;
              }
            }

            if (used[j]) {
              deps.add(allData[indices[0][0]][indices[0][1]]);
              dataline = new double[indices.length - 1];
              for (int k = 1; k < indices.length; k++) {
                dataline[k - 1] = Double.parseDouble(allData[indices[k][0]][indices[k][1]]);
              }
              indeps.add(dataline);
            }
          }

          count = 0;
          while (count < deps.size() && v.size() < 6) {
            if (!v.contains(deps.elementAt(count))) {
              v.add(deps.elementAt(count));
            }
            count++;
          }
          if (v.size() == 1) {
            log.reportError("Error - dependent variable has zero variation (" + v.elementAt(0)
                            + "); " + temp + " will not be calculated");
            used = new boolean[used.length];
          } else {
            logistic = (v.size() == 2);

            model =
                logistic ? new LogisticRegression(deps, indeps) : new LeastSquares(deps, indeps);
            dataline = temp.startsWith("$RESIDUAL") ? model.getResiduals() : model.getPredicteds();

            count = 0;
            for (int j = 0; j < idCollection.size(); j++) {
              allData = hash.get(idCollection.elementAt(j));
              if (used[j]) {
                allData[numFiles][i] = ext.formDeci(dataline[count++], 5);
              } else {
                allData[numFiles][i] = ".";
              }
            }
          }

        } else {
          doubleFun = false;
          for (int j = 0; j < idCollection.size(); j++) {
            allData = hash.get(idCollection.elementAt(j));
            if (indices.length == 1 && operator[0] == '#') {
              allData[numFiles][i] = allData[indices[0][0]][indices[0][1]];
            } else if (indices.length == 1 && operator[0] == '@') {
              allData[numFiles][i] = dateIt(idCollection.elementAt(j), subline[1],
                                            allData[indices[0][0]][indices[0][1]]);
            } else if (indices.length == 1 && operator[0] == '$'
                       && subline[1].startsWith("RECODE")) {
              line = subline[1].split("[,]+");
              if (line.length % 2 != 1) {
                log.reportError("Error - improper use of the RECODE function. Requires a variable name, a recoded value for each variable specified, as well as a final default value for non-missing but non-specified values");
                throw new Elision();
              }
              functionIndicies = new int[1][];
              functionIndicies[0] = getVector(line[1], expecteds, filenameKey);
              recodeFrom = new Vector<String>();
              recodeTo = new String[(line.length - 1) / 2];
              for (int l = 0; l < (line.length - 3) / 2; l++) {
                recodeFrom.add(line[2 + l * 2]);
                recodeTo[l + 1] = line[2 + l * 2 + 1];
              }
              recodeTo[0] = line[line.length - 1];
              if (!recodeFrom.contains(".")
                  && allData[functionIndicies[0][0]][functionIndicies[0][1]].equals(".")) {
                allData[numFiles][i] = ".";
              } else {
                allData[numFiles][i] =
                    recodeTo[recodeFrom.indexOf(allData[functionIndicies[0][0]][functionIndicies[0][1]])
                             + 1];
              }
            } else if (indices.length == 1 && operator[0] == '$' && subline[1].startsWith("HASH")) {
              line = subline[1].split("[,]+");
              if (line.length != 3) {
                log.reportError("Error - improper use of the HASH function. Requires a lookup filename (also listed at top), and the variable name that matches a key");
                throw new Elision();
              }
              functionIndicies = new int[1][];
              functionIndicies[0] = getVector(line[2], expecteds, filenameKey);
              index = ext.indexOfStr(line[1], hashNames);
              if (index == -1) {
                log.reportError("Error - HASH references a file that was not properly loaded");
                throw new Elision();
              }
              if (allData[functionIndicies[0][0]][functionIndicies[0][1]].equals(".")) {
                allData[numFiles][i] = ".";
              } else if (!hashes.elementAt(index)
                                .containsKey(allData[functionIndicies[0][0]][functionIndicies[0][1]])) {
                allData[numFiles][i] = "Hashtable did not have key '"
                                       + allData[functionIndicies[0][0]][functionIndicies[0][1]]
                                       + "'";
              } else {
                allData[numFiles][i] =
                    hashes.elementAt(index)
                          .get(allData[functionIndicies[0][0]][functionIndicies[0][1]]);
              }
            } else {
              score = 0;
              missing = 0;
              for (int k = 0; k < indices.length; k++) {
                if (checkData(allData[indices[k][0]][indices[k][1]]).equals(".")) {
                  missing++;
                }
                if (operator[k] != '$' && procDouble(allData[indices[k][0]][indices[k][1]])
                                          - procInt(allData[indices[k][0]][indices[k][1]]) > 0) {
                  doubleFun = true;
                }
                switch (operator[k]) {
                  case '+':
                    score += procDouble(allData[indices[k][0]][indices[k][1]]);
                    break;
                  case '-':
                    score -= procDouble(allData[indices[k][0]][indices[k][1]]);
                    break;
                  case '*':
                    score *= procDouble(allData[indices[k][0]][indices[k][1]]);
                    break;
                  case '/':
                    if (procDouble(allData[indices[k][0]][indices[k][1]]) == 0) {
                      if (countDivZero == 0) {
                        log.reportError("Error - Divide by zero");
                      }
                      score = 999;
                      countDivZero++;
                    } else {
                      score /= procDouble(allData[indices[k][0]][indices[k][1]]);
                    }
                    break;
                  case '!':
                    if (allData[indices[k][0]][indices[k][1]].equals("0")) {
                      score += 1;
                    } else if (procInt(allData[indices[k][0]][indices[k][1]]) != 0
                               && procInt(allData[indices[k][0]][indices[k][1]]) != 1) {
                      log.reportError("Error - '" + allData[indices[k][0]][indices[k][1]]
                                      + "' cannot be inversed for " + subline[k + 1]);
                    }
                    break;
                  case '$':
                    line = subline[k + 1].split("[,]+");
                    if (line[0].equals("EQUALS")) {
                      if (line.length != 3) {
                        log.reportError("Error - improper use of the EQUALS function. Requires two additional arguments: variable_name and target_value");
                        throw new Elision();
                      }
                      functionIndicies = new int[1][];
                      functionIndicies[0] = getVector(line[1], expecteds, filenameKey);
                      score +=
                          allData[functionIndicies[0][0]][functionIndicies[0][1]].equals(line[2]) ? 1
                                                                                                  : 0;
                    } else if (line[0].equals("AND")) {
                      if (line.length % 2 != 1) {
                        log.reportError("Error - improper use of the AND function. Requires an even number of additional arguments: (variable_name, target_value)*n");
                        throw new Elision();
                      }
                      fits = true;
                      functionIndicies = new int[1][];
                      for (int l = 0; l < (line.length - 1) / 2; l++) {
                        functionIndicies[0] = getVector(line[l * 2 + 1], expecteds, filenameKey);
                        if (!allData[functionIndicies[0][0]][functionIndicies[0][1]].equals(line[2])) {
                          fits = false;
                        }
                        if (allData[functionIndicies[0][0]][functionIndicies[0][1]].equals(".")) {
                          missing += line.length;
                        }
                      }
                      score += fits ? 1 : 0;
                    } else if (line[0].equals("OR")) {
                      if (line.length % 2 != 1) {
                        log.reportError("Error - improper use of the OR function. Requires an even number of additional arguments: (variable_name, target_value)*n");
                        throw new Elision();
                      }
                      fits = false;
                      functionIndicies = new int[1][];
                      for (int l = 0; l < (line.length - 1) / 2; l++) {
                        functionIndicies[0] = getVector(line[l * 2 + 1], expecteds, filenameKey);
                        if (allData[functionIndicies[0][0]][functionIndicies[0][1]].equals(line[2])) {
                          fits = true;
                        }
                        if (allData[functionIndicies[0][0]][functionIndicies[0][1]].equals(".")) {
                          missing += line.length;
                        }
                      }
                      score += fits ? 1 : 0;
                    } else if (line[0].equals("CONSTANT")) {
                      if (line.length != 2) {
                        log.reportError("Error - improper use of the CONSTANT function. Requires only the value of the constant");
                        throw new Elision();
                      }
                      score += procDouble(line[1]);
                    } else if (line[0].equals("GTE")) {
                      if (line.length != 3) {
                        log.reportError("Error - improper use of the GTE function. Requires two additional arguments: variable_name and cutoff_value");
                        throw new Elision();
                      }
                      functionIndicies = new int[2][];
                      functionIndicies[0] = getVector(line[1], expecteds, filenameKey);
                      if (line[2].indexOf(":") > 0) {
                        functionIndicies[1] = getVector(line[2], expecteds, filenameKey);
                        if (functionIndicies[1][0] == numFiles && functionIndicies[1][1] >= i) {
                          log.reportError("Error - cannot reference a variable (final:"
                                          + expecteds[numFiles][functionIndicies[1][1]]
                                          + ") before it is rendered. Reorder the final variables.");
                          throw new Elision();
                        }
                      }
                      if (allData[functionIndicies[0][0]][functionIndicies[0][1]].equals(".")) {
                        missing++;
                      } else if (functionIndicies[1] != null
                                 && allData[functionIndicies[1][0]][functionIndicies[1][1]].equals(".")) {
                        missing++;
                      } else {
                        score +=
                            procDouble(allData[functionIndicies[0][0]][functionIndicies[0][1]]) >= procDouble(line[2].indexOf(":") > 0 ? allData[functionIndicies[1][0]][functionIndicies[1][1]]
                                                                                                                                       : line[2]) ? 1
                                                                                                                                                  : 0;
                      }
                    } else if (line[0].equals("LTE")) {
                      if (line.length != 3) {
                        log.reportError("Error - improper use of the LTE function. Requires two additional arguments: variable_name and cutoff_value");
                        throw new Elision();
                      }
                      functionIndicies = new int[2][];
                      functionIndicies[0] = getVector(line[1], expecteds, filenameKey);
                      if (line[2].indexOf(":") > 0) {
                        functionIndicies[1] = getVector(line[2], expecteds, filenameKey);
                        if (functionIndicies[1][0] == numFiles && functionIndicies[1][1] >= i) {
                          log.reportError("Error - cannot reference a variable (final:"
                                          + expecteds[numFiles][functionIndicies[1][1]]
                                          + ") before it is rendered. Reorder the final variables.");
                          throw new Elision();
                        }
                      }
                      if (allData[functionIndicies[0][0]][functionIndicies[0][1]].equals(".")) {
                        missing++;
                      } else if (functionIndicies[1] != null
                                 && allData[functionIndicies[1][0]][functionIndicies[1][1]].equals(".")) {
                        missing++;
                      } else {
                        score +=
                            procDouble(allData[functionIndicies[0][0]][functionIndicies[0][1]]) <= procDouble(line[2].indexOf(":") > 0 ? allData[functionIndicies[1][0]][functionIndicies[1][1]]
                                                                                                                                       : line[2]) ? 1
                                                                                                                                                  : 0;
                      }
                    } else if (line[0].equals("FLOOR")) {
                      if (line.length != 2) {
                        log.reportError("Error - improper use of the FLOOR function. Requires only one additional argument: variable_name");
                        throw new Elision();
                      }
                      functionIndicies = new int[1][];
                      functionIndicies[0] = getVector(line[1], expecteds, filenameKey);
                      if (allData[functionIndicies[0][0]][functionIndicies[0][1]].equals(".")) {
                        missing++;
                      } else {
                        score +=
                            procDouble(allData[functionIndicies[0][0]][functionIndicies[0][1]].indexOf(".") > 0 ? allData[functionIndicies[0][0]][functionIndicies[0][1]].substring(0,
                                                                                                                                                                                    allData[functionIndicies[0][0]][functionIndicies[0][1]].indexOf("."))
                                                                                                                : allData[functionIndicies[0][0]][functionIndicies[0][1]]);
                      }
                    } else {
                      log.reportError("Error - Unrecognized function header: " + line[0]);
                      throw new Elision();
                    }
                    break;
                  default:
                }
              }
              if (missingThreshold == 0 ? missing == indices.length
                                        : (double) missing
                                          / (double) indices.length >= missingThreshold) {
                allData[numFiles][i] = ".";
              } else {
                allData[numFiles][i] =
                    (doubleFun ? ext.formDeci(score, 2, true) : ext.formDeci(score, 0));
              }
            }
            checkBounds(filenameKey[numFiles][0] + ":" + expecteds[numFiles][i],
                        idCollection.elementAt(j).length() == 8
                                                                ? formatUniqueID(idCollection.elementAt(j)) : idCollection.elementAt(j),
                        allData[numFiles][i], bounds[numFiles][i][0], bounds[numFiles][i][1]);
          }
        }

      }

      if (countDivZero > 0) {
        log.reportError("Divided by zero " + countDivZero + " times");
      }

      writer = new PrintWriter(new FileWriter(dir + output_filename));
      writer.print("UniqueID" + (filenameKey[0][1].split(":").length == 2 ? "\tFamID\tIndID" : ""));
      temp = "";
      for (int i = 0; i < specials.size(); i++) {
        writer.print(dropFromFinal[i] ? "" : "\t" + expecteds[numFiles][i]);
        temp +=
            dropFromFinal[i] ? ""
                             : (ext.indexOfStr(expecteds[numFiles][i], DONT_COMP) == -1
                                && ext.indexOfStr(expecteds[numFiles][i],
                                                  DONT_COMP_DEP) == -1 ? "\t"
                                                                         + expecteds[numFiles][i]
                                                                       : "");
      }
      writer.println();
      temp = temp.length() > 0 ? temp.substring(1) : temp;
      log.reportError(temp);

      for (int i = 0; i < idCollection.size(); i++) {
        unique_id = idCollection.elementAt(keys[i]);
        writer.print(unique_id + (filenameKey[0][1].split(":").length == 2
                                                                           ? "\t"
                                                                             + unique_id.substring(0,
                                                                                                   5)
                                                                             + "\t"
                                                                             + Integer.valueOf(unique_id.substring(5))
                                                                                      .intValue()
                                                                           : ""));
        allData = hash.get(idCollection.elementAt(keys[i]));
        for (int j = 0; j < allData[numFiles].length; j++) {
          writer.print(dropFromFinal[j] ? "" : "\t" + allData[numFiles][j]);
        }
        writer.println();
      }
      writer.close();

      writer = new PrintWriter(new FileWriter(dir + "example.ctl"));
      writer.println("\"" + new File(output_filename).getAbsolutePath() + "\" Affected=1 VPD=1");
      writer.println("AOO");
      writer.println(temp);
      writer.println();
      writer.println();
      writer.close();
    } catch (Elision e) {
      throw e;
    } catch (Exception e) {
      log.reportError(e.getMessage());
      log.reportException(e);
    }
  }

  public double[] parseBounds(String[] subline) {
    double[] bounds = new double[2];

    if (subline.length == 1) {
      bounds[0] = 0;
      bounds[1] = -1;
    } else if (subline.length == 2) {
      bounds[1] = Double.parseDouble(subline[1]);
      if (bounds[1] < 0) {
        log.reportError("Warning - upper bound for variable " + subline[0]
                        + " is below the default lower (zero); define new lower or else all will be flagged as invalid");
      }
    } else if (subline.length == 3) {
      bounds[0] = Double.parseDouble(subline[1]);
      bounds[1] = Double.parseDouble(subline[2]);
      if (bounds[1] < bounds[0]) {
        log.reportError("Error - upper bound " + bounds[1] + " for variable " + subline[0]
                        + " is less than the lower bound " + bounds[0]
                        + "; all will be flagged as invalid");
      }
    } else {
      log.reportError("Error - Don't know what to do with " + subline.length
                      + " arguments for variable '" + subline[0] + "'");
    }

    return bounds;
  }

  public void checkBounds(String variable, String id, String value, double lower, double upper) {
    if (!value.equals(".") && lower < upper) {
      try {
        if (Double.parseDouble(value) < lower) {
          log.reportError("Error - " + variable + " for " + id + " was below the lower bound ("
                          + lower + "): " + value);
        }
        if (Double.parseDouble(value) > upper) {
          log.reportError("Error - " + variable + " for " + id + " was above the upper bound ("
                          + upper + "): " + value);
        }
      } catch (NumberFormatException nfe) {
        log.reportError("Error - " + variable + " for " + id + " was not a number: " + value);
      }
    }

  }

  public static void dumpData(String dir, String id, String[][] data) throws IOException {
    PrintWriter writer = new PrintWriter(new FileWriter(dir + id + ".xln"));
    for (String[] element : data) {
      if (element != null) {
        for (int j = 0; j < element.length; j++) {
          writer.print(element[j] + "\t");
        }
      }
      writer.println();
    }
    writer.close();
  }

  public String formatUniqueID(String str) {
    return str.substring(0, 5) + "-" + Integer.valueOf(str.substring(5)).intValue();
  }

  public String dateIt(String unique_id, String variable, String date) {
    int month = 0, day = 0, year = 0;
    String[] errorChars = new String[] {"O", " ", "K"};

    if (date.equals("NA")) {
      date = "NNNN";
    }

    if (unique_id.length() == 8) {
      unique_id = formatUniqueID(unique_id);
    }

    if (date.equals(".")) {
      return ".";
    }

    if (date.length() == 4) {
      log.reportError("Error - " + variable + " for " + unique_id
                      + " is lacking missing data characters: '" + date + "'; should be 'UUUU"
                      + date + "'");
      date = "UUUU" + date;
    }
    if (date.length() == 6) {
      log.reportError("Error - " + variable + " for " + unique_id
                      + " is lacking missing data characters: '" + date + "'; should be 'UU" + date
                      + "'");
      date = "UU" + date;
    }
    if (date.length() != 8) {
      day = 0;
      for (int i = 0; i < date.length(); i++) {
        if (date.charAt(i) == 'N' || date.charAt(i) == 'U') {
          day++;
        }
      }
      if (date.length() == day) {
        log.reportError("Error - " + variable + " for " + unique_id
                        + " has the wrong number of characters: '" + date
                        + "'; should be 8 characters, not " + date.length());
        return ".";
      }

      log.reportError("Error - " + variable + " for " + unique_id
                      + " has a malformed date and is currently considered missing data: '" + date
                      + "'");
      return ".";
    }
    for (int i = 0; i < errorChars.length; i++) {
      if (date.indexOf(errorChars[i]) != -1) {
        log.reportError("Error - " + variable + " for " + unique_id + " has an invalid character: '"
                        + date + "'" + " -> replace letter '" + errorChars[i] + "' with "
                        + (i == 0 ? "number '0'" : "letter 'U'"));
        date = date.replaceAll(errorChars[i], (i == 0 ? "0" : "U"));
      }
    }
    date = date.replaceAll("N", "0");
    date = date.replaceAll("U", "0");
    date = date.replaceAll("A", "0"); // Ind 70730-9
    try {
      year = Integer.valueOf(date.substring(4)).intValue();
      if (year == 0) {
        return ".";
      }
      month = Integer.valueOf(date.substring(0, 2)).intValue();
      day = Integer.valueOf(date.substring(2, 4)).intValue();
    } catch (NumberFormatException nfe) {
      log.reportError("Error - " + variable + " for " + unique_id
                      + " has a malformed data and is currently considered missing data: '" + date
                      + "'");
      return ".";
    }
    if (month == 0) {
      month = 7;
      day = 4;
    } else if (day == 0) {
      day = 15;
    }
    month--;
    day--;

    return ext.formDeci(((double) ((new GregorianCalendar(year, month, day).getTimeInMillis()
                                    - new GregorianCalendar(year, 0, 0).getTimeInMillis())
                                   / 1000 / 60 / 60 / 24))
                        / 365 + year, 2, true);
  }

  public int[] getVector(String str, String[][] list, String[][] filenames) {
    int[] vec = new int[2];
    String[] subline = str.split(":");

    if (subline.length != 2) {
      log.reportError("Error - Variable needs to have 2 colon-delimited parts, filename:variable_name");
      log.reportError("        " + str + " has " + subline.length + " parts");
      return null;
    }

    vec[0] = -1;
    for (int i = 0; i < filenames.length; i++) {
      if (subline[0].equals(filenames[i][0])) {
        vec[0] = i;
      }
    }
    if (vec[0] == -1) {
      log.reportError("Error - malformed variable: " + subline[0]
                      + " is not a filename delineated at the top of the file");
      return null;
    }

    vec[1] = -1;
    for (int i = 0; i < list[vec[0]].length; i++) {
      if (subline[1].equals(list[vec[0]][i])) {
        vec[1] = i;
      }
    }

    if (vec[1] == -1) {
      log.reportError("Error - Could not find the variable '" + str + "'");
      return null;
    }

    return vec;
  }

  public static String fixMarks(String str) {
    String line = "";
    int on = 0;

    str = str.replaceAll("\"\"", "'");

    line += str.charAt(0);
    for (int i = 1; i < str.length(); i++) {
      if (str.charAt(i - 1) == ',' && str.charAt(i) == '\"' && on != 1) {
        on = 1;
      } else if (str.charAt(i) == '\"' && str.length() > (i + 1) && str.charAt(i + 1) == ','
                 && on == 1) {
        on = 2;
      } else if (str.charAt(i) == ',') {
        if (on == 1) {
          line += ';';
        } else {
          line += ',';
        }
      } else {
        line += str.charAt(i);
      }
    }
    if (on == 1) {
      return str;
    } else {
      return line;
    }
  }

  public double procDouble(String str) {
    double value = 0;

    if (checkData(str).equals(".") || checkData(str).equals("NA")) {
      value = 0;
    } else if (checkData(str).equals("Y")) {
      value = 1;
    } else {
      try {
        value = Double.parseDouble(str);
      } catch (Exception e) {
        log.reportError("Error - '" + str + "' is not a valid 'double' number");
      }
    }

    return value;
  }

  public int procInt(String str) {
    double value = 0;

    if (checkData(str).equals(".") || checkData(str).equals("NA")) {
      value = 0;
    } else if (checkData(str).equals("Y")) {
      value = 1;
    } else {
      try {
        value = Double.valueOf(str).doubleValue();
      } catch (Exception e) {
        log.reportError("Error - '" + str + "' is not a valid integer");
      }
    }

    return (int) value;
  }

  public static String checkData(String str) {
    if (str.equals(".") || str.equals("N") || str.equals("U") || str.equals("N.N")
        || str.equals("U.U") || str.equals("NN") || str.equals("UU") || str.equals("NN.NN")
        || str.equals("UU.UU") || str.equals("NNNN") || str.equals("UUUU") || str.equals("UUU")
        || str.equals("NNN") || str.equals("N.A") || str.equals("NA")) {
      return ".";
    }
    if (str.equals(". 0")) {
      return "0";
    }
    if (str.endsWith(".UU") || str.endsWith(".NN")) {
      return str.substring(0, str.lastIndexOf(".")) + ".00";
    }
    return str;
  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "crf_db_key.crf";

    boolean logging = true;
    boolean suppress = false;

    String usage = "\n" + "park.crfDB requires 0-1 arguments\n" + "   (1) filename (i.e. "
                   + filename + " (default)\n" + "   (2) log errors (i.e. log=" + logging
                   + " (default)\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("log=")) {
        logging = arg.split("=")[1].equalsIgnoreCase("true");
        numArgs--;
      } else if (arg.startsWith("-suppress")) {
        suppress = true;
        numArgs--;
      } else {
        filename = arg;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      new crfDB("", filename, logging);
    } catch (Exception e) {
      e.printStackTrace();
    }
    if (!suppress) {
      ext.waitForResponse("Press any key to close this window");
    }
  }
}
