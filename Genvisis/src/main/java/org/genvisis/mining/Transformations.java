package org.genvisis.mining;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.common.Array;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

import com.google.common.primitives.Doubles;

public class Transformations {
  private static final String[] LABELS =
      {"Identity", "Rank", "Log", "Inverse", "Square root", "Squared", "Cubed", "Box-Cox (maxLL)",
          "Box-Cox (minKurt)", "Normalized", "Standardized", "InverseNormalized", "NegativeLog10"};
  public static final int IDENTITY = 0;
  public static final int RANK = 1;
  public static final int LOG_NATURAL = 2;
  public static final int INVERSE = 3;
  public static final int SQUARE_ROOT = 4;
  public static final int SQUARED = 5;
  public static final int CUBED = 6;
  public static final int BOXCOX_LL = 7;
  public static final int BOXCOX_KURT = 8;
  public static final int NUM_TRANSFORMATIONS = 9;
  public static final int NORMALIZE = 9;
  public static final int STANDARDIZE_RANGE = 10;
  public static final int INVERSE_NORMALIZE = 11;
  public static final int QUANTILE = 12;
  public static final int INVERSE_TDIST_5DF = 13;
  public static final int NEGATIVE_LOG10 = 14;
  public static final int LOG10 = 15;
  public static final int MAD_SCALED = 16;
  public static final int X3 = 17;
  public static final int X5 = 18;

  public static double[] cubedTransform(double[] array) {
    double[] trans = new double[array.length];

    for (int i = 0; i < array.length; i++) {
      trans[i] = Math.pow(array[i], 3);
    }

    return trans;
  }

  public static void fromParameters(String filename, Logger log) {
    Vector<String> paramV;
    String[] defaults;
    String types;

    types = "";
    for (int i = 0; i < LABELS.length; i++) {
      types += " " + i + "=" + LABELS[i];
    }
    defaults = new String[] {"file=input.txt", "out=output.txt", "ignoreFirstLine=false", "col=0",
        "comma=false", "replace=false", "type=11", "#possible types:" + types};
    for (int i = 0; i < LABELS.length; i++) {
      Array.addStrToArray("# " + i + "=" + LABELS[i], defaults);
    }
    paramV = Files.parseControlFile(filename, "transform", defaults, log);
    if (paramV != null) {
      paramV.addElement("logfile=" + log.getFilename());
      main(Array.toStringArray(paramV));
    }
  }

  public static String getLabel(int type) {
    int max = 0;
    for (String element : LABELS) {
      if (element.length() > max) {
        max = element.length();
      }
    }
    if (type == -1) {
      return ext.formStr("", max, true);
    } else {
      return ext.formStr(LABELS[type], max, true);
    }
  }

  public static double[] inverseTransform(double[] array) {
    double[] trans = new double[array.length];

    for (int i = 0; i < array.length; i++) {
      trans[i] = 1 / array[i];
    }

    return trans;
  }

  public static double[] log10Transform(double[] array) {
    double[] trans = new double[array.length];

    for (int i = 0; i < array.length; i++) {
      trans[i] = Math.log10(array[i]);
    }

    return trans;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "input.dat";
    String outfile = null;
    boolean ignoreFirstLine = false;
    int col = 0;
    boolean commaDelimited = false;
    boolean replace = false;
    int type = NORMALIZE;
    String logfile = null;

    String usage = "\n" + "mining.Transformations requires 0-1 arguments\n"
        + "   (1) filename (i.e. file=" + filename + " (default))\n" + "   (2) outfile (i.e. out="
        + filename + "." + LABELS[type] + " (default))\n"
        + "   (3) ignore first [header] line (i.e. ignoreFirstLine=" + ignoreFirstLine
        + " (default))\n" + "   (4) column to be transformed (i.e. col=" + col + " (default))\n"
        + "   (5) comma delimited (i.e. comma=" + commaDelimited + " (default))\n"
        + "   (6) replace previous variables (i.e. replace=" + replace + " (default))\n"
        + "   (7) type of transformation (see below for options) (i.e. type=" + type
        + " (default))\n" + "";
    for (int i = 0; i < LABELS.length; i++) {
      usage += "       " + i + "=" + LABELS[i];
    }

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        outfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("ignoreFirstLine=")) {
        ignoreFirstLine = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("col=")) {
        col = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("comma=")) {
        commaDelimited = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("replace=")) {
        replace = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("type=")) {
        type = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("logfile=")) {
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
      // filename = "unr_phen.1";
      // commaDelimited = true;
      // col = 4;
      // replace = false;
      // ignoreFirstLine = true;
      // type = Transformations.INVERSE_NORMALIZE;
      //
      if (outfile == null) {
        outfile = filename + "." + LABELS[type];
      }
      transformFile(filename, outfile, ignoreFirstLine, col, commaDelimited, replace, type,
          new Logger(logfile));
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static double[] naturalLogTransform(double[] array) {
    double[] trans = new double[array.length];

    for (int i = 0; i < array.length; i++) {
      trans[i] = Math.log(array[i]);
    }

    return trans;
  }

  public static double[] negativeLog10Transform(double[] array) {
    double[] trans = new double[array.length];

    for (int i = 0; i < array.length; i++) {
      trans[i] = -1 * Math.log10(array[i]);
    }

    return trans;
  }

  public static double[] percentileTransform(double[] array) {
    double[] ranks = rankTransform(array);
    double[] newData = new double[array.length];

    for (int i = 0; i < newData.length; i++) {
      newData[i] = ranks[i] / array.length * 100;
    }

    return newData;
  }

  public static double[] rankTransform(double[] array) {
    int[] order = Sort.quicksort(array);
    double[] newData = new double[array.length];
    int count = 0, plus;

    while (count < array.length) {
      plus = 0;
      while (count + plus + 1 < array.length
          && array[order[count]] == array[order[count + plus + 1]]) {
        plus++;
      }
      for (int i = 0; i < plus + 1; i++) {
        newData[order[count + i]] = ((double) count * 2 + plus) / 2 + 1;
      }
      count += plus + 1;
    }

    return newData;
  }

  public static double[] sqrtTransform(double[] array) {
    double[] trans = new double[array.length];

    for (int i = 0; i < array.length; i++) {
      trans[i] = Math.sqrt(array[i]);
    }

    return trans;
  }

  public static double[] squaredTransform(double[] array) {
    double[] trans = new double[array.length];

    for (int i = 0; i < array.length; i++) {
      trans[i] = Math.pow(array[i], 2);
    }

    return trans;
  }

  public static double[] standardizeRange(double[] data) {
    return standardizeRange(data, 0, 1);
  }

  public static double[] standardizeRange(double[] data, double finalMin, double finalMax) {
    double[] newData = new double[data.length];
    int[] keys = Sort.quicksort(data);
    double min = data[keys[0]];
    double max = data[keys[data.length - 1]];

    for (int i = 0; i < newData.length; i++) {
      newData[i] = (data[i] - min) / (max - min) * (finalMax - finalMin) + finalMin;
    }

    return newData;
  }

  public static double[] transform(double[] array, int type) {
    return transform(array, type, new Logger());
  }

  public static double[] transform(double[] array, int type, Logger log) {
    switch (type) {
      case IDENTITY:
        return array;
      case RANK:
        return rankTransform(array);
      case NORMALIZE:
        return Array.normalize(array);
      case INVERSE_NORMALIZE:
        return Array.inverseNormalize(array);
      case QUANTILE:
        return Array.quantiles(array);
      case INVERSE_TDIST_5DF:
        return Array.inverseTdist(array, 5);
      case STANDARDIZE_RANGE:
        return standardizeRange(array);
      case LOG_NATURAL:
        return naturalLogTransform(array);
      case NEGATIVE_LOG10:
        return negativeLog10Transform(array);
      case LOG10:
        return log10Transform(array);
      case INVERSE:
        return inverseTransform(array);
      case SQUARE_ROOT:
        return sqrtTransform(array);
      case SQUARED:
        return squaredTransform(array);
      case CUBED:
        return cubedTransform(array);
      case BOXCOX_LL:
        return new BoxCox(array, log).getTransform_MaxLL();
      case BOXCOX_KURT:
        return new BoxCox(array, log).getTransform_MinKurt();
      case MAD_SCALED:
        BeastScore beastScore = new BeastScore(Array.toFloatArray(array), null, null, log);
        return Array.toDoubleArray(
            beastScore.getinverseTransformedDataScaleMAD(BeastScore.SCALE_FACTOR_MAD));
      case X3:
        return Array.multiply(array, 3);
      case X5:
        return Array.multiply(array, 5);

      default:
        log.reportError(
            "Error - '" + type + "' does not map to an implemented method; using NORMALIZE");
        return Array.normalize(array);
    }
  }

  public static double[][] transform(double[][] data, int type) {
    double[] array;

    for (int i = 0; i < data[0].length; i++) {
      array = new double[data.length];
      for (int j = 0; j < data.length; j++) {
        array[j] = data[j][i];
      }
      array = transform(array, type, new Logger());
      for (int j = 0; j < data.length; j++) {
        data[j][i] = array[j];
      }
    }

    return data;
  }

  public static float[] transform(float[] array, int type) {
    return Array.toFloatArray(transform(Array.toDoubleArray(array), type, new Logger()));
  }

  public static void transformFile(String filename, String outfile, boolean ignoreFirstLine,
      int column, boolean commaDelimited, boolean replace, int type, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, transformed;
    DoubleVector dv;
    IntVector duds;
    int count;

    line = HashVec.loadFileToStringArray(filename, false, ignoreFirstLine, new int[] {column}, true,
        false, commaDelimited ? "," : "[\\s]+");
    dv = new DoubleVector();
    duds = new IntVector();
    for (int i = 0; i < line.length; i++) {
      if (ext.isValidDouble(line[i])) {
        dv.add(Double.parseDouble(line[i]));
      } else {
        duds.add(i);
      }
    }
    transformed = Array.toStringArray(transform(Doubles.toArray(dv), type, new Logger()));
    for (int i = 0; i < duds.size(); i++) {
      Array.addStrToArray(".", transformed, duds.elementAt(i));
    }

    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = new PrintWriter(new FileWriter(outfile));
      if (ignoreFirstLine) {
        line = reader.readLine().trim().split(commaDelimited ? "," : "[\\s]+");
        if (!replace) {
          line = Array.insertStringAt(line[column] + "_" + Transformations.LABELS[type], line,
              column + 1);
        }
        writer.println(Array.toStr(line, commaDelimited ? "," : "\t"));
      }
      count = 0;
      while (reader.ready()) {
        line = reader.readLine().trim().split(commaDelimited ? "," : "[\\s]+");
        if (replace) {
          line[column] = transformed[count];
        } else {
          line = Array.insertStringAt(transformed[count], line, column + 1);
        }
        writer.println(Array.toStr(line, commaDelimited ? "," : "\t"));
        count++;
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

}
