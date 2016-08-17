package org.genvisis.parse;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Matrix;
import org.genvisis.common.Zip;
import org.genvisis.common.ext;
import org.genvisis.mining.Transformations;
import org.genvisis.stats.LeastSquares;

public class GEO_expression {
  public static final String DIR = "";

  public static void average(String filename) {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line, header1, header2, keys;
    Hashtable<String, Vector<String[]>> hashV;
    String temp;
    Vector<String[]> v;
    double sum;

    hashV = new Hashtable<String, Vector<String[]>>();
    try {
      reader = new BufferedReader(new FileReader(filename));
      header1 = reader.readLine().trim().split("[\\s]+");
      header2 = reader.readLine().trim().split("[\\s]+");
      while (reader.ready()) {
        temp = reader.readLine();
        line = temp.trim().split("[\\s]+");
        if (hashV.containsKey(line[0])) {
          v = hashV.get(line[0]);
        } else {
          hashV.put(line[0], v = new Vector<String[]>());
        }
        v.add(line);
      }
      reader.close();

      writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + "_avg.xln"));
      header1[1] = "N";
      writer.println(Array.toStr(header1));
      header2[1] = "N";
      writer.println(Array.toStr(header2));
      keys = HashVec.getKeys(hashV);
      for (String key : keys) {
        v = hashV.get(key);
        writer.print(key + "\t" + v.size());
        for (int j = 2; j < header1.length; j++) {
          sum = 0;
          for (int k = 0; k < v.size(); k++) {
            sum += Double.parseDouble(v.elementAt(k)[j]);
          }
          writer.print("\t" + sum / v.size());
        }
        writer.println();
      }
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  public static void checkFileContents(String dir) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, files, contents;
    String dest;

    files = Files.list(dir, ".gz", false);
    dest = dir + "unzipped/";
    new File(dest).mkdirs();
    try {
      writer = new PrintWriter(new FileWriter(dir + "list.xln"));
      for (String file : files) {
        Zip.unzipFile(dir + file, dir);
        contents = Files.list(dir, ".txt", false);
        writer.print(file);
        for (String content : contents) {
          try {
            reader = new BufferedReader(new FileReader(dir + content));
            line = reader.readLine().trim().split("[\\s]+");
            for (String element : line) {
              if (element.startsWith("NA")) {
                writer.print("\t" + element);
              }
            }
            reader.close();
          } catch (FileNotFoundException fnfe) {
            System.err
                .println("Error: file \"" + dir + content + "\" not found in current directory");
            System.exit(1);
          } catch (IOException ioe) {
            System.err.println("Error reading file \"" + dir + content + "\"");
            System.exit(2);
          }
          new File(dir + content).renameTo(new File(dest + content));

        }
        writer.println();
        writer.flush();
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + "list.xln");
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String filename = "GAK_region.dat";
    // String filename = "random1000b.txt";
    // String filename = "GAK_region_missing_meaned.dat";
    // String filename = "kids_random1000b.txt";
    String filename = "KARLN_TRIO.txt";


    String usage = "\n" + "parse.GEO_expression requires 0-1 arguments\n"
        + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Myron\\eQTLs\\GEO
    // DataSet GSE6536 new\\unzipped\\";
    // filename = dir + filename;
    try {
      parse(filename);
      // String dir = "C:\\Documents and Settings\\npankrat\\My
      // Documents\\gwas\\followUp\\GAK\\eQTL\\Duke\\";
      // filename = dir + filename;
      if (new File(ext.rootOf(filename, false) + "_transcripts.xln").exists()) {
        regress(ext.rootOf(filename, false) + "_transcripts.xln");
      } else {
        System.out.println("Error - file '" + ext.rootOf(filename, false) + "_transcripts.xln"
            + "' not found; nothing to parse");
      }
      average(ext.rootOf(filename, false) + "_transcripts_residuals.xln");
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parse(String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, files, data;
    String trav;
    String[] targets, names;
    int index, count;
    Hashtable<String, String> hash;

    files = Files.list(".", ".txt", false);
    targets = HashVec.loadFileToStringArray(filename, false, new int[] {0}, false);
    try {
      names = HashVec.loadFileToStringArray(filename, false, new int[] {1}, false);
    } catch (ArrayIndexOutOfBoundsException aioobe) {
      System.out.println(
          "No gene names in column 2 of the file (useful for generating easy to understand PLINK pheno files); copying transcript information...");
      names = targets;
    }

    try {
      writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + "_transcripts.xln"));
      writer.println("Sample\tfile\tbackground\t" + Array.toStr(targets));
      writer.println("Sample\tfile\tbackground\t" + Array.toStr(names));
      hash = new Hashtable<String, String>();
      count = 0;
      for (int i = 0; i < files.length; i++) {
        try {
          reader = new BufferedReader(new FileReader(files[i]));
          trav = reader.readLine().trim().split("[\\s]+")[0];
          data = Array.stringArray(targets.length + 1, ".");
          while (reader.ready()) {
            line = reader.readLine().trim().split("[\\s]+");
            if (i == 0) {
              index = ext.indexOfStr(line[0], targets);
              hash.put(line[0], index + "");
            } else {
              index = Integer.parseInt(hash.get(line[0]));
            }

            if (index != -1) {
              if (line.length < 3) {
                System.err.println("Error in '" + files[i] + "': " + Array.toStr(line));
              } else {
                data[index + 1] = line[1];
                data[0] = line[2];
              }
            }
            count++;
          }
          writer.println(trav + "\t" + files[i] + "\t" + Array.toStr(data));
          writer.flush();
          if (i == 0) {
            count = 0;
            for (int j = 1; j < data.length; j++) {
              if (!data[j].equals(".")) {
                count++;
              }
            }
            System.out.println(
                "Successfully identified " + count + " of " + (data.length - 1) + " transcripts");
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + files[i] + "\" not found in current directory");
          fnfe.printStackTrace();
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + files[i] + "\"");
          ioe.printStackTrace();
          System.exit(2);
        } catch (Exception e) {
          System.err.println("Error reading file \"" + files[i] + "\"");
          e.printStackTrace();
          System.exit(2);
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename) + "_transcripts.xln");
      e.printStackTrace();
    }
  }

  public static void regress(String filename) {
    PrintWriter writer;
    String[][] data;
    double[] indeps, deps;
    LeastSquares ls;
    String[] line;
    double[][] resids;

    data = HashVec.loadFileToStringMatrix(filename, false, null, false);
    indeps = Array.toDoubleArray(Array.subArray(Matrix.extractColumn(data, 2), 2));
    resids = new double[data[0].length - 3][];
    for (int i = 3; i < data[0].length; i++) {
      deps = Array.toDoubleArray(Array.subArray(Matrix.extractColumn(data, i), 2));
      ls = new LeastSquares(deps, Matrix.toMatrix(indeps));
      resids[i - 3] = ls.getResiduals();
      resids[i - 3] = Transformations.transform(resids[i - 3], Transformations.NORMALIZE);
    }

    try {
      writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + "_residuals.xln"));
      writer.println(Array.toStr(data[0]));
      writer.println(Array.toStr(data[1]));
      for (int i = 2; i < data.length; i++) {
        line = data[i][1].trim().split("_");
        writer.print(line[1] + "\t" + line[2] + "_" + line[3].charAt(0) + "\t" + data[i][2]);
        for (double[] resid : resids) {
          writer.print("\t" + resid[i - 2]);
        }
        writer.println();
      }

      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename, false) + "_residuals.xln");
      e.printStackTrace();
    }
  }
}
