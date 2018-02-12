// creates a table with 95% confidence intervals in the format requested by Lindsay Farrer
package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.stats.ProbDist;

public class FarrerParser {

  public static final String[] PLINK_FRQ_HEADER = {"CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"};
  public static final String[] LOGISTIC_HEADER_WITH_CI = {"CHR", "SNP", "BP", "A1", "TEST", "NMISS",
                                                          "OR", "SE", "L95", "U95", "STAT", "P"};
  public static final String[] PROBABEL_HEADER = {"name", "A1", "A2", "Freq1", "MAF", "Quality",
                                                  "Rsq", "n", "Mean_predictor_allele", "chrom"};

  public static void parse(String dir, String filename, String modelsFile, String probabelPrefix) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String model, trav;
    Hashtable<String, Hashtable<String, String>> hashes = new Hashtable<String, Hashtable<String, String>>();
    Hashtable<String, String> hash;
    Vector<String> v;
    String[] markers, models, files, finalSummary;
    String summary;
    double freq, beta, stderr;
    int numTerms, count;
    String[] terms;

    markers = HashVec.loadFileToStringArray(dir + filename, false, new int[] {0}, false);
    models = HashVec.loadFileToStringArray(dir + modelsFile, false, new int[] {0}, false);

    try {
      reader = new BufferedReader(new FileReader(dir + "plink.frq"));
      ext.checkHeader(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE), PLINK_FRQ_HEADER,
                      true);
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        if (ext.indexOfStr(line[1], markers) >= 0) {
          freq = Double.parseDouble(line[4]);
          summary = line[2] + "(" + ext.formDeci(freq, 4, true) + ")\t" + line[3] + "("
                    + ext.formDeci(1 - freq, 4, true) + ")";
          HashVec.addToHashHash(hashes, "freq", line[1], summary);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + "plink.frq"
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + "plink.frq" + "\"");
      System.exit(2);
    }

    files = Files.list(dir, "_add.out.txt");
    for (String file : files) {
      try {
        reader = new BufferedReader(new FileReader(dir + file));
        model = file.substring(0, file.indexOf("_chr"));
        if (model.startsWith(probabelPrefix)) {
          model = model.substring(probabelPrefix.length());
        }
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        for (int j = 0; j < PROBABEL_HEADER.length; j++) {
          if (!line[j].equalsIgnoreCase(PROBABEL_HEADER[j])) {
            System.err.println("Error - Expecting " + PROBABEL_HEADER[j] + " in column " + (j + 1)
                               + "; got " + line[j]);
          }
        }
        numTerms = (line.length - PROBABEL_HEADER.length - 1) / 2;
        terms = new String[numTerms];
        for (int j = 0; j < numTerms; j++) {
          if (line[PROBABEL_HEADER.length + j * 2 + 0].startsWith("beta_")
              && line[PROBABEL_HEADER.length + j * 2 + 1].startsWith("sebeta_")
              && line[PROBABEL_HEADER.length + j * 2
                      + 0].substring(("beta_").length())
                          .equals(line[PROBABEL_HEADER.length + j * 2
                                       + 1].substring(("sebeta_").length()))) {
            terms[j] = line[PROBABEL_HEADER.length + j * 2 + 0].substring(("beta_").length());
          } else {
            System.err.println("Error - expecting beta_ and sebeta_ prefixes on each of the "
                               + terms.length + " term(s), but found: ");
            System.err.println("        " + line[PROBABEL_HEADER.length + j * 2 + 0] + " and "
                               + line[PROBABEL_HEADER.length + j * 2 + 1]);
          }
        }
        HashVec.addToHashHash(hashes, model, "terms", ArrayUtils.toStr(terms));

        while (reader.ready()) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          if (ext.indexOfStr(line[0], markers) >= 0) {
            freq = Double.parseDouble(line[3]);
            for (int j = 0; j < numTerms; j++) {
              beta = Double.parseDouble(line[PROBABEL_HEADER.length + j * 2 + 0]);
              stderr = Double.parseDouble(line[PROBABEL_HEADER.length + j * 2 + 1]);
              if (freq > 0.50) {
                summary = line[2] + "(" + ext.formDeci(1 - freq, 4, true) + ")\t" + line[1] + "("
                          + ext.formDeci(freq, 4, true) + ")\t";
                if (terms[j].startsWith("SNP_")) {
                  beta *= -1;
                }
              } else {
                summary = line[1] + "(" + ext.formDeci(freq, 4, true) + ")\t" + line[2] + "("
                          + ext.formDeci(1 - freq, 4, true) + ")\t";
              }
              // summary +=
              // ext.formDeci(Math.exp(beta),3,true)+"\t("+ext.formDeci(Math.exp(beta-1.96*stderr),3,true)+",
              // "+ext.formDeci(Math.exp(beta+1.96*stderr),3,true)+")\t"+ProbDist.ChiDist(Double.parseDouble(line[interaction?14:12]),
              // 1);
              summary += ext.formDeci(Math.exp(beta), 3, true) + "\t("
                         + ext.formDeci(Math.exp(beta - 1.96 * stderr), 3, true) + ", "
                         + ext.formDeci(Math.exp(beta + 1.96 * stderr), 3, true) + ")\t"
                         + ProbDist.NormDist(beta / stderr);
              // summary +=
              // ext.formDeci(Math.exp(beta),3,true)+"\t("+ext.formDeci(Math.exp(beta-1.96*stderr),3,true)+",
              // "+ext.formDeci(Math.exp(beta+1.96*stderr),3,true)+")\t=NORMSDIST("+(-1*Math.abs(beta/stderr))+")*2";
              HashVec.addToHashHash(hashes, model, line[0] + "_" + terms[j], summary);
            }
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + dir + file + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + file + "\"");
        System.exit(2);
      }
    }

    files = Files.list(dir, ".assoc.logistic");
    for (String file : files) {
      try {
        reader = new BufferedReader(new FileReader(dir + file));
        model = file.substring(0, file.indexOf(".assoc.logistic"));
        ext.checkHeader(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE),
                        LOGISTIC_HEADER_WITH_CI, true);
        v = new Vector<String>();
        while (reader.ready()) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          if (ext.indexOfStr(line[1], markers) >= 0) {
            summary = ".\t.\t"
                      + (line[6].equals("NA") ? "NA"
                                              : ext.formDeci(Double.parseDouble(line[6]), 3, true))
                      + "\t("
                      + (line[6].equals("NA") ? "NA"
                                              : ext.formDeci(Double.parseDouble(line[8]), 3, true))
                      + ", "
                      + (line[9].equals("inf") ? "inf"
                                               : (line[6].equals("NA") ? "NA"
                                                                       : ext.formDeci(Double.parseDouble(line[9]),
                                                                                      3, true)))
                      + ")\t" + line[11];
            if (line[4].equals("ADD")) {
              trav = "SNP_add";
            } else if (line[4].startsWith("ADDx")) {
              trav = "SNP_" + line[4].substring(("ADDx").length());
            } else {
              trav = line[4];
            }
            HashVec.addToHashHash(hashes, model, line[1] + "_" + trav, summary);
            HashVec.addIfAbsent(trav, v);
          }
        }
        if (!hashes.get(model).containsKey("terms")) {
          HashVec.addToHashHash(hashes, model, "terms",
                                ArrayUtils.toStr(ArrayUtils.toStringArray(v)));
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + dir + file + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + file + "\"");
        System.exit(2);
      }
    }

    try {
      writer = Files.openAppropriateWriter(dir + ext.rootOf(filename) + "_parsed.xln");
      writer.print("Marker\tMinor allele\tMajor allele\tImputed minor allele\tImputed major allele");
      count = 0;
      for (String model2 : models) {
        terms = hashes.get(model2).get("terms").trim().split(PSF.Regex.GREEDY_WHITESPACE);
        count += terms.length;
        if (terms.length == 1) {
          writer.print("\t" + model2 + " OR\t95% CI\tp-value");
        } else {
          for (String term : terms) {
            writer.print("\t" + model2 + "_" + term + " OR\t95% CI\tp-value");
          }
        }
      }
      numTerms = count;
      writer.println();
      for (String marker : markers) {
        if (marker.equals("")) {
          writer.println();
        } else {
          writer.print(marker);
          hash = hashes.get("freq");
          if (hash != null && hash.containsKey(marker)) {
            writer.print("\t" + hash.get(marker));
          } else {
            writer.print("\t.\t.");
          }
          finalSummary = ArrayUtils.stringArray(2 + 3 * numTerms, "");
          count = 0;
          for (String model2 : models) {
            hash = hashes.get(model2);
            terms = hash.get("terms").trim().split(PSF.Regex.GREEDY_WHITESPACE);
            for (String term : terms) {
              if (hash != null && hash.containsKey(marker + "_" + term)) {
                line = hash.get(marker + "_" + term).split("\\t");
                if (!line[0].equals(".")) {
                  finalSummary[0] = line[0];
                  finalSummary[1] = line[1];
                }
                for (int m = 0; m < 3; m++) {
                  finalSummary[2 + count * 3 + m] = line[2 + m];
                }
              } else {
                for (int m = 0; m < 3; m++) {
                  finalSummary[2 + count * 3 + m] = ".";
                }
              }
              count++;
            }
          }
          writer.println("\t" + ArrayUtils.toStr(finalSummary));
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + ext.rootOf(filename) + "_parsed.xln");
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\candis\\Confidence
    // Intervals\\CI95_every_which_way\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\UMN\\Myron\\ExcisionPathway\\RareVariants\\summary\\";
    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\candis\\Confidence Intervals\\E4_interaction\\";
    String filename = "order.txt";
    String models = "models.txt";
    String prefix = "LOAD_Aff_";

    String usage = "\n" + "dead.FarrerParser requires 0-1 arguments\n"
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
    try {
      parse(dir, filename, models, prefix);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
