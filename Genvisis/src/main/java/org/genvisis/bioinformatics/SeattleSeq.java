package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class SeattleSeq {
  public static final String[][] NEEDS =
      {{"Chr"}, {"MapInfo", "Position"}, {"MarkerName", "SNP"}, {"RefStrand"}};
  public static final String[][] RELEVANTS =
      {{"chromosome"}, {"position"}, {"sampleAlleles"}, {"accession"}, {"functionGVS"},
          {"aminoAcids"}, {"geneList"}, {"inDBSNPOrNot"}, {"rsID"}, {"microRNAs"}};
  // "# inDBSNPOrNot", "", "position", "referenceBase", "sampleGenotype", "sampleAlleles",
  // "allelesDBSNP", "accession", "functionGVS", "functionDBSNP", "rsID", "aminoAcids",
  // "proteinPosition", "cDNAPosition", "polyPhen", "granthamScore", "scorePhastCons",
  // "consScoreGERP", "chimpAllele", "CNV", "geneList", "AfricanHapMapFreq", "EuropeanHapMapFreq",
  // "AsianHapMapFreq", "hasGenotypes", "dbSNPValidation", "repeatMasker", "tandemRepeat",
  // "clinicalAssociation", "distanceToSplice", "microRNAs", "proteinSequence"
  public static final String[] ORDER = {"frameshift-near-splice", "frameshift",
      "stop-gained-near-splice", "stop-gained", "stop-lost-near-splice", "stop-lost", "nonsense",
      "missense-near-splice", "missense", "coding-near-splice", "coding", "splice-donor",
      "splice-5", "splice-acceptor", "splice-3", "synonymous-near-splice", "intron-near-splice",
      "coding-unknown-near-splice", "synonymous", "coding-synonymous", "coding-unknown",
      "coding-notMod3", "non-coding-exon-near-splice", "non-coding-exon", "5-prime-UTR", "utr-5",
      "3-prime-UTR", "utr-3", "intron", "upstream-gene", "near-gene-5", "downstream-gene",
      "near-gene-3", "intergenic"};
  // public static final String[] BAD = {"missense", "stop-gained", "stop-lost",
  // "missense-near-splice", "splice-donor", "splice-acceptor"};
  public static final String[] BAD =
      {"frameshift-near-splice", "frameshift", "stop-gained-near-splice", "stop-gained",
          "stop-lost-near-splice", "stop-lost", "nonsense", "missense-near-splice", "missense",
          "coding", "splice-donor", "splice-5", "splice-acceptor", "splice-3"};
  public static final String[] NEUTRAL = {"intron-near-splice", "5-prime-UTR", "downstream-gene",
      "upstream-gene", "synonymous", "coding-synonymous", "intergenic", "non-coding-exon",
      "3-prime-UTR", "intron", "coding-notMod3"};

  public static char[] determineAlleles(String str) {
    String[] line;

    if (str.startsWith("[")) {
      str = str.substring(1);
    }

    if (str.startsWith("]")) {
      str = str.substring(0, str.length() - 1);
    }

    line = str.split("/");

    return new char[] {line[0].charAt(0), line[1].charAt(0)};
  }

  public static Hashtable<String, String[]> loadAllAnnotationInDir(String directory, Logger log) {
    BufferedReader reader;
    String[] files, line;
    Hashtable<String, String[]> finalHash;
    String markerName;
    String[] trav;
    String temp;
    Vector<String[]> v;
    int[] indices;
    String prev;
    boolean done;
    int type, worstType, worst;
    int linesSkipped;
    boolean problem;
    long time;

    problem = false;
    finalHash = new Hashtable<String, String[]>();

    if (directory == null) {
      log.reportError("The SeattleSeq annotation directory was null; returning an empty hashtable");
    } else if (!Files.exists(directory) || !Files.isDirectory(directory)) {
      log.reportError("Error - SeattleSeq annotation directory directory not found: " + directory);
      log.reportError("        returning an empty hashtable");
    } else {
      files = Files.list(directory, "SeattleSeqAnnotation", ".txt.gz", false, false);
      log.report(ext.getTime() + "\tFound " + files.length
          + " file(s) with a .SeattleSeq extension to include");
      for (int f = 0; f < files.length; f++) {
        time = new Date().getTime();
        log.report(
            "Processing SeattleSeq file " + (f + 1) + " of " + files.length + "\t" + files[f],
            false, true);
        try {
          reader = Files.getAppropriateReader(directory + files[f]);

          temp = reader.readLine().trim();
          temp = temp.substring(1).trim();
          line = temp.split("[\\s]+");
          indices = ext.indexFactors(RELEVANTS, line, false, true, true, log, true);

          prev = "";

          done = false;
          trav = null;
          linesSkipped = 0;
          v = new Vector<String[]>();
          while (!done) {
            if (reader.ready()) {
              line = reader.readLine().trim().split("[\\s]+");
              if (line[0].startsWith("#") || line.length < Array.max(indices)) {
                linesSkipped++;
                continue;
              }
            } else {
              done = true;
            }

            try {
              markerName = "chr" + line[1] + ":" + line[2] + "_" + line[3] + "_" + line[4];
              if (line[4].startsWith("I") || line[4].startsWith("D")) {
                markerName =
                    "chr" + line[1] + ":" + line[2] + "_" + line[3].charAt(0) + "_" + line[5];
              }
              if (done || (v.size() > 0 && !markerName.equals(prev))) {
                worst = -1;
                worstType = ORDER.length;
                type = -1;
                for (int i = 0; i < v.size(); i++) {
                  trav = v.elementAt(i);
                  type = ext.indexOfStr(trav[indices[4]], ORDER);
                  if (type == -1) {
                    log.reportError("unknown type: " + trav[indices[4]]);
                    log.reportError(Array.toStr(trav));
                    problem = true;
                  }
                  if (type < worstType) {
                    worstType = type;
                    worst = i;
                  }
                }
                trav = v.elementAt(worst);
                if (trav[indices[8]].equals("0")) {
                  trav[indices[8]] = ".";
                } else {
                  trav[indices[8]] = "rs" + trav[indices[8]];
                }

                trav = Array.subArray(trav, indices); // trim to the relevant columns
                trav = Array.addStrToArray(ext.indexOfStr(trav[4], BAD) >= 0 ? "1" : "0", trav, 0); // determine
                                                                                                    // if
                                                                                                    // mutation
                                                                                                    // is
                                                                                                    // missense/splice
                                                                                                    // or
                                                                                                    // worse
                finalHash.put(prev, trav);
                v = new Vector<String[]>();
              }
              v.add(line);
              prev = markerName;
            } catch (Exception e) {
              log.reportError("Error reading line: " + Array.toStr(line));
              log.reportException(e);
              problem = true;
              System.exit(1);
            }

          }

          log.report(
              " ...skipped " + linesSkipped + " line(s); finished in " + ext.getTimeElapsed(time));
          reader.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println(
              "Error: file \"" + directory + files[f] + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + directory + files[f] + "\"");
          System.exit(2);
        }
      }
    }

    if (problem) {
      return null;
    } else {
      return finalHash;
    }
  }

  public static Hashtable<String, String[]> loadAllBadAnnotationInDir(String directory,
      Logger log) {
    BufferedReader reader;
    String[] files, line;
    Hashtable<String, String[]> hash;
    String markerName;
    String function;

    hash = new Hashtable<String, String[]>();

    if (directory == null) {
      log.reportError("The SeattleSeq annotation directory was null; returning an empty hashtable");
    } else if (!Files.exists(directory) || !Files.isDirectory(directory)) {
      log.reportError("Error - SeattleSeq annotation directory directory not found: " + directory);
      log.reportError("        returning an empty hashtable");
    } else {
      files = Files.list(directory, "SeattleSeqAnnotation", ".txt.gz", false, false);
      log.report("Found " + files.length + " file(s) with a .SeattleSeq extension to include");
      for (String file : files) {
        try {
          reader = Files.getAppropriateReader(directory + file);
          while (reader.ready()) {
            line = reader.readLine().trim().split("\t", -1);
            if (line.length > 1) {
              markerName = "chr" + line[1] + ":" + line[2] + "_" + line[3] + "_" + line[4];
              if (!hash.containsKey(markerName) || hash.get(markerName) == null) {
                if (ext.indexOfStr(line[8], BAD) >= 0) {
                  function = line[8];
                  if (!line[11].equals("none")) {
                    function += " " + ext.replaceAllWith(line[11], ",",
                        line[12].substring(0, line[12].indexOf("/")));
                  }
                  function += "\t" + line[0];
                  hash.put(markerName, new String[] {function});
                } else {
                  hash.put(markerName, new String[0]);
                }
              }
            }

            // TODO
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          System.err
              .println("Error: file \"" + directory + file + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + directory + file + "\"");
          System.exit(2);
        }
      }
    }

    return hash;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "SeattleSeq.dat";
    String prep = null;
    String prepFreq = null; // "D:/home/npankrat/NCBI/ESP_exome_chip/EVS/freqInfo.out";
    String freqFilename = "D:/home/npankrat/NCBI/ESP_exome_chip/EVS/freqInfo_proc.dat";


    // Hashtable<String, String[]> hash = loadAllAnnotationInDir("D:/Logan/SuperNovo/SeattleSeq/",
    // new Logger());
    // String[] keys = HashVec.getKeys(hash);
    // for (int i = 0; i < keys.length; i++) {
    // System.out.println(keys[i]+"\t"+Array.toStr(hash.get(keys[i])));
    // }
    // System.exit(1);

    String usage = "\n" + "bioinformatics.SeattleSeq requires 0-1 arguments\n"
        + "   (1) filename of SeattleSeq annotation (i.e. file=" + filename + " (default))\n"
        + "   (2) filename of allele frequency data (i.e. freqFile=" + freqFilename
        + " (default))\n" + " OR:\n"
        + "   (1) name of file to prep for input to SeattleSeq (i.e. prep=filename.txt (not the default))\n"
        + " OR:\n"
        + "   (1) name of allele freq file to process (i.e. prepFreq=freqInfo.out (not the default))\n"
        + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("prep=")) {
        prep = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("prepFreq=")) {
        prepFreq = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("freqFile=")) {
        freqFilename = arg.split("=")[1];
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
      if (prep != null) {
        proc(filename);
      } else if (prepFreq != null) {
        parseFreq(prepFreq);
      } else {
        summarize(filename, freqFilename);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parseFreq(String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    double a1, a2, temp;

    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + "_proc.dat"));
      writer.println(reader.readLine());
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        for (int i = 1; i < line.length; i++) {
          a1 = ext.parseDoubleArg(line[i].split("/")[0]);
          a2 = ext.parseDoubleArg(line[i].split("/")[1]);
          if (a2 < a1) {
            temp = a1;
            a1 = a2;
            a2 = temp;
          }
          line[i] = a1 / (a1 + a2) + "";
        }
        writer.println(Array.toStr(line));
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

  public static void proc(String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    char temp;
    int[] indices;
    char[] alleles;

    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + ".input"));
      line = reader.readLine().trim().split("[\\s]+");
      indices = ext.indexFactors(NEEDS, line, false, true, true, true);
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (!line[indices[0]].startsWith("chr")) {
          line[indices[0]] = "chr" + line[indices[0]];
        }
        try {
          Integer.parseInt(line[indices[1]]);
        } catch (Exception e) {
          System.err.println("Error - invalid position ('" + line[indices[1]] + "')");
        }
        alleles = determineAlleles(line[indices[2]]);
        if (line[indices[3]].equals("-")) {
          temp = alleles[0];
          alleles[0] = alleles[1];
          alleles[1] = temp;
        } else if (!line[indices[3]].equals("+")) {
          System.err.println("Error - invalid strand ('" + line[indices[3]] + "')");
        }
        for (int i = 0; i < alleles.length; i++) {
          if (alleles[i] == 'I' || alleles[i] == 'D') {
            alleles[i] = 'A';
          }
        }
        writer.println(
            line[indices[0]] + "\t" + line[indices[1]] + "\t0\t" + alleles[0] + "\t" + alleles[1]);
      }
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  // the indels are not captured by SeattleSeq
  // there are triallelic markers on the exome chip, hence the Hashtable, but that doesn't
  // necessarily capture all, still missing a little over a hundred
  public static void summarize(String filename, String freqFilename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, trav;
    String temp;
    Hashtable<String, Vector<String[]>> hash = new Hashtable<String, Vector<String[]>>();
    Vector<String[]> v;
    int[] indices;
    Logger log;
    String prev;
    boolean done;
    int type, worstType, worst;
    String[] keys;
    Hashtable<String, String> hashFreq;
    double freq;
    int linesSkipped;

    hashFreq = HashVec.loadFileToHashString(freqFilename, new int[] {0}, new int[] {3}, false, "",
        true, false, false);

    try {
      reader = Files.getAppropriateReader(filename);
      writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + "_summary.out"));
      writer.println(
          Array.toStr(Matrix.extractColumn(RELEVANTS, 0)) + "\tAlleleFrequency\tMAF<1%\tMAF<5%");
      log = new Logger(ext.rootOf(filename, false) + ".log");
      temp = reader.readLine().trim();
      // if (temp.startsWith("#"));
      temp = temp.substring(1).trim();
      line = temp.split("[\\s]+");
      indices = ext.indexFactors(RELEVANTS, line, false, true, true, log, true);

      prev = "";

      done = false;
      linesSkipped = 0;
      while (!done) {
        if (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");
          if (line.length < Array.max(indices)) {
            linesSkipped++;
            continue;
          }
        } else {
          done = true;
        }
        if (hash.size() > 0 && !(line[indices[0]] + "\t" + line[indices[1]]).equals(prev) || done) {
          keys = HashVec.getKeys(hash);
          for (String key2 : keys) {
            v = hash.get(key2);
            worst = -1;
            worstType = ORDER.length;
            type = -1;
            for (int i = 0; i < v.size(); i++) {
              trav = v.elementAt(i);
              type = ext.indexOfStr(trav[indices[4]], ORDER);
              if (type == -1) {
                System.out.println("unknown type: " + trav[indices[4]]);
                System.out.println(Array.toStr(trav));
                System.exit(1);
              }
              if (type < worstType) {
                worstType = type;
                worst = i;
              }
            }
            trav = v.elementAt(worst);
            if (trav[indices[8]].equals("0")) {
              trav[indices[8]] = ".";
            } else {
              trav[indices[8]] = "rs" + trav[indices[8]];
            }
            for (int i = 0; i < indices.length; i++) {
              writer.print((i == 0 ? "" : "\t") + trav[indices[i]]);
            }
            if (hashFreq.containsKey(trav[indices[0]].toLowerCase() + "^" + trav[indices[1]])) {
              freq = Double.parseDouble(
                  hashFreq.get(trav[indices[0]].toLowerCase() + "^" + trav[indices[1]]));
              writer
                  .print("\t" + freq + "\t" + (freq < 0.01 ? 1 : 0) + "\t" + (freq < 0.05 ? 1 : 0));
            } else {
              writer.print("\t.\t.\t.");
            }
            writer.println();
            writer.flush();
          }
          hash = new Hashtable<String, Vector<String[]>>();
        }
        HashVec.addToHashArrayVec(hash, line[indices[2]], line);
        prev = line[indices[0]] + "\t" + line[indices[1]];
      }
      writer.close();
      reader.close();
      System.out.println("Skipped " + linesSkipped + " line(s)");
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }
}
