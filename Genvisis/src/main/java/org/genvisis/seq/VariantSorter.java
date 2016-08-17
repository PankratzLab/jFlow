// -Xms1024M -Xmx1024M
package org.genvisis.seq;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Collapsed;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;

public class VariantSorter {
  public static final double MAF_LOWER_BOUND = 0.01;
  public static final String[] SEATTLE_SEQ_REQS =
      {"geneList", "functionGVS", "# inDBSNPOrNot", "chromosome", "position", "referenceBase",
       "sampleGenotype", "AfricanHapMapFreq", "EuropeanHapMapFreq", "AsianHapMapFreq"};
  public static final String[] SIFT_REQS = {"Coordinates", "Prediction"};
  public static final String[] CATS =
      {"Total", "dbSNP>" + (MAF_LOWER_BOUND * 100) + "%", "dbSNP<" + (MAF_LOWER_BOUND * 100) + "%",
       "dbSNP_noFreq", "not_in_dbSNP", "deleterious", "deleteriousRare", "deleteriousNotInDBSNP"};
  public static final String[] FUNCS =
      {"missense", "nonsense", "coding-synonymous", "coding-notMod3", "splice-5", "splice-3",
       "utr-5", "utr-3", "near-gene-5", "near-gene-3", "intron", "intergenic"};
  public static final boolean[] FUNC_DISPS =
      {true, true, true, true, true, true, true, true, false, false, false, false};
  public static final String[] SENSE = {"A", "C", "G", "T", "I", "D"};
  public static final String[] ANTISENSE = {"T", "G", "C", "A", "I", "D"};
  public static final String[][] IUB_LOOKUPS =
      {{"A", "A"}, {"C", "C"}, {"G", "G"}, {"T", "T"}, {"R", "AG"}, {"Y", "CT"}, {"K", "GT"},
       {"M", "AC"}, {"S", "GC"}, {"W", "AT"}, {"B", "CGT"}, {"D", "AGT"}, {"H", "ACT"},
       {"V", "ACG"}, {"N", "AGCT"}};

  public static String convertIUBlookupCodes(String code, String ref) {
    char refAllele;
    String str, remaining;

    if (ref.length() != 1) {
      System.err.println("Error - reference allele must be a single nucleotide (not '" + ref
                         + "')");
      return null;
    }

    refAllele = ref.toUpperCase().charAt(0);
    if (refAllele != 'A' && refAllele != 'C' && refAllele != 'G' && refAllele != 'T') {
      System.err.println("Error - reference allele must be a single nucleotide (not '" + ref
                         + "')");
      return null;
    }

    str = "";
    for (String[] element : IUB_LOOKUPS) {
      if (code.equals(element[0])) {
        str = element[1];
      }
    }
    if (str.equals("")) {
      System.err.println("Error - '" + code + "' is not a recognized IUB code");
      return null;
    }

    remaining = "";
    for (int i = 0; i < str.length(); i++) {
      if (str.charAt(i) != refAllele) {
        remaining += str.charAt(i);
      }
    }

    if (remaining.length() > 1) {
      System.err.println("Error - '" + code
                         + "' is an ambiguous IUB code; codes for more than 2 alleles (" + str
                         + ") ");
      return null;
    }

    return remaining;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\SequencingProjectWithCIDR\\SeattleSeq_files\\OnTarget\\";
    // String suffix = "_snplist_filterbyBED_seattle_seq_output.txt";
    // String suffix = "_snplist_seattle_seq_output.txt";
    // String suffix = ".txt";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\SequencingProjectWithCIDR\\SNPlists\\older_variant_files\\by_bed\\filtered\\";
    // String suffix = ".snplist";

    String dir =
        "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\SequencingProjectWithCIDR\\SNPlists\\bfast_indels\\target\\";
    String suffix = ".snplist";

    String genes =
        "ESPNL,ARMC9,DNER,FAM132B,ALPPL2,DGKD,DIS3L2,EFHD1,MMRN1,MYEOV2,COL4A3,HRB,LRRK2,SH3BP4";
    boolean indels = true;

    String usage = "\n" + "seq.VariantSorter requires 0-1 arguments\n"
                   + "  ## run once to get SIFT input files; concatenate results into siftOutput.dat and re-run this program to incorporate data ##\n"
                   + "   (1) directory with SeattleSeq files (i.e. dir=" + dir + " (default))\n"
                   + "   (2) suffix with which to filter files (i.e. suffix=" + suffix
                   + " (default))\n" + "   (3) comma delimited list of favorite genes (i.e. genes="
                   + genes + " (default))\n"
                   + "   (4) files are snplists with indels (i.e. -indels (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("suffix=")) {
        suffix = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("genes=")) {
        suffix = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-indels")) {
        indels = ext.parseBooleanArg(arg);
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (indels) {
        parseIndels(dir, suffix, genes);
      } else {
        parse(dir, suffix, genes);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parse(String dir, String suffix, String favoriteGenes) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, genes, favGenes;
    String trav, gene, func, name, strand;
    Hashtable<String, Vector<String>> hash, favHits;
    Hashtable<String, String> siftInfo, geneNames, geneInfo, variants;
    Vector<String> v;
    String[] files;
    int[] indices;
    int[][] counts;
    int funcIndex, catIndex;
    boolean in_dbSNP;
    double freq;
    int orientation;

    favGenes = favoriteGenes.split(",");
    files = Collapsed.list(dir, suffix, false);
    System.out.println("Found " + files.length + " files to parse with the suffix: " + suffix);

    System.out.println("Loading aliases...");
    geneNames = Collapsed.loadFileToHashString(dir + "kgAlias.xln", 1, new int[] {0}, "\t", false);
    System.out.println("Loading strand information...");
    geneInfo = Collapsed.loadFileToHashString(dir + "knownGene.xln", 0, new int[] {2}, "\t", false);

    System.out.println("Parsing files...");
    hash = new Hashtable<String, Vector<String>>();
    variants = new Hashtable<String, String>();
    for (int i = 0; i < files.length; i++) {
      if (!files[i].startsWith("NA")) {
        System.out.print(".");
        try {
          reader = new BufferedReader(new FileReader(dir + files[i]));
          indices = Collapsed.indexFactors(SEATTLE_SEQ_REQS, reader.readLine().trim().split("\\t"),
                                           true, true);
          while (reader.ready()) {
            line = reader.readLine().trim().split("\\t");
            if (!line[0].startsWith("#")) {
              gene = line[indices[0]];
              func = line[indices[1]];
              if (Collapsed.indexOfStr(func, FUNCS) == -1) {
                System.err.println("Error - new func type: " + func);
                System.exit(1);
              }
              if (!gene.equals("none")) {
                trav = files[i];
                for (int j = 1; j < indices.length; j++) {
                  trav += "\t" + line[indices[j]];
                }
                Collapsed.addToHashVec(hash, gene, trav, true);
                if (func.equals("missense")) {
                  if (gene.indexOf(",") > 0 || !geneNames.containsKey(gene)) {
                    variants.put(line[indices[3]] + "," + line[indices[4]] + "," + "1" + ","
                                 + line[indices[5]] + "/"
                                 + convertIUBlookupCodes(line[indices[6]], line[indices[5]]), "");
                  } else {
                    name = geneNames.get(gene);
                    if (geneInfo.containsKey(name)) {
                      strand = geneInfo.get(name);
                      if (strand.equals("+")) {
                        orientation = 1;
                      } else if (strand.equals("-")) {
                        orientation = -1;
                      } else {
                        System.err.println("Error - '" + strand
                                           + "' is an invalid strand designation in knownGene");
                        orientation = -999;
                      }
                      if (orientation == 1) {
                        variants.put(line[indices[3]] + "," + line[indices[4]] + ",1,"
                                     + line[indices[5]] + "/"
                                     + convertIUBlookupCodes(line[indices[6]], line[indices[5]]),
                                     "");
                      } else if (orientation == -1) {
                        variants.put(line[indices[3]] + "," + line[indices[4]] + ",-1,"
                                     + ANTISENSE[Collapsed.indexOfStr(line[indices[5]], SENSE,
                                                                      false, true)]
                                     + "/"
                                     + ANTISENSE[Collapsed.indexOfStr(convertIUBlookupCodes(line[indices[6]],
                                                                                            line[indices[5]]),
                                                                      SENSE, false, true)],
                                     "");
                      }
                    } else {
                      System.err.println("Error - no gene name lookup for '" + gene
                                         + "' in kgAlias");
                      name = null;
                    }
                  }
                }
              }
            }
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + dir + files[i]
                             + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + dir + files[i] + "\"");
          System.exit(2);
        }
      }
    }

    Collapsed.writeList(Collapsed.getKeys(variants), dir + "siftInput.dat");

    siftInfo = new Hashtable<String, String>();
    if (Collapsed.exists(dir + "siftOutput.dat", false)) {
      try {
        reader = new BufferedReader(new FileReader(dir + "siftOutput.dat"));
        indices =
            Collapsed.indexFactors(SIFT_REQS, reader.readLine().trim().split("\\t"), true, true);
        while (reader.ready()) {
          line = reader.readLine().trim().split("\\t");
          if (line[indices[1]].startsWith("DAMAGING")) {
            siftInfo.put(line[indices[0]].split(",")[0] + "\t" + line[indices[0]].split(",")[1],
                         "");
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + dir + "siftOutput.dat"
                           + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + "siftOutput.dat" + "\"");
        System.exit(2);
      }
    }

    genes = Collapsed.getKeys(hash);
    favHits = new Hashtable<String, Vector<String>>();
    try {
      writer = new PrintWriter(new FileWriter(dir + "bins.xln"));
      for (String element : CATS) {
        for (int j = 0; j < FUNCS.length; j++) {
          if (FUNC_DISPS[j]) {
            writer.print("\t" + element);
          }
        }
      }
      writer.println();
      writer.print("Gene");
      for (String element : CATS) {
        for (int j = 0; j < FUNCS.length; j++) {
          if (FUNC_DISPS[j]) {
            writer.print("\t" + FUNCS[j]);
          }
        }
      }
      writer.println();
      for (String gene2 : genes) {
        v = hash.get(gene2);
        counts = new int[CATS.length][FUNCS.length];
        for (int j = 0; j < v.size(); j++) {
          line = v.elementAt(j).split("[\\s]+");
          funcIndex = Collapsed.indexOfStr(line[1], FUNCS);
          if (funcIndex == 0 && Collapsed.indexOfStr(gene2, favGenes) >= 0) {
            Collapsed.addToHashVec(favHits, gene2, Collapsed.toStr(line), false);
          }

          in_dbSNP = !line[2].equals("none");
          if (in_dbSNP) {
            freq = -1;
            for (int k = 0; k < 3; k++) {
              if (!line[7 + k].equals("NA")) {
                freq = Math.max(freq, Double.parseDouble(line[7 + k]));
              }
            }
            if (freq == -1) {
              catIndex = 3;
            } else if (freq < MAF_LOWER_BOUND) {
              catIndex = 2;
            } else {
              catIndex = 1;
            }
          } else {
            catIndex = 4;
          }

          // if (funcIndex == 0 && Collapsed.indexOfStr(genes[i], favGenes) >= 0) {
          if (funcIndex == 0 && catIndex == 4 && Collapsed.indexOfStr(gene2, favGenes) >= 0) {
            Collapsed.addToHashVec(favHits, gene2, Collapsed.toStr(line), false);
          }
          if (funcIndex <= 1 && catIndex == 4) {
            Collapsed.addToHashVec(favHits, "novelMissense", gene2 + "\t" + Collapsed.toStr(line),
                                   false);
          }

          counts[0][funcIndex]++;
          counts[catIndex][funcIndex]++;

          if (siftInfo.containsKey(line[3] + "\t" + line[4])) {
            counts[5][funcIndex]++;
            if (catIndex > 1) {
              counts[6][funcIndex]++;
            }
            if (catIndex == 4) {
              counts[7][funcIndex]++;
              // System.out.println(genes[i]+"\t"+Collapsed.toStr(line));
            }
          }

        }
        writer.print(gene2);
        for (int j = 0; j < CATS.length; j++) {
          for (int k = 0; k < FUNCS.length; k++) {
            if (FUNC_DISPS[k]) {
              writer.print("\t" + counts[j][k]);
            }
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + "bins.xln");
      e.printStackTrace();
    }

    genes = Collapsed.getKeys(favHits);
    for (String gene2 : genes) {
      Collapsed.writeList(Collapsed.toStringArray(favHits.get(gene2)), dir + gene2 + ".out");
    }
    Collapsed.writeList(Collapsed.toStringArray(favHits.get("novelMissense")),
                        dir + "novelMissense.out");

  }

  public static void parseIndels(String dir, String suffix, String favoriteGenes) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, genes, favGenes;
    Hashtable<String, Vector<String>> hash, favHits;
    Vector<String> v;
    String[] files, starts, stops;
    int[] loc;
    int[][][] geneLocs;
    String[] geneNames;
    int count;
    boolean exonic;
    Vector<String> plusHits, minusHits;

    favGenes = favoriteGenes.split(",");
    files = Collapsed.list(dir, suffix, false);
    System.out.println("Found " + files.length + " files to parse with the suffix: " + suffix);

    // System.out.println("Loading aliases...");
    // geneAliases = Collapsed.loadFileToHashString(dir+"kgAlias.xln", 0, new int[] {1}, "\t",
    // false);
    System.out.println("Loading gene/exon boundary information...");
    // geneInfo = Collapsed.loadFileToHashString(dir+"knownGene.xln", 0, new int[] {1,3,4,7,8,9,12},
    // "\t", false);

    count = Files.countLines(dir + "known_2q.xln", 1);
    geneLocs = new int[count][][];
    geneNames = new String[count];
    count = 0;
    try {
      reader = new BufferedReader(new FileReader(dir + "known_2q.xln"));
      reader.readLine(); // need to check header
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t", -1);
        geneNames[count] = line[12];
        geneLocs[count] = new int[Integer.parseInt(line[7]) + 2][];
        geneLocs[count][0] =
            new int[] {Positions.chromosomeNumber(line[1]), line[2].equals("+") ? 1 : 0};
        geneLocs[count][1] = new int[] {Integer.parseInt(line[3]), Integer.parseInt(line[4])};
        starts = line[8].split(",");
        stops = line[9].split(",");
        for (int i = 0; i < Integer.parseInt(line[7]); i++) {
          geneLocs[count][i + 2] =
              new int[] {Integer.parseInt(starts[i]), Integer.parseInt(stops[i])};
        }
        count++;
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + "known_2q.xln"
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + "known_2q.xln" + "\"");
      System.exit(2);
    }

    System.out.println("Parsing files...");
    hash = new Hashtable<String, Vector<String>>();
    for (int i = 0; i < files.length; i++) {
      if (!files[i].startsWith("NA")) {
        System.out.print(".");
        try {
          reader = new BufferedReader(new FileReader(dir + files[i]));
          while (reader.ready()) {
            line = reader.readLine().trim().split("\\t");
            if (!line[0].startsWith("#")) {
              loc = new int[] {Positions.chromosomeNumber(line[0]), Integer.parseInt(line[1])};
              plusHits = new Vector<String>();
              minusHits = new Vector<String>();
              for (int j = 0; j < geneLocs.length; j++) {
                if (loc[0] == geneLocs[j][0][0] && loc[1] >= geneLocs[j][1][0] - 5
                    && loc[1] <= geneLocs[j][1][1] + 5) {
                  exonic = false;
                  for (int k = 2; k < geneLocs[j].length; k++) {
                    if (loc[1] >= geneLocs[j][k][0] - 5 && loc[1] <= geneLocs[j][k][1] + 5) {
                      exonic = true;
                    }
                  }
                  if (exonic) {
                    if (geneLocs[j][0][1] == 1) {
                      HashVec.addIfAbsent(geneNames[j], plusHits);
                    } else {
                      HashVec.addIfAbsent(geneNames[j], minusHits);
                    }
                  }
                }
              }
              for (int j = 0; j < plusHits.size(); j++) {
                HashVec.addToHashVec(hash, plusHits.elementAt(j),
                                     files[i] + "\t" + Array.toStr(line), false);
              }
              for (int j = 0; j < minusHits.size(); j++) {
                HashVec.addToHashVec(hash, minusHits.elementAt(j),
                                     files[i] + "\t" + Array.toStr(line), false);
              }
            }
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + dir + "siftOutput.dat"
                             + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + dir + "siftOutput.dat" + "\"");
          System.exit(2);
        }
      }
    }

    genes = Collapsed.getKeys(hash);
    favHits = new Hashtable<String, Vector<String>>();
    try {
      writer = new PrintWriter(new FileWriter(dir + "bins.xln"));
      writer.println("Gene\tnumber of indels within an exon");
      for (String gene : genes) {
        v = hash.get(gene);
        writer.println(gene + "\t" + v.size());
        if (Collapsed.indexOfStr(gene, favGenes) >= 0) {
          for (int j = 0; j < v.size(); j++) {
            Collapsed.addToHashVec(favHits, gene, v.elementAt(j), false);
          }
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + "bins.xln");
      e.printStackTrace();
    }

    genes = Collapsed.getKeys(favHits);
    for (String gene : genes) {
      Collapsed.writeList(Collapsed.toStringArray(favHits.get(gene)), dir + gene + ".out");
    }
  }
}
