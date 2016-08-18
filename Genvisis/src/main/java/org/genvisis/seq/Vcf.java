// uses 1000G data from
// ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/supporting/ALL.wgs.project_consensus_vqsr2b.20101123.snps.low_coverage.sites.vcf.gz
// regenerating this file for each subpopulation using the chromosome specific genotype VCF files in
// the previous
package org.genvisis.seq;

import com.google.common.primitives.Ints;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.SerialHash;
import org.genvisis.filesys.VariantList;

public class Vcf {
  public static final String[] BITS_TO_KEEP = {"AC=", "AN=", "DB"};

  private static void lookup(String variantList, String vcfFile, String outfile,
                             String lookupReadyFile, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String temp;
    Hashtable<String, String[]> hash;
    String[][] variants;
    int[] indices;
    String[] alleles;
    int index;
    int[] alleleCounts, numAlleleCalled;
    double[] freqs;
    String[] info, splitInfo;
    boolean dbsnp;
    int missingSplitAlleleCounts;
    String infoToKeep;
    String currentChrom;
    IntVector iv;

    line = null;
    indices = null;
    variants = VariantList.parseList(variantList, log);
    // variants = HashVec.loadFileToStringMatrix(variantList, false, new int[] {0}, false);

    iv = new IntVector();
    for (int i = 0; i < Positions.CHR_CODES.length; i++) {
      if (Files.exists(vcfFile + "." + Positions.CHR_CODES[i] + ".serHash")) {
        iv.add(i);
      }
    }

    if (iv.size() > 0) {
      log.report("Found existing hashed databases from " + vcfFile
                 + " for the following chromomosomes: " + ext.listRanges(Ints.toArray(iv)));
    } else {
      log.report(ext.getTime() + "\tLoading data from : " + vcfFile);
      hash = new Hashtable<String, String[]>();
      try {
        reader = Files.getAppropriateReader(vcfFile);
        do {
          temp = reader.readLine();
        } while (reader.ready() && temp.startsWith("#") && !temp.startsWith("#CHROM"));

        indices = ext.indexFactors(new String[] {"#CHROM", "POS", "REF", "ALT", "ID", "INFO"},
                                   temp.trim().split("[\\s]+"), false, log, true, true);
        currentChrom = "";
        // while (reader.ready()) { // returning false at the same places in the file every time for
        // no reason
        do {
          temp = reader.readLine();
          if (temp != null) {
            line = temp.trim().split("[\\s]+");
            if (line[0].startsWith("#")) {
              log.reportError("Error - ran into some unexpected trailing comments: "
                              + Array.toStr(line));
            } else {
              infoToKeep = "";
              splitInfo = line[indices[5]].split(";");
              for (String element : splitInfo) {
                for (String element2 : BITS_TO_KEEP) {
                  if ((element2.endsWith("=") && element.startsWith(element2))
                      || (element.equals(element2))) {
                    infoToKeep += (infoToKeep.equals("") ? "" : ";") + element;
                  }
                }
              }
              hash.put(line[indices[0]] + "\t" + line[indices[1]],
                       new String[] {line[indices[2]], line[indices[3]], line[indices[4]],
                                     infoToKeep});
              // if (line[indices[1]].equals("193768186")) {
              // System.err.println("hola");
              // System.err.println(Array.toStr(line));
              // System.err.println(reader.readLine());
              // System.err.println(reader.readLine());
              // }

              // aliases = line[indices[4]].split(";");
              // for (int i = 0; i < aliases.length; i++) {
              // hash.put(aliases[i], new String[] {line[indices[2]], line[indices[3]],
              // line[indices[4]], line[indices[5]]});
              // }
            }
          }

          if (temp == null || !line[indices[0]].equals(currentChrom)) {
            log.report(ext.getTime() + "\tSaving hash (n=" + hash.size() + " entries) to " + vcfFile
                       + "." + currentChrom + ".serHash");
            SerialHash.createSerializedStringArrayHash(vcfFile + "." + currentChrom + ".serHash",
                                                       hash);
            hash.clear();
            if (temp != null) {
              currentChrom = line[indices[0]];
              System.out.println(ext.getTime() + "\tchr" + currentChrom);
            }
          }
        } while (temp != null);
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + vcfFile + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + vcfFile + "\"");
        System.exit(2);
      }
      log.report("Last element added: chr" + line[indices[0]] + "\t" + line[indices[1]]);
    }

    try {
      missingSplitAlleleCounts = 0;
      log.report(ext.getTime() + "\tAnnotating variants in " + variantList);
      writer = new PrintWriter(new FileWriter(outfile));
      writer.println("CHROM\tPOS\tREF\tALT\tAltCount\tNumTotalAlleles\tFreq\tNumAlts\tin_dbSNP\trsID");
      currentChrom = "";
      hash = null;
      for (int i = 0; i < variants.length; i++) {
        if (!variants[i][0].equals(currentChrom)) {
          currentChrom = variants[i][0];
          log.report("Loading hash for chr" + currentChrom);
          hash =
              SerialHash.loadSerializedStringArrayHash(vcfFile + "." + currentChrom + ".serHash");
        }
        if (hash.containsKey(variants[i][0] + "\t" + variants[i][1])) {
          line = hash.get(variants[i][0] + "\t" + variants[i][1]);
          // if (hash.containsKey(variants[i][0])) {
          // line = hash.get(variants[i][0]);
          if (!variants[i][2].equalsIgnoreCase(line[0])) {
            log.reportError("Error - mismatched reference allele at position chr" + variants[i][0]
                            + ":" + variants[0][1] + " (found " + variants[i][2] + ", expecting "
                            + line[0] + ")");
          }
          alleles = line[1].split(",");
          alleleCounts = Array.intArray(alleles.length, -1);
          numAlleleCalled = Array.intArray(alleles.length, -1);
          freqs = Array.doubleArray(alleles.length, -2);
          index = -1;
          dbsnp = false;
          for (int j = 0; j < alleles.length; j++) {
            info = line[3].split(";");
            for (String element : info) {
              if (element.startsWith("AC=")) {
                splitInfo = ext.parseStringArg(element, "").split(",");
                if (splitInfo.length != alleleCounts.length) {
                  alleleCounts[j] = Integer.parseInt(splitInfo[0]);
                  missingSplitAlleleCounts++;
                  if (missingSplitAlleleCounts == 1) {
                    log.reportError("Warning - at least one of the markers with multiple alleles listed does not have multiple allele counts listed");
                  }
                } else {
                  alleleCounts[j] = Integer.parseInt(splitInfo[j]);
                }
              }
              if (element.startsWith("AN=")) {
                splitInfo = ext.parseStringArg(element, "").split(",");
                if (splitInfo.length != alleleCounts.length) {
                  numAlleleCalled[j] = Integer.parseInt(splitInfo[0]);
                } else {
                  numAlleleCalled[j] = Integer.parseInt(splitInfo[j]);
                }
              }
              if (element.equals("DB")) {
                dbsnp = true;
              }
            }
            freqs[j] = (double) alleleCounts[j] / (double) numAlleleCalled[j];
            if (variants[i][3].equals(alleles[j])) {
              index = j;
            }
          }
          // writer.println(variants[i][0]+"\t"+Maths.min(freqs[0], 1-freqs[0])+"\t"+dbsnp);
          if (index == -1) {
            log.reportError("Warning - the alternate allele that was seen at position chr"
                            + variants[i][0] + ":" + variants[i][1] + " (" + variants[i][3]
                            + ") was not seen in the reference vcf");
            writer.println(variants[i][0] + "\t" + variants[i][1] + "\t" + variants[i][2] + "\t"
                           + variants[i][3] + "\t.\t.\t.\t" + alleles.length + "\t"
                           + (dbsnp ? "1" : "0") + "\t" + line[2]);
          } else {
            writer.println(variants[i][0] + "\t" + variants[i][1] + "\t" + variants[i][2] + "\t"
                           + variants[i][3] + "\t" + alleleCounts[index] + "\t"
                           + numAlleleCalled[index] + "\t" + freqs[index] + "\t" + alleles.length
                           + "\t" + (dbsnp ? "1" : "0") + "\t" + line[2]);
          }
        } else {
          writer.println(variants[i][0] + "\t" + variants[i][1] + "\t" + variants[i][2] + "\t"
                         + variants[i][3] + "\t.\t.\t.\t.\t.\t.");
        }
      }
      if (missingSplitAlleleCounts > 0) {
        log.reportError("       - there were " + missingSplitAlleleCounts
                        + " sites with multiple alleles that only listed a single allele count");
      }
      writer.close();

    } catch (Exception e) {
      System.err.println("Error writing to " + outfile);
      e.printStackTrace();
    }

    Files.makeVLookupReadyFile(outfile, lookupReadyFile, new int[] {0, 1, 2, 3},
                               new int[] {4, 5, 6, 7, 8, 9});
  }

  public static void createFromParameters(String filename, Logger log) {
    Vector<String> params;

    params =
        Files.parseControlFile(filename, "vcf",
                               new String[] {"list=snps.txt",
                                             "vcf=D:/home/npankrat/NCBI/1000G/ALL.wgs.project_consensus_vqsr2b.20101123.snps.low_coverage.sites.vcf.gz",
                                             "out1=finalProduct.out", "out2=vlookup_table.xln"},
                               log);

    if (params != null) {
      params.add("log=" + log.getFilename());
      main(Array.toStringArray(params));
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String variantList = "list.txt";
    // String vcfFile = "D:/home/npankrat/NCBI/1000G/CEU/CEU.low_coverage.2010_09.sites.vcf.gz";
    String vcfFile =
        "D:/home/npankrat/NCBI/1000G/ALL.wgs.project_consensus_vqsr2b.20101123.snps.low_coverage.sites.vcf.gz";
    String outfile = null;
    String lookupReadyFile = null;
    String logfile = null;
    Logger log;
    String errors = "";

    String usage = "\n" + "seq.Vcf requires 0-4 arguments\n"
                   + "   (1) filename of variants of interest (i.e. list=" + variantList
                   + " (default))\n" + "   (2) VCF file (i.e. vcf=" + vcfFile + " (default))\n"
                   + "   (3) tab-delimited outfile (i.e. out1=[listFile].out (default))\n"
                   + "   (4) ^-delimited VLOOKUP-ready outfile (i.e. out2=[listFile].xln (default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("list=")) {
        variantList = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("vcf=")) {
        vcfFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out1=")) {
        outfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out2=")) {
        lookupReadyFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        errors += "Error - invalid argument: " + arg + "\n";
      }
    }
    log = new Logger(logfile);
    if (numArgs != 0) {
      log.report(errors);
      log.report(usage);
      return;
    }
    try {
      if (outfile == null) {
        outfile = variantList + ".out";
      }
      if (lookupReadyFile == null) {
        lookupReadyFile = variantList + ".xln";
      }
      lookup(variantList, vcfFile, outfile, lookupReadyFile, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
