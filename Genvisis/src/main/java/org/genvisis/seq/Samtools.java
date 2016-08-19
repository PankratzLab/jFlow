package org.genvisis.seq;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class Samtools {
  public static final String[] COUNT_READS_PER_CHR =
                                                   {"echo `date`",
                                                    "echo \"Counting number of reads per chromosome for file $1\"",
                                                    "samtools view $1 | awk -F\"\t\" '{print $3}' | awk '",
                                                    "", "BEGIN {", "   pat1 = 0;", "   pat2 = 0;",
                                                    "   pat3 = 0;", "   pat4 = 0;", "   pat5 = 0;",
                                                    "   pat6 = 0;", "   pat7 = 0;", "   pat8 = 0;",
                                                    "   pat9 = 0;", "   pat10 = 0;",
                                                    "   pat11 = 0;", "   pat12 = 0;",
                                                    "   pat13 = 0;", "   pat14 = 0;",
                                                    "   pat15 = 0;", "   pat16 = 0;",
                                                    "   pat17 = 0;", "   pat18 = 0;",
                                                    "   pat19 = 0;", "   pat20 = 0;",
                                                    "   pat21 = 0;", "   pat22 = 0;",
                                                    "   patX = 0;", "   patY = 0;", "   patXY = 0;",
                                                    "   patM = 0;", "   patUnmapped = 0;", "}", "",
                                                    "/chr1$/ { pat1 += 1;}",
                                                    "/chr2$/ { pat2 += 1;}",
                                                    "/chr3$/ { pat3 += 1;}",
                                                    "/chr4$/ { pat4 += 1;}",
                                                    "/chr5$/ { pat5 += 1;}",
                                                    "/chr6$/ { pat6 += 1;}",
                                                    "/chr7$/ { pat7 += 1;}",
                                                    "/chr8$/ { pat8 += 1;}",
                                                    "/chr9$/ { pat9 += 1;}",
                                                    "/chr10$/ { pat10 += 1;}",
                                                    "/chr11$/ { pat11 += 1;}",
                                                    "/chr12$/ { pat12 += 1;}",
                                                    "/chr13$/ { pat13 += 1;}",
                                                    "/chr14$/ { pat14 += 1;}",
                                                    "/chr15$/ { pat15 += 1;}",
                                                    "/chr16$/ { pat16 += 1;}",
                                                    "/chr17$/ { pat17 += 1;}",
                                                    "/chr18$/ { pat18 += 1;}",
                                                    "/chr19$/ { pat19 += 1;}",
                                                    "/chr20$/ { pat20 += 1;}",
                                                    "/chr21$/ { pat21 += 1;}",
                                                    "/chr22$/ { pat22 += 1;}",
                                                    "/chrX$/ { patX += 1;}",
                                                    "/chrY$/ { patY += 1;}",
                                                    "/chrXY$/ { patXY += 1;}",
                                                    "/chrM$/ { patM += 1;}",
                                                    "/\\*$/ { patUnmapped += 1;}", "", "END {",
                                                    "   print \"chr1\", pat1;",
                                                    "   print \"chr2\", pat2;",
                                                    "   print \"chr3\", pat3;",
                                                    "   print \"chr4\", pat4;",
                                                    "   print \"chr5\", pat5;",
                                                    "   print \"chr6\", pat6;",
                                                    "   print \"chr7\", pat7;",
                                                    "   print \"chr8\", pat8;",
                                                    "   print \"chr9\", pat9;",
                                                    "   print \"chr10\", pat10;",
                                                    "   print \"chr11\", pat11;",
                                                    "   print \"chr12\", pat12;",
                                                    "   print \"chr13\", pat13;",
                                                    "   print \"chr14\", pat14;",
                                                    "   print \"chr15\", pat15;",
                                                    "   print \"chr16\", pat16;",
                                                    "   print \"chr17\", pat17;",
                                                    "   print \"chr18\", pat18;",
                                                    "   print \"chr19\", pat19;",
                                                    "   print \"chr20\", pat20;",
                                                    "   print \"chr21\", pat21;",
                                                    "   print \"chr22\", pat22;",
                                                    "   print \"chrX\", patX;",
                                                    "   print \"chrY\", patY;",
                                                    "   print \"chrXY\", patXY;",
                                                    "   print \"chrM\", patM;",
                                                    "   print \"*\", patUnmapped;", "}'",
                                                    "echo `date`",};

  public static void determineSex(String dir, int filesPerBatch, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String temp;
    Vector<String> v;
    String[] files;
    String pwd;
    String xCount, yCount;
    int autosomalCount;


    pwd = ext.pwd();
    dir = ext.verifyDirFormat(dir);

    if (!Files.exists(pwd + "samCountReadsPerChr")) {
      Files.writeList(COUNT_READS_PER_CHR, pwd + "samCountReadsPerChr");
      Files.chmod(pwd + "samCountReadsPerChr");
    }

    files = Files.list(dir, null, ".bam", false, false);
    if (files.length == 0) {
      log.reportError("There were no .bam files to be found in " + dir);
    }

    v = new Vector<String>();
    try {
      writer = new PrintWriter(new FileWriter(pwd + "sexPlot.xln"));
      writer.println("File\tnumXreads\tnumYreads\tnumAutosomalReads\txProp\tyProp\txExpRatio\tyExpRatio");
      for (String file : files) {
        if (Files.exists(pwd + file + ".chrCounts")) {
          try {
            xCount = "missing";
            yCount = "missing";
            autosomalCount = 0;
            reader = Files.getAppropriateReader(pwd + file + ".chrCounts");
            while (reader.ready()) {
              temp = reader.readLine();
              line = temp.split("[\\s]+");
              for (int chr = 1; chr <= 22; chr++) {
                if (line[0].equals("chr" + chr)) {
                  if (ext.isMissingValue(line[1])) {
                    log.reportError("Error - invalid read count for " + file + " for chr" + chr
                                    + ": " + line[1]);
                  } else {
                    autosomalCount += Integer.parseInt(line[1]);
                  }

                }
              }
              if (line[0].equals("chrX")) {
                if (ext.isValidInteger(line[1])) {
                  xCount = line[1];
                } else {
                  log.reportError("Error - invalid read count for " + file + " for chrX: "
                                  + line[1]);
                  xCount = "-1";
                }
              }
              if (line[0].equals("chrY")) {
                if (ext.isValidInteger(line[1])) {
                  yCount = line[1];
                } else {
                  log.reportError("Error - invalid read count for " + file + " for chrY: "
                                  + line[1]);
                  yCount = "-1";
                }
              }
              if (line.length != 2 && !temp.startsWith("Counting number of reads")) {
                if (line.length > 3 && line[3].indexOf(":") > 0) {
                  // timestamp
                } else {
                  log.reportError("Error message for '" + file + "': " + temp);
                }
              }
            }
            reader.close();
            writer.println(file + "\t" + xCount + "\t" + yCount + "\t" + autosomalCount + "\t"
                           + ext.formDeci(Double.parseDouble(xCount) / autosomalCount, 6) + "\t"
                           + ext.formDeci(Double.parseDouble(yCount) / autosomalCount, 6) + "\t"
                           + ext.formDeci(Double.parseDouble(xCount) / autosomalCount * 21.75, 4)
                           + "\t"
                           + ext.formDeci(Double.parseDouble(yCount) / autosomalCount * 327, 4));
          } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \"" + file + ".chrCounts" + "\" not found in directory"
                               + pwd);
            System.exit(1);
          } catch (IOException ioe) {
            System.err.println("Error reading file \"" + file + ".chrCounts" + "\"");
            System.exit(2);
          }
        } else {
          Files.write(pwd + "samCountReadsPerChr " + dir + file + " 1>> " + pwd + file
                      + ".chrCounts 2>> " + pwd + file + ".chrCounts",
                      "batches/" + file + ".batch");
          Files.chmod(pwd + "batches/" + file + ".batch");
          v.add(pwd + "batches/" + file + ".batch");
        }

      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + "sexPlot.xln");
      e.printStackTrace();
    }

    if (v.size() > 0) {
      log.report("There are " + v.size() + " .bam files remaining to be counted");
      Files.qsubMultiple(v, null, pwd + "batches/", "countPerChr", filesPerBatch, true, null, -1,
                         500, 0.5);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = "";
    String logfile = null;
    Logger log;
    boolean getSex = false;
    int filesPerBatch = 8;

    String usage = "\n" + "seq.Samtools requires 0-1 arguments\n" + "   (1) directory (i.e. dir="
                   + dir + " (default))\n"
                   + "   (2) determine sex for all bam files in directory (i.e. -getSex (not the default))\n"
                   + "   (3) number of files per batch (i.e. filesPerBatch=" + filesPerBatch
                   + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-getSex")) {
        getSex = true;
        numArgs--;
      } else if (arg.startsWith("filesPerBatch=")) {
        filesPerBatch = ext.parseIntArg(arg);
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
      if (getSex) {
        determineSex(dir, filesPerBatch, log);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
