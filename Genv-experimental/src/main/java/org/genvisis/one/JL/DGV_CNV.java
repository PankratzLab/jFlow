package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;

/**
 * @author lane0212 Parse file such as those listed here
 *         http://dgv.tcag.ca/dgv/app/downloads?ref=GRCh37/hg19
 */
public class DGV_CNV {
  private static final String[] HEADER_DGV = new String[] {"variantaccession", "chr", "start",
      "end", "varianttype", "variantsubtype", "reference", "pubmedid", "method", "platform",
      "mergedvariants", "supportingvariants", "mergedorsample", "frequency", "samplesize",
      "observedgains", "observedlosses", "cohortdescription", "genes", "samples"};
  private static final String[] COPY_NUMBER_VARIATION_MAP = new String[] {"chr", "start", "end",
      "state", "id", "type", "num_variants", "num_samples", "num_samples_multicounted",
      "num_studies", "variants", "samples", "studies", "African", "Asian", "European", "Mexican",
      "Middle_East", "Native_American", "Oceania", "South_American"};

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = "C:/bin/ref/1000GCNV/";

    String usage = "\n" + "one.JL.DGV_CNV requires 0-1 arguments\n" + "   (1) directory (i.e. dir="
        + dir + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
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
      Logger log = new Logger(dir + "dgvCNV_parser.log");
      parseToCNVs(dir, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parseToCNVs(String DGVDir, Logger log) {
    Hashtable<String, int[]> copyHash = new Hashtable<String, int[]>();
    copyHash.put("loss", new int[] {1});
    copyHash.put("duplication", new int[] {3});
    copyHash.put("gain", new int[] {3});
    copyHash.put("deletion", new int[] {1});
    copyHash.put("gain+loss", new int[] {1, 3});
    copyHash.put("insertion", new int[] {3});
    copyHash.put("Gain", new int[] {3});
    copyHash.put("Loss", new int[] {1});

    copyHash.put("complex", new int[] {2});

    copyHash.put("novel", new int[] {2});
    copyHash.put("inversion", new int[] {2});
    copyHash.put("mobile", new int[] {2});
    copyHash.put("sequence", new int[] {2});
    copyHash.put("tandem", new int[] {2});
    Hashtable<String, String> skipTypes = new Hashtable<String, String>();
    skipTypes.put("inversion", "inversion");
    skipTypes.put("tandem duplication", "tandem duplication");
    skipTypes.put("insertion", "insertion");
    skipTypes.put("mobile element insertion", "mobile element insertion");

    String[] filesToParse = Files.list(DGVDir, null, ".txt", true, false, true);
    for (String element : filesToParse) {
      String out = element + ".cnv";

      try {
        PrintWriter writer = new PrintWriter(new FileWriter(out));
        writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
        boolean dgv = true;
        BufferedReader reader = Files.getAppropriateReader(element);
        reader.readLine();

        int[] indices =
            ext.indexFactors(HEADER_DGV, Files.getHeaderOfFile(element, log), true, false);
        if (Array.countIf(indices, -1) > 0) {
          indices = ext.indexFactors(COPY_NUMBER_VARIATION_MAP, Files.getHeaderOfFile(element, log),
              true, false);
          dgv = false;
        }
        int lineNum = 0;
        // String[] header = Files.getHeaderOfFile(filesToParse[i], log);
        while (reader.ready()) {
          lineNum++;
          String[] line = reader.readLine().trim().split("\t");
          if (!dgv) {
            int chr = Positions.chromosomeNumber(line[indices[0]]);
            int start = Integer.parseInt(line[indices[1]]);
            int stop = Integer.parseInt(line[indices[2]]);
            String type = line[indices[5]];
            int[] cn = copyHash.get(type);
            if (cn.length != 1) {
              log.reportTimeError("what type " + type);
              System.exit(1);
            }

            int numVar = Integer.parseInt(line[indices[6]]);
            int doubleSamp = Integer.parseInt(line[indices[8]]);
            int numSamp = Integer.parseInt(line[indices[7]]);
            String[] samples = line[indices[11]].split(",");
            if (numSamp != samples.length) {
              System.out.println(numSamp + "\t" + samples.length);
              System.exit(1);
            }
            int effectiveNumVar = numVar - doubleSamp;
            for (int j = 0; j < effectiveNumVar; j++) {
              if (j >= samples.length) {
                break;
              } else {
                CNVariant tmp = new CNVariant(samples[j], samples[j], (byte) chr, start, stop,
                    cn[0], 99, -1, 99);
                writer.println(tmp.toPlinkFormat());
              }

            }
            // System.out.println(header[indices[6]] + "\t" + header[indices[7]] + "\t" +
            // header[indices[8]]);
            // System.out.println(numVar + "\t" + effectiveNumVar + "\t" + doubleSamp + "\t" +
            // numSamp + "\t" + samples.length);

          }

          else if (dgv) {

            int chr = Positions.chromosomeNumber(line[indices[1]]);
            int start = Integer.parseInt(line[indices[2]]);
            int stop = Integer.parseInt(line[indices[3]]);
            String type = line[indices[5]];
            String reference = line[indices[6]];
            if (indices[19] < line.length && !skipTypes.containsKey(type)
                && (reference.startsWith("1000_Genomes") || reference.startsWith("Conrad_"))) {

              try {
                String[] samples = line[indices[19]].split(",");
                int numDup = Integer.parseInt(line[indices[15]]);
                int numDel = Integer.parseInt(line[indices[16]]);
                if (numDel + numDup != samples.length && !type.equals("gain+loss")) {// gain plus
                                                                                     // loss does
                                                                                     // not equal
                                                                                     // sample size
                  System.out.println(Array.toStr(line));
                  System.out.println(lineNum + "\t" + line.length + "\t" + element);
                  log.reportTimeError("Del: " + numDel + " and Dup: " + numDup
                      + " does not add up to " + samples.length);
                  System.exit(1);
                } else {
                  int sampIndex = 0;
                  for (int j = 0; j < numDup; j++) {
                    CNVariant tmp = new CNVariant(samples[sampIndex], samples[sampIndex],
                        (byte) chr, start, stop, 3, 99, -1, 99);
                    writer.println(tmp.toPlinkFormat());
                    sampIndex++;
                    if (sampIndex >= samples.length) {
                      break;
                    }
                  }
                  if (sampIndex >= samples.length) {
                    sampIndex = Math.max(0, sampIndex - numDel);
                  }
                  for (int j = 0; j < numDel; j++) {
                    CNVariant tmp = new CNVariant(samples[sampIndex], samples[sampIndex],
                        (byte) chr, start, stop, 1, 99, -1, 99);
                    writer.println(tmp.toPlinkFormat());
                    sampIndex++;
                    if (sampIndex >= samples.length) {
                      break;
                    }
                  }
                }
              } catch (ArrayIndexOutOfBoundsException e) {
                System.out.println(Array.toStr(line));
                System.out.println(lineNum + "\t" + line.length + "\t" + element);
                log.reportException(e);
                e.printStackTrace();
                System.exit(1);
              }
            }
          }
        }

        reader.close();
        writer.close();

      } catch (FileNotFoundException fnfe) {
        log.reportError("Error: file \"" + element + "\" not found in current directory");
        return;
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + element + "\"");
        return;
      }
      System.out.println(element);

    }
  }

  public static void parseToCNVss(String DGVDir, Logger log) {
    Hashtable<String, int[]> copyHash = new Hashtable<String, int[]>();
    copyHash.put("loss", new int[] {1});
    copyHash.put("duplication", new int[] {3});
    copyHash.put("gain", new int[] {3});
    copyHash.put("deletion", new int[] {1});
    copyHash.put("gain+loss", new int[] {1, 3});
    copyHash.put("insertion", new int[] {3});
    copyHash.put("Gain", new int[] {3});
    copyHash.put("Loss", new int[] {1});

    copyHash.put("complex", new int[] {2});

    copyHash.put("novel", new int[] {2});
    copyHash.put("inversion", new int[] {2});
    copyHash.put("mobile", new int[] {2});
    copyHash.put("sequence", new int[] {2});
    copyHash.put("tandem", new int[] {2});

    String[] filesToParse = Files.list(DGVDir, null, ".txt", true, false, true);
    log.reportTimeInfo("Found " + filesToParse.length + " files to parse");
    for (String element : filesToParse) {
      String out = element + ".cnv";
      try {
        PrintWriter writer = new PrintWriter(new FileWriter(out));
        writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
        int[] indices =
            ext.indexFactors(HEADER_DGV, Files.getHeaderOfFile(element, log), true, false);
        try {
          BufferedReader reader = Files.getAppropriateReader(element);
          reader.readLine();
          int lineNum = 0;
          while (reader.ready()) {
            String[] line = reader.readLine().trim().split("[\\s]+");
            lineNum++;
            int chr = Positions.chromosomeNumber(line[indices[1]]);
            int start = Integer.parseInt(line[indices[2]]);
            int stop = Integer.parseInt(line[indices[3]]);
            String type = indices[5] < 0 ? line[indices[6]] : line[indices[5]];

            if (!copyHash.containsKey(type)) {

              log.reportTimeError("what type is " + type);
              System.exit(1);
            } else {
              int[] types = null;
              if (type.equals("gain+loss")) {
                int numDel = Integer.parseInt(line[indices[16]]);
                int numGain = Integer.parseInt(line[indices[17]]);
                types = new int[numDel + numGain];
                for (int j = 0; j < numDel; j++) {
                  types[j] = 1;
                }
                for (int j = numDel; j < numGain; j++) {
                  types[j] = 3;
                }

              } else {
                types = copyHash.get(type);
              }

              try {
                String[] samples = null;

                if (indices[20] < line.length) {
                  samples = line[indices[20]].split(",");

                  // else {
                  // samples = new String[] { "No_Sample_" + lineNum };
                  // }

                  for (String sample : samples) {
                    for (int type2 : types) {
                      CNVariant tmp =
                          new CNVariant(sample, sample, (byte) chr, start, stop, type2, 99, -1, -1);
                      writer.println(tmp.toPlinkFormat());
                    }
                  }
                }

              } catch (Exception e) {
                log.reportTimeError(Array.toStr(line));
                log.reportException(e);
                log.reportTimeError(element);
                log.reportTimeError(lineNum + "");
                System.exit(1);

              }

            }

          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          log.reportError("Error: file \"" + element + "\" not found in current directory");
          writer.close();
          return;
        } catch (IOException ioe) {
          log.reportError("Error reading file \"" + element + "\"");
          writer.close();
          return;
        }
        System.out.println(element);

        writer.close();
      } catch (Exception e) {
        log.reportError("Error writing to " + out);
        log.reportException(e);
      }
    }
  }

}
