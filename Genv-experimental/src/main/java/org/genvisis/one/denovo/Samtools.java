package org.genvisis.one.denovo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

public class Samtools {

  public static void extractRegions(String outputFilename, String fullPathToSuperNovoPhase2Output,
                                    String fullPathToTrioNameList, int windowInBp, Logger log) {
    Hashtable<String, Vector<String>> hash = new Hashtable<>();
    String[][] bamFilenamesByTrios;

    hash = loadFromFile(fullPathToSuperNovoPhase2Output, log);
    bamFilenamesByTrios = SuperNovo.loadNamesFromList(fullPathToTrioNameList);
    writerToFile(outputFilename, null, hash, bamFilenamesByTrios, windowInBp, log);
  }

  public static Hashtable<String, Vector<String>> loadFromFile(String filename, Logger log) {
    BufferedReader reader;
    String[] line;
    Hashtable<String, Vector<String>> hash;
    Vector<String> v;

    hash = new Hashtable<>();
    try {
      reader = new BufferedReader(new FileReader(filename));
      while (reader.ready()) {
        line = reader.readLine().split("\t");
        if (hash.containsKey(line[2])) {
          v = hash.get(line[2]);
        } else {
          v = new Vector<>();
        }
        v.add(line[4]);
      }
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }

    return hash;
  }

  public static void writerToFile(String filename, String xmlDir,
                                  Hashtable<String, Vector<String>> hash,
                                  String[][] samFilenamesByTrios, int windowInBp, Logger log) {
    PrintWriter writer;
    Vector<String> v;
    int index;
    Set<String> keySet;
    String chr;
    int pos;

    keySet = hash.keySet();
    try {
      for (String key : keySet) {
        chr = key.split(":")[0];
        pos = Integer.parseInt(key.split(":")[1]);
        v = hash.get(key);
        for (int i = 0; i < v.size(); i++) {
          index = getIndex(v.elementAt(i), samFilenamesByTrios, log);

          writer = Files.openAppropriateWriter(xmlDir + ".bamScript");
          writer.println("samtools view " + samFilenamesByTrios[index][1] + "-b " + chr + ":"
                         + (pos - windowInBp) + "-" + (pos + windowInBp) + " > "
                         + samFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos + "_C");
          writer.close();

          writer = Files.openAppropriateWriter(filename + ".bamScript");
          writer.println("samtools view " + samFilenamesByTrios[index][2] + "-b " + chr + ":"
                         + (pos - windowInBp) + "-" + (pos + windowInBp) + " > "
                         + samFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos + "_D");
          writer.close();

          writer = Files.openAppropriateWriter(filename + ".bamScript");
          writer.println("samtools view " + samFilenamesByTrios[index][3] + "-b " + chr + ":"
                         + (pos - windowInBp) + "-" + (pos + windowInBp) + " > "
                         + samFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos + "_M");
          writer.close();

          writer = Files.openAppropriateWriter(filename + ".IGV_Script");
          writer.println(getIgvXmlScript(xmlDir, chr, pos + "", samFilenamesByTrios[index]));
          writer.close();
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static int getIndex(String trioId, String[][] samFilenamesByTrios, Logger log) {
    int index;

    index = -1;
    for (int i = 0; i < samFilenamesByTrios.length; i++) {
      if (samFilenamesByTrios[i][0].equals(trioId)) {
        index = i;
        break;
      }
    }

    return index;
  }

  public static String getIgvXmlScript(String miniSamDir, String chr, String pos,
                                       String[] miniSamFilenamesOfOneTrio) {
    return "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>"
           + "\n<Session genome=\"hg19\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"chr"
           + chr + ":" + pos + "\" version=\"8\">" + "\n<Resources>" + "\n<Resource path=\""
           + miniSamDir + miniSamFilenamesOfOneTrio[1] + "\"/>" + "\n<Resource path=\"" + miniSamDir
           + miniSamFilenamesOfOneTrio[2] + "\"/>" + "\n<Resource path=\"" + miniSamDir
           + miniSamFilenamesOfOneTrio[3] + "\"/>" + "\n</Resources>" + "</Session>";
  }

  public static String getIgvLaunchScript(String fulPathToXml) {
    return "java -Xmx1200m -Dproduction=true -Djava.net.preferIPv4Stack=true -Dsun.java2d.noddraw=true -jar D:/logan/DeNovos/IGV/IGV_2.3.36/igv.jar "
           + fulPathToXml;
  }

  /**
   * @param scriptFileName
   * @param trioIdChrPos trioId + chr + pos
   * @param samFilenamesByTrios
   * @param windowInBp
   * @param log
   */
  public static void saveScriptsGeneratingMiniSamsOfOneGene(String scriptFileName,
                                                            String miniSamDir,
                                                            List<String> trioIdChrPos,
                                                            String[][] samFilenamesByTrios,
                                                            int windowInBp, Logger log) {
    PrintWriter writer;
    int index;
    String[] line;
    // String chr;
    int pos;

    try {
      writer = Files.openAppropriateWriter(scriptFileName);
      for (String key : trioIdChrPos) {
        line = key.split("\t");
        // chr = line[1];
        pos = Integer.parseInt(line[2]);
        index = getIndex(line[0], samFilenamesByTrios, log);

        if (index < 0) {
          log.reportError("Trio ID (" + line[0]
                          + ") from a phase 1 output file does not match any row of the Trio ID List file, - skipped generating some scripts for mini bam files");
        } else {
          writer.print("samtools view " + samFilenamesByTrios[index][1] + " -b chr" + line[1] + ":"
                       + Math.max(0, pos - windowInBp) + "-" + (pos + windowInBp) + " > "
                       + miniSamDir + samFilenamesByTrios[index][0] + "_chr" + line[1] + "_"
                       + line[2] + "_C.bam\n" + "samtools index " + miniSamDir
                       + samFilenamesByTrios[index][0] + "_chr" + line[1] + "_" + line[2]
                       + "_C.bam\n" + "samtools view " + samFilenamesByTrios[index][2] + " -b chr"
                       + line[1] + ":" + Math.max(0, pos - windowInBp) + "-" + (pos + windowInBp)
                       + " > " + miniSamDir + samFilenamesByTrios[index][0] + "_chr" + line[1] + "_"
                       + line[2] + "_D.bam\n" + "samtools index " + miniSamDir
                       + samFilenamesByTrios[index][0] + "_chr" + line[1] + "_" + line[2]
                       + "_D.bam\n" + "samtools view " + samFilenamesByTrios[index][3] + " -b chr"
                       + line[1] + ":" + Math.max(0, pos - windowInBp) + "-" + (pos + windowInBp)
                       + " > " + miniSamDir + samFilenamesByTrios[index][0] + "_chr" + line[1] + "_"
                       + line[2] + "_M.bam\n" + "samtools index " + miniSamDir
                       + samFilenamesByTrios[index][0] + "_chr" + line[1] + "_" + line[2]
                       + "_M.bam\n");
        }
      }

      writer.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /**
   * @param scriptFileName
   * @param trioIdChrPos trioId + chr + pos
   * @param samFilenamesByTrios
   * @param windowInBp
   * @param log
   */
  public static void saveScriptsGeneratingMiniSamsAllGenesAtOnce(String scriptFileName,
                                                                 String miniSamDir,
                                                                 Hashtable<String, Vector<String>> miniSamNeeded,
                                                                 String[][] samFilenamesByTrios,
                                                                 int windowInBp, Logger log) {
    PrintWriter writer;
    int index;
    Vector<String> trioIdChrPos;
    String[] line;
    // String chr;
    int pos;

    try {
      writer = Files.openAppropriateWriter(scriptFileName);
      for (String gene : miniSamNeeded.keySet()) {
        trioIdChrPos = miniSamNeeded.get(gene);
        for (String key : trioIdChrPos) {
          line = key.split("\t");
          // chr = line[1];
          pos = Integer.parseInt(line[2]);
          index = getIndex(line[0], samFilenamesByTrios, log);

          if (index < 0) {
            log.reportError("Trio ID (" + line[0]
                            + ") from a phase 1 output file does not match any row of the Trio ID List file, - skipped generating some scripts for mini bam files");
          } else {
            writer.print("samtools view " + samFilenamesByTrios[index][1] + " -b chr" + line[1]
                         + ":" + Math.max(0, pos - windowInBp) + "-" + (pos + windowInBp) + " > "
                         + miniSamDir + samFilenamesByTrios[index][0] + "_chr" + line[1] + "_"
                         + line[2] + "_C.bam\n" + "samtools index " + miniSamDir
                         + samFilenamesByTrios[index][0] + "_chr" + line[1] + "_" + line[2]
                         + "_C.bam\n" + "samtools view " + samFilenamesByTrios[index][2] + " -b chr"
                         + line[1] + ":" + Math.max(0, pos - windowInBp) + "-" + (pos + windowInBp)
                         + " > " + miniSamDir + samFilenamesByTrios[index][0] + "_chr" + line[1]
                         + "_" + line[2] + "_D.bam\n" + "samtools index " + miniSamDir
                         + samFilenamesByTrios[index][0] + "_chr" + line[1] + "_" + line[2]
                         + "_D.bam\n" + "samtools view " + samFilenamesByTrios[index][3] + " -b chr"
                         + line[1] + ":" + Math.max(0, pos - windowInBp) + "-" + (pos + windowInBp)
                         + " > " + miniSamDir + samFilenamesByTrios[index][0] + "_chr" + line[1]
                         + "_" + line[2] + "_M.bam\n" + "samtools index " + miniSamDir
                         + samFilenamesByTrios[index][0] + "_chr" + line[1] + "_" + line[2]
                         + "_M.bam\n");
          }
        }
      }
      writer.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void saveScriptsLaunchingIgv(String igvScriptDir, String miniSamDir,
                                             List<String> trioIdChrPos,
                                             String[][] samFilenamesByTrios, int windowInBp,
                                             Logger log) {
    int index;
    String[] line;

    for (String key : trioIdChrPos) {
      line = key.split("\t");
      index = getIndex(line[0], samFilenamesByTrios, log);

      if (index < 0) {
        log.reportError("Trio ID (" + line[0]
                        + ") from a phase 1 output file does not match any row of the Trio ID List file, - skipped generating some scripts for mini bam files");
      } else {
        if (!new File(igvScriptDir + line[0] + "_" + line[1] + "_" + line[2] + ".xml").exists()) {
          Files.write(getIgvXmlScript(miniSamDir, line[1], line[2],
                                      new String[] {line[0],
                                                    samFilenamesByTrios[index][0]
                                                             + "_chr" + line[1] + "_" + line[2]
                                                             + "_C.bam",
                                                    samFilenamesByTrios[index][0]
                                                                         + "_chr" + line[1] + "_"
                                                                         + line[2] + "_D.bam",
                                                    samFilenamesByTrios[index][0] + "_chr" + line[1]
                                                                                               + "_"
                                                                                               + line[2]
                                                                                               + "_M.bam"}),
                      igvScriptDir + line[0] + "_chr" + line[1] + "_" + line[2] + ".xml");
          Files.write(getIgvLaunchScript(igvScriptDir + line[0] + "_chr" + line[1] + "_" + line[2]
                                         + ".xml"),
                      igvScriptDir + line[0] + "_chr" + line[1] + "_" + line[2] + ".bat");
        }
      }
    }
  }

  public static HashSet<String> listFiles(String miniSamDir, Logger log) {
    String[] filenames;
    boolean[] isVisited;
    String[] line1;
    String[] line2;
    HashSet<String> result;
    Vector<String> currentTrio;

    result = new HashSet<>();
    filenames = Files.list(miniSamDir, ".bam");
    isVisited = new boolean[filenames.length];
    for (int i = 0; i < filenames.length; i++) {
      if (!isVisited[i]) {
        currentTrio = new Vector<>(3);
        currentTrio.add(filenames[i]);
        line1 = filenames[i].split("_");
        for (int j = 0; j < filenames.length; j++) {
          line2 = filenames[j].split("_");
          if (!isVisited[j] && j != i && line2[0].equals(line1[0]) && line2[1].equals(line1[1])
              && line2[2].equals(line1[2])) {
            isVisited[j] = true;
            currentTrio.add(filenames[j]);
          }
        }
        if (currentTrio.size() == 3) {
          result.add(line1[0] + "\t" + line1[1].substring(line1[1].indexOf("chr") + 3) + "\t"
                     + line1[2]);
          if (!new File(miniSamDir + line1[0] + "_" + line1[1] + "_" + line1[2]
                        + ".xml").exists()) {
            Files.write(getIgvXmlScript(miniSamDir, line1[1], line1[2],
                                        new String[] {line1[0], currentTrio.elementAt(0),
                                                      currentTrio.elementAt(1),
                                                      currentTrio.elementAt(2)}),
                        miniSamDir + line1[0] + "_" + line1[1] + "_" + line1[2] + ".xml");
          }
        }
      }
    }

    return result;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "Samtools.dat";
    String fullPathToSuperNovoPhase2Output = "N:/statgen/OS_Logan/SuperNovo/output/SuperNovo_summary.txt";
    String fullPathToTrioNameList = "Samtools.dat";
    int windowInBp = 5000;
    String logfile = null;
    Logger log;

    String usage = "\n" + "bioinformatics.Samtools requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("window=")) {
        windowInBp = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = ext.parseStringArg(arg, null);
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
      extractRegions(filename, fullPathToSuperNovoPhase2Output, fullPathToTrioNameList, windowInBp,
                     log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
