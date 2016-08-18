package org.genvisis.affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Vector;

import org.genvisis.cnv.manage.SourceFileParser;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class MergeChp implements Runnable {
  public static final String[][] SNP_HEADER_OPTIONS = {{"SNP", "rsID", "ProbeSetName"}};
  public static final String FILENAME_AS_ID_OPTION = "[FILENAME_ROOT]";

  private final String affyResultsDir;
  private final String[] files;
  // private long timeBegan;
  // private int threadId;
  private final String commonSubFolderPattern;
  private final String output;
  private final Logger log;

  public MergeChp(String affyResultsDir, String[] files, long timeBegan) {
    this(affyResultsDir, files, timeBegan, -1, "", "", new Logger());
  }

  public MergeChp(String affyResultsDir, String[] files, long timeBegan, int threadId,
                  String commonSubFolderPattern, String output, Logger log) {
    this.affyResultsDir = affyResultsDir;
    this.files = files;
    this.log = log;
    // this.timeBegan = timeBegan;
    // this.threadId = threadId;
    this.commonSubFolderPattern = commonSubFolderPattern;
    this.output = output;
  }

  @Override
  public void run() {

    BufferedReader reader;
    PrintWriter writer;
    String delimiter;
    // idHeader ,
    // idHeader = proj.getProperty(Project.ID_HEADER);
    delimiter = "\t";
    String[] line;
    String aline;
    // common folder output by apt-genotype, present in each subdirectory of Source
    String commonSubFolder = commonSubFolderPattern;
    // check source directory
    String[] dirList = Files.listDirectories(affyResultsDir, false);
    if (!affyResultsDir.equals("") && !new File(affyResultsDir).exists()) {
      log.reportTimeError("the Project source location is invalid: " + affyResultsDir);
      return;
    }

    int counts = 0;

    for (int j = 0; j < files.length; j++) {
      String outputSamp = output + files[j] + ".gz";
      if (!Files.exists(outputSamp)) {
        writer = Files.getAppropriateWriter(outputSamp);
        log.reportTimeInfo("merging files " + (j + 1) + " of " + files.length);

        for (int i = 0; i < dirList.length; i++) {
          try {
            reader = Files.getAppropriateReader(affyResultsDir + dirList[i] + commonSubFolder + "/"
                                                + files[j]);
            // filter comments
            do {
              line = reader.readLine().trim().split(delimiter, -1);

            } while (reader.ready() && (ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line,
                                                         false, true, false, false)[0] == -1));
            // if its the first directory, print the header

            if (i == 0) {
              writer.println(Array.toStr(line));
            }

            while (reader.ready()) {
              aline = reader.readLine().trim();
              writer.println(aline);
              counts++;
            }
            reader.close();

          } catch (FileNotFoundException fnfe) {
            log.reportTimeError("Error: file \"" + files[j] + "\" not found in " + affyResultsDir
                                + dirList[i]);
            return;
          } catch (IOException ioe) {
            log.reportTimeError("Error reading file \"" + affyResultsDir + dirList[i] + files[j]
                                + "\"");
            return;
          }
        }
        writer.close();
        log.reportTimeInfo("Merged file contains " + counts + " IDs");
        counts = 0;
      } else {
      }
    }
  }

  public static void combineChpFiles(String affyResultsDir, int numThreads,
                                     String commonSubFolderPattern, String commonFilename,
                                     String output, Logger log) {
    Vector<Vector<String>> fileCabinet;
    Thread[] threads;
    boolean complete;
    long timeBegan;
    // common folder output by apt-genotype, present in each subdirectory of Source
    String commonSubFolder = commonSubFolderPattern;
    // check source directory
    timeBegan = new Date().getTime();
    if (!affyResultsDir.equals("") && !new File(affyResultsDir).exists()) {
      System.err.println("Error - the Project source location is invalid: " + affyResultsDir);
      return;
    }

    String[] dirList = Files.listDirectories(affyResultsDir, false);
    String[] files =
        Files.list(affyResultsDir + dirList[0] + commonSubFolder, commonFilename, false);
    fileCabinet = new Vector<Vector<String>>();
    for (int i = 0; i < numThreads; i++) {
      fileCabinet.add(new Vector<String>());
    }
    for (int i = 0; i < files.length; i++) {
      fileCabinet.elementAt(i % numThreads).add(files[i]);

    }
    System.out.println("beginning to merge " + files.length + " files");
    threads = new Thread[numThreads];
    for (int i = 0; i < numThreads; i++) {
      threads[i] = new Thread(new MergeChp(affyResultsDir,
                                           fileCabinet.elementAt(i)
                                                      .toArray(new String[fileCabinet.elementAt(i)
                                                                                     .size()]),
                                           timeBegan, i, commonSubFolder, output, log));
      threads[i].start();
      try {
        Thread.sleep(100L);
      } catch (InterruptedException ex) {
      }
    }

    complete = false;
    while (!complete) {
      complete = true;
      for (int i = 0; i < numThreads; i++) {
        if (threads[i].isAlive()) {
          complete = false;
        }
      }
      if (!complete) {
        try {
          Thread.sleep(1000L);
        } catch (InterruptedException ex) {
        }
      }
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    System.out.println(numArgs);
    String output = "";
    int numThreads = 8;
    String commonSubFolderPattern = "";
    String commonFilename = ".txt";
    String affyResultsDir = "";
    String usage = "\n" + "affy.MergeChp requires 0-1 arguments\n"
                   + "   (1) Affy results directory (i.e. affyResultsDir= (default))\n"
                   + "   (2) number of threads to use (i.e. threads=" + numThreads + " (default))"
                   + "   (3) output directory (i.e. output=" + output + " (default))" + "\n";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        return;
      } else if (arg.startsWith("threads=")) {
        numThreads = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("out=")) {
        output = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("affyResultsDir=")) {
        affyResultsDir = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      return;
    }

    try {
      combineChpFiles(affyResultsDir, numThreads, commonSubFolderPattern, commonFilename, output,
                      new Logger());
    }

    catch (Exception e) {
      e.printStackTrace();
    }
    System.out.println("Done Merging Files");
  }
}
