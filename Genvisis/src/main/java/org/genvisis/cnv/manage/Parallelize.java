package org.genvisis.cnv.manage;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Files;
import org.genvisis.common.SerializedFiles;

public class Parallelize implements Runnable {
  public static final String[] OPTIONS = {};
  public static final int LRR_COMP = 1;
  public static final int DEFAULT_TYPE = LRR_COMP;

  private final Project proj;
  private final int type;
  private final int numThreads;
  private final int threadNumber;

  public Parallelize(Project proj, int type, int numThreads, int threadNumber) {
    this.proj = proj;
    this.type = type;
    this.numThreads = numThreads;
    this.threadNumber = threadNumber;
  }

  @Override
  public void run() {
    String[] samples;
    int count, total;

    switch (type) {
      case LRR_COMP:
        samples = proj.getSamples();
        new File(proj.PROJECT_DIRECTORY.getValue() + "comps/").mkdirs();
        count = 0;
        total = samples.length / numThreads;
        for (int i = 0; i < samples.length; i++) {
          if (i % numThreads == threadNumber) {
            System.out.println((++count) + " of " + total);
            proj.getFullSampleFromRandomAccessFile(samples[i]).compLRRs(proj);
          }
        }
        break;
      default:
        break;
    }
  }

  public static void startThreads(Project proj, int type, int numThreads) {
    PrintWriter writer;
    Thread[] threads;
    boolean wait;

    threads = new Thread[numThreads];
    for (int i = 0; i < numThreads; i++) {
      threads[i] = new Thread(new Parallelize(proj, type, numThreads, i));
      threads[i].start();
      try {
        Thread.sleep(100L);
      } catch (InterruptedException ex) {
      }
    }

    wait = true;
    while (wait) {
      wait = false;
      try {
        Thread.sleep(1000L);
      } catch (InterruptedException ex) {
      }
      for (int i = 0; i < numThreads; i++) {
        if (threads[i].isAlive()) {
          wait = true;
        }
      }
    }

    switch (type) {
      case LRR_COMP:
        String[] markerNames, files;
        double[] totalRaw, totalAbs;
        int[] counts;
        float[] diffs;

        markerNames = proj.getMarkerNames();
        files = Files.list(proj.PROJECT_DIRECTORY.getValue() + "comps/", ".comp", false);

        totalRaw = new double[markerNames.length];
        totalAbs = new double[markerNames.length];
        counts = new int[markerNames.length];
        for (int i = 0; i < files.length; i++) {
          System.out.println((i + 1) + " of " + files.length);
          diffs = (float[]) SerializedFiles.readSerial(proj.PROJECT_DIRECTORY.getValue() + "comps/"
                                                       + files[i]);
          for (int j = 0; j < markerNames.length; j++) {
            if (!Float.isNaN(diffs[j])) {
              totalRaw[j] += diffs[j];
              totalAbs[j] += Math.abs(diffs[j]);
              counts[j]++;
            }
          }
        }

        try {
          writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()
                                                  + "LRR_diff_list.xln"));
          writer.println("SNP\tAvg diff\tAvg abs(diff)\tcount");
          for (int i = 0; i < markerNames.length; i++) {
            writer.println(markerNames[i] + "\t" + (totalRaw[i] / counts[i]) + "\t"
                           + (totalAbs[i] / counts[i]) + "\t" + counts[i]);
          }
          writer.close();
        } catch (Exception e) {
          System.err.println("Error writing to " + "LRR_diff_list.xln");
          e.printStackTrace();
        }
        break;
      default:
        break;
    }
  }

  public static void waitForThreads(Thread[] threads) {
    boolean alive;

    alive = true;
    while (alive) {
      alive = false;
      for (Thread thread : threads) {
        if (thread.isAlive()) {
          alive = true;
        }
      }
      try {
        Thread.sleep(1000);
      } catch (InterruptedException ie) {
      }
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    Project proj;
    String filename = null;
    int numThreads = 1;
    int type = DEFAULT_TYPE;
    boolean tabulate = true;

    String usage = "\n" + "filesys.ParseIllumina requires 0-1 arguments\n"
                   + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + "   (2) number of threads to use (i.e. threads=" + numThreads + " (default))\n"
                   + "   (3) which algorithm to run (i.e. type=" + type + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("threads=")) {
        numThreads = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("type=")) {
        type = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    proj = new Project(filename, false);

    if (!proj.SOURCE_DIRECTORY.getValue(false, true).equals("")
        && !new File(proj.SOURCE_DIRECTORY.getValue(false, true)).exists()) {
      System.err.println("Error - the project source location is invalid: "
                         + proj.SOURCE_DIRECTORY.getValue(false, true));
      return;
    }

    try {
      if (tabulate) {
        // tabulate(proj);
      } else {
        startThreads(proj, type, numThreads);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
