package org.pankratzlab.common.qsub;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;

public class QsubWriter {

  public static void writeQsubHeader(
      PrintWriter writer,
      String filename,
      int totalMemoryRequestedInMb,
      double walltimeRequestedInHours,
      int numProcs,
      String nodeToUse) {
    writeQsubHeader(
        writer,
        filename,
        totalMemoryRequestedInMb,
        walltimeRequestedInHours,
        numProcs,
        nodeToUse,
        true,
        false);
  }

  public static void writeQsubHeaderGb(
      PrintWriter writer,
      String filename,
      int totalMemoryRequestedInGb,
      double walltimeRequestedInHours,
      int numProcs,
      String nodeToUse) {
    writeQsubHeader(
        writer,
        filename,
        totalMemoryRequestedInGb,
        walltimeRequestedInHours,
        numProcs,
        nodeToUse,
        true,
        true);
  }

  private static void writeQsubHeader(
      PrintWriter writer,
      String filename,
      int totalMemoryRequested,
      double walltimeRequestedInHours,
      int numProcs,
      String nodeToUse,
      boolean defualtMods,
      boolean useGb) {
    Vector<String> params;
    int hours, minutes;

    writer.println("#!/bin/bash");
    writer.println("#$ -cwd");
    writer.println("#$ -S /bin/bash");
    writer.println("#PBS -e $PBS_JOBNAME.$PBS_JOBID.e");
    writer.println("#PBS -o $PBS_JOBNAME.$PBS_JOBID.o");
    params = new Vector<>();
    if (totalMemoryRequested > 0) {
      params.add("mem=" + totalMemoryRequested + (useGb ? "gb" : "mb"));
    }
    if (walltimeRequestedInHours > 0) {
      hours = (int) Math.floor(walltimeRequestedInHours);
      minutes = (int) Math.ceil((walltimeRequestedInHours - hours) * 60.0);
      params.add("walltime=" + ext.formNum(hours, 2) + ":" + ext.formNum(minutes, 2) + ":00");
    }
    params.add("nodes=1:ppn=" + (numProcs <= 0 ? "1" : numProcs));
    if (params.size() > 0) {
      writer.println("#PBS -m ae"); // send mail when aborts or ends (add b, as in #PBS -m abe, for
      // begins as well)
      writer.println("#PBS -l " + ArrayUtils.toStr(ArrayUtils.toStringArray(params), ","));
    }

    if (nodeToUse != null) {
      writer.println("#$ -q *@" + nodeToUse);
    }
    writer.println();

    if (defualtMods) {
      writer.println(ArrayUtils.toStr(PSF.Load.getAllModules(), "\n"));
      writer.println(PSF.Cmd.getProfilePLWithOutput(ext.rootOf(filename, false) + ".profile"));
    }

    writer.println();
    writer.println("echo \"start " + ext.rootOf(filename) + " at: \" `date`");
    writer.println("/bin/hostname");
  }

  public static void qsub(
      String filename,
      String command,
      int totalMemoryRequestedInMb,
      double walltimeRequestedInHours,
      int numProcs) {
    PrintWriter writer;
    String[] lines;

    lines = command.split("\\n");
    writer = Files.getAppropriateWriter(filename);
    if (writer == null) {
      return;
    }

    writeQsubHeader(
        writer, filename, totalMemoryRequestedInMb, walltimeRequestedInHours, numProcs, null);

    boolean rewriteJavaCmd = (totalMemoryRequestedInMb / 1024) > 1; // default Java heap size is
    // min(1/4 mem avail, 1GB)

    for (String line : lines) {
      if (line.startsWith("java ")) {
        if (!line.contains("-Xmx") && rewriteJavaCmd) {
          int memG = (totalMemoryRequestedInMb / 1024);
          line = line.replace("java ", "java -Xmx" + memG + "G ");
        }
      }
      writer.println(line);
    }
    writer.println("echo \"end " + ext.rootOf(filename) + " at: \" `date`");
    writer.flush();
    writer.close();
    Files.chmod(filename, false);
  }


  public static void qsub(
      String root,
      int start,
      int stop,
      String commands,
      int memRequiredInMb,
      double walltimeRequestedInHours) {
    qsub("", root, start, stop, commands, memRequiredInMb, walltimeRequestedInHours, null);
  }

  public static void qsub(
      String dir,
      String root,
      int start,
      int stop,
      String commands,
      int memRequiredInMb,
      double walltimeRequestedInHours,
      String queue) {
    String[] lines;

    lines =
        qsub(
            dir,
            "chr#_" + root,
            start,
            stop,
            commands,
            null,
            memRequiredInMb,
            walltimeRequestedInHours,
            null);

    if (lines.length > 1) {
      Files.writeArray(lines, dir + "master." + (root == null ? "qsub" : root));
      Files.chmod(dir + "master." + (root == null ? "qsub" : root));
    }
  }

  public static String[] qsub(
      String dir,
      String filenameFormat,
      int start,
      int stop,
      String commands,
      String[] patterns,
      int totalMemoryRequestedInMb,
      double walltimeRequestedInHours,
      String nodeToUse) {
    return qsub(
        dir,
        filenameFormat,
        start,
        stop,
        commands,
        patterns,
        totalMemoryRequestedInMb,
        walltimeRequestedInHours,
        nodeToUse,
        null);
  }

  @Deprecated
  public static String[] qsub(
      String dir,
      String filenameFormat,
      int start,
      int stop,
      String commands,
      String[] patterns,
      int totalMemoryRequestedInMb,
      double walltimeRequestedInHours,
      String nodeToUse,
      String queueName) {
    PrintWriter writer;
    String filename;
    String[] lines;
    Vector<String> v;

    if (dir == null) {
      dir = "";
    } else if (!dir.equals("") && !new File(dir).exists()) {
      System.err.println(
          "Error - directory '" + dir + "' does not exist, cannot create batches there");
    }

    v = new Vector<>();
    for (int i = start; i <= stop; i++) {
      // filename = ext.parseDirectoryOfFile(prefix,
      // true)+(prefix==null?"":ext.removeDirectoryInfo(prefix))+(start!=stop||(start>0&&start<25)?(prefix==null||prefix.equals("")?i:(prefix.endsWith("chr")?"":".")+i):"")+(root==null?"":"_"+root)+".qsub";
      filename =
          ext.insertNumbers(filenameFormat, i) + (filenameFormat.endsWith(".qsub") ? "" : ".qsub");
      try {
        writer = Files.openAppropriateWriter(dir + filename);
        writeQsubHeader(
            writer, filename, totalMemoryRequestedInMb, walltimeRequestedInHours, 1, nodeToUse);
        if (patterns == null) {
          writer.println(ext.insertNumbers(commands, i));
        } else {
          writer.println("/bin/hostname");
          writer.println(">plug");
          writer.println("rep=0");
          writer.println("total_reps=0");
          writer.println("while [ -e \"plug\" ]; do ");
          writer.println(
              "    rep=$(java -jar "
                  + Files.ROOT_DIRECTORY
                  + PSF.Java.GENVISIS
                  + " common.Files -nextRep patterns="
                  + ArrayUtils.toStr(patterns, ",")
                  + " lastRep=$rep wait=1000)");
          writer.println("    echo \"Beginning replicate $rep\"");
          lines = commands.split("\n");
          for (String line : lines) {
            writer.println("    " + ext.replaceAllWith(line, "#", "$rep"));
          }
          writer.println("    rm $rep.taken");
          writer.println("    total_reps=$(expr $total_reps + 1)");
          writer.println("done");
          writer.println("echo \"Performed $total_reps replicate(s) within this thread\"");
        }
        writer.println("echo \"end " + ext.rootOf(filename) + " at: \" `date`");
        writer.close();
        Files.chmod(dir + filename, false);
        v.add("qsub " + (queueName == null ? "" : "-q " + queueName + " ") + filename);
      } catch (IOException ioe) {
        throw new RuntimeException("Problem creating " + dir + filename);
      }
    }
    return ArrayUtils.toStringArray(v);
  }

  public static void makeQsub(
      String filename,
      boolean multiple,
      int start,
      int stop,
      boolean separate,
      String[] patterns,
      boolean changeToCurrentWorkingDirectoryFirst) {
    String[] lines, qsubs;
    List<String> v;
    int numThreads;

    lines = HashVec.loadFileToStringArray(filename, false, null, false);

    if (multiple) { // could be more elegant and incorporated with those below if desired
      System.out.println("Creating a ScriptExecutor");
      numThreads = Math.min(24, Files.countLines(filename, 0));
      qsub(
          "scriptExecutorFor_" + ext.removeDirectoryInfo(filename),
          "cd "
              + ext.parseDirectoryOfFile(filename)
              + "\n"
              + Files.getRunString()
              + " one.ScriptExecutor file="
              + ext.removeDirectoryInfo(filename)
              + " threads="
              + numThreads,
          63000,
          12,
          numThreads);
    } else if (separate) {
      v = new ArrayList<>();
      for (int i = 0; i < lines.length; i++) {
        qsubs =
            qsub(
                "",
                ext.rootOf(filename) + (i + 1) + ".#",
                start,
                stop,
                (changeToCurrentWorkingDirectoryFirst ? "cd " + ext.pwd() + "\n" : "") + lines[i],
                patterns,
                Files.PBS_MEM,
                Files.PBS_PROC,
                null);
        v.add(qsubs[0]);
      }
      Files.writeArray(ArrayUtils.toStringArray(v), "master." + ext.rootOf(filename));
      Files.chmod("master." + ext.rootOf(filename));
    } else {
      if (changeToCurrentWorkingDirectoryFirst) {
        lines = ArrayUtils.addStrToArray("cd " + ext.pwd(), lines);
      }
      if (start == stop) {
        qsub(
            ext.rootOf(filename) + ".qsub",
            ArrayUtils.toStr(lines, "\n"),
            Files.PBS_MEM,
            Files.PBS_PROC,
            1);
      } else {
        qsubs =
            qsub(
                "",
                ext.rootOf(filename) + "#",
                start,
                stop,
                ArrayUtils.toStr(lines, "\n"),
                patterns,
                Files.PBS_MEM,
                Files.PBS_PROC,
                null);
        if (qsubs.length > 1) {
          Files.writeArray(qsubs, "master." + ext.rootOf(filename));
          Files.chmod("master." + ext.rootOf(filename));
        }
      }
    }
  }
}
