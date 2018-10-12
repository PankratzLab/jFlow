package org.pankratzlab.common;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Date;
import java.util.Enumeration;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;

public class Zip {

  public static void zip(String filename) {
    zip(new String[] {filename}, filename + ".zip", null, false);
  }

  public static void zip(String[] filenames, String zipfile, Logger log, boolean verbose) {
    byte[] buf = new byte[1024];
    boolean problem;

    if (log == null) {
      log = new Logger();
    }

    problem = false;
    for (int i = 0; i < filenames.length; i++) {
      if (!Files.exists(filenames[i])) {
        log.reportError("  Failed to find file '" + filenames[i] + "' destined for " + zipfile);
        problem = true;
      }
    }
    if (problem) {
      log.reportError("  aborting creation of " + zipfile);
      return;
    }

    try {
      ZipOutputStream out = new ZipOutputStream(new FileOutputStream(zipfile));

      for (String filename : filenames) {
        FileInputStream in = new FileInputStream(filename);
        if (verbose) {
          log.report(filename);
        }

        out.putNextEntry(new ZipEntry(filename));
        int len;
        while ((len = in.read(buf)) > 0) {
          out.write(buf, 0, len);
        }

        out.closeEntry();
        in.close();
      }

      out.close();
    } catch (IOException e) {
      log.reportError("Error creating zipfile '" + zipfile + "'");
      e.printStackTrace();
    }

  }

  public static void zipDirectory(String dir, String outfile) {
    try {
      ZipOutputStream zos = new ZipOutputStream(new FileOutputStream(outfile));
      zipDir(dir, zos);
      zos.close();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private static void zipDir(String dir2zip, ZipOutputStream zos) {
    try {
      File zipDir = new File(dir2zip);
      String[] dirList = zipDir.list();
      byte[] readBuffer = new byte[2156];
      int bytesIn = 0;

      for (String element : dirList) {
        File f = new File(zipDir, element);
        if (f.isDirectory()) {
          String filePath = f.getPath();
          zipDir(filePath, zos);
        } else {
          FileInputStream fis = new FileInputStream(f);
          String name = f.getPath();
          if (name.startsWith(".")) {
            name = name.substring(2);
          }
          ZipEntry anEntry = new ZipEntry(name);
          zos.putNextEntry(anEntry);
          while ((bytesIn = fis.read(readBuffer)) != -1) {
            zos.write(readBuffer, 0, bytesIn);
          }
          fis.close();
        }
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void unzipFile(String zipfile, String destination) {
    if (zipfile.endsWith(".gz")) {
      gunzip(zipfile, destination);
    } else {
      unzip(zipfile, destination);
    }
  }

  private static void unzip(String zipfile, String destination) {
    try {
      ZipInputStream in = new ZipInputStream(new FileInputStream(zipfile));

      if (!destination.endsWith("/") && !destination.endsWith("\\")) {
        destination += "/";
      }

      String[] files = list(zipfile);
      for (String file : files) {
        ZipEntry entry = in.getNextEntry();
        String filename = destination + entry.getName();

        if (!filename.endsWith("/") && !filename.endsWith("\\")) {
          int index = Math.max(filename.lastIndexOf("/"), filename.lastIndexOf("\\"));

          if (index >= 0) {
            new File(filename.substring(0, index)).mkdirs();
          }

          OutputStream out = new FileOutputStream(destination + entry.getName());

          byte[] buf = new byte[1024];
          int len;
          while ((len = in.read(buf)) > 0) {
            out.write(buf, 0, len);
          }

          out.close();
        }

      }

      in.close();
    } catch (IOException e) {
      System.err.println("Error unzipping '" + zipfile + "'");
      e.printStackTrace();
    }

  }

  public static void gzip(String filename) {
    gzip(filename, filename + ".gz", true);
  }

  public static void gzip(String filename, String outputFilename, boolean retainTimestamp) {
    GZIPOutputStream out;
    FileInputStream in;
    byte[] buf = new byte[1024];
    int len;

    try {
      out = new GZIPOutputStream(new FileOutputStream(outputFilename));
      in = new FileInputStream(filename);

      while ((len = in.read(buf)) > 0) {
        out.write(buf, 0, len);
      }
      in.close();
      out.close();
    } catch (IOException e) {
      System.err.println("Error creating zipfile '" + outputFilename + "'");
      e.printStackTrace();
    }

    if (retainTimestamp) {
      new File(outputFilename).setLastModified(new File(filename).lastModified());
    }
  }

  private static void gunzip(String zipfile, String destinationDirectory) {
    try {
      GZIPInputStream in = new GZIPInputStream(new FileInputStream(zipfile));

      if (!destinationDirectory.endsWith("/") && !destinationDirectory.endsWith("\\")) {
        destinationDirectory += "/";
      }

      OutputStream out = new FileOutputStream(destinationDirectory + ext.rootOf(zipfile, true));

      byte[] buf = new byte[1024];
      int len;
      while ((len = in.read(buf)) > 0) {
        out.write(buf, 0, len);
      }

      out.close();

      in.close();
    } catch (IOException e) {
      System.err.println("Error unzipping '" + zipfile + "'");
      e.printStackTrace();
    }

  }

  private static class GzipWorker implements Callable<Boolean> {

    private final String fileIn;
    private final String fileOut;
    private final Logger log;
    private final boolean overwrite;

    public GzipWorker(String fileIn, String fileOut, Logger log, boolean overwrite) {
      super();
      this.fileIn = fileIn;
      this.fileOut = fileOut;
      this.log = log;
      this.overwrite = overwrite;
    }

    @Override
    public Boolean call() throws Exception {
      if (Files.exists(fileOut) && !overwrite) {
        log.reportTimeWarning("The file " + fileOut + " already exists, skipping compression");
      } else {
        gzip(fileIn, fileOut, true);
      }
      return Files.exists(fileOut);
    }
  }

  private static class GzipProducer extends AbstractProducer<Boolean> {

    private final String[] filesToGzip;
    private final String outputDir;
    private final Logger log;
    private int index;
    private final boolean overwrite;

    public GzipProducer(String[] filesToGzip, String outputDir, Logger log, boolean overwrite) {
      super();
      this.filesToGzip = filesToGzip;
      this.outputDir = outputDir;
      this.log = log;
      index = 0;
      this.overwrite = overwrite;
    }

    @Override
    public boolean hasNext() {
      return index < filesToGzip.length;
    }

    @Override
    public Callable<Boolean> next() {
      log.report("compressing file #" + (index + 1) + " of " + filesToGzip.length + " ("
                 + filesToGzip[index] + ")");
      GzipWorker worker = new GzipWorker(filesToGzip[index],
                                         outputDir + ext.removeDirectoryInfo(filesToGzip[index])
                                                             + ".gz",
                                         log, overwrite);
      index++;
      return worker;
    }
  }

  public static void gzipAllTextFilesInDirectory(String dirin, String dirout, int numThreads,
                                                 Logger log) {
    String[] files;

    dirin = ext.verifyDirFormat(dirin);
    files = Files.list(dirin, null);

    if (files.length == 0) {
      log.reportError("Could not find any files in: " + dirin);
      return;
    }

    gzipMany(files, dirin, dirout, numThreads, log, false);
  }

  public static void gzipMany(String[] files, String dirin, String dirout, int numThreads,
                              Logger log, boolean overwrite) {
    if (dirin == null) {
      log.reportError("Error - dir cannot be null");
      return;
    }

    if (!Files.exists(dirin)) {
      log.reportError("Error - invalid directory: " + dirin);
      return;
    }

    dirin = ext.verifyDirFormat(dirin);
    if (dirout == null) {
      dirout = dirin + "compressed/";
    }

    new File(dirout).mkdirs();
    if (numThreads > 1) {
      GzipProducer producer = new GzipProducer(Files.toFullPaths(files, dirin), dirout, log,
                                               overwrite);
      WorkerTrain<Boolean> train = new WorkerTrain<Boolean>(producer, numThreads, numThreads, log);
      int index = 0;
      while (train.hasNext()) {
        if (!train.next()) {
          log.reportError("Could not compress file " + files[index]);
        }
        index++;
      }
    } else {
      for (int i = 0; i < files.length; i++) {
        log.report("compressing file #" + (i + 1) + " of " + files.length + " (" + files[i] + ")");
        gzip(dirin + files[i], dirout + files[i] + ".gz", true);
      }
    }

  }

  public static void unzipParticularFile(String target, String zipfile, String destination) {
    unzipParticularFiles(new String[] {target}, zipfile, destination);
  }

  public static void unzipParticularFiles(String[] targets, String zipfile, String destination) {
    try {
      ZipInputStream in = new ZipInputStream(new FileInputStream(zipfile));

      if (!destination.endsWith("/") && !destination.endsWith("\\")) {
        destination += "/";
      }

      String[] files = list(zipfile);
      int numErrors = 0;
      for (String target : targets) {
        // targets[i] = ext.replaceAllWith(targets[i], "\\", "/");
        if (ext.indexOfStr(target, files) == -1) {
          System.err.println("Error - '" + target + "' was not found in '" + zipfile + "'");
        }
      }
      if (numErrors > 0) {
        System.exit(1);
      }
      for (String file : files) {
        ZipEntry entry = in.getNextEntry();
        String filename = destination + entry.getName();

        if (ext.indexOfStr(entry.getName(), targets) != -1 && !filename.endsWith("/")
            && !filename.endsWith("\\")) {
          int index = Math.max(filename.lastIndexOf("/"), filename.lastIndexOf("\\"));

          if (index >= 0) {
            new File(filename.substring(0, index)).mkdirs();
          }

          OutputStream out = new FileOutputStream(destination + entry.getName());

          byte[] buf = new byte[1024];
          int len;
          while ((len = in.read(buf)) > 0) {
            out.write(buf, 0, len);
          }

          out.close();
        }

      }

      in.close();
    } catch (IOException e) {
      System.err.println("Error unzipping '" + zipfile + "'");
      e.printStackTrace();
    }

  }

  @SuppressWarnings("unchecked")
  public static String[] list(String filename) {
    Vector<String> v = new Vector<>();

    try {
      ZipFile zf = new ZipFile(filename);
      for (Enumeration<ZipEntry> entries = (Enumeration<ZipEntry>) zf.entries(); entries.hasMoreElements();) {
        v.add((entries.nextElement()).getName());
      }
      zf.close();
    } catch (IOException e) {
      System.err.println("Error listing contents of zipfile '" + filename + "'");
      e.printStackTrace();
    }

    return ArrayUtils.toStringArray(v);
  }

  public static void fromParameters(String filename, Logger log) {
    List<String> params;
    String[] files;
    int minutesToSleep;
    long time;
    String outputDirectory = "";
    boolean retainOriginalTimestamp = true;

    // get all files in the directory, excluding the crf itself and its corresponding log
    files = Files.list("./", ":" + ext.rootOf(filename), ":.crf", false);
    files = ArrayUtils.addStrToArray("# Output directory (i.e., sometimes more efficient if on a different drive); default is the same directory as the original file",
                                     files, 0);
    files = ArrayUtils.addStrToArray("outputDirectory=", files, 1);
    files = ArrayUtils.addStrToArray("# retain the timestamp of the original file", files, 2);
    files = ArrayUtils.addStrToArray("retainTimestamp=" + retainOriginalTimestamp, files, 3);
    files = ArrayUtils.addStrToArray("# Number of minutes to wait before beginning a step (can be used more than once)",
                                     files, 4);
    files = ArrayUtils.addStrToArray("delay_in_minutes=0", files, 5);

    params = Files.parseControlFile(filename, "gzip", files, log);

    time = new Date().getTime();
    if (params != null) {
      try {
        for (int i = 0; i < params.size(); i++) {
          filename = params.get(i);
          if (filename.startsWith("delay_in_minutes=")) {
            minutesToSleep = ext.parseIntArg(params.get(i));
            log.report("Going to sleep for "
                       + ext.getTimeElapsed(new Date().getTime() - minutesToSleep * 60 * 1000));
            try {
              Thread.sleep(minutesToSleep * 60 * 1000);
            } catch (InterruptedException ie) {}
          } else if (filename.startsWith("outputDirectory=")) {
            outputDirectory = ext.parseStringArg(params.get(i), "");
            log.report("Setting output directory to " + outputDirectory);
            if (!Files.isDirectory(outputDirectory)) {
              log.report("Output directory does not exist; attempting to create...");
              if (new File(outputDirectory).mkdirs()) {
                log.report("successful");
              } else {
                log.report("failed");
              }
            }
          } else if (filename.startsWith("retainTimestamp=")) {
            retainOriginalTimestamp = ext.parseBooleanArg(params.get(i));
            log.report("The timestamp from the original file will "
                       + (retainOriginalTimestamp ? "" : "NOT ") + "be retained");
          } else if (Files.exists(filename)) {
            log.report("\n" + ext.getTime() + "\tCompressing " + filename + " to " + outputDirectory
                       + filename + ".gz");
            gzip(filename, outputDirectory + filename + ".gz", retainOriginalTimestamp);
          } else {
            log.reportError("Error - could not find file " + filename);
          }

        }
        System.out.println(ext.getTime() + "\tfinished everything in " + ext.getTimeElapsed(time));
      } catch (Exception e) {
        System.err.println("Error at some point");
        e.printStackTrace();
      }
      ext.waitForResponse();
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    String dirin = "./";
    String dirout = null;
    String logfile = null;
    int numthreads = 1;
    Logger log;

    String usage = "\n" + "common.Zip requires 0-1 arguments\n"
                   + "   (1) name of file to be gzipped (i.e. file=in.csv (not the default))\n"
                   + " OR:\n" + "   (1) directory containing files to be gzipped (i.e. dirin="
                   + dirin + " (default))\n"
                   + "   (2) directory where the compressed files should be written to (i.e. dirout=[dirin]/compressed/ (default))\n"
                   + PSF.Ext.getNumThreadsCommand(3, numthreads) + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("dirin=")) {
        dirin = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("dirout=")) {
        dirout = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numthreads = ext.parseIntArg(arg);
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

      // dirin = "I:/GEDI/02_After_reclustering/00src/";
      // dirout = "C:/GEDI_compressed/";

      // dirin = "D:/data/GEDI/penn_data/";

      // dirin = "D:/data/genomestudio/GEDI_Exome/FinalReportsChargeTop/";
      // dirout = "C:/GEDI_exomeChargeTop_compressed/";

      // dirin = "C:/LoganReports/";
      // dirout = "I:/Logan/";

      // dirin = "C:/GEDI_Exome2/00src/";
      // dirout = "I:/GEDI/exome2_charge/";

      log = new Logger(logfile);
      if (filename != null) {
        gzip(filename);
      } else if (dirin != null) {
        gzipAllTextFilesInDirectory(dirin, dirout, numthreads, log);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
