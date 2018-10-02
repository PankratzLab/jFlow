package org.genvisis.cnv.manage;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.nio.file.FileSystems;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.core.CLI;
import com.google.common.collect.ImmutableMap;

public class TempFileTranspose {

  static class ListFileCheckoutSystem {

    static final int MAX_MINS_SLEPT = 10;

    public static void initFile(String file, String[] values) {
      if (!Files.exists(file)) {
        Files.writeArray(values, file);
      }
    }

    public static String[] checkout(String listFile, int pull, String label,
                                    Logger log) throws IOException {
      Path filePath = FileSystems.getDefault().getPath(listFile);
      String newFile = ext.rootOf(listFile, false) + "." + (label != null ? label + "." : "")
                       + Thread.currentThread().getId() + ".temp";
      Path newFilePath = FileSystems.getDefault().getPath(newFile);
      IOException e = null;
      double minsSlept = 0;
      do {
        try {
          // throws the caught exception if missing
          String[] list;
          synchronized (StandardCopyOption.ATOMIC_MOVE) {
            java.nio.file.Files.move(filePath, newFilePath, StandardCopyOption.ATOMIC_MOVE);
            list = HashVec.loadFileToStringArray(newFile, false, null, false);
          }
          if (list.length == 0) {
            return new String[0];
          }
          int p = Math.min(pull, list.length);
          int r = Math.max(0, list.length - p);
          String[] values = new String[p];
          String[] remain = new String[r];
          System.arraycopy(list, 0, values, 0, p);
          if (r > 0) {
            System.arraycopy(list, p, remain, 0, remain.length);
            Files.writeArray(remain, newFile);
            synchronized (StandardCopyOption.ATOMIC_MOVE) {
              java.nio.file.Files.move(newFilePath, filePath, StandardCopyOption.ATOMIC_MOVE);
            }
          }

          Files.writeArray(values, newFile);
          return values;
        } catch (NoSuchFileException e1) {
          e = e1;
          try {
            // 3 seconds
            Thread.sleep(1000 * 3);
          } catch (InterruptedException e2) {}
          minsSlept += .05;
        }
      } while (e != null && minsSlept < MAX_MINS_SLEPT);

      log.reportTimeWarning("Waited for " + minsSlept + " minutes, trying to checkout " + pull
                            + " values for " + label + ", + couldn't find " + listFile);

      return new String[0];
    }

    public static void checkin(String listFile, String[] values, int start, String label,
                               Logger log) throws IOException {
      Path fP = FileSystems.getDefault().getPath(listFile);
      String nF = ext.rootOf(listFile, false) + "." + label + "." + Thread.currentThread().getName()
                  + ".temp";
      Path nFP = FileSystems.getDefault().getPath(nF);
      IOException e = null;
      double minsSlept = 0;
      int maxMinsSlept = 20;
      do {
        try {
          java.nio.file.Files.move(fP, nFP, StandardCopyOption.ATOMIC_MOVE);

          int remain = values.length - start;
          String[] list = HashVec.loadFileToStringArray(nF, false, null, false);
          String[] total = new String[list.length + remain];
          System.arraycopy(values, start, total, 0, remain);
          System.arraycopy(list, 0, total, remain, list.length);
          Files.writeArray(total, nF);

          java.nio.file.Files.move(nFP, fP, StandardCopyOption.ATOMIC_MOVE);
          return;
        } catch (NoSuchFileException e1) {
          e = e1;
          try {
            // 30 seconds
            Thread.sleep(1000 * 30);
          } catch (InterruptedException e2) {}
          minsSlept += .5;
        }
      } while (e != null && minsSlept < maxMinsSlept);

      log.reportTimeWarning("Waited for " + minsSlept + " minutes, trying to readd "
                            + (values.length - start) + " values for " + label + ", couldn't find "
                            + listFile);

    }

  }

  Project proj;
  String tempDir;
  String label;
  private byte nullStatus = -1;
  private String qsubPath = null;

  private String[] discover() {
    LinkedHashSet<String> files = new LinkedHashSet<>();
    for (String s : proj.getMarkerNames()) {
      files.add(proj.MARKER_DATA_DIRECTORY.getValue()
                + proj.getMarkerLookup().get(s).split("\t")[0]);
    }
    return files.toArray(new String[files.size()]);
  }

  private byte getNullStatus() {
    String[] files = discover();
    return TransposeData.getNullstatusFromRandomAccessFile(files[0], false);
  }

  private int getMarkerCount(String mdRAF) throws IOException {
    byte[] p = readParameter(mdRAF);
    return Compression.bytesToInt(p, TransposeData.MARKERDATA_NUMMARKERS_START);
  }

  public String getTempFile(String markerFile) {
    return tempDir + ext.removeDirectoryInfo(markerFile) + ".tpd";
  }

  public void setupMarkerListFile() {
    String[] files = discover();
    List<String> toDo = new ArrayList<>();
    for (String file : files) {
      boolean add = true;
      String temp = getTempFile(file);
      if (Files.exists(temp)) {
        try {
          long size = Sample.getNBytesPerSampleMarker(TransposeData.getNullstatusFromRandomAccessFile(file,
                                                                                                      false))
                      * (long) Compression.bytesToInt(readParameter(file),
                                                      TransposeData.MARKERDATA_NUMMARKERS_START)
                      * proj.getSamples().length;
          long found = Files.getSize(getTempFile(file));
          if (found == size) {
            add = false;
          } else {
            proj.getLog()
                .reportTime("Re-processing " + file + "; found " + found + ", expected " + size);
          }
        } catch (IOException e) {
          // just recreate
        }
      }
      if (add) {
        toDo.add(file);
      }
    }
    proj.getLog().reportTime("Processing " + toDo.size() + " marker files.");
    ListFileCheckoutSystem.initFile(getMarkerListFile(), toDo.toArray(new String[toDo.size()]));
  }

  public String getMarkerListFile() {
    return tempDir + "temp.list";
  }

  public void runFirst() throws IOException {
    nullStatus = getNullStatus();
    final long f = MarkerSet.fingerprintForMarkers(proj);

    int proc = Runtime.getRuntime().availableProcessors();
    ExecutorService executor = Executors.newFixedThreadPool(proc);
    byte numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
    final String listFile = getMarkerListFile();

    for (int i = 0; i < proc; i++) {
      Runnable run = () -> {
        try {
          String[] files = null;
          while ((files = ListFileCheckoutSystem.checkout(listFile, 1, label,
                                                          proj.getLog())) != null
                 && files.length > 0) {
            processOneMDRAF(f, numBytesPerSampleMarker, files[0]);
          }
        } catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      };
      executor.submit(run);
    }
    executor.shutdown();
    try {
      executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }
    proj.getLog().reportTime("Creation of temporary files for reverse-transposing is complete!");
  }

  public void runFirstResub() throws IOException {
    nullStatus = getNullStatus();
    final long f = MarkerSet.fingerprintForMarkers(proj);

    int proc = Runtime.getRuntime().availableProcessors();
    ExecutorService executor = Executors.newFixedThreadPool(proc);
    final byte numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
    final String listFile = getMarkerListFile();

    final String[] preSel = ListFileCheckoutSystem.checkout(listFile, proc, label, proj.getLog());

    if (preSel.length > 0) {
      proj.getLog().reportTime("Checked out " + preSel.length + "; expected " + proc);
      for (int i = 0; i < preSel.length; i++) {
        final int ind = i;
        Runnable run = () -> {
          try {
            processOneMDRAF(f, numBytesPerSampleMarker, preSel[ind]);
          } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
          }
        };
        executor.submit(run);
        proj.getLog().report("Submitted process for mdRAF file: " + preSel[ind]);
      }
      executor.shutdown();
      try {
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
      } catch (InterruptedException e) {
        e.printStackTrace();
      }
      if (qsubPath != null) {
        CmdLine.run("qsub " + qsubPath, ext.pwd());
      }
    }
  }

  private void processOneMDRAF(final long f, byte numBytesPerSampleMarker,
                               String file) throws IOException, FileNotFoundException {
    byte[] parameter;
    byte[][] readBuffer;
    String out = getTempFile(file);
    parameter = readParameter(file);
    if (nullStatus != parameter[TransposeData.MARKERDATA_NULLSTATUS_START]) {
      throw new IllegalStateException("Error - null status was inconsistent between mdRAF files.  Found "
                                      + nullStatus + " and "
                                      + parameter[TransposeData.MARKERDATA_NULLSTATUS_START]);
    }
    // read entire mdRAF file into memory:
    readBuffer = MarkerDataLoader.loadFromMarkerDataRafWithoutDecompressRange(file, parameter, null,
                                                                              0, -1, f,
                                                                              proj.getLog());
    OutputStream os = new FileOutputStream(out);
    for (int s = 0, c = proj.getSamples().length; s < c; s++) {
      for (int m = 0; m < readBuffer.length; m++) { // should be all markers
        os.write(readBuffer[m], s * numBytesPerSampleMarker, numBytesPerSampleMarker);
      }
    }
    os.flush();
    os.close();

    os = null;
    readBuffer = null;
  }

  public String getSampleListFile() {
    return tempDir + "samp.list";
  }

  public void setupSampleListFile() {
    new File(proj.SAMPLE_DIRECTORY.getValue()).mkdirs();
    String[] samples = Arrays.copyOf(proj.getSamples(), proj.getSamples().length);
    for (int i = 0; i < samples.length; i++) {
      samples[i] = proj.SAMPLE_DIRECTORY.getValue() + samples[i] + ".sampRAF";
    }
    ListFileCheckoutSystem.initFile(getSampleListFile(), samples);
  }

  public void runSecond() throws IOException {
    nullStatus = getNullStatus();
    OutOfRangeValues outliers = OutOfRangeValues.construct(proj);
    String[] files = discover(); // deterministic / always the same order
    Map<String, Integer> markerCountMap = new HashMap<>();

    int threads = Runtime.getRuntime().availableProcessors();
    ExecutorService executor = Executors.newFixedThreadPool(threads);

    String listFile = getSampleListFile();

    int numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
    int numBytesPerSample = numBytesPerSampleMarker * proj.getMarkerNames().length;
    byte[] mkrCntBytes = Compression.intToBytes(proj.getMarkerNames().length);
    long fingerPrint = MarkerSet.fingerprintForSamples(proj);
    final ImmutableMap<String, Integer> sampleIndices = proj.getSampleIndices();

    HashMap<String, RandomAccessFile> readerMap = new HashMap<>();
    for (String file : files) {
      String out = getTempFile(file);
      // NOTE - file system must be capable of holding files.length number of file handles open at once
      readerMap.put(out, new RandomAccessFile(out, "r"));
      markerCountMap.put(out, getMarkerCount(file));
    }
    for (int t = 0; t < threads; t++) {
      Runnable run = () -> {
        String[] samps;
        try {
          while ((samps = ListFileCheckoutSystem.checkout(listFile, 1, label,
                                                          proj.getLog())) != null
                 && samps.length > 0) {
            processOneSAMPRAF(outliers, files, markerCountMap, numBytesPerSampleMarker,
                              numBytesPerSample, mkrCntBytes, fingerPrint, sampleIndices, readerMap,
                              samps[0]);
          }
        } catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        } catch (Exception e) {
          e.printStackTrace();
        }
      };
      executor.submit(run);
    }
    executor.shutdown();
    try {
      executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    for (RandomAccessFile f : readerMap.values()) {
      f.close();
    }
    proj.getLog().reportTime("Transposing marker files to sample files is complete!");
  }

  public void runSecondResub() throws IOException {
    nullStatus = getNullStatus();
    final OutOfRangeValues outliers = OutOfRangeValues.construct(proj);
    final String[] files = discover(); // deterministic / always the same order
    final Map<String, Integer> markerCountMap = new HashMap<>();

    final int threads = Runtime.getRuntime().availableProcessors();
    ExecutorService executor = Executors.newFixedThreadPool(threads);

    String listFile = getSampleListFile();
    final String[] preSel = ListFileCheckoutSystem.checkout(listFile, threads * 1000, label,
                                                            proj.getLog());

    int numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
    int numBytesPerSample = numBytesPerSampleMarker * proj.getMarkerNames().length;
    byte[] mkrCntBytes = Compression.intToBytes(proj.getMarkerNames().length);
    long fingerPrint = MarkerSet.fingerprintForSamples(proj);
    final ImmutableMap<String, Integer> sampleIndices = proj.getSampleIndices();

    if (preSel.length > 0) {
      proj.getLog()
          .reportTime("Checked out " + preSel.length + "; expected (max) " + threads * 1000);
      HashMap<String, RandomAccessFile> readerMap = new HashMap<>();
      for (String file : files) {
        String out = getTempFile(file);
        // NOTE - file system must be capable of holding files.length number of file handles open at once
        readerMap.put(out, new RandomAccessFile(out, "r"));
        markerCountMap.put(out, getMarkerCount(file));
      }
      AtomicInteger index = new AtomicInteger(0);
      for (int t = 0; t < threads; t++) {
        Runnable run = () -> {
          String samp;
          try {
            int ind = -1;
            while ((ind = index.getAndIncrement()) < preSel.length) {
              samp = preSel[ind];
              processOneSAMPRAF(outliers, files, markerCountMap, numBytesPerSampleMarker,
                                numBytesPerSample, mkrCntBytes, fingerPrint, sampleIndices,
                                readerMap, samp);
            }
          } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
          }
        };
        executor.submit(run);
      }
      executor.shutdown();
      try {
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
      } catch (InterruptedException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      for (RandomAccessFile f : readerMap.values()) {
        f.close();
      }
      if (qsubPath != null) {
        CmdLine.run("qsub " + qsubPath, ext.pwd());
      }
    }

  }

  private void processOneSAMPRAF(OutOfRangeValues outliers, String[] files,
                                 Map<String, Integer> markerCountMap, int numBytesPerSampleMarker,
                                 int numBytesPerSample, byte[] mkrCntBytes, long fingerPrint,
                                 final ImmutableMap<String, Integer> sampleIndices,
                                 HashMap<String, RandomAccessFile> readerMap,
                                 String samp) throws FileNotFoundException, IOException {
    Hashtable<String, Float> outs;
    byte[] buffer;
    RandomAccessFile sampFile;
    int sInd = -1;
    if (sampleIndices.containsKey(samp)) {
      sInd = sampleIndices.get(samp);
    } else if (sampleIndices.containsKey(ext.rootOf(samp, true))) {
      sInd = sampleIndices.get(ext.rootOf(samp, true));
    } else {
      throw new RuntimeException("Error - sample " + samp + " not found in the sample index map!");
    }
    sampFile = new RandomAccessFile(samp, "rw");

    sampFile.seek(0);
    sampFile.write(mkrCntBytes);
    sampFile.write(nullStatus);

    if (outliers.hasSample(samp)) {
      outs = outliers.getSampleOutliersForFile(proj, samp);
    } else if (outliers.hasSample(ext.rootOf(samp, true))) {
      outs = outliers.getSampleOutliersForFile(proj, ext.rootOf(samp, true));
    } else {
      proj.getLog().reportTimeWarning("No outliers found for sample " + samp);
      outs = new Hashtable<>();
    }
    byte[] outBytes = Compression.objToBytes(outs);
    if (outs.size() > 0) {
      sampFile.write(Compression.intToBytes(outBytes.length));
    } else {
      sampFile.write(Compression.intToBytes(0));
    }
    sampFile.write(Compression.longToBytes(fingerPrint));

    long aveR = 0;
    long aveW = 0;
    for (int i = 0; i < files.length; i++) {
      long tO = System.nanoTime();
      String out = getTempFile(files[i]);
      buffer = new byte[numBytesPerSampleMarker * markerCountMap.get(out)];
      long seekL = sInd * (long) buffer.length;
      RandomAccessFile tpRAF = readerMap.get(out);
      synchronized (tpRAF) {
        if (tpRAF.getFilePointer() != seekL) {
          tpRAF.seek(seekL);
        }
        tpRAF.read(buffer);
      }
      long tE = System.nanoTime() - tO;
      aveR += ((tE - aveR) / (i + 1));
      tO = System.nanoTime();
      sampFile.write(buffer);
      tE = System.nanoTime() - tO;
      aveW += ((tE - aveW) / (i + 1));
      buffer = null;
    }

    if (outs.size() > 0 && outBytes != null) {
      sampFile.write(outBytes);
    }
    sampFile.close();

    proj.getLog()
        .reportTime("Compiled sample " + samp + "; average read "
                    + ext.formatTimeElapsed(aveR, TimeUnit.NANOSECONDS) + ", average write "
                    + ext.formatTimeElapsed(aveW, TimeUnit.NANOSECONDS));
  }

  private byte[] readParameter(String file) throws IOException {
    RandomAccessFile fil = new RandomAccessFile(file, "r");
    byte[] parameterReadBuffer = new byte[TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN];
    fil.read(parameterReadBuffer);
    fil.close();
    return parameterReadBuffer;
  }

  public TempFileTranspose(Project p, String d, String l) {
    this.proj = p;
    this.tempDir = d;
    this.label = l;
    new File(tempDir).mkdirs();
  }

  public static void main(String[] args) throws IOException {
    String tempDir;
    String jobID;

    CLI cli = new CLI(TempFileTranspose.class);

    cli.addArg("proj", "Project properties file", true);
    cli.addArg("jobID", "jobID", true);
    cli.addArg("type", "Which type of files to run: 'M' for markers, 'S' for samples.", true);

    cli.addArg("temp", "Temporary file directory", false);
    cli.addArg("qsub",
               "QSUB file; if specified, the job will process N-Threads of files (either markers or samples) and then finish and submit the specified qsub file to the queue.",
               false);
    cli.addFlag("setup", "Create marker/sample lists");

    cli.parseWithExit(args);

    Project proj = new Project(cli.get("proj"));
    tempDir = cli.has("temp") ? cli.get("temp") : proj.PROJECT_DIRECTORY.getValue() + "temp/";
    jobID = cli.get("jobID");

    TempFileTranspose tft = new TempFileTranspose(proj, tempDir, jobID);

    switch (cli.get("type")) {
      case "M":
      case "m":
        if (cli.has("setup")) {
          tft.setupMarkerListFile();
        }
        if (cli.has("qsub")) {
          tft.qsubPath = cli.get("qsub");
          tft.runFirstResub();
        } else {
          tft.runFirst();
        }
        break;
      case "S":
      case "s":
        if (cli.has("setup")) {
          tft.setupSampleListFile();
        }
        if (cli.has("qsub")) {
          tft.qsubPath = cli.get("qsub");
          tft.runSecondResub();
        } else {
          tft.runSecond();
        }
        break;
      default:
        throw new IllegalArgumentException("Invalid type: " + cli.get("type"));
    }

  }

}
