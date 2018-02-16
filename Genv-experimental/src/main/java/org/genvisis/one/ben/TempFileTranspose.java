package org.genvisis.one.ben;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.nio.file.FileSystems;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import com.google.common.collect.ImmutableMap;

public class TempFileTranspose {

  static class ListFileCheckoutSystem {

    static final int MAX_MINS_SLEPT = 20;

    public static void initFile(String file, String[] values) {
      if (!Files.exists(file)) {
        Files.writeArray(values, file);
      }
    }

    public static String[] checkout(String listFile, int pull, String label,
                                    Logger log) throws IOException {
      Path fP = FileSystems.getDefault().getPath(listFile);
      String nF = ext.rootOf(listFile, false) + "." + label + "." + Thread.currentThread().getName()
                  + ".temp";
      Path nFP = FileSystems.getDefault().getPath(nF);
      IOException e = null;
      double minsSlept = 0;
      do {
        try {
          java.nio.file.Files.move(fP, nFP, StandardCopyOption.ATOMIC_MOVE);

          String[] list = HashVec.loadFileToStringArray(nF, false, null, false);
          if (list.length == 0) {
            return null;
          }
          String[] values = new String[pull];
          String[] remain = new String[list.length - pull];
          System.arraycopy(list, 0, values, 0, pull);
          System.arraycopy(list, pull, remain, 0, remain.length);
          Files.writeArray(remain, nF);

          java.nio.file.Files.move(nFP, fP, StandardCopyOption.ATOMIC_MOVE);

          Files.writeArray(values, nF);
          return values;
        } catch (NoSuchFileException e1) {
          e = e1;
          try {
            // 30 seconds
            Thread.sleep(1000 * 30);
          } catch (InterruptedException e2) {}
          minsSlept += .5;
        }
      } while (e != null && minsSlept < MAX_MINS_SLEPT);

      log.reportTimeWarning("Waited for " + minsSlept + " minutes, trying to checkout " + pull
                            + " values for " + label + ", + couldn't find " + listFile);

      return null;
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
                                                      TransposeData.MARKERDATA_NUMMARKERS_START);
          if (Files.getSize(getTempFile(file)) == size) {
            add = false;
          }
        } catch (IOException e) {
          // just recreate
        }
      }
      if (add) {
        toDo.add(temp);
      }
    }
    ListFileCheckoutSystem.initFile(tempDir + "temp.list", files);
  }

  public void runFirst() throws IOException {
    nullStatus = getNullStatus();
    final long f = MarkerSet.fingerprintForMarkers(proj);

    ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime()
                                                                   .availableProcessors());
    byte numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);

    final String listFile = tempDir + "temp.list";
    for (int i = 0, c1 = Runtime.getRuntime().availableProcessors(); i < c1; i++) {
      Runnable run = () -> {
        byte[] parameter;
        byte[][] readBuffer;
        try {
          String file = null;
          while ((file = ListFileCheckoutSystem.checkout(listFile, 1, label,
                                                         new Logger())[0]) != null) {
            String out = getTempFile(file);
            parameter = readParameter(file);
            if (nullStatus != parameter[TransposeData.MARKERDATA_NULLSTATUS_START]) {
              throw new IllegalStateException("Error - null status was inconsistent between mdRAF files.  Found "
                                              + nullStatus + " and "
                                              + parameter[TransposeData.MARKERDATA_NULLSTATUS_START]);
            }
            // read entire mdRAF file into memory:
            readBuffer = MarkerDataLoader.loadFromMarkerDataRafWithoutDecompressRange(file,
                                                                                      parameter,
                                                                                      null, 0, -1,
                                                                                      f,
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
  }

  static class Outliers {

    enum TYPE {
      X, Y, LRR;
    }

    public Map<String, Map<TYPE, Float>> getMarkerOutliers(String markerName) {
      return mkrMap.containsKey(markerName) ? mkrMap.get(markerName) : new HashMap<>();
    }

    public Map<String, Map<TYPE, Float>> getSampleOutliers(String sampleName) {
      return smpMap.containsKey(sampleName) ? smpMap.get(sampleName) : new HashMap<>();
    }

    public Hashtable<String, Float> getSampleOutliersForFile(Project proj, String sampleName) {
      Hashtable<String, Float> table = new Hashtable<>();
      Map<String, Map<TYPE, Float>> outs = getSampleOutliers(sampleName);
      for (Entry<String, Map<TYPE, Float>> e : outs.entrySet()) {
        for (Entry<TYPE, Float> e1 : e.getValue().entrySet()) {
          table.put(proj.getMarkerIndices().get(e.getKey()) + "\t"
                    + e1.getKey().name().toLowerCase(), e1.getValue());
        }
      }
      return table;
    }

    public static Outliers construct(Project proj, Hashtable<String, Float> set) {
      Outliers outs = new Outliers();
      String[] allMarkers = proj.getMarkerNames();
      for (Entry<String, Float> out : set.entrySet()) {
        String[] k = out.getKey().split("\t");
        int mInd = Integer.parseInt(k[0]);
        String sName = k[1];
        TYPE typ = TYPE.valueOf(k[2].toUpperCase());
        outs.add(allMarkers[mInd], sName, typ, out.getValue());
      }
      return outs;
    }

    private Outliers() {
      mkrMap = new HashMap<>();
      smpMap = new HashMap<>();
    }

    Map<String, Map<String, Map<TYPE, Float>>> mkrMap;
    Map<String, Map<String, Map<TYPE, Float>>> smpMap;

    private void add(String mkr, String samp, TYPE typ, Float value) {
      if (!mkrMap.containsKey(mkr)) {
        mkrMap.put(mkr, new HashMap<>());
      }
      Map<String, Map<TYPE, Float>> outMap = mkrMap.get(mkr);
      if (!outMap.containsKey(samp)) {
        outMap.put(samp, new HashMap<>());
      }
      outMap.get(samp).put(typ, value);
      if (!smpMap.containsKey(samp)) {
        smpMap.put(samp, new HashMap<>());
      }
      outMap = smpMap.get(samp);
      if (!outMap.containsKey(mkr)) {
        outMap.put(mkr, new HashMap<>());
      }
      outMap.get(mkr).put(typ, value);
    }

  }

  public void setupSampleListFile() {
    ListFileCheckoutSystem.initFile(tempDir + "samp.list", proj.getSamples());
  }

  public void runSecond() throws IOException {
    nullStatus = getNullStatus();
    Outliers outliers = Outliers.construct(proj, MarkerDataLoader.loadOutliers(proj));
    String[] files = discover(); // deterministic / always the same order
    Map<String, Integer> markerCountMap = new HashMap<>();

    HashMap<String, RandomAccessFile> readerMap = new HashMap<>();
    for (String file : files) {
      String out = tempDir + ext.removeDirectoryInfo(file) + ".tpd";
      readerMap.put(out, new RandomAccessFile(out, "r"));
      markerCountMap.put(out, getMarkerCount(file));
    }

    int numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
    int numBytesPerSample = numBytesPerSampleMarker * proj.getMarkerNames().length;
    byte[] mkrCntBytes = Compression.intToBytes(proj.getMarkerNames().length);
    long fingerPrint = MarkerSet.fingerprintForSamples(proj);
    final ImmutableMap<String, Integer> sampleIndices = proj.getSampleIndices();

    int threads = Runtime.getRuntime().availableProcessors();
    ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime()
                                                                   .availableProcessors());
    for (int t = 0; t < threads; t++) {
      Runnable run = () -> {
        Hashtable<String, Float> outs;
        byte[] buffer;
        String samp;
        RandomAccessFile sampFile;
        try {
          while ((samp = ListFileCheckoutSystem.checkout(tempDir + "samp.temp", 1, label,
                                                         new Logger())[0]) != null) {
            int sInd = sampleIndices.get(samp);
            sampFile = new RandomAccessFile(proj.SAMPLE_DIRECTORY.getValue() + samp
                                            + Sample.SAMPLE_FILE_EXTENSION, "rw");

            sampFile.write(mkrCntBytes);
            sampFile.write(nullStatus);

            outs = outliers.getSampleOutliersForFile(proj, samp);
            byte[] outBytes = null;
            if (outs.size() > 0) {
              outBytes = Compression.objToBytes(outs);

              sampFile.write(Compression.intToBytes(outBytes.length));

            } else {
              sampFile.write(Compression.intToBytes(0));
            }
            sampFile.write(Compression.longToBytes(fingerPrint));

            for (String file : files) {
              String out = tempDir + ext.removeDirectoryInfo(file) + ".tpd";
              buffer = new byte[numBytesPerSampleMarker * markerCountMap.get(out)];
              synchronized (out) {
                RandomAccessFile tpRAF = readerMap.get(out);
                long seekL = sInd * (long) numBytesPerSample;
                if (tpRAF.getFilePointer() != seekL) {
                  tpRAF.seek(seekL);
                }
                tpRAF.read(buffer);
              }
              sampFile.write(buffer);
              buffer = null;
            }

            if (outBytes != null) {
              sampFile.write(outBytes);
            }
            sampFile.close();
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
  }

}
