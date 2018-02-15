package org.genvisis.one.ben;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentLinkedQueue;
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
import org.genvisis.common.ext;

public class TempFileTranspose {

  Project proj;
  String tempDir;
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

  public void runFirst() throws IOException {
    nullStatus = getNullStatus();
    final long f = MarkerSet.fingerprint(proj.getSamples());
    String[] files = discover();
    ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime()
                                                                   .availableProcessors());
    for (String file : files) {
      String out = tempDir + ext.removeDirectoryInfo(file) + ".tpd";
      if (Files.exists(out)) continue;
      Runnable run = () -> {
        byte[] parameter;
        byte[][] readBuffer;
        byte numBytesPerSampleMarker;
        try {
          parameter = readParameter(file);
          if (nullStatus != parameter[TransposeData.MARKERDATA_NULLSTATUS_START]) {
            throw new IllegalStateException("Error - null status was inconsistent between mdRAF files.  Found "
                                            + nullStatus + " and "
                                            + parameter[TransposeData.MARKERDATA_NULLSTATUS_START]);
          }
          numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
          // read entire mdRAF file into memory:
          readBuffer = MarkerDataLoader.loadFromMarkerDataRafWithoutDecompressRange(file, parameter,
                                                                                    null, 0, -1, f,
                                                                                    proj.getLog());
          OutputStream os = new FileOutputStream(out);

          for (int s = 0, c = proj.getSamples().length; s < c; s++) {
            for (int m = 0; m < readBuffer.length; m++) { // should be all markers
              os.write(readBuffer[m], s * numBytesPerSampleMarker, numBytesPerSampleMarker);
            }
          }
          os.flush();
          os.close();

          readBuffer = null;
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
    long fingerPrint = MarkerSet.fingerprint(proj.getMarkerNames());

    ConcurrentLinkedQueue<String> sampleQueue = new ConcurrentLinkedQueue<>();
    Map<String, Integer> sampleIndexMap = proj.getSampleIndices();
    for (String s : proj.getSamples()) {
      sampleQueue.add(s);
    }
    int threads = Runtime.getRuntime().availableProcessors();
    ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime()
                                                                   .availableProcessors());
    for (int t = 0; t < threads; t++) {
      Runnable run = () -> {
        Hashtable<String, Float> outs;
        byte[] buffer;
        String samp;
        RandomAccessFile sampFile;
        while ((samp = sampleQueue.poll()) != null) {
          int sInd = sampleIndexMap.get(samp);
          try {
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
          } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
          }
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

  public TempFileTranspose(Project p, String d) {
    this.proj = p;
    this.tempDir = d;
  }

}
