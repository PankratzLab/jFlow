package org.genvisis.cnv.manage;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerLookup;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.Array;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;

import com.google.common.primitives.Ints;

public class MDL implements Iterator<MarkerData> {
  private final Project proj;
  private final String[] markerNames;
  private final MarkerSet markerSet;
  private final int numDecompressThreads;
  private WorkerTrain<MarkerData> decompTrain;
  private final MarkerLookup markerLookup;
  private final ArrayList<FileMatch> files;
  private int numLoaded, fileIndex, numLoadedForFile;
  private final int markerBuffer;
  private final Hashtable<String, String> missing;
  private BufferReader producer;
  private int reportEvery;
  private final long startTime;
  private boolean debugMode;
  private String currentFile;

  /**
   * @param proj
   * @param markerSet if null, will be loaded
   * @param markerNames them to load
   * @param numDecompressThreads number of threads used to decompress the marker
   * @param markerBuffer number of markers to hold in the queue for processing
   */
  public MDL(Project proj, MarkerSet markerSet, String[] markerNames, int numDecompressThreads,
             int markerBuffer) {
    this.proj = proj;
    missing = new Hashtable<String, String>();
    this.markerNames = markerNames;
    this.markerSet = markerSet == null ? proj.getMarkerSet() : markerSet;
    this.numDecompressThreads = numDecompressThreads;
    markerLookup = proj.getMarkerLookup();
    files = matchFileNames();
    this.markerBuffer = markerBuffer;
    numLoaded = 0;
    fileIndex = 0;
    reportEvery = 1000;
    startTime = new Date().getTime();
    debugMode = false;
  }

  /**
   * @param match starts a new loader for the file in {@link FileMatch}
   */
  private void initTrain(FileMatch match) {
    if (decompTrain != null) {
      decompTrain.shutdown();
    }
    producer = new BufferReader(proj, markerSet, match.getFileName(),
                                Ints.toArray(match.getFileIndices()),
                                Ints.toArray(match.getProjIndices()), debugMode);
    try {
      producer.init();
    } catch (IllegalStateException e) {
      proj.getLog().reportException(e);
      e.printStackTrace();
    } catch (IOException e) {
      proj.getLog().reportIOException(match.getFileName());
      e.printStackTrace();
    }
    currentFile = match.fileName;
    decompTrain = new WorkerTrain<MarkerData>(producer, numDecompressThreads, markerBuffer,
                                              proj.getLog());
  }

  @Override
  public boolean hasNext() {
    boolean hasNext = numLoaded < markerNames.length;
    if (!hasNext) {
      shutdown();
    }
    return hasNext;
  }

  public void shutdown() {
    if (decompTrain != null) {
      decompTrain.shutdown();
    }
  }

  @Override
  public MarkerData next() throws IllegalStateException {
    MarkerData markerData = null;
    try {
      if (numLoaded == 0 && fileIndex == 0) {
        initTrain(files.get(fileIndex));
        numLoadedForFile = 0;
      } else if (numLoadedForFile == files.get(fileIndex).getMarkersToLoad().size()) {
        fileIndex++;
        initTrain(files.get(fileIndex));
        numLoadedForFile = 0;
      }
      decompTrain.hasNext();
      if (!decompTrain.hasNext()) {
        String error = "Internal index error, halting";
        proj.getLog().reportTimeError(error);
        throw new IllegalStateException(error);
      }
      numLoadedForFile++;
      numLoaded++;
      if ((startTime - new Date().getTime()) % reportEvery == 0) {
        proj.getLog()
            .reportTimeInfo("Loaded " + numLoaded + " of " + markerNames.length + " markers");
      }
      markerData = decompTrain.next();
      if (!markerData.getMarkerName().equals(markerNames[numLoaded - 1])) {
        String error = "Internal index error - mismatched marker data found, halting";
        proj.getLog().reportTimeError(error);
        throw new IllegalStateException(error);
      }
    } catch (NullPointerException npe) {
      proj.getLog()
          .reportTimeError("Could not load " + markerNames[numLoaded - 1]
                           + ", this is usually caused by a corrupt or missing outlier file");

    }
    return markerData;
  }

  @Override
  public void remove() {
    // TODO Auto-generated method stub

  }

  /**
   *
   * @param reportEvery new reporting interval, in milliseconds
   */
  public void setReportEvery(int reportEvery) {
    this.reportEvery = reportEvery;
  }

  public void setDebugMode(boolean debugMode) {
    this.debugMode = debugMode;
  }

  public String getCurrentFile() {
    return currentFile;
  }

  /**
   * Stores the markers and indices for loading
   *
   */
  private static class FileMatch {
    private final String fileName;
    private final ArrayList<String> markersToLoad;
    private final ArrayList<Integer> projIndices;
    private final ArrayList<Integer> fileIndices;

    private FileMatch(String fileName) {
      super();
      this.fileName = fileName;
      markersToLoad = new ArrayList<String>(100000);
      projIndices = new ArrayList<Integer>(100000);
      fileIndices = new ArrayList<Integer>(100000);
    }

    private String getFileName() {
      return fileName;
    }

    private ArrayList<String> getMarkersToLoad() {
      return markersToLoad;
    }

    public ArrayList<Integer> getProjIndices() {
      return projIndices;
    }

    private ArrayList<Integer> getFileIndices() {
      return fileIndices;
    }

    private void add(String markerName, int projIndex, int fileIndex) {
      markersToLoad.add(markerName);
      projIndices.add(projIndex);
      fileIndices.add(fileIndex);
    }

  }

  /**
   * @return {@link FileMatch} objects for determining order and markers to scan from each file
   */
  private ArrayList<FileMatch> matchFileNames() {
    ArrayList<FileMatch> files = new ArrayList<MDL.FileMatch>();
    int[] indicesInProject = ext.indexLargeFactors(markerNames, proj.getMarkerNames(), true,
                                                   proj.getLog(), true, false);
    String currentFile = "";
    int currentIndex = -1;
    for (int i = 0; i < markerNames.length; i++) {
      if (markerLookup.contains(markerNames[i])) {
        String[] line = markerLookup.get(markerNames[i]).split("[\\s]+");
        if (!line[0].equals(currentFile)) {
          currentFile = line[0];
          files.add(new FileMatch(proj.MARKER_DATA_DIRECTORY.getValue(false, true) + currentFile));
          currentIndex++;
        }
        files.get(currentIndex).add(markerNames[i], indicesInProject[i], Integer.parseInt(line[1]));

      } else {
        missing.put(markerNames[i], markerNames[i]);
      }
    }
    if (missing.size() > 0) {
      proj.getLog()
          .reportTimeError("Could not find the following markers in the project"
                           + Array.toStr(missing.keySet().toArray(new String[missing.size()]),
                                         "\n"));
    }
    return files;
  }

  /**
   * {@link WorkerTrain.Producer} implementation to serve up marker data
   *
   */
  private static class BufferReader extends AbstractProducer<MarkerData> {
    private final byte[] parameterReadBuffer =
                                             new byte[TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN];
    private RandomAccessFile file;
    private byte nullStatus;
    private byte bytesPerSampleMarker;
    private int numBytesPerMarker, numSamplesObserved, numBytesMarkernamesSection, numLoaded;
    private final int[] markersIndicesInFile, markerIndicesInProject;
    private int[] positions;
    private byte[] chrs;
    private String[] names;
    private long fingerprint, sampleFingerprint;
    private final String currentMarkFilename;
    private final Project proj;
    private boolean isGcNull, isXNull, isYNull, isBafNull, isLrrNull, isGenotypeNull,
        isNegativeXYAllowed;
    private final boolean debugMode;
    private Hashtable<String, Float> outlierHash;
    private final MarkerSet markerSet;

    private BufferReader(Project proj, MarkerSet markerSet, String currentMarkFilename,
                         int[] markersIndicesInFile, int[] markerIndicesInProject,
                         boolean debugMode) {
      super();
      this.proj = proj;
      this.markerSet = markerSet;
      this.currentMarkFilename = currentMarkFilename;
      this.markersIndicesInFile = markersIndicesInFile;
      this.markerIndicesInProject = markerIndicesInProject;
      numLoaded = 0;
      this.debugMode = debugMode;
    }

    @SuppressWarnings("unchecked")
    private void init() throws IOException, IllegalStateException {
      sampleFingerprint = proj.getSampleList().getFingerprint();
      chrs = Array.subArray(markerSet.getChrs(), markerIndicesInProject);
      positions = Array.subArray(markerSet.getPositions(), markerIndicesInProject);
      names = Array.subArray(markerSet.getMarkerNames(), markerIndicesInProject);
      file = new RandomAccessFile(currentMarkFilename, "r");
      file.read(parameterReadBuffer);
      nullStatus = parameterReadBuffer[TransposeData.MARKERDATA_NULLSTATUS_START];
      bytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
      numBytesPerMarker = bytesPerSampleMarker * proj.getSamples().length;
      numSamplesObserved = Compression.bytesToInt(parameterReadBuffer,
                                                  TransposeData.MARKERDATA_NUMSAMPLES_START);
      if (numSamplesObserved != proj.getSamples().length) {
        String error = "mismatched number of samples between sample list (n="
                       + proj.getSamples().length + ") and file '" + currentMarkFilename + "' (n="
                       + numSamplesObserved + ")";
        proj.getLog().reportTimeError(error);
        throw new IllegalStateException(error);
      }
      fingerprint = Compression.bytesToLong(parameterReadBuffer,
                                            TransposeData.MARKERDATA_FINGERPRINT_START);
      if (fingerprint != sampleFingerprint) {
        String error = "mismatched sample fingerprints between sample list and file '"
                       + currentMarkFilename + "'";
        proj.getLog().reportError(error);
        throw new IllegalStateException(error);
      }
      numBytesMarkernamesSection =
                                 Compression.bytesToInt(parameterReadBuffer,
                                                        TransposeData.MARKERDATA_MARKERNAMELEN_START);
      isGcNull = Sample.isGcNull(nullStatus);
      isXNull = Sample.isXNull(nullStatus);
      isYNull = Sample.isYNull(nullStatus);
      isBafNull = Sample.isBafNull(nullStatus);
      isLrrNull = Sample.isLrrNull(nullStatus);
      isGenotypeNull = Sample.isAbAndForwardGenotypeNull(nullStatus);
      isNegativeXYAllowed = Sample.isNegativeXOrYAllowed(nullStatus);
      if (new File(proj.MARKER_DATA_DIRECTORY.getValue(false, true) + "outliers.ser").exists()) {
        outlierHash = (Hashtable<String, Float>) SerializedFiles.readSerial(
                                                                            proj.MARKER_DATA_DIRECTORY.getValue(false,
                                                                                                                true)
                                                                            + "outliers.ser");
        if (debugMode) {
          proj.getLog()
              .reportTimeInfo("Loading RAF: " + currentMarkFilename + " and outliers "
                              + proj.MARKER_DATA_DIRECTORY.getValue(false, true) + "outliers.ser");
        }
      } else {
        outlierHash = new Hashtable<String, Float>();
      }
    }

    @Override
    public boolean hasNext() {
      return numLoaded < markersIndicesInFile.length;
    }

    @Override
    public Callable<MarkerData> next() {
      long seekLocation = (long) TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN
                          + (long) numBytesMarkernamesSection
                          + markersIndicesInFile[numLoaded] * (long) numBytesPerMarker;
      byte[] buffer = new byte[numBytesPerMarker];
      try {
        file.seek(seekLocation);
        file.read(buffer);
        MDLWorker worker = new MDLWorker(buffer, bytesPerSampleMarker,
                                         markerIndicesInProject[numLoaded], proj.getSamples(),
                                         names[numLoaded], chrs[numLoaded], positions[numLoaded],
                                         isGcNull, isXNull, isYNull, isBafNull, isLrrNull,
                                         isGenotypeNull, isNegativeXYAllowed, outlierHash,
                                         sampleFingerprint, debugMode);
        numLoaded++;
        buffer = null;
        return worker;

      } catch (IOException e) {
        proj.getLog().reportIOException(currentMarkFilename);
        e.printStackTrace();
      } catch (Exception e) {
        proj.getLog().reportException(e);
      }
      numLoaded = markerIndicesInProject.length;
      return null;
    };

    @Override
    public void shutdown() {
      if (file != null) {
        try {
          file.close();
        } catch (IOException e) {
          proj.getLog().reportIOException(currentMarkFilename);
          e.printStackTrace();
        }
      }
    }
  }

  /**
   * Callable class to handle the decompression in another thread
   *
   */
  private static class MDLWorker implements Callable<MarkerData> {
    private final byte[] buffer;
    private final byte bytesPerSampleMarker;
    private final int markersIndexInProject;
    private final String[] allSampsInProj;
    private final String markerName;
    private final byte chr;
    private final int pos;
    private final boolean isGcNull, isXNull, isYNull, isBafNull, isLrrNull, isGenotypeNull,
        isNegativeXYAllowed, debugMode;
    private final Hashtable<String, Float> outOfRangeValues;
    private final long fingerprint;

    private MDLWorker(byte[] buffer, byte bytesPerSampleMarker, int markersIndexInProject,
                      String[] allSampsInProj, String markerName, byte chr, int pos,
                      boolean isGcNull, boolean isXNull, boolean isYNull, boolean isBafNull,
                      boolean isLrrNull, boolean isGenotypeNull, boolean isNegativeXYAllowed,
                      Hashtable<String, Float> outOfRangeValues, long fingerprint,
                      boolean debugMode) {
      super();
      this.buffer = buffer;
      this.bytesPerSampleMarker = bytesPerSampleMarker;
      this.markersIndexInProject = markersIndexInProject;
      this.allSampsInProj = allSampsInProj;
      this.markerName = markerName;
      this.chr = chr;
      this.pos = pos;
      this.isGcNull = isGcNull;
      this.isXNull = isXNull;
      this.isYNull = isYNull;
      this.isBafNull = isBafNull;
      this.isLrrNull = isLrrNull;
      this.isGenotypeNull = isGenotypeNull;
      this.isNegativeXYAllowed = isNegativeXYAllowed;
      this.outOfRangeValues = outOfRangeValues;
      this.fingerprint = fingerprint;
      this.debugMode = debugMode;
    }

    @Override
    public MarkerData call() throws Exception {
      int numSamplesProj = allSampsInProj.length;
      float[] gcs = null;
      float[] xs = null;
      float[] ys = null;
      float[] bafs = null;
      float[] lrrs = null;
      byte[] abGenotypes = null;
      byte[] forwardGenotypes = null;
      byte[] genotypeTmp;

      int indexReadBuffer = 0;

      int indexStart = 0;
      indexReadBuffer = indexStart;
      // time = new Date().getTime();
      if (!isGcNull) {
        gcs = new float[numSamplesProj];
        for (int j = 0; j < numSamplesProj; j++) {
          gcs[j] = Compression.gcBafDecompress(new byte[] {buffer[indexReadBuffer],
                                                           buffer[indexReadBuffer + 1]});
          indexReadBuffer += bytesPerSampleMarker;
        }

        indexStart += 2;
      }
      indexReadBuffer = indexStart;
      if (!isXNull) {
        xs = new float[numSamplesProj];
        for (int j = 0; j < numSamplesProj; j++) {
          if (isNegativeXYAllowed) {
            xs[j] = Compression.xyDecompressAllowNegative(new byte[] {buffer[indexReadBuffer],
                                                                      buffer[indexReadBuffer + 1]});
          } else {
            xs[j] = Compression.xyDecompressPositiveOnly(new byte[] {buffer[indexReadBuffer],
                                                                     buffer[indexReadBuffer + 1]});
          }
          if (xs[j] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLAG_FLOAT) {
            // xs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tx");
            String key = markersIndexInProject + "\t" + allSampsInProj[j] + "\tx";
            if (debugMode) {
              System.err.println("loading outlier " + key);
            }
            xs[j] = outOfRangeValues.get(key);

          }
          indexReadBuffer += bytesPerSampleMarker;
        }

        indexStart += 2;
      }
      indexReadBuffer = indexStart;
      if (!isYNull) {
        ys = new float[numSamplesProj];
        for (int j = 0; j < numSamplesProj; j++) {
          if (isNegativeXYAllowed) {
            ys[j] = Compression.xyDecompressAllowNegative(new byte[] {buffer[indexReadBuffer],
                                                                      buffer[indexReadBuffer + 1]});
          } else {
            ys[j] = Compression.xyDecompressPositiveOnly(new byte[] {buffer[indexReadBuffer],
                                                                     buffer[indexReadBuffer + 1]});
          }
          if (ys[j] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLAG_FLOAT) {
            // ys[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\ty");
            String key = markersIndexInProject + "\t" + allSampsInProj[j] + "\ty";
            if (debugMode) {
              System.err.println("loading outlier " + key);
            }
            ys[j] = outOfRangeValues.get(key);

          }
          indexReadBuffer += bytesPerSampleMarker;
        }

        indexStart += 2;
      }
      indexReadBuffer = indexStart;
      if (!isBafNull) {
        bafs = new float[numSamplesProj];
        for (int j = 0; j < numSamplesProj; j++) {
          bafs[j] = Compression.gcBafDecompress(new byte[] {buffer[indexReadBuffer],
                                                            buffer[indexReadBuffer + 1]});
          indexReadBuffer += bytesPerSampleMarker;
        }

        indexStart += 2;
      }
      indexReadBuffer = indexStart;
      if (!isLrrNull) {
        lrrs = new float[numSamplesProj];
        for (int j = 0; j < numSamplesProj; j++) {
          lrrs[j] = Compression.lrrDecompress(new byte[] {buffer[indexReadBuffer],
                                                          buffer[indexReadBuffer + 1],
                                                          buffer[indexReadBuffer + 2]});
          if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLAG_FLOAT) {
            // lrrs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tlrr");
            String key = markersIndexInProject + "\t" + allSampsInProj[j] + "\tlrr";
            if (debugMode) {
              System.err.println("loading outlier " + key);
            }
            lrrs[j] = outOfRangeValues.get(markersIndexInProject + "\t" + allSampsInProj[j]
                                           + "\tlrr");

          }
          indexReadBuffer += bytesPerSampleMarker;
        }

        indexStart += 3;
      }
      indexReadBuffer = indexStart;
      if (!isGenotypeNull) {
        abGenotypes = new byte[numSamplesProj];
        forwardGenotypes = new byte[numSamplesProj];
        for (int j = 0; j < numSamplesProj; j++) {
          genotypeTmp = Compression.genotypeDecompress(buffer[indexReadBuffer]);
          abGenotypes[j] = genotypeTmp[0];
          forwardGenotypes[j] = genotypeTmp[1];
          indexReadBuffer += bytesPerSampleMarker;
        }
      }
      return new MarkerData(markerName, chr, pos, fingerprint, gcs, null, null, xs, ys, null, null,
                            bafs, lrrs, abGenotypes, forwardGenotypes);
    }
  }

  public static void main(String[] args) {
    Project proj = new Project(null, false);
    System.out.println(proj.getMarkerNames().length);
    String[] markers = proj.getMarkerNames();
    long totalTime = System.currentTimeMillis();
    String totalTimeMDL = "";
    String totalTimemarkerDataLoader = "";

    int iter = 1;

    for (int i = 0; i < iter; i++) {
      int index = 0;
      MDL mdl = new MDL(proj, null, proj.getMarkerNames(), 3, 10);
      while (mdl.hasNext()) {
        try {
          MarkerData markerData = mdl.next();
          if (!markerData.getMarkerName().equals(markers[index])) {
            System.err.println("DSKLFJSD\t" + markerData.getMarkerName() + "\t" + markers[index]);
          }
          if (index % 10000 == 0) {
            proj.getLog().reportTimeInfo(index + " of " + markers.length);
          }
        } catch (NullPointerException npe) {
          System.out.println(markers[index]);
          System.exit(1);
        }
        index++;
      }

      mdl.shutdown();
    }
    totalTimeMDL = "TIME:MDL" + ext.getTimeElapsed(totalTime);

    totalTime = System.currentTimeMillis();
    for (int i = 0; i < iter; i++) {
      MarkerDataLoader markerDataLoader =
                                        MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj,
                                                                                                proj.getMarkerNames());
      for (int j = 0; j < proj.getMarkerNames().length; j++) {
        MarkerData markerData = markerDataLoader.requestMarkerData(j);
        if (!markerData.getMarkerName().equals(markers[j])) {
          System.err.println("DSKLFJSD\t" + markerData.getMarkerName() + "\t" + markers[j]);
        }
        if (j % 10000 == 0) {
          proj.getLog().reportTimeInfo(j + " of " + markers.length);
        }
        markerDataLoader.releaseIndex(j);
      }
    }
    totalTimemarkerDataLoader = "TIME:MarkerDataloader" + ext.getTimeElapsed(totalTime);
    proj.getLog().reportTimeInfo(totalTimeMDL);
    proj.getLog().reportTimeInfo(totalTimemarkerDataLoader);
  }

}
