package org.genvisis.cnv.manage;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.WorkerTrain;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;
import org.pankratzlab.common.ext;

import com.google.common.collect.ImmutableList;

public class MDL implements Iterator<MarkerData> {

  private final Project proj;
  private final List<Marker> markersToLoad;
  private final int numDecompressThreads;
  private WorkerTrain<MarkerData> decompTrain;
  private final MarkerLookup markerLookup;
  private final ArrayList<FileMatch> files;
  private int numLoaded, fileIndex, numLoadedForFile;
  private final int markerBuffer;
  private final Set<String> missing;
  private BufferReader producer;
  private int reportEvery;
  private final long startTime;
  private boolean debugMode;
  private String currentFile;

  /**
   * @param proj
   * @param markersToLoad them to load
   */
  public MDL(Project proj, List<Marker> markersToLoad) {
    this(proj, markersToLoad, 2, 100);
  }

  /**
   * @param proj
   * @param markerNames them to load
   */
  public MDL(Project proj, String[] markerNames) {
    this(proj, markerNames, 2, 100);
  }

  /**
   * @param proj
   * @param markerNames them to load
   * @param numDecompressThreads number of threads used to decompress the marker
   * @param markerBuffer number of markers to hold in the queue for processing
   */
  public MDL(Project proj, String[] markerNames, int numDecompressThreads, int markerBuffer) {
    this(proj, Arrays.stream(markerNames).map(proj.getMarkerSet().getMarkerNameMap()::get)
                     .collect(ImmutableList.toImmutableList()),
         numDecompressThreads, markerBuffer);
  }

  /**
   * @param proj
   * @param markersToLoad them to load
   * @param numDecompressThreads number of threads used to decompress the marker
   * @param markerBuffer number of markers to hold in the queue for processing
   */
  public MDL(Project proj, List<Marker> markersToLoad, int numDecompressThreads, int markerBuffer) {
    this.proj = proj;
    missing = new HashSet<>();
    this.markersToLoad = markersToLoad;
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
      decompTrain.close();
    }
    producer = new BufferReader(proj, match, debugMode);
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
    decompTrain = new WorkerTrain<>(producer, numDecompressThreads, markerBuffer, proj.getLog());
  }

  @Override
  public boolean hasNext() {
    boolean hasNext = numLoaded < markersToLoad.size();
    if (!hasNext) {
      shutdown();
    }
    return hasNext;
  }

  public void shutdown() {
    if (decompTrain != null) {
      decompTrain.close();
    }
  }

  @Override
  public MarkerData next() throws IllegalStateException {
    MarkerData markerData = null;
    try {
      if (numLoaded == 0 && fileIndex == 0) {
        initTrain(files.get(fileIndex));
        numLoadedForFile = 0;
      } else if (numLoadedForFile == files.get(fileIndex).getMatches().size()) {
        fileIndex++;
        initTrain(files.get(fileIndex));
        numLoadedForFile = 0;
      }
      decompTrain.hasNext();
      if (!decompTrain.hasNext()) {
        String error = "Internal index error, halting";
        proj.getLog().reportError(error);
        throw new IllegalStateException(error);
      }
      numLoadedForFile++;
      numLoaded++;
      if ((startTime - new Date().getTime()) % reportEvery == 0) {
        proj.getLog()
            .reportTimeInfo("Loaded " + numLoaded + " of " + markersToLoad.size() + " markers");
      }
      markerData = decompTrain.next();
      if (!markerData.getMarkerName().equals(markersToLoad.get(numLoaded - 1).getName())) {
        String error = "Internal index error - mismatched marker data found, halting";
        proj.getLog().reportError(error);
        throw new IllegalStateException(error);
      }
    } catch (NullPointerException npe) {
      proj.getLog().reportError("Could not load " + markersToLoad.get(numLoaded - 1).getName()
                                + ", this is usually caused by a corrupt or missing outlier file");

    }
    return markerData;
  }

  @Override
  public void remove() {
    // TODO Auto-generated method stub

  }

  /**
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
   */
  private static class FileMatch {

    private static class Match {

      private final Marker marker;
      private final int fileIndex;

      /**
       * @param marker
       * @param fileIndex
       */
      private Match(Marker marker, int fileIndex) {
        super();
        this.marker = marker;
        this.fileIndex = fileIndex;
      }

      /**
       * @return the marker
       */
      public Marker getMarker() {
        return marker;
      }

      /**
       * @return the fileIndex
       */
      public int getFileIndex() {
        return fileIndex;
      }

    }

    private final String fileName;
    private final List<Match> matches;

    private FileMatch(String fileName) {
      super();
      this.fileName = fileName;
      matches = new ArrayList<>(100000);
    }

    public String getFileName() {
      return fileName;
    }

    /**
     * @return the matches
     */
    public List<Match> getMatches() {
      return Collections.unmodifiableList(matches);
    }

    private void add(Marker marker, int fileIndex) {
      matches.add(new Match(marker, fileIndex));
    }

  }

  /**
   * @return {@link FileMatch} objects for determining order and markers to scan from each file
   */
  private ArrayList<FileMatch> matchFileNames() {
    ArrayList<FileMatch> files = new ArrayList<>();
    String currentFile = "";
    int currentIndex = -1;
    for (Marker marker : markersToLoad) {
      if (markerLookup.contains(marker.getName())) {
        String[] line = markerLookup.get(marker.getName()).split(PSF.Regex.GREEDY_WHITESPACE);
        if (!line[0].equals(currentFile)) {
          currentFile = line[0];
          files.add(new FileMatch(proj.MARKER_DATA_DIRECTORY.getValue(false, true) + currentFile));
          currentIndex++;
        }
        files.get(currentIndex).add(marker, Integer.parseInt(line[1]));

      } else {
        missing.add(marker.getName());
      }
    }
    if (!missing.isEmpty()) {
      proj.getLog().reportError("Could not find the following markers in the project"
                                + ArrayUtils.toStr(missing, "\n"));
    }
    return files;
  }

  /**
   * {@link WorkerTrain.Producer} implementation to serve up marker data
   */
  private static class BufferReader extends AbstractProducer<MarkerData> {

    private final byte[] parameterReadBuffer = new byte[TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN];
    private RandomAccessFile file;
    private byte nullStatus;
    private byte bytesPerSampleMarker;
    private int numBytesPerMarker, numSamplesObserved, numBytesMarkernamesSection, numLoaded;
    private long fingerprint, sampleFingerprint;
    private final FileMatch fileMatch;
    private final Project proj;
    private boolean isGcNull, isXNull, isYNull, isBafNull, isLrrNull, isGenotypeNull,
        isNegativeXYAllowed;
    private final boolean debugMode;
    private Hashtable<String, Float> outlierHash;
    private final MarkerDetailSet markerSet;

    private BufferReader(Project proj, FileMatch fileMatch, boolean debugMode) {
      super();
      this.proj = proj;
      this.markerSet = proj.getMarkerSet();
      this.fileMatch = fileMatch;
      numLoaded = 0;
      this.debugMode = debugMode;
    }

    @SuppressWarnings("unchecked")
    private void init() throws IOException, IllegalStateException {
      sampleFingerprint = proj.getSampleList().getFingerprint();
      file = new RandomAccessFile(fileMatch.getFileName(), "r");
      file.read(parameterReadBuffer);
      nullStatus = parameterReadBuffer[TransposeData.MARKERDATA_NULLSTATUS_START];
      bytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
      numBytesPerMarker = bytesPerSampleMarker * proj.getNumberOfParsedSamples();
      numSamplesObserved = Compression.bytesToInt(parameterReadBuffer,
                                                  TransposeData.MARKERDATA_NUMSAMPLES_START);
      if (numSamplesObserved != proj.getNumberOfParsedSamples()) {
        String error = "mismatched number of samples between sample list (n="
                       + proj.getNumberOfParsedSamples() + ") and file '" + fileMatch.getFileName()
                       + "' (n=" + numSamplesObserved + ")";
        proj.getLog().reportError(error);
        throw new IllegalStateException(error);
      }
      fingerprint = Compression.bytesToLong(parameterReadBuffer,
                                            TransposeData.MARKERDATA_FINGERPRINT_START);
      if (fingerprint != sampleFingerprint) {
        String error = "Error - mismatched sample fingerprints between sample list ("
                       + sampleFingerprint + ") and file '" + file + "' (" + fingerprint + ")";
        proj.getLog().reportError(error);
        throw new IllegalStateException(error);
      }
      numBytesMarkernamesSection = Compression.bytesToInt(parameterReadBuffer,
                                                          TransposeData.MARKERDATA_MARKERNAMELEN_START);
      isGcNull = Sample.isGcNull(nullStatus);
      isXNull = Sample.isXNull(nullStatus);
      isYNull = Sample.isYNull(nullStatus);
      isBafNull = Sample.isBafNull(nullStatus);
      isLrrNull = Sample.isLrrNull(nullStatus);
      isGenotypeNull = Sample.isAbAndForwardGenotypeNull(nullStatus);
      isNegativeXYAllowed = Sample.isNegativeXOrYAllowed(nullStatus);
      if (new File(proj.MARKER_DATA_DIRECTORY.getValue(false, true) + "outliers.ser").exists()) {
        outlierHash = (Hashtable<String, Float>) SerializedFiles.readSerial(proj.MARKER_DATA_DIRECTORY.getValue(false,
                                                                                                                true)
                                                                            + "outliers.ser");
        if (debugMode) {
          proj.getLog()
              .reportTimeInfo("Loading RAF: " + fileMatch.getFileName() + " and outliers "
                              + proj.MARKER_DATA_DIRECTORY.getValue(false, true) + "outliers.ser");
        }
      } else {
        outlierHash = new Hashtable<>();
      }
    }

    @Override
    public boolean hasNext() {
      return numLoaded < fileMatch.getMatches().size();
    }

    @Override
    public Callable<MarkerData> next() {
      FileMatch.Match nextMatch = fileMatch.getMatches().get(numLoaded);
      Marker nextMarker = nextMatch.getMarker();
      long seekLocation = (long) TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN
                          + (long) numBytesMarkernamesSection
                          + nextMatch.getFileIndex() * (long) numBytesPerMarker;
      byte[] buffer = new byte[numBytesPerMarker];
      try {
        file.seek(seekLocation);
        file.read(buffer);
        MDLWorker worker = new MDLWorker(buffer, bytesPerSampleMarker,
                                         markerSet.getMarkerIndexMap().get(nextMarker),
                                         proj.getSamples(), nextMarker.getName(),
                                         nextMarker.getChr(), nextMarker.getPosition(), isGcNull,
                                         isXNull, isYNull, isBafNull, isLrrNull, isGenotypeNull,
                                         isNegativeXYAllowed, outlierHash, sampleFingerprint,
                                         debugMode);
        numLoaded++;
        buffer = null;
        return worker;

      } catch (IOException e) {
        proj.getLog().reportIOException(fileMatch.getFileName());
        e.printStackTrace();
      } catch (Exception e) {
        proj.getLog().reportException(e);
      }
      numLoaded = fileMatch.getMatches().size();
      return null;
    };

    @Override
    public void shutdown() {
      if (file != null) {
        try {
          file.close();
        } catch (IOException e) {
          proj.getLog().reportIOException(fileMatch.getFileName());
          e.printStackTrace();
        }
      }
    }
  }

  /**
   * Callable class to handle the decompression in another thread
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
    Project proj = new Project(null);
    System.out.println(proj.getMarkerNames().length);
    String[] markers = proj.getMarkerNames();
    long totalTime = System.currentTimeMillis();
    String totalTimeMDL = "";
    String totalTimemarkerDataLoader = "";

    int iter = 1;

    for (int i = 0; i < iter; i++) {
      int index = 0;
      MDL mdl = new MDL(proj, proj.getMarkerNames(), 3, 10);
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
      MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj,
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
