package org.genvisis.one.JL.shadow;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.concurrent.ConcurrentHashMap;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CHROMOSOME_X_STRATEGY;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CORRECTION_TYPE;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.PcCorrectionProducer;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsResiduals;
import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerTrain;
import org.genvisis.stats.LeastSquares.LS_TYPE;
import com.google.common.collect.ImmutableList;

/**
 * Draft of a possibly a more friendly,less temporary version of shadowing Needs a few magic methods
 */
public class ShadowRework {

  private ShadowRework() {

  }

  static class ShadowMarkerDataWriter {

    int numInd = -1;
    long fingerprint;
    String outDir;

    HashMap<String, String> markerLookup = new HashMap<>();
    HashMap<String, Integer> markerIndexLocal = new HashMap<>();

    Map<String, Integer> mkrInds;
    String[] samples;

    Map<String, String[]> mkrNames = new ConcurrentHashMap<>();
    Map<String, RandomAccessFile> rafMap = new ConcurrentHashMap<>();
    Map<String, Byte> statMap = new ConcurrentHashMap<>();
    Map<String, Hashtable<String, Float>> oorTables = new ConcurrentHashMap<>();

    HashMap<Integer, Map<String, String>> chrFileMap = new HashMap<>();

    // proj is old proj
    public void setupMarkerFiles(Project proj) {
      samples = proj.getSamples();
      mkrInds = proj.getMarkerIndices();
      numInd = proj.getSamples().length;
      fingerprint = MarkerSet.fingerprintForMarkers(proj);
      NavigableMap<Byte, NavigableSet<Marker>> chrMap = proj.getMarkerSet().getChrMap();
      int numMarkers = 2500;
      for (Byte b : chrMap.keySet()) {
        List<Marker> mkrs = ImmutableList.copyOf(chrMap.get(b));
        int[] indLists = ArrayUtils.splitUp(mkrs.size(), (mkrs.size() / numMarkers) + 1);
        int total = 0;
        for (int cnt : indLists) {
          String mkrFile = getMDRAFName(b, total, total + cnt);
          String[] mkrNmArr = new String[cnt];
          for (int i = total; i < total + cnt; i++) {
            markerLookup.put(mkrs.get(i).getName(), mkrFile);
            markerIndexLocal.put(mkrs.get(i).getName(), i - total);
            mkrNmArr[i - total] = mkrs.get(i).getName();
          }
          mkrNames.put(mkrFile, mkrNmArr);
          oorTables.put(mkrFile, new Hashtable<>());
          Map<String, String> files = chrFileMap.get((int) b);
          if (files == null) {
            files = new HashMap<>();
            chrFileMap.put((int) b, files);
          }
          files.put(total + "\t" + (total + cnt), mkrFile);
          total = total + cnt;
        }
      }
    }

    public void write(MarkerData markerData) throws IOException, Elision {
      String mdrafName = markerLookup.get(markerData.getMarkerName());
      synchronized (mdrafName) {
        RandomAccessFile mdraf = rafMap.get(mdrafName);
        String[] mkrNmArr = mkrNames.get(mdrafName);
        if (mdraf == null) {
          byte nullStatus = Sample.updateNullStatus(markerData.getGCs(), null, null,
                                                    markerData.getBAFs(), markerData.getLRRs(),
                                                    markerData.getAbGenotypes(),
                                                    markerData.getForwardGenotypes(), false);
          mdraf = openMDRAF(mdrafName, numInd, nullStatus, fingerprint, mkrNmArr);
          rafMap.put(mdrafName, mdraf);
          statMap.put(mdrafName, nullStatus);
        }

        byte[] mkrBytes = Compression.objToBytes(mkrNames);
        int numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(statMap.get(mdrafName));
        int numBytesPerMarker = numBytesPerSampleMarker * numInd;
        long seek = TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN + mkrBytes.length
                    + markerIndexLocal.get(markerData.getMarkerName()) * numBytesPerMarker;
        // seek to location of marker in file, as we may be writing out of order
        mdraf.seek(seek);
        mdraf.write(markerData.compress(markerIndexLocal.get(markerData.getMarkerName()),
                                        statMap.get(mdrafName), oorTables.get(mdrafName), false));
      }
    }

    private RandomAccessFile openMDRAF(String filename, int nInd, byte nullStatus, long fingerprint,
                                       String[] mkrNames) throws IOException {
      byte[] mkrBytes = Compression.objToBytes(mkrNames);
      byte[] mdRAFHeader = TransposeData.getParameterSectionForMdRaf(nInd, mkrNames.length,
                                                                     nullStatus, fingerprint,
                                                                     mkrBytes);
      mkrBytes = null;

      RandomAccessFile mdRAF = new RandomAccessFile(filename, "rw");
      mdRAF.write(mdRAFHeader);
      mdRAFHeader = null;

      return mdRAF;
    }

    public void setOutputDirectory(String newTransposedDir) {
      this.outDir = newTransposedDir;
    }

    public void writeOutliers() throws IOException {
      Hashtable<String, Float> allOutliers = new Hashtable<>();
      for (Entry<String, Hashtable<String, Float>> oorEntry : oorTables.entrySet()) {
        RandomAccessFile mdRAF = rafMap.get(oorEntry.getKey());
        byte[] oorBytes = Compression.objToBytes(oorEntry.getValue());
        mdRAF.write(Compression.intToBytes(oorBytes.length));
        mdRAF.write(oorBytes);

        for (Entry<String, Float> entry : oorEntry.getValue().entrySet()) {
          String[] pts = entry.getKey().split("\t");
          int mkrInd = Integer.parseInt(pts[0]);
          int sampInd = Integer.parseInt(pts[1]);
          allOutliers.put(mkrInds.get(mkrNames.get(oorEntry.getKey())[mkrInd]) + "\t"
                          + samples[sampInd] + "\t" + pts[2], entry.getValue());
        }

      }

      SerializedFiles.writeSerial(allOutliers, outDir + "outliers.ser");
    }

    public void close() throws IOException {
      for (RandomAccessFile raf : rafMap.values()) {
        raf.close();
      }
    }
  }

  private static String getMDRAFName(int chr, int start, int end) {
    return "markers." + chr + "." + start + "." + end + MarkerData.MARKER_DATA_FILE_EXTENSION;
  }

  /**
   * @param proj Project to correct
   * @param principalComponentsResiduals PCs to do the correcting
   * @param preserveBafs preserve BAF values (NGS specific), you likely want false here
   * @param sampleSex for Sex specific clustering
   * @param samplesToUseCluster samples to seed correction
   * @param lType
   * @param numComponents number of PCs to correct for
   * @param numCorrectionThreads number of threads within a marker (max of 6 can be utilized)
   * @param numMarkerThreads number of markers corrected at once
   */
  public static void correctProject(Project proj, String newTransposedDir,
                                    PrincipalComponentsResiduals principalComponentsResiduals,
                                    boolean preserveBafs, int[] sampleSex,
                                    boolean[] samplesToUseCluster, CORRECTION_TYPE correctionType,
                                    CHROMOSOME_X_STRATEGY sexStrategy, int numComponents,
                                    int numCorrectionThreads, int numMarkerThreads) {

    Project shadowProject = new Project();
    // TODO update shadow project for new location of files,
    // transposed/samples dirs, etc
    ShadowMarkerDataWriter smdw = new ShadowMarkerDataWriter();
    smdw.setOutputDirectory(newTransposedDir);
    smdw.setupMarkerFiles(proj);

    String[] markers = proj.getMarkerNames(); // Correct the entire thing
    PcCorrectionProducer producer = new PcCorrectionProducer(principalComponentsResiduals,
                                                             numComponents, sampleSex,
                                                             samplesToUseCluster, LS_TYPE.REGULAR,
                                                             numCorrectionThreads, 1,
                                                             proj.getMarkerNames(), correctionType,
                                                             sexStrategy);
    ArrayList<String> notCorrected = new ArrayList<>();
    try (WorkerTrain<PrincipalComponentsIntensity> train = new WorkerTrain<>(producer,
                                                                             numMarkerThreads, 10,
                                                                             proj.getLog())) {
      int index = 0;
      while (train.hasNext()) {
        PrincipalComponentsIntensity principalComponentsIntensity = train.next();
        MarkerData markerData = principalComponentsIntensity.getCentroidCompute().getMarkerData();
        try {
          if (principalComponentsIntensity.isFail()) {
            notCorrected.add(markers[index]);
            /*
             * MDRAF requires knowing # of markers beforehand; this would require a double-pass (to
             * determine # successfully corrected) rather than streaming approach. Instead, either
             * write original data or write missing / dummy data.
             */
            smdw.write(markerData);
          } else {
            byte[] abGenotypes = principalComponentsIntensity.getGenotypesUsed();// for
            // now
            float[][] correctedXY = principalComponentsIntensity.getCorrectedIntensity(PrincipalComponentsIntensity.XY_RETURN,
                                                                                       true);
            float[][] correctedLRRBAF = principalComponentsIntensity.getCorrectedIntensity(PrincipalComponentsIntensity.BAF_LRR_RETURN,
                                                                                           true);
            markerData = new MarkerData(markerData.getMarkerName(), markerData.getChr(),
                                        markerData.getPosition(), markerData.getFingerprint(),
                                        markerData.getGCs(), null, null, correctedXY[0],
                                        correctedXY[1], null, null,
                                        preserveBafs ? markerData.getBAFs() : correctedLRRBAF[0],
                                        correctedLRRBAF[1], abGenotypes, abGenotypes);
            smdw.write(markerData);
          }
        } catch (IOException e) {
          proj.getLog().reportException(e);
          System.exit(1);
        } catch (Elision e) {
          proj.getLog().reportException(e);
          System.exit(1);
        }

        index++;
      }
    }

    try {
      smdw.writeOutliers();
      smdw.close();
    } catch (IOException e) {
      proj.getLog().reportException(e);
      System.exit(1);
    }

    if (!notCorrected.isEmpty()) {
      Files.writeArray(notCorrected.toArray(new String[notCorrected.size()]),
                       shadowProject.PROJECT_DIRECTORY.getValue() + notCorrected.size()
                                                                              + "_markersThatFailedCorrection.txt");
    }

    // call reverseTranspose next
  }
}
