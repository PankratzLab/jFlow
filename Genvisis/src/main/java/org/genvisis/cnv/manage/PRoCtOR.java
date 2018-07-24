package org.genvisis.cnv.manage;

import java.io.File;
import java.io.FilenameFilter;
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
import org.genvisis.cnv.analysis.PennCNVPrep;
import org.genvisis.cnv.analysis.pca.PCA;
import org.genvisis.cnv.analysis.pca.PCAPrep;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsApply;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsCompute.PRE_PROCESSING_METHOD;
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
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.stats.LeastSquares.LS_TYPE;
import com.google.common.collect.ImmutableList;

// PRincipal COmponents Residuals - PR [o] C[t]O R
public class PRoCtOR {

  private static final String SHADOW_PREP_DIR = "shadowPrep/";

  private static long getMaxSampleSize(Project proj) {
    File[] sampleFiles = (new File(proj.SAMPLE_DIRECTORY.getValue())).listFiles(new FilenameFilter() {

      @Override
      public boolean accept(File dir, String name) {
        return name.endsWith(Sample.SAMPLE_FILE_EXTENSION);
      }
    });
    long max = 0;
    for (File f : sampleFiles) {
      if (f.length() > max) {
        max = f.length();
      }
    }
    return max;
  }

  private static int getSampleChunks(Project proj, int numThreads) {
    long mem = Runtime.getRuntime().maxMemory();
    long samp = getMaxSampleSize(proj);
    // With the assumption that all samples in a given chunk can be open when a chunk is processed
    // by a thread.
    // Thus we want a number of chunks such that we will not run out of memory if each thread is
    // simultaneously
    // processing chunks of (chunk size * max_chunk_size) files.
    // A fraction of max memory is used to account for additional files which may need to be in
    // memory during
    // this process (e.g. marker data)
    double sampleChunks = (0.65 * mem) / (numThreads * samp);
    return (int) sampleChunks;
  }

  public static String shadow(Project proj, String tmpDir, String outputBase,
                              double markerCallRateFilter, boolean recomputeLRR_PCs,
                              CORRECTION_TYPE correctionType, CHROMOSOME_X_STRATEGY strategy,
                              int numComponents, int totalThreads,
                              boolean setupCNVCallingIfSuccessfull) {
    int numMarkerThreads = 1;
    int numThreads = (int) Math.ceil((double) totalThreads / (double) numMarkerThreads);
    boolean markerQC = true;
    String useFile = null;
    int sampleChunks = getSampleChunks(proj, numThreads);
    proj.getLog().report("Using " + sampleChunks + " sample chunks");

    int retCode = PCAPrep.prepPCA(proj, numThreads, outputBase, markerQC, markerCallRateFilter,
                                  useFile, proj.getSampleList(), proj.getLog());
    if (retCode != 42) {
      return PCAPrep.errorMessage(retCode);
    }
    PrincipalComponentsApply pcApply = PCA.generateFullPCA(proj, numComponents, outputBase,
                                                           recomputeLRR_PCs, true, null,
                                                           PRE_PROCESSING_METHOD.NONE,
                                                           proj.getLog());

    proj.getLog().reportTime("Setting PCs file: " + pcApply.getExtrapolatedPCsFile());
    proj.INTENSITY_PC_FILENAME.setValue(pcApply.getExtrapolatedPCsFile());

    if (correctionType == CORRECTION_TYPE.GENERATE_PCS_ONLY) {
      return "";
    }
    String projectDirectory = PennCNVPrep.getCorrectedProjectDirectory(proj, numComponents,
                                                                       correctionType, strategy);
    String markerDirectory = projectDirectory + proj.MARKER_DATA_DIRECTORY.getDefaultValue();

    PennCNVPrep.prepExport(proj, markerDirectory, numComponents, null, numThreads, numMarkerThreads,
                           LS_TYPE.REGULAR, false, correctionType, strategy);
    PennCNVPrep.exportSpecialPennCNV(proj, SHADOW_PREP_DIR, tmpDir, numComponents, null, numThreads,
                                     numMarkerThreads, true, LS_TYPE.REGULAR, sampleChunks, false,
                                     false, correctionType, strategy);
    if (setupCNVCallingIfSuccessfull) {
      GenvisisWorkflow.setupCNVCalling(projectDirectory);
    }
    return "";
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

    public void write(MarkerData markerData, boolean canXYBeNegative) throws IOException, Elision {
      String mdrafName = markerLookup.get(markerData.getMarkerName());
      synchronized (mdrafName) {
        RandomAccessFile mdraf = rafMap.get(mdrafName);
        String[] mkrNmArr = mkrNames.get(mdrafName);
        if (mdraf == null) {
          byte nullStatus = Sample.computeNullStatus(markerData.getGCs(), markerData.getXs(),
                                                     markerData.getYs(), markerData.getBAFs(),
                                                     markerData.getLRRs(),
                                                     markerData.getAbGenotypes(),
                                                     markerData.getForwardGenotypes(),
                                                     canXYBeNegative);
          mdraf = openMDRAF(outDir + mdrafName, numInd, nullStatus, fingerprint, mkrNmArr);
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
                                        statMap.get(mdrafName), oorTables.get(mdrafName),
                                        canXYBeNegative));
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
      new File(outDir).mkdirs();
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
                                    boolean[] samplesToUseCluster, String[] markers,
                                    CORRECTION_TYPE correctionType,
                                    CHROMOSOME_X_STRATEGY sexStrategy, int numComponents,
                                    int numCorrectionThreads, int numMarkerThreads) {

    Project shadowProject = new Project();
    // TODO update shadow project for new location of files,
    // transposed/samples dirs, etc
    ShadowMarkerDataWriter smdw = new ShadowMarkerDataWriter();
    smdw.setOutputDirectory(newTransposedDir);
    smdw.setupMarkerFiles(proj);
    Files.copyFile(proj.SAMPLELIST_FILENAME.getValue(),
                   shadowProject.SAMPLELIST_FILENAME.getValue());

    PcCorrectionProducer producer = new PcCorrectionProducer(principalComponentsResiduals,
                                                             numComponents, sampleSex,
                                                             samplesToUseCluster, LS_TYPE.REGULAR,
                                                             numCorrectionThreads, 1, markers,
                                                             correctionType, sexStrategy);
    try (WorkerTrain<PrincipalComponentsIntensity> train = new WorkerTrain<>(producer,
                                                                             numMarkerThreads, 10,
                                                                             proj.getLog())) {
      ArrayList<String> notCorrected = new ArrayList<>();
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
            smdw.write(markerData, proj.getArrayType().getCanXYBeNegative());
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
            smdw.write(markerData, proj.getArrayType().getCanXYBeNegative());
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
    }

    TransposeData.reverseTranspose(shadowProject);

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "D:/projects/Poynter.properties";
    String tempDir = null;
    String outputBase = MitoPipeline.FILE_BASE;
    double callrate = MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER;
    boolean recomputeLRR = false;
    int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
    CORRECTION_TYPE correctionType = CORRECTION_TYPE.XY;
    CHROMOSOME_X_STRATEGY strategy = CHROMOSOME_X_STRATEGY.BIOLOGICAL;
    int numThreads = Runtime.getRuntime().availableProcessors();
    boolean callCNVs = false;

    String usage = "\n" + "cnv.manage.PRoCtOR requires 0-1 arguments\n"
                   + "   (1) project properties filename (i.e. proj=" + filename + " (default))\n"
                   + "   (2) Number of principal components for correction (i.e. numComponents="
                   + numComponents + " (default))\n"
                   + "   (3) Output file full path and baseName for principal components correction files (i.e. outputBase="
                   + outputBase + " (default))\n"
                   + "   (4) Call-rate filter for determining high-quality markers (i.e. callrate="
                   + callrate + " (default))\n"
                   + "   (5) Flag specifying whether or not to re-compute Log-R Ratio values (usually false if LRRs already exist) (i.e. recomputeLRR="
                   + recomputeLRR + " (default))\n"
                   + "   (6) Type of correction.  Options include: "
                   + ArrayUtils.toStr(CORRECTION_TYPE.values(), ", ") + " (i.e. type="
                   + correctionType + " (default))\n"

                   + "   (7) Chromosome X correction strategy.  Options include: "
                   + ArrayUtils.toStr(CHROMOSOME_X_STRATEGY.values(), ", ") + " (i.e. sexStrategy="
                   + strategy + " (default))\n"

                   + "   (8) Total number of threads to use (i.e. numThreads=" + numThreads
                   + " (default))\n"

                   + "   (8) OPTIONAL: temp directory for intermediate files (which tend to be very large) (i.e. tmp="
                   + tempDir + " (default))\n"
                   + "   (9) OPTIONAL: Create a script with the commands required to process the corrected data and call CNVs (i.e. -callCNVs (not the default))\n"
                   + "" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("numComponents=")) {
        numComponents = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("outputBase=")) {
        outputBase = ext.parseStringArg(arg, outputBase);
        numArgs--;
      } else if (arg.startsWith("callrate=")) {
        callrate = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("recomputeLRR=")) {
        recomputeLRR = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("numThreads=")) {
        numThreads = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("type=")) {
        correctionType = CORRECTION_TYPE.valueOf(ext.parseStringArg(arg,
                                                                    correctionType.toString()));
        numArgs--;
      } else if (arg.startsWith("sexStrategy=")) {
        strategy = CHROMOSOME_X_STRATEGY.valueOf(ext.parseStringArg(arg, strategy.toString()));
        numArgs--;
      } else if (arg.startsWith("tmp=")) {
        tempDir = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("-callCNVs")) {
        callCNVs = true;
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
      Project proj = new Project(filename);
      String err = shadow(proj, tempDir, outputBase, callrate, recomputeLRR, correctionType,
                          strategy, numComponents, numThreads, callCNVs);
      if (!"".equals(err)) {
        System.err.println("Error - " + err);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
