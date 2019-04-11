package org.genvisis.cnv.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Timer;
import java.util.TimerTask;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;
import org.genvisis.cnv.filesys.AllelePair;
import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.seq.GenomeBuild;
import org.genvisis.seq.manage.mosdepth.AbstractParsingPipeline;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomicPosition;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.utils.gwas.DosageData;
import org.pankratzlab.utils.gwas.bgen.BGENBitMath;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.variant.variantcontext.Allele;

public class UKBBParsingPipeline extends AbstractParsingPipeline {

  int maxMDRAFThreadsPerChr = 8;
  int maxChrThreads = 4;
  int logEveryNMins = 5;

  // AB fields in RAF files are only 2 bits and can't hold full genotype values
  final boolean importGenoAsAB = true;

  Project proj;
  String sourceDir;
  HashMap<Integer, FileSet> fileSets;
  AtomicLong numberOfMarkersInProj = new AtomicLong(0);
  String famFile;
  String annotFile;
  String[][] famData;
  String[] allMarkers;
  Map<String, Marker> markerMap;

  public UKBBParsingPipeline() {
    super(2000, 5000);
  }

  public void setSourceDir(String sourceDir2) {
    sourceDir = sourceDir2;
  }

  public void setFamFile(String famFile2) {
    famFile = famFile2;
  }

  public void setAnnotationCSV(String csv) {
    annotFile = csv;
  }

  private void setThreadsPerChromosomeCount(int threadsPerChr) {
    maxMDRAFThreadsPerChr = threadsPerChr;
  }

  private void setChromosomeThreadCount(int chrThreads) {
    maxChrThreads = chrThreads;
  }

  public void runPipeline() {
    log.reportTime("Starting UKBioBank Parsing Pipeline...");
    createProject();
    createSampleList();
    createPED();
    createSampleData();
    discoverFilesets();
    if (fileSets.size() == 0) {
      log.reportError("No chrs found to parse!  Shutting down.");
      return;
    }
    createMarkerPositions();
    writeMarkerDetailSet();
    parseMarkerData();

    writeLookup();
    proj.writeMarkerSet();
    buildOutliers();
    createSampRAFsFromMDRAFs();

  }

  protected void setAdditionalProjectProperties() {
    proj.SOURCE_DIRECTORY.setValue(sourceDir);
    proj.SOURCE_FILENAME_EXTENSION.setValue("NULL");
    proj.SOURCE_FILENAME_EXTENSION.setValue("NULL");
    proj.GENOME_BUILD_VERSION.setValue(GenomeBuild.HG19);
    proj.ARRAY_TYPE.setValue(ARRAY.AFFY_AXIOM);
  }

  protected String[] parseSamples() {
    famData = HashVec.loadFileToStringMatrix(famFile, false, new int[] {0, 1, 2, 3, 4, 5});
    return Matrix.extractColumn(famData, 1);
  }

  protected int getNumSamples() {
    return famData.length;
  }

  protected int getNumMarkers() {
    return allMarkers.length;
  }

  protected void doWriteLookup() {
    Hashtable<String, String> lookup = new Hashtable<>();
    for (FileSet fs : fileSets.values()) {
      String mdRAF;
      List<String[]> mkrBins = fs.mkrBins;
      int start = 0;
      int end = 0;
      for (String[] bin : mkrBins) {
        end = start + bin.length;
        mdRAF = getMDRAFName(fs.chr, start, end);
        for (int m = 0; m < bin.length; m++) {
          lookup.put(bin[m], mdRAF + "\t" + m);
        }
        start = end;
      }
    }
    new MarkerLookup(lookup).serialize(proj.MARKERLOOKUP_FILENAME.getValue());
  }

  protected byte getNullStatus() {
    return Sample.computeNullStatus(new float[0], new float[0], new float[0], new float[0],
                                    new float[0], importGenoAsAB ? new byte[0] : null,
                                    importGenoAsAB ? null : new byte[0], false);
  }

  protected void createPED() {
    if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
      long t1 = System.nanoTime();
      PrintWriter writer = Files.getAppropriateWriter(proj.PEDIGREE_FILENAME.getValue());
      for (String[] ind : famData) {
        StringBuilder sb = new StringBuilder(ind[0]).append("\t");
        sb.append(ind[1]).append("\t");
        sb.append(ind[2]).append("\t");
        sb.append(ind[3]).append("\t");
        sb.append(ind[4]).append("\t");
        sb.append("0").append("\t");
        sb.append(ind[1]);
        writer.println(sb.toString());
      }
      writer.flush();
      writer.close();
      log.reportTime("Created ped file in " + ext.getTimeElapsedNanos(t1));
    } else {
      log.reportTime("Project sample list file already exists; skipping creation.");
    }
  }

  protected void createSampleData() {
    if (!Files.exists(proj.SAMPLE_DATA_FILENAME.getValue())) {
      try {
        SampleData.createSampleData(proj.PEDIGREE_FILENAME.getValue(), null, proj);
      } catch (Elision e) {
        e.printStackTrace();
        System.exit(1);
      }
    } else {
      log.reportTime("Project SampleData file already exists; skipping creation.");
    }
  }

  protected void discoverFilesets() {
    List<File> allSrcFiles = discoverAllSourceFiles(sourceDir, log);

    fileSets = new HashMap<>();

    for (File srcF : allSrcFiles) {
      String name = srcF.getName();
      String[] parts = name.split("_");
      int chr = Positions.chromosomeNumber(parts[2]);
      FileSet fs = fileSets.get(chr);
      if (fs == null) {
        fs = new FileSet();
        fs.log = this.log;
        fs.chr = chr;
        fs.famFile = famFile;
        fileSets.put(chr, fs);
      }

      switch (parts[1]) {
        case "l2r":
          if (fs.lrrFile != null) {
            log.reportError("Identified assumed duplicate " + "LRR" + " files: {" + fs.lrrFile
                            + " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
            System.exit(1);
          } else {
            fs.lrrFile = srcF.getAbsolutePath();
          }
          break;
        case "baf":
          if (fs.bafFile != null) {
            log.reportError("Identified assumed duplicate " + "BAF" + " files: {" + fs.bafFile
                            + " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
            System.exit(1);
          } else {
            fs.bafFile = srcF.getAbsolutePath();
          }
          break;
        case "con":
          if (fs.conFile != null) {
            log.reportError("Identified assumed duplicate " + "CON" + " files: {" + fs.conFile
                            + " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
            System.exit(1);
          } else {
            fs.conFile = srcF.getAbsolutePath();
          }
          break;
        case "int":
          if (fs.intFile != null) {
            log.reportError("Identified assumed duplicate " + "INT" + " files: {" + fs.intFile
                            + " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
            System.exit(1);
          } else {
            fs.intFile = srcF.getAbsolutePath();
          }
          break;
        case "snp":
          if (fs.bimFile != null) {
            log.reportError("Identified assumed duplicate " + "SNP" + " files: {" + fs.bimFile
                            + " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
            System.exit(1);
          } else {
            fs.bimFile = srcF.getAbsolutePath();
          }
          break;
        case "cal":
          if (fs.bedFile != null) {
            log.reportError("Identified assumed duplicate " + "CAL" + " files: {" + fs.bedFile
                            + " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
            System.exit(1);
          } else {
            fs.bedFile = srcF.getAbsolutePath();
          }
          break;
        default:
          log.reportTimeWarning("Parsing NOT YET IMPLEMENTED for the following data file: "
                                + srcF.getAbsolutePath());
          break;
      }
    }

    ArrayList<Integer> missingChrs = new ArrayList<>();
    ArrayList<Integer> incompleteChrs = new ArrayList<>();
    for (int i = 0; i < 27; i++) {
      if (!fileSets.containsKey(i)) {
        missingChrs.add(i);
      } else if (!fileSets.get(i).isComplete()) {
        incompleteChrs.add(i);
        fileSets.remove(i);
      }
    }

    if (!missingChrs.isEmpty()) {
      log.reportTimeWarning("Missing " + missingChrs.size() + " chrs: "
                            + ArrayUtils.toStr(missingChrs, ", "));
    }
    if (!incompleteChrs.isEmpty()) {
      log.reportTimeWarning("Removed chrs with incomplete files: "
                            + ArrayUtils.toStr(incompleteChrs, ", "));
    }
  }

  protected void writeMarkerDetailSet() {
    if (!Files.exists(proj.MARKER_DETAILS_FILENAME.getValue())) {
      try {
        long t1 = System.nanoTime();
        MarkerDetailSet mds = new MarkerDetailSet(loadMarkers());
        mds.serialize(proj.MARKER_DETAILS_FILENAME.getValue());
        markerMap = mds.getMarkerNameMap();
        log.reportTime("Created markerDetailSet in " + ext.getTimeElapsedNanos(t1));
      } catch (IOException e) {
        log.reportException(e);
      }
    } else {
      log.report("Project marker detail set file already exists at "
                 + proj.MARKER_DETAILS_FILENAME.getValue() + "; skipping creation.");
      markerMap = proj.getMarkerSet().getMarkerNameMap();
    }
  }

  protected void parseMarkerData() {
    long time = System.nanoTime();
    int numThreads = Math.min(fileSets.size(), maxChrThreads);
    ExecutorService executor = Executors.newFixedThreadPool(numThreads);

    final MarkerLogger mkrLog = new MarkerLogger(log, logEveryNMins);

    CountDownLatch cdl = new CountDownLatch(fileSets.size()) {

      int f = fileSets.size();

      @Override
      public void countDown() {
        super.countDown();
        log.reportTime("Marker Data Countdown: " + (--f) + " of " + fileSets.size());
      }
    };

    for (Integer v : fileSets.keySet()) {
      final FileSet fs = fileSets.get(v);
      Runnable runn = new Runnable() {

        @Override
        public void run() {
          try {
            fs.buildLookups(log);
            createMDRAF(fs, mkrLog);
            cdl.countDown();
            log.reportTime("Completed parsing chr" + fs.chr);
          } catch (IOException | Elision e) {
            log.reportError("FAILED TO PARSE DATA FILE FOR CHROMOSOME: " + fs.chr);
            e.printStackTrace();
            System.exit(1);
          } catch (Exception e) {
            log.reportError("FAILED TO PARSE DATA FILE FOR CHROMOSOME: " + fs.chr);
            e.printStackTrace();
            System.exit(1);
          }
        }
      };
      executor.submit(runn);
    }
    log.reportTime("Queued MarkerData parsing for " + fileSets.size() + " chrs with " + numThreads
                   + " threads.");

    executor.shutdown();
    try {
      cdl.await(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {
      /**/
    }
    mkrLog.stopLogging();

    log.reportTime("Completed parsing MarkerData in " + ext.getTimeElapsedNanos(time));
  }

  @Override
  protected void createMarkerPositions() {
    List<String> markerList = new ArrayList<>();
    String file = proj.MARKER_POSITION_FILENAME.getValue();
    if (!Files.exists(file)) {
      long time = System.nanoTime();
      PrintWriter writer = Files.getAppropriateWriter(file);
      writer.println("Marker\tChr\tPosition");

      for (int i = 0; i < 27; i++) {
        if (!fileSets.containsKey(i)) {
          continue;
        }
        FileSet fs = fileSets.get(i);
        fs.loadBIMData();
        numberOfMarkersInProj.addAndGet(fs.bimData.length);

        for (String[] parts : fs.bimData) {
          writer.println(parts[1] + "\t" + parts[0] + "\t" + parts[2]);
          markerList.add(parts[1]);
        }
      }
      writer.flush();
      writer.close();
      allMarkers = markerList.toArray(new String[markerList.size()]);

      log.reportTime("Completed markerPositions: " + proj.MARKER_POSITION_FILENAME.getValue()
                     + " in " + ext.getTimeElapsedNanos(time));
    } else {
      allMarkers = HashVec.loadFileToStringArray(file, true, new int[] {0}, false);
      numberOfMarkersInProj.set(getNumMarkers());

      for (int i = 0; i < 27; i++) {
        if (!fileSets.containsKey(i)) {
          continue;
        }
        FileSet fs = fileSets.get(i);
        fs.loadBIMData();
      }

      log.reportTime("Project markerPositions file already exists at " + file
                     + "; skipping creation (first marker found: " + allMarkers[0] + ").");
    }
  }

  private List<Marker> loadMarkers() throws IOException {
    String[] headerFactors = {"Affy SNP ID", "dbSNP RS ID", "Chromosome", "Physical Position",
                              "Allele A", "Allele B", "Ref Allele", "Flank"};
    int[] hdrInds = null;
    HashMap<String, Integer> markerIndices = new HashMap<>();
    ArrayList<Integer> allInds = new ArrayList<>(getNumMarkers());
    for (int i = 0; i < getNumMarkers(); i++) {
      markerIndices.put(allMarkers[i], i);
      allInds.add(i);
    }
    final Marker[] markers = new Marker[getNumMarkers()];
    String line = null;
    boolean foundHeader = false;
    String[] parts;
    String mkrName;
    int inAnnotNotProj = 0;
    int ind;
    int lines = 0;
    BufferedReader reader = Files.getAppropriateReader(annotFile);
    while ((line = reader.readLine()) != null) {
      lines++;
      if (line.charAt(0) == '#') continue;
      if (!foundHeader) {
        foundHeader = true; // first nonComment line
        hdrInds = ext.indexFactors(headerFactors, ext.splitCommasIntelligently(line, true, log),
                                   false);

        continue;
      }
      parts = ext.splitCommasIntelligently(line, true, log);
      mkrName = isMissing(parts[hdrInds[1]]) ? parts[hdrInds[0]] : parts[hdrInds[1]];

      if (!markerIndices.containsKey(mkrName)) {
        if (markerIndices.containsKey(parts[hdrInds[0]])) {
          mkrName = parts[hdrInds[0]];
          ind = markerIndices.get(mkrName);
        } else {
          inAnnotNotProj++;
          ind = -1;
        }
      } else {
        ind = markerIndices.get(mkrName);
      }
      GenomicPosition gPos;
      if (ind != -1) {
        allInds.remove(Integer.valueOf(ind));

        gPos = new GenomicPosition(!isMissing(parts[hdrInds[2]]) ? Positions.chromosomeNumber(parts[hdrInds[2]],
                                                                                              log)
                                                                 : (byte) 0,
                                   !isMissing(parts[hdrInds[3]]) ? Integer.parseInt(parts[hdrInds[3]])
                                                                 : 0);

        String a = parts[hdrInds[4]];
        String b = parts[hdrInds[5]];
        boolean aMiss = isMissing(a);
        boolean bMiss = isMissing(b);
        boolean aRef = parts[hdrInds[4]].equals(parts[hdrInds[6]]);

        if (aMiss && bMiss) {
          // skip
          a = "N";
          b = "N";
        } else if (aMiss) {
          a = Character.toString(parts[hdrInds[7]].charAt(parts[hdrInds[7]].indexOf('[') - 1));
          b = a + b;
        } else if (bMiss) {
          throw new IllegalStateException("B allele is an unexpected indel for marker " + mkrName);
        }

        markers[markerIndices.get(mkrName)] = new Marker(mkrName, gPos,
                                                         AllelePair.of(Allele.create(a, aRef),
                                                                       Allele.create(b, !aRef)));
      }
    }
    if (lines == 0 && annotFile.endsWith(".zip")) {
      throw new RuntimeException("Coudn't read the marker annotation file at " + annotFile
                                 + ".  The zipped annotation file contains two files: the annotation file itself and a README file.  Please unzip these files, set the annotation file argument to point to the unzipped annotation file, and rerun the pipeline.");
    }
    if (inAnnotNotProj > 0) {
      log.reportTimeWarning("Found " + inAnnotNotProj
                            + " markers in the annotation file that were not present in the project data.");
    }
    if (allInds.size() > 0) {
      log.reportError(allInds.size()
                      + " markers were not found in the annotation file! Cannot determine A/B alleles for these markers.  Genotypes may be flipped for these markers.  A file containing these markers and the source BIM data will be written to "
                      + proj.PROJECT_DIRECTORY.getValue() + "missingAB.snp");
      PrintWriter missing = Files.getAppropriateWriter(proj.PROJECT_DIRECTORY.getValue()
                                                       + "missingAB.snp");
      missing.println("SNP\tCHR\tPOS\tREF\tALT");
      for (int i : allInds) {
        String snp = allMarkers[i];
        for (FileSet fs : fileSets.values()) {
          if (fs.bimMap.containsKey(snp)) {
            String[] bimD = fs.bimData[fs.bimMap.get(snp)];
            markers[i] = new Marker(snp,
                                    new GenomicPosition((byte) fs.chr, Integer.parseInt(bimD[2])),
                                    AllelePair.of('N', 'N'));
            missing.println(bimD[1] + "\t" + bimD[0] + "\t" + bimD[2] + "\t" + bimD[3] + "\t"
                            + bimD[4]);
            break;
          }
        }
      }
      missing.close();
    }
    reader.close();
    markerIndices = null;
    return Arrays.asList(markers);
  }

  private boolean isMissing(String value) {
    return value.equals("---") || value.equals("\"---\"") || value.equals("-")
           || value.equals("\"-\"");
  }

  private void createMDRAF(FileSet fs, MarkerLogger mkrLog) throws IOException, Elision {
    final int nInd = getNumSamples();

    final byte nullStatus = getNullStatus();
    int numThreads = maxMDRAFThreadsPerChr;

    fs.loadBIMData(); // no-op if already loaded
    String[] mkrNames = Matrix.extractColumn(fs.bimData, 1);
    fs.mkrBins = ArrayUtils.splitUpArray(mkrNames, (mkrNames.length / numMarkersPerFile) + 1, log);
    mkrNames = null;

    String[] existing = readWritten();
    Map<Integer, String[]> toAdd = new TreeMap<>();
    int start = 0;
    int end = 0;
    for (int i = 0; i < fs.mkrBins.size(); i++) {
      end = start + fs.mkrBins.get(i).length;
      String mdRAFName = getMDRAFName(fs.chr, start, end);
      if (ext.indexOfStr(mdRAFName, existing) == -1) {
        toAdd.put(start, fs.mkrBins.get(i));
      }
      start = end;
    }

    Comparator<Entry<Integer, String[]>> entComp = new Comparator<Entry<Integer, String[]>>() {

      @Override
      public int compare(Entry<Integer, String[]> o1, Entry<Integer, String[]> o2) {
        return o1.getKey().compareTo(o2.getKey());
      }
    };
    List<List<Entry<Integer, String[]>>> binsOfBins = new ArrayList<>();
    int numPerBin = (toAdd.size() / numThreads) + 1;
    ArrayList<Entry<Integer, String[]>> binOfBins = new ArrayList<>();
    binsOfBins.add(binOfBins);
    int ind = 0;
    for (Entry<Integer, String[]> bin : toAdd.entrySet()) {
      binOfBins.add(bin);
      ind++;
      if (ind > numPerBin) {
        ind = 0;
        binOfBins = new ArrayList<>();
        binsOfBins.add(binOfBins);
      }
    }
    for (int i = binsOfBins.size() - 1; i >= 0; i--) {
      if (binsOfBins.get(i).isEmpty()) {
        binsOfBins.remove(i);
      } else {
        binsOfBins.get(i).sort(entComp);
      }
    }

    ArrayList<Runnable> runners = new ArrayList<>();
    int sz = binsOfBins.size();
    final CountDownLatch cdl = new CountDownLatch(binsOfBins.size()) {

      @Override
      public void countDown() {
        super.countDown();
        log.reportTime("COUNTED DOWN -- " + getCount() + " remaining of " + sz);
      }
    };

    for (List<Entry<Integer, String[]>> ent : binsOfBins) {
      runners.add(new Runnable() {

        @Override
        public void run() {

          String mdRAF;
          try {
            FileReaderSet frs = new FileReaderSet();
            frs.open(fs);

            byte[] magicBytes = new byte[3];
            frs.bedIn.read(magicBytes);
            if (magicBytes[2] == 0) {
              log.reportError("Error - .bed file is sample-dominant: " + fs.bedFile);
              System.exit(1);
            }

            int prevBatchEnd = 0;
            for (int i = 0; i < ent.size(); i++) {
              long t1 = System.nanoTime();
              Entry<Integer, String[]> ent2 = ent.get(i);
              mdRAF = write(fs, frs, prevBatchEnd, ent2.getKey(), ent2.getValue(), nullStatus, nInd,
                            mkrLog);
              success(mdRAF, t1);
              prevBatchEnd = ent2.getKey() + ent2.getValue().length;
            }

            frs.close();
            frs = null;
          } catch (Elision e) {
            e.printStackTrace();
          } catch (IOException e) {
            e.printStackTrace();
          } catch (Exception e) {
            e.printStackTrace();
          }
          cdl.countDown();
        }
      });
    }

    if (runners.size() > 0) {
      numThreads = Math.min(runners.size(), numThreads);
      ExecutorService executor = Executors.newFixedThreadPool(numThreads);
      for (int i = 0; i < runners.size(); i++) {
        executor.submit(runners.get(i));
        System.out.println("Queued bin " + (i + 1) + " of " + runners.size() + " (with "
                           + binsOfBins.get(i).size() + " bins) (found "
                           + (fs.mkrBins.size() - toAdd.size()) + " already parsed) for chr"
                           + fs.chr + " with " + numThreads + " threads.");
      }
      executor.shutdown();
      try {
        cdl.await(Long.MAX_VALUE, TimeUnit.DAYS);
      } catch (InterruptedException e) {
        /**/
      }
    }
  }

  private synchronized String[] readWritten() {
    String temp = proj.PROJECT_DIRECTORY.getValue() + "markersWritten.txt";
    return Files.exists(temp) ? HashVec.loadFileToStringArray(temp, false, null, false)
                              : new String[0];
  }

  private synchronized void success(String mkrFile, long t1) {
    String temp = proj.PROJECT_DIRECTORY.getValue() + "markersWritten.txt";
    Files.appendStringToFile(temp, mkrFile);
    log.reportTime("Wrote marker data file " + mkrFile + " in " + ext.getTimeElapsedNanos(t1));
  }

  private String getMDRAFName(int chr, int start, int end) {
    return "markers." + chr + "." + start + "." + end + MarkerData.MARKER_DATA_FILE_EXTENSION;
  }

  private String write(FileSet fs, FileReaderSet frs, int prevBatchEnd, int startBatchInd,
                       String[] mkrNames, byte nullStatus, int nInd,
                       MarkerLogger markerLogger) throws Elision, IOException {

    int lineSkips = startBatchInd - prevBatchEnd;
    if (lineSkips > 0) {
      long t1 = System.currentTimeMillis();
      frs.readAhead(fs, prevBatchEnd, lineSkips);
      log.reportTime("Chr" + fs.chr + " - skipped " + lineSkips + " lines from " + prevBatchEnd
                     + " to " + startBatchInd + " in " + ext.getTimeElapsed(t1));
    }

    RandomAccessFile mdRAF;

    Hashtable<String, Float> outOfRangeTable = new Hashtable<>();
    String mdRAFName = getMDRAFName(fs.chr, startBatchInd, (startBatchInd + mkrNames.length));

    mdRAF = openMDRAF(mdRAFName, getNumSamples(), nullStatus, mkrNames);

    int bedBlockSize = (int) Math.ceil(nInd / 4.0);
    int binBlockSize = nInd * 8;

    int bytesPerSamp = Sample.getNBytesPerSampleMarker(getNullStatus()); // 12
    int markerBlockSize = nInd * bytesPerSamp;
    byte[] mkrBuff;
    String line;
    String[] parts;
    int flipped = 0;
    // double log2 = Math.log(2);
    log.reportTime("Chr" + fs.chr + "; start=" + startBatchInd + "; reading markers...");
    for (int i = startBatchInd; i < startBatchInd + mkrNames.length; i++) {
      markerLogger.logMarker(fs.chr, startBatchInd, mkrNames.length, i);
      boolean flipGenos = false;
      if (markerMap.containsKey(mkrNames[i - startBatchInd])) {
        flipGenos = markerMap.get(mkrNames[i - startBatchInd]).getAlleleA().isNonReference();
        if (flipGenos) flipped++;
      } else {
        throw new IllegalStateException("Marker name " + mkrNames[i - startBatchInd]
                                        + " wasn't present in the MarkerDetailSet.");
      }

      mkrBuff = new byte[markerBlockSize];
      int buffInd = 0;
      int sampInd = 0;

      line = frs.conIn.readLine();
      parts = line.split(" ");
      line = null;
      if (parts.length != nInd) {
        log.reportError("Mismatched # of confidence scores {fnd: " + parts.length + ", exp: " + nInd
                        + "}.  File: " + fs.conFile);
        System.exit(1);
      }
      float gcV;
      for (String con : parts) {
        gcV = 1 - Float.parseFloat(con);
        if (gcV == 2) { // i.e. 1 - (-1)
          gcV = Float.NaN;
        }
        Compression.gcBafCompress(gcV, mkrBuff, sampInd);
        sampInd += bytesPerSamp;
      }
      buffInd += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
      sampInd = buffInd;
      parts = null;

      byte[] intensBytes = new byte[binBlockSize];
      frs.binIn.seek(((long) i) * ((long) binBlockSize));
      frs.binIn.read(intensBytes);

      float a, b, x, y;
      boolean oor;
      byte[] intA = new byte[4];
      byte[] intB = new byte[4];

      for (int bitInd = 0, binInd = 0; bitInd < nInd; bitInd++) {
        intA[0] = intensBytes[binInd++];
        intA[1] = intensBytes[binInd++];
        intA[2] = intensBytes[binInd++];
        intA[3] = intensBytes[binInd++];

        intB[0] = intensBytes[binInd++];
        intB[1] = intensBytes[binInd++];
        intB[2] = intensBytes[binInd++];
        intB[3] = intensBytes[binInd++];

        a = BGENBitMath.bytesToFloat(true, intA);
        b = BGENBitMath.bytesToFloat(true, intB);
        x = (float) (a / scaleFactor);
        y = (float) (b / scaleFactor);
        // x = (float) Maths.log2(a / b);
        // y = (float) (Maths.log2(a * b) / 2);

        oor = !Compression.xyCompressPositiveOnly(a == -1 ? Float.NaN : (float) (x), mkrBuff,
                                                  sampInd);
        if (oor) {
          outOfRangeTable.put((i - startBatchInd) + "\t" + bitInd + "\tx", x);
        }
        oor = !Compression.xyCompressPositiveOnly(b == -1 ? Float.NaN : (float) (y), mkrBuff,
                                                  sampInd + Compression.REDUCED_PRECISION_XY_NUM_BYTES);
        if (oor) {
          outOfRangeTable.put((i - startBatchInd) + "\t" + bitInd + "\ty", y);
        }
        sampInd += bytesPerSamp;
      }
      buffInd += Compression.REDUCED_PRECISION_XY_NUM_BYTES; // X
      buffInd += Compression.REDUCED_PRECISION_XY_NUM_BYTES; // Y
      sampInd = buffInd;
      intensBytes = null;

      line = frs.bafIn.readLine();
      parts = line.split(" ");
      line = null;
      if (parts.length != nInd) {
        log.reportError("Mismatched # of BAF values {fnd: " + parts.length + ", exp: " + nInd
                        + "}.  File: " + fs.bafFile);
        System.exit(1);
      }
      float baf = Float.NaN;
      for (String bafStr : parts) {
        try {
          baf = Float.parseFloat(bafStr);
        } catch (NumberFormatException e) {
          if (ext.isMissingValue(bafStr)) {
            baf = Float.NaN;
          } else {
            log.reportError("Malformed BAF value: " + bafStr + " for chr " + fs.chr);
            System.exit(1);
          }
        }
        Compression.gcBafCompress(baf, mkrBuff, sampInd);
        sampInd += bytesPerSamp;
      }
      buffInd += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
      sampInd = buffInd;
      parts = null;

      line = frs.lrrIn.readLine();
      parts = line.split(" ");
      line = null;
      if (parts.length != nInd) {
        log.reportError("Mismatched # of L2R values {fnd: " + parts.length + ", exp: " + nInd
                        + "}.  File: " + fs.lrrFile);
        System.exit(1);
      }
      float lrr = Float.NaN;
      for (int samp = 0; samp < parts.length; samp++) {
        try {
          lrr = Float.parseFloat(parts[samp]);
        } catch (NumberFormatException e) {
          if (ext.isMissingValue(parts[samp])) {
            lrr = Float.NaN;
          } else {
            log.reportError("Malformed LRR value: " + parts[samp] + " for chr " + fs.chr);
            System.exit(1);
          }
        }
        oor = -1 == Compression.lrrCompress(lrr, mkrBuff, sampInd);
        sampInd += bytesPerSamp;
        if (oor) {
          outOfRangeTable.put((i - startBatchInd) + "\t" + samp + "\tlrr", lrr);
        }
      }
      buffInd += Compression.REDUCED_PRECISION_LRR_NUM_BYTES;
      sampInd = buffInd;
      parts = null;

      byte[] markerBytes = new byte[bedBlockSize];
      frs.bedIn.seek(3L + ((long) i) * ((long) bedBlockSize)); // ALWAYS seek forward
      frs.bedIn.read(markerBytes);
      int tempInd = sampInd;
      for (int bitInd = 0; bitInd < markerBytes.length; bitInd++) {
        byte bedByte = markerBytes[bitInd];
        byte[] genotypes = DosageData.decodeBedByte(bedByte);
        for (int g = 0; g < genotypes.length; g++) {
          int idInd = bitInd * 4 + g;
          if (idInd == -1 || idInd >= nInd) {
            // tries to load three more samples than exist
            continue;
          }
          if (flipGenos && genotypes[g] != -1) {
            genotypes[g] = (byte) (2 - genotypes[g]);
          }
          mkrBuff[tempInd] = Compression.genotypeCompress(importGenoAsAB ? genotypes[g] : (byte) -1,
                                                          importGenoAsAB ? (byte) 0 : genotypes[g]);
          tempInd += bytesPerSamp;
        }
        genotypes = null;
      }
      buffInd += Compression.REDUCED_PRECISION_ABFORWARD_GENOTYPE_NUM_BYTES;
      sampInd = buffInd;
      markerBytes = null;

      mdRAF.write(mkrBuff);
      mkrBuff = null;
    }

    byte[] oorBytes = Compression.objToBytes(outOfRangeTable);
    mdRAF.write(Compression.intToBytes(oorBytes.length));
    mdRAF.write(oorBytes);
    log.reportTime("Wrote " + outOfRangeTable.size() + " outliers to " + mdRAFName);
    log.reportTime(flipped + " markers out of " + mkrNames.length + " had flipped genotypes.");

    mdRAF.close();
    // System.gc();

    oorBytes = null;
    outOfRangeTable = null;
    return mdRAFName;
  }

  private static void addDirs(File src, List<File> subDirs) {
    for (File f : src.listFiles(File::isDirectory)) {
      if (!subDirs.contains(f)) {
        subDirs.add(f);
        addDirs(f, subDirs);
      }
    }
  }

  private static List<File> discoverAllSourceFiles(String dir, Logger log) {
    File srcDir = new File(dir);
    ArrayList<File> subDirs = new ArrayList<>();
    subDirs.add(srcDir);
    addDirs(srcDir, subDirs);

    log.reportTime("Found " + subDirs.size() + " source directories to scan.");

    /**
     * Delim = "_"; <br />
     * At least 3 tokens; <br />
     * First token must be "ukb"; <br />
     * Second token must be from FileSet.EXP; <br />
     * Third token must start with "chr"; <br />
     * Will not accept any .fam or .gzi file; <br />
     */
    FilenameFilter ff = new FilenameFilter() {

      @Override
      public boolean accept(File dir, String name) {
        String[] parts = name.split("_");
        if (parts.length <= 2 || parts.length > 4) {
          return false;
        }
        if (name.endsWith(".fam")) {
          return false;
        }
        if (name.endsWith(".gzi")) {
          return false;
        }
        if (!parts[0].equals("ukb")) {
          return false;
        }
        if (ext.indexOfStr(parts[1], FileSet.EXP) == -1) {
          return false;
        }
        if (!parts[2].startsWith("chr")) {
          return false;
        }
        return true;
      }
    };

    ArrayList<File> allSrcFiles = new ArrayList<>();
    for (File f : subDirs) {
      for (File s : f.listFiles(ff)) {
        allSrcFiles.add(s);
      }
    }

    log.reportTime("Found " + allSrcFiles.size() + " potential source data files.");
    return allSrcFiles;
  }

  private class MarkerLogger {

    public MarkerLogger(Logger log, int minutes) {
      mkrMap = new ConcurrentHashMap<>();
      keyMap = new ConcurrentHashMap<>();
      this.log = log;
      timer = new Timer();
      startLogging(minutes);
    }

    Logger log;
    ConcurrentMap<String, Integer> mkrMap;
    ConcurrentMap<String, int[]> keyMap;

    public void logMarker(int chr, int bStrt, int bLen, int currMkr) {
      String keyStr = chr + "\t" + bStrt + "\t" + bLen;
      synchronized (keyStr) {
        Integer exists = mkrMap.put(keyStr, currMkr);
        if (exists == null) {
          keyMap.put(keyStr, new int[] {chr, bStrt, bLen});
        } else if (exists > currMkr) {
          mkrMap.put(keyStr, exists);
        }
      }
    }

    private void startLogging(int minutes) {
      // wait two minutes, then every N minutes
      timer.schedule(logging, 1000 * 60 * 2, 1000 * 60 * minutes);
    }

    public void stopLogging() {
      timer.cancel();
    }

    String prevReport = null;
    int sameReportStreak = 0;
    int streakReportAt = 3;
    Timer timer;
    TimerTask logging = new TimerTask() {

      @Override
      public void run() {
        if (keyMap.size() == 0) return;
        long t1 = System.nanoTime();
        List<String> sortedKeys = keyMap.keySet().stream().sorted((o1, o2) -> {
          String[] pts1 = o1.split("\t");
          String[] pts2 = o2.split("\t");
          for (int i = 0; i < pts1.length; i++) {
            int comp1 = Integer.parseInt(pts1[i]);
            int comp2 = Integer.parseInt(pts2[i]);
            int comp = Integer.compare(comp1, comp2);
            if (comp != 0) {
              return comp;
            }
          }
          return 0;
        }).collect(Collectors.toList());

        int max = 0;
        HashMap<Integer, ArrayList<String>> chrKeys = new HashMap<>();
        ArrayList<Integer> keys = new ArrayList<>();
        for (String k : sortedKeys) {
          int[] chrStrtLen = keyMap.get(k);
          int chr = chrStrtLen[0];
          double pct = (1 + (mkrMap.get(k) - chrStrtLen[1])) / (double) chrStrtLen[2];
          int v = (int) (100 * pct);
          if (!chrKeys.containsKey(chr)) {
            chrKeys.put(chr, new ArrayList<>());
            keys.add(chr);
          }
          chrKeys.get(chr).add(pct < 0.998 ? v + "%" : "C");
        }
        for (ArrayList<String> cK : chrKeys.values()) {
          if (cK.size() > max) {
            max = cK.size();
          }
        }
        StringBuilder header = new StringBuilder("CHR");
        for (Integer c : keys) {
          header.append("\t").append(c);
        }

        StringBuilder fullReport = new StringBuilder();
        fullReport.append(ext.getTime()).append(" - ").append(header.toString()).append("\n");
        for (int i = 0; i < max; i++) {
          StringBuilder line = new StringBuilder("BIN_").append(i);
          for (Integer c : keys) {
            line.append("\t");
            if (chrKeys.get(c).size() <= i) {
              line.append("-");
            } else {
              line.append(chrKeys.get(c).get(i));
            }
          }
          fullReport.append(ext.getTime()).append(" - ").append(line.toString()).append("\n");
        }

        fullReport.append("(wrote to log in " + ext.getTimeElapsedNanos(t1) + ")\n\n");
        if (fullReport.toString().equals(prevReport)) {
          sameReportStreak++;
        } else {
          sameReportStreak = 0;
        }
        if (sameReportStreak > 0 && sameReportStreak % streakReportAt == 0) {
          log.reportTimeWarning("Warning - log message has been the same " + sameReportStreak
                                + " times in a row");
          log.reportTimeWarning(ext.reportMemoryUsage());
          System.gc();
        }
        log.report(prevReport = fullReport.toString());
      }
    };

  }

  private interface ReaderProxy {

    public String readLine() throws IOException;

    public void close() throws IOException;

  }

  private class BGZipProxy implements ReaderProxy {

    BlockCompressedInputStream bcis;

    public BGZipProxy(BlockCompressedInputStream is) {
      this.bcis = is;
    }

    @Override
    public String readLine() throws IOException {
      return bcis.readLine();
    }

    @Override
    public void close() throws IOException {
      bcis.close();
    }

  }

  private class TextReaderProxy implements ReaderProxy {

    BufferedReader br;

    public TextReaderProxy(BufferedReader br1) {
      this.br = br1;
    }

    @Override
    public String readLine() throws IOException {
      return br.readLine();
    }

    @Override
    public void close() throws IOException {
      br.close();
    }

  }

  private class FileReaderSet {

    RandomAccessFile bedIn;
    RandomAccessFile binIn;
    ReaderProxy conIn;
    ReaderProxy lrrIn;
    ReaderProxy bafIn;

    public void open(FileSet fs) throws IOException {
      bedIn = new RandomAccessFile(fs.bedFile, "r");
      binIn = new RandomAccessFile(fs.intFile, "r");
      conIn = open(fs.conFile, 0);
      lrrIn = open(fs.lrrFile, 0);
      bafIn = open(fs.bafFile, 0);
    }

    public void readAhead(FileSet fs, int start, int binSize) throws IOException {
      int line = start + binSize;
      conIn = open(fs.conFile, fs.conInds.get(line));
      lrrIn = open(fs.lrrFile, fs.lrrInds.get(line));
      bafIn = open(fs.bafFile, fs.bafInds.get(line));
    }

    private ReaderProxy open(String file, long skip) throws IOException {
      InputStream is;
      if (file.endsWith(".gz")) {
        is = new BlockCompressedInputStream(new File(file));
        if (skip > 0) {
          ((BlockCompressedInputStream) is).seek(skip);
        }
        return new BGZipProxy((BlockCompressedInputStream) is);
        // return new BGZipReaderProxy(new BGZipReader(file));
      } else {
        is = new FileInputStream(file);
        if (skip > 0) {
          is.skip(skip);
        }
        return new TextReaderProxy(new BufferedReader(new InputStreamReader(is), 4000000));
      }
    }

    public void close() throws IOException {
      bedIn.close();
      binIn.close();
      conIn.close();
      lrrIn.close();
      bafIn.close();

      bedIn = null;
      binIn = null;
      conIn = null;
      lrrIn = null;
      bafIn = null;
    }
  }

  static class FileSet {

    int chr;
    String[][] bimData;
    HashMap<String, Integer> bimMap;

    String bimFile;
    String bedFile;
    String famFile;

    String bafFile;
    String lrrFile;
    String intFile;
    String conFile;

    /** genotype info: */
    String bgiFile;
    String bgenFile;

    HashMap<Integer, Long> conInds;
    HashMap<Integer, Long> lrrInds;
    HashMap<Integer, Long> bafInds;

    List<String[]> mkrBins;

    Logger log;

    /**
     * Expected files
     */
    final static String[] EXP = {"l2r", "baf", "con", "int", "snp", "cal",};

    public void buildLookups(Logger log) throws IOException {
      long t1 = System.nanoTime();
      ArrayList<Runnable> runners = new ArrayList<>();
      CountDownLatch cdl = new CountDownLatch(3);
      runners.add(new Runnable() {

        @Override
        public void run() {
          try {
            conInds = buildLookup(conFile, log);
          } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
          }
          cdl.countDown();
        }
      });
      runners.add(new Runnable() {

        @Override
        public void run() {
          try {
            lrrInds = buildLookup(lrrFile, log);
          } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
          }
          cdl.countDown();
        }
      });
      runners.add(new Runnable() {

        @Override
        public void run() {
          try {
            bafInds = buildLookup(bafFile, log);
          } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
          }
          cdl.countDown();
        }
      });

      ExecutorService executor = Executors.newFixedThreadPool(Math.min(runners.size(), 3));
      for (Runnable r : runners) {
        executor.submit(r);
      }
      executor.shutdown();
      try {
        cdl.await(Long.MAX_VALUE, TimeUnit.DAYS);
      } catch (InterruptedException e) {
        /**/
      }
      log.reportTime("Built chr" + chr + " lookups in " + ext.getTimeElapsedNanos(t1));
    }

    @SuppressWarnings("unchecked")
    private HashMap<Integer, Long> buildLookup(String file, Logger log) throws IOException {
      long t1 = System.nanoTime();
      String lookupFile = ext.rootOf(file, false) + "_lookup.dat";
      if (Files.exists(lookupFile)) {
        HashMap<Integer, Long> indMap = (HashMap<Integer, Long>) SerializedFiles.readSerial(lookupFile);
        return indMap;
      }
      log.reportTime("Building lookup file for " + ext.rootOf(lookupFile));
      HashMap<Integer, Long> lineIndices = new HashMap<>();

      int line = 0;
      if (file.endsWith(".gz")) {
        BlockCompressedInputStream bcis = new BlockCompressedInputStream(new File(file));

        while (bcis.readLine() != null) {
          line++;
          lineIndices.put(line, bcis.getFilePointer());
        }

        bcis.close();
      } else {
        FileInputStream fis = new FileInputStream(file);
        InputStreamReader isr = new InputStreamReader(fis);
        int bufferLength = 16384; // 16 Kb buffer
        ByteBuffer buffer = ByteBuffer.wrap(new byte[bufferLength]);
        FileChannel fc = fis.getChannel();

        int read = 0;
        while ((read = fc.read(buffer)) != -1) {
          byte[] array = buffer.array();
          if (read < bufferLength) {
            array = new byte[read];
            System.arraycopy(buffer.array(), 0, array, 0, read);
          }
          String bufferString = new String(array);
          for (int i = 0; i < bufferString.length(); i++) {
            if (bufferString.charAt(i) == '\n') {
              line++;
              long pos = fc.position() + 1;
              if (i < bufferString.length()) {
                pos -= bufferString.substring(i).getBytes().length;
              }
              lineIndices.put(line, pos);
            }
          }
          buffer.clear();
        }

        isr.close();
      }

      SerializedFiles.writeSerial(lineIndices, lookupFile);
      log.reportTime("Built {" + ext.rootOf(lookupFile) + "} in " + ext.getTimeElapsedNanos(t1));
      return lineIndices;
    }

    boolean isComplete() {
      boolean a = Files.exists(bimFile);
      boolean b = Files.exists(bedFile);
      boolean c = Files.exists(famFile);
      boolean d = Files.exists(bafFile);
      boolean e = Files.exists(lrrFile);
      boolean f = Files.exists(intFile);
      boolean g = Files.exists(conFile);
      if (!a) {
        log.reportTimeWarning("Missing BIM file for chr " + chr);
      }
      if (!b) {
        log.reportTimeWarning("Missing BED file for chr " + chr);
      }
      if (!c) {
        log.reportTimeWarning("Missing FAM file for chr " + chr);
      }
      if (!d) {
        log.reportTimeWarning("Missing BAF file for chr " + chr);
      }
      if (!e) {
        log.reportTimeWarning("Missing LRR file for chr " + chr);
      }
      if (!f) {
        log.reportTimeWarning("Missing INT file for chr " + chr);
      }
      if (!g) {
        log.reportTimeWarning("Missing CON file for chr " + chr);
      }
      return a && b && c && d && e && f && g;
    }

    void loadBIMData() {
      if (bimData != null) {
        return;
      }
      long t1 = System.nanoTime();
      bimData = HashVec.loadFileToStringMatrix(bimFile, false, new int[] {0, 1, 3, 4, 5});
      bimMap = new HashMap<>();
      for (int i = 0; i < bimData.length; i++) {
        bimMap.put(bimData[i][1], i);
      }
      log.reportTime("Loaded chr" + chr + " .bim data in " + ext.getTimeElapsedNanos(t1));
    }

  }

  private static final String ARG_SRC_DIR = "source";
  private static final String ARG_PROJ_DIR = "projDir";
  private static final String ARG_PROP_DIR = "propDir";
  private static final String ARG_PROJ_NAME = "projName";
  private static final String ARG_FAM_FILE = "fam";
  private static final String ARG_ANNOT_CSV = "annot";

  private static final String DESC_SRC_DIR = "Directory of UK BioBank source files";
  private static final String DESC_PROJ_DIR = "Directory in which to create parsed files";
  private static final String DESC_PROP_DIR = "Directory in which to create project properties file";
  private static final String DESC_PROJ_NAME = "Project name";
  private static final String DESC_FAM_FILE = "UKBioBank-provided, PLINK-formatted .fam file";
  private static final String DESC_ANNOT_CSV = "Affymetrix annotation .csv file";

  public static void main(String[] args) {
    String sourceDir = "";
    String projDir = "";
    String propFileDir = "";
    String projName = "";
    String famFile = "";
    String annotFile = "";
    int chrThreads;
    int threadsPerChr;
    int mkrsPerFile;
    float scale;
    String jobID;

    CLI cli = new CLI(UKBBParsingPipeline.class);

    cli.addArg(ARG_SRC_DIR, DESC_SRC_DIR, true);
    cli.addArg(ARG_PROJ_DIR, DESC_PROJ_DIR, true);
    cli.addArg(ARG_PROP_DIR, DESC_PROP_DIR, true);
    cli.addArgWithDefault(ARG_PROJ_NAME, DESC_PROJ_NAME, "UKBioBank");
    cli.addArg(ARG_FAM_FILE, DESC_FAM_FILE, true);
    cli.addArg(ARG_ANNOT_CSV, DESC_ANNOT_CSV, true);
    cli.addArg("chrThreads", "Number of chromosome threads", false);
    cli.addArg("threadsPerChr", "Number of threads PER chromosome", false);
    cli.addArg("markersPerFile", "Number of markers per mdRAF file", false);
    cli.addArg("scaleFactor", "Scaling factor for X/Y", false);
    cli.addArg("jobID", "Job Identifier");

    cli.parseWithExit(args);

    sourceDir = cli.get(ARG_SRC_DIR);
    projDir = cli.get(ARG_PROJ_DIR);
    propFileDir = cli.get(ARG_PROP_DIR);
    projName = cli.get(ARG_PROJ_NAME);
    famFile = cli.get(ARG_FAM_FILE);
    annotFile = cli.get(ARG_ANNOT_CSV);
    chrThreads = cli.has("chrThreads") ? cli.getI("chrThreads") : -1;
    threadsPerChr = cli.has("threadsPerChr") ? cli.getI("threadsPerChr") : -1;
    mkrsPerFile = cli.has("mkrsPerFile") ? cli.getI("mkrsPerFile") : -1;
    scale = (float) (cli.has("scaleFactor") ? cli.getD("scaleFactor") : -1);
    jobID = cli.has("jobID") ? cli.get("jobID") : null;

    UKBBParsingPipeline parser = new UKBBParsingPipeline();
    parser.setSourceDir(sourceDir);
    parser.setProjectDir(projDir);
    parser.setProjectPropertiesDir(propFileDir);
    parser.setProjectName(projName);
    parser.setFamFile(famFile);
    parser.setJobID(jobID);
    parser.setAnnotationCSV(annotFile);
    if (chrThreads != -1) {
      parser.setChromosomeThreadCount(chrThreads);
    }
    if (threadsPerChr != -1) {
      parser.setThreadsPerChromosomeCount(threadsPerChr);
    }
    if (mkrsPerFile != -1) {
      parser.setMarkersPerMDRAF(mkrsPerFile);
    }
    if (scale > 0) {
      parser.setScaleFactor(scale);
    }
    parser.runPipeline();
  }

}