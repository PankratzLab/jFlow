package org.genvisis.seq.manage.mosdepth;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.genvisis.cnv.filesys.AllelePair;
import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.seq.GenomeBuild;
import org.genvisis.seq.manage.AnnotatedBEDFeature;
import org.genvisis.seq.manage.BEDFileReader;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomicPosition;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.filesys.Segment;
import org.pankratzlab.common.stats.Maths;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.cache.RemovalListener;
import com.google.common.cache.RemovalNotification;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class MosdepthPipeline {

  private static final int MAX_FILES_OPEN = 5000;

  Logger log;

  String projDir;
  String propFileDir;
  String projName;
  double scaleFactor = 10;
  int numMarkersPerFile = 5000;
  int numThreads = Runtime.getRuntime().availableProcessors() / 2;

  public void setProjectName(String projName2) {
    projName = projName2;
  }

  public void setProjectPropertiesDir(String propFileDir2) {
    propFileDir = propFileDir2;
  }

  public void setProjectDir(String projDir2) {
    projDir = projDir2;
  }

  public void setNumThreads(int threads) {
    this.numThreads = threads;
  }

  Project proj;

  String useBed;
  String markerVCF;
  String genoVCF;

  public void setBinsToUseBED(String bedFile) {
    this.useBed = bedFile;
  }

  public void setSelectedMarkerVCF(String vcfFile) {
    this.markerVCF = vcfFile;
  }

  public void setGenotypeVCF(String vcfFile) {
    this.genoVCF = vcfFile;
  }

  public void setMosdepthDirectory(String dir, String ext) {
    this.mosdepthFiles = Files.list(dir, null, ext, false, true);
  }

  /**
   * Must be called after {@link MosdepthPipeline#setMosdepthDirectory(String, String)}
   * 
   * @param dir CRAM Read file directory, output from {@link CRAMSnpReader}
   */
  public void setCRAMReadDirectory(String dir) {
    for (String m : mosdepthFiles) {
      this.cramReadFiles.put(m, ext.rootOf(m, false) + CRAMSnpReader.CRAM_READS_EXT);
    }
  }

  Segment[] useBins;
  VCFFileReader genoReader;
  Map<String, String> genoIDLookup;

  LoadingCache<String, BEDFileReader> cramReadCache;
  LoadingCache<String, BEDFileReader> mosdepthCache;
  String[] mosdepthFiles;
  Map<String, String> cramReadFiles;
  Map<String, Double> medians;
  String[] samples;
  long fingerprintForMarkerFiles;

  MarkerDetailSet snpPosMDS;
  MarkerDetailSet binPosMDS;
  List<Marker> snpMarkers = new ArrayList<>();
  Map<Marker, Segment> binLookup = new HashMap<>();
  List<Marker> binMarkers = new ArrayList<>();
  Map<String, Marker> markerNameMap;
  String[][] binsOfMarkers;

  void run() throws IOException, Elision {
    if (log == null) {
      log = new Logger();
    }
    checkVars();
    createProject();

    loadBins();
    loadSNPs();
    createSampleList();
    loadMedians();
    mosdepthCache = buildCache();
    cramReadCache = buildCache();

    if (genoVCF != null) {
      log.reportTime("Opening genotype VCF file...");
      genoReader = new VCFFileReader(new File(genoVCF), true);
      List<String> genoSamples = genoReader.getFileHeader().getSampleNamesInOrder();
      genoIDLookup = new HashMap<>();
      for (String s : samples) {
        for (String g : genoSamples) {
          if (g.endsWith(s)) {
            genoIDLookup.put(s, g);
            break;
          }
        }
      }
    }

    read();

    cleanup();
  }

  private void checkVars() {
    if (projDir == null) {
      throw new IllegalArgumentException("No project directory set.");
    }
    if (propFileDir == null) {
      throw new IllegalArgumentException("No project property file directory set.");
    }
    if (projName == null) {
      throw new IllegalArgumentException("No project name set.");
    }
    if (mosdepthFiles == null) {
      throw new IllegalArgumentException("No mosdepth files set.");
    }
    if (useBed == null) {
      throw new IllegalArgumentException("No regions BED file set.");
    }
    if (markerVCF == null) {
      throw new IllegalArgumentException("No selected marker VCF file set.");
    }
    if (genoVCF == null) {
      throw new IllegalArgumentException("No genotype VCF file set.");
    }
  }

  private void cleanup() {
    if (genoReader != null) {
      genoReader.close();
      genoReader = null;
    }

    mosdepthCache.invalidateAll();
    mosdepthCache.cleanUp();
    mosdepthCache = null;

    binsOfMarkers = null;
    markerNameMap = null;
    binMarkers = null;
    binLookup = null;
    snpMarkers = null;
    binPosMDS = null;
    snpPosMDS = null;
    samples = null;
    mosdepthFiles = null;
    useBins = null;
  }

  private void createSampleList() {
    samples = new String[mosdepthFiles.length];
    for (int i = 0; i < mosdepthFiles.length; i++) {
      samples[i] = ext.rootRootOf(ext.removeDirectoryInfo((mosdepthFiles[i])));
    }
    fingerprintForMarkerFiles = MarkerSet.fingerprint(samples);
    if (!Files.exists(proj.SAMPLELIST_FILENAME.getValue())) {
      SampleList sl = new SampleList(samples);
      sl.serialize(proj.SAMPLELIST_FILENAME.getValue());
      sl = null;
      log.reportTime("Created sample list file: " + proj.SAMPLELIST_FILENAME.getValue());
    } else {
      log.reportTime("Project sample list file already exists; skipping creation.");
    }
  }

  private void loadMedians() {
    log.reportTime("Loading median depth values for autosomes...");
    long t1 = System.nanoTime();
    medians = new HashMap<>();
    for (int i = 0; i < mosdepthFiles.length; i++) {
      BEDFileReader reader = new BEDFileReader(mosdepthFiles[i], false);
      int sz = (int) reader.iterator().stream().filter((bf) -> {
        int c = (int) Positions.chromosomeNumber(bf.getContig());
        return c > 0 && c < 23;
      }).count();
      double median = ArrayUtils.median(reader.iterator().stream().filter((bf) -> {
        int c = (int) Positions.chromosomeNumber(bf.getContig());
        return c > 0 && c < 23;
      }).map(bf -> Double.parseDouble(bf.getName())), sz);
      reader.close();
      medians.put(samples[i], median);
    }
    log.reportTime("Computed autosomal median depth for " + medians.size() + " samples in "
                   + ext.getTimeElapsedNanos(t1));
  }

  private int getBaseIndex(String base) {
    switch (base) {
      case "A":
        return 0;
      case "G":
        return 1;
      case "C":
        return 2;
      case "T":
        return 3;
      case "N":
        return 4;
      default:
        return -1;
    }
  }

  private void read() throws IOException, Elision {
    long t2 = System.nanoTime();
    byte nullStatus = getNullStatus();
    int bytesPerSamp = Sample.getNBytesPerSampleMarker(nullStatus); // 12
    int markerBlockSize = samples.length * bytesPerSamp;

    ExecutorService exec = Executors.newFixedThreadPool(numThreads);
    // for each bin:
    for (int i = 0; i < binsOfMarkers.length; i++) {
      final int ind = i;
      exec.submit(new Runnable() {

        @Override
        public void run() {
          byte[] mkrBuff;
          long t1 = System.nanoTime();

          String[] markersInFile = binsOfMarkers[ind];
          String mdRAFName = "markers." + ind + MarkerData.MARKER_DATA_FILE_EXTENSION;
          log.reportTime("Parsing " + mdRAFName + "; " + (ind + 1) + " of " + binsOfMarkers.length
                         + " bins.");
          //    open and write header for mdRAF file containing current bin of markers:
          RandomAccessFile raf;
          try {
            raf = openMDRAF(mdRAFName, samples.length, nullStatus, fingerprintForMarkerFiles,
                            markersInFile);
          } catch (IOException e) {
            log.reportTime("CANNOT OPEN " + mdRAFName + ", skipping!");
            log.reportException(e);
            return;
          }
          Hashtable<String, Float> outOfRangeTable = new Hashtable<>();
          try {
            //    for each marker in bin:
            for (int m = 0; m < markersInFile.length; m++) {
              mkrBuff = new byte[markerBlockSize];

              Marker mkr = markerNameMap.get(markersInFile[m]);

              //        grab data from genotype vcf for marker:
              VariantContext match = null;
              if (genoReader != null) {
                List<VariantContext> genos;
                synchronized (genoReader) {
                  genos = genoReader.query(Positions.chromosomeNumberInverse(mkr.getChr()),
                                           mkr.getPosition() - 1, mkr.getPosition() + 1)
                                    .toList();
                }
                for (VariantContext vc : genos) {
                  if (vc.getID().equals(mkr.getName())) {
                    match = vc;
                    break;
                  } else {
                    Allele r = mkr.getRef();
                    Allele a = mkr.getAlt();
                    if (mkr.getGenomicPosition().getPosition() == vc.getStart()
                        && ((vc.getReference().basesMatch(r)
                             && vc.getAlternateAllele(0).basesMatch(a))
                            || (vc.getAlternateAllele(0).basesMatch(r)
                                && vc.getReference().basesMatch(a)))) {
                      match = vc;
                      break;
                    }
                  }
                }
                if (match == null) {
                  String msg = "Couldn't find an exact match between " + mkr.getName() + " and "
                               + genos.size() + " genotype records";

                  if (genos.size() > 0) {
                    msg += ": {";
                    for (int v = 0; v < genos.size(); v++) {
                      String gId = genos.get(v).getID();
                      if (gId.equals(".")) {
                        gId = genos.get(v).getContig() + ":" + genos.get(v).getStart();
                        if (genos.get(v).isIndel()) {
                          gId += "(I/D)";
                        }
                      }
                      msg += gId;
                      if (v < genos.size() - 1) {
                        msg += ", ";
                      }
                    }
                    msg += "}.";
                  } else {
                    msg += ".";
                  }
                  log.reportTimeWarning(msg);
                }
              }

              Map<String, Double> medianList = new HashMap<>();
              for (int s = 0; s < samples.length; s++) {
                double lrr;
                synchronized (samples[s]) {
                  lrr = loadMosdepth(match, mkr, s) / medians.get(samples[s]).doubleValue();
                }
                medianList.put(samples[s], lrr);
              }
              double markerMedian = ArrayUtils.median(medianList.values());

              //        for each sample in project,
              for (int s = 0; s < samples.length; s++) {

                double x = Double.NaN;
                double y = Double.NaN;

                CRAMRead reads = loadCRAMReads(mkr, s); // AGCTN
                if (!match.getReference().getBaseString().equals(reads.ref)) {
                  // TODO log mismatched ref allele between genotype vcf and selected-snp vcf!
                }
                int refInd = getBaseIndex(reads.ref);
                int altInd = getBaseIndex(reads.alt);

                x = reads.cnts[refInd] / scaleFactor;
                y = reads.cnts[altInd] / scaleFactor;

                double gc = 0.0;
                byte fg = 0;
                if (match != null) {

                  // TODO assign forward genotype and GC from allele counts 
                  Genotype g = match.getGenotype(genoIDLookup.get(samples[s]));
                  // gc
                  gc = g.getGQ() == -1 ? Double.NaN : ((double) g.getGQ()) / 100d;
                  // genotype
                  fg = (byte) ext.indexOfStr(g.getGenotypeString(), Sample.ALLELE_PAIRS);
                }

                // lrr
                double lrr = Maths.log2(medianList.get(samples[s]) / markerMedian);
                // baf
                double baf = 0.0;

                boolean noor;
                int sInd = s * bytesPerSamp;
                int bInd = sInd;
                Compression.gcBafCompress((float) gc, mkrBuff, bInd);
                bInd += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;

                noor = Compression.xyCompressPositiveOnly((float) x, mkrBuff, bInd);
                if (!noor) {
                  outOfRangeTable.put(m + "\t" + s + "\tx", (float) x);
                }
                bInd += Compression.REDUCED_PRECISION_XY_NUM_BYTES;

                noor = Compression.xyCompressPositiveOnly((float) y, mkrBuff, bInd);
                if (!noor) {
                  outOfRangeTable.put(m + "\t" + s + "\ty", (float) y);
                }
                bInd += Compression.REDUCED_PRECISION_XY_NUM_BYTES;

                Compression.gcBafCompress((float) baf, mkrBuff, bInd);
                bInd += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;

                noor = Compression.lrrCompress((float) lrr, mkrBuff, bInd) == 0;
                if (!noor) {
                  outOfRangeTable.put(m + "\t" + s + "\tlrr", (float) lrr);
                }
                bInd += Compression.REDUCED_PRECISION_LRR_NUM_BYTES;

                Compression.genotypeCompress((byte) -1, fg, mkrBuff, bInd);

                raf.write(mkrBuff);
              }
            }
          } catch (Elision e) {
            log.reportException(e);
          } catch (IOException e) {
            log.reportException(e);
          }

          try {
            //    write outliers and close file
            byte[] oorBytes = Compression.objToBytes(outOfRangeTable);
            raf.write(Compression.intToBytes(oorBytes.length));
            raf.write(oorBytes);
            raf.close();
          } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
          }
          log.reportTime("..... completed " + mdRAFName + " in " + ext.getTimeElapsedNanos(t1));
        }
      });
    }
    exec.shutdown();
    try {
      exec.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {
      log.reportException(e);
    }
    log.reportTime("Parsed " + binsOfMarkers.length + " marker data files in "
                   + ext.getTimeElapsedNanos(t2));
  }

  private LoadingCache<String, BEDFileReader> buildCache() {
    return CacheBuilder.newBuilder().softValues()
                       .removalListener(new RemovalListener<String, BEDFileReader>() {

                         @Override
                         public void onRemoval(RemovalNotification<String, BEDFileReader> notification) {
                           notification.getValue().close();
                         }
                       }).maximumSize(MAX_FILES_OPEN / 2).initialCapacity(MAX_FILES_OPEN / 2)
                       .build(new CacheLoader<String, BEDFileReader>() {

                         @Override
                         public BEDFileReader load(String key) throws Exception {
                           return BEDFileReader.createAnnotatedBEDFileReader(key, true);
                         }
                       });
  }

  private class CRAMRead {

    String ref;
    String alt;
    int[] cnts; // AGCTN
  }

  private CRAMRead loadCRAMReads(Marker mkr, int s) {
    String crF = cramReadFiles.get(mosdepthFiles[s]);
    BEDFileReader reader;
    boolean close = false;
    try {
      reader = cramReadCache.get(crF);
    } catch (ExecutionException e) {
      log.reportException(e);
      reader = BEDFileReader.createAnnotatedBEDFileReader(crF, true);
      close = true;
    }

    Segment seg = binLookup.get(mkr);
    List<BEDFeature> iter = reader.query(seg).toList();

    AnnotatedBEDFeature feat;
    if (iter.size() >= 1) {
      if (iter.size() > 1) {
        log.reportTimeWarning("Multiple depth scores found for " + samples[s] + ", bin "
                              + binLookup.get(mkr).getUCSClocation());
      }
      feat = (AnnotatedBEDFeature) iter.get(0);
    } else/* if (iter.size() == 0) */ {
      log.reportError("Missing depth score for " + samples[s] + ", bin "
                      + binLookup.get(mkr).getUCSClocation());
      feat = null;
    }
    CRAMRead read = null;
    if (feat != null) {
      read = new CRAMRead();
      read.ref = feat.getAnnotation(0);
      read.alt = feat.getAnnotation(1);
      read.cnts = new int[] {Integer.parseInt(feat.getAnnotation(2)),
                             Integer.parseInt(feat.getAnnotation(3)),
                             Integer.parseInt(feat.getAnnotation(4)),
                             Integer.parseInt(feat.getAnnotation(5)),};
    }

    if (close) {
      reader.close();
    }
    return read;
  }

  private double loadMosdepth(VariantContext match, Marker mkr, int s) {
    String file = mosdepthFiles[s];
    BEDFileReader reader;
    boolean close = false;
    try {
      reader = mosdepthCache.get(file);
    } catch (ExecutionException e) {
      log.reportException(e);
      reader = new BEDFileReader(file, true);
      close = true;
    }

    Segment seg = binLookup.get(mkr);
    List<BEDFeature> iter = reader.query(seg).toList();
    double v;
    if (iter.size() > 1) {
      log.reportTimeWarning("Multiple depth scores found for " + samples[s] + ", bin "
                            + binLookup.get(mkr).getUCSClocation());
      v = Double.parseDouble(iter.get(0).getName());
    } else if (iter.size() == 0) {
      log.reportError("Missing depth score for " + samples[s] + ", bin "
                      + binLookup.get(mkr).getUCSClocation());
      v = Double.NaN;
    } else {
      v = Double.parseDouble(iter.get(0).getName());
    }

    if (close) {
      reader.close();
    }

    return v;
  }

  private byte getNullStatus() {
    return Sample.computeNullStatus(new float[0], new float[0], new float[0], new float[0],
                                    new float[0], null, new byte[0], false);
  }

  private RandomAccessFile openMDRAF(String filename, int nInd, byte nullStatus, long fingerprint,
                                     String[] mkrNames) throws IOException {
    byte[] mkrBytes = Compression.objToBytes(mkrNames);
    byte[] mdRAFHeader = TransposeData.getParameterSectionForMdRaf(nInd, mkrNames.length,
                                                                   nullStatus, fingerprint,
                                                                   mkrBytes);
    mkrBytes = null;

    String file = proj.MARKER_DATA_DIRECTORY.getValue(true, true) + filename;
    if (Files.exists(file)) {
      new File(file).delete();
    }

    RandomAccessFile mdRAF = new RandomAccessFile(file, "rw");
    mdRAF.write(mdRAFHeader);
    mdRAFHeader = null;

    return mdRAF;
  }

  private void createProject() {
    String propFile = ext.verifyDirFormat(propFileDir)
                      + ext.replaceWithLinuxSafeCharacters(projName) + ".properties";
    if (!Files.exists(propFile)) {
      Files.write((new Project()).PROJECT_NAME.getName() + "=" + projName, propFile);
      proj = new Project(propFile);
      proj.PROJECT_NAME.setValue(projName);
      proj.PROJECT_DIRECTORY.setValue(projDir);
      proj.XY_SCALE_FACTOR.setValue(scaleFactor);
      proj.TARGET_MARKERS_FILENAMES.setValue(new String[] {});
      proj.ID_HEADER.setValue("NULL");
      proj.GENOME_BUILD_VERSION.setValue(GenomeBuild.HG38);

      proj.ARRAY_TYPE.setValue(ARRAY.NGS);

      proj.saveProperties();
      log.reportTime("Created project properties file: " + propFile);
    } else {
      log.reportTime("Project properties file already exists at " + propFile
                     + "; skipping creation.");
      proj = new Project(propFile);
    }
  }

  private void loadBins() {
    // 1-based indexing!
    BEDFileReader reader = new BEDFileReader(useBed, false);
    useBins = reader.loadAll(log).getStrictSegments();
    reader.close();
    log.reportTime("Loaded list of bins to be imported.");
  }

  private void loadSNPs() {
    VCFFileReader snpReader = new VCFFileReader(new File(markerVCF), true);
    for (Segment seg : useBins) {
      List<VariantContext> markerVC = snpReader.query(seg.getChromosomeUCSC(), seg.getStart(),
                                                      seg.getStop())
                                               .toList();
      if (markerVC.size() == 0) {
        log.reportTimeWarning("No snp found for bin " + seg.getUCSClocation()
                              + ". Bin will be skipped.");
        continue;
      }
      if (markerVC.size() > 1) {
        log.reportError("Multiple markers-to-use found for bin " + seg.getUCSClocation()
                        + ".  Bin will be skipped.");
        continue;
      }
      VariantContext chosenSnp = markerVC.get(0);

      Marker mSnpPos = new Marker(chosenSnp.getID(),
                                  new GenomicPosition(Positions.chromosomeNumber(chosenSnp.getContig()),
                                                      chosenSnp.getStart()),
                                  AllelePair.of(chosenSnp.getReference(),
                                                chosenSnp.getAlternateAlleles().get(0)));
      int stt = new Integer(chosenSnp.getAttribute("BINSTART").toString()).intValue();
      int stp = new Integer(chosenSnp.getAttribute("BINSTOP").toString()).intValue();
      int mid = stt + ((stp - stt) / 2) + 1;

      Marker mBinPos = new Marker(chosenSnp.getID(),
                                  new GenomicPosition(Positions.chromosomeNumber(chosenSnp.getContig()),
                                                      mid),
                                  AllelePair.of(chosenSnp.getReference(),
                                                chosenSnp.getAlternateAlleles().get(0)));
      binLookup.put(mSnpPos, seg);
      snpMarkers.add(mSnpPos);
      binMarkers.add(mBinPos);
    }
    snpReader.close();

    String f;
    if (!Files.exists(f = getMDSName(false))) {
      MarkerDetailSet mdsSnp = new MarkerDetailSet(snpMarkers);
      mdsSnp.serialize(f);
    } else {
      // TODO log
    }
    if (!Files.exists(f = getMDSName(true))) {
      MarkerDetailSet mdsBin = new MarkerDetailSet(binMarkers);
      mdsBin.serialize(f);
    } else {
      // TODO log
    }

    markerNameMap = new HashMap<>();
    int bins = (snpMarkers.size() / numMarkersPerFile);
    while ((bins * numMarkersPerFile) < snpMarkers.size()) {
      bins++;
    }
    int remainder = snpMarkers.size() - ((bins - 1) * numMarkersPerFile);
    binsOfMarkers = new String[bins][];
    for (int i = 0; i < binsOfMarkers.length - 1; i++) {
      binsOfMarkers[i] = new String[numMarkersPerFile];
    }
    binsOfMarkers[binsOfMarkers.length - 1] = new String[remainder];

    int currentBin = 0;
    int currentInd = 0;
    for (Marker m : snpMarkers) {
      markerNameMap.put(m.getName(), m);
      binsOfMarkers[currentBin][currentInd] = m.getName();
      currentInd++;
      if (currentInd == numMarkersPerFile) {
        currentBin++;
        currentInd = 0;
      }
    }

    log.reportTime("Found " + snpMarkers.size() + " markers to be parsed into " + bins
                   + " files with " + numMarkersPerFile + " markers per file.");
  }

  private String getMDSName(boolean binnedVersion) {
    return binnedVersion ? ext.rootOf(proj.MARKER_DETAILS_FILENAME.getValue(), false) + "_bins.ser"
                         : proj.MARKER_DETAILS_FILENAME.getValue();
  }

  public static void main(String[] args) throws IOException, Elision {
    MosdepthPipeline mi = new MosdepthPipeline();
    mi.setProjectDir("G:\\bamTesting\\topmed\\project\\");
    mi.setProjectName("TopmedMosdepth");
    mi.setProjectPropertiesDir("D:\\projects\\");
    mi.setNumThreads(Runtime.getRuntime().availableProcessors());

    mi.setBinsToUseBED("G:\\bamTesting\\snpSelection\\ReferenceGenomeBins.bed");
    //    mi.setGenotypeVCF("G:\\bamTesting\\EwingWGS\\ES_recalibrated_snps_indels.vcf.gz");
    mi.setSelectedMarkerVCF("G:\\bamTesting\\snpSelection\\selected.vcf");
    mi.setMosdepthDirectory("G:\\bamTesting\\topmed\\00src\\", ".bed.gz");
    mi.run();
  }

}
