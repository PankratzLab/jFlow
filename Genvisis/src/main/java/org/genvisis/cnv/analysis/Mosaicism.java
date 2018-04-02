// -Xms1024M -Xmx1024M
package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;
import java.util.concurrent.Callable;
import javax.swing.JOptionPane;
import org.genvisis.cnv.analysis.MosaicismDetect.MosaicBuilder;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerSet.PreparedMarkerSet;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.plots.MosaicPlot;
import org.genvisis.cnv.var.IndiPheno;
import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import com.google.common.primitives.Floats;

public class Mosaicism {

  public static final String[] HEADER = {"Sample", "Arm", "#CNVs", "Summed_Size", "%covered",
                                         "Custom_metric", "LRR_SD", "LRR_SD_flag",
                                         "Flagged_via_Mosaicism", "Mosaicism_level",
                                         "Mosaicism description"};
  public static final double LOWER_BOUND = 0.15;
  public static final double UPPER_BOUND = 0.85;

  public static void findOutliers(Project proj) {
    findOutliers(proj, -1);
  }

  public static void findOutliers(Project proj, int numthreads) {
    final PrintWriter writer;
    String[] samples;
    int chr;
    Hashtable<String, String> hash = new Hashtable<>();
    String[] markerNames;
    byte[] chrs;
    int[] positions;
    boolean[] snpDropped;
    int[][] chrBoundaries;
    MarkerDetailSet markerSet;

    hash = proj.getFilteredHash();

    chrBoundaries = new int[27][3];
    snpDropped = null;
    for (int i = 0; i < chrBoundaries.length; i++) {
      chrBoundaries[i][0] = chrBoundaries[i][1] = chrBoundaries[i][2] = -1;
    }
    markerSet = proj.getMarkerSet();
    markerNames = markerSet.getMarkerNames();
    chrs = markerSet.getChrs();
    positions = markerSet.getPositions();
    snpDropped = new boolean[markerNames.length];
    int[][] indicesByChr = markerSet.getIndicesByChr();
    System.out.println("Mosacism will be estimated using " + markerNames.length + " markers");
    chr = 0;
    for (int i = 0; i < markerNames.length; i++) {
      snpDropped[i] = hash.containsKey(markerNames[i]);
      if (positions[i] > Positions.CENTROMERE_MIDPOINTS[chr] && chrBoundaries[chr][1] == -1) {
        chrBoundaries[chr][1] = i;
      }
      if (chrs[i] > chr || i == markerNames.length - 1) {
        if (chr != 0) {
          chrBoundaries[chr][2] = i - 1;
        }
        chr = chrs[i];
        chrBoundaries[chr][0] = i;
      }
    }
    chrBoundaries[0][0] = 0;
    chrBoundaries[0][2] = markerNames.length - 1;

    for (int i = 0; i < chrBoundaries.length; i++) {
      if (chrBoundaries[i][0] == -1 || chrBoundaries[i][2] == -1) {
        System.err.println("Error - no data for chromosome '" + i + "'");
      }
    }

    PSF.checkInterrupted();
    samples = proj.getSamples();
    try {
      writer = Files.openAppropriateWriter(proj.MOSAIC_RESULTS_FILENAME.getValue(true, true));
      writer.println(ArrayUtils.toStr(MosaicPlot.MOSAICISM_HEADER));
      // samples = new String[] { "7355066051_R03C01", "7330686030_R02C01", "7159911135_R01C02" };
      // samples = new String[] { "7355066051_R03C01" };

      MosaicResultProducer producer = new MosaicResultProducer(proj, samples, snpDropped,
                                                               chrBoundaries, markerSet,
                                                               indicesByChr);
      try (WorkerTrain<String[]> train = new WorkerTrain<>(producer,
                                                           numthreads > 0 ? numthreads
                                                                          : proj.NUM_THREADS.getValue(),
                                                           2, proj.getLog())) {
        int index = 0;
        long timePer = System.currentTimeMillis();
        long time = System.currentTimeMillis();

        while (train.hasNext()) {
          PSF.checkInterrupted(new Runnable() {

            @Override
            public void run() {
              writer.close();
            }
          });
          try {
            String[] results = train.next();
            index++;
            if (index % numthreads == 0) {
              proj.getLog()
                  .reportTimeInfo((index) + " of " + samples.length + " in "
                                  + ext.getTimeElapsed(timePer) + ", total time at "
                                  + ext.getTimeElapsed(time));
              timePer = System.currentTimeMillis();
            }

            for (String result : results) {
              writer.println(result);
            }
          } catch (Exception e) {
            proj.getLog().reportException(e);
            producer.shutdown();
          }

        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to summary file");
      e.printStackTrace();
    }
  }

  private static class MosaicResultProducer extends AbstractProducer<String[]> {

    private final Project proj;
    private final String[] samples;
    private final boolean[] snpDropped;
    private final int[][] chrBoundaries;
    private final PreparedMarkerSet markerSet;
    private final int[][] indicesByChr;
    private int index;

    public MosaicResultProducer(Project proj, String[] samples, boolean[] snpDropped,
                                int[][] chrBoundaries, MarkerSetInfo markerSet,
                                int[][] indicesByChr) {
      super();
      this.proj = proj;
      this.samples = samples;
      this.snpDropped = snpDropped;
      this.chrBoundaries = chrBoundaries;
      this.indicesByChr = indicesByChr;
      this.markerSet = PreparedMarkerSet.getPreparedMarkerSet(markerSet);
      index = 0;
    }

    @Override
    public boolean hasNext() {
      return index < samples.length;
    }

    @Override
    public Callable<String[]> next() {
      final String currentSample = samples[index];
      Callable<String[]> callable = new Callable<String[]>() {

        @Override
        public String[] call() throws Exception {

          return getMosaicResults(proj, currentSample, snpDropped, chrBoundaries, markerSet,
                                  indicesByChr);
        }
      };
      index++;
      return callable;
    }

    @Override
    public void shutdown() {
      index = samples.length;
    }
  }

  private static String[] getMosaicResults(Project proj, String sample, boolean[] snpDropped,
                                           int[][] chrBoundaries, MarkerSetInfo markerSet,
                                           int[][] indicesByChr) {
    Sample samp;
    float baf;
    float[] lrrs;
    float[] bafs;
    samp = proj.getPartialSampleFromRandomAccessFile(sample);
    ArrayList<String> results = new ArrayList<>();
    if (samp.getFingerprint() != markerSet.getFingerprint()) {
      String error = "Error - cannot estimate mosaics if MarkerSet and Sample (" + sample
                     + ") don't use the same markers";
      throw new IllegalArgumentException(error);
    }
    lrrs = samp.getLRRs();
    bafs = samp.getBAFs();
    MosaicBuilder builder = new MosaicBuilder();
    builder.indicesByChr(indicesByChr);
    builder.verbose(false);
    builder.markerIndices(proj.getMarkerIndices());
    if (proj.getArrayType() == ARRAY.NGS) {
      proj.getLog().reportTimeWarning("Masking non-variant sites for project type " + ARRAY.NGS);
      boolean[] use = new boolean[markerSet.getMarkerNames().length];
      int numMasked = 0;
      for (int i = 0; i < use.length; i++) {
        boolean useit = !proj.getArrayType().isCNOnly(markerSet.getMarkerNames()[i]);
        if (!useit) {
          numMasked++;
        }
        use[i] = useit;
      }
      proj.getLog().reportTimeInfo(numMasked + " markers were masked");
      builder.use(use);
    }
    MosaicismDetect md = builder.build(proj, sample, markerSet, ArrayUtils.toDoubleArray(bafs));
    int[] positions = markerSet.getPositions();
    byte[] chrs = markerSet.getChrs();
    for (int j = 1; j <= 23; j++) {
      for (int arm = 0; arm < 2; arm++) {
        int startIndex = (arm == 0 ? chrBoundaries[j][0] : chrBoundaries[j][1]);
        int stopIndex = (arm == 0 ? chrBoundaries[j][1] : chrBoundaries[j][2] + 1);

        ArrayList<Float> lrrAl = new ArrayList<>(stopIndex + 10 - startIndex);
        ArrayList<Float> bafAl = new ArrayList<>(stopIndex + 10 - startIndex);
        for (int k = startIndex; k < stopIndex; k++) {
          if (!snpDropped[k] && !(md.getUse() != null && !md.getUse()[k])) {
            if (!Float.isNaN(lrrs[k])) {
              lrrAl.add(lrrs[k]);
            }
            baf = bafs[k];
            if (baf > LOWER_BOUND && baf < UPPER_BOUND) {
              bafAl.add(baf);
            }
          }
        }
        if (lrrAl.size() > 100) {
          if (chrs[startIndex] != (byte) j || chrs[stopIndex - 1] != (byte) j) {

            throw new IllegalStateException("Internal Error, mismatched chromosome indices, start ="
                                            + chrs[startIndex] + "\tstop = " + chrs[stopIndex]
                                            + "\tarm = " + arm);
          }
          Segment armSeg = new Segment((byte) j, positions[startIndex], positions[stopIndex - 1]);
          MosaicMetric mosaicMetrics = getMosiacMetric(md, armSeg, proj.getLog());

          int bafSize = bafAl.size();
          int lrrSize = lrrAl.size();
          float[] bafTmp = Floats.toArray(bafAl);
          String result = sample + "\t" + "chr" + j + (arm == 0 ? "p" : "q") + "\t" + lrrSize + "\t"
                          + ext.formDeci(ArrayUtils.mean(Floats.toArray(lrrAl)), 5) + "\t"
                          + bafAl.size()
                          + (bafSize > 10 ? "\t" + ext.formDeci(ArrayUtils.stdev(bafTmp, true), 5)
                                            + "\t"
                                            + ext.formDeci(ArrayUtils.iqrExclusive(bafTmp), 5)
                                          : "\t.\t.")
                          + "\t" + ext.formDeci((double) (lrrSize - bafSize) / (double) lrrSize, 5);
          result += "\t" + mosaicMetrics.getForcedCallproportionArmMosaic();
          result += "\t" + mosaicMetrics.getBpWeightedAverageArm();
          result += "\t" + mosaicMetrics.getBpWeightedAverageCalled();
          result += "\t" + mosaicMetrics.getMosaicRegionsDetected();
          result += "\t" + mosaicMetrics.getBpMosiac();
          result += "\t" + mosaicMetrics.getBpArm();
          result += "\t" + mosaicMetrics.getProportionBpCalledMosaic();
          results.add(result);
        }
      }
    }
    return ArrayUtils.toStringArray(results);
  }

  private static class MosaicMetric {

    private double forcedCallproportionArmMosaic;
    private double bpWeightedAverageArm;
    private double bpWeightedAverageCalled;
    private int bpMosiac;
    private int bpArm;
    private double proportionBpCalledMosaic;
    private int mosaicRegionsDetected;

    public MosaicMetric(double forcedCallproportionArmMosaic, double bpWeightedAverageArm,
                        double bpWeightedAverageCalled, int mosaicRegionsDetected, int bpMosiac,
                        int bpArm, double proportionBpCalledMosaic) {
      super();
      this.forcedCallproportionArmMosaic = forcedCallproportionArmMosaic;
      this.bpWeightedAverageArm = bpWeightedAverageArm;
      this.bpWeightedAverageCalled = bpWeightedAverageCalled;
      this.bpMosiac = bpMosiac;
      this.bpArm = bpArm;
      this.proportionBpCalledMosaic = proportionBpCalledMosaic;
      this.mosaicRegionsDetected = mosaicRegionsDetected;

    }

    private int getMosaicRegionsDetected() {
      return mosaicRegionsDetected;
    }

    private double getBpWeightedAverageCalled() {
      return bpWeightedAverageCalled;
    }

    private void setBpWeightedAverageCalled(double bpWeightedAverageCalled) {
      this.bpWeightedAverageCalled = bpWeightedAverageCalled;
    }

    private void setMosaicRegionsDetected(int mosaicRegionsDetected) {
      this.mosaicRegionsDetected = mosaicRegionsDetected;
    }

    private int getBpMosiac() {
      return bpMosiac;
    }

    private void setBpMosiac(int bpMosiac) {
      this.bpMosiac = bpMosiac;
    }

    private int getBpArm() {
      return bpArm;
    }

    private void setBpArm(int bpArm) {
      this.bpArm = bpArm;
    }

    private double getProportionBpCalledMosaic() {
      return proportionBpCalledMosaic;
    }

    private void setProportionBpCalledMosaic(double proportionBpCalledMosaic) {
      this.proportionBpCalledMosaic = proportionBpCalledMosaic;
    }

    private double getForcedCallproportionArmMosaic() {
      return forcedCallproportionArmMosaic;
    }

    private void setForcedCallproportionArmMosaic(double percentArmMosaic) {
      forcedCallproportionArmMosaic = percentArmMosaic;
    }

    private double getBpWeightedAverageArm() {
      return bpWeightedAverageArm;
    }

    private void setBpWeightedAverageArm(double bpWeightedAverage) {
      bpWeightedAverageArm = bpWeightedAverage;
    }

  }

  private static MosaicMetric getMosiacMetric(MosaicismDetect md, Segment seg, Logger log) {
    LocusSet<MosaicRegion> mosSet = md.callMosaic(seg, true);
    MosaicMetric mosaicMetric = new MosaicMetric(-1, -1, -1, -1, -1, -1, -1);
    int bpArm = seg.getSize();
    mosaicMetric.setBpArm(bpArm);
    if (mosSet.getLoci().length != 1 || !seg.equals(mosSet.getLoci()[0])) {
      log.reportError("Mosaic caller not in force call mode");
      log.reportError(seg.getUCSClocation() + " went in, and "
                      + mosSet.getLoci()[0].getUCSClocation() + " came out");
    } else if (seg.getChr() < 23) {// can't call chr23 yet
      mosaicMetric.setForcedCallproportionArmMosaic(mosSet.getLoci()[0].getScore());
    }

    LocusSet<MosaicRegion> tmp = md.callMosaic(seg, false);

    mosaicMetric.setMosaicRegionsDetected(tmp.getLoci().length);
    if (tmp.getLoci().length < 1 || seg.getChr() >= 23) {

    } else {
      int bpCalledMosaic = 0;
      double bpWeightedSum = 0;
      int numRegions = 0;
      for (int i = 0; i < tmp.getLoci().length; i++) {
        if (tmp.getLoci()[i].getNumMarkers() > 2 * md.getMovingFactor()) {// a healthy filter
          bpWeightedSum += tmp.getLoci()[i].getSize() * tmp.getLoci()[i].getScore();
          numRegions++;
          bpCalledMosaic += tmp.getLoci()[i].getSize();
        }
      }
      mosaicMetric.setBpMosiac(bpCalledMosaic);
      double proportionCalledMosaic = (double) bpCalledMosaic / bpArm;
      mosaicMetric.setProportionBpCalledMosaic(proportionCalledMosaic);
      double bpWeightedAverageArm = 0;
      double bpWeightedAverageCalled = 0;
      bpWeightedAverageCalled = bpWeightedSum / bpCalledMosaic;
      bpWeightedAverageArm = bpWeightedSum / bpArm;
      mosaicMetric.setMosaicRegionsDetected(numRegions);
      mosaicMetric.setBpWeightedAverageArm(bpWeightedAverageArm);
      mosaicMetric.setBpWeightedAverageCalled(bpWeightedAverageCalled);
    }
    return mosaicMetric;
  }

  public static void checkForOverlap(Project proj, String listOfMosaicArms) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Vector<String> v;
    Vector<String[]> list;
    int count;
    SampleData sampleData;
    int chr, sum;
    Segment arm;
    IndiPheno indiPheno;
    CNVariant[] cnvs;
    String[] sampleList;
    String[][] listOfArms;
    Hashtable<String, String> lrrsdHash, mosaicHash;
    double proportion;
    long time;
    String[] cnvFiles;

    time = new Date().getTime();
    cnvFiles = proj.CNV_FILENAMES.getValue();
    if (cnvFiles.length == 0) {
      System.err.println("Error - need to specify the name of a CNV file in the project properties file before running Mosaicism.checkForOverlap()");
      return;
    }
    sampleData = proj.getSampleData(new String[] {cnvFiles[0]});
    if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + "lrr_sd.xln")) {
      lrrsdHash = HashVec.loadFileToHashString(proj.PROJECT_DIRECTORY.getValue() + "lrr_sd.xln",
                                               false);
    } else {
      System.err.println("Warning - could not find 'lrr_sd.xln' in project directory; no flags will be generated");
      lrrsdHash = new Hashtable<>();
    }
    if (Files.exists(proj.MOSAIC_COLOR_CODES_FILENAME.getValue(false, false))) {
      mosaicHash = HashVec.loadFileToHashString(proj.MOSAIC_COLOR_CODES_FILENAME.getValue(false,
                                                                                          false),
                                                new int[] {0, 1}, new int[] {2, 3}, false, "\t",
                                                true, true);
    } else {
      System.err.println("Warning - could not find "
                         + proj.MOSAIC_COLOR_CODES_FILENAME.getValue(false, false)
                         + "; no annotation possible");
      mosaicHash = new Hashtable<>();
    }

    if (listOfMosaicArms.toLowerCase().endsWith("/all")) {
      sampleList = sampleData.getListOfSamples();
      listOfArms = new String[sampleList.length * 22 * 2][2];
      for (int i = 0; i < sampleList.length; i++) {
        for (int j = 0; j < 22; j++) {
          listOfArms[i * 22 * 2 + j * 2 + 0][0] = sampleList[i];
          listOfArms[i * 22 * 2 + j * 2 + 0][1] = "chr" + (j + 1) + "p";
          listOfArms[i * 22 * 2 + j * 2 + 1][0] = sampleList[i];
          listOfArms[i * 22 * 2 + j * 2 + 1][1] = "chr" + (j + 1) + "q";
        }
      }
    } else {
      list = new Vector<>();
      try {
        reader = new BufferedReader(new FileReader(listOfMosaicArms));
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        if (!ext.checkHeader(line, new String[] {"Sample", "Arm"}, false)) {
          reader.close();
          return;
        }

        while (reader.ready()) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          list.add(new String[] {line[0], line[1]});
        }
        reader.close();

      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + listOfMosaicArms
                           + "\" not found in current directory");
        return;
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + listOfMosaicArms + "\"");
        return;
      }
      listOfArms = Matrix.toStringArrays(list);
    }

    v = new Vector<>();
    try {
      writer = Files.openAppropriateWriter(ext.rootOf(listOfMosaicArms, false) + "_counts.xln");
      writer.println(ArrayUtils.toStr(HEADER));
      for (String[] listOfArm : listOfArms) {
        indiPheno = sampleData.getIndiPheno(listOfArm[0]);

        if (indiPheno == null) {
          HashVec.addIfAbsent(listOfArm[0], v);
        } else {
          chr = Integer.parseInt(listOfArm[1].substring(3, listOfArm[1].length() - 1));
          if (listOfArm[1].charAt(listOfArm[1].length() - 1) == 'p') {
            arm = new Segment((byte) chr, 0, Positions.CENTROMERE_MIDPOINTS[chr]);
          } else {
            arm = new Segment((byte) chr, Positions.CENTROMERE_MIDPOINTS[chr],
                              Positions.CHROMOSOME_LENGTHS_B36_HG18[chr]);
          }

          cnvs = indiPheno.getCNVs(0, chr);
          if (cnvs == null) {
            cnvs = new CNVariant[0];
          }

          count = 0;
          sum = 0;
          for (CNVariant cnv : cnvs) {
            if (cnv.overlaps(arm)) {
              count++;
              sum += cnv.amountOfOverlapInBasepairs(arm);
            }
          }
          proportion = (double) sum / (double) arm.getSize();
          writer.print(listOfArm[0] + "\t" + listOfArm[1]);
          writer.print("\t" + count + "\t" + sum + "\t" + proportion + "\t"
                       + (400 * proportion + count));
          if (lrrsdHash.containsKey(listOfArm[0])) {
            writer.print("\t" + lrrsdHash.get(listOfArm[0]) + "\t"
                         + (Double.parseDouble(lrrsdHash.get(listOfArm[0])) > 0.28 ? 1 : 0));
          } else {
            writer.print("\t.\t.");
          }
          if (mosaicHash.containsKey(listOfArm[0] + "\t" + listOfArm[1])) {
            writer.print("\t1\t" + mosaicHash.get(listOfArm[0] + "\t" + listOfArm[1]));
          } else {
            writer.print("\t0\t.\t.");
          }
          writer.println();
        }
      }
      writer.close();

      if (v.size() > 0) {
        writer = Files.openAppropriateWriter(proj.PROJECT_DIRECTORY.getValue()
                                             + "SAMPLES_IN_CNVFILE_NOT_IN_SAMPLE_DATA.txt");
        for (int i = 0; i < v.size(); i++) {
          writer.println(v.elementAt(i));
        }
        writer.close();
        JOptionPane.showMessageDialog(null,
                                      "There were " + v.size()
                                            + " samples not present in the SampleData file; check file in the project directory for a list",
                                      "Error", JOptionPane.ERROR_MESSAGE);

      }
    } catch (IOException ioe) {
      System.err.println("Error writing to file \"" + ext.rootOf(listOfMosaicArms, false)
                         + "_counts.xln" + "\"");
      return;
    }
    System.out.println("Finished in " + ext.getTimeElapsed(time));
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    String arms = "MosaicArms.txt";
    Project proj;
    boolean check = false;
    int numthreads = 24;

    String usage = "\n" + "filesys.ParseIllumina requires 0-1 arguments\n"
                   + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + "   (2) check for overlap between mosaic arms and CNV calls in the first CNV file listed in the project file (i.e. -check (not the default))\n"
                   + "   (3) mosaic arms file (i.e. arms=MosaicArms.txt (default))\n"
                   + PSF.Ext.getNumThreadsCommand(4, numthreads) +

                   "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        return;
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("arms=")) {
        arms = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-check")) {
        check = true;
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numthreads = ext.parseIntArg(arg);
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      return;
    }

    try {
      proj = new Project(filename);
      if (check) {
        checkForOverlap(proj, arms);
      } else {
        findOutliers(proj);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
