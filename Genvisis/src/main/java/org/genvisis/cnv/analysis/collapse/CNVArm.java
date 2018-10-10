/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import org.genvisis.cnv.Resources.GENOME_BUILD;
import org.genvisis.cnv.filesys.CNVariant;
import org.genvisis.cnv.filesys.Project;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerTrain;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.WorkerTrain.Producer;
import org.pankratzlab.common.CLI;
import org.pankratzlab.shared.filesys.LocusSet;
import org.pankratzlab.shared.filesys.Segment;
import org.pankratzlab.shared.filesys.LocusSet.TO_STRING_TYPE;

/**
 * Methods to replace CNVs with a forced call across an entire arm - e.g. when an entire chromosomal
 * arm is duplicated but PennCNV only called portion
 */
public class CNVArm extends CNVariant {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  /**
   * 
   */

  private static final String[] INPUT_FILE_HEADER = new String[] {"DNA", "ARM", "ANNOTATION"};

  private int numDupsOverlapped;
  private int numDelsOverlapped;
  private int lengthDupsOverlapped;
  private int lengthDelsOverlapped;

  private static String[] getArmHeader() {
    return new String[] {"NumDupsOverlapped", "NumDelsOverlapped", "lengthDupsOverlapped",
                         "lengthDelsOverlapped", "propDup", "propDel", "PredictedCN"};

  }

  private CNVArm(CNVBuilder builder, int numDupsReplaced, int numDelsReplaced,
                 int lengthDupsOverlapped, int lengthDelsOverlapped) {
    super(builder);
    this.numDupsOverlapped = numDupsReplaced;
    this.numDelsOverlapped = numDelsReplaced;
    this.lengthDupsOverlapped = lengthDupsOverlapped;
    this.lengthDelsOverlapped = lengthDelsOverlapped;
  }

  private int getPredictedCN() {
    if (lengthDelsOverlapped > 2 * lengthDupsOverlapped) {
      return 1;
    } else if (lengthDupsOverlapped > 2 * lengthDelsOverlapped) {
      return 3;
    } else {
      return 2;
    }
  }

  private double getPropDup() {
    return (double) lengthDupsOverlapped / getSize();
  }

  private double getPropDel() {
    return (double) lengthDelsOverlapped / getSize();
  }

  @Override
  public String toAnalysisString() {
    return super.toAnalysisString() + "\t" + numDupsOverlapped + "\t" + numDelsOverlapped + "\t"
           + lengthDupsOverlapped + "\t" + lengthDelsOverlapped + "\t" + getPropDup() + "\t"
           + getPropDel() + "\t" + getPredictedCN();
  }

  @Override
  public String[] getHeader() {
    return ArrayUtils.concatAll(super.getHeader(), getArmHeader());
  }

  private static <T extends CNVariant> LocusSet<CNVArm> replace(LocusSet<T> set,
                                                                LocusSet<T> replacers, Logger log) {
    if (replacers.getLoci().length > 0) {
      List<CNVArm> replaced = new ArrayList<>();
      for (T t : set.getLoci()) {
        T[] olaps = replacers.getOverLappingLoci(t);
        if (olaps == null || olaps.length == 0) {// does not overlap any of the force call regions
          replaced.add(new CNVArm(new CNVBuilder(t), t.getCN() < 2 ? 1 : 0, t.getCN() > 2 ? 1 : 0,
                                  0, 0));
        } else if (olaps.length > 0) {
          for (T olap : olaps) {
            if (olap.amountOfOverlapInBasepairs(t) != t.getSize()) {// completely contained in olap
              LocusSet<Segment> clean = t.removeAll(replacers.getLoci(), log);
              log.reportTimeInfo("converting overlap for " + t.toAnalysisString() + " from "
                                 + olap.toAnalysisString() + " to " + clean.getLoci().length
                                 + " separate segments");
              for (Segment seg : clean.getLoci()) {
                CNVBuilder builder = new CNVBuilder(t);
                builder.chr(seg.getChr());
                builder.start(seg.getStart());
                builder.stop(seg.getStop());
                replaced.add(new CNVArm(builder, t.getCN() < 2 ? 1 : 0, t.getCN() > 2 ? 1 : 0, 0,
                                        0));
              }
            }
          }
        }
      }
      for (T t : replacers.getLoci()) {
        T[] olaps = set.getOverLappingLoci(t);
        int numDupsOverlapped = 0;
        int numDelsOverlapped = 0;
        int lengthDupsOverlapped = 0;
        int lengthDelsOverlapped = 0;

        if (olaps != null) {
          for (T to : olaps) {
            if (to.getCN() < 2) {
              numDelsOverlapped++;
              lengthDelsOverlapped += to.getSize();
            } else if (to.getCN() > 2) {
              numDupsOverlapped++;
              lengthDupsOverlapped += to.getSize();
            }
          }
        }
        replaced.add(new CNVArm(new CNVBuilder(t).score(1000), numDupsOverlapped, numDelsOverlapped,
                                lengthDupsOverlapped, lengthDelsOverlapped));
      }
      return new LocusSet<>(replaced, true, log);
    } else {
      return null;
    }
  }

  private static LocusSet<CNVariant> loadArmVariants(String file, GENOME_BUILD build,
                                                     boolean includeMos, Logger log) {
    List<CNVariant> tmps = new ArrayList<>();

    if (!Files.headerOfFileContainsAll(file, INPUT_FILE_HEADER, log)) {
      throw new IllegalArgumentException("Invalid header for file " + file + ", need "
                                         + ArrayUtils.toStr(INPUT_FILE_HEADER));
    }
    try {
      BufferedReader reader = Files.getAppropriateReader(file);
      int[] indices = ext.indexFactors(INPUT_FILE_HEADER, reader.readLine().trim().split("\t"),
                                       true);
      while (reader.ready()) {
        String[] line = reader.readLine().trim().split("\t");

        Segment seg = getLoc(build, log, line[indices[1]]);
        int CN = 2;
        if (line[indices[2]].equals("DEL")) {
          CN = -1;
        } else if (line[indices[2]].equals("DUP")) {
          CN = 3;
        }
        System.out.println(line[indices[2]] + "\t" + CN);
        if (CN != 2 | includeMos) {
          CNVariant tmp = new CNVariant(line[indices[0]], line[indices[0]], seg.getChr(),
                                        seg.getStart(), seg.getStop(), CN, Double.NaN, -1, -1);
          tmps.add(tmp);

        }
      }
      reader.close();

    } catch (IOException e) {
      log.reportException(e);
    }
    if (tmps.isEmpty()) {
      log.reportError("no regions found to force call");
      return null;
    }
    return new LocusSet<>(tmps, true, log);
  }

  private static Segment getLoc(GENOME_BUILD build, Logger log, String chrarm) {
    //    int[] loc = Positions.parseUCSClocation(chrarm, Centromeres.getCentromereMidPoints(build, log),
    //                                            Centromeres.getChrLength(build, log));
    //    if (build != null) {
    throw new IllegalArgumentException("Fix centromere lookup first");
    //    }
    //    return new Segment((byte) loc[0], loc[1], loc[2]);
    //    return null;
  }

  private static LocusSet<CNVariant> convert(LocusSet<CNVariant> set, String dna, Logger log) {

    ArrayList<CNVariant> conv = new ArrayList<>();
    HashSet<String> added = new HashSet<>();
    for (CNVariant cnv : set.getLoci()) {
      CNVariant ctmp = new CNVBuilder(cnv).familyID(dna).individualID(dna).build();
      if (!added.contains(ctmp.toAnalysisString())) {
        added.add(ctmp.toAnalysisString());
        conv.add(ctmp);
      }
    }
    return new LocusSet<>(conv, true, log);
  }

  private static LocusSet<CNVariant> getAllArms(GENOME_BUILD build, Logger log) {
    List<CNVariant> cnvs = new ArrayList<>();
    for (int i = 1; i < 23; i++) {
      cnvs.add(new CNVariant(getLoc(build, log, "chr" + i + "p").getUCSClocation()));
      cnvs.add(new CNVariant(getLoc(build, log, "chr" + i + "q").getUCSClocation()));
    }
    return new LocusSet<>(cnvs, true, log);
  }

  private static class CNVForceCallerProducer implements Producer<CNVForceCaller> {

    private final String[] samples;
    private final boolean[] markersToUse;
    private int index;
    private final Project proj;
    private final GENOME_BUILD build;
    private final int[][] indicesByChr;

    private CNVForceCallerProducer(String[] samples, Project proj, GENOME_BUILD build,
                                   boolean[] markersToUse, int[][] indicesByChr) {
      super();
      this.samples = samples;
      this.proj = proj;
      this.build = build;
      this.markersToUse = markersToUse;
      this.indicesByChr = indicesByChr;
    }

    /*
     * (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    @Override
    public boolean hasNext() {
      return index < samples.length;
    }

    /*
     * (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    @Override
    public Callable<CNVForceCaller> next() {
      final String sample = samples[index];
      index++;
      return new Callable<CNVForceCaller>() {

        @Override
        public CNVForceCaller call() throws Exception {
          proj.getLog().reportTimeInfo(sample);
          return CNVForceCaller.call(proj, proj.getFullSampleFromRandomAccessFile(sample),
                                     convert(getAllArms(build, proj.getLog()), sample,
                                             proj.getLog()),
                                     CNVForceCaller.VALID_METHODS, markersToUse, null,
                                     indicesByChr);
        }
      };
    }

    /*
     * (non-Javadoc)
     * @see org.genvisis.common.WorkerTrain.Producer#shutdown()
     */
    @Override
    public void shutdown() {

    }

  }

  private static boolean[] getMarkersToUse(Project proj) {
    String[] names = proj.getMarkerNames();
    boolean[] use = ArrayUtils.booleanArray(names.length, true);
    HashSet<String> excludes = new HashSet<>();
    if (proj.FILTERED_MARKERS_FILENAME.getValue().equals("")) {
      excludes = new HashSet<>();
    } else if (Files.exists(proj.FILTERED_MARKERS_FILENAME.getValue())) {
      excludes = HashVec.loadFileToHashSet(proj.FILTERED_MARKERS_FILENAME.getValue(), false);
    } else {
      excludes = new HashSet<>();
      ;
    }
    proj.getLog().reportTimeInfo(excludes.size() + " markers loaded from "
                                 + proj.FILTERED_MARKERS_FILENAME.getValue());
    for (int i = 0; i < use.length; i++) {
      if (excludes.contains(names[i])) {
        use[i] = false;
      }
    }
    return use;
  }

  private static void run(Project proj, String cnvFile, String armFile, GENOME_BUILD build,
                          boolean includeMos, boolean annotate) {

    LocusSet<CNVariant> all = CNVariant.loadLocSet(cnvFile, proj.getLog());
    LocusSet<CNVariant> replaces = loadArmVariants(armFile, build, includeMos, proj.getLog());
    String outDir = proj.PROJECT_DIRECTORY.getValue() + "armForceCall/";
    new File(outDir).mkdirs();
    replaces.writeRegions(outDir + ext.removeDirectoryInfo(armFile) + ".cnv",
                          TO_STRING_TYPE.REGULAR, true, proj.getLog());
    Map<String, LocusSet<CNVariant>> indsAll = CNVariant.breakIntoInds(all, proj.getLog());
    Map<String, LocusSet<CNVariant>> indsReplace = CNVariant.breakIntoInds(replaces, proj.getLog());
    proj.getLog().reportTimeInfo("Loaded " + all.getLoci().length + " cnvs from "
                                 + indsAll.keySet().size() + " samples");
    proj.getLog().reportTimeInfo("Loaded " + replaces.getLoci().length + " arms to force call from "
                                 + indsReplace.keySet().size() + " samples");
    String outFile = outDir + ext.removeDirectoryInfo(armFile) + ".replaced.cnv";

    String outFileAllReplaced = outDir + ext.removeDirectoryInfo(armFile) + ".all.replaced.cnv";
    String outCallFile = outDir + ext.removeDirectoryInfo(armFile) + ".all.forceCall.cnv";

    proj.getLog().reportTimeInfo(getAllArms(build, proj.getLog()).getLoci().length + " total arms");
    if (annotate) {
      PrintWriter writer = Files.getAppropriateWriter(outCallFile);
      boolean header = true;

      boolean[] markersToUse = getMarkersToUse(proj);

      CNVForceCallerProducer producer = new CNVForceCallerProducer(proj.getSamples(), proj, build,
                                                                   markersToUse,
                                                                   proj.getMarkerSet()
                                                                       .getIndicesByChr());
      WorkerTrain<CNVForceCaller> train = new WorkerTrain<>(producer, 4, 4, proj.getLog());
      int num = 0;
      while (train.hasNext()) {
        CNVForceCaller called = train.next();
        called.dumpToFile(writer, header, proj.getLog());
        header = false;
        proj.getLog().reportTimeInfo("On sample n=" + num);
        num++;
      }
      writer.close();
    }

    List<CNVariant> replacedCNVs = new ArrayList<>();

    for (String sample : indsAll.keySet()) {
      if (indsReplace.containsKey(sample) && indsReplace.get(sample).getLoci().length > 0) {

        LocusSet<CNVariant> current = indsAll.containsKey(sample) ? indsAll.get(sample)
                                                                  : new LocusSet<>(new CNVariant[0],
                                                                                   true,
                                                                                   proj.getLog());
        LocusSet<CNVArm> tmp = replace(current, indsReplace.get(sample), proj.getLog());
        for (CNVArm c : tmp.getLoci()) {
          replacedCNVs.add(c);
        }
      } else {
        for (CNVariant c : indsAll.get(sample).getLoci()) {
          replacedCNVs.add(new CNVArm(new CNVBuilder(c), 0, 0, 0, 0));
        }
      }
    }
    new LocusSet<>(replacedCNVs, true, proj.getLog()).writeRegions(outFile, TO_STRING_TYPE.REGULAR,
                                                                   true, proj.getLog());

    List<CNVariant> allreplacedCNVs = new ArrayList<>();

    for (String dna : proj.getSamples()) {
      String sample = dna + "\t" + dna;
      LocusSet<CNVariant> replaceHere = convert(getAllArms(GENOME_BUILD.HG19, proj.getLog()), dna,
                                                proj.getLog());
      if (indsAll.containsKey(sample) && indsAll.get(sample).getLoci().length > 0) {

        LocusSet<CNVariant> current = indsAll.containsKey(sample) ? indsAll.get(sample)
                                                                  : new LocusSet<>(new CNVariant[0],
                                                                                   true,
                                                                                   proj.getLog());
        LocusSet<CNVArm> tmp = replace(current, replaceHere, proj.getLog());
        for (CNVArm c : tmp.getLoci()) {
          allreplacedCNVs.add(c);
        }
      } else {
        for (CNVariant c : replaceHere.getLoci()) {
          allreplacedCNVs.add(new CNVArm(new CNVBuilder(c), 0, 0, 0, 0));
        }
      }
    }
    new LocusSet<>(allreplacedCNVs, true, proj.getLog()).writeRegions(outFileAllReplaced,
                                                                      TO_STRING_TYPE.REGULAR, true,
                                                                      proj.getLog());

  }

  public static void main(String[] args) {
    CLI c = new CLI(CNVArm.class);
    c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ);
    c.addArg("cnvFile", "cnvFile to armitize");
    c.addArg("armFile", "file of arms to force call");
    c.addArg("build", "genome build, options are " + ArrayUtils.toStr(GENOME_BUILD.values()), ",");

    c.parseWithExit(args);
    boolean annotate = false;

    run(new Project(c.get(CLI.ARG_PROJ)), c.get("cnvFile"), c.get("armFile"),
        GENOME_BUILD.valueOf(c.get("build")), false, true);
    Files.copyFileUsingFileChannels(c.get("armFile"), c.get("armFile") + ".mos", new Logger());
    run(new Project(c.get(CLI.ARG_PROJ)), c.get("cnvFile"), c.get("armFile") + ".mos",
        GENOME_BUILD.valueOf(c.get("build")), true, annotate);

  }
}
