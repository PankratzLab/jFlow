package org.genvisis.seq.analysis;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.StringJoiner;
import java.util.concurrent.Callable;

import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerHive;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Segment;

public class SimpleTallyGene {

  private static class Params implements Callable<Params> {

    private final String vcf;
    private final Collection<Segment> segs;
    private final String name;
    private final String vpopFile;
    private final String omimDir;
    private final VariantContextFilter filter;
    private final double maf;

    /**
     * @param vcf vcf file to use
     * @param seg seqment to tally
     * @param name name of this segment
     * @param vpopFile
     * @param omimDir
     * @param filter
     * @param maf max maf
     */
    private Params(String vcf, Collection<Segment> segs, String name, String vpopFile,
                   String omimDir, VariantContextFilter filter, double maf) {
      super();
      this.vcf = vcf;
      this.segs = segs;
      this.name = name;
      this.vpopFile = vpopFile;
      this.omimDir = omimDir;
      this.filter = filter;
      this.maf = maf;

    }

    @Override
    public Params call() throws Exception {
      String dir = ext.rootOf(vpopFile, false) + "_" + name + "/";
      new File(dir).mkdirs();
      Logger log = new Logger(dir + "log.log");

      String newVpop = dir + ext.removeDirectoryInfo(vpopFile);
      String segFile = dir + name + ".segs";
      StringJoiner segJoiner = new StringJoiner("\n");
      for (Segment seg : segs) {
        segJoiner.add(seg.getChr() + "\t" + seg.getStart() + "\t" + seg.getStop());
      }
      Files.write(segJoiner.toString(), segFile);
      String subVcf = VCFOps.extractSegments(vcf, segFile, 100, null,
                                             ext.rootOf(vpopFile, false) + "_extractedVcfs/", false,
                                             true, true, false, null, 1, log);
      Files.copyFile(vpopFile, newVpop);
      if (Files.exists(ext.rootOf(vpopFile, false) + ".lowerQualitySamples.txt")) {
        Files.copyFile(ext.rootOf(vpopFile, false) + ".lowerQualitySamples.txt",
                       dir + ext.rootOf(vpopFile, true) + ".lowerQualitySamples.txt");
      }
      VCFSimpleTally.test(subVcf, new String[] {newVpop}, omimDir, null, null, maf, true, true,
                          filter, false, new HashSet<>(), 24);
      return this;
    }

  }

  /**
   * @param vcf
   * @param vpop
   * @param maf
   * @param seg
   * @param name
   * @param omimDir
   */
  public static void run(String vcf, String vpop, double maf, Collection<Segment> segs, String name,
                         String omimDir) {

    WorkerHive<Params> hive = new WorkerHive<>(1, 1, new Logger());
    hive.addCallable(new Params(vcf, segs, name, vpop, omimDir, null, maf));
    hive.execute(true);

  }

  @Deprecated
  public static void run(String vcf, String vpop, double[] mafs) {

    String[] names = new String[] {"Mito"};
    Segment[] segs = new Segment[] {new Segment("chrM:1-20000")};
    mafs = mafs == null ? new double[] {1.2} : mafs;
    String tnVpop = "/Volumes/Beta/data/Cushings/mito/CUSHINGS_TUMOR.vpop";

    WorkerHive<Params> hive = new WorkerHive<>(1, 1, new Logger());
    String omimDir = "/Volumes/Beta/ref/OMIM/";
    for (int i = 0; i < names.length; i++) {
      for (double maf : mafs) {
        hive.addCallable(new Params(vcf, Arrays.asList(segs[i]), names[i], vpop, omimDir, null,
                                    maf));
        if (maf == 0) {
          hive.addCallable(new Params(vcf, Arrays.asList(segs[i]), names[i], tnVpop, omimDir, null,
                                      maf));
        }
      }
    }
    hive.execute(true);

  }

  /**
   * @param args
   */
  public static void main(String[] args) {
    CLI c = new CLI(SimpleTallyGene.class);

    c.addArgWithDefault("vcf", "vcf to tally", "a.vcf");
    c.addArgWithDefault("vpop", "vpop to use", "a.vpop");
    c.addArgWithDefault("mafs", "mafs to use,comma delimited", "0.0,1.2");
    c.addArgWithDefault("segment", "UCSC segments , ; delimited",
                        "chr18:20714428-20840534;chr17:7571720-7590868");
    c.addArgWithDefault("name", "typically gene name, comma delimited", "CABLES1,TP53");
    c.addArgWithDefault("omimDir", "omim directory", "/Volumes/Beta/ref/OMIM/");
    c.addFlag("combine", "If multiple segments are provided, combine for a single analysis");

    c.parseWithExit(args);

    String vcf = c.get("vcf");
    String vpop = c.get("vpop");
    double[] mafs = ArrayUtils.toDoubleArray(c.get("mafs").split(","));

    Segment[] segs = Segment.getSegments(c.get("segment").split(";"));
    String[] names = c.get("name").split(",");
    String omimDir = c.get("omimDir");

    for (double maf : mafs) {
      if (c.has("combine")) {
        if (names.length != 1) {
          throw new IllegalArgumentException("Must have a single name when combining");
        }
        run(vcf, vpop, maf, Arrays.asList(segs), names[0], omimDir);
      } else {
        if (names.length != segs.length) {
          throw new IllegalArgumentException("Must have a name for each segment");
        }
        for (int i = 0; i < names.length; i++) {
          run(vcf, vpop, maf, Arrays.asList(segs[i]), names[i], omimDir);
        }
      }
    }

  }

}
