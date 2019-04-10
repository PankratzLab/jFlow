/**
 * 
 */
package org.genvisis.cnv.hmm;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Stream;
import org.genvisis.cnv.filesys.CNVariant;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.filesys.LocusSet;
import com.google.common.math.Stats;
import com.google.common.math.StatsAccumulator;

/**
 * Generates estimates of LRR/BAF distributions using existing calls, which can come from an
 * orthogonal technology, an un-optimized .hmm, etc
 */
public class HMMParameterizer {

  /**
   * 
   */
  private static final String MAX_CNV = "maxCNV";
  /**
   * 
   */
  private static final String MIN_P = "minP";
  /**
   * 
   */
  private static final String MIN_BP = "minBP";
  /**
   * 
   */
  private static final String CNV_FILE = "cnvFile";

  static void run(CLI c) {
    Project proj = new Project(c.get(CLI.ARG_PROJ));
    String cnvFile = c.get(CNV_FILE);
    double minbp = c.getD(MIN_BP);
    double minP = c.getD(MIN_P);
    int maxCalls = c.getI(MAX_CNV);
    LocusSet<CNVariant> tmp = CNVariant.loadLocSet(cnvFile, proj.getLog());
    Map<String, LocusSet<CNVariant>> set = CNVariant.breakIntoInds(tmp, proj.getLog());
    proj.getLog().reportTimeInfo("Found " + tmp.getLoci().length + " cnvs across " + set.size()
                                 + " samples");

    List<StatsAccumulator> lrrstats = new ArrayList<>();
    for (int i = 0; i < 6; i++) {
      lrrstats.add(new StatsAccumulator());
    }
    List<StatsAccumulator> bafstats = new ArrayList<>();
    for (int i = 0; i < 5; i++) {
      bafstats.add(new StatsAccumulator());
    }
    MarkerDetailSet markerDetailSet = proj.getMarkerSet();
    Map<Marker, Integer> markerIndexMap = markerDetailSet.getMarkerIndexMap();
    for (String sample : set.keySet()) {
      LocusSet<CNVariant> sampSet = set.get(sample);
      if (sampSet.getLoci().length < maxCalls) {
        HashSet<Marker> notAccountedFor = new HashSet<>();
        notAccountedFor.addAll(markerDetailSet.markersAsList());
        String dna = sample.split("\t")[0];
        Sample samp = proj.getFullSampleFromRandomAccessFile(dna);
        if (samp == null) {
          String error = "Sample " + dna + " not found, ensure DNA is used in cnv file " + cnvFile;
          proj.getLog().reportError(error);
          throw new IllegalArgumentException(error);
        }

        proj.getLog().reportTimeInfo("Processing sample " + dna);
        double[] lrrs = ArrayUtils.toDoubleArray(samp.getLRRs());
        double[] bafs = ArrayUtils.toDoubleArray(samp.getBAFs());
        byte[] genos = samp.getAB_Genotypes();

        for (CNVariant cnv : sampSet.getLoci()) {
          Stream<Marker> cnvMarks = markerDetailSet.viewMarkersInSeg(cnv);
          cnvMarks.forEach(cnvMark -> {
            double lrr = lrrs[markerIndexMap.get(cnvMark)];
            if (Double.isFinite(lrr) && cnv.getSize() > minbp && cnv.getNumMarkers() > minP) {
              lrrstats.get(cnv.getCN()).add(lrr);
            }

            double baf = bafs[markerIndexMap.get(cnvMark)];
            byte geno = genos[markerIndexMap.get(cnvMark)];

            if (Double.isFinite(baf) && geno >= 0 && cnv.getSize() > minbp
                && cnv.getNumMarkers() > minP) {
              // TODO needs review
              switch (cnv.getCN()) {
                case 0:
                  if (geno != 1) {
                    bafstats.get(4).add(baf);
                  }
                  break;
                case 1:
                  if (geno != 1) {
                    bafstats.get(0).add(Math.min(baf, 1 - baf));
                  }
                  break;
                case 2:
                  if (geno == 1) {
                    bafstats.get(3).add(baf);
                  } else {
                    bafstats.get(0).add(Math.min(baf, 1 - baf));
                  }
                  break;
                case PennHmm.LOH_FLAG:
                  if (geno != 1) {
                    bafstats.get(0).add(Math.min(baf, 1 - baf));
                  }
                  break;
                case 3:

                  if (geno != 1) {
                    bafstats.get(0).add(Math.min(baf, 1 - baf));
                  } else {
                    bafstats.get(1).add(baf);
                  }
                  break;

                case 4:
                  if (geno != 1) {
                    bafstats.get(0).add(Math.min(baf, 1 - baf));
                  } else {
                    bafstats.get(2).add(baf);
                  }
                  break;

                default:
                  throw new IllegalArgumentException("Invalid CN state " + cnv.getCN());
              }
            }

            notAccountedFor.remove(cnvMark);
          });
        }

        for (Marker notAccountedForMarker : notAccountedFor) {
          double lrr = lrrs[markerIndexMap.get(notAccountedForMarker)];
          if (Double.isFinite(lrr)) {
            lrrstats.get(2).add(lrr);
          }
        }
        proj.getLog().reportTimeInfo("LRR:\n" + getSnapShot(lrrstats));
        proj.getLog().reportTimeInfo("BAF:\n" + getSnapShot(bafstats));

      }
    }

  }

  private static String getSnapShot(List<StatsAccumulator> stats) {
    StringJoiner joiner = new StringJoiner("\n");
    for (int i = 0; i < stats.size(); i++) {
      if (stats.get(i).count() > 1) {
        Stats s = stats.get(i).snapshot();
        joiner.add("STATE=" + i + ";MEAN=" + s.mean() + ";SD=" + s.sampleStandardDeviation());
      }
    }

    return joiner.toString();

  }

  public static void main(String[] args) {
    CLI c = new CLI(HMMParameterizer.class);
    c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ, true);
    c.addArg(CNV_FILE, "use these cnvs to parameterize the .hmm file", true);
    c.addArgWithDefault(MIN_BP, "minimum size in bps", "10000");
    c.addArgWithDefault(MIN_P, "minimum size in probes", "10");
    c.addArgWithDefault(MAX_CNV, "maximum number cnvs per sample", "150");

    c.parseWithExit(args);
    run(c);
  }

}
