package org.genvisis.cnv.qc;

import java.util.List;
import java.util.Set;

import org.genvisis.cnv.Resources;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.filesys.Segment;
import org.pankratzlab.utils.gwas.Qc;
import org.pankratzlab.utils.gwas.RelationAncestryQc;

import com.google.common.collect.Sets;

public class IntensityMarkers {

  Project proj;
  List<String> oneHits;
  Set<String> pruneKeeps;
  Set<String> problematicDrops;
  Set<String> markers;

  private IntensityMarkers(Project proj) {
    this.proj = proj;
  }

  public static Set<String> getIntensityMarkers(Project proj) {
    IntensityMarkers im = new IntensityMarkers(proj);
    im.getOneHitWonders();
    im.loadProblemRegions();
    im.loadPruneList();
    im.filter();
    return im.markers;
  }

  private static List<String> blastHits(Project proj) {
    return MarkerBlastQC.getOneHitWonders(proj, proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                          MarkerBlastQC.DEFAULT_CROSS_HYBE_THRESHOLD,
                                          proj.getLog());
  }

  private void getOneHitWonders() {
    if (Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue())) {
      proj.getLog().report("Loading one-hit wonders markers from blast vcf.");
      oneHits = blastHits(proj);
    } else {
      switch (proj.ARRAY_TYPE.getValue()) {
        case AFFY_AXIOM:
          if (!Files.exists(proj.INTENSITY_PC_MARKERS_FILENAME.getValue())) {
            // TODO change this when Genvisis can import Axiom arrays as RSIDs
            Files.copyFile(Resources.axiomTx(proj.getLog())
                                    .genome(proj.GENOME_BUILD_VERSION.getValue())
                                    .getIntensityPCMarkers_ProbeSets().get(),
                           proj.INTENSITY_PC_MARKERS_FILENAME.getValue());
          }
          oneHits = HashVec.loadFileToVec(proj.INTENSITY_PC_MARKERS_FILENAME.getValue(), false,
                                          false, false);
          break;
        default:
          throw new IllegalArgumentException("BLAST VCF is required for one-hit wonders file creation.");
      }
    }
  }

  private void loadProblemRegions() {
    String probRgns = Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), proj.getLog())
                               .getProblematicRegions().get();
    String[] rgns = HashVec.loadFileToStringArray(probRgns, false, new int[] {0}, false);
    problematicDrops = Sets.newHashSet();
    for (String r : rgns) {
      Segment seg = new Segment(r);
      proj.getMarkerSet().viewMarkersInSeg(seg).map(Marker::getName).forEach(problematicDrops::add);
    }
  }

  private void loadPruneList() {
    // TODO determine "plink" dir name elsewhere?
    String file = proj.PROJECT_DIRECTORY.getValue() + "plink/" + Qc.QC_SUBDIR
                  + RelationAncestryQc.LD_PRUNING_DIR + "plink.prune.out";
    if (Files.exists(file)) {
      proj.getLog().reportTime("Retaining only plink-pruned markers found in " + file + ".");
      pruneKeeps = HashVec.loadFileToHashSet(file, new int[] {0}, "", false);
    }
  }

  private void filter() {
    markers = Sets.newHashSet(oneHits);
    markers.removeAll(problematicDrops);
    if (pruneKeeps != null) {
      markers.retainAll(pruneKeeps);
    }
  }

}