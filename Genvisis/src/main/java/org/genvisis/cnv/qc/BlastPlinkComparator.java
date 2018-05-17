package org.genvisis.cnv.qc;

import java.util.HashSet;
import java.util.Set;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;
import com.google.common.collect.Sets;

public class BlastPlinkComparator {

  Project proj;

  public BlastPlinkComparator(Project project) {
    this.proj = proj;
  }

  public void run() {
    Set<String> oneHits = loadOneHitters();
    Set<String> plinkQC = loadPlinkQCPassingMarkers();

    int qcPassClean = 0;
    int qcPassNoClean = 0;
    int qcFailClean = 0;
    int qcFailNoClean = 0;

    String[] allMkrs = proj.getMarkerNames();
    for (String s : allMkrs) {
      boolean hit = oneHits.contains(s);
      boolean pqc = plinkQC.contains(s);

      if (hit && pqc) {
        qcPassClean++;
      } else if (hit && !pqc) {
        qcFailClean++;
      } else if (!hit && pqc) {
        qcPassNoClean++;
      } else if (!hit && !pqc) {
        qcFailNoClean++;
      }

    }

    int maxLen = 0;
    int len = 0;
    if ((len = Integer.toString(qcPassClean).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcPassNoClean).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcFailClean).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcFailNoClean).length()) > maxLen) {
      maxLen = len;
    }

    int allMkrsLen = allMkrs.length;

    String rep1 = ext.formStr(Integer.toString(qcPassClean), maxLen) + " ("
                  + ext.formDeci(((double) qcPassClean * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that pass marker QC that have a clean exact BLAST hit";
    String rep2 = ext.formStr(Integer.toString(qcFailClean), maxLen) + " ("
                  + ext.formDeci(((double) qcFailClean * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that do not pass marker QC that have a clean exact BLAST hit";
    String rep3 = ext.formStr(Integer.toString(qcPassNoClean), maxLen) + " ("
                  + ext.formDeci(((double) qcPassNoClean * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that pass marker QC that do not have an exact BLAST hit";
    String rep4 = ext.formStr(Integer.toString(qcFailNoClean), maxLen) + " ("
                  + ext.formDeci(((double) qcFailNoClean * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that do not pass marker QC that do not have an exact BLAST hit";

    System.out.println(rep1);
    System.out.println(rep2);
    System.out.println(rep3);
    System.out.println(rep4);
  }

  private Set<String> loadOneHitters() {
    return Sets.newHashSet(MarkerBlastQC.getOneHitWonders(proj,
                                                          proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                                          MarkerBlastQC.DEFAULT_CROSS_HYBE_THRESHOLD,
                                                          proj.getLog()));
  }

  private Set<String> loadPlinkQCPassingMarkers() {
    String[] allMarkers = proj.getMarkerNames();
    String[] droppedMkrs = HashVec.loadFileToStringArray(proj.PROJECT_DIRECTORY.getValue()
                                                         + "plink/marker_qc/mind_drops.dat", false,
                                                         new int[] {0}, false);
    Set<String> drops = Sets.newHashSet(droppedMkrs);
    Set<String> returnMkrs = new HashSet<>();
    for (String s : allMarkers) {
      if (!drops.contains(s)) {
        returnMkrs.add(s);
      }
    }
    return returnMkrs;
  }

  public static void main(String[] args) {
    CLI cli = new CLI(BlastPlinkComparator.class);

    cli.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ);

    cli.parseWithExit(args);

    new BlastPlinkComparator(new Project(cli.get(CLI.ARG_PROJ))).run();

  }

}
