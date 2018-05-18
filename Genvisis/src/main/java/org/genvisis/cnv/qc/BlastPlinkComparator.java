package org.genvisis.cnv.qc;

import java.util.HashSet;
import java.util.Set;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.MarkerBlastQC.QCResults;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;
import com.google.common.collect.Sets;

public class BlastPlinkComparator {

  Project proj;

  public BlastPlinkComparator(Project project) {
    this.proj = project;
  }

  public void run() {
    QCResults oneHits = loadOneHitters();
    Set<String> plinkQC = loadPlinkQCPassingMarkers();

    int qcPassPerf = 0;
    int qcFailPerf = 0;
    int qcPassClose = 0;
    int qcFailClose = 0;
    int qcPassAmbig = 0;
    int qcFailAmbig = 0;
    int qcPassBad = 0;
    int qcFailBad = 0;
    int qcPassMiss = 0;
    int qcFailMiss = 0;

    String[] allMkrs = proj.getMarkerNames();
    for (String s : allMkrs) {
      boolean hit = oneHits.getPerfect().contains(s);
      boolean clo = oneHits.getMatches().contains(s);
      boolean amb = oneHits.getAmbig().contains(s);
      boolean bad = oneHits.getBad().contains(s);
      boolean miss = oneHits.getMissing().contains(s);
      boolean pqc = plinkQC.contains(s);

      if (pqc && hit) {
        qcPassPerf++;
      } else if (!pqc && hit) {
        qcFailPerf++;
      } else if (pqc && clo) {
        qcPassClose++;
      } else if (!pqc && clo) {
        qcFailClose++;
      } else if (pqc && amb) {
        qcPassAmbig++;
      } else if (!pqc && amb) {
        qcFailAmbig++;
      } else if (pqc && bad) {
        qcPassBad++;
      } else if (!pqc && bad) {
        qcFailBad++;
      } else if (pqc && miss) {
        qcPassMiss++;
      } else if (!pqc && miss) {
        qcFailMiss++;
      }

    }

    int maxLen = 0;
    int len = 0;
    if ((len = Integer.toString(qcPassPerf).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcFailPerf).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcPassClose).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcFailClose).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcPassAmbig).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcFailAmbig).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcPassBad).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcFailBad).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcPassMiss).length()) > maxLen) {
      maxLen = len;
    }
    if ((len = Integer.toString(qcFailMiss).length()) > maxLen) {
      maxLen = len;
    }

    int allMkrsLen = allMkrs.length;

    String rep1 = ext.formStr(Integer.toString(qcPassPerf), maxLen) + " ("
                  + ext.formDeci(((double) qcPassPerf * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that pass marker QC that have a clean exact BLAST hit";
    String rep2 = ext.formStr(Integer.toString(qcFailPerf), maxLen) + " ("
                  + ext.formDeci(((double) qcFailPerf * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that do not pass marker QC that have a clean exact BLAST hit";
    String rep3 = ext.formStr(Integer.toString(qcPassClose), maxLen) + " ("
                  + ext.formDeci(((double) qcPassClose * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that pass marker QC that have a single close, but not exact, BLAST hit";
    String rep4 = ext.formStr(Integer.toString(qcFailClose), maxLen) + " ("
                  + ext.formDeci(((double) qcFailClose * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that do not pass marker QC that have a single close, but not exact, BLAST hit";
    String rep5 = ext.formStr(Integer.toString(qcPassAmbig), maxLen) + " ("
                  + ext.formDeci(((double) qcPassAmbig * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that pass marker QC that have a multiple close, but not exact, BLAST hits";
    String rep6 = ext.formStr(Integer.toString(qcFailAmbig), maxLen) + " ("
                  + ext.formDeci(((double) qcFailAmbig * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that do not pass marker QC that have a multiple close, but not exact, BLAST hits";
    String rep7 = ext.formStr(Integer.toString(qcPassBad), maxLen) + " ("
                  + ext.formDeci(((double) qcPassBad * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that pass marker QC that do not have any close BLAST hits";
    String rep8 = ext.formStr(Integer.toString(qcFailBad), maxLen) + " ("
                  + ext.formDeci(((double) qcFailBad * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that do not pass marker QC that do not have any close BLAST hits";
    String rep9 = ext.formStr(Integer.toString(qcPassMiss), maxLen) + " ("
                  + ext.formDeci(((double) qcPassMiss * 100) / allMkrsLen, 2, 2, true)
                  + "%) Number of markers that pass marker QC that were not in the BLAST results.";
    String rep10 = ext.formStr(Integer.toString(qcFailMiss), maxLen) + " ("
                   + ext.formDeci(((double) qcFailMiss * 100) / allMkrsLen, 2, 2, true)
                   + "%) Number of markers that do not pass marker QC that were not in the BLAST results.";

    System.out.println(rep1);
    System.out.println(rep2);
    System.out.println(rep3);
    System.out.println(rep4);
    System.out.println(rep5);
    System.out.println(rep6);
    System.out.println(rep7);
    System.out.println(rep8);
    System.out.println(rep9);
    System.out.println(rep10);
  }

  private QCResults loadOneHitters() {
    return MarkerBlastQC.getSingleHitMarkers(proj, proj.BLAST_ANNOTATION_FILENAME.getValue());
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
