package org.genvisis.one.JL;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Map;

import org.genvisis.cnv.filesys.CNVariant;
import org.genvisis.cnv.seq.qc.Mappability;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.LocusSet;
import org.pankratzlab.common.filesys.Segment;
import org.pankratzlab.common.stats.Rscript.RScatter;
import org.pankratzlab.common.stats.Rscript.SCATTER_TYPE;

public class CompExomeDepthConcordance {

  public static void comp(String[] files, Hashtable<String, String> iidMatch, String mapFile,
                          String callSubset, String problematicRegionFIle) {
    ArrayList<LocusSet<CNVariant>> sets = new ArrayList<>();
    String resultDir = ext.parseDirectoryOfFile(files[0]) + "results/";
    new File(resultDir).mkdirs();
    Logger log = new Logger();
    for (String file : files) {
      log.reportTimeInfo("Loading " + file);
      LocusSet<CNVariant> tmp = CNVariant.loadLocSet(file, log);
      Map<String, LocusSet<CNVariant>> inds = CNVariant.breakIntoInds(tmp, log);

      ArrayList<CNVariant> indsForThisFile = new ArrayList<>();
      for (String ind : iidMatch.keySet()) {
        if (inds.containsKey(ind)) {
          inds.get(ind).addAll(indsForThisFile);
        }
      }
      LocusSet<CNVariant> indSet = new LocusSet<CNVariant>(indsForThisFile.toArray(new CNVariant[indsForThisFile.size()]),
                                                           true, log) {

        /**
         *
         */
        private static final long serialVersionUID = 1L;

      };
      log.reportTimeInfo(file + " had " + indsForThisFile.size() + " cnvs matching inds");
      sets.add(indSet);
    }

    for (int i = 0; i < sets.size(); i++) {
      LocusSet<CNVariant> current = sets.get(i);
      ArrayList<CNVariant> match = new ArrayList<>();
      ArrayList<CNVariant> noMatch = new ArrayList<>();

      for (int j = 0; j < sets.size(); j++) {

        if (j != i) {
          LocusSet<CNVariant> comp = sets.get(j);
          for (int k = 0; k < current.getLoci().length; k++) {
            CNVariant currentCompare = current.getLoci()[k];
            CNVariant[] olaps = comp.getOverLappingLoci(currentCompare);
            boolean hasMatch = false;
            if (olaps != null && olaps.length > 0) {
              for (CNVariant olap : olaps) {
                if (olap.getCN() == currentCompare.getCN()) {
                  hasMatch = true;
                } else {
                  // noMatch.add(currentCompare);
                }
              }
            }
            if (hasMatch) {
              match.add(currentCompare);
            } else {
              noMatch.add(currentCompare);

            }

          }
        }
      }

      if (match.size() + noMatch.size() != current.getLoci().length) {
        throw new IllegalArgumentException("ERROR");
      }
      String outMatch = resultDir + ext.addToRoot(ext.removeDirectoryInfo(files[i]), ".matches");
      for (int j = 0; j < match.size(); j++) {
        // System.out.println(match.get(j).toPlinkFormat());
      }
      // System.exit(1);
      LocusSet<CNVariant> matchSet = new LocusSet<CNVariant>(match.toArray(new CNVariant[match.size()]),
                                                             true, log) {

        /**
         *
         */
        private static final long serialVersionUID = 1L;

      };
      dumpAndSummarize(matchSet, outMatch, mapFile, callSubset, problematicRegionFIle, log);

      String outNOMatch = resultDir + ext.addToRoot(ext.removeDirectoryInfo(files[i]), ".noMatch");

      LocusSet<CNVariant> noMatchSet = new LocusSet<CNVariant>(noMatch.toArray(new CNVariant[noMatch.size()]),
                                                               true, log) {

        /**
         *
         */
        private static final long serialVersionUID = 1L;

      };
      dumpAndSummarize(noMatchSet, outNOMatch, mapFile, callSubset, problematicRegionFIle, log);
      System.out.println(match.size());
      System.out.println(noMatch.size());

      // matchSet.writeRegions(outMatch, TO_STRING_TYPE.REGULAR, true, log);
      //
      // noMatchSet.writeRegions(outNOMatch, TO_STRING_TYPE.REGULAR, true, log);

    }
  }

  private static void dumpAndSummarize(LocusSet<CNVariant> out, String outFile,
                                       String mappabilityFile, String callSubsetBed,
                                       String problematicRegionFIle, Logger log) {
    Mappability<CNVariant> cnMappability = new Mappability<>(out, mappabilityFile, callSubsetBed,
                                                             log);
    cnMappability.computeMappability();
    LocusSet<Segment> pSet = LocusSet.loadSegmentSetFromFile(problematicRegionFIle, 0, 1, 2, 0,
                                                             true, true, 0, log);

    try {
      PrintWriter writer = Files.openAppropriateWriter(outFile);
      writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER)
                     + "\tmappability\tInProblematicRegion");
      for (int i = 0; i < out.getLoci().length; i++) {
        CNVariant tmp = out.getLoci()[i];
        writer.println(tmp.toPlinkFormat() + "\t"
                       + (cnMappability.getMappabilityResults().get(i).getAverageMapScore() * 100)
                       + "\t" + (pSet.getOverLappingLoci(tmp) != null));
      }

      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + outFile);
      log.reportException(e);
    }
    String[] ys = new String[] {"SCORE", "SITES", "mappability"};
    RScatter rsScatter = new RScatter(outFile, outFile + ".rscript",
                                      ext.removeDirectoryInfo(outFile), outFile + ".jpg", "TYPE",
                                      ys, SCATTER_TYPE.BOX, log);
    rsScatter.setOverWriteExisting(true);
    rsScatter.execute();
  }

  public static void main(String[] args) {
    String dir = "D:/data/Project_Tsai_21_25_26_spector/cnvs/compConcord/";
    String[] cnvFiles = Files.listFullPaths(dir, ".cnvs");
    Hashtable<String, String> iidMatch = new Hashtable<>();
    iidMatch.put("D100\tD100",
                 "HapMap_Control_CAGAGAGG-CTCTCTAT\tHapMap_Control_CAGAGAGG-CTCTCTAT");
    iidMatch.put("HapMap_Control_CAGAGAGG-CTCTCTAT\tHapMap_Control_CAGAGAGG-CTCTCTAT",
                 "D100\tD100");
    String mapFile = "C:/bin/ref/mappability/hg19/wgEncodeCrgMapabilityAlign100mer.bedGraph";
    String callSubset = "C:/bin/ExomeDepth/exons.hg19.chr.bed";
    String problematicRegionFIle = dir + "problematicRegions.txt";
    comp(cnvFiles, iidMatch, mapFile, callSubset, problematicRegionFIle);
  }

}
