package org.genvisis.one.JL.ks;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.genvisis.cnv.analysis.pod.PODAnalysis;
import org.genvisis.cnv.analysis.pod.PODAnalysis.PODResults;
import org.genvisis.cnv.analysis.pod.PODAnalysis.SEARCH_SPACE_TYPE;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.filesys.FamilyStructure;
import org.pankratzlab.common.filesys.LocusSet;
import org.pankratzlab.common.filesys.Segment;

public class testKSPOD {

  static List<Segment> searchSpace = Arrays.asList(new Segment[] {new Segment("chrX:142000000-155000000"),
                                                                  new Segment("chrX:51000000-97000000"),
                                                                  new Segment("chrX:18000000-38000000"),
                                                                  new Segment("chrX:41000000-62000000"),
                                                                  new Segment("chrX:115000000-146000000"),
                                                                  new Segment("chrX:31000000-135000000")});

  enum TYPE {
    FULL, RESTRICTED;
  }

  public static void main(String[] args) {
    Project proj = new Project("/Users/Kitty/.genvisis/projects/Poynter1.properties");
    String ped = "/Volumes/Beta/data/Poynter/Poynter1/pedigree_use.dat";

    Project proj2 = new Project("/Users/Kitty/.genvisis/projects/Poynter2.properties");
    String ped2 = "/Volumes/Beta/data/Poynter/Poynter2/pedigree.dat";

    boolean runMos = false;

    if (runMos) {
      //
      // String outDir = proj.PROJECT_DIRECTORY.getValue() + "mos/";
      // new File(outDir).mkdirs();
      // MosaicismDetect.callMosaicRegions(proj, outDir + "mos.mos", 2);
      //
      // String outDir2 = proj2.PROJECT_DIRECTORY.getValue() + "mos/";
      // new File(outDir2).mkdirs();
      // MosaicismDetect.callMosaicRegions(proj2, outDir2 + "mos.mos", 2);

    }

    else {
      List<Segment> fullX = Arrays.asList(new Segment[] {new Segment("chrX:0-2147483600")});
      LocusSet<Segment> set = new LocusSet<>(fullX, true, proj.getLog());

      List<Segment> restricted = Arrays.asList(set.removeThese(new LocusSet<>(searchSpace, true,
                                                                              proj.getLog()),
                                                               0)
                                                  .getLoci());

      for (TYPE type : TYPE.values()) {
        runKS(proj, proj.getMarkerSet(), ped, type == TYPE.FULL ? fullX : restricted, type,
              type == TYPE.FULL ? SEARCH_SPACE_TYPE.INDIVIDUAL : SEARCH_SPACE_TYPE.INCLUSIVE);

        runKS(proj2, proj2.getMarkerSet(), ped2, type == TYPE.FULL ? fullX : restricted, type,
              type == TYPE.FULL ? SEARCH_SPACE_TYPE.INDIVIDUAL : SEARCH_SPACE_TYPE.INCLUSIVE);
      }
    }

  }

  private static void runKS(Project proj, MarkerDetailSet markerDetailSet, String ped,
                            List<Segment> segments, TYPE type, SEARCH_SPACE_TYPE searchType) {

    String outDir = proj.PROJECT_DIRECTORY.getValue() + "KS/";
    new File(outDir).mkdirs();
    String outFile = outDir + type + "ks.summary.txt";
    proj.getLog().reportTimeInfo(outFile);
    // if (!Files.exists(outFile)) {
    PrintWriter writer = Files.getAppropriateWriter(outFile);

    Pedigree pedigree = new Pedigree(proj, ped);
    FamilyStructure.PedigreeUtils.loadPOPairs(pedigree, false, null, null, true);
    FamilyStructure.PedigreeUtils.loadCompleteTrios(pedigree, null, null, true);

    // Pedigree.l
    ArrayList<int[]> completePairs = pedigree.cached_complete_trios;
    ArrayList<String[]> poPairs = pedigree.cached_poPairsIDs;
    HashMap<String, Integer> map = pedigree.getfidiidToIndexMap();

    //
    writer.println(ArrayUtils.toStr(PODResults.getHeader()) + "\tSEX\tPROJECT\tHAS_FULL");
    HashSet<String> completed = new HashSet<>();
    for (int[] pair : completePairs) {
      if (pair != null) {
        // pedigree.
        String offDNA = pedigree.getiDNA(pair[0]);
        String faDNA = pedigree.getiDNA(pair[1]);
        String moDNA = pedigree.getiDNA(pair[2]);
        completed.add(offDNA);
        List<PODResults> results = PODAnalysis.analyze(proj, markerDetailSet, segments, offDNA,
                                                       moDNA, faDNA, searchType, 0);

        writer.flush();
        String stat = "NA";
        if (pedigree.getGender(pair[0]) == (byte) 1) {
          stat = "MALE";
        } else if (pedigree.getGender(pair[0]) == (byte) 2) {
          stat = "FEMALE";
        }
        for (PODResults poResults : results) {
          writer.println(poResults.toString() + "\t" + stat + "\t" + proj.PROJECT_NAME.getValue()
                         + "\tTRUE");
        }
      }
    }

    for (String[] pair : poPairs) {
      if (pair != null) {
        System.out.println(ArrayUtils.toStr(pair));
        if (map.containsKey(pair[0])) {
          int offIndex = map.get(pair[1]);
          String offDNA = pedigree.getiDNA(offIndex);
          System.out.println(offDNA + '\t' + pair[0] + "\t" + pedigree.getFaDNAIndex(offIndex)
                             + "\t" + pedigree.getMoDNAIndex(offIndex));
          System.out.println(proj.getSamples().length + "\t" + pedigree.getDnas().length);
          if (!completed.contains(offDNA)) {
            // pedigree.
            String faDNA = null;
            String moDNA = null;
            if (pedigree.getFaDNAIndex(offIndex) > 0
                && pedigree.getiDNA(map.get(pair[0])) == pedigree.getiDNA(pedigree.getIndexOfFaInIDs(offIndex))) {
              faDNA = pedigree.getiDNA(pedigree.getIndexOfFaInIDs(offIndex));
            }
            if (pedigree.getMoDNAIndex(offIndex) > 0
                && pedigree.getiDNA(map.get(pair[0])) == pedigree.getiDNA(pedigree.getIndexOfMoInIDs(offIndex))) {
              moDNA = pedigree.getiDNA(pedigree.getIndexOfMoInIDs(offIndex));
            }
            if (moDNA != null && faDNA != null) {
              throw new IllegalArgumentException("should have been completed");
            }
            if (moDNA != null || faDNA != null) {

              List<PODResults> results = PODAnalysis.analyze(proj, markerDetailSet, segments,
                                                             offDNA, moDNA, faDNA, searchType, 0);

              writer.flush();
              String stat = "NA";
              if (pedigree.getGender(offIndex) == (byte) 1) {
                stat = "MALE";
              } else if (pedigree.getGender(offIndex) == (byte) 2) {
                stat = "FEMALE";
              }
              for (PODResults poResults : results) {
                writer.println(poResults.toString() + "\t" + stat + "\t"
                               + proj.PROJECT_NAME.getValue() + "\t" + completed.contains(offDNA));
              }
            }
          }
        }
      }
    }
    writer.close();

    System.out.println(proj.PEDIGREE_FILENAME.getValue());
    // }
  }

}
