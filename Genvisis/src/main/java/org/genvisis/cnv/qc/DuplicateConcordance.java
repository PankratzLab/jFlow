package org.genvisis.cnv.qc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class DuplicateConcordance {

  private final double projectConcordance;
  private final double projectNonMissingConcordance;
  private final int markersChecked;
  private final int duplicatePairsChecked;


  /**
   * @param discordantCalls number of marker calls that did not match
   * @param nonMissingDiscordantCalls number of marker calls that did not match, excluding cases
   *        where one call was missing
   * @param markersChecked number of markers checked
   * @param duplicatePairsChecked number of duplicate pairs checked
   */
  private DuplicateConcordance(int discordantCalls, int nonMissingDiscordantCalls,
                               int markersChecked, int duplicatePairsChecked) {
    super();
    int totalChecks = markersChecked * duplicatePairsChecked;
    projectConcordance = (double) (totalChecks - discordantCalls) / totalChecks;
    projectNonMissingConcordance = (double) (totalChecks - nonMissingDiscordantCalls) / totalChecks;
    this.markersChecked = markersChecked;
    this.duplicatePairsChecked = duplicatePairsChecked;
  }

  public double getProjectNonMissingConcordance() {
    return projectNonMissingConcordance;
  }

  public double getProjectConcordance() {
    return projectConcordance;
  }

  public int getMarkersChecked() {
    return markersChecked;
  }

  public int getDuplicatePairsChecked() {
    return duplicatePairsChecked;
  }

  public String getConcordanceString() {
    return "Duplicate Concordance was calculated to be " + projectConcordance
           + " (including missing calls) and " + projectNonMissingConcordance
           + " (excluding missing calls) using " + duplicatePairsChecked
           + " pairs of duplicates at " + markersChecked + " markers.";
  }

  /**
   * 
   * @param proj Project to calculate duplicate concordance for
   * @param targetMarkers Markers to use in concordance checks or null to check all markers
   * @return
   */
  public static DuplicateConcordance calculateDuplicateConcordances(Project proj,
                                                                    String[] targetMarkers) {
    Logger log = proj.getLog();
    ClusterFilterCollection clusterFilterCollection = proj.getClusterFilterCollection();
    String[] markerNames;
    int[] markerIndices;
    if (targetMarkers == null) {
      markerNames = proj.getMarkerNames();
      markerIndices = null;
    } else {
      markerNames = targetMarkers;
      markerIndices =
          ext.indexLargeFactors(markerNames, proj.getMarkerNames(), true, log, false, false);
      for (int i = 0; i < markerIndices.length; i++) {
        if (markerIndices[i] == -1) {
          log.reportTimeError("Marker " + markerNames[i] + " could not be found in project");
          return null;
        }
      }
    }
    String sampleData = proj.SAMPLE_DATA_FILENAME.getValue();
    if (sampleData == null) {
      log.reportTimeError("Project Sample Data file is not defined, cannot determine duplicates");
      return null;
    }
    if (!Files.exists(sampleData)) {
      log.reportTimeError("Project Sample Data file, " + sampleData
                          + " could not be found, cannot determine duplicates");
      return null;
    }
    String[] sampleDataHeader = Files.getHeaderOfFile(sampleData, log);
    String[] sampleDataCols = new String[] {"DNA", "CLASS=Exclude", "DuplicateId"};
    int[] sampleDataIndices =
        ext.indexFactors(sampleDataCols, sampleDataHeader, false, log, false, false);
    for (int i = 0; i < sampleDataIndices.length; i++) {
      if (sampleDataIndices[i] == -1) {
        log.reportTimeError("Could not find " + sampleDataCols[i]
                            + " in Sample Data file, cannot determine duplicates");
        return null;
      }
    }
    String[][] sampleInfo = HashVec.loadFileToStringMatrix(sampleData, true, sampleDataIndices,
                                                           proj.JAR_STATUS.getValue());
    HashVec.loadFileToStringArray(sampleData, true, sampleDataIndices, false);
    HashMap<String, HashSet<String>> duplicateSets = new HashMap<String, HashSet<String>>();
    for (String[] sampleLine : sampleInfo) {
      String dna = sampleLine[0];
      boolean exclude = sampleLine[1].equals("1");
      String duplicateID = sampleLine[2];

      if (!duplicateID.equals(".") && !exclude) {
        HashSet<String> duplicateSet = duplicateSets.get(duplicateID);
        if (duplicateSet == null) {
          duplicateSet = new HashSet<String>();
          duplicateSets.put(duplicateID, duplicateSet);
        }
        duplicateSet.add(dna);
      }
    }
    int discordantCalls = 0;
    int nonMissingDiscordantCalls = 0;
    int pairsChecked = 0;

    for (HashSet<String> duplicateSet : duplicateSets.values()) {
      HashSet<String> loopDuplicateSet = new HashSet<String>(duplicateSet);
      for (String dna1 : loopDuplicateSet) {
        duplicateSet.remove(dna1);
        if (!duplicateSet.isEmpty()) {
          Sample sample1 = proj.getFullSampleFromRandomAccessFile(dna1);
          if (sample1 == null) {
            log.reportTimeError("Could not find data for Sample " + dna1
                                + ", will not be used to calculate concordance");
            continue;
          }
          for (String dna2 : duplicateSet) {
            Sample sample2 = proj.getFullSampleFromRandomAccessFile(dna2);
            if (sample2 == null) {
              log.reportTimeError("Could not find data for Sample " + dna2
                                  + ", will not be used to calculate concordance");
              continue;
            }
            pairsChecked++;
            byte[] s1Genotypes, s2Genotypes;
            if (clusterFilterCollection == null) {
              s1Genotypes = sample1.getAB_Genotypes(markerIndices);
              s2Genotypes = sample2.getAB_Genotypes(markerIndices);
            } else {
              s1Genotypes = sample1.getAB_GenotypesAfterFilters(markerNames, markerIndices,
                                                                clusterFilterCollection, 0.0f);
              s2Genotypes = sample2.getAB_GenotypesAfterFilters(markerNames, markerIndices,
                                                                clusterFilterCollection, 0.0f);
            }
            for (int i = 0; i < s1Genotypes.length; i++) {
              if (s1Genotypes[i] != s2Genotypes[i]) {
                discordantCalls++;
                if (s1Genotypes[i] != -1 && s2Genotypes[i] != -1) {
                  nonMissingDiscordantCalls++;
                }
              }
            }
          }
        }
      }
    }

    if (pairsChecked == 0) {
      log.reportTimeError("No duplicates could be compared, duplicate concordance cannot be calculated");
      return null;
    }

    return new DuplicateConcordance(discordantCalls, nonMissingDiscordantCalls, markerNames.length,
                                    pairsChecked);

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    Project proj = null;
    String markerKeeps = null;
    String markerDrops = null;

    String usage = "\n" + "cnv.qc.DuplicateConcordance requires 1-2 arguments\n"
                   + "   (1) Project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false)
                   + " (not the default))\n" + "AND\n"
                   + "   (2) File of markers to use (i.e. markerKeeps=keeps.txt (not the default))\n"
                   + "OR\n"
                   + "   (2) File of markers to not use (i.e. markerDrops=drops.txt (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        proj = new Project(arg.split("=")[1], false);
        numArgs--;
      } else if (arg.startsWith("markerKeeps=")) {
        markerKeeps = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("markerDrops=")) {
        markerDrops = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (proj == null) {
        System.err.println("Project must be defined");
        System.err.println(usage);
        System.exit(1);
      }
      if (markerKeeps != null && markerDrops != null) {
        System.err.println("Include a marker keeps or drops file but not both");
        System.err.println(usage);
        System.exit(1);
      }
      String[] targetMarkers;
      if (markerKeeps != null) {
        targetMarkers = proj.getTargetMarkers(markerKeeps);
      } else if (markerDrops != null) {
        Set<String> excludes = HashVec.loadFileToHashSet(markerDrops, false);
        ArrayList<String> markers = new ArrayList<String>();
        for (String marker : proj.getMarkerNames()) {
          if (!excludes.contains(marker)) {
            markers.add(marker);
          }
        }
        targetMarkers = Array.toStringArray(markers);
      } else {
        targetMarkers = null;
      }

      DuplicateConcordance duplicateConcordance =
          calculateDuplicateConcordances(proj, targetMarkers);
      if (duplicateConcordance != null) {
        proj.getLog().report(duplicateConcordance.getConcordanceString());
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
