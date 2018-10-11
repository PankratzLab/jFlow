package org.pankratzlab.internal.seq.manage;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import org.genvisis.seq.manage.VCFOps.PLINK_SET_MODE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.gwas.MatchSamples;
import org.pankratzlab.gwas.MatchesVisualized;
import org.pankratzlab.gwas.MergeDatasets;
import org.pankratzlab.shared.filesys.LocusSet;
import org.pankratzlab.shared.filesys.Segment;
import org.pankratzlab.shared.qsub.Qsub;
import org.pankratzlab.utils.gwas.RelationAncestryQc;

/**
 * Class for common actions on a VCF
 */
public class VCFOps {

  public enum UTILITY_TYPE {
    /**
     * Convert a vcf to plink and run gwas QC
     */
    CONVERT_PLINK,

    /**
     * Determine Homogeneity for populations in a vcf
     */
    HOMOGENEITY
  }

  /**
   * @param vcf a vcf file to convert to plink
   * @param rootOut the root output for the plink files
   * @param log
   */

  public static String[] convertToPlinkSet(String outputDir, String vcf, String rootOut,
                                           PLINK_SET_MODE mode, Logger log) {
    String[] plinkCommand = null;
    String dir = outputDir == null ? ext.parseDirectoryOfFile(vcf) + "plink" + ext.rootOf(vcf) + "/"
                                   : outputDir;

    new File(dir).mkdirs();

    rootOut = dir + rootOut;
    String[] outFiles = PSF.Plink.getPlinkBedBimFam(rootOut);
    if (!Files.exists("", outFiles)) {
      plinkCommand = PSF.Plink.getPlinkVCFCommand(vcf, rootOut);
      if (CmdLine.runCommandWithFileChecks(plinkCommand, "", new String[] {vcf}, outFiles, true,
                                           true, false, log)) {

      }
    } else {
      log.reportTimeWarning("Detected that the following files already exist "
                            + ArrayUtils.toStr(outFiles));
    }
    if (Files.exists("", outFiles)) {
      // Hashtable<String, String> newIDS = new Hashtable<String, String>();
      fixFamFile(log, outFiles[2]);
      log.reportTimeInfo("MODE=" + mode);
      if (mode == PLINK_SET_MODE.GWAS_QC) {

        RelationAncestryQc.fullGamut(dir, rootOut, false,
                                     new Logger(dir + "fullGamutOfMarkerAndSampleQC.log"));
        String mdsFile = dir + RelationAncestryQc.GENOME_DIR + "mds20.mds";
        if (Files.exists(mdsFile)) {
          // fixMdsFile(log, dir, newIDS, mdsFile);
          // CmdLine.run("runEigenstratWoHapMap", dir + Qc.GENOME_DIR);
          // CmdLine.run("runEigenstrat2", dir + Qc.GENOME_DIR);
          // fixMdsFile(log, dir + Qc.GENOME_DIR, newIDS, combo_fancy_postnormed_eigens.xln);
        }
      } else if (mode == PLINK_SET_MODE.HOMOGENEITY) {

        if (!Files.exists(ext.parseDirectoryOfFile(outFiles[2]) + "hardy.hwe")) {
          CmdLine.run("plink --bfile " + rootOut
                      + " --maf 0 --geno 1 --mind 1 --hardy --out hardy --noweb",
                      ext.parseDirectoryOfFile(outFiles[2]));
        } else {
          log.reportTimeInfo("Found file " + ext.parseDirectoryOfFile(outFiles[2]) + "hardy.hwe");
        }
      }
    }
    if (!Files.exists(dir + ".qc.pbs")) {
      String gwasQC = ArrayUtils.toStr(PSF.Load.getAllModules(), "\n") + "\n" + Files.getRunString()
                      + " gwas.Qc dir=" + dir;
      Qsub.qsub(dir + "qc.pbs", gwasQC, 62000, 24, 16);
    }
    return outFiles;
  }

  private static Hashtable<String, String> fixFamFile(Logger log, String famFile) {
    Hashtable<String, String> changedIds = new Hashtable<>();
    Files.copyFile(famFile, famFile + ".bak");
    String[][] fam = HashVec.loadFileToStringMatrix(famFile, false, new int[] {0, 1, 2, 3, 4, 5});
    String[][] newfam = new String[fam.length][fam[0].length];
    boolean newSex = false;
    // String[][] fam = HashVec.loadFileToStringMatrix(, false, new int[]{1,2,3,4,5,6},
    // PSF.Regex.GREEDY_WHITESPACE,
    // false, 1000, false);
    int noSexcount = 0;

    for (String[] element : fam) {
      if (element[4].equals("0")) {
        noSexcount++;
      }
    }
    if (noSexcount == fam.length) {
      newSex = true;
      log.reportTimeWarning("Assigning alternating sex specifications");
    }
    Hashtable<String, String> uniqIds = new Hashtable<>();
    boolean swit = true;
    for (int i = 0; i < fam.length; i++) {
      String FidIid = fam[i][0];
      if (FidIid.length() >= 37) {
        String newID = FidIid.substring(0, 36);
        log.reportTimeWarning("Changing " + FidIid + " to " + newID);
        changedIds.put(newID, FidIid);
        // uniqIds.put(FidIid, newID);
        FidIid = newID;
      }
      uniqIds.put(FidIid, FidIid);

      newfam[i][0] = FidIid;
      newfam[i][1] = FidIid;
      if (newSex && swit) {
        newfam[i][4] = "1";
        swit = false;
      } else if (newSex && !swit) {
        newfam[i][4] = "2";
        swit = true;
      } else {
        newfam[i][4] = fam[i][4];
      }
      for (int j = 0; j < newfam[i].length; j++) {
        if (j != 4 && j != 0 && j != 1) {
          newfam[i][j] = fam[i][j];
        }
      }

    }
    if (uniqIds.size() != fam.length) {
      log.reportError("Could not remedy fam file");
    } else {
      log.reportTimeInfo("fixed fam file");
      Files.writeMatrix(newfam, famFile, "\t");
      if (changedIds.size() > 0) {
        try {
          PrintWriter writer = Files.openAppropriateWriter(famFile + ".changedIds");
          for (String newID : changedIds.keySet()) {
            writer.print(newID + "\t" + changedIds.get(newID));
          }
          writer.close();
        } catch (Exception e) {
          log.reportError("Error writing to " + famFile + ".changedIds");
          log.reportException(e);
        }
      }
    }
    return changedIds;
  }

  /**
   * @param vcf runs {@link org.pankratzlab.gwas.RelationAncestryQc#fullGamut(String, boolean)}
   *          after converting to plink* files if neccesary
   * @param log
   */
  public static void vcfGwasQC(String vcf, Logger log) {
    if (Files.exists(vcf)) {
      String dir = ext.parseDirectoryOfFile(vcf);
      String[] plinkFiles = PSF.Plink.getPlinkBedBimFam("plink");
      if (!Files.exists(dir, plinkFiles)) {
        log.reportTimeInfo("Generating plink files for " + vcf + " in " + dir);
        convertToPlinkSet(null, vcf, "plink", PLINK_SET_MODE.GWAS_QC, log);
      }
      log.reportTimeInfo("Running gwas.qc on the following files in " + dir + ":");
      log.reportTimeInfo("\t" + ArrayUtils.toStr(plinkFiles, "\n"));
      // gwas.Qc.fullGamut(dir, false, new Logger(dir + "fullGamutOfMarkerAndSampleQC.log"));
    } else {
      log.reportFileNotFound(vcf);
    }
  }

  /**
   * @param vcf
   * @param fullPathToPopFiles split the vcf by the definitions in this file and run homogeneity
   *          tests
   * @param log
   */
  public static void runHomoGeneity(String vcf, String[] fullPathToPopFiles, Logger log) {
    HashSet<String> toRemoveHash = new HashSet<>();
    String[] toRemove = new String[] {};
    Segment[] toRemoveSeg = new Segment[] {};
    double callRate = 0.80;
    double hwe = .00001;
    int numBarnsPerSample = 5;
    String finalSamples = vcf + ".finalSamples";
    String[] matchUpVpops = Files.listFullPaths(ext.parseDirectoryOfFile(vcf), ".homogeneity.vpop");
    if (matchUpVpops.length < 1) {
      log.reportError("Required file(s) ending with .homogeneity.vpop in directory "
                      + ext.parseDirectoryOfFile(vcf) + " were not found");
      return;
    }
    if (!Files.exists(finalSamples)) {
      log.reportError("Required file " + finalSamples + " is missing");
      return;
    }
    String finalDir = ext.parseDirectoryOfFile(vcf) + "homogeneity/";
    new File(finalDir).mkdirs();
    String idFile = finalDir + "variants_Removed.txt";
    String mdsFile = finalDir + "/genome/mds20.mds";

    if (!Files.exists(idFile) || !Files.exists(mdsFile)) {

      String[] samples = HashVec.loadFileToStringArray(finalSamples, false, new int[] {0}, false);
      log.reportTimeInfo("Found " + samples.length + " samples for the final analysis");
      for (String fullPathToPopFile : fullPathToPopFiles) {
        String[] splits = VcfPopulation.splitVcfByPopulation(vcf, fullPathToPopFile, false, false,
                                                             false, log);
        String[] dirs = new String[splits.length];
        String dir = ext.parseDirectoryOfFile(vcf) + ext.rootOf(fullPathToPopFile) + "/";
        for (int j = 0; j < splits.length; j++) {

          String export = dir + "plink_" + ext.rootOf(splits[j]) + "/";
          dirs[j] = ext.parseDirectoryOfFile(convertToPlinkSet(export, splits[j], "plink",
                                                               PLINK_SET_MODE.HOMOGENEITY, log)[0]);
          if (org.genvisis.seq.manage.VCFOps.getSamplesInFile(splits[j]).length > 50) {
            String callRateFiltered = dir + ext.rootOf(splits[j]) + ".CR." + callRate + ".hwe."
                                      + hwe + ".txt";
            if (!Files.exists(callRateFiltered)) {
              org.genvisis.seq.manage.VCFOps.reportCallRateHWEFiltered(splits[j], callRateFiltered,
                                                                       callRate, hwe, log);
            }
            // LocusSet<Segment> segs =
            // LocusSet.loadSegmentSetFromFile(ext.addToRoot(callRateFiltered, ".segment"), 0, 1, 2,
            // 0, true, true, 0, log);
            // toRemoveSeg = Array.concatAll(toRemoveSeg, segs.getLoci());
            String[] callRateRemove = HashVec.loadFileToStringArray(callRateFiltered, false,
                                                                    new int[] {0}, true);
            for (String element : callRateRemove) {
              toRemoveHash.add(element);
            }
            log.reportTimeInfo(callRateRemove.length + " variants removed from " + splits[j]
                               + " at callrate " + callRateFiltered);
            toRemove = ArrayUtils.unique(ArrayUtils.concatAll(toRemove, callRateRemove));
          }
        }
        String lackOfHomoGeneity = dir + MergeDatasets.CHI_SQUARE_DROPS_FILENAME;
        String problems = dir + "problematic.dat";

        if (!Files.exists(lackOfHomoGeneity)) {
          MergeDatasets.checkForHomogeneity(null, dirs, dir, "ALL", log);

        } else {
          log.reportTimeInfo("Found " + lackOfHomoGeneity
                             + ", assuming this has run to completion");
        }

        String[] lackOfHomoGeneityIDs = HashVec.loadFileToStringArray(lackOfHomoGeneity, false,
                                                                      new int[] {0}, true);
        log.reportTimeInfo(lackOfHomoGeneityIDs.length + " markers lacking homogeneity from "
                           + fullPathToPopFile);
        toRemove = ArrayUtils.unique(ArrayUtils.concatAll(toRemove, lackOfHomoGeneityIDs));
        for (String lackOfHomoGeneityID : lackOfHomoGeneityIDs) {
          toRemoveHash.add(lackOfHomoGeneityID);
        }
        if (Files.exists(problems)) {
          String[] problematic = HashVec.loadFileToStringArray(problems, false, new int[] {0},
                                                               true);
          log.reportTimeInfo(problematic.length + " markers with problems from "
                             + fullPathToPopFile);
          // toRemove = Array.unique(Array.concatAll(toRemove, problematic));
        }
      }

      Files.writeArray(toRemoveHash.toArray(new String[toRemoveHash.size()]), idFile);
      HashSet<String> sampleHash = new HashSet<>();
      for (String sample : samples) {
        sampleHash.add(sample);
      }
      LocusSet<Segment> sort = new LocusSet<Segment>(toRemoveSeg, true, log) {

        /**
        *
        */
        private static final long serialVersionUID = 1L;

      };

      log.reportTimeInfo("Removing " + toRemove.length + " variants from " + vcf);
      log.reportTimeInfo("Removing " + toRemoveSeg.length + " segments as well");
      log.reportTimeInfo("Subsetting to " + sampleHash.size() + " samples");
      String extractVCF = org.genvisis.seq.manage.VCFOps.extractIDs(vcf, idFile, finalDir, true,
                                                                    true, sampleHash,
                                                                    sort.getLoci(), true, false,
                                                                    log);

      if (!Files.exists(mdsFile)) {
        convertToPlinkSet(finalDir, extractVCF, "plink", PLINK_SET_MODE.GWAS_QC, log);
      }
    } else {
      log.reportTimeWarning("found file " + idFile + " and " + mdsFile
                            + " , assuming processing up to this point has been completed");
    }
    org.pankratzlab.utils.widgets.TabVersion.make(mdsFile);
    mdsFile = mdsFile + ".xln";

    for (int i = 1; i < matchUpVpops.length; i++) {
      VcfPopulation vpop = VcfPopulation.load(matchUpVpops[i], POPULATION_TYPE.ANCHOR_BARNACLE,
                                              log);
      vpop.report();
      String matchDir = finalDir + "match_" + ext.rootOf(matchUpVpops[i]) + "/";
      new File(matchDir).mkdirs();

      String factorFile = matchDir + "plink.mds";
      Files.copyFileUsingFileChannels(new File(mdsFile), new File(factorFile), log);
      String[] barnacleIdsPresent = HashVec.loadFileToStringArray(factorFile, true, new int[] {0},
                                                                  true);

      Set<String> anchors = vpop.getSuperPop().get(VcfPopulation.ANCHOR);
      Set<String> barnacles = vpop.getSuperPop().get(VcfPopulation.BARNACLE);

      ArrayList<String> currentBarns = new ArrayList<>();

      for (int j = 0; j < numBarnsPerSample; j++) {
        HashSet<String> barnaclesPresent = new HashSet<>();
        String[] currentBarnsA = currentBarns.toArray(new String[currentBarns.size()]);
        System.out.println("SIZE" + ArrayUtils.unique(currentBarnsA).length);
        for (String barn : barnacles) {
          if (ext.indexOfStr(barn, barnacleIdsPresent) >= 0
              && ext.indexOfStr(barn, currentBarnsA) < 0) {
            barnaclesPresent.add(barn);
          } else {
            // log.reportTimeWarning("Missing sample " + barn + " in file " + factorFile);
          }
        }

        System.out.println(barnaclesPresent.size());
        try {
          Thread.sleep(100);
        } catch (InterruptedException ie) {}
        String anchorList = matchDir + "anchors.txt";
        String barnacleList = matchDir + j + "barnacles.txt";

        Files.writeArray(anchors.toArray(new String[anchors.size()]), anchorList);
        Files.writeArray(barnaclesPresent.toArray(new String[barnaclesPresent.size()]),
                         barnacleList);
        String[] run = new String[] {"C1", "C3"};
        System.out.println("RUNNING match1");
        String matchFile = MatchSamples.matchMaker(matchDir, ext.removeDirectoryInfo(anchorList),
                                                   ext.removeDirectoryInfo(barnacleList),
                                                   ext.removeDirectoryInfo(factorFile), run,
                                                   new double[] {1, 1}, false);
        System.out.println("RUNNING match2");

        matchFile = MatchSamples.normalizeDistances(matchDir, matchFile, 0, 100);
        System.out.println("RUNNING match3");

        String pairs = matchDir + MatchSamples.matchPairs(matchDir, matchFile, true);
        System.out.println("RUNNING match4");

        Files.copyFileUsingFileChannels(new File(pairs), new File(pairs + j + ".selection"), log);

        String[] barnesPicked = HashVec.loadFileToStringArray(pairs, true, new int[] {1}, true);
        String[] deletes = Files.listFullPaths(matchDir, ".xln");
        new MatchesVisualized(matchDir, ext.removeDirectoryInfo(anchorList),
                              ext.removeDirectoryInfo(barnacleList),
                              ext.removeDirectoryInfo(factorFile),
                              ext.indexFactors(run, Files.getHeaderOfFile(factorFile, log), true),
                              ext.removeDirectoryInfo(pairs));

        for (String delete : deletes) {
          new File(delete).delete();
        }
        System.out.println(ArrayUtils.toStr(barnesPicked));
        for (String element : barnesPicked) {
          currentBarns.add(element);
        }

      }
      String finalVpop = matchDir + "barnacle.vpop";
      try {
        PrintWriter writer = Files.openAppropriateWriter(finalVpop);
        writer.println(ArrayUtils.toStr(VcfPopulation.HEADER));
        for (int j = 0; j < currentBarns.size(); j++) {
          writer.println(currentBarns.get(j) + "\t" + VcfPopulation.CONTROL + "\t"
                         + VcfPopulation.CONTROL);
        }
        for (String anchor : anchors) {
          writer.println(anchor + "\t" + VcfPopulation.CASE + "\t" + VcfPopulation.CASE);
        }
        writer.close();
      } catch (Exception e) {
        log.reportError("Error writing to " + finalVpop);
        log.reportException(e);
      }
    }

  }

  private static final String VCF_COMMAND = "vcf=";
  private static final String UTILITY_COMMAND = "utility=";
  public static final String COMMAND_VCF_OPS_EXTRACT = "vcfExtract";
  public static final String COMMAND_VCF_EXTRACT_DESCRIPTION = "extract a file of segments from a vcf and bam files";

  public static void main(String[] args) {
    int numArgs = args.length;
    String vcf = "Avcf.vcf";
    String populationFile = null;
    String logfile = null;
    UTILITY_TYPE type = UTILITY_TYPE.HOMOGENEITY;

    String usage = "\n" + "seq.analysis.VCFUtils requires 0-1 arguments\n";
    usage += "   (1) full path to a vcf file (i.e. " + VCF_COMMAND + vcf + " (default))\n" + "";
    usage += "   (2) utility type (i.e. " + UTILITY_COMMAND + type + " (default))\n" + "";
    usage += "   (3) full path to a file (can be comma delimited for homogeneity utility) defining a population for the vcf (i.e. vpopFile= (no default))\n"
             + "";
    usage += "   (4) full path to a file name with chr,start,stop or *.bim to extract (i.e. segs= (no default))\n"
             + "";

    usage += "   NOTE: available utilities are:\n";

    for (int i = 0; i < UTILITY_TYPE.values().length; i++) {
      usage += UTILITY_TYPE.values()[i] + "\n";
    }
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith(VCF_COMMAND)) {
        vcf = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith(UTILITY_COMMAND)) {
        type = UTILITY_TYPE.valueOf(ext.parseStringArg(arg, ""));
        numArgs--;
      } else if (arg.startsWith("vpopFile=")) {
        populationFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
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
      Logger log = new Logger(logfile);
      log.reportTimeInfo("Running utiltity type: " + type);

      switch (type) {
        case CONVERT_PLINK:
          convertToPlinkSet(null, vcf, "plink", PLINK_SET_MODE.GWAS_QC, log);
        case HOMOGENEITY:
          runHomoGeneity(vcf, populationFile.split(","), log);
          break;
        default:
          System.err.println("Invalid utility type: Available are ->");
          for (int i = 0; i < UTILITY_TYPE.values().length; i++) {
            usage += UTILITY_TYPE.values()[i] + "\n";
          }
          break;
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
