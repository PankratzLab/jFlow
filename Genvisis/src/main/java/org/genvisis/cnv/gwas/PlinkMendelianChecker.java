package org.genvisis.cnv.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.cnv.qc.MendelErrors.MendelErrorCheck;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.var.SampleData;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.utils.gwas.Plink;
import com.google.common.collect.ImmutableSet;

public class PlinkMendelianChecker {

  public static final String TRIOS_FILENAME = "trios.xln";
  public static final String REL_CHECKS_FILENAME = "relationshipChecks.xln";
  public static final Set<String> OUTPUTS = ImmutableSet.of(TRIOS_FILENAME, REL_CHECKS_FILENAME);

  static class Pair {

    String fid1;
    String fid2;
    String iid1;
    String iid2;

    public Pair(String fid1, String iid1, String fid2, String iid2) {
      this.fid1 = fid1;
      this.iid1 = iid1;
      this.fid2 = fid2;
      this.iid2 = iid2;
    }

    public Pair(String fid, String iid1, String iid2) {
      fid1 = fid;
      this.iid1 = iid1;
      fid2 = fid;
      this.iid2 = iid2;
    }

  }

  static class GenomeLoader {

    private static final int FID1_INDEX = 0;
    private static final int IID1_INDEX = 1;
    private static final int FID2_INDEX = 2;
    private static final int IID2_INDEX = 3;
    private static final int RT_INDEX = 4;
    private static final int EZ_INDEX = 5;
    private static final int Z0_INDEX = 6;
    private static final int Z1_INDEX = 7;
    private static final int Z2_INDEX = 8;
    private static final int PI_HAT_INDEX = 9;
    private static final int PHE_INDEX = 10;
    private static final int DST_INDEX = 11;
    private static final int PPC_INDEX = 12;
    private static final int RATIO_INDEX = 13;

    HashMap<String, HashMap<String, String>> pairData;
    ArrayList<String> unrelLines;

    private GenomeLoader() {
      pairData = new HashMap<>();
      unrelLines = new ArrayList<>();
    }

    // TODO only goes one direction?
    static GenomeLoader run(String genomeFile, List<Pair> pairs) {
      BufferedReader reader;
      String line;

      GenomeLoader gl = new GenomeLoader();

      HashSet<String> famSet = new HashSet<>();

      HashMap<String, HashSet<String>> pairSets = new HashMap<>();

      for (Pair p : pairs) {
        famSet.add(p.fid1);
        // famSet.add(p.fid2); // uncomment for bi-directional
        HashSet<String> indiPairs = pairSets.get(p.fid1 + "\t" + p.iid1);
        if (indiPairs == null) {
          indiPairs = new HashSet<>();
          pairSets.put(p.fid1 + "\t" + p.iid1, indiPairs);
        }
        indiPairs.add(p.fid2 + "\t" + p.iid2);
        // not set up for bi-directional
      }

      // temp[1] // FID1
      // temp[2] // IID1
      // temp[3] // FID2
      // temp[4] // IID2
      // temp[5] // RT
      // temp[6] // EZ
      // temp[7] // Z0
      // temp[8] // Z1
      // temp[9] // Z2
      // temp[10] // PI_HAT
      // temp[11] // PHE
      // temp[12] // DST
      // temp[13] // PPC
      // temp[14] // RATIO

      try {
        reader = Files.getAppropriateReader(genomeFile);
        line = reader.readLine(); // read header
        String fid1, iid1, fid2, iid2;
        while ((line = reader.readLine()) != null) {
          line = line.trim();
          // String[] temp = line.split(PSF.Regex.GREEDY_WHITESPACE); // slow for 32mil lines
          int breakInd = line.indexOf(" ");
          fid1 = line.substring(0, breakInd);
          if (!famSet.contains(fid1)) {
            continue; // comment out for bi-directional
          }
          while (line.charAt(breakInd) == ' ') {
            breakInd++;
          }
          int break2Ind = line.indexOf(" ", breakInd);
          iid1 = line.substring(breakInd, break2Ind);
          if (!pairSets.containsKey(fid1 + "\t" + iid1)) {
            continue; // comment out for bi-directional
          }
          while (line.charAt(break2Ind) == ' ') {
            break2Ind++;
          }
          breakInd = break2Ind;
          break2Ind = line.indexOf(" ", breakInd);
          fid2 = line.substring(breakInd, break2Ind);
          // if (!famSet.contains(fid1) && !famSet.contains(fid2)) continue; // uncomment for
          // bi-directional
          while (line.charAt(break2Ind) == ' ') {
            break2Ind++;
          }
          breakInd = break2Ind;
          break2Ind = line.indexOf(" ", breakInd);
          iid2 = line.substring(breakInd, break2Ind);
          // if (!pairSets.containsKey(fid1 + "\t" + iid1) && !pairSets.containsKey(fid2 + "\t" +
          // iid2)) continue; // uncomment for bi-directional

          if (pairSets.get(fid1 + "\t" + iid1).contains(fid2 + "\t" + iid2)) {
            HashMap<String, String> indiMap = gl.pairData.get(fid1 + "\t" + iid1);
            if (indiMap == null) {
              indiMap = new HashMap<>();
              gl.pairData.put(fid1 + "\t" + iid1, indiMap);
            }
            indiMap.put(fid2 + "\t" + iid2, line);
          } else {
            while (line.charAt(break2Ind) == ' ') {
              break2Ind++;
            }
            breakInd = break2Ind;
            break2Ind = line.indexOf(" ", breakInd);
            // RT
            while (line.charAt(break2Ind) == ' ') {
              break2Ind++;
            }
            breakInd = break2Ind;
            break2Ind = line.indexOf(" ", breakInd);
            // EZ
            while (line.charAt(break2Ind) == ' ') {
              break2Ind++;
            }
            breakInd = break2Ind;
            break2Ind = line.indexOf(" ", breakInd);
            // Z0
            while (line.charAt(break2Ind) == ' ') {
              break2Ind++;
            }
            breakInd = break2Ind;
            break2Ind = line.indexOf(" ", breakInd);
            // Z1
            while (line.charAt(break2Ind) == ' ') {
              break2Ind++;
            }
            breakInd = break2Ind;
            break2Ind = line.indexOf(" ", breakInd);
            // Z2
            while (line.charAt(break2Ind) == ' ') {
              break2Ind++;
            }
            breakInd = break2Ind;
            break2Ind = line.indexOf(" ", breakInd);
            String piHatStr = line.substring(breakInd, break2Ind);
            if (Double.parseDouble(piHatStr) > 0.2) {
              gl.unrelLines.add(line);
            }
          }
          // not set up for bi-directional
        }
        reader.close();
        reader = null;
      } catch (IOException e) {
        e.printStackTrace();
      }

      return gl;
    }

  }

  static class MendelLoader {

    HashMap<String, ArrayList<String>> errorMarkersMapFather;
    HashMap<String, ArrayList<String>> errorMarkersMapMother;

    private MendelLoader() {
      errorMarkersMapFather = new HashMap<>();
      errorMarkersMapMother = new HashMap<>();
    }

    static MendelLoader run(String mendelFile) {
      MendelLoader ml = new MendelLoader();
      BufferedReader reader;
      String line;
      String[] temp;
      try {
        reader = Files.getAppropriateReader(mendelFile);
        line = reader.readLine().trim();
        if (ext.checkHeader(line.split(PSF.Regex.GREEDY_WHITESPACE),
                            new String[] {"FID", "KID", "CHR", "SNP", "CODE", "ERROR"}, false)) {
          while ((line = reader.readLine()) != null) {
            line = line.trim();
            temp = line.split(PSF.Regex.GREEDY_WHITESPACE);

            // FID
            // KID
            // CHR
            // SNP
            // CODE
            // ERROR (could be split into 5 when from PLINK due to spaces)

            try {
              MendelErrorCheck mec = new MendelErrorCheck(Integer.parseInt(temp[4]));

              String fidiid = temp[0] + "\t" + temp[1];

              if (mec.hasFaMendelError()) {
                ArrayList<String> mkrs = ml.errorMarkersMapFather.get(fidiid);
                if (mkrs == null) {
                  mkrs = new ArrayList<>();
                  ml.errorMarkersMapFather.put(fidiid, mkrs);
                }
                mkrs.add(temp[3]);
              }
              if (mec.hasMoMendelError()) {
                ArrayList<String> mkrs = ml.errorMarkersMapMother.get(fidiid);
                if (mkrs == null) {
                  mkrs = new ArrayList<>();
                  ml.errorMarkersMapMother.put(fidiid, mkrs);
                }
                mkrs.add(temp[3]);
              }
            } catch (NumberFormatException e) {
              System.err.println("Warning - Non-integer mendelian error code (" + temp[4]
                                 + ") ignored!");
            }

          }
        }
        reader.close();
        reader = null;
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      return ml;
    }

  }

  final Project project;
  final Pedigree ped;
  final String genomeFile; // from property - plink.genome
  final boolean genomeDNA; // is the genome file FID/IID or DNA/DNA?
  final String mendelFile; // from MarkerMetrics.fullQC - outputs a property file with an appended
                          // name
  final String outDir;
  final Logger log;

  public PlinkMendelianChecker(Project project) {
    this.project = project;
    ped = new Pedigree(project, project.PEDIGREE_FILENAME.getValue(false, false));
    mendelFile = ext.rootOf(project.MARKER_METRICS_FILENAME.getValue(true, false), false)
                 + MarkerMetrics.DEFAULT_MENDEL_FILE_SUFFIX;
    genomeFile = project.GENOME_CLUSTER_FILENAME.getValue();
    log = project.getLog();
    PlinkData.ExportIDScheme idScheme = PlinkData.detectExportIDScheme(project, genomeFile);
    switch (idScheme) {
      case DNA_DNA:
        genomeDNA = true;
        break;
      case FID_IID:
        genomeDNA = false;
        break;
      default:
        log.reportError(this.getClass().getSimpleName() + " does not support "
                        + PlinkData.ExportIDScheme.class.getSimpleName() + idScheme.toString());
        throw new IllegalArgumentException("Invalid " + PlinkData.ExportIDScheme.class.getName());
    }
    outDir = parseOutputDirectory(project);

  }

  public PlinkMendelianChecker(String pedFile, String mendelFile, String genomeFile,
                               boolean genomeDNA, String outDir, Logger log) {
    project = null;
    this.ped = new Pedigree(project, pedFile);
    this.mendelFile = mendelFile;
    this.genomeFile = genomeFile;
    this.genomeDNA = genomeDNA;
    this.outDir = ext.verifyDirFormat(outDir);
    this.log = log;
  }

  public void run() {
    SampleData sampleData = null;
    SampleQC sampQC = null;
    String[] samples = null;
    HashSet<String> pedDNA;
    HashMap<String, String> dnaLookup, idLookup;
    HashMap<String, String[]> pedToFAMO;
    HashMap<String, ArrayList<String>> childrenMap;
    HashMap<String, Integer> qcIndexMap = null;
    PrintWriter writer;
    GenomeLoader gl = null;
    MendelLoader ml = null;

    if (project != null) {
      System.out.println(ext.getTime() + "]\tLoading Project data...");
      samples = project.getSamples();
      sampleData = project.getSampleData(false);
      sampQC = SampleQC.loadSampleQC(project);
      qcIndexMap = new HashMap<>();
      if (sampQC != null) {
        for (int i = 0; i < sampQC.getSamples().length; i++) {
          qcIndexMap.put(sampQC.getSamples()[i], i);
        }
      }
    }

    pedDNA = new HashSet<>();
    pedToFAMO = new HashMap<>();
    childrenMap = new HashMap<>();
    dnaLookup = new HashMap<>();
    idLookup = new HashMap<>();

    for (int i = 0; i < ped.getIDs().length; i++) {
      // if (pe == null) {
      // continue;
      // }
      if (ext.isMissingValue(ped.getDnas()[i])) {
        continue;
      }
      if (sampleData != null) {
        if (sampleData.individualShouldBeExcluded(ped.getDnas()[i])) {
          continue;
        }
      }
      dnaLookup.put(ped.getFID(i) + "\t" + ped.getIID(i), ped.getDnas()[i]);
      idLookup.put(ped.getDnas()[i], ped.getFID(i) + "\t" + ped.getIID(i));
      pedDNA.add(ped.getDnas()[i]);
      pedToFAMO.put(ped.getFID(i) + "\t" + ped.getIID(i),
                    new String[] {"0".equals(ped.getFA(i)) ? "."
                                                           : ped.getFID(i) + "\t" + ped.getFA(i),
                                  "0".equals(ped.getMO(i)) ? "."
                                                           : ped.getFID(i) + "\t" + ped.getMO(i)});
      if (!"0".equals(ped.getFA(i))) {
        ArrayList<String> children = childrenMap.get(ped.getFID(i) + "\t" + ped.getFA(i));
        if (children == null) {
          children = new ArrayList<>();
          childrenMap.put(ped.getFID(i) + "\t" + ped.getFA(i), children);
        }
        children.add(ped.getFID(i) + "\t" + ped.getIID(i));
      }
      if (!"0".equals(ped.getMO(i))) {
        ArrayList<String> children = childrenMap.get(ped.getFID(i) + "\t" + ped.getMO(i));
        if (children == null) {
          children = new ArrayList<>();
          childrenMap.put(ped.getFID(i) + "\t" + ped.getMO(i), children);
        }
        children.add(ped.getFID(i) + "\t" + ped.getIID(i));
      }
    }

    ArrayList<Pair> pairs = new ArrayList<>();
    for (java.util.Map.Entry<String, ArrayList<String>> childrenList : childrenMap.entrySet()) {
      String fidiid = childrenList.getKey();

      for (String childFIDIID : childrenList.getValue()) {
        String[] spl1 = childFIDIID.split("\t");
        String[] spl2 = fidiid.split("\t");
        String childDNA = dnaLookup.get(childFIDIID); // shouldn't be nul, unless child DNA value
                                                     // was missing
        String parentDNA = dnaLookup.get(fidiid);
        if (childDNA == null) {
          continue;
        }
        if (parentDNA != null) {
          if (genomeDNA) {
            pairs.add(new Pair(childDNA, childDNA, parentDNA, parentDNA)); // for DNA/DNA genome
                                                                          // file
            pairs.add(new Pair(parentDNA, parentDNA, childDNA, childDNA)); // either bi-directional
                                                                          // here or in genome
                                                                          // loader
          } else {
            pairs.add(new Pair(spl1[0], spl1[1], spl2[0], spl2[1])); // for fid/iid genome file
            pairs.add(new Pair(spl2[0], spl2[1], spl1[0], spl1[1])); // either bi-directional here
                                                                    // or in genomeloader
          }
        }

        String[] famo = pedToFAMO.get(childFIDIID);
        ArrayList<String> spouseChildren = new ArrayList<>();
        if (famo[0].equals(fidiid) && !".".equals(famo[1])) {
          spouseChildren = childrenMap.get(famo[1]);
        } else if (famo[1].equals(fidiid) && !".".equals(famo[0])) {
          spouseChildren = childrenMap.get(famo[0]);
        }

        HashSet<String> sameParentSibs = new HashSet<>(childrenList.getValue());
        for (String otherChild : spouseChildren) {
          if (otherChild.equals(childFIDIID)) {
            continue;
          }
          spl2 = otherChild.split("\t");
          if (sameParentSibs.contains(otherChild)) {
            sameParentSibs.remove(otherChild);
          }
          String otherChildDNA = dnaLookup.get(otherChild);
          if (otherChildDNA != null) {
            if (genomeDNA) {
              pairs.add(new Pair(childDNA, childDNA, otherChildDNA, otherChildDNA)); // sib1 -> sib2
              pairs.add(new Pair(otherChildDNA, otherChildDNA, childDNA, childDNA)); // sib2 -> sib1
            } else {
              pairs.add(new Pair(spl1[0], spl1[1], spl2[0], spl2[1])); // sib1 -> sib2
              pairs.add(new Pair(spl2[0], spl2[1], spl1[0], spl1[1])); // sib2 -> sib1
            }
          }
        }
        for (String halfSib : sameParentSibs) {
          if (halfSib.equals(childFIDIID)) {
            continue;
          }
          String halfSibDNA = dnaLookup.get(halfSib);
          if (halfSibDNA != null) {
            if (genomeDNA) {
              pairs.add(new Pair(childDNA, childDNA, halfSibDNA, halfSibDNA));
              pairs.add(new Pair(halfSibDNA, halfSibDNA, childDNA, childDNA));
            } else {
              spl2 = halfSib.split("\t");
              pairs.add(new Pair(spl1[0], spl1[1], spl2[0], spl2[1])); // sib1 -> sib2
              pairs.add(new Pair(spl2[0], spl2[1], spl1[0], spl1[1])); // sib2 -> sib1
            }
          }
        }
      }
    }

    if (genomeFile != null && (new File(genomeFile)).exists()) {
      System.out.println(ext.getTime() + "]\tLoading GenomicCluster data...");
      gl = GenomeLoader.run(genomeFile, pairs);
    } else {
      log.reportTimeWarning("No Genome data found - output will be limited. Is "
                            + project.GENOME_CLUSTER_FILENAME.getName() + " set appropriately?");
    }
    if (mendelFile != null && (new File(mendelFile)).exists()) {
      System.out.println(ext.getTime() + "]\tLoading MendelianError data...");
      ml = MendelLoader.run(mendelFile);
    } else {
      log.reportTimeWarning("Mendelian data not found"
                            + (mendelFile == null ? "" : " at " + mendelFile)
                            + " - output will be limited");
    }

    StringBuilder sb;
    sb = new StringBuilder().append(org.pankratzlab.common.PSF.Plink.FAMID).append("\t")
                            .append(org.pankratzlab.common.PSF.Plink.INDID).append("\t")
                            .append(org.pankratzlab.common.PSF.Plink.FATHER_ID).append("\t")
                            .append(org.pankratzlab.common.PSF.Plink.MOTHER_ID).append("\t")
                            .append("SEX").append("\t")
                            .append(org.pankratzlab.common.PSF.Plink.IND_DNA).append("\t")
                            .append(org.pankratzlab.common.PSF.Plink.FATHER_DNA).append("\t")
                            .append(org.pankratzlab.common.PSF.Plink.MOTHER_DNA).append("\t");
    if (gl != null) {
      sb.append("IBD0_FATHER").append("\t").append("IBD1_FATHER").append("\t").append("IBD2_FATHER")
        .append("\t").append("PIHAT_FATHER").append("\t").append("IBD0_MOTHER").append("\t")
        .append("IBD1_MOTHER").append("\t").append("IBD2_MOTHER").append("\t")
        .append("PIHAT_MOTHER").append("\t");
    }
    if (ml != null) {
      sb.append("MendelianErrorsFather").append("\t").append("MendelianErrorsMother").append("\t");
    }

    if (sampQC != null) {
      sb.append("LRR_SD_CHILD").append("\t").append("LRR_SD_FATHER").append("\t")
        .append("LRR_SD_MOTHER").append("\t").append("CALLRATE_CHILD").append("\t")
        .append("CALLRATE_FATHER").append("\t").append("CALLRATE_MOTHER").append("\t");
    }

    if (sampleData != null) {
      sb.append("EXCLUDE_CHILD").append("\t").append("EXCLUDE_FATHER").append("\t")
        .append("EXCLUDE_MOTHER").append("\t");
    }

    sb.append("COMPLETE_TRIO").append("\t");

    System.out.println(ext.getTime() + "]\tWriting result data to " + outDir + TRIOS_FILENAME);

    int missingCount = 0;

    String outputHeader = sb.toString();
    writer = Files.getAppropriateWriter(outDir + TRIOS_FILENAME);
    writer.println(outputHeader);
    for (int i = 0; i < ped.getIDs().length; i++) {
      // PedigreeEntry pe = ped.getPedigreeEntries()[i];
      // if (pe == null) continue;
      // if (null == ped.getMO() || "0".equals(ped.getMO()) || null == ped.getFA() ||
      // "0".equals(ped.getFA())) continue;
      sb = new StringBuilder();
      sb.append(ped.getFID(i)).append("\t").append(ped.getIID(i)).append("\t").append(ped.getFA(i))
        .append("\t").append(ped.getMO(i)).append("\t").append(ped.getGender(i)).append("\t")
        .append(ped.getiDNA(i)).append("\t");

      boolean missingChildDNA = ext.isMissingValue(ped.getiDNA(i));
      boolean missingFADNA = false;
      boolean missingMODNA = false;
      String faDNA;
      if (ped.getFaDNAIndex(i) == Pedigree.MISSING_DNA_INDEX || samples == null) {
        if ("0".equals(ped.getFA(i))
            || dnaLookup.get(ped.getFID(i) + "\t" + ped.getFA(i)) == null) {
          missingFADNA = true;
          faDNA = ".";
        } else {
          faDNA = dnaLookup.get(ped.getFID(i) + "\t" + ped.getFA(i));
        }
      } else {
        faDNA = samples[ped.getFaDNAIndex(i)];
      }
      sb.append(faDNA).append("\t");
      String moDNA;
      if (ped.getMoDNAIndex(i) == -1 || samples == null) {
        if ("0".equals(ped.getMO(i))
            || dnaLookup.get(ped.getFID(i) + "\t" + ped.getMO(i)) == null) {
          missingMODNA = true;
          moDNA = ".";
        } else {
          moDNA = dnaLookup.get(ped.getFID(i) + "\t" + ped.getMO(i));
        }
      } else {
        moDNA = samples[ped.getMoDNAIndex(i)];
      }
      sb.append(moDNA).append("\t");
      if (missingFADNA && missingMODNA) {
        // Founder; remove
        continue;
      }
      boolean missingDNA = missingChildDNA || missingFADNA || missingMODNA;

      if (gl != null) {

        HashMap<String, String> genoLines = gl.pairData.get(getGenomeFIDIID(i));

        String key = getGenomeFatherFIDIID(i);
        if (genoLines == null || (!".".equals(faDNA) && !genoLines.containsKey(key))) {
          genoLines = gl.pairData.get(key);
          key = getGenomeFIDIID(i);
        }

        if (genoLines == null || ".".equals(faDNA) || !genoLines.containsKey(key)) {
          missingCount++;
          sb.append(".").append("\t").append(".").append("\t").append(".").append("\t").append(".")
            .append("\t");
        } else {
          String[] tmpGL = genoLines.get(key).split(PSF.Regex.GREEDY_WHITESPACE);
          sb.append(tmpGL[6]).append("\t").append(tmpGL[7]).append("\t").append(tmpGL[8])
            .append("\t").append(tmpGL[9]).append("\t");
        }

        genoLines = gl.pairData.get(getGenomeFIDIID(i));

        key = getGenomeMotherFIDIID(i);
        if (genoLines == null || (!".".equals(moDNA) && !genoLines.containsKey(key))) {
          genoLines = gl.pairData.get(key);
          key = getGenomeFIDIID(i);
        }

        if (genoLines == null || ".".equals(moDNA) || !genoLines.containsKey(key)) {
          missingCount++;
          sb.append(".").append("\t").append(".").append("\t").append(".").append("\t").append(".")
            .append("\t");
        } else {
          String[] tmpGL = genoLines.get(key).split(PSF.Regex.GREEDY_WHITESPACE);
          sb.append(tmpGL[6]).append("\t").append(tmpGL[7]).append("\t").append(tmpGL[8])
            .append("\t").append(tmpGL[9]).append("\t");
        }
      }

      if (ml != null) {
        ArrayList<String> errors = ml.errorMarkersMapFather.get(idLookup.get(ped.getiDNA(i)));
        if (errors == null) {
          sb.append(0).append("\t");
        } else {
          sb.append(errors.size()).append("\t");
        }
        errors = ml.errorMarkersMapMother.get(idLookup.get(ped.getiDNA(i)));
        if (errors == null) {
          sb.append(0).append("\t");
        } else {
          sb.append(errors.size()).append("\t");
        }
      }

      boolean highLRRSD = false;
      boolean lowCallrate = false;
      if (sampQC != null) {
        double[] data = sampQC.getDataFor("LRR_SD");
        Integer indexInt = qcIndexMap.get(ped.getiDNA(i));
        Integer faIndex = null;
        if (ped.getFaDNAIndex(i) == -1
            && dnaLookup.get(ped.getFID(i) + "\t" + ped.getFA(i)) != null) {
          faIndex = qcIndexMap.get(dnaLookup.get(ped.getFID(i) + "\t" + ped.getFA(i)));
        } else {
          faIndex = qcIndexMap.get(Integer.toString(ped.getFaDNAIndex(i)));
        }
        Integer moIndex = null;
        if (ped.getMoDNAIndex(i) == -1) {
          if (dnaLookup.get(ped.getFID(i) + "\t" + ped.getMO(i)) != null) {
            moIndex = qcIndexMap.get(dnaLookup.get(ped.getFID(i) + "\t" + ped.getMO(i)));
          }
        } else {
          moIndex = qcIndexMap.get(ped.getMO(i));
        }
        if (indexInt == null) {
          sb.append(".").append("\t");
        } else {
          if (data[indexInt.intValue()] > 0.5) {
            highLRRSD = true;
          }
          sb.append(data[indexInt.intValue()]).append("\t");
        }
        if (faIndex == null) {
          sb.append(".").append("\t");
        } else {
          if (data[indexInt.intValue()] > 0.5) {
            highLRRSD = true;
          }
          sb.append(data[faIndex.intValue()]).append("\t");
        }
        if (moIndex == null) {
          sb.append(".").append("\t");
        } else {
          if (data[indexInt.intValue()] > 0.5) {
            highLRRSD = true;
          }
          sb.append(data[moIndex.intValue()]).append("\t");
        }
        data = sampQC.getDataFor("Genotype_callrate");
        if (indexInt == null) {
          sb.append(".").append("\t");
        } else {
          if (data[indexInt.intValue()] < 0.95) {
            lowCallrate = true;
          }
          sb.append(data[indexInt.intValue()]).append("\t");
        }
        if (faIndex == null) {
          sb.append(".").append("\t");
        } else {
          if (data[indexInt.intValue()] < 0.95) {
            lowCallrate = true;
          }
          sb.append(data[faIndex.intValue()]).append("\t");
        }
        if (moIndex == null) {
          sb.append(".").append("\t");
        } else {
          if (data[indexInt.intValue()] < 0.95) {
            lowCallrate = true;
          }
          sb.append(data[moIndex.intValue()]).append("\t");
        }
      }

      boolean excluded = false;
      if (sampleData != null) {
        boolean ex1 = sampleData.individualShouldBeExcluded(ped.getiDNA(i));
        if (ex1) {
          excluded = true;
        }
        sb.append(ex1 ? "1" : "0").append("\t");
        if (ped.getFaDNAIndex(i) == -1 || samples == null) {
          if (dnaLookup.get(ped.getFID(i) + "\t" + ped.getFA(i)) != null) {
            boolean ex = sampleData.individualShouldBeExcluded(dnaLookup.get(ped.getFID(i) + "\t"
                                                                             + ped.getFA(i)));
            if (ex) {
              excluded = true;
            }
            sb.append(ex ? "1" : "0").append("\t");
          } else {
            sb.append(".").append("\t");
          }
        } else {
          boolean ex = sampleData.individualShouldBeExcluded(samples[ped.getFaDNAIndex(i)]);
          if (ex) {
            excluded = true;
          }
          sb.append(ex ? "1" : "0").append("\t");
        }
        if (ped.getMoDNAIndex(i) == -1 || samples == null) {
          if (dnaLookup.get(ped.getFID(i) + "\t" + ped.getMO(i)) != null) {
            boolean ex = sampleData.individualShouldBeExcluded(dnaLookup.get(ped.getFID(i) + "\t"
                                                                             + ped.getMO(i)));
            if (ex) {
              excluded = true;
            }
            sb.append(ex ? "1" : "0").append("\t");
          } else {
            sb.append(".").append("\t");
          }
        } else {
          boolean ex = sampleData.individualShouldBeExcluded(samples[ped.getMoDNAIndex(i)]);
          if (ex) {
            excluded = true;
          }
          sb.append(ex ? "1" : "0").append("\t");
        }
      }

      sb.append((missingDNA || highLRRSD || lowCallrate || excluded) ? "0" : "1");

      writer.println(sb.toString());
    }

    writer.flush();
    writer.close();

    System.out.println(ext.getTime() + "]\tMissing genome pair data for " + missingCount
                       + " individuals");

    writeFamily(gl, pedToFAMO, childrenMap, dnaLookup, idLookup);
  }

  private String getGenomeFID(int i) {
    return genomeDNA ? ped.getiDNA(i) : ped.getFID(i);
  }

  private String getGenomeIID(int i) {
    return genomeDNA ? ped.getiDNA(i) : ped.getIID(i);
  }

  private String getGenomeFIDIID(int i) {
    return genomeDNA ? ped.getiDNA(i) + "\t" + ped.getiDNA(i)
                     : ped.getFID(i) + "\t" + ped.getIID(i);
  }

  private String formGenomeFIDIID(String fidiid) {
    if (genomeDNA) {
      int i = ped.getIndIndex(fidiid);
      if (i < 0) throw new IllegalArgumentException("Invalid fidiid: " + fidiid);
      String dna = ped.getiDNA(i);
      return dna + "\t" + dna;
    } else {
      return fidiid;
    }
  }

  private String getGenomeFatherFIDIID(int i) {
    String fid;
    String iid;
    if (genomeDNA) {
      int faIndex = ped.getIndexOfFaInIDs(i);
      iid = faIndex < 0 ? "." : ped.getiDNA(faIndex);
      fid = iid;
    } else {
      fid = ped.getFID(i);
      iid = ped.getFA(i);
    }
    return fid + "\t" + iid;
  }

  private String getGenomeMotherFIDIID(int i) {
    String fid;
    String iid;
    if (genomeDNA) {
      int moIndex = ped.getIndexOfMoInIDs(i);
      iid = moIndex < 0 ? "." : ped.getiDNA(moIndex);
      fid = iid;
    } else {
      fid = ped.getFID(i);
      iid = ped.getMO(i);
    }
    return fid + "\t" + iid;

  }

  private void writeFamily(GenomeLoader gl, HashMap<String, String[]> pedToFAMO,
                           HashMap<String, ArrayList<String>> childrenMap,
                           HashMap<String, String> dnaLookup, HashMap<String, String> idLookup) {
    PrintWriter writer;
    StringBuilder sb;

    HashSet<String> writtenPairs = new HashSet<>();

    writer = Files.getAppropriateWriter(outDir + REL_CHECKS_FILENAME);

    String header = "FID1\tIID1\tDNA1\tFID2\tIID2\tDNA2\tPUTATIVE_REL\t";
    if (gl != null) {
      header += "IBD0" + "\t" + "IBD1" + "\t" + "IBD2" + "\t" + "PIHAT" + "\t" + "DERIVED_REL"
                + "\t" + "REL_MATCH" + "\t";
    }
    writer.println(header);

    for (java.util.Map.Entry<String, ArrayList<String>> childrenList : childrenMap.entrySet()) {
      String fidiid = childrenList.getKey();
      String parentDNA = dnaLookup.get(fidiid);

      // if (parentDNA != null) {
      // continue;
      // }

      for (String childFIDIID : childrenList.getValue()) {
        String childDNA = dnaLookup.get(childFIDIID);
        if (childDNA == null) {
          continue;
        }
        String[] famo = pedToFAMO.get(childFIDIID);
        ArrayList<String> spouseChildren = new ArrayList<>();
        if (famo[0].equals(fidiid) && !".".equals(famo[1])) {
          spouseChildren = childrenMap.get(famo[1]);
        } else if (famo[1].equals(fidiid) && !".".equals(famo[0])) {
          spouseChildren = childrenMap.get(famo[0]);
        }
        HashSet<String> sameParentSibs = new HashSet<>(childrenList.getValue());

        if (parentDNA != null && !writtenPairs.contains(parentDNA + "\t" + childDNA)) {
          writtenPairs.add(parentDNA + "\t" + childDNA);
          writtenPairs.add(childDNA + "\t" + parentDNA);
          sb = new StringBuilder();
          sb.append(fidiid).append("\t");
          sb.append(parentDNA).append("\t");
          sb.append(childFIDIID).append("\t");

          sb.append(childDNA).append("\t");

          sb.append("PO").append("\t");

          if (gl != null) {
            if (parentDNA == null || childDNA == null) {
              sb.append(".\t.\t.\t.\t.\t.\t");
            } else {
              String genoFIDIID = formGenomeFIDIID(fidiid);
              HashMap<String, String> genoData = gl.pairData.get(genoFIDIID);
              String key = formGenomeFIDIID(childFIDIID);
              if (genoData == null || !genoData.containsKey(key)) {
                genoData = gl.pairData.get(key);
                key = genoFIDIID;
              }

              if (genoData != null) {
                String genoDataLine = genoData.get(key);
                if (genoDataLine != null) {
                  String[] tmpGL = genoDataLine.split(PSF.Regex.GREEDY_WHITESPACE);
                  sb.append(tmpGL[6]).append("\t").append(tmpGL[7]).append("\t").append(tmpGL[8])
                    .append("\t").append(tmpGL[9]).append("\t");
                  String rel = deriveRelationship(tmpGL[6], tmpGL[7], tmpGL[8], tmpGL[9]);
                  sb.append(rel).append("\t");
                  if (!rel.equals("parent-offspring")) {
                    sb.append("0").append("\t");
                  } else {
                    sb.append("1").append("\t");
                  }
                } else {
                  sb.append(".\t.\t.\t.\t.\t.\t");
                }
              } else {
                sb.append(".\t.\t.\t.\t.\t.\t");
              }
            }
          }
          writer.println(sb.toString());
        }

        for (String otherChild : spouseChildren) {
          if (otherChild.equals(childFIDIID)) {
            continue;
          }

          String otherChildDNA = dnaLookup.get(otherChild);
          if (otherChildDNA == null || writtenPairs.contains(otherChildDNA + "\t" + childDNA)) {
            continue;
          }
          writtenPairs.add(childDNA + "\t" + otherChildDNA);
          writtenPairs.add(otherChildDNA + "\t" + childDNA);

          sb = new StringBuilder();
          sb.append(childFIDIID).append("\t");

          if (childDNA != null) {
            sb.append(childDNA).append("\t");
          } else {
            sb.append(".\t");
          }

          sb.append(otherChild).append("\t");

          sb.append(otherChildDNA).append("\t");

          String expRel;
          if (sameParentSibs.contains(otherChild)) {
            sameParentSibs.remove(otherChild);
            sb.append("SIB").append("\t");
            expRel = "SIB";
          } else {
            sb.append("HALFSIB").append("\t");
            expRel = "HALFSIB";
          }
          if (gl != null) {

            if (childDNA == null || otherChildDNA == null) {
              sb.append(".\t.\t.\t.\t.\t.\t");
            } else {
              String genoChildFIDIID = formGenomeFIDIID(childFIDIID);
              HashMap<String, String> genoData = gl.pairData.get(genoChildFIDIID);

              String key = formGenomeFIDIID(otherChild);
              if (genoData == null || !genoData.containsKey(key)) {
                genoData = gl.pairData.get(key);
                key = genoChildFIDIID;
              }

              if (genoData != null) {
                String genoDataLine = genoData.get(key);
                if (genoDataLine != null) {
                  String[] tmpGL = genoDataLine.split(PSF.Regex.GREEDY_WHITESPACE);
                  sb.append(tmpGL[6]).append("\t").append(tmpGL[7]).append("\t").append(tmpGL[8])
                    .append("\t").append(tmpGL[9]).append("\t");
                  String rel = deriveRelationship(tmpGL[6], tmpGL[7], tmpGL[8], tmpGL[9]);
                  sb.append(rel).append("\t");
                  if ((expRel.equals("SIB") && rel.equals("sibling"))
                      || (expRel.equals("HALFSIB") && rel.equals("first cousins,halfsibs"))) {
                    sb.append("1").append("\t");
                  } else {
                    sb.append("0").append("\t");
                  }
                } else {
                  sb.append(".\t.\t.\t.\t.\t.\t");
                }
              } else {
                sb.append(".\t.\t.\t.\t.\t.\t");
              }
            }
          }
          writer.println(sb.toString());
        }
        for (String halfSib : sameParentSibs) {
          if (halfSib.equals(childFIDIID)) {
            continue;
          }
          String halfSibDNA = dnaLookup.get(halfSib);
          if (halfSibDNA == null || writtenPairs.contains(childDNA + "\t" + halfSibDNA)) {
            continue;
          }
          writtenPairs.add(childDNA + "\t" + halfSibDNA);
          writtenPairs.add(halfSibDNA + "\t" + childDNA);
          sb = new StringBuilder();
          sb.append(childFIDIID).append("\t");
          if (childDNA != null) {
            sb.append(childDNA).append("\t");
          } else {
            sb.append(".\t");
          }
          sb.append(halfSib).append("\t");
          sb.append(halfSibDNA).append("\t");

          sb.append("HALFSIB").append("\t");

          if (gl != null) {
            if (childDNA == null || halfSibDNA == null) {
              sb.append(".\t.\t.\t.\t.\t.\t");
            } else {
              String genoChildFIDIID = formGenomeFIDIID(childFIDIID);
              HashMap<String, String> genoData = gl.pairData.get(genoChildFIDIID);

              String key = formGenomeFIDIID(halfSib);
              if (genoData == null || !genoData.containsKey(key)) {
                genoData = gl.pairData.get(key);
                key = genoChildFIDIID;
              }

              if (genoData != null) {
                String genoDataLine = genoData.get(key);
                if (genoDataLine != null) {
                  String[] tmpGL = genoDataLine.split(PSF.Regex.GREEDY_WHITESPACE);
                  sb.append(tmpGL[6]).append("\t").append(tmpGL[7]).append("\t").append(tmpGL[8])
                    .append("\t").append(tmpGL[9]).append("\t");
                  String rel = deriveRelationship(tmpGL[6], tmpGL[7], tmpGL[8], tmpGL[9]);
                  sb.append(rel).append("\t");
                  if (rel.equals("first cousins,halfsibs")) {
                    sb.append("1").append("\t");
                  } else {
                    sb.append("0").append("\t");
                  }
                } else {
                  sb.append(".\t.\t.\t.\t.\t.\t");
                }
              } else {
                sb.append(".\t.\t.\t.\t.\t.\t");
              }
            }
          }
          writer.println(sb.toString());
        }

      }

    }

    if (gl != null) {
      for (String unrelLine : gl.unrelLines) {
        sb = new StringBuilder();
        String[] parts = unrelLine.split(PSF.Regex.GREEDY_WHITESPACE);

        String fidiid1;
        String dna1;
        String fidiid2;
        String dna2;
        if (genomeDNA) {
          dna1 = parts[0];
          dna2 = parts[2];
          fidiid1 = idLookup.get(dna1);
          fidiid2 = idLookup.get(dna2);
        } else {
          fidiid1 = parts[0] + "\t" + parts[1];
          fidiid2 = parts[2] + "\t" + parts[3];
          dna1 = ped.getiDNA(ped.getIndIndex(fidiid1));
          dna2 = ped.getiDNA(ped.getIndIndex(fidiid2));
        }
        sb.append(fidiid1).append("\t").append(dna1).append("\t");
        sb.append(fidiid2).append("\t").append(dna2).append("\t");
        sb.append("UN").append("\t").append(parts[6]).append("\t").append(parts[7]).append("\t")
          .append(parts[8]).append("\t").append(parts[9]).append("\t");
        String rel = deriveRelationship(parts[6], parts[7], parts[8], parts[9]);
        sb.append(rel).append("\t");
        if (!rel.equals("UN")) {
          sb.append("0").append("\t");
        } else {
          sb.append("1").append("\t");
        }
        writer.println(sb.toString());
      }
    }

    writer.flush();
    writer.close();

  }

  private String deriveRelationship(String z0Str, String z1Str, String z2Str, String piHatStr) {
    double z0 = Double.parseDouble(z0Str);
    double z1 = Double.parseDouble(z1Str);
    double z2 = Double.parseDouble(z2Str);
    double piHat = Double.parseDouble(piHatStr);

    for (int i = 0; i < Plink.THRESHOLDS.length; i++) {
      if (z0 >= Plink.THRESHOLDS[i][0] && z1 >= Plink.THRESHOLDS[i][1]
          && z2 >= Plink.THRESHOLDS[i][2] && piHat >= Plink.THRESHOLDS[i][3]) {
        return Plink.FLAGS[i];
      }
    }
    return "UN";
  }

  public static String parseOutputDirectory(Project proj) {
    return proj.PROJECT_DIRECTORY.getValue();
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String projFile = null;
    String ped = null;
    String mendel = null;
    String genome = null;
    String out = null;
    boolean genomeDNA = false;

    String usage = "\n" + "gwas.PlinkMendelianChecker requires 0-1 arguments\n"
                   + "   (1) Project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + "  OR \n"
                   + "   (1) File with pedigree data (i.e. pedigree=pedigree.dat (not the default))\n"
                   + "   (2) (optional) File with Mendelian Error data (i.e. mendel=markerQualityChecks.mendel (not the default))\n"
                   + "   (3) (optional) File with genomic cluster data (i.e. genome=plink.genome (not the default))\n"
                   + "   (4) (optional) If a genomeic cluster file is specified, specify if the id columns are FID/IID or DNA/DNA (i.e. genomeDNA=TRUE (not the default)) \n"
                   + "   (5) (optional) Directory of output (i.e. out=/path/to/dir/ (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        projFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("pedigree=")) {
        ped = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("mendel=")) {
        mendel = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("genome=")) {
        genome = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        out = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("genomeDNA=")) {
        genomeDNA = ext.parseBooleanArg(arg);
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
      if (projFile != null) {
        (new PlinkMendelianChecker(new Project(projFile))).run();
      } else if (ped != null) {
        if (out == null) {
          out = ext.parseDirectoryOfFile(ped);
        }
        (new PlinkMendelianChecker(ped, mendel, genome, genomeDNA, out,
                                   new Logger("plinkMendelianChecks.log"))).run();
      } else {
        System.err.println(usage);
        System.exit(1);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
