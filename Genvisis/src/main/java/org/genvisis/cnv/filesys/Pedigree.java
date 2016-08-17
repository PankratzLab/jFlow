package org.genvisis.cnv.filesys;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.cnv.qc.MendelErrors;
import org.genvisis.cnv.qc.MendelErrors.MendelErrorCheck;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.FamilyStructure;



public class Pedigree extends FamilyStructure {

  public static class PedigreeUtils {

    public static MendelErrorCheck[] checkMendelErrors(Pedigree pedigree, MarkerData markerData,
        boolean[] samplesToCheck, String[] sex, ClusterFilterCollection clusterFilters,
        float gcThreshold, Logger log) {
      if (pedigree.getProject() == null) {
        log.reportError(
            ext.getTime() + "]\t Error - cannot run checkMendelErrors without a Project");
        return null;
      }
      MendelErrorCheck[] mendelErrorChecks =
          new MendelErrorCheck[pedigree.getProject().getSamples().length];
      byte[] genotypes = markerData.getAbGenotypesAfterFilters(clusterFilters,
          markerData.getMarkerName(), gcThreshold, log);
      if (!pedigree.isProjectOrder()) {
        log.reportTimeError("Pedigree file must be in project order, internal error");
        return null;
      } else {
        for (int i = 0; i < pedigree.getIDs().length; i++) {
          if (/* this check isn't valid *//* pedigreeEntries[i] != null && */(samplesToCheck == null
              || samplesToCheck[i]) && pedigree.getIDNAIndex(i) >= 0) {
            int sampleIndex = pedigree.getIDNAIndex(i);
            int faDNAIndex = pedigree.getFaDNAIndex(i);
            int moDNAIndex = pedigree.getMoDNAIndex(i);

            byte faGenotype = -1;
            if (faDNAIndex >= 0 && (samplesToCheck == null || samplesToCheck[faDNAIndex])) {
              faGenotype = genotypes[faDNAIndex];
            }
            byte moGenotype = -1;
            if (moDNAIndex >= 0 && (samplesToCheck == null || samplesToCheck[moDNAIndex])) {
              moGenotype = genotypes[moDNAIndex];
            }
            int sampleSex = -1;
            try {
              if (sex != null) {
                sampleSex = Integer.parseInt(sex[pedigree.getIDNAIndex(i)]);
              }
            } catch (NumberFormatException nfe) {

            }
            // System.out.println(faGenotype+"\t"+moGenotype);
            MendelErrors mendelErrors = new MendelErrors(markerData.getChr(), sampleSex,
                genotypes[sampleIndex], faGenotype, moGenotype);
            mendelErrorChecks[i] = mendelErrors.checkMendelError();
          } else {
            mendelErrorChecks[i] =
                new MendelErrors(markerData.getChr(), -1, (byte) -1, (byte) -1, (byte) -1)
                    .checkMendelError();
          }
        }
      }
      return mendelErrorChecks;
    }

    public static ArrayList<int[]> loadCompleteTrios(FamilyStructure ped,
        HashSet<String> excludedFIDIIDs, HashSet<String> includedFIDIIDs, boolean cache) {
      if (ped.cached_complete_trios != null) {
        return ped.cached_complete_trios;
      }
      ArrayList<String[]> allTrios = new ArrayList<String[]>();
      ArrayList<int[]> trios = new ArrayList<int[]>();
      for (int i = 0; i < ped.getIDs().length; i++) {
        int faInd = -1;
        int moInd = -1;
        if (!FamilyStructure.MISSING_ID_STR.equals(ped.getFA(i))
            && (excludedFIDIIDs == null
                || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getFA(i))
                    && (includedFIDIIDs == null
                        || includedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getFA(i))))
            && !FamilyStructure.MISSING_ID_STR.equals(ped.getMO(i))
            && (excludedFIDIIDs == null
                || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getMO(i))
                    && (includedFIDIIDs == null
                        || includedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getMO(i))))) {
          if (cache) {
            allTrios.add(new String[] {ped.getFID(i) + "\t" + ped.getIID(i),
                ped.getFID(i) + "\t" + ped.getFA(i), ped.getFID(i) + "\t" + ped.getMO(i)});
          }
          if ((faInd = ped.getIndexOfFaInIDs(i)) >= 0 && (moInd = ped.getIndexOfMoInIDs(i)) >= 0) {
            trios.add(new int[] {i, faInd, moInd});
          }
        }
      }
      if (cache) {
        ped.cached_complete_trios = trios;
        ped.cached_all_trios = allTrios;
      }
      return trios;
    }

    public static HashMap<String, ArrayList<String>> loadParentToChildrenMap(FamilyStructure ped,
        boolean completeOnly, HashSet<String> excludedFIDIIDs, HashSet<String> includedFIDIIDs,
        boolean cache) {
      if (ped.cached_parentToChildrenMap != null) {
        return ped.cached_parentToChildrenMap;
      }
      HashMap<String, ArrayList<String>> parentMap = new HashMap<String, ArrayList<String>>();

      for (int i = 0; i < ped.getIDs().length; i++) {
        if (!FamilyStructure.MISSING_ID_STR.equals(ped.getFA(i))
            && (excludedFIDIIDs == null
                || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getFA(i))
                    && (includedFIDIIDs == null
                        || includedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getFA(i))))
            && (!completeOnly || (/* faInd = */ped.getIndexOfFaInIDs(i)) >= 0)) {
          ArrayList<String> children = parentMap.get(ped.getFID(i) + "\t" + ped.getFA(i));
          if (children == null) {
            children = new ArrayList<String>();
            parentMap.put(ped.getFID(i) + "\t" + ped.getFA(i), children);
          }
          children.add(ped.getFID(i) + "\t" + ped.getIID(i));
        }
        if (!FamilyStructure.MISSING_ID_STR.equals(ped.getMO(i))
            && (excludedFIDIIDs == null
                || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getMO(i))
                    && (includedFIDIIDs == null
                        || includedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getMO(i))))
            && (!completeOnly || (/* moInd = */ped.getIndexOfMoInIDs(i)) >= 0)) {
          ArrayList<String> children = parentMap.get(ped.getFID(i) + "\t" + ped.getMO(i));
          if (children == null) {
            children = new ArrayList<String>();
            parentMap.put(ped.getFID(i) + "\t" + ped.getMO(i), children);
          }
          children.add(ped.getFID(i) + "\t" + ped.getIID(i));
        }
      }
      if (cache) {
        ped.cached_parentToChildrenMap = parentMap;
        ped.cached_parentMapIsCompleteOnly = completeOnly;
      }
      return parentMap;
    }

    public static ArrayList<String[]> loadPOPairs(FamilyStructure ped, boolean completeOnly,
        HashSet<String> excludedFIDIIDs, HashSet<String> includedFIDIIDs, boolean cache) {
      if (ped.cached_poPairsIDs != null) {
        return ped.cached_poPairsIDs;
      }
      ArrayList<String[]> pairs = new ArrayList<String[]>();
      ArrayList<int[]> completePairs = new ArrayList<int[]>();
      for (int i = 0; i < ped.getIDs().length; i++) {
        int faInd = -1;
        int moInd = -1;
        if (!FamilyStructure.MISSING_ID_STR.equals(ped.getFA(i))
            && (excludedFIDIIDs == null
                || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getFA(i))
                    && (includedFIDIIDs == null
                        || includedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getFA(i))))
            && (!completeOnly || (faInd = ped.getIndexOfFaInIDs(i)) >= 0)) {
          pairs.add(new String[] {ped.getFID(i) + "\t" + ped.getFA(i),
              ped.getFID(i) + "\t" + ped.getIID(i)});
          if (completeOnly) {
            completePairs.add(new int[] {faInd, i});
          }
        }
        if (!FamilyStructure.MISSING_ID_STR.equals(ped.getMO(i))
            && (excludedFIDIIDs == null
                || !excludedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getMO(i))
                    && (includedFIDIIDs == null
                        || includedFIDIIDs.contains(ped.getFID(i) + "\t" + ped.getMO(i))))
            && (!completeOnly || (moInd = ped.getIndexOfMoInIDs(i)) >= 0)) {
          pairs.add(new String[] {ped.getFID(i) + "\t" + ped.getMO(i),
              ped.getFID(i) + "\t" + ped.getIID(i)});
          if (completeOnly) {
            completePairs.add(new int[] {moInd, i});
          }
        }
      }
      if (cache) {
        ped.cached_poPairsIDs = pairs;
        if (completeOnly) {
          ped.cached_poPairsCompleteOnly = completePairs;
        }
      }
      return pairs;
    }

    // also loads trios and p-o hash
    public static ArrayList<String[]> loadSibs(FamilyStructure ped, boolean completeOnly,
        HashSet<String> excludedFIDIIDs, HashSet<String> includedFIDIIDs, boolean cache) {
      if (ped.cached_all_trios == null) {
        loadCompleteTrios(ped, excludedFIDIIDs, includedFIDIIDs, true); // will also create
                                                                        // all_trios
      }
      if (ped.cached_all_trios == null) {
        // ERROR
      }
      if (ped.cached_sib_pairs != null) {
        return ped.cached_sib_pairs;
      }
      HashMap<String, ArrayList<String>> parentToChildren =
          loadParentToChildrenMap(ped, completeOnly, excludedFIDIIDs, includedFIDIIDs, cache);

      // at this point, only non-excluded IDs are present in all_trios and parentToChildren
      ArrayList<String[]> sibPairs = new ArrayList<String[]>();
      for (String[] trio : ped.cached_all_trios) {
        ArrayList<String> faChildren = parentToChildren.get(trio[1]);
        if (faChildren == null) {
          // Error!
        } else if (faChildren.size() == 1 && !trio[0].equals(faChildren.get(0))) {
          // Error!
        }
        ArrayList<String> moChildren = parentToChildren.get(trio[2]);
        if (moChildren == null) {
          // Error!
        } else if (moChildren.size() == 1 && !trio[0].equals(faChildren.get(0))) {
          // Error!
        }
        HashSet<String> unionSet = new HashSet<String>();
        unionSet.addAll(faChildren);
        unionSet.retainAll(moChildren);
        if (unionSet.size() == 0) {
          continue; // no sibs
        } else {
          for (String sib : unionSet) {
            if (!sib.equals(trio[0])) {
              sibPairs.add(new String[] {sib, trio[0]});
              sibPairs.add(new String[] {trio[0], sib});
            }
          }
        }
      }
      if (cache) {
        ped.cached_sib_pairs = sibPairs;
      }
      return sibPairs;
    }


  }

  public static final int MISSING_DNA_INDEX = -1;

  public static void build(Project proj, String newPedFile, String[] samples, boolean overwrite) {
    // FID IID FAID MOID SEX PHENO DNA MZTWINID
    PrintWriter writer;
    SampleData sd = null;
    String file;
    Logger log;
    Hashtable<String, String> sexDict;

    log = proj.getLog();
    file = newPedFile == null ? proj.PEDIGREE_FILENAME.getValue() : newPedFile;
    if (Files.exists(file) && !overwrite) {
      log.reportError("Error - file " + file
          + " already exists; overwrite flag must be set or a different output file must be specified.");
      return;
    }

    String sexFile = proj.SEXCHECK_RESULTS_FILENAME.getValue();

    if (Files.exists(sexFile)) {
      log.report("Found sex check file.");
      // sexDict = HashVec.loadFileToHashString(sexFile, new int[]{1, 2}, new int[]{3, 4}, false,
      // "\t", true, false, true);
      sexDict = HashVec.loadFileToHashString(sexFile, new int[] {0}, new int[] {3, 4}, false, "\t",
          true, false, true);
    } else {
      log.report(
          "Warning - no sex check file found, sex will be set to '0' for all individuals.  Otherwise, first run SexChecks and then re-run.");
      sexDict = new Hashtable<String, String>();
    }

    samples = samples == null ? proj.getSamples() : samples;
    if (Files.exists(proj.SAMPLE_DATA_FILENAME.getValue())) {
      sd = proj.getSampleData(0, false);
    }

    writer = Files.getAppropriateWriter(file);
    for (String sample : samples) {
      String[] ids = sd == null ? new String[] {sample, sample + "\t" + sample} : sd.lookup(sample);
      String sexCode = "0";
      if (sexDict.size() > 0 && sexDict.containsKey(sample)) {
        sexCode = sexDict.get(sample).split("\t")[1]; // est. sex
      } else if (sd != null) {
        sexCode = "" + sd.getSexForIndividual(sample); // putative sex
      }
      sexCode = "" + SexChecks.mapEstimatedSexToSex(sexCode);
      writer.println(ids[1] + "\t0\t0\t" + sexCode + "\t1\t" + ids[0] + "\t.");
    }
    writer.flush();
    writer.close();

  }

  private static int getSampleIndex(String sample, SampleData sampleData, String[] projectSamples) {
    int sampleIndex = MISSING_DNA_INDEX;
    if (sample != null && !sample.equals(FamilyStructure.MISSING_ID_STR) && sampleData != null
        && projectSamples != null && sampleData.lookup(sample) != null) {
      sampleIndex = ext.indexOfStr(sampleData.lookup(sample)[0], projectSamples);
    }
    return sampleIndex;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    String out = null;
    boolean overwrite = false;
    Project proj;

    String usage = "\n" + "cnv.filesys.Pedigree requires 1-3 arguments\n"
        + "   (1) project properties filename (i.e. proj="
        + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
        + "   (2) OPTIONAL: Pedigree output filename (if omitted, output file will be project property) (i.e. out="
        + out + " (default))\n"
        + "   (1) OPTIONAL: overwrite flag, will overwrite the output file if it already exists (i.e. -overwrite (not the default))\n"
        + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        out = arg.split("=")[1];
        numArgs--;
      } else if (arg.equals("-overwrite")) {
        overwrite = true;
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
      proj = new Project(filename, false);
      build(proj, out, null, overwrite);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private final Project project;
  private final boolean nullProject;

  private final int[][] dnaIndicesInProject;

  boolean projectOrder;

  private final String pedigreeFile;

  public Pedigree(Project proj) {
    this(proj, proj.PEDIGREE_FILENAME.getValue(), true);
  }

  /**
   * 
   * @param proj
   * @param projectOrder (Used for ProjectUtils.checkMendelErrors()) Indicates if the order of
   *        entries in the project's pedigree file (from the PEDIGREE_FILENAME property) matches the
   *        internal order of the project samples.
   */
  public Pedigree(Project proj, boolean projectOrder) {
    this(proj, proj.PEDIGREE_FILENAME.getValue(), projectOrder);
  }

  /**
   * 
   * @param proj
   * @param pedigreeFile
   * @param projectOrder (Used for ProjectUtils.checkMendelErrors()) Indicates if the order of
   *        entries in the given pedigree file matches the internal order of the project samples.
   */
  public Pedigree(Project proj, String pedigreeFile, boolean projectOrder) {
    super(pedigreeFile, true, proj == null ? new Logger() : proj.getLog());
    project = proj;
    nullProject = proj == null;
    dnaIndicesInProject = new int[ids.length][];
    this.projectOrder = projectOrder;
    SampleData sampleData = nullProject ? null : proj.getSampleData(0, false);
    String[] samples = nullProject ? null : proj.getSamples();
    this.pedigreeFile = pedigreeFile;
    for (int i = 0; i < ids.length; i++) {
      int iDNAIndex = MISSING_DNA_INDEX;
      int faDNAIndex = MISSING_DNA_INDEX;
      int moDNAIndex = MISSING_DNA_INDEX;

      if (proj != null) {
        iDNAIndex = getSampleIndex(dnas[i], sampleData, samples);

        int faIDIndex = getIndexOfFaInIDs(i);
        if (faIDIndex != -1) {
          faDNAIndex = getSampleIndex(dnas[faIDIndex], sampleData, samples);
        }

        int moIDIndex = getIndexOfMoInIDs(i);
        if (moIDIndex != -1) {
          moDNAIndex = getSampleIndex(dnas[moIDIndex], sampleData, samples);
        }
      }

      dnaIndicesInProject[i] = new int[] {iDNAIndex, faDNAIndex, moDNAIndex};
    }

  }

  public int getFaDNAIndex(int index) {
    return dnaIndicesInProject[index][1];
  }

  public int getIDNAIndex(int index) {
    return dnaIndicesInProject[index][0];
  }

  public int getMoDNAIndex(int index) {
    return dnaIndicesInProject[index][2];
  }

  public String getPedigreeFile() {
    return pedigreeFile;
  }

  public Project getProject() {
    return project;
  }

  public boolean isProjectOrder() {
    return projectOrder;
  }

}
