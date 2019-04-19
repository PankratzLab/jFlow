package org.genvisis.cnv.filesys;

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;

import org.genvisis.cnv.qc.MendelErrors;
import org.genvisis.cnv.qc.MendelErrors.MendelErrorCheck;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.cnv.var.SampleData;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.FamilyStructure;

public class Pedigree extends FamilyStructure {

  private final Project project;
  private final boolean nullProject;
  private final int[][] dnaIndicesInProject;
  public static final int MISSING_DNA_INDEX = -1;
  boolean projectOrder;
  private final String pedigreeFile;

  public Pedigree(Project proj) {
    this(proj, proj.PEDIGREE_FILENAME.getValue());
  }

  /**
   * @param proj
   * @param pedigreeFile
   */
  public Pedigree(Project proj, String pedigreeFile) {
    super(pedigreeFile, true, proj == null ? new Logger() : proj.getLog());
    project = proj;
    nullProject = proj == null;
    dnaIndicesInProject = new int[ids.length][];
    SampleData sampleData = nullProject ? null : proj.getSampleData(false);
    String[] samples = nullProject ? null : proj.getSamples();
    if (dnas != null && samples != null && samples.length == dnas.length) {
      projectOrder = true;
      // Compare project samples and dnas. If we find a mismatch, set projectOrder to false
      for (int i = 0; i < samples.length && projectOrder; i++) {
        projectOrder = samples[i].equals(dnas[i]);
      }
    }

    Map<String, Integer> sampleIndices = samples == null ? null : ArrayUtils.indexMap(samples);
    this.pedigreeFile = pedigreeFile;
    for (int i = 0; i < ids.length; i++) {
      int iDNAIndex = MISSING_DNA_INDEX;
      int faDNAIndex = MISSING_DNA_INDEX;
      int moDNAIndex = MISSING_DNA_INDEX;

      if (proj != null) {
        iDNAIndex = getSampleIndex(dnas[i], sampleData, sampleIndices);

        int faIDIndex = getIndexOfFaInIDs(i);
        if (faIDIndex != -1) {
          faDNAIndex = getSampleIndex(dnas[faIDIndex], sampleData, sampleIndices);
        }

        int moIDIndex = getIndexOfMoInIDs(i);
        if (moIDIndex != -1) {
          moDNAIndex = getSampleIndex(dnas[moIDIndex], sampleData, sampleIndices);
        }
      }

      dnaIndicesInProject[i] = new int[] {iDNAIndex, faDNAIndex, moDNAIndex};
    }
  }

  private static int getSampleIndex(String sample, SampleData sampleData,
                                    Map<String, Integer> projectSamples) {
    int sampleIndex = MISSING_DNA_INDEX;
    if (sample != null && !sample.equals(FamilyStructure.MISSING_ID_STR) && sampleData != null
        && projectSamples != null && sampleData.lookup(sample) != null) {
      sampleIndex = projectSamples.get(sampleData.lookup(sample)[0]);
    }
    return sampleIndex;
  }

  /**
   * @return true iff the order and number of samples in this pedigree matches those of the
   *         (non-null) accompanying project.
   */
  public boolean isProjectOrder() {
    return projectOrder;
  }

  public Project getProject() {
    return project;
  }

  public int getIDNAIndex(int index) {
    return dnaIndicesInProject[index][0];
  }

  public int getFaDNAIndex(int index) {
    return dnaIndicesInProject[index][1];
  }

  public int getMoDNAIndex(int index) {
    return dnaIndicesInProject[index][2];
  }

  public String getPedigreeFile() {
    return pedigreeFile;
  }

  /**
   * @return Mapping of sample ID to {@link MendelErrorCheck} instance for that individual. Note
   *         that the keys are always for children, with parents being reachable through the error
   *         check.
   */
  public static Map<String, MendelErrorCheck> checkMendelErrors(Pedigree pedigree,
                                                                MarkerData markerData,
                                                                boolean[] samplesToCheck,
                                                                String[] sex,
                                                                ClusterFilterCollection clusterFilters,
                                                                float gcThreshold, Logger log) {
    Project proj = pedigree.getProject();
    if (proj == null) {
      log.reportError(ext.getTime() + "]\t Error - cannot run checkMendelErrors without a Project");
      return null;
    }
    Map<String, MendelErrorCheck> mendelErrorChecks = new HashMap<>();
    byte[] genotypes = markerData.getAbGenotypesAfterFilters(clusterFilters,
                                                             markerData.getMarkerName(),
                                                             gcThreshold, log);
    for (int i = 0; i < pedigree.getIDs().length; i++) {
      int sampleIndex = pedigree.getIDNAIndex(i);
      MendelErrors mendelErrors = null;
      if (sampleIndex >= 0 && (samplesToCheck == null || samplesToCheck[sampleIndex])) {
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
            sampleSex = Integer.parseInt(sex[sampleIndex]);
          }
        } catch (NumberFormatException nfe) {

        }
        // System.out.println(faGenotype+"\t"+moGenotype);
        mendelErrors = new MendelErrors(markerData.getChr(), sampleSex, genotypes[sampleIndex],
                                        faGenotype, moGenotype);
      } else {
        mendelErrors = new MendelErrors(markerData.getChr(), -1, (byte) -1, (byte) -1, (byte) -1);
      }
      mendelErrorChecks.put(pedigree.getiDNA(i), mendelErrors.checkMendelError());
    }
    return mendelErrorChecks;
  }

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
                                             true, true);
    } else {
      log.report("Warning - no sex check file found, sex will be set to '0' for all individuals.  Otherwise, first run SexChecks and then re-run.");
      sexDict = new Hashtable<>();
    }

    samples = samples == null ? proj.getSamples() : samples;
    if (Files.exists(proj.SAMPLE_DATA_FILENAME.getValue())) {
      sd = proj.getSampleData(false);
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
      proj = new Project(filename);
      build(proj, out, null, overwrite);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
