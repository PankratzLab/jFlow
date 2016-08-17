package org.genvisis.cnv.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.cnv.analysis.pca.PrincipalComponentsResiduals;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.stats.Quantiles;

/**
 * Class that automates adding sample QC metrics to sample data, and also parses the qc metrics to
 * the desired quantile class
 * <p>
 * Currently requires that all samples are present
 *
 */
public class SampleQC {
  private static double CONTAM_LRR_SD_THRESHOLD = 0.35;
  private static double CONTAM_CALLRATE_THRESHOLD = 0.9;

  private static double EXCLUDE_LRR_SD_THRESHOLD = 0.5;

  private static void appendToSampleData(Project proj, Hashtable<String, String> hashtable,
      String[] header, int numQ, boolean justQuantiles) {
    SampleData sampledata = proj.getSampleData(0, false);
    proj.getLog()
        .reportTimeInfo("Adding " + header.length + " columns to sample data based on sample QC");
    sampledata.addData(hashtable, "DNA", header, ".", "\t", proj.getLog());
  }

  private static String[] developMetricsHeader(Quantiles[] quantiles, String[] titles, int numQ,
      boolean justQuantiles) {
    String[] qcHeader = new String[justQuantiles ? quantiles.length : quantiles.length * 2];
    int curIndex = 0;
    for (int i = 0; i < quantiles.length; i++) {
      if (!justQuantiles) {
        qcHeader[curIndex] = titles[i];
        curIndex++;
      }
      qcHeader[curIndex] = getClassForQuantile(quantiles[i], titles[i], numQ);
      curIndex++;
    }
    return qcHeader;
  }

  private static String developMetricsLine(int sampleIndex, Quantiles[] quantiles,
      double[][] variableDomMatrix, boolean justQuantiles) {
    String line = "";
    for (int i = 0; i < quantiles.length; i++) {
      line += (i == 0 ? "" : "\t");
      if (!justQuantiles) {
        line += variableDomMatrix[i][sampleIndex] + "\t";
      }
      line += quantiles[i].getQuantileMembershipAsRoundedInt()[sampleIndex];
    }
    return line;
  }

  /**
   * Creates sample data friendly class titles for column headers
   */
  private static String getClassForQuantile(Quantiles quantiles, String qcTitle, int numQ) {
    String thisClass = "CLASS=QUANTILE_" + numQ + "_" + qcTitle;
    int[] uniqLabels = Array.toIntArray(
        Array.unique(Array.toStringArray(quantiles.getQuantileMembershipAsRoundedInt())));
    double[] uniqQs =
        Array.toDoubleArray(Array.unique(Array.toStringArray(quantiles.getQuantileMembership())));
    int[] orderLabels = Sort.quicksort(uniqLabels);
    int[] orderQs = Sort.quicksort(uniqQs);
    for (int i = 0; i < orderLabels.length; i++) {
      thisClass += ";" + uniqLabels[orderLabels[i]] + "=q_" + ext.formDeci(uniqQs[orderQs[i]], 3);
    }
    return thisClass;
  }

  public static SampleQC loadSampleQC(Project proj) {
    return loadSampleQC(proj, LrrSd.SAMPLE_COLUMN, LrrSd.NUMERIC_COLUMNS, false, false, null);
  }

  public static SampleQC loadSampleQC(Project proj, boolean generate, boolean gcCorrectedLrrSd,
      String duplicatesSetFile) {
    return loadSampleQC(proj, LrrSd.SAMPLE_COLUMN, LrrSd.NUMERIC_COLUMNS, generate,
        gcCorrectedLrrSd, duplicatesSetFile);
  }

  /**
   * @param proj
   * @param sampleColumnName header of the column containing sample names
   * @param qcTitlesToLoad qc titles to load from the sample QC file
   * @param generateSampleQC generate sampleQC if missing
   * @param duplicatesSetFile filename for file with 3 columns: FID IID DuplicateID, null to not
   *        check duplicates
   * @return
   */
  public static SampleQC loadSampleQC(Project proj, String sampleColumnName,
      String[] qcTitlesToLoad, boolean generateSampleQC, boolean gcCorrectedLrrSd,
      String duplicatesSetFile) {
    // String lrrSdToLoad = proj.getFilename(proj.SAMPLE_QC_FILENAME);
    String lrrSdToLoad = proj.SAMPLE_QC_FILENAME.getValue();
    SampleQC sampleQC = null;
    if (!Files.exists(lrrSdToLoad) && !generateSampleQC) {
      proj.getLog().reportTimeError("Could not find sample QC file " + lrrSdToLoad);
    } else {
      if (!Files.exists(lrrSdToLoad) && generateSampleQC) {
        proj.getLog().reportTimeInfo("Attempting to generate sample QC file " + lrrSdToLoad);
        LrrSd.init(proj, null, null, proj.NUM_THREADS.getValue());
        if (!Files.exists(lrrSdToLoad)) {
          proj.getLog().reportTimeError("Could not generate sample QC file " + lrrSdToLoad);
        }
      }
      proj.getLog().reportTimeInfo("Loading qc data from " + lrrSdToLoad);
      try {
        BufferedReader reader = Files.getAppropriateReader(lrrSdToLoad);
        String[] header = reader.readLine().trim().split("[\\s]+");
        int[] indicesToLoad =
            ext.indexFactors(qcTitlesToLoad, header, true, proj.getLog(), true, false);
        int sampleColumn = ext.indexOfStr(sampleColumnName, header);
        if (Array.countIf(indicesToLoad, -1) > 0 || sampleColumn < 0) {
          proj.getLog()
              .reportTimeError("Could not find all desired columns in qc file " + lrrSdToLoad);
          proj.getLog().reportTimeError(
              "Consider re-creating " + lrrSdToLoad + " if sample qc has been updated");

        } else if (!verifyAllProjectSamples(proj, lrrSdToLoad, sampleColumn)) {
          proj.getLog().reportTimeError(
              "Could not find all of the projects samples in qc file " + lrrSdToLoad);
        } else {
          sampleQC = new SampleQC(proj, qcTitlesToLoad, gcCorrectedLrrSd);

          while (reader.ready()) {
            String[] line = reader.readLine().trim().split("[\\s]+");
            String sample = line[sampleColumn];
            if (ext.indexOfStr(sampleColumnName, line) < 0) {
              for (int i = 0; i < indicesToLoad.length; i++) {
                double data = Double.NaN;
                try {
                  data = Double.parseDouble(line[indicesToLoad[i]]);
                } catch (NumberFormatException e) {
                  proj.getLog().reportTimeWarning("line " + Array.toStr(line)
                      + " contained an invalid number for qc column " + qcTitlesToLoad[i]);
                }
                sampleQC.addToMatrix(sample, i, data);
              }
            }
          }
          if (!sampleQC.verify()) {
            proj.getLog().reportTimeError("Could not verify that all data has been added");
            return null;
          }
        }
        reader.close();
        proj.getLog().reportTimeInfo("Finished loading qc data from " + lrrSdToLoad);

      } catch (FileNotFoundException fnfe) {
        proj.getLog()
            .reportTimeError("file \"" + lrrSdToLoad + "\" not found in current directory");
        return null;
      } catch (IOException ioe) {
        proj.getLog().reportTimeError("Error reading file \"" + lrrSdToLoad + "\"");
        return null;
      }

      proj.getLog().reportTimeInfo("Filtering empty columns from " + lrrSdToLoad);
      int numFiltered = sampleQC.removeEmptyMetrics();
      if (numFiltered == -1) {
        return null;
      }
      proj.getLog()
          .reportTimeInfo("Filtered " + numFiltered + " empty columns from " + lrrSdToLoad);

      proj.getLog().reportTimeInfo("Finding samples to exclude");
      int numExcluded = sampleQC.addExcludes();
      proj.getLog().reportTimeInfo("Found " + numExcluded + " samples to exclude");

      if (sampleQC.addPedigreeData()) {
        if (duplicatesSetFile != null) {
          sampleQC.addDuplicatesData(duplicatesSetFile);
        }
      }
      sampleQC.addUse_cnvs();
    }
    return sampleQC;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    int numQ = 5;
    int numPCs = 0;
    boolean justQuantiles = true;
    boolean gcCorrectedLrrSd = false;
    String duplicatesSetFile = null;
    boolean correctFidIids = false;

    String usage = "\n" + "cnv.qc.SampleQC requires 0-1 arguments\n";
    usage += "   (1) filename (i.e. proj=" + filename + " (default))\n" + "";
    usage +=
        "   (2) number of quantiles to divide the sample QC metrics into (5 = quintiles, 100 = percentiles) (i.e. numQ="
            + numQ + " (default))\n" + "";
    usage += "   (3) add only class (quantiled) qc data to sample data  (i.e. justQuantiles="
        + justQuantiles + " (default))\n" + "";
    usage +=
        "   (4) if a pc file is available, add this many pcs to the sample data file , must be set to >=1 to be added  (i.e. numPCs="
            + numPCs + " (default,no addition))\n" + "";
    usage += "   (5) use GC corrected LRR SD for filtering (i.e. gcCorrectedLrrSd="
        + gcCorrectedLrrSd + " (default))\n" + "";
    usage +=
        "   (6) duplicates set file with 3 columns (FID IID DuplicateID) to ID duplicates, also adds use columns  (i.e. duplicatesSetFile=\\quality_control\\genome\\plink.genome_duplicatesSet.dat (not the default))\n"
            + "";
    usage += "   (7) replace fid and iid with pedigree values (i.e. correctFidIids="
        + correctFidIids + " (default))\n" + "";

    usage +=
        "   NOTE: the projects sample qc file must be present, for the qc metrics to be added to sample data";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("numQ=")) {
        numQ = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("numPCs=")) {
        numPCs = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("justQuantiles=")) {
        justQuantiles = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("gcCorrectedLrrSd=")) {
        gcCorrectedLrrSd = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("duplicatesSetFile=")) {
        duplicatesSetFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("correctFidIids=")) {
        correctFidIids = ext.parseBooleanArg(arg);
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
      Project proj = new Project(filename, false);
      parseAndAddToSampleData(proj, numQ, numPCs, justQuantiles, gcCorrectedLrrSd,
          duplicatesSetFile, correctFidIids);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parseAndAddToSampleData(Project proj, int numQ, int numPCs,
      boolean justQuantiles, boolean gcCorrectedLrrSd, String duplicatesSetFile,
      boolean correctFidIids) {
    // TODO Make gcCorrectedLrrSd functional, put FID/IID in appropriate columns (2&3?)
    SampleQC sampleQC = loadSampleQC(proj, false, gcCorrectedLrrSd, duplicatesSetFile);
    sampleQC.addQCsToSampleData(numQ, numPCs, justQuantiles);
    if (numPCs > 0) {
      sampleQC.addPCsToSampleData(numQ, numPCs, justQuantiles);
    }
    if (correctFidIids) {
      sampleQC.correctFidIidCols();
    }
  }

  private static boolean verifyAllProjectSamples(Project proj, String lrrSdToLoad,
      int sampleColumn) {
    String[] projSamples = proj.getSamples();
    String[] fileSamples =
        HashVec.loadFileToStringArray(lrrSdToLoad, false, new int[] {sampleColumn}, false);
    int[] indices =
        ext.indexLargeFactors(fileSamples, projSamples, true, proj.getLog(), false, false);
    return fileSamples.length - Array.countIf(indices, -1) == projSamples.length;
  }

  private final Project proj;
  private final String[] samples;
  private double[][] qcMatrix;
  private String[] qctitles;

  private int qcsAdded;
  private String[] fidiids;
  private HashMap<String, Integer> fidiidToIndex;

  private final boolean gcCorrectedLrrSd;

  private final boolean[] excludes;

  private final String[] excludeNotes;

  private String[] mzTwinIds;

  private boolean checkDuplicates;

  private String[] duplicateIds;

  private boolean[] uses;

  private String[] useNotes;

  private boolean[] use_cnvs;

  private String[] use_cnvNotes;

  private SampleQC(Project proj, String[] qctitles, boolean gcCorrectedLrrSd) {
    super();
    this.proj = proj;
    samples = proj.getSamples();

    qcMatrix = new double[qctitles.length][samples.length];
    this.qctitles = qctitles;
    qcsAdded = 0;

    this.gcCorrectedLrrSd = gcCorrectedLrrSd;

    excludes = Array.booleanArray(samples.length, false);
    excludeNotes = Array.stringArray(samples.length, ".");

    checkDuplicates = false;
  }

  private void addDuplicatesData(String duplicatesSetFile) {
    proj.getLog().reportTimeInfo("Identifying duplicate samples");
    try {
      BufferedReader reader = Files.getAppropriateReader(duplicatesSetFile);
      duplicateIds = Array.stringArray(samples.length, ".");
      HashMap<String, HashSet<Integer>> duplicateSets = new HashMap<String, HashSet<Integer>>();
      String[] line;
      int duplicatesFound = 0;
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (line.length != 3) {
          proj.getLog().reportTimeError("file \"" + duplicatesSetFile
              + "\" contains at least one line that is not 3 columns: " + Array.toStr(line));
          return;
        }
        String fidiid = line[0] + "\t" + line[1];
        String duplicateId = line[2];
        try {
          int index = fidiidToIndex.get(fidiid);
          duplicateIds[index] = duplicateId;
          duplicatesFound++;
          HashSet<Integer> dupeSet = duplicateSets.get(duplicateId);
          if (dupeSet == null) {
            dupeSet = new HashSet<Integer>();
            duplicateSets.put(duplicateId, dupeSet);
          }
          dupeSet.add(index);
        } catch (NullPointerException npe) {
          proj.getLog().reportTimeError(
              "FID/IID pair (" + fidiid + ") in duplicates missing from project samples");
        }
      }
      checkDuplicates = true;
      proj.getLog().reportTimeInfo("Identified " + duplicatesFound + " duplicate samples");
      removeDuplicatesFromUse(duplicateSets);
    } catch (FileNotFoundException fnfe) {
      proj.getLog()
          .reportTimeError("file \"" + duplicatesSetFile + "\" not found in current directory");
    } catch (IOException ioe) {
      proj.getLog().reportTimeError("Error reading file \"" + duplicatesSetFile + "\"");
    }
  }

  private int addExcludes() {
    int numExcluded = 0;
    int callrateIndex = getCallrateIndex();
    int lrr_sdIndex = getLrr_sdIndex();

    for (int i = 0; i < samples.length; i++) {
      boolean callrate = callrateIndex != -1 && !Double.isNaN(qcMatrix[callrateIndex][i]);
      boolean lrr_sd = lrr_sdIndex != -1 && !Double.isNaN(qcMatrix[lrr_sdIndex][i]);
      if (callrate && qcMatrix[callrateIndex][i] < proj.SAMPLE_CALLRATE_THRESHOLD.getValue()) {
        if (addToExcludes(i, "Callrate < " + proj.SAMPLE_CALLRATE_THRESHOLD.getValue())) {
          numExcluded++;
        }
      }
      if (lrr_sd && qcMatrix[lrr_sdIndex][i] > EXCLUDE_LRR_SD_THRESHOLD) {
        if (addToExcludes(i, "LRR_SD > " + EXCLUDE_LRR_SD_THRESHOLD)) {
          numExcluded++;
        }
      }
      if (callrate && lrr_sd && qcMatrix[callrateIndex][i] < CONTAM_CALLRATE_THRESHOLD
          && qcMatrix[lrr_sdIndex][i] < CONTAM_LRR_SD_THRESHOLD) {
        if (addToExcludes(i, "Possible Contamination (LRR_SD < " + CONTAM_LRR_SD_THRESHOLD
            + ", CALLRATE < " + CONTAM_CALLRATE_THRESHOLD + ")")) {
          numExcluded++;
        }
      }
    }
    uses = Array.booleanNegative(excludes);
    useNotes = new String[samples.length];
    for (int i = 0; i < useNotes.length; i++) {
      useNotes[i] = uses[i] ? "." : "Excluded";
    }
    return numExcluded;
  }

  public void addPCsToSampleData(int numQ, int numPCs, boolean justQuantiles) {
    proj.getLog().reportTimeInfo("Adding " + numPCs + " to sample data");
    PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
    double[][] pcBasisSubset = new double[numPCs][];
    String[] pcTitles = new String[numPCs];
    if (pcResiduals.isSortedByProject()) {
      for (int i = 0; i < pcTitles.length; i++) {
        pcTitles[i] = "PC" + (i + 1);
        pcBasisSubset[i] = pcResiduals.getBasisAt((i + 1));
      }
      Quantiles[] quantiles =
          Quantiles.qetQuantilesFor(numQ, pcBasisSubset, pcTitles, proj.getLog());
      Hashtable<String, String> hashtable = new Hashtable<String, String>();
      for (int i = 0; i < proj.getSamples().length; i++) {
        String qcInfo = developMetricsLine(i, quantiles, pcBasisSubset, justQuantiles);
        hashtable.put(proj.getSamples()[i], qcInfo);
      }
      String[] header = developMetricsHeader(quantiles, pcTitles, numQ, justQuantiles);
      appendToSampleData(proj, hashtable, header, numQ, justQuantiles);
    } else {
      proj.getLog()
          .reportTimeError("PCs are not sorted by project, currently this is not supported");
    }
  }

  private boolean addPedigreeData() {
    Pedigree ped = proj.loadPedigree();
    if (ped == null) {
      proj.getLog().reportTimeError("No pedigree file found - cannot map samples to FID/IID");
      return false;
    }
    fidiids = new String[samples.length];
    fidiidToIndex = new HashMap<String, Integer>();
    mzTwinIds = Array.stringArray(samples.length, ".");
    HashMap<String, Integer> sampleIndices = new HashMap<String, Integer>();
    for (int i = 0; i < samples.length; i++) {
      sampleIndices.put(samples[i], i);
    }
    for (int i = 0; i < ped.getIDs().length; i++) {
      try {
        String fidiid = ped.getFID(i) + "\t" + ped.getIID(i);
        String mzTwinId = ped.getMzTwinId(i);
        int sampleIndex = sampleIndices.get(ped.getiDNA(i));

        fidiids[sampleIndex] = fidiid;
        fidiidToIndex.put(fidiid, sampleIndex);
        if (mzTwinId != null) {
          mzTwinIds[sampleIndex] = mzTwinId;
        }
      } catch (NullPointerException npe) {
        proj.getLog().reportTimeError(
            "Sample " + ped.getiDNA(i) + " exists in Pedigree but not in Project samples");
      }
    }
    return true;
  }

  public void addQCsToSampleData(int numQ, boolean justQuantiles) {
    addQCsToSampleData(numQ, 0, justQuantiles);
  }

  public void addQCsToSampleData(int numQ, int numPCs, boolean justQuantiles) {
    Quantiles[] quantiles = Quantiles.qetQuantilesFor(numQ, qcMatrix, qctitles, proj.getLog());
    Hashtable<String, String> hashtable = developHash(quantiles, justQuantiles);
    String[] header = developHeader(quantiles, numQ, justQuantiles);
    appendToSampleData(proj, hashtable, header, numQ, justQuantiles);
  }


  private boolean addToExcludes(int sampleIndex, String reason) {
    if (!excludes[sampleIndex]) {
      excludes[sampleIndex] = true;
      excludeNotes[sampleIndex] = reason;
      return true;
    }
    excludeNotes[sampleIndex] += ("; " + reason);
    return false;
  }

  private void addToMatrix(int sampleIndex, int qcTitleIndex, double data) {
    qcMatrix[qcTitleIndex][sampleIndex] = data;
    qcsAdded++;
    if (qcsAdded > qctitles.length * samples.length) {
      proj.getLog().reportTimeError("Internal Error; too many QC metrics have been added");
    }
  }

  private void addToMatrix(String sample, int qcTitleIndex, double data) {
    int sampleIndex = ext.indexOfStr(sample, samples);
    addToMatrix(sampleIndex, qcTitleIndex, data);
  }

  private void addUse_cnvs() {
    int lrr_sdIndex = getLrr_sdIndex();
    proj.getLog().reportTime("Calculating use_cnvs column");
    int numLowIntensity = 0;
    use_cnvs = uses.clone();
    use_cnvNotes = new String[use_cnvs.length];
    for (int i = 0; i < use_cnvs.length; i++) {
      if (!use_cnvs[i]) {
        use_cnvNotes[i] = "Not in Use";
      } else if (qcMatrix[lrr_sdIndex][i] > proj.LRRSD_CUTOFF.getValue()) {
        if (removeFromUse_cnvs(i, "Low Quality Intensity")) {
          numLowIntensity++;
        }
      } else {
        use_cnvNotes[i] = ".";
      }
    }

    proj.getLog().reportTime("Calculated use_cnvs column, dropping " + numLowIntensity
        + " low quality intensity samples");
  }

  public boolean correctFidIidCols() {
    if (fidiids == null) {
      return false;
    }
    Hashtable<String, String> fidiidhash = new Hashtable<String, String>();
    for (int i = 0; i < samples.length; i++) {
      fidiidhash.put(samples[i], fidiids[i]);
    }
    SampleData sampledata = proj.getSampleData(0, false);
    String[] header = new String[] {"FID", "IID"};
    proj.getLog()
        .reportTimeInfo("Replacing FID/IID columns in Sample Data with FID/IID from Pedigree");
    if (sampledata.replaceData(fidiidhash, "DNA", header, "\t", proj.getLog())) {
      proj.getLog()
          .reportTimeInfo("Replaced FID/IID columns in Sample Data with FID/IID from Pedigree");
      HashSet<String> sampleDataFidIids = new HashSet<String>();
      int nonUnique = 0;
      for (String fidiid : fidiidhash.values()) {
        if (!sampleDataFidIids.add(fidiid)) {
          nonUnique++;
        }
      }
      if (nonUnique > 0) {
        proj.getLog()
            .reportTimeWarning(nonUnique + " non-unique FID/IID pairs now exist in Sample Data!");
      }
      return true;
    }
    proj.getLog().reportTimeError("Replacing FID/IID columns in Sample Data failed");
    return false;
  }

  /**
   * Prepares a hash to be used in the addition to sample data
   */
  private Hashtable<String, String> developHash(Quantiles[] quantiles, boolean justQuantiles) {
    Hashtable<String, String> hashtable = new Hashtable<String, String>();
    for (int i = 0; i < samples.length; i++) {
      String qcInfo = (excludes[i] ? "1" : "0") + "\t" + excludeNotes[i];
      if (checkDuplicates) {
        qcInfo += "\t" + duplicateIds[i];
        qcInfo += "\t" + (uses[i] ? "1" : "0") + "\t" + useNotes[i];
        qcInfo += "\t" + (use_cnvs[i] ? "1" : "0") + "\t" + use_cnvNotes[i];
      }
      if (fidiids != null) {
        qcInfo += "\t" + mzTwinIds[i];
      }
      qcInfo += "\t" + developMetricsLine(i, quantiles, qcMatrix, justQuantiles);
      hashtable.put(samples[i], qcInfo);
    }
    return hashtable;
  }

  private String[] developHeader(Quantiles[] quantiles, int numQ, boolean justQuantiles) {
    ArrayList<String> header = new ArrayList<String>();
    header.add("CLASS=Exclude");
    header.add("ExcludeNote");
    if (checkDuplicates) {
      header.add("DuplicateId");
      header.add("Use");
      header.add("UseNote");
      header.add("Use_cnv");
      header.add("Use_cnvNote");
    }
    if (fidiids != null) {
      header.add("mzTwinID");
    }
    return Array.combine(header.toArray(new String[] {}),
        developMetricsHeader(quantiles, qctitles, numQ, justQuantiles));
  }

  private int getCallrateIndex() {
    return ext.indexOfStr("Genotype_callrate", qctitles);
  }

  public double[] getDataFor(String qctitle) {
    int indexToExtract = ext.indexOfStr(qctitle, qctitles);
    if (indexToExtract < 0) {
      proj.getLog().reportTimeError("Invalid title " + qctitle + " returning null");
      return null;
    } else {
      return qcMatrix[indexToExtract];
    }

  }

  private int getLrr_sdIndex() {
    return gcCorrectedLrrSd ? ext.indexOfStr("LRR_SD_Post_Correction", getQctitles())
        : ext.indexOfStr("LRR_SD", getQctitles());
  }

  public double[][] getQcMatrix() {
    return qcMatrix;
  }

  public String[] getQctitles() {
    return qctitles;
  }

  public String[] getSamples() {
    return samples;
  }

  private void removeDuplicatesFromUse(HashMap<String, HashSet<Integer>> duplicateSets) {
    proj.getLog().reportTimeInfo("Choosing duplicates to drop");

    int callrateIndex = getCallrateIndex();
    int lrr_sdIndex = getLrr_sdIndex();

    if (callrateIndex == -1 && lrr_sdIndex == -1) {
      proj.getLog().reportTimeWarning(
          "No callrates or LRR SDs found, duplicates will not be preferentially selected");
    }

    for (HashSet<Integer> dupeSet : duplicateSets.values()) {
      for (int index : dupeSet) {
        if (!uses[index] || !mzTwinIds[index].equals(".")) {
          continue;
        }
        for (int compIndex : dupeSet) {
          if (compIndex == index || !uses[compIndex]) {
            continue;
          }
          if (!mzTwinIds[compIndex].equals(".")) {
            // Always keep MZ Twins
            if (removeFromUses(index, "Duplicate of sample marked as MZ Twin")) {
              proj.getLog().reportTimeWarning(
                  "Sample " + samples[index] + " is duplicate of MZ Twin and was dropped");
            }
          } else if (callrateIndex != -1 && !Double.isNaN(qcMatrix[callrateIndex][index])
              && !Double.isNaN(qcMatrix[callrateIndex][compIndex])
              && qcMatrix[callrateIndex][index] != qcMatrix[callrateIndex][compIndex]) {
            // Otherwise, if samples have non-equal callrates, keep higher callrate
            if (qcMatrix[callrateIndex][index] < qcMatrix[callrateIndex][compIndex]) {
              removeFromUses(index, "Low callrate Duplicate");
            } else {
              removeFromUses(compIndex, "Low callrate Duplicate");
            }
          } else if (lrr_sdIndex != -1 && !Double.isNaN(qcMatrix[lrr_sdIndex][index])
              && !Double.isNaN(qcMatrix[lrr_sdIndex][compIndex])
              && qcMatrix[lrr_sdIndex][index] != qcMatrix[lrr_sdIndex][compIndex]) {
            // Otherwise, if samples have non-equal LRR_SDs, keep lower LRR_SD
            if (qcMatrix[lrr_sdIndex][index] > qcMatrix[lrr_sdIndex][compIndex]) {
              removeFromUses(index, "High LRR_SD Duplicate");
            } else {
              removeFromUses(compIndex, "High LRR_SD Duplicate");
            }
          } else {
            removeFromUses(compIndex, "Arbitrarily selected duplicate");
          }
        }
      }
    }
    proj.getLog().reportTimeInfo("Finished choosing duplicates to drop");
  }

  private int removeEmptyMetrics() {
    if (!verify()) {
      proj.getLog().reportTimeError("Could not verify that all data was added");
      return -1;
    }

    boolean[] useQCs = new boolean[getQcMatrix().length];
    Arrays.fill(useQCs, false);
    int numFiltered = getQcMatrix().length;
    for (int i = 0; i < getQcMatrix().length; i++) {
      for (int j = 0; j < getQcMatrix()[i].length; j++) {
        if (!Double.isNaN(getQcMatrix()[i][j])) {
          useQCs[i] = true;
          numFiltered--;
          break;
        }
      }
    }

    qcMatrix = Array.subArray(getQcMatrix(), useQCs);
    qctitles = Array.subArray(getQctitles(), useQCs);
    qcsAdded = qctitles.length * samples.length;

    if (gcCorrectedLrrSd && getLrr_sdIndex() == -1) {
      proj.getLog().reportTimeError(
          "GC corrected LRR SD is missing, run again without selecting GC corrected LRR SD or after regenerating Sample QC file with GC correction");
      return -1;
    }

    return numFiltered;
  }

  private boolean removeFromUse_cnvs(int sampleIndex, String reason) {
    if (use_cnvs[sampleIndex]) {
      use_cnvs[sampleIndex] = false;
      use_cnvNotes[sampleIndex] = reason;
      return true;
    }
    use_cnvNotes[sampleIndex] += ("; " + reason);
    return false;
  }

  private boolean removeFromUses(int sampleIndex, String reason) {
    if (uses[sampleIndex]) {
      uses[sampleIndex] = false;
      useNotes[sampleIndex] = reason;
      return true;
    }
    useNotes[sampleIndex] += ("; " + reason);
    return false;
  }

  public boolean verify() {
    return qcsAdded == qctitles.length * samples.length;
  }
}
