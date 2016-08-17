package org.genvisis.cnv.var;

import java.awt.GraphicsEnvironment;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.SwingUtilities;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.common.Aliases;
import org.genvisis.common.Array;
import org.genvisis.common.CountVector;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariantHash;
import org.genvisis.filesys.Segment;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

public class SampleData {
  /**
   * A helper class to facilitate creating sample data from either .ped or Sample_Map.csv formats
   * 
   */
  public static class Individual {
    private final String fid;
    private final String iid;
    private final String fa;

    private final String mo;
    private final String sex;
    private final String aff;
    private String dna;
    private String[] sampleMapLine;
    private String[] sampleMapHeader;

    /**
     * @param indices of the required columns
     * @param sampleMapLine a line from the sampleMap file
     * @param header
     */
    public Individual(int[] indices, String[] sampleMapLine, String[] header) {
      fid = sampleMapLine[indices[1]];
      iid = sampleMapLine[indices[2]];
      fa = "NA";
      mo = "NA";
      sex = indices[3] == -1 ? "NA" : parseSex(sampleMapLine[indices[3]]);
      dna = sampleMapLine[indices[2]];
      aff = "NA";
      sampleMapHeader = header;
      this.sampleMapLine = sampleMapLine;
    }

    public Individual(String dna, String fid, String iid, String fa, String mo, String sex,
                      String aff) {
      this.fid = fid;
      this.iid = iid;
      this.fa = fa;
      this.mo = mo;
      this.sex = sex;
      this.aff = aff;
      this.dna = dna;
    }

    public String getDna() {
      return dna;
    }

    public String getSampDataFormat() {
      if (sampleMapHeader == null) {
        return dna + "\t" + fid + "\t" + iid + "\t" + fa + "\t" + mo + "\t" + sex + "\t" + aff;
      } else {
        return dna + "\t" + fid + "\t" + iid + "\t" + fa + "\t" + mo + "\t" + sex + "\t" + aff
               + "\t" + Array.toStr(sampleMapLine);
      }
    }

    public String[] getSampleMapHeader() {
      return sampleMapHeader;
    }

    public String parseSex(String sex) {
      String s = "-1";
      if (sex.toLowerCase().equals(MitoPipeline.SEX[0])) {
        s = "2";
      } else if (sex.toLowerCase().equals(MitoPipeline.SEX[1])) {
        s = "1";
      }
      return s;
    }

    public void setDna(String dna) {
      this.dna = dna;
    }

  }

  public static final String HEATMAP = "Heat map";
  public static final String GENOTYPE = "Genotype";
  public static final String[] BASIC_CLASSES = {"All", HEATMAP, GENOTYPE};
  public static final String[] MINIMAL_SAMPLE_DATA_HEADER = {"DNA", "FID", "IID"};
  public static final String[] EUPHEMISMS = {"CleanedSex", "Sex", "Gender",};
  public static final String[] EXCLUDE_ALIASES = {"Exclude"};
  public static final String[][][] KEYS_FOR_BASIC_CLASSES =
      {{{"0", "All"}}, {{"1", "A/A"}, {"2", "A/B"}, {"3", "B/B"}},
       {{"1", "A/A"}, {"2", "A/B"}, {"3", "B/B"}}, {{"1", "A/A"}, {"2", "A/B"}, {"3", "B/B"}}};

  // public static final String[][] LINKERS = {
  // //TODO - Rohit: Removed Sample from first Linker. Confirm with Nathan if this is okay.
  // {"IndividualID", "ID", "IID", "UID", "UniqueID", "IndID"},
  // {"Family ID", "FamID", "FID"},
  // {"DNA/Sample", "DNA", "DNA#", "Sample", "LabID"},
  // {"MarkerName", "Marker", "SNP", "Variant", "VariantName"}, // will link to Scatter Plot
  // {"Region", "UCSC", "Band", "Arm"}, // will link to Trailer
  // {"Chromosome", "Chr"}, // secondary link to Trailer
  // {"Position", "Pos", "Start", "Begin"}, // secondary link to Trailer
  // {"Stop Position", "Stop", "End"} // secondary link to Trailer
  // };
  public static final String[][] LINKERS =
      {Aliases.INDIVIDUAL_ID, Aliases.FAMILY_ID, Aliases.DNA, Aliases.MARKER_NAMES, Aliases.REGION,
       Aliases.CHRS, Array.combine(Aliases.POSITIONS, Aliases.POSITIONS_START),
       Aliases.POSITIONS_STOP};
  // Hashtable<String, Integer> linkKeyIndex;
  // Hashtable<String, ArrayList<Integer>> colorKeyIndex;
  public static final int IID_INDEX_IN_LINKERS = 0;
  public static final int FID_INDEX_IN_LINKERS = 1;
  public static final int DNA_INDEX_IN_LINKERS = 2;
  private static final String NO_VALUE_FOUND = ".";

  // public static final String[] BASIC_FILTERS = {"GC"};

  public static final String PLINK_CLASS_PREFIX = "PLINK: ";

  private static String backupSampleData(Project proj, Logger log) {
    String sampleDatafilename = proj.SAMPLE_DATA_FILENAME.getValue();
    String bakDir = proj.BACKUP_DIRECTORY.getValue(true, true);
    try {
      Thread.sleep(100);
    } catch (InterruptedException e1) {
      log.reportException(e1);
    }

    String bakFile =
        bakDir + ext.removeDirectoryInfo(sampleDatafilename) + "_" + System.currentTimeMillis();
    // fail if can't back up file
    if (!Files.copyFile(sampleDatafilename, bakFile) || !new File(bakFile).exists()) {
      log.reportError("Error - could not back up sample data file " + sampleDatafilename
                      + ", will not add data");
      return null;
    }
    log.report("Info - backed up sample data file (" + sampleDatafilename + ") to " + bakFile);

    return bakFile;
  }

  /**
   * Creates a minimal sample data file using sample names. Sample names are assigned for all
   * entries in {@link SampleData#MINIMAL_SAMPLE_DATA_HEADER} (currently DNA,FID,IID)
   * <p>
   * Warning - if a sample data file already exists for the project, a new one will not be created
   * <p>
   * Warning - samples must be parsed prior to using this method as sample names are retrieved with
   * {@link Project#getSamples}
   * 
   * @return true if a minimal sample data file was created, false if not
   */
  public static boolean createMinimalSampleData(Project proj) {
    boolean created = false;
    String sampleDatafilename = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
    Logger log = proj.getLog();

    if (Files.exists(sampleDatafilename)) {
      log.reportError("Error - Sample data file " + sampleDatafilename
                      + " already exists, will not create a new one");
    } else {
      String[] samples = proj.getSamples();
      if (samples == null || samples.length < 1) {
        log.reportError("Error - Could not generate sample data file " + sampleDatafilename
                        + ", samples have not been parsed");
      } else {
        try {
          PrintWriter writer = new PrintWriter(new FileWriter(sampleDatafilename));
          writer.println(Array.toStr(MINIMAL_SAMPLE_DATA_HEADER) + "\tCLASS=" + EUPHEMISMS[1]);
          for (String sample : samples) {
            for (int j = 0; j < MINIMAL_SAMPLE_DATA_HEADER.length; j++) {
              writer.print(sample + (j < (MINIMAL_SAMPLE_DATA_HEADER.length - 1) ? "\t" : ""));
            }
            writer.print("\t-1");
            writer.println();
          }
          writer.close();
          created = true;
        } catch (Exception e) {
          log.reportError("Error writing to " + sampleDatafilename);
          log.reportException(e);
          created = false;
        }
      }
    }
    return created;
  }

  /**
   * @param pedFile ped format file to create sample data,
   * @param sampleMapCsv sample_Map.csv format file to create sample data
   * @param proj an existing, or newly created project
   * @param log Note: if the pedFile and sampleMapCsv file are both null, we create a minimal sample
   *        data instead Note: if sample data already exists, we leave it alone
   * @throws Elision
   */
  public static int createSampleData(String pedFile, String sampleMapCsv,
                                     Project proj) throws Elision {
    String sampleDataFilename = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
    if ((sampleMapCsv == null || "".equals(sampleMapCsv) || !Files.exists(sampleMapCsv))
        && (pedFile == null || "".equals(pedFile) || !Files.exists(pedFile))
        && !Files.exists(sampleDataFilename)) {
      proj.getLog().report(
                           "Neither a sample manifest nor a sample map file was provided; generating sample data file at: "
                           + sampleDataFilename);
      createMinimalSampleData(proj);
      return 0;
    } else if (!Files.exists(sampleDataFilename)) {
      if (pedFile != null && !"".equals(pedFile) && Files.exists(pedFile)) {
        generateSampleDataPed(proj, pedFile);
        return 0;
      } else {
        generateSampleDataMap(proj, sampleMapCsv);
        return 0;
      }
    } else {
      proj.getLog().report("Detected that a sampleData file already exists at " + sampleDataFilename
                           + ", skipping sampleData creation");
      return -1;
    }
  }

  public static int[] determineKeyIndices(String filename) {
    String[] header;
    header = Files.getHeaderOfFile(filename, null);
    int[] linkKeyIndices = ext.indexFactors(LINKERS, header, false, true, false, null, false);

    if (linkKeyIndices[0] == -1) {
      System.out.println("ID linker not automatically identified for file '" + filename
                         + "'; assuming the first column.");
      linkKeyIndices[0] = 0;
    }
    return linkKeyIndices;
  }

  /**
   * We use the Individual class as input so that we only need one method to generate the sample
   * data
   */
  private static void generateSampleData(Project proj, Individual[] inds) {
    // String sampleDataFile = proj.PROJECT_DIRECTORY.getValue() +
    // proj.getProperty(proj.SAMPLE_DATA_FILENAME);
    String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
    Logger log = proj.getLog();

    try {
      PrintWriter writer = new PrintWriter(new FileWriter(sampleDataFile));
      String[] classed = MitoPipeline.PED_INPUT;
      classed[5] = "Class=Sex";
      writer.println(Array.toStr(classed)
                     + (inds[0].getSampleMapHeader() == null ? ""
                                                             : "\t"
                                                               + Array.toStr(inds[0].getSampleMapHeader())));
      for (Individual ind : inds) {
        writer.println(ind.getSampDataFormat());
      }
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + sampleDataFile
                      + "\" could not be written to (it's probably open)");
      log.reportException(fnfe);
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + sampleDataFile + "\"");
      log.reportException(ioe);
    }
  }

  private static void generateSampleDataMap(Project proj, String sampleMapCsv) throws Elision {
    SampleData.generateSampleData(proj, SampleData.loadSampleMapFile(proj.getSamples(),
                                                                     sampleMapCsv, proj.getLog()));
  }

  private static void generateSampleDataPed(Project proj, String pedFile) throws Elision {
    SampleData.generateSampleData(proj, SampleData.loadPedInputFile(proj.getSamples(), pedFile,
                                                                    proj.getLog()));
  }

  private static Individual[] loadPedInputFile(String[] samples, String pedFile,
                                               Logger log) throws Elision {
    String[] line;
    ArrayList<Individual> al = new ArrayList<Individual>();

    HashSet<String> sampSet = new HashSet<String>();
    for (String sa : samples) {
      sampSet.add(sa);
    }

    ArrayList<String> foundMissing = new ArrayList<String>();

    try {
      BufferedReader reader = Files.getReader(pedFile, false, true, false);
      String temp = reader.readLine().trim();
      String delim = ext.determineDelimiter(temp);
      line = temp.split(delim);
      int[] indices = ext.indexFactors(MitoPipeline.PED_INPUT, line, true, false);
      boolean allElementsMissing = true;
      for (int ind : indices) {
        if (ind >= 0) {
          allElementsMissing = false;
          break;
        }
      }
      Individual indi;
      if (!allElementsMissing) {
        for (int indice : indices) {
          if (indice < 0) {
            log.reportError("Error - Improper formatting of the pedigree file, cannot generate sampleData");
            log.reportError("Warning - Parsing can proceed, but a sample data file is needed to generate principal components");
            throw new Elision("Improper File Format");
          }
        }
      } else {
        log.report("No header elements found, assuming no header and standard pedigree.dat format...");
        indices = new int[] {6, 0, 1, 2, 3, 4, 5};
        indi =
            new Individual(line[indices[0]], line[indices[1]], line[indices[2]], line[indices[3]],
                           line[indices[4]], line[indices[5]], line[indices[6]]);
        al.add(indi);
        if (!sampSet.remove(indi.iid)) {
          foundMissing.add(indi.iid);
        }
      }
      temp = null;
      while ((temp = reader.readLine()) != null) {
        line = temp.trim().split(delim);
        indi =
            new Individual(line[indices[0]], line[indices[1]], line[indices[2]], line[indices[3]],
                           line[indices[4]], line[indices[5]], line[indices[6]]);
        al.add(indi);
        if (!sampSet.remove(indi.iid)) {
          foundMissing.add(indi.iid);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + pedFile + "\" not found in current directory");
      throw new Elision("Specified File is Missing");
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + pedFile + "\"");
      throw new Elision("Generic Error Reading File");
    }
    if (foundMissing.size() > 0) {
      log.reportError("MAJOR ERROR - found " + foundMissing.size()
                      + " IIDs in Pedigree file that do not have data files.");
      log.reportError("\tMissing IIDs: " + foundMissing.toString());
    }
    if (sampSet.size() > 0) {
      log.reportError("MAJOR ERROR - found " + sampSet.size()
                      + " IIDs with data files not listed in Pedigree file.");
      log.reportError("\tMissing IIDs: " + sampSet.toString());
    }
    if (foundMissing.size() > 0 || sampSet.size() > 0) {
      log.reportError("MAJOR ERROR - please fix your Pedigree file to represent all individuals with data files.");
      throw new Elision("Missing Individuals Either in File or Data");
    }

    return al.toArray(new Individual[al.size()]);
  }

  private static Individual[] loadSampleMapFile(String[] samples, String sampleMapCsv,
                                                Logger log) throws Elision {
    String[] line;
    ArrayList<Individual> al = new ArrayList<Individual>();
    log.report("Using Sample Map file " + sampleMapCsv);

    HashSet<String> sampSet = new HashSet<String>();
    for (String sa : samples) {
      sampSet.add(sa);
    }

    String delim = ",";
    ArrayList<String> foundMissing = new ArrayList<String>();
    try {
      BufferedReader reader = Files.getReader(sampleMapCsv, false, true, false);
      line = reader.readLine().trim().split(delim);
      String[] header = line;
      int[] indices = ext.indexFactors(MitoPipeline.SAMPLEMAP_INPUT, line, true, false);
      if (indices[1] == -1 || indices[2] == -1) {
        log.reportError("Error - Columns \"" + MitoPipeline.SAMPLEMAP_INPUT[1] + "\" and \""
                        + MitoPipeline.SAMPLEMAP_INPUT[2] + "\" must be provided in .csv format "
                        + sampleMapCsv);
        throw new Elision("Improper File Format");
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split(delim);
        Individual indi = new Individual(indices, line, header);
        al.add(indi);
        if (!sampSet.remove(indi.iid)) {
          foundMissing.add(indi.iid);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + sampleMapCsv + "\" not found in current directory");
      throw new Elision("Specified File is Missing");
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + sampleMapCsv + "\"");
      throw new Elision("Generic Error Reading File");
    }
    if (foundMissing.size() > 0) {
      log.reportError("MAJOR ERROR - found " + foundMissing.size()
                      + " IIDs in Sample_Map file that do not have data files.");
      log.reportError("\tMissing IIDs: " + foundMissing.toString());
    }
    if (sampSet.size() > 0) {
      log.reportError("MAJOR ERROR - found " + sampSet.size()
                      + " IIDs with data files not listed in Sample_Map file.");
      log.reportError("\tMissing IIDs: " + sampSet.toString());
    }
    if (foundMissing.size() > 0 || sampSet.size() > 0) {
      log.reportError("MAJOR ERROR - please fix your Sample_Map file to represent all individuals with data files.");
      throw new Elision("Missing Individuals Either in File or Data");
    }

    return al.toArray(new Individual[al.size()]);
  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = null;
    String ped = null;
    String samp = null;

    String usage = "\n" + "cnv.var.SampleData requires 1-3 arguments\n"
                   + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + " OR\n" + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + "   (2a) OPTIONAL: .ped file (i.e. ped=data.ped (not the default))\n"
                   + "      OR\n"
                   + "   (2b) OPTIONAL: sampleMap.csv file (i.e. sampleMap=sampleMap.csv (not the default))\n"
                   + "\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("ped=")) {
        ped = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("sampleMap=")) {
        samp = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    Project thisProject = new Project(filename, false);

    if (Files.exists(thisProject.SAMPLE_DATA_FILENAME.getValue(false, false),
                     thisProject.JAR_STATUS.getValue())) {
      thisProject.getSampleData(2, false);
    } else {
      try {
        createSampleData(ped, samp, thisProject);
      } catch (Elision e) {
        // ignore, as log output has already been written if an error is thrown
      }
    }
  }

  /**
   * Popped this warning into a new Thread for gui reporting
   */
  private static void warnMissingFiles(final Project proj, final Vector<String> missingFiles) {
    final String message =
        "The following CNV " + (missingFiles.size() > 1 ? "files do not" : "file does not")
                           + " exist and therefore " + (missingFiles.size() > 1 ? "were" : "was")
                           + " not loaded:\n     "
                           + Array.toStr(Array.toStringArray(missingFiles), "\n     ")
                           + "\n\nTo prevent this message in the future, either find "
                           + (missingFiles.size() > 1 ? "these files \n     or remove them"
                                                      : "this file \n     or remove it")
                           + " from CNV_FILENAMES under CNV Files in Project Properties";
    if (GraphicsEnvironment.isHeadless()) {
      proj.getLog().reportTimeWarning(message);
    } else {
      Runnable mR = new Runnable() {
        @Override
        public void run() {
          proj.message(message);
        }
      };
      SwingUtilities.invokeLater(mR);
    }
  }

  private Project proj;
  private String[] basicClasses;
  private String[] filters;
  private String[] covars;
  private String[] classes;
  private String[][][] classColorKeys;

  private String[] cnvClasses;

  // public Hashtable<String, Integer> getLinkKeyIndex() {
  // return linkKeyIndex;
  // }
  //
  // public void setLinkKeyIndex(Hashtable<String, Integer> linkKeyIndex) {
  // this.linkKeyIndex = linkKeyIndex;
  // }

  private String[] plinkClasses;

  // private Hashtable<String,String> sampleLookup;
  // private Hashtable<String,String> famIndLookup;
  // private Hashtable<String,String> indLookup;
  private Hashtable<String, String[]> lookup;

  private Hashtable<String, IndiPheno> sampleHash;

  private boolean failedToLoad;

  private int sexClassIndex;

  private int excludeClassIndex;

  private boolean containsDNA;

  private boolean containsFID;

  private boolean containsIID;

  private volatile boolean loadedCNVs = false;

  public SampleData(final Project proj, int numberOfBasicClassesToUse,
                    final String[] cnvFilenamesVar) {
    BufferedReader reader;
    String[] line, header;
    IntVector filterIs = new IntVector();
    IntVector covarIs = new IntVector();
    IntVector classIs = new IntVector();
    IntVector iv;
    DoubleVector dv;
    IndiPheno indi;
    int dnaIndex, famIndex, indIndex;
    // Hashtable<String,IndiPheno> sampleHash;
    String filename;
    CountVector sexCountHash;
    int[] sexValues;
    String[] ids;
    Logger log;

    this.proj = proj;
    log = proj.getLog();

    failedToLoad = true;
    final String[] cnvFilenames;
    if (cnvFilenamesVar == null) {
      cnvFilenames = new String[0];
    } else {
      cnvFilenames = cnvFilenamesVar;
    }

    containsDNA = containsFID = containsIID = true;
    // linkKeyIndex = new Hashtable<String, Integer>();
    // colorKeyIndex = new Hashtable<String, ArrayList<Integer>>();

    if (numberOfBasicClassesToUse > BASIC_CLASSES.length) {
      log.reportError("Error - selected number of basic classes to use exceeds the number defined");
      numberOfBasicClassesToUse = BASIC_CLASSES.length;
    } else if (numberOfBasicClassesToUse > 0) {
      numberOfBasicClassesToUse = BASIC_CLASSES.length;
    }
    basicClasses = new String[numberOfBasicClassesToUse];
    for (int i = 0; i < basicClasses.length; i++) {
      basicClasses[i] = BASIC_CLASSES[i];
    }

    try {
      // filename = proj.getFilename(proj.SAMPLE_DATA_FILENAME);
      filename = proj.SAMPLE_DATA_FILENAME.getValue();
      if (!Files.exists(filename, proj.JAR_STATUS.getValue())) {
        proj.message("SampleData file does not exist: " + filename);
        SampleData.createMinimalSampleData(proj);
        if (!Files.exists(filename, proj.JAR_STATUS.getValue())) {
          proj.message("Could not create a minimal sample data file at: " + filename);
          return;
        } else {
          proj.message("Generated a temporary sample data file at " + filename);
        }
      }
      reader = Files.getReader(filename, proj.JAR_STATUS.getValue(), true, true); // to do, don't
                                                                                  // kill?
      header = reader.readLine().split("\t");
      dnaIndex = ext.indexOfStr("DNA", header);
      famIndex = ext.indexOfStr("FID", header);
      indIndex = ext.indexOfStr("IID", header);
      if (dnaIndex == -1) {
        log.reportError("Error - 'DNA' was not a header in the SampleData file; assuming lookup with first column");
        dnaIndex = 0;
        containsDNA = false;
      }
      if (cnvFilenames.length > 0 && famIndex == -1) {
        log.reportError("Error - 'FID' was not a header in the SampleData file; lookup for cnv data may be inaccurate");
        // cnvFilesnames = new String[0];
      }
      if (cnvFilenames.length > 0 && indIndex == -1) {
        log.reportError("Error - 'IID' was not a header in the SampleData file; lookup for cnv data may be inaccurate");
        // cnvFilesnames = new String[0];
      }
      if (famIndex == -1) {
        log.reportError("Error - 'FID' was not a header in the SampleData file; assuming family ID is in the second column");
        famIndex = 1;
        containsFID = false;
      }
      if (indIndex == -1) {
        log.reportError("Error - 'IID' was not a header in the SampleData file; assuming individual ID is in the third column");
        indIndex = 2;
        containsIID = false;
      }
      for (int i = 1; i < header.length; i++) {
        if (header[i].toUpperCase().startsWith("FILTER=")) {
          filterIs.add(i);
        } else if (header[i].toUpperCase().startsWith("CLASS=")) {
          classIs.add(i);
        } else if (header[i].toUpperCase().startsWith("COVAR=")) {
          covarIs.add(i);
        }
      }
      filters = new String[filterIs.size()];
      for (int i = 0; i < filters.length; i++) {
        filters[i] = header[filterIs.elementAt(i)].split("=")[1];
      }
      covars = new String[covarIs.size()];
      for (int i = 0; i < covars.length; i++) {
        covars[i] = header[covarIs.elementAt(i)].split("=")[1];
      }
      classes = new String[classIs.size()];
      classColorKeys = new String[classIs.size()][][];
      for (int i = 0; i < classes.length; i++) {
        line = header[classIs.elementAt(i)].split(";");
        classes[i] = line[0].split("=")[1];
        classColorKeys[i] = new String[line.length - 1][];
        for (int j = 1; j < line.length; j++) {
          classColorKeys[i][j - 1] = line[j].split("=");
          if (classColorKeys[i][j - 1].length != 2) {
            log.reportError("Error - invalid key for class '" + classes[i]
                            + "'; must use format #=Key (not '" + line[j]
                            + "'), separated by semicolons");
            classColorKeys[i][j - 1] = new String[0];
          }
        }
      }

      sexClassIndex = ext.indexFactors(new String[][] {EUPHEMISMS}, classes, true, false, true,
                                       log.getLevel() >= 1 ? true : false, log, false)[0];
      excludeClassIndex = ext.indexFactors(new String[][] {EXCLUDE_ALIASES}, classes, false, false,
                                           true, log.getLevel() >= 1 ? true : false, log, false)[0];
      log.report("Class list: " + Array.toStr(classes), true, true, 1);

      sexCountHash = new CountVector();
      sampleHash = new Hashtable<String, IndiPheno>();
      lookup = new Hashtable<String, String[]>();
      while (reader.ready()) {
        line = reader.readLine().split("\t", -1);
        indi = new IndiPheno();

        ids = new String[] {line[dnaIndex], line[famIndex] + "\t" + line[indIndex], line[indIndex]};
        lookup.put(line[dnaIndex].toLowerCase(), ids);
        lookup.put(line[famIndex].toLowerCase() + "\t" + line[indIndex].toLowerCase(), ids);
        lookup.put(line[indIndex].toLowerCase(), ids);

        dv = new DoubleVector();
        for (int i = 0; i < filterIs.size(); i++) {
          dv.add(Double.parseDouble(line[filterIs.elementAt(i)]));
        }
        indi.setFilters(Doubles.toArray(dv));

        dv = new DoubleVector();
        for (int i = 0; i < covarIs.size(); i++) {
          dv.add(ext.isMissingValue(line[covarIs.elementAt(i)]) ? Double.NaN
                                                                : Double.parseDouble(line[covarIs.elementAt(i)]));
        }
        indi.setCovars(Doubles.toArray(dv));

        iv = new IntVector();
        for (int i = 0; i < classIs.size(); i++) {
          iv.add(ext.isMissingValue(line[classIs.elementAt(i)])
                 || Integer.parseInt(line[classIs.elementAt(i)]) < 0 ? Integer.MIN_VALUE
                                                                     : Integer.parseInt(line[classIs.elementAt(i)]));
        }
        indi.setClasses(Ints.toArray(iv));
        if (sexClassIndex != -1) {
          sexCountHash.add(indi.getClasses()[sexClassIndex] + "");
        }

        sampleHash.put(line[dnaIndex].toLowerCase(), indi);
      }
      reader.close();

      if (log.getLevel() >= 1) {
        if (sexClassIndex != -1) {
          sexValues = Array.toIntArray(sexCountHash.getValues());
          sexValues = Sort.putInOrder(sexValues, Sort.quicksort(sexValues, Sort.DESCENDING));
          if (Array.countIf(sexValues, 2) == 0) {
            log.reportError("Warning - no females listed in SampleData file; make sure 1=male and 2=female in the coding");
            // proj.message("descending "+ Array.toStr(sexValues, " ")+"\tError - warning no females
            // listed in SampleData file; make sure 1=male and 2=female in the coding");
          }
        } else {
          proj.message("Error - variable names '" + Array.toStr(EUPHEMISMS, "/")
                       + "' was not found in the SampleData file; make sure 1=male and 2=female in the coding");
        }
      }
    } catch (FileNotFoundException fnfe) {
      // System.err.println("Error: file \""+proj.getFilename(proj.SAMPLE_DATA_FILENAME)+"\" not
      // found in current directory");
      System.err.println("Error: file \"" + proj.SAMPLE_DATA_FILENAME.getValue()
                         + "\" not found in current directory");
      // System.exit(1);
    } catch (IOException ioe) {
      // System.err.println("Error reading file
      // \""+proj.getFilename(proj.SAMPLE_DATA_FILENAME)+"\"");
      System.err.println("Error reading file \"" + proj.SAMPLE_DATA_FILENAME.getValue() + "\"");
      // System.exit(2);
    }

    cnvClasses = new String[0];
    if (cnvFilenames.length > 0) {
      loadedCNVs = false;
      Runnable cnvLoadingRunnable = new Runnable() {
        @Override
        public void run() {
          loadCNVs(cnvFilenames, proj.JAR_STATUS.getValue(), proj.getLog());
          IndiPheno.setCNVsLoaded();
          loadedCNVs = true;
        }
      };
      Thread cnvLoadingThread = new Thread(cnvLoadingRunnable);
      cnvLoadingThread.setName("CNVLoadingThread");
      cnvLoadingThread.setDaemon(true);
      cnvLoadingThread.start();
    } else {
      loadedCNVs = true;
      IndiPheno.setCNVsLoaded();

    }

    String[] plinkFilenames = proj.PLINK_DIR_FILEROOTS.getValue();
    if (plinkFilenames.length > 0) {
      loadPlinkFiles(plinkFilenames, proj.JAR_STATUS.getValue());
    } else {
      plinkClasses = new String[0];
    }

    failedToLoad = false;
  }

  /**
   * Add column(s) to sample data according to linker (linker must be a single column i.e "DNA"
   * ,"FID","IID",and not "FID\tIID") that corresponds to a single column in the sample data file,
   * new columns will be appended, not inserted.
   * <p>
   * Warning - existing matching columns are left in place with their header edited for posterity
   * <p>
   * Warning - alternate names for linker are not determined, please do that prior to use
   * 
   * @param linkData Hashtable with the linker entry as the key (DNA, , FID, IID, plate, sex etc),
   *        and data to add as value. If adding more than one column, the data to add can be
   *        delimited
   * @param linker The linker column header in sample Data, the linker should be verified prior to
   *        being used in this function i.e check that the linker is present before hand and
   *        determine the appropriate linker. Data will not be added if linker is not present in
   *        header.
   * @param columnHeaders The headers to be added. If greater than 1, the length of the header must
   *        equal the length of the (delimited) data for every addition, missing keys will be
   *        replaced with columnHeaders.length missingData values
   * @param missingData What to replace missing data with (when a key is not found)
   * @param linkDataDelimiter The delimiter to use for the data in linkData, or <code>null</code> if
   *        only adding one column of data
   * @param log a log!
   * 
   * @return true if data was added successfully, false if not
   */
  public boolean addData(Hashtable<String, String> linkData, String linker, String[] columnHeaders,
                         String missingData, String linkDataDelimiter, Logger log) {
    boolean add = true;
    boolean writerWasOpened = false;
    BufferedReader reader;
    PrintWriter writer;
    // String sampleDatafilename = proj.getFilename(proj.SAMPLE_DATA_FILENAME);
    String sampleDatafilename = proj.SAMPLE_DATA_FILENAME.getValue();
    // in memory backup
    String[][] sampleDataMatrix =
        HashVec.loadFileToStringMatrix(sampleDatafilename, false, null, false);
    String[] sampleDataHeader = Files.getHeaderOfFile(sampleDatafilename, log);
    int linkerIndex = ext.indexOfStr(linker, sampleDataHeader);
    // fail if can't find linker
    if (linkerIndex == -1) {
      log.reportError("Error - could not find linker " + linker + " in sample data file "
                      + sampleDatafilename + ", will not add data");
      add = false;
    } else {
      String bakFile = backupSampleData(proj, log);
      if (bakFile == null) {
        add = false;
      } else {
        log.report("Info - backed up sample data file (" + sampleDatafilename + ") to " + bakFile);
        int numMissing = 0;
        int numAdded = 0;
        try {
          reader = new BufferedReader(new FileReader(bakFile));
          String[] tmpHeader = reader.readLine().trim().split("\t");
          // fail if header of original and backup do not match
          if (!Array.equals(tmpHeader, sampleDataHeader, true)) {
            log.reportError("Error - backup sample data " + bakFile
                            + " does not contain the same header as " + sampleDatafilename
                            + ", will not add data");
            add = false;
          } else {
            writerWasOpened = true;
            writer = Files.getAppropriateWriter(sampleDatafilename);
            String dateTime = ext.getDate() + "_" + ext.getTime();
            for (String columnHeader : columnHeaders) {
              int[] indices = ext.indicesOfStr(columnHeader, sampleDataHeader, false, true);
              for (int index : indices) {
                sampleDataHeader[index] = "x" + sampleDataHeader[index] + "_replaced_" + dateTime;
              }
            }
            writer.println(Array.toStr(sampleDataHeader) + "\t" + Array.toStr(columnHeaders));
            String[] blanks = new String[columnHeaders.length];
            Arrays.fill(blanks, missingData);
            String[] line;
            while (reader.ready() && add) {
              line = reader.readLine().trim().split("\t");
              if (linkData.containsKey(line[linkerIndex])) {
                if (columnHeaders.length > 1) {
                  String[] data = linkData.get(line[linkerIndex]).split(linkDataDelimiter);
                  if (data.length == columnHeaders.length) {
                    writer.println(Array.toStr(line) + "\t" + Array.toStr(data));
                    numAdded++;
                  } // fail if added data for a linker is the wrong length after delimiting
                  else {
                    log.reportError("Error - the number of entries in the column header "
                                    + Array.toStr(columnHeaders) + " (" + columnHeaders.length
                                    + ") does not equal the number of entries in the data "
                                    + Array.toStr(data) + " (" + data.length + ")");
                    log.reportError("Cancelling the addition and replacing sample data with backup");
                    add = false;
                  }
                } else {
                  String data = linkData.get(line[linkerIndex]);
                  writer.println(Array.toStr(line) + "\t" + data);
                  numAdded++;
                }
              } else {
                writer.println(Array.toStr(line) + "\t" + Array.toStr(blanks));
                numMissing++;
              }
            }
            writer.close();
          }
          reader.close();
          if (add) {
            // TODO check for new classes/filters/etc
            for (String header : columnHeaders) {
              if (header.toUpperCase().startsWith("CLASS=")) {
                classes = Array.addStrToArray(header.split("=")[1], classes);
              }
            }

            log.report("Info - added new data to sample data for " + numAdded
                       + (numAdded > 1 ? " samples " : " samples"));
            if (numMissing > 0) {
              log.report("Warning - " + numMissing + " " + (numMissing > 1 ? "samples" : "samples")
                         + " had missing data ");
            }
          } else if (!add && writerWasOpened) {
            if (!Files.copyFile(bakFile, sampleDatafilename)) {
              Files.writeMatrix(sampleDataMatrix, sampleDatafilename, "\t");
            }
          }
        } catch (FileNotFoundException e) {
          log.reportError("Error: Sample Data backup file " + bakFile + " not found");
          log.reportException(e);
          add = false;
        } catch (IOException e) {
          log.reportError("Error: unable to read sample data backup file" + bakFile);
          log.reportException(e);
          add = false;
        }
      }
    }
    return add;
  }

  public void addToSampleData(Hashtable<String, String> colorKeyValue, String recentSelectionFile,
                              int selectedColorKey) {
    String sampleDatafilename;
    BufferedReader reader;
    BufferedWriter writer;
    String[] inLineArry;
    String bakFile;
    String inLine;
    String colorKeyHeader;
    String[] keys;
    boolean covar, negativeValues, largerThanByte;
    String trav;

    sampleDatafilename = proj.SAMPLE_DATA_FILENAME.getValue(false, false);

    if (!Files.exists(sampleDatafilename, proj.JAR_STATUS.getValue())) {
      // JOptionPane.showMessageDialog(null, "Cannot add as a color key without an existing
      // SampleData file", "Error", JOptionPane.ERROR_MESSAGE);
      System.err.println("Cannot add as a color key without an existing SampleData file");
      return;
    }

    reader = null;
    writer = null;

    System.out.println("Sample data: " + sampleDatafilename);
    bakFile = proj.archiveFile(sampleDatafilename); // create backup of sample data file
    colorKeyHeader = Files.getHeaderOfFile(recentSelectionFile, null)[selectedColorKey - 1];

    covar = false;
    negativeValues = false;
    largerThanByte = false;
    keys = HashVec.getKeys(colorKeyValue, false, false);
    for (String key : keys) {
      trav = colorKeyValue.get(key);
      if (!ext.isMissingValue(trav) && !ext.isValidInteger(trav)) {
        covar = true;
      }
      if (ext.isValidDouble(trav) && Double.parseDouble(trav) < 0) {
        negativeValues = true;
      }
      if (ext.isValidDouble(trav) && Double.parseDouble(trav) > Byte.MAX_VALUE) {
        largerThanByte = true;
      }
    }

    if (covar) {
      // TODO: make these thing show up
      // JOptionPane.showMessageDialog(null, "Variable '"+colorKeyHeader+"' contains a quantitative
      // meaure and will be added as a COVAR in SampleData and not as a Class", "Warning",
      // JOptionPane.ERROR_MESSAGE);
      System.out.println("Variable '" + colorKeyHeader
                         + "' contains a quantitative measure and will be added as a COVAR in SampleData and not as a Class");
    } else if (negativeValues) {
      // JOptionPane.showMessageDialog(null, "Variable '"+colorKeyHeader+"' contains negative values
      // and will be added as a COVAR in SampleData and not as a Class", "Warning",
      // JOptionPane.ERROR_MESSAGE);
      System.out.println("Variable '" + colorKeyHeader
                         + "' contains negative values and will be added as a COVAR in SampleData and not as a Class");
      covar = true;
    } else if (largerThanByte) {
      // JOptionPane.showMessageDialog(null, "Variable '"+colorKeyHeader+"' contains values larger
      // than 128 and will be added as a COVAR in SampleData and not as a Class", "Warning",
      // JOptionPane.ERROR_MESSAGE);
      System.out.println("Variable '" + colorKeyHeader
                         + "' contains values larger than 128 and will be added as a COVAR in SampleData and not as a Class");
      covar = true;
    }

    try {
      reader = new BufferedReader(new FileReader(bakFile));
      writer = new BufferedWriter(new FileWriter(sampleDatafilename));
      inLine = reader.readLine();
      int samDataIndex = getSampleDataHeaders(inLine)[DNA_INDEX_IN_LINKERS];
      inLine = inLine + "\t" + (covar ? "Covar=" : "Class=") + colorKeyHeader;
      writer.write(inLine); // write the headers
      while (reader.ready()) {
        writer.newLine();
        inLine = reader.readLine();
        if (inLine.contains("\t")) {
          inLineArry = inLine.trim().split("\t", -1);
        } else {
          inLineArry = inLine.trim().split("[\\s]+");
        }
        if (colorKeyValue.containsKey(inLineArry[samDataIndex])) {
          inLine = inLine + "\t" + colorKeyValue.get(inLineArry[samDataIndex]);
        } else {
          inLine = inLine + "\t" + NO_VALUE_FOUND;
        }
        writer.write(inLine);
      }
    } catch (FileNotFoundException e) {
      System.out.println("Error: Sample Data backup file not found");
    } catch (IOException e) {
      System.out.println("Error: unable to read sample data backup file");
    } finally {
      closeStream(reader);
      closeStream(writer);
      // twoDPanel.paintAgain();
    }
    // reloadSampleDataUI();
    System.out.println(colorKeyHeader.split(";")[0] + " set as color key and added to Sample Data");
    // JOptionPane.showMessageDialog(null, colorKeyHeader.split(";")[0] + " set as color key and
    // added to Sample Data", "Information", JOptionPane.INFORMATION_MESSAGE);
  }

  public void closeStream(Closeable s) {
    try {
      if (s != null) {
        s.close();
      }
    } catch (IOException e) {
      // Log or rethrow as unchecked (like RuntimException) ;)
    }
  }

  public boolean containsDNA() {
    return containsDNA;
  }

  public boolean containsFID() {
    return containsFID;
  }

  public boolean containsIID() {
    return containsIID;
  }

  public Hashtable<String, String> createHashWithSampleID(Hashtable<String, String> colorKeyValue) {
    Hashtable<String, String> colorKeyValueHash;

    colorKeyValueHash = new Hashtable<String, String>();
    for (String key : colorKeyValue.keySet()) {
      colorKeyValueHash.put(lookup(key)[0], colorKeyValue.get(key));
    }

    return colorKeyValueHash;
  }

  public byte determineCodeFromClass(int currentClass, byte alleleCount, IndiPheno indi, byte chr,
                                     int position) {
    int[] classes, indices;
    CNVariant[] segs;
    int index;

    indices = getClassCategoryAndIndex(currentClass);
    switch (indices[0]) {
      case 0:
        if (basicClasses[indices[1]].equals("All")) {
          return 0;
        } else if (basicClasses[indices[1]].equals(GENOTYPE)) {
          return (byte) (alleleCount + 1);
        } else if (basicClasses[indices[1]].equals(HEATMAP)) {
          return 0;
        } else {
          System.err.println("Error - codeFromClass not defined for type: " + indices[1]);
          System.err.println("        (" + basicClasses[indices[1]] + ")");
          return 0;
        }
        // break; // DEAD CODE
      case 1:
        classes = indi.getClasses();
        if (classes[indices[1]] == Integer.MIN_VALUE) {
          return -1;
        } else {
          return (byte) classes[indices[1]];
        }
        // break; // DEAD CODE
      case 2:
        segs = indi.getCNVs(indices[1], chr);
        if (segs == null) {
          return 0;
        } else {
          index =
              Segment.binarySearchForOverlap(new Segment(chr, position - 2, position + 2), segs);
          if (index == -1) {
            return 0;
          } else {
            return (byte) (segs[index].getCN() + 1);
          }
        }
        // break; // DEAD CODE
      case 3:
        return (byte) (alleleCount + 1); // fake as genotype - we'll load PLINK genotypes elsewhere
                                         // as needed
      default:
        System.err.println("Error - invalid class index");
        return 0;
    }
  }

  public boolean failedToLoad() {
    return failedToLoad;
  }

  public String[][] getActualClassColorKey(int index) {
    return classColorKeys[index];
  }

  public String getActualClassName(int index) {
    return classes[index];
  }

  public String[] getBasicClasses() {
    return basicClasses;
  }

  public int[] getClassCategoryAndIndex(int index) {
    int[] indices = new int[2];

    if (index < basicClasses.length) {
      indices[0] = 0;
      indices[1] = index;
    } else if (index < basicClasses.length + classes.length) {
      indices[0] = 1;
      indices[1] = index - basicClasses.length;
    } else if (index < basicClasses.length + classes.length + cnvClasses.length) {
      indices[0] = 2;
      indices[1] = index - basicClasses.length - classes.length;
    } else if (index < basicClasses.length + classes.length + cnvClasses.length
                       + plinkClasses.length) {
      indices[0] = 3;
      indices[1] = index - basicClasses.length - classes.length - cnvClasses.length;
    } else {
      System.err.println("Error - invalid class index");
    }

    return indices;
  }

  public String[] getClasses() {
    return getClasses(false);
  }

  public String[] getClasses(boolean includeBasicClasses) {
    String[] result;

    if (!includeBasicClasses && classes == null) {
      result = new String[0];
    } else if (!includeBasicClasses && classes != null) {
      result = classes;
    } else if (includeBasicClasses && classes == null) {
      result = basicClasses;
    } else {
      result = new String[basicClasses.length + classes.length];
      for (int i = 0; i < basicClasses.length; i++) {
        result[i] = basicClasses[i];
      }
      for (int i = 0; i < classes.length; i++) {
        result[i + basicClasses.length] = classes[i];
      }
    }

    return result;
  }

  public byte getClassForInd(String sample, int currentClass) {
    int[] classes;
    IndiPheno indi;

    indi = sampleHash.get(sample.toLowerCase());
    if (indi == null) {
      return 0;
    }
    classes = indi.getClasses();
    if (currentClass == -1 || classes[currentClass] == Integer.MIN_VALUE) {
      return -1;
    } else {
      return (byte) classes[currentClass];
    }
  }

  public String getClassName(int index) {
    int[] indices = getClassCategoryAndIndex(index);

    switch (indices[0]) {
      case 0:
        return basicClasses[indices[1]];
      case 1:
        return classes[indices[1]];
      case 2:
        return cnvClasses[indices[1]];
      case 3:
        return plinkClasses[indices[1]];
      default:
        return null;
    }
  }

  public String[] getCnvClasses() {
    return cnvClasses;
  }

  public boolean getCNVsLoaded() {
    return loadedCNVs;
  }

  public String[] getCovars() {
    if (covars == null) {
      return new String[0];
    } else {
      return covars;
    }
  }

  public String[] getFilters() {
    if (filters == null) {
      return new String[0];
    } else {
      return filters;
    }
  }

  public IndiPheno getIndiFromSampleHash(String sampleID) {
    return sampleHash.get(sampleID.toLowerCase());
  }

  public IndiPheno getIndiPheno(String sample) {
    return sampleHash.get(sample.toLowerCase());
  }

  public String[] getListOfSamples() {
    return HashVec.getKeys(sampleHash);
  }

  public int getNumActualClasses() {
    if (classes == null) {
      return -1;
    }
    return classes.length;
  }

  public int getNumClasses() {
    return basicClasses.length + classes.length + cnvClasses.length + plinkClasses.length;
  }

  // public void initLinkKey(String filename) {
  // int[] linkKeyColumnLabels = determineKeyIndices(filename);
  // if (linkKeyColumnLabels[DNA_INDEX_IN_LINKERS] >= 0) {
  // // {"DNA/Sample", "DNA", "DNA#", "Sample", "LabID"} exists
  // linkKeyIndex.put(filename, DNA_INDEX_IN_LINKERS);
  // // JOptionPane.showMessageDialog(null, "Link is set to: " +
  // Arrays.toString(LINKERS[DNA_INDEX_IN_LINKERS]), "Information",
  // JOptionPane.INFORMATION_MESSAGE);
  // System.out.println("Link key set to: " + Arrays.toString(LINKERS[DNA_INDEX_IN_LINKERS]));
  // } else if (linkKeyColumnLabels[FID_INDEX_IN_LINKERS] >= 0) {
  // linkKeyIndex.put(filename, FID_INDEX_IN_LINKERS);
  // // JOptionPane.showMessageDialog(null, "Link is set to: " +
  // Arrays.toString(LINKERS[FID_INDEX_IN_LINKERS]), "Information",
  // JOptionPane.INFORMATION_MESSAGE);
  // System.out.println("Link key set to: " + Arrays.toString(LINKERS[FID_INDEX_IN_LINKERS]));
  // } else if (linkKeyColumnLabels[IID_INDEX_IN_LINKERS] >= 0) {
  // linkKeyIndex.put(filename, IID_INDEX_IN_LINKERS);
  // // JOptionPane.showMessageDialog(null, "Link is set to: " +
  // Arrays.toString(LINKERS[IID_INDEX_IN_LINKERS]), "Information",
  // JOptionPane.INFORMATION_MESSAGE);
  // System.out.println("Link key set to: " + Arrays.toString(LINKERS[IID_INDEX_IN_LINKERS]));
  // } else {
  // // JOptionPane.showMessageDialog(null, "Unable to initialize the link key. Please select a link
  // key manually.", "Error", JOptionPane.ERROR_MESSAGE);
  // System.out.println("Unable to initialize the link key.");
  // }
  // }

  // public void setLinkKey(String filename, int selectedLinkKey) {
  // int[] linkKeyColumnLabels;
  //
  // linkKeyColumnLabels = determineKeyIndices(filename);
  //
  // for (int i = 0; i < linkKeyColumnLabels.length; i++) {
  // if ((linkKeyColumnLabels[i] + 1) == selectedLinkKey) {
  // linkKeyIndex.put(filename, i);
  // System.out.println("Link Key set to: " + Arrays.toString(LINKERS[i]));
  // // createLinkKeyToDataHash(selectedNodes[0][0], linkKeyColumnLabels);
  // // JOptionPane.showMessageDialog(null, "Link is set to: " + Arrays.toString(LINKERS[i]),
  // "Information", JOptionPane.INFORMATION_MESSAGE);
  // System.err.println("Link is set to:" + Arrays.toString(LINKERS[i]));
  // return;
  // }
  // }
  // proj.message("Unable to set link key. Please make sure you are selecting a valid key");
  // }



  // public void setColorKey(String dataFile, int selectedColorKey) {
  // ArrayList<Integer> colorKeys;
  //
  // if (colorKeyIndex.containsKey(dataFile)) {
  // colorKeys = colorKeyIndex.get(dataFile);
  // } else {
  // colorKeyIndex.put(dataFile, colorKeys = new ArrayList<Integer>());
  // }
  // for (Integer colorKey : colorKeys) {
  // if (colorKey == selectedColorKey) {
  // System.out.println("Error: Already set as color key");
  // // TODO: Make these thing to display
  // // JOptionPane.showMessageDialog(null, "Error: Already sey as color key", "Error",
  // // JOptionPane.ERROR_MESSAGE);
  // System.err.println("Error: Already set as color key");
  // return;
  // }
  // }
  // colorKeys.add(selectedColorKey); // add to colorKeys
  // setColorKeyHandler(dataFile, selectedColorKey);
  // }

  // public void setColorKeyHandler(String filename, int selectedColorKey) {
  // Hashtable<String, String> colorKeyValue;
  // int[] linkKeyColumnLabels;
  //
  // linkKeyColumnLabels = determineKeyIndices(filename);
  // colorKeyValue = new Hashtable<String, String>();
  // if (linkKeyIndex.containsKey(filename)) {
  // switch (linkKeyIndex.get(filename)) {
  // case DNA_INDEX_IN_LINKERS:
  // colorKeyValue = HashVec.loadFileToHashString(filename, new int[] {
  // linkKeyColumnLabels[linkKeyIndex.get(filename)] }, new int[] { selectedColorKey - 1 }, false,
  // "", true, false, false);
  // break;
  // case FID_INDEX_IN_LINKERS:
  // colorKeyValue = HashVec.loadFileToHashString(filename, new int[] {
  // linkKeyColumnLabels[linkKeyIndex.get(filename)], linkKeyColumnLabels[IID_INDEX_IN_LINKERS] },
  // new int[] { selectedColorKey - 1 }, false, "", true, false, false);
  // colorKeyValue = createHashWithSampleID(colorKeyValue); // colorkey value hash with key as
  // sampleID
  // break;
  // case IID_INDEX_IN_LINKERS:
  // colorKeyValue = HashVec.loadFileToHashString(filename, new int[] {
  // linkKeyColumnLabels[linkKeyIndex.get(filename)] }, new int[] { selectedColorKey - 1 }, false,
  // "", true, false, false);
  // colorKeyValue = createHashWithSampleID(colorKeyValue); // colorkey value hash with key as
  // sampleID
  // break;
  // default:
  // System.out.println("Error: Unable to read color key values. Invalid link key.");
  // // TODO: display this
  // // JOptionPane.showMessageDialog(null, "Error: Unable to read color key values. Invalid link
  // key.", "Error", JOptionPane.ERROR_MESSAGE);
  // System.err.println("Error: Unable to read color key values. Invalid link key.");
  // break;
  // }
  // addToSampleData(colorKeyValue, filename, selectedColorKey);
  // } else {
  // System.err.println("Error: Unable to find link key index for file: " + filename + "\n Failed to
  // set the color key. Please make sure link key is set before setting color key");
  // }
  // }

  public int getNumCNVClasses() {
    if (cnvClasses == null) {
      return 0;
    }
    return cnvClasses.length;
  }

  public int getNumPLINKClasses() {
    if (plinkClasses == null) {
      return 0;
    }
    return plinkClasses.length;
  }

  public String[] getPlinkClasses() {
    return plinkClasses;
  }

  /**
   * Function to indentify the headers in the sample data file
   *
   * @param header : a string containing all the headers read as string
   */
  public int[] getSampleDataHeaders(String header) {
    String[] headersArray;
    int[] indices;

    if (header.contains("\t")) {
      headersArray = header.trim().split("\t", -1);
    } else {
      headersArray = header.trim().split("[\\s]+");
    }
    indices = ext.indexFactors(LINKERS, headersArray, false, true, false, null, false);

    if (indices[0] == -1) {
      System.err.println("ID linker not automatically identified for Sample Data. Assuming the first column.");
      indices[0] = 0;
    }
    System.out.println("The header indices in Sample data are: " + Arrays.toString(headersArray));

    return indices;
  }

  public Hashtable<String, IndiPheno> getSampleHash() {
    return sampleHash;
  }

  public int getSexClassIndex() {
    return sexClassIndex;
  }

  /**
   * @return the sex designation (1=male and 2=female) of all samples in the project, relies on
   *         {@link SampleData#getSexForIndividual(String)}
   */
  public int[] getSexForAllIndividuals() {
    String[] samples = proj.getSamples();
    int[] sexes = new int[samples.length];

    if (sexClassIndex < 0) {
      proj.getLog()
          .reportTimeWarning("Variable names '" + Array.toStr(EUPHEMISMS, "/")
                             + "' were not found in the SampleData file; make sure 1=male and 2=female in the coding to retrive sample sex");
    }
    for (int i = 0; i < samples.length; i++) {
      sexes[i] = getSexForIndividual(samples[i]);
    }
    return sexes;
  }

  public int getSexForIndividual(String id) {
    IndiPheno indi;
    String[] ids;

    indi = sampleHash.get(id.toLowerCase());
    // indi = sampleHash.get("S_"+id);
    // sampleLookup.put("FI_"+line[famIndex]+"\t"+line[indIndex], "S_"+line[dnaIndex]);
    // sampleLookup.put("I_"+line[indIndex], "S_"+line[dnaIndex]);
    // famIndLookup.put("I_"+line[indIndex], "FI_"+line[famIndex]+"\t"+line[indIndex]);
    // famIndLookup.put("S_"+line[dnaIndex], "FI_"+line[famIndex]+"\t"+line[indIndex]);
    // indLookup.put("FI_"+line[famIndex]+"\t"+line[dnaIndex], "I_"+line[indIndex]);
    // indLookup.put("S_"+line[dnaIndex], "I_"+line[indIndex]);
    if (indi == null) {
      ids = lookup.get(id.toLowerCase());
      if (ids != null) {
        indi = sampleHash.get(ids[0]);
      }
    }

    if (indi == null) {
      System.err.println("Error - id '" + id + "' was not present in the SampleData");
      return -1;
    }

    if (sexClassIndex == -1 || indi.getClasses()[sexClassIndex] == Integer.MIN_VALUE) {
      return -1;
    } else {
      return indi.getClasses()[sexClassIndex];
    }
  }

  public boolean hasExcludedIndividuals() {
    return excludeClassIndex != -1;
  }

  public boolean individualShouldBeExcluded(String id) {
    IndiPheno indi;
    String[] ids;

    if (excludeClassIndex == -1) {
      return false;
    }

    indi = sampleHash.get(id.toLowerCase());
    if (indi == null) {
      ids = lookup.get(id.toLowerCase());
      if (ids != null) {
        indi = sampleHash.get(ids[0]);
      }
    }

    if (indi == null) {
      // System.err.println("Error - id '"+id+"' was not present in the SampleData");
      return false;
    }

    if (indi.getClasses()[excludeClassIndex] == Integer.MIN_VALUE) {
      return false;
    } else {
      return indi.getClasses()[excludeClassIndex] == 1;
    }
  }

  public void loadCNVs(String[] files, boolean jar, Logger log) {
    Vector<Hashtable<String, CNVariant[]>> finalHashes;
    CNVariantHash[] cnvhs;
    IndiPheno indi;
    String[] inds;
    String trav;
    long time;
    Vector<String> missingFiles = new Vector<String>();

    for (int i = 0; i < files.length; i++) {
      if (!Files.exists(files[i], proj.JAR_STATUS.getValue())) {
        missingFiles.add(files[i]);
        files = Array.removeFromArray(files, i);
        i--;
      }
    }
    if (missingFiles.size() > 0) {
      warnMissingFiles(proj, missingFiles);
    }

    time = new Date().getTime();
    cnvClasses = new String[files.length];
    cnvhs = new CNVariantHash[files.length];
    for (int i = 0; i < files.length; i++) {
      cnvClasses[i] = ext.rootOf(files[i]);
      System.out.println(i + "\t" + cnvClasses[i]);
      cnvhs[i] = CNVariantHash.load(files[i], CNVariantHash.CONSTRUCT_BY_IND, jar, log);
    }
    log.report("Read in CNV data in " + ext.getTimeElapsed(time));

    time = new Date().getTime();
    inds = HashVec.getKeys(sampleHash);
    for (String ind : inds) {
      indi = sampleHash.get(ind);
      trav = lookup.get(ind.toLowerCase())[1];

      finalHashes = new Vector<Hashtable<String, CNVariant[]>>();
      for (int j = 0; j < files.length; j++) {
        finalHashes.add(cnvhs[j].getDataFor(trav));
      }
      indi.setCNVclasses(finalHashes);
    }
    log.report("Added CNV data in " + ext.getTimeElapsed(time));
  }

  public void loadPlinkFiles(String[] files, boolean jar) {
    String[] tempPlinkClasses = new String[files.length];
    int cnt = 0;
    for (String file : files) {
      if ("".equals(file)) {
        continue;
      }
      tempPlinkClasses[cnt] = PLINK_CLASS_PREFIX + ext.rootOf(file);
      cnt++;
    }
    plinkClasses = Array.trimArray(tempPlinkClasses);
  }

  public String[] lookup(String str) {
    return lookup.get(str.toLowerCase());
  }

  public void removeColorKey(String colorKey) {

    String sampleDatafilename = proj.SAMPLE_DATA_FILENAME.getValue();

    System.out.println("Sample data: " + sampleDatafilename);

    String[] sampeleDataHeader = Files.getHeaderOfFile(sampleDatafilename, null); // header of
                                                                                  // sample data
    int i;
    for (i = 0; i < sampeleDataHeader.length; i++) {
      String[] splitOnEquals = sampeleDataHeader[i].split("=", 2);
      if (splitOnEquals.length > 0 && splitOnEquals[0].equalsIgnoreCase("CLASS")) {
        if (splitOnEquals[1].split(";", 2)[0].equalsIgnoreCase(colorKey)) {
          // color key found at position i in header columns
          break;
        }
      }
    }
    if (i == sampeleDataHeader.length) {
      // column to be deleted was not foung in sample data
      // JOptionPane.showMessageDialog(null, "Error: Unable to find the specified column in Sample
      // Data for deletion", "Error", JOptionPane.ERROR_MESSAGE);
      System.out.println("Error: Unable to find the specified column in Sample Data for deletion");
    } else {
      // the column at i is to be deleted
      int[] colToLoad = new int[sampeleDataHeader.length - 1];
      int index = 0, col = 0;
      while (index < (sampeleDataHeader.length - 1)) {
        if (col != i) {
          colToLoad[index++] = col;
        }
        col++;
      }

      // load the sample data without the color key column which has to bed deleted
      String[][] sampleDataMatrix =
          HashVec.loadFileToStringMatrix(sampleDatafilename, false, colToLoad, false);

      String sampleDataDelimiter = Files.determineDelimiter(sampleDatafilename, null);

      String bakFile = proj.archiveFile(sampleDatafilename); // create backup of sample data file
      System.out.println("Deleting color key " + colorKey
                         + " from sample data. Sample data backup: " + bakFile);

      // write the new sample data which does not have the removed color key column

      Files.writeMatrix(sampleDataMatrix, sampleDatafilename, sampleDataDelimiter);
      // JOptionPane.showMessageDialog(null, colorKey + "deleted in sample data", "Information",
      // JOptionPane.INFORMATION_MESSAGE);
      System.err.println(colorKey + "deleted in sample data");
    }
  }

  /**
   * Replace column(s) in sample data according to linker (linker must be a single column i.e "DNA"
   * ,"FID","IID",and not "FID\tIID") that corresponds to a single column in the sample data file,
   * new columns will be merged with existing columns
   * <p>
   * Warning - alternate names for linker are not determined, please do that prior to use
   * <p>
   * Warning - missing values are left unchanged
   * 
   * @param linkData Hashtable with the linker entry as the key (DNA, , FID, IID, plate, sex etc),
   *        and data to set as value. If replacing more than one column, the data to replace can be
   *        delimited If there are entries in Sample Data missing, the existing data is added to
   *        this Hashtable and left unchanged in the Sample Data file
   * @param linker The linker column header in sample Data, the linker should be verified prior to
   *        being used in this function i.e check that the linker is present before hand and
   *        determine the appropriate linker. Data will not be replaced if linker is not present in
   *        header.
   * @param columnHeaders The columns to be replaced. If greater than 1, the length of the header
   *        must equal the length of the (delimited) data for every replacement
   * @param linkDataDelimiter The delimiter to use for the data in linkData, or <code>null</code> if
   *        only replacing one column of data
   * @param log a log!
   * 
   * @return true if data was replaced successfully, false if not
   */
  public boolean replaceData(Hashtable<String, String> linkData, String linker,
                             String[] columnHeaders, String linkDataDelimiter, Logger log) {
    boolean replace = true;
    boolean writerWasOpened = false;
    BufferedReader reader;
    PrintWriter writer;
    String sampleDatafilename = proj.SAMPLE_DATA_FILENAME.getValue();
    // in memory backup
    String[][] sampleDataMatrix =
        HashVec.loadFileToStringMatrix(sampleDatafilename, false, null, false);
    String[] sampleDataHeader = Files.getHeaderOfFile(sampleDatafilename, log);
    int linkerIndex = ext.indexOfStr(linker, sampleDataHeader);
    // fail if can't find linker
    if (linkerIndex == -1) {
      log.reportError("Error - could not find linker " + linker + " in sample data file "
                      + sampleDatafilename + ", will not replace data");
      replace = false;
    } else {
      String bakFile = backupSampleData(proj, log);
      if (bakFile == null) {
        replace = false;
      } else {
        log.report("Info - backed up sample data file (" + sampleDatafilename + ") to " + bakFile);
        int numUnchanged = 0;
        int numReplaced = 0;

        try {
          reader = new BufferedReader(new FileReader(bakFile));
          String[] tmpHeader = reader.readLine().trim().split("\t");
          // fail if header of original and backup do not match
          if (!Array.equals(tmpHeader, sampleDataHeader, true)) {
            log.reportError("Error - backup sample data " + bakFile
                            + " does not contain the same header as " + sampleDatafilename
                            + ", will not replace data");
            replace = false;
          } else {
            writerWasOpened = true;
            writer = Files.getAppropriateWriter(sampleDatafilename);
            String dateTime = ext.getDate() + "_" + ext.getTime();
            int[] existingIndices = new int[columnHeaders.length];
            for (int i = 0; i < columnHeaders.length; i++) {
              int[] indices = ext.indicesOfStr(columnHeaders[i], sampleDataHeader, false, true);
              if (indices.length > 0) {
                existingIndices[i] = indices[0];
                for (int index = 1; index < indices.length; index++) {
                  // If more than one instance of the header exists, deactivate the others
                  sampleDataHeader[indices[index]] =
                      "x" + sampleDataHeader[indices[index]] + "_replaced_" + dateTime;
                  log.report("Warning - Additional instance of column \"" + columnHeaders[i]
                             + "\" found, header changed to: " + sampleDataHeader[indices[index]]);
                }
              } else {
                replace = false;
                log.reportError("Error - column \"" + columnHeaders[i] + "\" does not exist in "
                                + sampleDatafilename + ", will not replace data");
              }
            }
            if (replace) {
              writer.println(Array.toStr(sampleDataHeader));
              String[] line;
              HashSet<String> linkers = new HashSet<String>();
              while (reader.ready() && replace) {
                line = reader.readLine().trim().split("\t");
                if (!linkers.add(line[linkerIndex])) {
                  log.reportError("Error - Specified linker \"" + linker
                                  + "\" is not unique, the value \"" + line[linkerIndex]
                                  + "\" appears more than once");
                  log.reportError("Cancelling the replacement and replacing sample data with backup");
                  replace = false;
                } else if (linkData.containsKey(line[linkerIndex])) {
                  String[] data;
                  if (columnHeaders.length > 1) {
                    data = linkData.get(line[linkerIndex]).split(linkDataDelimiter);
                  } else {
                    data = new String[] {linkData.get(line[linkerIndex])};
                  }

                  if (data.length != columnHeaders.length) {
                    // fail if replaced data for a linker is the wrong length after delimiting
                    log.reportError("Error - the number of entries in the column header "
                                    + Array.toStr(columnHeaders) + " (" + columnHeaders.length
                                    + ") does not equal the number of entries in the data "
                                    + Array.toStr(data) + " (" + data.length + ")");
                    log.reportError("Cancelling the replacement and replacing sample data with backup");
                    replace = false;
                  } else {
                    for (int i = 0; i < data.length; i++) {
                      line[existingIndices[i]] = data[i];
                    }
                    writer.println(Array.toStr(line));
                    numReplaced++;
                  }
                } else {
                  writer.println(Array.toStr(line));
                  String[] existingVals = new String[existingIndices.length];
                  for (int i = 0; i < existingIndices.length; i++) {
                    existingVals[i] = line[existingIndices[i]];
                  }
                  linkData.put(line[linkerIndex],
                               Array.toStr(existingVals,
                                           linkDataDelimiter != null ? linkDataDelimiter : "\t"));
                  numUnchanged++;
                }
              }
              writer.close();
            }
          }
          reader.close();
          if (replace) {
            log.report("Info - replaced data in sample data for " + numReplaced
                       + (numReplaced > 1 ? " samples " : " samples"));
            if (numUnchanged > 0) {
              log.report("Warning - " + numUnchanged + " "
                         + (numUnchanged > 1 ? "samples" : "samples") + " were left unchanged ");
            }
          } else if (!replace && writerWasOpened) {
            if (!Files.copyFile(bakFile, sampleDatafilename)) {
              Files.writeMatrix(sampleDataMatrix, sampleDatafilename, "\t");
            }
          }
        } catch (FileNotFoundException e) {
          log.reportError("Error: Sample Data backup file " + bakFile + " not found");
          log.reportException(e);
          replace = false;
        } catch (IOException e) {
          log.reportError("Error: unable to read sample data backup file" + bakFile);
          log.reportException(e);
          replace = false;
        }
      }
    }
    return replace;
  }

  // public void initLinkKey(String filename) {
  // int[] linkKeyColumnLabels = determineKeyIndices(filename);
  // if (linkKeyColumnLabels[DNA_INDEX_IN_LINKERS] >= 0) {
  // // {"DNA/Sample", "DNA", "DNA#", "Sample", "LabID"} exists
  // linkKeyIndex.put(filename, DNA_INDEX_IN_LINKERS);
  // // JOptionPane.showMessageDialog(null, "Link is set to: " +
  // Arrays.toString(LINKERS[DNA_INDEX_IN_LINKERS]), "Information",
  // JOptionPane.INFORMATION_MESSAGE);
  // System.out.println("Link key set to: " + Arrays.toString(LINKERS[DNA_INDEX_IN_LINKERS]));
  // } else if (linkKeyColumnLabels[FID_INDEX_IN_LINKERS] >= 0) {
  // linkKeyIndex.put(filename, FID_INDEX_IN_LINKERS);
  // // JOptionPane.showMessageDialog(null, "Link is set to: " +
  // Arrays.toString(LINKERS[FID_INDEX_IN_LINKERS]), "Information",
  // JOptionPane.INFORMATION_MESSAGE);
  // System.out.println("Link key set to: " + Arrays.toString(LINKERS[FID_INDEX_IN_LINKERS]));
  // } else if (linkKeyColumnLabels[IID_INDEX_IN_LINKERS] >= 0) {
  // linkKeyIndex.put(filename, IID_INDEX_IN_LINKERS);
  // // JOptionPane.showMessageDialog(null, "Link is set to: " +
  // Arrays.toString(LINKERS[IID_INDEX_IN_LINKERS]), "Information",
  // JOptionPane.INFORMATION_MESSAGE);
  // System.out.println("Link key set to: " + Arrays.toString(LINKERS[IID_INDEX_IN_LINKERS]));
  // } else {
  // // JOptionPane.showMessageDialog(null, "Unable to initialize the link key. Please select a link
  // key manually.", "Error", JOptionPane.ERROR_MESSAGE);
  // System.out.println("Unable to initialize the link key.");
  // }
  // }

  // public void setLinkKey(String filename, int selectedLinkKey) {
  // int[] linkKeyColumnLabels;
  //
  // linkKeyColumnLabels = determineKeyIndices(filename);
  //
  // for (int i = 0; i < linkKeyColumnLabels.length; i++) {
  // if ((linkKeyColumnLabels[i] + 1) == selectedLinkKey) {
  // linkKeyIndex.put(filename, i);
  // System.out.println("Link Key set to: " + Arrays.toString(LINKERS[i]));
  // // createLinkKeyToDataHash(selectedNodes[0][0], linkKeyColumnLabels);
  // // JOptionPane.showMessageDialog(null, "Link is set to: " + Arrays.toString(LINKERS[i]),
  // "Information", JOptionPane.INFORMATION_MESSAGE);
  // System.err.println("Link is set to:" + Arrays.toString(LINKERS[i]));
  // return;
  // }
  // }
  // proj.message("Unable to set link key. Please make sure you are selecting a valid key");
  // }

  public void setCnvClasses(String[] cnvClasses) {
    this.cnvClasses = cnvClasses;
  }
}


