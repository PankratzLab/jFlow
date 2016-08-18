package org.genvisis.cyto;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.FileChooser;
import org.genvisis.cnv.manage.SourceFileParser;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.Segment;

public class CytoPanel extends JPanel implements ActionListener {
  private static final long serialVersionUID = 1L;
  public static final String CNV_DIR = "CYTO_VARIANTS/";

  private static final String[] SAMP_DATA_TO_ADD = {"DNA", "FID", "IID", "CLASS=SEX"};
  private static final String SELECT_TITLE =
      "Please select the directory containing workbench and raw data files";
  private static final String WELCOME_TITLE = "Welcome to the Genvisis parser for Agilent files";
  private static final String WORKBENCH_TITLE = "Workbench Aberration File(s)";
  private static final String SAMP_DATA_TITLE = "Raw Data File(s)";
  private static final String SELECT_NEW = "Select a new";
  private static final String WORKBENCH = " workbench file";
  private static final String REMOVE = "REMOVE File";
  private static final String SAMPLES = " feature extraction file";
  private static final String NO_FILE_SELECTED = " No file selected";
  private static final String WORKBENCH_EXT = "workbench.xls";
  private static final String SAMPLES_EXT = ".txt";
  private static final String PANEL_LOG = "parser.log";
  private static final String FINAL_CNV = "ALL_Variants" + CytoCNVariant.CNV_EXT;

  private static final String SER = ".ser";
  private static final String ALL = ".all";
  private static final String IND = ".ind";

  private static final int DEFAULT_MAX_NUM_WORKBENCH = 4;
  private static final int DEFAULT_MAX_NUM_SAMPLES = 16;

  private final CytoGUI cytoGUI;
  private final Project proj;
  private String startDir;
  private Thread thread;
  private String[] workBenchFiles, importFiles, sampleNames, sampleIndices;
  private JLabel jLabelMain, jLabelBench, jLabelSamp;
  private JButton selectFilesButton, parseButton, cancelButton;
  private FileBox[] benchfileBoxes, sampfileBoxes;
  private final Logger log;

  public CytoPanel(Project proj, CytoGUI cytoGUI, String startDir, Logger log) {
    this.proj = proj;
    this.startDir = startDir;
    this.log = (log == null ? new Logger(proj.PROJECT_DIRECTORY.getValue()
                                         + ext.getTimestampForFilename() + PANEL_LOG)
                            : log);
    this.cytoGUI = cytoGUI;
    setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
    initLabels();
    add(jLabelMain);
    initButtons();
    setVisible(true);
  }

  private void initLabels() {
    jLabelMain = getLabel(WELCOME_TITLE, true);
    jLabelBench = getLabel(WORKBENCH_TITLE, false);
    jLabelSamp = getLabel(SAMP_DATA_TITLE, false);
  }

  private JLabel getLabel(String title, boolean visible) {
    JLabel jLabel = new JLabel();
    jLabel.setText(title);
    jLabel.setVisible(visible);
    return jLabel;
  }

  /**
   * add select file button, init parse and cancel buttonsb but do not show
   */
  public void initButtons() {
    selectFilesButton =
        getButton("Select Directory", "Select the directory containing data to parse", true);
    parseButton = getButton("PARSE", "Import data to Genvisis", true);
    cancelButton = getButton("Cancel", "", true);
    add(selectFilesButton);
    // add(importButton);
  }

  private JButton getButton(String title, String tip, boolean enabled) {
    JButton jbutton = new JButton(title);
    jbutton.setToolTipText(tip);
    jbutton.addActionListener(this);
    jbutton.setEnabled(enabled);
    return jbutton;
  }

  /**
   * Reset to the welcome type screen
   */
  private void reset() {
    jLabelMain.setVisible(true);
    jLabelBench.setVisible(false);
    jLabelSamp.setVisible(false);
    parseButton.setVisible(false);
    cancelButton.setVisible(false);
    removeAllBoxes(benchfileBoxes);
    removeAllBoxes(sampfileBoxes);
    benchfileBoxes = new FileBox[0];
    sampfileBoxes = new FileBox[0];
    workBenchFiles = new String[0];
    selectFilesButton.setVisible(true);
    cytoGUI.pack();
  }

  @Override
  public void actionPerformed(ActionEvent e) {

    JComponent source = (JComponent) e.getSource();
    // gather the files from the directory navigated to
    if (source.equals(selectFilesButton)) {
      FileChooser fileChooser =
          new FileChooser(this, proj.PROJECT_DIRECTORY.getValue(), false, true, SELECT_TITLE, log);
      if (fileChooser.isSelected()) {
        startDir = fileChooser.getNavDir();
        log.report(startDir);
        workBenchFiles = Files.list(startDir, WORKBENCH_EXT, false);
        checkFiles(workBenchFiles, WORKBENCH_EXT, true);
        if (!errorFindingFromWorkBench()) {
          // assign fileBoxes and update the panel to reflect found files
          updateFiles();
        } else {
          reset();
        }

      }
    } else if (source.equals(cancelButton)) {
      reset();
    } else if (source.equals(parseButton) && checkThread()) {
      boolean keepGoing = addToSampleData(proj, sampfileBoxes, log);
      if (keepGoing) {
        thread = new Thread(new Runnable() {
          @Override
          public void run() {
            // run the heavy lifting in a separate thread so we can update with messages
            importFiles = getFinalFullPath(sampfileBoxes);

            if (checkAlreadyParsed()) {
              parseAndImport();
            }

            if (keepGoing()) {
              parseWorkbench();
            }
            writeVariants();
            promptComplete();
            reset();
          }
        });
        thread.start();
      }
    }

  }

  // TODO, check add new files, or remove it -> likely can just remove it since we are going to
  // auto-detect everything
  /**
   * If we find too many files for a sample, that is bad. If we do not find one, we let it slide.
   * Currently we only allow .txt extenstions (SAMPLES_EXT)
   */
  private boolean errorFindingFromWorkBench() {
    String[] allFiles = Files.list(startDir, SAMPLES_EXT, false);
    boolean error = false;
    ArrayList<String> tmpFile = new ArrayList<String>();
    ArrayList<String> tmpName = new ArrayList<String>();
    ArrayList<String> tmpIndex = new ArrayList<String>();
    for (int i = 0; i < workBenchFiles.length; i++) {
      String[] headers = CytoCNVariant.getSampleNames(startDir + workBenchFiles[i], log);
      for (int j = 0; j < headers.length; j++) {
        boolean found = false;
        for (String allFile : allFiles) {
          if (!found && match(allFile, headers[j])) {
            tmpFile.add(allFile);
            tmpName.add(ext.replaceWithLinuxSafeCharacters(headers[j], true));
            tmpIndex.add("" + (j + 1));
            found = true;
          } else if (found && match(allFile, headers[j])) {
            multiMatch(headers[j], workBenchFiles[i]);
            return true;
          }
        }
        if (!found) {
          if (!NoMatchfor(parseToMatch(headers[j]), workBenchFiles[i])) {
            return false;
          }
        }
      }
    }
    importFiles = tmpFile.toArray(new String[tmpFile.size()]);
    sampleNames = tmpName.toArray(new String[tmpName.size()]);
    sampleIndices = tmpIndex.toArray(new String[tmpIndex.size()]);
    return error;
  }

  private boolean match(String file, String header) {
    return parseToMatch(file).equals(parseToMatch(header));
  }

  /**
   * Parses a string to the match format serving as a link between workbench files and feature
   * extraction
   */
  private static String parseToMatch(String original) {
    String forMatch = original;
    if (forMatch.endsWith(".xls")) {
      forMatch = "";
    } else {
      forMatch = ext.rootOf(forMatch);
      forMatch = ext.replaceWithLinuxSafeCharacters(forMatch, true);
      forMatch = forMatch.substring(forMatch.length() - 3, forMatch.length());
    }
    return forMatch;
  }

  /**
   * To prevent running the same thing twice (and concurrently)
   * 
   * @return true if the thread is available
   */
  private boolean checkThread() {
    if (thread != null && thread.isAlive()) {
      String[] options = new String[] {"OK"};
      String dialoge =
          "Warning - processing files in the background, please wait until completion\n";
      JOptionPane.showOptionDialog(null, dialoge, "Error", JOptionPane.DEFAULT_OPTION,
                                   JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
      return false;
    }
    return true;
  }

  /**
   * When everything finished
   */
  private void promptComplete() {
    String[] options = new String[] {"OK"};
    String dialoge = "Parsing complete\n";
    JOptionPane.showOptionDialog(null, dialoge, "Done...", JOptionPane.DEFAULT_OPTION,
                                 JOptionPane.QUESTION_MESSAGE, null, options, options[0]);

  }

  /**
   * If there are multiple files matching a particular sample from the workbench
   */
  private void multiMatch(String sampleExt, String workbenchFile) {
    String[] options = new String[] {"OK"};
    String dialoge = "Error: multiple files ended with the extension " + parseToMatch(sampleExt)
                     + SAMPLES_EXT + " corresponding to sample found in " + workbenchFile + " \n";
    dialoge += "Please make sure there is only one file ending with " + parseToMatch(sampleExt)
               + SAMPLES_EXT;

    JOptionPane.showOptionDialog(null, dialoge, "Error", JOptionPane.DEFAULT_OPTION,
                                 JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
  }

  /**
   * If there are no files matching a particular sample from the workbench
   */
  private boolean NoMatchfor(String sampleExt, String workbenchFile) {
    String[] options = new String[] {"Cancel", "OK"};
    String dialoge = "Warning: Could not find a file ending with extension " + sampleExt
                     + SAMPLES_EXT + " corresponding to sample found in " + workbenchFile + " \n";
    int response =
        JOptionPane.showOptionDialog(null, dialoge, "Error", JOptionPane.DEFAULT_OPTION,
                                     JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
    return (response == 0 ? false : true);
  }

  /**
   * If files to import are present, imports the files to Genvisis. We also delete all samRAF and
   * mkRAF files
   */
  private void parseAndImport() {
    if (checkFiles(importFiles, SAMPLES, true)) {
      Files.write("", proj.SAMPLE_DIRECTORY.getValue(true, true)
                      + SourceFileParser.OVERWRITE_OPTION_FILE);
      SourceFileParser.deleteAllFilesInSampleDirectory(proj);
      SourceFileParser.deleteAllFilesInMarkerDataDirectory(proj);

      log.report(ext.getTime() + " Info - found " + importFiles.length + " files");
      CytoAgilentParse.parseCytoToGenvisis(proj, importFiles, log);
      // cnv.manage.ParseIllumina.createFiles(proj, proj.getInt(proj.NUM_THREADS));
      org.genvisis.cnv.manage.SourceFileParser.createFiles(proj, proj.NUM_THREADS.getValue());
      deleteSampleList();
      TransposeData.transposeData(proj, 2000000000, false);

      new File(proj.SAMPLE_DIRECTORY.getValue(true, true)
               + SourceFileParser.OVERWRITE_OPTION_FILE).delete();
    }
  }

  /**
   * Parse the *workbench.xls files, if something failed during import to Genvisis, this should
   * still work, but Beast scores will just be NA
   */
  private void parseWorkbench() {
    try {
      workBenchFiles = getFinalFullPath(benchfileBoxes);
      if (checkFiles(workBenchFiles, WORKBENCH_TITLE, false)) {
        // CytoCompare.compare(proj, workBenchFiles, proj.getFilename(proj.COMMON_CNP_FILENAME),
        // proj.getFilename(proj.REPORTED_CNP_FILENAME),
        // proj.getFilename(proj.UNREPORTED_CNP_FILENAME), CytoCompare.DEFAULT_SCORE_THRESHOLD,
        // startDir, true, new Logger(startDir + ext.getTimestampForFilename() + "_" +
        // CytoCompare.DEFAULT_LOG));
        CytoCompare.compare(proj, workBenchFiles, proj.COMMON_CNP_FILENAME.getValue(),
                            proj.REPORTED_CNP_FILENAME.getValue(),
                            proj.UNREPORTED_CNP_FILENAME.getValue(),
                            CytoCompare.DEFAULT_SCORE_THRESHOLD, startDir, true,
                            new Logger(startDir + ext.getTimestampForFilename() + "_"
                                       + CytoCompare.DEFAULT_LOG));
      }
    } catch (Exception exc) {
      JOptionPane.showMessageDialog(this,
                                    "An error occured while parsing the following "
                                          + (workBenchFiles.length == 1 ? " file: " : " files: ")
                                          + Array.toStr(workBenchFiles, "\n"));
      log.reportError("An error occured while parsing the following "
                      + (workBenchFiles.length == 1 ? " file: " : " files: ")
                      + Array.toStr(workBenchFiles, "\n"));
      log.reportException(exc);
    }
  }

  /**
   * We are assuming that the workbench files are the final version (valid, and desired to be used).
   * getFinalFullPath(benchfileBoxes) should be called prior Cats
   *
   */
  private void writeVariants() {
    String dir = proj.PROJECT_DIRECTORY.getValue() + CNV_DIR;
    if (!Files.exists(dir)) {
      new File(dir).mkdir();
    }
    if (Files.exists(dir + FINAL_CNV)) {
      new File(dir + FINAL_CNV).delete();
    }

    if (Files.exists(dir + FINAL_CNV + ALL + SER)) {
      new File(dir + FINAL_CNV + ALL + SER).delete();
    }
    if (Files.exists(dir + FINAL_CNV + IND + SER)) {
      new File(dir + FINAL_CNV + IND + SER).delete();
    }
    for (String workBenchFile : workBenchFiles) {
      CytoCNVariant.writeIndCNVariantFiles(CytoCNVariant.directToInds(workBenchFile, log), dir,
                                           log);
    }
    String[] cnvsToCat = Files.toFullPaths(Files.list(dir, CytoCNVariant.CNV_EXT, false), dir);
    cnvsToCat = filterCNP(proj, cnvsToCat);
    Files.cat(cnvsToCat, dir + FINAL_CNV, new int[0], log);
    if (!projectHasCNVFile(proj, CNV_DIR + FINAL_CNV, log)) {
      addCNVFile(proj, CNV_DIR + FINAL_CNV);
    }
    convertCNPfiles(proj, CNV_DIR, log);
    // TODO document and regions list etc
  }


  /**
   * Filter out the CNP files (cyto data base list of variants) to avoid writing them to all
   * variants
   */
  private static String[] filterCNP(Project proj, String[] files) {
    ArrayList<String> filter = new ArrayList<String>();
    for (int i = 0; i < files.length; i++) {
      // if (!linuxSafeRootEquals(proj.getFilename(proj.COMMON_CNP_FILENAME),files[i]) &&
      // !linuxSafeRootEquals(proj.getFilename(proj.REPORTED_CNP_FILENAME),files[i]) &&
      // !linuxSafeRootEquals(proj.getFilename(proj.UNREPORTED_CNP_FILENAME),files[i])) {
      if (!linuxSafeRootEquals(proj.COMMON_CNP_FILENAME.getValue(), files[i])
          && !linuxSafeRootEquals(proj.REPORTED_CNP_FILENAME.getValue(), files[i])
          && !linuxSafeRootEquals(proj.UNREPORTED_CNP_FILENAME.getValue(), files[i])) {
        filter.add(files[i]);
      }

    }
    return filter.toArray(new String[filter.size()]);
  }

  /**
   * Check to see if the file matches after replacing by linux and taking root
   */
  private static boolean linuxSafeRootEquals(String toLinux, String toEqual) {
    if (!ext.replaceWithLinuxSafeCharacters(ext.rootOf(toLinux), true).equals(ext.rootOf(toEqual))
        && !ext.rootOf(toLinux).equals(ext.rootOf(toEqual))) {
      return false;
    } else {
      return true;
    }
  }

  private static void convertCNPfiles(Project proj, String cnvdir, Logger log) {
    // convertCNPfiles(proj, proj.getFilename(proj.COMMON_CNP_FILENAME), cnvdir +
    // ext.replaceWithLinuxSafeCharacters(ext.rootOf(proj.getFilename(proj.COMMON_CNP_FILENAME)),
    // true) + CytoCNVariant.CNV_EXT, log);
    // convertCNPfiles(proj, proj.getFilename(proj.REPORTED_CNP_FILENAME), cnvdir +
    // ext.replaceWithLinuxSafeCharacters(ext.rootOf(proj.getFilename(proj.REPORTED_CNP_FILENAME)),
    // true) + CytoCNVariant.CNV_EXT, log);
    // convertCNPfiles(proj, proj.getFilename(proj.UNREPORTED_CNP_FILENAME), cnvdir +
    // ext.replaceWithLinuxSafeCharacters(ext.rootOf(proj.getFilename(proj.UNREPORTED_CNP_FILENAME)),
    // true) + CytoCNVariant.CNV_EXT, log);
    convertCNPfiles(proj,
                    proj.COMMON_CNP_FILENAME.getValue(), cnvdir
                                                         + ext.replaceWithLinuxSafeCharacters(ext.rootOf(proj.COMMON_CNP_FILENAME.getValue()),
                                                                                              true)
                                                         + CytoCNVariant.CNV_EXT,
                    log);
    convertCNPfiles(proj,
                    proj.REPORTED_CNP_FILENAME.getValue(), cnvdir
                                                           + ext.replaceWithLinuxSafeCharacters(ext.rootOf(proj.REPORTED_CNP_FILENAME.getValue()),
                                                                                                true)
                                                           + CytoCNVariant.CNV_EXT,
                    log);
    convertCNPfiles(proj,
                    proj.UNREPORTED_CNP_FILENAME.getValue(), cnvdir
                                                             + ext.replaceWithLinuxSafeCharacters(ext.rootOf(proj.UNREPORTED_CNP_FILENAME.getValue()),
                                                                                                  true)
                                                             + CytoCNVariant.CNV_EXT,
                    log);
  }

  /**
   * We do not convert if the cnp file is already present, we do not add to properties if it is
   * already present
   */
  private static void convertCNPfiles(Project proj, String cnpFile, String output, Logger log) {
    if (!projectHasCNVFile(proj, output, log)) {
      addCNVFile(proj, output);
    }
    if (!Files.exists(proj.PROJECT_DIRECTORY.getValue() + output)) {
      convertCNPFile(cnpFile, proj.PROJECT_DIRECTORY.getValue() + output, log);
    }
  }

  /**
   * add a cnv file to the projects cnv list
   */
  private static void addCNVFile(Project proj, String cnvFile) {
    // proj.setProperty(proj.CNV_FILENAMES, proj.getProperty(proj.CNV_FILENAMES) +
    // (proj.getProperty(proj.CNV_FILENAMES).length() > 0 ? ";" : "") + cnvFile);
    String[] vals = proj.getProperty(proj.CNV_FILENAMES);
    proj.setProperty(proj.CNV_FILENAMES, Array.addStrToArray(cnvFile, vals));
    proj.saveProperties();
  }

  /**
   * If the cnvFile is listed in the projects cnv list, return true
   */
  private static boolean projectHasCNVFile(Project proj, String cnvFile, Logger log) {
    boolean has = false;
    // String[] current = proj.getProperty(proj.CNV_FILENAMES).split(";");
    String[] current = proj.getProperty(proj.CNV_FILENAMES);
    if (current != null && current.length > 0
        && ext.indexOfStr(cnvFile, current, true, true, log, false) >= 0) {
      has = true;
    }
    return has;
  }

  /**
   * Generates a mock file of CNVariants from the cnpFile, FID = filename, IID = location
   */
  private static boolean convertCNPFile(String cnpFile, String output, Logger log) {
    boolean converted = false;
    if (!Files.exists(output)) {
      org.genvisis.filesys.Segment[] segs = CytoCompare.loadsegs(cnpFile, log);
      if (segs != null) {
        try {
          PrintWriter writer = new PrintWriter(new FileWriter(output));
          writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
          for (Segment seg : segs) {
            writer.println(new CNVariant(ext.rootOf(output), seg.getUCSClocation(), seg.getChr(),
                                         seg.getStart(), seg.getStop(), 2, 10f, 10,
                                         0).toPlinkFormat());
          }
          writer.close();
          converted = true;
        } catch (Exception e) {
          log.reportError("Error writing to " + output);
          log.reportException(e);
        }
      }
    } else {
      converted = true;
    }
    return converted;
  }

  /**
   * We delete the sample list so that it is created only from available files in the transpose step
   * (after sample import), otherwise things can get goofy
   */
  private void deleteSampleList() {

    // if (Files.exists(proj.getFilename(proj.SAMPLELIST_FILENAME)) && proj.getSampleList() != null
    // && proj.getSampleList().getSamples().length > 0) {
    // new File(proj.getFilename(proj.SAMPLELIST_FILENAME)).delete();
    if (Files.exists(proj.SAMPLELIST_FILENAME.getValue()) && proj.getSampleList() != null
        && proj.getSampleList().getSamples().length > 0) {
      new File(proj.SAMPLELIST_FILENAME.getValue()).delete();
    }
  }

  /**
   * @return default true, false if files have already been parsed and want to cancel
   */
  private boolean checkAlreadyParsed() {
    ArrayList<String> alreadyParsed = new ArrayList<String>();
    if (importFiles != null && importFiles.length > 0) {
      for (String importFile : importFiles) {
        String file = proj.SOURCE_DIRECTORY.getValue(false, true) + ext.rootOf(importFile, true)
                      + proj.getProperty(proj.SOURCE_FILENAME_EXTENSION);
        if (Files.exists(file)) {
          alreadyParsed.add(file);
        }
      }
    }
    if (alreadyParsed.size() > 0) {
      String[] overwriteOptions = new String[] {"Overwrite All", "Cancel parser"};
      String samples =
          alreadyParsed.size() + " " + (alreadyParsed.size() > 1 ? " of the selected samples "
                                                                 : " the selected sample");
      int response =
          JOptionPane.showOptionDialog(null,
                                       "These data (at least  " + samples
                                             + ") have already been parsed to the project source directory "
                                             + proj.SOURCE_DIRECTORY.getValue(false, true) + ".\n"
                                             + "This happens if you inadvertently restarted the parser or if the parser was interrupted and manually restarted, or if there are non-unique filenames.\n"
                                             + "Select \"" + overwriteOptions[0]
                                             + "\" earlier files to overwrite the " + samples
                                             + " already processed.\n"
                                             + "Otherwise, cancel and manually remove the files.\n"
                                             + "What would you like to do?",
                                       "Samples already exist", JOptionPane.DEFAULT_OPTION,
                                       JOptionPane.QUESTION_MESSAGE, null, overwriteOptions,
                                       overwriteOptions[1]);

      switch (response) {
        case 0:
          for (int i = 0; i < alreadyParsed.size(); i++) {
            new File(alreadyParsed.get(i)).delete();
          }
          return true;
        case 1:
          return false;
        default:
          return false;
      }
    }
    return true;
  }

  /**
   * Here we create FileBox objects corresponding to the files discovered when listing the
   * directory, feature extraction files get special treatment with index and sample name in the
   * workbench file
   */
  private void updateFiles() {
    int numWorkBoxes = Math.min(workBenchFiles.length, DEFAULT_MAX_NUM_WORKBENCH);
    int numSampBoxes = Math.min(importFiles.length, DEFAULT_MAX_NUM_SAMPLES);
    checkTooMany(workBenchFiles, DEFAULT_MAX_NUM_WORKBENCH, WORKBENCH_EXT);
    checkTooMany(importFiles, DEFAULT_MAX_NUM_SAMPLES, SAMPLES_EXT);

    jLabelMain.setVisible(false);
    selectFilesButton.setVisible(false);
    jLabelBench.setVisible(true);
    add(jLabelBench);
    benchfileBoxes = new FileBox[numWorkBoxes];
    for (int i = 0; i < benchfileBoxes.length; i++) {
      log.report("" + i);
      if (i < workBenchFiles.length) {
        benchfileBoxes[i] = new FileBox(cytoGUI, startDir, workBenchFiles[i],
                                        SELECT_NEW + WORKBENCH, REMOVE, "-1", null);
      } else {
        benchfileBoxes[i] = new FileBox(cytoGUI, NO_FILE_SELECTED, null, SELECT_NEW + WORKBENCH,
                                        REMOVE, "-1", null);
      }
      if (i < DEFAULT_MAX_NUM_WORKBENCH) {
        benchfileBoxes[i].setVisible(true);
        add(benchfileBoxes[i]);
      }
    }
    jLabelSamp.setVisible(true);
    add(jLabelSamp);

    sampfileBoxes = new FileBox[numSampBoxes];
    for (int i = 0; i < sampfileBoxes.length; i++) {
      if (i < importFiles.length) {
        sampfileBoxes[i] = new FileBox(cytoGUI, startDir, importFiles[i], SELECT_NEW + SAMPLES,
                                       REMOVE, sampleIndices[i], sampleNames[i]);
      } else {
        sampfileBoxes[i] =
            new FileBox(cytoGUI, NO_FILE_SELECTED, null, SELECT_NEW + SAMPLES, REMOVE, "-1", null);
      }
      // we only add up to a max amount to the Gui
      if (i < DEFAULT_MAX_NUM_SAMPLES) {
        sampfileBoxes[i].setVisible(true);
        add(sampfileBoxes[i]);
      }

    }
    // ScrollFileBox workerScrollBox = new ScrollFileBox(benchfileBoxes);
    cancelButton.setVisible(true);
    parseButton.setVisible(true);

    add(cancelButton, BorderLayout.EAST);
    add(parseButton, BorderLayout.WEST);

    cytoGUI.setSize(800, 600);
    cytoGUI.pack();
  }

  /**
   * @param files the files found
   * @param tooMany how many we will display in the gui
   * @param type the type of file found
   */
  public void checkTooMany(String[] files, int tooMany, String type) {
    if (files.length > tooMany) {
      String[] options = new String[] {"OK", "See all files found"};
      String dialoge =
          "Warning - detected " + files.length + " files of type " + type + ", only " + tooMany
                       + " will be displayed. Please make sure that the number detected is correct.";
      int response =
          JOptionPane.showOptionDialog(null, dialoge, "Warning", JOptionPane.DEFAULT_OPTION,
                                       JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
      if (response == 1) {
        JOptionPane.showOptionDialog(null, Array.toStr(files, "\n"), "Warning",
                                     JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null,
                                     new String[] {options[0]}, options[0]);
      }
    }
  }

  /**
   * Alert the user if ids are already present, if that is ok by them we skip the addition We check
   * both DNA, and FID/IID
   */
  private static boolean verifySampleDataIssues(String[] issues, String sampleDataFile) {
    boolean keepGoing = true;
    if (issues != null && issues.length > 0) {
      String[] options = new String[] {"Continue", "Cancel"};
      String dialoge =
          "Warning - the following samples were already present in the sample data file \""
                       + sampleDataFile + "\".\n";
      dialoge +=
          "If you are re-parsing the same samples, select \"Continue\", otherwise please ensure that samples have unique names\n";
      dialoge += "Please note that a FID/IID combination must be unique\n";
      dialoge += Array.toStr(issues, "\n");
      int response =
          JOptionPane.showOptionDialog(null, dialoge, "Error", JOptionPane.DEFAULT_OPTION,
                                       JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
      if (response != 0) {
        keepGoing = false;
      }
    }
    return keepGoing;

  }

  /**
   * Check if the individuals are already present
   */
  private static boolean checkSampleData(Project proj, FileBox[] fileBoxes, Logger log) {
    ArrayList<String> issues = new ArrayList<String>();
    SampleData sampleData = proj.getSampleData(0, false);
    for (FileBox fileBoxe : fileBoxes) {
      try {
        String[] id = sampleData.lookup(fileBoxe.getDNA());
        id = sampleData.lookup(fileBoxe.getFID() + "\t" + fileBoxe.getIID());
        log.report("Previous Individual: " + id[0] + "\t" + id[1]);
        issues.add("   DNA:  " + fileBoxe.getDNA() + "\t   FID:  " + fileBoxe.getFID()
                   + "\t  IID:  " + fileBoxe.getIID());
        fileBoxe.setAddToSampleData(false);
      } catch (NullPointerException npe) {
        log.report("New Individual: " + fileBoxe.getDNA() + "\t" + fileBoxe.getFID() + "\t"
                   + fileBoxe.getIID());
      }
    }
    proj.resetSampleData();
    // return verifySampleDataIssues(issues.toArray(new String[issues.size()]),
    // proj.getFilename(proj.SAMPLE_DATA_FILENAME));
    return verifySampleDataIssues(issues.toArray(new String[issues.size()]),
                                  proj.SAMPLE_DATA_FILENAME.getValue());
  }

  /**
   * Two branches, if sample data does not exist, we create it. Else we check and append
   */
  private static boolean addToSampleData(Project proj, FileBox[] fileBoxes, Logger log) {
    boolean keepGoing = false;
    fileBoxes = getBoxesToUse(fileBoxes);
    // String sampleDataFile = proj.getFilename(proj.SAMPLE_DATA_FILENAME);
    String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue();
    if (!Files.exists(sampleDataFile)) {
      keepGoing = writeNewSampleData(fileBoxes, sampleDataFile);
    } else {

      keepGoing = checkSampleData(proj, fileBoxes, log);
      if (keepGoing) {
        String[] header = Files.getHeaderOfFile(sampleDataFile, log);
        int DNAIndex = ext.indexOfStr(SAMP_DATA_TO_ADD[0], header);
        keepGoing = warnIndex(DNAIndex, SAMP_DATA_TO_ADD[0], sampleDataFile);
        int FIDIndex = ext.indexOfStr(SAMP_DATA_TO_ADD[1], header);
        keepGoing =
            (keepGoing ? warnIndex(DNAIndex, SAMP_DATA_TO_ADD[1], sampleDataFile) : keepGoing);
        int IIDIndex = ext.indexOfStr(SAMP_DATA_TO_ADD[2], header);
        keepGoing =
            (keepGoing ? warnIndex(DNAIndex, SAMP_DATA_TO_ADD[2], sampleDataFile) : keepGoing);
        int SEXINDEX = ext.indexOfStr(SAMP_DATA_TO_ADD[3], header);
        keepGoing =
            (keepGoing ? warnIndex(DNAIndex, SAMP_DATA_TO_ADD[3], sampleDataFile) : keepGoing);
        if (keepGoing) {
          try {
            PrintWriter writer = new PrintWriter(new FileWriter(sampleDataFile, true));
            for (FileBox fileBoxe : fileBoxes) {
              if (fileBoxe.isAddToSampleData()) {
                String print = "";
                for (int j = 0; j < header.length; j++) {
                  if (j == DNAIndex) {
                    print += fileBoxe.getDNA() + (j < (header.length - 1) ? "\t" : "");
                  } else if (j == FIDIndex) {
                    print += fileBoxe.getFID() + (j < (header.length - 1) ? "\t" : "");
                  } else if (j == IIDIndex) {
                    print += fileBoxe.getIID() + (j < (header.length - 1) ? "\t" : "");
                  } else if (j == SEXINDEX) {
                    print += "0" + (j < (header.length - 1) ? "\t" : "");
                  } else {
                    print += "." + (j < (header.length - 1) ? "\t" : "");
                  }
                }
                writer.println(print);
              }
            }
            writer.close();

          } catch (IOException e) {
            log.reportError("Error - could not write to sample data file " + sampleDataFile);
            keepGoing = false;
            e.printStackTrace();
          }
        }
      }
    }
    return keepGoing;
  }

  private static boolean writeNewSampleData(FileBox[] fileBoxes, String sampleDataFile) {
    boolean keepGoing;
    PrintWriter writer = Files.getAppropriateWriter(sampleDataFile);
    writer.println(Array.toStr(SAMP_DATA_TO_ADD));
    for (FileBox fileBoxe : fileBoxes) {
      writer.println(fileBoxe.getDNA() + "\t" + fileBoxe.getFID() + "\t" + fileBoxe.getIID()
                     + "\t0");
    }
    writer.close();
    keepGoing = true;
    return keepGoing;
  }

  /**
   *
   * We require certain elements in the sample data, if the are not present (ext returned -1) we
   * warn and cancel
   */
  public static boolean warnIndex(int index, String type, String sampleDataFile) {
    boolean keepGoing = true;
    if (index < 0) {
      String[] options = new String[] {"OK",};
      String dialoge =
          "Error: could not find column " + type + " in sample data file " + sampleDataFile;
      JOptionPane.showOptionDialog(null, dialoge, "Error", JOptionPane.DEFAULT_OPTION,
                                   JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
      keepGoing = false;
    }
    return keepGoing;
  }

  /**
   * @param fileBoxes set all of these to invisible
   */
  public void removeAllBoxes(FileBox[] fileBoxes) {
    if (fileBoxes != null && fileBoxes.length > 0) {
      for (FileBox fileBoxe : fileBoxes) {
        fileBoxe.setVisible(false);
      }
    }
  }

  /**
   * We retrieve the text from the FileBox's textfields, and pass it on if it was not removed by the
   * user
   */
  private static String[] getFinalFullPath(FileBox[] fileBoxes) {
    String[] finalFiles;
    fileBoxes = getBoxesToUse(fileBoxes);
    finalFiles = new String[fileBoxes.length];
    for (int i = 0; i < fileBoxes.length; i++) {
      if (fileBoxes[i].isUsed()) {
        finalFiles[i] = fileBoxes[i].getfullPath();
      }
    }
    return finalFiles;
  }

  private static FileBox[] getBoxesToUse(FileBox[] fileBoxes) {
    ArrayList<FileBox> boxesToUse = new ArrayList<FileBox>();
    for (FileBox fileBoxe : fileBoxes) {
      if (fileBoxe.isUsed()) {
        boxesToUse.add(fileBoxe);
      }
    }
    return boxesToUse.toArray(new FileBox[boxesToUse.size()]);
  }

  /**
   * @param array array of files, check for null and length greater than 1
   * @param category for the files
   * @return
   */
  private static boolean checkFiles(String[] array, String category, boolean message) {
    if (array == null || array.length < 1) {
      if (message) {
        String[] options = new String[] {"OK"};
        String dialoge = "There seems to be no files to parse for type\"" + category
                         + "\", skipping this step\n";
        JOptionPane.showOptionDialog(null, dialoge, "Error", JOptionPane.DEFAULT_OPTION,
                                     JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
      }
      return false;
    }
    return true;
  }

  /**
   * Here we check for the files to compare overlap against the aberrations in the workbench files.
   * We allow them not to exist, but let everyone know first
   */
  private boolean keepGoing() {
    int response = 1;
    String[] options = new String[] {"Cancel", "Continue without"};
    // if (!Files.exists(proj.getFilename(proj.COMMON_CNP_FILENAME))) {
    if (!Files.exists(proj.COMMON_CNP_FILENAME.getValue())) {
      // String dialoge = "Could not find the Common CNP file " +
      // proj.getFilename(proj.COMMON_CNP_FILENAME);
      String dialoge = "Could not find the Common CNP file " + proj.COMMON_CNP_FILENAME.getValue();
      response =
          JOptionPane.showOptionDialog(null, dialoge, "Error", JOptionPane.DEFAULT_OPTION,
                                       JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
    }
    // if (parseResponse(response) && !Files.exists(proj.getFilename(proj.REPORTED_CNP_FILENAME))) {
    // String dialoge = "Could not find the Reported CNP file " +
    // proj.getFilename(proj.REPORTED_CNP_FILENAME);
    if (parseResponse(response) && !Files.exists(proj.REPORTED_CNP_FILENAME.getValue())) {
      String dialoge =
          "Could not find the Reported CNP file " + proj.REPORTED_CNP_FILENAME.getValue();
      response =
          JOptionPane.showOptionDialog(this, dialoge, "Error", JOptionPane.DEFAULT_OPTION,
                                       JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
    }
    // if (parseResponse(response) && !Files.exists(proj.getFilename(proj.COMMON_CNP_FILENAME))) {
    // String dialoge = "Could not find the Un - Reported CNP file " +
    // proj.getFilename(proj.UNREPORTED_CNP_FILENAME);
    if (parseResponse(response) && !Files.exists(proj.COMMON_CNP_FILENAME.getValue())) {
      String dialoge =
          "Could not find the  Un - Reported CNP file " + proj.UNREPORTED_CNP_FILENAME.getValue();
      response =
          JOptionPane.showOptionDialog(this, dialoge, "Error", JOptionPane.DEFAULT_OPTION,
                                       JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
    }
    return parseResponse(response);
  }

  private boolean parseResponse(int response) {
    if (response == 1) {
      return true;
    } else {
      return false;
    }
  }

  /**
   * A helper class to remove/track files discovered
   *
   */
  public class FileBox extends JPanel implements ActionListener {

    private static final long serialVersionUID = 1L;
    private final String file;
    private final String startDir;
    private final JButton selectButton, removeButton;
    private final JTextField jFileNameTextField;
    private JTextField jindexTextField;
    private JTextField jSampleNameTextField;
    private final JFrame jframe;
    // once removed, or if text field is blank, we set this to false
    private boolean used, addToSampleData;

    /**
     * @param startDir starting directory
     * @param file file to put as text
     * @param titleSelect title of select button
     * @param titleRemove title of the remove button
     */

    public FileBox(JFrame jframe, String startDir, String file, String titleSelect,
                   String titleRemove, String workBenchIndex, String sampleName) {
      this.startDir = startDir;
      this.file = file;
      jFileNameTextField = new JTextField();
      selectButton = new JButton(titleSelect);
      removeButton = new JButton(titleRemove);
      this.jframe = jframe;
      used = false;
      addToSampleData = true;
      setText(null);
      jFileNameTextField.addActionListener(this);
      selectButton.addActionListener(this);
      removeButton.addActionListener(this);

      add(selectButton, BorderLayout.EAST);
      add(jFileNameTextField, BorderLayout.WEST);
      setIndexText(workBenchIndex);
      setSampleNametext(sampleName);
      add(removeButton, BorderLayout.WEST);

    }

    public void setIndexText(String index) {
      if (Integer.parseInt(index) >= 0) {
        jindexTextField = new JTextField(index);
        jindexTextField.setToolTipText("This is the sample index in the workbench.xls file");
        jindexTextField.addActionListener(this);
        add(jindexTextField, BorderLayout.EAST);
      }
    }

    public void setSampleNametext(String sampleName) {
      if (sampleName != null) {
        jSampleNameTextField = new JTextField(sampleName);
        jSampleNameTextField.setToolTipText("This is the sample name in the workbench.xls file");
        jSampleNameTextField.addActionListener(this);
        add(jSampleNameTextField, BorderLayout.EAST);
      }
    }


    public boolean isUsed() {
      return used;
    }

    public boolean isAddToSampleData() {
      return addToSampleData;
    }

    public void setAddToSampleData(boolean addToSampleData) {
      this.addToSampleData = addToSampleData;
    }

    public String getText() {
      return jFileNameTextField.getText();
    }

    public String getfullPath() {
      return startDir + jFileNameTextField.getText();
    }

    private void setText(String newText) {
      String text = "";

      if (newText != null) {
        text = newText;
      } else if (file != null) {
        text += file;
        used = true;
      }

      jFileNameTextField.setText(text);
    }

    public String getDNA() {
      return ext.replaceWithLinuxSafeCharacters(ext.rootOf(getText()), true);
    }

    public String getIID() {
      return jSampleNameTextField.getText();
    }

    public String getFID() {
      return jindexTextField.getText();
    }

    @Override
    public void actionPerformed(ActionEvent e) {
      JComponent source = (JComponent) e.getSource();
      if (source.equals(selectButton)) {
        FileChooser fileChooser = new FileChooser(this, proj.PROJECT_DIRECTORY.getValue(), false,
                                                  false, SELECT_TITLE, log);
        if (fileChooser.getFiles() != null && fileChooser.getFiles().length > 0) {
          jFileNameTextField.setText(fileChooser.getFiles()[0]);
        }
      } else if (source.equals(removeButton)) {
        int response = 1;
        String[] overwriteOptions = new String[] {"Remove", "Cancel"};
        response =
            JOptionPane.showOptionDialog(null, "Remove " + jFileNameTextField.getText() + "?",
                                         "Cancel", JOptionPane.DEFAULT_OPTION,
                                         JOptionPane.QUESTION_MESSAGE, null, overwriteOptions,
                                         overwriteOptions[1]);
        if (response == 0) {
          setVisible(false);
          used = false;
          jframe.pack();
        }
      }
    }
  }
}
