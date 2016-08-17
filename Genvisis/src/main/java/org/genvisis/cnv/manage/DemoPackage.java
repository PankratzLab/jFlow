package org.genvisis.cnv.manage;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.concurrent.Callable;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import org.genvisis.cnv.LaunchProperties;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.FileChooser;
import org.genvisis.cnv.manage.DemoProject.DEMO_TYPE;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

public class DemoPackage {
  private static class DemoPackageWorker implements Callable<Boolean>, Runnable {
    private final Project proj;
    private final String markersFile;
    private final String samplesFile;
    String demoDirectory;
    private final int numThreads;
    private final boolean overwriteExisting;

    public DemoPackageWorker(Project proj, String markersFile, String samplesFile,
                             String demoDirectory, int numThreads, boolean overwriteExisting) {
      super();
      this.proj = proj;
      this.markersFile = markersFile;
      this.samplesFile = samplesFile;
      this.demoDirectory = demoDirectory;
      this.numThreads = numThreads;
      this.overwriteExisting = overwriteExisting;
    }

    @Override
    public Boolean call() throws Exception {
      generateDemoPackage(proj, demoDirectory, markersFile, samplesFile, numThreads,
                          overwriteExisting);
      return true;
      // TODO Auto-generated method stub
    }

    @Override
    public void run() {
      generateDemoPackage(proj, demoDirectory, markersFile, samplesFile, numThreads,
                          overwriteExisting);
      // TODO Auto-generated method stub

    }

  }

  private static class SelectButton extends JButton implements ActionListener {
    private static final long serialVersionUID = 1L;
    private final JTextField text;
    private final Project proj;
    private final JLabel exists;
    private final JFrame jFrame;

    public SelectButton(JFrame jFrame, Project proj, String title, JTextField text, JLabel exists) {
      super(title);
      addActionListener(this);
      this.text = text;
      this.text.getDocument().addDocumentListener(new DocumentListener() {

        @Override
        public void changedUpdate(DocumentEvent e) {
          updateText();
        }

        @Override
        public void insertUpdate(DocumentEvent e) {
          updateText();
        }

        @Override
        public void removeUpdate(DocumentEvent e) {
          updateText();
        }
      });
      this.proj = proj;
      this.exists = exists;
      this.jFrame = jFrame;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
      FileChooser fileChooser = new FileChooser(null, proj.PROJECT_DIRECTORY.getValue(), true,
                                                false, getText(), proj.getLog());
      fileChooser.setApproveButtonToolTipText("Note this file will only be applied to the "
                                              + DemoProject.DEMO_TYPE.SAMPLE_FOCUS + " demo set");
      if (fileChooser.isSelected()) {
        text.setText(fileChooser.getFiles()[0]);
        updateText();
        jFrame.pack();
      }

    }

    public String getCurrentFile() {
      String file = null;
      if (text.getText() != null && !text.getText().equals("")) {
        file = text.getText();
        if (!Files.exists(file)) {
          file = null;
        }
      }
      return file;
    }

    private void updateText() {
      if (Files.exists(text.getText())) {
        exists.setText("File exists");
      } else {
        exists.setText("File does not exist");
      }
    }
  }

  public static final String GENERATE_DEMO_PACKAGE = "Generate a demo package";

  public static void generateDemoPackage(Project proj, String demoDirectory, String markersFile,
                                         String samplesFile, int numThreads,
                                         boolean overwriteExisting) {

    demoDirectory =
        demoDirectory == null ? proj.DEMO_DIRECTORY.getValue(true, false) : demoDirectory;
    proj.getLog().reportTimeInfo("Using demo directory " + demoDirectory);
    DemoPackage demoPackage = new DemoPackage(demoDirectory, proj.getLog());
    if (!demoPackage.isFail()) {
      demoPackage.generateDemoProject(proj, markersFile, samplesFile, numThreads,
                                      overwriteExisting);
    }
    proj.getLog().reportTimeInfo("Finished Generating a demo in " + demoDirectory);

  }

  public static void guiAccess(Project proj) {
    String fileExists = "File exists";
    String fileDoesNotExist = "File does not exist";

    // String samplesFile = !Files.exists(proj.getFilename(proj.SAMPLE_SUBSET_FILENAME)) ? null :
    // proj.getFilename(proj.SAMPLE_SUBSET_FILENAME);
    // String markersFile = !Files.exists(proj.getFilename(proj.TARGET_MARKERS_FILENAME)) ? null :
    // proj.getFilename(proj.TARGET_MARKERS_FILENAME);
    String samplesFile =
        !Files.exists(proj.SAMPLE_SUBSET_FILENAME.getValue()) ? null
                                                              : proj.SAMPLE_SUBSET_FILENAME.getValue();
    String markersFile =
        !Files.exists(proj.TARGET_MARKERS_FILENAMES.getDefaultValueString()) ? null
                                                                             : proj.TARGET_MARKERS_FILENAMES.getDefaultValueString();
    JLabel sampleExists = new JLabel(samplesFile == null ? fileDoesNotExist : fileExists);
    JLabel markerExists = new JLabel(markersFile == null ? fileDoesNotExist : fileExists);

    int width = 100;
    final class MainPanel extends JPanel implements ActionListener {
      private static final long serialVersionUID = 1L;

      @Override
      public void actionPerformed(ActionEvent e) {
        // TODO Auto-generated method stub

      }

    }

    final class GoButton extends JButton implements ActionListener {
      private static final long serialVersionUID = 1L;
      private final SelectButton markerButton;
      private final SelectButton sampButton;
      private final Project proj;
      private final JFrame jFrame;

      public GoButton(Project proj, SelectButton markerButton, SelectButton sampButton,
                      JFrame jFrame) {
        super("Create Demo");
        this.markerButton = markerButton;
        this.sampButton = sampButton;
        this.proj = proj;
        this.jFrame = jFrame;
        addActionListener(this);
      }

      @Override
      public void actionPerformed(ActionEvent e) {
        // new Thread(new DemoPackageWorker(proj, markerButton.getCurrentFile(),
        // sampButton.getCurrentFile(), null, proj.getInt(proj.NUM_THREADS), true)).start();
        new Thread(new DemoPackageWorker(proj, markerButton.getCurrentFile(),
                                         sampButton.getCurrentFile(), null,
                                         proj.NUM_THREADS.getValue(), true)).start();
        jFrame.pack();

      }
    }
    JFrame jFrame = new JFrame();
    Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
    jFrame.setLocation(dim.width / 2 - jFrame.getSize().width / 2,
                       dim.height / 2 - jFrame.getSize().height / 2);
    jFrame.setSize(width, 200);
    jFrame.setVisible(true);
    jFrame.setTitle("Generate Demo Project");

    MainPanel mainPanel = new MainPanel();
    mainPanel.setLayout(new BorderLayout());

    JPanel filePanel = new JPanel(new GridLayout(2, 0));
    // JTextField sampFileText = new
    // JTextField(proj.getFilename(proj.SAMPLE_SUBSET_FILENAME).replaceAll("\"", ""));
    JTextField sampFileText =
        new JTextField(proj.SAMPLE_SUBSET_FILENAME.getValue().replaceAll("\"", ""));
    sampFileText.setSize(width, 30);
    // JTextField markFileText = new
    // JTextField(proj.getFilename(proj.TARGET_MARKERS_FILENAME).replaceAll("\"", ""));
    JTextField markFileText = new JTextField(proj.TARGET_MARKERS_FILENAMES.getValue()[0]);
    filePanel.add(markFileText, BorderLayout.NORTH);
    filePanel.add(sampFileText, BorderLayout.SOUTH);

    JPanel selectPanel = new JPanel(new GridLayout(2, 2));
    SelectButton markerButton =
        new SelectButton(jFrame, proj, "Select Marker File", markFileText, markerExists);
    selectPanel.add(markerButton, BorderLayout.NORTH);
    selectPanel.add(markerExists);
    SelectButton sampButton =
        new SelectButton(jFrame, proj, "Select Sample File", sampFileText, sampleExists);
    selectPanel.add(sampButton, BorderLayout.SOUTH);
    selectPanel.add(sampleExists);
    GoButton goButton = new GoButton(proj, markerButton, sampButton, jFrame);

    mainPanel.add(filePanel, BorderLayout.WEST);
    mainPanel.add(selectPanel, BorderLayout.EAST);
    mainPanel.add(goButton, BorderLayout.SOUTH);

    jFrame.add(mainPanel);
    jFrame.pack();

    // FileChooser fileChooser = new FileChooser(null, proj.getProjectDir(), true, false, "Select a
    // file with a subset of samples to export", proj.getLog());
    // fileChooser.setApproveButtonToolTipText("Note this file will only be applied to the " +
    // DemoProject.DEMO_TYPE.SAMPLE_FOCUS + " demo set");
    // String samplesFile = fileChooser.getSelectedFile().toString();
    // System.out.println(samplesFile);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    String markersFile = null;
    String samplesFile = null;
    String demoDirectory = null;
    int numThreads = 0;
    boolean overwriteExisting = false;
    String logfile = null;
    Project proj;

    String usage = "\n" + "cnv.manage.GenerateDemoPackage requires 0-1 arguments\n";
    usage += "   (1) project properties filename (i.e. proj="
             + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n" + "";
    usage += "   OPTIONAL";
    usage +=
        "   (2) full path to a file with markers to export for marker focused project (i.e. markersFile= ( default exports \"TARGET_MARKERS_FILENAME\", and all markers if targets does not exist))\n"
             + "";
    usage +=
        "   (3) full path to a file of samples to export for sample focused project (i.e. samplesFile= ( default exports \"Project.SAMPLE_SUBSET_FILENAME\" and all samples if subset file does not exist))\n"
             + "";
    usage += "   (4) number of threads (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads
             + " ( defaults to \"NUM_THREADS\" ))\n" + "";
    usage +=
        "   (5) force overwrite option for existing demos (i.e. -overwriteExisting (not the default))\n"
             + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("markersFile=")) {
        markersFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("samplesFile=")) {
        samplesFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("demoDir=")) {
        demoDirectory = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numThreads = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("-overwriteExisting")) {
        overwriteExisting = true;
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
      proj = new Project(filename, logfile, false);

      // guiAccess(proj);
      generateDemoPackage(proj, demoDirectory, markersFile, samplesFile, numThreads,
                          overwriteExisting);

      // makeUpScatter(proj, proj.getFilename(Project.DISPLAY_MARKERS_FILENAME));
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void makeUpScatter(Project proj, String filenameOfMarkerList) {
    if (new File("demo/").exists()) {
      new File("demo/").renameTo(new File(Files.backup("demo", "", "")));
      // new File("projects/").renameTo(new File(Files.backup("demo/projects", "", "")));
      new File("projects/").delete();
    }
    new File("demo/data/").mkdirs();
    new File("demo/transposed/").mkdirs();
    new File("projects/").mkdirs();
    Files.copyFile(filenameOfMarkerList, "demo/data/test.txt");
    // Files.copyFile(proj.getFilename(proj.SAMPLE_DATA_FILENAME), "demo/data/" +
    // ext.removeDirectoryInfo(proj.getFilename(proj.SAMPLE_DATA_FILENAME)));
    Files.copyFile(proj.SAMPLE_DATA_FILENAME.getValue(),
                   "demo/data/" + ext.removeDirectoryInfo(proj.SAMPLE_DATA_FILENAME.getValue()));

    // proj.setJarStatus(true);
    proj.JAR_STATUS.setValue(true);
    proj.setProperty(proj.PROJECT_DIRECTORY, "demo/");
    proj.setProperty(proj.SOURCE_DIRECTORY, "demo/");
    proj.saveProperties("projects/" + ext.removeDirectoryInfo(proj.getPropertyFilename()));
  }

  private final String demoDirectory;

  private final Logger log;

  private LaunchProperties newLaunchProperties;

  private boolean fail;

  public DemoPackage(String demoDirectory, Logger log) {
    super();
    this.log = log;
    this.demoDirectory = demoDirectory;
    fail = false;
    initDemoPackage();
  }

  /**
   * @param proj
   * @param markersFile currently only used for {@link DemoProject.DEMO_TYPE#MARKER_FOCUS}
   * @param samplesFile currently only used for {@link DemoProject.DEMO_TYPE#SAMPLE_FOCUS}
   * @param numThreads number of threads for export and subsets
   * @param overwriteExisting overwrite existing demo projects Note:
   */
  public void generateDemoProject(Project proj, String markersFile, String samplesFile,
                                  int numThreads, boolean overwriteExisting) {
    // numThreads = numThreads > 0 ? numThreads : proj.getInt(proj.NUM_THREADS);
    numThreads = numThreads > 0 ? numThreads : proj.NUM_THREADS.getValue();

    DemoProject demoProjectMarkerFocus =
        new DemoProject(proj, demoDirectory, overwriteExisting, DEMO_TYPE.MARKER_FOCUS);
    demoProjectMarkerFocus.createProjectDemo(markersFile, null, numThreads);// all samples for these
                                                                            // markers...
    saveDemoProperties(proj, demoProjectMarkerFocus, false);

    DemoProject demoProjectSampleFocus =
        new DemoProject(proj, demoDirectory, overwriteExisting, DEMO_TYPE.SAMPLE_FOCUS);
    demoProjectSampleFocus.createProjectDemo(null, samplesFile, numThreads);// all markers for these
                                                                            // samples...
    saveDemoProperties(proj, demoProjectSampleFocus, true);
  }

  private void initDemoPackage() {
    new File(demoDirectory).mkdirs();
    String demoPark = demoDirectory + org.genvisis.common.PSF.Java.GENVISIS;
    String demoVis = demoDirectory + org.genvisis.common.PSF.Java.GENVISIS;
    String runningJar = getClass().getProtectionDomain().getCodeSource().getLocation().getFile();
    if (runningJar.endsWith(org.genvisis.common.PSF.Java.GENVISIS)) {
      String copyRunning;
      String copyOther;
      if (runningJar.endsWith(org.genvisis.common.PSF.Java.GENVISIS)) {
        copyRunning = demoPark;
        copyOther = demoVis;
      } else {
        copyRunning = demoVis;
        copyOther = demoPark;
      }
      if (!Files.exists(demoVis)) {
        if (!Files.exists(copyRunning)
            && !runningJar.endsWith(org.genvisis.common.PSF.Java.GENVISIS)) {
          log.reportTimeInfo("Detected " + runningJar + ", copying to " + copyRunning
                             + "\n\t (this takes a while due to byte by byte copying)");
          if (Files.copyFileUsingFileChannels(runningJar, copyRunning, log)) {
            log.reportTimeInfo("Finished copying " + runningJar + ", to " + copyRunning);
          } else {
            log.reportTimeError("Could not copy " + runningJar + " to " + copyRunning);
            fail = true;
          }
        }

        String other =
            ext.parseDirectoryOfFile(runningJar, false) + ext.removeDirectoryInfo(copyOther);
        if (Files.exists(other) && !Files.exists(copyOther)
            && !other.endsWith(org.genvisis.common.PSF.Java.GENVISIS)) {
          if (Files.exists(other)) {
            log.reportTimeInfo("Detected " + other + ", copying to " + copyOther
                               + "\n\t (this takes a while due to byte by byte copying)");
            if (Files.copyFileUsingFileChannels(other, copyOther, log)) {
              log.reportTimeInfo("Finished copying " + other + ", to " + copyOther);
            } else {
              log.reportTimeError("Could not copy " + other + " to " + copyOther);
              fail = true;
            }
          } else {
            log.reportTimeError("Did not detect " + other + " , halting ");
            fail = true;
          }
        }
      }
    } else {
      log.reportTimeError("Could not detect proper jar file, found " + runningJar
                          + " and it should have ended with "
                          + org.genvisis.common.PSF.Java.GENVISIS);
      log.reportTimeError("This could be because you are running from eclipse without a jar file");
      // fail = true;
    }
    String launchProperties = demoDirectory + LaunchProperties.DEFAULT_PROPERTIES_FILE;
    org.genvisis.cnv.Launch.initLaunchProperties(launchProperties, true, true);

    newLaunchProperties = new LaunchProperties(launchProperties);
    System.out.println(newLaunchProperties.getProperty(LaunchProperties.PROJECTS_DIR));
  }

  public boolean isFail() {
    return fail;
  }

  private void saveDemoProperties(Project proj, DemoProject demoProject, boolean setToDefault) {
    String projectsDir = newLaunchProperties.getDirectory();
    if (Files.isRelativePath(projectsDir)) {
      projectsDir = ext.parseDirectoryOfFile(newLaunchProperties.getFilename()) + projectsDir;
    }
    demoProject.setProperty(demoProject.PROJECT_NAME,
                            proj.PROJECT_NAME.getValue() + "_" + demoProject.getdType());
    String newProjectFile = projectsDir + demoProject.PROJECT_NAME.getValue() + ".properties";
    demoProject.setProperty(demoProject.PROJECT_DIRECTORY,
                            demoProject.PROJECT_NAME.getValue() + "/");
    Files.copyFile(proj.getPropertyFilename(), newProjectFile);
    Files.copyFile(proj.getPropertyFilename(), projectsDir + "example.properties");

    if (!demoProject.isFail()) {
      try {
        // PrintWriter writer = new PrintWriter(new FileWriter(newProjectFile));
        // writer.println("##");
        // writer.close();

        // demoProject.store(writer, demoProject.getNameOfProject());
        demoProject.setPropertyFilename(newProjectFile);
        demoProject.saveProperties();
      } catch (Exception e) {
        log.reportError("Error writing to " + newProjectFile);
        log.reportException(e);
      }
    }

    if (setToDefault) {
      newLaunchProperties.setProperty(LaunchProperties.LAST_PROJECT_OPENED,
                                      demoProject.PROJECT_NAME.getValue() + ".properties");
      newLaunchProperties.save();
    }
  }

  // Files.copyFileFromJar(other, copyOther);
  // Files.copyFile(other, copyOther);
  //
  // FilesNIO.copyFile(other, copyOther, log);

  // Files.copyFile(other, copyOther);
  // Files.copyFile(from, to)
  // Files.copyFileFromJar(runningJar, copyRunning);
  //
  // String bat = "Launch.bat";
  // String sh = "../cnv/" + LAUNCH + SH;
  // System.out.println(bat + "\t" + sh);
  // //InputStream is = cnv.getClass().getResourceAsStream(bat);
  // // try {
  // // //System.out.println(is.read());
  // // } catch (IOException e) {
  // // // TODO Auto-generated catch block
  // // e.printStackTrace();
  // // }
  //
  // String DEFAULT_LAUNCH_SH = "/cnv/Launch.txt";
  // InputStream is =this.getClass().getResourceAsStream(DEFAULT_LAUNCH_SH);
  // InputStreamReader isr = new InputStreamReader(is);
  // BufferedReader reader = new BufferedReader(isr);
  // try {
  // while (reader.ready()){
  // System.out.println(reader.readLine());
  // }
  // } catch (IOException e) {
  // // TODO Auto-generated catch block
  // e.printStackTrace();
  // }
  //
  // Files.copyFileFromJar(DEFAULT_LAUNCH_SH, demoDirectory + LAUNCH + SH);
  // //Files.copyFileFromJar(sh, demoDirectory + LAUNCH + SH);
}
