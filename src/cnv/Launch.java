package cnv;

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import common.*;
import cyto.CytoGUI;
import cnv.analysis.DeNovoCNV;
import cnv.analysis.Mosaicism;
import cnv.filesys.*;
//import cnv.gui.GuiManager;
//import cnv.gui.PropertyEditor;
import cnv.manage.*;
import cnv.plots.*;
import cnv.qc.MarkerMetrics;

public class Launch extends JFrame implements ActionListener, WindowListener, ItemListener {
	public static final long serialVersionUID = 1L;
	
	public static final String VERSION = "0.60";

	public static final String EXIT = "Exit";
	public static final String EDIT = "Project Properties";
	public static final String REFRESH = "Refresh";

	public static final String MAP_FILES = "Map .csv files to IDs";
	public static final String GENERATE_MARKER_POSITIONS = "Generate marker positions file";
	public static final String PARSE_FILES_CSV = "Parse .csv files";
	public static final String TRANSPOSE_DATA = "Transpose data";
	public static final String KITANDKABOODLE = "Kit and Kaboodle";

	public static final String CHECK_SEX = "Check sex";
	public static final String LRR_SD = "LRR Stdevs";
	public static final String CNP_SCAN = "Scan for CNPs";
	public static final String MOSAICISM = "Determine mosaic arms";
	public static final String MARKER_METRICS = "Full QC marker metrics";
	public static final String FILTER_MARKER_METRICS = "Filter marker metrics";
	public static final String TALLY_MARKER_ANNOTATIONS = "Tally marker annotations";
	public static final String TALLY_WITHOUT_DETERMINING_DROPS = "Tally without determining dropped markers (much faster)";
	
	public static final String SCATTER = "Scatter module";
	public static final String QQ = "QQ module";
	public static final String STRAT = "Stratify module";
	public static final String MOSAIC_PLOT = "Mosaic plot module";
	public static final String SEX_PLOT = "Sex module";
	public static final String TRAILER = "Trailer module";
	public static final String TWOD = "2D Plot";
	public static final String LINE_PLOT = "Line Plot";
	public static final String COMP = "Comp module";

	public static final String GENERATE_ABLOOKUP = "Generate AB Lookup";
	public static final String GENERATE_PLINK_FILES = "Generate PLINK files";
	public static final String GENERATE_PLINK_BINARY_FILES = "Generate PLINK binary files";
	public static final String GENERATE_PENNCNV_FILES = "Generate PennCNV files";
	public static final String PARSE_RAW_PENNCNV_RESULTS = "Parse raw PennCNV results files";
	public static final String POPULATIONBAF = "Compute Population BAF file";
	public static final String GCMODEL = "Compute GC model file";
	public static final String DENOVO_CNV = "De Novo CNV";
	public static final String EXPORT_CNVS = "Export CNVs to Pedfile format";
	public static final String CYTO_WORKBENCH = "Parse workbench files";
	public static final String TEST = "Test new program";
	
	public static String[][] MENUS = {{"File", "Select Project", EDIT, "Preferences", EXIT},
			{"Data", MAP_FILES, GENERATE_MARKER_POSITIONS, PARSE_FILES_CSV, TRANSPOSE_DATA, KITANDKABOODLE},
			{"Quality", CHECK_SEX, LRR_SD, CNP_SCAN, MOSAICISM, MARKER_METRICS, FILTER_MARKER_METRICS, TALLY_MARKER_ANNOTATIONS, TALLY_WITHOUT_DETERMINING_DROPS},
			{"Plots", SCATTER, QQ, STRAT, MOSAIC_PLOT, SEX_PLOT, TRAILER, TWOD, LINE_PLOT, COMP},
			{"Tools", GENERATE_ABLOOKUP, GENERATE_PLINK_FILES, GENERATE_PLINK_BINARY_FILES, GENERATE_PENNCNV_FILES, PARSE_RAW_PENNCNV_RESULTS, POPULATIONBAF, GCMODEL, DENOVO_CNV, EXPORT_CNVS, CYTO_WORKBENCH, TEST},
			{"Help", "Contents", "Search", "About"}};

	
	private Project proj;
	private boolean jar;
	private JComboBox<String> projectsBox;
	private String[] projects;
    private LaunchProperties launchProperties;
    private String launchPropertiesFile;
    private JTextArea output;
    private JScrollPane scrollPane;
    private Vector<Thread> threadsRunning;
    private int indexOfCurrentProj;
    private long timestampOfPropertiesFile;
    private long timestampOfSampleDataFile;
    private Logger log;

	public Launch(String launchPropertiesFile, boolean jar) {
		super("Genvisis");
		this.jar = jar;
		this.launchPropertiesFile = launchPropertiesFile;
		this.timestampOfPropertiesFile = -1;
		this.timestampOfSampleDataFile = -1;
		this.threadsRunning = new Vector<Thread>();
	}

	public void loadProjects() {
		String[] projectNames;
		
		projects = Files.list(launchProperties.getDirectory(), ".properties", false);
		projectNames = new String[projects.length];
		for (int i = 0; i<projectNames.length; i++) {
			projectNames[i] = ext.rootOf(projects[i], true);
        }
		projectsBox.setModel(new DefaultComboBoxModel<String>(projectNames));
	}
	

	public void loadProject() {
		proj = new Project(launchProperties.getDirectory() + projects[indexOfCurrentProj], jar);
		proj.setGuiState(true);
		timestampOfPropertiesFile = new Date().getTime();
		timestampOfSampleDataFile = new Date().getTime();
		if (!Files.exists(proj.getProjectDir(), proj.getJarStatus())) {
			JOptionPane.showMessageDialog(null, "Error - the directory ('"+proj.getProjectDir()+"') for project '"+proj.getNameOfProject()+"' did not exist; creating now. If this was in error, please edit the property file.", "Error", JOptionPane.ERROR_MESSAGE);
		}
		
		log = proj.getLog();
	    log.linkTextArea(output);
	}

	public void setIndexOfCurrentProject(String projPropertiesFileName) {
		indexOfCurrentProj = 0;
		for (int i = 0; i < projects.length; i++) {
			if (projects[i].equals(projPropertiesFileName)) {
				indexOfCurrentProj = i;
			}
		}
	}
	
	@SuppressWarnings("rawtypes")
	public static void setUIFont (Font newFont){
		Enumeration keys = UIManager.getDefaults().keys();
		while (keys.hasMoreElements()) {
			Object key = keys.nextElement();
			Object value = UIManager.get (key);
			if (value != null && value instanceof javax.swing.plaf.FontUIResource) {
				Font oldFont = UIManager.getFont(key);
				UIManager.put(key, newFont.deriveFont(oldFont.getStyle(), oldFont.getSize2D()));
			}
		}
	}
	
	private static void createAndShowGUI() {
    	String launchPropertiesFile;
        Launch frame;
    	String path;
    	
//		try {
//			UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
//		} catch (Exception e) {
//			System.err.println("Failed loading LookandFeel: com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
//			System.err.println(e);
//			try {
//				UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
//			} catch (Exception e2) {
//				System.err.println("Failed loading CrossPlatformLookAndFeel");
//				System.err.println(e2);
//			}
//		}
//		
//    	try {
//    		setUIFont(Font.createFont(Font.TRUETYPE_FONT, ClassLoader.getSystemResourceAsStream("cnv/gui/LiberationSansRegular.ttf")));
//		} catch (Exception e) {
//			System.err.println("Error - failed to load custom font 'LiberationSans'");
//			e.printStackTrace();
//		}

		//Create and set up the content pane.
    	launchPropertiesFile = LaunchProperties.DEFAULT_PROPERTIES_FILE;
    	if (!new File(launchPropertiesFile).exists()) {
//			frame.output.append("Could not find file '"+launchPropertiesFile+"'; generating a blank one that you can populate with your own project links");
			try {
				path = ext.parseDirectoryOfFile(new File(launchPropertiesFile).getCanonicalPath());
			} catch (IOException ioe) {
				path = "";
			}
			new File(path+"projects/").mkdirs();
			new File(path+"example/").mkdirs();
			Files.writeList(new String[] {"LAST_PROJECT_OPENED=example.properties", "PROJECTS_DIR="+path+"projects/"}, launchPropertiesFile);
	    	if (!new File(path+"projects/example.properties").exists()) {
	    		Files.writeList(new String[] {"PROJECT_NAME=Example", "PROJECT_DIRECTORY=example/", "SOURCE_DIRECTORY=sourceFiles/"}, path+"projects/example.properties");
	    	}
    	}
    	frame = new Launch(launchPropertiesFile, false);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		frame.setJMenuBar(frame.topMenuBar());
        frame.setContentPane(frame.createContentPane());
        frame.createPopupMenu();

		frame.setSize(650, 550);
		frame.setLocation(300,200);

//		frame.pack();
//		frame.setBounds(20, 20, 1000, 600);
		frame.addWindowListener(frame);
//		frame.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH);
		frame.setVisible(true);

		// TODO only instantiate when used
		frame.setIndexOfCurrentProject(frame.launchProperties.getProperty(LaunchProperties.LAST_PROJECT_OPENED));
		frame.loadProject();
    }

    public Container createContentPane() {
	    //Create the content-pane-to-be.
	    JPanel contentPane = new JPanel(new BorderLayout());
	    contentPane.setOpaque(true);
	
	    contentPane.add(topIconBar(), BorderLayout.NORTH);
	
	    //Create a scrolled text area.
	    output = new JTextArea(5, 30);
	    output.setEditable(false);
	    scrollPane = new JScrollPane(output);
//	    log.linkTextArea(output);
//	    log.report("Genvisis, v0.60\n(c)2012 Nathan Pankratz, GNU General Public License, v2\n\n"+(new Date()));
//	    output.append("Genvisis, v0.60\n(c)2012 Nathan Pankratz, GNU General Public License, v2\n\n"+(new Date()));
	
	    //Add the text area to the content pane.
	    contentPane.add(scrollPane, BorderLayout.CENTER);
	
	    return contentPane;
	}

	private JMenuBar topMenuBar() {
		JMenuBar menuBar;
		JMenu menu, submenu;
		JMenuItem menuItem;
		Hashtable<Character,String> hash;
					
        launchProperties = new LaunchProperties(launchPropertiesFile);
		menuBar = new JMenuBar();
		for (int i=0; i<MENUS.length; i++) {
			menu = new JMenu(MENUS[i][0]);
	        menu.setMnemonic((int)MENUS[i][0].charAt(0));
			menuBar.add(menu);
			hash = new Hashtable<Character, String>();
			for (int j=1; j<MENUS[i].length; j++) {
				if (MENUS[i][j]=="") {
					break;
				}
				if (MENUS[i][j]=="1") {
					menu.addSeparator();
				} else if (MENUS[i][j].equals("Select Project")) {
					submenu = new JMenu(MENUS[i][j]);
//			        submenu.setMnemonic(KeyEvent.VK_S);

			        menuItem = new JMenuItem("New");
//			        menuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_2, ActionEvent.ALT_MASK));
			        menuItem.addActionListener(this);
			        submenu.add(menuItem);
					projects = Files.list(launchProperties.getDirectory(), ".properties", false);
					for (int k = 0; k<projects.length; k++) {
				        menuItem = new JMenuItem(ext.rootOf(projects[k], true)+" ");
//				        menuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_2, ActionEvent.ALT_MASK));
				        menuItem.addActionListener(this);
				        submenu.add(menuItem);
			        }
//					projectsBox.setModel(new DefaultComboBoxModel<String>(projectNames));
					menu.add(submenu);
				} else {
					menuItem = new JMenuItem(MENUS[i][j]);
					for (int k = 0; k < MENUS[i][j].length(); k++) {
						if (!hash.containsKey(MENUS[i][j].toLowerCase().charAt(k))) {
							menuItem.setMnemonic((int)MENUS[i][j].toLowerCase().charAt(k));
							hash.put(MENUS[i][j].toLowerCase().charAt(k), "");
							k = MENUS[i][j].length();
						}
					}
//					menuItem.setAccelerator(KeyStroke.getKeyStroke((int)(j+"").charAt(0), ActionEvent.ALT_MASK));
					menuItem.addActionListener(this);
					//TODO What is this?
					menuItem.getAccessibleContext().setAccessibleDescription("Under development");
	//				menuItem = new JMenuItem(menus[i][j], KeyEvent.VK_O);
					menu.add(menuItem);
				}
			}
		}
		return menuBar;
	}

	private JPanel topIconBar() {
		JPanel iconBar;
		JButton button;
		String[] icons = {"images/save1.png", "images/edit1.png", "images/refresh.gif", "images/scatterPlot2.png", "images/trailerPlot2.png", "images/qqplot.gif", "images/recluster1.png", "images/twoDPlot1.jpg"};
		String[] commands = {"", EDIT, REFRESH, SCATTER, TRAILER, QQ, LINE_PLOT, TWOD};
		
		iconBar = new JPanel();
		iconBar.setLayout(new FlowLayout(FlowLayout.LEFT));
		for (int i=0; i<icons.length; i++) {
			button = new JButton(Grafik.getImageIcon(icons[i], true));
//			button.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif", true));
//			button.addActionListener(this);
			button.setActionCommand(commands[i]);
	    	button.addActionListener(this);
	    	button.setToolTipText(commands[i]);
			button.setPreferredSize(new Dimension(25, 25));
	        button.setBorder(null);
			iconBar.add(button);
		}
        addComponentsToPane(iconBar);

		return iconBar;
	}

	public void addComponentsToPane(final Container pane) {
//		String lastProjectOpened;
        
        launchProperties = new LaunchProperties(launchPropertiesFile);
        projectsBox = new JComboBox<String>();
        loadProjects();
        // In JDK1.4 this prevents action events from being fired when the  up/down arrow keys are used on the dropdown menu
        projectsBox.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
		
		setIndexOfCurrentProject(launchProperties.getProperty(LaunchProperties.LAST_PROJECT_OPENED));
		projectsBox.setSelectedIndex(indexOfCurrentProj);
		projectsBox.addItemListener(this);
        pane.add(projectsBox);
    }

	public void createPopupMenu() {
        JMenuItem menuItem;
 
        //Create the popup menu.
        JPopupMenu popup = new JPopupMenu();
        menuItem = new JMenuItem("A popup menu item");
        menuItem.addActionListener(this);
        popup.add(menuItem);
        menuItem = new JMenuItem("Another popup menu item");
        menuItem.addActionListener(this);
        popup.add(menuItem);
 
        //Add listener to the text area so the popup menu can come up.
//        MouseListener popupListener = new PopupListener(popup);
//        output.addMouseListener(popupListener);
    }

	public class IndependentThread implements Runnable {
		private Project proj;
		private String command;
		
		public IndependentThread(Project proj, String command) {
			this.proj = proj;
			this.command = command;
		}

		@Override
		public void run() {
			if (command.equals(MAP_FILES)) {
				cnv.manage.ParseIllumina.mapFilenamesToSamples(proj, "filenamesMappedToSamples.txt");
			} else if (command.equals(GENERATE_MARKER_POSITIONS)) {
				cnv.manage.Markers.generateMarkerPositions(proj, proj.getLocationOfSNP_Map());
			} else if (command.equals(PARSE_FILES_CSV)) {
				cnv.manage.ParseIllumina.createFiles(proj, proj.getInt(Project.NUM_THREADS));
			} else if (command.equals(CHECK_SEX)) {
				cnv.qc.SexChecks.sexCheck(proj);
			} else if (command.equals(TRANSPOSE_DATA)) {
				TransposeData.transposeData(proj, 2000000000, false);
			} else if (command.equals(GENERATE_ABLOOKUP)) {
				ABLookup abLookup;
				String filename;
				
				filename = proj.getProjectDir()+ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
				if (!Files.exists(filename)) {
					abLookup = new ABLookup();
					abLookup.parseFromOriginalGenotypes(proj);
					abLookup.writeToFile(filename);
				}
				
				ABLookup.fillInMissingAlleles(proj, filename, proj.getLocationOfSNP_Map());
			} else if (command.equals(GENERATE_PLINK_FILES)) {
				String filename;
				boolean success;

				filename = ClusterFilterCollection.getClusterFilterFilenameSelection(proj);
//				System.out.println("using "+filename);
				if (filename == null) {
					log.report("No ClusterFilterCollection will be used");
				} else {
					log.report("The ClusterFilterCollection in '"+proj.getProperty(Project.DATA_DIRECTORY)+"/"+filename+"' will be used");
				}
				if ( filename==null || (!filename.equals("cancel")) ) {
//						String lookupTable = ClusterFilterCollection.getGenotypeLookupTableSelection(proj);
//						if () {
//							
//						}
					//success = PlinkData.saveGenvisisToPlinkBedSet(proj, "plinkZack", filename, -1, true, log);
					success = cnv.manage.PlinkFormat.createPlink(proj, "gwas", filename);
					if (success) {
						try {
							log.report("Converting ped/map files to binary PLINK files...", false, true);
							if (CmdLine.run("plink --file gwas --make-bed --out plink", proj.getProjectDir())) {
								log.report("complete!");
							} else {
								log.report("PLINK conversion failed");
							}
//							CmdLine.run("plink --bfile plink --recode --out gwas_plink_reverse", proj.getProjectDir());
							new File(proj.getProjectDir()+"genome/").mkdirs();
							CmdLine.run("plink --bfile ../plink --freq", proj.getProjectDir()+"genome/");
							CmdLine.run("plink --bfile ../plink --missing", proj.getProjectDir()+"genome/");
//							CmdLine.run("plink --bfile plink --mind 0.1 --geno 0.9 --make-bed --out plinkSlim", proj.getProjectDir());
						} catch (Exception e) {}
					}
		//			vis cnv.manage.PlinkFormat root=../plink genome=6
				}

			} else if (command.equals(GENERATE_PLINK_BINARY_FILES)) {
				String filename;

				filename = ClusterFilterCollection.getClusterFilterFilenameSelection(proj);
				if (filename == null) {
					log.report("No ClusterFilterCollection will be used");
				} else {
					log.report("The ClusterFilterCollection in '"+proj.getProperty(Project.DATA_DIRECTORY)+"/"+filename+"' will be used");
				}

				if (PlinkData.saveGenvisisToPlinkBedSet(proj, "plinkZack", filename, -1, true)) {
					log.report("Success!");
				}
			} else if (command.equals(GENERATE_PENNCNV_FILES)) {
				cnv.analysis.AnalysisFormats.penncnv(proj, proj.getSampleList().getSamples(), null);
			} else if (command.equals(PARSE_RAW_PENNCNV_RESULTS)) {
				// TODO make dialog to ask for filenames with a JCheckBox for denovo parsing
				cnv.analysis.PennCNV.parseWarnings(proj, "penncnv.log");
				cnv.analysis.PennCNV.parseResults(proj, "penncnv.rawcnv", false);
			} else if (command.equals(LRR_SD)) {
				cnv.qc.LrrSd.init(proj, null, null, Integer.parseInt(proj.getProperty(Project.NUM_THREADS)));
			} else if (command.equals(CNP_SCAN)) {
//					TODO Genotyping
//					new ScanForCnp(proj, "CNPScanResult.txt");
			} else if (command.equals(DENOVO_CNV)) {
				DeNovoCNV.main(null);
			} else if (command.equals(SCATTER)) {
				ScatterPlot.createAndShowGUI(proj, null, null, false);
			} else if (command.equals(QQ)) {
				QQPlot.loadPvals(proj.getFilenames(Project.QQ_FILENAMES, true), "Q-Q Plot", Boolean.valueOf(proj.getProperty(Project.DISPLAY_QUANTILES)), Boolean.valueOf(proj.getProperty(Project.DISPLAY_STANDARD_QQ)), Boolean.valueOf(proj.getProperty(Project.DISPLAY_ROTATED_QQ)), -1, false, proj.getFloat(Project.QQ_MAX_NEG_LOG10_PVALUE), proj.getLog());
			} else if (command.equals(STRAT)) {
				StratPlot.loadStratificationResults(proj);
			} else if (command.equals(MOSAICISM)) {
				Mosaicism.findOutliers(proj);
			} else if (command.equals(MOSAIC_PLOT)) {
				MosaicPlot.loadMosaicismResults(proj);
			} else if (command.equals(SEX_PLOT)) {
				SexPlot.loadGenderResults(proj);
			} else if (command.equals(TRAILER)) {
				new Trailer(proj, null, proj.getFilenames(Project.CNV_FILENAMES), Trailer.DEFAULT_LOCATION);
			} else if (command.equals(TWOD)) {
//				TwoDPlot twoDP = 
				TwoDPlot.createAndShowGUI(proj);
				//TODO: Sample call to test this functionality. Should be removed when seems to work fine.
				// replace the filename according to the path on local machine
//				twoDP.showSpecificFile(proj, "/Users/rohitsinha/Documents/development/ra/practice/sexCheck.xln", 5, 9, proj.getLog());
			} else if (command.equals(LINE_PLOT)) {
				LinePlot.createAndShowGUI(proj);
			} else if (command.equals(COMP)) {
				new CompPlot(proj);
			} else if (command.equals(POPULATIONBAF)) {
				cnv.analysis.PennCNV.populationBAF(proj);
			} else if (command.equals(EXPORT_CNVS)) {
				cnv.manage.ExportCNVsToPedFormat.main(null);
			} else if (command.equals(CYTO_WORKBENCH)) {
				new CytoGUI(proj, proj.getProjectDir(), null);
			} else if (command.equals(TEST)) {
//				log.report("No new program to test");
//				ScatterPlot.createAndShowGUI(proj, null, null, false);

				cnv.qc.SexChecks.sexCheck(proj);
				cnv.qc.LrrSd.init(proj, null, null, Integer.parseInt(proj.getProperty(Project.NUM_THREADS)));
				Mosaicism.findOutliers(proj);

				cnv.manage.PlinkFormat.createPlink(proj, "gwas", null);
				CmdLine.run("plink --file gwas --make-bed --out plink", proj.getProjectDir());
				new File(proj.getProjectDir()+"genome/").mkdirs();
				CmdLine.run("plink --bfile ../plink --freq", proj.getProjectDir()+"genome/");
				CmdLine.run("plink --bfile ../plink --missing", proj.getProjectDir()+"genome/");

				
			} else if (command.equals(GCMODEL)) {
				cnv.analysis.PennCNV.gcModel(proj, "/projects/gcModel/gc5Base.txt", "/projects/gcModel/ourResult.gcModel", 100);
			} else if (command.equals(MARKER_METRICS)) {
				cnv.qc.MarkerMetrics.fullQC(proj, null, null);
			} else if (command.equals(FILTER_MARKER_METRICS)) {
				cnv.qc.MarkerMetrics.filterMetrics(proj);
			} else if (command.equals(TALLY_MARKER_ANNOTATIONS)) {
				MarkerMetrics.tallyFlaggedReviewedChangedAndDropped(proj, true);
			} else if (command.equals(TALLY_WITHOUT_DETERMINING_DROPS)) {
				MarkerMetrics.tallyFlaggedReviewedChangedAndDropped(proj, false);
			} else if (command.equals(KITANDKABOODLE)) {
				String filename;

				if (!Files.exists(proj.getFilename(Project.MARKER_POSITION_FILENAME, false, false))) {
					log.reportError("Could not find required file "+proj.getFilename(Project.MARKER_POSITION_FILENAME, false, false)+"\n    attempting to generate one for you...");
					filename = proj.getLocationOfSNP_Map();
					if (filename == null) {
						return;
					}
					log.report("Generating from "+filename);
					cnv.manage.Markers.generateMarkerPositions(proj, filename);
				}
				cnv.manage.ParseIllumina.createFiles(proj, proj.getInt(Project.NUM_THREADS));

				TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
				cnv.qc.LrrSd.init(proj, null, null, Integer.parseInt(proj.getProperty(Project.NUM_THREADS)));
				cnv.qc.SexChecks.sexCheck(proj);

				cnv.manage.PlinkFormat.createPlink(proj, "gwas", null);
				CmdLine.run("plink --file gwas --make-bed --out plink", proj.getProjectDir());
				new File(proj.getProjectDir()+"genome/").mkdirs();
				CmdLine.run("plink --bfile ../plink --freq", proj.getProjectDir()+"genome/");
				CmdLine.run("plink --bfile ../plink --missing", proj.getProjectDir()+"genome/");



			} else {
				log.reportError("Error - unknown command: "+command);
			}
		}

	}
	
	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();
		Thread thread;
		
//		output.append("Action performed: " + command + "\n");
//		output.setCaretPosition(output.getDocument().getLength());
		if (log == null) {
			log = new Logger();
		}
		log.report("Action performed: " + command + "\n");
		
		if (proj == null) {
			log.report("Trying again to load project");
			loadProject();
		} else if (timestampOfPropertiesFile < new File(proj.getPropertyFilename()).lastModified()) {
			log.report("Detected a change in the project properties file; reloading from '"+proj.getPropertyFilename()+"'");
			proj = null;
			loadProject();
		} else {
//			log.report("No change in properties file");
		}

		if (proj != null && timestampOfSampleDataFile < new File(proj.getFilename(Project.SAMPLE_DATA_FILENAME, false, false)).lastModified()) {
			log.report("Detected a change in the sampleData file; reloading sample data");
			proj.resetSampleData();
		}	
		
		if (command.equals(EXIT)) {
			System.exit(0);
		} else if (command.equals(EDIT)) {
//			new PropertyEditor(proj);
			
			int index = projectsBox.getSelectedIndex();
			String dir = launchProperties.getDirectory();
			try {
				if (System.getProperty("os.name").startsWith("Windows")) {
					Runtime.getRuntime().exec("C:\\Windows\\System32\\Notepad.exe \""+dir+projects[index]+"\"");
					if (!new File(dir+projects[index]).exists()) {
						log.report("Tried to open "+projects[index]+" which does not exist");
					}
				} else {
					JOptionPane.showMessageDialog(null, "This button currently only works for the Windows operating system; a full feature property editor will be arriving in due course", "Sorry", JOptionPane.ERROR_MESSAGE);
				}
			} catch (IOException ioe) {
				log.reportError("Error - failed to open Notepad");
			}
		} else if (command.equals(REFRESH)) {
	        loadProjects();
			log.report("Refreshed list of projects");
		} else if (command.endsWith(" ")) {
			for (int i=0; i<projects.length; i++) {
				if (command.equals(ext.rootOf(projects[i])+" ")) {
					projectsBox.setSelectedIndex(i);
					log.report("Selecting: "+projects[i]);
				}
			}
		} else {
			thread = new Thread(new IndependentThread(proj, command));
			thread.start();
			threadsRunning.add(thread);
		}
	}
    
    public void windowOpened(WindowEvent we) {}
    
    public void windowClosing(WindowEvent we) {
//    	boolean alive;
//    	
////    	notify all threads that they need to close
//    	alive = false;
//    	for (int i = 0; i < threadsRunning.size(); i++) {
//    		if (threadsRunning.elementAt(i) != null) {
//    			threadsRunning.elementAt(i).interrupt();
//    		}
//		}
//
//    	alive = false;
//    	for (int i = 0; i < threadsRunning.size(); i++) {
//    		if (threadsRunning.elementAt(i) != null && threadsRunning.elementAt(i).isAlive()) {
//    			alive = true;
//    		}
//		}
//    	
//    	if (!alive) {
//    		GuiManager.disposeOfParentFrame(this);
//    	}
    }
    
    public void windowClosed(WindowEvent we) {}
    
    public void windowIconified(WindowEvent we) {}
    
    public void windowDeiconified(WindowEvent we) {}
    
    public void windowActivated(WindowEvent we) {}
    
    public void windowDeactivated(WindowEvent we) {}

    public void itemStateChanged(ItemEvent e) {
//    	byte i;
    	if (e.getStateChange()==ItemEvent.SELECTED) {
    		indexOfCurrentProj = projectsBox.getSelectedIndex();
//			proj = new Project(launchProperties.getDirectory() + projects[indexOfCurrentProject], jar);
			loadProject();
//			output.append("\nCurrent project: " + ext.rootOf(projects[indexOfCurrentProj]) + "\n");
//			output.setCaretPosition(output.getDocument().getLength());
			log.report("\nCurrent project: " + ext.rootOf(projects[indexOfCurrentProj]) + "\n");

			launchProperties.setProperty(LaunchProperties.LAST_PROJECT_OPENED, projects[projectsBox.getSelectedIndex()]);
			launchProperties.save();
    	}
    }
    
	public static String getDefaultDebugProjectFile(boolean verbose) {
		LaunchProperties launchProperties;
		String dir, filename;
		
		launchProperties = new LaunchProperties(LaunchProperties.DEFAULT_PROPERTIES_FILE);
		dir = launchProperties.getProperty(LaunchProperties.PROJECTS_DIR);
		filename = launchProperties.getProperty(LaunchProperties.DEBUG_PROJECT_FILENAME);
		if (dir == null || filename == null) {
			if (verbose) {
				System.err.println("Warning - you are trying to access the default debug project properties file, but there is no '"+LaunchProperties.DEBUG_PROJECT_FILENAME+"=' property listed in '"+LaunchProperties.DEFAULT_PROPERTIES_FILE+"'. The default filename is being set to \"default.properties\" in the current directory. However, if that does not exist either, then the program will likely end in an error.");
			}
			dir = "./";
			filename = "default.properties";
		} else if (!Files.exists(dir) || !Files.exists(dir+filename)) {
			if (verbose) {
				System.err.println("Error - default debug project properties file does not exist: "+dir+filename);
			}
		} else {
			if (verbose) {
				System.out.println("The default debug project properties file is currently set to '"+dir+filename+"'");
			}
		}
		
		return dir+filename;
	}

	public static void main(String[] args) {
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
            	createAndShowGUI();
            }
        });
	}
}
