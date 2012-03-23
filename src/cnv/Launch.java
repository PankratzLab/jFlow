package cnv;

import java.io.*;
//import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import common.*;
//import cnv.analysis.Mosaicism;
import cnv.analysis.Mosaicism;
import cnv.filesys.*;
import cnv.manage.*;
import cnv.plots.*;

public class Launch extends JFrame implements ActionListener, WindowListener {
	public static final long serialVersionUID = 1L;
	public static String filename = LaunchProperties.DEFAULT_PROPERTIES_FILE;

	public static final String MAP_FILES = "Map .csv files to IDs";
	public static final String GENERATE_MARKER_POSITIONS = "Generate marker positions file";
	public static final String PARSE_FILES = "Parse .csv files";
	public static final String EXTRACT_PLOTS = "Extract plots";
	public static final String SLIM_PLOTS = "Slim plots";
	public static final String GENERATE_PLINK_FILES = "Generate PLINK files";
	public static final String LRR_SD = "Check LRR stdevs";
	public static final String SCATTER = "Scatter module";
	public static final String QQ = "QQ module";
	public static final String STRAT = "Stratify module";
	public static final String MOSAICISM = "Determine mosaic arms";
	public static final String MOSAIC_PLOT = "Mosaic module";
	public static final String TRAILER = "Trailer module";
	public static final String GENDER_PLOT = "Gender module";
	public static final String EDIT = "Edit";
	public static final String REFRESH = "Refresh";
	public static final String TEST = "Test";

	public static final String[] BUTTONS = {MAP_FILES, GENERATE_MARKER_POSITIONS, PARSE_FILES, EXTRACT_PLOTS, SLIM_PLOTS, GENERATE_PLINK_FILES, LRR_SD, SCATTER, QQ, STRAT, MOSAICISM, MOSAIC_PLOT, GENDER_PLOT, TRAILER, TEST}; 

	private boolean jar;
	private JComboBox projectsBox;
	private String[] projects;
    private LaunchProperties props;
    private String launchPropertiesFile;

	public Launch(String launchPropertiesFile, boolean jar) {
		super("CNVis");
		this.jar = jar;
		this.launchPropertiesFile = launchPropertiesFile;
	}

	public void loadProjects() {
		String[] projectNames;
		
		projects = Files.list(props.getDirectory(), ".properties", false);
		projectNames = new String[projects.length];
		for (int i = 0; i<projectNames.length; i++) {
			projectNames[i] = ext.rootOf(projects[i], true);
        }
		projectsBox.setModel(new DefaultComboBoxModel(projectNames));
	}
	
    public void addComponentsToPane(final Container pane) {
        JPanel projectPanel, buttonPanel;
		GridLayout layout;
		JButton button;
		String lastProjectOpened;
        
        props = new LaunchProperties(launchPropertiesFile);
        projectsBox = new JComboBox();
        loadProjects();
        // In JDK1.4 this prevents action events from being fired when the  up/down arrow keys are used on the dropdown menu
        projectsBox.putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
		
        projectPanel = new JPanel();
        projectPanel.add(projectsBox);
        
        button = new JButton(Grafik.getImageIcon("images/edit.gif", true));
        button.setPreferredSize(new Dimension(button.getIcon().getIconWidth(), button.getIcon().getIconHeight()));
        button.setBorder(null);
        button.setActionCommand(EDIT);
    	button.addActionListener(this);
        projectPanel.add(button);
        button = new JButton(Grafik.getImageIcon("images/refresh.gif", true));
        button.setPreferredSize(new Dimension(button.getIcon().getIconWidth(), button.getIcon().getIconHeight()));
        button.setBorder(null);
        button.setActionCommand(REFRESH);
    	button.addActionListener(this);
        projectPanel.add(button);
        
        projectPanel.setBackground(Color.WHITE);
        
		lastProjectOpened = props.getProperty(LaunchProperties.LAST_PROJECT_OPENED);
		for (int i = 0; i<projects.length; i++) {
			if (projects[i].equals(lastProjectOpened)) {
				projectsBox.setSelectedIndex(i);
			}
        }

        buttonPanel = new JPanel();
        layout = new GridLayout(0, 1);
        buttonPanel.setLayout(layout);
//        buttonPanel.setPreferredSize(new Dimension(400, 160));
        layout.setHgap(10);
        layout.setVgap(10);
        
        for (int i = 0; i<BUTTONS.length; i++) {
        	button = new JButton(BUTTONS[i]);
        	button.addActionListener(this);
        	buttonPanel.add(button);
        }
        
        pane.add(projectPanel, BorderLayout.NORTH);
        pane.add(new JSeparator(), BorderLayout.CENTER);
        pane.add(buttonPanel, BorderLayout.SOUTH);
    }
    
    public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();
		Project proj;
		String dir;		
		int index;
		
		index = projectsBox.getSelectedIndex();
		dir = props.getDirectory();
		proj = new Project(dir+projects[index], jar);
		
		try {
			if (command.equals(MAP_FILES)) {
	//			vis cnv.manage.ParseIllumina -mapFiles
			} else if (command.equals(GENERATE_MARKER_POSITIONS)) {
				cnv.manage.Markers.generateMarkerPositions(proj, "SNP_Map.csv");
			} else if (command.equals(PARSE_FILES)) {
				cnv.manage.ParseIllumina.createFiles(proj, 2);
	//			nohup vis -Xmx1g cnv.manage.ParseIllumina threads=6 proj=current.proj
			} else if (command.equals(EXTRACT_PLOTS)) {
				ExtractPlots.extractAll(proj, 0, !ExtractPlots.containsLRR(proj));
				//ExtractPlots.extractAll(proj, 0, true); // compact if no LRR was provided
			} else if (command.equals(SLIM_PLOTS)) {
				ExtractPlots.breakUpMarkerCollections(proj, Integer.parseInt(proj.getProperty(Project.NUM_MARKERS_PER_FILE)));
	//			nohup vis -Xmx15g -d64 cnv.manage.ExtractPlots per=250 proj=current.proj
			} else if (command.equals(GENERATE_PLINK_FILES)) {
				String filename = ClusterFilterCollection.getClusterFilterFilenameSelection(proj);
				if ( filename==null || (!filename.equals("cancel")) ) {
					cnv.manage.PlinkFormat.createPlink(proj, filename);
					CmdLine.run("plink --file gwas --make-bed --out plink", proj.getProjectDir());
					new File(proj.getProjectDir()+"genome/").mkdirs();
					CmdLine.run("plink --bfile ../plink --freq", proj.getProjectDir()+"genome/");
					CmdLine.run("plink --bfile ../plink --missing", proj.getProjectDir()+"genome/");
		//			vis cnv.manage.PlinkFormat root=../plink genome=6
				}
			} else if (command.equals(LRR_SD)) {
				cnv.qc.LrrSd.init(proj, null, Integer.parseInt(proj.getProperty(Project.NUM_THREADS)));
			} else if (command.equals(SCATTER)) {
				new ScatterPlot(proj);
			} else if (command.equals(QQ)) {
				QQPlot.loadPvals(proj.getFilenames(Project.QQ_FILENAMES, true), Boolean.valueOf(proj.getProperty(Project.DISPLAY_QUANTILES)), Boolean.valueOf(proj.getProperty(Project.DISPLAY_STANDARD_QQ)), Boolean.valueOf(proj.getProperty(Project.DISPLAY_ROTATED_QQ)));
			} else if (command.equals(STRAT)) {
				StratPlot.loadStratificationResults(proj);
			} else if (command.equals(MOSAICISM)) {
				Mosaicism.findOutliers(proj);
			} else if (command.equals(MOSAIC_PLOT)) {
				MosaicPlot.loadMosaicismResults(proj);
			} else if (command.equals(TRAILER)) {
				new Trailer(proj, null, proj.getFilenames(Project.CNV_FILENAMES), Trailer.DEFAULT_LOCATION);
			} else if (command.equals(GENDER_PLOT)) {
				JOptionPane.showMessageDialog(null, "Currently under construction. Please check back later.");
				GenderPlot.loadGenderResults(proj);
			} else if (command.equals(EDIT)) {
				try {
	//				Runtime.getRuntime().exec("C:\\Program Files\\Windows NT\\Accessories\\wordpad.exe \"C:"+ext.replaceAllWith(projects[index], "/", "\\")+"\"");
					Runtime.getRuntime().exec("C:\\Windows\\System32\\Notepad.exe \""+dir+projects[index]+"\"");
					System.out.println("tried to open "+projects[index]+" which "+(new File(dir+projects[index]).exists()?"does":"does not")+" exist");
				} catch (IOException ioe) {
					System.err.println("Error - failed to open WordPad");
				}
			} else if (command.equals(REFRESH)) {
		        loadProjects();
				System.out.println("Refreshed list of projects");
			} else if (command.equals(TEST)) {
				System.out.println("Testing latest subroutine");
//				Mosaicism.checkForOverlap(proj);
				ParseIllumina.parseAlleleLookup(proj);
			} else {
				System.err.println("Error - unknown command '"+command+"'");
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
    
    public void windowOpened(WindowEvent we) {}
    
    public void windowClosing(WindowEvent we) {
    	props.setProperty(LaunchProperties.LAST_PROJECT_OPENED, projects[projectsBox.getSelectedIndex()]);
    	props.save();
    }
    
    public void windowClosed(WindowEvent we) {}
    
    public void windowIconified(WindowEvent we) {}
    
    public void windowDeiconified(WindowEvent we) {}
    
    public void windowActivated(WindowEvent we) {}
    
    public void windowDeactivated(WindowEvent we) {}    

    private static void createAndShowGUI(String launchPropertiesFile) {
    	Launch frame;
    	String path;
    	
    	if (!new File(launchPropertiesFile).exists()) {
    		JOptionPane.showMessageDialog(null, "Could not find file '"+launchPropertiesFile+"'; generating a blank one that you can populate with your own project links", "FYI", JOptionPane.ERROR_MESSAGE);
    		try {
    			path = ext.parseDirectoryOfFile(new File(launchPropertiesFile).getCanonicalPath());
    		} catch (IOException ioe) {
    			path = "";
    		}
    		new File(path+"projects/").mkdirs();
    		Files.writeList(new String[] {"LAST_PROJECT_OPENED=example.properties", "PROJECTS_DIR="+path+"projects/"}, launchPropertiesFile);
        	if (!new File(path+"projects/example.properties").exists()) {
        		Files.writeList(new String[] {"PROJECT_NAME=Example", "PROJECT_DIRECTORY=example/", "SOURCE_DIRECTORY=example/source/"}, path+"projects/example.properties");
        	}
    	}
    	frame = new Launch(launchPropertiesFile, false);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.addComponentsToPane(frame.getContentPane());
        frame.addWindowListener(frame);
        frame.pack();
        frame.setVisible(true);
    }

	public static void main(String[] args) {
		int numArgs = args.length;
		
		String usage = "\n"+
		"cnv.Launch requires 1 argument\n"+
		"   (1) name of file with project locations (i.e. file="+filename+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
            	createAndShowGUI(filename);
            }
        });
	}
}
