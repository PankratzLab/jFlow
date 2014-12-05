// replace all exits with a proper disposal of all resources
package cnv.plots;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;

import stats.CTable;
import stats.ContingencyTable;
import stats.ProbDist;

import java.awt.*;
import java.awt.event.*;

import cnv.analysis.pca.PrincipalComponentsIntensity;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.*;
import cnv.gui.AnnotationAction;
import cnv.gui.AutoSaveForScatterPlot;
import cnv.gui.ColorKeyPanel;
import cnv.gui.CycleRadio;
import cnv.manage.MarkerDataLoader;
import common.*;
import cnv.var.*;

public class ScatterPlot extends JPanel implements ActionListener, WindowListener {
	public static final long serialVersionUID = 1L;
	public static final byte DEFAULT_SIZE = 8;
	public static final int DEFAULT_GC_THRESHOLD = 15;
	private static final String ALT_UP = "ALT UP";
	private static final String ALT_DOWN = "ALT DOWN";
	private static final String ALT_LEFT = "ALT LEFT";
	private static final String ALT_RIGHT = "ALT RIGHT";
	private static final String FIRST = "First";
	private static final String PREVIOUS = "Previous";
	private static final String NEXT = "Next";
	private static final String LAST = "Last";
	private static final String SYMMETRY = "Symmetry";
	private static final String CORRECTION = "Correction";
	private static final String CLUSTER_FILTER_BACKWARD = "Backward";
	private static final String CLUSTER_FILTER_FORWARD = "Forward";
	private static final String CLUSTER_FILTER_DELETE = "Delete";
	private static final String CAPTURE = "Screen capture";
	private static final String DUMP = "Dump raw data";
	private static final String TRAVERSE_ALL = "Traverse all";
	private static final String TRAVERSE_ANNOTATED_ONLY = "Traverse annotated only";
	private static final String TRAVERSE_UNANNOTATED_ONLY = "Traverse unannotated only";
	private static final String[] RADIOBUTTON_TEXTS = new String[] {TRAVERSE_ALL, TRAVERSE_ANNOTATED_ONLY, TRAVERSE_UNANNOTATED_ONLY};
	public static final String MASK_MISSING = "Mask missing values";
	public static final String UNMASK_MISSING = "Unmask missing values";
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String DEFAULT_MESSAGE = "enter new annotation here";
	public static final String[] GENOTYPE_OPTIONS = new String[] {"-","A/A","A/B","B/B"};
	public static final int NUM_MARKERS_TO_SAVE_IN_HISTORY = 10;


	private JButton first, previous, next, last;
	private JTextField navigationField;
	private ScatterPanel scatPanel;
	private JLabel sizeLabel;
	private JLabel gcLabel;
	private JPanel qcPanel;
	private JLabel pcLabel, nStageStDevLabel, correctionRatioLabel;
	private JScrollPane annotationScrollPane;
	private JPanel annotationPanel;
	private JPanel annotationPanelLowerPart;
	private JCheckBox[] annotationCheckBoxes;
	private JRadioButton[] fileterRadioButtons;
	private JTextField addAnnotationField;
	private boolean showAnnotationShortcuts;
	private boolean annotationAutoAdv;
	private boolean showAllMarkersOrNot;
	private boolean showAnnotatedOrUnannotated;
	private boolean isInitilizing;
	private char[] annotationKeys;
	private JComboBox<String> newGenotype;
	private SpringLayout annotationPanelLowerPartLayout;


	private Project proj;
	private String[] masterMarkerList;
	private String[] masterCommentList;
	private String[] markerList;
	private String[] commentList;
	private float[][][][] cents;
	private String[] centList;
	private JCheckBox[] centBoxes;
	private boolean[] displayCents;
	private JLabel[] centLabels;
	private int markerIndex;
	private int markerIndexBak;
	private int[] markerIndexHistory;
	private int previousMarkerIndex;
	private JTextField markerName, commentLabel;
	private String[] samples;
	private int plot_type;
	private byte size;
	private float gcThreshold;
	private long sampleListFingerprint;
	private MarkerLookup markerLookup;
	private boolean jar;
	private SampleData sampleData;
	private JCheckBox symmetryBox;
	private JCheckBox correctionBox;
	private boolean maskMissing;
	private ClusterFilterCollection clusterFilterCollection;
	private AnnotationCollection annotationCollection;
	private byte currentClusterFilter;
	private JTextField clusterFilterNavigation;
	private boolean isClusterFilterUpdated;
	private boolean isAnnotationUpdated;
	private String sessionID;
	private AutoSaveForScatterPlot autoSave;
	private JRadioButton[] typeRadioButtons;
	private MarkerDataLoader markerDataLoader;
    private ColorKeyPanel colorKeyPanel;
	private Color[] colorScheme;
	private int indexOfAnnotationUsedAsMarkerList; 
	private boolean fail;
	private boolean exitOnClose;
	private Logger log;
	private double stdevFilter, correctionRatio;
	private int numComponents ;
	private PrincipalComponentsResiduals pcResids;
	
	public ScatterPlot(Project project, String[] initMarkerList, String[] initCommentList, boolean exitOnClose) {
		SampleList sampleList;
		long time;
		
		time = new Date().getTime();

		proj = project;
		jar = proj.getJarStatus();
		log = proj.getLog();
		size = DEFAULT_SIZE;
		this.exitOnClose = exitOnClose;
		gcThreshold = (float)DEFAULT_GC_THRESHOLD/100f;
		stdevFilter =0;
		correctionRatio =PrincipalComponentsIntensity.DEFAULT_CORRECTION_RATIO;
		markerIndexHistory = Array.intArray(NUM_MARKERS_TO_SAVE_IN_HISTORY, -1);

		if (!Files.exists(proj.getDir(Project.MARKER_DATA_DIRECTORY), proj.getJarStatus())) {
			JOptionPane.showMessageDialog(null, "Directory "+proj.getProperty(Project.MARKER_DATA_DIRECTORY)+" does not exist; the raw data needs to be parsed and transposed before it can be visualized", "Error", JOptionPane.ERROR_MESSAGE);
			fail = true;
			return;
		}

		if (Files.list(proj.getDir(Project.MARKER_DATA_DIRECTORY), "marker", MarkerData.MARKER_DATA_FILE_EXTENSION, false, proj.getJarStatus()).length==0) {
			JOptionPane.showMessageDialog(null, "There is no data in directory "+proj.getProperty(Project.MARKER_DATA_DIRECTORY)+"; the raw data needs to be parsed and transposed before it can be visualized", "Error", JOptionPane.ERROR_MESSAGE);
			fail = true;
			return;
		}

		markerLookup = proj.getMarkerLookup();
		if (markerLookup == null) {
			fail = true;
			return;
		}

		sampleList = proj.getSampleList();
		samples = sampleList.getSamples();
		sampleListFingerprint = sampleList.getFingerprint();
		sampleData = proj.getSampleData(3, true);
		
		fail = sampleData.failedToLoad();
		if (fail) {
			proj.getLog().reportError("Without a SampleData file, ScatterPlot will not start");
			JOptionPane.showMessageDialog(null, "Without a SampleData file, ScatterPlot will not start", "Error", JOptionPane.ERROR_MESSAGE);
			return;
		}
		
		masterMarkerList = initMarkerList;
		masterCommentList = initCommentList;
		if (masterMarkerList == null) {
			loadMarkerListFromFile();
			if (masterMarkerList == null) {
				fail = true;
				return;
			}
			if (masterMarkerList.length == 0) {
				JOptionPane.showMessageDialog(null, "Error - file '"+proj.getFilename(Project.DISPLAY_MARKERS_FILENAME)+"' was devoid of any valid markers", "Error", JOptionPane.ERROR_MESSAGE);
				fail = true;
				return;
			}
		}
		if (masterCommentList == null) {
			masterCommentList = Array.stringArray(masterMarkerList.length, "");
		}
		
		markerList = masterMarkerList;
		commentList = masterCommentList;
		showAllMarkersOrNot = true;
		showAnnotatedOrUnannotated = true;
		showAnnotationShortcuts = true;
		isInitilizing = true;
		indexOfAnnotationUsedAsMarkerList = -1;
		
		loadMarkerDataFromList(0);
		pcResids = loadPcResids();// returns null if not found, marker data should return original x/y if null
		numComponents = 0;//initialize to 0 PCs
		log.reportError("3\t"+ext.getTimeElapsed(time));
		loadCentroids();
		sessionID = (new Date().getTime()+"").substring(5);
		isClusterFilterUpdated = false;
		isAnnotationUpdated = false;
		autoSave = null;
		fail = !loadClusterFilterFiles();
		if (fail) {
			proj.getLog().reportError("Chose to ignore prompt for autosaved cluster filters; ScatterPlot will not start");
			return;
		}
		log.reportError("4\t"+ext.getTimeElapsed(time));
		
		annotationAutoAdv = true;

		scatPanel = new ScatterPanel(this);
		loadAnnotationCollection();
		colorScheme = scatPanel.getColorScheme();
		
		setLayout(new BorderLayout());
		add(scatPanel, BorderLayout.CENTER);
		add(markerPanel(), BorderLayout.NORTH);
		add(eastPanel(), BorderLayout.EAST);
		colorKeyPanel = new ColorKeyPanel(sampleData, scatPanel);
		add(colorKeyPanel, BorderLayout.SOUTH);

		scatPanel.setPointsGeneratable(true);
		scatPanel.setQcPanelUpdatable(true);
		scatPanel.setExtraLayersVisible(new byte[] {99});
		updateGUI();
		displayIndex(navigationField);
//		clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+" of "+clusterFilterCollection.getSize(getMarkerName()));
//		currentClusterFilter=0;
//		setCurrentClusterFilter((byte) (clusterFilterCollection.getSize(getMarkerName())-1));
//		setCurrentClusterFilter();
		if (clusterFilterCollection==null) {
			currentClusterFilter = 0;
		} else {
			currentClusterFilter = (byte)(clusterFilterCollection.getSize(getMarkerName())-1);
		}
//		scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
		displayClusterFilterIndex();
		activateAllAnnotationMaps();

//    	newGenotype.setSelectedIndex(clusterFilterCollection.getGenotype(getMarkerName(), currentClusterFilter)+1);
		symmetryBox.setSelected(true);
		if (centList.length > 0) {
			centBoxes[0].setSelected(true);
		}
		
		inputMapAndActionMap(this);
//		inputMapAndActionMap(annotationPanel);
//		inputMapAndActionMap(annotationScrollPane);

//		next.getInputMap().put(KeyStroke.getKeyStroke("space"), NEXT);
//		next.setActionMap(actionMap);
//		previous.setActionMap(actionMap);
		scatPanel.grabFocus();

		
		setBounds(20, 20, 1000, 720);
		setVisible(true);
	}
	
	public boolean failed() {
		return fail;
	}
	
	private JPanel markerPanel() {
		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new GridLayout(3, 1));
		//markerName = new JLabel("", JLabel.CENTER);
		markerName = new JTextField("");
		markerName.setBorder(null);
		markerName.setEditable(false);
		markerName.setBackground(BACKGROUND_COLOR);
		markerName.setHorizontalAlignment(JTextField.CENTER);
		markerName.setFont(new Font("Arial", 0, 20));
		descrPanel.add(markerName);

		//commentLabel = new JLabel("", JLabel.CENTER);
		commentLabel = new JTextField("");
		commentLabel.setBorder(null);
		commentLabel.setEditable(false);
		commentLabel.setBackground(BACKGROUND_COLOR);
		commentLabel.setHorizontalAlignment(JTextField.CENTER);
		commentLabel.setFont(new Font("Arial", 0, 14));
		descrPanel.add(commentLabel);

		JPanel navigationPanel = new JPanel();
		first = new JButton(Grafik.getImageIcon("images/firstLast/First.gif", true));
		first.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif", true));
		first.addActionListener(this);
		first.setActionCommand(FIRST);
		first.setPreferredSize(new Dimension(20, 20));
		previous = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif", true));
		previous.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		previous.addActionListener(this);
		previous.setActionCommand(PREVIOUS);
		previous.setPreferredSize(new Dimension(20, 20));
		navigationField = new JTextField("", 8);
		navigationField.setHorizontalAlignment(JTextField.CENTER);
		navigationField.setFont(new Font("Arial", 0, 14));
		//navigationField.setEditable(false);
		navigationField.setBackground(BACKGROUND_COLOR);
		navigationField.addActionListener(new ActionListener() {
	        public void actionPerformed(ActionEvent e) {
				try {
					int trav = Integer.valueOf(((JTextField)e.getSource()).getText().split("[\\s]+")[0]).intValue()-1;
					if (trav >=0 && trav < markerList.length) {
						markerIndex = trav;
						updateMarkerIndexHistory();
					}
				} catch (NumberFormatException nfe) {}
				displayIndex((JTextField)e.getSource());
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);
				setCurrentClusterFilter();
				updateGUI();
				displayClusterFilterIndex();
	        }
		});
		
		
		next = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif", true));
		next.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		next.addActionListener(this);
		next.setActionCommand(NEXT);
		next.setPreferredSize(new Dimension(20, 20));
		last = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif", true));
		last.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif", true));
		last.addActionListener(this);
		last.setActionCommand(LAST);
		last.setPreferredSize(new Dimension(20, 20));
		navigationPanel.add(first);
		navigationPanel.add(previous);
		navigationPanel.add(navigationField);
		navigationPanel.add(next);
		navigationPanel.add(last);

		navigationPanel.setBackground(BACKGROUND_COLOR);
		descrPanel.add(navigationPanel);
		descrPanel.setBackground(BACKGROUND_COLOR);
		return descrPanel;
    }

	private JPanel eastPanel() {
		JPanel eastPanel;
		JTabbedPane tabbedPane;

		tabbedPane = new JTabbedPane();
		tabbedPane.setBackground(BACKGROUND_COLOR);

		tabbedPane.addTab("Control", null, controlPanel(), "Manipulate the chart");
		tabbedPane.setMnemonicAt(0, KeyEvent.VK_1);

		tabbedPane.addTab("Centroid", null, centroidPanel(), "Displays the centroid");
		tabbedPane.setMnemonicAt(1, KeyEvent.VK_2);

//		annotationPanel = new JPanel();
//		annotationPanel.add(createAnnotationPanel());
//		annotationPanelComment = new JPanel();
//		annotationPanel.add(annotationPanelComment);
		annotationPanel();
//		annotationScrollPane = new JScrollPane(annotationPanel);
//		annotationPanel.addKeyListener(this);
//		annotationPanel.setFocusable(true);
//		annotationPanel.requestFocus();
		tabbedPane.addTab("Annotation", null, annotationScrollPane, "Annotation to this marker");
		tabbedPane.setMnemonicAt(2, KeyEvent.VK_3);

		qcPanel = new JPanel();
		qcPanel.setLayout(new GridLayout(8, 1));
		qcPanel.setBackground(BACKGROUND_COLOR);

		eastPanel = new JPanel();
		eastPanel.setLayout(new BoxLayout(eastPanel, BoxLayout.PAGE_AXIS));
		eastPanel.setBackground(BACKGROUND_COLOR);
		eastPanel.add(tabbedPane);
		eastPanel.add(clusterFilterPanel());
//		tabbedPanel.add(plotTypePanel());
		eastPanel.add(qcPanel);

		return eastPanel;
    }

	private JPanel plotTypePanel() {
		JPanel plotTypePanel = new JPanel();
//		plotTypePanel.setLayout(new BoxLayout(plotTypePanel, BoxLayout.Y_AXIS));
		plotTypePanel.setLayout(new GridLayout(4, 1));
		plotTypePanel.setBackground(BACKGROUND_COLOR);

		ItemListener typeListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				JRadioButton jrb = (JRadioButton)ie.getItem();
				if (jrb.isSelected()) {
					for (int i = 0; i<MarkerData.TYPES.length; i++) {
						if (jrb.getText().equals(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1])) {
							plot_type = i;
							scatPanel.setPointsGeneratable(true);
							scatPanel.setQcPanelUpdatable(true);
							updateGUI();
						}
					}
				}
			}
		};
		// --- Beginning of the original block ---
		ButtonGroup typeRadio = new ButtonGroup();
		typeRadioButtons = new JRadioButton[MarkerData.TYPES.length];
		for (int i = 0; i<MarkerData.TYPES.length; i++) {
			typeRadioButtons[i] = new JRadioButton(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1], false);
			typeRadioButtons[i].setFont(new Font("Arial", 0, 14));
			typeRadio.add(typeRadioButtons[i]);
			typeRadioButtons[i].addItemListener(typeListener);
			typeRadioButtons[i].setBackground(BACKGROUND_COLOR);
//			tabPanel.add(typeRadioButtons[i], gbc);
			plotTypePanel.add(typeRadioButtons[i]);
		}
		typeRadioButtons[1].setSelected(true);

		return plotTypePanel;
	}

	private JPanel sizeSliderPanel() {
		JPanel sizeSliderPanel = new JPanel();
		sizeSliderPanel.setLayout(new BoxLayout(sizeSliderPanel, BoxLayout.Y_AXIS));
		sizeSliderPanel.setBackground(BACKGROUND_COLOR);

		JSlider slider = new JSlider(JSlider.HORIZONTAL, 2, 20, DEFAULT_SIZE);
//		slider.setSize(new Dimension(150, 20));
		slider.setBackground(BACKGROUND_COLOR);
		sizeLabel = new JLabel("Size = "+DEFAULT_SIZE, JLabel.CENTER);
		sizeLabel.setFont(new Font("Arial", Font.PLAIN, 16));
//		tabPanel.add(sizeLabel, gbc);
		sizeSliderPanel.add(sizeLabel);

		slider.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent ce) {
				JSlider slider = (JSlider)ce.getSource();
				sizeLabel.setText("Size = "+slider.getValue());
				size = (byte)slider.getValue();
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(false);
				scatPanel.paintAgain();
			}
		});
//		tabPanel.add(slider, gbc);
		sizeSliderPanel.add(slider);

		return sizeSliderPanel;
	}

	private JPanel gcSliderPanel() {
		JPanel gcSliderPanel = new JPanel();
		gcSliderPanel.setLayout(new BoxLayout(gcSliderPanel, BoxLayout.Y_AXIS));
		gcSliderPanel.setBackground(BACKGROUND_COLOR);

		JSlider slider = new JSlider(JSlider.HORIZONTAL, 2, 20, DEFAULT_SIZE);
//		slider.setSize(new Dimension(150, 20));
		slider.setBackground(BACKGROUND_COLOR);
		slider = new JSlider(JSlider.HORIZONTAL, 0, 100, DEFAULT_SIZE);
		slider.setBackground(BACKGROUND_COLOR);
		gcLabel = new JLabel("GC > 0."+DEFAULT_GC_THRESHOLD, JLabel.CENTER);
		gcLabel.setFont(new Font("Arial", Font.PLAIN, 16));
//		tabPanel.add(gcLabel, gbc);
		gcSliderPanel.add(gcLabel);

		slider.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent ce) {
				JSlider slider = (JSlider)ce.getSource();
				gcThreshold = (float)slider.getValue()/100f;
				gcLabel.setText("GC > "+ext.formDeci(gcThreshold, 2, true));
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);
				scatPanel.paintAgain();
				//qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");
			}
		});
//		tabPanel.add(slider, gbc);
		gcSliderPanel.add(slider);

		return gcSliderPanel;
	}
	
	private JPanel pcSliderPanel() {
		JPanel pcSliderPanel = new JPanel();
		pcSliderPanel.setLayout(new BoxLayout(pcSliderPanel, BoxLayout.Y_AXIS));
		pcSliderPanel.setBackground(BACKGROUND_COLOR);

		JSlider slider = new JSlider(JSlider.HORIZONTAL, 2, 20, DEFAULT_SIZE);
		// slider.setSize(new Dimension(150, 20));
		slider.setBackground(BACKGROUND_COLOR);
		slider = new JSlider(JSlider.HORIZONTAL, 0,  Integer.parseInt(proj.getProperty(Project.INTENSITY_PC_NUM_COMPONENTS)), 0);
		slider.setValue(0);
		slider.setBackground(BACKGROUND_COLOR);
		if (pcResids == null) {
			pcLabel = new JLabel("PC file not detected", JLabel.CENTER);
		} else {
			pcLabel = new JLabel("Total PCs Loaded = " + pcResids.getNumComponents(), JLabel.CENTER);
		}
		pcLabel.setFont(new Font("Arial", Font.PLAIN, 16));
		pcSliderPanel.add(pcLabel);

		slider.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent ce) {
				JSlider slider = (JSlider) ce.getSource();
				numComponents = slider.getValue();
				pcLabel.setText("NumPCs = "+numComponents);
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);
				scatPanel.paintAgain();
			}
		});
		pcSliderPanel.add(slider);

		return pcSliderPanel;
	}

	private JPanel nstageStdevSliderPanel() {
		JPanel pcSliderPanel = new JPanel();
		pcSliderPanel.setLayout(new BoxLayout(pcSliderPanel, BoxLayout.Y_AXIS));
		pcSliderPanel.setBackground(BACKGROUND_COLOR);

		JSlider slider = new JSlider(JSlider.HORIZONTAL, 0, 100, DEFAULT_SIZE);
		// slider.setSize(new Dimension(150, 20));
		slider.setBackground(BACKGROUND_COLOR);
		slider = new JSlider(JSlider.HORIZONTAL, 0,  100, 0);
		slider.setValue(0);
		slider.setBackground(BACKGROUND_COLOR);
		nStageStDevLabel = new JLabel("Standard deviation filter "+stdevFilter, JLabel.CENTER);
		nStageStDevLabel.setFont(new Font("Arial", Font.PLAIN, 16));
		// tabPanel.add(gcLabel, gbc);
		pcSliderPanel.add(nStageStDevLabel);

		slider.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent ce) {
				JSlider slider = (JSlider) ce.getSource();
				stdevFilter = (double)slider.getValue()/25;
				nStageStDevLabel.setText("Alternate Genotype stdev Filter "+stdevFilter);
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);
				scatPanel.paintAgain();
				// qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");
			}
		});
		// tabPanel.add(slider, gbc);
		pcSliderPanel.add(slider);

		return pcSliderPanel;
}

	private JPanel nstageCorrectionRatio() {
		JPanel pcSliderPanel = new JPanel();
		pcSliderPanel.setLayout(new BoxLayout(pcSliderPanel, BoxLayout.Y_AXIS));
		pcSliderPanel.setBackground(BACKGROUND_COLOR);

		JSlider slider = new JSlider(JSlider.HORIZONTAL, 0, 100, DEFAULT_SIZE);
		// slider.setSize(new Dimension(150, 20));
		slider.setBackground(BACKGROUND_COLOR);
		slider = new JSlider(JSlider.HORIZONTAL, 0, 100, 0);
		slider.setValue(0);
		slider.setBackground(BACKGROUND_COLOR);
		correctionRatioLabel = new JLabel("Correction Ratio " + correctionRatio, JLabel.CENTER);
		correctionRatioLabel.setFont(new Font("Arial", Font.PLAIN, 16));
		String usage = "This limits the ratio of PCs to the number of samples in a cluster...";

		correctionRatioLabel.setToolTipText(usage);
		// tabPanel.add(gcLabel, gbc);
		pcSliderPanel.add(correctionRatioLabel);

		slider.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent ce) {
				JSlider slider = (JSlider) ce.getSource();
				correctionRatio = (double) slider.getValue() / 100;
				correctionRatioLabel.setText("Correction Ratio " + correctionRatio);
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);
				scatPanel.paintAgain();
				// qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");
			}
		});
		// tabPanel.add(slider, gbc);
		pcSliderPanel.add(slider);

		return pcSliderPanel;
	}
	private JPanel clusterFilterPanel() {
		JPanel clusterFilterPanel = new JPanel();
		clusterFilterPanel.setBackground(BACKGROUND_COLOR);
//		clusterFilterPanel.setLayout(new BoxLayout(clusterFilterPanel, BoxLayout.PAGE_AXIS));
		clusterFilterPanel.setSize(50, 10);
		
//		JButton first1 = new JButton(Grafik.getImageIcon("images/firstLast/First.gif", true));
//		first1.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif", true));
//		first1.addActionListener(this);
//		first1.setActionCommand(FIRST);
//		first1.setPreferredSize(new Dimension(20, 20));
		JButton backward = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif", true));
		backward.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		backward.addActionListener(this);
		backward.setActionCommand(CLUSTER_FILTER_BACKWARD);
		backward.setPreferredSize(new Dimension(20, 20));
		clusterFilterNavigation = new JTextField("", 5);
		clusterFilterNavigation.setHorizontalAlignment(JTextField.CENTER);
		clusterFilterNavigation.setFont(new Font("Arial", 0, 14));
		//navigationField.setEditable(false);
		clusterFilterNavigation.setBackground(BACKGROUND_COLOR);
		clusterFilterNavigation.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent focusevent) {}

			public void focusLost(FocusEvent fe) {
				try {
					int trav = Integer.valueOf(((JTextField)fe.getSource()).getText().split("[\\s]+")[0]).intValue()-1;
					if (trav >= 0 && trav < clusterFilterCollection.getSize(getMarkerName())) {
//						currentClusterFilter = (byte) trav;
						setCurrentClusterFilter((byte) trav);
					}
				} catch (NumberFormatException nfe) {}
//				clusterFilterNavigation.setText((currentClusterFilter+1)+" of "+clusterFilterCollection.getSize(getMarkerName()));
//				displayClusterFilterIndex();
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);
//				scatPanel.generateRectangles();
				updateGUI();
				displayClusterFilterIndex();
			}
		});

		JButton forward = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif", true));
		forward.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		forward.addActionListener(this);
		forward.setActionCommand(CLUSTER_FILTER_FORWARD);
		forward.setPreferredSize(new Dimension(20, 20));
//		JButton last1 = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif", true));
//		last1.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif", true));
//		last1.addActionListener(this);
//		last1.setActionCommand(LAST);
//		last1.setPreferredSize(new Dimension(20, 20));
		JButton delete = new JButton(Grafik.getImageIcon("images/delete9.png", true));
		delete.setDisabledIcon(Grafik.getImageIcon("images/delete8.png", true));
		delete.addActionListener(this);
		delete.setActionCommand(CLUSTER_FILTER_DELETE);
		delete.setPreferredSize(new Dimension(20, 20));
//		clusterFilterPanel.add(first1);
		clusterFilterPanel.add(backward);
		clusterFilterPanel.add(clusterFilterNavigation);
		clusterFilterPanel.add(forward);
//		clusterFilterPanel.add(last1);
		clusterFilterPanel.add(delete);

    	newGenotype = new JComboBox<String>(GENOTYPE_OPTIONS);
		ActionListener newGenotypeListener = new ActionListener() {
			@SuppressWarnings("unchecked")
			public void actionPerformed(ActionEvent e) {
				byte newGenotypeSelected;
				String newGenoSelected;

				newGenoSelected = (String) ((JComboBox<String>) e.getSource()).getSelectedItem();
				newGenotypeSelected = -2;
				for (int i = 0; i < GENOTYPE_OPTIONS.length; i++) {
					if (newGenoSelected.equals(GENOTYPE_OPTIONS[i])) {
						newGenotypeSelected = (byte) (i - 1);
						break;
					}
				}
				if ((clusterFilterCollection.getGenotype(getMarkerName(), currentClusterFilter) != newGenotypeSelected)) {
					updateCurrentClusterFilterGenotype(newGenotypeSelected, false);
				}
			}
		};
    	newGenotype.addActionListener(newGenotypeListener);
//		newGenotype.setActionCommand(NEW_GENOTYPE);//???
    	clusterFilterPanel.add(newGenotype);
//		typePanel.add(newGenotype, gbc);

		clusterFilterPanel.setBackground(BACKGROUND_COLOR);

		return clusterFilterPanel;
	}

	private JPanel controlPanel() {
		JPanel controlPanel;

		controlPanel = new JPanel();
		controlPanel.setLayout(new GridBagLayout());
		controlPanel.setSize(50, 100);
		controlPanel.setBackground(BACKGROUND_COLOR);

		GridBagConstraints gbc = new GridBagConstraints();
		gbc.insets = new Insets(1,3,0,30);
		gbc.weightx = 1.0;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.gridwidth = GridBagConstraints.REMAINDER;

		controlPanel.add(sizeSliderPanel(), gbc);
		controlPanel.add(gcSliderPanel(), gbc);
		controlPanel.add(pcSliderPanel(), gbc);
		controlPanel.add(nstageStdevSliderPanel(), gbc);
		controlPanel.add(nstageCorrectionRatio(), gbc);
		
		ItemListener symmetryListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
//				scatPanel.setPointsGenerated(true); ??? Why not true?
//				scatPanel.setUpdateQcPanel(false); ??? Why cannot set to false?
				updateGUI();
			}
		};
		
		symmetryBox = new JCheckBox("Symmetric axes");
		symmetryBox.setFont(new Font("Arial", 0, 14));
		symmetryBox.addItemListener(symmetryListener);
		symmetryBox.setBackground(BACKGROUND_COLOR);

		
		ItemListener correctionListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);
				scatPanel.paintAgain();
				updateGUI();
			}
		};
		correctionBox = new JCheckBox("Correct Data");
		correctionBox.setToolTipText("The keyboard shortcut for this function is ctrl+D");
		correctionBox.setFont(new Font("Arial", 0, 14));
		correctionBox.addItemListener(correctionListener);
		correctionBox.setBackground(BACKGROUND_COLOR);
		
		
		controlPanel.add(symmetryBox, gbc);
		controlPanel.add(correctionBox, gbc);

		JButton button = new JButton(CAPTURE);
		button.addActionListener(this);
		button.setActionCommand(CAPTURE);
		controlPanel.add(button, gbc);

		button = new JButton(DUMP);
		button.addActionListener(this);
		button.setActionCommand(DUMP);
		controlPanel.add(button, gbc);
		
		button = new JButton(MASK_MISSING);
		button.addActionListener(this);
		button.setActionCommand(MASK_MISSING);
		controlPanel.add(button, gbc);
		
		controlPanel.add(plotTypePanel(), gbc);

		return controlPanel;
	}

//	private JPanel controlPanel() {
//		JPanel controlPanel;
//		JPanel panel;
//		JButton button;
//
//		controlPanel = new JPanel();
//		controlPanel.setSize(100, 500);
//		controlPanel.setBackground(BACKGROUND_COLOR);
//
////		GridBagConstraints gbc = new GridBagConstraints();
////		gbc.insets = new Insets(1,3,0,30);
////		gbc.weightx = 1.0;
////		gbc.fill = GridBagConstraints.HORIZONTAL;
////		gbc.gridwidth = GridBagConstraints.REMAINDER;
////		controlPanel.setLayout(new BoxLayout(controlPanel, BoxLayout.Y_AXIS));
//		SpringLayout springLayout = new SpringLayout();
//		controlPanel.setLayout(springLayout);
//
//		panel = sizeSliderPanel();
//		controlPanel.add(panel);
//		springLayout.putConstraint(SpringLayout.WEST, panel, 5, SpringLayout.WEST, controlPanel);
//		springLayout.putConstraint(SpringLayout.NORTH, panel, 5, SpringLayout.NORTH, controlPanel);
//		panel = gcSliderPanel();
//		controlPanel.add(panel);
//		springLayout.putConstraint(SpringLayout.WEST, panel, 5, SpringLayout.WEST, controlPanel);
//		springLayout.putConstraint(SpringLayout.NORTH, panel, 45, SpringLayout.NORTH, controlPanel);
//
//		ItemListener symmetryListener = new ItemListener() {
//			public void itemStateChanged(ItemEvent ie) {
////				scatPanel.setPointsGenerated(true); ??? Why not true?
////				scatPanel.setUpdateQcPanel(false); ??? Why cannot set to false?
//				updateGUI();
//			}
//		};
//		
//		symmetryBox = new JCheckBox("Symmetric axes");
//		symmetryBox.setFont(new Font("Arial", 0, 14));
//		symmetryBox.addItemListener(symmetryListener);
//		symmetryBox.setBackground(BACKGROUND_COLOR);
//
//		controlPanel.add(symmetryBox);
//		springLayout.putConstraint(SpringLayout.WEST, symmetryBox, 5, SpringLayout.WEST, controlPanel);
//		springLayout.putConstraint(SpringLayout.NORTH, symmetryBox, 85, SpringLayout.NORTH, controlPanel);
////		tabPanel.add(symmetryBox);
//		
//		button = new JButton(CAPTURE);
//		button.addActionListener(this);
//		button.setActionCommand(CAPTURE);
//		controlPanel.add(button);
//		springLayout.putConstraint(SpringLayout.WEST, button, 5, SpringLayout.WEST, controlPanel);
//		springLayout.putConstraint(SpringLayout.NORTH, button, 115, SpringLayout.NORTH, controlPanel);
////		tabPanel.add(button);
//
//		button = new JButton(DUMP);
//		button.addActionListener(this);
//		button.setActionCommand(DUMP);
//		controlPanel.add(button);
//		springLayout.putConstraint(SpringLayout.WEST, button, 5, SpringLayout.WEST, controlPanel);
//		springLayout.putConstraint(SpringLayout.NORTH, button, 145, SpringLayout.NORTH, controlPanel);
////		tabPanel.add(button);
//		
//		button = new JButton(MASK_MISSING);
//		button.addActionListener(this);
//		button.setActionCommand(MASK_MISSING);
//		controlPanel.add(button);
//		springLayout.putConstraint(SpringLayout.WEST, button, 5, SpringLayout.WEST, controlPanel);
//		springLayout.putConstraint(SpringLayout.NORTH, button, 175, SpringLayout.NORTH, controlPanel);
//		
//		panel = plotTypePanel();
//		controlPanel.add(panel);
//		springLayout.putConstraint(SpringLayout.WEST, panel, 5, SpringLayout.WEST, controlPanel);
//		springLayout.putConstraint(SpringLayout.NORTH, panel, 205, SpringLayout.NORTH, controlPanel);
//		
//
//		return controlPanel;
//	}


//	private JPanel controlPanel() {
//		JPanel controlPanel;
//		JPanel panel;
//		JButton button;
//
//		controlPanel = new JPanel();
//		controlPanel.setSize(100, 500);
//		controlPanel.setBackground(BACKGROUND_COLOR);
//
//		controlPanel.setLayout(new BoxLayout(controlPanel, BoxLayout.Y_AXIS));
////		SpringLayout springLayout = new SpringLayout();
////		controlPanel.setLayout(springLayout);
//
//		panel = sizeSliderPanel();
//		controlPanel.add(panel);
////		springLayout.putConstraint(SpringLayout.WEST, panel, 5, SpringLayout.WEST, controlPanel);
////		springLayout.putConstraint(SpringLayout.NORTH, panel, 5, SpringLayout.NORTH, controlPanel);
//		panel = gcSliderPanel();
//		controlPanel.add(panel);
////		springLayout.putConstraint(SpringLayout.WEST, panel, 5, SpringLayout.WEST, controlPanel);
////		springLayout.putConstraint(SpringLayout.NORTH, panel, 45, SpringLayout.NORTH, controlPanel);
//
//		ItemListener symmetryListener = new ItemListener() {
//			public void itemStateChanged(ItemEvent ie) {
////				scatPanel.setPointsGenerated(true); ??? Why not true?
////				scatPanel.setUpdateQcPanel(false); ??? Why cannot set to false?
//				updateGUI();
//			}
//		};
//		
//		symmetryBox = new JCheckBox("Symmetric axes");
//		symmetryBox.setFont(new Font("Arial", 0, 14));
//		symmetryBox.addItemListener(symmetryListener);
//		symmetryBox.setBackground(BACKGROUND_COLOR);
//
//		controlPanel.add(symmetryBox);
////		springLayout.putConstraint(SpringLayout.WEST, symmetryBox, 5, SpringLayout.WEST, controlPanel);
////		springLayout.putConstraint(SpringLayout.NORTH, symmetryBox, 85, SpringLayout.NORTH, controlPanel);
////		tabPanel.add(symmetryBox);
//		
//		button = new JButton(CAPTURE);
//		button.addActionListener(this);
//		button.setActionCommand(CAPTURE);
//		controlPanel.add(button);
////		springLayout.putConstraint(SpringLayout.WEST, button, 5, SpringLayout.WEST, controlPanel);
////		springLayout.putConstraint(SpringLayout.NORTH, button, 115, SpringLayout.NORTH, controlPanel);
////		tabPanel.add(button);
//
//		button = new JButton(DUMP);
//		button.addActionListener(this);
//		button.setActionCommand(DUMP);
//		controlPanel.add(button);
////		springLayout.putConstraint(SpringLayout.WEST, button, 5, SpringLayout.WEST, controlPanel);
////		springLayout.putConstraint(SpringLayout.NORTH, button, 145, SpringLayout.NORTH, controlPanel);
////		tabPanel.add(button);
//		
//		button = new JButton(MASK_MISSING);
//		button.addActionListener(this);
//		button.setActionCommand(MASK_MISSING);
//		controlPanel.add(button);
////		springLayout.putConstraint(SpringLayout.WEST, button, 5, SpringLayout.WEST, controlPanel);
////		springLayout.putConstraint(SpringLayout.NORTH, button, 175, SpringLayout.NORTH, controlPanel);
//		
//		panel = plotTypePanel();
//		controlPanel.add(panel);
////		springLayout.putConstraint(SpringLayout.WEST, panel, 5, SpringLayout.WEST, controlPanel);
////		springLayout.putConstraint(SpringLayout.NORTH, panel, 205, SpringLayout.NORTH, controlPanel);
//		
//
//		return controlPanel;
//	}


	private JPanel centroidPanel() {
		JPanel centroidPanel;

		centroidPanel = new JPanel();
//		tabPanel.setLayout(new BoxLayout(tabPanel, BoxLayout.Y_AXIS));
		centroidPanel.setLayout(new GridBagLayout());
		centroidPanel.setSize(50, 100);
		centroidPanel.setBackground(BACKGROUND_COLOR);

		GridBagConstraints gbc = new GridBagConstraints();
		gbc.insets = new Insets(1,3,0,30);
		gbc.weightx = 1.0;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.gridwidth = GridBagConstraints.REMAINDER;

		ItemListener centListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				int index= ext.indexOfStr(((JCheckBox)ie.getSource()).getText(), centList);
				displayCents[index] = ((JCheckBox)ie.getSource()).isSelected();
				centLabels[index].setVisible(displayCents[index]);
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(false);
				updateGUI();
			}
		};
		JLabel label = new JLabel("  ");
		label.setFont(new Font("Arial", 0, 20));
		centroidPanel.add(label, gbc);
//		tabPanel.add(label);
		label = new JLabel("Centroids:");
		label.setFont(new Font("Arial", 0, 20));
		label.setHorizontalAlignment(JLabel.CENTER);
		centroidPanel.add(label, gbc);
//		tabPanel.add(label);
		centBoxes = new JCheckBox[centList.length];
		displayCents = new boolean[centList.length];
		centLabels = new JLabel[centList.length];
		//log.report(centList.length+"\n");
		for (int i = 0; i<centList.length; i++) {
			centBoxes[i] = new JCheckBox(centList[i]);
			centBoxes[i].setFont(new Font("Arial", 0, 14));
			centBoxes[i].setSelected(displayCents[i]);
			centBoxes[i].addItemListener(centListener);
//			centBoxes[i].setBorder(BorderFactory.createLineBorder(ColorKeyPanel.DEFAULT_COLORS[5+i], 5));
			centBoxes[i].setBorder(BorderFactory.createLineBorder(colorScheme[5+i], 5));
			centBoxes[i].setBorderPainted(true);
			centBoxes[i].setBackground(BACKGROUND_COLOR);
			centroidPanel.add(centBoxes[i], gbc);
//			tabPanel.add(centBoxes[i]);
			
			centLabels[i] = new JLabel("LRR correlation not performed");
			centLabels[i].setVisible(displayCents[i]);
			centroidPanel.add(centLabels[i], gbc);
//			tabPanel.add(centLabels[i]);
		}

		return centroidPanel;
	}

	private void annotationPanel() {
		annotationPanel = new JPanel();
//		annotationScrollPane = new JScrollPane(annotationPanel, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
//		annotationScrollPane.setVisible(true);
//		annotationScrollPane.createVerticalScrollBar();

		annotationPanel.setBackground(Color.WHITE);
		annotationPanel.setLayout(new BoxLayout(annotationPanel, BoxLayout.Y_AXIS));
//		annotationPanelLayout = new SpringLayout();
//		annotationPanel.setLayout(annotationPanelLayout);
		annotationPanel.add(annotationPanelUpper());
//		annotationPanelLayout.putConstraint(SpringLayout.WEST, temp1, 5, SpringLayout.WEST, annotationPanel);
//		annotationPanelLayout.putConstraint(SpringLayout.NORTH, temp1, 5, SpringLayout.NORTH, annotationPanel);
		annotationPanelLowerPart = new JPanel();
		annotationPanelLowerPartLayout = new SpringLayout();
		annotationPanelLowerPart.setLayout(annotationPanelLowerPartLayout);
		annotationPanelLowerPart.setBackground(BACKGROUND_COLOR);
		annotationPanel.add(annotationPanelLower());
//		annotationPanelLayout.putConstraint(SpringLayout.WEST, annotationPanelLowerPart, 5, SpringLayout.WEST, annotationPanel);
//		annotationPanelLayout.putConstraint(SpringLayout.NORTH, annotationPanelLowerPart, 25, SpringLayout.NORTH, annotationPanel);
		annotationScrollPane = new JScrollPane(annotationPanel);
		annotationScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		annotationScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
	}

	private JPanel annotationPanelUpper() {
		JPanel annotationPanelUpper;
		JCheckBox checkBox;
		JButton button;

		annotationPanelUpper = new JPanel();
		annotationPanelUpper.setLayout(new BoxLayout(annotationPanelUpper, BoxLayout.X_AXIS));
		annotationPanelUpper.setBackground(BACKGROUND_COLOR);

		checkBox = new JCheckBox("Shortcut");
		checkBox.setBackground(Color.WHITE);
		checkBox.setHorizontalTextPosition(SwingConstants.LEFT);
		checkBox.setMnemonic(KeyEvent.VK_S);
		checkBox.setSelected(true);
		checkBox.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				JCheckBox checkBox = (JCheckBox) e.getSource();
				if (checkBox.isSelected()) {
					showAnnotationShortcuts = true;
					activateAllAnnotationMaps();
				} else {
					showAnnotationShortcuts = false;
					deactivateAllAnnotationMaps();
				}
			}});
		annotationPanelUpper.add(checkBox);

		checkBox = new JCheckBox("Auto Advance");
		checkBox.setBackground(Color.WHITE);
		checkBox.setHorizontalTextPosition(SwingConstants.LEFT);
		checkBox.setMnemonic(KeyEvent.VK_A);
		checkBox.setSelected(true);
		checkBox.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				JCheckBox checkBox = (JCheckBox) e.getSource();
				if (checkBox.isSelected()) {
					annotationAutoAdv = true;
				} else {
					annotationAutoAdv = false;
				}
			}});
		annotationPanelUpper.add(checkBox);

//		button = new JButton("Export");
		button = new JButton(Grafik.getImageIcon("images/export-icon.png", true));
		button.setToolTipText("Export");
		button.setBackground(Color.WHITE);
		button.setBorder(null);
		button.setSize(10, 11);
		button.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if (e.getActionCommand() != null) {
					annotationCollection.dumpLists(proj);
				}
			}
		});
		annotationPanelUpper.add(button);

//		button = new JButton("Import");
		button = new JButton(Grafik.getImageIcon("images/import-icon.png", true));
		button.setToolTipText("Import");
		button.setBackground(Color.WHITE);
		button.setBorder(null);
		button.setSize(10, 11);
		button.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if (e.getActionCommand() != null) {
					String filename, filenameBak;
					JFileChooser fileChooser;
					int returnVal;
					String[] options;
					int choice;
					
					fileChooser = new JFileChooser(new File(proj.getProjectDir()));
					returnVal = fileChooser.showOpenDialog(null);

			        if (returnVal == JFileChooser.APPROVE_OPTION) {
			        	filename = fileChooser.getSelectedFile().getAbsolutePath();
			    		if (isAnnotationUpdated) {
			    			options = new String[] {"Yes, overwrite", "No"};
			    			choice = JOptionPane.showOptionDialog(null, "New Annotations have been generated. Do you want to save them to the permanent file?", "Overwrite permanent file?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
			    			if (choice == 0) {
			    				filenameBak = proj.getFilename(proj.getProjectDir(), false, false) + "annotations_bak_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".ser";
			    				log.report("Writing to " + filenameBak);
			    				annotationCollection.serialize(filenameBak);
			    			}
			    		}
			    		deactivateAllAnnotationMaps();
			    		if (annotationCollection == null) {
			    			annotationCollection = AnnotationCollection.loadFromLists(filename, null);
							startAutoSaveToTempFile();
			    		} else {
			    			AnnotationCollection.appendFromLists(filename, annotationCollection, proj.getMarkerNames(), null);
			    		}
			    		setAnnotationUpdated(true);
						annotationKeys = annotationCollection.getKeys();
			    		activateAllAnnotationMaps();
						annotationPanelLower();
			        }
				}
			}
		});
		annotationPanelUpper.add(button);

//		annotationPanelControl.setVisible(true);
		return annotationPanelUpper;
	}

	public JPanel annotationPanelLower() {
		JButton removeButton;
		ButtonGroup radioButtonGroup;
		int horizontalMargin;
		int componentHeight;
		int currentHorizontalPos;
		int maxWidth;
		
		annotationPanelLowerPart.removeAll();
		annotationPanelLowerPart.repaint();

		horizontalMargin = 2;
		componentHeight = 20;
		currentHorizontalPos = 0;
		radioButtonGroup = new ButtonGroup();
		fileterRadioButtons = new JRadioButton[RADIOBUTTON_TEXTS.length];
		for (int i = 0; i < RADIOBUTTON_TEXTS.length; i++) {
			fileterRadioButtons[i] = new JRadioButton(RADIOBUTTON_TEXTS[i] + " (n=placeholder)");
			fileterRadioButtons[i].setBackground(Color.WHITE);
			fileterRadioButtons[i].addActionListener(new ActionListener(){
				@Override
				public void actionPerformed(ActionEvent e) {
					String commandText;
					commandText = e.getActionCommand();
					if (commandText.startsWith(TRAVERSE_ALL + " (n=")) {
						showAllMarkersOrNot = true;
					} else if (commandText.startsWith(TRAVERSE_ANNOTATED_ONLY + " (n=")) {
						showAllMarkersOrNot = false;
						showAnnotatedOrUnannotated = true;
					} else {
						showAllMarkersOrNot = false;
						showAnnotatedOrUnannotated = false;
					}
				}
			});
			radioButtonGroup.add(fileterRadioButtons[i]);
			annotationPanelLowerPart.add(fileterRadioButtons[i]);
			annotationPanelLowerPartLayout.putConstraint(SpringLayout.WEST, fileterRadioButtons[i], 5, SpringLayout.WEST, annotationPanelLowerPart);
			annotationPanelLowerPartLayout.putConstraint(SpringLayout.NORTH, fileterRadioButtons[i], currentHorizontalPos, SpringLayout.NORTH, annotationPanelLowerPart);
//			horizontalPos += (horizontalMargin + fileterRadioButtons[i].getSize().height);
			currentHorizontalPos += (horizontalMargin + componentHeight);
		}
		fileterRadioButtons[0].setSelected(true);
		
		MouseListener mouseListenerForAnnotationCheckBoxes = new MouseListener() {
			public void mouseReleased(MouseEvent e) {}
			public void mousePressed(MouseEvent e) {}
			public void mouseExited(MouseEvent e) {}
			public void mouseEntered(MouseEvent e) {}

//			public void mouseClicked(MouseEvent e) {
//				JCheckBox source;
//				String annotation;
//				int annotationIndex;
//
//				if (e.getButton()==MouseEvent.BUTTON3) {
//					annotationIndex = -1;
//					source = (JCheckBox)e.getSource();
//					annotation = source.getText();
//					for (int i = 0; i < annotationKeys.length; i ++) {
//						if (annotationCollection.getDescriptionForComment(annotationKeys[i], showAnnotationShortcuts, true).equals(annotation)) {
//							annotationIndex = i;
//							break;
//						}
//					}
//					
//					if (indexOfAnnotationUsedAsMarkerList == -1) {
//						for (int i = 0; i < annotationIndex; i ++) {
//							annotationCheckBoxes[i].setEnabled(false);
//						}
//						for (int i = annotationIndex + 1; i < annotationCheckBoxes.length; i ++) {
//							annotationCheckBoxes[i].setEnabled(false);
//						}
//						indexOfAnnotationUsedAsMarkerList = annotationIndex;
//					} else if (indexOfAnnotationUsedAsMarkerList == annotationIndex) {
//						for (int i = 0; i < annotationCheckBoxes.length; i ++) {
//							annotationCheckBoxes[i].setEnabled(true);
//						}
//						indexOfAnnotationUsedAsMarkerList = -1;
//					}
//				}

			@Override
			public void mouseClicked(MouseEvent e) {
				JCheckBox source;
				JPopupMenu menu;
				String annotation;
				int annotationIndex;

				source = (JCheckBox)e.getSource();
				annotation = source.getText();
				
				annotationIndex = -1;
				for (int i = 0; i < annotationKeys.length; i ++) {
					if (annotationCollection.getDescriptionForComment(annotationKeys[i], showAnnotationShortcuts, true).equals(annotation)) {
						annotationIndex = i;
						break;
					}
				}

				if (annotation.charAt(0) == '\'' && annotation.charAt(2) == '\'') {
					annotation = annotation.substring(4);
				}
				annotation = annotation.substring(0, annotation.lastIndexOf(" (n="));

				if (e.getButton()==MouseEvent.BUTTON3) {
					
					menu = new JPopupMenu();
					menu.setName(annotationIndex + "");

					menu.add(new AbstractAction("Display all markers annotated with " + annotation) {
						private static final long serialVersionUID = 1L;

						@Override
						public void actionPerformed(ActionEvent e1) {
							String annotation;
							int annotationIndex;
							
							annotationIndex = -1;
							annotation = e1.getActionCommand();
							annotation = annotation.substring("Display all markers annotated with ".length());
							for (int i = 0; i < annotationKeys.length; i ++) {
								if (annotationCollection.getDescriptionForComment(annotationKeys[i], false, false).equals(annotation)) {
									annotationIndex = i;
									break;
								}
							}

							if (!saveClusterFilterAndAnnotationCollection()) {
								return;
							}

							indexOfAnnotationUsedAsMarkerList = -2;
							markerIndexBak = markerIndex;
							markerList = annotationCollection.getMarkerLists(annotationKeys[annotationIndex]);
							commentList = new String[markerList.length];
							loadMarkerDataFromList(0);
							displayIndex(navigationField);
							updateGUI();
						}

						public boolean isEnabled() {
							return indexOfAnnotationUsedAsMarkerList != -2;
						}
					});

					menu.add(new AbstractAction("Revert marker list to original list") {
						private static final long serialVersionUID = 1L;

						@Override
						public void actionPerformed(ActionEvent e2) {
							if (!saveClusterFilterAndAnnotationCollection()) {
								return;
							}

							for (int i = 0; i < annotationCheckBoxes.length; i ++) {
								annotationCheckBoxes[i].setEnabled(true);
							}
							indexOfAnnotationUsedAsMarkerList = -1;
							markerList = masterMarkerList;
							commentList = masterCommentList;
							loadMarkerDataFromList(markerIndexBak);
							displayIndex(navigationField);
							updateGUI();
						}
						
						@Override
						public boolean isEnabled() {
							return indexOfAnnotationUsedAsMarkerList != -1;
						}
					});

					menu.add(new AbstractAction((indexOfAnnotationUsedAsMarkerList == -1?"Limit":"Unlimit")+" Current List to Those with Annotation " + annotation) {
						private static final long serialVersionUID = 1L;

						@Override
						public void actionPerformed(ActionEvent e3) {
							String annotation;
							int annotationIndex;
							
							annotationIndex = -1;
							annotation = e3.getActionCommand();
							annotation = annotation.substring(((indexOfAnnotationUsedAsMarkerList == -1?"Limit":"Unlimit")+" Current List to Those with Annotation ").length());
							for (int i = 0; i < annotationKeys.length; i ++) {
								if (annotationCollection.getDescriptionForComment(annotationKeys[i], false, false).toLowerCase().equals(annotation)) {
									annotationIndex = i;
									break;
								}
							}

							if (indexOfAnnotationUsedAsMarkerList == -1) {
								for (int i = 0; i < annotationIndex; i ++) {
									annotationCheckBoxes[i].setEnabled(false);
								}
								for (int i = annotationIndex + 1; i < annotationCheckBoxes.length; i ++) {
									annotationCheckBoxes[i].setEnabled(false);
								}
								indexOfAnnotationUsedAsMarkerList = annotationIndex;
							} else if (indexOfAnnotationUsedAsMarkerList == annotationIndex) {
								for (int i = 0; i < annotationCheckBoxes.length; i ++) {
									annotationCheckBoxes[i].setEnabled(true);
								}
								indexOfAnnotationUsedAsMarkerList = -1;
							}
						}
					});
					
					menu.show(source, e.getX(), e.getY());
				}
			}
		};

		maxWidth = 200;
		annotationCheckBoxes = new JCheckBox[annotationKeys.length];
		for (int i=0; annotationKeys != null && i < annotationKeys.length; i++) {
			annotationCheckBoxes[i] = new JCheckBox(annotationCollection.getDescriptionForComment(annotationKeys[i], showAnnotationShortcuts, true));
			maxWidth = Math.max(annotationCheckBoxes[i].getPreferredSize().width, maxWidth);
			annotationCheckBoxes[i].setBackground(Color.WHITE);
			if (annotationCollection.markerHasAnnotation(markerList[markerIndex], annotationKeys[i])) {
				annotationCheckBoxes[i].setSelected(true);
			}
			annotationCheckBoxes[i].addItemListener(new ItemListener() {
				public void itemStateChanged(ItemEvent itemEvent) {
					if (! isInitilizing) {
						JCheckBox checkBox;
						char currentKey;
						
						checkBox = (JCheckBox)itemEvent.getSource();
						currentKey = checkBox.getText().charAt(1);
						
				        if (checkBox.getModel().isSelected()) {
				        	annotationCollection.addAnnotationForMarker(markerList[markerIndex], currentKey);
				        } else {
				        	annotationCollection.removeAnnotationForMarker(markerList[markerIndex], currentKey);
				        }
				        checkBox.setText(annotationCollection.getDescriptionForComment(currentKey, showAnnotationShortcuts, true));
				        updateAnnotationPanelFilterRadioButtons();
//				        isAnnotationUpdated = true;
				        setAnnotationUpdated(true);
					}
				}
			});
			annotationCheckBoxes[i].addMouseListener(mouseListenerForAnnotationCheckBoxes);
			
			annotationPanelLowerPart.add(annotationCheckBoxes[i]);
			annotationPanelLowerPartLayout.putConstraint(SpringLayout.WEST, annotationCheckBoxes[i], 5, SpringLayout.WEST, annotationPanelLowerPart);
			annotationPanelLowerPartLayout.putConstraint(SpringLayout.NORTH, annotationCheckBoxes[i], currentHorizontalPos, SpringLayout.NORTH, annotationPanelLowerPart);
//			horizontalPos += (horizontalMargin + annotationCheckBoxes[i].getSize().height);
			currentHorizontalPos += (horizontalMargin + componentHeight);
			removeButton = new JButton(Grafik.getImageIcon("images/delete2sm.png", true));
			removeButton.setActionCommand("removeButton" + i);
			removeButton.setBorder(null);
			removeButton.setSize(6, 6);
			removeButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					String commandText;
					int index;

					commandText = e.getActionCommand();
					index = Integer.parseInt(commandText.substring("removeButton".length(), commandText.length()));
					annotationCollection.removeAnnotation(proj, annotationKeys[index]);
					removeAnnotationFromMaps(annotationKeys[index]);
					annotationKeys = annotationCollection.getKeys();
//					isAnnotationUpdated = true;
					setAnnotationUpdated(true);
					annotationPanelLower();
				}
			});
			annotationPanelLowerPart.add(removeButton);
			annotationPanelLowerPartLayout.putConstraint(SpringLayout.WEST, removeButton, 5, SpringLayout.EAST, annotationCheckBoxes[i]);
			annotationPanelLowerPartLayout.putConstraint(SpringLayout.NORTH, removeButton, 7, SpringLayout.NORTH, annotationCheckBoxes[i]);
		}
		
		addAnnotationField = new JTextField(DEFAULT_MESSAGE);
		
		addAnnotationField.addFocusListener(new FocusListener() {
			public void focusLost(FocusEvent e) {
				if (showAnnotationShortcuts) {
					activateAllAnnotationMaps();
				}
				if (addAnnotationField.getText().equals("")) {
					addAnnotationField.setText(DEFAULT_MESSAGE);
					Font font = addAnnotationField.getFont();
					addAnnotationField.setFont(new Font(font.getFontName(), Font.ITALIC, font.getSize()));
				}
				
			}
			public void focusGained(FocusEvent e) {
				if (addAnnotationField.getText().equals(DEFAULT_MESSAGE)) {
					addAnnotationField.setText("");
					Font font = addAnnotationField.getFont();
					addAnnotationField.setFont(new Font(font.getFontName(), Font.PLAIN, font.getSize()));
				}
				deactivateAllAnnotationMaps();
			}
		});
		
		addAnnotationField.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				String str;
				char c;
				
				str = ((JTextField)e.getSource()).getText();
//				c = str.toLowerCase().charAt(0);
//				str = str.substring(1);
				c = AnnotationCollection.assignKey(str, annotationCollection);
				if (c == 0) {
					log.report("Failed to add the new annotation key '" + str + "', because cannot automatically assign a shortcut for it.");
				} else {
					if (annotationCollection == null) {
						annotationCollection = new AnnotationCollection();
						startAutoSaveToTempFile();
					}
					annotationCollection.addAnnotation(c, str);
					addAnnotationToMaps(c);
					annotationKeys = annotationCollection.getKeys();
					annotationPanelLower();
//					isAnnotationUpdated = true;
					setAnnotationUpdated(true);
				}
			}
		});

		annotationPanelLowerPart.add(addAnnotationField);
		annotationPanelLowerPartLayout.putConstraint(SpringLayout.WEST, addAnnotationField, 5, SpringLayout.WEST, annotationPanelLowerPart);
		annotationPanelLowerPartLayout.putConstraint(SpringLayout.NORTH, addAnnotationField, currentHorizontalPos + 3, SpringLayout.NORTH, annotationPanelLowerPart);
//		horizontalPos += (horizontalMargin + addAnnotationField.getSize().height);
//		horizontalPos += (horizontalMargin + componentHeight);

		// TODO This could be fixed such that annotationPanel is always a fixed width and we have an inner panel with horizontal scroll bar, but this has proved too complex for the moment
//		annotationPanelLowerPartLayout.getConstraints(annotationPanelLowerPart).setConstraint(SpringLayout.EAST, Spring.constant(250, 250, 250));
//		annotationPanel.setPreferredSize(new Dimension(200, 500));
		annotationPanelLowerPart.setPreferredSize(new Dimension(maxWidth+42, ((annotationKeys == null?0:annotationKeys.length)+4)*22));

		annotationPanelLowerPart.validate();
		
		return annotationPanelLowerPart;
	}

	private void addAnnotationToMaps(char c) {
		scatPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke((int)(c+"").toUpperCase().charAt(0), 0), "annotation\t"+c);
		scatPanel.getActionMap().put("annotation\t"+c, new AnnotationAction(this, c, true));

		scatPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke((int)(c+"").toUpperCase().charAt(0), InputEvent.SHIFT_MASK), "annotation\tShift_"+c);
		scatPanel.getActionMap().put("annotation\tShift_"+c, new AnnotationAction(this, c, false));
	}
	
	private void activateAllAnnotationMaps() {
		for (int i = 0; i < annotationKeys.length; i++) {
			addAnnotationToMaps(annotationKeys[i]);
		}
	}
	
	private void deactivateAllAnnotationMaps() {
		for (int i = 0; i < annotationKeys.length; i++) {
			removeAnnotationFromMaps(annotationKeys[i]);
		}
	}

	private void removeAnnotationFromMaps(char c) {
		scatPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).remove(KeyStroke.getKeyStroke((int)(c+"").toUpperCase().charAt(0), 0));
		scatPanel.getActionMap().remove("annotation\t"+c);
	}	
	
	private void inputMapAndActionMap(JComponent comp) {
		InputMap inputMap;
		ActionMap actionMap;
		
		inputMap = comp.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		actionMap = comp.getActionMap();
		
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.ALT_MASK), ALT_UP);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.ALT_MASK), ALT_DOWN);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), ALT_LEFT);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), ALT_RIGHT);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_HOME, InputEvent.CTRL_MASK), FIRST);
//		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, InputEvent.CTRL_MASK), PREVIOUS);
//		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, InputEvent.CTRL_MASK), NEXT);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, 0), PREVIOUS);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, 0), NEXT);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_END, InputEvent.CTRL_MASK), LAST);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_MASK), SYMMETRY);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.CTRL_MASK), CORRECTION);

		actionMap.put(ALT_UP, new CycleRadio(typeRadioButtons, -1));
		actionMap.put(ALT_DOWN, new CycleRadio(typeRadioButtons, 1));
		actionMap.put(ALT_LEFT, new CycleRadio(colorKeyPanel.getClassRadioButtons(), -1));
		actionMap.put(ALT_RIGHT, new CycleRadio(colorKeyPanel.getClassRadioButtons(), 1));
		actionMap.put(FIRST, new AbstractAction() {
			public static final long serialVersionUID = 4L;

			public void actionPerformed(ActionEvent e) {
				first();
			}
		});
		actionMap.put(PREVIOUS, new AbstractAction() {
			public static final long serialVersionUID = 5L;

			public void actionPerformed(ActionEvent e) {
				previous();
			}
		});
		actionMap.put(NEXT, new AbstractAction() {
			public static final long serialVersionUID = 6L;

			public void actionPerformed(ActionEvent e) {
				next();
			}
		});
		actionMap.put(LAST, new AbstractAction() {
			public static final long serialVersionUID = 7L;

			public void actionPerformed(ActionEvent e) {
				last();
			}
		});
		actionMap.put(SYMMETRY, new AbstractAction() {
			public static final long serialVersionUID = 7L;

			public void actionPerformed(ActionEvent e) {
				symmetryBox.setSelected(!symmetryBox.isSelected());
			}
		});
		actionMap.put(CORRECTION, new AbstractAction() {
			public static final long serialVersionUID = 7L;

			public void actionPerformed(ActionEvent e) {
				correctionBox.setSelected(!correctionBox.isSelected());
			}
		});

//		panel.setActionMap(actionMap);
	}
	
	public void first() {
//		markerIndex = 0;
		markerIndex = getAvailableMarker(false, false, showAllMarkersOrNot, showAnnotatedOrUnannotated);
		finishProcessing();
	}
	
	public void previous() {
		markerIndex = getAvailableMarker(false, true, showAllMarkersOrNot, showAnnotatedOrUnannotated);
		finishProcessing();
	}
	
	
	public void next() {
		markerIndex = getAvailableMarker(true, true, showAllMarkersOrNot, showAnnotatedOrUnannotated);
		finishProcessing();
	}

	public void last() {
//		markerIndex = markerList.length - 1;
		markerIndex = getAvailableMarker(true, false, showAllMarkersOrNot, showAnnotatedOrUnannotated);
		finishProcessing();
	}
	
	public void finishProcessing() {
		updateMarkerIndexHistory();
		displayIndex(navigationField);
		scatPanel.setPointsGeneratable(true);
		scatPanel.setQcPanelUpdatable(true);
		setCurrentClusterFilter();
		updateGUI();
		displayClusterFilterIndex();
	}

	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();
		String filename;
		int count;

		if (command.equals(FIRST)) {
			first();
		} else if (command.equals(PREVIOUS)) {
			previous();
		} else if (command.equals(NEXT)) {
			next();
		} else if (command.equals(LAST)) {
			last();
		} else if (command.equals(CLUSTER_FILTER_BACKWARD)) {
			if (clusterFilterCollection.getSize(getMarkerName())>0) {
				scatPanel.rectangles[currentClusterFilter].setColor((byte)7);
				//scatPanel.setRectangles();
				setCurrentClusterFilter((byte) Math.max(currentClusterFilter-1, 0));
//				scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
				//clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+" of "+clusterFilterCollection.getSize(getMarkerName()));
				displayClusterFilterIndex();
//				annotationPanelLowerPart();
				scatPanel.setPointsGeneratable(true);	// why not scatPanel.setPointsGeneratable(false);
				//scatPanel.setUpdateQcPanel(true);
				updateGUI();
			}
		} else if (command.equals(CLUSTER_FILTER_FORWARD)) {
			if (clusterFilterCollection.getSize(getMarkerName())>0) {
				scatPanel.rectangles[currentClusterFilter].setColor((byte)7);
				//scatPanel.setRectangles();
				setCurrentClusterFilter((byte) Math.min(currentClusterFilter+1, (clusterFilterCollection.getSize(getMarkerName())==0?0:(clusterFilterCollection.getSize(getMarkerName())-1))));
//				scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
				//clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+" of "+clusterFilterCollection.getSize(getMarkerName()));
				displayClusterFilterIndex();
//				annotationPanelLowerPart();
				scatPanel.setPointsGeneratable(true);	// why not scatPanel.setPointsGeneratable(false);
				//scatPanel.setUpdateQcPanel(true);
				updateGUI();
			}
		} else if (command.equals(CLUSTER_FILTER_DELETE)) {
			if (clusterFilterCollection.getSize(getMarkerName())>0) {
				//scatPanel.clearRectangles(currentClusterFilter);
				clusterFilterCollection.deleteClusterFilter(getMarkerName(), currentClusterFilter);
				setCurrentClusterFilter((byte) Math.min(currentClusterFilter, clusterFilterCollection.getSize(getMarkerName())-1));
//				startAutoSave();
				setClusterFilterUpdated(true);
				//clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+" of "+clusterFilterCollection.getSize(getMarkerName()));
				displayClusterFilterIndex();
//				annotationPanelLowerPart();
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);
				//scatPanel.repaint();
				//scatPanel.paintAgain();
				scatPanel.generateRectangles();
//				if (clusterFilterCollection.getSize(getMarkerName())>0) {
//					scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
//				}
				updateGUI();
			}
		} else if (command.equals(CAPTURE)) {
			count = 1;
			do {
				filename = markerList[markerIndex]+"_"+MarkerData.TYPES[plot_type][0]+"-"+MarkerData.TYPES[plot_type][1]+"_"+null+(count==1?"":" v"+count);
				filename = ext.replaceWithLinuxSafeCharacters(filename, true);
				count++;
			} while (new File(proj.getProjectDir()+filename+".png").exists());
			scatPanel.screenCapture(proj.getProjectDir()+filename+".png");
		} else if (command.equals(DUMP)) {
			count = 1;
			do {
				filename = markerList[markerIndex]+"_dump"+(count==1?"":" v"+count);
				filename = ext.replaceWithLinuxSafeCharacters(filename, true);
				count++;
			} while (new File(proj.getProjectDir()+filename+".xln").exists());
//			markerData[markerIndex].writeToFile(samples, proj.getProjectDir()+filename+".xln");
			getCurrentMarkerData().dump(sampleData, proj.getProjectDir()+filename+".xln", samples, false);
		} else if (command.equals(MASK_MISSING) || command.equals(UNMASK_MISSING)) {
			maskMissing = !maskMissing;
			((JButton)ae.getSource()).setText(maskMissing?UNMASK_MISSING:MASK_MISSING);
			scatPanel.setPointsGeneratable(true);
			scatPanel.setQcPanelUpdatable(true);
			updateGUI();
		} else {
			log.reportError("Error - unknown command '"+command+"'");
		}
	}

	private int getAvailableMarker(boolean forwardOrBackward, boolean firstOrLastInThatDirection, boolean allMarkersOrOnlyThoseAnnotatedOrUnannotated, boolean annotatedOrUnannotated) {
		int result;
		int begin;
		int end;
		int step;


		result = markerIndex;
		if (forwardOrBackward) {
			if (firstOrLastInThatDirection) {
				begin = markerIndex + 1;
				end = markerList.length;
				step = 1;
			} else {
				begin = markerList.length - 1;
				end = markerIndex;
				step = -1;
			}
		} else {
			if (firstOrLastInThatDirection) {
				begin = markerIndex - 1;
				end = -1;
				step = -1;
			} else {
				begin = 0;
				end = markerIndex;
				step = 1;
			}
		}
		for (int i = begin; i != end; i = i + step) {
			if (indexOfAnnotationUsedAsMarkerList < 0) {
				if (allMarkersOrOnlyThoseAnnotatedOrUnannotated || !(annotatedOrUnannotated ^ annotationCollection.markerHasAnyAnnotation(markerList[i]))) {
					result = i;
					break;
				}
			} else if (annotationCollection.markerHasAnnotation(markerList[i], annotationKeys[indexOfAnnotationUsedAsMarkerList])) {
				result = i;
				break;
			}
		}

		return result;
	}


//	private int getAvailableMarker(boolean forwardOrBackward, boolean firstOrLastInThatDirection, boolean allMarkersOrOnlyThoseAnnotatedOrUnannotated, boolean annotatedOrUnannotated) {
//		int result;
//
//		result = markerIndex;
//		if (allMarkersOrOnlyThoseAnnotatedOrUnannotated && indexOfAnnotationControllingMarkerList == -1) {
//			if (firstOrLastInThatDirection) {
//				if (forwardOrBackward) {
//					result = Math.min(markerIndex + 1, markerList.length - 1);
//				} else {
//					result = Math.max(markerIndex - 1, 0);
//				}
//			} else {
//				if (forwardOrBackward) {
//					result = markerList.length - 1;
//				} else {
//					result = 0;
//				}
//			}
//		} else {
//			if (firstOrLastInThatDirection) {
//				if (forwardOrBackward) {
//					if (indexOfAnnotationControllingMarkerList == -1) {
//						if (annotatedOrUnannotated) {
//							for (int i = markerIndex + 1; i < isAnnotated.length; i++) {
//								if (isAnnotated[i]) {
//									result = i;
//									break;
//								}
//							}
//						} else {
//							for (int i = markerIndex + 1; i < isAnnotated.length; i++) {
//								if (! isAnnotated[i]) {
//									result = i;
//									break;
//								}
//							}
//						}
//					} else {
//						for (int i = markerIndex + 1; i < isAnnotated.length; i++) {
//							if (annotationCollection.markerHasAnnotation(markerList[i], annotationKeys[indexOfAnnotationControllingMarkerList])) {
//								result = i;
//								break;
//							}
//						}
//					}
//				} else {
//					if (indexOfAnnotationControllingMarkerList == -1) {
//						if (annotatedOrUnannotated) {
//							for (int i = markerIndex - 1; i >= 0; i--) {
//								if (isAnnotated[i]) {
//									result = i;
//									break;
//								}
//							}
//						} else {
//							for (int i = markerIndex - 1; i >= 0; i--) {
//								if (! isAnnotated[i]) {
//									result = i;
//									break;
//								}
//							}
//						}
//					} else {
//						for (int i = markerIndex - 1; i >= 0; i--) {
//							if (annotationCollection.markerHasAnnotation(markerList[i], annotationKeys[indexOfAnnotationControllingMarkerList])) {
//								result = i;
//								break;
//							}
//						}
//					}
//				}
//			} else {
//				if (forwardOrBackward) {
//					if (indexOfAnnotationControllingMarkerList == -1) {
//						if (annotatedOrUnannotated) {
//							for (int i = isAnnotated.length - 1; i >= 0; i--) {
//								if (isAnnotated[i]) {
//									result = i;
//									break;
//								}
//							}
//						} else {
//							for (int i = isAnnotated.length - 1; i >= 0; i--) {
//								if (! isAnnotated[i]) {
//									result = i;
//									break;
//								}
//							}
//						}
//					} else {
//						for (int i = isAnnotated.length - 1; i >= 0; i--) {
//							if (annotationCollection.markerHasAnnotation(markerList[i], annotationKeys[indexOfAnnotationControllingMarkerList])) {
//								result = i;
//								break;
//							}
//						}
//					}
//				} else {
//					if (indexOfAnnotationControllingMarkerList == -1) {
//						if (annotatedOrUnannotated) {
//							for (int i = 0; i < isAnnotated.length; i++) {
//								if (isAnnotated[i]) {
//									result = i;
//									break;
//								}
//							}
//						} else {
//							for (int i = 0; i < isAnnotated.length; i++) {
//								if (! isAnnotated[i]) {
//									result = i;
//									break;
//								}
//							}
//						}
//					} else {
//						for (int i = 0; i < isAnnotated.length; i++) {
//							if (annotationCollection.markerHasAnnotation(markerList[i], annotationKeys[indexOfAnnotationControllingMarkerList])) {
//								result = i;
//								break;
//							}
//						}
//					}
//				}
//			}
//		}
//
//		return result;
//	}

	private boolean loadClusterFilterFiles() {
		String[] otherClusterFilterFiles;
		int choice;
		String[] options = new String[] {"Yes, load and delete old file", "No, just delete old file", "Cancel and close ScatterPlot"};
		String clusterFilterFilename;
		
		clusterFilterFilename = proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME);
		otherClusterFilterFiles = Files.list(proj.getDir(Project.DATA_DIRECTORY), ".tempClusterFilters.ser", jar);
		if (otherClusterFilterFiles.length > 0) {
			choice = JOptionPane.showOptionDialog(null, "Error - either multiple instances of ScatterPlot are running or ScatterPlot failed to close properly\n" +
														"last time. The ability to generate new ClusterFilters will be disabled until this file has been\n" +
														"removed. Do you want to load the contents of the temporary file into memory before it is deleted?",
												  "Error", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
			if (choice == 0) {
				// load the last one in otherClusterFilerFiles[]
				clusterFilterCollection = ClusterFilterCollection.load(proj.getDir(Project.DATA_DIRECTORY) + otherClusterFilterFiles[otherClusterFilterFiles.length-1], jar);
				for (int i=0; i<otherClusterFilterFiles.length; i++) {
					(new File(proj.getDir(Project.DATA_DIRECTORY)+otherClusterFilterFiles[i])).delete();
				}
				startAutoSaveToTempFile();
				setClusterFilterUpdated(true);
			} else if (choice == 1) {
				// load permanent
				if (Files.exists(clusterFilterFilename, jar) ) {
					clusterFilterCollection = ClusterFilterCollection.load(clusterFilterFilename, jar);
				} else {
					clusterFilterCollection = new ClusterFilterCollection();
				}
				for (int i=0; i<otherClusterFilterFiles.length; i++) {
					(new File(proj.getDir(Project.DATA_DIRECTORY)+otherClusterFilterFiles[i])).delete();
				}
			} else {
				return false;
			}
		} else if (Files.exists(clusterFilterFilename, jar) ) {
			clusterFilterCollection = ClusterFilterCollection.load(clusterFilterFilename, jar);
		} else {
			clusterFilterCollection = new ClusterFilterCollection();
		}

		startAutoSaveToTempFile();

		return true;
	}

	public void loadAnnotationCollection() {
		String filename;
		
		filename = proj.getFilename(Project.ANNOTATION_FILENAME, false, false);
		if (new File(filename).exists()) {
			log.report("Loading annotation from: "+filename);
			annotationCollection = (AnnotationCollection) Files.readSerial(filename);
		} else {
			log.report("Could not find annotation file: "+filename);
			annotationCollection = new AnnotationCollection();
			annotationCollection.addAnnotation('e', "Extra heterozygote clusters");
			annotationCollection.addAnnotation('m', "Monomorphic");
		}
		startAutoSaveToTempFile();

		annotationKeys = annotationCollection.getKeys();
	}

	public long getSampleFingerprint() {
		return sampleListFingerprint;
	}

	public Project getProject() {
		return proj;
	}

	public String[] getSamples() {
		return samples;
	}

	public SampleData getSampleData() {
		return sampleData;
	}

//	public MarkerData[] getMarkerData1() {
//		return markerData;
//	}
//
	public int getMarkerIndex() {
		return markerIndex;
	}

	public String getMarkerName() {
		return markerList[markerIndex];
	}

//	public void setCurrentClass(byte newCurrentClass) {
//		currentClass = newCurrentClass;
//	}

	public int getCurrentClass() {
		return colorKeyPanel.getCurrentClass();
	}

	public int getPlotType() {
		return plot_type;
	}

	public byte getPointSize() {
		return size;
	}

	public float getGCthreshold() {
		return gcThreshold;
	}

	public boolean[] getDisplayCents() {
		return displayCents;
	}

	public ClusterFilterCollection getClusterFilterCollection() {
		return clusterFilterCollection;
	}

	public byte getCurrentClusterFilter() {
		while (currentClusterFilter >= clusterFilterCollection.getSize(getMarkerName())) {
			currentClusterFilter--;
		}
		return currentClusterFilter;
	}
	
	public float[][][][] getCents() {
		return cents;
	}
	
	public void setCurrentClusterFilter(byte currentClusterFilter) {
		ClusterFilter currentFilter; 
		
		this.currentClusterFilter = currentClusterFilter;
		
//		if (currentClusterFilter >= 0) {
//		log.report("Before: image=="+(scatPanel.image==null?"null":"BImg")+"\t");
		if (clusterFilterCollection.getSize(getMarkerName())>0) {
			newGenotype.setSelectedIndex(clusterFilterCollection.getGenotype(markerList[markerIndex], currentClusterFilter)+1);
			if (scatPanel.getRectangles()!=null) {
				scatPanel.generateRectangles();
				scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
			}
//		}
//
//		if (clusterFilterCollection.getSize(getMarkerName()) > 0) {
			currentFilter = clusterFilterCollection.getClusterFilters(getMarkerName()).get(currentClusterFilter);
			scatPanel.highlightPoints(getCurrentMarkerData().getHighlightStatus(currentFilter));
		}
//		log.report("After: image=="+(scatPanel.image==null?"null":"BImg"));

	}
	
	public void setCurrentClusterFilter() {
		setCurrentClusterFilter((byte) (clusterFilterCollection.getSize(getMarkerName())-1));
	}
	
	public void updateCentLabels() {
		double[] comp;
		String str;
		MarkerData markerData;

		if (markerDataLoader != null) {
			markerData = getCurrentMarkerData();
			if (markerData.getLRRs() != null) {
				for (int i = 0; i<centList.length; i++) {
					comp = markerData.compareLRRs(cents[i][markerIndex]);
					str = comp[0]+"";
		
					for (int j = 2; j<str.length(); j++) {
						if (str.charAt(j) != '9') {
							str = str.substring(0, Math.min(j+2, str.length()));
						}
			        }
					if (centLabels != null) {
						centLabels[i].setText("LRR corr="+str+", err="+ext.formDeci(comp[1], 3));
					}
		        }
			}
		}
	}

	public boolean maskMissing() {
		return maskMissing;
	}
	
	public void loadMarkerListFromFile() {
		BufferedReader reader;
		Vector<String> markerNames = new Vector<String>();
		Vector<String> markerComments = new Vector<String>();
		String[] line;
		String filename;
		Vector<String> missingMarkers;
		
		missingMarkers = new Vector<String>();
		filename = proj.getFilename(Project.DISPLAY_MARKERS_FILENAME);
		//log.report("filename: "+filename);
		try {
			try {
				reader = Files.getReader(filename, jar, true, false);
				if (reader == null) {
					JOptionPane.showMessageDialog(null, "Failed to load '"+filename+"'; this is the designated filename in the project properties file. You will need to create a list of the markers that you want to review and place them in this file.", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				}
				while (reader.ready()) {
					line = reader.readLine().trim().split("\t", -1);
					if (markerLookup.contains(line[0])) {
						markerNames.add(line[0]);
						markerComments.add(line.length>1?line[1]:"");
					} else {
//						log.reportError("Error - could not find "+line[0]+" in the lookup table");
						missingMarkers.add(line[0]);
					}
				}
				reader.close();
				if (missingMarkers.size() > 0) {
					JOptionPane.showMessageDialog(null, "Error - the following markers were not found in the MarkerSet:\n"+Array.toStr(Array.toStringArray(missingMarkers), "\n"), "Error", JOptionPane.ERROR_MESSAGE);
				}
			} catch (FileNotFoundException fnfe) {
				JOptionPane.showMessageDialog(null, "Error - could not find \""+filename+"\"", "Error", JOptionPane.ERROR_MESSAGE);
				// handle error by closing window
			} catch (Exception e) {
				log.reportError("Error reading file \""+filename+"\"");
				e.printStackTrace();
				return;
			}

			masterMarkerList = Array.toStringArray(markerNames);
			masterCommentList = Array.toStringArray(markerComments);
			
			markerIndex = 0;
			markerIndexHistory = new int[NUM_MARKERS_TO_SAVE_IN_HISTORY];
		} catch (Exception e) {
			log.reportError("Error loading: "+filename);
			e.printStackTrace();
		}
	}
	
	public void loadMarkerDataFromList(int newMarkerIndex) {
		markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerList);
//		markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSameThread(proj, markerList);

		markerIndex = newMarkerIndex;
		updateMarkerIndexHistory();
		previousMarkerIndex = -1;
//		navigationField.getActionListeners()[0].actionPerformed(e);
	}
	
	public boolean markerDataIsActive() {
		return markerDataLoader != null;
	}
	
	public MarkerData getCurrentMarkerData() {
		MarkerData markerData;
		int count;
		
		count = 0;
		while((markerData = markerDataLoader.getMarkerData(markerIndex)) == null) {
			markerDataLoader.requestIndexBeTheNextFileToLoad(markerIndex);
			try {
				Thread.sleep(250);
				count++;
			} catch (InterruptedException ie) {
			}
			if (count > 8 && count % 8 == 0) {
				log.reportError("Error - ScatterPlot has been waiting on markerDataLoader to load " + markerList[markerIndex] + " for " + (count/4) + " secounds");
			}
		}
		
		return markerData;
	}

	public void loadCentroids() {
		String[] markerNames, files;
		Centroids trav;
		float[][][] travCents, targetCents;
		MarkerSet set;
		Vector<float[][][]> v;
		Vector<String> fileList;
		Hashtable<String,String> hash;
		int[] indices;
		
		files = Files.list(proj.getDir(Project.DATA_DIRECTORY), ".cent", jar);
		//log.report("Found "+files.length+" .cent files in "+proj.getDir(Project.DATA_DIRECTORY));
		
		hash = new Hashtable<String,String>();
		for (int i = 0; i<markerList.length; i++) {
			hash.put(markerList[i], i+"");
        }
		set = proj.getMarkerSet();
		markerNames = set.getMarkerNames();
		indices = new int[markerList.length];
		for (int i = 0; i<markerNames.length; i++) {
			if (hash.containsKey(markerNames[i])) {
				indices[Integer.parseInt(hash.get(markerNames[i]))] = i;
			}
        }
		
		v = new Vector<float[][][]>();
		fileList = new Vector<String>(); 
		for (int i = 0; i<files.length; i++) {
			log.report("<", false, true);
			trav = Centroids.load(proj.getDir(Project.DATA_DIRECTORY)+files[i], jar);			
			log.report(">", false, true);
			if (trav.getFingerprint() != set.getFingerprint()) {
				log.reportError("Error - Centroids file '"+proj.getDir(Project.DATA_DIRECTORY) + files[i]+"' does not match up with the fingerprint of the current marker set and therefore will not load; if you don't want to see this error message again, then remove the .cent extension from this file.");
			} else {
				travCents = trav.getCentroids();
				targetCents = new float[markerList.length][][];
				for (int j = 0; j<markerList.length; j++) {
					targetCents[j] = travCents[indices[j]];
                }
				v.add(targetCents);
				fileList.add(ext.rootOf(files[i]));
			}
        }
		log.report("");
		cents = new float[v.size()][][][];
		for (int i = 0; i<cents.length; i++) {
			cents[i] = v.elementAt(i);
        }
		centList = Array.toStringArray(fileList);
	}

	public void updateGUI() {
//		log.report("Entering updateGUI()");
		if (markerDataLoader == null) {
			return;
		}
		
//		if (markerData[markerIndex]==null) {
//			log.report("Marker data is null at index "+markerIndex);
////			create visual that we're loading data, (write to scatter panel "Loading data from file")
//			loadMarkerData(markerIndex);
//		}
		
		if (markerList.length==0) {
			markerName.setText("Error: marker data was not successfully loaded");
			commentLabel.setText("Check to make sure MarkerLookup is synchronized with the current data");

		} else {
			markerName.setText(markerList[markerIndex]);
			commentLabel.setText(commentList[markerIndex]);
		}
		if (plot_type >= 2) {
			symmetryBox.setEnabled(false);
			scatPanel.setSymmetricAxes(false);
		} else {
			symmetryBox.setEnabled(true);
			scatPanel.setSymmetricAxes(symmetryBox.isSelected());
		}
		if (plot_type == 3) {
			boolean recomputed = false;
			for (int i = 0; i<centBoxes.length; i++) {
				if (centBoxes[i].isSelected()) {
					if (!recomputed) {
						getCurrentMarkerData().recompute(cents[i][markerIndex]);
						recomputed = true;
					} else {
						centBoxes[i].setSelected(false);
					}
				}

            }
		}
		if(plot_type==5){
			symmetryBox.setEnabled(true);
			scatPanel.setSymmetricAxes(symmetryBox.isSelected());
		}
		if (fileterRadioButtons != null) {
			updateAnnotationPanelFilterRadioButtons();
			updateAnnotationPanelAnnotationCheckBoxes();
		}

		scatPanel.paintAgain();
		if (markerIndex != previousMarkerIndex) {
			updateCentLabels();
		}
		previousMarkerIndex = markerIndex;
//		log.report("    exiting updateGUI()");
	}

	public void updateColorKey(Hashtable<String,String> hash) {
		colorKeyPanel.updateColorKey(hash);
	}

	public Hashtable<String, String> getDisabledClassValues() {
		return colorKeyPanel.getDisabledClassValues();
	}

	public void updateQcPanel(byte chr, int[] genotype, String[] sex, String[] otherClass) {
		int numCalledGenotypes;
		double callrate;
		JLabel qcPanelLabel;
		//JLabel qcCallRateLabel;u
		//JLabel qcHwePvalueLabel;u
		int[] alleleCounts;
		double hweP;
		double minorAlleleFrequency;
		//int[][] sexContingecyTable = new int[2][2];
		CTable classCount;
		String[] called;
		int currentClass;
		int numIncludedSamples;
		
		numIncludedSamples = 0;
		numCalledGenotypes = 0;
		alleleCounts = new int[3];
		called = new String[genotype.length];
		for (int i=0; i<genotype.length; i++) {
			numIncludedSamples += genotype[i]==-3?0:1;
		}
		if (chr==23) {
			for (int i=0; i<genotype.length; i++) {
				if (genotype[i]<=0) {
					called[i] = "-1";
				} else {
					called[i] = "1";
					if (sex[i].equals("2")) {
						alleleCounts[genotype[i]-1]++;
					}
					numCalledGenotypes++;
				}
			}
		} else if (chr==24) {
			for (int i=0; i<genotype.length; i++) {
				if (genotype[i]<=0) {
					called[i] = "-1";
				} else {
					called[i] = "1";
					numCalledGenotypes++;
				}
			}
		} else {
			for (int i=0; i<genotype.length; i++) {
				if (genotype[i]<=0) {
					called[i] = "-1";
				} else {
					called[i] = "1";
					alleleCounts[genotype[i]-1] ++;
					numCalledGenotypes++;
				}
				
			}
		}
		hweP = AlleleFreq.HWEsig(alleleCounts);
		callrate=(double)numCalledGenotypes/(double)numIncludedSamples*100.0;
		if (alleleCounts[0] > alleleCounts[2]) {
			minorAlleleFrequency = (double)(alleleCounts[1] + alleleCounts[2]) / (alleleCounts[0] + 2 * alleleCounts[1] + alleleCounts[2]);
		} else {
			minorAlleleFrequency = (double)(alleleCounts[0] + alleleCounts[1]) / (alleleCounts[0] + 2 * alleleCounts[1] + alleleCounts[2]);
		}
		
		qcPanel.removeAll();
		qcPanel.repaint();
		
        qcPanelLabel = new JLabel("", JLabel.CENTER);
        qcPanel.add(qcPanelLabel);
        qcPanelLabel = new JLabel("QC Metrics", JLabel.CENTER);
        qcPanelLabel.setFont(new Font("Arial", 0, 20));
        qcPanel.add(qcPanelLabel);
		
        qcPanelLabel = new JLabel("Chromosome: " + chr, JLabel.LEFT);
        qcPanelLabel.setFont(new Font("Arial", 0, 14));
        qcPanel.add(qcPanelLabel);

        qcPanelLabel = new JLabel("Callrate: "+ext.formDeci(callrate, 4)+"%"+"                           ", JLabel.LEFT);
        qcPanelLabel.setFont(new Font("Arial", 0, 14));
        qcPanel.add(qcPanelLabel);

        qcPanelLabel = new JLabel("Minor Allele Frequency: " + ext.prettyP(minorAlleleFrequency), JLabel.LEFT);
        qcPanelLabel.setFont(new Font("Arial", 0, 14));
		qcPanel.add(qcPanelLabel);

        qcPanelLabel = new JLabel("HWE p-value: "+ext.prettyP(hweP), JLabel.LEFT);
        qcPanelLabel.setFont(new Font("Arial", 0, 14));
		qcPanel.add(qcPanelLabel);

		ToolTipManager.sharedInstance().setDismissDelay(100000);

		classCount = new CTable(called, sex);//This is the problem.
		classCount.setCustomNullValues(Array.addStrToArray("-1", CTable.DEFAULT_NULL_VALUES));
		classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"}, {"1","Genotype NOT missing"}}, sampleData.getActualClassColorKey(0));
		qcPanelLabel = new JLabel("Callrate by sex: "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(), false), 1) ), JLabel.LEFT);
		//classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"}, {"1","Genotype NOT missing"}}, Matrix.addRow(sampleData.getActualClassColorKey(0), new String[] {null, "missing"}));
		qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
        qcPanelLabel.setFont(new Font("Arial", 0, 14));
		qcPanel.add(qcPanelLabel);

		classCount = new CTable(CTable.extrapolateCounts(sex, genotype));
		classCount.setCustomNullValues(Array.addStrToArray("-1", CTable.DEFAULT_NULL_VALUES));
		classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(0), new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}});
		qcPanelLabel = new JLabel("Allele Freq by sex: "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(), false), 1)), JLabel.LEFT);
		//classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(0), new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}, {".","Missing"}});
		qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
        qcPanelLabel.setFont(new Font("Arial", 0, 14));
		qcPanel.add(qcPanelLabel);

		currentClass = getCurrentClass();
		if (currentClass>SampleData.BASIC_CLASSES.length && currentClass<sampleData.getBasicClasses().length+sampleData.getNumActualClasses()) {
			classCount = new CTable(called, otherClass);//This is the problem.
			classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"}, {"1","Genotype NOT missing"}}, sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length));
			qcPanelLabel = new JLabel("Callrate by "+sampleData.getClassName(currentClass)+": "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(), false), 1) ), JLabel.LEFT);
			classCount.setCustomNullValues(Array.addStrToArray("-1", Array.addStrToArray("0", CTable.DEFAULT_NULL_VALUES)));
			//classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"}, {"1","Genotype NOT missing"}}, Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length), new String[] {null, "missing"}));
			qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
	        qcPanelLabel.setFont(new Font("Arial", 0, 14));
			qcPanel.add(qcPanelLabel);
	
			classCount = new CTable(CTable.extrapolateCounts(otherClass, genotype));
			classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length), new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}});
			//classCount.replaceIdWithLabel(SampleData.KEYS_FOR_BASIC_CLASSES[1],sampleData.getActualClassColorKey(0));
			qcPanelLabel = new JLabel("Allele Freq by "+sampleData.getClassName(currentClass)+": "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(), false), 1)), JLabel.LEFT);
			classCount.setCustomNullValues(Array.addStrToArray("-1", Array.addStrToArray("0", CTable.DEFAULT_NULL_VALUES)));
			//classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length), new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}});
			qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
	        qcPanelLabel.setFont(new Font("Arial", 0, 14));
			qcPanel.add(qcPanelLabel);
		}

//		qcPanelLabel = new JLabel("Minor Allele Freq: " + (new DecimalFormat("#.####").format(classCount.getMinorAlleleFrequency())), JLabel.LEFT);
//      qcPanelLabel.setFont(new Font("Arial", 0, 14));
//		qcPanel.add(qcPanelLabel);

//		AlleleFreq.calcFrequency(genotypes);
		
		
		/*
		qcPanelLabel = new JLabel(""+sexContingecyTable[0][1], JLabel.LEFT);
        qcPanelLabel.setFont(new Font("Arial", 0, 14));
		qcPanel.add(qcPanelLabel);
		qcPanelLabel = new JLabel(""+sexContingecyTable[1][0], JLabel.LEFT);
        qcPanelLabel.setFont(new Font("Arial", 0, 14));
		qcPanel.add(qcPanelLabel);
		qcPanelLabel = new JLabel(""+sexContingecyTable[1][1], JLabel.LEFT);
        qcPanelLabel.setFont(new Font("Arial", 0, 14));
		qcPanel.add(qcPanelLabel);
		*/

		qcPanel.validate();
	}

	public void updateCurrentClusterFiltersGenotype() {
		byte currentGenotype, newGenotype;

		currentGenotype = clusterFilterCollection.getGenotype(getMarkerName(), currentClusterFilter);
		if (currentGenotype == (ScatterPlot.GENOTYPE_OPTIONS.length - 2)) {
			newGenotype = -1;
		} else {
			newGenotype = (byte) (currentGenotype + 1);
		}
		updateCurrentClusterFilterGenotype(newGenotype, true);
	}

	public void updateCurrentClusterFilterGenotype(byte newGenotypeSelected, boolean updateGenotypeComboBox) {
		clusterFilterCollection.updateGenotype(getMarkerName(), currentClusterFilter, newGenotypeSelected);
		setClusterFilterUpdated(true);
		scatPanel.setPointsGeneratable(true);
		scatPanel.setQcPanelUpdatable(true);
		scatPanel.paintAgain();
		if (updateGenotypeComboBox) {
			newGenotype.setSelectedIndex(newGenotypeSelected + 1);
		}
	}

	public void updateAnnotationPanelFilterRadioButtons() {
		int numAnnotated;
		
		numAnnotated = 0;
		for (int i = 0; i < markerList.length; i++) {
			if (annotationCollection.markerHasAnyAnnotation(markerList[i])) {
				numAnnotated++;
			}
		}
		
        for (int i=0; i<fileterRadioButtons.length; i++) {
        	fileterRadioButtons[i].setText(RADIOBUTTON_TEXTS[i] + " (n=" + (i==0? markerList.length : (i==1? numAnnotated : (markerList.length - numAnnotated))) + ")");
        }
	}

	public void updateAnnotationPanelAnnotationCheckBoxes() {
		isInitilizing = true;
		for (int i=0; i<annotationCheckBoxes.length; i++) {
        	annotationCheckBoxes[i].setText(annotationCollection.getDescriptionForComment(annotationKeys[i], showAnnotationShortcuts, true));
			if (annotationCollection.markerHasAnnotation(markerList[markerIndex], annotationKeys[i])) {
				annotationCheckBoxes[i].setSelected(true);
			} else {
				annotationCheckBoxes[i].setSelected(false);
			}
        }
		isInitilizing = false;
	}

	public void checkAnnotationBox(char c) {
		int index;
		
		index = ext.indexOfChar(c, annotationKeys);
		if (! annotationCheckBoxes[index].isSelected()) {
			annotationCheckBoxes[index].setSelected(true);
			if (annotationAutoAdv) {
				actionPerformed(new ActionEvent(next, 0, NEXT));
			}
		}
	}

	public void uncheckAnnotationBox(char c) {
		int index;
		
		index = ext.indexOfChar(c, annotationKeys);
		if (annotationCheckBoxes[index].isSelected()) {
			annotationCheckBoxes[index].setSelected(false);
		}
	}

	public void displayIndex(JTextField field) {
		field.setText((markerIndex + 1) + " of " + markerList.length);
	}
	
	public void displayClusterFilterIndex() {
		clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))
										+" of "
										+clusterFilterCollection.getSize(getMarkerName()));
	}

	public void startAutoSaveToTempFile() {
		if (autoSave == null) {
			autoSave = new AutoSaveForScatterPlot(clusterFilterCollection, proj.getDir(Project.DATA_DIRECTORY) + sessionID + ".tempClusterFilters.ser", annotationCollection, proj.getDir(Project.DATA_DIRECTORY) + sessionID + ".tempAnnotation.ser", 30);
			new Thread(autoSave).start();
		} else if (clusterFilterCollection != null && autoSave.isClusterFilterNull()) {
			autoSave.addToAutoSave(clusterFilterCollection, proj.getDir(Project.DATA_DIRECTORY) + sessionID + ".tempClusterFilters.ser");
		} else if (annotationCollection != null && autoSave.isAnnotationNull()) {
			autoSave.addToAutoSave(annotationCollection, proj.getDir(Project.DATA_DIRECTORY) + sessionID + ".tempAnnotation.ser");
		}
		autoSave.saveNow();
	}

	public void updateMarkerIndexHistory() {
		// Oldest history is stored in the end of the array
		for (int i = 0; i < markerIndexHistory.length - 1 ; i++) {
			markerIndexHistory[i+1] = markerIndexHistory[i];
		}
		markerIndexHistory[0] = markerIndex;
	}

	public void goBackToIndexHistory(int indexInHistory) {
		if (markerIndexHistory[indexInHistory] == -1 ) {
			log.report("Cannot roll back to more than the history has recorded.");
		} else {
			markerIndex = markerIndexHistory[indexInHistory];
			updateMarkerIndexHistory();
		}
	}

	public void setClusterFilterUpdated(boolean status) {
		isClusterFilterUpdated = status;
		autoSave.setClusterFilterUpdated(status);
	}
	
	public void setAnnotationUpdated (boolean status) {
		isAnnotationUpdated = status;
		autoSave.setAnnotationUpdated(status);
	}
	
	public boolean saveClusterFilterAndAnnotationCollection() {
		String[] options;
		int choice;
		String filename;
		String clusterFilterFilename;

		clusterFilterFilename = proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, false, false);
		options = new String[] {"Yes, overwrite", "No"};
		if (isClusterFilterUpdated) {
			choice = JOptionPane.showOptionDialog(null, "New ClusterFilters have been generated. Do you want to save them to the permanent file?", "Overwrite permanent file?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
			if (choice == 0) {
				clusterFilterCollection.serialize(clusterFilterFilename);
				proj.archiveFile(clusterFilterFilename);
				clusterFilterCollection.serialize(clusterFilterFilename);
				log.report("clusterFilters should be backed up now");
			} else if (choice == -1) {
				return false;
			} else {
				//TODO As a double security, move sessionID + ".tempClusterFilters.ser" to BACKUP_DIRECTORY. But then need to delete it from BACKUP_DIRECTORY at some point of time.
				// need a rule for that. also need the code for deletion.
			}
			setClusterFilterUpdated(false);
		}

		if (isAnnotationUpdated) {
			choice = JOptionPane.showOptionDialog(null, "New Annotations have been generated. Do you want to save them to the permanent file?", "Overwrite permanent file?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
			if (choice == 0) {
				filename = proj.getFilename(Project.ANNOTATION_FILENAME, false, false);
				annotationCollection.serialize(filename);
				proj.archiveFile(filename);
				annotationCollection.serialize(filename);
			} else if (choice == -1) {
				return false;
			}
			setAnnotationUpdated(false);
		}
		
		return true;
	}

	public void windowActivated(WindowEvent e) {}

	public void windowClosed(WindowEvent e) {}

	public void windowClosing(WindowEvent e) {
		int count;
		
		if (!saveClusterFilterAndAnnotationCollection()) {
			return;
		}

		//TODO notify all threads (e.g., MarkerDataLoader) that they need to close

		if (autoSave != null) {
			autoSave.kill();
		}
		if (markerDataLoader != null) {
			markerDataLoader.kill();
			count = 0;
			while (!markerDataLoader.killComplete()) {
				try {
					Thread.sleep(250);
				} catch (InterruptedException ie) {}
				count++;
				if (count > 0) {
					log.reportError("Waiting for markerDataLoader to wind down...");
				}
			}
		}
		
		if (exitOnClose) {
			System.exit(1); // keep this exit, of course
		} else {
			((JFrame)(this.getParent().getParent().getParent())).dispose();
		}
	}

	public void windowDeactivated(WindowEvent e) {}

	public void windowDeiconified(WindowEvent e) {}

	public void windowIconified(WindowEvent e) {}

	public void windowOpened(WindowEvent e) {}

    public static void createAndShowGUI(Project proj, String[] markerList, String[] commentList, boolean exitOnClose) {
        JFrame frame;
    	ScatterPlot scatterPlot;

    	scatterPlot = new ScatterPlot(proj, markerList, commentList, exitOnClose);
    	
    	if (!scatterPlot.failed()) {
	    	frame = new JFrame("Genvisis - ScatterPlot - " + proj.getNameOfProject());
			frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
	        frame.setContentPane(scatterPlot);
			frame.addWindowListener(scatterPlot);
	
	        frame.pack();
			frame.setSize(1200, 870);
			frame.setVisible(true);
    	}
    }

    public static void main(String[] args) {
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
            	createAndShowGUI(new Project(cnv.Launch.getDefaultDebugProjectFile(true), false), null, null, true);
            }
        });
	}
    
    public PrincipalComponentsResiduals getPcResids() {
		return pcResids;
	}

	private  PrincipalComponentsResiduals loadPcResids() {
		String pcFile = proj.getFilename(Project.INTENSITY_PC_FILENAME);
		PrincipalComponentsResiduals pcResids;
		if (Files.exists(proj.getProjectDir() + ext.removeDirectoryInfo(pcFile))) {
			proj.getLog().report("Info - loading " + ext.removeDirectoryInfo(pcFile));
			pcResids = new PrincipalComponentsResiduals(proj, ext.removeDirectoryInfo(pcFile), null, Integer.parseInt(proj.getProperty(Project.INTENSITY_PC_NUM_COMPONENTS)), false, 0, false, false, null);
			setNumComponents(Math.min(numComponents, pcResids.getTotalNumComponents()));
		} else {
			proj.getLog().report("Info - did not find " + proj.getProjectDir() + ext.removeDirectoryInfo(pcFile));
			pcResids = null;
		}
		return pcResids;
	}

	public int getNumComponents() {
		return numComponents;
	}

	public void setNumComponents(int numComponents) {
		this.numComponents = numComponents;
	}

	public double getstdevFilter() {
		return stdevFilter;
	}
	

	public double getCorrectionRatio() {
		return correctionRatio;
	}

	public JCheckBox getCorrectionBox() {
		return correctionBox;
	}
	
	
}
