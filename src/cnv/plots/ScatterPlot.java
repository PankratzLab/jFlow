package cnv.plots;

import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;

import stats.CTable;
import stats.ContingencyTable;
import stats.ProbDist;

import java.awt.*;
import java.awt.event.*;

import cnv.filesys.*;
import cnv.gui.AnnotationAction;
import cnv.gui.AutoSaveClusterFilterCollection;
import cnv.gui.CycleRadio;
import cnv.manage.MarkerDataLoader;
import common.*;
import cnv.var.*;

//public class ScatterPlot extends JFrame implements ActionListener, WindowListener, KeyListener {
public class ScatterPlot extends JFrame implements ActionListener, WindowListener {
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
	private static final String CLUSTER_FILTER_BACKWARD = "Backward";
	private static final String CLUSTER_FILTER_FORWARD = "Forward";
	private static final String CLUSTER_FILTER_DELETE = "Delete";
	private static final String CAPTURE = "Screen capture";
	private static final String DUMP = "Dump raw data";
	private static final String SHOW_ALL = "Show all";
	private static final String SHOW_ANNOTATED_ONLY = "Show annotated only";
	private static final String SHOW_UNANNOTATED_ONLY = "Show unannotated only";
	private static final String[] RADIOBUTTON_TEXTS = new String[] {SHOW_ALL, SHOW_ANNOTATED_ONLY, SHOW_UNANNOTATED_ONLY};
	public static final String MASK_MISSING = "Mask missing values";
	public static final String UNMASK_MISSING = "Unmask missing values";
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String DEFAULT_MESSAGE = "enter new annotation here";


	private JButton first, previous, next, last;
	private JTextField navigationField;
	private JPanel classPanel;
//	private JPanel legendPanel;
//	private JPanel bottomPanel;
	private ScatterPanel scatPanel;
	private JLabel sizeLabel;
	private JLabel gcLabel;
//	private JPanel typePanel;
//	private JTabbedPane tabbedPane;
	private JPanel qcPanel;
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
	private boolean[] isAnnotated;
	private int annotated;
	private char[] annotationKeys;
	private JComboBox<String> newGenotype;
	private boolean isInitilizing;

	private Project proj;
//	private MarkerData[] markerData;
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
	private int previousMarkerIndex;
	//private JLabel markerName, commentLabel;
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
	private boolean maskMissing;
	private ClusterFilterCollection clusterFilterCollection;
	private AnnotationCollection annotationCollection;
	private byte currentClusterFilter;
	private JTextField clusterFilterNavigation;
	private boolean clusterFilterCollectionUpdated;
	private boolean annotationUpdated;
	private String sessionID;
	private AutoSaveClusterFilterCollection autoSaveCFC;
	private JRadioButton[] typeRadioButtons;
	private MarkerDataLoader markerDataLoader;
	private Thread thread2;
    private ColorKeyPanel colorKeyPanel;
	private Color[] colorScheme;
	
	
	public ScatterPlot(Project project, String[] initMarkerList, String[] initCommentList) {
		super("ScatterPlot");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		long time = new Date().getTime();
		
		proj = project;
		jar = proj.getJarStatus();
		size = DEFAULT_SIZE;
		gcThreshold = (float)DEFAULT_GC_THRESHOLD/100f;
//		clusterFilterCollectionUpdated = false;

		SampleList list = proj.getSampleList();
//		SampleList list2 = proj.getSampleList2();
		samples = list.getSamples();
		sampleListFingerprint = list.getFingerprint();
		sampleData = proj.getSampleData(2, true);
		markerLookup = proj.getMarkerLookup();
		
		masterMarkerList = initMarkerList;
		masterCommentList = initCommentList;
		if (masterMarkerList == null) {
			loadMarkerListFromFile();
			if (masterMarkerList == null) {
				return;
			}
		}
		if (masterCommentList == null) {
			masterCommentList = Array.stringArray(masterMarkerList.length, "");
		}
		
		markerList = masterMarkerList;
		commentList = masterCommentList;
		isAnnotated = new boolean[markerList.length];
		annotated = 0;
		showAllMarkersOrNot = true;
		showAnnotatedOrUnannotated = true;
		showAnnotationShortcuts = true;
		isInitilizing = true;
		
		loadMarkerDataFromList();
		
		System.err.println("3\t"+ext.getTimeElapsed(time));
		loadCentroids();
		addWindowListener(this);
		sessionID = (new Date().getTime()+"").substring(5);
		loadClusterFilterFiles();
		autoSaveCFC = null;
		System.err.println("4\t"+ext.getTimeElapsed(time));
		
//		annotationShortcuts = true;
		annotationAutoAdv = true;

		scatPanel = new ScatterPanel(this);
		loadAnnotationCollection();
		colorScheme = scatPanel.getColorScheme();
		getContentPane().add(scatPanel, BorderLayout.CENTER);
		getContentPane().add(markerPanel(), BorderLayout.NORTH);
		getContentPane().add(eastPanel(), BorderLayout.EAST);
//		getContentPane().add(colorPanel(), BorderLayout.SOUTH);
		colorKeyPanel = new ColorKeyPanel(sampleData, scatPanel);
		getContentPane().add(colorKeyPanel, BorderLayout.SOUTH);

		inputMapAndActionMap();

		scatPanel.setPointsGeneratable(true);//zx
		scatPanel.setQcPanelUpdatable(true);//zx???
//		scatPanel.generateRectangles();
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
		updateAnnotationPanel();
		activateAllAnnotationMaps();
		clusterFilterCollectionUpdated = false;
		annotationUpdated = false;

//    	newGenotype.setSelectedIndex(clusterFilterCollection.getGenotype(getMarkerName(), currentClusterFilter)+1);
		symmetryBox.setSelected(true);
		if (centList.length > 0) {
			centBoxes[0].setSelected(true);
		}
		

		next.getInputMap().put(KeyStroke.getKeyStroke("space"), NEXT);
//		next.setActionMap(actionMap);
//		previous.setActionMap(actionMap);
		scatPanel.grabFocus();

		setBounds(20, 20, 1000, 720);
		setVisible(true);
	}
	
	private JPanel markerPanel() {
		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new GridLayout(3, 1));
		//markerName = new JLabel("", JLabel.CENTER);
		markerName = new JTextField("");//zx
		markerName.setBorder(null);//zx
		markerName.setEditable(false);//zx
		markerName.setBackground(BACKGROUND_COLOR);//zx
		markerName.setHorizontalAlignment(JTextField.CENTER);
		markerName.setFont(new Font("Arial", 0, 20));
		descrPanel.add(markerName);

		//commentLabel = new JLabel("", JLabel.CENTER);
		commentLabel = new JTextField("");//zx
		commentLabel.setBorder(null);//zx
		commentLabel.setEditable(false);//zx
		commentLabel.setBackground(BACKGROUND_COLOR);//zx
		commentLabel.setHorizontalAlignment(JTextField.CENTER);//zx
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
		//navigationField.setEditable(false);//zx
		navigationField.setBackground(BACKGROUND_COLOR);//zx
		navigationField.addActionListener(new ActionListener() {
	        public void actionPerformed(ActionEvent e) {
				try {
					int trav = Integer.valueOf(((JTextField)e.getSource()).getText().split("[\\s]+")[0]).intValue()-1;
					if (trav>=0&&trav<markerList.length) {
						markerIndex = trav;
					}
				} catch (NumberFormatException nfe) {}
				displayIndex((JTextField)e.getSource());
				scatPanel.setPointsGeneratable(true);//zx
				scatPanel.setQcPanelUpdatable(true);//zx???
				setCurrentClusterFilter();
				updateGUI();
				displayClusterFilterIndex();
//				annotationPanelLowerPart();
				updateAnnotationPanel();
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
		JPanel tabbedPanel;
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
		createAnnotationPanel();
		annotationScrollPane = new JScrollPane(annotationPanel);
//		annotationPanel.addKeyListener(this);
//		annotationPanel.setFocusable(true);
//		annotationPanel.requestFocus();
		tabbedPane.addTab("Annotation", null, annotationScrollPane, "Annotation to this marker");
		tabbedPane.setMnemonicAt(2, KeyEvent.VK_3);

		qcPanel = new JPanel();//zx
		qcPanel.setLayout(new GridLayout(8, 1));//zx
		qcPanel.setBackground(BACKGROUND_COLOR);//zx

		tabbedPanel = new JPanel();
		tabbedPanel.setLayout(new BoxLayout(tabbedPanel, BoxLayout.PAGE_AXIS));
		tabbedPanel.setBackground(BACKGROUND_COLOR);
		tabbedPanel.add(tabbedPane);
		tabbedPanel.add(clusterFilterPanel());
//		tabbedPanel.add(plotTypePanel());
		tabbedPanel.add(qcPanel);//zx

		return tabbedPanel;
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
				scatPanel.setPointsGeneratable(true);//zx
				scatPanel.setQcPanelUpdatable(false);//zx???
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
				scatPanel.setPointsGeneratable(true);//zx
				scatPanel.setQcPanelUpdatable(true);//zx
				scatPanel.paintAgain();
				//qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");//zx
			}
		});
//		tabPanel.add(slider, gbc);
		gcSliderPanel.add(slider);

		return gcSliderPanel;
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
		//navigationField.setEditable(false);//zx
		clusterFilterNavigation.setBackground(BACKGROUND_COLOR);//zx
		clusterFilterNavigation.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent focusevent) {}

			public void focusLost(FocusEvent fe) {
				try {
					int trav = Integer.valueOf(((JTextField)fe.getSource()).getText().split("[\\s]+")[0]).intValue()-1;
					if (trav>=0&&trav<clusterFilterCollection.getSize(getMarkerName())) {
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
//				updateAnnotationPanel();
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

    	newGenotype = new JComboBox<String>(new String[] {"-","A/A","A/B","B/B"});
		ActionListener newGenotypeListener = new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				byte newGenotypeSelected;
				if ((String)newGenotype.getSelectedItem()=="-") {
					newGenotypeSelected=(byte)-1;
				} else if ((String)newGenotype.getSelectedItem()=="A/A") {
					newGenotypeSelected=(byte)0;
				} else if ((String)newGenotype.getSelectedItem()=="A/B") {
					newGenotypeSelected=(byte)1;
				} else if ((String)newGenotype.getSelectedItem()=="B/B") {
					newGenotypeSelected=(byte)2;
				} else {
					newGenotypeSelected=(byte)-9;
				}
//				byte index = (byte) (clusterFilterCollection.getSize(getMarkerName())-1);
				if ((clusterFilterCollection.getGenotype(getMarkerName(), currentClusterFilter)!=newGenotypeSelected)) {					//??????
					saveClusterFilterCollection();
				}
				clusterFilterCollection.updateGenotype(getMarkerName(), currentClusterFilter, newGenotypeSelected);
//				System.out.println("e.getActionCommand: "+e.getActionCommand()+"\te.getWhen: "+e.getWhen()+"\te.getSource: "+e.getSource());
//				if (e.getActionCommand().equals("comboBoxChanged")) {					//??????
//				if (e.getWhen()>1000) {					//??????
//					saveClusterFilterCollection();
//				}
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);//zx???
				//updateGUI();
				scatPanel.paintAgain();
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
//		tabPanel.setLayout(new BoxLayout(tabPanel, BoxLayout.Y_AXIS));
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

		ItemListener symmetryListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
//				scatPanel.setPointsGenerated(true);//zx ??? Why not true?
//				scatPanel.setUpdateQcPanel(false);//zx ??? Why cannot set to false?
				updateGUI();
			}
		};
		
		symmetryBox = new JCheckBox("Symmetric axes");
		symmetryBox.setFont(new Font("Arial", 0, 14));
		symmetryBox.addItemListener(symmetryListener);
		symmetryBox.setBackground(BACKGROUND_COLOR);

		controlPanel.add(symmetryBox, gbc);
//		tabPanel.add(symmetryBox);
		
		JButton button = new JButton(CAPTURE);
		button.addActionListener(this);
		button.setActionCommand(CAPTURE);
		controlPanel.add(button, gbc);
//		tabPanel.add(button);

		button = new JButton(DUMP);
		button.addActionListener(this);
		button.setActionCommand(DUMP);
		controlPanel.add(button, gbc);
//		tabPanel.add(button);
		
		button = new JButton(MASK_MISSING);
		button.addActionListener(this);
		button.setActionCommand(MASK_MISSING);
		controlPanel.add(button, gbc);
		
		controlPanel.add(plotTypePanel(), gbc);

		return controlPanel;
	}

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
				scatPanel.setPointsGeneratable(true);//zx
				scatPanel.setQcPanelUpdatable(false);//zx
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
		//System.out.println(centList.length+"\n");//zx
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

	private void createAnnotationPanel() {
		annotationPanel = new JPanel();
		annotationPanel.setBackground(Color.WHITE);
		annotationPanel.setLayout(new BoxLayout(annotationPanel, BoxLayout.Y_AXIS));
		annotationPanel.add(annotationPanelUpperPart());
		annotationPanelLowerPart = new JPanel();
		annotationPanelLowerPart.setLayout(new BoxLayout(annotationPanelLowerPart, BoxLayout.Y_AXIS));
		annotationPanelLowerPart.setBackground(BACKGROUND_COLOR);
		annotationPanel.add(annotationPanelLowerPart);
		annotationPanelLowerPart();
		annotationScrollPane = new JScrollPane(annotationPanel);
	}

	private JPanel annotationPanelUpperPart() {
		JPanel annotationPanelControl;
		JPanel subPanel;
		JCheckBox checkBox;

		annotationPanelControl = new JPanel();
		annotationPanelControl.setLayout(new BoxLayout(annotationPanelControl, BoxLayout.Y_AXIS));
		annotationPanelControl.setBackground(BACKGROUND_COLOR);

		subPanel = new JPanel();
		subPanel.setLayout(new BoxLayout(subPanel, BoxLayout.X_AXIS));
		subPanel.setBackground(Color.WHITE);

		checkBox = new JCheckBox("Shortcut");
		checkBox.setBackground(Color.WHITE);
		checkBox.setHorizontalTextPosition(SwingConstants.LEFT);
		checkBox.setSelected(true);
		checkBox.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				JCheckBox checkBox = (JCheckBox) e.getSource();
				if (checkBox.isSelected()) {
					showAnnotationShortcuts = true;
					activateAllAnnotationMaps();
//					annotationPanelLowerPart();
					updateAnnotationPanelAnnotationCheckBoxes();
				} else {
					showAnnotationShortcuts = false;
					deactivateAllAnnotationMaps();
//					annotationPanelLowerPart();
					updateAnnotationPanelAnnotationCheckBoxes();
				}
			}});
//		annotationPanelControl.add(checkBox);
		subPanel.add(checkBox);

		checkBox = new JCheckBox("Auto Advance");
		checkBox.setBackground(Color.WHITE);
		checkBox.setHorizontalTextPosition(SwingConstants.LEFT);
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
//		annotationPanelControl.add(checkBox);
		subPanel.add(checkBox);
		annotationPanelControl.add(subPanel);


//		annotationPanelControl.setVisible(true);
		return annotationPanelControl;
	}

	public void annotationPanelLowerPart() {
		JPanel panel;
		JButton removeButton;
		JLabel label;
		ButtonGroup radioButtonGroup;
		
		annotationPanelLowerPart.removeAll();
		annotationPanelLowerPart.repaint();

		label = new JLabel("Filtering");
		panel = new JPanel();
		panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
//		BoxLayout tmp = new BoxLayout(panel, BoxLayout.Y_AXIS);
//		tmp.setAlignmentX();
//		panel.setAlignmentX(panel.LEFT_ALIGNMENT);
		panel.setBackground(Color.WHITE);
		panel.add(label);
		radioButtonGroup = new ButtonGroup();
		fileterRadioButtons = new JRadioButton[RADIOBUTTON_TEXTS.length];
		for (int i = 0; i < RADIOBUTTON_TEXTS.length; i++) {
			fileterRadioButtons[i] = new JRadioButton(RADIOBUTTON_TEXTS[i] + " (n=" + (i==0? isAnnotated.length : (i==1? annotated : (isAnnotated.length - annotated))) + ")");
			fileterRadioButtons[i].setBackground(Color.WHITE);
			fileterRadioButtons[i].addActionListener(new ActionListener(){
				@Override
				public void actionPerformed(ActionEvent e) {
					String commandText;
					commandText = e.getActionCommand();
					if (commandText.startsWith(SHOW_ALL + " (n=")) {
						showAllMarkersOrNot = true;
					} else if (commandText.startsWith(SHOW_ANNOTATED_ONLY + " (n=")) {
						showAllMarkersOrNot = false;
						showAnnotatedOrUnannotated = true;
					} else {
						showAllMarkersOrNot = false;
						showAnnotatedOrUnannotated = false;
					}
				}
			});
			radioButtonGroup.add(fileterRadioButtons[i]);
			panel.add(fileterRadioButtons[i]);
		}
		fileterRadioButtons[0].setSelected(true);
//		annotationPanel.add(radioButtonGroup);
		annotationPanelLowerPart.add(panel);

		annotationCheckBoxes = new JCheckBox[annotationKeys.length];
		for (int i=0; annotationKeys != null && i < annotationKeys.length; i++) {
			annotationCheckBoxes[i] = new JCheckBox(annotationCollection.getDescriptionForComment(annotationKeys[i], showAnnotationShortcuts));
//			annotationCheckBoxes[i] = new JCheckBox();
			annotationCheckBoxes[i].setBackground(Color.WHITE);
			if (annotationCollection.markerHasAnnotation(markerList[markerIndex], annotationKeys[i])) {
				annotationCheckBoxes[i].setSelected(true);
			}
			annotationCheckBoxes[i].addItemListener(new ItemListener() {
				public void itemStateChanged(ItemEvent itemEvent) {
					if (! isInitilizing) {
						JCheckBox checkBox;
						char currentKey;
						boolean prev;
						
						checkBox = (JCheckBox)itemEvent.getSource();
						currentKey = checkBox.getText().charAt(1);
						
				        if (checkBox.getModel().isSelected()) {
				        	annotationCollection.addAnnotationForMarker(markerList[markerIndex], currentKey);
				        	if (! isAnnotated[markerIndex]) {
				        		annotated ++;
				        	}
				        	isAnnotated[markerIndex] = true;
				        } else {
				        	annotationCollection.removeAnnotationForMarker(markerList[markerIndex], currentKey);
				        	prev = isAnnotated[markerIndex];
				        	isAnnotated[markerIndex] = (annotationCollection.annotationsForMarker(markerList[markerIndex]).length != 0);
				        	if (isAnnotated[markerIndex] != prev) {
				        		annotated --;
				        	}
				        }
				        checkBox.setText(annotationCollection.getDescriptionForComment(currentKey, showAnnotationShortcuts));
				        updateAnnotationPanelFilterRadioButtons();
				        annotationUpdated = true;
					}
				}
			});
			panel = new JPanel();
			panel.setBackground(Color.WHITE);
//			panel.setLayout(new BoxLayout(panel, BoxLayout.X_AXIS));
			panel.add(annotationCheckBoxes[i]);
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
					annotationCollection.removeAnnotation(annotationKeys[index]);
					removeAnnotationFromMaps(annotationKeys[index]);
					annotationKeys = annotationCollection.getKeys();
					annotationUpdated = true;
					annotationPanelLowerPart();
				}
			});
			panel.add(removeButton);
			panel.setAlignmentY(Component.LEFT_ALIGNMENT); // TODO not working yet
			annotationPanelLowerPart.add(panel);
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
				c = str.toLowerCase().charAt(0);
				str = str.substring(1);
				annotationCollection.addAnnotation(c, str);
				addAnnotationToMaps(c);
				annotationKeys = annotationCollection.getKeys();
				annotationPanelLowerPart();
				annotationUpdated = true;
			}
		});

//		addAnnotationKey = new JCheckBox("Add annotation key");
//		addAnnotationKey.setBackground(Color.WHITE);
//		addAnnotationKey.addItemListener(new ItemListener() {
//			public void itemStateChanged(ItemEvent itemEvent) {
//				JCheckBox checkBox;
//				
//				checkBox = (JCheckBox)itemEvent.getSource();
//				
//		        if (checkBox.getModel().isSelected()) {
//		        	String[][] text = Editor.getInput();
//		        	for (int i=0; i<text.length; i++) {
//			        	annotationCollection.addAnnotation(text[i][0].charAt(0), text[i][1]);
//		        	}
//		        	annotationKeys = annotationCollection.getKeys();
//		        }
//			}
//		});
//		annotationPanelComment.add(addAnnotationKey);
		annotationPanelLowerPart.add(addAnnotationField);

		annotationPanelLowerPart.validate();
	}

	private void addAnnotationToMaps(char c) {
		scatPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke((int)(c+"").toUpperCase().charAt(0), 0), "annotation\t"+c);
		scatPanel.getActionMap().put("annotation\t"+c, new AnnotationAction(this, c));
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
	
	private void inputMapAndActionMap() {
		InputMap inputMap;
		ActionMap actionMap;
		
		inputMap = scatPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
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

		actionMap = scatPanel.getActionMap();
		actionMap.put(ALT_UP, new CycleRadio(typeRadioButtons, -1));
		actionMap.put(ALT_DOWN, new CycleRadio(typeRadioButtons, 1));
		actionMap.put(ALT_LEFT, new CycleRadio(colorKeyPanel.getClassRadioButtons(), -1));
		actionMap.put(ALT_RIGHT, new CycleRadio(colorKeyPanel.getClassRadioButtons(), 1));
		actionMap.put(FIRST, new AbstractAction() {
			public static final long serialVersionUID = 4L;

			public void actionPerformed(ActionEvent e) {
				markerIndex = 0;
				displayIndex(navigationField);
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);
				updateGUI();
			}
		});
		actionMap.put(PREVIOUS, new AbstractAction() {
			public static final long serialVersionUID = 5L;

			public void actionPerformed(ActionEvent e) {
				markerIndex = Math.max(markerIndex-1, 0);
				displayIndex(navigationField);
				scatPanel.setPointsGeneratable(true);
				scatPanel.setQcPanelUpdatable(true);
				updateGUI();
			}
		});
		actionMap.put(NEXT, new AbstractAction() {
			public static final long serialVersionUID = 6L;

			public void actionPerformed(ActionEvent e) {
				markerIndex = Math.min(markerIndex+1, markerList.length-1);
				displayIndex(navigationField);
				scatPanel.setPointsGeneratable(true);//zx
				scatPanel.setQcPanelUpdatable(true);//zx???
				updateGUI();
			}
		});
		actionMap.put(LAST, new AbstractAction() {
			public static final long serialVersionUID = 7L;

			public void actionPerformed(ActionEvent e) {
				markerIndex = markerList.length-1;
				displayIndex(navigationField);
				scatPanel.setPointsGeneratable(true);//zx
				scatPanel.setQcPanelUpdatable(true);//zx???
				updateGUI();
			}
		});
		scatPanel.setActionMap(actionMap);
	}

	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();
		String filename;
		int count;

		if (command.equals(FIRST)) {
//			markerIndex = 0;
			markerIndex = getAvailableMarker(false, false, showAllMarkersOrNot, showAnnotatedOrUnannotated);
			displayIndex(navigationField);
			scatPanel.setPointsGeneratable(true);//zx
			scatPanel.setQcPanelUpdatable(true);
//			currentClusterFilter=(byte) (clusterFilterCollection.getSize(getMarkerName())-1);
			setCurrentClusterFilter();
//			scatPanel.generateRectangles();
//			if (clusterFilterCollection.getSize(getMarkerName())>0) {
//				scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
//			}
			updateGUI();
			displayClusterFilterIndex();
			isInitilizing = true;
			updateAnnotationPanelAnnotationCheckBoxes();
		} else if (command.equals(PREVIOUS)) {
//			markerIndex = Math.max(markerIndex-1, 0);
			markerIndex = getAvailableMarker(false, true, showAllMarkersOrNot, showAnnotatedOrUnannotated);
//			if (markerData[markerIndex] == null) {
//			if (!loaded[markerIndex]) {
//				loadMarkerData(markerIndex);
//				loaded[markerIndex] = true;
////				scatPanel.updateMarkerData(markerData);
//			}
			displayIndex(navigationField);
			scatPanel.setPointsGeneratable(true);//zx
			scatPanel.setQcPanelUpdatable(true);
//			currentClusterFilter=(byte) (clusterFilterCollection.getSize(getMarkerName())-1);
			setCurrentClusterFilter();
//			scatPanel.generateRectangles();
//			if (clusterFilterCollection.getSize(getMarkerName())>0) {
//				scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
//			}
			updateGUI();
			displayClusterFilterIndex();
			isInitilizing = true;
			updateAnnotationPanelAnnotationCheckBoxes();
		} else if (command.equals(NEXT)) {
//			markerIndex = Math.min(markerIndex+1, markerList.length-1);
			markerIndex = getAvailableMarker(true, true, showAllMarkersOrNot, showAnnotatedOrUnannotated);
//			if (markerData[markerIndex] == null) {
//			if (!loaded[markerIndex]) {
//				loadMarkerData(markerIndex);
//				loaded[markerIndex] = true;
////				scatPanel.updateMarkerData(markerData);
//			}
			displayIndex(navigationField);
			scatPanel.setPointsGeneratable(true);//zx
			scatPanel.setQcPanelUpdatable(true);
//			currentClusterFilter=(byte) (clusterFilterCollection.getSize(getMarkerName())-1);
			setCurrentClusterFilter();
//			scatPanel.generateRectangles();
//			if (clusterFilterCollection.getSize(getMarkerName())>0) {
//				scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
//			}
			updateGUI();
			displayClusterFilterIndex();
			isInitilizing = true;
			updateAnnotationPanelAnnotationCheckBoxes();
		} else if (command.equals(LAST)) {
//			markerIndex = markerList.length-1;
			markerIndex = getAvailableMarker(true, false, showAllMarkersOrNot, showAnnotatedOrUnannotated);
//			if (markerData[markerIndex] == null) {
//			if (!loaded[markerIndex]) {
//				loadMarkerData(markerIndex);
//				loaded[markerIndex] = true;
////				scatPanel.updateMarkerData(markerData);
//			}
			displayIndex(navigationField);
			scatPanel.setPointsGeneratable(true);//zx
			scatPanel.setQcPanelUpdatable(true);
//			currentClusterFilter=(byte) (clusterFilterCollection.getSize(getMarkerName())-1);
			setCurrentClusterFilter();
//			scatPanel.generateRectangles();
//			if (clusterFilterCollection.getSize(getMarkerName())>0) {
//				scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
//			}
			updateGUI();
			displayClusterFilterIndex();
			isInitilizing = true;
			updateAnnotationPanelAnnotationCheckBoxes();
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
				saveClusterFilterCollection();
				//clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+" of "+clusterFilterCollection.getSize(getMarkerName()));
				displayClusterFilterIndex();
//				annotationPanelLowerPart();
				scatPanel.setPointsGeneratable(true);//zx
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
				filename = markerList[markerIndex]+"_"+MarkerData.TYPES[plot_type][0]+"-"+MarkerData.TYPES[plot_type][1]+(count==1?"":" v"+count);
				count++;
			} while (new File(proj.getProjectDir()+filename+".png").exists());
			scatPanel.screenCapture(proj.getProjectDir()+filename+".png");
		} else if (command.equals(DUMP)) {
			count = 1;
			do {
				filename = markerList[markerIndex]+"_dump"+(count==1?"":" v"+count);
				count++;
			} while (new File(proj.getProjectDir()+filename+".xln").exists());
//			markerData[markerIndex].writeToFile(samples, proj.getProjectDir()+filename+".xln");
			getCurrentMarkerData().dump(proj.getProjectDir()+filename+".xln", samples);
		} else if (command.equals(MASK_MISSING) || command.equals(UNMASK_MISSING)) {
			maskMissing = !maskMissing;
			((JButton)ae.getSource()).setText(maskMissing?UNMASK_MISSING:MASK_MISSING);
			scatPanel.setPointsGeneratable(true);//zx
			scatPanel.setQcPanelUpdatable(true);//zx
			updateGUI();
		} else {
			System.err.println("Error - unknown command '"+command+"'");
		}
	}

	private int getAvailableMarker(boolean forwardOrBackward, boolean firstOrLastInThatDirection, boolean allMarkersOrOnlyThoseAnnotatedOrUnannotated, boolean annotatedOrUnannotated) {
		int result;

		result = markerIndex;
		if (allMarkersOrOnlyThoseAnnotatedOrUnannotated) {
			if (firstOrLastInThatDirection) {
				if (forwardOrBackward) {
					result = Math.min(markerIndex + 1, markerList.length - 1);
				} else {
					result = Math.max(markerIndex - 1, 0);
				}
			} else {
				if (forwardOrBackward) {
					result = markerList.length - 1;
				} else {
					result = 0;
				}
			}
		} else {
			if (firstOrLastInThatDirection) {
				if (forwardOrBackward) {
					if (annotatedOrUnannotated) {
						for (int i = markerIndex + 1; i < isAnnotated.length; i++) {
							if (isAnnotated[i]) {
								result = i;
								break;
							}
						}
					} else {
						for (int i = markerIndex + 1; i < isAnnotated.length; i++) {
							if (! isAnnotated[i]) {
								result = i;
								break;
							}
						}
					}
				} else {
					if (annotatedOrUnannotated) {
						for (int i = markerIndex - 1; i >= 0; i--) {
							if (isAnnotated[i]) {
								result = i;
								break;
							}
						}
					} else {
						for (int i = markerIndex - 1; i >= 0; i--) {
							if (! isAnnotated[i]) {
								result = i;
								break;
							}
						}
					}
				}
			} else {
				if (forwardOrBackward) {
					if (annotatedOrUnannotated) {
						for (int i = isAnnotated.length - 1; i >= 0; i--) {
							if (isAnnotated[i]) {
								result = i;
								break;
							}
						}
					} else {
						for (int i = isAnnotated.length - 1; i >= 0; i--) {
							if (! isAnnotated[i]) {
								result = i;
								break;
							}
						}
					}
				} else {
					if (annotatedOrUnannotated) {
						for (int i = 0; i < isAnnotated.length; i++) {
							if (isAnnotated[i]) {
								result = i;
								break;
							}
						}
					} else {
						for (int i = 0; i < isAnnotated.length; i++) {
							if (! isAnnotated[i]) {
								result = i;
								break;
							}
						}
					}
				}
			}
		}

		return result;
	}

	private void loadClusterFilterFiles() {
		String[] otherClusterFilTerFiles;
		int choice;
		String[] options = new String[] {"Yes, load and delete old file", "No, delete old file", "Cancel and close ScatterPlot"};
		
		otherClusterFilTerFiles = Files.list(proj.getDir(Project.DATA_DIRECTORY), ".tempClusterFilters.ser", jar);
		if (otherClusterFilTerFiles.length > 0) {
			choice = JOptionPane.showOptionDialog(null, "Error - either multiple instances of ScatterPlot are running or ScatterPlot failed to close properly\n" +
														"last time. The ability to generate new ClusterFilters will be disabled until this file has been\n" +
														"removed. Do you want to delete it and load the contents of the temporary file into memory?",
												  "Error", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
			if (choice == 0) {
				// load the last one in otherClusterFilerFiles[]
				clusterFilterCollection = ClusterFilterCollection.load(proj.getDir(Project.DATA_DIRECTORY)+otherClusterFilTerFiles[otherClusterFilTerFiles.length-1], jar);
				for (int i=0; i<otherClusterFilTerFiles.length; i++) {
					(new File(proj.getDir(Project.DATA_DIRECTORY)+otherClusterFilTerFiles[i])).delete();
				}
			} else if (choice == 1) {
				// load permanent
				if (Files.exists(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, true), jar) ) {
					clusterFilterCollection = ClusterFilterCollection.load(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, true), jar);
				} else {
					clusterFilterCollection = new ClusterFilterCollection();
				}
				for (int i=0; i<otherClusterFilTerFiles.length; i++) {
					(new File(proj.getDir(Project.DATA_DIRECTORY)+otherClusterFilTerFiles[i])).delete();
				}
			} else {
				System.exit(1);
			}
		} else if (Files.exists(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, true), jar) ) {
			// load
			clusterFilterCollection = ClusterFilterCollection.load(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, true), jar);
		} else {
			// create new
			clusterFilterCollection = new ClusterFilterCollection();
		}
	}

	public void loadAnnotationCollection() {
		String[] annotatedMarkers;

		if (new File(proj.getFilename(Project.ANNOTATION_FILENAME, false, false)).exists()) {
			annotationCollection = (AnnotationCollection) Files.readSerial(proj.getFilename(Project.ANNOTATION_FILENAME, false, false));

		} else {
			annotationCollection = new AnnotationCollection();
			annotationCollection.addAnnotation('u', "Ugly");
			annotationCollection.addAnnotation('m', "Monomorphic");
		}

		annotatedMarkers = annotationCollection.getMarkerLists();
		for (int i = 0; i < isAnnotated.length; i++) {
			for (int j = 0; j < annotatedMarkers.length; j++) {
				if (markerList[i].equals(annotatedMarkers[j])) {
					isAnnotated[i] = true;
					annotated ++;
					break;
				}
			}
		}

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
//		System.out.print("Before: image=="+(scatPanel.image==null?"null":"BImg")+"\t");
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
//		System.out.println("After: image=="+(scatPanel.image==null?"null":"BImg"));

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
					centLabels[i].setText("LRR corr="+str+", err="+ext.formDeci(comp[1], 3));
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
		//System.out.println("filename: "+filename);//zx
		try {
			try {
				reader = Files.getReader(filename, jar, true, false);
				if (reader == null) {
					JOptionPane.showMessageDialog(null, "Failed to load '"+filename+"'", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				}
				while (reader.ready()) {
					line = reader.readLine().trim().split("\t", -1);
					if (markerLookup.contains(line[0])) {
						markerNames.add(line[0]);
						markerComments.add(line.length>1?line[1]:"");
					} else {
//						System.err.println("Error - could not find "+line[0]+" in the lookup table");
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
				System.err.println("Error reading file \""+filename+"\"");
				e.printStackTrace();
				System.exit(2);
			}

			masterMarkerList = Array.toStringArray(markerNames);
			masterCommentList = Array.toStringArray(markerComments);
		} catch (Exception e) {
			System.err.println("Error loading: "+filename);
			e.printStackTrace();
		}
	}
	
	public void loadMarkerDataFromList() {
		markerDataLoader = new MarkerDataLoader(proj, markerList, -1);
		thread2 = new Thread(markerDataLoader);
		thread2.start();

		markerIndex = 0;
		previousMarkerIndex = -1;
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
				System.err.println("Error - have been waiting on markerDataLoader to load "+markerList[markerIndex]+" for "+(count/4)+" secounds");
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
		//System.out.println("Found "+files.length+" .cent files in "+proj.getDir(Project.DATA_DIRECTORY));
		
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
			System.out.print("<");
			trav = Centroids.load(proj.getDir(Project.DATA_DIRECTORY)+files[i], jar);
			System.out.print(">");
			if (trav.getFingerprint() != set.getFingerprint()) {
				System.err.println("Error - Centroids file '"+files[i]+"' does not match up with the fingerprint of the current marker set");
				System.exit(1);
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
		System.out.println();
		cents = new float[v.size()][][][];
		for (int i = 0; i<cents.length; i++) {
			cents[i] = v.elementAt(i);
        }
		centList = Array.toStringArray(fileList);
	}

	public void updateGUI() {
		if (markerDataLoader == null) {
			return;
		}
		
//		if (markerData[markerIndex]==null) {
//			System.out.println("Marker data is null at index "+markerIndex);
////			create visual that we're loading data, (write to scatter panel "Loading data from file")
//			loadMarkerData(markerIndex);
//		}
		
		if (markerList.length==0) {
			markerName.setText("Error: marker data was not successfully loaded");
			commentLabel.setText("Check to make sure MarkerLookup is synchronized with the current data");
			classPanel.setEnabled(false);
		} else {
			markerName.setText(markerList[markerIndex]);
			commentLabel.setText(commentList[markerIndex]);

			if (classPanel!=null) {
				classPanel.setEnabled(getCurrentMarkerData().getFingerprint()==sampleListFingerprint);
			}
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
		
		scatPanel.paintAgain();
		if (markerIndex != previousMarkerIndex) {
			updateCentLabels();
		}
		previousMarkerIndex = markerIndex;
	}

	public void updateColorKey(Hashtable<String,String> hash) {
		colorKeyPanel.updateColorKey(hash);
	}

	public void updateQcPanel(byte chr, int[] genotype, String[] sex, String[] otherClass) {
		float callRate=0;//zx
		JLabel qcPanelLabel;//zx
		//JLabel qcCallRateLabel;//zxu
		//JLabel qcHwePvalueLabel;//zxu
		int[] alleleCounts;
		double hweP;
		//int[][] sexContingecyTable = new int[2][2];
		CTable classCount;
		String[] called;
		int currentClass;
		
		alleleCounts = new int[3];
		called = new String[genotype.length];
		if (chr==23) {
			for (int i=0; i<genotype.length; i++) {
				if (genotype[i]<=0) {
					called[i] = "-1";
					callRate++;//zx
				} else {
					called[i] = "1";
					if (sex[i].equals("2")) {
						alleleCounts[genotype[i]-1]++;
					}
				}
			}
		} else if (chr==24) {
			for (int i=0; i<genotype.length; i++) {
				if (genotype[i]<=0) {
					called[i] = "-1";
					callRate++;//zx
				} else {
					called[i] = "1";
				}
			}
		} else {
			for (int i=0; i<genotype.length; i++) {
				if (genotype[i]<=0) {
					called[i] = "-1";
					callRate++;//zx
				} else {
					called[i] = "1";
					alleleCounts[genotype[i]-1]++;
				}
				
			}
		}
		hweP = AlleleFreq.HWEsig(alleleCounts);
		callRate=(genotype.length-callRate)*100/genotype.length;//zx
		
		qcPanel.removeAll();
		qcPanel.repaint();
		
        qcPanelLabel = new JLabel("", JLabel.CENTER);//zx
        qcPanel.add(qcPanelLabel);//zx
        qcPanelLabel = new JLabel("QC Metrics", JLabel.CENTER);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 20));//zx
        qcPanel.add(qcPanelLabel);//zx
		
        qcPanelLabel = new JLabel("Callrate: "+callRate+"%"+"                           ", JLabel.LEFT);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
        qcPanel.add(qcPanelLabel);//zx

        qcPanelLabel = new JLabel("HWE p-value: "+ext.prettyP(hweP), JLabel.LEFT);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcPanelLabel);//zx

		ToolTipManager.sharedInstance().setDismissDelay(100000);

		classCount = new CTable(called, sex);//This is the problem.
		classCount.setCustomNullValues(Array.addStrToArray("-1", CTable.DEFAULT_NULL_VALUES));
		classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"}, {"1","Genotype NOT missing"}}, sampleData.getActualClassColorKey(0));
		qcPanelLabel = new JLabel("Callrate by sex: "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(), false), 1) ), JLabel.LEFT);//zx
		//classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"}, {"1","Genotype NOT missing"}}, Matrix.addRow(sampleData.getActualClassColorKey(0), new String[] {null, "missing"}));
		qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcPanelLabel);//zx

		classCount = new CTable(CTable.extrapolateCounts(sex, genotype));
		classCount.setCustomNullValues(Array.addStrToArray("-1", CTable.DEFAULT_NULL_VALUES));
		classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(0), new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}});
		qcPanelLabel = new JLabel("Allele Freq by sex: "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(), false), 1)), JLabel.LEFT);//zx
		//classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(0), new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}, {".","Missing"}});
		qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcPanelLabel);//zx

		currentClass = getCurrentClass();
		if (currentClass>=(SampleData.BASIC_CLASSES.length+1)) {
			classCount = new CTable(called, otherClass);//This is the problem.
			classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"}, {"1","Genotype NOT missing"}}, sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length));
			qcPanelLabel = new JLabel("Callrate by "+sampleData.getClassName(currentClass)+": "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(), false), 1) ), JLabel.LEFT);//zx
			classCount.setCustomNullValues(Array.addStrToArray("-1", Array.addStrToArray("0", CTable.DEFAULT_NULL_VALUES)));
			//classCount.setCustomLabelsAndOrder(new String[][] {{"-1","Genotype missing"}, {"1","Genotype NOT missing"}}, Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length), new String[] {null, "missing"}));
			qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
	        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
			qcPanel.add(qcPanelLabel);//zx
	
			classCount = new CTable(CTable.extrapolateCounts(otherClass, genotype));
			classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length), new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}});
			//classCount.replaceIdWithLabel(SampleData.KEYS_FOR_BASIC_CLASSES[1],sampleData.getActualClassColorKey(0));
			qcPanelLabel = new JLabel("Allele Freq by "+sampleData.getClassName(currentClass)+": "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTable(), false), 1)), JLabel.LEFT);//zx
			classCount.setCustomNullValues(Array.addStrToArray("-1", Array.addStrToArray("0", CTable.DEFAULT_NULL_VALUES)));
			//classCount.setCustomLabelsAndOrder(Matrix.addRow(sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length), new String[] {null, "missing"}), new String[][] {{"A","Allele A"}, {"B","Allele B"}});
			qcPanelLabel.setToolTipText(classCount.getCTableInHtml());
	        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
			qcPanel.add(qcPanelLabel);//zx
		}

//		qcPanelLabel = new JLabel("Minor Allele Freq: " + (new DecimalFormat("#.####").format(classCount.getMinorAlleleFrequency())), JLabel.LEFT);//zx
//      qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
//		qcPanel.add(qcPanelLabel);//zx

//		AlleleFreq.calcFrequency(genotypes);
		
		
		/*
		qcPanelLabel = new JLabel(""+sexContingecyTable[0][1], JLabel.LEFT);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcPanelLabel);//zx
		qcPanelLabel = new JLabel(""+sexContingecyTable[1][0], JLabel.LEFT);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcPanelLabel);//zx
		qcPanelLabel = new JLabel(""+sexContingecyTable[1][1], JLabel.LEFT);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcPanelLabel);//zx
		*/

		qcPanel.validate();
	}


//	public void updateAnnotationPanel() {
//		JPanel panel;
//		JButton removeButton;
//		BoxLayout boxLayout;
//		
//		annotationPanelComment.removeAll();
//		annotationPanelComment.repaint();
//
////		annotationPanel = new JPanel();
//		boxLayout = new BoxLayout(annotationPanelComment, BoxLayout.Y_AXIS);
//		annotationPanelComment.setLayout(boxLayout);
//		annotationPanelComment.setBackground(BACKGROUND_COLOR);
//
//		annotationCheckBoxes = new JCheckBox[annotationKeys.length];
//		for (int i=0; annotationKeys != null && i < annotationKeys.length; i++) {
//			annotationCheckBoxes[i] = new JCheckBox(annotationCollection.getDescriptionForComment(annotationKeys[i], showAnnotationShortcuts));
//			annotationCheckBoxes[i].setBackground(Color.WHITE);
//			if (annotationCollection.markerHasAnnotation(markerList[markerIndex], annotationKeys[i])) {
//				annotationCheckBoxes[i].setSelected(true);
//			}
//			annotationCheckBoxes[i].addItemListener(new ItemListener() {
//				public void itemStateChanged(ItemEvent itemEvent) {
//					JCheckBox checkBox;
//					char currentKey;
//					
//					checkBox = (JCheckBox)itemEvent.getSource();
//					currentKey = checkBox.getText().charAt(1);
//					
//			        if (checkBox.getModel().isSelected()) {
//			        	annotationCollection.addAnnotationForMarker(markerList[markerIndex], currentKey);
//			        	isAnnotated[markerIndex] = true;
//			        } else {
//			        	annotationCollection.removeAnnotationForMarker(markerList[markerIndex], currentKey);
//			        	isAnnotated[markerIndex] = (annotationCollection.annotationsForMarker(markerList[markerIndex]).length != 0);
//			        }
//			        checkBox.setText(annotationCollection.getDescriptionForComment(currentKey, showAnnotationShortcuts));
//			        annotationUpdated = true;
//				}
//			});
//			panel = new JPanel();
//			panel.setBackground(Color.WHITE);
////			panel.setLayout(new BoxLayout(panel, BoxLayout.X_AXIS));
//			panel.add(annotationCheckBoxes[i]);
//			removeButton = new JButton(Grafik.getImageIcon("images/delete2sm.png", true));
//			removeButton.setActionCommand("removeButton" + i);
//			removeButton.setBorder(null);
//			removeButton.setSize(6, 6);
//			removeButton.addActionListener(new ActionListener() {
//				public void actionPerformed(ActionEvent e) {
//					String commandText;
//					int index;
//
//					commandText = e.getActionCommand();
//					index = Integer.parseInt(commandText.substring("removeButton".length(), commandText.length()));
//					annotationCollection.removeAnnotation(annotationKeys[index]);
//					removeAnnotationFromMaps(annotationKeys[index]);
//					annotationKeys = annotationCollection.getKeys();
//					annotationUpdated = true;
//				}
//			});
//			panel.add(removeButton);
//			panel.setAlignmentY(Component.LEFT_ALIGNMENT); // TODO not working yet
//			annotationPanelComment.add(panel);
//		}
//		
//		addAnnotationField = new JTextField("default");
//		
//		addAnnotationField.addFocusListener(new FocusListener() {
//			public void focusLost(FocusEvent e) {
//				if (showAnnotationShortcuts) {
//					activateAllAnnotationMaps();
//				}
//			}
//			public void focusGained(FocusEvent e) {
//				deactivateAllAnnotationMaps();
//			}
//		});
//		
//		addAnnotationField.addActionListener(new ActionListener() {
//			public void actionPerformed(ActionEvent e) {
//				String str;
//				char c;
//				
//				str = ((JTextField)e.getSource()).getText();
//				c = str.toLowerCase().charAt(0);
//				str = str.substring(1);
//				annotationCollection.addAnnotation(c, str);
//				addAnnotationToMaps(c);
//				annotationKeys = annotationCollection.getKeys();
//				updateAnnotationPanel();
//				annotationUpdated = true;
//			}
//		});
//
////		addAnnotationKey = new JCheckBox("Add annotation key");
////		addAnnotationKey.setBackground(Color.WHITE);
////		addAnnotationKey.addItemListener(new ItemListener() {
////			public void itemStateChanged(ItemEvent itemEvent) {
////				JCheckBox checkBox;
////				
////				checkBox = (JCheckBox)itemEvent.getSource();
////				
////		        if (checkBox.getModel().isSelected()) {
////		        	String[][] text = Editor.getInput();
////		        	for (int i=0; i<text.length; i++) {
////			        	annotationCollection.addAnnotation(text[i][0].charAt(0), text[i][1]);
////		        	}
////		        	annotationKeys = annotationCollection.getKeys();
////		        }
////			}
////		});
////		annotationPanelComment.add(addAnnotationKey);
//		annotationPanelComment.add(addAnnotationField);
//
//		annotationPanelComment.validate();
//	}


	public void updateAnnotationPanel() {
		updateAnnotationPanelFilterRadioButtons();
		updateAnnotationPanelAnnotationCheckBoxes();
	}

	public void updateAnnotationPanelFilterRadioButtons() {
        for (int i=0; i<fileterRadioButtons.length; i++) {
        	fileterRadioButtons[i].setText(RADIOBUTTON_TEXTS[i] + " (n=" + (i==0? isAnnotated.length : (i==1? annotated : (isAnnotated.length - annotated))) + ")");
        }
	}

	public void updateAnnotationPanelAnnotationCheckBoxes() {
		for (int i=0; i<annotationCheckBoxes.length; i++) {
        	annotationCheckBoxes[i].setText(annotationCollection.getDescriptionForComment(annotationKeys[i], showAnnotationShortcuts));
			if (annotationCollection.markerHasAnnotation(markerList[markerIndex], annotationKeys[i])) {
				annotationCheckBoxes[i].setSelected(true);
			} else {
				annotationCheckBoxes[i].setSelected(false);
			}
        }
		isInitilizing = false;
	}

	public void toggleAnnotationBox(char c) {
		int index;
		
		index = ext.indexOfChar(c, annotationKeys);
		annotationCheckBoxes[index].setSelected(!annotationCheckBoxes[index].isSelected());
//		isAnnotated[markerIndex] = true;
//		annotationUpdated = true;
		if (annotationAutoAdv) {
			actionPerformed(new ActionEvent(next, 0, NEXT));
		}
	}

	public void displayIndex(JTextField field) {
		field.setText((markerIndex+1)+" of "+markerList.length);
	}
	
	public void displayClusterFilterIndex() {
		clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))
										+" of "
										+clusterFilterCollection.getSize(getMarkerName()));
	}

	public void saveClusterFilterCollection() {
		clusterFilterCollectionUpdated(true);
		autoSaveCFC.saveNow();
	}
	
	public void clusterFilterCollectionUpdated(boolean status) {
		if (status && autoSaveCFC == null) {
			autoSaveCFC = new AutoSaveClusterFilterCollection(clusterFilterCollection, proj.getDir(Project.DATA_DIRECTORY)+sessionID+".tempClusterFilters.ser", 30);
			new Thread(autoSaveCFC).start();
		}
		clusterFilterCollectionUpdated = status;
	}

	public static void main(String[] args) {
		String filename = Project.DEFAULT_SCATTER_PROJECT;
//		boolean jar = args.length>0&&args[0].equals("-notJar")?false:true;
		boolean jar = false;

		try {
			new ScatterPlot(new Project(filename, jar), null, null);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void windowActivated(WindowEvent e) {}

	public void windowClosed(WindowEvent e) {}

	public void windowClosing(WindowEvent e) {
		String[] options;
		int choice;
		String filename;
		
		options = new String[] {"Yes, overwrite", "No"};
		if (clusterFilterCollectionUpdated) {
			choice = JOptionPane.showOptionDialog(null, "New ClusterFilters have been generated. Do you want to save them to the permanent file?", "Overwrite permanent file?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
			if (choice == 0) {
				clusterFilterCollection.serialize(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, false));
			}
			autoSaveCFC.kill();
		}

		if (annotationUpdated) {
			choice = JOptionPane.showOptionDialog(null, "New Annotations have been generated. Do you want to save them to the permanent file?", "Overwrite permanent file?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
			if (choice == 0) {
				filename = proj.getFilename(Project.ANNOTATION_FILENAME, false, false);
				System.out.println("Writing to "+filename);
				Files.writeSerial(annotationCollection, filename);
			}
		}

		//TODO notify all threads (e.g., MarkerDataLoader) that they need to close
		markerDataLoader.kill();
		
	}

	public void windowDeactivated(WindowEvent e) {}

	public void windowDeiconified(WindowEvent e) {}

	public void windowIconified(WindowEvent e) {}

	public void windowOpened(WindowEvent e) {}

//	@Override
//	public void keyTyped(KeyEvent e) {
//		String shortcut;
//		String temp;
//
//		System.err.println("key typed");
//		
//		shortcut = e.getKeyChar() + "";
//		for (int i=0; i<commentShortcuts.length; i++) {
//			if (commentShortcuts[i].equals(shortcut)) {
//				temp = annotationHash.get(getMarkerName());
//				if (temp == null || ! temp.contains((comments[i]))) {
//					annotationHash.put(getMarkerName(), (temp == null? "" : temp + "\t") + comments[i]);
////					updateAnnotationPanel();
//					actionPerformed(new ActionEvent(next, 0, NEXT));
//					System.out.println(annotationHash.get(getMarkerName()));
//				}
//			}
//		}
//	}

//	@Override
//	public void keyPressed(KeyEvent e) {
//	}

//	@Override
//	public void keyReleased(KeyEvent e) {
//	}

}
