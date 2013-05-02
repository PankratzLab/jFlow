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
import cnv.gui.AutoSaveClusterFilterCollection;
import cnv.gui.CycleRadio;
import cnv.manage.MarkerDataLoader;
import common.*;
import cnv.var.*;

// TODO needs major cleanup
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
	public static final String MASK_MISSING = "Mask missing values";
	public static final String UNMASK_MISSING = "Unmask missing values";
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	private JButton first, previous, next, last;
	private JTextField navigationField;
	private JPanel classPanel;
//	private JPanel legendPanel;
//	private JPanel bottomPanel;
	private ScatterPanel scatPanel;
	private JLabel sizeLabel;
	private JLabel gcLabel;
//	private JPanel typePanel;
	private JTabbedPane tabbedPane;
	private JPanel qcPanel;//zx
	//private JLabel qcPanelLabel;//zx
	//private JLabel qcCallRateLabel;//zxu
	//private JLabel qcHwePvalueLabel;//zxu
	//private JPael typePanelRadioButton;//
	private JComboBox<String> newGenotype;

	private Project proj;
//	private MarkerData[] markerData;
	private String[] markerList;
	private float[][][][] cents;
	private String[] centList;
	private JCheckBox[] centBoxes;
	private boolean[] displayCents;
	private JLabel[] centLabels;
	private String[] commentList;
	private int markerIndex;
	private int previousMarkerIndex;
	//private JLabel markerName, commentLabel;
	private JTextField markerName, commentLabel;
	private String[] samples;
//	private byte currentClass;
	private int plot_type;
	private byte size;
	private float gcThreshold;
//	private Hashtable<String,IndiPheno> sampleData;
	private long sampleListFingerprint;
	private MarkerLookup markerLookup;
	private boolean jar;
	private SampleData sampleData;
	private JCheckBox symmetryBox;
	private boolean maskMissing;
	private ClusterFilterCollection clusterFilterCollection;
	private byte currentClusterFilter;//zx
	private JTextField clusterFilterNavigation;//zx
	private boolean clusterFilterCollectionUpdated;
	private String sessionID;
	private AutoSaveClusterFilterCollection autoSaveCFC;
	private JRadioButton[] typeRadioButtons;
	private JRadioButton[] classRadioButtons;
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
		
		markerList = initMarkerList;
		commentList = initCommentList;
		if (markerList == null) {
			loadMarkerListFromFile();
			if (markerList == null) {
				return;
			}
		}
		loadMarkerDataFromList();
		if (commentList == null) {
			commentList = Array.stringArray(markerList.length, "");
		}
		
		System.err.println("3\t"+ext.getTimeElapsed(time));
		loadCentroids();
		addWindowListener(this);
		sessionID = (new Date().getTime()+"").substring(5);
		loadClusterFilterFiles();
		autoSaveCFC = null;
		System.err.println("4\t"+ext.getTimeElapsed(time));
		
		scatPanel = new ScatterPanel(this);
		colorScheme = scatPanel.getColorScheme();
		getContentPane().add(scatPanel, BorderLayout.CENTER);
		getContentPane().add(markerPanel(), BorderLayout.NORTH);
		getContentPane().add(controlPanel(), BorderLayout.EAST);
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
		clusterFilterCollectionUpdated = false;
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

	private JComponent markerPanel() {
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
		navigationField.addFocusListener(new FocusListener() {
			public void focusGained(FocusEvent focusevent) {}

			public void focusLost(FocusEvent fe) {
				try {
					int trav = Integer.valueOf(((JTextField)fe.getSource()).getText().split("[\\s]+")[0]).intValue()-1;
					if (trav>=0&&trav<markerList.length) {
						markerIndex = trav;
//						scatPanel.setPointsGenerated(false);//zx
//						scatPanel.setUpdateQcPanel(true);//zx???
//						updateGUI();
					}
				} catch (NumberFormatException nfe) {}
				displayIndex((JTextField)fe.getSource());
				scatPanel.setPointsGeneratable(true);//zx
				scatPanel.setQcPanelUpdatable(true);//zx???
//				currentClusterFilter=(byte) (clusterFilterCollection.getSize(getMarkerName())-1);
				setCurrentClusterFilter();
//				scatPanel.generateRectangles();
//				if (clusterFilterCollection.getSize(getMarkerName())>0) {
//					scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
//				}
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

	private JComponent controlPanel() {
		JPanel controlPanel = new JPanel();
		controlPanel.setLayout(new BoxLayout(controlPanel, BoxLayout.PAGE_AXIS));
		//typePanel.setLayout(new GridLayout(20, 1));
//		typePanel.setLayout(new GridBagLayout());
		controlPanel.setBackground(BACKGROUND_COLOR);
//		typePanel.setBorder();
//		typePanel.setSize(50, 100);
		
        GridBagConstraints gbc = new GridBagConstraints();   
        gbc.insets = new Insets(1,3,0,30);   
        gbc.weightx = 1.0;   
        gbc.fill = GridBagConstraints.HORIZONTAL;   
        gbc.gridwidth = GridBagConstraints.REMAINDER;   

		JPanel tabPanel;
		tabbedPane = new JTabbedPane();
		tabbedPane.setBackground(BACKGROUND_COLOR);
		tabPanel = new JPanel();
//		tabPanel.setLayout(new BoxLayout(tabPanel, BoxLayout.Y_AXIS));
		tabPanel.setLayout(new GridBagLayout());
		tabPanel.setSize(50, 100);
		tabPanel.setBackground(BACKGROUND_COLOR);

		tabPanel.add(sizeSliderPanel(), gbc);
		tabPanel.add(gcSliderPanel(), gbc);

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

		tabPanel.add(symmetryBox, gbc);
//		tabPanel.add(symmetryBox);
		
		JButton button = new JButton(CAPTURE);
		button.addActionListener(this);
		button.setActionCommand(CAPTURE);
		tabPanel.add(button, gbc);
//		tabPanel.add(button);

		button = new JButton(DUMP);
		button.addActionListener(this);
		button.setActionCommand(DUMP);
		tabPanel.add(button, gbc);
//		tabPanel.add(button);
		
		button = new JButton(MASK_MISSING);
		button.addActionListener(this);
		button.setActionCommand(MASK_MISSING);
		tabPanel.add(button, gbc);
//		tabPanel.add(button);

//		typePanel.addTab("Control", null, tabPanel, "Manipulate the chart");
//		typePanel.setMnemonicAt(0, KeyEvent.VK_1);
//
//		tabPanel = new JPanel();

		tabPanel.add(clusterFilterPanel(), gbc);
//		tabPanel.add(clusterFilterPanel());

		tabbedPane.addTab("Control", null, tabPanel, "Manipulate the chart");
		tabbedPane.setMnemonicAt(0, KeyEvent.VK_1);

		tabPanel = new JPanel();
		
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
		tabPanel.add(label, gbc);
//		tabPanel.add(label);
		label = new JLabel("Centroids:");
		label.setFont(new Font("Arial", 0, 20));
		label.setHorizontalAlignment(JLabel.CENTER);
		tabPanel.add(label, gbc);
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
			tabPanel.add(centBoxes[i], gbc);
//			tabPanel.add(centBoxes[i]);
			
			centLabels[i] = new JLabel("LRR correlation not performed");
			centLabels[i].setVisible(displayCents[i]);
			tabPanel.add(centLabels[i], gbc);
//			tabPanel.add(centLabels[i]);
		}
		
        //JLabel padding = new JLabel();//np
        //gbc.weighty = 1.0;
        //typePanel.add(padding, gbc);
		//private JLabel qcPanelLabel;//zx
		//private JLabel qcCallRateLabel;//zxu
		//private JLabel qcHwePvalueLabel;//zxu

		tabbedPane.addTab("Centroid", null, tabPanel, "Displays the centroid");
		tabbedPane.setMnemonicAt(1, KeyEvent.VK_2);

//		tabPanel = new JPanel();

		qcPanel = new JPanel();//zx
		qcPanel.setLayout(new GridLayout(8, 1));//zx
		qcPanel.setBackground(BACKGROUND_COLOR);//zx
		/*
        qcPanelLabel = new JLabel("                      ", JLabel.CENTER);//zx
        qcPanel.add(qcPanelLabel);//zx
        qcPanelLabel = new JLabel("QC Metrics", JLabel.CENTER);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 20));//zx
        qcPanel.add(qcPanelLabel);//zx
		
        qcCallRateLabel = new JLabel("Call Rate: "+ScatterPanel.getCallRate()+"%", JLabel.LEFT);//zx
        qcCallRateLabel.setFont(new Font("Arial", 0, 14));//zx
        qcPanel.add(qcCallRateLabel);//zx

        qcHwePvalueLabel = new JLabel("HWE p-value: ", JLabel.LEFT);//zx
        qcHwePvalueLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcHwePvalueLabel);//zx
		*/

//		tabbedPane.addTab("QC", null, tabPanel, "Displays the QC result");
//		tabbedPane.setMnemonicAt(2, KeyEvent.VK_3);

		controlPanel.add(tabbedPane);

		controlPanel.add(plotTypePanel());
		controlPanel.add(qcPanel);//zx
        return controlPanel;
    }

	private JComponent plotTypePanel() {
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
							scatPanel.setPointsGeneratable(true);//zx
							scatPanel.setQcPanelUpdatable(true);//zx???
//							scatPanel.generateRectangles();//zx???
//							if (clusterFilterCollection.getSize(getMarkerName())>0) {
//								scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
//							}
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
		// --- End of the original block ---
		
		/*
		// --- Beginning of the new block ---
		typePanelRadioButton = new JPanel();//
		typePanelRadioButton.setLayout(new GridLayout(MarkerData.TYPES.length,2));//
		ButtonGroup typeRadio = new ButtonGroup();
		JRadioButton[] typeRadioButtons = new JRadioButton[MarkerData.TYPES.length];
		JTextField[] typeRadioTexts = new JTextField[MarkerData.TYPES.length];//
		for (int i = 0; i<MarkerData.TYPES.length; i++) {
			//typeRadioButtons[i] = new JRadioButton(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1], false);
			typeRadioButtons[i] = new JRadioButton(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1], false);//
			//typeRadioButtons[i].setFont(new Font("Arial", 0, 14));
			typeRadio.add(typeRadioButtons[i]);
			typeRadioButtons[i].addItemListener(typeListener);
			typeRadioButtons[i].setBackground(BACKGROUND_COLOR);
			//typePanel.add(typeRadioButtons[i], gbc);
			typePanelRadioButton.add(typeRadioButtons[i]);//
			typeRadioTexts[i]=new JTextField(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1]);//
			typeRadioTexts[i].setBackground(BACKGROUND_COLOR);//
			typeRadioTexts[i].setFont(new Font("Arial",0,14));//
			typeRadioTexts[i].setBorder(null);//
			typePanelRadioButton.add(typeRadioTexts[i]);//
		}
		typePanel.add(typePanelRadioButton, gbc);//
		// --- End of the new block ---
		 */
	}

	private JComponent sizeSliderPanel() {
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

	private JComponent gcSliderPanel() {
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

//	private JComponent colorPanel() {
//		JPanel bottomPanel = new JPanel();
//		bottomPanel.setLayout(new GridLayout(2, 1));
//        //bottomPanel.setBackground(BACKGROUND_COLOR);
//		classPanel = new JPanel();
//		JLabel label = new JLabel("Color code by:");
//		label.setFont(new Font("Arial", 0, 14));
//		classPanel.add(label);
//
//		ItemListener classListener = new ItemListener() {
//			public void itemStateChanged(ItemEvent ie) {
//				JRadioButton jrb = (JRadioButton)ie.getItem();
//				if (jrb.isSelected()) {
//					for (int i = 0; i<sampleData.getNumClasses(); i++) {
//						if (jrb.getText().equals(sampleData.getClassName(i))) {
//							currentClass = i;
////							scatPanel.setPointsGenerated(true);//zx Why should be false?
//							scatPanel.setPointsGeneratable(true);//zx
//							scatPanel.setQcPanelUpdatable(true);//zx
//							updateGUI();
//						}
//					}
//				}
//			}
//		};
//		ButtonGroup classRadio = new ButtonGroup();
//		classRadioButtons = new JRadioButton[sampleData.getNumClasses()];
//		for (int i = 0; i<sampleData.getNumClasses(); i++) {
//			classRadioButtons[i] = new JRadioButton(sampleData.getClassName(i), false);
//			classRadioButtons[i].setFont(new Font("Arial", 0, 14));
//			classRadio.add(classRadioButtons[i]);
//			classRadioButtons[i].addItemListener(classListener);
//			classRadioButtons[i].setBackground(BACKGROUND_COLOR);
//			classPanel.add(classRadioButtons[i]);
//		}
//		classPanel.setBackground(BACKGROUND_COLOR);
//		bottomPanel.add(classPanel);
//		
//		legendPanel = new JPanel();
//        legendPanel.setBackground(BACKGROUND_COLOR);
//
//        //JLabel legend1 = new JLabel("Color Key: ");
//		//legend1.setFont(new Font("Arial", 0, 14));
//		//legendPanel.add(legend1);
//		
////		for (int i=0; i<sampleData.getActualClassColorKey(currentClass).length; i++){
////			legendPanel.add(new JLabel(new ColorIcon(12,12,ScatterPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[i][0])])));
////			legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[i][1]));
////		}
////		legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[0][0])])));
////		legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[0][1]));
////		legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[1][0])])));
////		legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[1][1]));
//		
//		bottomPanel.add(legendPanel);
//		classRadioButtons[1].setSelected(true);
//
//		return bottomPanel;
//	}

	private JComponent clusterFilterPanel() {
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
				displayClusterFilterIndex();
				scatPanel.setPointsGeneratable(true);//zx
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
				scatPanel.setPointsGeneratable(true);//zx
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

	private void inputMapAndActionMap() {
		InputMap inputMap = scatPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
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
		ActionMap actionMap = scatPanel.getActionMap();
		actionMap.put(ALT_UP, new CycleRadio(typeRadioButtons, -1));
		actionMap.put(ALT_DOWN, new CycleRadio(typeRadioButtons, 1));
		actionMap.put(ALT_LEFT, new CycleRadio(classRadioButtons, -1));
		actionMap.put(ALT_RIGHT, new CycleRadio(classRadioButtons, 1));
		actionMap.put(FIRST, new AbstractAction() {
			public static final long serialVersionUID = 4L;

			public void actionPerformed(ActionEvent e) {
				markerIndex = 0;
				displayIndex(navigationField);
				scatPanel.setPointsGeneratable(true);//zx
				scatPanel.setQcPanelUpdatable(true);//zx???
				updateGUI();
			}
		});
		actionMap.put(PREVIOUS, new AbstractAction() {
			public static final long serialVersionUID = 5L;

			public void actionPerformed(ActionEvent e) {
				markerIndex = Math.max(markerIndex-1, 0);
				displayIndex(navigationField);
				scatPanel.setPointsGeneratable(true);//zx
				scatPanel.setQcPanelUpdatable(true);//zx???
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
			markerIndex = 0;
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
		} else if (command.equals(PREVIOUS)) {
			markerIndex = Math.max(markerIndex-1, 0);
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
		} else if (command.equals(NEXT)) {
			markerIndex = Math.min(markerIndex+1, markerList.length-1);
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
		} else if (command.equals(LAST)) {
			markerIndex = markerList.length-1;
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
		} else if (command.equals(CLUSTER_FILTER_BACKWARD)) {
			if (clusterFilterCollection.getSize(getMarkerName())>0) {
				scatPanel.rectangles[currentClusterFilter].setColor((byte)7);
				//scatPanel.setRectangles();
				setCurrentClusterFilter((byte) Math.max(currentClusterFilter-1, 0));
//				scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
				//clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+" of "+clusterFilterCollection.getSize(getMarkerName()));
				displayClusterFilterIndex();
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

			markerList = Array.toStringArray(markerNames);
			commentList = Array.toStringArray(markerComments);
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
	
	
//	public synchronized void loadMarkerData(int markerIndex) {
////		thread2.suspend();
//		markerDataLoader.requestIndexBeTheNextFileToLoad(markerIndex);
////		markerData[markerIndex] = MarkerSet.loadFromList(proj, new String[] {markerList[markerIndex]})[0];
////		loaded[markerIndex] = true;
////		thread2.resume();
//
////		try {
//////			Thread.sleep(5000);
//////			thread2.sleep(5000); //TODO
//////			thread2.wait(5000);
//////			thread2.interrupt();
////			thread2.wait();
////			markerData[markerIndex] = MarkerSet.loadFromList(proj, new String[] {markerList[markerIndex]})[0];
////			loaded[markerIndex] = true;
////			thread2.notify();
////		} catch (InterruptedException e) {
////			// TODO Auto-generated catch block
////			e.printStackTrace();
////		}
//	}

//	public void loadMarkerData(int markerIndex) {
////		thread2.suspend();
////		markerData[markerIndex] = MarkerSet.loadFromList(proj, new String[] {markerList[markerIndex]})[0];
////		loaded[markerIndex] = true;
////		thread2.resume();
//
//		try {
////			Thread.sleep(5000);
//			thread2.sleep(5000); //TODO
////			thread2.wait(5000);
////			thread2.interrupt();
//			markerData[markerIndex] = MarkerSet.loadFromList(proj, new String[] {markerList[markerIndex]})[0];
//			loaded[markerIndex] = true;
//		} catch (InterruptedException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//	}
//
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


//	public void updateColorKey(Hashtable<String,String> hash) {
//		JLabel label, block;
//		String[][] colorKeys;
//		String[] keys;
//		legendPanel.removeAll();
//		legendPanel.repaint();
//		//JLabel legend = new JLabel("Color Key: ");
//		//legend.setFont(new Font("Arial", 0, 14));
//		//legendPanel.add(legend);
//		
//		//legendPanel.removeAll();
//		//legendPanel.repaint();
//		label = new JLabel("Color key:");
//		label.setFont(new Font("Arial", 0, 14));
//		legendPanel.add(label);
//		if (currentClass < SampleData.BASIC_CLASSES.length) {
//			colorKeys = SampleData.KEYS_FOR_BASIC_CLASSES[currentClass];
//		} else {		
//			colorKeys = sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length);
//		}
//		for (int i = 0; i<colorKeys.length; i++) {
//			block = new JLabel(new ColorIcon(12, 12, ScatterPanel.DEFAULT_COLORS[Integer.parseInt(colorKeys[i][0])]));
//			label = new JLabel(colorKeys[i][1]+" (n="+(hash.containsKey(colorKeys[i][0])?hash.get(colorKeys[i][0]):"0")+")");
//			hash.remove(colorKeys[i][0]);
//			label.setFont(new Font("Arial", 0, 14));
//			legendPanel.add(block);
//			legendPanel.add(label);
//		}
//		keys = HashVec.getKeys(hash);
//		for (int i = 0; i<keys.length; i++) {
//			if (!keys[i].equals("-1")) {
//				block = new JLabel(new ColorIcon(12, 12, ScatterPanel.DEFAULT_COLORS[Integer.parseInt(keys[i])]));
//				label = new JLabel((keys[i].equals("0")?"missing":keys[i])+" (n="+hash.get(keys[i])+")");
//				label.setFont(new Font("Arial", 0, 14));
//				legendPanel.add(block);
//				legendPanel.add(label);
//			}
//		}
//
//		
///**		
//		//colorKeys = sampleData.getActualClassColorKey(currentClass);
//		for (int i=0; i<sampleData.getActualClassColorKey(currentClass).length; i++){
//			legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[i][0])])));
//			legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[i][1]));
//			hash.remove(arg0)
//		}
////		for (int i = 0; i<colorKeys.length; i++) {
////			label = new JLabel(colorKeys[i][1]+" (n="+(hash.containsKey(colorKeys[i][0])?hash.get(colorKeys[i][0]):"0")+")");
////		}
//		/*
//		keys = HashVec.getKeys(hash);
//		for (int i = 0; i<keys.length; i++) {
//			if (!keys[i].equals("-1")) {
//				block = new JLabel(new ColorIcon(12, 12, ScatterPanel.DEFAULT_COLORS[Integer.parseInt(keys[i])]));
//				label = new JLabel((keys[i].equals("0")?"missing":keys[i])+" (n="+hash.get(keys[i])+")");
//				label.setFont(new Font("Arial", 0, 14));
//				legendPanel.add(block);
//				legendPanel.add(label);
//			}
//		}
//		*/
////*/
//
//		legendPanel.validate();
//		
//		/*
//		legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[0][0])])));
//		legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[0][1]));
//		legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[1][0])])));
//		legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[1][1]));
//		*/
//	}


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
		
		options = new String[] {"Yes, overwrite", "No"};
		if (clusterFilterCollectionUpdated) {
			choice = JOptionPane.showOptionDialog(null, "New ClusterFilters have been generated. Do you want to save them to the permanent file?", "Overwrite permanent file?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
			if (choice == 0) {
				clusterFilterCollection.serialize(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, false));
			}
			autoSaveCFC.kill();
		} else {
		}

		//TODO notify all threads (e.g., MarkerDataLoader) that they need to close
		markerDataLoader.kill();
		
	}

	public void windowDeactivated(WindowEvent e) {}

	public void windowDeiconified(WindowEvent e) {}

	public void windowIconified(WindowEvent e) {}

	public void windowOpened(WindowEvent e) {}

//	private void generateTypePanel () {
//		typePanel = new JPanel();
//		//typePanel.setLayout(new BoxLayout(typePanel, BoxLayout.PAGE_AXIS));
//		//typePanel.setLayout(new GridLayout(20, 1));
//		typePanel.setLayout(new GridBagLayout());
//		
//        GridBagConstraints gbc = new GridBagConstraints();   
//        gbc.insets = new Insets(1,3,0,30);   
//        gbc.weightx = 1.0;   
//        gbc.fill = GridBagConstraints.HORIZONTAL;   
//        gbc.gridwidth = GridBagConstraints.REMAINDER;   
//
//		ItemListener typeListener = new ItemListener() {
//			public void itemStateChanged(ItemEvent ie) {
//				JRadioButton jrb = (JRadioButton)ie.getItem();
//				if (jrb.isSelected()) {
//					for (int i = 0; i<MarkerData.TYPES.length; i++) {
//						if (jrb.getText().equals(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1])) {
//							plot_type = i;
//							scatPanel.setPointsGenerated(false);//zx
//							scatPanel.setUpdateQcPanel(true);//zx???
////							scatPanel.generateRectangles();//zx???
////							if (clusterFilterCollection.getSize(getMarkerName())>0) {
////								scatPanel.rectangles[currentClusterFilter].setColor((byte)0);
////							}
//							updateGUI();
//						}
//					}
//				}
//			}
//		};
//		
//		// --- Beginning of the original block ---
//		ButtonGroup typeRadio = new ButtonGroup();
//		JRadioButton[] typeRadioButtons = new JRadioButton[MarkerData.TYPES.length];
//		for (int i = 0; i<MarkerData.TYPES.length; i++) {
//			typeRadioButtons[i] = new JRadioButton(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1], false);
//			typeRadioButtons[i].setFont(new Font("Arial", 0, 14));
//			typeRadio.add(typeRadioButtons[i]);
//			typeRadioButtons[i].addItemListener(typeListener);
//			typeRadioButtons[i].setBackground(BACKGROUND_COLOR);
//			typePanel.add(typeRadioButtons[i], gbc);
////			tab1Panel.add(typeRadioButtons[i]);
//		}
//		// --- End of the original block ---
//		
//		/*
//		// --- Beginning of the new block ---
//		typePanelRadioButton = new JPanel();//
//		typePanelRadioButton.setLayout(new GridLayout(MarkerData.TYPES.length,2));//
//		ButtonGroup typeRadio = new ButtonGroup();
//		JRadioButton[] typeRadioButtons = new JRadioButton[MarkerData.TYPES.length];
//		JTextField[] typeRadioTexts = new JTextField[MarkerData.TYPES.length];//
//		for (int i = 0; i<MarkerData.TYPES.length; i++) {
//			//typeRadioButtons[i] = new JRadioButton(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1], false);
//			typeRadioButtons[i] = new JRadioButton(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1], false);//
//			//typeRadioButtons[i].setFont(new Font("Arial", 0, 14));
//			typeRadio.add(typeRadioButtons[i]);
//			typeRadioButtons[i].addItemListener(typeListener);
//			typeRadioButtons[i].setBackground(BACKGROUND_COLOR);
//			//typePanel.add(typeRadioButtons[i], gbc);
//			typePanelRadioButton.add(typeRadioButtons[i]);//
//			typeRadioTexts[i]=new JTextField(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1]);//
//			typeRadioTexts[i].setBackground(BACKGROUND_COLOR);//
//			typeRadioTexts[i].setFont(new Font("Arial",0,14));//
//			typeRadioTexts[i].setBorder(null);//
//			typePanelRadioButton.add(typeRadioTexts[i]);//
//		}
//		typePanel.add(typePanelRadioButton, gbc);//
//		// --- End of the new block ---
//		 */
//
//		JSlider slider = new JSlider(JSlider.HORIZONTAL, 2, 20, DEFAULT_SIZE);
//		slider.setSize(new Dimension(250, 20));
//		slider.setBackground(BACKGROUND_COLOR);
//		sizeLabel = new JLabel("Size = "+DEFAULT_SIZE, JLabel.CENTER);
//		sizeLabel.setFont(new Font("Arial", Font.PLAIN, 16));
//		typePanel.add(sizeLabel, gbc);
////		tab1Panel.add(sizeLabel);
//
//		slider.addChangeListener(new ChangeListener() {
//			public void stateChanged(ChangeEvent ce) {
//				JSlider slider = (JSlider)ce.getSource();
//				sizeLabel.setText("Size = "+slider.getValue());
//				size = (byte)slider.getValue();
//				scatPanel.setPointsGenerated(true);//zx
//				scatPanel.setUpdateQcPanel(false);//zx???
//				scatPanel.paintAgain();
//			}
//		});
//		typePanel.add(slider, gbc);
////		tab1Panel.add(slider);
//
//		slider = new JSlider(JSlider.HORIZONTAL, 0, 100, DEFAULT_SIZE);
//		slider.setBackground(BACKGROUND_COLOR);
//		gcLabel = new JLabel("GC > 0."+DEFAULT_GC_THRESHOLD, JLabel.CENTER);
//		gcLabel.setFont(new Font("Arial", Font.PLAIN, 16));
//		typePanel.add(gcLabel, gbc);
////		tab1Panel.add(gcLabel);
//		slider.addChangeListener(new ChangeListener() {
//			public void stateChanged(ChangeEvent ce) {
//				JSlider slider = (JSlider)ce.getSource();
//				gcThreshold = (float)slider.getValue()/100f;
//				gcLabel.setText("GC > "+ext.formDeci(gcThreshold, 2, true));
//				scatPanel.setPointsGenerated(false);//zx
//				scatPanel.setUpdateQcPanel(true);//zx
//				scatPanel.paintAgain();
//				//qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");//zx
//			}
//		});
//		typePanel.add(slider, gbc);
////		tab1Panel.add(slider);
//
////		typePanel.addTab("Tab 1", null, tab1Panel, "Does nothing");
////		typePanel.setMnemonicAt(0, KeyEvent.VK_1);
////		typePanel.addTab("Tab 2", null, sizeLabel, "Does nothing");
////		typePanel.setMnemonicAt(1, KeyEvent.VK_2);
//
//		ItemListener symmetryListener = new ItemListener() {
//			public void itemStateChanged(ItemEvent ie) {
////				scatPanel.setPointsGenerated(true);//zx ??? Why not true?
////				scatPanel.setUpdateQcPanel(false);//zx ??? Why cannot set to false?
//				updateGUI();
//			}
//		};
//		
//		symmetryBox = new JCheckBox("Symmetric axes");
//		symmetryBox.setFont(new Font("Arial", 0, 14));
//		symmetryBox.addItemListener(symmetryListener);
//		symmetryBox.setBackground(BACKGROUND_COLOR);
//
//		typePanel.add(symmetryBox, gbc);
//		
//		JButton button = new JButton(CAPTURE);
//		button.addActionListener(this);
//		button.setActionCommand(CAPTURE);
//		typePanel.add(button, gbc);
//
//		button = new JButton(DUMP);
//		button.addActionListener(this);
//		button.setActionCommand(DUMP);
//		typePanel.add(button, gbc);
//		
//		button = new JButton(MASK_MISSING);
//		button.addActionListener(this);
//		button.setActionCommand(MASK_MISSING);
//		typePanel.add(button, gbc);
//
//		
//		//----------Begin of Cluster Filter--------------
//		JPanel clusterFilterPanel = new JPanel();
////		JButton first1 = new JButton(Grafik.getImageIcon("images/firstLast/First.gif", true));
////		first1.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif", true));
////		first1.addActionListener(this);
////		first1.setActionCommand(FIRST);
////		first1.setPreferredSize(new Dimension(20, 20));
//		JButton backward = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif", true));
//		backward.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
//		backward.addActionListener(this);
//		backward.setActionCommand(CLUSTER_FILTER_BACKWARD);
//		backward.setPreferredSize(new Dimension(20, 20));
//		clusterFilterNavigation = new JTextField("", 5);
//		clusterFilterNavigation.setHorizontalAlignment(JTextField.CENTER);
//		clusterFilterNavigation.setFont(new Font("Arial", 0, 14));
//		//navigationField.setEditable(false);//zx
//		clusterFilterNavigation.setBackground(BACKGROUND_COLOR);//zx
//		clusterFilterNavigation.addFocusListener(new FocusListener() {
//			public void focusGained(FocusEvent focusevent) {}
//
//			public void focusLost(FocusEvent fe) {
//				try {
//					int trav = Integer.valueOf(((JTextField)fe.getSource()).getText().split("[\\s]+")[0]).intValue()-1;
//					if (trav>=0&&trav<clusterFilterCollection.getSize(getMarkerName())) {
////						currentClusterFilter = (byte) trav;
//						setCurrentClusterFilter((byte) trav);
//					}
//				} catch (NumberFormatException nfe) {}
////				clusterFilterNavigation.setText((currentClusterFilter+1)+" of "+clusterFilterCollection.getSize(getMarkerName()));
//				displayClusterFilterIndex();
//				scatPanel.setPointsGenerated(true);//zx
////				scatPanel.setPointsGenerated(false);//zx
//				scatPanel.setUpdateQcPanel(true);
////				scatPanel.generateRectangles();
//				updateGUI();
//				displayClusterFilterIndex();
//			}
//		});
//
//		JButton forward = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif", true));
//		forward.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
//		forward.addActionListener(this);
//		forward.setActionCommand(CLUSTER_FILTER_FORWARD);
//		forward.setPreferredSize(new Dimension(20, 20));
////		JButton last1 = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif", true));
////		last1.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif", true));
////		last1.addActionListener(this);
////		last1.setActionCommand(LAST);
////		last1.setPreferredSize(new Dimension(20, 20));
//		JButton delete = new JButton(Grafik.getImageIcon("images/delete9.png", true));
//		delete.setDisabledIcon(Grafik.getImageIcon("images/delete8.png", true));
//		delete.addActionListener(this);
//		delete.setActionCommand(CLUSTER_FILTER_DELETE);
//		delete.setPreferredSize(new Dimension(20, 20));
////		clusterFilterPanel.add(first1);
//		clusterFilterPanel.add(backward);
//		clusterFilterPanel.add(clusterFilterNavigation);
//		clusterFilterPanel.add(forward);
////		clusterFilterPanel.add(last1);
//		clusterFilterPanel.add(delete);
//
//    	newGenotype = new JComboBox(new String[] {"-","A/A","A/B","B/B"});
//		ActionListener newGenotypeListener = new ActionListener() {
//			public void actionPerformed(ActionEvent e) {
//				byte newGenotypeSelected;
//				if ((String)newGenotype.getSelectedItem()=="-") {
//					newGenotypeSelected=(byte)-1;
//				} else if ((String)newGenotype.getSelectedItem()=="A/A") {
//					newGenotypeSelected=(byte)0;
//				} else if ((String)newGenotype.getSelectedItem()=="A/B") {
//					newGenotypeSelected=(byte)1;
//				} else if ((String)newGenotype.getSelectedItem()=="B/B") {
//					newGenotypeSelected=(byte)2;
//				} else {
//					newGenotypeSelected=(byte)-9;
//				}
////				byte index = (byte) (clusterFilterCollection.getSize(getMarkerName())-1);
//				clusterFilterCollection.updateGenotype(getMarkerName(), currentClusterFilter, newGenotypeSelected);
//				saveClusterFilterCollection();
////				scatPanel.setPointsGenerated(true);//zx
//				scatPanel.setPointsGenerated(false);//zx
//				scatPanel.setUpdateQcPanel(true);//zx???
//				//updateGUI();
//				scatPanel.paintAgain();
//			}
//			public void comboBoxChanged(ActionEvent e) {
//				System.out.println("Why cannot I see this comboBox changed?");
//			}
//		};
//    	newGenotype.addActionListener(newGenotypeListener);
////		newGenotype.setActionCommand(NEW_GENOTYPE);//???
//    	clusterFilterPanel.add(newGenotype);
////		typePanel.add(newGenotype, gbc);
//
//		clusterFilterPanel.setBackground(BACKGROUND_COLOR);
//		typePanel.add(clusterFilterPanel, gbc);
//		//-----------End of Cluster Filter------------------
//		
//		ItemListener centListener = new ItemListener() {
//			public void itemStateChanged(ItemEvent ie) {
//				int index= ext.indexOfStr(((JCheckBox)ie.getSource()).getText(), centList);
//				displayCents[index] = ((JCheckBox)ie.getSource()).isSelected();
//				centLabels[index].setVisible(displayCents[index]);
//				scatPanel.setPointsGenerated(false);//zx
//				scatPanel.setUpdateQcPanel(false);//zx
//				updateGUI();
//			}
//		};
//		JLabel label = new JLabel("  ");
//		label.setFont(new Font("Arial", 0, 20));
//		typePanel.add(label, gbc);
//		label = new JLabel("Centroids:");
//		label.setFont(new Font("Arial", 0, 20));
//		label.setHorizontalAlignment(JLabel.CENTER);
//		typePanel.add(label, gbc);
//		centBoxes = new JCheckBox[centList.length];
//		displayCents = new boolean[centList.length];
//		centLabels = new JLabel[centList.length];
//		//System.out.println(centList.length+"\n");//zx
//		for (int i = 0; i<centList.length; i++) {
//			centBoxes[i] = new JCheckBox(centList[i]);
//			centBoxes[i].setFont(new Font("Arial", 0, 14));
//			centBoxes[i].setSelected(displayCents[i]);
//			centBoxes[i].addItemListener(centListener);
//			centBoxes[i].setBorder(BorderFactory.createLineBorder(ScatterPanel.DEFAULT_COLORS[5+i], 5));
//			centBoxes[i].setBorderPainted(true);
//			centBoxes[i].setBackground(BACKGROUND_COLOR);
//			typePanel.add(centBoxes[i], gbc);
//			
//			centLabels[i] = new JLabel("LRR correlation not performed");
//			centLabels[i].setVisible(displayCents[i]);
//			typePanel.add(centLabels[i], gbc);
//		}
//		
//        //JLabel padding = new JLabel();//np
//        //gbc.weighty = 1.0;
//        //typePanel.add(padding, gbc);
//		//private JLabel qcPanelLabel;//zx
//		//private JLabel qcCallRateLabel;//zxu
//		//private JLabel qcHwePvalueLabel;//zxu
//
//		qcPanel = new JPanel();//zx
//		qcPanel.setLayout(new GridLayout(8, 1));//zx
//		qcPanel.setBackground(BACKGROUND_COLOR);//zx
//		/*
//        qcPanelLabel = new JLabel("                      ", JLabel.CENTER);//zx
//        qcPanel.add(qcPanelLabel);//zx
//        qcPanelLabel = new JLabel("QC Metrics", JLabel.CENTER);//zx
//        qcPanelLabel.setFont(new Font("Arial", 0, 20));//zx
//        qcPanel.add(qcPanelLabel);//zx
//		
//        qcCallRateLabel = new JLabel("Call Rate: "+ScatterPanel.getCallRate()+"%", JLabel.LEFT);//zx
//        qcCallRateLabel.setFont(new Font("Arial", 0, 14));//zx
//        qcPanel.add(qcCallRateLabel);//zx
//
//        qcHwePvalueLabel = new JLabel("HWE p-value: ", JLabel.LEFT);//zx
//        qcHwePvalueLabel.setFont(new Font("Arial", 0, 14));//zx
//		qcPanel.add(qcHwePvalueLabel);//zx
//		*/
//		typePanel.add(qcPanel);//zx
//
//        typePanel.setBackground(BACKGROUND_COLOR);
//    }

}
