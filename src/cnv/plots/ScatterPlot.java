package cnv.plots;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;

import stats.ContingencyTable;
import stats.ProbDist;
import stats.CTable;

import java.awt.*;
import java.awt.event.*;

import cnv.filesys.*;
import cnv.gui.ColorIcon;
import cnv.gui.CycleRadio;
import common.*;
import cnv.var.*;

public class ScatterPlot extends JFrame implements ActionListener {
	public static final long serialVersionUID = 1L;
	public static final byte DEFAULT_SIZE = 8;
	public static final int DEFAULT_GC_THRESHOLD = 25;
	private static final String ALT_UP = "ALT UP";
	private static final String ALT_DOWN = "ALT DOWN";
	private static final String ALT_LEFT = "ALT LEFT";
	private static final String ALT_RIGHT = "ALT RIGHT";
	private static final String FIRST = "First";
	private static final String PREVIOUS = "Previous";
	private static final String NEXT = "Next";
	private static final String LAST = "Last";
	private static final String CAPTURE = "Screen capture";
	private static final String DUMP = "Dump raw data";
	public static final String MASK_MISSING = "Mask missing values";
	public static final String UNMASK_MISSING = "Unmask missing values";
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	private JButton first, previous, next, last;
	private JTextField navigationField;
	private JPanel classPanel;
	private JPanel legendPanel;
	private JPanel bottomPanel;
	private ScatterPanel scatPanel;
	private JLabel sizeLabel;
	private JLabel gcLabel;
	private JPanel typePanel;
	private JPanel qcPanel;//zx
	//private JLabel qcPanelLabel;//zx
	//private JLabel qcCallRateLabel;//zxu
	//private JLabel qcHwePvalueLabel;//zxu
	//private JPanel typePanelRadioButton;//

	private Project proj;
	private MarkerData[] markerData;
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
	private int currentClass;
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

	public ScatterPlot(Project project) {
		super("ScatterPlot");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		proj = project;
		jar = proj.getJarStatus();
		size = DEFAULT_SIZE;
		gcThreshold = (float)DEFAULT_GC_THRESHOLD/100f;

		SampleList list = proj.getSampleList();
		samples = list.getSamples();
		sampleListFingerprint = list.getFingerprint();
		sampleData = new SampleData(proj, true);
		markerLookup = proj.getMarkerLookup();
		loadFile();
		if (markerList == null) {
			return;
		}
		loadCentroids();

		scatPanel = new ScatterPanel(this);
		getContentPane().add(scatPanel, BorderLayout.CENTER);

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
						updateGUI();
					}
				} catch (NumberFormatException nfe) {}
				displayIndex((JTextField)fe.getSource());
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
		getContentPane().add(descrPanel, BorderLayout.NORTH);

		typePanel = new JPanel();
		//typePanel.setLayout(new BoxLayout(typePanel, BoxLayout.PAGE_AXIS));
		//typePanel.setLayout(new GridLayout(20, 1));
		typePanel.setLayout(new GridBagLayout());
		
        GridBagConstraints gbc = new GridBagConstraints();   
        gbc.insets = new Insets(1,3,0,30);   
        gbc.weightx = 1.0;   
        gbc.fill = GridBagConstraints.HORIZONTAL;   
        gbc.gridwidth = GridBagConstraints.REMAINDER;   

		ItemListener typeListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				JRadioButton jrb = (JRadioButton)ie.getItem();
				if (jrb.isSelected()) {
					for (int i = 0; i<MarkerData.TYPES.length; i++) {
						if (jrb.getText().equals(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1])) {
							plot_type = i;
							updateGUI();
						}
					}
				}
			}
		};
		// --- Beginning of the original block ---
		ButtonGroup typeRadio = new ButtonGroup();
		JRadioButton[] typeRadioButtons = new JRadioButton[MarkerData.TYPES.length];
		for (int i = 0; i<MarkerData.TYPES.length; i++) {
			typeRadioButtons[i] = new JRadioButton(MarkerData.TYPES[i][0]+"/"+MarkerData.TYPES[i][1], false);
			typeRadioButtons[i].setFont(new Font("Arial", 0, 14));
			typeRadio.add(typeRadioButtons[i]);
			typeRadioButtons[i].addItemListener(typeListener);
			typeRadioButtons[i].setBackground(BACKGROUND_COLOR);
			typePanel.add(typeRadioButtons[i], gbc);
		}
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

		JSlider slider = new JSlider(JSlider.HORIZONTAL, 2, 20, DEFAULT_SIZE);
		slider.setSize(new Dimension(250, 20));
		slider.setBackground(BACKGROUND_COLOR);
		sizeLabel = new JLabel("Size = "+DEFAULT_SIZE, JLabel.CENTER);
		sizeLabel.setFont(new Font("Arial", Font.PLAIN, 16));
		typePanel.add(sizeLabel, gbc);
		slider.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent ce) {
				JSlider slider = (JSlider)ce.getSource();
				sizeLabel.setText("Size = "+slider.getValue());
				size = (byte)slider.getValue();
				scatPanel.paintAgain();
			}
		});
		typePanel.add(slider, gbc);

		slider = new JSlider(JSlider.HORIZONTAL, 0, 100, DEFAULT_SIZE);
		slider.setBackground(BACKGROUND_COLOR);
		gcLabel = new JLabel("GC > 0."+DEFAULT_GC_THRESHOLD, JLabel.CENTER);
		gcLabel.setFont(new Font("Arial", Font.PLAIN, 16));
		typePanel.add(gcLabel, gbc);
		slider.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent ce) {
				JSlider slider = (JSlider)ce.getSource();
				gcThreshold = (float)slider.getValue()/100f;
				gcLabel.setText("GC > "+ext.formDeci(gcThreshold, 2, true));
				scatPanel.paintAgain();
				//qcCallRateLabel.setText("Call Rate: "+ScatterPanel.getCallRate()+"%");//zx
			}
		});
		typePanel.add(slider, gbc);


		ItemListener symmetryListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				updateGUI();
			}
		};
		
		symmetryBox = new JCheckBox("Symmetric axes");
		symmetryBox.setFont(new Font("Arial", 0, 14));
		symmetryBox.addItemListener(symmetryListener);
		symmetryBox.setBackground(BACKGROUND_COLOR);

		typePanel.add(symmetryBox, gbc);
		
		JButton button = new JButton(CAPTURE);
		button.addActionListener(this);
		button.setActionCommand(CAPTURE);
		typePanel.add(button, gbc);

		button = new JButton(DUMP);
		button.addActionListener(this);
		button.setActionCommand(DUMP);
		typePanel.add(button, gbc);
		
		button = new JButton(MASK_MISSING);
		button.addActionListener(this);
		button.setActionCommand(MASK_MISSING);
		typePanel.add(button, gbc);
		
		ItemListener centListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				int index= ext.indexOfStr(((JCheckBox)ie.getSource()).getText(), centList);
				displayCents[index] = ((JCheckBox)ie.getSource()).isSelected();
				centLabels[index].setVisible(displayCents[index]);
				updateGUI();
			}
		};
		JLabel label = new JLabel("  ");
		label.setFont(new Font("Arial", 0, 20));
		typePanel.add(label, gbc);
		label = new JLabel("Centroids:");
		label.setFont(new Font("Arial", 0, 20));
		label.setHorizontalAlignment(JLabel.CENTER);
		typePanel.add(label, gbc);
		centBoxes = new JCheckBox[centList.length];
		displayCents = new boolean[centList.length];
		centLabels = new JLabel[centList.length];
		//System.out.println(centList.length+"\n");//zx
		for (int i = 0; i<centList.length; i++) {
			centBoxes[i] = new JCheckBox(centList[i]);
			centBoxes[i].setFont(new Font("Arial", 0, 14));
			centBoxes[i].setSelected(displayCents[i]);
			centBoxes[i].addItemListener(centListener);
			centBoxes[i].setBorder(BorderFactory.createLineBorder(ScatterPanel.DEFAULT_COLORS[5+i], 5));
			centBoxes[i].setBorderPainted(true);
			centBoxes[i].setBackground(BACKGROUND_COLOR);
			typePanel.add(centBoxes[i], gbc);
			
			centLabels[i] = new JLabel("LRR correlation not performed");
			centLabels[i].setVisible(displayCents[i]);
			typePanel.add(centLabels[i], gbc);
		}
		
        //JLabel padding = new JLabel();//np
        //gbc.weighty = 1.0;
        //typePanel.add(padding, gbc);
		//private JLabel qcPanelLabel;//zx
		//private JLabel qcCallRateLabel;//zxu
		//private JLabel qcHwePvalueLabel;//zxu

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
		typePanel.add(qcPanel);//zx

        typePanel.setBackground(BACKGROUND_COLOR);
		getContentPane().add(typePanel, BorderLayout.EAST);

		bottomPanel = new JPanel();
		bottomPanel.setLayout(new GridLayout(2, 1));
        //bottomPanel.setBackground(BACKGROUND_COLOR);
		classPanel = new JPanel();
		label = new JLabel("Color code by:");
		label.setFont(new Font("Arial", 0, 14));
		classPanel.add(label);

		ItemListener classListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				JRadioButton jrb = (JRadioButton)ie.getItem();
				if (jrb.isSelected()) {
					for (int i = 0; i<sampleData.getNumClasses(); i++) {
						if (jrb.getText().equals(sampleData.getClassName(i))) {
							currentClass = i;
							updateGUI();
						}
					}
				}
			}
		};
		ButtonGroup classRadio = new ButtonGroup();
		JRadioButton[] classRadioButtons = new JRadioButton[sampleData.getNumClasses()];
		for (int i = 0; i<sampleData.getNumClasses(); i++) {
			classRadioButtons[i] = new JRadioButton(sampleData.getClassName(i), false);
			classRadioButtons[i].setFont(new Font("Arial", 0, 14));
			classRadio.add(classRadioButtons[i]);
			classRadioButtons[i].addItemListener(classListener);
			classRadioButtons[i].setBackground(BACKGROUND_COLOR);
			classPanel.add(classRadioButtons[i]);
		}
		classPanel.setBackground(BACKGROUND_COLOR);
		bottomPanel.add(classPanel);
		
		legendPanel = new JPanel();
        legendPanel.setBackground(BACKGROUND_COLOR);
		//JLabel legend1 = new JLabel("Color Key: ");
		//legend1.setFont(new Font("Arial", 0, 14));
		//legendPanel.add(legend1);
		
		/*
		for (int i=0; i<sampleData.getActualClassColorKey(currentClass).length; i++){
			legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[i][0])])));
			legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[i][1]));
		}
		legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[0][0])])));
		legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[0][1]));
		legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[1][0])])));
		legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[1][1]));
		
		// test point 1
		System.out.println("Length of the two dimisional Array: \t"+sampleData.getActualClassColorKey(currentClass).length);
		System.out.println("currentClass\t dimension2\t dimension3\t Value");
		System.out.println(currentClass+"\t 0\t 0\t"+sampleData.getActualClassColorKey(currentClass)[0][0]);
		System.out.println(currentClass+"\t 0\t 1\t"+sampleData.getActualClassColorKey(currentClass)[0][1]);
		System.out.println(currentClass+"\t 1\t 1\t"+sampleData.getActualClassColorKey(currentClass)[1][0]);
		System.out.println(currentClass+"\t 1\t 1\t"+sampleData.getActualClassColorKey(currentClass)[1][1]);
		*/
			
		bottomPanel.add(legendPanel);
		
		getContentPane().add(bottomPanel, BorderLayout.SOUTH);
		classRadioButtons[1].setSelected(true);

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
				updateGUI();
			}
		});
		actionMap.put(PREVIOUS, new AbstractAction() {
			public static final long serialVersionUID = 5L;

			public void actionPerformed(ActionEvent e) {
				markerIndex = Math.max(markerIndex-1, 0);
				displayIndex(navigationField);
				updateGUI();
			}
		});
		actionMap.put(NEXT, new AbstractAction() {
			public static final long serialVersionUID = 6L;

			public void actionPerformed(ActionEvent e) {
				markerIndex = Math.min(markerIndex+1, markerList.length-1);
				displayIndex(navigationField);
				updateGUI();
			}
		});
		actionMap.put(LAST, new AbstractAction() {
			public static final long serialVersionUID = 7L;

			public void actionPerformed(ActionEvent e) {
				markerIndex = markerList.length-1;
				displayIndex(navigationField);
				updateGUI();
			}
		});
		scatPanel.setActionMap(actionMap);

		updateGUI();
		displayIndex(navigationField);
		typeRadioButtons[1].setSelected(true);
		symmetryBox.setSelected(true);
		if (centList.length > 0) {
			centBoxes[0].setSelected(true);
		}
		
		next.getInputMap().put(KeyStroke.getKeyStroke("space"), NEXT);
		next.setActionMap(actionMap);
		previous.setActionMap(actionMap);
		scatPanel.grabFocus();

		setBounds(20, 20, 1000, 720);
		setVisible(true);
	}

	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();
		String filename;
		int count;

		if (command.equals(FIRST)) {
			markerIndex = 0;
			displayIndex(navigationField);
			updateGUI();
		} else if (command.equals(PREVIOUS)) {
			markerIndex = Math.max(markerIndex-1, 0);
			displayIndex(navigationField);
			updateGUI();
		} else if (command.equals(NEXT)) {
			markerIndex = Math.min(markerIndex+1, markerList.length-1);
			displayIndex(navigationField);
			updateGUI();
		} else if (command.equals(LAST)) {
			markerIndex = markerList.length-1;
			displayIndex(navigationField);
			updateGUI();
		} else if (command.equals(CAPTURE)) {
			System.out.println("command.equals(CAPTURE)");//zx
			count = 1;
			do {
				filename = markerList[markerIndex]+"_"+MarkerData.TYPES[plot_type][0]+"-"+MarkerData.TYPES[plot_type][1]+(count==1?"":" v"+count);
				System.out.println(filename);//zx
				count++;
			} while (new File(proj.getProjectDir()+filename+".png").exists());
			scatPanel.screenCapture(proj.getProjectDir()+filename+".png");
		} else if (command.equals(DUMP)) {
			count = 1;
			do {
				filename = markerList[markerIndex]+"_dump"+(count==1?"":" v"+count);
				count++;
			} while (new File(proj.getProjectDir()+filename+".xln").exists());
			markerData[markerIndex].writeToFile(samples, proj.getProjectDir()+filename+".xln");
		} else if (command.equals(MASK_MISSING) || command.equals(UNMASK_MISSING)) {
			maskMissing = !maskMissing;
			((JButton)ae.getSource()).setText(maskMissing?UNMASK_MISSING:MASK_MISSING);
			updateGUI();
		} else {
			System.err.println("Error - unknown command '"+command+"'");
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

	public MarkerData[] getMarkerData() {
		return markerData;
	}

	public int getMarkerIndex() {
		return markerIndex;
	}

	public int getCurrentClass() {
		return currentClass;
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

	public float[][][][] getCents() {
		return cents;
	}
	
	public void updateCentLabels() {
		double[] comp;
		String str;

		if (markerData[markerIndex].getLRRs() != null) {
			for (int i = 0; i<centList.length; i++) {
				comp = markerData[markerIndex].compareLRRs(cents[i][markerIndex]);
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

	public boolean maskMissing() {
		return maskMissing;
	}
	
	public void loadFile() {
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
				System.err.println("Error: file \""+filename+"\" not found in current directory");
				System.exit(1);
			} catch (Exception e) {
				System.err.println("Error reading file \""+filename+"\"");
				e.printStackTrace();
				System.exit(2);
			}

			markerList = Array.toStringArray(markerNames);
			commentList = Array.toStringArray(markerComments);
			markerData = MarkerSet.loadFromList(proj, markerList);
			//System.out.println("markerName: "+markerData.+);//zx
			markerIndex = 0;
			previousMarkerIndex = -1;
		} catch (Exception e) {
			System.err.println("Error loading: "+filename);
			e.printStackTrace();
		}
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
		if (markerList.length==0) {
			markerName.setText("Error: marker data was not successfully loaded");
			commentLabel.setText("Check to make sure MarkerLookup is synchronized with the current data");
			classPanel.setEnabled(false);
		} else {
			markerName.setText(markerList[markerIndex]);
			commentLabel.setText(commentList[markerIndex]);

			if (classPanel!=null) {
				classPanel.setEnabled(markerData[markerIndex].getFingerprint()==sampleListFingerprint);
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
						markerData[markerIndex].recompute(cents[i][markerIndex]);
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
		JLabel label, block;
		String[][] colorKeys;
		String[] keys;
		legendPanel.removeAll();
		legendPanel.repaint();
		//JLabel legend = new JLabel("Color Key: ");
		//legend.setFont(new Font("Arial", 0, 14));
		//legendPanel.add(legend);
		
		//legendPanel.removeAll();
		//legendPanel.repaint();
		label = new JLabel("Color key:");
		label.setFont(new Font("Arial", 0, 14));
		legendPanel.add(label);
		if (currentClass < SampleData.BASIC_CLASSES.length) {
			colorKeys = SampleData.KEYS_FOR_BASIC_CLASSES[currentClass];
		} else {		
			colorKeys = sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length);
		}
		for (int i = 0; i<colorKeys.length; i++) {
			block = new JLabel(new ColorIcon(12, 12, StratPanel.DEFAULT_COLORS[Integer.parseInt(colorKeys[i][0])]));
			label = new JLabel(colorKeys[i][1]+" (n="+(hash.containsKey(colorKeys[i][0])?hash.get(colorKeys[i][0]):"0")+")");
			hash.remove(colorKeys[i][0]);
			label.setFont(new Font("Arial", 0, 14));
			legendPanel.add(block);
			legendPanel.add(label);
		}
		keys = HashVec.getKeys(hash);
		for (int i = 0; i<keys.length; i++) {
			if (!keys[i].equals("-1")) {
				block = new JLabel(new ColorIcon(12, 12, StratPanel.DEFAULT_COLORS[Integer.parseInt(keys[i])]));
				label = new JLabel((keys[i].equals("0")?"missing":keys[i])+" (n="+hash.get(keys[i])+")");
				label.setFont(new Font("Arial", 0, 14));
				legendPanel.add(block);
				legendPanel.add(label);
			}
		}

		
/**		
		//colorKeys = sampleData.getActualClassColorKey(currentClass);
		for (int i=0; i<sampleData.getActualClassColorKey(currentClass).length; i++){
			legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[i][0])])));
			legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[i][1]));
			hash.remove(arg0)
		}
//		for (int i = 0; i<colorKeys.length; i++) {
//			label = new JLabel(colorKeys[i][1]+" (n="+(hash.containsKey(colorKeys[i][0])?hash.get(colorKeys[i][0]):"0")+")");
//		}
		/*
		keys = HashVec.getKeys(hash);
		for (int i = 0; i<keys.length; i++) {
			if (!keys[i].equals("-1")) {
				block = new JLabel(new ColorIcon(12, 12, StratPanel.DEFAULT_COLORS[Integer.parseInt(keys[i])]));
				label = new JLabel((keys[i].equals("0")?"missing":keys[i])+" (n="+hash.get(keys[i])+")");
				label.setFont(new Font("Arial", 0, 14));
				legendPanel.add(block);
				legendPanel.add(label);
			}
		}
		*/
//*/

		legendPanel.validate();
		
		/*
		legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[0][0])])));
		legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[0][1]));
		legendPanel.add(new JLabel(new ColorIcon(12,12,scatPanel.DEFAULT_COLORS[Integer.parseInt(sampleData.getActualClassColorKey(currentClass)[1][0])])));
		legendPanel.add(new JLabel(sampleData.getActualClassColorKey(currentClass)[1][1]));
		*/
	}


	public void updateQcPanel(int[][] dataForQc) {
		float callRate=0;//zx
		JLabel qcPanelLabel;//zx
		//JLabel qcCallRateLabel;//zxu
		//JLabel qcHwePvalueLabel;//zxu
		int[] alleleCounts;
		double hweP;
		//int[][] sexContingecyTable = new int[2][2];
		CTable classCount;
		
		alleleCounts = new int[3]; 
		for (int i=0; i<dataForQc[0].length; i++){//zx
			if (dataForQc[0][i]<=0) {
				callRate++;//zx
			} else {
				alleleCounts[dataForQc[0][i]-1]++;
			}
		}
		hweP = AlleleFreq.HWEsig(alleleCounts);
		callRate=(dataForQc[0].length-callRate)*100/dataForQc[0].length;//zx
		
		qcPanel.removeAll();
		qcPanel.repaint();
		
        qcPanelLabel = new JLabel("", JLabel.CENTER);//zx
        qcPanel.add(qcPanelLabel);//zx
        qcPanelLabel = new JLabel("QC Metrics", JLabel.CENTER);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 20));//zx
        qcPanel.add(qcPanelLabel);//zx
		
        //qcCallRateLabel = new JLabel("Call Rate: "+(dataForQc.length-missing)*100/dataForQc.length+"%", JLabel.LEFT);//zx
        qcPanelLabel = new JLabel("Callrate: "+callRate+"%"+"                           ", JLabel.LEFT);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
        qcPanel.add(qcPanelLabel);//zx

        qcPanelLabel = new JLabel("HWE p-value: "+ext.prettyP(hweP), JLabel.LEFT);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcPanelLabel);//zx

		ToolTipManager.sharedInstance().setDismissDelay(100000);

		classCount = new CTable(dataForQc[0], dataForQc[1], SampleData.KEYS_FOR_BASIC_CLASSES[1],
				sampleData.getActualClassColorKey(0));
		qcPanelLabel = new JLabel("Callrate by "+sampleData.getClassName(2)+": "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTableForCallRate()), 1) ), JLabel.LEFT);//zx
		qcPanelLabel.setToolTipText(classCount.generateToolTipTextForCallRate());
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcPanelLabel);//zx

		if (currentClass>2) {
			classCount = new CTable(dataForQc[0], dataForQc[2], SampleData.KEYS_FOR_BASIC_CLASSES[1],
					currentClass<SampleData.BASIC_CLASSES.length?SampleData.KEYS_FOR_BASIC_CLASSES[currentClass]:sampleData.getActualClassColorKey(currentClass-SampleData.BASIC_CLASSES.length));
			qcPanelLabel = new JLabel("Callrate by "+sampleData.getClassName(currentClass)+": "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getContingencyTableForCallRate()), 1) ), JLabel.LEFT);//zx
			qcPanelLabel.setToolTipText(classCount.generateToolTipTextForCallRate());
	        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
			qcPanel.add(qcPanelLabel);//zx
		}

		classCount = new CTable(dataForQc[0], dataForQc[1], SampleData.KEYS_FOR_BASIC_CLASSES[1],
				sampleData.getActualClassColorKey(0));
		qcPanelLabel = new JLabel("Allele Freq by Sex: "+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(classCount.getAlleleFreqBySex()), 1)), JLabel.LEFT);//zx
		qcPanelLabel.setToolTipText(classCount.generateToolTipTextForAllelFreq());
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcPanelLabel);//zx

		qcPanelLabel = new JLabel("Minor Allele Freq: " + (new DecimalFormat("#.####").format(classCount.getMinorAlleleFrequency())), JLabel.LEFT);//zx
        qcPanelLabel.setFont(new Font("Arial", 0, 14));//zx
		qcPanel.add(qcPanelLabel);//zx

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

	public static void main(String[] args) {
		String filename = Project.DEFAULT_SCATTER_PROJECT;
		boolean jar = args.length>0&&args[0].equals("-notJar")?false:true;

		try {
			new ScatterPlot(new Project(filename, jar));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
