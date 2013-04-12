package cnv.plots;

import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

import cnv.filesys.*;
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
	private ScatterPanel scatPanel;
	private JLabel sizeLabel;
	private JLabel gcLabel;
	private JPanel typePanel;

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
	private JLabel markerName, commentLabel;
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
	private Color[] colorScheme;

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
		colorScheme = scatPanel.getColorScheme();
		getContentPane().add(scatPanel, BorderLayout.CENTER);

		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new GridLayout(3, 1));
		markerName = new JLabel("", JLabel.CENTER);
		markerName.setFont(new Font("Arial", 0, 20));
		descrPanel.add(markerName);

		commentLabel = new JLabel("", JLabel.CENTER);
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
		// typePanel.setLayout(new BoxLayout(typePanel, BoxLayout.PAGE_AXIS));
//		typePanel.setLayout(new GridLayout(20, 1));
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
		for (int i = 0; i<centList.length; i++) {
			centBoxes[i] = new JCheckBox(centList[i]);
			centBoxes[i].setFont(new Font("Arial", 0, 14));
			centBoxes[i].setSelected(displayCents[i]);
			centBoxes[i].addItemListener(centListener);
			centBoxes[i].setBorder(BorderFactory.createLineBorder(colorScheme[5+i], 5));
			centBoxes[i].setBorderPainted(true);
			centBoxes[i].setBackground(BACKGROUND_COLOR);
			typePanel.add(centBoxes[i], gbc);
			
			centLabels[i] = new JLabel("LRR correlation not performed");
			centLabels[i].setVisible(displayCents[i]);
			typePanel.add(centLabels[i], gbc);
		}
		
        JLabel padding = new JLabel();   
        gbc.weighty = 1.0;   
        typePanel.add(padding, gbc);   

        typePanel.setBackground(BACKGROUND_COLOR);
		getContentPane().add(typePanel, BorderLayout.EAST);

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
		getContentPane().add(classPanel, BorderLayout.SOUTH);
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
		
		filename = proj.getFilename(Project.DISPLAY_MARKERS_FILENAME);
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
						System.err.println("Error - could not find "+line[0]+" in the lookup table");
					}
				}
				reader.close();
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
		System.out.println("Found "+files.length+" .cent files in "+proj.getDir(Project.DATA_DIRECTORY));
		
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
