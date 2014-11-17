package cnv.plots;

import java.io.*;
import java.util.*;
import common.*;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.event.*;

import cnv.filesys.Project;
import cnv.var.*;
import cnv.gui.*;

public class StratPlot extends JFrame implements ActionListener, TreeSelectionListener {
	public static final long serialVersionUID = 1L;
	public static final String DEFAULT_FILENAME = "samplesCliskedInStratify.xln";
	public static final String[] MOSAICISM_HEADER = {"Sample", "Band", "LRR N", "mean LRR", "BAF N", "SD of BAF (0.15-0.85)", "IQR of BAF (0.15-0.85)", "%Homo"};
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String SWAP_AXES = "Swap Axes";
	public static final String INVERT_X = "Invert X axis";
	public static final String INVERT_Y = "Invert Y axis";
	public static final String MASK_MISSING = "Mask missing values";
	public static final String UNMASK_MISSING = "Unmask missing values";
	private static final String REFRESH_SAMPLE_DATA = "Refresh SampleData";
	public static final String[] BUTTONS = {SWAP_AXES, INVERT_X, INVERT_Y, MASK_MISSING, REFRESH_SAMPLE_DATA};

	private Project proj;
	private CheckBoxTree tree;
	private StratPanel stratPanel;
	private JLabel descriptor;
	private JPanel legendPanel;
	private SampleData sampleData;
	private int currentVariable;
	private boolean swapAxes;
	private boolean maskMissing;
	private Logger log;
//	private boolean fail; // TODO Use this
	
	// TODO need to move frame, etc out of constructor to fail as the others do
	public StratPlot(Project project, String[][] names, Hashtable<String,float[][]> hash) {
		super("Genvisis - Stratify - " + project.getNameOfProject());
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		JPanel treePanel, classPanel, classPanelTop;
		
		proj = project;
		log = proj.getLog();

		treePanel = new JPanel();
		treePanel.setBackground(BACKGROUND_COLOR);
		treePanel.setLayout(new BorderLayout());
		
		JPanel buttonPanel = new JPanel();
		buttonPanel.setBackground(BACKGROUND_COLOR);
		buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.Y_AXIS));
		
		JButton button;
		for (int i = 0; i<BUTTONS.length; i++) {
			button = new JButton(BUTTONS[i]);
			button.setMinimumSize(new Dimension(200, 20));
			button.setMaximumSize(new Dimension(200, 20));
	        button.setAlignmentX(Component.CENTER_ALIGNMENT);
			button.addActionListener(this);
			buttonPanel.add(button);
        }
		treePanel.add(buttonPanel, BorderLayout.NORTH);
		
		tree = new CheckBoxTree(names, 2);
		treePanel.add(new JScrollPane(tree), BorderLayout.CENTER);
		getContentPane().add(treePanel, BorderLayout.WEST);
		tree.addTreeSelectionListener(this);
		treePanel.setPreferredSize(new Dimension(200,500));

		
		sampleData = proj.getSampleData(3, false);
		if (!sampleData.containsDNA()) {
			log.reportError("Without a DNA column in the SampleData file, ScatterPlot will not start");
			JOptionPane.showMessageDialog(null, "Error - Without a DNA column in the SampleData file, ScatterPlot will not start", "Error",
					JOptionPane.ERROR_MESSAGE);
			return;
		}
		
		if (sampleData.getNumActualClasses() == -1) {
			log.reportError("Error - Failed to load SampleData... closing");
			return;
		}
		if (sampleData.getNumActualClasses() == 0) {
			log.reportError("Error - this isn't going to work; you don't have any classes defined in SampleData.txt");
			return;
		}
		
		stratPanel = new StratPanel(this, names, hash);
		// panel.setToolTipText("");
		getContentPane().add(stratPanel, BorderLayout.CENTER);

		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new GridLayout(2, 1));

		JLabel label = new JLabel("Stratification Plot", JLabel.CENTER);
		label.setFont(new Font("Arial", 0, 20));
		descrPanel.add(label);

		descriptor = new JLabel("", JLabel.CENTER);
		descriptor.setFont(new Font("Arial", 0, 14));
		descrPanel.add(descriptor);
		descrPanel.setBackground(BACKGROUND_COLOR);

		getContentPane().add(descrPanel, BorderLayout.NORTH);

		classPanel = new JPanel(new BorderLayout());
		classPanelTop = new JPanel();
		legendPanel = new JPanel();
		label = new JLabel("Color code by:");
		label.setFont(new Font("Arial", 0, 14));
		classPanelTop.add(label);

		ItemListener classListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				JRadioButton jrb = (JRadioButton)ie.getItem();
				if (jrb.isSelected()) {
					for (int i = 0; i<sampleData.getNumActualClasses(); i++) {
						if (jrb.getText().equals(sampleData.getActualClassName(i))) {
							currentVariable = i;
							updateGUI();
						}
					}
				}
			}
		};
		ButtonGroup classRadio = new ButtonGroup();
		JRadioButton[] classRadioButtons = new JRadioButton[sampleData.getNumActualClasses()];
		for (int i = 0; i<sampleData.getNumActualClasses(); i++) {
			classRadioButtons[i] = new JRadioButton(sampleData.getActualClassName(i), false);
			classRadioButtons[i].setFont(new Font("Arial", 0, 14));
			classRadio.add(classRadioButtons[i]);
			classRadioButtons[i].addItemListener(classListener);
			classRadioButtons[i].setBackground(BACKGROUND_COLOR);
			classPanelTop.add(classRadioButtons[i]);
		}
		classPanelTop.setBackground(BACKGROUND_COLOR);
		classPanelTop.setLayout(new FlowLayout());
		classPanelTop.addComponentListener(new JPanelFlowLayoutComponentListener());
		classRadioButtons[0].setSelected(true);

		legendPanel.setBackground(BACKGROUND_COLOR);
		legendPanel.setLayout(new FlowLayout());
		legendPanel.add(new JLabel("Place holder"));
//		classPanelBottom.addComponentListener(new JPanelFlowLayoutComponentListener());
		
		classPanel.setBackground(BACKGROUND_COLOR);
		classPanel.add(classPanelTop, BorderLayout.NORTH);
		classPanel.add(legendPanel, BorderLayout.SOUTH);
		getContentPane().add(classPanel, BorderLayout.SOUTH);

//		tree.selectFirstTwo();
		
		repaint();

		setBounds(20, 20, 1000, 720);
//		setBounds(-10, 170, 990, 788);
		setVisible(true);

//		unnecessary leads to a double rendering
		stratPanel.createImage();
		stratPanel.repaint();
	}
	
	public Project getProject() {
		return proj;
	}
	
	public void updateGUI() {
		stratPanel.paintAgain();
	}
	
	public void updateColorKey(Hashtable<String,String> hash) {
		JLabel label, block;
		String[][] colorKeys;
		String[] keys;
		
		legendPanel.removeAll();
		legendPanel.repaint();
		label = new JLabel("Color key:");
		label.setFont(new Font("Arial", 0, 14));
		legendPanel.add(label);
		colorKeys = sampleData.getActualClassColorKey(currentVariable);
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
		
		legendPanel.validate();
	}

	public SampleData getSampleData() {
		return sampleData;
	}
	
	public void setDescriptor(String str) {
		descriptor.setText(str);
	}

	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();
		
		if (command.equals(SWAP_AXES)) {
			swapAxes = !swapAxes;
			updateGUI();
		} else if (command.equals(INVERT_X)) {
			stratPanel.toggleXinversion();
			updateGUI();
		} else if (command.equals(INVERT_Y)) {
			stratPanel.toggleYinversion();
			updateGUI();
		} else if (command.equals(MASK_MISSING) || command.equals(UNMASK_MISSING)) {
			maskMissing = !maskMissing;
			((JButton)ae.getSource()).setText(maskMissing?UNMASK_MISSING:MASK_MISSING);
			updateGUI();
		} else if (command.equals(REFRESH_SAMPLE_DATA)) {
			proj.resetSampleData();
			sampleData = proj.getSampleData(2, false);
			stratPanel.pushSampleData();
			updateGUI();
		} else {
			log.reportError("Error - unknown command '"+command+"'");
		}
	}
	
	public void valueChanged(TreeSelectionEvent tse) {
		updateGUI();
	}
	
	public int getCurrentVariable() {
		return currentVariable;
	}

	public int[][] getCurrentPair() {
		int[][] currentPair = tree.getSelectionIndices();
//		System.out.println(Array.toStr(currentPair[1]) +"\t"+ Array.toStr(currentPair[0])); // TODO this being called twice each time 
		return swapAxes?new int[][] {currentPair[1], currentPair[0]}:currentPair;
	}
	
	public boolean maskMissing() {
		return maskMissing;
	}

	public static void loadStratificationResults(Project proj) {
		Hashtable<String, float[][]> hash;
		Vector<String> stratFiles;
		BufferedReader reader;
		String[][] names;
		float[][] data;
		String[] line;
		boolean sol;
		String trav;
		int n;
		Logger log;

		log = proj.getLog();
		stratFiles = proj.getStratResults();
		names = new String[stratFiles.size()][];
		hash = new Hashtable<String,float[][]>();
		for (int i = 0; i<names.length; i++) {
			try {
	            reader = Files.getReader(stratFiles.elementAt(i), proj.getJarStatus(), true, false);
	            line = reader.readLine().trim().split("[\\s]+");
	            if (!line[0].equals("FID") || !line[1].equals("IID")) {
	            	log.reportError("Error - different format than expected; first two columns should be FID and IID");
	            	throw new IOException();
	            }
	            sol = line[2].equals("SOL");
	            n = line.length-(sol?3:2);
	            names[i] = new String[n+1];
	            for (int j = 0; j<n; j++) {
	            	names[i][j+1] = line[(sol?3:2)+j];
                }
	            while (reader.ready()) {
	            	line = reader.readLine().trim().split("[\\s]+");
	            	trav = line[0]+"\t"+line[1];
	            	if (hash.containsKey(trav)) {
	            		data = hash.get(trav);
	            	} else {
	            		hash.put(trav, data = new float[stratFiles.size()][]);
	            	}
	            	data[i] = new float[n];
	            	for (int j = 0; j<n; j++) {
	            		data[i][j] = line[(sol?3:2)+j].equals(".")?Float.NaN:Float.parseFloat(line[(sol?3:2)+j]);
                    }
	            }
	            reader.close();
            } catch (FileNotFoundException fnfe) {
            	log.reportError("Error: file \""+stratFiles.elementAt(i)+"\" not found in current directory");
	            names[i] = new String[1];
            } catch (IOException ioe) {
            	log.reportError("Error reading file \""+stratFiles.elementAt(i)+"\"");
	            names[i] = new String[1];
            }
            names[i][0] = ext.rootOf(stratFiles.elementAt(i), true);
        }

		new StratPlot(proj, names, hash);
	}

	public static void main(String[] args) {
		try {
			loadStratificationResults(new Project(cnv.Launch.getDefaultDebugProjectFile(true), false));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
