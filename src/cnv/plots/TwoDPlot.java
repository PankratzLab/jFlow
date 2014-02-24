package cnv.plots;

import java.io.*;
import java.util.*;

import javax.swing.*;
import javax.swing.event.*;

import java.awt.*;
import java.awt.event.*;

import cnv.filesys.*;
import cnv.gui.CheckBoxTree;
import cnv.gui.ColorKeyPanel;
import cnv.gui.GuiManager;
import common.*;
import cnv.var.*;

public class TwoDPlot extends JPanel implements WindowListener, ActionListener, TreeSelectionListener { 
	public static final long serialVersionUID = 1L;
	public static final byte DEFAULT_SIZE = 8;
	public static final int DEFAULT_GC_THRESHOLD = 25;
	private static final String ALT_UP = "ALT UP";
	private static final String ALT_DOWN = "ALT DOWN";
	private static final String ALT_LEFT = "ALT LEFT";
	private static final String ALT_RIGHT = "ALT RIGHT";
//	public static final String MASK_MISSING = "Mask missing values";
//	public static final String UNMASK_MISSING = "Unmask missing values";
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String ADD_DATA_FILE = "Add Data File";
	public static final String REMOVE_DATA_FILE = "Remove Data File";
	public static final String SET_AS_COLORKEY = "Set as Color Key";
	public static final String SET_AS_LINKKEY = "Set as Link Key";
	public static final String[] BUTTONS = {ADD_DATA_FILE, REMOVE_DATA_FILE, SET_AS_COLORKEY, SET_AS_LINKKEY};
	public static final String[][] LINKERS = {
			//TODO - Rohit: Removed Sample from first Linker. Confirm with Nathan if this is okay.
		{"IndividualID", "ID", "IID", "UID", "UniqueID", "IndID"},
		{"Family ID", "FamID", "FID"}, 
		{"DNA/Sample", "DNA", "DNA#", "Sample", "LabID"}, 
		{"MarkerName", "Marker", "SNP", "Variant", "VariantName"}, // will link to Scatter Plot
		{"Region", "UCSC", "Band", "Arm"},	// will link to Trailer
		{"Chromosome", "Chr"},	// secondary link to Trailer
		{"Position", "Pos", "Start", "Begin"}, // secondary link to Trailer
		{"Stop Position", "Stop", "End"} // secondary link to Trailer
	};
	public static final int IID_INDEX_IN_LINKERS = 0;
	public static final int FID_INDEX_IN_LINKERS = 1;
	public static final int DNA_INDEX_IN_LINKERS = 2;
	public static int MARKER_INDEX_IN_LINKERS = 3;
	public static int REGION_INDEX_IN_LINKERS = 4;
	public static int CHR_INDEX_IN_LINKERS = 5;
	public static int POS_INDEX_IN_LINKERS = 6;
	public static int STOP_POS_INDEX_IN_LINKERS = 7;

	//	private JPanel bottomPanel;
	private TwoDPanel twoDPanel;
	private JLayeredPane layeredPane;
	private ColorKeyPanel colorKeyPanel;

	private Project proj;
//	private int currentClass;
//	private int plot_type;
	private byte size;
	private SampleData sampleData;
//	private boolean maskMissing;
//	private JRadioButton[] typeRadioButtons;
//	private JRadioButton[] classRadioButtons;
	private JButton flipButton, invXButton, invYButton;
	private boolean flipStatus, xInvStatus, yInvStatus;
	private CheckBoxTree tree;
	private Vector<String> treeFilenameLookup;
//	Hashtable<String, String[][]> dataHash;
	Hashtable<String, Vector<String[]>> dataHash;
//	Hashtable<String, Hashtable<String, String[]>> dataHash;
	Hashtable<String, String[]> namesHash;
	Hashtable<String, boolean[]> numericHash;
	String[][] treeFileVariableNameLookup;
	Hashtable<String, int[]> keyIndices;
	Hashtable<String, Integer> linkKeyIndex;
	Hashtable<String, Hashtable<String, String[]>> linkKeyToDataHash;
	Hashtable<String, ArrayList<Integer>> colorKeyIndex;
//	Vector<String> colorKeyVariables;
//	String[][] colorKeyUniqueValues;
	Logger log;
//	Vector<Integer> rowsSelected;
//	Vector<String[]> linkKeyValues;
	
	public TwoDPlot() {
		this(null, new Logger());
	}

	public TwoDPlot(Project project, Logger log) {
		String[] previouslyLoadedFiles;
		
		this.log = log;
		proj = project;
		size = DEFAULT_SIZE;

		sampleData = proj.getSampleData(2, false);
		treeFilenameLookup = new Vector<String>();
		dataHash = new Hashtable<String, Vector<String[]>>();
		namesHash = new Hashtable<String, String[]>();
		numericHash = new Hashtable<String, boolean[]>();
		keyIndices = new Hashtable<String, int[]>();
		linkKeyIndex = new Hashtable<String, Integer>();
		linkKeyToDataHash = new Hashtable<String, Hashtable<String, String[]>>();
		colorKeyIndex = new Hashtable<String, ArrayList<Integer>>();

		previouslyLoadedFiles = proj.getFilenames(Project.TWOD_LOADED_FILENAMES);
		for (int i = 0; i < previouslyLoadedFiles.length; i++) {
			loadFile(previouslyLoadedFiles[i]);
		}
		
		
		setLayout(new BorderLayout());
		
//		SpringLayout layout = new SpringLayout();
//        setLayout(layout);

		twoDPanel = new TwoDPanel(this, log);
//		twoDPanel.setBounds(0,0,1000,600);
//		twoDPanel.setPreferredSize(new Dimension(1000, 600));

//		UIManager.put("PopupMenuUI", "CustomPopupMenuUI");

		layeredPane = new JLayeredPane();
		layeredPane.setLayout(new BorderLayout());
//		SpringLayout layout = new SpringLayout();
//		layeredPane.setLayout(layout);

		generateFlipButton();
		layeredPane.add(flipButton);
		generateInvXButton();
		layeredPane.add(invXButton);
		generateInvYButton();
		layeredPane.add(invYButton);
		
		layeredPane.add(twoDPanel);
		layeredPane.setPreferredSize(new Dimension(1000, 600));

		JPanel treePanel = new JPanel();
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

		initializeTree();
		updateTree();

		treePanel.add(new JScrollPane(tree), BorderLayout.CENTER);
		tree.addTreeSelectionListener(this);
		JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, treePanel, layeredPane);
		splitPane.setBackground(Color.WHITE);
		splitPane.setOneTouchExpandable(true);
		splitPane.setDividerLocation(150);

		Dimension minimumSize = new Dimension(100, 50);
		layeredPane.setMinimumSize(minimumSize);

		add(splitPane, BorderLayout.CENTER);


//		add(twoDPanel, BorderLayout.CENTER);
//		add(layeredPane, BorderLayout.CENTER);
		
		colorKeyPanel = new ColorKeyPanel(sampleData, twoDPanel);
		add(colorKeyPanel, BorderLayout.SOUTH);

		inputMapAndActionMap();

		twoDPanel.setPointsGeneratable(true);//zx
		twoDPanel.setExtraLayersVisible(new byte[] {99});
		updateGUI();
		
		twoDPanel.grabFocus();
		
		layeredPane.addComponentListener(new ComponentListener() {
			public void componentShown(ComponentEvent e) {}
			public void componentResized(ComponentEvent e) {
//				layeredPane.setSize(new Dimension(getWidth()-20, getHeight()-50));	//???zx
				flipButton.setBounds(70,layeredPane.getHeight()-75,38,38);
				invXButton.setBounds(70,layeredPane.getHeight()-35,38,13);
				invYButton.setBounds(55,layeredPane.getHeight()-75,13,38);
//				layeredPane.repaint();
//				twoDPanel.paintAgain();
				
			}
			public void componentMoved(ComponentEvent e) {}
			public void componentHidden(ComponentEvent e) {}
		});

		setVisible(true);
	}
	
	public void refreshOtherButtons() {
//		twoDPanel.repaint();
		flipButton.repaint();
		invXButton.repaint();
		invYButton.repaint();
	}

	private void generateFlipButton() {
		flipButton = new JButton(Grafik.getImageIcon("images/flip_and_invert/flip_10p.jpg", false));
		flipButton.setRolloverIcon(Grafik.getImageIcon("images/flip_and_invert/flip_10p_blue.jpg", false));
		flipButton.setToolTipText("Inverts axes");
		flipButton.setBorder(null);
		flipButton.setVisible(true);
		flipStatus = true;
		flipButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
//				twoDPanel.setSwapable(flipStatus);
				twoDPanel.setPointsGeneratable(true);
				twoDPanel.setSwapAxes(flipStatus);
				twoDPanel.paintAgain();
				if (flipStatus) {
					flipStatus = false;
				} else {
					flipStatus = true;
				}
			}
		});
	}

	private void generateInvXButton() {
		invXButton = new JButton(Grafik.getImageIcon("images/flip_and_invert/right_10.gif", true));
		invXButton.setBorder(null);
		invXButton.setVisible(true);
		xInvStatus = true;
		invXButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				twoDPanel.setPointsGeneratable(true);
				twoDPanel.setXinversion(xInvStatus);
				twoDPanel.paintAgain();
				if (xInvStatus) {
					invXButton.setIcon(Grafik.getImageIcon("images/flip_and_invert/left_10.gif", true));
				} else {
					invXButton.setIcon(Grafik.getImageIcon("images/flip_and_invert/right_10.gif", true));
				}
				xInvStatus = !xInvStatus;
			}
		});
	}

	private void generateInvYButton() {
		invYButton = new JButton(Grafik.getImageIcon("images/flip_and_invert/up_10.gif", true));
		invYButton.setBorder(null);
		invYButton.setVisible(true);
		yInvStatus = true;
		invYButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				twoDPanel.setPointsGeneratable(true);
				twoDPanel.setYinversion(yInvStatus);
				twoDPanel.paintAgain();
				if (yInvStatus) {
					invYButton.setIcon(Grafik.getImageIcon("images/flip_and_invert/down_10.gif", true));
				} else {
					invYButton.setIcon(Grafik.getImageIcon("images/flip_and_invert/up_10.gif", true));
				}
				yInvStatus = !yInvStatus;
			}
		});
	}

	private void inputMapAndActionMap() {
		InputMap inputMap = twoDPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.ALT_MASK), ALT_UP);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.ALT_MASK), ALT_DOWN);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), ALT_LEFT);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), ALT_RIGHT);
		ActionMap actionMap = twoDPanel.getActionMap();
		twoDPanel.setActionMap(actionMap);
	}

	private JMenuBar menuBar() {
		JMenuBar menuBar;
		JMenu menu;
		JMenuItem menuItemExit, menuItemOpen, menuItemSave;
		
		menuBar = new JMenuBar();
		menu = new JMenu("File");
        menu.setMnemonic(KeyEvent.VK_F);
		menuBar.add(menu);
		menuItemOpen = new JMenuItem("Open", KeyEvent.VK_O);
		menuItemOpen.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fileChooser = new JFileChooser(proj.getProjectDir());
				int fileOpenActionSelected = fileChooser.showOpenDialog(null);
		        if (fileOpenActionSelected == JFileChooser.APPROVE_OPTION) {
		    		loadFile(ext.replaceAllWith(fileChooser.getSelectedFile().toString(), "\\", "/"));
		    		updateTree();
		        }
			}
		});
		menu.add(menuItemOpen);
		menuItemSave = new JMenuItem("Save Image", KeyEvent.VK_S);
		menuItemSave.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fileChooser = new JFileChooser();
				int fileOpenActionSelected = fileChooser.showOpenDialog(null);
		        if (fileOpenActionSelected == JFileChooser.APPROVE_OPTION) {
		            File fileToOpen = fileChooser.getSelectedFile();
		            twoDPanel.screenCapture(fileToOpen.toString()+".png");	//??? zx: How to avoid twoDPanel being static?
		        }
			}
		});
		menu.add(menuItemSave);
		menuItemExit = new JMenuItem("Close", KeyEvent.VK_C);
		menuItemExit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				System.exit(0);
			}
		});
		menu.add(menuItemExit);
		menu = new JMenu("Edit");
        menu.setMnemonic(KeyEvent.VK_E);
		menuBar.add(menu);
		menu.add(new JMenuItem("Cut"));
		menu.add(new JMenuItem("Copy"));
		menu.add(new JMenuItem("Paste"));
		menu.add(new JMenuItem("Paste Image"));
		menu.add(new JMenuItem("Find"));
		menu = new JMenu("Help");
		menuBar.add(menu);
		menu.add(new JMenuItem("Contents"));
		menu.add(new JMenuItem("Search"));
		menu.add(new JMenuItem("About"));
		return menuBar;
	}
	
	public void actionPerformed(ActionEvent ae) {
		boolean found = false;
		String command = ae.getActionCommand();
		byte numberOfSelectedNodes;
		String[] keys;

		if (command.equals(ADD_DATA_FILE)) {
			JFileChooser fileChooser = new JFileChooser(proj.getProjectDir());
			int fileOpenActionSelected = fileChooser.showOpenDialog(null);
	        if (fileOpenActionSelected == JFileChooser.APPROVE_OPTION) {
        		for (int i=0; tree!=null && i<tree.getModel().getChildCount(tree.getModel().getRoot()); i++) {
        			if (ext.replaceAllWith(fileChooser.getSelectedFile().toString(), "\\", "/").equals(tree.getModel().getChild(tree.getModel().getRoot(),i).toString())) {
        				found = true;
        				break;
        			}
        		}
	        	if (found) {
	        		JOptionPane.showMessageDialog(null,"The data file has already been opened.");
	        	} else {
	        		loadFile(ext.replaceAllWith(fileChooser.getSelectedFile().toString(), "\\", "/"));
	        		updateTree();
	        	}
	        }
		} else if (command.equals(REMOVE_DATA_FILE)) {
			numberOfSelectedNodes = (byte) tree.getSelectedPathComponent();
			if (numberOfSelectedNodes != -1) {
				keys = HashVec.getKeys(dataHash); // it is better for keys to be a local variable rather than a global variable, otherwise it needs to be updated every time something is added or deleted
				tree.deleteSelectedNode();
				dataHash.remove(keys[numberOfSelectedNodes]);//TODO tree.getSelectionValues()[0][0] is not the branch to delete.
				namesHash.remove(keys[numberOfSelectedNodes]);
				keyIndices.remove(keys[numberOfSelectedNodes]);
				numericHash.remove(keys[numberOfSelectedNodes]);
				treeFilenameLookup.remove(keys[numberOfSelectedNodes]);
			}
//			System.out.println(dataHash.size()+"\t"+namesHash.size()+"\t"+numericHash.size()+"\t"+treeFilenameLookup.size());
//			System.out.println("\n==== End =====\n");
		} else if (command.equals(SET_AS_COLORKEY)) {
//			colorKeyVariables.add(tree.getSelectedPathComponentName());
//			colorKeyVariables.add(getNamesSelected()[0]);
			System.out.println("getSelectedPathComponent: " + tree.getSelectedPathComponent() + "\t + getNamesSelected: " + Arrays.toString(getNamesSelected()));
			int[] selectedLinkKey = tree.getSelectionRows();
			setColorKey(selectedLinkKey);
		} else if (command.equals(SET_AS_LINKKEY)) {
			int[] selectedLinkKey = tree.getSelectionRows();
			setLinkKey(selectedLinkKey);
		} else {
			System.err.println("Error - unknown command '"+command+"'");
		}
	}

	public void setColorKey(int[] selectedColorKey) {
		String[][] selectedNodes = tree.getSelectionValues();
		ArrayList<Integer> colorKeys;
		if (selectedColorKey.length == 1) {
			if (colorKeyIndex.containsKey(selectedNodes[0][0])) {
				colorKeys = colorKeyIndex.get(selectedNodes[0][0]);
			} else {
				colorKeyIndex.put(selectedNodes[0][0], colorKeys = new ArrayList<Integer>());
			}
			for (int i = 0; i < colorKeys.size(); i++) {
				if (colorKeys.get(i) == selectedColorKey[0])
					return;
			}
			colorKeys.add(selectedColorKey[0]);
		}
	}

	public void setLinkKey(int[] selectedLinkKey) {
		String[][] selectedNodes = tree.getSelectionValues();
		int[] linkKeyColumnLabels;
		if (selectedLinkKey.length == 1) {
			if (keyIndices.containsKey(selectedNodes[0][0])) {
				linkKeyColumnLabels = keyIndices.get(selectedNodes[0][0]);
			} else {
				JOptionPane.showMessageDialog(null, "There was a problem in your selection. Please select again", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			System.out.println("LinkColumn: " + Arrays.toString(linkKeyColumnLabels));
			System.out.println("SelectedNode: " + Arrays.toString(selectedLinkKey));

			for (int i = 0; i < linkKeyColumnLabels.length; i++) {
				if ((linkKeyColumnLabels[i] + 1) == selectedLinkKey[0]) {
					linkKeyIndex.put(selectedNodes[0][0], i);
					System.out.println("Link Key set to: " + Arrays.toString(LINKERS[i]));
					createLinkKeyToDataHash(selectedNodes[0][0], linkKeyColumnLabels);
					JOptionPane.showMessageDialog(null, "Link is set to: " + Arrays.toString(LINKERS[i]), "Information", JOptionPane.INFORMATION_MESSAGE);
					return;
				}
			}
			JOptionPane.showMessageDialog(null, "Unable to set link key. Please make sure you are selecting valid key", "Error", JOptionPane.ERROR_MESSAGE);
		}
	}

	public void createLinkKeyToDataHash(String filename, int[] linkKeyColumnLabels) {
		String inLine;
		String[] inLineArray;
		Hashtable<String, String[]> curFileLinkKeyDataHash = new Hashtable<String, String[]>();
		try {
			BufferedReader fileReader = new BufferedReader(new FileReader(filename));
			inLine = fileReader.readLine(); // headers: skip this line

			System.out.println("link column labels:" + Arrays.toString(linkKeyColumnLabels));
			System.out.println("link key index: " + linkKeyIndex.get(filename));
			while (fileReader.ready()) {
				if (inLine.contains("\t")) {
					inLineArray = fileReader.readLine().trim().split("\t", -1);
				} else {
					inLineArray = fileReader.readLine().trim().split("[\\s]+");
				}
				if (linkKeyIndex.containsKey(filename)) {
					switch (linkKeyIndex.get(filename)) {
					case DNA_INDEX_IN_LINKERS:
						curFileLinkKeyDataHash.put(inLineArray[linkKeyColumnLabels[linkKeyIndex.get(filename)]], inLineArray);
						break;
					case FID_INDEX_IN_LINKERS:
						curFileLinkKeyDataHash.put(inLineArray[linkKeyColumnLabels[linkKeyIndex.get(filename)]] + "\t" + inLineArray[linkKeyColumnLabels[IID_INDEX_IN_LINKERS]], inLineArray);
						break;
					case IID_INDEX_IN_LINKERS:
						curFileLinkKeyDataHash.put(inLineArray[linkKeyColumnLabels[linkKeyIndex.get(filename)]], inLineArray);
						break;
					default:
						System.out.println("Error: Invalid link key.");
						break;
					}
				}
			}
			fileReader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
		}
		if (!curFileLinkKeyDataHash.isEmpty())
			linkKeyToDataHash.put(filename, curFileLinkKeyDataHash);
		System.out.println("LinkKeyDataHash:" + linkKeyToDataHash.toString());
	}

	public void initLinkKey(int[] linkKeyColumnLabels, String filename) {
		if (linkKeyColumnLabels[DNA_INDEX_IN_LINKERS] >= 0) {
			// {"DNA/Sample", "DNA", "DNA#", "Sample", "LabID"} exists
			linkKeyIndex.put(filename, DNA_INDEX_IN_LINKERS);
			JOptionPane.showMessageDialog(null, "Link is set to: " + Arrays.toString(LINKERS[DNA_INDEX_IN_LINKERS]), "Information", JOptionPane.INFORMATION_MESSAGE);
			System.out.println("Link key set to: " + Arrays.toString(LINKERS[DNA_INDEX_IN_LINKERS]));
		} else if (linkKeyColumnLabels[FID_INDEX_IN_LINKERS] >= 0) {
			linkKeyIndex.put(filename, FID_INDEX_IN_LINKERS);
			JOptionPane.showMessageDialog(null, "Link is set to: " + Arrays.toString(LINKERS[FID_INDEX_IN_LINKERS]), "Information", JOptionPane.INFORMATION_MESSAGE);
			System.out.println("Link key set to: " + Arrays.toString(LINKERS[FID_INDEX_IN_LINKERS]));
		} else if (linkKeyColumnLabels[IID_INDEX_IN_LINKERS] >= 0) {
			linkKeyIndex.put(filename, IID_INDEX_IN_LINKERS);
			JOptionPane.showMessageDialog(null, "Link is set to: " + Arrays.toString(LINKERS[IID_INDEX_IN_LINKERS]), "Information", JOptionPane.INFORMATION_MESSAGE);
			System.out.println("Link key set to: " + Arrays.toString(LINKERS[IID_INDEX_IN_LINKERS]));
		} else {
			JOptionPane.showMessageDialog(null, "Unable to initialize the link key. Please select a link key manually.", "Error", JOptionPane.ERROR_MESSAGE);
			System.out.println("Unable to initialize the link key.");
		}
	}

	public byte getPointSize() {
		return size;
	}

	public String[] getNamesSelected() {
		String[] result = new String[2];
		String[][] selectionValues = tree.getSelectionValues();
		
		if (selectionValues[0][0] == null || selectionValues[1][0] == null) {
			result[0] = "";
			result[1] = "";
		} else {
			result[0] = ext.removeDirectoryInfo(selectionValues[0][0]) + " _ " + namesHash.get(selectionValues[0][0])[Integer.parseInt(selectionValues[0][1])];
			result[1] = ext.removeDirectoryInfo(selectionValues[1][0]) + " _ " + namesHash.get(selectionValues[1][0])[Integer.parseInt(selectionValues[1][1])];
		}
		return result;
	}
	
	public int[] getCurrentLinkKeyColumnLabels() {
		String[][] selectedValues;
		
		selectedValues = tree.getSelectionValues();
		if (selectedValues == null ||  selectedValues[0][0] == null) {
			return null;
		}
		return keyIndices.get(selectedValues[0][0]);
	}

//	public Vector<String[]> getCurrentLinkKeyValues() {
//		return linkKeyValues;
//	}

//	public String[][] getCurrentLinkKeyValues() {
//		String[][] linkKeyValues;
//		int[] linkKeyColumnLabels;
//		String[][] selectedValues;
//		Vector<String[]> currentData;
//		
//		selectedValues = tree.getSelectionValues();
//		if (selectedValues[0][0] == null ||  selectedValues[0][1] == null) {
//			return null;
//		}
//		currentData = dataHash.get(selectedValues[0][0]);
//		linkKeyColumnLabels = keyIndices.get(selectedValues[0][0]);
//		linkKeyValues = new String[rowsSelected.size()][linkKeyColumnLabels.length];	//TODO Problem is here.
//		for (int i=0; i<linkKeyValues.length; i++) {
//			for (int j=0; j<linkKeyColumnLabels.length; j++) {
//				if (linkKeyColumnLabels[j] >= 0) {
//					linkKeyValues[i][j] = currentData.elementAt(rowsSelected.elementAt(i))[linkKeyColumnLabels[j]];
//				}
//			}
//		}
//		return linkKeyValues;
//	}


//	public String[][] getCurrentLinkKeyValues() {
//		String[][] linkKeyValues;
//		int[] linkKeyColumnLabels;
//		Vector<String[]> dataOfSelectedFile;
//		String[][] selectedValues;
//		
//		selectedValues = tree.getSelectionValues();
//		if (selectedValues[0][0] == null ||  selectedValues[0][1] == null) {
//			return null;
//		}
//		linkKeyColumnLabels = keyIndices.get(selectedValues[0][0]);
//		dataOfSelectedFile = dataHash.get(tree.getSelectionValues()[0][0]);
//		linkKeyValues = new String[dataOfSelectedFile.size()][linkKeyColumnLabels.length];	//TODO Problem is here.
//		for (int i=0; i<linkKeyValues.length; i++) {
//			for (int j=0; j<linkKeyColumnLabels.length; j++) {
//				if (linkKeyColumnLabels[j] >= 0) {
//					linkKeyValues[i][j] = dataOfSelectedFile.elementAt(i)[linkKeyColumnLabels[j]];
//				}
//			}
//		}
//		return linkKeyValues;
//	}

	public Project getProject() {
		return proj;
	}

//	public int[] getCurrentDataLabels() {
//		int currentClass = colorKeyPanel.getCurrentClass();
//		int currentDataLabels;
//
//		currentDataLabels = new int[currentData.length];
//		for (int i=0; )
//	}



//	public float[][] getDataSelected() {
//		float[][] result;
//		String[][] selectedNodes = tree.getSelectionValues();
//
////		System.out.println("tree.getSelectionValues():\t"+temp[0][0]+", "+temp[0][1]+", "+temp[1][0]+", "+temp[1][1]);
//
//		if (selectedNodes[0][0] == null || selectedNodes[1][0] == null) {
//			result = new float[0][0];
////			result = null;
//		} else {
//			result = new float[dataHash.get(selectedNodes[0][0]).size()][2];
//			for (int i=0; i<result.length; i++) {
//				result[i][0] = Float.valueOf( dataHash.get(selectedNodes[0][0]).elementAt(i)[Integer.parseInt(selectedNodes[0][1])]);
//				result[i][1] = Float.valueOf( dataHash.get(selectedNodes[1][0]).elementAt(i)[Integer.parseInt(selectedNodes[1][1])]);
//			}
//		}
//		return result;
//	}



//	public float[][] getDataSelected() {
//		float[][] result;
//		String[][] selectedNodes;
//		Vector<String[]> dataOfSelectedFile;
//		Hashtable<String, String> xHash, yHash;
//		String[] line;
//		int selectedColumn;
//		String[] keys;
//		Vector<String[]> v;
//
//		selectedNodes = tree.getSelectionValues();
//		if (selectedNodes[0][0] == null || selectedNodes[1][0] == null) {
//			result = new float[0][0];
////			result = null;
//		} else {
//			selectedColumn = Integer.parseInt(selectedNodes[0][1]);
//			dataOfSelectedFile = dataHash.get(selectedNodes[0][0]);
//			if (keyIndices.get(selectedNodes[0][0]) == null) {
//				System.out.println("Please set link key for the file '" + selectedNodes[0][0] + "'. Alternatively, you could change the corresponding column lable to IID or ID and restart the program.");
//				return new float[0][0];
//			}
//			xHash = new Hashtable<String, String>();
//			for (int i=0; i<dataOfSelectedFile.size(); i++) {
//				line = dataOfSelectedFile.elementAt(i);
//				xHash.put(line[keyIndices.get(selectedNodes[0][0])[0]], line[selectedColumn]);
//			}
//
//			selectedColumn = Integer.parseInt(selectedNodes[1][1]);
//			dataOfSelectedFile = dataHash.get(selectedNodes[1][0]);
//			if (keyIndices.get(selectedNodes[1][0]) == null) {
//				System.out.println("Please set link key for the file '" + selectedNodes[1][0] + "'. Alternatively, you could change the corresponding column lable to IID or ID and restart the program.");
//				return new float[0][0];
//			}
//			yHash = new Hashtable<String, String>();
//			for (int i=0; i<dataOfSelectedFile.size(); i++) {
//				line = dataOfSelectedFile.elementAt(i);
//				yHash.put(line[keyIndices.get(selectedNodes[1][0])[0]], line[selectedColumn]);
//			}
//
//			keys = HashVec.getKeys(xHash, false, false);
//			v = new Vector<String[]>();
//			for (int i=0; i<keys.length; i++) {
//				if (yHash.containsKey(keys[i])) {
//					v.add(new String[] {keys[i], xHash.get(keys[i]), yHash.get(keys[i])});
//				}
//			}
//			result = new float[v.size()][2];
//			for (int i=0; i<result.length; i++) {
//				result[i][0] = Float.parseFloat(v.elementAt(i)[1]);
//				result[i][1] = Float.parseFloat(v.elementAt(i)[2]);
//			}
//		}
//		return result;
//	}



	public Vector<String[]> getDataSelected() {
		return getDataSelected(false);
	}

	public Vector<String[]> getDataSelected(boolean includeColorKeyValue) {
//		String[][] result;
		String[][] selectedNodes;
		Vector<String[]> dataOfSelectedFile;
		Hashtable<String, String[]> xHash;
		Hashtable<String, String> yHash;
		String[] inLine, outLine;
		int selectedColumn;
		String[] keys;
		Vector<String[]> v;
		int currentClass;
		String[] ids;
		byte colorCode;
		int[] linkKeyColumnLabels;
		byte index;

		selectedNodes = tree.getSelectionValues();
		v = new Vector<String[]>();
		if (selectedNodes[0][0] != null && selectedNodes[1][0] != null && keyIndices.get(selectedNodes[0][0]) != null && keyIndices.get(selectedNodes[1][0]) != null) {
//			rowsSelected = new Vector<Integer>();
//			linkKeyValues = new Vector<String[]>();
			selectedColumn = Integer.parseInt(selectedNodes[0][1]);
			dataOfSelectedFile = dataHash.get(selectedNodes[0][0]);
			currentClass = colorKeyPanel.getCurrentClass();
			if (currentClass < SampleData.BASIC_CLASSES.length && SampleData.BASIC_CLASSES[currentClass].equals(SampleData.HEATMAP)) {
				twoDPanel.setChartType(AbstractPanel.HEAT_MAP_TYPE);
			} else {
				twoDPanel.setChartType(AbstractPanel.SCATTER_PLOT_TYPE);
			}
			
			linkKeyColumnLabels = keyIndices.get(selectedNodes[0][0]);
			index = (byte) (includeColorKeyValue? 4 : 3);
			xHash = new Hashtable<String, String[]>();
			for (int i=0; i<dataOfSelectedFile.size(); i++) {
				inLine = dataOfSelectedFile.elementAt(i);
				outLine = new String[linkKeyColumnLabels.length + index];
				outLine[0] = inLine[0];
				outLine[1] = inLine[selectedColumn];
				for (int j=0; j<linkKeyColumnLabels.length; j++) {
					if (linkKeyColumnLabels[j] >= 0) {
						outLine[j + index] = inLine[linkKeyColumnLabels[j]];
					}
				}
				xHash.put(inLine[0], outLine);
			}

			selectedColumn = Integer.parseInt(selectedNodes[1][1]);
			dataOfSelectedFile = dataHash.get(selectedNodes[1][0]);
			yHash = new Hashtable<String, String>();
			for (int i=0; i<dataOfSelectedFile.size(); i++) {
				inLine = dataOfSelectedFile.elementAt(i);
				yHash.put(inLine[keyIndices.get(selectedNodes[1][0])[0]], inLine[selectedColumn]);
			}

			keys = HashVec.getKeys(xHash, false, false);
			v = new Vector<String[]>();
			if (includeColorKeyValue) {
				for (int i=0; i<keys.length; i++) {
					if (yHash.containsKey(keys[i])) {
						ids = sampleData.lookup(keys[i]);
						if (ids == null) {
							colorCode = 0;
						} else {
							colorCode = sampleData.determineCodeFromClass(currentClass, (byte)0, sampleData.getIndiFromSampleHash(ids[0]), (byte)0, 0);
						}
						inLine = xHash.get(keys[i]);
						inLine[2] = yHash.get(keys[i]);
						inLine[3] = colorCode + "";
//						v.add(new String[] {keys[i], line[0], yHash.get(keys[i]), colorCode + ""});
						v.add(inLine);
					}
				}

			} else {
				for (int i=0; i<keys.length; i++) {
					if (yHash.containsKey(keys[i])) {
						inLine = xHash.get(keys[i]);
						inLine[2] = yHash.get(keys[i]);
//						v.add(new String[] {keys[i], xHash.get(keys[i]), yHash.get(keys[i]), "all other keys"});
						v.add(inLine);
					}
				}
			}
		}
//		return result;
		return v;
	}


//	public float[][] getDataSelected() {
//		return getDataSelected(false);
//	}

//	public float[][] getDataSelected(boolean includeColorKeyValue) {
//		float[][] result;
//		String[][] selectedNodes;
//		Vector<String[]> dataOfSelectedFile;
//		Hashtable<String, String> xHash, yHash;
//		String[] line;
//		int selectedColumn;
//		String[] keys;
//		Vector<String[]> v;
//		int currentClass;
//		String[] ids;
//		byte colorCode;
//
//		selectedNodes = tree.getSelectionValues();
//		if (selectedNodes[0][0] == null || selectedNodes[1][0] == null) {
//			result = new float[0][0];
////			result = null;
//		} else {
//			selectedColumn = Integer.parseInt(selectedNodes[0][1]);
//			dataOfSelectedFile = dataHash.get(selectedNodes[0][0]);
//			currentClass = colorKeyPanel.getCurrentClass();
//			if (keyIndices.get(selectedNodes[0][0]) == null) {
//				System.out.println("Please set link key for the file '" + selectedNodes[0][0] + "'. Alternatively, you could change the corresponding column lable to IID or ID and restart the program.");
//				return new float[0][0];
//			}
//			xHash = new Hashtable<String, String>();
//			for (int i=0; i<dataOfSelectedFile.size(); i++) {
//				line = dataOfSelectedFile.elementAt(i);
//				xHash.put(line[keyIndices.get(selectedNodes[0][0])[0]], line[selectedColumn]);
//			}
//
//			selectedColumn = Integer.parseInt(selectedNodes[1][1]);
//			dataOfSelectedFile = dataHash.get(selectedNodes[1][0]);
//			if (keyIndices.get(selectedNodes[1][0]) == null) {
//				System.out.println("Please set link key for the file '" + selectedNodes[1][0] + "'. Alternatively, you could change the corresponding column lable to IID or ID and restart the program.");
//				return new float[0][0];
//			}
//			yHash = new Hashtable<String, String>();
//			for (int i=0; i<dataOfSelectedFile.size(); i++) {
//				line = dataOfSelectedFile.elementAt(i);
//				yHash.put(line[keyIndices.get(selectedNodes[1][0])[0]], line[selectedColumn]);
//			}
//
//			keys = HashVec.getKeys(xHash, false, false);
//			v = new Vector<String[]>();
//			if (includeColorKeyValue) {
//				for (int i=0; i<keys.length; i++) {
//					if (yHash.containsKey(keys[i])) {
//							ids = sampleData.lookup(keys[i]);
//							if (ids == null) {
//								colorCode = 0;
//							} else {
//								colorCode = sampleData.determineCodeFromClass(currentClass, (byte)0, sampleData.getIndiFromSampleHash(ids[0]), (byte)0, 0);
//							}
//							v.add(new String[] {keys[i], xHash.get(keys[i]), yHash.get(keys[i]), colorCode + ""});
//					}
//				}
//				result = new float[v.size()][3];
//				for (int i=0; i<result.length; i++) {
//					result[i][0] = Float.parseFloat(v.elementAt(i)[1]);
//					result[i][1] = Float.parseFloat(v.elementAt(i)[2]);
//					result[i][2] = Integer.parseInt(v.elementAt(i)[3]);
//				}
//
//			} else {
//				for (int i=0; i<keys.length; i++) {
//					if (yHash.containsKey(keys[i])) {
//						v.add(new String[] {keys[i], xHash.get(keys[i]), yHash.get(keys[i])});
//					}
//				}
//				result = new float[v.size()][2];
//				for (int i=0; i<result.length; i++) {
//					result[i][0] = Float.parseFloat(v.elementAt(i)[1]);
//					result[i][1] = Float.parseFloat(v.elementAt(i)[2]);
//				}
//			}
//		}
//		return result;
//	}


	//	public int[][] getCurrentPair() {
//		if (tree.getSelectionIndices()[0][0]<0 || tree.getSelectionIndices()[0][1]<0 || tree.getSelectionIndices()[1][0]<0 || tree.getSelectionIndices()[1][1]<0) {
//			return null;
//		} else {
//			return tree.getSelectionIndices();
//		}
//	}

//	public boolean maskMissing() {
//		return maskMissing;
//	}

	public void updateGUI() {
//		if (markerList.length==0) {
//			markerName.setText("Error: marker data was not successfully loaded");
//			commentLabel.setText("Check to make sure MarkerLookup is synchronized with the current data");
//			classPanel.setEnabled(false);
//		} else {
//			markerName.setText(markerList[markerIndex]);
//			commentLabel.setText(commentList[markerIndex]);
//
//			if (classPanel!=null) {
//				classPanel.setEnabled(markerData[markerIndex].getFingerprint()==sampleListFingerprint);
//			}
//		}
		
		
		
		/*
		if (plot_type >= 2) {
			symmetryBox.setEnabled(false);
			twoDPanel.setSymmetricAxes(false);
		} else {
			symmetryBox.setEnabled(true);
			twoDPanel.setSymmetricAxes(symmetryBox.isSelected());
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
		*/
		
		twoDPanel.paintAgain();
	}

	public void updateColorKey(Hashtable<String,String> hash) {
		colorKeyPanel.updateColorKey(hash);
	}


//	public void displayIndex(JTextField field) {
//		field.setText((markerIndex+1)+" of "+markerList.length);
//	}
	

	public void windowActivated(WindowEvent e) {}

	public void windowClosed(WindowEvent e) {}

	@Override
	public void windowClosing(WindowEvent e) {
		String[] options;
		int choice;
		String filenames, selections, message;

		filenames = Array.toStr(Array.toStringArray(treeFilenameLookup), ";");
		selections = "";
		
		options = new String[] {"Yes", "No", "Cancel"};
		message = "";
		
		if (!proj.getProperty(Project.TWOD_LOADED_FILENAMES).equals(filenames)) {
			if (filenames.equals("")) {
				message = "All files have been unloaded from 2D Plot.";
			} else {
				message = "A different set of files have been loaded/unloaded from 2D Plot.";
			}
		}
		
		if (!proj.getProperty(Project.TWOD_LOADED_VARIABLES).equals(selections)) {
			if (!message.equals("")) {
				message = message.substring(0, message.length()-1)+", and new variables have been selected.";
			} else {
				message = "New variables have been selected.";
			}
		}
		
		if (message.equals("")) {
			GuiManager.disposeOfParentFrame(this);
		} else {
			System.out.println("message: '"+message+"'");
			choice = JOptionPane.showOptionDialog(null, message+" Would you like to keep this configuration for the next time 2D Plot is loaded?", "Preserve 2D Plot workspace?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
			if (choice == 0) {
				proj.setProperty(Project.TWOD_LOADED_FILENAMES, filenames);
				proj.setProperty(Project.TWOD_LOADED_VARIABLES, selections);
				proj.saveProperties();
				GuiManager.disposeOfParentFrame(this);
			} else if (choice == 1) {
				GuiManager.disposeOfParentFrame(this);
			} else if (choice == -1 || choice == 2) {
				// keep frame alive
			}
		}
	}

	public void windowDeactivated(WindowEvent e) {}

	public void windowDeiconified(WindowEvent e) {}

	public void windowIconified(WindowEvent e) {}

	public void windowOpened(WindowEvent e) {}

	public void valueChanged(TreeSelectionEvent e) {
//		System.out.println("Tree value changed listener"+"\ttreeSelectionIndices: "+tree.getSelectionIndices()[0][0]+","+tree.getSelectionIndices()[0][1]+"\t"+tree.getSelectionIndices()[1][0]+","+tree.getSelectionIndices()[1][1]);
		twoDPanel.setPointsGeneratable(true);
		twoDPanel.createImage();
		twoDPanel.paintAgain();
	}

	public void loadFile(String filename) {
		BufferedReader reader;
		String[] header, line;
		String readBuffer;
		int[] linkKeyIndices;

		if (treeFilenameLookup.contains(filename)) {
			return;
		}
		try {
//	        reader = Files.getReader(new FileReader(filename), proj.getJarStatus(), true, false);
			reader = new BufferedReader(new FileReader(filename));
			treeFilenameLookup.add(filename);
//			filename = ext.removeDirectoryInfo(filename); // Strip "/" from filename

			readBuffer = reader.readLine();
			if (readBuffer.contains("\t")) {
				header = readBuffer.trim().split("\t",-1);
			} else {
				header = readBuffer.trim().split("[\\s]+");
			}
			namesHash.put(filename, header);
			
			
			linkKeyIndices = ext.indexFactors(LINKERS, header, false, true, false, log, false);
			
			if (linkKeyIndices[0] == -1) {
				log.report("ID linker not automatically identified for file '" + filename + "'; assuming the first column.");
				linkKeyIndices[0] = 0;
			}
			
			keyIndices.put(filename, linkKeyIndices);
			
			numericHash.put(filename, new boolean[namesHash.get(filename).length]);
        	for (int i=0; i<numericHash.get(filename).length; i++) {
        		numericHash.get(filename)[i] = true;
        	}

			initLinkKey(linkKeyIndices, filename);	// initialize the link key
			createLinkKeyToDataHash(filename, linkKeyIndices);
        	dataHash.put(filename, new Vector<String[]>());
            while (reader.ready()) {
				if (readBuffer.contains("\t")) {
					line = reader.readLine().trim().split("\t",-1);
				} else {
					line = reader.readLine().trim().split("[\\s]+");
				}
            	dataHash.get(filename).add(line);
            	for (int i=0; i<header.length; i++) {
            		if (!ext.isValidDouble(line[i])) {
            			numericHash.get(filename)[i] = false;
            		}
            	}
            }
            reader.close();
        } catch (FileNotFoundException fnfe) {
        	System.err.println("Error: file \""+filename+"\" not found in current directory");
        } catch (IOException ioe) {
            System.err.println("Error reading file \""+filename+"\"");
        }
	}

	private void initializeTree() {
		tree = new CheckBoxTree(new String[0], new String[0], new String[0][], new boolean[0], 2);
	}
	
	private void updateTree() {
		String[] namesOfBranches;
		String[] branchHandles;
		String[][] namesOfNodes;

		
//        TreePath[] prev;
//        
//    	prev = tree.getSelectionPaths();
//
        treeFileVariableNameLookup = new String[treeFilenameLookup.size()][];
        for (int i=0; i<treeFilenameLookup.size(); i++) {
    		treeFileVariableNameLookup[i] = namesHash.get(treeFilenameLookup.elementAt(i));
        	if (tree==null) {
        		namesOfBranches = new String[1];
        		branchHandles = new String[1];
        		namesOfNodes = new String[1][];
        		namesOfBranches[0]=ext.removeDirectoryInfo(treeFilenameLookup.elementAt(i));
        		branchHandles[0]=treeFilenameLookup.elementAt(i);
        		namesOfNodes[0] = treeFileVariableNameLookup[i];
        		tree = new CheckBoxTree(namesOfBranches, branchHandles, namesOfNodes, numericHash.get(treeFilenameLookup.elementAt(0)), 2);
        	} else {
//        		if(i>=tree.getModel().getChildCount(tree.getModel().getRoot())) {
//        			tree.addNode(treeFilenameLookup.elementAt(i), i, treeFileVariableNameLookup[i]);
//        		}
        		boolean found = false;
        		for (int j=0; j<tree.getModel().getChildCount(tree.getModel().getRoot()); j++) {
        			if(tree.getModel().getChild(tree.getModel().getRoot(), j).toString().equals(ext.removeDirectoryInfo(treeFilenameLookup.elementAt(i)))) {
        				found = true;
        			}
        		}
        		if (!found) {
        			tree.addNode(ext.removeDirectoryInfo(treeFilenameLookup.elementAt(i)), treeFilenameLookup.elementAt(i), treeFileVariableNameLookup[i], numericHash.get(treeFilenameLookup.elementAt(i)));
        		}
        	}
       }
	}

    /**
	     * Create the GUI and show it.  For thread safety,
	     * this method should be invoked from the
	     * event-dispatching thread.
	     */
	public static void createAndShowGUI(Project proj, Logger log) {

		//Create and set up the window.
		JFrame frame = new JFrame("2D Plot");
//        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);

        //Create and set up the content pane.
        TwoDPlot twoDPlot = new TwoDPlot(proj, log);
        frame.setJMenuBar(twoDPlot.menuBar());
        twoDPlot.setOpaque(true); //content panes must be opaque
        frame.setContentPane(twoDPlot);
        frame.addWindowListener(twoDPlot);
        frame.setBounds(20, 20, 1000, 600);
//		frame.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }

	public static void main(String[] args) {
		ext.verifyDirFormat("");

        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
//                createAndShowGUI(new Project("C:/workspace/Genvisis/projects/twodplot.properties", false), new Logger());
//                createAndShowGUI(new Project("C:/workspace/Genvisis/projects/GEDI_exome.properties", false), new Logger());
                createAndShowGUI(new Project(Project.DEFAULT_PROJECT, false), new Logger());
            }
        });
		
	}
}
