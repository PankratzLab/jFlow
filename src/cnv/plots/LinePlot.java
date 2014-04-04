package cnv.plots;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import javax.swing.*;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;

import cnv.filesys.Project;
import cnv.gui.CheckBoxTree;
import cnv.gui.ColorIcon;

import common.Grafik;
import common.HashVec;
import common.Logger;
import common.ext;

public class LinePlot extends JPanel implements WindowListener, ActionListener, TreeSelectionListener {
	public static final long serialVersionUID = 1L;
	public static final byte DEFAULT_SIZE = 8;

	public static final String LINE_PLOT_DIR_NAME = "LinePlot";
	public static final String COMMON_FILE_EXTENSION = ".xln";
	public static final String LINE_PLOT_INPUT_FILENAME = "input" + COMMON_FILE_EXTENSION;
	public static final String LINE_PLOT_COLOR_FILENAME = "ColorMap" + COMMON_FILE_EXTENSION;
	private static final String ALT_UP = "ALT UP";
	private static final String ALT_DOWN = "ALT DOWN";
	private static final String ALT_LEFT = "ALT LEFT";
	private static final String ALT_RIGHT = "ALT RIGHT";
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String ADD_DATA_FILE = "Add Data Files";
	public static final String REMOVE_DATA_FILE = "Remove Data File";
	public static final String SET_AS_COLORKEY = "Set as Color Key";
	public static final String SET_AS_LINKKEY = "Set as Link Key";
	public static final String[] BUTTONS = { ADD_DATA_FILE, REMOVE_DATA_FILE, SET_AS_COLORKEY, SET_AS_LINKKEY };
	public static final String[][] LINKERS = { { "IndividualID", "ID", "IID", "UID", "UniqueID", "IndID", "Sample" }, { "Family ID", "FamID", "FID" }, { "DNA/Sample", "DNA", "DNA#", "Sample", "LabID" }, { "MarkerName", "Marker", "SNP", "Variant", "VariantName" }, // will link to Scatter Plot
	{ "Region", "UCSC", "Band", "Arm" }, // will link to Trailer
	{ "Chromosome", "Chr" }, // secondary link to Trailer
	{ "Position", "Pos", "Start", "Begin" }, // secondary link to Trailer
	{ "Stop Position", "Stop", "End" } // secondary link to Trailer
	};

	private LinePanel linePanel;
	private JLayeredPane layeredPane;
	private JPanel legendPanel;
	private Project proj;
	private JButton flipButton, invXButton, invYButton;
	private boolean flipStatus, xInvStatus, yInvStatus;
	private CheckBoxTree tree;
	private Vector<String> treeFilenameLookup;
	Hashtable<String, Vector<String[]>> dataHash;
	Hashtable<String, Vector<String[]>> commentHash;
	Hashtable<String, String[]> namesHash;
	Hashtable<String, boolean[]> numericHash;
	String[][] treeFileVariableNameLookup;
	Hashtable<String, int[]> keyIndices;
	Hashtable<String, Byte> groupToColorHash;
	Hashtable<String, Byte> tempGroupToColorHash;
	Hashtable<String, Boolean> enabledGroups;
	Hashtable<String, Hashtable<String, String>> colorKeyHash;
	Hashtable<String, JLabel[]> groupToColorLabelHash;
	Logger log;

	public LinePlot() {
		this(null, new Logger());
	}

	public LinePlot(Project project, Logger log) {
		String[] previouslyLoadedFiles;
		this.log = log;
		proj = project;
		treeFilenameLookup = new Vector<String>();
		// TODO Need to save the previously loaded files in other location.
		previouslyLoadedFiles = proj.getFilenames(Project.TWOD_LOADED_FILENAMES);
		dataHash = new Hashtable<String, Vector<String[]>>();
		commentHash = new Hashtable<String, Vector<String[]>>();
		namesHash = new Hashtable<String, String[]>();
		numericHash = new Hashtable<String, boolean[]>();
		keyIndices = new Hashtable<String, int[]>();
		enabledGroups = new Hashtable<String, Boolean>();
		colorKeyHash = new Hashtable<String, Hashtable<String, String>>();
		tempGroupToColorHash = new Hashtable<String, Byte>();
		groupToColorHash = new Hashtable<String, Byte>();
		groupToColorLabelHash = new Hashtable<String, JLabel[]>();
		for (int i = 0; i < previouslyLoadedFiles.length; i++) {
			loadFile(previouslyLoadedFiles[i]);
		}
		setLayout(new BorderLayout());

		linePanel = new LinePanel(this, log);

		layeredPane = new JLayeredPane();
		layeredPane.setLayout(new BorderLayout());

		generateFlipButton();
		layeredPane.add(flipButton);
		generateInvXButton();
		layeredPane.add(invXButton);
		generateInvYButton();
		layeredPane.add(invYButton);

		layeredPane.add(linePanel);
		layeredPane.setPreferredSize(new Dimension(1000, 600)); // ???zx

		// ******* New code starts here ************
		JPanel treePanel = new JPanel();
		treePanel.setBackground(BACKGROUND_COLOR);
		treePanel.setLayout(new BorderLayout());

		JPanel buttonPanel = new JPanel();
		buttonPanel.setBackground(BACKGROUND_COLOR);
		buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.Y_AXIS));

		JButton button;
		for (int i = 0; i < BUTTONS.length; i++) {
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

		// Provide minimum sizes for the two components in the split pane
		Dimension minimumSize = new Dimension(100, 50);
		layeredPane.setMinimumSize(minimumSize);

		add(splitPane, BorderLayout.CENTER);

		legendPanel = new JPanel();

		add(legendPanel, BorderLayout.SOUTH);

		inputMapAndActionMap();

		linePanel.setPointsGeneratable(true);// zx
		linePanel.setExtraLayersVisible(new byte[] { 99 });
		updateGUI();

		linePanel.grabFocus();

		layeredPane.addComponentListener(new ComponentListener() {
			public void componentShown(ComponentEvent e) {
			}

			public void componentResized(ComponentEvent e) {
				flipButton.setBounds(70, layeredPane.getHeight() - 75, 38, 38);
				invXButton.setBounds(70, layeredPane.getHeight() - 35, 38, 13);
				invYButton.setBounds(55, layeredPane.getHeight() - 75, 13, 38);
			}

			public void componentMoved(ComponentEvent e) {
			}

			public void componentHidden(ComponentEvent e) {
			}
		});

		setVisible(true);
		readColorKeyFile(proj.getProjectDir() + File.separator + LINE_PLOT_DIR_NAME + File.separator + LINE_PLOT_COLOR_FILENAME);
		updateColorKey();
	}

	public void refreshOtherButtons() {
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
				linePanel.setPointsGeneratable(true);
				linePanel.setSwapAxes(flipStatus);
				linePanel.paintAgain();
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
				linePanel.setPointsGeneratable(true);
				linePanel.setXinversion(xInvStatus);
				linePanel.paintAgain();
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
				linePanel.setPointsGeneratable(true);
				linePanel.setYinversion(yInvStatus);
				linePanel.paintAgain();
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
		InputMap inputMap = linePanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.ALT_MASK), ALT_UP);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.ALT_MASK), ALT_DOWN);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), ALT_LEFT);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), ALT_RIGHT);
		ActionMap actionMap = linePanel.getActionMap();
		linePanel.setActionMap(actionMap);
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
					displayAll();
					tree.expandRow(0);
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
					linePanel.screenCapture(fileToOpen.toString() + ".png"); // ??? zx: How to avoid LinePanel being static?
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

	/**
	 * Function to get the color for a filename (celltype)
	 * 
	 * @param filename
	 *            the filename of the celltype
	 * @return the index of the color from DEFAULT_COLOR in {@link LinePanel} as byte else index of black color as byte
	 */
	public byte getColorFromFilename(String filename) {

		Set<String> keySet = colorKeyHash.keySet();
		Iterator<String> it = keySet.iterator();
		while (it.hasNext()) {
			String thisGroup = it.next();
			for (String thisFilename : colorKeyHash.get(thisGroup).keySet()) {
				if (thisFilename.equals(filename)) {
					if (groupToColorHash.containsKey(thisGroup)) {
						return (groupToColorHash.get(thisGroup));
					}
				}
			}
		}
		return (byte) LinePanel.DEFAULT_COLORS_BLACK_INDEX;
	}

	/**
	 * Function to read the color key file
	 * 
	 * @param filename
	 *            the color key filename
	 */
	public void readColorKeyFile(String filename) {
		String curLine;
		String[] curLineParams;
		Hashtable<String, String> colorKeyHashValue;

		try {
			BufferedReader reader = new BufferedReader(new FileReader(filename));
			// skip the headers in file
			reader.readLine();
			while (reader.ready()) {
				curLine = reader.readLine();
				curLineParams = curLine.trim().split("[\\s]+");
				if (colorKeyHash.containsKey(curLineParams[2])) { // if this group is already there
					colorKeyHashValue = colorKeyHash.get(curLineParams[2]);
				} else {
					colorKeyHash.put(curLineParams[2], colorKeyHashValue = new Hashtable<String, String>());
				}

				colorKeyHashValue.put(curLineParams[0], curLineParams[1]); // put filename as key and human readable
																			// celltype as value

			}
			reader.close();
		} catch (IOException e) {
			log.report("Unable to read color key file" + e.getMessage());
			System.exit(1);
		}
	}

	/**
	 * Function to display all the files by selecting all the checkboxes initially
	 */
	public void displayAll() {
		if (namesHash != null) {
			for (String thisProject : namesHash.keySet()) {
				String[] thisProjectFiles = namesHash.get(thisProject);
				for (int i = 0; i < thisProjectFiles.length; i++) {
					tree.performCheckBoxAction(thisProjectFiles[i], ItemEvent.SELECTED);
				}
			}
			for (String thisGroup : colorKeyHash.keySet()) {
				enabledGroups.put(thisGroup, false);
			}
		}
	}

	/**
	 * Function to perform a given action on a all the files (checkboxes) in a group
	 * 
	 * @param groupName
	 *            the name of the group
	 * @param action
	 *            the action to be performed
	 */
	public void performGroupCheckboxAction(String groupName, int action) {
		Hashtable<String, String> groupValuesHash = colorKeyHash.get(groupName);
		for (String thisFile : groupValuesHash.keySet()) {
			tree.performCheckBoxAction(thisFile, action);
		}
	}

	/**
	 * Function to draw and update the ColorKeyPanel and also defined the mouse click action
	 */
	public void updateColorKey() {
		JLabel label, block;
		String[] keys;
		MouseListener mouseListenerForColorKey;

		mouseListenerForColorKey = new MouseListener() {
			public void mouseReleased(MouseEvent e) {
			}

			public void mousePressed(MouseEvent e) {
			}

			public void mouseExited(MouseEvent e) {
			}

			public void mouseEntered(MouseEvent e) {
			}

			public void mouseClicked(MouseEvent e) {

				String classValue;

				if (e.getButton() == MouseEvent.BUTTON1) {

					classValue = ((JLabel) e.getSource()).getName();
					// if this checkbox is already selected
					if (enabledGroups.containsKey(classValue)) {
						// if selected to display as without color
						if (enabledGroups.get(classValue) == true) {
							enabledGroups.remove(classValue);
							// get the the color for this group from temporary hash
							if (tempGroupToColorHash.containsKey(classValue)) {
								groupToColorHash.put(classValue, tempGroupToColorHash.get(classValue));
							}
							// set the Group color block in ColorKeyPanel to white
							if (groupToColorLabelHash.containsKey(classValue)) {
								JLabel thisColorLabel[] = groupToColorLabelHash.get(classValue);
								thisColorLabel[0] = new JLabel(new ColorIcon(12, 12, LinePanel.DEFAULT_COLORS[LinePanel.DEFAULT_COLORS.length - 1]));
								generateColorKeyPanel();
							}
							performGroupCheckboxAction(classValue, ItemEvent.DESELECTED);
						} else { // if selected to show with color
							enabledGroups.put(classValue, true); // reverse the flag
							// Store this group to a temp hash
							if (groupToColorHash.containsKey(classValue)) {
								tempGroupToColorHash.put(classValue, groupToColorHash.get(classValue));
							}
							// set this group color to black
							groupToColorHash.put(classValue, (byte) LinePanel.DEFAULT_COLORS_BLACK_INDEX);
							// set this group color block in ColorKeyPanel to black
							if (groupToColorLabelHash.containsKey(classValue)) {
								JLabel thisColorLabel[] = groupToColorLabelHash.get(classValue);
								thisColorLabel[0] = new JLabel(new ColorIcon(12, 12, LinePanel.DEFAULT_COLORS[LinePanel.DEFAULT_COLORS_BLACK_INDEX]));
								generateColorKeyPanel();
							}
							performGroupCheckboxAction(classValue, ItemEvent.SELECTED);
						}
					} else { // if this checkbox is deselected
						enabledGroups.put(classValue, false); // add to enabledGroups
						// set the color block to this group's color in the ColorKeyPanel
						if (groupToColorLabelHash.containsKey(classValue)) {
							JLabel thisColorLabel[] = groupToColorLabelHash.get(classValue);
							thisColorLabel[0] = new JLabel(new ColorIcon(12, 12, LinePanel.DEFAULT_COLORS[groupToColorHash.get(classValue)]));
							generateColorKeyPanel();
						}
						performGroupCheckboxAction(classValue, ItemEvent.SELECTED);
					}
					log.report("Selected files for this group: " + enabledGroups.toString());
					linePanel.paintAgain();
				}
			}
		};

		// ColorKeyPanel
		label = new JLabel("Color key:");
		label.setFont(new Font("Arial", 0, 14));
		legendPanel.add(label);
		keys = HashVec.getKeys(colorKeyHash);
		for (int i = 0; i < keys.length; i++) {
			if (!keys[i].equals("-1")) {
				groupToColorHash.put(keys[i], (byte) (i + 1));
				// the color block
				block = new JLabel(new ColorIcon(12, 12, LinePanel.DEFAULT_COLORS[i + 1]));
				block.setName(keys[i] + "_ColorBlock");
				block.addMouseListener(mouseListenerForColorKey);

				// the color label
				label = new JLabel(keys[i]);
				label.setName(keys[i]);
				label.setFont(new Font("Arial", 0, 14));
				label.addMouseListener(mouseListenerForColorKey);

				groupToColorLabelHash.put(label.getName(), new JLabel[] { block, label });
			}
		}
		generateColorKeyPanel();
	}

	/**
	 * Function to generate the color blocks and their labels for ColorKeyPanel
	 */
	public void generateColorKeyPanel() {

		legendPanel.removeAll();
		legendPanel.repaint();

		for (String thisGroup : groupToColorLabelHash.keySet()) {
			legendPanel.add(groupToColorLabelHash.get(thisGroup)[0]); // color black
			legendPanel.add(groupToColorLabelHash.get(thisGroup)[1]); // color label
		}
		legendPanel.validate();
	}

	public void windowActivated(WindowEvent e) {
	}

	public void windowClosed(WindowEvent e) {
	}

	public void windowClosing(WindowEvent e) {
	}

	public void windowDeactivated(WindowEvent e) {
	}

	public void windowDeiconified(WindowEvent e) {
	}

	public void windowIconified(WindowEvent e) {
	}

	public void windowOpened(WindowEvent e) {
	}

	public void valueChanged(TreeSelectionEvent e) {
		linePanel.setPointsGeneratable(true);
		linePanel.createImage();
		linePanel.paintAgain();
	}

	/**
	 * Function to load the data for a given filename
	 * 
	 * @param filename
	 *            the filename
	 */
	public void loadFile(String filename) {

		BufferedReader reader;
		String[] line;
		String readBuffer;
		int[] linkKeyIndices;

		// if this file is already added
		if (treeFilenameLookup.contains(filename)) {
			return;
		}

		try {
			reader = new BufferedReader(new FileReader(filename));
			treeFilenameLookup.add(filename);
			readBuffer = reader.readLine();

			// Split the line on whitespaces and get all the headers in the file
			if (readBuffer.contains("\t")) {
				line = readBuffer.trim().split("\t", -1);
			} else {
				line = readBuffer.trim().split("[\\s]+");
			}

			// Adding headers of the file
			namesHash.put(filename, line);

			linkKeyIndices = ext.indexFactors(LINKERS, line, false, true, false, log, false);

			if (linkKeyIndices[0] == -1) {
				log.report("ID linker not automatically identified for file '" + filename + "'; assuming the first column.");
				linkKeyIndices[0] = 0;
			}

			keyIndices.put(filename, linkKeyIndices);

			numericHash.put(filename, new boolean[namesHash.get(filename).length]);
			for (int i = 0; i < numericHash.get(filename).length; i++) {
				numericHash.get(filename)[i] = true;
			}

			dataHash.put(filename, new Vector<String[]>());
			commentHash.put(filename, new Vector<String[]>());

			while (reader.ready()) {
				if (readBuffer.contains("\t")) {
					line = reader.readLine().trim().split("\t", -1);
				} else {
					line = reader.readLine().trim().split("[\\s]+");
				}
				ArrayList<String[]> sepRecords = getComments(line);
				commentHash.get(filename).add(sepRecords.get(1));
				dataHash.get(filename).add(sepRecords.get(0));
				for (int i = 0; i < sepRecords.get(0).length; i++) {
					if (!ext.isValidDouble(sepRecords.get(0)[i])) {
						numericHash.get(filename)[i] = false;
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.report("Error: file \"" + filename + "\" not found in current directory");
		} catch (IOException ioe) {
			log.report("Error reading file \"" + filename + "\"");
		}
	}

	/**
	 * Function to get comments from a given line
	 * 
	 * @param line
	 *            String[] which contains the line
	 * @return {@link ArrayList} of String[] in which data is at 0th position and comment is at 1st position
	 */
	private ArrayList<String[]> getComments(String[] line) {
		String[] comment = new String[line.length];
		if (line != null) {
			for (int i = 0; i < line.length; i++) {
				if (line[i].contains(":")) {
					// The data has comment separated by :
					// split and get comments and update data
					String[] temp = line[i].split(":");
					// the first record is data and rest all is comment
					// update line to have just data
					line[i] = temp[0];
					// get the comment; if the comment has : inside it then it will be
					// splitted into multiple parts so marge it
					String thisComment = null;
					for (int j = 1; j < temp.length; j++) {
						thisComment += temp[j];
					}
					// Store the complete comment
					comment[i] = thisComment;
				} else {
					// the data does not have comments
					comment[i] = "NA";
				}
			}
		}
		ArrayList<String[]> result = new ArrayList<String[]>();
		result.add(line);
		result.add(comment);
		return result;
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
				for (int i = 0; tree != null && i < tree.getModel().getChildCount(tree.getModel().getRoot()); i++) {
					if (ext.replaceAllWith(fileChooser.getSelectedFile().toString(), "\\", "/").equals(tree.getModel().getChild(tree.getModel().getRoot(), i).toString())) {
						found = true;
						break;
					}
				}
				if (found) {
					JOptionPane.showMessageDialog(null, "The data file has already been opened.");
				} else {
					loadFile(ext.replaceAllWith(fileChooser.getSelectedFile().toString(), "\\", "/"));
					updateTree();
					displayAll();
					tree.expandRow(0);
				}
			}
		} else if (command.equals(REMOVE_DATA_FILE)) {
			numberOfSelectedNodes = (byte) tree.getSelectedPathComponent();
			if (numberOfSelectedNodes != -1) {
				keys = HashVec.getKeys(dataHash); // keys is better to be block variable than a class variable. Otherwise, keys need to be updated every time there is an adding or deleting.
				tree.deleteSelectedNode();
				dataHash.remove(keys[numberOfSelectedNodes]);// TODO tree.getSelectionValues()[0][0] is not the branch to delete.
				namesHash.remove(keys[numberOfSelectedNodes]);
				keyIndices.remove(keys[numberOfSelectedNodes]);
				numericHash.remove(keys[numberOfSelectedNodes]);
				treeFilenameLookup.remove(keys[numberOfSelectedNodes]);
			}
		} else if (command.equals(SET_AS_COLORKEY)) {
		} else if (command.equals(SET_AS_LINKKEY)) {
			tree.getSelectionRows();
		} else {
			System.err.println("Error - unknown command '" + command + "'");
		}
	}

	public String[] getNamesSelected() {
		String[][] selectionValues = tree.getSelectionValues();
		String[] result = new String[selectionValues.length];
		for (int i = 0; i < selectionValues.length; i++) {
			if (selectionValues[i][0] == null) {
				result[i] = "";
			} else {
				result[i] = ext.removeDirectoryInfo(selectionValues[i][0]) + " _ " + namesHash.get(selectionValues[i][0])[Integer.parseInt(selectionValues[i][1])];
			}
		}
		return result;
	}

	public ArrayList<Vector<String[]>> getDataSelected(boolean includeColorKeyValue) {
		String[][] selectedNodes;
		Vector<String[]> dataOfSelectedFile;
		Hashtable<String, String[]> xHash;
		Hashtable<String, String> yHash;
		String[] inLine;
		int selectedColumn;
		String[] keys;
		Vector<String[]> v;
		byte colorCode;

		selectedNodes = tree.getSelectionValues();
		v = new Vector<String[]>();
		ArrayList<Vector<String[]>> result = new ArrayList<Vector<String[]>>();
		for (int a = 0; a < selectedNodes.length; a++) {
			if (selectedNodes[a][0] != null && keyIndices.get(selectedNodes[a][0]) != null) {
				xHash = new Hashtable<String, String[]>();
				selectedColumn = Integer.parseInt(selectedNodes[a][1]);
				dataOfSelectedFile = dataHash.get(selectedNodes[a][0]);
				yHash = new Hashtable<String, String>();
				for (int i = 0; i < dataOfSelectedFile.size(); i++) {
					inLine = dataOfSelectedFile.elementAt(i);
					yHash.put(String.valueOf(i), inLine[selectedColumn]);
					xHash.put(String.valueOf(i), new String[] { String.valueOf(i) });
				}
				String recordName = namesHash.get(selectedNodes[a][0])[selectedColumn];
				keys = HashVec.getKeys(xHash, false, false);
				v = new Vector<String[]>();
				if (includeColorKeyValue) {
					for (int i = 0; i < keys.length; i++) {
						if (yHash.containsKey(keys[i])) {
							colorCode = getColorFromFilename(recordName);
							v.add(new String[] { xHash.get(keys[i])[0], xHash.get(keys[i])[0], yHash.get(keys[i]), colorCode + "", xHash.get(keys[i])[0] });
						}
					}
				} else {
					for (int i = 0; i < keys.length; i++) {
						if (yHash.containsKey(keys[i])) {
							v.add(new String[] { xHash.get(keys[i])[0], xHash.get(keys[i])[0], yHash.get(keys[i]), "all other keys", xHash.get(keys[i])[0] });
						}
					}
				}
			}
			result.add(v);
		}
		return result;
	}

	private void initializeTree() {
		tree = new CheckBoxTree(new String[0], new String[0], new String[0][], new boolean[0], 2);
	}

	private void updateTree() {
		String[] namesOfBranches;
		String[] branchHandles;
		String[][] namesOfNodes;

		treeFileVariableNameLookup = new String[treeFilenameLookup.size()][];
		int maxSelectable = 0;
		for (int i = 0; i < treeFilenameLookup.size(); i++) {
			treeFileVariableNameLookup[i] = namesHash.get(treeFilenameLookup.elementAt(i));
			maxSelectable += treeFileVariableNameLookup[i].length;
			tree.updateSelections(maxSelectable);

			if (tree == null) {
				namesOfBranches = new String[1];
				branchHandles = new String[1];
				namesOfNodes = new String[1][];
				namesOfBranches[0] = ext.removeDirectoryInfo(treeFilenameLookup.elementAt(i));
				branchHandles[0] = treeFilenameLookup.elementAt(i);
				namesOfNodes[0] = treeFileVariableNameLookup[i];
				tree = new CheckBoxTree(namesOfBranches, branchHandles, namesOfNodes, numericHash.get(treeFilenameLookup.elementAt(0)), 2);
			} else {
				boolean found = false;
				for (int j = 0; j < tree.getModel().getChildCount(tree.getModel().getRoot()); j++) {
					if (tree.getModel().getChild(tree.getModel().getRoot(), j).toString().equals(ext.removeDirectoryInfo(treeFilenameLookup.elementAt(i)))) {
						found = true;
					}
				}
				if (!found) {
					tree.addNode(ext.removeDirectoryInfo(treeFilenameLookup.elementAt(i)), treeFilenameLookup.elementAt(i), treeFileVariableNameLookup[i], numericHash.get(treeFilenameLookup.elementAt(i)), null);
				}
			}
		}
	}

	/**
	 * Create the GUI and show it. For thread safety, this method should be invoked from the event-dispatching thread.
	 * 
	 * @param proj
	 *            {@link Project} the project
	 * @param log
	 *            {@link Logger} a logger
	 */
	public static void createAndShowGUI(Project proj, Logger log) {

		// Create and set up the window.
		JFrame frame = new JFrame("Line Plot");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		// Create and set up the content pane.
		log.report("Creating new Line Plot object");
		LinePlot twoDPlot = new LinePlot(proj, log);
		frame.setJMenuBar(twoDPlot.menuBar());
		twoDPlot.setOpaque(true); // content panes must be opaque
		frame.setContentPane(twoDPlot);
		frame.addWindowListener(twoDPlot);
		frame.setBounds(20, 20, 1000, 600);

		// Display the window.
		frame.pack();
		frame.setVisible(true);
	}

	public void updateGUI() {
		linePanel.paintAgain();
	}

	public int[] getCurrentLinkKeyColumnLabels() {
		String[][] selectedValues;

		selectedValues = tree.getSelectionValues();
		if (selectedValues == null || selectedValues[0][0] == null) {
			return null;
		}
		return keyIndices.get(selectedValues[0][0]);
	}

	public Project getProject() {
		return proj;
	}

	public byte getPointSize() {
		return DEFAULT_SIZE;
	}

	public static void main(String[] args) {

		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				createAndShowGUI(new Project(Project.DEFAULT_PROJECT, false), new Logger());
			}
		});

	}
}