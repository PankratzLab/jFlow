package org.genvisis.cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
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
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.ImageIcon;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLayeredPane;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.CheckBoxTree;
import org.genvisis.cnv.gui.ColorKeyPanel;
import org.genvisis.cnv.gui.GuiManager;
import org.genvisis.cnv.gui.UITools;
import org.genvisis.cnv.plots.data.DataPipe;
import org.genvisis.cnv.var.IndiPheno;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Grafik;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Numbers;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;

import com.google.common.base.Strings;

public class TwoDPlot extends JPanel implements WindowListener, ActionListener, TreeSelectionListener {
	public static final long serialVersionUID = 1L;
	public static final byte DEFAULT_SIZE = 8;
	public static final int DEFAULT_GC_THRESHOLD = 25;
	public static final String TWO_D_STICH_COMMAND = "twoDScreenshots";
	public static final String TWO_D_STICH_COMMAND_DESCRIPTION = "Have two-d plot create a series of images";
	private static final String ALT_UP = "ALT UP";
	private static final String ALT_DOWN = "ALT DOWN";
	private static final String ALT_LEFT = "ALT LEFT";
	private static final String ALT_RIGHT = "ALT RIGHT";
	// public static final String MASK_MISSING = "Mask missing values";
	// public static final String UNMASK_MISSING = "Unmask missing values";
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String ADD_DATA_FILE = "Add Data File";
	public static final String REMOVE_DATA_FILE = "Remove Data";
	public static final String REMOVE_ALL = "Clear Data";
	public static final String CREATE_SCREENS = "Create Screenshots";
	// public static final String SET_AS_COLORKEY = "Set as Color Key";
	// public static final String SET_AS_LINKKEY = "Set as Link Key";

	public static final String[] COMMENT_FIELD_ALIASES = {
																												"COMMENT",
																												"NOTE"
	};
	public static final String[][] LINKERS = {Aliases.INDIVIDUAL_ID, Aliases.FAMILY_ID, Aliases.DNA,
																						Aliases.MARKER_NAMES, Aliases.REGION, Aliases.CHRS,
																						ArrayUtils.combine(Aliases.POSITIONS,
																															 Aliases.POSITIONS_START),
																						Aliases.POSITIONS_STOP,
																						COMMENT_FIELD_ALIASES};
	public static final int IID_INDEX_IN_LINKERS = 0;
	public static final int FID_INDEX_IN_LINKERS = 1;
	public static final int DNA_INDEX_IN_LINKERS = 2;
	public static final int MARKER_INDEX_IN_LINKERS = 3;
	public static final int REGION_INDEX_IN_LINKERS = 4;
	public static final int CHR_INDEX_IN_LINKERS = 5;
	public static final int POS_INDEX_IN_LINKERS = 6;
	public static final int STOP_POS_INDEX_IN_LINKERS = 7;
	public static final int COMMENT_INDEX_IN_LINKERS = 8;

	public static final String[] MISSING_VALUES = {".", "NA"};

	/*
	 * regex to match and pull out column title, chromosome, and position
	 */
	private static String DATA_TITLE_REGEX = "(.+)_chr([0-2[XYM]]?[0-9[Y]]*):([\\d,]+-[\\d,]+).*";

	// private JPanel bottomPanel;
	private TwoDPanel twoDPanel;
	private JLayeredPane layeredPane;
	protected ColorKeyPanel colorKeyPanel;

	private Project proj;
	private byte size;
	private SampleData sampleData;
	private JButton flipButton, invXButton, invYButton;
	private volatile boolean flipStatus, xInvStatus, yInvStatus, hideExcludes;
	private volatile boolean isHistPlot;
	private CheckBoxTree tree;
	private ArrayList<String> dataKeys;
	// Map file short names to absolute path
	private Map<String, String> filenameMap;
	HashMap<String, ArrayList<String[]>> dataHash;
	HashMap<String, String[]> dataColumnsHash;
	HashMap<String, boolean[]> validColumnsHash;
	HashMap<String, HashMap<String, DataPipe>> fileColumnDataPipes;
	HashSet<String> validExts;

	class ColumnMetaData {
		String lbl;
		String chr;
		String region;
		String start;
		String stop;

		public ColumnMetaData(String l, String c, String r, String s1, String s2) {
			this.lbl = l;
			this.chr = c;
			this.region = r;
			this.start = s1;
			this.stop = s2;
		}
	}

	HashMap<String, HashMap<Integer, ColumnMetaData>> columnMetaData;

	HashMap<String, int[]> linkerIndices;
	Logger log;
	Pattern headerPattern = Pattern.compile(DATA_TITLE_REGEX);
	private final boolean promptOnClose;

	private volatile boolean generatingScreenshots;
	private HashMap<String, Integer> colorData;
	private JPanel treePanel;
	private JScrollPane scrollPane;
	private final HashSet<String> selectedDataHash = new HashSet<String>();
	private final ImageIcon flip1_1 = Grafik.getImageIcon("images/flip_and_invert/flip_10p.jpg");
	private final ImageIcon flip1_2 = Grafik.getImageIcon("images/flip_and_invert/flip_10p_blue_3.jpg");
	private final ImageIcon flip2_1 = Grafik.getImageIcon("images/flip_and_invert/flip_10p_inv.jpg");
	private final ImageIcon flip2_2 = Grafik.getImageIcon("images/flip_and_invert/flip_10p_inv_blue_3.jpg");
	private final ImageIcon flipX1_1 = Grafik.getImageIcon("images/flip_and_invert/right_10.gif");
	private final ImageIcon flipX1_2 = Grafik.getImageIcon("images/flip_and_invert/right_10_blue.jpg");
	private final ImageIcon flipX2_1 = Grafik.getImageIcon("images/flip_and_invert/left_10.gif");
	private final ImageIcon flipX2_2 = Grafik.getImageIcon("images/flip_and_invert/left_10_blue.jpg");
	private final ImageIcon flipY1_1 = Grafik.getImageIcon("images/flip_and_invert/up_10.gif");
	private final ImageIcon flipY1_2 = Grafik.getImageIcon("images/flip_and_invert/up_10_blue.jpg");
	private final ImageIcon flipY2_1 = Grafik.getImageIcon("images/flip_and_invert/down_10.gif");
	private final ImageIcon flipY2_2 = Grafik.getImageIcon("images/flip_and_invert/down_10_blue.jpg");

	private TwoDPlot() {
		this(new Project(), false, null);
	}

	private TwoDPlot(Project proj, boolean promptOnClose, String[] fileExts) {
		this.promptOnClose = promptOnClose;
		String[] previouslyLoadedFiles;

		this.proj = proj;
		log = proj == null ? new Logger() : proj.getLog();
		size = DEFAULT_SIZE;

		if (proj != null && Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false))) {
			sampleData = proj.getSampleData(false);
		} else {
			sampleData = null;
		}
		validExts = new HashSet<String>();
		dataKeys = new ArrayList<String>();
		dataHash = new HashMap<String, ArrayList<String[]>>();
		filenameMap = new HashMap<>();
		dataColumnsHash = new HashMap<String, String[]>();
		validColumnsHash = new HashMap<String, boolean[]>();
		fileColumnDataPipes = new HashMap<>();
		linkerIndices = new HashMap<String, int[]>();
		columnMetaData = new HashMap<String, HashMap<Integer, ColumnMetaData>>();

		if (fileExts != null) {
			for (String ext : fileExts) {
				validExts.add(ext);
			}
		}

		setLayout(new BorderLayout());

		twoDPanel = new TwoDPanel(this);

		// UIManager.put("PopupMenuUI", "CustomPopupMenuUI");

		layeredPane = new JLayeredPane() {
			private static final long serialVersionUID = 1L;

			@Override
			public boolean isOptimizedDrawingEnabled() {
				return false; // override to force proper z-ordered painting
			}
		};
		layeredPane.setLayout(new BorderLayout());

		generateFlipButton();
		layeredPane.add(flipButton);
		generateInvXButton();
		layeredPane.add(invXButton);
		generateInvYButton();
		layeredPane.add(invYButton);

		layeredPane.setComponentZOrder(flipButton, 0);
		layeredPane.setComponentZOrder(invXButton, 0);
		layeredPane.setComponentZOrder(invYButton, 0);

		layeredPane.add(twoDPanel);
		layeredPane.setComponentZOrder(twoDPanel, 3);
		layeredPane.setPreferredSize(new Dimension(1000, 600));

		treePanel = new JPanel();
		treePanel.setBackground(BACKGROUND_COLOR);
		treePanel.setLayout(new BorderLayout());

		tree = new CheckBoxTree(new String[0], new String[0], new String[0][], new boolean[0], 2);
		updateTreeForNewData();

		scrollPane = new JScrollPane(tree);
		treePanel.add(scrollPane, BorderLayout.CENTER);
		tree.addTreeSelectionListener(this);
		JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, treePanel, layeredPane);
		splitPane.setBackground(Color.WHITE);
		splitPane.setOneTouchExpandable(true);
		splitPane.setDividerLocation(150);

		Dimension minimumSize = new Dimension(100, 50);
		layeredPane.setMinimumSize(minimumSize);

		add(splitPane, BorderLayout.CENTER);

		colorKeyPanel = new ColorKeyPanel(sampleData, twoDPanel, TwoDPanel.DEFAULT_COLORS);
		layeredPane.add(colorKeyPanel, BorderLayout.SOUTH);
		layeredPane.setComponentZOrder(colorKeyPanel, 3);

		inputMapAndActionMap();

		twoDPanel.setPointsGeneratable(true);// zx
		twoDPanel.setExtraLayersVisible(new byte[] {99});
		updateGUI();

		twoDPanel.grabFocus();

		layeredPane.addComponentListener(new ComponentListener() {
			@Override
			public void componentShown(ComponentEvent e) {}

			@Override
			public void componentResized(ComponentEvent e) {
				flipButton.setBounds(70, layeredPane.getHeight() - 75, 38, 38);
				invXButton.setBounds(70, layeredPane.getHeight() - 35, 38, 13);
				invYButton.setBounds(55, layeredPane.getHeight() - 75, 13, 38);
			}

			@Override
			public void componentMoved(ComponentEvent e) {}

			@Override
			public void componentHidden(ComponentEvent e) {}
		});

		if (proj != null && promptOnClose) {
			previouslyLoadedFiles = proj.TWOD_LOADED_FILENAMES.getValue();
			String errMsg = "";
			Set<String> failed = new HashSet<String>();
			for (String previouslyLoadedFile : previouslyLoadedFiles) {
				String result = loadFile(previouslyLoadedFile);
				if (!result.isEmpty()) {
					failed.add(previouslyLoadedFile);
					errMsg += result + "\n";
				}
			}

			if (!failed.isEmpty()) {
				proj.message("The following file(s) failed to load and will no longer be displayed by default:\n"
										 + errMsg);

				String[] loadedSuccessfully = new String[previouslyLoadedFiles.length - failed.size()];
				int i = 0;
				for (String s : previouslyLoadedFiles) {
					if (!failed.contains(s)) {
						loadedSuccessfully[i++] = s;
					}
				}
				proj.TWOD_LOADED_FILENAMES.setValue(loadedSuccessfully);
			}

			updateTreeForNewData();
			String[] sel = proj.TWOD_LOADED_VARIABLES.getValue();
			List<String> passed = new ArrayList<String>();
			// Determine which variable(s) of this file were selected
			for (String s : sel) {
				// Selection syntax is "filename|variable_index"
				String[] fileSel = s.split("\\|");

				if (fileSel.length != 2) {
					// No variables were selected
					continue;
				}

				if (Files.isRelativePath(fileSel[0])) {
					fileSel[0] = proj.PROJECT_DIRECTORY.getValue() + fileSel[0];
				}

				// Skip selections from files that were removed
				if (!failed.contains(fileSel[0])) {
					passed.add(s);
					tree.performCheckBoxAction(dataColumnsHash.get(fileSel[0])[Integer.parseInt(fileSel[1])],
																		 ItemEvent.SELECTED);
				}
			}

			// Update the project property if necessary
			if (passed.size() != sel.length) {
				proj.TWOD_LOADED_VARIABLES.setValue(passed.toArray(new String[passed.size()]));
			}
		}

		setVisible(true);
		generateShortcutMenus();
	}

	public void refreshOtherButtons() {
		flipButton.repaint();
		invXButton.repaint();
		invYButton.repaint();
	}

	public TwoDPanel getPanel() {
		return twoDPanel;
	}

	// Mouse listener for the checkbox tree
	MouseListener checkBoxMouseListener = new MouseListener() {

		@Override
		public void mouseClicked(MouseEvent e) {
			String[][] selectionValues = tree.getSelectionValues();
			selectedDataHash.clear();
			for (String[] selectionValue : selectionValues) {
				selectedDataHash.add(selectionValue[0] + "***" + selectionValue[1]);
			}
		}

		@Override
		public void mousePressed(MouseEvent e) {}

		@Override
		public void mouseReleased(MouseEvent e) {

		}

		@Override
		public void mouseEntered(MouseEvent e) {}

		@Override
		public void mouseExited(MouseEvent e) {}
	};

	private void generateFlipButton() {
		flipButton = new JButton(flip1_1);
		flipButton.setRolloverIcon(flip1_2);
		flipButton.setToolTipText("Inverts X & Y Axes");
		flipButton.setBorder(null);
		flipButton.setVisible(true);
		flipStatus = true;
		flipButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				twoDPanel.setPointsGeneratable(true);
				twoDPanel.setSwapAxes(flipStatus);
				twoDPanel.paintAgain();
				if (flipStatus) {
					flipStatus = false;
					flipButton.setIcon(flip2_1);
					flipButton.setRolloverIcon(flip2_2);
				} else {
					flipStatus = true;
					flipButton.setIcon(flip1_1);
					flipButton.setRolloverIcon(flip1_2);
				}
			}
		});
		flipButton.setDoubleBuffered(true);
		flipButton.setContentAreaFilled(false);
	}

	private void generateInvXButton() {
		invXButton = new JButton(flipX1_1);
		invXButton.setRolloverIcon(flipX1_2);
		invXButton.setBorder(null);
		invXButton.setVisible(true);
		invXButton.setToolTipText("Inverts X Axis");
		xInvStatus = true;
		invXButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				twoDPanel.setPointsGeneratable(true);
				twoDPanel.setXinversion(xInvStatus);
				twoDPanel.paintAgain();
				if (xInvStatus) {
					invXButton.setIcon(flipX2_1);
					invXButton.setRolloverIcon(flipX2_2);
				} else {
					invXButton.setIcon(flipX1_1);
					invXButton.setRolloverIcon(flipX1_2);
				}
				xInvStatus = !xInvStatus;
			}
		});
		invXButton.setDoubleBuffered(true);
		invXButton.setContentAreaFilled(false);
	}

	private void generateInvYButton() {
		invYButton = new JButton(flipY1_1);
		invYButton.setRolloverIcon(flipY1_2);
		invYButton.setBorder(null);
		invYButton.setVisible(true);
		invYButton.setToolTipText("Inverts Y Axis");
		yInvStatus = true;
		invYButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				twoDPanel.setPointsGeneratable(true);
				twoDPanel.setYinversion(yInvStatus);
				twoDPanel.paintAgain();
				if (yInvStatus) {
					invYButton.setIcon(flipY2_1);
					invYButton.setRolloverIcon(flipY2_2);
				} else {
					invYButton.setIcon(flipY1_1);
					invYButton.setRolloverIcon(flipY1_2);
				}
				yInvStatus = !yInvStatus;
			}
		});
		invYButton.setDoubleBuffered(true);
		invYButton.setContentAreaFilled(false);
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
		JMenuItem menuItemExit, menuItemOpen, menuItemRemove, menuItemRemoveAll, menuItemScreens, menuItemSave;
		final JCheckBoxMenuItem menuItemExclude, menuItemHist;

		menuBar = new JMenuBar();
		menu = new JMenu("File");
		menu.setMnemonic(KeyEvent.VK_F);
		menuBar.add(menu);
		menuItemOpen = new JMenuItem("Open File", KeyEvent.VK_O);
		menuItemOpen.addActionListener(e -> {
			addFile();
		});
		menu.add(menuItemOpen);
		menuItemRemove = new JMenuItem("Remove Selected Data", KeyEvent.VK_O);
		menuItemRemove.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				removeSelectedData();
			}
		});
		menu.add(menuItemRemove);
		menuItemRemoveAll = new JMenuItem("Remove All Data", KeyEvent.VK_O);
		menuItemRemoveAll.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				removeAllData();
			}
		});
		menu.add(menuItemRemoveAll);
		menuItemSave = new JMenuItem("Save Current Image", KeyEvent.VK_S);
		menuItemSave.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				JFileChooser fileChooser = new JFileChooser(proj != null
																																? proj.PROJECT_DIRECTORY.getValue()
																																: ".");
				int fileOpenActionSelected = fileChooser.showSaveDialog(TwoDPlot.this);
				if (fileOpenActionSelected == JFileChooser.APPROVE_OPTION) {
					File fileToOpen = fileChooser.getSelectedFile();
					twoDPanel.screenCapture(fileToOpen.toString() + ".png");
				}
			}
		});
		menu.add(menuItemSave);
		menuItemScreens = new JMenuItem("Create Screenshots from File", KeyEvent.VK_S);
		menuItemScreens.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				createScreenshotsFromFile();
			}
		});
		menu.add(menuItemScreens);
		menuItemExit = new JMenuItem("Close", KeyEvent.VK_C);
		menuItemExit.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				SwingUtilities.invokeLater(new Runnable() {
					@Override
					public void run() {
						setVisible(false);
						Component parent = getParentComponent();
						parent.dispatchEvent(new WindowEvent((JFrame) parent, WindowEvent.WINDOW_CLOSING));
					}
				});
			}
		});
		menu.add(menuItemExit);
		menu = new JMenu("View");
		menu.setMnemonic(KeyEvent.VK_V);
		menuBar.add(menu);

		hideExcludes = true;
		menuItemExclude = new JCheckBoxMenuItem("Hide Excluded", true);
		menuItemExclude.setMnemonic(KeyEvent.VK_E);
		menuItemExclude.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				hideExcludes = menuItemExclude.isSelected();
				twoDPanel.paintAgain();
			}
		});
		menu.add(menuItemExclude);

		menuItemHist = new JCheckBoxMenuItem("Display as Histogram", false);
		menuItemHist.setMnemonic(KeyEvent.VK_H);
		menuItemHist.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				setHistogram(menuItemHist.isSelected());

				tree.setMaxSelections(isHistogram() ? 1 : 2);

				twoDPanel.paintAgain();
			}
		});
		menu.add(menuItemHist);

		return menuBar;
	}


	private Component getParentComponent() {
		Component c = TwoDPlot.this.getParent();
		while (!(c instanceof JFrame)) {
			c = c.getParent();
		}
		return c;
	}

	@Override
	public void actionPerformed(ActionEvent ae) {
		String command;

		command = ae.getActionCommand();
		if (command.equals(ADD_DATA_FILE)) {
			addFile();
		} else if (command.equals(REMOVE_DATA_FILE)) {
			removeSelectedData();
		} else if (command.equals(REMOVE_ALL)) {
			removeAllData();
		} else if (command.equals(CREATE_SCREENS)) {
			createScreenshotsFromFile();
		} else {
			System.err.println("Error - unknown command '" + command + "'");
		}
	}

	private void createScreenshotsFromFile() {
		JFileChooser jfc = new JFileChooser(proj != null ? proj.PROJECT_DIRECTORY.getValue()
																										: (new File(".")).toString());
		jfc.setMultiSelectionEnabled(true);
		if (jfc.showOpenDialog(TwoDPlot.this) == JFileChooser.APPROVE_OPTION) {
			File selFile = jfc.getSelectedFile();
			List<String> params = parseControlFile(selFile.getAbsolutePath(), log);
			// final String projFile = params.get(0).split("=")[1];
			final String baseDir = params.get(1).split("=")[1];
			final ArrayList<ScreenToCapture> screens = condenseCtrlFile(params.subList(2, params.size()),
																																	true);
			createScreenshots(baseDir, screens);
		}
	}

	private void addFile() {
		JFileChooser fileChooser;
		boolean found = false;
		fileChooser = new JFileChooser(proj != null ? proj.PROJECT_DIRECTORY.getValue() : ".");
		if (fileChooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
			for (int i = 0; tree != null
											&& i < tree.getModel().getChildCount(tree.getModel().getRoot()); i++) {
				if (ext.replaceAllWith(fileChooser.getSelectedFile().toString(), "\\", "/")
							 .equals(tree.getModel().getChild(tree.getModel().getRoot(), i).toString())) {
					found = true;
					break;
				}
			}
			if (found) {
				JOptionPane.showMessageDialog(null, "The data file has already been opened.");
			} else {
				loadFile(ext.replaceAllWith(fileChooser.getSelectedFile().toString(), "\\", "/"));
				updateTreeForNewData();
			}
		}
	}

	public void removeAllData() {
		dataKeys.clear();
		dataHash.clear();
		dataColumnsHash.clear();
		validColumnsHash.clear();
		linkerIndices.clear();
		columnMetaData.clear();
		fileColumnDataPipes.clear();
		tree.deleteAllNodes();
		updateGUI();
	}

	public void removeSelectedData() {
		String toRemove = tree.getSelectedPathComponentName();
		if (toRemove != null) {
			// When removing items from the GUI we may just have the file name and not the absolute path
			// Ensure both are purged.
			String path = filenameMap.get(toRemove);
			tree.deleteSelectedNode();
			for (String s : new String[] {toRemove, path}) {
				dataHash.remove(s);
				dataColumnsHash.remove(s);
				linkerIndices.remove(s);
				validColumnsHash.remove(s);
				fileColumnDataPipes.remove(s);
				dataKeys.remove(s);
			}
		}
	}


	public void reloadSampleDataUI() {
		if (proj != null) {
			// Note: This line resets the sample data by setting it to null which will cause linkKeyIndex
			// to
			// be lost.
			// Make sure to set it back after resetting sample data
			// Hashtable<String, Integer> linkKeyIndexCopy = new Hashtable<String,
			// Integer>(sampleData.getLinkKeyIndex());
			proj.resetSampleData();
			if (Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false))) {
				sampleData = proj.getSampleData(false);
			}
			sampleData = proj.getSampleData(false);
			// sampleData.setLinkKeyIndex(linkKeyIndexCopy);
			colorKeyPanel.updateSampleData(sampleData);
		}
		colorKeyPanel.updateColorVariablePanel();

		twoDPanel.paintAgain();
		generateShortcutMenus();
	}

	// public Hashtable<String, String> createHashWithSampleID(Hashtable<String, String>
	// colorKeyValue) {
	// Hashtable<String, String> colorKeyValueHash;
	//
	// colorKeyValueHash = new Hashtable<String, String>();
	// for (String key : colorKeyValue.keySet()) {
	// colorKeyValueHash.put(twoDPanel.getSampleData().lookup(key)[0], colorKeyValue.get(key));
	// }
	//
	// return colorKeyValueHash;
	// }

	// public void setLinkKeyHandler(int[] selectedLinkKey) {
	// String selectedCol;
	//
	// selectedCol = tree.getSelectionValues()[0][0];
	//
	// if (selectedLinkKey.length == 1) {
	// if (!linkerIndices.containsKey(selectedCol)) {
	// JOptionPane.showMessageDialog(null, "There was a problem in your selection. Please select
	// again", "Error", JOptionPane.ERROR_MESSAGE);
	// return;
	// }
	// sampleData.setLinkKey(selectedCol, selectedLinkKey[0]);
	// }
	// }

	// public void createLinkKeyToDataHash(String filename, int[] linkKeyColumnLabels) {
	// String inLine;
	// String[] inLineArray;
	// Hashtable<String, String[]> curFileLinkKeyDataHash = new Hashtable<String, String[]>();
	// try {
	// BufferedReader fileReader = new BufferedReader(new FileReader(filename));
	// inLine = fileReader.readLine(); // headers: skip this line
	//
	// System.out.println("link column labels:" + Arrays.toString(linkKeyColumnLabels));
	// System.out.println("link key index: " + linkKeyIndex.get(filename));
	// while (fileReader.ready()) {
	// if (inLine.contains("\t")) {
	// inLineArray = fileReader.readLine().trim().split("\t", -1);
	// } else {
	// inLineArray = fileReader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
	// }
	// if (linkKeyIndex.containsKey(filename)) {
	// switch (linkKeyIndex.get(filename)) {
	// case DNA_INDEX_IN_LINKERS:
	// curFileLinkKeyDataHash.put(inLineArray[linkKeyColumnLabels[linkKeyIndex.get(filename)]],
	// inLineArray);
	// break;
	// case FID_INDEX_IN_LINKERS:
	// curFileLinkKeyDataHash.put(inLineArray[linkKeyColumnLabels[linkKeyIndex.get(filename)]] + "\t"
	// + inLineArray[linkKeyColumnLabels[IID_INDEX_IN_LINKERS]], inLineArray);
	// break;
	// case IID_INDEX_IN_LINKERS:
	// curFileLinkKeyDataHash.put(inLineArray[linkKeyColumnLabels[linkKeyIndex.get(filename)]],
	// inLineArray);
	// break;
	// default:
	// System.out.println("Error: Invalid link key.");
	// break;
	// }
	// }
	// }
	// fileReader.close();
	// } catch (FileNotFoundException fnfe) {
	// System.err.println("Error: file \"" + filename + "\" not found in current directory");
	// } catch (IOException ioe) {
	// System.err.println("Error reading file \"" + filename + "\"");
	// }
	// if (!curFileLinkKeyDataHash.isEmpty())
	// linkKeyToDataHash.put(filename, curFileLinkKeyDataHash);
	// System.out.println("LinkKeyDataHash:" + linkKeyToDataHash.toString());
	// }

	public byte getPointSize() {
		return size;
	}

	public void setPointSize(byte sz) {
		size = sz;
	}

	public String[] getNamesSelected() {
		ArrayList<String> resultList = new ArrayList<String>();
		String[][] metaData = getCurrentColumnMetaData();
		String[][] selectedValues = tree.getSelectionValues();
		for (int i = 0; i < metaData.length; i++) {
			resultList.add(metaData[i] == null
																				? dataColumnsHash.get(selectedValues[i][0])[Integer.parseInt(selectedValues[i][1])]
																				: metaData[i][4] + " chr" + metaData[i][0] + ":"
																					+ metaData[i][1]);
		}
		String[] result = resultList.toArray(new String[resultList.size()]);
		return result;
	}

	public int[] getCurrentLinkKeyColumnIndices() {
		String[][] selectedValues;

		selectedValues = tree.getSelectionValues();
		if (selectedValues == null || selectedValues[0][0] == null) {
			return null;
		}
		return linkerIndices.get(selectedValues[0][0]);
	}

	public Project getProject() {
		return proj;
	}

	public ArrayList<String[]> getDataSelected(boolean includeColorKeyValue) {
		String[][] selectedNodes;
		ArrayList<String[]> dataOfSelectedFile;
		HashMap<String, String[]> xHash;
		String[] inLine, outLine;
		int selectedColumn;
		String selectedFile;
		Set<String> keys;
		ArrayList<String[]> v;
		int currentClass;
		String[] ids;
		byte colorCode;
		int[] linkKeyColumnIndices;
		byte index;

		DataPipe pipeX;
		DataPipe pipeY;

		// file/collection name and column number of selected data collection
		selectedNodes = tree.getSelectionValues();
		v = new ArrayList<String[]>();
		xHash = new HashMap<String, String[]>();
		if (selectedNodes[0][0] != null && linkerIndices.get(selectedNodes[0][0]) != null) {
			selectedColumn = Integer.parseInt(selectedNodes[0][1]);
			selectedFile = selectedNodes[0][0];
			pipeX = fileColumnDataPipes.get(selectedFile)
																 .get(dataColumnsHash.get(selectedFile)[selectedColumn]);
			dataOfSelectedFile = dataHash.get(selectedFile);
			currentClass = colorKeyPanel.getCurrentClass();
			if (currentClass < SampleData.BASIC_CLASSES.length
					&& SampleData.BASIC_CLASSES[currentClass].equals(SampleData.HEATMAP)) {
				twoDPanel.setChartType(AbstractPanel.HEAT_MAP_TYPE);
			} else {
				twoDPanel.setChartType(AbstractPanel.SCATTER_PLOT_TYPE);
			}

			linkKeyColumnIndices = linkerIndices.get(selectedFile);
			index = (byte) (includeColorKeyValue ? 4 : 3);
			for (int i = 0; i < dataOfSelectedFile.size(); i++) {
				inLine = dataOfSelectedFile.get(i);

				if (sampleData != null
						&& (hideExcludes || colorKeyPanel.getDisabledClassValues().size() > 0)) {

					// Try DNA first
					IndiPheno indi = null;
					String curSample = null;
					if (linkKeyColumnIndices[DNA_INDEX_IN_LINKERS] >= 0) {
						curSample = inLine[linkKeyColumnIndices[DNA_INDEX_IN_LINKERS]];
						indi = sampleData.getIndiFromSampleHash(curSample);
					}
					if (indi == null) {
						// Try FID/IID
						if (linkKeyColumnIndices[FID_INDEX_IN_LINKERS] >= 0
								&& linkKeyColumnIndices[IID_INDEX_IN_LINKERS] >= 0) {
							curSample = inLine[linkKeyColumnIndices[FID_INDEX_IN_LINKERS]] + "\t"
													+ inLine[linkKeyColumnIndices[IID_INDEX_IN_LINKERS]];
							curSample = sampleData.lookup(curSample)[0];
							indi = sampleData.getIndiFromSampleHash(curSample);
						}
					}

					// Check if hidden
					if (curSample != null && hideExcludes
							&& sampleData.individualShouldBeExcluded(curSample)) {
						continue;
					}

					if (indi != null) {
						int sampColorKey = sampleData.determineCodeFromClass(currentClass, (byte) 0, indi,
																																 (byte) 0, 0);

						if (colorKeyPanel.getDisabledClassValues()
														 .containsKey(currentClass + "\t" + sampColorKey)) {
							continue;
						}
					} else {
						// TODO assign color based on value, check if values are hidden
					}
				}
				outLine = new String[linkKeyColumnIndices.length + index];
				outLine[0] = pipeX == null ? inLine[0] : pipeX.pipe(inLine[0]);
				if (outLine[0] == null) { // rejected value from pipe
					continue;
				}
				outLine[1] = inLine[selectedColumn];
				for (int j = 0; j < linkKeyColumnIndices.length; j++) {
					if (linkKeyColumnIndices[j] >= 0) {
						outLine[j + index] = inLine[linkKeyColumnIndices[j]];
					}
				}
				String key = linkerIndices.get(selectedFile)[DNA_INDEX_IN_LINKERS] == -1
																																								? i + ""
																																								: inLine[linkerIndices.get(selectedFile)[DNA_INDEX_IN_LINKERS]];
				xHash.put(key, outLine);
			}
			if (selectedNodes.length > 1 && selectedNodes[1] != null && selectedNodes[1][0] != null
					&& linkerIndices.get(selectedNodes[1][0]) != null) {
				selectedColumn = Integer.parseInt(selectedNodes[1][1]);
				selectedFile = selectedNodes[1][0];
				pipeY = fileColumnDataPipes.get(selectedFile)
																	 .get(dataColumnsHash.get(selectedFile)[selectedColumn]);
				dataOfSelectedFile = dataHash.get(selectedFile);
				for (int i = 0; i < dataOfSelectedFile.size(); i++) {
					inLine = dataOfSelectedFile.get(i);
					String key = linkerIndices.get(selectedFile)[DNA_INDEX_IN_LINKERS] == -1
																																									? i + ""
																																									: inLine[linkerIndices.get(selectedFile)[DNA_INDEX_IN_LINKERS]];
					if (xHash.containsKey(key)) {
						String val = pipeY == null ? inLine[selectedColumn]
																			: pipeY.pipe(inLine[selectedColumn]);
						if (val == null) { // rejected by pipe
							xHash.remove(key);
							continue;
						}
						xHash.get(key)[2] = val;
						if (includeColorKeyValue) {
							if (sampleData != null) {
								ids = sampleData.lookup(key);
								if (ids == null) {
									colorCode = 0;
								} else if (generatingScreenshots) {
									colorCode = getColorForScreenshot(ids[0]);
								} else {
									String[][] metaData = getCurrentColumnMetaData();
									int chr = 0, start = 0, stop = 0;
									if (metaData != null && metaData.length != 0 && metaData[0] != null) {
										chr = Positions.chromosomeNumber(metaData[0][0]);
										start = Numbers.parseWithLocale(metaData[0][2]);
										stop = Numbers.parseWithLocale(metaData[0][3]);
									}
									colorCode = sampleData.determineCodeFromClass(currentClass,
																																(byte) 0,
																																sampleData.getIndiFromSampleHash(ids[0]),
																																(byte) chr,
																																start + ((stop - start) / 2));
								}
							} else {
								colorCode = 0;
								// TODO color by value
							}
							xHash.get(key)[3] = colorCode + "";
						}
					}
				}
			}
			keys = xHash.keySet();
			v = new ArrayList<String[]>();
			for (String key : keys) {
				inLine = xHash.get(key);
				if (includeColorKeyValue) {
					if (sampleData != null) {
						ids = sampleData.lookup(key);
						if (ids == null) {
							colorCode = 0;
						} else if (generatingScreenshots) {
							colorCode = getColorForScreenshot(ids[0]);
						} else {
							String[][] metaData = getCurrentColumnMetaData();
							int chr = 0, start = 0, stop = 0;
							if (metaData != null && metaData.length != 0 && metaData[0] != null) {
								chr = Positions.chromosomeNumber(metaData[0][0]);
								start = Numbers.parseWithLocale(metaData[0][2]);
								stop = Numbers.parseWithLocale(metaData[0][3]);
							}
							colorCode = sampleData.determineCodeFromClass(currentClass,
																														(byte) 0,
																														sampleData.getIndiFromSampleHash(ids[0]),
																														(byte) chr,
																														start + ((stop - start) / 2));
						}
					} else {
						colorCode = 0;
						// TODO color by value
					}
					inLine[3] = colorCode + "";
				}
				if (inLine[2] != null || isHistogram()) {
					v.add(inLine);
				}
			}
		}
		return v;
	}

	public void setHistogram(boolean b) {
		isHistPlot = b;
		size = b ? 0 : (size == 0 ? DEFAULT_SIZE : size);
	}

	public boolean isHistogram() {
		return isHistPlot;
	}

	public int getSelectedDataHash() {
		return selectedDataHash.hashCode();
	}

	public void updateGUI() {
		twoDPanel.paintAgain();
	}

	public void updateColorKey(Hashtable<String, String> hash) {
		colorKeyPanel.updateColorKey(hash);
	}


	public void generateShortcutMenus() {
		JPanel colorClass = colorKeyPanel.getClassVariablesPanel();
		Component[] components = colorClass.getComponents();
		MouseListener mouseListenerForRadios = new MouseListener() {
			@Override
			public void mouseReleased(MouseEvent e) {}

			@Override
			public void mousePressed(MouseEvent e) {}

			@Override
			public void mouseExited(MouseEvent e) {}

			@Override
			public void mouseEntered(MouseEvent e) {}

			@Override
			public void mouseClicked(MouseEvent e) {
				JRadioButton source;
				JPopupMenu menu;
				String annotation;
				// int annotationIndex;

				source = (JRadioButton) e.getSource();
				annotation = source.getText();

				if (e.getButton() == MouseEvent.BUTTON3) {

					menu = new JPopupMenu();
					menu.setName("Action");

					menu.add(new AbstractAction("Delete: " + annotation) {
						private static final long serialVersionUID = 1L;

						@Override
						public void actionPerformed(ActionEvent e1) {
							String annotation = e1.getActionCommand();
							annotation = annotation.substring("Delete: ".length());
							removeColorKeyHandler(annotation);
						}

					});
					menu.show(source, e.getX(), e.getY());
				}
			}
		};
		for (Component comp : components) {
			if (comp instanceof JRadioButton) {
				comp.addMouseListener(mouseListenerForRadios);
			}
		}
	}

	public void removeColorKeyHandler(String colorKeyName) {
		sampleData.removeColorKey(colorKeyName);
		reloadSampleDataUI();
	}

	public void showSpecificFile(String filename, int colForX, int colForY) {
		if (proj != null && promptOnClose) {
			String[] prevFiles = proj.TWOD_LOADED_FILENAMES.getValue();
			if (Arrays.binarySearch(prevFiles, filename) < 0) {
				// the supplied file was not found so load it
				loadFile(filename); // load the file
			}
		}

		updateTreeForNewData();

		if (colForX >= 0) {
			// select the x axis
			tree.performCheckBoxAction(dataColumnsHash.get(filename)[colForX], ItemEvent.SELECTED);
		}

		if (colForY >= 0) {
			// select the x axis
			tree.performCheckBoxAction(dataColumnsHash.get(filename)[colForY], ItemEvent.SELECTED);
		}

		updateGUI();
		tree.expandRow(dataKeys.indexOf(filename));
		twoDPanel.paintAgain();
	}

	@Override
	public void windowActivated(WindowEvent e) {}

	@Override
	public void windowClosed(WindowEvent e) {}

	@Override
	public void windowClosing(WindowEvent e) {
		if (e == null) {
			GuiManager.disposeOfParentFrame(this);
			return;
		}
		if (proj == null || !promptOnClose) { // TODO add flag for dispose or hide on close
			GuiManager.hideParentFrame(this);
			return;
		}

		String[] options;
		int choice;
		String filenames, selections = "", message;

		filenames = ArrayUtils.toStr(ArrayUtils.toStringArray(dataKeys), ";");

		// find the selected nodes in the plot and create a string from them delimited by ;
		String[][] selectedNodes = tree.getSelectionValues();
		for (int i = 0; i < selectedNodes.length; i++) {
			if (selectedNodes[i] != null && selectedNodes[i].length == 2 && selectedNodes[i][0] != null
					&& selectedNodes[i][1] != null) {
				selections += selectedNodes[i][0] + "|" + selectedNodes[i][1];
				if (i < selectedNodes.length - 1) {
					selections += ";";
				}
			}
		}
		selections = selections.replaceAll(proj.PROJECT_DIRECTORY.getValue(), "");

		options = new String[] {"Yes", "No", "Cancel"};
		message = "";

		if (!ArrayUtils.toStr(proj.TWOD_LOADED_FILENAMES.getValue(), ";").equals(filenames)) {
			if (filenames.equals("")) {
				message = "All files have been unloaded from 2D Plot.";
			} else {
				message = "A different set of files have been loaded/unloaded from 2D Plot.";
			}
		}

		if (!ArrayUtils.toStr(proj.TWOD_LOADED_VARIABLES.getValue(), ";").equals(selections)) {
			if (!Strings.isNullOrEmpty(message)) {
				if (!Strings.isNullOrEmpty(selections)) {
					message = message.substring(0, message.length() - 1)
										+ ", and new variables have been selected.";
				} else {
					message = message.substring(0, message.length() - 1)
										+ ", and all variables have been unselected.";
				}
			} else {
				if (!Strings.isNullOrEmpty(selections)) {
					message = "New variables have been selected.";
				} else {
					message = "All variables have been unselected from 2D Plot.";
				}
			}
		}

		if (message.equals("")) {
			GuiManager.disposeOfParentFrame(this);
		} else {
			System.out.println("message: '" + message + "'");
			choice = JOptionPane.showOptionDialog(null,
																						message
																								+ " Would you like to keep this configuration for the next time 2D Plot is loaded?",
																						"Preserve 2D Plot workspace?",
																						JOptionPane.YES_NO_CANCEL_OPTION,
																						JOptionPane.QUESTION_MESSAGE, null, options,
																						options[0]);
			if (choice == JOptionPane.YES_OPTION) {
				proj.TWOD_LOADED_FILENAMES.setValue(ArrayUtils.toStringArray(dataKeys));
				proj.TWOD_LOADED_VARIABLES.setValue(selections.split(";"));
				proj.saveProperties();
				GuiManager.disposeOfParentFrame(this);
			} else if (choice == 1) {
				GuiManager.disposeOfParentFrame(this);
			} else if (choice == -1 || choice == 2) {
				// keep frame alive
			}
		}
	}

	@Override
	public void windowDeactivated(WindowEvent e) {}

	@Override
	public void windowDeiconified(WindowEvent e) {}

	@Override
	public void windowIconified(WindowEvent e) {}

	@Override
	public void windowOpened(WindowEvent e) {}

	int valueChangedCounter = 0;

	@Override
	public void valueChanged(TreeSelectionEvent e) {
		// if (valueChangedCounter < 3) {
		// valueChangedCounter++;
		// return;
		// }
		// valueChangedCounter = 0;
		twoDPanel.setPointsGeneratable(true);
		// twoDPanel.createImage(); // calling paintAgain sets image == null, so why call 'createImage'
		// right before doing so?
		twoDPanel.paintAgain();
	}


	public static class ScreenToCapture {
		String outputName;
		String dataXFile;
		String dataYFile;
		String colorFile;
		String title = "TEST";
		int xDataIndex, yDataIndex, colorIndex;
		int xIDIndex, yIDIndex, colorIDIndex;
		float minX, minY, maxX, maxY;
		boolean hideExcluded;
		boolean isHistogram;
		boolean createColorKey;
		boolean includeColorKey;

		public ScreenToCapture(String[] files, int[] dataIndices, int[] idIndices,
													 float[] displayWindow, boolean excluded, boolean colorKey,
													 boolean appendColorKey, boolean hist, String outputFilename) {
			dataXFile = files[0];
			dataYFile = files[1];
			colorFile = files[2];
			xDataIndex = dataIndices[0];
			yDataIndex = dataIndices[1];
			colorIndex = dataIndices[2];
			xIDIndex = idIndices[0];
			yIDIndex = idIndices[1];
			colorIDIndex = idIndices[2];
			minX = displayWindow[0];
			maxX = displayWindow[1];
			minY = displayWindow[2];
			maxY = displayWindow[3];
			hideExcluded = excluded;
			isHistogram = hist;
			createColorKey = colorKey;
			includeColorKey = appendColorKey;
			outputName = outputFilename;
		}
	}

	public static final String COMMAND_TWO_D_SCREENSHOTS = "twoDScreenshots";

	private static List<String> parseControlFile(String filename, Logger log) {
		List<String> params;

		params = Files.parseControlFile(filename,
																		COMMAND_TWO_D_SCREENSHOTS,
																		new String[] {
																									"proj=",
																									"dir=",
																									"# Per-line:",
																									"# after being set once, all tags can be duplicated in a subsequent line by omitting them; the only necessary tags are {xDataColumn=<> and yDataColumn=<>}. If the output tag is not specified for each line, the created file will be named as <xFilename>_<xDataColumn>_<yFilename>_<yDataColumn>",
																									"# fileX=<> fileY=<> xDataColumn=<> yDataColumn=<> xIDColumn=<> yIDColumn=<> colorFile=<> colorDataColumn=<> colorIDColumn=<> minX=<> minY=<> maxX=<> maxY=<> hideExcluded=True/False colorKey=True/False includeColorKey=True/False isHistogram=True/False output=<>",},
																		log);

		return params;
	}

	private static final String SCRN_TAG_OUTPUT = "output";
	private static final String SCRN_TAG_TITLE = "title";

	private static final String SCRN_TAG_COLOR_FILE = "colorFile";
	private static final String SCRN_TAG_COLOR_DATA_COL = "colorDataColumn";
	private static final String SCRN_TAG_COLOR_ID_COL = "colorIDColumn";

	private static ArrayList<ScreenToCapture> condenseCtrlFile(java.util.List<String> ctrlLines,
																														 boolean fail) {
		// TODO refactor these into an enum or objects, defining:
		// - is required?
		// - use prev if missing?

		HashSet<String> tagSet = new HashSet<String>();
		tagSet.add("fileX");
		tagSet.add("fileY");
		tagSet.add("xDataColumn");
		tagSet.add("yDataColumn");
		tagSet.add("xIDColumn");
		tagSet.add("yIDColumn");
		tagSet.add("minX");
		tagSet.add("minY");
		tagSet.add("maxX");
		tagSet.add("maxY");
		tagSet.add(SCRN_TAG_COLOR_FILE);
		tagSet.add(SCRN_TAG_COLOR_DATA_COL);
		tagSet.add(SCRN_TAG_COLOR_ID_COL);
		tagSet.add("hideExcluded");
		tagSet.add("isHistogram");
		tagSet.add("colorKey");
		tagSet.add("includeColorKey");
		tagSet.add(SCRN_TAG_OUTPUT);
		tagSet.add(SCRN_TAG_TITLE);

		HashSet<String> notRqrd = new HashSet<>();
		notRqrd.add(SCRN_TAG_OUTPUT);
		notRqrd.add(SCRN_TAG_TITLE);
		notRqrd.add(SCRN_TAG_COLOR_ID_COL);
		notRqrd.add(SCRN_TAG_COLOR_DATA_COL);
		notRqrd.add(SCRN_TAG_COLOR_FILE);

		HashMap<String, ArrayList<String>> tagValues = new HashMap<String, ArrayList<String>>();
		for (String tagKey : tagSet) {
			tagValues.put(tagKey, new ArrayList<String>());
		}

		HashSet<String> lineTagEntries = new HashSet<String>();
		for (String line : ctrlLines) {
			String[] lineTags = line.trim().split(ext.REGEX_TO_SPLIT_SPACES_NOT_IN_QUOTES);
			for (String lineTag : lineTags) {
				String tagKey = lineTag.split("=")[0];
				String tagValue = lineTag.split("=").length > 1 ? lineTag.split("=")[1] : "-1";

				if (!tagSet.contains(tagKey)) {
					System.err.println("Error - argument \"" + tagKey + "\" is not valid.");
					if (fail) {
						throw new IllegalArgumentException("Error - argument \"" + tagKey + "\" is not valid.");
					}
				} else {
					if (lineTagEntries.contains(tagKey)) {
						System.err.println("Error - duplicate argument \"" + tagKey + "\" detected.");
						if (fail) {
							throw new IllegalArgumentException("Error - duplicate argument \"" + tagKey
																								 + "\" detected.");
						}
					} else {
						lineTagEntries.add(tagKey);
						tagValues.get(tagKey).add(tagValue);
					}
				}

			}

			for (String tag : tagSet) {
				if (!lineTagEntries.contains(tag)) {
					ArrayList<String> values = tagValues.get(tag);
					if (values.isEmpty()) {
						if (notRqrd.contains(tag)) {
							values.add(null);
						} else {
							System.err.println("Error - missing argument \"" + tag + "\" detected.");
							if (fail) {
								throw new IllegalArgumentException("Error - missing argument \"" + tag
																									 + "\" detected.");
							}
						}
					} else {
						if (tag.equals(SCRN_TAG_OUTPUT) || tag.equalsIgnoreCase(SCRN_TAG_TITLE)) {
							values.add(null);
						}
						values.add(values.get(values.size() - 1));
					}
				}
			}

			lineTagEntries.clear();
		}

		ArrayList<ScreenToCapture> caps = new ArrayList<TwoDPlot.ScreenToCapture>();
		int capCnt = ctrlLines.size();
		for (int i = 0; i < capCnt; i++) {
			String[] files = new String[3];
			int[] idCols = new int[3];
			int[] dataCols = new int[3];
			float[] window = new float[4];
			boolean hideExcludes = true;
			boolean isHistogram = false;
			boolean colorKey = false;
			boolean inclKey = false;
			String output = null;

			files[0] = tagValues.get("fileX").get(i);
			files[1] = tagValues.get("fileY").get(i);
			files[2] = tagValues.get("colorFile").get(i);

			dataCols[0] = Integer.parseInt(tagValues.get("xDataColumn").get(i));
			dataCols[1] = Integer.parseInt(tagValues.get("yDataColumn").get(i));
			dataCols[2] = tagValues.get("colorDataColumn").get(i) == null
																																	 ? -1
																																	 : Integer.parseInt(tagValues.get("colorDataColumn")
																																															 .get(i));

			idCols[0] = Integer.parseInt(tagValues.get("xIDColumn").get(i));
			idCols[1] = Integer.parseInt(tagValues.get("yIDColumn").get(i));
			idCols[2] = tagValues.get("colorIDColumn").get(i) == null
																															 ? -1
																															 : Integer.parseInt(tagValues.get("colorIDColumn")
																																													 .get(i));

			window[0] = Float.parseFloat(tagValues.get("minX").get(i));
			window[1] = Float.parseFloat(tagValues.get("maxX").get(i));
			window[2] = Float.parseFloat(tagValues.get("minY").get(i));
			window[3] = Float.parseFloat(tagValues.get("maxY").get(i));

			hideExcludes = Boolean.parseBoolean(tagValues.get("hideExcluded").get(i));
			isHistogram = Boolean.parseBoolean(tagValues.get("isHistogram").get(i));
			colorKey = Boolean.parseBoolean(tagValues.get("colorKey").get(i));
			inclKey = Boolean.parseBoolean(tagValues.get("includeColorKey").get(i));

			output = tagValues.get("output").get(i);

			ScreenToCapture sc = new ScreenToCapture(files, dataCols, idCols, window, hideExcludes,
																							 colorKey, inclKey, isHistogram, output);
			caps.add(sc);
		}

		return caps;
	}

	public void createScreenshots(String baseDir, List<ScreenToCapture> screens) {
		generatingScreenshots = true;
		HashSet<String> dataFiles = new HashSet<String>();

		for (ScreenToCapture cap : screens) {
			if (cap.dataXFile != null && Files.exists(baseDir + cap.dataXFile)) {
				dataFiles.add(cap.dataXFile);
			}
			if (cap.dataYFile != null && Files.exists(baseDir + cap.dataYFile)) {
				dataFiles.add(cap.dataYFile);
			}
			if (cap.colorFile != null && Files.exists(baseDir + cap.colorFile)) {
				dataFiles.add(cap.colorFile);
			}
		}

		for (String file : dataFiles) {
			loadFile(baseDir + file);
		}

		twoDPanel.setPointsGeneratable(true);
		twoDPanel.paintAgain();
		updateTreeForNewData();

		for (ScreenToCapture screencap : screens) {
			hideExcludes = screencap.hideExcluded;
			setHistogram(screencap.isHistogram);

			twoDPanel.setTitle(screencap.title);
			twoDPanel.setDisplayTitle(screencap.title != null);

			twoDPanel.setForcePlotXmin(screencap.minX);
			twoDPanel.setForcePlotXmax(screencap.maxX);
			twoDPanel.setForcePlotYmin(screencap.minY);
			twoDPanel.setForcePlotYmax(screencap.maxY);

			boolean colorLoaded = false;
			if (screencap.colorFile != null && Files.exists(baseDir + screencap.colorFile)) {
				loadColor(baseDir, screencap);
				colorLoaded = true; // vulnerable to issues actually loading color file, but good enough for
														// now? TODO need some way to specify HeatMap/Genotype coloration
			}

			if (screencap.dataXFile != null) {
				tree.performCheckBoxAction(screencap.dataXFile,
																	 dataColumnsHash.get(baseDir
																											 + screencap.dataXFile)[screencap.xDataIndex],
																	 ItemEvent.SELECTED);
			}
			if (screencap.dataYFile != null) {
				// 07/23/15 BRC -- this was screencap.dataXFile, which doesn't make sense (changed it to
				// dataYFile). If it breaks, put it back.
				tree.performCheckBoxAction(screencap.dataYFile,
																	 dataColumnsHash.get(baseDir
																											 + screencap.dataYFile)[screencap.yDataIndex],
																	 ItemEvent.SELECTED);
			}

			twoDPanel.setChartType(AbstractPanel.SCATTER_PLOT_TYPE);
			colorKeyPanel.getClassRadioButtons()[colorLoaded
																											? (colorKeyPanel.getClassRadioButtons().length
																											- 1)
																											: 0].setSelected(true);

			twoDPanel.createImage();

			int count = 1;
			String basename = screencap.outputName;
			boolean addPng = !screencap.outputName.endsWith(".png");
			if (basename == null) {
				if (screencap.dataXFile != null) {
					basename += ext.rootOf(screencap.dataXFile, true);
					basename += "_";
					basename += dataColumnsHash.get(baseDir + screencap.dataXFile)[screencap.xDataIndex];
				}
				if (screencap.dataYFile != null) {
					if (!basename.equals("")) {
						basename += "_";
					}
					basename += ext.rootOf(screencap.dataYFile, true);
					basename += "_";
					basename += dataColumnsHash.get(baseDir + screencap.dataYFile)[screencap.yDataIndex];
				}
			}
			String screenname = basename;
			while ((new File(baseDir + ext.replaceWithLinuxSafeCharacters(screenname, true)
											 + (addPng ? ".png" : ""))).exists()) {
				screenname = basename + "_v" + count;
				count++;
			}

			screenname = baseDir + ext.replaceWithLinuxSafeCharacters(screenname, true)
									 + (addPng ? ".png" : "");

			if (screencap.createColorKey) {
				// Use a JFrame for it's 'pack()' method - this shrinks colorKeyPanel to the minimum
				// required dimensions
				JFrame frame = new JFrame();
				// then change the layout (briefly) to disable WrapLayout's line-wrapping
				colorKeyPanel.classValuesPanel.setLayout(new FlowLayout(FlowLayout.CENTER, 0, 0));
				frame.add(colorKeyPanel.classValuesPanel);
				frame.pack();
				BufferedImage bi = new BufferedImage(colorKeyPanel.classValuesPanel.getWidth(),
																						 colorKeyPanel.classValuesPanel.getHeight(),
																						 BufferedImage.TYPE_INT_ARGB);
				colorKeyPanel.classValuesPanel.paint(bi.createGraphics());
				// then reset and dispose of extra resources
				colorKeyPanel.classValuesPanel.setLayout(new org.genvisis.cnv.gui.WrapLayout(
																																										 FlowLayout.CENTER,
																																										 0, 0));
				frame.removeAll();
				frame.dispose();
				frame = null;
				if (screencap.includeColorKey) {
					int totW, totH;
					totW = Math.max(bi.getWidth(), twoDPanel.getImage().getWidth());
					totH = bi.getHeight() + twoDPanel.getImage().getHeight();
					BufferedImage img = new BufferedImage(totW, totH, BufferedImage.TYPE_INT_ARGB);
					Graphics2D g = img.createGraphics();
					g.drawImage(twoDPanel.getImage(), 0, 0, null);
					int x = (int) ((.5 * twoDPanel.getWidth()) - (.5 * bi.getWidth()));
					g.drawImage(bi, x, twoDPanel.getHeight(), null);
					// g.drawImage(bi, 0, twoDPanel.getHeight(), null);
					try {
						ImageIO.write(img, "png", new File(screenname));
					} catch (IOException ie) {
						JOptionPane.showMessageDialog(null, "Error while trying to save the plot");
					}
				} else {
					twoDPanel.screenCapture(screenname);

					String legName = baseDir + ext.replaceWithLinuxSafeCharacters(basename, true)
													 + "_legend.png";
					count = 1;
					while ((new File(legName).exists())) {
						legName = baseDir + ext.replaceWithLinuxSafeCharacters(basename, true) + "_legend_v"
											+ count + ".png";
						count++;
					}
					try {
						ImageIO.write(bi, "png", new File(legName));
					} catch (IOException ie) {
						JOptionPane.showMessageDialog(null, "Error while trying to save the plot legend");
					}
				}
			} else {
				twoDPanel.screenCapture(screenname);
			}

		}

	}

	private void loadColor(String baseDir, ScreenToCapture screencap) {
		Hashtable<String, String> colorData = HashVec.loadFileToHashString(baseDir
																																					 + screencap.colorFile,
																																			 screencap.colorIDIndex,
																																			 new int[] {screencap.colorIndex},
																																			 "\t", true);
		this.colorData = new HashMap<String, Integer>();
		for (java.util.Map.Entry<String, String> entry : colorData.entrySet()) {
			this.colorData.put(entry.getKey(), Integer.valueOf(entry.getValue()));
		}
	}

	private byte getColorForScreenshot(String id) {
		return (byte) (colorData == null || colorData.get(id) == null ? 0 : colorData.get(id)
																																								 .intValue());
	}



	/**
	 * @return Empty string if successful, otherwise an error message.
	 */
	public String loadFile(String filename) {
		BufferedReader reader;
		String[] header, line;
		String readBuffer;
		int lineLength;

		if (dataKeys.contains(filename)) {
			return "";
		}
		if (!validExts.isEmpty()) {
			boolean found = false;
			for (String s : validExts) {
				if (filename.endsWith(s)) {
					found = true;
					break;
				}
			}
			if (!found) {
				String e = "Error - extension of file {" + filename
									 + "} is invalid.  Valid extensions are: {" + validExts.toString() + "}";
				log.reportError(e);
				return e;
			}
		}
		try {
			long t1 = System.currentTimeMillis();
			reader = new BufferedReader(new FileReader(filename));

			readBuffer = reader.readLine();
			String delim = ext.determineDelimiter(readBuffer);
			header = readBuffer.trim().split(delim);

			// TODO add column selection and data filtering / transformations
			// don't load invalid data (or discard once loaded and display warning)

			lineLength = header.length;

			ArrayList<String[]> data = new ArrayList<String[]>();
			boolean[] valid = ArrayUtils.booleanArray(header.length, true);

			String tempLine = "";
			int cnt = 0;
			while ((tempLine = reader.readLine()) != null) {
				if ("".equals(tempLine)) {
					continue;
				}
				line = tempLine.trim().split(delim);
				if (line.length < lineLength) {
					String e = "Error on line " + cnt + " of file " + filename + ": expected "
										 + header.length + " columns, but only found " + line.length;
					log.reportError(e);
					reader.close();
					return e;
				}
				valid = ArrayUtils.booleanArrayAnd(valid, validate(line));
				data.add(line);
				cnt++;
			}
			reader.close();

			filenameMap.put(new File(filename).getName(), filename);
			addDataHeader(filename, header);
			dataHash.put(filename, data);
			validColumnsHash.put(filename, valid);
			dataKeys.add(filename);

			log.reportTimeElapsed("Read file " + filename + " in: ", t1);

		} catch (FileNotFoundException fnfe) {
			String e = "Error: file \"" + filename + "\" not found in current directory";
			log.reportError(e);
			return e;
		} catch (IOException ioe) {
			String e = "Error reading file \"" + filename + "\"";
			log.reportError(e);
			return e;
		}

		return "";
	}

	private void addDataHeader(String filename, String[] header) {
		int[] linkKeyIndices;

		dataColumnsHash.put(filename, header);

		linkKeyIndices = ext.indexFactors(LINKERS, header, false, true, false, log);

		if (linkKeyIndices[DNA_INDEX_IN_LINKERS] == -1) {
			log.report("Sample ID not automatically identified for file '" + filename
								 + "'.");
			linkKeyIndices[DNA_INDEX_IN_LINKERS] = -1;
		}

		columnMetaData.put(filename, new HashMap<Integer, ColumnMetaData>());
		fileColumnDataPipes.put(filename, new HashMap<>());

		String region, chr, start, stop, lbl;
		String[] regSplit;
		// REGION
		if (linkKeyIndices[REGION_INDEX_IN_LINKERS] == -1 || linkKeyIndices[CHR_INDEX_IN_LINKERS] == -1
				|| linkKeyIndices[POS_INDEX_IN_LINKERS] == -1
				|| linkKeyIndices[STOP_POS_INDEX_IN_LINKERS] == -1) {
			// columnMetaData
			for (int i = 0; i < header.length; i++) {
				Matcher matcher = headerPattern.matcher(header[i]);
				if (matcher.matches()) {
					lbl = matcher.group(1);
					chr = matcher.group(2);
					region = matcher.group(3);
					regSplit = region.split("-");
					start = regSplit[0];
					stop = regSplit[1];
					columnMetaData.get(filename).put(i, new ColumnMetaData(lbl, chr, region, start, stop));
				}
			}
		}
		linkerIndices.put(filename, linkKeyIndices);
	}

	private boolean[] validate(String[] data) {
		boolean[] result = ArrayUtils.booleanArray(data.length, true);
		for (int i = 0; i < data.length; i++) {
			if (!ext.isValidDouble(data[i])) {
				boolean res = false;
				for (String miss : MISSING_VALUES) {
					if (miss.equals(data[i])) {
						res = true;
						break;
					}
				}
				result[i] = res;
			}
		}
		return result;
	}

	private void updateTreeForNewData() {
		String[] namesOfBranches;
		String[] branchHandles;
		String[][] namesOfNodes;

		String[][] treeFileVariableNameLookup;
		treeFileVariableNameLookup = new String[dataKeys.size()][];
		for (int i = 0; i < dataKeys.size(); i++) {
			String key = dataKeys.get(i);
			treeFileVariableNameLookup[i] = dataColumnsHash.get(key);
			if (tree == null) {
				System.err.println("Error - CheckBoxTree was null!  TODO check if not setting the treeListener causes issues...");
				namesOfBranches = new String[1];
				branchHandles = new String[1];
				namesOfNodes = new String[1][];
				namesOfBranches[0] = ext.removeDirectoryInfo(key);
				branchHandles[0] = key;
				namesOfNodes[0] = treeFileVariableNameLookup[i];
				tree = new CheckBoxTree(namesOfBranches, branchHandles, namesOfNodes,
																validColumnsHash.get(key), /* isHistPlot ? 1 : */2);
			} else {
				boolean found = false;
				for (int j = 0; j < tree.getModel().getChildCount(tree.getModel().getRoot()); j++) {
					if (tree.getModel().getChild(tree.getModel().getRoot(), j).toString()
									.equals(ext.removeDirectoryInfo(key))) {
						found = true;
					}
				}
				if (!found) {
					tree.addNode(ext.removeDirectoryInfo(key), key, treeFileVariableNameLookup[i],
											 validColumnsHash.get(key), checkBoxMouseListener);
				}
			}
		}
	}

	/**
	 * Create the GUI and show it. For thread safety, this method should be invoked from the
	 * event-dispatching thread.
	 *
	 * @param proj Project
	 * @param show setVisible
	 * @param promptOnClose prompt to save files/vars on close
	 * @param fileExts allowed file extensions for files (null for all files)
	 * @return
	 */
	public static TwoDPlot createGUI(Project proj, boolean show, boolean promptOnClose,
																	 String... fileExts) {
		boolean headless = GraphicsEnvironment.isHeadless();
		JFrame frame = null;
		TwoDPlot twoDPlot = new TwoDPlot(proj, promptOnClose, fileExts);
		twoDPlot.setOpaque(true); // content panes must be opaque
		if (!headless) {
			frame = new JFrame("Genvisis - 2D Plot"
												 + (proj != null ? " - " + proj.PROJECT_NAME.getValue() : ""));
			frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
			frame.setJMenuBar(twoDPlot.menuBar());
			frame.setContentPane(twoDPlot);
			frame.addWindowListener(twoDPlot);
			frame.setMinimumSize(new Dimension(20, 20));
			UITools.setSize(frame, 0.75, 0.75);
			// frame.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH);

			// Display the window.
			frame.pack();
			UITools.centerComponent(frame);
			frame.setVisible(show);
		}
		return twoDPlot;
	}

	/*
	 * CHR POS START STOP LABEL
	 */
	public String[][] getCurrentColumnMetaData() {
		ArrayList<String[]> resultList = new ArrayList<String[]>();
		HashSet<String> inclSet = new HashSet<String>();
		String[][] selectionValues = tree.getSelectionValues();
		for (String[] selectionValue : selectionValues) {
			String[] result = null;
			if (selectionValue != null && selectionValue[0] != null) {
				String file = selectionValue[0];
				int col = Integer.parseInt(selectionValue[1]);
				ColumnMetaData metaData = getColumnMetaData(file, col);
				if (metaData != null) {
					result = new String[] {metaData.chr, metaData.region, metaData.start, metaData.stop,
																 metaData.lbl};
				} else {
					String colName = dataColumnsHash.get(file)[col];
					Matcher m = headerPattern.matcher(colName);
					if (m.matches()) {
						result = new String[] {m.group(2), m.group(3), m.group(3).split("-")[0],
																	 m.group(3).split("-")[1], m.group(1)};
					} else {
						result = null;
					}
				}
				if (result != null) {
					if (inclSet.contains(result[0] + result[1] + result[2] + result[3] + result[4])) {
						result = null;
					} else {
						inclSet.add(result[0] + result[1] + result[2] + result[3] + result[4]);
					}
				}
				resultList.add(result);
			}
		}
		String[][] returnResults = new String[resultList.size()][];
		for (int i = 0; i < returnResults.length; i++) {
			returnResults[i] = resultList.get(i);
		}
		return returnResults;
	}

	public ColumnMetaData getColumnMetaData(String filename, int index) {
		return columnMetaData.get(filename).get(index);
	}

	public static void fromParameters(String filename, Logger log) {
		List<String> params = parseControlFile(filename, log);

		if (params != null) {
			String[] parts = params.get(0).split("=");
			final String projFile = parts.length == 1 ? null : parts[1];
			final String baseDir = params.get(1).split("=")[1];
			final ArrayList<ScreenToCapture> screens = condenseCtrlFile(params.subList(2, params.size()),
																																	true);
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				@Override
				public void run() {
					Project proj = projFile == null ? null : new Project(projFile);
					TwoDPlot tdp = createGUI(proj, false, false);
					tdp.createScreenshots(baseDir, screens);
					tdp.windowClosing(null);
				}
			});
		}
	}

	public static void main(String[] args) {
		// HashSet<String> tagSet = new HashSet<String>();
		// tagSet.add("fileXs");
		// tagSet.add("fileYs");
		// tagSet.add("xDataColumns");
		// tagSet.add("yDataColumns");
		// tagSet.add("xIDColumns");
		// tagSet.add("yIDColumns");
		// tagSet.add("colorFiles");
		// tagSet.add("colorDataColumns");
		// tagSet.add("colorIDColumns");
		// tagSet.add("excludeFiles");
		// tagSet.add("excludeDataColumns");
		// tagSet.add("excludeIDColumns");
		// tagSet.add("minXs");
		// tagSet.add("minYs");
		// tagSet.add("maxXs");
		// tagSet.add("maxYs");
		//
		// HashMap<String, String> tagValues = new HashMap<String, String>();
		//
		// String usage = "";
		//
		// System.out.println();
		// for (int i = 0; i<args.length; i++) {
		// if
		// (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help"))
		// {
		// System.err.println(usage);
		// System.exit(1);
		// } else if (tagSet.contains(args[i].split("=")[0])) {
		// tagValues.put(args[i].split("=")[0], args[i].split("=")[1]);
		// }
		// }
		//
		//
		//
		// fromParameters("D:/data/FarrarReparse/classification/twoDscreenshots.dat", new Logger());
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				createGUI(new Project(), true, false, null);
				// createAndShowGUI(new Project(cnv.Launch.getDefaultDebugProjectFile(true), false));
			}
		});
	}
}
