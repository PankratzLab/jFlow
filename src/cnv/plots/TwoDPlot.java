package cnv.plots;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.tree.TreePath;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;

import cnv.filesys.*;
import cnv.filesys.Project.StringListProperty;
import cnv.gui.CheckBoxTree;
import cnv.gui.ColorKeyPanel;
import cnv.gui.GuiManager;
import cnv.var.*;
import common.*;

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
//	public static final String MASK_MISSING = "Mask missing values";
//	public static final String UNMASK_MISSING = "Unmask missing values";
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String ADD_DATA_FILE = "Add Data File";
	public static final String REMOVE_DATA_FILE = "Remove Data";
	public static final String REMOVE_ALL = "Clear Data";
	public static final String CREATE_SCREENS = "Create Screenshots";
//	public static final String SET_AS_COLORKEY = "Set as Color Key";
//	public static final String SET_AS_LINKKEY = "Set as Link Key";
//	private static final String NO_VALUE_FOUND = ".";
//	public static final String[] BUTTONS = {ADD_DATA_FILE, REMOVE_DATA_FILE, REMOVE_ALL, CREATE_SCREENS};
//	public static final String[][] LINKERS = {
//			//TODO - Rohit: Removed Sample from first Linker. Confirm with Nathan if this is okay.
//		{"IndividualID", "ID", "IID", "UID", "UniqueID", "IndID"},
//		{"Family ID", "FamID", "FID"}, 
//		{"DNA/Sample", "DNA", "DNA#", "Sample", "LabID"}, 
//		{"MarkerName", "Marker", "SNP", "Variant", "VariantName"}, // will link to Scatter Plot
//		{"Region", "UCSC", "Band", "Arm"},	// will link to Trailer
//		{"Chromosome", "Chr"},	// secondary link to Trailer
//		{"Position", "Pos", "Start", "Begin"}, // secondary link to Trailer
//		{"Stop Position", "Stop", "End"} // secondary link to Trailer
//	};
	public static final String[][] LINKERS = {
	    Aliases.INDIVIDUAL_ID,
	    Aliases.FAMILY_ID,
	    Aliases.DNA,
	    Aliases.MARKER_NAMES,
	    Aliases.REGION,
	    Aliases.CHRS,
	    Array.combine(Aliases.POSITIONS, Aliases.POSITIONS_START),
	    Aliases.POSITIONS_STOP
	};
	public static final int IID_INDEX_IN_LINKERS = 0;
	public static final int FID_INDEX_IN_LINKERS = 1;
	public static final int DNA_INDEX_IN_LINKERS = 2;
	public static int MARKER_INDEX_IN_LINKERS = 3;
	public static int REGION_INDEX_IN_LINKERS = 4;
	public static int CHR_INDEX_IN_LINKERS = 5;
	public static int POS_INDEX_IN_LINKERS = 6;
	public static int STOP_POS_INDEX_IN_LINKERS = 7;
	
	public static final String[] MISSING_VALUES = {"."};
	
	/*
	 * regex to match and pull out column title, chromosome, and position
	 */
	private static String DATA_TITLE_REGEX = "(.+)_chr([0-2[XYM]]?[0-9[Y]]*):([\\d,]+-[\\d,]+).*";
	
	//	private JPanel bottomPanel;
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
	HashMap<String, ArrayList<String[]>> dataHash;
	HashMap<String, String[]> dataColumnsHash;
	HashMap<String, boolean[]> validColumnsHash;
	HashSet<String> validExts;
//	lbl, chr, region, start, stop
	HashMap<String, HashMap<Integer, String[]>> columnMetaData;
	String[][] treeFileVariableNameLookup;
	HashMap<String, int[]> linkerIndices;
	Logger log;
	Pattern headerPattern = Pattern.compile(DATA_TITLE_REGEX);
	private final boolean promptOnClose;
    private StringListProperty filesProperty;

    private volatile boolean generatingScreenshots;
    private HashMap<String, Integer> colorData;
    private JPanel treePanel;
    private JScrollPane scrollPane;
    private HashSet<String> selectedDataHash = new HashSet<String>();
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
		this(new Project(), true, null, null);
	}

	private TwoDPlot(Project proj, boolean promptOnClose, String[] fileExts, StringListProperty filenames) {
	    this.promptOnClose = promptOnClose;
		String[] previouslyLoadedFiles;
		
		this.proj = proj;
		log = proj.getLog();
		size = DEFAULT_SIZE;

		if (Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false), proj.JAR_STATUS.getValue())) {
			sampleData = proj.getSampleData(2, false);
		} else {
			sampleData = null;
		}
		validExts = new HashSet<String>();
		dataKeys = new ArrayList<String>();
		dataHash = new HashMap<String, ArrayList<String[]>>();
		dataColumnsHash = new HashMap<String, String[]>();
		validColumnsHash = new HashMap<String, boolean[]>();
		linkerIndices = new HashMap<String, int[]>();
		columnMetaData = new HashMap<String, HashMap<Integer,String[]>>();
		
		if (fileExts != null) {
    		for (String ext : fileExts) {
    		    validExts.add(ext);
    		}
		}
		
		this.filesProperty = filenames;
		previouslyLoadedFiles = filenames == null ? new String[0] : filenames.getValue();
		for (String previouslyLoadedFile : previouslyLoadedFiles) {
			loadFile(previouslyLoadedFile);
		}
		setLayout(new BorderLayout());
		
		twoDPanel = new TwoDPanel(this);

//		UIManager.put("PopupMenuUI", "CustomPopupMenuUI");

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

		twoDPanel.setPointsGeneratable(true);//zx
		twoDPanel.setExtraLayersVisible(new byte[] {99});
		updateGUI();
		
		twoDPanel.grabFocus();
		
		layeredPane.addComponentListener(new ComponentListener() {
			public void componentShown(ComponentEvent e) {}
			public void componentResized(ComponentEvent e) {
				flipButton.setBounds(70, layeredPane.getHeight() - 75, 38, 38);
				invXButton.setBounds(70, layeredPane.getHeight() - 35, 38, 13);
				invYButton.setBounds(55, layeredPane.getHeight() - 75, 13, 38);
			}
			public void componentMoved(ComponentEvent e) {}
			public void componentHidden(ComponentEvent e) {}
		});
		
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
    		    for (int i = 0; i < selectionValues.length; i++) {
    //                String file = selectionValues[i][0];
    //                int col = Integer.parseInt(selectionValues[i][1]);
                    selectedDataHash.add(selectionValues[i][0] + "***" + selectionValues[i][1]);
    		    }
    //			JCheckBox source;
    //			JPopupMenu menu;
    //			source = (JCheckBox) e.getSource();
    //			if (e.getButton() == MouseEvent.BUTTON3) {
    //
    //				menu = new JPopupMenu();
    //				menu.setName("Actions");
    //				menu.add(new AbstractAction("Set As Color Key") {
    //					private static final long serialVersionUID = 1L;
    //
    //					@Override
    //					public void actionPerformed(ActionEvent e1) {
    //						setColorKeyHandler(tree.getSelectionRows());
    //					}
    //				});
    //				menu.add(new AbstractAction("Set As Link Key") {
    //					private static final long serialVersionUID = 1L;
    //
    //					@Override
    //					public void actionPerformed(ActionEvent e1) {
    //						setLinkKeyHandler(tree.getSelectionRows());
    //					}
    //				});
    //				menu.show(source, e.getX(), e.getY());
    //			}
    		}
    
    		@Override
    		public void mousePressed(MouseEvent e) {
    		}
    
    		@Override
    		public void mouseReleased(MouseEvent e) {
    
    		}
    
    		@Override
    		public void mouseEntered(MouseEvent e) {
    		}
    
    		@Override
    		public void mouseExited(MouseEvent e) {
    		}
    	};

    private void generateFlipButton() {
		flipButton = new JButton(flip1_1);
		flipButton.setRolloverIcon(flip1_2);
		flipButton.setToolTipText("Inverts X & Y Axes");
		flipButton.setBorder(null);
		flipButton.setVisible(true);
		flipStatus = true;
		flipButton.addActionListener(new ActionListener() {
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
		menuItemOpen.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			    addFile();
			}
		});
		menu.add(menuItemOpen);
        menuItemRemove = new JMenuItem("Remove Selected Data", KeyEvent.VK_O);
        menuItemRemove.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                removeSelectedData();
            }
        });
        menu.add(menuItemRemove);
        menuItemRemoveAll = new JMenuItem("Remove All Data", KeyEvent.VK_O);
        menuItemRemoveAll.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                removeAllData();
            }
        });
        menu.add(menuItemRemoveAll);
		menuItemSave = new JMenuItem("Save Current Image", KeyEvent.VK_S);
		menuItemSave.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fileChooser = new JFileChooser(proj != null ? proj.PROJECT_DIRECTORY.getValue() : ".");
				int fileOpenActionSelected = fileChooser.showSaveDialog(TwoDPlot.this);
		        if (fileOpenActionSelected == JFileChooser.APPROVE_OPTION) {
		            File fileToOpen = fileChooser.getSelectedFile();
		            twoDPanel.screenCapture(fileToOpen.toString()+".png");
		        }
			}
		});
		menu.add(menuItemSave);
		menuItemScreens = new JMenuItem("Create Screenshots from File", KeyEvent.VK_S);
		menuItemScreens.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent e) {
		        createScreenshotsFromFile();
		    }
		});
		menu.add(menuItemScreens);
		menuItemExit = new JMenuItem("Close", KeyEvent.VK_C);
		menuItemExit.addActionListener(new ActionListener() {
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
				
//				tree = new CheckBoxTree(new String[0], new String[0], new String[0][], new boolean[0], /*isHistPlot ? 1 : */2);
//                updateTree();
//                
//                treePanel.remove(scrollPane);
//                scrollPane = new JScrollPane(tree);
//                treePanel.add(scrollPane, BorderLayout.CENTER);
//                tree.addTreeSelectionListener(TwoDPlot.this);
//                TwoDPlot.this.revalidate();
//                updateGUI();
				
				twoDPanel.paintAgain();
			}
		});
		menu.add(menuItemHist);
		
		
//		menu.add(new JMenuItem("Cut"));
//		menu.add(new JMenuItem("Copy"));
//		menu.add(new JMenuItem("Paste"));
//		menu.add(new JMenuItem("Paste Image"));
//		menu.add(new JMenuItem("Find"));
//		menu = new JMenu("Help");
//		menuBar.add(menu);
//		menu.add(new JMenuItem("Contents"));
//		menu.add(new JMenuItem("Search"));
//		menu.add(new JMenuItem("About"));
		return menuBar;
	}
	

    private Component getParentComponent() {
        Component c = TwoDPlot.this.getParent();
        while (!(c instanceof JFrame)) {
            c = c.getParent();
        }
        return c;
    }
	
	public void actionPerformed(ActionEvent ae) {
//		byte numberOfSelectedNodes;
		String command;
//		String[] keys;
		
		command = ae.getActionCommand();
		if (command.equals(ADD_DATA_FILE)) {
			addFile();
		} else if (command.equals(REMOVE_DATA_FILE)) {
		    removeSelectedData();
//			System.out.println(dataHash.size()+"\t"+namesHash.size()+"\t"+numericHash.size()+"\t"+treeFilenameLookup.size());
//			System.out.println("\n==== End =====\n");
//		} else if (command.equals(SET_AS_COLORKEY)) {
////			colorKeyVariables.add(tree.getSelectedPathComponentName());
////			colorKeyVariables.add(getNamesSelected()[0]);
//			System.out.println("getSelectedPathComponent: " + tree.getSelectedPathComponent() + "\t + getNamesSelected: " + Arrays.toString(getNamesSelected()));
//			setColorKey(tree.getSelectionRows());
//		} else if (command.equals(SET_AS_LINKKEY)) {
//			setLinkKeyHandler(tree.getSelectionRows());
		} else if (command.equals(REMOVE_ALL)) {
		    removeAllData();
		} else if (command.equals(CREATE_SCREENS)) {
		    createScreenshotsFromFile();
		} else {
			System.err.println("Error - unknown command '"+command+"'");
		}
	}
	
	private void createScreenshotsFromFile() {
        JFileChooser jfc = new JFileChooser(proj != null ? proj.PROJECT_DIRECTORY.getValue() : (new File(".")).toString());
        jfc.setMultiSelectionEnabled(true);
        if (jfc.showOpenDialog(TwoDPlot.this) == JFileChooser.APPROVE_OPTION) {
            File selFile = jfc.getSelectedFile();
            Vector<String> params = parseControlFile(selFile.getAbsolutePath(), log);
//          final String projFile = params.get(0).split("=")[1];
            final String baseDir = params.get(1).split("=")[1];
            final ArrayList<ScreenToCapture> screens = condenseCtrlFile(params.subList(2, params.size()), true);
            createScreenshots(baseDir, screens);
        }
    }

    private void addFile() {
        JFileChooser fileChooser;
        boolean found = false;
	    fileChooser = new JFileChooser(proj != null ? proj.PROJECT_DIRECTORY.getValue() : ".");
        if (fileChooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
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
        tree.deleteAllNodes();
        updateGUI();
	}
	
	public void removeSelectedData() {
        byte numberOfSelectedNodes;
        String[] keys;
        numberOfSelectedNodes = (byte) tree.getSelectedPathComponent();
        if (numberOfSelectedNodes != -1) {
            keys = HashVec.getKeys(dataHash); // it is better for keys to be a local variable rather than a global variable, otherwise it needs to be updated every time something is added or deleted
            tree.deleteSelectedNode();
            dataHash.remove(keys[numberOfSelectedNodes]);//TODO tree.getSelectionValues()[0][0] is not the branch to delete.
            dataColumnsHash.remove(keys[numberOfSelectedNodes]);
            linkerIndices.remove(keys[numberOfSelectedNodes]);
            validColumnsHash.remove(keys[numberOfSelectedNodes]);
            dataKeys.remove(keys[numberOfSelectedNodes]);
        }
	}
	
//	public void setColorKeyHandler(int[] selectedColorKey){
//		String[][] selectedNodes;
//
//		selectedNodes = tree.getSelectionValues();
//		sampleData.setColorKey(selectedNodes[1][0], selectedColorKey[0]);
//		reloadSampleDataUI();
//	}

//	public void setColorKey(int[] selectedColorKey) {
//		String[][] selectedNodes;
//		ArrayList<Integer> colorKeys;
//
//		selectedNodes = tree.getSelectionValues();
//		if (selectedColorKey.length == 1) {
//			if (colorKeyIndex.containsKey(selectedNodes[1][0])) {
//				colorKeys = colorKeyIndex.get(selectedNodes[1][0]);
//			} else {
//				colorKeyIndex.put(selectedNodes[1][0], colorKeys = new ArrayList<Integer>());
//			}
//			for (Integer colorKey : colorKeys) {
//				if (colorKey == selectedColorKey[0]) {
//					System.out.println("Error: Already sey as color key");
//					JOptionPane.showMessageDialog(null, "Error: Already sey as color key", "Error",
//							JOptionPane.ERROR_MESSAGE);
//					return;
//				}
//			}
//			colorKeys.add(selectedColorKey[0]);	// add to colorKeys
//			setColorKeyHandler(selectedNodes[1][0], selectedColorKey[0]);
//		}
//	}
//
//	public void setColorKeyHandler(String recentSelectionFile, int selectedColorKey) {
//		Hashtable<String, String> colorKeyValue;
//		int[] linkKeyColumnLabels;
//
//		linkKeyColumnLabels = keyIndices.get(recentSelectionFile);
//		colorKeyValue = new Hashtable<String, String>();
//		if (linkKeyIndex.containsKey(recentSelectionFile)) {
//			switch (linkKeyIndex.get(recentSelectionFile)) {
//			case DNA_INDEX_IN_LINKERS:
//				colorKeyValue = HashVec.loadFileToHashString(recentSelectionFile, new int[]{linkKeyColumnLabels[linkKeyIndex.get(recentSelectionFile)]}, new int[]{selectedColorKey-1}, false, "",true, false, false);
//				break;
//			case FID_INDEX_IN_LINKERS:
//				colorKeyValue = HashVec.loadFileToHashString(recentSelectionFile, new int[]{linkKeyColumnLabels[linkKeyIndex.get(recentSelectionFile)], linkKeyColumnLabels[IID_INDEX_IN_LINKERS]}, new int[]{selectedColorKey-1}, false, "",true, false, false);
//				colorKeyValue = createHashWithSampleID(colorKeyValue);	// colorkey value hash with key as sampleID
//				break;
//			case IID_INDEX_IN_LINKERS:
//				colorKeyValue = HashVec.loadFileToHashString(recentSelectionFile, new int[]{linkKeyColumnLabels[linkKeyIndex.get(recentSelectionFile)]}, new int[]{selectedColorKey-1}, false, "",true, false, false);
//				colorKeyValue = createHashWithSampleID(colorKeyValue);	// colorkey value hash with key as sampleID
//				break;
//			default:
//				System.out.println("Error: Unable to read color key values. Invalid link key.");
//				JOptionPane.showMessageDialog(null, "Error: Unable to read color key values. Invalid link key.", "Error", JOptionPane.ERROR_MESSAGE);
//				break;
//			}
//		}
//		addToSampleData(colorKeyValue, recentSelectionFile, selectedColorKey);
//	}


//	/**
//	 * Function to indentify the headers in the sample data file
//	 * @param header: a string containing all the headers read as string
//	 */
//	public int[] getSampleDataHeaders(String header) {
//		String[] headersArray;
//		int[] indices;
//
//		if (header.contains("\t")) {
//			headersArray = header.trim().split("\t",-1);
//		} else {
//			headersArray = header.trim().split("[\\s]+");
//		}
//		indices = ext.indexFactors(LINKERS, headersArray, false, true, false, log, false);
//
//		if (indices[0] == -1) {
//			log.report("ID linker not automatically identified for Sample Data. Assuming the first column.");
//			indices[0] = 0;
//		}
//		System.out.println("The header indices in Sample data are: " + Arrays.toString(headersArray));
//
//		return indices;
//	}

//	public void addToSampleData(Hashtable<String, String> colorKeyValue, String recentSelectionFile, int selectedColorKey) {
//		String sampleDatafilename;
//		BufferedReader reader;
//		BufferedWriter writer;
//		String[] inLineArry;
//		String bakFile;
//		String inLine;
//		String colorKeyHeader;
//		String[] keys;
//		boolean covar, negativeValues, largerThanByte;
//		String trav;
//
//		sampleDatafilename = proj.getFilename(Project.SAMPLE_DATA_FILENAME, false, false);
//
//		if (!Files.exists(sampleDatafilename, proj.getJarStatus())) {
//			JOptionPane.showMessageDialog(null, "Cannot add as a color key without an existing SampleData file", "Error", JOptionPane.ERROR_MESSAGE);
//			log.reportError("Cannot add as a color key without an existing SampleData file");
//			return;
//		}
//
//		reader = null;
//		writer = null;
//
//		System.out.println("Sample data: " + sampleDatafilename);
//		bakFile = proj.archiveFile(sampleDatafilename);	// create backup of sample data file
//		colorKeyHeader = namesHash.get(recentSelectionFile)[selectedColorKey-1];
//
//		covar = false;
//		negativeValues = false;
//		largerThanByte = false;
//		keys = HashVec.getKeys(colorKeyValue, false, false);
//		for (String key : keys) {
//			trav = colorKeyValue.get(key);
//			if (!ext.isMissingValue(trav) && !ext.isValidInteger(trav)) {
//				covar = true;
//			}
//			if (ext.isValidDouble(trav) && Double.parseDouble(trav) < 0) {
//				negativeValues = true;
//			}
//			if (ext.isValidDouble(trav) && Double.parseDouble(trav) > Byte.MAX_VALUE) {
//				largerThanByte = true;
//			}
//		}
//
//		if (covar) {
//			JOptionPane.showMessageDialog(null, "Variable '"+colorKeyHeader+"' contains a quantitative meaure and will be added as a COVAR in SampleData and not as a Class", "Warning", JOptionPane.ERROR_MESSAGE);
//		} else if (negativeValues) {
//			JOptionPane.showMessageDialog(null, "Variable '"+colorKeyHeader+"' contains negative values and will be added as a COVAR in SampleData and not as a Class", "Warning", JOptionPane.ERROR_MESSAGE);
//			covar = true;
//		} else if (largerThanByte) {
//			JOptionPane.showMessageDialog(null, "Variable '"+colorKeyHeader+"' contains values larger than 128 and will be added as a COVAR in SampleData and not as a Class", "Warning", JOptionPane.ERROR_MESSAGE);
//			covar = true;
//		}
//
//		try {
//			reader = new BufferedReader(new FileReader(bakFile));
//			writer = new BufferedWriter(new FileWriter(sampleDatafilename));
//			inLine = reader.readLine();
//			int samDataIndex = getSampleDataHeaders(inLine)[DNA_INDEX_IN_LINKERS];
//			inLine = inLine + "\t"+(covar?"Covar=":"Class=") + colorKeyHeader;
//			writer.write(inLine);	// write the headers
//			while(reader.ready()) {
//				writer.newLine();
//				inLine = reader.readLine();
//				if (inLine.contains("\t")) {
//					inLineArry = inLine.trim().split("\t",-1);
//				} else {
//					inLineArry = inLine.trim().split("[\\s]+");
//				}
//				if (colorKeyValue.containsKey(inLineArry[samDataIndex])) {
//					inLine = inLine + "\t" + colorKeyValue.get(inLineArry[samDataIndex]);
//				} else {
//					inLine = inLine + "\t" + NO_VALUE_FOUND;
//				}
//				writer.write(inLine);
//			}
//		} catch (FileNotFoundException e) {
//			System.out.println("Error: Sample Data backup file not found");
//		} catch (IOException e) {
//			System.out.println("Error: unable to read sample data backup file");
//		} finally {
//			closeStream(reader);
//			closeStream(writer);
////			twoDPanel.paintAgain();
//		}
//		reloadSampleDataUI();
//		System.out.println(colorKeyHeader.split(";")[0] + " set as color key and added to Sample Data");
//		JOptionPane.showMessageDialog(null, colorKeyHeader.split(";")[0] + " set as color key and added to Sample Data", "Information", JOptionPane.INFORMATION_MESSAGE);
//	}

	public void reloadSampleDataUI(){
		//Note: This line resets the sample data by setting it to null which will cause linkKeyIndex to be lost.
		//Make sure to set it back after resetting sample data
//		Hashtable<String, Integer> linkKeyIndexCopy = new Hashtable<String, Integer>(sampleData.getLinkKeyIndex());
		proj.resetSampleData();
		if (Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false), proj.JAR_STATUS.getValue())) {
			sampleData = proj.getSampleData(2, false);
		}
		sampleData = proj.getSampleData(2, false);
//		sampleData.setLinkKeyIndex(linkKeyIndexCopy);

		colorKeyPanel.updateSampleData(sampleData);
		colorKeyPanel.updateColorVariablePanel();

		twoDPanel.paintAgain();
		generateShortcutMenus();
	}
	
//	public Hashtable<String, String> createHashWithSampleID(Hashtable<String, String> colorKeyValue) {
//		Hashtable<String, String> colorKeyValueHash;
//
//		colorKeyValueHash = new Hashtable<String, String>();
//		for (String key : colorKeyValue.keySet()) {
//			 colorKeyValueHash.put(twoDPanel.getSampleData().lookup(key)[0], colorKeyValue.get(key));
//		}
//
//		return colorKeyValueHash;
//	}

//	public void setLinkKeyHandler(int[] selectedLinkKey) {
//		String selectedCol;
//
//		selectedCol = tree.getSelectionValues()[0][0];
//
//		if (selectedLinkKey.length == 1) {
//			if (!linkerIndices.containsKey(selectedCol)) {
//				JOptionPane.showMessageDialog(null, "There was a problem in your selection. Please select again", "Error", JOptionPane.ERROR_MESSAGE);
//				return;
//			}
//			sampleData.setLinkKey(selectedCol, selectedLinkKey[0]);
//		}
//	}

//	public void createLinkKeyToDataHash(String filename, int[] linkKeyColumnLabels) {
//		String inLine;
//		String[] inLineArray;
//		Hashtable<String, String[]> curFileLinkKeyDataHash = new Hashtable<String, String[]>();
//		try {
//			BufferedReader fileReader = new BufferedReader(new FileReader(filename));
//			inLine = fileReader.readLine(); // headers: skip this line
//
//			System.out.println("link column labels:" + Arrays.toString(linkKeyColumnLabels));
//			System.out.println("link key index: " + linkKeyIndex.get(filename));
//			while (fileReader.ready()) {
//				if (inLine.contains("\t")) {
//					inLineArray = fileReader.readLine().trim().split("\t", -1);
//				} else {
//					inLineArray = fileReader.readLine().trim().split("[\\s]+");
//				}
//				if (linkKeyIndex.containsKey(filename)) {
//					switch (linkKeyIndex.get(filename)) {
//					case DNA_INDEX_IN_LINKERS:
//						curFileLinkKeyDataHash.put(inLineArray[linkKeyColumnLabels[linkKeyIndex.get(filename)]], inLineArray);
//						break;
//					case FID_INDEX_IN_LINKERS:
//						curFileLinkKeyDataHash.put(inLineArray[linkKeyColumnLabels[linkKeyIndex.get(filename)]] + "\t" + inLineArray[linkKeyColumnLabels[IID_INDEX_IN_LINKERS]], inLineArray);
//						break;
//					case IID_INDEX_IN_LINKERS:
//						curFileLinkKeyDataHash.put(inLineArray[linkKeyColumnLabels[linkKeyIndex.get(filename)]], inLineArray);
//						break;
//					default:
//						System.out.println("Error: Invalid link key.");
//						break;
//					}
//				}
//			}
//			fileReader.close();
//		} catch (FileNotFoundException fnfe) {
//			System.err.println("Error: file \"" + filename + "\" not found in current directory");
//		} catch (IOException ioe) {
//			System.err.println("Error reading file \"" + filename + "\"");
//		}
//		if (!curFileLinkKeyDataHash.isEmpty())
//			linkKeyToDataHash.put(filename, curFileLinkKeyDataHash);
//		System.out.println("LinkKeyDataHash:" + linkKeyToDataHash.toString());
//	}

	public byte getPointSize() {
		return size;
	}
	
	public void setPointSize(byte sz) {
	    this.size = sz;
	}
	
	public String[] getNamesSelected() {
	    ArrayList<String> resultList = new ArrayList<String>();
	    String[][] metaData = getCurrentColumnMetaData();
	    String[][] selectedValues = tree.getSelectionValues();
	    for (int i = 0; i < metaData.length; i++) {
	        resultList.add(metaData[i] == null ? dataColumnsHash.get(selectedValues[i][0])[Integer.parseInt(selectedValues[i][1])] : metaData[i][4] + " chr" + metaData[i][0] + ":" + metaData[i][1]);
	    }
	    String[] result = resultList.toArray(new String[resultList.size()]);
	    return result;
	}

	public int[] getCurrentLinkKeyColumnIndices() {
		String[][] selectedValues;
		
		selectedValues = tree.getSelectionValues();
		if (selectedValues == null ||  selectedValues[0][0] == null) {
			return null;
		}
		return linkerIndices.get(selectedValues[0][0]);
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

	public ArrayList<String[]> getDataSelected(boolean includeColorKeyValue) {
		String[][] selectedNodes;
		ArrayList<String[]> dataOfSelectedFile;
		HashMap<String, String[]> xHash;
		HashMap<String, String> yHash;
		String[] inLine, outLine;
		int selectedColumn;
		String selectedFile;
		String[] keys;
		ArrayList<String[]> v;
		int currentClass;
		String[] ids;
		byte colorCode;
		int[] linkKeyColumnIndices;
		byte index;
		
		// file/collection name and column number of selected data collection
		selectedNodes = tree.getSelectionValues();
		v = new ArrayList<String[]>();
        xHash = new HashMap<String, String[]>();
        yHash = new HashMap<String, String>();
		if (selectedNodes[0][0] != null && linkerIndices.get(selectedNodes[0][0]) != null) {
			selectedColumn = Integer.parseInt(selectedNodes[0][1]);
			selectedFile = selectedNodes[0][0];
			dataOfSelectedFile = dataHash.get(selectedFile);
			currentClass = colorKeyPanel.getCurrentClass();
			if (currentClass < SampleData.BASIC_CLASSES.length && SampleData.BASIC_CLASSES[currentClass].equals(SampleData.HEATMAP)) {
				twoDPanel.setChartType(AbstractPanel.HEAT_MAP_TYPE);
			} else {
				twoDPanel.setChartType(AbstractPanel.SCATTER_PLOT_TYPE);
			}
			
			linkKeyColumnIndices = linkerIndices.get(selectedFile);
			index = (byte) (includeColorKeyValue? 4 : 3);
			outer : for (int i = 0; i < dataOfSelectedFile.size(); i++) {
				inLine = dataOfSelectedFile.get(i);
				
				if (sampleData != null && (hideExcludes || colorKeyPanel.getDisabledClassValues().size() > 0)) {
					for (int k = 0; k < 3; k++) {
						if (linkKeyColumnIndices[k] >= 0) {
							String curSample = inLine[linkKeyColumnIndices[k]];

							IndiPheno indi = sampleData.getIndiFromSampleHash(curSample);

							int sampColorKey = 0;
							if (indi != null) {
								sampColorKey = sampleData.determineCodeFromClass(currentClass, (byte) 0, indi, (byte) 0, 0);
							}
							if (hideExcludes && sampleData.individualShouldBeExcluded(curSample)) {
								continue outer;
							}
							if (colorKeyPanel.getDisabledClassValues().containsKey(currentClass + "\t" + sampColorKey)) {
								continue outer;
							}
						}
					}
				}
				outLine = new String[linkKeyColumnIndices.length + index];
				outLine[0] = inLine[0];
				outLine[1] = inLine[selectedColumn];
				for (int j = 0; j < linkKeyColumnIndices.length; j++) {
					if (linkKeyColumnIndices[j] >= 0) {
						outLine[j + index] = inLine[linkKeyColumnIndices[j]];
					}
				}
				xHash.put(inLine[linkerIndices.get(selectedFile)[IID_INDEX_IN_LINKERS]], outLine);
			}
    		if (selectedNodes.length > 1 && selectedNodes[1] != null && selectedNodes[1][0] != null && linkerIndices.get(selectedNodes[1][0]) != null) {
    			selectedColumn = Integer.parseInt(selectedNodes[1][1]);
    			selectedFile = selectedNodes[1][0];
    			dataOfSelectedFile = dataHash.get(selectedFile);
    			for (int i = 0; i < dataOfSelectedFile.size(); i++) {
    				inLine = dataOfSelectedFile.get(i);
    				yHash.put(inLine[linkerIndices.get(selectedFile)[IID_INDEX_IN_LINKERS]], inLine[selectedColumn]);
    			}
    		}
			keys = HashVec.getKeys(xHash, false, false);
			v = new ArrayList<String[]>();
			for (String key : keys) {
			    inLine = xHash.get(key);
			    boolean hasY = false;
				if (yHash.containsKey(key)) {
					inLine[2] = yHash.get(key);
					hasY = true;
				}
				if (hasY || isHistogram()) {
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
        				            chr = Integer.parseInt(metaData[0][0]);
        				            start = Integer.parseInt(metaData[0][2]);
        				            stop = Integer.parseInt(metaData[0][3]);
    				            }
    				            colorCode = sampleData.determineCodeFromClass(currentClass, (byte) 0, sampleData.getIndiFromSampleHash(ids[0]), (byte) chr, start + ((stop - start) / 2));
    				        }
    				    } else {
    				        colorCode = 0;
    				    }
    				    inLine[3] = colorCode + "";
    				}
    //					if (!inLine[3].equals("0"))
    				v.add(inLine);
				}
    		}
		}
		return v;
	}
	
	public void setHistogram(boolean b) {
        this.isHistPlot = b;
        this.size = b ? 0 : (size == 0 ? DEFAULT_SIZE : size);
    }

	public boolean isHistogram() {
	    return this.isHistPlot;
	}
	
    public int getSelectedDataHash() {
	    return selectedDataHash.hashCode();
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
//		generateShortcutMenus();
	}


	public void generateShortcutMenus() {
		JPanel colorClass = colorKeyPanel.getClassVariablesPanel();
		Component[] components = colorClass.getComponents();
		MouseListener mouseListenerForRadios = new MouseListener() {
			public void mouseReleased(MouseEvent e) {
			}

			public void mousePressed(MouseEvent e) {
			}

			public void mouseExited(MouseEvent e) {
			}

			public void mouseEntered(MouseEvent e) {
			}

			@Override
			public void mouseClicked(MouseEvent e) {
				JRadioButton source;
				JPopupMenu menu;
				String annotation;
//				int annotationIndex;

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

	public void removeColorKeyHandler(String colorKeyName){
		sampleData.removeColorKey(colorKeyName);
		reloadSampleDataUI();
	}

	public void showSpecificFile(Project proj, String filename, int colForX, int colForY) {
		if (filesProperty != null) {
    	    String[] prevFiles = filesProperty.getValue();
    		if(Arrays.binarySearch(prevFiles, filename) < 0) {
    			// the supplied file was not found so load it
    			loadFile(filename);	// load the file
    		}
		}
		if(colForX != Integer.MIN_VALUE){
			// select the x axis
			tree.performCheckBoxAction(dataColumnsHash.get(filename)[colForX], ItemEvent.SELECTED);
		}

		if(colForY != Integer.MIN_VALUE){
			// select the x axis
			tree.performCheckBoxAction(dataColumnsHash.get(filename)[colForY], ItemEvent.SELECTED);
		}

		updateTreeForNewData();
		updateGUI();
		tree.expandRow(dataKeys.indexOf(filename));
		twoDPanel.paintAgain();
		
	}
	public void windowActivated(WindowEvent e) {}

	public void windowClosed(WindowEvent e) {}

	@Override
	public void windowClosing(WindowEvent e) {
		if (e == null) {
			GuiManager.disposeOfParentFrame(this);
			return;
		}
		if (!promptOnClose) { // TODO add flag for dispose or hide on close
		    GuiManager.hideParentFrame(this);
		    return;
		}
		
	    String[] options;
		int choice;
		String filenames, selections = "", message;

		filenames = Array.toStr(Array.toStringArray(dataKeys), ";");

		// find the selected nodes in the plot and create a string from them delimited by ;
		String[][] selectedNodes = tree.getSelectionValues();
		for (int i = 0; i < selectedNodes.length; i++) {
		    selections += Array.toStr(selectedNodes[i]) + ";";
		}
		if (selections.endsWith(";")) { selections = selections.substring(0, selections.length() - 1); }
//		if(selectedNodes.length != 0 && selectedNodes[0][0] != null && selectedNodes[1][0] != null) {
//			selections = Array.toStr(selectedNodes[0], ",") + ";" + Array.toStr(selectedNodes[1], ",");
//		}

		options = new String[] {"Yes", "No", "Cancel"};
		message = "";
		
		if (!Array.toStr(filesProperty.getValue()).equals(filenames)) {
			if (filenames.equals("")) {
				message = "All files have been unloaded from 2D Plot.";
			} else {
				message = "A different set of files have been loaded/unloaded from 2D Plot.";
			}
		}
		
		if (!filesProperty.getValue().equals(selections)) {
			if (!message.equals("")) {
				if(!selections.equals(""))
					message = message.substring(0, message.length()-1)+", and new variables have been selected.";
				else
					message = message.substring(0, message.length()-1)+", and all variables have been unselected.";
			} else {
				if(!selections.equals(""))
					message = "New variables have been selected.";
				else
					message = "All variables have been unselected from 2D Plot.";
			}
		}
		
		if (message.equals("")) {
			GuiManager.disposeOfParentFrame(this);
		} else {
			System.out.println("message: '"+message+"'");
			choice = JOptionPane.showOptionDialog(null, message+" Would you like to keep this configuration for the next time 2D Plot is loaded?", "Preserve 2D Plot workspace?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
			if (choice == 0) {
			    filesProperty.setValue(Array.toStringArray(dataKeys));
				proj.setProperty(proj.TWOD_LOADED_VARIABLES, selections); // TODO add vars property field
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

	int valueChangedCounter = 0;
	public void valueChanged(TreeSelectionEvent e) {
//	    if (valueChangedCounter < 3) {
//	        valueChangedCounter++;
//	        return;
//	    }
//	    valueChangedCounter = 0;
		twoDPanel.setPointsGeneratable(true);
//		twoDPanel.createImage(); // calling paintAgain sets image == null, so why call 'createImage' right before doing so?
		twoDPanel.paintAgain();
	}
	
	
	public static class ScreenToCapture {
		String dataXFile;
		String dataYFile;
		String colorFile;
		int xDataIndex, yDataIndex, colorIndex;
		int xIDIndex, yIDIndex, colorIDIndex;
		float minX, minY, maxX, maxY;
		boolean hideExcluded;
		boolean isHistogram;
		boolean createColorKey;
		boolean includeColorKey;
		
		public ScreenToCapture(String[] files, int[] dataIndices, int[] idIndices, float[] displayWindow, boolean excluded, boolean colorKey, boolean appendColorKey, boolean hist) {
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
		}
	}
	
	private static Vector<String> parseControlFile(String filename, Logger log) {
		Vector<String> params;
		
		params = Files.parseControlFile(filename, "twoDScreenshots", new String[] {
				"proj=",
				"dir=",
				"# Per-line:",
				"# after being set once, all tags can be duplicated in a subsequent line by omitting them; the only necessary tags are {xDataColumn=<> and yDataColumn=<>}",
				"# fileX=<> fileY=<> xDataColumn=<> yDataColumn=<> xIDColumn=<> yIDColumn=<> colorFile=<> colorDataColumn=<> colorIDColumn=<> minX=<> minY=<> maxX=<> maxY=<> hideExcluded=True/False colorKey=True/False includeColorKey=True/False isHistogram=True/False",
		}, log);
		
		return params;
	}

	private static ArrayList<ScreenToCapture> condenseCtrlFile(java.util.List<String> ctrlLines, boolean fail) {
		HashSet<String> tagSet = new HashSet<String>();
		tagSet.add("fileX");
		tagSet.add("fileY");
		tagSet.add("colorFile");
		tagSet.add("xDataColumn");
		tagSet.add("yDataColumn");
		tagSet.add("colorDataColumn");
		tagSet.add("xIDColumn");
		tagSet.add("yIDColumn");
		tagSet.add("colorIDColumn");
		tagSet.add("minX");
		tagSet.add("minY");
		tagSet.add("maxX");
		tagSet.add("maxY");
		tagSet.add("hideExcluded");
		tagSet.add("isHistogram");
		tagSet.add("colorKey");
		tagSet.add("includeColorKey");
		
		HashMap<String, ArrayList<String>> tagValues = new HashMap<String, ArrayList<String>>();
		for (String tagKey : tagSet) tagValues.put(tagKey, new ArrayList<String>());
		
		HashSet<String> lineTagEntries = new HashSet<String>();
		for (String line : ctrlLines) {
			String[] lineTags = line.split("\\s+");
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
					        throw new IllegalArgumentException("Error - duplicate argument \"" + tagKey + "\" detected.");
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
                        System.err.println("Error - missing argument \"" + tag + "\" detected.");
                        if (fail) {
                            throw new IllegalArgumentException("Error - missing argument \"" + tag + "\" detected.");
                        }
					} else {
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
		
			files[0] = tagValues.get("fileX").get(i);
			files[1] = tagValues.get("fileY").get(i);
			files[2] = tagValues.get("colorFile").get(i);
			
			dataCols[0] = Integer.parseInt(tagValues.get("xDataColumn").get(i));
			dataCols[1] = Integer.parseInt(tagValues.get("yDataColumn").get(i));
			dataCols[2] = Integer.parseInt(tagValues.get("colorDataColumn").get(i));
			
			idCols[0] = Integer.parseInt(tagValues.get("xIDColumn").get(i));
			idCols[1] = Integer.parseInt(tagValues.get("yIDColumn").get(i));
			idCols[2] = Integer.parseInt(tagValues.get("colorIDColumn").get(i));

			window[0] = Float.parseFloat(tagValues.get("minX").get(i));
			window[1] = Float.parseFloat(tagValues.get("maxX").get(i));
			window[2] = Float.parseFloat(tagValues.get("minY").get(i));
			window[3] = Float.parseFloat(tagValues.get("maxY").get(i));
			
			hideExcludes = Boolean.parseBoolean(tagValues.get("hideExcluded").get(i));
			isHistogram = Boolean.parseBoolean(tagValues.get("isHistogram").get(i));
			colorKey = Boolean.parseBoolean(tagValues.get("colorKey").get(i));
			inclKey = Boolean.parseBoolean(tagValues.get("includeColorKey").get(i));
			
			ScreenToCapture sc = new ScreenToCapture(files, dataCols, idCols, window, hideExcludes, colorKey, inclKey, isHistogram);
			caps.add(sc);
		}
		
		return caps;
	}
	
	public void createScreenshots(String baseDir, ArrayList<ScreenToCapture> screens) {
		this.generatingScreenshots = true;
		HashSet<String> dataFiles = new HashSet<String>();

		for (ScreenToCapture cap : screens) {
			if (cap.dataXFile != null && Files.exists(baseDir + cap.dataXFile)) dataFiles.add(cap.dataXFile);
			if (cap.dataYFile != null && Files.exists(baseDir + cap.dataYFile)) dataFiles.add(cap.dataYFile);
			if (cap.colorFile != null && Files.exists(baseDir + cap.colorFile)) dataFiles.add(cap.colorFile);
		}
		
		for (String file : dataFiles) {
			loadFile(baseDir + file);
		}
		
		twoDPanel.setPointsGeneratable(true);	
		twoDPanel.paintAgain();
		updateTreeForNewData();
		
		for (ScreenToCapture screencap : screens) {			
		    this.hideExcludes = screencap.hideExcluded;
		    this.setHistogram(screencap.isHistogram);
		    
		    twoDPanel.setForcePlotXmin(screencap.minX);
		    twoDPanel.setForcePlotXmax(screencap.maxX);
		    twoDPanel.setForcePlotYmin(screencap.minY);
		    twoDPanel.setForcePlotYmax(screencap.maxY);
		    
		    boolean colorLoaded = false;
			if (screencap.colorFile != null && Files.exists(baseDir + screencap.colorFile)) {
				loadColor(baseDir, screencap);
				colorLoaded = true; // vulnerable to issues actually loading color file, but good enough for now? TODO need some way to specify HeatMap/Genotype coloration
			}
			
			if (screencap.dataXFile != null) {
			    tree.performCheckBoxAction(screencap.dataXFile, dataColumnsHash.get(baseDir + screencap.dataXFile)[screencap.xDataIndex], ItemEvent.SELECTED);
			}
			if (screencap.dataYFile != null) {
			    // 07/23/15 BRC -- this was screencap.dataXFile, which doesn't make sense (changed it to dataYFile).  If it breaks, put it back.
			    tree.performCheckBoxAction(screencap.dataYFile, dataColumnsHash.get(baseDir + screencap.dataYFile)[screencap.yDataIndex], ItemEvent.SELECTED);
			}

			twoDPanel.setChartType(AbstractPanel.SCATTER_PLOT_TYPE);
			this.colorKeyPanel.getClassRadioButtons()[colorLoaded ? (this.colorKeyPanel.getClassRadioButtons().length - 1) : 0].setSelected(true);
			
			twoDPanel.createImage();
			
			int count = 1;
			String basename = "";
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
			String screenname = basename;
			while((new File(baseDir + ext.replaceWithLinuxSafeCharacters(screenname, true) + ".png")).exists()) {
				screenname = basename + "_v" + count;
				count++;
			}
			
			screenname = baseDir + ext.replaceWithLinuxSafeCharacters(screenname, true) + ".png";
			
			if (screencap.createColorKey) {
				// Use a JFrame for it's 'pack()' method - this shrinks colorKeyPanel to the minimum required dimensions
				JFrame frame = new JFrame();
				// then change the layout (briefly) to disable WrapLayout's line-wrapping
				colorKeyPanel.classValuesPanel.setLayout(new FlowLayout(FlowLayout.CENTER, 0, 0));
				frame.add(colorKeyPanel.classValuesPanel);
				frame.pack();
				BufferedImage bi = new BufferedImage(colorKeyPanel.classValuesPanel.getWidth(), colorKeyPanel.classValuesPanel.getHeight(), BufferedImage.TYPE_INT_ARGB);
				colorKeyPanel.classValuesPanel.paint(bi.createGraphics());
				// then reset and dispose of extra resources
				colorKeyPanel.classValuesPanel.setLayout(new cnv.gui.WrapLayout(FlowLayout.CENTER, 0, 0));
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
//					g.drawImage(bi, 0, twoDPanel.getHeight(), null);
					try {
						ImageIO.write(img, "png", new File(screenname));
					} catch (IOException ie) {
						JOptionPane.showMessageDialog(null, "Error while trying to save the plot");
					}
				} else {
					twoDPanel.screenCapture(screenname);
					
					String legName = baseDir + ext.replaceWithLinuxSafeCharacters(basename, true) + "_legend.png";
					count = 1;
					while((new File(legName).exists())) {
						legName = baseDir + ext.replaceWithLinuxSafeCharacters(basename, true) + "_legend_v" + count + ".png";
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
		Hashtable<String, String> colorData = HashVec.loadFileToHashString(baseDir + screencap.colorFile, screencap.colorIDIndex, new int[]{screencap.colorIndex}, "\t", true);
		this.colorData = new HashMap<String, Integer>();
		for (java.util.Map.Entry<String, String> entry : colorData.entrySet()) {
			this.colorData.put(entry.getKey(), Integer.valueOf(entry.getValue()));
		}
	}

	private byte getColorForScreenshot(String id) {
		return (byte) (this.colorData.get(id) == null ? 0 : this.colorData.get(id).intValue());
	}

	public void addDataSource(String dataName, String[][] data, String[] columns) {
	    if (dataKeys.contains(dataName)) {
	        return;
	    }
	    if (dataName == null ||
	            columns == null ||
	            data == null) 
	        return;
	    dataKeys.add(dataName);
	    addDataHeader(dataName, columns);

//        if (sampleData != null) {
//            sampleData.initLinkKey(dataName);   // initialize the link key
//        }
        //createLinkKeyToDataHash(filename, linkKeyIndices);
        dataHash.put(dataName, new ArrayList<String[]>());
        
        for (String[] line : data) {
            if (line == null) continue;
            dataHash.get(dataName).add(line);
            validColumnsHash.put(dataName, validate(line));
        }
	    
        updateTreeForNewData();
	}
	
	public void loadFile(String filename) {
		BufferedReader reader;
		String[] header, line;
		String readBuffer;

		if (dataKeys.contains(filename)) {
			return;
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
    		    log.reportError("Error - extension of file {" + filename + "} is invalid.  Valid extensions are: {" + validExts.toString() + "}");
    		    return;
    		}
		}
		try {
		    long t1 = System.currentTimeMillis();
			reader = new BufferedReader(new FileReader(filename));
			dataKeys.add(filename);

			readBuffer = reader.readLine();
			String delim = ext.determineDelimiter(readBuffer);
			header = readBuffer.trim().split(delim);
			
			addDataHeader(filename, header);
			
//			validColumnsHash.put(filename, Array.booleanArray(dataColumnsHash.get(filename).length, true));
        	
//        	if (sampleData != null) {
//        	    sampleData.initLinkKey(filename);	// initialize the link key
//        	}
			//createLinkKeyToDataHash(filename, linkKeyIndices);
        	dataHash.put(filename, new ArrayList<String[]>());
        	
        	String tempLine = "";
        	while ((tempLine = reader.readLine()) != null) {
        		if ("".equals(tempLine)) continue;
				line = tempLine.trim().split(delim);
            	dataHash.get(filename).add(line);
            	validColumnsHash.put(filename, validate(line));
            }
            reader.close();
            
            log.reportTimeElapsed("Read file " + filename + " in: ", t1);
		
		} catch (FileNotFoundException fnfe) {
            log.reportError("Error: file \""+filename+"\" not found in current directory");
        } catch (IOException ioe) {
            log.reportError("Error reading file \""+filename+"\"");
        }
	}
	
	private void addDataHeader(String filename, String[] header) {
        int[] linkKeyIndices;
        
        dataColumnsHash.put(filename, header);
        
        linkKeyIndices = ext.indexFactors(LINKERS, header, false, true, false, log, false);
        
        if (linkKeyIndices[IID_INDEX_IN_LINKERS] == -1) {
            log.report("ID linker not automatically identified for file '" + filename + "'; assuming the first column.");
            linkKeyIndices[IID_INDEX_IN_LINKERS] = 0;
        }
        
        columnMetaData.put(filename, new HashMap<Integer, String[]>());
        
        String region, chr, start, stop, lbl;
        String[] regSplit;
        // REGION
        if (linkKeyIndices[REGION_INDEX_IN_LINKERS] == -1 || linkKeyIndices[CHR_INDEX_IN_LINKERS] == -1 || linkKeyIndices[POS_INDEX_IN_LINKERS] == -1 || linkKeyIndices[STOP_POS_INDEX_IN_LINKERS] == -1) {
            //columnMetaData
            for (int i = 0; i < header.length; i++) {
                Matcher matcher = headerPattern.matcher(header[i]);
                if (matcher.matches()) {
                    lbl = matcher.group(1);
                    chr = matcher.group(2);
                    region = matcher.group(3);
                    regSplit = region.split("-");
                    start = regSplit[0];
                    stop = regSplit[1];
                    columnMetaData.get(filename).put(i, new String[]{lbl, chr, region, start, stop});
                }
            }
        }
        linkerIndices.put(filename, linkKeyIndices);
    }

    private boolean[] validate(String[] data) {
	    boolean[] result = Array.booleanArray(data.length, true);
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
        		tree = new CheckBoxTree(namesOfBranches, branchHandles, namesOfNodes, validColumnsHash.get(key), /*isHistPlot ? 1 :*/ 2);
        	} else {
        		boolean found = false;
        		for (int j=0; j<tree.getModel().getChildCount(tree.getModel().getRoot()); j++) {
        			if (tree.getModel().getChild(tree.getModel().getRoot(), j).toString().equals(ext.removeDirectoryInfo(key))) {
        				found = true;
        			}
        		}
        		if (!found) {
        			tree.addNode(ext.removeDirectoryInfo(key), key, treeFileVariableNameLookup[i], validColumnsHash.get(key), checkBoxMouseListener);
        		}
        	}
       }
	}

	/**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
	 * @param proj Project
	 * @param show setVisible
	 * @param promptOnClose prompt to save files/vars on close
	 * @param fileExts allowed file extensions for files (null for all files)
	 * @return
	 */
	public static TwoDPlot createGUI(Project proj, boolean show, boolean promptOnClose, StringListProperty filenamesProperty, String... fileExts) {
		JFrame frame = new JFrame("Genvisis - 2D Plot - " + proj.PROJECT_NAME.getValue());
		frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);

        //Create and set up the content pane.
        TwoDPlot twoDPlot = new TwoDPlot(proj, promptOnClose, fileExts, filenamesProperty);
        frame.setJMenuBar(twoDPlot.menuBar());
        twoDPlot.setOpaque(true); //content panes must be opaque
        frame.setContentPane(twoDPlot);
        frame.addWindowListener(twoDPlot);
        frame.setBounds(20, 20, 1000, 600);
//		frame.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH);

        //Display the window.
        frame.pack();
        frame.setVisible(show);

		return twoDPlot;
    }
	
	/*
	 * CHR POS START STOP LABEL
	 */
	public String[][] getCurrentColumnMetaData() {
	    ArrayList<String[]> resultList = new ArrayList<String[]>();
        HashSet<String> inclSet = new HashSet<String>();
        String[][] selectionValues = tree.getSelectionValues();
        for (int i = 0; i < selectionValues.length; i++) {
            String[] result = null;
            if (selectionValues[i] != null && selectionValues[i][0] != null) {
                String file = selectionValues[i][0];
                int col = Integer.parseInt(selectionValues[i][1]);
                String[] metaData = getColumnMetaData(file, col);
                if (metaData != null) {
                    result = new String[]{metaData[1], metaData[2], metaData[3], metaData[4], metaData[0]};
                } else {
                    String colName = dataColumnsHash.get(file)[col]; 
                    Matcher m = headerPattern.matcher(colName);
                    if (m.matches()) {
                        result = new String[]{m.group(2), m.group(3), m.group(3).split("-")[0], m.group(3).split("-")[1], m.group(1)};
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
	
	public String[] getColumnMetaData(String filename, int index) {
		return columnMetaData.get(filename).get(index);
	}
	
	public static void fromParameters(String filename, Logger log) {
		Vector<String> params = parseControlFile(filename, log);
		
		if (params != null) {
			final String projFile = params.get(0).split("=")[1];
			final String baseDir = params.get(1).split("=")[1];
			final ArrayList<ScreenToCapture> screens = condenseCtrlFile(params.subList(2, params.size()), true);
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
	            public void run() {
	                Project proj = new Project(projFile, false);
	                TwoDPlot tdp = createGUI(proj, false, false, proj.TWOD_LOADED_FILENAMES);
	                tdp.createScreenshots(baseDir, screens);
	                tdp.windowClosing(null);
	            }
	        });
		}
	}

	public static void main(String[] args) {
//		HashSet<String> tagSet = new HashSet<String>();
//		tagSet.add("fileXs");
//		tagSet.add("fileYs");
//		tagSet.add("xDataColumns");
//		tagSet.add("yDataColumns");
//		tagSet.add("xIDColumns");
//		tagSet.add("yIDColumns");
//		tagSet.add("colorFiles");
//		tagSet.add("colorDataColumns");
//		tagSet.add("colorIDColumns");
//		tagSet.add("excludeFiles");
//		tagSet.add("excludeDataColumns");
//		tagSet.add("excludeIDColumns");
//		tagSet.add("minXs");
//		tagSet.add("minYs");
//		tagSet.add("maxXs");
//		tagSet.add("maxYs");
//		
//		HashMap<String, String> tagValues = new HashMap<String, String>();
//		
//		String usage = "";
//		
//		System.out.println();
//		for (int i = 0; i<args.length; i++) {
//			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
//				System.err.println(usage);
//				System.exit(1);
//			} else if (tagSet.contains(args[i].split("=")[0])) {
//				tagValues.put(args[i].split("=")[0], args[i].split("=")[1]);
//			}
//		}
//		
//		
//		
//		fromParameters("D:/data/FarrarReparse/classification/twoDscreenshots.dat", new Logger());
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createGUI(new Project(), true, false, null, null);
//                createAndShowGUI(new Project(cnv.Launch.getDefaultDebugProjectFile(true), false));
            }
        });
	}
}


