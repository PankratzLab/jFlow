package cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.LineNumberReader;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.TreeMap;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLayeredPane;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;

import cnv.filesys.Project;
import common.Aliases;
import common.Array;
import common.Files;
import common.Grafik;
import common.Logger;
import common.ext;

class MetaStudy {
	private final ArrayList<StudyData> studies;
	private ArrayList<StudyData> sorted;
	private final HashMap<String, StudyData> nameMap;
	private final float metaBeta;
	private final float metaStderr;
	private final float[] metaConf = new float[2];
    private ArrayList<String> sortOrder = null;
    private boolean shouldSort;
    private boolean currentSortIsNaturalSort = false;
	
	public MetaStudy(float metaBeta, float metaStderr) {
		studies = new ArrayList<StudyData>();
		nameMap = new HashMap<String, StudyData>();
		this.metaBeta = metaBeta;
		this.metaStderr = metaStderr;
		this.metaConf[0] = (float) (metaBeta - 1.96 * metaStderr);
		this.metaConf[1] = (float) (metaBeta + 1.96 * metaStderr);
	}
	
	public void addStudy(StudyData studyData) {
		studies.add(studyData);
		nameMap.put(studyData.getLabel(), studyData);
	}

	String findLongestStudyName() {
		String longest = "";
		for(StudyData ft : studies){
			longest = longest.length() < ft.getDisplayLabel().length() ? ft.getDisplayLabel() : longest;
		}
		return longest;
	}
	
	float calcSumZScore() {
		float sum = 0;
		for	(StudyData ft : studies){
			sum += ft.getZScore();
		}
		return sum;
	}
	
	float findMaxZScore() {
		float max = Float.MIN_VALUE;
		for (StudyData data: studies) {
			max = Math.max(max, data.getZScore());
		}
		return max;
	}
	
	private ArrayList<StudyData> getSorted(ArrayList<String> order) {
	    if (this.sorted == null || this.sorted.isEmpty() || currentSortIsNaturalSort) {
	        this.sorted = new ArrayList<StudyData>();
	        for (int i = order.size() - 1; i >= 0; i--) {
	            String name = ext.replaceAllWith(order.get(i), ForestPlot.REPLACEMENTS_FOOLISHLY_HARD_CODED);
	            String repl = null;
	            if (name.equals("")) {
	                this.sorted.add(new StudyBreak());
	            } else {
	                if (name.split("\t").length > 1) {
	                    repl = name.split("\t")[1];
	                    name = name.split("\t")[0];
	                }
	                StudyData sd = nameMap.get(name);
	                if (sd == null) {
	                    sd = new StudyBreak();
	                } else {
    	                if (repl != null) {
    	                    sd.setReplacementLabel(repl);
    	                }
	                }
	                sorted.add(sd);
	            }
	        }
	    }
        currentSortIsNaturalSort = false;
	    return this.sorted;
	}
	
	private ArrayList<StudyData> getSorted() {
		if (this.sorted == null || !currentSortIsNaturalSort) {
			this.sorted = new ArrayList<StudyData>();
			
			TreeMap<String, String> zeroStudyMap = new TreeMap<String, String>();
			TreeMap<Float, String> betaStudyMap = new TreeMap<Float, String>();
			for (StudyData study : studies) {
				if (study.getBeta(false) == 0.0f) {
					zeroStudyMap.put(study.getLabel(), study.getLabel());
				} else {
					betaStudyMap.put(study.getBeta(false), study.getLabel());
				}
			}
			ArrayList<StudyData> desc = new ArrayList<StudyData>();
			for (java.util.Map.Entry<String, String> entry : zeroStudyMap.entrySet()) {
				desc.add(nameMap.get(entry.getValue()));
			}
			for (java.util.Map.Entry<Float, String> entry : betaStudyMap.entrySet()) {
				desc.add(nameMap.get(entry.getValue()));
			}
			for (int i = desc.size() - 1; i >= 0; i--) {
				sorted.add(desc.get(i));
			}
		}
		currentSortIsNaturalSort = true;
		return this.sorted;
	}

	public ArrayList<StudyData> getStudies() {
	    return this.shouldSort ? (this.sortOrder == null || this.sortOrder.isEmpty() ? getSorted() : getSorted(this.sortOrder)) : this.studies;
	}

	public float[] getMetaConf(boolean odds) {
		return odds ? new float[]{(float) Math.exp(metaConf[0]), (float) Math.exp(metaConf[1])} : metaConf;
	}

//	public float getMetaBeta() {
//		return metaBeta;
//	}
//
//	public float getMetaStderr() {
//	    return metaStderr;
//	}
	public float getMetaBeta(boolean odds) {
	    return (float) (odds ? Math.exp(metaBeta) : metaBeta);
	}
	
	public float getMetaStderr(boolean odds) {
	    return (float) (odds ? Math.exp(metaStderr) : metaStderr);
	}

    public void setSort(boolean sortedDisplay, ArrayList<String> sortOrder) {
        this.shouldSort = sortedDisplay;
        this.sortOrder = sortOrder;
        this.sorted = null;
    }
	
}

class StudyBreak extends StudyData {
    // placeholder class for visual breaks
    public StudyBreak() {
        this ("", 0f, 0f, 0, (byte) 0);
    }
    private StudyBreak(String label, float beta, float stderr, int color, byte shape) {
        super(label, beta, stderr, color, shape);
    }
}

class StudyData {
	private final String label;
	private String replLabel = null;
	private final float beta;
	private final float stderr;
	private int color;
	private byte shape;
	private final float[] confInterval;
	private final float zScore;

	public StudyData(String label, float beta, float stderr, int color, byte shape) {
		this.label = label;
		this.beta = beta;
		this.stderr = stderr;
		this.color = color;
		this.shape = shape;
		this.confInterval = new float[2];
		this.confInterval[0] = (float) (beta - 1.96 * stderr);
		this.confInterval[1] = (float) (beta + 1.96 * stderr);
		this.zScore = stderr == 0.0f ? 0.0f : Math.abs(beta / stderr);
	}

	public void setReplacementLabel(String repl) {
       this.replLabel = repl;
    }
	
	public String getDisplayLabel() {
	    if (replLabel != null) {
	        return replLabel;
	    }
	    return label;
	}
	
    public String getLabel() {
		return label;
	}

	public float getBeta(boolean odds) {
		return (float) (odds ? Math.exp(beta) : beta);
	}

	public float getStderr(boolean odds) {
		return (float) (odds ? Math.exp(stderr) : stderr);
	}

	public int getColor() {
		return color;
	}

	public byte getShape() {
		return shape;
	}

	public float[] getConfInterval(boolean odds) {
		return odds ? new float[]{(float) Math.exp(confInterval[0]), (float) Math.exp(confInterval[1])} : confInterval;
	}

	public float getZScore() {
		return zScore;
	}
}

class ForestInput {
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((comment == null) ? 0 : comment.hashCode());
		result = prime * result + ((file == null) ? 0 : file.hashCode());
		result = prime * result + ((marker == null) ? 0 : marker.hashCode());
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ForestInput other = (ForestInput) obj;
		if (comment == null) {
			if (other.comment != null)
				return false;
		} else if (!comment.equals(other.comment))
			return false;
		if (file == null) {
			if (other.file != null)
				return false;
		} else if (!file.equals(other.file))
			return false;
		if (marker == null) {
			if (other.marker != null)
				return false;
		} else if (!marker.equals(other.marker))
			return false;
		return true;
	}
	
	final String marker;
	final String file;
	final String comment;
	int[] metaIndicies;
	HashMap<String, Integer> studyToColIndexMap;
	ArrayList<String> studyList;
	MetaStudy ms;
	
	public ForestInput(String marker, String file, String comment) {
		this.marker = marker;
		this.file = file;
		this.comment = comment;
		this.studyToColIndexMap = new HashMap<String, Integer>();
		this.studyList = new ArrayList<String>();
	}

	public void addStudy(String string, int i) {
		studyList.add(string);
		studyToColIndexMap.put(string, i);
	}
	public MetaStudy getMetaStudy() { return ms; }
	public void setMetaStudy(MetaStudy ms) { this.ms = ms; }
}

public class ForestPlot extends JFrame implements WindowListener {
	private static final long serialVersionUID = 1L;
	
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String ADD_DATA_FILE = "Add Data File";
	public static final String REMOVE_DATA_FILE = "Remove Data File";
//	private static final String ALT_UP = "ALT UP";
//	private static final String ALT_DOWN = "ALT DOWN";
//	private static final String ALT_LEFT = "ALT LEFT";
//	private static final String ALT_RIGHT = "ALT RIGHT";
	private static final String[] BETA_META_HEADERS = {"beta", "effect"};
	private static final String[] SE_META_HEADERS = {"se", "stderr"};
	private static final String BETA_PREFIX = "beta.";
	private static final String SE_PREFIX = "se.";
	private static final String FIRST = "First";
	private static final String PREVIOUS = "Previous";
	private static final String NEXT = "Next";
	private static final String LAST = "Last";
	private static final String MARKER_LIST_NEW_FILE = "Add New File(s)...";
	private static final String MARKER_LIST_PLACEHOLDER = "Select Input File...";
//	private LinkedHashSet<ForestInput> data = new LinkedHashSet<ForestInput>();
	private ArrayList<ForestInput> dataIndices = new ArrayList<ForestInput>();
//	private HashMap<ForestInput, MetaStudy> dataToMetaMap = new HashMap<ForestInput, MetaStudy>();
	private MetaStudy currMetaStudy;
	private String plotLabel;
	private Logger log;
	private ForestPanel forestPanel;
	private float maxZScore;
	private float sumZScore;
	private String longestStudyName;
//	private JCheckBox chkSortStudies;
	private JLayeredPane layeredPane;
	private JButton first, previous, next, last;
//	private JButton btnScreen, btnScreenAll;
	private JTextField navigationField;
	private JComboBox<String> markerFileList;
	private HashMap<String, String> markerFileNameLoc = new HashMap<String, String>();
//	int curMarkerIndex;
	private int currentDataIndex;
	private boolean atleastOneStudy;
    private boolean sortedDisplay = false;

	private String markerFileName;

	protected Project proj;

	private volatile boolean loadingFile;
	private Thread loadingThread;
	public static final  String[][] REPLACEMENTS_FOOLISHLY_HARD_CODED = new String[][] {
			{"_WBC_TOTAL", ""},
			{"_WBC_NEUTRO", ""},
			{"_", " "},
	};
	
	public ForestPlot(Project proj) {
		super("Genvisis - Forest Plot - " + proj.PROJECT_NAME.getValue());
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		this.proj = proj;
		this.log = proj.getLog();
		setup();
		setBounds(20, 20, 1000, 600);
		pack();
		setVisible(true);
		this.addWindowListener(this);
	}
	
	public ForestPlot(String markerFile, Logger log) {
		super("Genvisis - Forest Plot");
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		this.log = log;
		this.markerFileName = markerFile;
		setup();
		setBounds(20, 20, 1000, 600);
		pack();
		setVisible(true);
		this.addWindowListener(this);
	}

	private void setup() {
		setLayout(new BorderLayout());

		forestPanel = new ForestPanel(this, log);
		forestPanel.setLayout(new BorderLayout());

		layeredPane = new JLayeredPane();
		layeredPane.setLayout(new BorderLayout());

		layeredPane.add(forestPanel);
		layeredPane.setPreferredSize(new Dimension(1000, 600));

		JPanel treePanel = new JPanel();
		treePanel.setBackground(BACKGROUND_COLOR);
		treePanel.setLayout(new BorderLayout());
		
		forestPanel.add(createControlPanel(), BorderLayout.NORTH);
		
		Dimension minimumSize = new Dimension(100, 50);
		layeredPane.setMinimumSize(minimumSize);

		add(layeredPane, BorderLayout.CENTER);
		inputMapAndActionMap();

		forestPanel.grabFocus();

		progressBar = new JProgressBar();
		progressBar.setPreferredSize(new Dimension(progressBar.getPreferredSize().width, 22));
		progressBar.setMinimumSize(new Dimension(progressBar.getPreferredSize().width, 22));
		
		JButton btnCancelProg = new JButton(Grafik.getImageIcon("images/delete2sm.png", true));
		btnCancelProg.setMinimumSize(new Dimension(30, 20));
		btnCancelProg.setPreferredSize(new Dimension(30, 20));
		btnCancelProg.setMaximumSize(new Dimension(30, 20));
		btnCancelProg.setAlignmentX(Component.RIGHT_ALIGNMENT);
		btnCancelProg.addActionListener(new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent e) {
				ForestPlot.this.loadingThread.interrupt();
			}
		});
		
		progressBar.setLayout(new BorderLayout());
		progressBar.add(btnCancelProg, BorderLayout.EAST);
		
		add(progressBar, BorderLayout.SOUTH);
		
		setJMenuBar(createMenuBar());
		
		setVisible(true);
		
		loadMarkerFile();
		// generateShortcutMenus();
	}
	
	private JMenuBar createMenuBar() {
		JMenuBar bar = new JMenuBar();
		
		JMenu actions = new JMenu("Actions", true);
		actions.setMnemonic(KeyEvent.VK_A);
		
		final JCheckBoxMenuItem oddsRatioButton = new JCheckBoxMenuItem();
		AbstractAction oddsRatioAction = new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                ForestPlot.this.setOddsRatioDisplay(oddsRatioButton.isSelected());
                ForestPlot.this.forestPanel.paintAgain();                
            }
        };
		oddsRatioButton.setAction(oddsRatioAction);
		oddsRatioButton.setText("Display as Odds Ratios");
		oddsRatioButton.setMnemonic(KeyEvent.VK_O);
		actions.add(oddsRatioButton);
		
		actions.add(new JSeparator());
        
		JMenuItem createSortFile = new JMenuItem();
		createSortFile.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                generateStudyNameFile();
            }
		});
		createSortFile.setText("Create Sort File");
		actions.add(createSortFile);
		
		final JRadioButtonMenuItem noSortButton = new JRadioButtonMenuItem();
		AbstractAction noSortAction = new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                ForestPlot.this.setSortedDisplay(false);
                ForestPlot.this.forestPanel.paintAgain();                
            }
        };
        noSortButton.setAction(noSortAction);
        noSortButton.setText("No Sorting");
        noSortButton.setMnemonic(KeyEvent.VK_N);
        actions.add(noSortButton);
		
		final JRadioButtonMenuItem sortStudies = new JRadioButtonMenuItem();
		AbstractAction sortAction = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent arg0) {
				ForestPlot.this.setSortedDisplay(true);
				ForestPlot.this.forestPanel.paintAgain();
			}
		};
		sortStudies.setAction(sortAction);
		sortStudies.setText("Sort Studies");
		sortStudies.setMnemonic(KeyEvent.VK_S);
		actions.add(sortStudies);
		
		final JRadioButtonMenuItem sortFileStudies = new JRadioButtonMenuItem();
		AbstractAction sortFileAction = new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                
                JFileChooser jfc = new JFileChooser(new File(proj == null ? ext.parseDirectoryOfFile(markerFileName) : proj.PROJECT_DIRECTORY.getValue()));
                int returnVal = jfc.showOpenDialog(null);
                
                if (returnVal == JFileChooser.APPROVE_OPTION) {
                    loadOrderFile(jfc.getSelectedFile().getAbsolutePath(), true);
                    ForestPlot.this.forestPanel.paintAgain();
                } else {
                    noSortButton.setSelected(true);
                    ForestPlot.this.setSortedDisplay(false);
                    ForestPlot.this.forestPanel.paintAgain();
                }
            }
        };
    	sortFileStudies.setAction(sortFileAction);
    	sortFileStudies.setText("Sort Studies with File");
    	sortFileStudies.setMnemonic(KeyEvent.VK_S);
    	actions.add(sortFileStudies);
        
        ButtonGroup sortGroup = new ButtonGroup();
        sortGroup.add(noSortButton);
        sortGroup.add(sortStudies);
        sortGroup.add(sortFileStudies);
        
        noSortButton.setSelected(true);
		
		AbstractAction screenAction1 = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent e) {
				if (!loadingFile && (markerFileName != null && !"".equals(markerFileName))) {
					ForestPlot.this.screenCap();
				}
			}
		};
		AbstractAction screenAction2 = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent e) {
				if (!loadingFile && (markerFileName != null && !"".equals(markerFileName))) {
					ForestPlot.this.screenCapAll();
				}
			}
		};

		actions.add(new JSeparator());
		
		JMenuItem screenCap1 = new JMenuItem();
		screenCap1.setAction(screenAction1);
		screenCap1.setText("Screen Capture");
		screenCap1.setMnemonic(KeyEvent.VK_C);
		actions.add(screenCap1);
		
		JMenuItem screenCap2 = new JMenuItem();
		screenCap2.setAction(screenAction2);
		screenCap2.setText("Screen Capture All");
		screenCap2.setMnemonic(KeyEvent.VK_A);
		actions.add(screenCap2);
		
		bar.add(actions);
		
		return bar;
	}

	public void setOddsRatioDisplay(boolean selected) {
        this.forestPanel.oddsDisplay = selected;
    }

    private JPanel createControlPanel() {
		first = new JButton(Grafik.getImageIcon("images/firstLast/First.gif", true));
		first.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif", true));
		first.addActionListener(navFirstAction);
		first.setActionCommand(FIRST);
		first.setPreferredSize(new Dimension(20, 20));
		previous = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif", true));
		previous.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		previous.addActionListener(navPrevAction);
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
					if (trav >=0 && trav < getDataIndices().size()) {
						setCurrentData(trav);
						updateForestPlot();
					}
				} catch (NumberFormatException nfe) {
					System.out.println("Please enter a valid integer which is in range");
				}
				//displayIndex((JTextField) e.getSource());
				forestPanel.setPointsGeneratable(true);
				updateForestPlot();
			}
		});
	
		next = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif", true));
		next.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		next.addActionListener(navNextAction);
		next.setActionCommand(NEXT);
		next.setPreferredSize(new Dimension(20, 20));
		last = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif", true));
		last.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif", true));
		last.addActionListener(navLastAction);
		last.setActionCommand(LAST);
		last.setPreferredSize(new Dimension(20, 20));
		
		JPanel subNavPanel = new JPanel();
		subNavPanel.setBackground(BACKGROUND_COLOR);
		subNavPanel.add(first);
		subNavPanel.add(previous);
		subNavPanel.add(navigationField);
		subNavPanel.add(next);
		subNavPanel.add(last);
		
		Vector<String> items = new Vector<String>();
		items.add(MARKER_LIST_PLACEHOLDER);
		if (proj != null) {
			String[] files = proj.FOREST_PLOT_FILENAMES.getValue();
			String name;
			for (String file : files) {
				name = ext.rootOf(file);
				markerFileNameLoc.put(name, file);
				items.add(name);
			}
		}
		items.add(MARKER_LIST_NEW_FILE);
		markerFileList = new JComboBox<String>(items);
		markerFileList.setFont(new Font("Arial", 0, 12));
		markerFileList.setMinimumSize(new Dimension(200, 20));
		markerFileList.setPreferredSize(new Dimension(200, 20));
		markerFileList.setMaximumSize(new Dimension(200, 20));
		AbstractAction markerFileSelectAction = new AbstractAction() {
			private static final long serialVersionUID = 1L;
	
			@SuppressWarnings("unchecked")
			@Override
			public void actionPerformed(ActionEvent e) {
				String shortName = (String) ((JComboBox<String>)e.getSource()).getSelectedItem();
				if (!loadingFile && !MARKER_LIST_NEW_FILE.equals(shortName) && !MARKER_LIST_PLACEHOLDER.equals(shortName)) {
					String file = markerFileNameLoc.get(shortName);
					if (file != null && file.equals(ForestPlot.this.markerFileName)) {
						return;
					}
					ForestPlot.this.markerFileName = file;
					loadMarkerFile();
				} else if (loadingFile && MARKER_LIST_PLACEHOLDER.equals(shortName)) {
					// do nothing
				} else if (loadingFile || MARKER_LIST_PLACEHOLDER.equals(shortName)) {
					// leave as currently selected marker
					if (ForestPlot.this.markerFileName != "" && ForestPlot.this.markerFileName != null) {
						((JComboBox<String>)e.getSource()).setSelectedItem(ext.rootOf(ForestPlot.this.markerFileName));
					}
					return;
				} else if (MARKER_LIST_NEW_FILE.equals(shortName)) {
					chooseNewFiles();
					if (ForestPlot.this.markerFileName != null && !"".equals(ForestPlot.this.markerFileName)) {
						((JComboBox<String>)e.getSource()).setSelectedItem(ext.rootOf(ForestPlot.this.markerFileName));
					}
				}
			}
		};
		markerFileList.addActionListener(markerFileSelectAction);
		
		JPanel navigationPanel = new JPanel();
		navigationPanel.setLayout(new BoxLayout(navigationPanel, BoxLayout.Y_AXIS));
		navigationPanel.setBackground(BACKGROUND_COLOR);
		navigationPanel.add(subNavPanel);
		navigationPanel.add(markerFileList);
		navigationPanel.add(Box.createVerticalGlue());
		
//		final JPanel subOptPanel = new JPanel();
//		subOptPanel.setLayout(new BoxLayout(subOptPanel, BoxLayout.Y_AXIS));
//		subOptPanel.setBackground(BACKGROUND_COLOR);
//		subOptPanel.setAlignmentX(Component.RIGHT_ALIGNMENT);
//		subOptPanel.add(btnScreen);
//		subOptPanel.add(Box.createRigidArea(new Dimension(0,5)));
//		subOptPanel.add(btnScreenAll);
		
//		JPanel optPanel = new JPanel() {
//			private static final long serialVersionUID = 1L;
//			@Override
//			public Dimension getMaximumSize() {
//				return super.getPreferredSize();
//			}
//		};
//		optPanel.setBackground(BACKGROUND_COLOR);
//		optPanel.add(chkSortStudies);
//		optPanel.add(subOptPanel);
	
		
		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new BoxLayout(descrPanel, BoxLayout.X_AXIS));
		// All of these are needed (or at least, have some effect)
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
//		descrPanel.add(Box.createHorizontalGlue());
//		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(navigationPanel);
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
//		descrPanel.add(optPanel);
		descrPanel.setBackground(BACKGROUND_COLOR);
		return descrPanel;
	}

    public boolean waitForLoad() {
        while(this.loadingFile) {
            Thread.yield();
        }
        return true;
    }
    
	private void loadMarkerFile() {
		if (!this.loadingFile) {
			this.loadingFile = true;
			this.loadingThread = new Thread(new Runnable() {
				@Override
				public void run() {
					reloadData();
					forestPanel.setPointsGeneratable(ForestPlot.this.markerFileName != null);
					forestPanel.setRectangleGeneratable(ForestPlot.this.markerFileName != null);
					forestPanel.setExtraLayersVisible(new byte[] { 99 });
					try {
						SwingUtilities.invokeAndWait(new Runnable() {
							@Override
							public void run() {
								updateGUI();
								displayIndex(navigationField);
								progressBar.setStringPainted(false);
								progressBar.setValue(0);
								progressBar.setString(null);
								progressBar.setIndeterminate(false);
								progressBar.repaint();
							}
						});
					} catch (InvocationTargetException e) {
						// TODO Auto-generated catch block
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
					}
					loadingFile = false;
				}
			}, "ForestPlot_loadMarkerFile");
			this.loadingThread.start();
		}
	}

	private void reloadData() {
		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				@Override
				public void run() {
					progressBar.setStringPainted(true);
					progressBar.setString("Loading marker list...");
					progressBar.setValue(0);
					progressBar.setIndeterminate(true);
					progressBar.repaint();
				}
			});
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
		}
		atleastOneStudy = false;
		
		this.setDataIndices(new ArrayList<ForestInput>());
		if (this.markerFileName != null) {
			this.getDataIndices().addAll(readMarkerFile(this.markerFileName));
		} 
		
//		dataToMetaMap = new HashMap<ForestInput, MetaStudy>();

		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				@Override
				public void run() {
					progressBar.setString("Calculating datafile size(s)...");
					progressBar.repaint();
				}
			});
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
		}
		
		if (!this.getDataIndices().isEmpty() && !Thread.interrupted()) {
			loadStudyData();
			setCurrentData(0);
		} else {
			clearCurrentData();
		}
		
	}

	private void generateStudyNameFile() {
        String fileName = ext.verifyDirFormat(ext.parseDirectoryOfFile(markerFileName)) + "sort.txt";
        ArrayList<StudyData> currentData = getCurrentMetaStudy().getStudies();
        ArrayList<String> names = new ArrayList<String>();
        for (StudyData sd : currentData) {
            names.add(sd.getLabel());
        }
        Files.writeArrayList(names, fileName);
    }

    private void updateForestPlot(){
		forestPanel.setPointsGeneratable(true);
		forestPanel.setRectangleGeneratable(true);
		forestPanel.setExtraLayersVisible(new byte[] { 99 });
		displayIndex(navigationField);
		updateGUI();
	}
	
	private void chooseNewFiles() {
		JFileChooser jfc = new JFileChooser((proj != null ? proj.PROJECT_DIRECTORY.getValue() : ext.parseDirectoryOfFile(markerFileName)));
		jfc.setMultiSelectionEnabled(true);
		if (jfc.showOpenDialog(ForestPlot.this) == JFileChooser.APPROVE_OPTION) {
			File[] files = jfc.getSelectedFiles();
			if (files.length > 0) {
				boolean[] keep = Array.booleanArray(files.length, true);
				for (int i = 0; i < files.length; i++) {
					for (String fileName : markerFileNameLoc.keySet()) {
						if (ext.rootOf(files[i].toString()).equals(fileName)) {
							keep[i] = false;
						}
					}
				}
				File[] keptFiles = Array.subArray(files, keep);
				File[] discards = Array.subArray(files, Array.booleanNegative(keep));
				
				if (discards.length > 0) {
					StringBuilder msg = new StringBuilder("The following data file(s) are already present:");
					for (File disc : discards) {
						msg.append("\n").append(disc.getName());
					}
					JOptionPane.showMessageDialog(ForestPlot.this, msg.toString()); 
				}
				
				for (File kept : keptFiles) {
					addFileToList(kept.getAbsolutePath());
				}
			} else {
				File file = jfc.getSelectedFile();
				boolean keep = true;
				for (String fileName : markerFileNameLoc.keySet()) {
					if (ext.rootOf(file.toString()).equals(fileName)) {
						keep = false;
					}
				}
				
				if (!keep) {
					StringBuilder msg = new StringBuilder("The following data file is already present:\n").append(file.getName());
					JOptionPane.showMessageDialog(ForestPlot.this, msg.toString()); 
				} else {
					addFileToList(file.getAbsolutePath());
				}
				
			}
			
		}
	}
	
	private void addFileToList(String rawfile) {
		String file = ext.verifyDirFormat(rawfile);
		file = file.substring(0, file.length() - 1);
		int num = markerFileList.getModel().getSize() - 2;
		Vector<String> currFiles = new Vector<String>();
		if (num > 0) {
			for (int i = 1; i < num + 1; i++) {
				currFiles.add(markerFileList.getModel().getElementAt(i));
			}
		}
		String name = ext.rootOf(file);
		markerFileNameLoc.put(name, file);
		currFiles.add(name);
		currFiles.add(0, MARKER_LIST_PLACEHOLDER);
		currFiles.add(MARKER_LIST_NEW_FILE);
		
		markerFileList.setModel(new DefaultComboBoxModel<String>(currFiles));
	}
	
	/**
	 * NOTE: be sure to call 'waitForLoad' before calling screenCap, otherwise data might not be loaded in time.
	 * @param subdir
	 * @param odds
	 */
	public void screenCapAll(String subdir, boolean odds, boolean versionIfExists) {
	    setOddsRatioDisplay(odds);
        int currentSelection = getCurrentDataIndex();
        ArrayList<ForestInput> data = getDataIndices();
        for (int i = 0; i < data.size(); i++) {
            setCurrentData(i);
            forestPanel.createImage();
            String marker = "", filename = "", dataFile = "";
            int count = 1;
            String root = (proj == null ? ext.parseDirectoryOfFile(markerFileName) : proj.PROJECT_DIRECTORY.getValue());
            root = ext.verifyDirFormat(root);
            if (subdir != null && !subdir.equals("")) {
                root += subdir;
                root = ext.verifyDirFormat(root);
            }
            marker = getDataIndices().get(getCurrentDataIndex()).marker;
            dataFile = ext.rootOf(getDataIndices().get(getCurrentDataIndex()).file, true);
            filename = marker + "_" + dataFile;
            filename = ext.replaceWithLinuxSafeCharacters(filename, true);
            if (new File(root + filename + ".png").exists()) {
                if (versionIfExists) {
                    while (new File(root+filename+".png").exists()) {
                        filename = marker + "_" + dataFile + "_v" + count;
                        filename = ext.replaceWithLinuxSafeCharacters(filename, true);
                        count++;
                    }
                }
            }
            if (log != null) {
                log.report("Writing screenshot to file " + root + filename + ".png");
            } else {
                System.out.println("Writing screenshot to file " + root + filename + ".png");
            }
            ForestPlot.this.forestPanel.screenCapture(root+filename+".png");
        }
        setCurrentData(currentSelection);
        updateForestPlot();
	}
	
	private void screenCapAll() {
		int currentSelection = getCurrentDataIndex();
		ArrayList<ForestInput> data = getDataIndices();
		for (int i = 0; i < data.size(); i++) {
			setCurrentData(i);
			forestPanel.createImage();
			String marker = "", filename = "", ctrlFile = "";
			int count = 1;
			String root = (proj == null ? ext.parseDirectoryOfFile(markerFileName) : proj.PROJECT_DIRECTORY.getValue());
			marker = getDataIndices().get(getCurrentDataIndex()).marker;
			ctrlFile = ext.rootOf(markerFileName);
			filename = marker + "_" + ctrlFile;
			filename = ext.replaceWithLinuxSafeCharacters(filename, true);
			while (new File(root+filename+".png").exists()) {
				filename = marker + "_" + ctrlFile + "_v" + count;
				filename = ext.replaceWithLinuxSafeCharacters(filename, true);
				count++;
			}
            if (log != null) {
                log.report("Writing screenshot to file " + root + filename + ".png");
            } else {
                System.out.println("Writing screenshot to file " + root + filename + ".png");
            }
			ForestPlot.this.forestPanel.screenCapture(root+filename+".png");
		}
		setCurrentData(currentSelection);
		updateForestPlot();
	}
	
	private void screenCap() {
		String marker = "", filename = "", ctrlFile = "";
		int count = 1;
		String root = (proj == null ? ext.parseDirectoryOfFile(markerFileName) : proj.PROJECT_DIRECTORY.getValue());
		marker = getDataIndices().get(getCurrentDataIndex()).marker;
		ctrlFile = ext.rootOf(markerFileName);
		filename = marker + "_" + ctrlFile;
		filename = ext.replaceWithLinuxSafeCharacters(filename, true);
		while (new File(root+filename+".png").exists()) {
			filename = marker + "_" + ctrlFile + "_v" + count;
			filename = ext.replaceWithLinuxSafeCharacters(filename, true);
			count++;
		}
        if (log != null) {
            log.report("Writing screenshot to file " + root + filename + ".png");
        } else {
            System.out.println("Writing screenshot to file " + root + filename + ".png");
        }
		ForestPlot.this.forestPanel.screenCapture(root+filename+".png");
	}
	
	private void displayIndex(JTextField field) {
		if (getDataIndices().size() == 0) {
			field.setText("0 of 0");
		} else {
			field.setText((getCurrentDataIndex() + 1) + " of " + getDataIndices().size());
		}
	}
	
	private void setCurrentData(int index) {
		if (this.getDataIndices().size() == 0 || index < 0 || index > this.getDataIndices().size()) return;
		setCurrentDataIndex(index);
//		setCurrentMetaStudy(dataToMetaMap.get(getDataIndices().get(index)));
		setCurrentMetaStudy(getDataIndices().get(index).getMetaStudy());
		getCurrentMetaStudy().setSort(isSortedDisplay(), getSortOrder());
		if (getCurrentMetaStudy() == null) {
		    String msg = "Error - could not set index to "+index+" since the data did not load properly; check to see if any results files are missing";
		    if (log != null) {
		        log.reportError(msg);
		    } else {
		        System.err.println(msg);
		    }
			return;
		}
		maxZScore = getCurrentMetaStudy().findMaxZScore();
		maxZScore = getCurrentMetaStudy().findMaxZScore();
		sumZScore = getCurrentMetaStudy().calcSumZScore();
		longestStudyName = getCurrentMetaStudy().findLongestStudyName();
		setPlotLabel(getDataIndices().get(index).marker);
	}
	
	private void clearCurrentData() {
		setCurrentDataIndex(-1);
		setCurrentMetaStudy(null);
		maxZScore = 0;
		maxZScore = 0;
		sumZScore = 0;
		longestStudyName = "";
		setPlotLabel("");
	}
	
	public String getLongestStudyName() {
		return longestStudyName;
	}

	private void loadStudyData() throws RuntimeException {
		HashMap<String, ArrayList<ForestInput>> files = new HashMap<String, ArrayList<ForestInput>>();
		for (ForestInput fi : getDataIndices()) {
			ArrayList<ForestInput> inputList = files.get(fi.file);
			if (inputList == null) {
				inputList = new ArrayList<ForestInput>();
				files.put(fi.file, inputList);
			}
			inputList.add(fi);
			if (Thread.interrupted()) {
				interruptLoading();
				return;
			}
		}
		long tempSize = 0;
		final HashMap<String, Integer> progSteps = new HashMap<String, Integer>();
		for (String file : files.keySet()) {
			int sz = Files.getSize(file, false);
			progSteps.put(file, sz);
			tempSize += sz;
			if (Thread.interrupted()) {
				interruptLoading();
				return;
			}
		}
		final long totalSize = tempSize;
		
		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				@Override
				public void run() {
					progressBar.setValue(0);
					progressBar.setStringPainted(true);
					progressBar.setString(null);
					progressBar.setIndeterminate(false);
					progressBar.repaint();
				}
			});
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
		}
		
		for (final java.util.Map.Entry<String, ArrayList<ForestInput>> fileMap : files.entrySet()) {
			if (Thread.interrupted()) {
				interruptLoading();
				return;
			}
			final int progStep = (int) (((double)progSteps.get(fileMap.getKey()) / (double)totalSize) * 100);
			
			int lineStep = 0;
			int lines = 0;
			try {
				LineNumberReader lnr = new LineNumberReader(new java.io.FileReader(fileMap.getKey()));
				lnr.skip(Long.MAX_VALUE);
				lines = lnr.getLineNumber();
				lnr.close();
				lineStep = lines > progStep ? 1 : progStep / lines;
			} catch (IOException e1) {
				// TODO Auto-generated catch block
			}
			if (Thread.interrupted()) {
				interruptLoading();
				return;
			}
			final int lineProg = lineStep;
			final int totalLines = lines;
			
			BufferedReader dataReader;
			dataReader = Files.getReader(fileMap.getKey(), 
											false, // not a jar file
											true, // verbose mode on 
											log, 
											false // don't kill the whole process, esp. if we're running a GUI
											);
			if (dataReader == null) {
				continue;
			}
			if (Thread.interrupted()) {
				interruptLoading();
				return;
			}
			String delimiter = Files.determineDelimiter(fileMap.getKey(), log);
			String header;
			int lineCnt = 1;
			try {
				header = dataReader.readLine();	// skip header
				String[] hdr;
				if (delimiter.startsWith(",")) {
					hdr = ext.splitCommasIntelligently(header, true, log);
				} else {
					hdr = header.trim().split(delimiter);
				}
				int idIndex = ext.indexFactors(new String[][]{Aliases.MARKER_NAMES}, hdr, false, true, false, false)[0];
				
				for (ForestInput inputData : fileMap.getValue()) {
					if (Thread.interrupted()) {
						interruptLoading();
						return;
					}
					try {
						mapMarkersToCol(inputData, header);
					} catch (InterruptedException e) {
						inputData.metaIndicies = null;
						inputData.studyList.clear();
						inputData.studyToColIndexMap.clear();
						interruptLoading();
						return;
					}
				}
				while (dataReader.ready() && !Thread.interrupted()) {
					String readLine = dataReader.readLine();
					lineCnt++;
					if (lineProg > 1 || lineCnt % (totalLines / progStep > 0 ? totalLines / progStep : 1) == 0 /*lineProg > 0*/) {
						try {
							SwingUtilities.invokeAndWait(new Runnable() {
								@Override
								public void run() {
									progressBar.setValue(progressBar.getValue() + lineProg);
									progressBar.repaint();
								}
							});
						} catch (InvocationTargetException e) {
							// TODO Auto-generated catch block
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
						}
					}
					String readData[] = delimiter.equals(",") ? ext.splitCommasIntelligently(readLine, true, log) : readLine.split(delimiter);
					String markerName = readData[idIndex];
					for (ForestInput inputData : fileMap.getValue()) {
						if (Thread.interrupted()) {
							dataReader.close();
							interruptLoading();
							return;
						}
						if(inputData.marker.equals(markerName)) {
//							dataToMetaMap.put(inputData, getMetaStudy(inputData, readData));
							getMetaStudy(inputData, readData);
							atleastOneStudy = true;
						}
					}
				}
				dataReader.close();
				
			} catch (IOException e) {
				log.reportException(e);
			}
	
			if(!atleastOneStudy) {
				log.reportError("Not able to find data for file '"+fileMap.getKey()+"'. Please make sure the given markers are correct and included in data file.");
			}
			
			try {
				SwingUtilities.invokeAndWait(new Runnable() {
					@Override
					public void run() {
						progressBar.setValue(progressBar.getValue() + progStep);
						progressBar.repaint();
					}
				});
			} catch (InvocationTargetException e) {
				// TODO Auto-generated catch block
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
			}

		}
		
		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				@Override
				public void run() {
					progressBar.setValue(0);
					progressBar.setStringPainted(true);
					progressBar.setString("Displaying data...");
					progressBar.setIndeterminate(true);
					progressBar.repaint();
				}
			});
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
		}
	}
	
	private void interruptLoading() {
		this.atleastOneStudy = false;
		this.markerFileName = "";
		this.getDataIndices().clear();
//		this.dataToMetaMap.clear();
		this.clearCurrentData();
		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				@Override
				public void run() {
					progressBar.setValue(0);
					progressBar.setStringPainted(false);
					progressBar.setString("");
					progressBar.setIndeterminate(false);
					progressBar.repaint();
					markerFileList.setSelectedItem(ForestPlot.MARKER_LIST_PLACEHOLDER);
					forestPanel.paintAgain();
					ForestPlot.this.repaint();
				}
			});
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
		}
	}
	
	private void getMetaStudy(ForestInput input, String[] readData) {
		String metaB, metaS;
		float metaBeta, metaStderr;
		
		metaB = readData[input.metaIndicies[0]];
		metaS = readData[input.metaIndicies[1]];
		metaBeta = ext.isValidDouble(metaB) ? Float.parseFloat(metaB) : 0.0f;
		metaStderr = ext.isValidDouble(metaS) ? Float.parseFloat(metaS) : 0.0f;
		
		MetaStudy ms = new MetaStudy(metaBeta, metaStderr);
		ArrayList<StudyData> studies = getStudyEntries(input, readData);
		for (int i = 0; i < studies.size(); i++) {
			ms.addStudy(studies.get(i));
		}
		
		input.setMetaStudy(ms);
//		return ms;
	}

	private ArrayList<StudyData> getStudyEntries(ForestInput input, String[] readData) {
		ArrayList<StudyData> studies = new ArrayList<StudyData>();
		String betaVal, seVal;
		for (int i = input.studyList.size() - 1; i >= 0; i--) {
			String studyName = input.studyList.get(i);
			betaVal = readData[input.studyToColIndexMap.get(studyName)];
			seVal = readData[input.studyToColIndexMap.get(studyName) + 1];
			float beta = ext.isValidDouble(betaVal) ? Float.parseFloat(betaVal) : 0.0f;
			float stderr = ext.isValidDouble(seVal) ? Float.parseFloat(seVal) : 0.0f;
			studies.add(new StudyData(ext.replaceAllWith(studyName, REPLACEMENTS_FOOLISHLY_HARD_CODED), beta, stderr, 0, PlotPoint.FILLED_CIRCLE));
		}
		return studies;
	}
	
	private void mapMarkersToCol(ForestInput data, String hdr) throws RuntimeException, InterruptedException {
		String delim = data.file.toLowerCase().endsWith(".csv") ? ",!" : ext.determineDelimiter(hdr);
		String[] dataFileHeaders = delim.startsWith(",") ? ext.splitCommasIntelligently(hdr, delim.endsWith("!"), log) : hdr.trim().split(delim);
		for (int i = 0; i < dataFileHeaders.length; i++) {
			for (int j = 0; j < BETA_META_HEADERS.length; j++) {
				if (dataFileHeaders[i].toLowerCase().equals(BETA_META_HEADERS[j])) {
					if (dataFileHeaders[i + 1].toLowerCase().startsWith(SE_META_HEADERS[j])) {
						data.metaIndicies = new int[]{i, i+1};
					}
				}
			}
			if (dataFileHeaders[i].toLowerCase().startsWith(BETA_PREFIX)) {
				if (dataFileHeaders[i + 1].toLowerCase().startsWith(SE_PREFIX)) {
					if (data.studyToColIndexMap.containsKey(dataFileHeaders[i].split("\\.")[1])) {
						throw new RuntimeException("Malformed data file: Duplicate study name found in file");
					} else {
						data.addStudy(dataFileHeaders[i].split("\\.")[1], i);
					}
				} else {
					throw new RuntimeException("Malformed data file: SE is not present after Beta for: " + dataFileHeaders[i]);
				}
			}
		}
		if (data.metaIndicies == null) {
			log.reportError("Error - no overall beta/se pairing or effect/stderr pairing was found in file "+data.file);
		}
		dataFileHeaders = null;
	}
	
	private LinkedHashSet<ForestInput> readMarkerFile(String markerFile) {
	    String file;
		LinkedHashSet<ForestInput> markerNames = new LinkedHashSet<ForestInput>();
		BufferedReader markerReader = Files.getReader(markerFile, false, true, false);
		
		if (markerReader != null) {
    		try {
    			while (markerReader.ready() && !Thread.interrupted()) {
    				String[] line = markerReader.readLine().trim().split("\\t");
				    if (line.length >= 2) {
				        file = line[1];
				        if (!file.contains(":") && !file.startsWith("/") && !Files.exists(file)) {
				            if (Files.exists(ext.verifyDirFormat(ext.parseDirectoryOfFile(markerFile)) + file)) {
				                file = ext.verifyDirFormat(ext.parseDirectoryOfFile(markerFile)) + file;
				            } else {
				                if (log != null) {
				                    log.reportError("Error - file " + file + " not found!");
				                } else {
				                    System.err.println("Error - file " + file + " not found!");
				                }
				            }
				        }
    					markerNames.add(new ForestInput(line[0], file, line.length > 2 ? line[2] : ""));
    				} else if (line.length == 1) {
    					markerNames.add(new ForestInput(line[0], "", ""));
    				}
    			}
    		} catch (IOException e) {
    		    if (log != null) {
    		        log.reportException(e);
    		    } else {
    		        e.printStackTrace();
    		    }
    		}
		}
		
		return markerReader == null || Thread.interrupted() ? new LinkedHashSet<ForestInput>() : markerNames;
	}

	public float getMaxZScore() {
		return maxZScore;
	}

	public float getSumZScore() {
		return sumZScore;
	}

	public void updateGUI() {
		forestPanel.paintAgain();
	}

	private void inputMapAndActionMap() {
		InputMap inputMap = forestPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.ALT_MASK), FIRST);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.ALT_MASK), LAST);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), PREVIOUS);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), NEXT);
		ActionMap actionMap = forestPanel.getActionMap();
		actionMap.put(FIRST, navFirstAction);
		actionMap.put(LAST, navLastAction);
		actionMap.put(PREVIOUS, navPrevAction);
		actionMap.put(NEXT, navNextAction);
		forestPanel.setActionMap(actionMap);
	}

	AbstractAction navFirstAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;
		@Override
		public void actionPerformed(ActionEvent e) {
			if(getCurrentDataIndex() != 0){
				setCurrentData(0);
				updateForestPlot();
			}
		}
	};

	AbstractAction navPrevAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;
		@Override
		public void actionPerformed(ActionEvent e) {
			if(getCurrentDataIndex() != 0){
				setCurrentData(getCurrentDataIndex() - 1);
				updateForestPlot();
			}
		}
	};

	AbstractAction navNextAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;
		@Override
		public void actionPerformed(ActionEvent e) {
			if(getCurrentDataIndex() < getDataIndices().size() - 1){
				setCurrentData(getCurrentDataIndex() + 1);
				updateForestPlot();
			}
		}
	};

	AbstractAction navLastAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;
		@Override
		public void actionPerformed(ActionEvent e) {
			if(getCurrentDataIndex() < getDataIndices().size() - 1){
				setCurrentData(getDataIndices().size() - 1);
				updateForestPlot();
			}
		}
	};

	private JProgressBar progressBar;
	private ArrayList<String> sortOrder = null;
	
	public void loadOrderFile(String filename, boolean shouldSort) {
	    if (!Files.exists(filename)) {
	        String msg = "Error - study order file (" + filename + ") not found!";
	        if (log != null) {
	            log.reportError(msg);
	        } else {
	            System.err.println(msg);
	        }
	        return;
	    }
	    ArrayList<String> order = new ArrayList<String>();
	    try {
    	    BufferedReader reader = Files.getAppropriateReader(filename);
    	    String line = null;
            while((line = reader.readLine()) != null) {
                order.add(line.trim());
            }
            reader.close();
        } catch (IOException e) {
            if (proj != null) {
                proj.message("Error occurred while loading study order file: " + e.getMessage());
            }
            if (log != null) {
                log.reportException(e);
            } else {
                e.printStackTrace();
            }
        }
	    if (log != null) {
	        log.report("Loaded Study Order File: " + filename);
	    } else {
	        System.out.println("Loaded Study Order File: " + filename);
	    }
	    this.setSortOrder(order);
	    this.setSortedDisplay(shouldSort);
        if (currMetaStudy != null) {
            currMetaStudy.setSort(this.isSortedDisplay(), this.getSortOrder());
        }
	}
	
	public MetaStudy getCurrentMetaStudy() {
	    if (currMetaStudy != null) {
	        currMetaStudy.setSort(this.isSortedDisplay(), this.getSortOrder());
	    }
		return currMetaStudy;
	}


	private void setCurrentMetaStudy(MetaStudy currMetaStudy) {
		this.currMetaStudy = currMetaStudy;
	}


	public int getCurrentDataIndex() {
		return currentDataIndex;
	}


	private void setCurrentDataIndex(int currentDataIndex) {
		this.currentDataIndex = currentDataIndex;
	}

	public ArrayList<ForestInput> getDataIndices() {
		return dataIndices;
	}


	private void setDataIndices(ArrayList<ForestInput> dataIndices) {
		this.dataIndices = dataIndices;
	}


	public String getPlotLabel() {
		return plotLabel;
	}


	private void setPlotLabel(String plotLabel) {
		this.plotLabel = plotLabel;
	}

    public void setSortedDisplay(boolean sorted) {
        this.sortedDisplay = sorted;
    }
    
	@Override
	public void windowOpened(WindowEvent e) {/**/}

	@Override
	public void windowClosing(WindowEvent e) {
		if (proj == null) {
			this.setVisible(false);
			this.dispose();
			return;
		}
		
//		String[] projFiles = proj.getFilenames(Project.FOREST_PLOT_FILENAMES);
		String[] projFiles = proj.FOREST_PLOT_FILENAMES.getValue();
		String[] currFiles = markerFileNameLoc.values().toArray(new String[]{});
		
		ArrayList<String> newFiles = new ArrayList<String>();
		for (String file : currFiles) {
			boolean found = false;
			for (String oldFile : projFiles) {
				if (oldFile.equals(file)) {
					found = true;
				}
			}
			if (!found) {
				newFiles.add(file);
			}
		}
		
		if (newFiles.size() == 0) {
			this.setVisible(false);
			this.dispose();
			return;
		}
//		String newProp = Array.toStr(currFiles, ";");
		
		String message = newFiles.size() + " files have been added.  ";
		int choice = JOptionPane.showOptionDialog(null, message+" Would you like to keep this configuration for the next time Forest Plot is loaded?", "Preserve Forest Plot workspace?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
		if (choice == 0) {
//			proj.setProperty(Project.FOREST_PLOT_FILENAMES, newProp);
//			proj.setProperty(Project.FOREST_PLOT_FILENAMES, currFiles);
			proj.FOREST_PLOT_FILENAMES.setValue(currFiles);
			proj.saveProperties();
		} else if (choice == -1 || choice == 2) {
			return;
		}

		this.setVisible(false);
		this.dispose();
	}

	@Override
	public void windowClosed(WindowEvent e) {/**/}

	@Override
	public void windowIconified(WindowEvent e) {/**/}

	@Override
	public void windowDeiconified(WindowEvent e) {/**/}

	@Override
	public void windowActivated(WindowEvent e) {/**/}

	@Override
	public void windowDeactivated(WindowEvent e) {/**/}

	public static void main(String[] args) {
				int numArgs = args.length;
		//		String betaSource = "SeqMeta_results.csv";
				String markerList = "markersToDisplay.txt";
	//			String sourceType = "SeqMeta";
				String filename = null;
				String logfile = null;
				final Logger log;
		
				String usage = "\n" + "cnv.plots.ForestPlot requires 1 arguments\n" +
						"  (1) Name of the file with the list of markers (SeqMeta format), files, and comments, to display (i.e. markerList="+ markerList +" (default))\n" +
						"OR\n" +
						"  (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n";
		
				for (int i = 0; i < args.length; i++) {
					if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
						System.err.println(usage);
						System.exit(1);
		//			} else if (args[i].startsWith("betaSource=")) {
		//				betaSource = args[i].split("=")[1];
		//				numArgs--;
	//				} else if (args[i].startsWith("type=")) {
	//					sourceType = args[i].split("=")[1];
	//					numArgs--;
					} else if (args[i].startsWith("proj=")) {
						filename = args[i].split("=")[1];
						numArgs--;
					} else if (args[i].startsWith("markerList=") || args[i].startsWith("file=")) {
						markerList = args[i].split("=")[1];
						numArgs--;
					} else if (args[i].startsWith("log=")) {
						logfile = args[i].split("=")[1];
						numArgs--;
					} else {
						System.err.println("Error - invalid argument: " + args[i]);
					}
				}
				if (numArgs != 0) {
					System.err.println(usage);
					System.exit(1);
				}
				try {
					final Project proj = (filename != null ? new Project(filename, false) : null);
					log = new Logger(logfile);
		
//					markerList = "list.txt";
					
		//			final String finalDataFile = betaSource;
					final String finalMarkerFile = markerList;
					SwingUtilities.invokeLater(new Runnable() {
						public void run() {
							if (proj != null) {
								new ForestPlot(proj);
							} else {
								new ForestPlot(finalMarkerFile, log);
							}
						}
					});
				} catch (Exception e) {
					e.printStackTrace();
				}
			}

    public boolean isSortedDisplay() {
        return sortedDisplay;
    }

    public ArrayList<String> getSortOrder() {
        return sortOrder;
    }

    public void setSortOrder(ArrayList<String> sortOrder) {
        this.sortOrder = sortOrder;
    }
}
