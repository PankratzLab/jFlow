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
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.TreeMap;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLayeredPane;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;

import cnv.filesys.Project;
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
	
	public MetaStudy(float metaBeta, float metaStderr) {
		studies = new ArrayList<StudyData>();
		nameMap = new HashMap<String, StudyData>();
		this.metaBeta = metaBeta;
		this.metaStderr = metaStderr;
		this.getMetaConf()[0] = (float) (metaBeta - 1.96 * metaStderr);
		this.getMetaConf()[1] = (float) (metaBeta + 1.96 * metaStderr);
	}
	
	public void addStudy(StudyData studyData) {
		getStudies().add(studyData);
		nameMap.put(studyData.getLabel(), studyData);
	}

	String findLongestStudyName() {
		String longest = "";
		for(StudyData ft : getStudies()){
			longest = longest.length() < ft.getLabel().length() ? ft.getLabel() : longest;
		}
		return longest;
	}
	
	float calcSumZScore() {
		float sum = 0;
		for	(StudyData ft : getStudies()){
			sum += ft.getZScore();
		}
		return sum;
	}
	
	float findMaxZScore() {
		float max = Float.MIN_VALUE;
		for (StudyData data: getStudies()) {
			max = Math.max(max, data.getZScore());
		}
		return max;
	}

	public ArrayList<StudyData> getStudies(boolean sorted) {
		return sorted ? getSorted() : this.studies;
	}
	
	private ArrayList<StudyData> getSorted() {
		if (this.sorted == null) {
			this.sorted = new ArrayList<StudyData>();
			
			TreeMap<String, String> zeroStudyMap = new TreeMap<String, String>();
			TreeMap<Float, String> betaStudyMap = new TreeMap<Float, String>();
			for (StudyData study : getStudies()) {
				if (study.getBeta() == 0.0f) {
					zeroStudyMap.put(study.getLabel(), study.getLabel());
				} else {
					betaStudyMap.put(study.getBeta(), study.getLabel());
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
		
		return this.sorted;
	}

	public ArrayList<StudyData> getStudies() {
		return studies;
	}

	public float[] getMetaConf() {
		return metaConf;
	}

	public float getMetaBeta() {
		return metaBeta;
	}

	public float getMetaStderr() {
		return metaStderr;
	}
	
}

class StudyData {
	private final String label;
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

	public String getLabel() {
		return label;
	}

	public float getBeta() {
		return beta;
	}

	public float getStderr() {
		return stderr;
	}

	public int getColor() {
		return color;
	}

	public byte getShape() {
		return shape;
	}

	public float[] getConfInterval() {
		return confInterval;
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
	
	public ForestInput(String marker, String file, String comment) {
		this.marker = marker;
		this.file = file;
		this.comment = comment;
		this.studyToColIndexMap = new HashMap<String, Integer>();
	}
	
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
	private static final String BETA_META_HEADER = "beta";
	private static final String SE_META_HEADER = "se";
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
	private HashMap<ForestInput, MetaStudy> dataToMetaMap = new HashMap<ForestInput, MetaStudy>();
	private MetaStudy currMetaStudy;
	private String plotLabel;
	private Logger log;
	private ForestPanel forestPanel;
	private float maxZScore;
	private float sumZScore;
	private String longestStudyName;
	private JCheckBox chkSortStudies;
	private JLayeredPane layeredPane;
	private JButton first, previous, next, last;
	private JButton btnScreen, btnScreenAll;
	private JTextField navigationField;
	private JComboBox<String> markerFileList;
	private HashMap<String, String> markerFileNameLoc = new HashMap<String, String>();
//	int curMarkerIndex;
	private int currentDataIndex;
	private boolean atleastOneStudy;

	private String markerFileName;

	protected Project proj;

	private volatile boolean loadingFile;
	private Thread loadingThread;
	
	public ForestPlot(Project proj) {
		super("Genvisis - Forest Plot - " + proj.getNameOfProject());
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
		super("Forest Plot");
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
		
		setVisible(true);
		
		loadMarkerFile();
		// generateShortcutMenus();
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
		
		AbstractAction sortAction = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent arg0) {
				 ForestPlot.this.forestPanel.setSortedDisplay(chkSortStudies.isSelected());
				 ForestPlot.this.forestPanel.paintAgain();
			}
		};
		chkSortStudies = new JCheckBox(sortAction);
		chkSortStudies.setText("Sort Studies");
		chkSortStudies.setBackground(BACKGROUND_COLOR);
		chkSortStudies.setAlignmentX(Component.RIGHT_ALIGNMENT);
		
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
		btnScreen = new JButton(screenAction1);
		btnScreen.setText("Screen Capture");
		btnScreen.setAlignmentX(Component.RIGHT_ALIGNMENT);
		btnScreen.setMinimumSize(new Dimension(140, 26));
		btnScreen.setPreferredSize(new Dimension(140, 26));
		btnScreen.setMaximumSize(new Dimension(140, 26));
		
		btnScreenAll = new JButton(screenAction2);
		btnScreenAll.setText("Screen Capture All");
		btnScreenAll.setAlignmentX(Component.RIGHT_ALIGNMENT);
		
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
			String[] files = proj.getFilenames(Project.FOREST_PLOT_FILES);
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
		
		final JPanel subOptPanel = new JPanel();
		subOptPanel.setLayout(new BoxLayout(subOptPanel, BoxLayout.Y_AXIS));
		subOptPanel.setBackground(BACKGROUND_COLOR);
		subOptPanel.setAlignmentX(Component.RIGHT_ALIGNMENT);
		subOptPanel.add(btnScreen);
		subOptPanel.add(Box.createRigidArea(new Dimension(0,5)));
		subOptPanel.add(btnScreenAll);
		
		JPanel optPanel = new JPanel() {
			private static final long serialVersionUID = 1L;
			@Override
			public Dimension getMaximumSize() {
				return super.getPreferredSize();
			}
		};
		optPanel.setBackground(BACKGROUND_COLOR);
		optPanel.add(chkSortStudies);
		optPanel.add(subOptPanel);
	
		
		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new BoxLayout(descrPanel, BoxLayout.X_AXIS));
		// All of these are needed (or at least, have some effect)
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(navigationPanel);
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(optPanel);
		descrPanel.setBackground(BACKGROUND_COLOR);
		return descrPanel;
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
		
		dataToMetaMap = new HashMap<ForestInput, MetaStudy>();

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

	private void updateForestPlot(){
		forestPanel.setPointsGeneratable(true);
		forestPanel.setRectangleGeneratable(true);
		forestPanel.setExtraLayersVisible(new byte[] { 99 });
		displayIndex(navigationField);
		updateGUI();
	}
	
	private void chooseNewFiles() {
		JFileChooser jfc = new JFileChooser((proj != null ? proj.getProjectDir() : ext.parseDirectoryOfFile(markerFileName)));
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
	
	private void screenCapAll() {
		int currentSelection = getCurrentDataIndex();
		ArrayList<ForestInput> data = getDataIndices();
		for (int i = 0; i < data.size(); i++) {
			setCurrentData(i);
			forestPanel.createImage();
			String marker = "", filename = "", ctrlFile = "";
			int count = 1;
			String root = (proj == null ? ext.parseDirectoryOfFile(markerFileName) : proj.getProjectDir());
			marker = getDataIndices().get(getCurrentDataIndex()).marker;
			ctrlFile = ext.rootOf(markerFileName);
			filename = marker + "_" + ctrlFile;
			filename = ext.replaceWithLinuxSafeCharacters(filename, true);
			while (new File(root+filename+".png").exists()) {
				filename = marker + "_" + ctrlFile + "_v" + count;
				filename = ext.replaceWithLinuxSafeCharacters(filename, true);
				count++;
			}
			log.report("Writing screenshot to file " + root + filename + ".png");
			ForestPlot.this.forestPanel.screenCapture(root+filename+".png");
		}
		setCurrentData(currentSelection);
		updateForestPlot();
	}
	
	private void screenCap() {
		String marker = "", filename = "", ctrlFile = "";
		int count = 1;
		String root = (proj == null ? ext.parseDirectoryOfFile(markerFileName) : proj.getProjectDir());
		marker = getDataIndices().get(getCurrentDataIndex()).marker;
		ctrlFile = ext.rootOf(markerFileName);
		filename = marker + "_" + ctrlFile;
		filename = ext.replaceWithLinuxSafeCharacters(filename, true);
		while (new File(root+filename+".png").exists()) {
			filename = marker + "_" + ctrlFile + "_v" + count;
			filename = ext.replaceWithLinuxSafeCharacters(filename, true);
			count++;
		}
		log.report("Writing screenshot to file " + root + filename + ".png");
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
		setCurrentMetaStudy(dataToMetaMap.get(getDataIndices().get(index)));
		if (getCurrentMetaStudy() == null) {
			log.reportError("Error - could not set index to "+index+" since the data did not load properly; check to see if any results files are missing");
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
			
			BufferedReader dataReader = Files.getReader(fileMap.getKey(), 
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
				
				for (ForestInput inputData : fileMap.getValue()) {
					if (Thread.interrupted()) {
						interruptLoading();
						return;
					}
					try {
						mapMarkersToCol(inputData, header);
					} catch (InterruptedException e) {
						inputData.metaIndicies = null;
						inputData.studyToColIndexMap.clear();
						interruptLoading();
						return;
					}
				}
				while (dataReader.ready() && !Thread.interrupted()) {
					String readLine = dataReader.readLine();
					lineCnt++;
					if (lineProg > 1 || lineCnt % (totalLines / progStep) == 0 /*lineProg > 0*/) {
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
					String markerName = readData[1];
					for (ForestInput inputData : fileMap.getValue()) {
						if (Thread.interrupted()) {
							dataReader.close();
							interruptLoading();
							return;
						}
						if(inputData.marker.equals(markerName)) {
							dataToMetaMap.put(inputData, getMetaStudy(inputData, readData));
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
		this.dataToMetaMap.clear();
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
	
	private MetaStudy getMetaStudy(ForestInput input, String[] readData) {
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
		
		return ms;
	}

	private ArrayList<StudyData> getStudyEntries(ForestInput input, String[] readData) {
		ArrayList<StudyData> studies = new ArrayList<StudyData>();
		String betaVal, seVal;
		for(String studyName : input.studyToColIndexMap.keySet()){
			betaVal = readData[input.studyToColIndexMap.get(studyName)];
			seVal = readData[input.studyToColIndexMap.get(studyName) + 1];
			float beta = ext.isValidDouble(betaVal) ? Float.parseFloat(betaVal) : 0.0f;
			float stderr = ext.isValidDouble(seVal) ? Float.parseFloat(seVal) : 0.0f;
			studies.add(new StudyData(studyName, beta, stderr, 0, PlotPoint.FILLED_CIRCLE));
		}
		return studies;
	}
	
	private void mapMarkersToCol(ForestInput data, String hdr) throws RuntimeException, InterruptedException {
		String delim = data.file.toLowerCase().endsWith(".csv") ? ",!" : ext.determineDelimiter(hdr);
		String[] dataFileHeaders = delim.startsWith(",") ? ext.splitCommasIntelligently(hdr, delim.endsWith("!"), log) : hdr.trim().split(delim);
		for (int i = 0; i < dataFileHeaders.length; i++) {
			if (dataFileHeaders[i].toLowerCase().equals(BETA_META_HEADER)) {
				if (dataFileHeaders[i + 1].toLowerCase().startsWith(SE_META_HEADER)) {
					data.metaIndicies = new int[]{i, i+1};
				}
			} else if (dataFileHeaders[i].toLowerCase().startsWith(BETA_PREFIX)) {
				if (dataFileHeaders[i + 1].toLowerCase().startsWith(SE_PREFIX)) {
					if (data.studyToColIndexMap.containsKey(dataFileHeaders[i].split("\\.")[1])) {
						throw new RuntimeException("Malformed data file: Duplicate study name found in file");
					} else {
						data.studyToColIndexMap.put(dataFileHeaders[i].split("\\.")[1], i);
					}
				} else {
					throw new RuntimeException("Malformed data file: SE is not present after Beta for: " + dataFileHeaders[i]);
				}
			}
		}
		if (data.metaIndicies == null) {
			log.reportError("Error - never found the beta/stderr pairing for file "+data.file);
		}
		dataFileHeaders = null;
	}
	
	private LinkedHashSet<ForestInput> readMarkerFile(String markerFile) {
		LinkedHashSet<ForestInput> markerNames = new LinkedHashSet<ForestInput>();
		BufferedReader markerReader = Files.getReader(markerFile, false, true, true);
		
		try {
			while (markerReader.ready() && !Thread.interrupted()) {
				String[] line = markerReader.readLine().trim().split("\\t");
				if (line.length >= 2) {
					markerNames.add(new ForestInput(line[0], line[1], line.length > 2 ? line[2] : ""));
				} else if (line.length == 1) {
					markerNames.add(new ForestInput(line[0], "", ""));
				}
			}
		} catch (IOException e) {
			log.reportException(e);
		}
		
		return Thread.interrupted() ? new LinkedHashSet<ForestInput>() : markerNames;
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
	
	public MetaStudy getCurrentMetaStudy() {
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


	@Override
	public void windowOpened(WindowEvent e) {/**/}

	@Override
	public void windowClosing(WindowEvent e) {
		if (proj == null) {
			this.setVisible(false);
			this.dispose();
			return;
		}
		
		String[] projFiles = proj.getFilenames(Project.FOREST_PLOT_FILES);
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
		String newProp = Array.toStr(currFiles, ";");
		
		String message = newFiles.size() + " files have been added.  ";
		int choice = JOptionPane.showOptionDialog(null, message+" Would you like to keep this configuration for the next time Forest Plot is loaded?", "Preserve Forest Plot workspace?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
		if (choice == 0) {
			proj.setProperty(Project.FOREST_PLOT_FILES, newProp);
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
					} else if (args[i].startsWith("markerList=")) {
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
}
