package cnv.plots;

import java.awt.*;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

import javax.swing.*;

import common.Files;
import common.Grafik;
import common.Logger;
import common.ext;

class MetaStudy {
	final ArrayList<StudyData> studies;
	ArrayList<StudyData> sorted;
	private final HashMap<String, StudyData> nameMap;
	final float metaBeta;
	final float metaStderr;
	final float[] metaConf = new float[2];
	
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

	String findLongestStudyNameSize() {
		String longest = "";
		for(StudyData ft : studies){
			longest = longest.length() < ft.getLabel().length() ? ft.getLabel() : longest;
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

	public ArrayList<StudyData> getStudies(boolean sorted) {
		return sorted ? getSorted() : this.studies;
	}
	
	private ArrayList<StudyData> getSorted() {
		if (this.sorted == null) {
			this.sorted = new ArrayList<StudyData>();
			
			TreeMap<String, String> zeroStudyMap = new TreeMap<String, String>();
			TreeMap<Float, String> betaStudyMap = new TreeMap<Float, String>();
			for (StudyData study : studies) {
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

public class ForestPlot extends JPanel implements ActionListener{
	private static final long serialVersionUID = 1L;
	
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String ADD_DATA_FILE = "Add Data File";
	public static final String REMOVE_DATA_FILE = "Remove Data File";
	private static final String ALT_UP = "ALT UP";
	private static final String ALT_DOWN = "ALT DOWN";
	private static final String ALT_LEFT = "ALT LEFT";
	private static final String ALT_RIGHT = "ALT RIGHT";
	private static final String BETA_META_HEADER = "beta";
	private static final String SE_META_HEADER = "se";
	private static final String BETA_PREFIX = "beta.";
	private static final String SE_PREFIX = "se.";
	private static final String FIRST = "First";
	private static final String PREVIOUS = "Previous";
	private static final String NEXT = "Next";
	private static final String LAST = "Last";
	private LinkedHashSet<ForestInput> data;
	ArrayList<ForestInput> dataIndices;
	HashMap<ForestInput, MetaStudy> dataToMetaMap;
	MetaStudy currMetaStudy;
	String plotLabel;
	Logger log;
	ForestPanel forestPanel;
	float maxZScore;
	float sumZScore;
	String longestStudyNameSize;
	private JCheckBox btnSortStudies;
	private JButton flipButton, invXButton, invYButton;
	private JLayeredPane layeredPane;
	private JButton first, previous, next, last;
	private JTextField navigationField;
//	int curMarkerIndex;
	int currentDataIndex;
	private boolean atleastOneStudy;


	public ForestPlot(/*String dataFile, */String markerFile, Logger log) {
		this.log = log;
		atleastOneStudy = false;
		this.data = readMarkerFile(markerFile);
		this.dataIndices = new ArrayList<ForestInput>();
		this.dataIndices.addAll(data);
		dataToMetaMap = new HashMap<ForestInput, MetaStudy>();
		try{
			mapMarkersToCol();
		} catch(RuntimeException rte){
			// TODO handle runtime exceptions!
			log.reportException(rte);
		}

		Iterator<ForestInput> iter = this.data.iterator();
		while(iter.hasNext()) {
			ForestInput data = iter.next();
			try{
				loadStudyData(data);
			} catch (RuntimeException re){
				// TODO handle runtime exceptions!
				log.reportException(re);
			}
		}
		
		setCurrentData(0);
		
		setLayout(new BorderLayout());

		forestPanel = new ForestPanel(this, log);

		layeredPane = new JLayeredPane();
		layeredPane.setLayout(new BorderLayout());

		layeredPane.add(forestPanel);
		layeredPane.setPreferredSize(new Dimension(1000, 600));

		JPanel treePanel = new JPanel();
		treePanel.setBackground(BACKGROUND_COLOR);
		treePanel.setLayout(new BorderLayout());

		//JPanel infoPanel = new JPanel();
		//infoPanel.setBackground(BACKGROUND_COLOR);
		//infoPanel.setLayout(new BoxLayout(infoPanel, BoxLayout.Y_AXIS));

		//JLabel header = new JLabel("Information");
		//header.setMinimumSize(new Dimension(200, 20));
		//header.setMaximumSize(new Dimension(200, 20));
		//header.setAlignmentX(Component.CENTER_ALIGNMENT);
		// button.addActionListener(this);
		//infoPanel.add(header);

		forestPanel.add(markerPanel(), BorderLayout.SOUTH);
		//treePanel.add(infoPanel, BorderLayout.NORTH);


		// initializeTree();
		// updateTree();
		//
		// treePanel.add(new JScrollPane(tree), BorderLayout.CENTER);
		// tree.addTreeSelectionListener(this);
		JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, layeredPane, treePanel);
		splitPane.setBackground(Color.WHITE);
		splitPane.setOneTouchExpandable(true);
		splitPane.setDividerLocation(getWidth() - 150);

		Dimension minimumSize = new Dimension(100, 50);
		layeredPane.setMinimumSize(minimumSize);

		add(splitPane, BorderLayout.CENTER);
		inputMapAndActionMap();

		forestPanel.setPointsGeneratable(true);// zx
		forestPanel.setRectangleGeneratable(true);// zx
		forestPanel.setExtraLayersVisible(new byte[] { 99 });
		updateGUI();
		displayIndex(navigationField);

		forestPanel.grabFocus();

		layeredPane.addComponentListener(new ComponentListener() {
			public void componentShown(ComponentEvent e) {
			}

			public void componentResized(ComponentEvent e) {
				// layeredPane.setSize(new Dimension(getWidth()-20, getHeight()-50)); //???zx
//				flipButton.setBounds(70, layeredPane.getHeight() - 75, 38, 38);
//				invXButton.setBounds(70, layeredPane.getHeight() - 35, 38, 13);
//				invYButton.setBounds(55, layeredPane.getHeight() - 75, 13, 38);
				// layeredPane.repaint();
				// twoDPanel.paintAgain();

			}

			public void componentMoved(ComponentEvent e) {
			}

			public void componentHidden(ComponentEvent e) {
			}
		});

		setVisible(true);
		// generateShortcutMenus();
	}


	public void updateForestPlot(){
		forestPanel.setPointsGeneratable(true);// zx
		forestPanel.setRectangleGeneratable(true);// zx
		forestPanel.setExtraLayersVisible(new byte[] { 99 });
		displayIndex(navigationField);
		updateGUI();
	}
	
	private JPanel markerPanel() {
		JPanel descrPanel = new JPanel();
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
		//navigationField.setEditable(false);
		navigationField.setBackground(BACKGROUND_COLOR);
		navigationField.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					int trav = Integer.valueOf(((JTextField)e.getSource()).getText().split("[\\s]+")[0]).intValue()-1;
					if (trav >=0 && trav < dataIndices.size()) {
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
		next.addActionListener(this);
		next.setActionCommand(NEXT);
		next.setPreferredSize(new Dimension(20, 20));
		last = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif", true));
		last.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif", true));
		last.addActionListener(this);
		last.setActionCommand(LAST);
		last.setPreferredSize(new Dimension(20, 20));
		
		AbstractAction sortAction = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent arg0) {
				 ForestPlot.this.forestPanel.sortedDisplay = btnSortStudies.isSelected();
				 ForestPlot.this.forestPanel.paintAgain();
			}
		};
		btnSortStudies = new JCheckBox(sortAction);
		btnSortStudies.setText("Sort Studies");
		btnSortStudies.setBackground(BACKGROUND_COLOR);
		
		navigationPanel.add(first);
		navigationPanel.add(previous);
		navigationPanel.add(navigationField);
		navigationPanel.add(next);
		navigationPanel.add(last);
		navigationPanel.add(btnSortStudies);

		navigationPanel.setBackground(BACKGROUND_COLOR);
		descrPanel.add(navigationPanel);
		descrPanel.setBackground(BACKGROUND_COLOR);
		return descrPanel;
	}
	
	public void displayIndex(JTextField field) {
		field.setText((currentDataIndex + 1) + " of " + dataIndices.size());
	}
	
	private void setCurrentData(int index) {
		currentDataIndex = index;
		currMetaStudy = dataToMetaMap.get(dataIndices.get(index));
		maxZScore = currMetaStudy.findMaxZScore();
		maxZScore = currMetaStudy.findMaxZScore();
		sumZScore = currMetaStudy.calcSumZScore();
		longestStudyNameSize = currMetaStudy.findLongestStudyNameSize();
		plotLabel = dataIndices.get(index).marker;
	}
	
	public String getLongestStudyNameSize() {
		return longestStudyNameSize;
	}

	private void loadStudyData(ForestInput data) throws RuntimeException{
		BufferedReader dataReader = Files.getReader(data.file, 
														false, // not a jar file
														true, // verbose mode on 
														log, 
														false // DON'T KILL THE WHOLE DAMN SYSTEM (esp. if we're running a GUI)
														);
		String delimiter = Files.determineDelimiter(data.file, log);
		try {
			dataReader.readLine();	// skip header
			
			while(dataReader.ready()) {
				String readLine = dataReader.readLine();
				String readData[] = readLine.split(delimiter);
				String markerName = readData[1];
				if(data.marker.equals(markerName)) {
					dataToMetaMap.put(data, getMetaStudy(data, readData));
					atleastOneStudy = true;
				}
			}
		} catch (IOException e) {
			log.reportException(e);
		}

		if(!atleastOneStudy) {
			throw new RuntimeException("Not able to find data for the given markers. Please make sure the given markers are correct and included in the data file.");
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
	
	private void mapMarkersToCol() throws RuntimeException {
		Iterator<ForestInput> iter = this.data.iterator();
		while(iter.hasNext()) {
			ForestInput data = iter.next();
			
			String[] dataFileHeaders = Files.getHeaderOfFile(data.file, this.log);
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
			dataFileHeaders = null;
			
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
	
	private LinkedHashSet<ForestInput> readMarkerFile(String markerFile) {
		LinkedHashSet<ForestInput> markerNames = new LinkedHashSet<ForestInput>();
		BufferedReader markerReader = Files.getReader(markerFile, false, true, true);
		
		try {
			while (markerReader.ready()) {
				String[] line = markerReader.readLine().trim().split("\\t");
				markerNames.add(new ForestInput(line[0], line[1], line.length > 2 ? line[2] : ""));
			}
		} catch (IOException e) {
			log.reportException(e);
		}
		
		return markerNames;
	}

	public static void createAndShowGUI(String markerFile, Logger log) {

		// Create and set up the window.
		JFrame frame = new JFrame("Forest Plot");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		// Create and set up the content pane.
		log.report("Creating new Forest Plot object");
		ForestPlot forestPlot = new ForestPlot(markerFile, log);
		// frame.setJMenuBar(twoDPlot.menuBar());
		forestPlot.setOpaque(true); // content panes must be opaque
		frame.setContentPane(forestPlot);
		// frame.addWindowListener(twoDPlot);
		frame.setBounds(20, 20, 1000, 600);

		// Display the window.
		frame.pack();
		frame.setVisible(true);
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
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.ALT_MASK), ALT_UP);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.ALT_MASK), ALT_DOWN);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), ALT_LEFT);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), ALT_RIGHT);
		ActionMap actionMap = forestPanel.getActionMap();
		forestPanel.setActionMap(actionMap);
	}

	@Override
	public void actionPerformed(ActionEvent ae) {

		String command = ae.getActionCommand();
//		String filename;
//		int count;

		if (command.equals(FIRST)) {
			first();
		} else if (command.equals(PREVIOUS)) {
			previous();
		} else if (command.equals(NEXT)) {
			next();
		} else if (command.equals(LAST)) {
			last();
		} else {
			log.reportError("Error - unknown command '"+command+"'");
		}

	}

	public void first() {
		if(currentDataIndex != 0){
			setCurrentData(0);
			updateForestPlot();
		}
	}

	public void previous() {
		if(currentDataIndex != 0){
			setCurrentData(currentDataIndex - 1);
			updateForestPlot();
		}
	}

	public void next() {
		if(currentDataIndex != dataIndices.size() - 1){
			setCurrentData(currentDataIndex + 1);
			updateForestPlot();
		}
	}

	public void last() {
		if(currentDataIndex < dataIndices.size() - 1){
			setCurrentData( dataIndices.size() - 1);
			updateForestPlot();
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
//		String betaSource = "SeqMeta_results.csv";
		String markerList = "markersToDisplay.txt";
		String sourceType = "SeqMeta";
		String logfile = null;
		final Logger log;

		String usage = "\n" + "cnv.plots.ForestPlot requires 2 arguments\n" +
//				"(1) Name of the file with betas and standard errors (i.e. betaSource=" + betaSource +"(default))\n" +
				"(1) File type (i.e. type=" + sourceType +" (default))\n" +
				"(2) Name of the file with the list of markers to display (i.e. markerList="+ markerList +" (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
//			} else if (args[i].startsWith("betaSource=")) {
//				betaSource = args[i].split("=")[1];
//				numArgs--;
			} else if (args[i].startsWith("type=")) {
				sourceType = args[i].split("=")[1];
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
			log = new Logger(logfile);

//			final String finalDataFile = betaSource;
			final String finalMarkerFile = markerList;
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					createAndShowGUI(finalMarkerFile, log);
				}
			});
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
