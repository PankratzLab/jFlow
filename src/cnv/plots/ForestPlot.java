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

/**
 * Forest Plot class
 *
 * @author Rohit Sinha
 */
public class ForestPlot extends JPanel implements ActionListener{

	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String ADD_DATA_FILE = "Add Data File";
	public static final String REMOVE_DATA_FILE = "Remove Data File";
	private static final String ALT_UP = "ALT UP";
	private static final String ALT_DOWN = "ALT DOWN";
	private static final String ALT_LEFT = "ALT LEFT";
	private static final String ALT_RIGHT = "ALT RIGHT";
	private static final String BETA_PREFIX = "beta.";
	private static final String SE_PREFIX = "se.";
	private static final String FIRST = "First";
	private static final String PREVIOUS = "Previous";
	private static final String NEXT = "Next";
	private static final String LAST = "Last";
	private  String dataFile;
	private LinkedHashSet<String> markers;
	private ArrayList<String> markersIndexes;
	private HashMap<String, Integer> markerToColMap;
	HashMap<String, ArrayList<ForestTree>> markersToTreesMap;
	ArrayList<ForestTree> curTrees;
	String plotLabel;
	Logger log;
	ForestPanel forestPanel;
	float maxZScore;
	float sumZScore;
	String longestStudyNameSize;
	private JButton flipButton, invXButton, invYButton;
	private boolean flipStatus, xInvStatus, yInvStatus;
	private JLayeredPane layeredPane;
	String[] dataFileHeaders;
	private JButton first, previous, next, last;
	private JTextField navigationField;
	int curMarkerIndex;
	private boolean atleastOneTree;


	public ForestPlot(String dataFile, String markerFile, Logger log) {
		this.log = log;
		atleastOneTree = false;
		this.dataFile = dataFile;
		markersIndexes = new ArrayList<String>();
		this.markers = readMarkerNames(markerFile);
		this.dataFileHeaders = Files.getHeaderOfFile(dataFile, this.log);
		markerToColMap = new HashMap<String, Integer>();
		markersToTreesMap = new HashMap<String, ArrayList<ForestTree>>();
		try{
			mapMarkersToCol();
		} catch(RuntimeException rte){
			log.reportException(rte);
		}

		System.out.println(markerToColMap.toString());

		try{
			loadTrees();
		} catch (RuntimeException re){
			log.reportException(re);
			System.exit(1);
		}
		curMarkerIndex = 0;
		setCurTree(markersIndexes.get(0));

		setLayout(new BorderLayout());

		forestPanel = new ForestPanel(this, log);

		layeredPane = new JLayeredPane();
		layeredPane.setLayout(new BorderLayout());
		generateFlipButton();
		layeredPane.add(flipButton);
		generateInvXButton();
		layeredPane.add(invXButton);
		generateInvYButton();
		layeredPane.add(invYButton);

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
				flipButton.setBounds(70, layeredPane.getHeight() - 75, 38, 38);
				invXButton.setBounds(70, layeredPane.getHeight() - 35, 38, 13);
				invYButton.setBounds(55, layeredPane.getHeight() - 75, 13, 38);
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


	public void  updateForestPlot(){
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
					if (trav >=0 && trav < markersIndexes.size()) {
						curMarkerIndex = trav;
						setCurTree(markersIndexes.get(curMarkerIndex));
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
		navigationPanel.add(first);
		navigationPanel.add(previous);
		navigationPanel.add(navigationField);
		navigationPanel.add(next);
		navigationPanel.add(last);

		navigationPanel.setBackground(BACKGROUND_COLOR);
		descrPanel.add(navigationPanel);
		descrPanel.setBackground(BACKGROUND_COLOR);
		return descrPanel;
	}

	public void displayIndex(JTextField field) {
		field.setText((curMarkerIndex + 1) + " of " + markersIndexes.size());
	}

	private void setCurTree(String markerName) {
		curTrees = markersToTreesMap.get(markerName);
		maxZScore = findMaxZScore();
		sumZScore = calcSumZScore();
		longestStudyNameSize = findLongestStudyNameSize();
		plotLabel = markerName;
	}

	private String findLongestStudyNameSize() {
		String longest = "";
		for(ForestTree ft : curTrees){
			longest = longest.length() < ft.getLabel().length() ? ft.getLabel() : longest;
		}
		return longest;
	}

	public String getLongestStudyNameSize() {
		return longestStudyNameSize;
	}

	private float calcSumZScore() {
		float sum=0;
		for	(ForestTree ft : curTrees){
			sum += ft.getzScore();
		}
		return sum;
	}

	private void loadTrees() throws RuntimeException{
		BufferedReader dataReader = Files.getReader(dataFile, false,true, log, true);
		String delimiter = Files.determineDelimiter(dataFile, log);
		try {
			dataReader.readLine();	// skip headers

			while(dataReader.ready()){
				String readLine = dataReader.readLine();
				String readData[] = readLine.split(delimiter);
				if(markers.contains(readData[1])){
					getAllTrees(readData);
					atleastOneTree = true;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		if(!atleastOneTree)
			throw new RuntimeException("Not able to find trees with the given markers. Please make sure that markers are right");
	}

	private void getAllTrees(String[] readData) {
		if(markersToTreesMap.containsKey(readData[1])){
			throw new RuntimeException("Malformed Data File: Duplicate marker name found");
		} else {
			markersToTreesMap.put(readData[1], getTrees(readData));
		}
	}

	private ArrayList<ForestTree> getTrees(String[] readData) {
		ArrayList<ForestTree> trees = new ArrayList<ForestTree>();
		for(String studyName : markerToColMap.keySet()){
			float beta = ext.isValidDouble(readData[markerToColMap.get(studyName)]) ? Float.parseFloat(readData[markerToColMap.get(studyName)]) : 0.0f;
			float stderr =  ext.isValidDouble(readData[markerToColMap.get(studyName) + 1]) ? Float.parseFloat(readData[markerToColMap.get(studyName) + 1]) : 0.0f;
			trees.add(new ForestTree(studyName, beta, stderr,  0, PlotPoint.FILLED_CIRCLE));
		}
		return trees;
	}

	private void mapMarkersToCol() throws RuntimeException {

		for (int i = 0; i < dataFileHeaders.length; i++) {
			if (dataFileHeaders[i].toLowerCase().startsWith(BETA_PREFIX)) {
				if (dataFileHeaders[i + 1].toLowerCase().startsWith(SE_PREFIX)) {
					if (markerToColMap.containsKey(dataFileHeaders[i].split("\\.")[1])) {
						throw new RuntimeException("Malformed data file: Duplicate study name found in file");
					} else {
						markerToColMap.put(dataFileHeaders[i].split("\\.")[1], i);
					}
				} else {
					throw new RuntimeException("Malformed data file: SE is not present after Beta for: " + dataFileHeaders[i]);
				}
			}
		}
	}

	private LinkedHashSet<String> readMarkerNames(String markerFile) {
		LinkedHashSet<String> markerNames = new LinkedHashSet<String>();
		BufferedReader markerReader = Files.getReader(markerFile, false, true, true);

		try {
			while (markerReader.ready()) {
				markerNames.add(markerReader.readLine().trim());

			}
		} catch (IOException e) {
			log.reportException(e);
		}
		markersIndexes.addAll(markerNames);
		return markerNames;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String betaSource = "SeqMeta_results.csv";
		String markerList = "markersToDisplay.txt";
		String sourceType = "SeqMeta";
		String logfile = null;
		final Logger log;

		String usage = "\n" + "cnv.plots.ForestPlot requires 3 arguments\n" +
				"(1) Name of the file with betas and standard errors (i.e. betaSource=" + betaSource +"(default))\n" +
				"(2) File type (i.e. type=" + sourceType +" (default))" +
				"(3) Name of the file with the list of markers to display (i.e. markerList="+ markerList +" (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("betaSource=")) {
				betaSource = args[i].split("=")[1];
				numArgs--;
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

			final String finalDataFile = betaSource;
			final String finalMarkerFile = markerList;
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					createAndShowGUI(finalDataFile, finalMarkerFile, log);
				}
			});
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void createAndShowGUI(String dataFile, String markerFile, Logger log) {

		// Create and set up the window.
		JFrame frame = new JFrame("Forest Plot");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		// Create and set up the content pane.
		log.report("Creating new Forest Plot object");
		ForestPlot forestPlot = new ForestPlot(dataFile, markerFile, log);
		// frame.setJMenuBar(twoDPlot.menuBar());
		forestPlot.setOpaque(true); // content panes must be opaque
		frame.setContentPane(forestPlot);
		// frame.addWindowListener(twoDPlot);
		frame.setBounds(20, 20, 1000, 600);

		// Display the window.
		frame.pack();
		frame.setVisible(true);
	}

//	public static ArrayList<ForestTree> getTrees() {
//		ArrayList<ForestTree> trees = new ArrayList<ForestTree>();
//		trees.add(new ForestTree("PROGENI/GenePD", 0.296154f, 0.0834038f, 0, PlotPoint.FILLED_CIRCLE));
//		trees.add(new ForestTree("NGRC", 0.105856f, 0.0559677f, 0, PlotPoint.FILLED_CIRCLE));
//		trees.add(new ForestTree("23andMe", 0.213202f, 0.027064f, 0, PlotPoint.FILLED_CIRCLE));
//		trees.add(new ForestTree("Summary", 0.156191625f, 0.022027313f, 0, PlotPoint.FILLED_CIRCLE));
//		return trees;
//	}

	public float getMaxZScore() {
		return maxZScore;
	}

	public float getSumZScore() {
		return sumZScore;
	}

	private float findMaxZScore() {
		float max = Float.MIN_VALUE;
		for (ForestTree tree : curTrees) {
			max = Math.max(max, tree.getzScore());
		}
		return max;
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

	private void generateFlipButton() {
		flipButton = new JButton(Grafik.getImageIcon("images/flip_and_invert/flip_10p.jpg", false));
		flipButton.setRolloverIcon(Grafik.getImageIcon("images/flip_and_invert/flip_10p_blue.jpg", false));
		flipButton.setToolTipText("Inverts axes");
		flipButton.setBorder(null);
		flipButton.setVisible(true);
		flipStatus = true;
		flipButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				// twoDPanel.setSwapable(flipStatus);
				forestPanel.setPointsGeneratable(true);
				forestPanel.setRectangleGeneratable(true);// zx
				forestPanel.setSwapAxes(flipStatus);
				forestPanel.paintAgain();
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
				forestPanel.setPointsGeneratable(true);
				forestPanel.setRectangleGeneratable(true);// zx
				forestPanel.setXinversion(xInvStatus);
				forestPanel.paintAgain();
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
				forestPanel.setPointsGeneratable(true);
				forestPanel.setRectangleGeneratable(true);// zx
				forestPanel.setYinversion(yInvStatus);
				forestPanel.paintAgain();
				if (yInvStatus) {
					invYButton.setIcon(Grafik.getImageIcon("images/flip_and_invert/down_10.gif", true));
				} else {
					invYButton.setIcon(Grafik.getImageIcon("images/flip_and_invert/up_10.gif", true));
				}
				yInvStatus = !yInvStatus;
			}
		});
	}

	@Override
	public void actionPerformed(ActionEvent ae) {

		String command = ae.getActionCommand();
		String filename;
		int count;

		if (command.equals(FIRST)) {
			first();
		} else if (command.equals(PREVIOUS)) {
			previous();
		} else if (command.equals(NEXT)) {
			next();
		} else if (command.equals(LAST)) {
			last();
		}	else {
			log.reportError("Error - unknown command '"+command+"'");
		}

	}

	public void first() {
		if(curMarkerIndex != 0){
			curMarkerIndex = 0;
			setCurTree(markersIndexes.get(curMarkerIndex));
			updateForestPlot();
		}
	}

	public void previous() {
		if(curMarkerIndex != 0){
			curMarkerIndex--;
			setCurTree(markersIndexes.get(curMarkerIndex));
			updateForestPlot();
		}
	}

	public void next() {
		if(curMarkerIndex != markersIndexes.size()-1){
			curMarkerIndex++;
			setCurTree(markersIndexes.get(curMarkerIndex));
			updateForestPlot();
		}
	}

	public void last() {
		if(curMarkerIndex < markersIndexes.size()-1){
			curMarkerIndex = markersIndexes.size()-1;
			setCurTree(markersIndexes.get(curMarkerIndex));
			updateForestPlot();
		}
	}
}

class ForestTree {
	public String label;
	public float beta;
	public float stderr;
	public int color;
	public byte shape;
	public float[] confInterval;
	public float zScore;

	public ForestTree(String label, float beta, float stderr, int color, byte shape) {
		this.label = label;
		this.beta = beta;
		this.stderr = stderr;
		this.color = color;
		this.shape = shape;
		this.confInterval = new float[2];
		this.confInterval[0] = (float) (beta - 1.96 * stderr);
		this.confInterval[1] = (float) (beta + 1.96 * stderr);
		this.zScore = stderr == 0.0f ? 0.0f : beta / stderr;
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

	public float getzScore() {
		return zScore;
	}
}
