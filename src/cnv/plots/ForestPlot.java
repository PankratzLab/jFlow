package cnv.plots;

import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;

import javax.swing.*;

import common.Grafik;
import common.Logger;

/**
 * Forest Plot class
 * 
 * @author Rohit Sinha
 */
public class ForestPlot extends JPanel {

	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String ADD_DATA_FILE = "Add Data File";
	public static final String REMOVE_DATA_FILE = "Remove Data File";
	private static final String ALT_UP = "ALT UP";
	private static final String ALT_DOWN = "ALT DOWN";
	private static final String ALT_LEFT = "ALT LEFT";
	private static final String ALT_RIGHT = "ALT RIGHT";
	ArrayList<ForestTree> trees;
	String plotLabel;
	Logger log;
	ForestPanel forestPanel;
	float maxZScore;
	private JButton flipButton, invXButton, invYButton;
	private boolean flipStatus, xInvStatus, yInvStatus;
	private JLayeredPane layeredPane;

	public ForestPlot(ArrayList<ForestTree> trees, String plotLabel, Logger log) {
		this.trees = trees;
		this.plotLabel = plotLabel;
		this.log = log;
		this.maxZScore = findMaxZScore();

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

		JPanel infoPanel = new JPanel();
		infoPanel.setBackground(BACKGROUND_COLOR);
		infoPanel.setLayout(new BoxLayout(infoPanel, BoxLayout.Y_AXIS));

		JLabel header = new JLabel("Information");
		header.setMinimumSize(new Dimension(200, 20));
		header.setMaximumSize(new Dimension(200, 20));
		header.setAlignmentX(Component.CENTER_ALIGNMENT);
		// button.addActionListener(this);
		infoPanel.add(header);

		treePanel.add(infoPanel, BorderLayout.NORTH);

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

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "ForestPlot.dat";
		String logfile = null;
		final Logger log;

		String usage = "\n" + "cnv.plots.ForestPlot requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
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
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					createAndShowGUI(getTrees(), "rs12338", log);
				}
			});
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void createAndShowGUI(ArrayList<ForestTree> trees, String plotLabel, Logger log) {

		// Create and set up the window.
		JFrame frame = new JFrame("Forest Plot");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		// Create and set up the content pane.
		log.report("Creating new Forest Plot object");
		ForestPlot forestPlot = new ForestPlot(getTrees(), plotLabel, log);
		// frame.setJMenuBar(twoDPlot.menuBar());
		forestPlot.setOpaque(true); // content panes must be opaque
		frame.setContentPane(forestPlot);
		// frame.addWindowListener(twoDPlot);
		frame.setBounds(20, 20, 1000, 600);

		// Display the window.
		frame.pack();
		frame.setVisible(true);
	}

	public static ArrayList<ForestTree> getTrees() {
		ArrayList<ForestTree> trees = new ArrayList<ForestTree>();
		trees.add(new ForestTree("PROGENI/GenePD", 0.296154f, 0.0834038f, 0, PlotPoint.FILLED_CIRCLE));
		trees.add(new ForestTree("NGRC", 0.105856f, 0.0559677f, 0, PlotPoint.FILLED_CIRCLE));
		trees.add(new ForestTree("23andMe", 0.213202f, 0.027064f, 0, PlotPoint.FILLED_CIRCLE));
		trees.add(new ForestTree("Summary", 0.156191625f, 0.022027313f, 0, PlotPoint.FILLED_CIRCLE));
		return trees;
	}

	public float getMaxZScore() {
		return maxZScore;
	}

	private float findMaxZScore() {
		float max = Float.MIN_VALUE;
		for (ForestTree tree : trees) {
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
		this.zScore = beta / stderr;
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
