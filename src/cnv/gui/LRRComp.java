package cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.Border;

import common.Logger;
import common.ext;
import cnv.analysis.MedianLRRWorker;
import cnv.filesys.Project;
import cnv.manage.Transforms;
import cnv.plots.LinePlot;
import cnv.plots.TwoDPlot;

public class LRRComp extends JFrame implements Runnable {
	private static final long serialVersionUID = 1L;
	// once a job has been started, used to track completion
	private volatile int computeComplete = 42;
	public static final String BLANKREGION = "Paste new Region here...";
	public static final String[] REGION_TEXT_FIELD_LABELS = { "Input UCSC or probeset-based regions of Interest (one per Line):", "Progress..." };
	public static final String[] CLASSES_TO_DUMP = { "IID" };

	private int transformationType;
	private int scope;
	private String initRegion;
	private MedianLRRWorker medianLRRWorker;
	private Project proj;
	private String outputBase;

	public LRRComp(Project proj, String initRegion) {
		this.proj = proj;
		this.initRegion = initRegion;
		this.transformationType = 0;
		this.scope = 0;
	}

	public void run() {
		createAndShowGUI();
	}

	public void createAndShowGUI() {
		setSize(500, 360);
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		setTitle("Median Log R Ratio Settings");
		WindowListener exitListener = new WindowAdapter() {
			// check for running job before exit
			@Override
			public void windowClosing(WindowEvent e) {
				String ObjButtons[] = { "Exit Anyway", "Cancel" };
				if (computeComplete != 42 && !medianLRRWorker.isDone()) {
					int promptResult = JOptionPane.showOptionDialog(null, "A Thread is computing Median Log R Ratios\nThread will continue if you exit", "Warning - Running thread", JOptionPane.DEFAULT_OPTION, JOptionPane.WARNING_MESSAGE, null, ObjButtons, ObjButtons[1]);
					System.out.println(promptResult);
					if (promptResult == 0) {
						dispose();
					}
				} else {
					dispose();
				}
			}
		};
		this.addWindowListener(exitListener);
		TransformationPanel transformationPanel = new TransformationPanel();
		transformationPanel.setLayout(new BoxLayout(transformationPanel, BoxLayout.Y_AXIS));
		add(transformationPanel);
		this.setLocation(dim.width / 2 - this.getSize().width / 2, dim.height / 2 - this.getSize().height / 2);
		setVisible(true);
	}

	private class TransformationPanel extends JPanel implements ActionListener {
		private static final long serialVersionUID = 1L;
		private RegionTextField regionTextField;
		private ComputeButton computeButton;
		private TwoDPlotButton twoDPlotButton;
		private JProgressBar progressBar;

		private TransformationPanel() {
			outputBase = Transforms.TRANFORMATIONS[transformationType];
			this.setLayout(new BorderLayout());
			this.regionTextField = new RegionTextField(initRegion, 10, 10);
			progressBar = new JProgressBar(0, 100);
			computeButton = new ComputeButton(this);
			twoDPlotButton = new TwoDPlotButton(this);
			addLabel("Select a Log R Ratio transformation: ");
			ActionListener actionListener = getradioListener();
			addTransformButtons(actionListener, transformationType);
			addLabel("Transform by: ");
			addScopeButtons(actionListener, scope);
			addLabel(REGION_TEXT_FIELD_LABELS[0]);
			add(regionTextField, BorderLayout.CENTER);
			JScrollPane scroll = new JScrollPane(regionTextField);
			add(scroll);
			add(computeButton, BorderLayout.EAST);
			// TODO add action to launch 2D plot with created file
			add(twoDPlotButton, BorderLayout.WEST);
		}

		@Override
		public void actionPerformed(ActionEvent actionEvent) {
			JComponent source = (JComponent) actionEvent.getSource();
			if (source.equals(computeButton)) {
				if (computeComplete != 42 && !medianLRRWorker.isDone()) {
					JOptionPane.showMessageDialog(this, "Thread is busy computing median Log R Ratios");
				} else {
					startJob();
				}
			} else if (source.equals(twoDPlotButton)) {
				javax.swing.SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						LinePlot linePlot = new LinePlot(proj, new Logger());
					}
				});
			}
		}

		private void startJob() {
			String ObjButtons[] = { "OK", "Cancel" };
			int promptResult = JOptionPane.showOptionDialog(this, "Compute median using " + Transforms.TRANFORMATIONS[transformationType] + "?", "Log R Ratio " + Transforms.TRANFORMATIONS[transformationType], JOptionPane.DEFAULT_OPTION, JOptionPane.WARNING_MESSAGE, null, ObjButtons, ObjButtons[1]);
			if (promptResult == 0) {
				outputBase = outputBase.replaceAll(" ", "_");
				outputBase = outputBase.replaceAll("-", "_");
				String promptName = JOptionPane.showInputDialog(this, "Add an Analysis Name to " + outputBase + "?");
				outputBase = promptName + "_" + outputBase;
				System.out.println(outputBase + proj.getProjectDir());
				this.add(progressBar);
				progressBar.setVisible(true);
				progressBar.setStringPainted(true);
				computeComplete = 0;
				medianLRRWorker = new MedianLRRWorker(proj, regionTextField.getText().split("\n"), transformationType, scope, outputBase, progressBar, null);
				medianLRRWorker.execute();
				this.revalidate();
			}
		}

		private JLabel addLabel(String text) {
			JLabel label = new JLabel(text);
			label.setFont(new Font("Arial", 0, 14));
			this.add(label);
			return label;
		}

		// buttons for scope of transformation
		private void addScopeButtons(ActionListener actionListener, int initTransformationType) {
			ButtonGroup scopeRadio = new ButtonGroup();
			JRadioButton[] scopeRadioButtons = new JRadioButton[Transforms.SCOPES.length];
			for (int i = 0; i < Transforms.SCOPES.length; i++) {
				scopeRadioButtons[i] = new JRadioButton(Transforms.SCOPES[i], false);
				scopeRadioButtons[i].setFont(new Font("Arial", 0, 14));
				scopeRadio.add(scopeRadioButtons[i]);
				scopeRadioButtons[i].addActionListener(actionListener);
				this.add(scopeRadioButtons[i]);
			}
			// not initiating scope since not valid for raw values
			// scopeRadioButtons[initTransformationType].setSelected(true);
		}

		// button for type of transformation
		private void addTransformButtons(ActionListener actionListener, int initScope) {
			ButtonGroup typeRadio = new ButtonGroup();
			JRadioButton[] transformationRadioButtons = new JRadioButton[Transforms.TRANFORMATIONS.length];
			for (int i = 0; i < Transforms.TRANFORMATIONS.length; i++) {
				transformationRadioButtons[i] = new JRadioButton(Transforms.TRANFORMATIONS[i], false);
				transformationRadioButtons[i].setFont(new Font("Arial", 0, 14));
				typeRadio.add(transformationRadioButtons[i]);
				transformationRadioButtons[i].addActionListener(actionListener);
				this.add(transformationRadioButtons[i]);
			}
			transformationRadioButtons[initScope].setSelected(true);
		}
	}

	// input area
	private class RegionTextField extends JTextArea {
		private static final long serialVersionUID = 1L;

		private RegionTextField(String region, int size1, int size2) {
			this.setFont(new Font("Tahoma", Font.PLAIN, 14));
			this.setText(region);
			this.setSize(size1, size2);
			Border border = BorderFactory.createLineBorder(Color.BLACK);
			this.setBorder(border);
			this.setMaximumSize(this.getPreferredSize());
			// this.addActionListener(actionListener);
		}
	}

	private class ComputeButton extends JButton {
		private static final long serialVersionUID = 1L;

		private ComputeButton(ActionListener actionListener) {
			this.setFont(new Font("Tahoma", Font.PLAIN, 14));
			this.setText("Compute");
			this.setToolTipText("Compute Median Log R Ratios");
			this.addActionListener(actionListener);
		}
	}

	private class TwoDPlotButton extends JButton {
		private static final long serialVersionUID = 1L;

		private TwoDPlotButton(ActionListener actionListener) {
			this.setFont(new Font("Tahoma", Font.PLAIN, 14));
			this.setText("2D Plot");
			this.setToolTipText("Launch 2D plot to visualize results");
			this.addActionListener(actionListener);
		}
	}

	public ActionListener getradioListener() {
		return new ActionListener() {
			public void actionPerformed(ActionEvent actionEvent) {
				if (ext.indexOfStr(actionEvent.getActionCommand(), Transforms.TRANFORMATIONS, true, true, new Logger(), false) >= 0) {
					transformationType = ext.indexOfStr(actionEvent.getActionCommand(), Transforms.TRANFORMATIONS, true, true, new Logger(), false);
					outputBase = Transforms.TRANFORMATIONS[transformationType];
				} else if (ext.indexOfStr(actionEvent.getActionCommand(), Transforms.SCOPES, true, true, new Logger(), false) >= 0 && transformationType == 0) {
					JOptionPane.showMessageDialog(null, "Transform by Chromosome or Genome not valid for Raw Values");
				} else if (ext.indexOfStr(actionEvent.getActionCommand(), Transforms.SCOPES, true, true, new Logger(), false) >= 0) {
					scope = ext.indexOfStr(actionEvent.getActionCommand(), Transforms.SCOPES, true, true, new Logger(), false);
					outputBase = Transforms.TRANFORMATIONS[transformationType] + "_" + Transforms.SCOPES[scope];
				} else {
					System.err.println("Error - could not find transformation type");
				}
				System.out.println(actionEvent.getActionCommand() + "\t" + transformationType + "\t" + scope);
			}
		};
	}
}

// String[] classes = { "IID" };
// // public static void dump(Project proj, String[] phenotypes, String mlrrSetFile, String regionToDumpOrNullForAll, int transformationToUse, Logger log) {
// MedianLRR.dump(proj, classes, proj.getProjectDir() + outputBase + ".mlrr", null, 0, log);
// TwoDPlot twoDPlot = new TwoDPlot(proj, log);
// twoDPlot.loadFile(proj.getProjectDir() + outputBase + "_dump.xln");

// private class ValidateButton extends JButton {
// private static final long serialVersionUID = 1L;
//
// private ValidateButton(ActionListener actionListener) {
// this.setFont(new Font("Tahoma", Font.PLAIN, 14));
// this.setText("Validate Regions");
// this.setToolTipText("Check UCSC formats");
// this.addActionListener(actionListener);
// }
//
// }

// private class Compute implements Runnable {
// private volatile boolean isRunning;
//
// private Project proj;
// private Segment[] segs;
// private int transformationType;
// private int scope;
// private Logger log;
//
// public Compute(Project proj, Segment[] segs, int transformationType, int scope, Logger log) {
// this.proj = proj;
// this.log = log;
// this.segs = segs;
// this.transformationType = transformationType;
// this.scope = scope;
// this.isRunning = true;
// }
//
// public Compute() {
// this.isRunning = false;
// }
//
// public void kill() {
// log.report("Interupting current Thread...");
// while (computethread.isAlive()) {
// computethread.interrupt();
// Thread.currentThread().interrupt();
// isRunning = false;
// }
// log.report("Thread killed :(");
// }
//
// public boolean isRunning() {
// return isRunning;
// }
//
// public void run() {
//
// while (isRunning && !Thread.currentThread().isInterrupted()) {
// if (transformationType == 0) {
// MedianLRR.createFilesFromMarkerData(proj, segs, outputBase, log);
// MedianLRR.dump(proj, CLASSES_TO_DUMP, proj.getProjectDir() + outputBase + ".mlrr", null, 0, log);
// isRunning = false;
//
// }
// }
// }
// }

// class ComputeWorker extends SwingWorker<Void, String> {
// String computeOutputBase = outputBase;
// Logger computelog = new Logger(proj.getProjectDir() + outputBase + ".log");
// Project computeProject = proj;
// ProgressBarDialog progressBarDialog;
//
// public ComputeWorker(ProgressBarDialog progressBarDialog) {
// this.progressBarDialog = progressBarDialog;
// }
//
// protected Void doInBackground() throws Exception {
// computeComplete = 0;
// // public ProgressBarDialog(String frameText, int min, int max, int width, int height) {
//
// MedianLRRWorker.createFilesFromMarkerData(computeProject, segs, computeOutputBase, computelog);
// Thread.sleep(100);
// publish("Step 2");
// MedianLRRWorker.dump(computeProject, CLASSES_TO_DUMP, computeProject.getProjectDir() + computeOutputBase + ".mlrr", null, 0, computelog);
// Thread.sleep(1000);
// computeComplete = 42;
// return null;
// }
//
// protected void process(ArrayList<String> chunks) {
// jlabel.setText(chunks.get(chunks.size() - 1)); // The last value in this array is all we care about.
// System.out.println(chunks.get(chunks.size() - 1));
// }
//
// protected void done() {
// try {
// get();
// JOptionPane.showMessageDialog(null, "Log R Ratio Summarization Complete");
// } catch (ExecutionException e) {
// computelog.reportError("Error - could not compute median Log R ratio values");
// } catch (InterruptedException e) {
// computelog.reportError("Error - Median Log R ratio computation was interupted ");
// }
// }
// }

// private boolean validateRegions() {
// String notValid = "\n";
// String[] regions = regionTextField.getText().split("\n");
// boolean allvalid = true;
// if (regions == null) {
// JOptionPane.showMessageDialog(this, "No regions Present");
// }
// int numBad = 0;
// for (int i = 0; i < regions.length; i++) {
// if (!checkRegions(regions[i])) {
// notValid += regions[i] + "\n";
// allvalid = false;
// numBad++;
// }
// }
// if (!allvalid) {
// JOptionPane.showMessageDialog(this, (numBad > 1 ? "Invalid formats: " : "Invalid format: ") + notValid);
// }
// return allvalid;
// }
//
// private boolean checkRegions(String region) {
// boolean valid = false;
// int[] newLocation = Positions.parseUCSClocation(region);
// if ((newLocation == null || newLocation.length != 3 || newLocation[0] < 0) || (newLocation[1] < 0) || (newLocation[2] < 0)) {
// valid = false;
// } else {
// valid = true;
// }
// return valid;
// }
