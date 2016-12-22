package org.genvisis.cnv.gui;

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
import java.util.Arrays;
import java.util.concurrent.ExecutionException;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.border.Border;

import org.genvisis.cnv.analysis.MedianLRRWorker;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.SEX_CHROMOSOME_STRATEGY;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.Transforms;
import org.genvisis.cnv.plots.TwoDPlot;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class LRRComp extends JFrame implements Runnable {
	private static final long serialVersionUID = 1L;
	// once a job has been started, used to track completion
	private volatile int computeComplete = 42;
	public static final String FILENAME = "Enter Analysis Name";

	public static final String[] CLASSES_TO_DUMP = {"IID"};
	public static final String[] BASIC_CORRECTION = {"None", "Recompute LRR"};
	public static final String[] BASIC_CORRECTION_TIPS = {"No Correction Procedure will be performed, select this option if you wish to transform the data",
																												"Recompute Log R Ratios"};

	public static final String[] EXTRA_CORRECTION = {"Correct LRR", "Correct XY"};
	public static final String[] REGION_TEXT_FIELD_LABELS =
																												{	"Input UCSC or probeset-based regions of Interest (one per Line):",
																													"Progress...", "Enter Analysis Name Here",
																													"Transform by: ",
																													"Select a Log R Ratio transformation: ",
																													"Select a correction method",
																													"Select for homozygous markers only",
																													"If " + Array.toStr(EXTRA_CORRECTION,
																																							", or ") + " are selected, choose sex-specific correction strategy for chrX (if present)"};

	private int transformationType;
	private int scope;
	private SEX_CHROMOSOME_STRATEGY strategy = SEX_CHROMOSOME_STRATEGY.ARTIFICIAL;
	private final String initRegion;
	private MedianLRRWorker medianLRRWorker;
	private final Project proj;
	private String outputBase;
	private final boolean[] correctionParams;
	private JCheckBox homozygousCheckBox;// this is fairly specific for mitochondrial markers

	public LRRComp(Project proj, String initRegion) {
		this.proj = proj;
		this.initRegion = initRegion;
		transformationType = 0;
		scope = 0;
		correctionParams = new boolean[BASIC_CORRECTION.length + EXTRA_CORRECTION.length];
		Arrays.fill(correctionParams, false);
		correctionParams[0] = true;// defaults to no correction
	}

	@Override
	public void run() {
		createAndShowGUI();
	}

	public void createAndShowGUI() {
		setSize(500, 500);
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		setTitle("Median Log R Ratio Settings");
		WindowListener exitListener = new WindowAdapter() {
			// check for running job before exit
			@Override
			public void windowClosing(WindowEvent e) {
				String ObjButtons[] = {"Exit Anyway", "Cancel"};
				if (computeComplete != 42 && !medianLRRWorker.isDone()) {
					int promptResult = JOptionPane.showOptionDialog(null,
																													"A Thread is computing Median Log R Ratios\nThread will continue if you exit",
																													"Warning - Running thread",
																													JOptionPane.DEFAULT_OPTION,
																													JOptionPane.WARNING_MESSAGE, null,
																													ObjButtons, ObjButtons[1]);
					System.out.println(promptResult);
					if (promptResult == 0) {
						dispose();
					}
				} else {
					dispose();
				}
			}
		};
		addWindowListener(exitListener);
		TransformationPanel transformationPanel = new TransformationPanel();
		transformationPanel.setLayout(new BoxLayout(transformationPanel, BoxLayout.Y_AXIS));
		add(transformationPanel);
		this.setLocation(dim.width / 2	- this.getSize().width / 2,
											dim.height / 2 - this.getSize().height / 2);
		setVisible(true);
	}

	private class TransformationPanel extends JPanel implements ActionListener {
		private static final long serialVersionUID = 1L;
		private final RegionTextField regionTextField;
		private final ComputeButton computeButton;
		private final TwoDPlotButton twoDPlotButton;
		private final JProgressBar progressBar;
		private final FileInputArea fileInputArea;

		private TransformationPanel() {
			outputBase = Transforms.TRANFORMATIONS[transformationType];
			setLayout(new BorderLayout());
			regionTextField = new RegionTextField(initRegion, 10, 10);
			progressBar = new JProgressBar(0, 100);
			computeButton = new ComputeButton(this);
			twoDPlotButton = new TwoDPlotButton(this);
			ActionListener actionListener = getradioListener();

			homozygousCheckBox = new JCheckBox(REGION_TEXT_FIELD_LABELS[6]);
			homozygousCheckBox.setSelected(false);
			homozygousCheckBox.setVisible(true);
			homozygousCheckBox.setToolTipText("Only Homozygous calls will be included in the median value computation (does not change other metrics)");
			homozygousCheckBox.addActionListener(actionListener);

			addLabel(REGION_TEXT_FIELD_LABELS[5]);
			addCorrectionButtons(actionListener, 0);

			addLabel(REGION_TEXT_FIELD_LABELS[4]);
			addTransformButtons(actionListener, transformationType);
			addLabel(REGION_TEXT_FIELD_LABELS[3]);
			addScopeButtons(actionListener, scope);
			addLabel(REGION_TEXT_FIELD_LABELS[2]);
			fileInputArea = new FileInputArea(initRegion, 10, this);
			add(fileInputArea);
			addLabel(REGION_TEXT_FIELD_LABELS[7]);
			addSexStrategyButtons(actionListener, 1);
			addLabel(REGION_TEXT_FIELD_LABELS[0]);
			add(regionTextField, BorderLayout.CENTER);
			JScrollPane scroll = new JScrollPane(regionTextField);
			add(scroll);
			add(homozygousCheckBox);
			add(computeButton, BorderLayout.EAST);
			add(twoDPlotButton, BorderLayout.WEST);
		}

		@Override
		public void actionPerformed(ActionEvent actionEvent) {
			JComponent source = (JComponent) actionEvent.getSource();
			if (computeComplete != 42 && !medianLRRWorker.isDone()) {
				JOptionPane.showMessageDialog(this, "Thread is busy computing median Log R Ratios");
			} else if (source.equals(computeButton)) {
				if (validateFileName()) {
					startJob();
				}
			} else if (source.equals(twoDPlotButton)) {
				if (computeComplete == 42) {
					JOptionPane.showMessageDialog(this, "Please compute Median values before visualizing");
					resetOutputBase();
				} else {
					String fileNameToVisualize;
					revalidate();
					try {
						fileNameToVisualize = medianLRRWorker.get();
						TwoDPlot twoDplot = TwoDPlot.createGUI(proj, true, true);
						twoDplot.showSpecificFile(fileNameToVisualize, 2, 3);
						twoDplot.updateGUI();
						// twoDPlot;

					} catch (InterruptedException e) {
						JOptionPane.showMessageDialog(this,
																					"Thread was interupted when computing Median Log R Ratios");
						e.printStackTrace();
					} catch (ExecutionException e) {
						JOptionPane.showMessageDialog(this, "There was an error computing Median Log R Ratios");
						e.printStackTrace();
					}
					resetOutputBase();
				}
			}

		}

		private boolean validateFileName() {
			boolean valid = true;
			String customName = fileInputArea.getText();
			customName = ext.replaceWithLinuxSafeCharacters(customName, true);
			outputBase = customName + "_" + outputBase;
			outputBase = ext.replaceWithLinuxSafeCharacters(outputBase, true);
			if (!MedianLRRWorker.checkExists(proj, outputBase)) {
				valid = true;
			} else {
				valid = false;
				String ObjButtons[] = {"Overwrite", "Cancel"};
				int promptResult = JOptionPane.showOptionDialog(this,
																												"The Files for Analysis "	+ outputBase
																															+ " Exist",
																												"Warning - Analysis Files Exist",
																												JOptionPane.DEFAULT_OPTION,
																												JOptionPane.WARNING_MESSAGE, null,
																												ObjButtons, ObjButtons[1]);
				if (promptResult == 0) {
					valid = true;
				} else {
					outputBase = outputBase.replaceFirst(customName, "");
				}
			}
			return valid;
		}

		private void startJob() {
			String ObjButtons[] = {"OK", "Cancel"};
			int promptResult = JOptionPane.showOptionDialog(this,
																											"Compute median using "
																															+ Transforms.TRANFORMATIONS[transformationType]
																														+ "?",
																											"Log R Ratio " + Transforms.TRANFORMATIONS[transformationType],
																											JOptionPane.DEFAULT_OPTION,
																											JOptionPane.WARNING_MESSAGE, null, ObjButtons,
																											ObjButtons[1]);
			if (promptResult == 0) {
				add(progressBar);
				progressBar.setVisible(true);
				progressBar.setStringPainted(true);
				computeComplete = 0;
				medianLRRWorker = new MedianLRRWorker(proj, regionTextField.getText().split("\n"),
																							transformationType, scope, outputBase, progressBar,
																							correctionParams[1], correctionParams[2],
																							correctionParams[3], homozygousCheckBox.isSelected(),
																							strategy, proj.getLog());// TODO homozygous
																															// box
				medianLRRWorker.execute();
				revalidate();
			}
			resetOutputBase();
		}

		private void resetOutputBase() {
			outputBase = outputBase.replaceFirst(	ext.replaceWithLinuxSafeCharacters(fileInputArea.getText()
																																							+ "_", true),
																						"");
		}

		private JLabel addLabel(String text) {
			JLabel label = new JLabel(text);
			label.setFont(new Font("Arial", 0, 14));
			add(label);
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
		ButtonGroup typeRadio = new ButtonGroup();

		private void addTransformButtons(ActionListener actionListener, int initScope) {
			JRadioButton[] transformationRadioButtons =
																								new JRadioButton[Transforms.TRANFORMATIONS.length];
			for (int i = 0; i < Transforms.TRANFORMATIONS.length; i++) {
				transformationRadioButtons[i] = new JRadioButton(Transforms.TRANFORMATIONS[i], false);
				transformationRadioButtons[i].setFont(new Font("Arial", 0, 14));
				typeRadio.add(transformationRadioButtons[i]);
				transformationRadioButtons[i].addActionListener(actionListener);
				this.add(transformationRadioButtons[i]);
			}
			transformationRadioButtons[initScope].setSelected(true);
		}

		ButtonGroup sexStrategy = new ButtonGroup();

		private void addSexStrategyButtons(ActionListener actionListener, int initScope) {
			JRadioButton[] sexStrategyButtons = new JRadioButton[SEX_CHROMOSOME_STRATEGY.values().length];
			for (int i = 0; i < SEX_CHROMOSOME_STRATEGY.values().length; i++) {
				sexStrategyButtons[i] = new JRadioButton(	SEX_CHROMOSOME_STRATEGY.values()[i].toString(),
																									false);
				sexStrategyButtons[i].setFont(new Font("Arial", 0, 14));
				sexStrategyButtons[i].setToolTipText(SEX_CHROMOSOME_STRATEGY.values()[i].getToolTip());
				typeRadio.add(sexStrategyButtons[i]);
				sexStrategyButtons[i].addActionListener(actionListener);
				this.add(sexStrategyButtons[i]);
			}
			sexStrategyButtons[initScope].setSelected(true);
		}



		private void addCorrectionButtons(ActionListener actionListener, int initCorrection) {
			ButtonGroup correctionRadio = new ButtonGroup();
			boolean extra = false;
			int numButtons = BASIC_CORRECTION.length;
			// if (Files.exists(proj.getFilename(proj.INTENSITY_PC_FILENAME))) {
			if (Files.exists(proj.INTENSITY_PC_FILENAME.getValue())) {
				extra = true;
				numButtons = BASIC_CORRECTION.length + EXTRA_CORRECTION.length;
			}
			JRadioButton[] correctionRadioButtons = new JRadioButton[numButtons];
			for (int i = 0; i < BASIC_CORRECTION.length; i++) {
				correctionRadioButtons[i] = new JRadioButton(BASIC_CORRECTION[i], false);
				correctionRadioButtons[i].setFont(new Font("Arial", 0, 14));
				correctionRadio.add(correctionRadioButtons[i]);
				correctionRadioButtons[i].addActionListener(actionListener);
				this.add(correctionRadioButtons[i]);
			}
			if (extra) {
				for (int i = BASIC_CORRECTION.length; i < numButtons; i++) {
					correctionRadioButtons[i] =
																		new JRadioButton(	EXTRA_CORRECTION[i - BASIC_CORRECTION.length],
																											false);
					correctionRadioButtons[i].setFont(new Font("Arial", 0, 14));
					correctionRadio.add(correctionRadioButtons[i]);
					correctionRadioButtons[i].addActionListener(actionListener);
					this.add(correctionRadioButtons[i]);
				}
			}
			correctionRadioButtons[initCorrection].setSelected(true);
		}
	}

	// input area
	private class RegionTextField extends JTextArea {
		private static final long serialVersionUID = 1L;

		private RegionTextField(String region, int width, int height) {
			setFont(new Font("Tahoma", Font.PLAIN, 14));
			setText(region);
			setSize(width, height);
			Border border = BorderFactory.createLineBorder(Color.BLACK);
			setBorder(border);
			setMaximumSize(getPreferredSize());

		}
	}

	private class FileInputArea extends JTextField {
		private static final long serialVersionUID = 1L;

		private FileInputArea(String init, int width, ActionListener actionListener) {
			setFont(new Font("Tahoma", Font.PLAIN, 14));
			setText(init);
			Border border = BorderFactory.createLineBorder(Color.BLACK);
			setBorder(border);
			setMaximumSize(new Dimension(Integer.MAX_VALUE, getMinimumSize().height));
			// setMaximumSize(this.getPreferredSize());
			addActionListener(actionListener);
			setVisible(true);
		}
	}

	private class ComputeButton extends JButton {
		private static final long serialVersionUID = 1L;

		private ComputeButton(ActionListener actionListener) {
			setFont(new Font("Tahoma", Font.PLAIN, 14));
			setText("Compute");
			setToolTipText("Compute Median Log R Ratios");
			addActionListener(actionListener);
		}
	}

	private class TwoDPlotButton extends JButton {
		private static final long serialVersionUID = 1L;

		private TwoDPlotButton(ActionListener actionListener) {
			setFont(new Font("Tahoma", Font.PLAIN, 14));
			setText("2D Plot");
			setToolTipText("Launch 2D plot to visualize results");
			addActionListener(actionListener);
		}
	}

	public ActionListener getradioListener() {
		return new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent actionEvent) {

				if (ext.indexOfStr(	actionEvent.getActionCommand(), Transforms.TRANFORMATIONS, true, true,
														proj.getLog(), false) >= 0) {
					transformationType = ext.indexOfStr(actionEvent.getActionCommand(),
																							Transforms.TRANFORMATIONS, true, true, proj.getLog(),
																							false);
					outputBase = Transforms.TRANFORMATIONS[transformationType];
				} else if (ext.indexOfStr(actionEvent.getActionCommand(), Transforms.SCOPES, true, true,
																	proj.getLog(), false) >= 0
										&& transformationType == 0) {
					JOptionPane.showMessageDialog(null,
																				"Transform by Chromosome or Genome not valid for Raw Values");
				} else if (ext.indexOfStr(actionEvent.getActionCommand(), Transforms.SCOPES, true, true,
																	proj.getLog(), false) >= 0) {
					scope = ext.indexOfStr(	actionEvent.getActionCommand(), Transforms.SCOPES, true, true,
																	proj.getLog(), false);
					outputBase = Transforms.TRANFORMATIONS[transformationType]	+ "_"
												+ Transforms.SCOPES[scope];
				} else if (ext.indexOfStr(actionEvent.getActionCommand(), BASIC_CORRECTION, true, true,
																	proj.getLog(), false) >= 0) {
					int index = ext.indexOfStr(	actionEvent.getActionCommand(), BASIC_CORRECTION, true, true,
																			proj.getLog(), false);
					Arrays.fill(correctionParams, false);
					correctionParams[index] = true;
					outputBase = BASIC_CORRECTION[index];

				} else if (ext.indexOfStr(actionEvent.getActionCommand(), EXTRA_CORRECTION, true, true,
																	proj.getLog(), false) >= 0) {
					int index = ext.indexOfStr(	actionEvent.getActionCommand(), EXTRA_CORRECTION, true, true,
																			proj.getLog(), false);
					Arrays.fill(correctionParams, false);
					correctionParams[(BASIC_CORRECTION.length + index)] = true;
					outputBase = EXTRA_CORRECTION[index];
				} else {
					try {
					strategy = SEX_CHROMOSOME_STRATEGY.valueOf(actionEvent.getActionCommand());
					} catch (IllegalArgumentException ile) {

					}
				}
				if (transformationType != 0 && !correctionParams[0]) {
					JOptionPane.showMessageDialog(null,
																				"Intensity Correction is currently not valid for Transformed Values");
				}
				if (transformationType != 0 && homozygousCheckBox.isSelected()) {
					JOptionPane.showMessageDialog(null,
																				"Selecting homozygous markers is currently not valid for Transformed Values");
					homozygousCheckBox.setSelected(false);
				}
			}
		};
	}
}
