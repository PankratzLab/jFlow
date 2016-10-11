package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import org.genvisis.cnv.Launch;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.GenvisisWorkflow;
import org.genvisis.cnv.manage.GenvisisWorkflow.RequirementInputType;
import org.genvisis.cnv.manage.GenvisisWorkflow.STEP;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Grafik;
import org.genvisis.common.ext;

import net.miginfocom.swing.MigLayout;

public class GenvisisWorkflowGUI extends JDialog {

	private static final long serialVersionUID = 1L;

	private final JPanel contentPanel = new JPanel();


	private final ConcurrentHashMap<STEP, JCheckBox> checkBoxes =
																															new ConcurrentHashMap<STEP, JCheckBox>();
	private final ConcurrentHashMap<STEP, JLabel> descLabels = new ConcurrentHashMap<STEP, JLabel>();
	private final ConcurrentHashMap<STEP, ArrayList<JLabel>> requirementsLabels =
																																							new ConcurrentHashMap<STEP, ArrayList<JLabel>>();
	private final ConcurrentHashMap<STEP, JAccordionPanel> panels =
																																new ConcurrentHashMap<GenvisisWorkflow.STEP, JAccordionPanel>();
	public ConcurrentHashMap<STEP, ArrayList<? extends JComponent>> varFields =
																																						new ConcurrentHashMap<GenvisisWorkflow.STEP, ArrayList<? extends JComponent>>();
	public ConcurrentHashMap<STEP, JProgressBar> progBars =
																												new ConcurrentHashMap<GenvisisWorkflow.STEP, JProgressBar>();
	public ConcurrentHashMap<STEP, ArrayList<JButton>> fileBtns =
																															new ConcurrentHashMap<GenvisisWorkflow.STEP, ArrayList<JButton>>();
	public ConcurrentHashMap<STEP, JLabel> alreadyRunLbls =
																												new ConcurrentHashMap<GenvisisWorkflow.STEP, JLabel>();
	public ConcurrentHashMap<STEP, JButton> cancelStepBtns = new ConcurrentHashMap<GenvisisWorkflow.STEP, JButton>();

	Project proj;

	private static final String TOP_LABEL = "Genvisis Project Workflow:";

	volatile boolean cancelled = false;
	volatile boolean[] selected;
	STEP[] steps;
	int DEFAULT_SCROLL_SPEED = 16;

	/**
	 * Create the dialog.
	 */
	public GenvisisWorkflowGUI(Project proj2, final Launch launch) {
		if (proj2 == null) {
			proj = createNewProject(launch.getLaunchProperties().getListOfProjectNames());
		} else {
			proj = proj2;
		}
		if (proj == null) {
			doClose();
			return;
		} else {
			launch.loadProjects();
			launch.setIndexOfCurrentProject(ext.removeDirectoryInfo(proj.getPropertyFilename()));
			proj = launch.loadProject();
		}
		proj.getLog().report("Launching Genvisis Project Pipeline");
		steps = GenvisisWorkflow.getStepsForProject(proj);
		selected = Array.booleanArray(steps.length, true);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		JPanel optionPanel = new JPanel();
		JScrollPane scrollPane = new JScrollPane(optionPanel);
		scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		scrollPane.getVerticalScrollBar().setUnitIncrement(DEFAULT_SCROLL_SPEED);
		contentPanel.setLayout(new MigLayout("", "[grow]", "[][][]"));
		contentPanel.add(scrollPane, "cell 0 2, grow");
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		optionPanel.setLayout(new MigLayout("", "[grow]", "[][][]"));
		{
			JLabel lblKitAndKaboodle = new JLabel(TOP_LABEL);
			lblKitAndKaboodle.setFont(new Font("Arial", Font.BOLD, 16));
			contentPanel.add(lblKitAndKaboodle, "cell 0 0,alignx center");
		}
		{
			for (int i = 0; i < steps.length; i++) {
				final int index = i;
				JAccordionPanel panel = createPanel(index);
				optionPanel.add(panel, "cell 0 " + (i) + ", alignx left, growx");
				panels.put(steps[i], panel);
			}
		}
		Insets btnInsets = new Insets(0, 1, 0, 1);
		{
			JPanel buttonPane = new JPanel();
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			buttonPane.setLayout(new MigLayout("hidemode 2", "[][][][][][][][grow][][][]", "[]"));

			JLabel lblSelect = new JLabel("Select:");
			buttonPane.add(lblSelect, "flowx,cell 0 0");

			btnSelectAll = new JButton("All");
			btnSelectAll.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					for (Entry<STEP, JCheckBox> entry : checkBoxes.entrySet()) {
						entry.getValue().setSelected(true);
					}
					for (int i = 0; i < selected.length; i++) {
						selected[i] = true;
					}
					refreshLabels(GenvisisWorkflowGUI.this, Arrays.asList(steps));
				}
			});
			btnSelectAll.setMargin(btnInsets);
			buttonPane.add(btnSelectAll, "cell 0 0");

			btnDeselectAll = new JButton("None");
			btnDeselectAll.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					for (Entry<STEP, JCheckBox> entry : checkBoxes.entrySet()) {
						entry.getValue().setSelected(false);
					}
					for (int i = 0; i < selected.length; i++) {
						selected[i] = false;
					}
					refreshLabels(GenvisisWorkflowGUI.this, Arrays.asList(steps));
				}
			});
			btnDeselectAll.setMargin(btnInsets);
			buttonPane.add(btnDeselectAll, "cell 1 0");

			btnSelectValid = new JButton("Valid");
			btnSelectValid.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					int stepIndex = -1;
					HashMap<STEP, Boolean> selectedSteps = new HashMap<GenvisisWorkflow.STEP, Boolean>();
					for (Entry<STEP, JCheckBox> entry : checkBoxes.entrySet()) {
						selectedSteps.put(entry.getKey(), true); // pretend everything is selected
					}
					HashMap<STEP, ArrayList<String>> variables = getVariables();
					for (final STEP step : steps) {
						stepIndex++;
						if (step == null || checkBoxes.get(step) == null || varFields.get(step) == null) {
							continue;
						}
						if (!step.checkIfOutputExists(proj, variables)) {
							boolean check = step.hasRequirements(proj, selectedSteps, variables);
							checkBoxes.get(step).setSelected(check);
							selected[stepIndex] = check;
							selectedSteps.put(step, check);
						} else {
							selectedSteps.put(step, false);
							selected[stepIndex] = false;
							checkBoxes.get(step).setSelected(false);
						}
					}
					refreshLabels(GenvisisWorkflowGUI.this, Arrays.asList(steps));
				}
			});
			btnSelectValid.setMargin(btnInsets);
			buttonPane.add(btnSelectValid, "cell 2 0");

			JSeparator separator = new JSeparator();
			separator.setOrientation(SwingConstants.VERTICAL);
			buttonPane.add(separator, "cell 3 0,growy");

			JLabel lblCollapse = new JLabel("Collapse:");
			buttonPane.add(lblCollapse, "cell 4 0");

			JButton btnAll = new JButton("All");
			btnAll.setMargin(btnInsets);
			btnAll.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					for (Entry<STEP, JAccordionPanel> entry : panels.entrySet()) {
						entry.getValue().shrink();
					}
				}
			});
			buttonPane.add(btnAll, "cell 5 0");

			JButton btnNone = new JButton("None");
			btnNone.setMargin(btnInsets);
			btnNone.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					for (Entry<STEP, JAccordionPanel> entry : panels.entrySet()) {
						entry.getValue().expand();
					}
				}
			});
			buttonPane.add(btnNone, "cell 6 0");

			progVal = new JProgressBar(0, steps.length);
			progVal.setValue(0);
			progVal.setStringPainted(true);
			progVal.setString("Validating Steps...");
			progVal.setVisible(false);
			buttonPane.add(progVal, "growx, cell 7 0");

			AbstractAction listener = new AbstractAction() {
				private static final long serialVersionUID = 1L;

				@Override
				public void actionPerformed(ActionEvent e) {
					if (running) {
						return;
					}
					if (e.getActionCommand().equals("Close")) {
						cancelled = true;
						doClose();
					} else if (e.getActionCommand().equals("Export")) {
						exportToText();
					} else {
						run();
					}
				}
			};
			btnExportToText = new JButton("Export To Text");
			btnExportToText.setActionCommand("Export");
			btnExportToText.addActionListener(listener);
			btnExportToText.setMargin(btnInsets);
			buttonPane.add(btnExportToText, "cell 8 0, alignx right");
			btnOk = new JButton("Run");
			btnOk.setMargin(btnInsets);
			btnOk.setActionCommand("Run");
			btnOk.setMnemonic(KeyEvent.VK_O);
			buttonPane.add(btnOk, "cell 9 0,alignx right");
			getRootPane().setDefaultButton(btnOk);
			btnCancel = new JButton("Close");
			btnCancel.setMargin(btnInsets);
			btnCancel.setActionCommand("Close");
			btnCancel.setMnemonic(KeyEvent.VK_C);
			buttonPane.add(btnCancel, "cell 10 0,alignx left");

			btnOk.addActionListener(listener);
			btnCancel.addActionListener(listener);

		}
		pack();
		setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				if (running) {
					return;
				}
				doClose();
			}
		});
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				HashMap<STEP, ArrayList<String>> variables = getVariables();
				for (int i = 0; i < steps.length; i++) {
					STEP step = steps[i];
					if (step.checkIfOutputExists(proj, variables)) {
						checkBoxes.get(step).setSelected(false);
						alreadyRunLbls.get(step).setVisible(true);
						panels.get(step).shrink();
						selected[i] = false;
					}
				}
			}
		});
		refreshLabels(this, Arrays.asList(steps));
		setBounds(100, 100, 750, 850);
		setTitle(TOP_LABEL);
		addWindowFocusListener(new WindowFocusListener() {
			@Override
			public void windowLostFocus(WindowEvent e) {}

			@Override
			public void windowGainedFocus(WindowEvent e) {
				launch.toFront();
				GenvisisWorkflowGUI.this.toFront();
			}
		});

		InputMap inMap = getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		ActionMap actMap = getRootPane().getActionMap();

		KeyStroke escKey = KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0);
		AbstractAction escapeAction = new AbstractAction() {
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent e) {
				doClose();
			}
		};
		inMap.put(escKey, "Action.escape");
		actMap.put("Action.escape", escapeAction);

		// panels start out visible to help with spacing (otherwise the containing jscrollpane is too
		// small)
		// for (JPanel panel : panels.values()) {
		// ((JAccordionPanel) panel).shrink();
		// }
	}

	protected void doClose() {
		// TODO check if running, may not be able to stop current step. Ask if should interrupt. Maybe
		// do cleanup if desired?
		cancelled = true;
		setVisible(false);
		dispose();
	}

	AbstractAction interruptAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;

		@Override
		public void actionPerformed(ActionEvent arg0) {
			if (runThread != null) {
				runThread.interrupt();
			}
		}
	};

	AbstractAction fileSelectAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;

		@Override
		public void actionPerformed(ActionEvent e) {
			String[] parts = e.getActionCommand().split(":");
			int stepIndex = Integer.parseInt(parts[0]);
			int fieldIndex = Integer.parseInt(parts[1]);
			JTextField fileField = (JTextField) varFields.get(steps[stepIndex]).get(fieldIndex);

			String current = fileField.getText();

			String dir = current.equals("")	? proj.PROJECT_DIRECTORY.getValue(false, false)
																			: ext.parseDirectoryOfFile(current);
			JFileChooser chooser = new JFileChooser(dir);
			chooser.setMultiSelectionEnabled(false);
			final STEP step = steps[stepIndex];
			RequirementInputType[][] reqs = step.reqTypes;
			int accum = 0;
			RequirementInputType type = null;
			outer: for (int i = 0; i < reqs.length; i++) {
				for (int j = 0; j < reqs[i].length; j++) {
					if (reqs[i][j] == RequirementInputType.NONE) {
						continue;
					}
					if (accum == fieldIndex) {
						type = reqs[i][j];
						break outer;
					}
					accum++;
				}
			}

			if (type == RequirementInputType.FILE) {
				chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
				chooser.setDialogTitle("Select File");
			} else if (type == RequirementInputType.DIR) {
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				chooser.setDialogTitle("Select Directory");
			}

			int retValue = chooser.showDialog(GenvisisWorkflowGUI.this, "Select");

			if (retValue == JFileChooser.CANCEL_OPTION) {
				refreshLabels(GenvisisWorkflowGUI.this, step.getRelatedSteps());
				return;
			} else {
				File newFile = chooser.getSelectedFile();
				String txt = ext.verifyDirFormat(newFile.getAbsolutePath());
				if (type == RequirementInputType.FILE) {
					txt = txt.substring(0, txt.length() - 1);
				}
				fileField.setText(txt);
			}
			refreshLabels(GenvisisWorkflowGUI.this, step.getRelatedSteps());
		}
	};

	private JAccordionPanel createPanel(final int index) {
		final STEP step = steps[index];
		final JAccordionPanel panel = new JAccordionPanel();

		final JCheckBox chckbx = new JCheckBox();
		chckbx.setAction(new StepRefresher(GenvisisWorkflowGUI.this, step) {
			/**
			*
			*/
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent e) {
				selected[index] = chckbx.isSelected();
				super.actionPerformed(e);
			}
		});
		chckbx.setFont(chckbx.getFont().deriveFont(Font.PLAIN, 14));
		Grafik.scaleCheckBoxIcon(chckbx);
		chckbx.setSelected(selected[index]);
		chckbx.setToolTipText(step.stepDesc);
		chckbx.setText((index + 1) + ": " + step.stepName);

		checkBoxes.put(step, chckbx);
		panel.topPanel.add(chckbx, "cell 0 0");

		JProgressBar stepProgBar = new JProgressBar();
		progBars.put(step, stepProgBar);
		stepProgBar.setVisible(false);
		panel.topPanel.add(stepProgBar, "cell 1 0, alignx right, hidemode 3, split 2");
		JLabel alreadyRanLbl = new JLabel("Output Already Exists!");
		alreadyRanLbl.setVisible(false);
		alreadyRunLbls.put(step, alreadyRanLbl);
		panel.topPanel.add(alreadyRanLbl, "cell 1 0, alignx right, hidemode 3");
		JButton cancelStepButton = new JButton(interruptAction);
		cancelStepButton.setIcon(Grafik.getImageIcon("images/redx.png"));
		cancelStepButton.setVisible(false);
		cancelStepButton.setIconTextGap(0);
		cancelStepButton.setMargin(new Insets(1, 1, 0, 1));
		cancelStepBtns.put(step, cancelStepButton);
		panel.topPanel.add(cancelStepButton, "cell 1 0, alignx right, hidemode 3");

		JLabel descLbl = new JLabel("<html><center><p>" + step.stepDesc + "</p></center></html>");
		descLbl.setVerticalAlignment(SwingConstants.TOP);
		descLbl.setHorizontalAlignment(SwingConstants.LEFT);
		descLbl.setFont(descLbl.getFont().deriveFont(Font.PLAIN));
		descLabels.put(step, descLbl);

		panel.setBorder(new LineBorder(Color.GRAY.brighter(), 1, true));
		String rows = "[][]";
		for (int i = 0; i < step.getRequirements().length; i++) {
			rows = rows + "[]";
		}
		panel.contentPanel.setLayout(new MigLayout("", "[200px,grow]push[200px,grow]", rows));
		panel.contentPanel.add(descLbl, "cell 0 0");

		String[][] reqs = step.getRequirements();

		String reqSteps = "abcdefghijklmnop"; // probably won't need that many
		if (reqs.length > 0) {
			JLabel reqLbl = new JLabel("Requires:");
			panel.contentPanel.add(reqLbl, "cell 0 1");

			ArrayList<JLabel> reqLbls = new ArrayList<JLabel>();
			int rowIndex = 2;
			for (int i = 0; i < reqs.length; i++) {
				// AND

				for (int j = 0; j < reqs[i].length; j++) {
					// OR
					JLabel indLbl = new JLabel(""	+ (i + 1) + (reqs[i].length > 1 ? reqSteps.charAt(j) : "")
																			+ ". ");
					indLbl.setFont(indLbl.getFont().deriveFont(Font.PLAIN, 9));
					panel.contentPanel.add(indLbl, "gapleft 25, aligny top, split 1, cell 0 " + rowIndex);
					JLabel requirementLbl = new JLabel("<html><p>" + reqs[i][j] + "</p></html>");
					requirementLbl.setFont(requirementLbl.getFont().deriveFont(Font.PLAIN, 9));
					panel.contentPanel.add(requirementLbl, "aligny top, cell 0 " + rowIndex);
					reqLbls.add(requirementLbl);
					rowIndex++;
					if (j < reqs[i].length - 1) {
						JLabel orLbl = new JLabel("OR");
						orLbl.setFont(orLbl.getFont().deriveFont(Font.PLAIN, 10));
						panel.contentPanel.add(orLbl, "gapleft 18, cell 0 " + rowIndex);
						rowIndex++;
					}

				}
				requirementsLabels.put(step, reqLbls);

				if (i < reqs.length - 1) {
					JLabel andLbl = new JLabel("AND");
					andLbl.setFont(andLbl.getFont().deriveFont(Font.PLAIN, 10));
					panel.contentPanel.add(andLbl, "gapleft 7, cell 0 " + rowIndex);
					rowIndex++;
				}

			}

			RequirementInputType[][] inputTypes = step.getRequirementInputTypes();
			ArrayList<JComponent> reqInputFields = new ArrayList<JComponent>();
			int reqIndex = 0;
			rowIndex = 2;
			for (int i = 0; i < inputTypes.length; i++) {
				for (int j = 0; j < inputTypes[i].length; j++) {

					if (inputTypes[i][j] == RequirementInputType.BOOL) {
						JCheckBox checkBox = new JCheckBox();
						checkBox.setAction(new StepRefresher(GenvisisWorkflowGUI.this, step));
						checkBox.setFont(checkBox.getFont().deriveFont(14));
						Grafik.scaleCheckBoxIcon(checkBox);
						checkBox.setVerticalAlignment(SwingConstants.TOP);
						checkBox.setHorizontalAlignment(SwingConstants.RIGHT);
						boolean sel = false;
						sel = Boolean.valueOf(step.getRequirementDefaults(proj)[reqIndex].toString());
						checkBox.setSelected(sel);
						reqIndex++;
						reqInputFields.add(checkBox);
						panel.contentPanel.add(	checkBox,
																		"alignx right, aligny center, growx, gapleft 20, cell 1 "
																							+ rowIndex);
					} else if (inputTypes[i][j] != RequirementInputType.NONE) {
						JTextField textField = new JTextField();
						// textField.setHorizontalAlignment(JTextField.RIGHT);
						// textField.addKeyListener(new KeyAdapter() {
						// @Override
						// public void keyTyped(KeyEvent e) {
						// super.keyTyped(e);
						// refreshLabels();
						// }
						// });
						textField.getDocument().addDocumentListener(new TextChangedListener() {
							@Override
							public void changedUpdate(DocumentEvent e) {
								refreshLabels(GenvisisWorkflowGUI.this, step.getRelatedSteps());
							}
						});
						textField.setText(step.getRequirementDefaults(proj)[reqIndex].toString());
						reqIndex++;
						reqInputFields.add(textField);
						panel.contentPanel.add(	textField,
																		"alignx right, aligny center, growx, gapleft 20, split 1, cell 1 "
																								+ rowIndex);
						if (inputTypes[i][j] == RequirementInputType.FILE
								|| inputTypes[i][j] == RequirementInputType.DIR) {
							JButton fileBtn = new JButton();
							fileBtn.setAction(fileSelectAction);
							fileBtn.setText("...");
							fileBtn.setActionCommand((index) + ":" + (reqInputFields.size() - 1));
							fileBtn.setMargin(new Insets(1, 4, 0, 3));
							panel.contentPanel.add(fileBtn, "cell 1 " + rowIndex);
							ArrayList<JButton> list = fileBtns.get(step);
							if (list == null) {
								list = new ArrayList<JButton>();
								fileBtns.put(step, list);
							}
							list.add(fileBtn);
						}
					}
					rowIndex++;
					if (j < reqs[i].length - 1) {
						rowIndex++;
					}
				}
				rowIndex++;
			}
			varFields.put(step, reqInputFields);

		}

		return panel;

	}

	public boolean[] getSelectedOptions() {
		return selected;
	}

	public HashMap<STEP, Boolean> getSelectedOptionsMap() {
		HashMap<STEP, Boolean> selectedSteps = new HashMap<GenvisisWorkflow.STEP, Boolean>();
		for (Entry<STEP, JCheckBox> entry : checkBoxes.entrySet()) {
			selectedSteps.put(entry.getKey(), entry.getValue().isSelected());
		}
		return selectedSteps;
	}

	public boolean getCancelled() {
		return cancelled;
	}

	public void startStep(STEP step) {
		if (alreadyRunLbls.get(step).isVisible()) {
			alreadyRunLbls.get(step).setVisible(false);
		}
		progBars.get(step).setString("Running...");
		progBars.get(step).setStringPainted(true);
		progBars.get(step).setIndeterminate(true);
		progBars.get(step).setVisible(true);
		cancelStepBtns.get(step).setVisible(true);
		cancelStepBtns.get(step).setEnabled(true);
	}

	public void endStep(STEP step, FINAL_CODE code) {
		String resp = code.getMessage();
		progBars.get(step).setString(resp);
		progBars.get(step).setIndeterminate(false);
		cancelStepBtns.get(step).setEnabled(false);
	}


	/**
	 * Helper {@link DocumentListener} that redirects both {@link #insertUpdate(DocumentEvent)} and
	 * {@link #removeUpdate(DocumentEvent)} to {@link #changedUpdate(DocumentEvent)}.
	 */
	public abstract static class TextChangedListener implements DocumentListener {
		@Override
		public void insertUpdate(DocumentEvent e) {
			changedUpdate(e);
		}

		@Override
		public void removeUpdate(DocumentEvent e) {
			changedUpdate(e);
		}
	}
	/**
	 * Helper {@link ActionListener} to refresh one or more UI steps.
	 */
	public static class StepRefresher extends AbstractAction {

		/**
		*
		*/
		private static final long serialVersionUID = 1L;
		private final transient Set<STEP> stepsToRefresh;
		private final GenvisisWorkflowGUI refrenceGUI;

		/**
		 * @param gui UI which is displaying the given STEPs.
		 * @param steps STEPs to validate by this {@link ActionListener}.
		 */
		public StepRefresher(final GenvisisWorkflowGUI gui, final STEP... steps) {
			refrenceGUI = gui;
			stepsToRefresh = new HashSet<STEP>();
			for (final STEP s : steps) {
				stepsToRefresh.addAll(s.getRelatedSteps());
			}

		}

		@Override
		public void actionPerformed(ActionEvent e) {
			refreshLabels(refrenceGUI, stepsToRefresh);
		}
	}

	/**
	 * Validate all elements of the given {@link STEP}s and refresh the specified UI.
	 *
	 * @param stepsToRefresh
	 */
	public static void refreshLabels(	final GenvisisWorkflowGUI gui,
																		final Collection<STEP> stepsToRefresh) {
		new Thread(new Runnable() {
			@Override
			public void run() {
				try {
					SwingUtilities.invokeAndWait(new Runnable() {
						@Override
						public void run() {
							gui.progVal.setValue(0);
							gui.progVal.setVisible(true);
						}
					});
				} catch (InvocationTargetException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				final Color greenDark = Color.GREEN.darker();
				final Color dark = Color.GRAY;
				HashMap<STEP, Boolean> selectedSteps = new HashMap<GenvisisWorkflow.STEP, Boolean>();
				for (Entry<STEP, JCheckBox> entry : gui.checkBoxes.entrySet()) {
					selectedSteps.put(entry.getKey(), entry.getValue().isSelected());
				}
				int i = 0;
				for (final STEP step : stepsToRefresh) {
					if (step == null || gui.checkBoxes.get(step) == null || gui.varFields.get(step) == null) {
						continue;
					}
					HashMap<STEP, ArrayList<String>> variables = gui.getVariables();
					final int update = ++i;
					if (!step.checkIfOutputExists(gui.proj, variables)
							|| gui.checkBoxes.get(step).isSelected()) {
						boolean check = step.hasRequirements(gui.proj, selectedSteps, variables);
						gui.descLabels.get(step).setForeground(check ? greenDark : Color.RED);
						gui.checkBoxes.get(step).setForeground(check ? greenDark : Color.RED);
						final ArrayList<JLabel> reqLbls = gui.requirementsLabels.get(step);
						final boolean[][] reqVals = step.checkRequirements(gui.proj, selectedSteps, variables);
						// SwingUtilities.invokeLater(new Runnable() {
						// @Override
						// public void run() {
						// int lblIndex = 0;
						// for (int i = 0; i < reqVals.length; i++) {
						// boolean hasAny = Array.booleanArraySum(reqVals[i]) > 0;
						// for (int j = 0; j < reqVals[i].length; j++) {
						// reqLbls.get(lblIndex).setForeground(reqVals[i][j] ? greenDark : hasAny ? Color.GRAY :
						// Color.RED);
						// lblIndex++;
						// }
						// }
						// }
						// });
						try {
							SwingUtilities.invokeAndWait(new Runnable() {
								@Override
								public void run() {
									int lblIndex = 0;
									for (boolean[] reqVal : reqVals) {
										boolean hasAny = Array.booleanArraySum(reqVal) > 0;
										for (boolean element : reqVal) {
											reqLbls	.get(lblIndex)
															.setForeground(element ? greenDark : hasAny ? Color.GRAY : Color.RED);
											lblIndex++;
										}
									}
									gui.progVal.setValue(update);
								}
							});
						} catch (InvocationTargetException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					} else {
						final boolean[][] reqVals = step.checkRequirements(gui.proj, selectedSteps, variables);
						final ArrayList<JLabel> reqLbls = gui.requirementsLabels.get(step);
						// SwingUtilities.invokeLater(new Runnable() {
						// @Override
						// public void run() {
						// int lblIndex = 0;
						// checkBoxes.get(step).setSelected(false);
						// alreadyRunLbls.get(step).setVisible(true);
						// descLabels.get(step).setForeground(dark);
						// checkBoxes.get(step).setForeground(dark);
						// for (int i = 0; i < reqVals.length; i++) {
						// for (int j = 0; j < reqVals[i].length; j++) {
						// reqLbls.get(lblIndex).setForeground(dark);
						// lblIndex++;
						// }
						// }
						// }
						// });
						try {
							SwingUtilities.invokeAndWait(new Runnable() {
								@Override
								public void run() {
									int lblIndex = 0;
									gui.checkBoxes.get(step).setSelected(false);
									gui.alreadyRunLbls.get(step).setVisible(true);
									gui.descLabels.get(step).setForeground(dark);
									gui.checkBoxes.get(step).setForeground(dark);
									for (boolean[] reqVal : reqVals) {
										for (boolean element : reqVal) {
											reqLbls.get(lblIndex).setForeground(dark);
											lblIndex++;
										}
									}
									gui.progVal.setValue(update);
								}
							});
						} catch (InvocationTargetException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
				try {
					SwingUtilities.invokeAndWait(new Runnable() {
						@Override
						public void run() {
							gui.progVal.setValue(0);
							gui.progVal.setVisible(false);
						}
					});
				} catch (InvocationTargetException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

			}
		}).start();

	}

	private volatile boolean running = false;

	private JButton btnSelectAll;
	private JButton btnDeselectAll;
	private JButton btnSelectValid;
	private JProgressBar progVal;

	private Thread runThread;

	private JButton btnExportToText;
	private JButton btnOk;
	private JButton btnCancel;

	private void lockup(final boolean lock) {
		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				@Override
				public void run() {
					for (JCheckBox box : checkBoxes.values()) {
						box.setEnabled(!lock);
					}
					for (ArrayList<? extends JComponent> flds : varFields.values()) {
						for (JComponent fld : flds) {
							fld.setEnabled(!lock);
						}
					}
					for (ArrayList<JButton> btns : fileBtns.values()) {
						for (JButton btn : btns) {
							btn.setEnabled(!lock);
						}
					}
					btnSelectAll.setEnabled(!lock);
					btnDeselectAll.setEnabled(!lock);
					btnSelectValid.setEnabled(!lock);
					btnCancel.setEnabled(!lock);
					btnOk.setEnabled(!lock);
					btnExportToText.setEnabled(!lock);
				}
			});
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	private void exportToText() {
		running = true;
		new Thread(new Runnable() {
			@Override
			public void run() {
				lockup(true);

				boolean[] options = getSelectedOptions();
				HashMap<STEP, Boolean> selectedSteps = getSelectedOptionsMap();
				HashMap<STEP, ArrayList<String>> variables = getVariables();
				if (checkRequirementsAndNotify(selectedSteps, variables)) {
					StringBuilder output = new StringBuilder("## Genvisis Project Pipeline - Stepwise Commands\n\n");
					for (int i = 0; i < options.length; i++) {
						if (options[i]) {
							String cmd = steps[i].getCommandLine(proj, variables);
							output.append("## ").append(steps[i].stepName).append("\n");
							output.append("echo \" start ").append(steps[i].stepName).append(" at: \" `date`")
										.append("\n");
							output.append(cmd).append("\n");
							output.append("echo \" end ").append(steps[i].stepName).append(" at: \" `date`")
										.append("\n");
							output.append("\n\n");
						}
					}
					Files.write(output.toString(),
											proj.PROJECT_DIRECTORY.getValue() + "GenvisisPipeline.run");
					Files.qsub(proj.PROJECT_DIRECTORY.getValue()	+ "GenvisisPipeline."
											+ ext.getTimestampForFilename() + ".pbs", output.toString(), 50000, 48, 24);
					proj.message("GenvisisPipeline commands written to "	+ proj.PROJECT_DIRECTORY.getValue()
												+ "GenvisisPipeline.run", "Command File Written",
												JOptionPane.INFORMATION_MESSAGE);
				}

				lockup(false);
				running = false;
			}
		}).start();
	}

	enum FINAL_CODE {
										COMPLETE("Complete"), FAILED("Failed"), CANCELLED("Cancelled");

		private String message;

		FINAL_CODE(String msg) {
			message = msg;
		}

		public String getMessage() {
			return message;
		}
	}

	private void run() {
		running = true;
		runThread = new Thread(new Runnable() {
			@Override
			public void run() {
				lockup(true);

				boolean[] options = getSelectedOptions();
				HashMap<STEP, Boolean> selectedSteps = getSelectedOptionsMap();
				HashMap<STEP, ArrayList<String>> variables = getVariables();
				if (checkRequirementsAndNotify(selectedSteps, variables)) {
					for (int i = 0; i < options.length; i++) {
						if (options[i]) {
							startStep(steps[i]);
							Throwable e = null;
							FINAL_CODE code = FINAL_CODE.COMPLETE;
							try {
								steps[i].setNecessaryPreRunProperties(proj, variables);
								steps[i].run(proj, variables);
							} catch (Throwable e1) {
								if (e1 instanceof RuntimeException
										&& ((RuntimeException) e1).getCause() instanceof InterruptedException) {
									code = FINAL_CODE.CANCELLED;
								}
								e = e1;
								System.gc();
							}
							if (code != FINAL_CODE.CANCELLED
									&& (e != null	|| steps[i].getFailed()
											|| !steps[i].checkIfOutputExists(proj, variables))) {
								code = FINAL_CODE.FAILED;
							}
							endStep(steps[i], code);
							if (code == FINAL_CODE.FAILED) {
								StringBuilder failureMessage = new StringBuilder("Error Occurred on Step ")	.append(i + 1)
																																														.append(":");
								if (e != null) {
									proj.getLog().reportException(e, 0);
									failureMessage.append("\n").append(e.getMessage());
								} else if (steps[i].getFailed()) {
									for (String msg : steps[i].getFailureMessages()) {
										failureMessage.append("\n").append(msg);
									}
								} else if (!steps[i].checkIfOutputExists(proj, variables)) {
									failureMessage.append("\nUnknown error occurred.");
								}
								failureMessage.append("\nPlease check project log for more details.");
								String[] opts = {"Continue", "Retry", "Cancel"};
								int opt = JOptionPane.showOptionDialog(	GenvisisWorkflowGUI.this,
																												failureMessage.toString(), "Error!",
																												JOptionPane.YES_NO_OPTION,
																												JOptionPane.ERROR_MESSAGE, null, opts,
																												opts[2]);
								if (opt == JOptionPane.CLOSED_OPTION || opt == 2) { // closed or cancel
									for (STEP step : steps) {
										step.resetRun();
									}
									break;
								} else if (opt == 1) { // retry
									steps[i].resetRun();
									i--;
								} else {
									// continue
								}
							} else if (code == FINAL_CODE.CANCELLED) {
								steps[i].gracefulDeath(proj);
								// TODO remove message when gracefulDeath is implemented for each step
								JOptionPane.showMessageDialog(GenvisisWorkflowGUI.this,
																							"Error - cleanup of cancelled steps is not implemented.  Please clean or remove any generated files and try again.",
																							"Error", JOptionPane.ERROR_MESSAGE);
								boolean foundMore = false;
								for (int k = i + 1; k < steps.length; k++) {
									if (options[k]) {
										foundMore = true;
										break;
									}
								}
								boolean continueExec = false;
								if (foundMore) {
									int opt = JOptionPane.showConfirmDialog(GenvisisWorkflowGUI.this,
																													"A step was cancelled.  Do you wish to continue?",
																													"Step Cancelled",
																													JOptionPane.YES_NO_OPTION);
									if (opt == JOptionPane.YES_OPTION) {
										continueExec = true;
									}
								}
								if (!continueExec) {
									for (STEP step : steps) {
										step.resetRun();
									}
									break;
								}
							}

						}
					}
				}

				lockup(false);
				running = false;
			}
		});
		runThread.start();

	}

	private boolean checkRequirementsAndNotify(	HashMap<STEP, Boolean> selectedSteps,
																							HashMap<STEP, ArrayList<String>> variables) {
		boolean[] options = getSelectedOptions();

		ArrayList<String> reqMsgs = new ArrayList<String>();
		for (int i = 0; i < options.length; i++) {
			if (options[i]) {
				STEP step = steps[i];
				if (!step.hasRequirements(proj, selectedSteps, variables)) {
					reqMsgs.add((i + 1) + ". " + step.stepName);
				}
			}
		}
		boolean retVal = true;
		if (reqMsgs.size() > 0) {
			StringBuilder msg =
												new StringBuilder("Prerequisite requirements are unmet for the following steps:");
			for (String str : reqMsgs) {
				msg.append("\n").append(str);
			}
			proj.message(msg.toString());
			retVal = false;
		}
		return retVal;
	}

	private HashMap<STEP, ArrayList<String>> getVariables() {
		HashMap<STEP, ArrayList<String>> returnVars =
																								new HashMap<GenvisisWorkflow.STEP, ArrayList<String>>();
		for (Entry<STEP, ArrayList<? extends JComponent>> entry : varFields.entrySet()) {
			ArrayList<String> values = new ArrayList<String>();
			returnVars.put(entry.getKey(), values);
			for (JComponent j : entry.getValue()) {
				String val = "";
				if (j instanceof JTextField) {
					val = ((JTextField) j).getText().trim();
				} else if (j instanceof JCheckBox) {
					val = ((JCheckBox) j).isSelected() + "";
				} else if (j instanceof JSpinner) {
					val = ((JSpinner) j).getValue().toString();
				}
				values.add(val);
			}
		}
		return returnVars;
	}

	private Project createNewProject(String[] existing) {
		ProjectCreationGUI createGUI = new ProjectCreationGUI(existing);
		createGUI.setModal(true);
		createGUI.setVisible(true);

		if (createGUI.wasCancelled()) {
			return null;
		} else {
			return createGUI.getCreatedProject();
		}
	}


}
