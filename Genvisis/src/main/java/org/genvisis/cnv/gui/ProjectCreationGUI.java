package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import org.genvisis.cnv.LaunchProperties;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import org.genvisis.cnv.filesys.SourceFileHeaderData;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.manage.SourceFileParser;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Grafik;
import org.genvisis.common.ext;

import net.miginfocom.swing.MigLayout;

public class ProjectCreationGUI extends JDialog {

	private static final long serialVersionUID = 1L;

	private static final String PROJECT_DIR_TOOLTIP = "<html>Directory in which to create, store, and manage all project files.</html>";
	private static final String SOURCE_DIR_TOOLTIP = "<html>Directory of source (e.g. FinalReport.txt.gz) files; this can be different than the Project Directory.</html>";
	private static final String SOURCE_EXT_TOOLTIP = "<html>Extension of source files (e.g. for \"FinalReport.txt.gz\", the extension would be \".txt.gz\".</html>";
	private static final String XY_TOOLTIP = "<html>Suggested values for X/Y correction:<br />Illumina: use the default of 1.<br />Affymetrix: use 100.<br />DBGAP: use 2000.</html>";
	private final JPanel contentPane;
	private final JTextField txtFldProjName;
	private final JTextField txtFldProjDir;
	private final JTextField txtFldSrcDir;
	private final JTextField txtFldSrcExt;
	private final JTextField txtFldTgtMkrs;
	// private JSpinner spinnerLrrSd;
	private Project proj;
	private volatile boolean cancelled = false;
	private final Action fileSelectAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;

		@Override
		public void actionPerformed(ActionEvent e) {
			JFileChooser jfc = new JFileChooser();
			String curr = ext.pwd();
			JTextField fld = null;
			if (e.getActionCommand().equals("PROJECT")) {
				String txt = txtFldProjDir.getText().trim();
				if (!txt.equals("") && !txt.equals("./")) {
					curr = txt;
				}
				jfc.setCurrentDirectory(new File(curr));
				jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				jfc.setDialogTitle("Select Project Directory");
				jfc.setMultiSelectionEnabled(false);
				fld = txtFldProjDir;
			} else if (e.getActionCommand().equals("SOURCE")) {
				String txt = txtFldSrcDir.getText().trim();
				if (!txt.equals("") && !txt.equals("./")) {
					curr = txt;
				} else {
					txt = txtFldProjDir.getText().trim();
					if (!txt.equals("") && !txt.equals("./")) {
						curr = txt;
					}
				}
				jfc.setCurrentDirectory(new File(curr));
				jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				jfc.setDialogTitle("Select Source Directory");
				jfc.setMultiSelectionEnabled(false);
				fld = txtFldSrcDir;
			} else if (e.getActionCommand().equals("TARGET")) {
				String txt = txtFldProjDir.getText().trim();
				if (!txt.equals("") && !txt.equals("./")) {
					curr = txt;
				}
				jfc.setCurrentDirectory(new File(curr));
				jfc.setDialogTitle("Select Target Markers File");
				jfc.setMultiSelectionEnabled(false);
				fld = txtFldTgtMkrs;
			}

			int val = jfc.showDialog(ProjectCreationGUI.this, "Select");
			if (val == JFileChooser.APPROVE_OPTION) {
				File f = jfc.getSelectedFile();
				if (fld != null) {
					String txt = ext.verifyDirFormat(f.getAbsolutePath());
					try {
						txt = ext.verifyDirFormat(f.getCanonicalPath());
					} catch (IOException e1) {
					}
					if (jfc.getFileSelectionMode() == FileChooser.FILES_ONLY && txt.length() > 1
							&& txt.endsWith("/")) {
						txt = txt.substring(0, txt.length() - 1);
					}
					fld.setText(txt);
				}

			}

		}
	};

	private final JComboBox comboBoxArrayType;
	// private JComboBox<String> comboBoxArrayType;
	private final JLabel lblSrcFileStatus;
	private final JSpinner spinnerXY;
	private final JProgressBar progressBar;
	private final String[] existingNames;

	private boolean multipleExts = false;
	private String currentDir = "";

	public static Project runCreationGUI() {
		ProjectCreationGUI createGUI = new ProjectCreationGUI(LaunchProperties.getListOfProjectNames());
		createGUI.setModal(true);
		createGUI.setVisible(true);
		UITools.setSize(createGUI, new Dimension(550, 500));
		createGUI.pack();

		if (createGUI.wasCancelled()) {
			return null;
		} else {
			return createGUI.getCreatedProject();
		}
	}

	/**
	 * Create the frame.
	 */
	public ProjectCreationGUI(String[] existingProjectNames) {
		setLayout(new BorderLayout());
		setTitle("Genvisis: Create New Project");
		existingNames = existingProjectNames;
		Project proj = new Project();

		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(new MigLayout("", "[grow][10px:10px:10px][grow][grow]",
																				"[grow][][grow][][][grow][][][][][grow][][][][grow][]"));

		JLabel lblGenvisisProjectCreation = new JLabel("Genvisis Project Creation");
		lblGenvisisProjectCreation.setFont(new Font("Arial", Font.BOLD, 16));
		contentPane.add(lblGenvisisProjectCreation, "cell 0 0 4 1,alignx center");

		JSeparator separator = new JSeparator();
		contentPane.add(separator, "cell 0 1 4 1,growx");

		JLabel lblProjectName = new JLabel("Project Name:");
		contentPane.add(lblProjectName, "cell 0 3,alignx trailing");

		txtFldProjName = new JTextField(proj.PROJECT_NAME.getDefaultValueString());
		contentPane.add(txtFldProjName, "cell 2 3,growx");
		txtFldProjName.setColumns(10);

		JLabel lblProjectDirectory = new JLabel("Project Directory:");
		contentPane.add(lblProjectDirectory, "cell 0 4,split 2, alignx right");
		contentPane.add(Grafik.getToolTipIconLabel(PROJECT_DIR_TOOLTIP), "cell 0 4, alignx right");

		txtFldProjDir = new JTextField(proj.PROJECT_DIRECTORY.getDefaultValueString());
		contentPane.add(txtFldProjDir, "flowx,cell 2 4,growx");
		txtFldProjDir.setColumns(10);

		JLabel lblSourceFileDirectory = new JLabel("Source File Directory:");
		contentPane.add(lblSourceFileDirectory, "cell 0 6,split 2, alignx right");
		contentPane.add(Grafik.getToolTipIconLabel(SOURCE_DIR_TOOLTIP), "cell 0 6, alignx right");

		CaretListener checkSource = new CaretListener() {
			@Override
			public void caretUpdate(CaretEvent e) {
				updateSourceFileNotice();
			}
		};



		txtFldSrcDir = new JTextField(proj.SOURCE_DIRECTORY.getDefaultValueString());
		txtFldSrcDir.addCaretListener(checkSource);
		contentPane.add(txtFldSrcDir, "flowx,cell 2 6,growx");
		txtFldSrcDir.setColumns(10);

		JLabel lblSourceFileExtension = new JLabel("Source File Extension:");
		contentPane.add(lblSourceFileExtension, "cell 0 7,split 2,alignx right");
		contentPane.add(Grafik.getToolTipIconLabel(SOURCE_EXT_TOOLTIP), "cell 0 7, alignx right");

		txtFldSrcExt = new JTextField(proj.SOURCE_FILENAME_EXTENSION.getDefaultValueString());
		txtFldSrcExt.addCaretListener(checkSource);
		contentPane.add(txtFldSrcExt, "cell 2 7,growx");
		txtFldSrcExt.setColumns(10);

		JLabel lblArrayType = new JLabel("Array Type:");
		contentPane.add(lblArrayType, "cell 0 8,alignx trailing");

		comboBoxArrayType = new JComboBox(Project.ARRAY.values());
		// comboBoxArrayType = new JComboBox<String>();
		comboBoxArrayType.setFont(comboBoxArrayType.getFont().deriveFont(Font.PLAIN));
		contentPane.add(comboBoxArrayType, "cell 2 8,growx");

		lblSrcFileStatus = new JLabel("");
		contentPane.add(lblSrcFileStatus, "cell 2 10,alignx right,aligny top");

		JLabel lblXyCorrectionRatio = new JLabel("X/Y Correction Ratio:");
		contentPane.add(lblXyCorrectionRatio, "cell 0 11,split 2, alignx right");
		contentPane.add(Grafik.getToolTipIconLabel(XY_TOOLTIP), "cell 0 11, alignx right");

		spinnerXY = new JSpinner();
		spinnerXY.setModel(new SpinnerNumberModel(proj.XY_SCALE_FACTOR.getDefaultValue().doubleValue(),
																							0.001, 10000000000d, 1.0));
		contentPane.add(spinnerXY, "cell 2 11,growx");

		// JLabel lblLogrRatioStddev = new JLabel("Log-R Ratio Std.Dev. Cut-off Threshold:");
		// contentPane.add(lblLogrRatioStddev, "cell 0 12,alignx trailing");

		// spinnerLrrSd = new JSpinner();
		// spinnerLrrSd.setModel(new
		// SpinnerNumberModel(proj.LRRSD_CUTOFF.getDefaultValue().doubleValue(), 0.0, 3.0, 0.1));
		// contentPane.add(spinnerLrrSd, "cell 2 12,growx");

		JLabel lblTargetMarkersFile = new JLabel("[Optional] Target Markers File:");
		contentPane.add(lblTargetMarkersFile, "cell 0 12,alignx trailing");

		txtFldTgtMkrs = new JTextField(proj.TARGET_MARKERS_FILENAMES.getDefaultValueString());
		contentPane.add(txtFldTgtMkrs, "flowx,cell 2 12,growx");
		txtFldTgtMkrs.setColumns(10);

		Insets fileBtnInsets = new Insets(0, 3, 0, 3);
		JButton btnProjDir = new JButton(fileSelectAction);
		btnProjDir.setMargin(fileBtnInsets);
		btnProjDir.setText("...");
		btnProjDir.setActionCommand("PROJECT");
		contentPane.add(btnProjDir, "cell 2 4");

		JButton btnSrcDir = new JButton(fileSelectAction);
		btnSrcDir.setMargin(fileBtnInsets);
		btnSrcDir.setText("...");
		btnSrcDir.setActionCommand("SOURCE");
		contentPane.add(btnSrcDir, "cell 2 6");

		JButton btnTgtMkrs = new JButton(fileSelectAction);
		btnTgtMkrs.setMargin(fileBtnInsets);
		btnTgtMkrs.setText("...");
		btnTgtMkrs.setActionCommand("TARGET");
		contentPane.add(btnTgtMkrs, "cell 2 12");

		JSeparator separator_1 = new JSeparator();
		contentPane.add(separator_1, "cell 0 15 4 1,growx");

		JPanel panel = new JPanel();
		contentPane.add(panel, "south");
		panel.setLayout(new MigLayout("", "[grow][]", "[]"));

		JButton btnCreate = new JButton("Create [Skip Validation]");
		btnCreate.setToolTipText("<html>Create new project, skipping source file validation.<br />  ONLY select if all source files are guaranteed <br />to be correct, valid, and uniform in structure.</html>");
		btnCreate.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if (checkValues()) {
					int resp = JOptionPane.showConfirmDialog(ProjectCreationGUI.this,
																									 "<html>You are waiving the opporunity to review your project structure.<br />Are you sure that all source files are valid, correct, and uniform in structure?<br /><br />[If not, select 'Validate and Create' to interactively review project structure]</html>",
																									 "Confirm File Validity",
																									 JOptionPane.YES_NO_OPTION);
					if (resp == JOptionPane.YES_OPTION) {
						// Don't care if createProject doesn't work nicely
						if (!createProject(false)) {
							JOptionPane.showMessageDialog(ProjectCreationGUI.this,
																						"Could not create project - please check your inputs and try again.",
																						"Project Creation Failed", JOptionPane.ERROR_MESSAGE);
						} else {
							doClose(false);
						}
					}
				}
			}
		});

		progressBar = new JProgressBar();
		progressBar.setVisible(false);
		panel.add(progressBar, "hidemode 2, cell 0 0");

		btnCreate.setMnemonic(KeyEvent.VK_F);
		panel.add(btnCreate, "flowx,cell 1 0");

		JButton btnCreateAndValidate = new JButton("Validate and Create");
		btnCreateAndValidate.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if (checkValues()) {
					if (createProject(true)) {
						doClose(false);
					}
				}
			}
		});
		btnCreateAndValidate.setMnemonic(KeyEvent.VK_V);
		panel.add(btnCreateAndValidate, "cell 1 0");

		JButton btnCancel = new JButton("Cancel");
		btnCancel.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				doClose(true);
			}
		});
		btnCancel.setMnemonic(KeyEvent.VK_C);
		panel.add(btnCancel, "cell 1 0");

		updateSourceFileNotice();
	}

	private boolean checkValues() {
		String name = txtFldProjName.getText().trim();
		String projDir = txtFldProjDir.getText().trim();
		String srcDir = txtFldSrcDir.getText().trim();
		String srcExt = txtFldSrcExt.getText().trim();
		// String idHdr = txtFldIDHdr.getText().trim();
		// double lrrSd = ((Double)spinnerLrrSd.getValue()).doubleValue();
		// double xy = ((Double)spinnerXY.getValue()).doubleValue();
		String tgtMkrs = txtFldTgtMkrs.getText().trim();

		boolean validProjDir = false;
		try {
			File f = (new File(projDir));
			f.getCanonicalPath();
			validProjDir = true; // f.exists();
		} catch (IOException e) {
		}
		boolean validSrcDir = false;
		try {
			File f = (new File(srcDir));
			f.getCanonicalPath();
			validSrcDir = f.exists();
		} catch (IOException e) {
		}
		boolean validTgtMkrs = false;
		try {
			File f = (new File(tgtMkrs));
			f.getCanonicalPath();
			validTgtMkrs = true;
		} catch (IOException e) {
		}

		boolean nameCheck = !name.equals("")
												&& !name.equals((new Project()).PROJECT_NAME.getDefaultValueString());
		if (nameCheck) {
			for (String s : existingNames) {
				if (s.equals(name)) {
					nameCheck = false;
					break;
				}
			}
		}

		boolean[] checks = {nameCheck, validProjDir, validSrcDir, !srcExt.equals(""),
												// !idHdr.equals("") ,
												validTgtMkrs};// Files.exists(tgtMkrs)};

		if (ArrayUtils.booleanArraySum(checks) < checks.length) {
			StringBuilder errorMsg = new StringBuilder();
			for (int i = 0; i < checks.length; i++) {
				if (!checks[i]) {
					errorMsg.append("\n").append(getErrorFor(i));
				}
			}

			JOptionPane.showMessageDialog(null, errorMsg.toString(), "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;

	}

	private void updateSourceFileNotice() {
		int cnt = getValidSourceFileCount();
		if (cnt <= 0) {
			lblSrcFileStatus.setText("No source files found");
			lblSrcFileStatus.setForeground(Color.RED);
		} else if (!multipleExts) {
			lblSrcFileStatus.setText(cnt + " source file(s) found");
			lblSrcFileStatus.setForeground(Color.GREEN.darker());
		} else {
			lblSrcFileStatus.setText("WARNING: Multiple file types detected. \n" + cnt
															 + " source file(s) found");
			lblSrcFileStatus.setForeground(Color.ORANGE.darker());
		}
		lblSrcFileStatus.setVisible(true);
	}

	private int getValidSourceFileCount() {
		String srcDir = txtFldSrcDir.getText().trim();
		String srcExt = txtFldSrcExt.getText().trim();

		boolean validSrcDir = false;
		try {
			File f = (new File(srcDir));
			f.getCanonicalPath();
			validSrcDir = f.exists();
		} catch (IOException e) {
		}

		if (!validSrcDir)
			return -1;

		if (!currentDir.equals(srcDir)) {
			HashMap<String, Integer> extensions = new HashMap<String, Integer>();
			// look for the most likely extension
			for (String s : (new File(srcDir).list())) {
				String[] split = s.split("\\.", 2); // only split on the first . to capture things like
																						// .tar.gz

				if (split.length <= 1)
					continue;

				String ext = "." + split[1];

				if (!extensions.containsKey(ext))
					extensions.put(ext, 1);
				else
					extensions.replace(ext, extensions.get(ext) + 1);
			}

			if (extensions.isEmpty())
				return -1;

			multipleExts = extensions.size() > 1;
			int mostCommonExt = -1;

			for (String key : extensions.keySet()) {
				int freq = extensions.get(key);
				if (freq > mostCommonExt) {
					srcExt = key;
					mostCommonExt = freq;
				}
			}

			currentDir = srcDir;
		}

		if ("".equals(srcExt)) {
			return -1;
		}

		final String srcExtFinal = srcExt;

		if (!srcExt.equals(txtFldSrcExt.getText().trim())) {

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					txtFldSrcExt.setText(srcExtFinal);
				}
			});
		}
		String[] files = (new File(srcDir)).list(new FilenameFilter() {
			@Override
			public boolean accept(File arg0, String arg1) {
				return arg1.endsWith(srcExtFinal);
			}
		});

		return files.length;
		// TODO check if files are valid, check if file headers are valid, check if file headers contain
		// ID field
		// TODO use FinalReport header class (to be constructed)
	}

	private String getErrorFor(int index) {
		switch (index) {
			case 0: // Project Name
				return "Project name must have a value, must not be 'New Project', and may not already exist";
			case 1: // Project dir
				return "Project directory must be a valid directory path";
			case 2: // src dir
				return "Source directory must be a valid and existing directory path";
			case 3: // src ext
				return "Source file extension cannot be empty";
			case 4: // id hdr
				return "ID header cannot be empty";
			case 5: // tgt mkrs
				return "Target markers must be a valid file";
			default:
				return "Unknown error - please check all options and try again";
		}
	}

	private boolean createProject(boolean actuallyValidate) {
		String name = txtFldProjName.getText().trim();
		String projDir = txtFldProjDir.getText().trim();
		String srcDir = txtFldSrcDir.getText().trim();
		String srcExt = txtFldSrcExt.getText().trim();
		// String idHdr = txtFldIDHdr.getText().trim();
		// double lrrSd = ((Double)spinnerLrrSd.getValue()).doubleValue();
		double xy = ((Double) spinnerXY.getValue()).doubleValue();
		String tgtMkrs = txtFldTgtMkrs.getText().trim();

		HashMap<String, SourceFileHeaderData> headers = SourceFileHeaderData.validate(srcDir, srcExt,
																																									actuallyValidate,
																																									new org.genvisis.common.Logger(),
																																									progressBar);
		if (headers == null) {
			// errors found in headers - check output and retry?
			return false;
		}
		SourceFileHeaderData reportHdr = null;
		for (SourceFileHeaderData d : headers.values()) {
			if (reportHdr == null) {
				reportHdr = d;
				break;
			}
		}
		// TODO do column assignment
		SourceFileHeaderGUI gui = new SourceFileHeaderGUI(reportHdr);
		if (actuallyValidate) {
			gui.setModal(true);
			gui.setVisible(true);
			if (gui.wasCancelled()) {
				return false;
			}
		}
		String sourceDelim = null;
		int sampCol = -100;
		String[] cols = null;
		for (SourceFileHeaderData d : headers.values()) {
			if (sourceDelim == null) {
				sourceDelim = d.getSourceFileDelimiter();
			}
			if (cols == null) {
				cols = d.cols;
			}
			d.colSnpIdent = gui.getSelectedSNPIndex();
			d.colSampleIdent = gui.getSelectedSampleID();
			if (sampCol == -100) {
				sampCol = d.colSampleIdent;
			}
			d.colGeno1 = gui.getSelectedGeno1();
			d.colGeno2 = gui.getSelectedGeno2();
			d.colGenoAB1 = gui.getSelectedAB1();
			d.colGenoAB2 = gui.getSelectedAB2();
			d.colBAF = gui.getSelectedBAF();
			d.colLRR = gui.getSelectedLRR();
			d.colGC = gui.getSelectedGC();
			d.colR = gui.getSelectedR();
			d.colTheta = gui.getSelectedTheta();
			d.colX = gui.getSelectedX();
			d.colY = gui.getSelectedY();
			d.colXRaw = gui.getSelectedXRaw();
			d.colYRaw = gui.getSelectedYRaw();
		}

		String path = MitoPipeline.initGenvisisProject();
		File file = new File(projDir);

		Project dummy = new Project();
		dummy.setGuiState(true);
		if (!file.exists()) {
			if (file.mkdirs()) {
				dummy.getLog().report("\n" + ext.getTime() + " Created directory " + projDir);
			} else {
				dummy.getLog().reportError("Error - failed to create  " + projDir
																	 + ", please manually create it unless it already exists");
				dummy.message("Error - failed to create  " + projDir
											+ ", please manually create it unless it already exists");
				return false;
			}
		}
		String filename = path + name + MitoPipeline.PROJECT_EXT;
		if (Files.exists(filename)) {
			dummy.getLog().reportError("Project " + name + " already exists");
			dummy.message("Error - Project " + name + " already exists");
			return false;
		} else {
			Files.write((new Project()).PROJECT_NAME.getName() + "=" + name, filename);
		}
		Project actualProj = new Project(filename, false);
		actualProj.PROJECT_NAME.setValue(name);
		actualProj.PROJECT_DIRECTORY.setValue(projDir);
		actualProj.SOURCE_DIRECTORY.setValue(srcDir);
		actualProj.SOURCE_FILENAME_EXTENSION.setValue(srcExt);
		// actualProj.LRRSD_CUTOFF.setValue(lrrSd);
		actualProj.XY_SCALE_FACTOR.setValue(xy);
		actualProj.TARGET_MARKERS_FILENAMES.setValue(new String[] {ext.removeDirectoryInfo(tgtMkrs)});
		actualProj.ARRAY_TYPE.setValue((ARRAY) comboBoxArrayType.getSelectedItem());
		// if (abLookup != null && Files.exists(projectDirectory + abLookup)) {
		// proj.setProperty(proj.AB_LOOKUP_FILENAME, ext.removeDirectoryInfo(abLookup));
		// }
		actualProj.ID_HEADER.setValue(sampCol == SourceFileHeaderGUI.FILENAME_IND ? SourceFileParser.FILENAME_AS_ID_OPTION
																																							: cols[sampCol]);
		actualProj.SOURCE_FILE_DELIMITER.setValue(SOURCE_FILE_DELIMITERS.getDelimiter(sourceDelim));
		actualProj.saveProperties();
		actualProj.setSourceFileHeaders(headers);
		proj = actualProj;
		return true;
	}

	public Project getCreatedProject() {
		return proj;
	}

	private void doClose(boolean cancel) {
		cancelled = cancel;
		setVisible(false);
	}

	public boolean wasCancelled() {
		return cancelled;
	}
}
