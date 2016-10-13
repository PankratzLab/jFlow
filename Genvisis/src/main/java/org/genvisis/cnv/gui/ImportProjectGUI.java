package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.GROUP;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.prop.Property;
import org.genvisis.common.Files;
import org.genvisis.common.Grafik;
import org.genvisis.common.ext;

import net.miginfocom.swing.MigLayout;

/**
 * Dialog for creating a project properties file from an existing directory of data and files.
 */
public class ImportProjectGUI extends JDialog {

	private static final long serialVersionUID = 1L;

	private final JPanel contentPanel = new JPanel();
	private JTextField txtFldProjName;
	private JTextField txtFldProjDir;
	private JTextField txtFldSampleDataFile;
	private JLabel lblFoundProjectStatus;
	// private JLabel lblFoundSampleDataStatus;
	private JLabel lblFoundSamplesStatus;
	private JLabel lblFoundSampleListStatus;
	private JLabel lblFoundMarkerLookupStatus;
	private JLabel lblFoundDataStatus;
	private JLabel lblFoundMarkerListStatus;
	private JLabel lblFoundTransposedStatus;

	private static final ImageIcon redX = Grafik.getImageIcon("images/redx.png");
	private static final ImageIcon tick = Grafik.getImageIcon("images/tick.png");
	private static final String PROJECT = "PROJECT";
	private static final String SAMPLEDATA = "SAMPLEDATA";
	private static final String DEFAULT_PROJ_NAME = "New Project";

	volatile boolean cancelled = false;

	private final String propertyFilePath = MitoPipeline.initGenvisisProject();


	private final Action fileSelectAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;

		@Override
		public void actionPerformed(ActionEvent e) {
			JFileChooser jfc = new JFileChooser();
			String curr = ext.pwd();
			JTextField fld = null;
			if (PROJECT.equals(e.getActionCommand())) {
				String txt = txtFldProjDir.getText().trim();
				if (!txt.equals("") && !txt.equals("./")) {
					curr = txt;
				}
				jfc.setCurrentDirectory(new File(curr));
				jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				jfc.setDialogTitle("Select Project Directory");
				jfc.setMultiSelectionEnabled(false);
				fld = txtFldProjDir;
			} else if (SAMPLEDATA.equals(e.getActionCommand())) {
				String txt = txtFldSampleDataFile.getText().trim();
				if (!txt.equals("") && !txt.equals("./")) {
					curr = txt;
				} else {
					txt = txtFldProjDir.getText().trim();
					if (!txt.equals("") && !txt.equals("./")) {
						curr = txt;
					}
				}
				jfc.setCurrentDirectory(new File(curr));
				jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
				jfc.setDialogTitle("Select SampleData File");
				jfc.setMultiSelectionEnabled(false);
				fld = txtFldSampleDataFile;
			}
			int val = jfc.showDialog(ImportProjectGUI.this, "Select");
			if (val == JFileChooser.APPROVE_OPTION) {
				File f = jfc.getSelectedFile();
				if (fld != null) {
					String txt = ext.verifyDirFormat(f.getAbsolutePath());
					try {
						txt = ext.verifyDirFormat(f.getCanonicalPath());
					} catch (IOException e1) {
					}
					if (jfc.getFileSelectionMode() == FileChooser.FILES_ONLY	&& txt.length() > 1
							&& txt.endsWith("/")) {
						txt = txt.substring(0, txt.length() - 1);
					}
					// if a project directory is selected and the current project name is "" or "New project"
					// update the project name too
					if (PROJECT.equals(e.getActionCommand())
							&& (txtFldProjName.getText().isEmpty()
									|| txtFldProjName.getText().equals(DEFAULT_PROJ_NAME))) {
						txtFldProjName.setText(new File(txt).getName());
					}
					fld.setText(txt);
				}

			}

		}
	};

	/**
	 * Create the dialog.
	 */
	public ImportProjectGUI() {
		setTitle("Genvisis - Create Project from existing data");

		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		setBounds(100, 100, 350, 600);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(new MigLayout(	"", "[grow][10px:10px:10px][grow 70][grow 40]",
																					"[grow][][grow][][][][][][][][][][][][]"));
		CaretListener caretListener = new CaretListener() {
			@Override
			public void caretUpdate(CaretEvent arg0) {
				updateFound(getStatuses());
			}
		};
		Insets btnInsets = new Insets(-1, 3, -2, 3);
		{
			JLabel lblGenvisisProjectImport = new JLabel("Genvisis Project Import");
			lblGenvisisProjectImport.setFont(new Font("Arial", Font.BOLD, 16));
			contentPanel.add(lblGenvisisProjectImport, "cell 0 0 5 1,alignx center");
		}
		{
			JSeparator separator = new JSeparator();
			contentPanel.add(separator, "cell 0 1 5 1,growx");
		}
		{
			JLabel lblProjectName = new JLabel("New Project Name:");
			contentPanel.add(lblProjectName, "cel   l 0 3,alignx right");
		}
		{
			txtFldProjName = new JTextField(DEFAULT_PROJ_NAME);
			txtFldProjName.setColumns(10);
			contentPanel.add(txtFldProjName, "cell 2 3,growx");
		}
		{
			JLabel lblProjectDirectory = new JLabel("Project Directory:");
			contentPanel.add(lblProjectDirectory, "cell 0 4,alignx right");
		}
		{
			txtFldProjDir = new JTextField("./");
			txtFldProjDir.setColumns(10);
			txtFldProjDir.addCaretListener(caretListener);
			contentPanel.add(txtFldProjDir, "flowx,cell 2 4,growx");
		}
		{
			JButton btnSelectProjDir = new JButton(fileSelectAction);
			btnSelectProjDir.setText("...");
			btnSelectProjDir.setMargin(btnInsets);
			btnSelectProjDir.setActionCommand(PROJECT);
			contentPanel.add(btnSelectProjDir, "cell 2 4");
		}
		// {
		// JLabel lblSampleDataFile = new JLabel("<html><code>SampleData</code> File:</html>");
		// contentPanel.add(lblSampleDataFile, "cell 0 5,alignx right");
		// }
		// {
		// txtFldSampleDataFile = new JTextField("./");
		// txtFldSampleDataFile.setColumns(10);
		// txtFldSampleDataFile.addCaretListener(caretListener);
		// contentPanel.add(txtFldSampleDataFile, "flowx,cell 2 5,growx");
		// }
		// {
		// JButton btnSelectSampleDataFile = new JButton(fileSelectAction);
		// btnSelectSampleDataFile.setText("...");
		// btnSelectSampleDataFile.setMargin(btnInsets);
		// btnSelectSampleDataFile.setActionCommand("SAMPLEDATA");
		// contentPanel.add(btnSelectSampleDataFile, "cell 2 5");
		// }
		{
			JSeparator separator = new JSeparator();
			contentPanel.add(separator, "cell 0 6 4 1,growx");
		}
		{
			JPanel panel = new JPanel();
			contentPanel.add(panel, "cell 0 7 4 7,grow");
			panel.setLayout(new MigLayout("", "[grow][grow][][grow]", "[][][][][][][][][][][][]"));
			{
				JLabel lblRequired = new JLabel("<html><u>Required:</u><html>");
				panel.add(lblRequired, "cell 1 0,growy, alignx left,growy");
			}
			{
				JLabel lblValidProjectDirectory = new JLabel("Valid Project Directory:");
				panel.add(lblValidProjectDirectory, "cell 1 1,alignx right,growy");
				lblFoundProjectStatus = new JLabel(redX);
				panel.add(lblFoundProjectStatus, "cell 2 1,alignx left,growy");
			}
			JLabel lblFoundDataDirectory = new JLabel("<html>Found <code>data/</code> Directory:</html>");
			panel.add(lblFoundDataDirectory, "cell 1 2,alignx right,growy");
			// {
			// JLabel lblValidSampledataFile = new JLabel("<html>Valid <code>SampleData</code>
			// File:</html>");
			// panel.add(lblValidSampledataFile, "cell 4 4,alignx right");
			// }
			// lblFoundSampleDataStatus = new JLabel(redX);
			// panel.add(lblFoundSampleDataStatus, "cell 5 4,alignx left");
			{
				JLabel lblFoundSampleList = new JLabel("<html>Found <code>SampleList</code> File:</html>");
				panel.add(lblFoundSampleList, "cell 1 3,alignx right,growy");
			}
			{
				lblFoundSampleListStatus = new JLabel(redX);
				panel.add(lblFoundSampleListStatus, "cell 2 3,alignx left,growy");
			}
			{
				JLabel lblFoundMarkerList = new JLabel("<html>Found <code>MarkerSet</code> File:</html>");
				panel.add(lblFoundMarkerList, "cell 1 4,alignx right,growy");
			}
			{
				lblFoundDataStatus = new JLabel(redX);
				panel.add(lblFoundDataStatus, "cell 2 2,alignx left,growy");
			}
			lblFoundMarkerListStatus = new JLabel(redX);
			panel.add(lblFoundMarkerListStatus, "cell 2 4,alignx left,growy");
			{
				JLabel lblSample = new JLabel("<html><u>Sample Import:</u></html>");
				panel.add(lblSample, "cell 1 5,alignx left,growy");
			}
			JLabel lblFoundSamplesDirectory =
																			new JLabel("<html>Found <code>samples/</code> Directory:</html>");
			panel.add(lblFoundSamplesDirectory, "cell 1 6,alignx right,growy");
			lblFoundSamplesStatus = new JLabel(redX);
			panel.add(lblFoundSamplesStatus, "cell 2 6,alignx left,growy");
			{
				JLabel lblMarkerDataImport = new JLabel("<html><u>Marker Import:</u></html>");
				panel.add(lblMarkerDataImport, "cell 1 7,alignx left,growy");
			}
			{
				JLabel lblFoundTransposedDirectory = new JLabel("<html>Found <code>transposed/</code> Directory:</html>");
				panel.add(lblFoundTransposedDirectory, "cell 1 8,alignx right,growy");
			}
			lblFoundTransposedStatus = new JLabel(redX);
			panel.add(lblFoundTransposedStatus, "cell 2 8,alignx left,growy");
			JLabel lblFoundMarkerlookupFile =
																			new JLabel("<html>Found <code>MarkerLookup</code> File:</html>");
			panel.add(lblFoundMarkerlookupFile, "cell 1 9,alignx right,growy");
			{
				lblFoundMarkerLookupStatus = new JLabel(redX);
				panel.add(lblFoundMarkerLookupStatus, "cell 2 9,alignx left,growy");
			}
			{
				JLabel lblpropertyImport = new JLabel("<html><u>Property Import:</u></html>");
				panel.add(lblpropertyImport, "cell 1 10,growy");
			}
			{
				JLabel lblFoundImportMeta =
																	new JLabel("<html>Found <code>import.ser</code> meta File:</html>");
				panel.add(lblFoundImportMeta, "cell 1 11,alignx right,growy");
			}
			{
				lblFoundImportMetaStatus = new JLabel(redX);
				panel.add(lblFoundImportMetaStatus, "cell 2 11");
			}
		}
		{
			JLabel lblNewLabel = new JLabel("<html><center>Note: Remember to set other properties as needed!</center></html>");
			lblNewLabel.setFont(new Font("Tahoma", Font.BOLD, 14));
			contentPanel.add(lblNewLabel, "cell 0 14 4 1,alignx center");
		}
		{
			JPanel buttonPane = new JPanel();
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			buttonPane.setLayout(new MigLayout("", "[grow][47px][65px]", "[][23px]"));
			{
				JSeparator separator = new JSeparator();
				buttonPane.add(separator, "cell 0 0 3 1,growx");
			}
			{
				JButton okButton = new JButton("OK");
				okButton.addActionListener(new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent arg0) {
						close(false);
					}
				});
				okButton.setActionCommand("OK");
				buttonPane.add(okButton, "cell 1 1,alignx left,aligny top");
				getRootPane().setDefaultButton(okButton);
			}
			{
				JButton cancelButton = new JButton("Cancel");
				cancelButton.addActionListener(new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent arg0) {
						close(true);
					}
				});
				cancelButton.setActionCommand("Cancel");
				buttonPane.add(cancelButton, "cell 2 1,alignx left,aligny top");
			}
		}
		{
			addWindowListener(new WindowAdapter() {
				@Override
				public void windowClosing(WindowEvent e) {
					close(true);
					super.windowClosing(e);
				}
			});
		}
	}

	private void close(boolean cancelled) {
		if (!cancelled) {
			boolean[] statuses = getStatuses();
			if (statuses[0] && statuses[1] && statuses[2] && statuses[3]) {
				// good
				// if (statuses[4] && statuses[5]) {
				// // samples
				// }
				// if (statuses[6] && statuses[7]) {
				// // markers
				// }
			} else {
				return;
			}
		}
		this.cancelled = cancelled;
		setVisible(false);
	}

	private boolean checkProjectName() {
		String name = txtFldProjName.getText().trim();
		if (name.isEmpty()	|| DEFAULT_PROJ_NAME.equals(name) || name.length() > 23
				|| new File(propertyFilePath + name + MitoPipeline.PROJECT_EXT).exists()) {
			JOptionPane.showMessageDialog(null,
																		"Project name must be 1-23 characters in length, must not be \""
																						+ DEFAULT_PROJ_NAME
																					+ "\", and must not clash with an existing project.",
																		"Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;
	}

	// TODO these should ideally refer to Project/Property defaults
	private final String DEFAULT_SAMPLE_DIR = "samples/";
	private final String DEFAULT_DATA_DIR = "data/";
	private final String DEFAULT_TRANSPOSE_DIR = "transposed/";
	private final String DEFAULT_SAMPLELIST = DEFAULT_DATA_DIR + "samples.bis";
	private final String DEFAULT_SAMPLELIST_ALT = DEFAULT_DATA_DIR + "samples.ser";
	private final String DEFAULT_MARKERLIST = DEFAULT_DATA_DIR + "markers.bim";
	private final String DEFAULT_MARKERLIST_ALT = DEFAULT_DATA_DIR + "markers.ser";
	private final String DEFAULT_MARKERLOOKUP = DEFAULT_DATA_DIR + "markerLookup.bml";
	private final String DEFAULT_MARKERLOOKUP_ALT = DEFAULT_DATA_DIR + "markerLookup.ser";

	private JLabel lblFoundImportMetaStatus;

	private boolean[] getStatuses() {
		String baseDir = ext.verifyDirFormat(txtFldProjDir.getText().trim());
		boolean foundProject = Files.exists(baseDir);

		boolean foundSamples = Files.exists(baseDir + DEFAULT_SAMPLE_DIR);
		boolean foundData = Files.exists(baseDir + DEFAULT_DATA_DIR);
		boolean foundTransposed = Files.exists(baseDir + DEFAULT_TRANSPOSE_DIR);
		boolean foundSampleList = Files.exists(baseDir + DEFAULT_SAMPLELIST)
															|| Files.exists(baseDir + DEFAULT_SAMPLELIST_ALT);
		boolean foundMarkerList = Files.exists(baseDir + DEFAULT_MARKERLIST)
															|| Files.exists(baseDir + DEFAULT_MARKERLIST_ALT);
		boolean foundMarkerLookup = Files.exists(baseDir + DEFAULT_MARKERLOOKUP)
																|| Files.exists(baseDir + DEFAULT_MARKERLOOKUP_ALT);
		boolean foundImportMeta = Files.exists(baseDir + DEFAULT_DATA_DIR + Project.IMPORT_FILE);
		return new boolean[] {/* 0 */ foundProject, /* 1 */ foundData, /* 2 */ foundSampleList,
													/* 3 */ foundMarkerList,
													// /*4*/ foundSampleData,
													/* 5 */ foundSamples, /* 6 */ foundTransposed, /* 7 */ foundMarkerLookup,
													/* 8 */ foundImportMeta};
	}

	private void updateFound(final boolean[] statuses) {

		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				lblFoundProjectStatus.setIcon(statuses[0] ? tick : redX);
				lblFoundDataStatus.setIcon(statuses[1] ? tick : redX);
				lblFoundSampleListStatus.setIcon(statuses[2] ? tick : redX);
				lblFoundMarkerListStatus.setIcon(statuses[3] ? tick : redX);
				// lblFoundSampleDataStatus.setIcon(statuses[4] ? tick : redX);
				lblFoundSamplesStatus.setIcon(statuses[4] ? tick : redX);
				lblFoundTransposedStatus.setIcon(statuses[5] ? tick : redX);
				lblFoundMarkerLookupStatus.setIcon(statuses[6] ? tick : redX);
				lblFoundImportMetaStatus.setIcon(statuses[7] ? tick : redX);
			}
		});
	}

	public boolean getCancelled() {
		return cancelled;
	}

	public String getNewProjectFilename() {
		String name = txtFldProjName.getText().trim();
		String filename = propertyFilePath + name + MitoPipeline.PROJECT_EXT;
		return filename;
	}

	public boolean run() {
		String name = txtFldProjName.getText().trim();
		String filename = propertyFilePath + name + MitoPipeline.PROJECT_EXT;
		if (!checkProjectName() || Files.exists(filename)) {
			return false;
		} else {
			Files.write((new Project()).PROJECT_NAME.getName() + "=" + name, filename);
		}
		
		String projDir = txtFldProjDir.getText().trim();
		Project actualProj = new Project(filename, false);
		actualProj.PROJECT_NAME.setValue(name);
		actualProj.PROJECT_DIRECTORY.setValue(projDir);

		Map<String, String> importProps = actualProj.loadImportMetaFile();
		if (importProps != null && !importProps.isEmpty()) {
			List<Property<?>> props = actualProj.getProperties(GROUP.IMPORT);
			for (Property<?> p : props) {
				if (importProps.containsKey(p.getName())) {
					p.parseValue(importProps.get(p.getName()));
				}
			}
		}
		
        boolean foundOldSampleList = Files.exists(projDir + DEFAULT_SAMPLELIST) && !Files.exists(projDir + DEFAULT_SAMPLELIST_ALT);
        boolean foundOldMarkerList = Files.exists(projDir + DEFAULT_MARKERLIST) && !Files.exists(projDir + DEFAULT_MARKERLIST_ALT);
        boolean foundOldMarkerLookup = Files.exists(projDir + DEFAULT_MARKERLOOKUP) && !Files.exists(projDir + DEFAULT_MARKERLOOKUP_ALT);
        
        if (foundOldMarkerList) {
          actualProj.MARKERSET_FILENAME.setValue(projDir + DEFAULT_MARKERLIST);
        }
        if (foundOldMarkerLookup) {
          actualProj.MARKERLOOKUP_FILENAME.setValue(projDir + DEFAULT_MARKERLOOKUP);
        }
        if (foundOldSampleList) {
          actualProj.SAMPLELIST_FILENAME.setValue(projDir + DEFAULT_SAMPLELIST);
        }

		// actualProj.SOURCE_DIRECTORY.setValue(srcDir);
		// actualProj.SOURCE_FILENAME_EXTENSION.setValue(srcExt);
		// actualProj.LRRSD_CUTOFF.setValue(lrrSd);
		// actualProj.XY_SCALE_FACTOR.setValue(xy);
		// TODO should set the following two?
		// actualProj.TARGET_MARKERS_FILENAMES.setValue(new String[]{ext.removeDirectoryInfo(tgtMkrs)});
		// actualProj.ARRAY_TYPE.setValue((ARRAY) comboBoxArrayType.getSelectedItem());
		// actualProj.ID_HEADER.setValue(sampCol == SourceFileHeaderGUI.FILENAME_IND ?
		// SourceFileParser.FILENAME_AS_ID_OPTION : cols[sampCol]);
		// actualProj.SOURCE_FILE_DELIMITER.setValue(SOURCE_FILE_DELIMITERS.getDelimiter(sourceDelim));

		// actualProj.SAMPLE_DATA_FILENAME.setValue(txtFldSampleDataFile.getText().trim());
		actualProj.saveProperties();
		return true;
	}



}
