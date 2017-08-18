package org.genvisis.cnv.plots;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.ListCellRenderer;
import javax.swing.ListSelectionModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;

import net.miginfocom.swing.MigLayout;

import org.genvisis.cnv.plots.data.DataFile;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;

public class ManhattanLoadGUI extends JDialog {

	private final String TEMP_DIR = "D:/data/ny_choanal/omni2.5v1.2/plink/";

	private final JPanel contentPanel = new JPanel();
	private JTextField textField;

	private Logger log = new Logger(); // TODO
	private JList<JCheckBox> list;
	private JList<JCheckBox> chrList;

	protected static Border noFocusBorder = new EmptyBorder(1, 1, 1, 1);

	private JPanel filterPanel;
	private int closeCode = JOptionPane.CANCEL_OPTION;

	private DataFile dataFile;

	HashMap<Integer, JCheckBox> chrBoxes = new HashMap<>();
	HashMap<String, JCheckBox> columnBoxes = new HashMap<>();

	private JSpinner spinner;

	private JCheckBox chkPValFilter;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		try {
			ManhattanLoadGUI dialog = new ManhattanLoadGUI();
			dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
			dialog.setVisible(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Create the dialog.
	 */
	public ManhattanLoadGUI() {
		setBounds(100, 100, 450, 456);
		setTitle("Load File");
		setModal(true);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(new MigLayout("ins 0", "[424px,grow][][][424px,grow]", "[][][][grow]"));
		{
			JLabel lblSelectFile = new JLabel("Select File:");
			contentPanel.add(lblSelectFile, "cell 0 0 4 1");
		}
		{
			textField = new JTextField();
			contentPanel.add(textField, "flowx,cell 0 1 4 1,growx");
			textField.setColumns(10);
		}
		{
			JLabel lblSelectColumnsTo = new JLabel("Select Columns to Load:");
			contentPanel.add(lblSelectColumnsTo, "cell 0 2");
		}
		{
			JScrollPane scrollPane = new JScrollPane();
			contentPanel.add(scrollPane, "cell 0 3,grow");
			{
				list = new JList<>();
				scrollPane.setViewportView(list);
				setupList(list);
			}
		}
		{
			JSeparator separator = new JSeparator();
			separator.setOrientation(SwingConstants.VERTICAL);
			contentPanel.add(separator, "cell 1 2 1 2,growy");
		}
		filterPanel = new JPanel();
		contentPanel.add(filterPanel, "cell 3 2 1 2,grow");
		filterPanel.setLayout(new MigLayout("ins 5 0 0 5", "[grow]", "[][][grow][]"));
		{
			chkPValFilter = new JCheckBox("P-Value Filter:");
			chkPValFilter.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					spinner.setEnabled(chkPValFilter.isSelected());
				}
			});
			chkPValFilter.setSelected(true);
			chkPValFilter.setFont(new Font("Tahoma", Font.PLAIN, 11));
			filterPanel.add(chkPValFilter, "flowx,cell 0 0,alignx left");
		}
		{
			JLabel lblSelectChrs = new JLabel("Select Chrs:");
			lblSelectChrs.setFont(new Font("Tahoma", Font.PLAIN, 11));
			filterPanel.add(lblSelectChrs, "cell 0 1,alignx center");
		}
		{
			JScrollPane scrollPane = new JScrollPane();
			filterPanel.add(scrollPane, "cell 0 2,grow");
			chrList = new JList<JCheckBox>();
			scrollPane.setViewportView(chrList);
			setupList(chrList);
			DefaultListModel<JCheckBox> dlm = new DefaultListModel<>();
			for (int i = 0; i < 27; i++) {
				String txt;
				switch (i) {
					case 23:
						txt = "X";
						break;
					case 24:
						txt = "Y";
						break;
					case 25:
						txt = "XY";
						break;
					case 26:
						txt = "MT";
						break;
					default:
						txt = i + "";
						break;
				}
				JCheckBox box = new JCheckBox(txt);
				chrBoxes.put(i, box);
				box.setSelected(i > 0);
				dlm.addElement(box);
			}
			chrList.setModel(dlm);
		}
		{
			JLabel lblLessthan = new JLabel("<=");
			lblLessthan.setHorizontalAlignment(SwingConstants.TRAILING);
			filterPanel.add(lblLessthan, "cell 0 0,alignx center");
		}
		{
			spinner = new JSpinner();
			spinner.setModel(new SpinnerNumberModel(0.05, 0.0, 1.0, 0.01));
			filterPanel.add(spinner, "cell 0 0,growx");
		}
		{
			JButton btnSelectAll = new JButton("Select All");
			btnSelectAll.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					for (JCheckBox box : chrBoxes.values()) {
						box.setSelected(true);
					}
					chrList.repaint();
				}
			});
			filterPanel.add(btnSelectAll, "flowx,cell 0 3,growx");
		}
		{
			JButton btnSelectNone = new JButton("Select None");
			btnSelectNone.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					for (JCheckBox box : chrBoxes.values()) {
						box.setSelected(false);
					}
					chrList.repaint();
				}
			});
			filterPanel.add(btnSelectNone, "cell 0 3,growx");
		}
		{
			JButton btnLoad = new JButton(">");
			btnLoad.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					selectFile();
				}
			});
			contentPanel.add(btnLoad, "cell 0 1 4 1");
		}
		{
			JPanel buttonPane = new JPanel();
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				JButton okButton = new JButton("OK");
				okButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						if (textField.getText().trim().equals("")) {
							return;
						}
						if (!Files.exists(textField.getText().trim())) {
							JOptionPane.showMessageDialog(ManhattanLoadGUI.this, "Error - file not found!",
																						"Error - Missing File!",
																						JOptionPane.ERROR_MESSAGE);
							return;
						}
						close(false);
					}
				});
				okButton.setActionCommand("OK");
				buttonPane.add(okButton);
				getRootPane().setDefaultButton(okButton);
			}
			{
				JButton cancelButton = new JButton("Cancel");
				cancelButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						close(true);
					}
				});
				cancelButton.setActionCommand("Cancel");
				buttonPane.add(cancelButton);
			}
		}
	}

	private void setupList(JList<JCheckBox> list) {
		list.setCellRenderer(new CellRenderer());
		list.addMouseListener(new MouseAdapter() {
			int prevIndex = -1;

			public void mousePressed(MouseEvent e) {
				int index = list.locationToIndex(e.getPoint());
				if (index != -1) {
					JCheckBox checkbox = (JCheckBox) list.getModel().getElementAt(index);
					if (checkbox.isEnabled()) {
						boolean cl = e.getClickCount() >= 2;
						boolean li = prevIndex != -1 && index == prevIndex;
						if (cl || li) {
							checkbox.setSelected(!checkbox.isSelected());
						}
						prevIndex = index;
					}
					repaint();
				}
			}
		});
		list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
	}

	private void close(boolean cancel) {
		if (!cancel) {
			this.closeCode = JFileChooser.APPROVE_OPTION;
		}
		setVisible(false);
	}

	private void selectFile() {
		JFileChooser jfc = new JFileChooser(TEMP_DIR);
		jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
		jfc.setMultiSelectionEnabled(false);
		jfc.setDialogTitle("Select File...");
		int code = jfc.showOpenDialog(this);
		if (code == JFileChooser.APPROVE_OPTION) {
			File f = jfc.getSelectedFile();
			String fil = f.getAbsolutePath();
			try {
				fil = f.getCanonicalPath();
			} catch (IOException e) {
			}
			textField.setText(fil);
			update();
		}
	}

	private void update() {
		String fil = textField.getText().trim();
		if ("".equals(fil) || !Files.exists(fil)) {
			// TODO
			return;
		}

		DataFile dataFile = new DataFile(fil, log, ManhattanPlot.LINKERS, ManhattanPlot.REQ_LINKS);
		if (!dataFile.hasRequiredData()) {
			boolean missChr = dataFile.getLinkedColumnName(ManhattanPlot.CHR_LINKER) == null;
			boolean missPos = dataFile.getLinkedColumnName(ManhattanPlot.POS_LINKER) == null;
			boolean missPVal = dataFile.getLinkedColumnName(ManhattanPlot.PVAL_LINKER) == null;
			String miss = "";
			if (missChr) {
				miss += "Chr";
			}
			if (missPos) {
				if (miss.length() > 0) {
					miss += ", ";
				}
				miss += "Pos";
			}
			if (missPVal) {
				if (miss.length() > 0) {
					miss += ", ";
				}
				miss += "P-val";
			}
			JOptionPane.showMessageDialog(this, "Error - missing " + miss + " data!",
																		"Error - Missing Data!",
																		JOptionPane.ERROR_MESSAGE);
			return;
		}
		this.dataFile = dataFile;

		String[] columns = dataFile.getColumns();
		int[] links = dataFile.getLinkIndices();

		DefaultListModel<JCheckBox> listModel = new DefaultListModel<>();
		String chrStr = dataFile.getLinkedColumnName(ManhattanPlot.CHR_LINKER);
		JCheckBox chrBox = new JCheckBox(chrStr);
		chrBox.setSelected(true);
		chrBox.setEnabled(false);
		columnBoxes.put(chrStr, chrBox);
		String posStr = dataFile.getLinkedColumnName(ManhattanPlot.POS_LINKER);
		JCheckBox posBox = new JCheckBox(posStr);
		posBox.setSelected(true);
		posBox.setEnabled(false);
		columnBoxes.put(posStr, posBox);
		String pvalStr = dataFile.getLinkedColumnName(ManhattanPlot.PVAL_LINKER);
		JCheckBox pvalBox = new JCheckBox(pvalStr);
		pvalBox.setSelected(true);
		pvalBox.setEnabled(false);
		columnBoxes.put(pvalStr, pvalBox);

		listModel.addElement(chrBox);
		listModel.addElement(posBox);
		listModel.addElement(pvalBox);

		for (int i = 0; i < columns.length; i++) {
			if (i == links[ManhattanPlot.CHR_LINKER])
				continue;
			if (i == links[ManhattanPlot.POS_LINKER])
				continue;
			if (i == links[ManhattanPlot.PVAL_LINKER])
				continue;

			JCheckBox box = new JCheckBox(columns[i]);
			if (i == links[ManhattanPlot.SNP_LINKER]) {
				box.setSelected(true);
			}
			columnBoxes.put(columns[i], box);
			listModel.addElement(box);
		}

		list.setModel(listModel);
		repaint();
	}

	protected class CellRenderer implements ListCellRenderer<JCheckBox> {
		public Component getListCellRendererComponent(
																									JList<? extends JCheckBox> list, JCheckBox value,
																									int index,
																									boolean isSelected, boolean cellHasFocus) {
			JCheckBox checkbox = value;

			checkbox.setBackground(isSelected ? list.getSelectionBackground()
																				: list.getBackground());
			checkbox.setForeground(isSelected ? list.getSelectionForeground()
																				: list.getForeground());
			checkbox.setEnabled(value.isEnabled());
			checkbox.setFont(getFont());
			checkbox.setFocusPainted(false);
			checkbox.setBorderPainted(true);
			checkbox.setBorder(isSelected ? UIManager
																							 .getBorder("List.focusCellHighlightBorder")
																		: noFocusBorder);
			return checkbox;
		}
	}

	public int getCloseCode() {
		return closeCode;
	}

	public boolean[] getSelectedChrs() {
		boolean[] sel = new boolean[27];
		for (int i = 0; i < sel.length; i++) {
			sel[i] = chrBoxes.get(i).isSelected();
		}
		return sel;
	}

	public double getPValueThreshold() {
		return chkPValFilter.isSelected() ? (double) spinner.getValue() : Double.NaN;
	}

	public String[] getSelectedColumns() {
		ArrayList<String> cols = new ArrayList<>();
		for (Entry<String, JCheckBox> ent : columnBoxes.entrySet()) {
			if (ent.getValue().isSelected()) {
				cols.add(ent.getKey());
			}
		}
		return cols.toArray(new String[cols.size()]);
	}

	public DataFile getDataFile() {
		return dataFile;
	}

}
