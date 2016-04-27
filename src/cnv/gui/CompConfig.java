package cnv.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import net.miginfocom.swing.MigLayout;
import cnv.filesys.Project;
import cnv.plots.CompPlot;
import cnv.plots.Trailer;
import cnv.var.SampleData;
import filesys.CNVariant;

public class CompConfig extends JPanel implements ChangeListener, ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private String displayMode = "Full";
	private int probes = 0; // Default to 0 probes
	private int minSize = 0; // Default to no minimum size
	private int qualityScore = 0; // Default to a score of 0
	private int rectangleHeight = 10; // Default rectangle height
	CompPlot compPlot;

	JSlider probesSlider;
	JLabel lblProbes;

	JSlider minSizeSlider;
	JLabel lblMinimumSizekb;

	JSlider qualityScoreSlider;
	JLabel lblQualityScore;

	JSlider rectangleHeightSlider;
	JLabel lblRectangleHeight;

	CNVPanel cnvPanel;

	//

	/**
	 * Create the panel.
	 */
	public CompConfig(CompPlot cp) {
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		compPlot = cp;

		cnvPanel = new CNVPanel(cp);
		cnvPanel.setDisplayMode(displayMode);
		add(cnvPanel);

		add(Box.createGlue());

		JPanel configPanel = new JPanel();
		configPanel.setLayout(new BoxLayout(configPanel, BoxLayout.PAGE_AXIS));

		// Configure whether the display mode is Full (variants displayed are separated and color-coded by file)
		// or abbreviated (CNVs with multiple copies appear as only one block)
		JPanel dmPanel = new JPanel(new FlowLayout());
		JLabel lblDisplayMode = new JLabel("Display Mode");
		dmPanel.add(lblDisplayMode);

		JComboBox<String> comboBox = new JComboBox<String>();
		comboBox.setModel(new DefaultComboBoxModel<String>(new String[] { "Full", "Pack", "Collapsed" }));
		comboBox.addActionListener(this);
		dmPanel.add(comboBox);
		configPanel.add(dmPanel);

		// Filter by the number of probes
		lblProbes = new JLabel("Probes (" + probes + ")");
		lblProbes.setAlignmentX(Component.CENTER_ALIGNMENT);
		configPanel.add(lblProbes);

		add(Box.createGlue());

		probesSlider = new JSlider(JSlider.HORIZONTAL, 0, 50, probes);
		probesSlider.setMajorTickSpacing(10);
		probesSlider.setMinorTickSpacing(5);
		probesSlider.setPaintTicks(true);
		probesSlider.setPaintLabels(true);
		probesSlider.addChangeListener(this);
		configPanel.add(probesSlider);

		add(Box.createGlue());

		// Filter by the minimum size in kilobases
		lblMinimumSizekb = new JLabel("Minimum Size in kb (" + minSize + ")");
		lblMinimumSizekb.setAlignmentX(Component.CENTER_ALIGNMENT);
		configPanel.add(lblMinimumSizekb);

		add(Box.createGlue());

		minSizeSlider = new JSlider(JSlider.HORIZONTAL, 0, 1000, minSize);
		minSizeSlider.setMajorTickSpacing(200);
		minSizeSlider.setMinorTickSpacing(100);
		minSizeSlider.setPaintTicks(true);
		minSizeSlider.setPaintLabels(true);
		minSizeSlider.addChangeListener(this);
		configPanel.add(minSizeSlider);

		add(Box.createGlue());

		// Filter by the minimum quality score
		lblQualityScore = new JLabel("Quality Score (" + qualityScore + ")");
		lblQualityScore.setAlignmentX(Component.CENTER_ALIGNMENT);
		configPanel.add(lblQualityScore);

		add(Box.createGlue());

		qualityScoreSlider = new JSlider(0, 100, qualityScore);
		qualityScoreSlider.setMajorTickSpacing(20);
		qualityScoreSlider.setMinorTickSpacing(10);
		qualityScoreSlider.setPaintTicks(true);
		qualityScoreSlider.setPaintLabels(true);
		qualityScoreSlider.addChangeListener(this);
		configPanel.add(qualityScoreSlider);

		add(Box.createGlue());

		// Scale the height of the rectangles
		lblRectangleHeight = new JLabel("Rectangle Height (" + rectangleHeight + ")");
		lblRectangleHeight.setAlignmentX(Component.CENTER_ALIGNMENT);
		configPanel.add(lblRectangleHeight);

		add(Box.createGlue());

		rectangleHeightSlider = new JSlider(1, 25, rectangleHeight);
		rectangleHeightSlider.setMajorTickSpacing(24);
		rectangleHeightSlider.setMinorTickSpacing(1);
		rectangleHeightSlider.setPaintTicks(true);
		rectangleHeightSlider.setPaintLabels(true);
		rectangleHeightSlider.addChangeListener(this);
		configPanel.add(rectangleHeightSlider);

		add(Box.createGlue());

		add(configPanel);
	}

	// Monitor the sliders for changes
	@Override
	public void stateChanged(ChangeEvent arg0) {
		JSlider source = (JSlider) arg0.getSource();
		if (source.equals(probesSlider)) {
			int p = source.getValue();
			this.firePropertyChange("probes", probes, p);
			probes = p;
			lblProbes.setText("Probes (" + probes + ")");
		} else if (source.equals(minSizeSlider)) {
			int ms = source.getValue();
			this.firePropertyChange("minSize", minSize, ms);
			minSize = ms;
			lblMinimumSizekb.setText("Minimum Size in kb (" + minSize + ")");
		} else if (source.equals(qualityScoreSlider)) {
			int qs = source.getValue();
			this.firePropertyChange("qualityScore", qualityScore, qs);
			qualityScore = qs;
			lblQualityScore.setText("Quality Score (" + qualityScore + ")");
		} else if (source.equals(rectangleHeightSlider)) {
			int height = source.getValue();
			this.firePropertyChange("rectangleHeight", rectangleHeight, height);
			rectangleHeight = height;
			lblRectangleHeight.setText("Rectangle Height (" + rectangleHeight + ")");
		}
	}

	// Monitor the combobox for changes
	@Override
	public void actionPerformed(ActionEvent arg0) {
		@SuppressWarnings("unchecked")
		JComboBox<String> cb = (JComboBox<String>) arg0.getSource();
		String mode = (String) cb.getSelectedItem();
		firePropertyChange("displayMode", displayMode, mode);
		displayMode = mode;
		cnvPanel.setDisplayMode(displayMode);
	}

	public String getDisplayMode() {
		return displayMode;
	}

	public int getProbes() {
		return probes;
	}

	public int getMinSize() {
		return minSize;
	}

	public int getQualityScore() {
		return qualityScore;
	}

	public int getRectangleHeight() {
		return rectangleHeight;
	}

	public void setSelectedCNVs(ArrayList<CNVariant> cnvs) {
		cnvPanel.setCNVs(cnvs);
	}

	public CompPlot getPlot() {
		return compPlot;
	}

}

class CNVPanel extends JPanel implements ActionListener {
	private static final long serialVersionUID = 1L;
	ArrayList<CNVariant> selectedCNVs;
	CNVariant selectedCNV;
	CNVariant oldCNV;
	JPanel cnvPanel;
	JLabel iid; // Individual ID
	JLabel fid; // Family ID
	JLabel length; // CNV length
	JLabel copies; // CNV copy number
	JLabel probes; // Number of probes
	JLabel score; // Quality score
	String displayMode;
	JScrollPane cnvPane;
	JComboBox<CNVariant> cnvList;
	JScrollPane cnvScroll;

	JLabel cnvListLabel;
	CNVCheckList checkList;
	JButton selectAll;
	JButton selectNone;
	JButton trailerButton;
	LaunchAction launchTrailer;
	CompPlot compPlot;

	public CNVPanel(CompPlot cp) {
		compPlot = cp;
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
//		setLayout(new MigLayout("debug", "[]", "[][][]"));
		add(new JLabel("Selected CNV:"));
		iid = new JLabel();
		fid = new JLabel();
		length = new JLabel();
		copies = new JLabel();
		probes = new JLabel();
		score = new JLabel();
		selectAll = new JButton("Select All");
		selectNone = new JButton("Select None");
		trailerButton = new JButton("To Trailer");
		trailerButton.setToolTipText("Launch Trailer to examine selected CNVs");

		// Panel for CNV information
		cnvPanel = new JPanel();
		cnvPanel.setLayout(new GridLayout(7, 2));
		cnvPanel.setPreferredSize(new Dimension(250, 120));
		cnvPanel.setMinimumSize(new Dimension(250, 120));
		cnvPanel.setMaximumSize(new Dimension(250, 120));

		// Individual ID
		cnvPanel.add(new JLabel("IID:"));
		cnvPanel.add(iid);

		// Family ID
		cnvPanel.add(new JLabel("FID:"));
		cnvPanel.add(fid);

		// Length
		cnvPanel.add(new JLabel("Length:"));
		cnvPanel.add(length);

		// Copies
		cnvPanel.add(new JLabel("Copies"));
		cnvPanel.add(copies);

		// Probes
		cnvPanel.add(new JLabel("Probes:"));
		cnvPanel.add(probes);

		// Score
		cnvPanel.add(new JLabel("Score:"));
		cnvPanel.add(score);

		cnvListLabel = new JLabel("Select CNVs:");
		cnvPanel.add(cnvListLabel);

//		add(cnvPanel, "cell 0 0");
		add(cnvPanel);

		cnvScroll = new JScrollPane();
//		add(cnvScroll, "cell 0 1, hidemode 3");
		add(cnvScroll);
		cnvScroll.setVisible(false);
		
		JPanel btnPanel = new JPanel();
		btnPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
//		add(selectAll, "cell 0 2, split 3");
		btnPanel.add(selectAll);
		selectAll.setVisible(false);
		selectAll.addActionListener(this);

//		add(selectNone,  "cell 0 2");
		btnPanel.add(selectNone);
		selectNone.setVisible(false);
		selectNone.addActionListener(this);

		// Link off to Trailer
//		add(trailerButton, "cell 0 2");
		btnPanel.add(trailerButton);
		trailerButton.setEnabled(false);
		trailerButton.addActionListener(this);
		
		add(btnPanel);
	}

	// Update the fields
	public void setCNVs(ArrayList<CNVariant> cnvs) {
		selectedCNVs = new ArrayList<CNVariant>(cnvs);

		// In collapsed mode, if there are multiple CNVs associated, add a combo box that lets you select which CNV to look at
		if (displayMode.equals("Collapsed")) {
			if (selectedCNVs.size() > 1) {
				cnvScroll.setPreferredSize(new Dimension(100, 100));
				JPanel checkPanel = new JPanel(new GridLayout(selectedCNVs.size(), 1));
				checkList = new CNVCheckList(selectedCNVs);
				for (JCheckBox checkBox : checkList.checkList) {
					checkPanel.add(checkBox);
				}
				cnvScroll.add(checkPanel);
				cnvScroll.setViewportView(checkPanel);
				cnvScroll.setVisible(true);
				selectAll.setVisible(true);
				selectNone.setVisible(true);
			} else {
				selectAll.setVisible(false);
				selectNone.setVisible(false);
				cnvScroll.setVisible(false);
			}

			selectedCNV = selectedCNVs.get(0);

			setCNVText();
		} else {
			// There's only one CNV, so select it
			selectedCNV = selectedCNVs.get(0);

			// Clear the combo box
			if (cnvList != null) {
				cnvPanel.remove(cnvList);
				cnvListLabel.setText("");
				cnvPanel.validate();
			}

			setCNVText();
		}

		// Don't enable the button if there aren't any CNVs selected
		if (selectedCNVs.size() > 0) {
			trailerButton.setText("To Trailer");
			trailerButton.setIcon(null);
			trailerButton.setEnabled(true);
		} else {
			trailerButton.setEnabled(false);
		}

		cnvPanel.repaint();
	}

	// Update the text with the currently selected CNV
	public void setCNVText() {
		iid.setText(selectedCNV.getIndividualID());
		fid.setText(selectedCNV.getFamilyID());
		length.setText("" + selectedCNV.getSize());
		copies.setText("" + selectedCNV.getCN());
		probes.setText("" + selectedCNV.getNumMarkers());
		score.setText("" + selectedCNV.getScore());
	}

	// Clear the CNV panel
	public void clearCNVText() {
		iid.setText("");
		fid.setText("");
		length.setText("");
		copies.setText("");
		probes.setText("");
		score.setText("");
		trailerButton.setEnabled(false);
		selectAll.setVisible(false);
		selectNone.setVisible(false);
		cnvScroll.setVisible(false);

		if (cnvList != null) {
			cnvPanel.remove(cnvList);
			cnvListLabel.setText("");
			cnvPanel.validate();
		}
	}

	public void setDisplayMode(String mode) {
		displayMode = mode;
		clearCNVText();
		selectedCNV = null;
		selectedCNVs = null;
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		if (arg0.getSource().equals(cnvList)) {
			oldCNV = selectedCNV;
			setCNVText();
		} else if (arg0.getSource().equals(selectAll)) {
			checkList.selectAll();
		} else if (arg0.getSource().equals(selectNone)) {
			checkList.selectNone();
		} else if (arg0.getSource().equals(trailerButton)) {
			// Launch 1 or more instances of Trailer
			Project proj = compPlot.getProject();
			SampleData sampleData = compPlot.getProject().getSampleData(2, true);
			int[] location = compPlot.getCPLocation();
//			int window = Integer.parseInt(compPlot.getProject().getProperty(Project.WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER));
			int window = compPlot.getProject().getProperty(proj.WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER);

			if (selectedCNVs.size() > 1) {
				// More than 4 seems to run out of heap space
				if (checkList.getSelected().size() > 4) {
					String[] options = { "Yes", "No" };
					int answer = JOptionPane.showOptionDialog(null, "Warning - this will launch " + checkList.getSelected().size() + " instances of Trailer\nProceed?", "Warning", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[1]);
					if (answer == JOptionPane.YES_OPTION) {
						for (CNVariant cnv : checkList.getSelected()) {
							String markerPosition = "chr" + location[0] + ":" + (cnv.getStart() - window) + "-" + (cnv.getStop() + window);

							// Strip p or q from the end
							if (markerPosition.endsWith("p") || markerPosition.endsWith("q")) {
								markerPosition = markerPosition.substring(0, markerPosition.length() - 1);
							}

							String trailerID = cnv.getFamilyID() + "\t" + cnv.getIndividualID();

							new Trailer(proj, sampleData.lookup(trailerID)[0], proj.CNV_FILENAMES.getValue(), markerPosition);
						}
					}
				} else if (checkList.getSelected().size() == 0) {
					JOptionPane.showMessageDialog(null, "No CNVs selected");
				} else {
					for (CNVariant cnv : checkList.getSelected()) {
						String markerPosition = "chr" + location[0] + ":" + (cnv.getStart() - window) + "-" + (cnv.getStop() + window);
						// Strip p or q from the end
						if (markerPosition.endsWith("p") || markerPosition.endsWith("q")) {
							markerPosition = markerPosition.substring(0, markerPosition.length() - 1);
						}

						String trailerID = cnv.getFamilyID() + "\t" + cnv.getIndividualID();

						new Trailer(proj, sampleData.lookup(trailerID)[0], proj.CNV_FILENAMES.getValue(), markerPosition);
					}
				}
			} else {
				String markerPosition = "chr" + location[0] + ":" + (selectedCNV.getStart() - window) + "-" + (selectedCNV.getStop() + window);

				// Strip p or q from the end
				if (markerPosition.endsWith("p") || markerPosition.endsWith("q")) {
					markerPosition = markerPosition.substring(0, markerPosition.length() - 1);
				}

				String trailerID = selectedCNV.getFamilyID() + "\t" + selectedCNV.getIndividualID();
				// TODO errors if IDs is null
				String[] ids = sampleData.lookup(trailerID);
				if (ids == null) {
					proj.message("Error - could not find a lookup for individual "+selectedCNV.getFamilyID()+"-"+selectedCNV.getIndividualID()+" in the SampleData file; cannot launch Trailer without knowing which DNA sample to open");
				} else {
					String id = ids[0];
					String[] fns = proj.CNV_FILENAMES.getValue();
					new Trailer(proj, id, fns, markerPosition);
				}
			}
		}
	}
}

/**
 * Create a checklist of variants
 * 
 * @author Michael Vieths
 * 
 */
class CNVCheckList {
	ArrayList<JCheckBox> checkList;
	ArrayList<CNVariant> variants;

	public CNVCheckList(ArrayList<CNVariant> variants) {
		this.variants = variants;
		checkList = new ArrayList<JCheckBox>();
		for (CNVariant cnv : variants) {
			checkList.add(new JCheckBox(cnv.getIndividualID()));
		}
	}

	public void selectAll() {
		for (JCheckBox box : checkList) {
			box.setSelected(true);
		}
	}

	public void selectNone() {
		for (JCheckBox box : checkList) {
			box.setSelected(false);
		}
	}

	ArrayList<CNVariant> getSelected() {
		ArrayList<CNVariant> selected = new ArrayList<CNVariant>();
		for (int i = 0; i < checkList.size(); i++) {
			if (checkList.get(i).isSelected()) {
				selected.add(variants.get(i));
			}
		}
		return selected;
	}
}
