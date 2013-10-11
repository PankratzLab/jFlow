package cnv.gui;

import java.awt.Color;
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
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import cnv.filesys.Project;
import cnv.plots.CompPlot;
import cnv.var.CNVariant;
import cnv.var.SampleData;

public class CompConfig extends JPanel implements ChangeListener, ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private String displayMode = "Full";
	private int probes = 5; // Default to 5 probes
	private int minSize = 0; // Default to no minimum size
	private int qualityScore = 10; // Default to a score of 10
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
	JComboBox<CNVariant> cnvList;
	JLabel cnvListLabel;
	JButton trailerButton;
	LaunchAction launchTrailer;
	CompPlot compPlot;

	public CNVPanel(CompPlot cp) {
		compPlot = cp;
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		add(new JLabel("Selected CNV:"));
		iid = new JLabel();
		fid = new JLabel();
		length = new JLabel();
		copies = new JLabel();
		probes = new JLabel();
		score = new JLabel();
		trailerButton = new JButton("To Trailer");

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

		cnvListLabel = new JLabel("");
		cnvPanel.add(cnvListLabel);

		add(cnvPanel);

		// Link off to Trailer
		add(trailerButton);
		trailerButton.setEnabled(false);
	}

	// Update the fields
	public void setCNVs(ArrayList<CNVariant> cnvs) {
		selectedCNVs = cnvs;

		// In collapsed mode, if there are multiple CNVs associated, add a combo box that lets you select which CNV to look at
		if (displayMode.equals("Collapsed")) {
			// Clear the combo box
			if (cnvList != null) {
				cnvPanel.remove(cnvList);
				cnvListLabel.setText("");
				cnvPanel.validate();
			}

			// Repopulate and add back the combo box if needed
			if (selectedCNVs.size() > 1) {
				cnvList = new JComboBox<CNVariant>();
				cnvList.addActionListener(this);
				CNVariant[] cnvArray = selectedCNVs.toArray(new CNVariant[selectedCNVs.size()]);
				cnvList.setModel(new DefaultComboBoxModel<CNVariant>(cnvArray));
				cnvListLabel.setText("CNVs:");
				cnvPanel.add(cnvListLabel);
				cnvPanel.add(cnvList);
				selectedCNV = (CNVariant) cnvList.getSelectedItem();
			} else {
				selectedCNV = selectedCNVs.get(0);
			}

			setCNVText();
		} else {
			// There's only one CNV, so select it
			selectedCNV = selectedCNVs.get(0);
			setCNVText();
		}

		// Don't enable the button if there aren't any CNVs selected
		if (selectedCNVs.size() > 0) {
			int[] location = compPlot.getCPLocation();
			int window = Integer.parseInt(compPlot.getProject().getProperty(Project.WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER));
			String markerPosition = "chr" + location[0] + ":" + (selectedCNV.getStart() - window) + "-" + (selectedCNV.getStop() + window);
			// System.out.println("Window is " + window + " markerPosition is " + markerPosition);
			String trailerID = selectedCNV.getFamilyID() + "\t" + selectedCNV.getIndividualID();
			SampleData sampleData = compPlot.getProject().getSampleData(2, true);

			launchTrailer = new LaunchAction(compPlot.getProject(), sampleData.lookup(trailerID)[0], markerPosition, Color.BLACK);
			trailerButton.setAction(launchTrailer);
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

	public void setDisplayMode(String mode) {
		displayMode = mode;
	}

	// Monitor the combobox for changes
	public void actionPerformed(ActionEvent arg0) {
		if (arg0.getSource().equals(cnvList)) {
			oldCNV = selectedCNV;
			selectedCNV = (CNVariant) cnvList.getSelectedItem();
			setCNVText();
		}
	}
}