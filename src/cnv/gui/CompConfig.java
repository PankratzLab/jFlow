package cnv.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import cnv.var.CNVariant;

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
	public CompConfig() {
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));

		cnvPanel = new CNVPanel();
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

		probesSlider = new JSlider(JSlider.HORIZONTAL, 0, 50, probes);
		probesSlider.setMajorTickSpacing(10);
		probesSlider.setMinorTickSpacing(5);
		probesSlider.setPaintTicks(true);
		probesSlider.setPaintLabels(true);
		probesSlider.addChangeListener(this);
		configPanel.add(probesSlider);

		// Filter by the minimum size in kilobases
		lblMinimumSizekb = new JLabel("Minimum Size in kb (" + minSize + ")");
		lblMinimumSizekb.setAlignmentX(Component.CENTER_ALIGNMENT);
		configPanel.add(lblMinimumSizekb);

		minSizeSlider = new JSlider(JSlider.HORIZONTAL, 0, 1000, minSize);
		minSizeSlider.setMajorTickSpacing(200);
		minSizeSlider.setMinorTickSpacing(100);
		minSizeSlider.setPaintTicks(true);
		minSizeSlider.setPaintLabels(true);
		minSizeSlider.addChangeListener(this);
		configPanel.add(minSizeSlider);

		// Filter by the minimum quality score
		lblQualityScore = new JLabel("Quality Score (" + qualityScore + ")");
		lblQualityScore.setAlignmentX(Component.CENTER_ALIGNMENT);
		configPanel.add(lblQualityScore);

		qualityScoreSlider = new JSlider(0, 100, qualityScore);
		qualityScoreSlider.setMajorTickSpacing(20);
		qualityScoreSlider.setMinorTickSpacing(10);
		qualityScoreSlider.setPaintTicks(true);
		qualityScoreSlider.setPaintLabels(true);
		qualityScoreSlider.addChangeListener(this);
		configPanel.add(qualityScoreSlider);

		// Scale the height of the rectangles
		lblRectangleHeight = new JLabel("Rectangle Height (" + rectangleHeight + ")");
		lblRectangleHeight.setAlignmentX(Component.CENTER_ALIGNMENT);
		configPanel.add(lblRectangleHeight);

		rectangleHeightSlider = new JSlider(1, 25, rectangleHeight);
		rectangleHeightSlider.setMajorTickSpacing(24);
		rectangleHeightSlider.setMinorTickSpacing(1);
		rectangleHeightSlider.setPaintTicks(true);
		rectangleHeightSlider.setPaintLabels(true);
		rectangleHeightSlider.addChangeListener(this);
		configPanel.add(rectangleHeightSlider);

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

	public void setSelectedCNVs(Vector<CNVariant> cnvs) {
		cnvPanel.setCNVs(cnvs);
	}

}

class CNVPanel extends JPanel {
	private static final long serialVersionUID = 1L;
	Vector<CNVariant> selectedCNVs;
	JPanel cnvPanel;
	JLabel iid; // Individual ID
	JLabel fid; // Family ID
	JLabel length; // CNV length
	JLabel copies; // CNV copy number
	JLabel probes; // Number of probes
	JLabel score; // Quality score
	String displayMode;

	public CNVPanel() {
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		add(new JLabel("Selected CNV:"));
		iid = new JLabel();
		fid = new JLabel();
		length = new JLabel();
		copies = new JLabel();
		probes = new JLabel();
		score = new JLabel();

		// Panel for CNV information
		cnvPanel = new JPanel();
		cnvPanel.setLayout(new GridLayout(6, 2));
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

		add(cnvPanel);
	}

	// Update the fields
	public void setCNVs(Vector<CNVariant> cnvs) {
		selectedCNVs = cnvs;

		if (displayMode.equals("Collapsed")) {
			int maxLength = 0;
			int minLength = Integer.MAX_VALUE;
			int maxCopies = 0;
			int minCopies = Integer.MAX_VALUE;
			int maxProbes = 0;
			int minProbes = Integer.MAX_VALUE;
			double maxScore = 0;
			double minScore = Double.MAX_VALUE;
			Vector<String> iids = new Vector<String>();
			Vector<String> fids = new Vector<String>();

			for (CNVariant cnv : selectedCNVs) {
				String iid = cnv.getIndividualID();
				String fid = cnv.getFamilyID();
				int cnvLength = cnv.getSize();
				int cnvCopies = cnv.getCN();
				int cnvProbes = cnv.getNumMarkers();
				double cnvScore = cnv.getScore();

				// Count iids
				if (!iids.contains(iid)) {
					iids.add(iid);
				}

				// Count fids
				if (!fids.contains(fid)) {
					fids.add(fid);
				}

				// Set min/max length
				if (cnvLength < minLength) {
					minLength = cnvLength;
				}
				if (cnvLength > maxLength) {
					maxLength = cnvLength;
				}

				// Set min/max copies
				if (cnvCopies < minCopies) {
					minCopies = cnvCopies;
				}
				if (cnvCopies > maxCopies) {
					maxCopies = cnvCopies;
				}

				// Set min/max probes
				if (cnvProbes < minProbes) {
					minProbes = cnvProbes;
				}
				if (cnvProbes > maxProbes) {
					maxProbes = cnvProbes;
				}

				// Set min/max score
				if (cnvScore < minScore) {
					minScore = cnvScore;
				}
				if (cnvScore > maxScore) {
					maxScore = cnvScore;
				}
			}

			if (iids.size() > 1) {
				iid.setText("Multiple (" + iids.size() + ")");
			} else {
				iid.setText("" + iids.get(0));
			}

			if (fids.size() > 1) {
				fid.setText("Multiple (" + fids.size() + ")");
			} else {
				fid.setText("" + fids.get(0));
			}

			if (minLength != maxLength) {
				length.setText(minLength + "-" + maxLength);
			} else {
				length.setText("" + minLength);
			}

			if (minCopies != maxCopies) {
				copies.setText(minCopies + "-" + maxCopies);
			} else {
				copies.setText("" + minCopies);
			}

			if (minProbes != maxProbes) {
				probes.setText(minProbes + "-" + maxProbes);
			} else {
				probes.setText("" + minProbes);
			}

			if (minScore != maxScore) {
				score.setText(minScore + "-" + maxScore);
			} else {
				score.setText("" + minScore);
			}
		} else {
			CNVariant selectedCNV = selectedCNVs.get(0);
			iid.setText(selectedCNV.getIndividualID());
			fid.setText(selectedCNV.getFamilyID());
			length.setText("" + selectedCNV.getSize());
			copies.setText("" + selectedCNV.getCN());
			probes.setText("" + selectedCNV.getNumMarkers());
			score.setText("" + selectedCNV.getScore());
		}
		cnvPanel.repaint();

	}

	public void setDisplayMode(String mode) {
		displayMode = mode;
	}
}