package cnv.gui;

import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class CompConfig extends JPanel implements ChangeListener, ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private String displayMode = "Full";
	private int probes = 5; // Default to 5 probes
	private int minSize = 0; // Default to no minimum size
	private int qualityScore = 10; // Default to a score of 10

	JSlider probesSlider;
	JLabel lblProbes;

	JSlider minSizeSlider;
	JLabel lblMinimumSizekb;

	JSlider qualityScoreSlider;
	JLabel lblQualityScore;

	/**
	 * Create the panel.
	 */
	public CompConfig() {
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));

		JPanel dmPanel = new JPanel(new FlowLayout());

		// Configure whether the display mode is Full (variants displayed are separated and color-coded by file)
		// or abbreviated (CNVs with multiple copies appear as only one block)
		JLabel lblDisplayMode = new JLabel("Display Mode");
		dmPanel.add(lblDisplayMode);

		JComboBox<String> comboBox = new JComboBox<String>();
		comboBox.setModel(new DefaultComboBoxModel<String>(new String[] { "Full", "Abbreviated", "Compressed" }));
		comboBox.addActionListener(this);
		dmPanel.add(comboBox);
		add(dmPanel);

		// Filter by the number of probes
		lblProbes = new JLabel("Probes (" + probes + ")");
		lblProbes.setAlignmentX(Component.CENTER_ALIGNMENT);
		add(lblProbes);

		probesSlider = new JSlider(JSlider.HORIZONTAL, 0, 50, probes);

		probesSlider.setMajorTickSpacing(10);
		probesSlider.setMinorTickSpacing(5);
		probesSlider.setPaintTicks(true);
		probesSlider.setPaintLabels(true);
		probesSlider.addChangeListener(this);
		add(probesSlider);

		// Filter by the minimum size in kilobases
		lblMinimumSizekb = new JLabel("Minimum Size in kb (" + minSize + ")");
		lblMinimumSizekb.setAlignmentX(Component.CENTER_ALIGNMENT);
		add(lblMinimumSizekb);

		minSizeSlider = new JSlider(JSlider.HORIZONTAL, 0, 1000, minSize);

		minSizeSlider.setMajorTickSpacing(200);
		minSizeSlider.setMinorTickSpacing(100);
		minSizeSlider.setPaintTicks(true);
		minSizeSlider.setPaintLabels(true);
		minSizeSlider.addChangeListener(this);
		add(minSizeSlider);

		// Filter by the minimum quality score
		lblQualityScore = new JLabel("Quality Score (" + qualityScore + ")");
		lblQualityScore.setAlignmentX(Component.CENTER_ALIGNMENT);
		add(lblQualityScore);

		qualityScoreSlider = new JSlider(0, 100, qualityScore);

		qualityScoreSlider.setMajorTickSpacing(20);
		qualityScoreSlider.setMinorTickSpacing(10);
		qualityScoreSlider.setPaintTicks(true);
		qualityScoreSlider.setPaintLabels(true);
		qualityScoreSlider.addChangeListener(this);
		add(qualityScoreSlider);
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
		}
	}

	// Monitor the combobox for changes
	public void actionPerformed(ActionEvent arg0) {
		@SuppressWarnings("unchecked")
		JComboBox<String> cb = (JComboBox<String>) arg0.getSource();
		String mode = (String) cb.getSelectedItem();
		firePropertyChange("displayMode", displayMode, mode);
		displayMode = mode;
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

}
