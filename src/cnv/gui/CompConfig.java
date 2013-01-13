package cnv.gui;

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
	private int probes = 25; // Default to 25 probes
	private int minSize = 100;
	private int qualityScore = 50;

	JSlider probesSlider;
	JSlider minSizeSlider;
	JSlider qualityScoreSlider;

	/**
	 * Create the panel.
	 */
	public CompConfig() {
		// setLayout(null);
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

		JPanel dmPanel = new JPanel(new FlowLayout());

		// Configure whether the display mode is Full (variants displayed are separated and color-coded by file)
		// or abbreviated (CNVs with multiple copies appear as only one block)
		JLabel lblDisplayMode = new JLabel("Display Mode");
		dmPanel.add(lblDisplayMode);

		JComboBox<String> comboBox = new JComboBox<String>();
		comboBox.setModel(new DefaultComboBoxModel<String>(new String[] { "Full", "Abbreviated" }));
		comboBox.setBounds(121, 8, 84, 20);
		comboBox.addActionListener(this);
		dmPanel.add(comboBox);
		add(dmPanel);

		// Filter by the number of probes
		JLabel lblProbes = new JLabel("Probes");
		add(lblProbes);

		probesSlider = new JSlider();
		probesSlider.setSnapToTicks(true);
		probesSlider.setPaintTicks(true);
		probesSlider.setPaintLabels(true);
		probesSlider.setMaximum(50);
		probesSlider.setValue(probes);
		probesSlider.addChangeListener(this);
		add(probesSlider);

		// Filter by the minimum size in kilobases
		JLabel lblMinimumSizekb = new JLabel("Minimum Size (kb)");
		add(lblMinimumSizekb);

		minSizeSlider = new JSlider();
		minSizeSlider.setMaximum(1000);
		minSizeSlider.setValue(minSize);
		minSizeSlider.setPaintTicks(true);
		minSizeSlider.setPaintLabels(true);
		minSizeSlider.addChangeListener(this);
		add(minSizeSlider);

		// Filter by the minimum quality score
		JLabel lblQualityScore = new JLabel("Quality Score");
		add(lblQualityScore);

		qualityScoreSlider = new JSlider();
		qualityScoreSlider.setMaximum(100);
		qualityScoreSlider.setValue(qualityScore);
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
		} else if (source.equals(minSizeSlider)) {
			int ms = source.getValue();
			this.firePropertyChange("minSize", minSize, ms);
			minSize = ms;
		} else if (source.equals(qualityScoreSlider)) {
			int qs = source.getValue();
			this.firePropertyChange("qualityScore", qualityScore, qs);
			qualityScore = qs;
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
