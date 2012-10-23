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

	private String displayMode;
	private int probes;
	private int minSize;
	private int qualityScore;

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
		// lblDisplayMode.setBounds(29, 11, 110, 14);
		dmPanel.add(lblDisplayMode);

		JComboBox<String> comboBox = new JComboBox<String>();
		comboBox.setModel(new DefaultComboBoxModel<String>(new String[] { "Full", "Abbreviated" }));
		comboBox.setBounds(121, 8, 84, 20);
		comboBox.addActionListener(this);
		dmPanel.add(comboBox);
		add(dmPanel);

		// Filter by the number of probes
		JLabel lblProbes = new JLabel("Probes");
		// lblProbes.setBounds(101, 39, 200, 14);
		add(lblProbes);

		probesSlider = new JSlider();
		probesSlider.setSnapToTicks(true);
		probesSlider.setPaintTicks(true);
		probesSlider.setPaintLabels(true);
		probesSlider.setMaximum(50);
		// probesSlider.setBounds(17, 59, 200, 23);
		probesSlider.addChangeListener(this);
		add(probesSlider);

		// Filter by the minimum size in kilobases
		JLabel lblMinimumSizekb = new JLabel("Minimum Size (kb)");
		// lblMinimumSizekb.setBounds(75, 93, 200, 14);
		add(lblMinimumSizekb);

		minSizeSlider = new JSlider();
		minSizeSlider.setMaximum(1000);
		minSizeSlider.setPaintTicks(true);
		minSizeSlider.setPaintLabels(true);
		// minSizeSlider.setBounds(17, 118, 200, 23);
		minSizeSlider.addChangeListener(this);
		add(minSizeSlider);

		// Filter by the minimum quality score
		JLabel lblQualityScore = new JLabel("Quality Score");
		// lblQualityScore.setBounds(85, 152, 200, 14);
		add(lblQualityScore);

		qualityScoreSlider = new JSlider();
		qualityScoreSlider.setPaintTicks(true);
		qualityScoreSlider.setPaintLabels(true);
		// qualityScoreSlider.setBounds(17, 177, 200, 23);
		qualityScoreSlider.addChangeListener(this);
		add(qualityScoreSlider);
	}

	// Monitor the sliders for changes
	public void stateChanged(ChangeEvent arg0) {
		JSlider source = (JSlider) arg0.getSource();
		if (source.equals(probesSlider)) {
			setProbes(source.getValue());
		} else if (source.equals(minSizeSlider)) {
			setMinSize(source.getValue());
		} else if (source.equals(qualityScoreSlider)) {
			setQualityScore(source.getValue());
		}
	}

	// Monitor the combobox for changes
	public void actionPerformed(ActionEvent arg0) {
		@SuppressWarnings("unchecked")
		JComboBox<String> cb = (JComboBox<String>) arg0.getSource();
		setDisplayMode((String) cb.getSelectedItem());
	}

	public String getDisplayMode() {
		return displayMode;
	}

	public void setDisplayMode(String displayMode) {
		this.displayMode = displayMode;
	}

	public int getProbes() {
		return probes;
	}

	public void setProbes(int probes) {
		this.probes = probes;
	}

	public int getMinSize() {
		return minSize;
	}

	public void setMinSize(int minSize) {
		this.minSize = minSize;
	}

	public int getQualityScore() {
		return qualityScore;
	}

	public void setQualityScore(int qualityScore) {
		this.qualityScore = qualityScore;
	}
}
