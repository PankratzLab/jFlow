package cnv.gui;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;

public class CompConfig extends JPanel {

	/**
	 * Create the panel.
	 */
	public CompConfig() {
		setLayout(null);

		JLabel lblDisplayMode = new JLabel("Display Mode");
		lblDisplayMode.setBounds(29, 11, 63, 14);
		add(lblDisplayMode);

		JComboBox comboBox = new JComboBox();
		comboBox.setModel(new DefaultComboBoxModel(new String[] { "Full", "Abbreviated" }));
		comboBox.setBounds(121, 8, 84, 20);
		add(comboBox);

		JLabel lblProbes = new JLabel("Probes");
		lblProbes.setBounds(101, 39, 33, 14);
		add(lblProbes);

		JSlider slider = new JSlider();
		slider.setSnapToTicks(true);
		slider.setPaintTicks(true);
		slider.setPaintLabels(true);
		slider.setMaximum(50);
		slider.setBounds(17, 59, 200, 23);
		add(slider);

		JLabel lblMinimumSizekb = new JLabel("Minimum Size (kb)");
		lblMinimumSizekb.setBounds(75, 93, 84, 14);
		add(lblMinimumSizekb);

		JSlider slider_1 = new JSlider();
		slider_1.setMaximum(1000);
		slider_1.setPaintTicks(true);
		slider_1.setPaintLabels(true);
		slider_1.setBounds(17, 118, 200, 23);
		add(slider_1);

		JLabel lblQualityScore = new JLabel("Quality Score");
		lblQualityScore.setBounds(85, 152, 64, 14);
		add(lblQualityScore);

		JSlider slider_2 = new JSlider();
		slider_2.setPaintTicks(true);
		slider_2.setPaintLabels(true);
		slider_2.setBounds(17, 177, 200, 23);
		add(slider_2);

	}
}
