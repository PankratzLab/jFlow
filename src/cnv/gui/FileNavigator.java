package cnv.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.ListCellRenderer;

public class FileNavigator extends JPanel implements ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Create the panel.
	 */
	public FileNavigator(String[] files, Color[] colors) {

		// JButton leftButton = new JButton(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		// add(leftButton);

		JComboBox<String> comboBox = new JComboBox<String>();
		comboBox.setModel(new DefaultComboBoxModel<String>(files));
		comboBox.setFont(new Font("Tahoma", Font.PLAIN, 14));
		add(comboBox);

		ComboBoxRenderer renderer = new ComboBoxRenderer(comboBox);

		renderer.setColors(colors);
		renderer.setStrings(files);

		comboBox.setRenderer(renderer);

		// JButton rightButton = new JButton(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		// add(rightButton);

	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub

	}
}

class ComboBoxRenderer extends JPanel implements ListCellRenderer<Object> {

	private static final long serialVersionUID = -1L;
	private Color[] colors;
	private String[] strings;

	JPanel textPanel;
	JLabel text;

	public ComboBoxRenderer(JComboBox<String> combo) {

		textPanel = new JPanel();
		textPanel.add(this);
		text = new JLabel();
		text.setOpaque(true);
		text.setFont(combo.getFont());
		textPanel.add(text);
	}

	public void setColors(Color[] col) {
		colors = col;
	}

	public void setStrings(String[] str) {
		strings = str;
	}

	public Color[] getColors() {
		return colors;
	}

	public String[] getStrings() {
		return strings;
	}

	@Override
	public Component getListCellRendererComponent(JList<?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {

		if (isSelected) {
			setBackground(list.getSelectionBackground());
		} else {
			setBackground(Color.WHITE);
		}

		text.setBackground(getBackground());

		text.setText(value.toString());
		if (index > -1) {
			text.setForeground(colors[index]);
		}
		return text;
	}
}