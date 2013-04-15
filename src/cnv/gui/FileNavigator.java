package cnv.gui;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JPanel;

import common.Grafik;

public class FileNavigator extends JPanel implements ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Create the panel.
	 */
	public FileNavigator(String[] files) {

		JButton leftButton = new JButton(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		add(leftButton);

		JComboBox<String> comboBox = new JComboBox<String>();
		comboBox.setModel(new DefaultComboBoxModel<String>(files));
		comboBox.setFont(new Font("Tahoma", Font.PLAIN, 14));
		add(comboBox);

		JButton rightButton = new JButton(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		add(rightButton);

	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub

	}
}
