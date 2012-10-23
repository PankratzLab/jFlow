package cnv.gui;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JTextField;

import common.Grafik;

public class RegionNavigator extends JPanel implements ActionListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private JTextField textField;
	private String region;
	JButton firstButton, leftButton, rightButton, lastButton;

	/**
	 * Create the panel.
	 */
	public RegionNavigator(String myRegion) {
		region = myRegion;

		firstButton = new JButton(Grafik.getImageIcon("images/firstLast/dFirst.gif", true));
		firstButton.addActionListener(this);
		add(firstButton);

		leftButton = new JButton(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		leftButton.addActionListener(this);
		add(leftButton);

		textField = new JTextField(region);
		textField.setFont(new Font("Tahoma", Font.PLAIN, 14));
		textField.addActionListener(this);
		add(textField);
		textField.setColumns(20);

		rightButton = new JButton(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		rightButton.addActionListener(this);
		add(rightButton);

		lastButton = new JButton(Grafik.getImageIcon("images/firstLast/dLast.gif", true));
		lastButton.addActionListener(this);
		add(lastButton);
	}

	public String getRegion() {
		return textField.getText();
	}

	public void setRegion(String myRegion) {
		region = myRegion;
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		JComponent source = (JComponent) arg0.getSource();
		if (source.equals(firstButton)) {
			System.out.println("First button pressed");
		} else if (source.equals(leftButton)) {
			System.out.println("Left button pressed");
		} else if (source.equals(rightButton)) {
			System.out.println("Right button pressed");
		} else if (source.equals(lastButton)) {
			System.out.println("Last button pressed");
		} else if (source.equals(textField)) {
			System.out.println("Changed region to " + textField.getText());
			region = textField.getText();
			this.getParent().dispatchEvent(arg0);
		}
	}
}
