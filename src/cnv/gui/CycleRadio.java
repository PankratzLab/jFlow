package cnv.gui;

import javax.swing.*;
import java.awt.event.*;

public class CycleRadio extends AbstractAction {
	public static final long serialVersionUID = 4L;

	private JRadioButton[] radioButtons;
	private int direction;

	public CycleRadio(JRadioButton[] radioButtons, int direction) {
		super("ChangeType");
		this.radioButtons = radioButtons;
		this.direction = direction;
	}

	public void actionPerformed(ActionEvent e) {
		int index = 0;
		for (int i = 0; i<radioButtons.length; i++) {
			if (radioButtons[i].isSelected()) {
				index = i;
			}
		}
		index += direction;
		if (index<0) {
			index = 0;
		}
		if (index>radioButtons.length-1) {
			index = radioButtons.length-1;
		}
		radioButtons[index].setSelected(true);
	}
}
