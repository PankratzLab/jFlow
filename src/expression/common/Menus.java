package common;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.event.*;

public class Menus {
	public static JMenu createMenu(JMenuBar bar, String name, char mnemonic, MenuListener listener) {
		JMenu menu = (JMenu)bar.add(new JMenu(name));

		if (mnemonic!=-1) {
			menu.setMnemonic(mnemonic);
		}
		if (listener!=null) {
			menu.addMenuListener(listener);
		}
		return menu;
	}

	public static JMenuItem createMenuItem(JMenu menu, String name, int mnemonic, String tooltip, int accelerator, boolean checkbox, ActionListener listener) {
		JMenuItem menuItem = checkbox?(JCheckBoxMenuItem)menu.add(new JCheckBoxMenuItem(name)):(JMenuItem)menu.add(new JMenuItem(name));
		if (mnemonic!=-1) {
			menuItem.setMnemonic(mnemonic);
		}
		if (!tooltip.equals("")) {
			menuItem.setToolTipText(tooltip);
		}
		if (accelerator!=-1) {
			menuItem.setAccelerator(KeyStroke.getKeyStroke(accelerator, ActionEvent.CTRL_MASK));
		}
		if (listener!=null) {
			menuItem.addActionListener(listener);
		}
		if (checkbox) {
			((JCheckBoxMenuItem)menuItem).setState(false);
		}
		return menuItem;
	}

}
