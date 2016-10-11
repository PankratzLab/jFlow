package org.genvisis.common;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;
import javax.swing.event.MenuListener;

public class Menus {
	public static JMenu createMenu(JMenuBar bar, String name, char mnemonic, MenuListener listener) {
		JMenu menu = bar.add(new JMenu(name));

		if (mnemonic != (char) -1) {
			menu.setMnemonic(mnemonic);
		}
		if (listener != null) {
			menu.addMenuListener(listener);
		}
		return menu;
	}

	public static JMenuItem createMenuItem(	JMenu menu, String name, int mnemonic, String tooltip,
																					int accelerator, boolean checkbox,
																					ActionListener listener) {
		JMenuItem menuItem = checkbox	? (JCheckBoxMenuItem) menu.add(new JCheckBoxMenuItem(name))
																	: (JMenuItem) menu.add(new JMenuItem(name));
		if (mnemonic != -1) {
			menuItem.setMnemonic(mnemonic);
		}
		if (!tooltip.equals("")) {
			menuItem.setToolTipText(tooltip);
		}
		if (accelerator != -1) {
			menuItem.setAccelerator(KeyStroke.getKeyStroke(accelerator, ActionEvent.CTRL_MASK));
		}
		if (listener != null) {
			menuItem.addActionListener(listener);
		}
		if (checkbox) {
			((JCheckBoxMenuItem) menuItem).setState(false);
		}
		return menuItem;
	}

}
