package cnv.gui;

import java.awt.Font;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenu;

import common.Files;
import common.ext;

/**
 * @author lane0212 For creating a little menu that has actions associated with each existing file
 */
public abstract class FileActionMenu {

	private String[] actionFiles;
	private JMenu actionMenu;
	private ArrayList<String> existingFiles;
	private HashMap<String, JCheckBoxMenuItem> buttonMap;
	private Hashtable<String, Integer> namePathMap;
	private String menuTitle;
	private Font font = new Font("Arial", 0, 12);

	public FileActionMenu(String menuTitle, String[] actionFiles) {
		super();
		this.menuTitle = menuTitle;
		this.actionFiles = actionFiles;
		init();
	}

	public ItemListener getListener() {
		return null;
	}

	private void init() {
		this.existingFiles = new ArrayList<String>();
		this.buttonMap = new HashMap<String, JCheckBoxMenuItem>();
		this.namePathMap = new Hashtable<String, Integer>();
		if (actionFiles != null) {
			for (int i = 0; i < actionFiles.length; i++) {
				if (Files.exists(actionFiles[i])) {
					existingFiles.add(actionFiles[i]);
					String name = ext.rootOf(actionFiles[i]);
					namePathMap.put(name, i);
					JCheckBoxMenuItem boxMenuItem = new JCheckBoxMenuItem(name);
					buttonMap.put(name, boxMenuItem);
				}
			}
		}
		
	}

	public void developMenu() {
		this.actionMenu = new JMenu(menuTitle);
		ButtonGroup buttonGroup = new ButtonGroup();
		JCheckBoxMenuItem blankBox = new JCheckBoxMenuItem("None");
		blankBox.addItemListener(getListener());
		blankBox.setFont(font);
		blankBox.setSelected(true);
		buttonGroup.add(blankBox);
		buttonMap.put("None", blankBox);
		actionMenu.setFont(font);
		actionMenu.add(blankBox);
		for (String key : namePathMap.keySet()) {
			buttonMap.get(key).addItemListener(getListener());
			buttonMap.get(key).setFont(font);
			buttonGroup.add(buttonMap.get(key));
			actionMenu.add(buttonMap.get(key));
		}
	}

	public String[] getActionFiles() {
		return actionFiles;
	}

	public JMenu getActionMenu() {
		return actionMenu;
	}


	public ArrayList<String> getExistingFiles() {
		return existingFiles;
	}

	public HashMap<String, JCheckBoxMenuItem> getButtonMap() {
		return buttonMap;
	}

	public Hashtable<String, Integer> getNamePathMap() {
		return namePathMap;
	}

	public String getMenuTitle() {
		return menuTitle;
	}

	public Font getFont() {
		return font;
	}


}
