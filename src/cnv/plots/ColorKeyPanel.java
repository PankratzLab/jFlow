/**
 * 
 */
package cnv.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import java.util.*;
import javax.swing.*;

import common.HashVec;

import cnv.gui.ColorIcon;
import cnv.var.SampleData;

/**
 * @author zxu
 *
 */
public class ColorKeyPanel extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static final Color[] DEFAULT_COLORS = {new Color(33, 31, 53), // dark dark
												  new Color(23, 58, 172), // dark blue
												  new Color(201, 30, 10), // deep red
												  new Color(140, 20, 180), // deep purple
												  new Color(33, 87, 0), // dark green
												  new Color(55, 129, 252), // light blue
												  new Color(94, 88, 214), // light purple
												  new Color(189, 243, 61), // light green
												  new Color(217, 109, 194), // pink
												  new Color(0, 0, 128), // ALL KINDS OF BLUES
												  new Color(100, 149, 237),
												  new Color(72, 61, 139),
												  new Color(106, 90, 205),
												  new Color(123, 104, 238),
												  new Color(132, 112, 255),
												  new Color(0, 0, 205),
												  new Color(65, 105, 225),
												  new Color(0, 0, 255),
												  new Color(30, 144, 255),
												  new Color(0, 191, 255),
												  new Color(135, 206, 250),
												  new Color(135, 206, 250),
												  new Color(70, 130, 180),
												  new Color(176, 196, 222),
												  new Color(173, 216, 230),
												  new Color(176, 224, 230),
												  new Color(175, 238, 238),
												  new Color(0, 206, 209),
												  new Color(72, 209, 204),
												  new Color(64, 224, 208),
												  new Color(0, 255, 255),
												  new Color(224, 255, 255)
												 };

	private static final Color BACKGROUND_COLOR = Color.WHITE;

	private Color[] colorScheme;
	private JPanel classVariablesPanel;
	private JPanel classValuesPanel;
	private int currentClass;
	private JRadioButton[] classRadioButtons;
	private SampleData sampleData;
	private AbstractPanel sisterPanel;
	private Hashtable<String, String> disabledClassValues;

	public ColorKeyPanel(SampleData newSampleData, AbstractPanel newSisterPanel) {
		this(newSampleData, newSisterPanel, DEFAULT_COLORS);
	}

	public ColorKeyPanel(SampleData newSampleData, AbstractPanel newSisterPanel, Color[] newColorScheme) {
		this.sampleData = newSampleData;
		this.sisterPanel = newSisterPanel;
		this.colorScheme = newColorScheme;

		setLayout(new GridLayout(2, 1));
		classVariablesPanel = new JPanel();
		JLabel label = new JLabel("Color code by:");
		label.setFont(new Font("Arial", 0, 14));
		classVariablesPanel.add(label);

		ItemListener classListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				JRadioButton jrb = (JRadioButton)ie.getItem();
				if (jrb.isSelected()) {
					for (byte i = 0; i<sampleData.getNumClasses(); i++) {
						if (jrb.getText().equals(sampleData.getClassName(i))) {
							currentClass = i;
							sisterPanel.paintAgain();
						}
					}
				}
			}
		};
		ButtonGroup classRadio = new ButtonGroup();
		classRadioButtons = new JRadioButton[sampleData.getNumClasses()];
		for (int i = 0; i<sampleData.getNumClasses(); i++) {
			classRadioButtons[i] = new JRadioButton(sampleData.getClassName(i), false);
			classRadioButtons[i].setFont(new Font("Arial", 0, 14));
			classRadio.add(classRadioButtons[i]);
			classRadioButtons[i].addItemListener(classListener);
			classRadioButtons[i].setBackground(BACKGROUND_COLOR);
			classVariablesPanel.add(classRadioButtons[i]);
		}
		classVariablesPanel.setBackground(BACKGROUND_COLOR);
		add(classVariablesPanel);
		
		classValuesPanel = new JPanel();
        classValuesPanel.setBackground(BACKGROUND_COLOR);

		add(classValuesPanel);
		classRadioButtons[1].setSelected(true);
		disabledClassValues = new Hashtable<String, String>();
	}

	public void updateColorKey(Hashtable<String,String> currentClassUniqueValues) {
		JLabel label, block;
		String[][] colorKeys;
		String[] keys;
		int numBasicClasses;
		MouseListener mouseListenerForColorKey;
		
		mouseListenerForColorKey = new MouseListener() {
			public void mouseReleased(MouseEvent e) {}
			public void mousePressed(MouseEvent e) {}
			public void mouseExited(MouseEvent e) {}
			public void mouseEntered(MouseEvent e) {}

			public void mouseClicked(MouseEvent e) {
				String classValue;

				if (e.getButton()==MouseEvent.BUTTON1) {
					classValue = ((JLabel)e.getSource()).getName();
					if (disabledClassValues.containsKey(classValue)) {
						disabledClassValues.remove(classValue);
					} else {
						disabledClassValues.put(classValue, "");
					}
					sisterPanel.paintAgain();
				}
			}
		};


		
		numBasicClasses = sampleData.getBasicClasses().length;
		
		classValuesPanel.removeAll();
		classValuesPanel.repaint();
		
		label = new JLabel("Color key:");
		label.setFont(new Font("Arial", 0, 14));
		classValuesPanel.add(label);
		if (currentClass < numBasicClasses) {
			colorKeys = SampleData.KEYS_FOR_BASIC_CLASSES[currentClass];
		} else {		
			colorKeys = sampleData.getActualClassColorKey(currentClass - numBasicClasses);
		}
		for (int i = 0; i<colorKeys.length; i++) {
			block = new JLabel(new ColorIcon(12, 12, colorScheme[Integer.parseInt(colorKeys[i][0])]));
			block.setName(currentClass+"\t"+colorKeys[i][0]);
			block.addMouseListener(mouseListenerForColorKey);
			label = new JLabel(colorKeys[i][1]+" (n="+(currentClassUniqueValues.containsKey(colorKeys[i][0])?currentClassUniqueValues.get(colorKeys[i][0]):"0")+")");
			label.setName(currentClass+"\t"+colorKeys[i][0]);
			label.addMouseListener(mouseListenerForColorKey);
			currentClassUniqueValues.remove(colorKeys[i][0]);
			if (disabledClassValues.containsKey(currentClass+"\t"+colorKeys[i][0])) {
				label.setForeground(Color.RED);
				label.setFont(new Font("Arial", Font.ITALIC, 14));
			} else {
				label.setForeground(Color.BLACK);
				label.setFont(new Font("Arial", 0, 14));
			}
			classValuesPanel.add(block);
			classValuesPanel.add(label);
		}

		keys = HashVec.getKeys(currentClassUniqueValues);

		for (int i = 0; i<keys.length; i++) {
			if (!keys[i].equals("-1")) {
				block = new JLabel(new ColorIcon(12, 12, colorScheme[Integer.parseInt(keys[i])]));
				block.setName(currentClass+"\t"+keys[i]);
				block.addMouseListener(mouseListenerForColorKey);
				label = new JLabel((keys[i].equals("0")?"missing":keys[i])+" (n="+currentClassUniqueValues.get(keys[i])+")");
				label.setName(currentClass+"\t"+keys[i]);
				label.setFont(new Font("Arial", 0, 14));
				label.addMouseListener(mouseListenerForColorKey);
				classValuesPanel.add(block);
				classValuesPanel.add(label);
			}
		}

		classValuesPanel.validate();
	}

	public int getCurrentClass() {
		return currentClass;
	}

	public Hashtable<String, String> getDisabledClassValues() {
		return disabledClassValues;
	}
	
	public JRadioButton[] getClassRadioButtons() {
		return classRadioButtons;
	}
}
