/**
 *
 */
package org.genvisis.cnv.gui;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Hashtable;

import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import org.genvisis.cnv.plots.AbstractPanel;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.HashVec;

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
      new Color(100, 149, 237), new Color(72, 61, 139), new Color(106, 90, 205),
      new Color(123, 104, 238), new Color(132, 112, 255), new Color(0, 0, 205),
      new Color(65, 105, 225), new Color(0, 0, 255), new Color(30, 144, 255),
      new Color(0, 191, 255), new Color(135, 206, 250), new Color(135, 206, 250),
      new Color(70, 130, 180), new Color(176, 196, 222), new Color(173, 216, 230),
      new Color(176, 224, 230), new Color(175, 238, 238), new Color(0, 206, 209),
      new Color(72, 209, 204), new Color(64, 224, 208), new Color(0, 255, 255),
      new Color(224, 255, 255)};

  private static final Color BACKGROUND_COLOR = Color.WHITE;

  private final Color[] colorScheme;
  private final JPanel classVariablesPanel;
  public JPanel classValuesPanel;
  private int currentClass;
  private JRadioButton[] classRadioButtons;

  private SampleData sampleData;

  private AbstractPanel sisterPanel;

  private Hashtable<String, String> disabledClassValues;
  final ItemListener defaultClassListener = new ItemListener() {
    @Override
    public void itemStateChanged(ItemEvent ie) {
      JRadioButton jrb = (JRadioButton) ie.getItem();
      if (jrb.isSelected()) {
        if (sampleData != null) {
          for (byte i = 0; i < sampleData.getNumClasses(); i++) {
            if (jrb.getText().equals(sampleData.getClassName(i))) {
              currentClass = i;
              getSisterPanel().paintAgain();
            }
          }
        } else if (classRadioButtons != null && classRadioButtons.length > 0) {
          for (int i = 0; i < classRadioButtons.length; i++) {
            if (classRadioButtons[i] != null && classRadioButtons[i].equals(jrb)) {
              currentClass = i;
              getSisterPanel().paintAgain();
            }
          }
        }
      }
    }
  };

  public ColorKeyPanel(SampleData newSampleData, AbstractPanel newSisterPanel) {
    this(newSampleData, newSisterPanel, DEFAULT_COLORS);
  }

  public ColorKeyPanel(SampleData newSampleData, AbstractPanel newSisterPanel,
      Color[] newColorScheme) {
    this(newSampleData, newSisterPanel, newColorScheme, null, 0);
  }

  public ColorKeyPanel(SampleData newSampleData, AbstractPanel newSisterPanel,
      Color[] newColorScheme, ItemListener listener, int clazz) {
    sampleData = newSampleData;
    setSisterPanel(newSisterPanel);
    colorScheme = newColorScheme;
    currentClass = clazz;

    setLayout(new GridLayout(2, 1));
    classVariablesPanel = new JPanel(new WrapLayout(FlowLayout.CENTER, 0, 0));
    classVariablesPanel.setBackground(BACKGROUND_COLOR);
    updateColorVariablePanel(listener == null ? defaultClassListener : listener);
    add(classVariablesPanel);
    classValuesPanel = new JPanel(new WrapLayout(FlowLayout.CENTER, 0, 0));
    classValuesPanel.setBackground(BACKGROUND_COLOR);
    add(classValuesPanel);
  }

  /**
   * @param label label to color if the key is disabled
   * @param key
   */
  private void colorDisabled(JLabel label, String key) {
    if (disabledClassValues.containsKey(currentClass + "\t" + key)) {
      label.setForeground(Color.RED);
      label.setFont(new Font("Arial", Font.ITALIC, 14));
    } else {
      label.setForeground(Color.BLACK);
      label.setFont(new Font("Arial", 0, 14));
    }
  }

  public JRadioButton[] getClassRadioButtons() {
    return classRadioButtons;
  }

  public JPanel getClassValuesPanel() {
    return classValuesPanel;
  }

  public JPanel getClassVariablesPanel() {
    return classVariablesPanel;
  }

  public int getCurrentClass() {
    return currentClass;
  }

  public Hashtable<String, String> getDisabledClassValues() {
    return disabledClassValues;
  }

  public AbstractPanel getSisterPanel() {
    return sisterPanel;
  }

  public void setCurrentClass(int clazz) {
    currentClass = clazz;
  }

  public void setSisterPanel(AbstractPanel sisterPanel) {
    this.sisterPanel = sisterPanel;
  }

  public void updateColorKey(Hashtable<String, String> currentClassUniqueValues) {
    JLabel label, block;
    String[][] colorKeys;
    String[] keys;
    int numBasicClasses, numRegularClasses, numCNVClasses, numPLINKClasses;
    MouseListener mouseListenerForColorKey;

    mouseListenerForColorKey = new MouseListener() {
      @Override
      public void mouseClicked(MouseEvent e) {
        String classValue;

        if (e.getButton() == MouseEvent.BUTTON1) {
          classValue = ((JLabel) e.getSource()).getName();
          if (disabledClassValues.containsKey(classValue)) {
            disabledClassValues.remove(classValue);
          } else {
            disabledClassValues.put(classValue, "");
          }
          getSisterPanel().paintAgain();
        }
      }

      @Override
      public void mouseEntered(MouseEvent e) {}

      @Override
      public void mouseExited(MouseEvent e) {}

      @Override
      public void mousePressed(MouseEvent e) {}

      @Override
      public void mouseReleased(MouseEvent e) {}
    };

    if (sampleData == null) {
      numBasicClasses = 2;
      numRegularClasses = 0;
      numCNVClasses = 0;
      numPLINKClasses = 0;
    } else {
      numBasicClasses = sampleData.getBasicClasses().length;
      numRegularClasses = sampleData.getNumActualClasses();
      numCNVClasses = sampleData.getNumCNVClasses();
      numPLINKClasses = sampleData.getNumPLINKClasses();
    }

    classValuesPanel.removeAll();
    classValuesPanel.repaint();

    label = new JLabel("Color key:");
    label.setFont(new Font("Arial", 0, 14));
    classValuesPanel.add(label);
    if (currentClass < numBasicClasses) {
      colorKeys = SampleData.KEYS_FOR_BASIC_CLASSES[currentClass];
    } else if (currentClass < numBasicClasses + numRegularClasses) {
      colorKeys = sampleData.getActualClassColorKey(currentClass - numBasicClasses);
    } else if (currentClass < numBasicClasses + numRegularClasses + numCNVClasses) {
      colorKeys = new String[][] {{"0", "Not in a CNV"}, {"1", "Homozygous deletion"},
          {"2", "Heterozygous deletion"}, {"4", "Duplication"}, {"5", "Triplication"},};
    } else if (currentClass < numBasicClasses + numRegularClasses + numCNVClasses
        + numPLINKClasses) {
      // colorKeys = new String[][] {{"-1", "Missing"}, {"0", "AA"}, {"1", "AB"}, {"2", "BB"}};
      colorKeys = new String[][] {{"0", "Missing"}, {"1", "AA"}, {"2", "AB"}, {"3", "BB"}};
    } else {
      colorKeys = new String[][] {{"0", "Not in a CNV"}, {"1", "Homozygous deletion"},
          {"2", "Heterozygous deletion"}, {"4", "Duplication"}, {"5", "Triplication"},};
    }
    for (String disabled : disabledClassValues.keySet()) {// so that labels are added even if the
                                                          // class value was not seen
      if (!currentClassUniqueValues.containsKey(disabled)) {
        currentClassUniqueValues.put(disabled.split("\t")[1], "0");
      }
    }

    for (String[] colorKey : colorKeys) {
      block = new JLabel(new ColorIcon(12, 12, colorScheme[Integer.parseInt(colorKey[0])]));
      block.setName(currentClass + "\t" + colorKey[0]);
      block.addMouseListener(mouseListenerForColorKey);
      label = new JLabel(colorKey[1] + " (n=" + (currentClassUniqueValues.containsKey(colorKey[0])
          ? currentClassUniqueValues.get(colorKey[0]) : "0") + ")");
      label.setToolTipText("Remove " + colorKey[1] + " from the display by selecting here");
      label.setName(currentClass + "\t" + colorKey[0]);
      label.addMouseListener(mouseListenerForColorKey);
      currentClassUniqueValues.remove(colorKey[0]);
      colorDisabled(label, colorKey[0]);

      JPanel labelEnclosure = new JPanel();
      labelEnclosure.setBackground(BACKGROUND_COLOR);
      labelEnclosure.add(block);
      labelEnclosure.add(label);
      classValuesPanel.add(labelEnclosure);
    }

    keys = HashVec.getKeys(currentClassUniqueValues);

    for (int i = 0; i < keys.length; i++) {
      if (!keys[i].equals("-1")) {
        block = new JLabel(new ColorIcon(12, 12,
            colorScheme[Math.min(Integer.parseInt(keys[i]), colorScheme.length - 1)]));
        block.setName(currentClass + "\t" + keys[i]);
        block.addMouseListener(mouseListenerForColorKey);
        label = new JLabel((keys[i].equals("0") ? "missing" : keys[i]) + " (n="
            + currentClassUniqueValues.get(keys[i]) + ")");
        label.setName(currentClass + "\t" + keys[i]);
        label.setFont(new Font("Arial", 0, 14));
        label.addMouseListener(mouseListenerForColorKey);
        colorDisabled(label, keys[i]);

        JPanel labelEnclosure = new JPanel();
        labelEnclosure.setBackground(BACKGROUND_COLOR);
        labelEnclosure.add(block);
        labelEnclosure.add(label);
        classValuesPanel.add(labelEnclosure);
      }
    }

    classValuesPanel.validate();
    if (classVariablesPanel != null) {
      classVariablesPanel.validate();
    }
  }

  public void updateColorVariablePanel() {
    updateColorVariablePanel(defaultClassListener);
  }

  public void updateColorVariablePanel(ItemListener classListener) {
    classVariablesPanel.removeAll();

    JLabel label = new JLabel("Color code by:");
    label.setFont(new Font("Arial", 0, 14));
    classVariablesPanel.add(label);


    ButtonGroup classRadio = new ButtonGroup();
    if (sampleData == null) {
      classRadioButtons = new JRadioButton[2];
      classRadioButtons[0] =
          new JRadioButton("All points (link to a SampleData file for more options)", false);
      classRadioButtons[0].setFont(new Font("Arial", 0, 14));
      classRadio.add(classRadioButtons[0]);
      classRadioButtons[0].addItemListener(classListener);
      classRadioButtons[0].setBackground(BACKGROUND_COLOR);
      classVariablesPanel.add(classRadioButtons[0]);
      classRadioButtons[0].setSelected(true);

      classRadioButtons[1] = new JRadioButton("Heat map", false);
      classRadioButtons[1].setFont(new Font("Arial", 0, 14));
      classRadio.add(classRadioButtons[1]);
      classRadioButtons[1].addItemListener(classListener);
      classRadioButtons[1].setBackground(BACKGROUND_COLOR);
      classVariablesPanel.add(classRadioButtons[1]);
      // classRadioButtons[0].setSelected(true);
    } else {
      classRadioButtons = new JRadioButton[sampleData.getNumClasses()];
      for (int i = 0; i < sampleData.getNumClasses(); i++) {
        classRadioButtons[i] = new JRadioButton(sampleData.getClassName(i), false);
        classRadioButtons[i].setFont(new Font("Arial", 0, 14));
        classRadio.add(classRadioButtons[i]);
        classRadioButtons[i].addItemListener(classListener);
        classRadioButtons[i].setBackground(BACKGROUND_COLOR);
        classVariablesPanel.add(classRadioButtons[i]);
      }
      // classRadioButtons[Math.max(sampleData.getBasicClasses().length-1, 0)].setSelected(true);
      classRadioButtons[getCurrentClass()].setSelected(true);
    }

    disabledClassValues = new Hashtable<String, String>();

    classVariablesPanel.validate();
    if (classValuesPanel != null) {
      classValuesPanel.validate();
    }
  }

  public void updateSampleData(SampleData sampleData) {
    this.sampleData = sampleData;
  }

}
