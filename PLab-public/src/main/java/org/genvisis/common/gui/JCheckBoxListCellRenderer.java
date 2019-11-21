package org.genvisis.common.gui;

import java.awt.Component;
import java.awt.Font;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JCheckBox;
import javax.swing.JList;
import javax.swing.ListCellRenderer;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;

public class JCheckBoxListCellRenderer implements ListCellRenderer<JCheckBox> {

  public static void setupList(JList<JCheckBox> list) {
    list.addMouseListener(new MouseAdapter() {

      int prevIndex = -1;

      public void mousePressed(MouseEvent e) {
        int index = list.locationToIndex(e.getPoint());
        if (index != -1) {
          JCheckBox checkbox = (JCheckBox) list.getModel().getElementAt(index);
          if (checkbox.isEnabled()) {
            boolean cl = e.getClickCount() >= 2;
            boolean li = prevIndex != -1 && index == prevIndex;
            if (cl || li) {
              checkbox.setSelected(!checkbox.isSelected());
            }
            prevIndex = index;
          }
          list.repaint();
        }
      }
    });
  }

  private static final Border BORDER = new EmptyBorder(1, 1, 1, 1);

  private final Font font;

  public JCheckBoxListCellRenderer() {
    this.font = null;
  }

  public JCheckBoxListCellRenderer(Font font) {
    this.font = font;
  }

  public Component getListCellRendererComponent(JList<? extends JCheckBox> list, JCheckBox value,
                                                int index, boolean isSelected,
                                                boolean cellHasFocus) {
    JCheckBox checkbox = value;

    checkbox.setBackground(isSelected ? list.getSelectionBackground() : list.getBackground());
    checkbox.setForeground(isSelected ? list.getSelectionForeground() : list.getForeground());
    checkbox.setEnabled(value.isEnabled());
    if (font != null) {
      checkbox.setFont(font);
    }
    checkbox.setFocusPainted(false);
    checkbox.setBorderPainted(true);
    checkbox.setBorder(isSelected ? UIManager.getBorder("List.focusCellHighlightBorder") : BORDER);
    return checkbox;
  }

}
