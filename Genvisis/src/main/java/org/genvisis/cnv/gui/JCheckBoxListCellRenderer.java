package org.genvisis.cnv.gui;

import java.awt.Component;
import java.awt.Font;

import javax.swing.JCheckBox;
import javax.swing.JList;
import javax.swing.ListCellRenderer;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;

public class JCheckBoxListCellRenderer implements ListCellRenderer<JCheckBox> {

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
