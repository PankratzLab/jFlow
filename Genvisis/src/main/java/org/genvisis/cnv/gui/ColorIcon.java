package org.genvisis.cnv.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Insets;

import javax.swing.Icon;

public class ColorIcon extends Component implements Icon {
  public static final long serialVersionUID = 1L;
  private final int width;
  private final int height;
  private Color color;
  private Color border;
  private final Insets insets;

  public ColorIcon() {
    this(32, 16);
  }

  public ColorIcon(int width, int height) {
    this(width, height, Color.black);
  }

  public ColorIcon(int width, int height, Color c) {
    this.width = width;
    this.height = height;

    color = c;
    border = Color.black;
    insets = new Insets(1, 1, 1, 1);
  }

  public Color getColor() {
    return color;
  }

  @Override
  public int getIconHeight() {
    return height;
  }

  @Override
  public int getIconWidth() {
    return width;
  }

  @Override
  public void paintIcon(Component c, Graphics g, int x, int y) {
    g.setColor(border);
    g.drawRect(x, y, width - 1, height - 2);

    x += insets.left;
    y += insets.top;

    int w = width - insets.left - insets.right;
    int h = height - insets.top - insets.bottom - 1;

    g.setColor(color);
    g.fillRect(x, y, w, h);
  }

  public void setBorderColor(Color c) {
    border = c;
  }

  public void setColor(Color c) {
    color = c;
  }
}
