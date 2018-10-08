package org.pankratzlab.shared.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GraphicsConfiguration;
import java.awt.Insets;
import java.awt.Toolkit;

/**
 * Static utility class to help with UI tasks
 */
public final class UITools {

  /**
   * Calls {@link Component#setLocation(java.awt.Point)} to a position that will center the given
   * component on the screen, based the component's dimensions and the size of
   * {@link Toolkit#getScreenSize()}
   */
  public static void centerComponent(Component comp) {
    Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
    Insets scnMax = Toolkit.getDefaultToolkit().getScreenInsets(comp.getGraphicsConfiguration());
    int taskBarSize = scnMax.bottom;
    int totalHeight = dim.height - (taskBarSize > 0 ? taskBarSize : (int) (dim.height * 0.075));

    comp.setLocation(dim.width / 2 - comp.getSize().width / 2,
                     totalHeight / 2 - comp.getSize().height / 2);
  }

  /**
   * Sets the {@link Component#preferredSize()} to the given width and height as a percent of the
   * available {@link Toolkit#getScreenSize()}.
   */
  public static void setSize(Component c, double widthPct, double heightPct) {
    Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
    setSize(c, (int) Math.round(widthPct * dim.getWidth()),
            (int) Math.round(heightPct * dim.getHeight()));
  }

  /**
   * Sets the {@link Component#preferredSize()} to the given width and height, with ranges from
   * [0.1, 0.9] of the {@link Toolkit#getScreenSize()}.
   */
  public static void setSize(Component c, int width, int height) {
    Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
    setSize(c, new Dimension(capValue(width, dim.getWidth()), capValue(height, dim.getHeight())));
  }

  /**
   * Sets the {@link Component#preferredSize()} to the given dimension, reducing height to prevent
   * the window from spawning behind the taskbar.
   */
  public static void setSize(Component c, Dimension d) {
    Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
    GraphicsConfiguration configuration = c.getGraphicsConfiguration();
    Insets scnMax = configuration == null ? null : Toolkit.getDefaultToolkit()
                                                          .getScreenInsets(configuration);
    int taskBarSize = scnMax == null ? 0 : scnMax.bottom;
    int totalHeight = dim.height - (taskBarSize > 0 ? taskBarSize : (int) (dim.height * 0.075));

    int height = d.height >= totalHeight ? totalHeight : d.height;

    c.setPreferredSize(new Dimension(d.width, height));
  }

  /**
   * @return the capped value of {@code dim} such that it is falls in the range of [0.1 * cap, 0.9 *
   *         cap].
   */
  private static int capValue(int dim, double cap) {
    int min = (int) (0.1 * cap);
    int max = (int) (0.9 * cap);
    return dim < min ? min : Math.min(dim, max);
  }
}
