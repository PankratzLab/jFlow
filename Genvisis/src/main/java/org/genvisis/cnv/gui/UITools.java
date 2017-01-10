package org.genvisis.cnv.gui;

import java.awt.Component;
import java.awt.Dimension;
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
		comp.setLocation(dim.width / 2 - comp.getSize().width / 2,
		                 dim.height / 2 - comp.getSize().height / 2);
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
		c.setPreferredSize(new Dimension(capValue(width, dim.getWidth()), capValue(height, dim.getHeight())));
	}

	/**
	 * @return the capped value of {@code dim} such that it is falls in the range of [0.1 * cap, 0.9 *
	 *         cap].
	 */
	private static int capValue(int dim, double cap) {
		int min = (int) (0.1 * cap);
		int max = (int) (0.9 * cap);
		return dim < min ? min : Math.min(dim,  max);
	}
}
