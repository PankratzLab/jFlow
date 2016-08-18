/**
 * 
 */
package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.JPanel;

import org.genvisis.cnv.gui.ChromosomeViewer;
import org.genvisis.cnv.var.CNVRectangle;
import org.genvisis.cnv.var.CNVRectangles;
import org.genvisis.common.Positions;
import org.genvisis.filesys.CNVariant;

/**
 * @author Michael Vieths
 * 
 */
public class CompPanel extends JPanel implements MouseListener, MouseMotionListener, MouseWheelListener {
	public static final long serialVersionUID = 1L;

	CompPlot plot;
	float scalingFactor;
	int startBase, endBase;
	ArrayList<CNVariant> selectedCNVs;
	int rectangleHeight = 10;
	int oldMaxY = 0;
	String displayMode;
	ChromosomeViewer chrViewer;

	private CNVRectangles cnvRectangles;
	ArrayList<org.genvisis.cnv.var.CNVRectangle> rectangles;

	public CompPanel(CompPlot cp) {
		plot = cp;
		cnvRectangles = new CNVRectangles();
		addMouseListener(this);
		addMouseMotionListener(this);
		addMouseWheelListener(this);
	}

	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);

		// Find the maximum Y so we can properly set the dimensions
		int maxY = 0;
		for (CNVRectangle cnvRect : rectangles) {
			Rectangle rect = cnvRect.getRect();
			int y = (int) rect.getY();
			if (y > maxY) {
				maxY = y;
			}
		}

		// Only revalidate/change dimensions if they've changed
		if (maxY != oldMaxY) {
			oldMaxY = maxY;
			Dimension d = getPreferredSize();
			d.height = maxY;
			d.width = chrViewer.getWidth();

			// Set the maximum Y so we see all of the rectangles
			setPreferredSize(d);
			// Revalidate so the scroll bar will be updated properly
			revalidate();
		}

		// Render all of the rectangles
		for (CNVRectangle cnvRect : rectangles) {
			Rectangle rect = cnvRect.getRect();
			int x = (int) rect.getX();
			int y = (int) rect.getY();
			int width = (int) rect.getWidth();
			int height = (int) rect.getHeight();

			g.setColor(cnvRect.getCNVColor());

			if (cnvRect.isSelected()) {
				g.fillRect(x, y, width, height);
				// Draw a black border around the selected CNV
				g.setColor(Color.BLACK);
				g.setPaintMode();
				g.drawRect(x, y, width, height);
			} else {
				g.fillRect(x, y, width, height);
			}

			// In Collapsed mode, draw a 'xN' to indicate that there are N CNVs associated with this rectangle
			if (displayMode.equals("Collapsed") && cnvRect.getCNVs().size() > 1) {
				String numCNVs = "x" + cnvRect.getCNVs().size();
				// If there are 0 copies it may end up black, so don't draw black-on-black
				if (cnvRect.getCNV().getCN() == 0) {
					g.setColor(Color.WHITE);
				} else {
					g.setColor(Color.BLACK);
				}
				g.setPaintMode();
				Font font = getFont().deriveFont((float) rectangleHeight);
				g.setFont(font);
				g.drawString(numCNVs, x, y + rectangleHeight);
			}
		}
	}

	/**
	 * Use the appropriately arranged array of rectangles for our current display mode
	 * 
	 * @param cnvRects
	 */
	void setCNVRectangles(CNVRectangles cnvRects) {
		cnvRectangles = cnvRects;
		if (displayMode.equals("Pack")) {
			setRectangles(cnvRectangles.getPackedRectangles());
		} else if (displayMode.equals("Collapsed")) {
			setRectangles(cnvRectangles.getCollapsedRectangles());
		} else {
			setRectangles(cnvRectangles.getFullRectangles());
		}
	}

	/**
	 * Set the current rectangles and repaint the window
	 * 
	 * @param rects
	 */
	void setRectangles(ArrayList<CNVRectangle> rects) {
		rectangles = rects;
		repaint();
	}

	/**
	 * Need to know how big the visible window is in bases so we can scale it down to the panel width
	 * 
	 * @param window
	 */
	void setWindow(int start, int end) {
		startBase = start;
		endBase = end;
		int window = endBase - startBase;
		scalingFactor = (float) getWidth() / window;
		// scalingFactor = 1;
		cnvRectangles.setScalingFactor(scalingFactor);
	}

	/**
	 * Return the current scaling factor
	 * 
	 * @return
	 */
	public float getScalingFactor() {
		return scalingFactor;
	}

	/**
	 * Set the height of the rectangles (Configured in CompConfig)
	 * 
	 * @param height
	 */
	void setRectangleHeight(int height) {
		rectangleHeight = height;
		cnvRectangles.setRectangleHeight(rectangleHeight);
	}

	public void setChromosomeViewer(ChromosomeViewer viewer) {
		chrViewer = viewer;
	}

	/**
	 * Set the display mode (Configured in CompConfig)
	 * 
	 * Also retrieves the correct arrangement of rectangles
	 * 
	 * @param dm
	 */
	public void setDisplayMode(String dm) {
		displayMode = dm;
		if (displayMode.equals("Pack")) {
			setRectangles(cnvRectangles.getPackedRectangles());
		} else if (displayMode.equals("Collapsed")) {
			setRectangles(cnvRectangles.getCollapsedRectangles());
		} else {
			setRectangles(cnvRectangles.getFullRectangles());
		}
		selectedCNVs = null;
	}

	/**
	 * When clicking on the panel, mark any rectangles under the mouse as selected
	 */
	@Override
	public void mouseClicked(MouseEvent e) {
		for (CNVRectangle cnvRect : rectangles) {
			Rectangle rect = cnvRect.getRect();
			if (rect.contains(e.getPoint())) {
				ArrayList<CNVariant> currentCNVs = cnvRect.getCNVs();
				firePropertyChange("selectedCNV", selectedCNVs, currentCNVs);
				selectedCNVs = currentCNVs;
				cnvRect.setSelected(true);
			} else {
				cnvRect.setSelected(false);
			}

		}
		repaint();
	}

	/**
	 * Check to see if the mouse cursor is over a rectangle and display an informational popup
	 */
	@Override
	public void mouseMoved(MouseEvent e) {
		setToolTipText(null);

		for (CNVRectangle cnvRect : rectangles) {
			Rectangle rect = cnvRect.getRect();
			if (rect.contains(e.getPoint())) {
				ArrayList<CNVariant> currentCNVs = cnvRect.getCNVs();
				CNVariant cnv = cnvRect.getCNV();
				String toolTipText = "<html>IID: " + cnv.getIndividualID() + "<br/>FID: " + cnv.getFamilyID() + "<br/>Length: " + cnv.getSize() + "<br/>Copies: " + cnv.getCN() + "<br/>Probes: " + cnv.getNumMarkers() + "<br/>Score: " + cnv.getScore();
				// Add an indication of how many other CNVs are in this collapsed view
				if (currentCNVs.size() > 1) {
					toolTipText += "<br/>Plus " + (currentCNVs.size() - 1) + " others</html>";
				} else {
					toolTipText += "</html>";
				}
				setToolTipText(toolTipText);
				break;
			}
		}
	}

	double clickStart;

	@Override
	public void mousePressed(MouseEvent e) {
		// Not doing anything on mouse press
		clickStart = e.getPoint().getX();
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		// Not doing anything on mouse release
	}

	@Override
	public void mouseEntered(MouseEvent e) {
		// Not doing anything on mouse enter
	}

	@Override
	public void mouseExited(MouseEvent e) {
		// Not doing anything on mouse exit
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		int[] newLocation = plot.getCPLocation();
		int chromosomeLength = Positions.CHROMOSOME_LENGTHS_B36_HG18[newLocation[0]];
		int x = e.getPoint().x;
		double diff = (int) (clickStart - x);

		double moveRatio = (diff / getWidth());
		double window = (endBase - startBase);
		double offset = window * moveRatio;

		clickStart = x;

		int loc1 = (int) (newLocation[1] + offset);
		int loc2 = (int) (newLocation[2] + offset);

		// Only change the window if we're still in a valid range
		if (offset < 0) {
			if (loc1 >= 0) {
				newLocation[1] = loc1;
				newLocation[2] = loc2;
			}
		} else {
			if (loc2 <= chromosomeLength) {
				newLocation[1] = loc1;
				newLocation[2] = loc2;
			}
		}
		plot.setCPLocation(newLocation);
	}

	public void mouseWheelMoved(MouseWheelEvent e) {
		int[] newLocation = plot.getCPLocation();
		int[] oldLocation = Arrays.copyOf(newLocation, newLocation.length);
//		int chromosomeLength = Positions.CHROMOSOME_LENGTHS_B36_HG18[newLocation[0]]; // TODO make this build specific
		int chromosomeLength = Positions.CHROMOSOME_LENGTHS_B37_HG19[newLocation[0]];

		int rotation = e.getWheelRotation();
		double width = endBase - startBase;

		// Get the X position of the mouse
		int mouseX = e.getX();

		// Figure out the relative position
		double percentage = (double) mouseX / (double) getWidth();

		// Figure out which base we're moused over
		int mouseBase = startBase + (int) ((double) mouseX / scalingFactor);

		if (rotation > 0) {
			// Zoom out 10%
			double newWidth = width / 0.9;
			// Figure out how much of this width comes before mouseX
			double leftBases = (newWidth - 1) * percentage;
			// Figure out how much of this width comes after mouseX
			double rightBases = newWidth - leftBases - 1;

			int loc1 = (int) (mouseBase - leftBases);
			int loc2 = (int) (mouseBase + rightBases);

			// Zoom out
			if (loc1 < 0) {
				newLocation[1] = 0;
			} else {
				newLocation[1] = loc1;
			}

			if (loc2 > chromosomeLength) {
				newLocation[2] = chromosomeLength;
			} else {
				newLocation[2] = loc2;
			}
		} else if (rotation < 0) {
			// Zoom in 10%
			double newWidth = width * 0.9;
			// Figure out how much of this width comes before mouseX
			double leftBases = (newWidth - 1) * percentage;
			// Figure out how much of this width comes after mouseX
			double rightBases = newWidth - leftBases - 1;

			int loc1 = (int) (mouseBase - leftBases);
			int loc2 = (int) (mouseBase + rightBases);

			// Zoom in
			if ((loc2 - loc1) < 100) {
				// Don't zoom in further than a 100bp window.
				// Otherwise we can get caught in a position where the width never changes and we can't zoom back out again.
			} else {
				newLocation[1] = loc1;
				newLocation[2] = loc2;
			}
		}

		// Don't reset the location if it hasn't changed
		if (!Arrays.equals(oldLocation, newLocation)) {
			plot.setCPLocation(newLocation);
		}
	}
}
