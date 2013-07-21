/**
 * 
 */
package cnv.plots;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.Vector;

import javax.swing.JPanel;

import cnv.var.CNVariant;

/**
 * @author Michael Vieths
 * 
 */
public class CompPanel extends JPanel implements MouseListener, MouseMotionListener {
	public static final long serialVersionUID = 1L;
	CNVRectangle[] rectangles;
	float scalingFactor;
	int startBase, endBase;
	Vector<CNVariant> selectedCNVs;
	int rectangleHeight = 10;
	int lowestStart = 0; // Represents the lowest startX of all rectangles in the window
	String displayMode;
	ArrayList<ArrayList<CNVRectangle>> lines = null;

	public CompPanel(CompPlot cp) {
		lines = new ArrayList<ArrayList<CNVRectangle>>();
		addMouseListener(this);
		addMouseMotionListener(this);
	}

	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		if (displayMode.equals("Pack")) {
			// Pack the CNVs so they fill a row
			packRectangles();

			setPreferredSize(new Dimension(800, (lines.size() * rectangleHeight) + lines.size()));

			for (int i = 0; i < lines.size(); i++) {
				ArrayList<CNVRectangle> cnvRects = lines.get(i);
				int y = (i * rectangleHeight) + i;
				for (int j = 0; j < cnvRects.size(); j++) {
					CNVRectangle cnvRect = cnvRects.get(j);
					int width = Math.round(((int) cnvRect.getStopXValue() - (int) cnvRect.getStartXValue()) * scalingFactor);
					int x = Math.round((int) cnvRect.getStartXValue() * scalingFactor);

					// Store the rectangle for later bounds checking
					cnvRect.setRect(x, y, width, rectangleHeight);

					// CNV color is based on number of copies; set during CompPlot.loadCNVs
					g.setColor(cnvRect.getCNVColor());

					if (cnvRect.isSelected()) {
						g.fillRect(x, y, width, rectangleHeight);
						// Draw a 1-pixel wide black rectangle around the selected CNV
						g.setColor(Color.BLACK);
						g.setPaintMode();
						// Subtract 1 from y to ensure we're drawing the bottom line of the highlight rectangle
						g.drawRect(x, y - 1, width, rectangleHeight);
					} else {
						g.fillRect(x, y, width, rectangleHeight);
					}
				}
			}
		} else {
			// We'll need to scale the relative base to the window size
			setPreferredSize(new Dimension(800, (rectangles.length * rectangleHeight) + rectangles.length));

			for (int i = 0; i < rectangles.length; i++) {
				int width = Math.round(((int) rectangles[i].getStopXValue() - (int) rectangles[i].getStartXValue()) * scalingFactor);
				int x = Math.round((int) rectangles[i].getStartXValue() * scalingFactor);

				// Ensure a 1 pixel gap between CNVs
				int y = (i * rectangleHeight) + i;

				// Store the rectangle for later bounds checking
				rectangles[i].setRect(x, y, width, rectangleHeight);

				// CNV color is based on number of copies; set during CompPlot.loadCNVs
				g.setColor(rectangles[i].getCNVColor());

				if (rectangles[i].isSelected()) {
					g.fillRect(x, y, width, rectangleHeight);
					// Draw a 1-pixel wide black rectangle around the selected CNV
					g.setColor(Color.BLACK);
					g.setPaintMode();
					// Subtract 1 from y to ensure we're drawing the bottom line of the highlight rectangle
					g.drawRect(x, y - 1, width, rectangleHeight);
				} else {
					g.fillRect(x, y, width, rectangleHeight);
				}
			}
		}
	}

	private void packRectangles() {
		while (hasUnused()) {
			CNVRectangle cnvRect = getLeftMost(lowestStart);
			if (cnvRect != null) {
				ArrayList<CNVRectangle> newLine = new ArrayList<CNVRectangle>();
				newLine.add(cnvRect);
				cnvRect.setUsed(true);
				do {
					cnvRect = getLeftMost((int) cnvRect.getStopXValue());
					if (cnvRect == null) {
						break;
					} else {
						newLine.add(cnvRect);
						cnvRect.setUsed(true);
					}
				} while (cnvRect != null);
				lines.add(newLine);
			} else {
				break;
			}
		}
	}

	private int getLowestRect() {
		int lowest = 0;
		for (int i = 0; i < rectangles.length; i++) {
			int startX = (int) rectangles[i].getStartXValue();
			if (startX < lowest) {
				lowest = startX;
			}
		}
		return lowest;
	}

	private boolean hasUnused() {
		boolean unused = false;
		for (int i = 0; i < rectangles.length; i++) {
			if (!rectangles[i].isInUse()) {
				unused = true;
				break;
			}
		}
		return unused;
	}

	private void clearUsed() {
		for (int i = 0; i < rectangles.length; i++) {
			rectangles[i].setUsed(false);
		}
	}

	private CNVRectangle getLeftMost(int startX) {
		int lastPosition = startX + 2; // Want a 2-pixel buffer between CNVs
		if (startX == lowestStart) {
			// The first time through we'll be starting at lowestStart so we don't want the offset
			lastPosition = lowestStart;
		}

		CNVRectangle leftMostCNV = null;
		for (int i = 0; i < rectangles.length; i++) {
			if (((int) rectangles[i].getStartXValue() >= lastPosition) && (!rectangles[i].isInUse())) {
				if (leftMostCNV == null) {
					leftMostCNV = rectangles[i];
				} else {
					int oldDifference = (int) leftMostCNV.getStartXValue() - lastPosition;
					int newDifference = (int) rectangles[i].getStartXValue() - lastPosition;
					if (newDifference < oldDifference) {
						leftMostCNV = rectangles[i];
					}
				}
			}
		}
		return leftMostCNV;
	}

	/**
	 * Provide the list of rectangles to render
	 * 
	 * @param rects
	 */
	void setRectangles(CNVRectangle[] rects) {
		rectangles = rects;
		lowestStart = getLowestRect();
		lines.clear();
		clearUsed();
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
	}

	/**
	 * Set the height of the rectangles (Configured in CompConfig)
	 * 
	 * @param height
	 */
	void setRectangleHeight(int height) {
		rectangleHeight = height;
	}

	/**
	 * Set the display mode (Configured in CompConfig)
	 * 
	 * @param dm
	 */
	public void setDisplayMode(String dm) {
		displayMode = dm;
	}

	/*
	 * Implemented interfaces
	 */
	@Override
	public void mouseClicked(MouseEvent e) {
		for (int i = 0; i < rectangles.length; i++) {
			Rectangle rect = rectangles[i].getRect();
			if (rect.contains(e.getPoint())) {
				Vector<CNVariant> currentCNVs = rectangles[i].getCNVs();
				firePropertyChange("selectedCNV", selectedCNVs, currentCNVs);
				selectedCNVs = currentCNVs;
				rectangles[i].setSelected(true);
				// TODO link off to Trailer
			} else {
				rectangles[i].setSelected(false);
			}
		}
		repaint();
	}

	@Override
	public void mouseMoved(MouseEvent e) {
		setToolTipText(null);
		for (int i = 0; i < rectangles.length; i++) {
			Rectangle rect = rectangles[i].getRect();
			if (rect.contains(e.getPoint())) {
				Vector<CNVariant> currentCNVs = rectangles[i].getCNVs();
				CNVariant cnv = rectangles[i].getCNV();
				String toolTipText = "<html>IID: " + cnv.getIndividualID() + "<br/>FID: " + cnv.getFamilyID() + "<br/>Length: " + cnv.getSize() + "<br/>Copies: " + cnv.getCN() + "<br/>Probes: " + cnv.getNumMarkers() + "<br/>Score: " + cnv.getScore();
				// Add an indication of how many other CNVs are in this collapsed view
				if (currentCNVs.size() > 1) {
					toolTipText += "<br/>Plus " + currentCNVs.size() + " others</html>";
				} else {
					toolTipText += "</html>";
				}
				setToolTipText(toolTipText);
				break;
			}
		}
	}

	@Override
	public void mousePressed(MouseEvent e) {
		// Not doing anything on mouse press
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
		// Not doing anything on mouse drag
	}
}
