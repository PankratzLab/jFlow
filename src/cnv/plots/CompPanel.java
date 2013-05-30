/**
 * 
 */
package cnv.plots;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

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
	CNVariant selectedCNV;
	int rectangleHeight = 10;

	public CompPanel(CompPlot cp) {
		addMouseListener(this);
		addMouseMotionListener(this);
	}

	void setRectangles(CNVRectangle[] rects) {
		rectangles = rects;
		repaint();
	}

	void setRectangleHeight(int height) {
		rectangleHeight = height;
	}

	/**
	 * Need to know how big the visible window is in bases so we can scale it down to the panel width
	 * 
	 * @param window
	 */
	void setWindow(int window) {
		scalingFactor = (float) getWidth() / window;
	}

	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		// We'll need to scale the relative base to the window size
		for (int i = 0; i < rectangles.length; i++) {
			int width = Math.round(((int) rectangles[i].getStopXValue() - (int) rectangles[i].getStartXValue()) * scalingFactor);
			int height = rectangleHeight;
			int x = Math.round((int) rectangles[i].getStartXValue() * scalingFactor);
			int y = (i * height);

			// Store the rectangle for later bounds checking
			rectangles[i].setRect(x, y, width, height);

			g.setColor(rectangles[i].getCNVColor());

			if (rectangles[i].isSelected()) {
				g.fillRect(x, y, width, height);
				g.setColor(Color.BLACK);
				g.setPaintMode();
				g.drawRect(x, y - 1, width, height);
			} else {
				g.fillRect(x, y, width, height);
			}
		}
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		for (int i = 0; i < rectangles.length; i++) {
			Rectangle rect = rectangles[i].getRect();
			if (rect.contains(e.getPoint())) {
				CNVariant currentCNV = rectangles[i].getCNV();
				firePropertyChange("selectedCNV", selectedCNV, currentCNV);
				selectedCNV = currentCNV;
				rectangles[i].setSelected(true);
				// TODO link off to Trailer
			} else {
				rectangles[i].setSelected(false);
			}
		}
		repaint();
	}

	@Override
	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseDragged(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseMoved(MouseEvent e) {
		setToolTipText(null);
		for (int i = 0; i < rectangles.length; i++) {
			Rectangle rect = rectangles[i].getRect();
			if (rect.contains(e.getPoint())) {
				CNVariant cnv = rectangles[i].getCNV();
				String toolTipText = "<html>IID: " + cnv.getIndividualID() + "<br/>FID: " + cnv.getFamilyID() + "<br/>Length: " + cnv.getSize() + "<br/>Copies: " + cnv.getCN() + "<br/>Probes: " + cnv.getNumMarkers() + "<br/>Score: " + cnv.getScore() + "</html>";
				setToolTipText(toolTipText);
				break;
			}
		}
	}
}
