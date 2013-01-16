/**
 * 
 */
package cnv.plots;

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
public class CompPanel extends JPanel implements MouseListener, MouseMotionListener {// extends AbstractPanel {
	public static final long serialVersionUID = 1L;
	CNVRectangle[] rectangles;
	float scalingFactor;

	public CompPanel(CompPlot cp) {
		// super();
		// displayXaxis = false;
		// displayYaxis = false;
		// colorScheme = new Color[3];
		// colorScheme[0] = Color.blue;
		// colorScheme[1] = Color.red;
		// colorScheme[2] = Color.green;
		addMouseListener(this);
		addMouseMotionListener(this);
	}

	void setRectangles(CNVRectangle[] rects) {
		rectangles = rects;
		repaint();
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
			int width = Math.round(((int) rectangles[i].getStopX() - (int) rectangles[i].getStartX()) * scalingFactor);
			int height = 10;
			int x = Math.round((int) rectangles[i].getStartX() * scalingFactor);
			int y = (i * height);

			// Store the rectangle for later bounds checking
			rectangles[i].setRect(x, y, width, height);

			g.setColor(rectangles[i].getCNVColor());
			g.fillRect(x, y, width, height);
		}
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		for (int i = 0; i < rectangles.length; i++) {
			Rectangle rect = rectangles[i].getRect();
			if (rect.contains(e.getPoint())) {
				// TODO link off to Trailer
				break;
			}
		}
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
