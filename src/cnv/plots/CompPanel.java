/**
 * 
 */
package cnv.plots;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JPanel;

/**
 * @author Michael Vieths
 * 
 */
public class CompPanel extends JPanel implements MouseListener {// extends AbstractPanel {
	public static final long serialVersionUID = 1L;
	CNVRectangle[] rectangles;
	float scalingFactor;

	private Color[] colorScheme = { Color.RED, Color.GREEN, Color.BLUE };

	public CompPanel(CompPlot cp) {
		// super();
		// displayXaxis = false;
		// displayYaxis = false;
		// colorScheme = new Color[3];
		// colorScheme[0] = Color.blue;
		// colorScheme[1] = Color.red;
		// colorScheme[2] = Color.green;
		addMouseListener(this);
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

			g.setColor(rectangles[i].getCNVColor());
			g.fillRect(x, y, width, height);
		}
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub
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
}
