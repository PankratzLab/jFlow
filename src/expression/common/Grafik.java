package common;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import javax.swing.ImageIcon;

public class Grafik {
	public static void drawThickLine(Graphics g, int x1, int y1, int x2, int y2, int thickness, Color c) {
		int dX = x2-x1;
		int dY = y2-y1;
		double lineLength = Math.sqrt(dX*dX+dY*dY);
		g.setColor(c);

		double scale = (double)(thickness)/(2*lineLength);

		// The x,y increments from an endpoint needed to create a rectangle...
		double ddx = -scale*(double)dY;
		double ddy = scale*(double)dX;
		ddx += (ddx>0)?0.5:-0.5;
		ddy += (ddy>0)?0.5:-0.5;
		int dx = (int)ddx;
		int dy = (int)ddy;

		// Now we can compute the corner points...
		int xPoints[] = new int[4];
		int yPoints[] = new int[4];

		xPoints[0] = x1+dx;
		yPoints[0] = y1+dy;
		xPoints[1] = x1-dx;
		yPoints[1] = y1-dy;
		xPoints[2] = x2-dx;
		yPoints[2] = y2-dy;
		xPoints[3] = x2+dx;
		yPoints[3] = y2+dy;

		g.fillPolygon(xPoints, yPoints, 4);
	}
	
	public static void drawCircle(Graphics g, int x, int y, int size, boolean filled, Color c) {
		g.setColor(c);
		drawCircle(g, x, y, size, filled);
	}

	public static void drawCircle(Graphics g, int x, int y, int size, boolean filled) {
		if (filled) {
			g.fillOval(x-size/2, y-size/2, size, size);
		} else {
			g.drawOval(x-size/2, y-size/2, size, size);
		}
	}

	public static ImageIcon getImageIcon(String filename, boolean jar) {
		ImageIcon iicon;

		if (jar) {
			try {
				iicon = new ImageIcon(ClassLoader.getSystemResource(filename));
			} catch (NullPointerException npe) {
				iicon = null;
			}
			if (iicon == null) { // can be null either way given the JVM
				iicon = new ImageIcon(filename);
			}
		} else {
			iicon = new ImageIcon(filename);
		}

		return iicon;		
	}

	public static BufferedImage rotateImage(BufferedImage image, boolean leftNotRight) {
		BufferedImage rotated;
		AffineTransform at;
		Graphics2D g2;
		int h, w;

		h = image.getHeight();
		w = image.getWidth();
		rotated = new BufferedImage(h, w, image.getType());
		g2 = rotated.createGraphics();

		at = AffineTransform.getTranslateInstance((h-w)/2.0, (w-h)/2.0);
		at.rotate(Math.toRadians(leftNotRight?270:90), w/2.0, h/2.0);
		g2.drawRenderedImage(image, at);
		g2.dispose();

		return rotated;
	}
	
	public static int getTextWidth(String text, Graphics g) {
		return g.getFontMetrics(g.getFont()).stringWidth(text);
	}
}
