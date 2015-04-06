package common;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.UIManager;

public class Grafik {

	public static void scaleCheckBoxIcon(JCheckBox checkbox){
		// from http://stackoverflow.com/a/26995048/875496 
	    boolean previousState = checkbox.isSelected();
	    checkbox.setSelected(false);
	    FontMetrics boxFontMetrics =  checkbox.getFontMetrics(checkbox.getFont());
	    Icon boxIcon = UIManager.getIcon("CheckBox.icon");
	    BufferedImage boxImage = new BufferedImage(
	        boxIcon.getIconWidth(), boxIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB
	    );
	    Graphics graphics = boxImage.createGraphics();
	    try{
	        boxIcon.paintIcon(checkbox, graphics, 0, 0);
	    }finally{
	        graphics.dispose();
	    }
	    ImageIcon newBoxImage = new ImageIcon(boxImage);
	    Image finalBoxImage = newBoxImage.getImage().getScaledInstance(
	        boxFontMetrics.getHeight(), boxFontMetrics.getHeight(), Image.SCALE_SMOOTH
	    );
	    checkbox.setIcon(new ImageIcon(finalBoxImage));
	
	    checkbox.setSelected(true);
	    Icon checkedBoxIcon = UIManager.getIcon("CheckBox.icon");
	    BufferedImage checkedBoxImage = new BufferedImage(
	        checkedBoxIcon.getIconWidth(),  checkedBoxIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB
	    );
	    Graphics checkedGraphics = checkedBoxImage.createGraphics();
	    try{
	        checkedBoxIcon.paintIcon(checkbox, checkedGraphics, 0, 0);
	    }finally{
	        checkedGraphics.dispose();
	    }
	    ImageIcon newCheckedBoxImage = new ImageIcon(checkedBoxImage);
	    Image finalCheckedBoxImage = newCheckedBoxImage.getImage().getScaledInstance(
	        boxFontMetrics.getHeight(), boxFontMetrics.getHeight(), Image.SCALE_SMOOTH
	    );
	    checkbox.setSelectedIcon(new ImageIcon(finalCheckedBoxImage));
	    checkbox.setSelected(false);
	    checkbox.setSelected(previousState);
	}
	
	
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

	public static void drawTriangle(Graphics g, int x, int y, int size, boolean filled) {
		g.drawPolygon(new int[] {x-size/2, x, x+size/2},
					  new int[] {y+size/2, y-size/2, y+size/2},
					  3);

		if (filled) {
			for (int i = 0; i < size; i++) {
				g.drawPolygon(new int[] {x-(size/2-i), x, x+(size/2-i)},
						  new int[] {y+(size/2), y-size/2, y+(size/2)},
						  3);
			}

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
	
	public static int[] getHeatmapColor(double value) {
		int[] color;
		
		if (value < 0 || value > 1) {
			System.err.println("Error - heatmap value needs to be between zero and one, inclusive");
		}
		
		color = new int[] {0,0,0};
		if (value < 0.25) {
			color[0] = 0;
			color[1] = (int) (255 * value / 0.25);
			color[2] = 255;
		} else if (value < 0.50) {
			color[0] = 0;
			color[1] = 255;
			color[2] = (int) (255 - 255 * (value - .25) / .25);
		} else if (value < 0.75) {
			color[0] = (int) (255 * (value - .5 ) / .25);
			color[1] = 255;
			color[2] = 0;
		} else {
			color[0] = 255;
			color[1] = (int) (255 - 255 * (value - .75) / .25);
			color[2] = 0;
		}
		
		return color;
	}
}
