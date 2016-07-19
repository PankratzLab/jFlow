package common;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.ToolTipManager;
import javax.swing.UIManager;

public class Grafik {

	public static void scaleCheckBoxIcon(JCheckBox checkbox) {
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
	    drawThickLine(g, x1, y1, x2, y2, thickness, c, 0);
	}
	
	public static void drawThickLine(Graphics g, int x1, int y1, int x2, int y2, int thickness, Color c, int direction) {
	    drawThickLine(g, x1, y1, x2, y2, thickness, c, direction, true);
	}
 
	public static void drawThickLine(Graphics g, int x1, int y1, int x2, int y2, int thickness, Color c, int direction, boolean scaleThickness) {
	    int dX = x2-x1;
	    int dY = y2-y1;
	    g.setColor(c);

	    int dx, dy;
	    
	    if (scaleThickness) {
    	    double lineLength = Math.sqrt(dX*dX+dY*dY);
    	    double scale = (double)(thickness)/(2*lineLength);
    	    // The x,y increments from an endpoint needed to create a rectangle...
    	    double ddx = -scale*(double)dY;
    	    double ddy = scale*(double)dX;
    	    // TODO what are these two alterations doing?  Preliminary removal doesn't seem to alter visual results...
    	    ddx += (ddx>0)?0.5:-0.5;
    	    ddy += (ddy>0)?0.5:-0.5;
    	    dx = (int)ddx;
    	    dy = (int)ddy;
	    } else {
            double radO = Math.atan2(dX, dY);
            double radT = (Math.PI / 2d) - radO;
            dx = (int) Math.round(-(thickness / 2d) * Math.sin(radT));
            dy = (int) Math.round((thickness / 2d) * Math.cos(radT));
	    }
	    
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
	    
	    if (direction != 0) {
    	    int xH = x1 + (dX / 2);
    	    int yH = y1 + (dY / 2);
    	    int len = thickness * 4;
    	    
    	    double radO = Math.atan2(dX, dY);
    	    double radT = (Math.PI / 2d) - radO;
    	    double radT45 = radT + (Math.PI / 4d); // add 45deg
    	    double xD = -(len) * Math.sin(radT45);
    	    double yD = (len) * Math.cos(radT45);
    	    
    	    drawThickLine(g, xH, yH, (int)(direction == 1 ? xH + xD : xH - xD), (int)(yH + yD), 2, c, 0, scaleThickness);
    	    

            xH = x1 + (dX / 3);
            yH = y1 + (dY / 3);
            len = thickness * 4;
            
            radO = Math.atan2(dX, dY);
            radT = (Math.PI / 2d) - radO;
            radT45 = radT + (Math.PI / 4d); // add 45deg
            xD = -(len) * Math.sin(radT45);
            yD = (len) * Math.cos(radT45);
            
            drawThickLine(g, xH, yH, (int)(direction == 1 ? xH + xD : xH - xD), (int)(yH + yD), 2, c, 0, scaleThickness);

            xH = x1 + 2 * (dX / 3);
            yH = y1 + 2 * (dY / 3);
            len = thickness * 4;
            
            radO = Math.atan2(dX, dY);
            radT = (Math.PI / 2d) - radO;
            radT45 = radT + (Math.PI / 4d); // add 45deg
            xD = -(len) * Math.sin(radT45);
            yD = (len) * Math.cos(radT45);
            
            drawThickLine(g, xH, yH, (int)(direction == 1 ? xH + xD : xH - xD), (int)(yH + yD), 2, c, 0, scaleThickness);
            

    	    
    	    
	    }
	    
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

	public static ImageIcon getImageIcon(String filename) {
		ImageIcon iicon;

//		if (jar) {
			try {
				iicon = new ImageIcon(ClassLoader.getSystemResource(filename));
			} catch (NullPointerException npe) {
				iicon = null;
			}
			if (iicon == null) { // can be null either way given the JVM
				iicon = new ImageIcon(filename);
			}
//		} else {
//			iicon = new ImageIcon(filename);
//		}

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
		
		double[] gates = {.25, .5, .75};
		color = new int[] {0,0,0};
		if (value < gates[0]) {
			color[0] = 0;
			color[1] = (int) (255 * value / 0.25);
			color[2] = 255;
		} else if (value < gates[1]) {
			color[0] = 0;
			color[1] = 255;
			color[2] = (int) (255 - 255 * (value - gates[0]) / .25);
		} else if (value < gates[2]) {
			color[0] = (int) (255 * (value - gates[1] ) / .25);
			color[1] = 255;
			color[2] = 0;
		} else {
			color[0] = 255;
			color[1] = (int) (255 - 255 * (value - gates[2]) / .25);
			color[2] = 0;
		}
		
		return color;
	}


    public static JLabel getToolTipIconLabel(String tooltip) {
        JLabel tooltipLbl = new JLabel("");
        tooltipLbl.setIcon(Grafik.getImageIcon("images/question-mark.png"));
        tooltipLbl.setToolTipText(tooltip);
        tooltipLbl.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                super.mouseClicked(e);
                ToolTipManager.sharedInstance().mouseMoved(e);
            }
        });
        return tooltipLbl;
    }
}
