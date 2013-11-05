package cnv.var;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;

import cnv.plots.GenericRectangle;

/**
 * CNVRectangle describes a rectangle to be rendered in CompPanel
 * 
 * @author Michael Vieths
 * 
 *         Contains the CNVariant and a color based on the file from which it came
 */
public class CNVRectangle extends GenericRectangle implements Comparable<CNVRectangle> {
    private ArrayList<CNVariant> cnvs;
    private Color                CNVColor;
    private Rectangle            rect;
    private boolean              selected;
    private int                  quantity; // How many CNVs are represented by this rectangle
    private boolean              inUse;

    public CNVRectangle(float startX, float stopX, byte thickness, boolean fill, boolean roundedCorners, byte color, byte layer) {
        // Y coord doesn't matter, that'll get set based on context
        super(startX, 0, stopX, 0, thickness, fill, roundedCorners, color, layer);
        quantity = 1;
        selected = false;
        inUse = false;
        cnvs = new ArrayList<CNVariant>();
    }

    public CNVRectangle(CNVariant variant, int offset) {
        this(((int) variant.getStart() - offset), ((int) variant.getStop() - offset), (byte) 2, true, true, (byte) 2, (byte) 1);
    }

    public ArrayList<CNVariant> getCNVs() {
        return cnvs;
    }

    public CNVariant getCNV() {
        return cnvs.get(0);
    }

    public void setCNVs(ArrayList<CNVariant> variants) {
        cnvs = variants;
    }

    public void addCNV(CNVariant variant) {
        cnvs.add(variant);
    }

    public Color getCNVColor() {
        return CNVColor;
    }

    public int getQuantity() {
        return quantity;
    }

    public void setQuantity(int newQuantity) {
        quantity = newQuantity;
    }

    public void setCNVColor(Color color) {
        CNVColor = color;
    }

    /**
     * Sets the color to render for the CNV. Adjusts the brightness based on number of copies.
     * 
     * @param CNVColor
     */
    public void setCNVColor(Color CNVColor, String displayMode) {
        // Default to 2 copies
        int copies = 2;

        // Only change the brightness if we're in Full or Compressed mode.
        // There's only one CNV in this CNVRectangle for Full and Pack, could be multiple if it's Compressed
        if (cnvs.size() > 0) {
            if (displayMode.equals("Full") || displayMode.equals("Pack")) {
                copies = cnvs.get(0).getCN();
            }
        }

        // Need to adjust the brightness
        float[] hsbVals = Color.RGBtoHSB(CNVColor.getRed(), CNVColor.getGreen(), CNVColor.getBlue(), null);
        float newBrightness = hsbVals[2];

        if (copies > 2) {
            // It's a duplication, make it brighter
            newBrightness *= (copies - 2);
        } else if (copies == 1) {
            // It's a deletion, make it darker
            newBrightness *= 0.5f;
        } else if (copies == 0) {
            // No copies, make it much darker
            newBrightness *= 0.2f;
        } else {
            // Normal number of copies, no change in brightness
        }
        this.CNVColor = Color.getHSBColor(hsbVals[0], hsbVals[1], newBrightness);
    }

    public Rectangle getRect() {
        return rect;
    }

    public void setRect(int x, int y, int width, int height) {
        rect = new Rectangle(x, y, width, height);
    }

    public void setSelected(boolean sel) {
        selected = sel;
    }

    public boolean isSelected() {
        return selected;
    }

    public void setUsed(boolean used) {
        inUse = used;
    }

    public boolean isInUse() {
        return inUse;
    }

    @Override
    // Allow sorting the entire list based first on start position, then on length
    public int compareTo(CNVRectangle o) {
        int retValue = 0;
        float start1 = getStartXValue();
        float start2 = o.getStartXValue();
        float length1 = getStopXValue() - start1;
        float length2 = o.getStopXValue() - start2;

        if (start1 > start2) {
            retValue = 1;
        } else if (start1 < start2) {
            retValue = -1;
        } else {
            // They start at the same spot, but their lengths may not be the same
            if (length1 > length2) {
                retValue = 1;
            } else if (length1 < length2) {
                retValue = -1;
            } else {
                retValue = 0;
            }
        }

        return retValue;
    }
}
