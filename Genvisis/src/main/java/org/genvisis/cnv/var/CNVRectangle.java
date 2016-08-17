package org.genvisis.cnv.var;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;

import org.genvisis.cnv.plots.GenericRectangle;
import org.genvisis.filesys.CNVariant;

/**
 * CNVRectangle describes a rectangle to be rendered in CompPanel
 *
 * @author Michael Vieths
 *
 *         Contains the CNVariant and a color based on the file from which it came
 */
public class CNVRectangle extends GenericRectangle implements Comparable<CNVRectangle> {
  private ArrayList<CNVariant> cnvs;
  private Color CNVColor;
  private Rectangle rect;
  private boolean selected;
  private int quantity; // How many CNVs are represented by this rectangle
  private boolean inUse;
  private String filename;

  public CNVRectangle(CNVariant variant, int offset) {
    this((variant.getStart() - offset), (variant.getStop() - offset), (byte) 2, true, true,
        (byte) 2, (byte) 1);
  }

  public CNVRectangle(float startX, float stopX, byte thickness, boolean fill,
      boolean roundedCorners, byte color, byte layer) {
    // Y coord doesn't matter, that'll get set based on context
    super(startX, 0, stopX, 0, thickness, fill, roundedCorners, color, layer, false);
    quantity = 1;
    selected = false;
    inUse = false;
    cnvs = new ArrayList<CNVariant>();
  }

  /**
   * Associates a CNV with this rectangle
   * 
   * @param variant The CNV to be added
   */
  public void addCNV(CNVariant variant) {
    cnvs.add(variant);
  }

  /**
   * Allow sorting the entire list based first on start position, then on length
   */
  @Override
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

  /**
   * 
   * @return The first CNV in the list (the only CNV in non-Collapsed display modes)
   */
  public CNVariant getCNV() {
    return cnvs.get(0);
  }

  /**
   * @return The Color of the rectangle
   */
  public Color getCNVColor() {
    return CNVColor;
  }

  /**
   * 
   * @return All CNVs associated with this rectangle
   */
  public ArrayList<CNVariant> getCNVs() {
    return cnvs;
  }

  /**
   * @return The associated filename for this CNV
   */
  public String getFilename() {
    return filename;
  }

  /**
   * @return The number of CNVs associated with this rectangle
   */
  public int getQuantity() {
    return quantity;
  }

  /**
   * @return The Rectangle object associated with this CNV
   */
  public Rectangle getRect() {
    return rect;
  }

  /**
   * Indicates the rectangle has already been sorted
   * 
   * @return true if it's been sorted, false otherwise
   */
  public boolean isInUse() {
    return inUse;
  }

  /**
   * 
   * @return true if the rectangle is selected, false otherwise
   */
  public boolean isSelected() {
    return selected;
  }

  /**
   * @param color The Color in which the rectangle should be rendered
   */
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
    // There's only one CNV in this CNVRectangle for Full and Pack, could be multiple if it's
    // Compressed
    if (cnvs.size() > 0) {
      if (displayMode.equals("Full") || displayMode.equals("Pack")) {
        copies = cnvs.get(0).getCN();
      }
    }

    // Need to adjust the brightness
    float[] hsbVals =
        Color.RGBtoHSB(CNVColor.getRed(), CNVColor.getGreen(), CNVColor.getBlue(), null);
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

  /**
   * 
   * @param variants A list of CNVs to be associated with this rectangle
   */
  public void setCNVs(ArrayList<CNVariant> variants) {
    cnvs = variants;
  }

  /**
   * 
   * @param filename The name of the file from which this CNV originated
   */
  public void setFilename(String filename) {
    this.filename = filename;
  }

  /**
   * @param newQuantity The number of CNVs associated with this rectangle
   */
  public void setQuantity(int newQuantity) {
    quantity = newQuantity;
  }

  /**
   * Defines a rectangle for this CNV. X and Y are not direct coordinates, they are base positions.
   * They get translated to relative coordinates in CompPanel.
   * 
   * @param x Start X coordinate for this rectangle
   * @param y Start Y coordinate for this rectangle
   * @param width Width of this rectangle
   * @param height Height of this rectangle
   */
  public void setRect(int x, int y, int width, int height) {
    rect = new Rectangle(x, y, width, height);
  }

  /**
   * Indicates whether this rectangle is selected. Selected rectangles have a black border drawn
   * around them
   * 
   * @param sel True to select it, false otherwise
   */
  public void setSelected(boolean sel) {
    selected = sel;
  }

  /**
   * Flag for CNVRectangles to use when deciding whether the rectangle has been sorted
   * 
   * @param used
   */
  public void setUsed(boolean used) {
    inUse = used;
  }
}
