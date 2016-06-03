package cnv.plots;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Rectangle2D.Float;

public class GenericRectangle {
	private float startXValue;
	private float startYValue;
	private float stopXValue;
	private float stopYValue;
	private byte thickness;
	private boolean fill;
	private boolean roundedCorners;
	private byte color;
	private byte fillColor;
	private byte layer;
	private boolean editable;
//	private Rectangle2D myRect;

	public GenericRectangle(float startX, float startY, float stopX, float stopY, byte thickness, boolean fill, boolean roundedCorners, byte color, byte layer, boolean editable) {
		this.startXValue = startX;
		this.startYValue = startY;
		this.stopXValue = stopX;
		this.stopYValue = stopY;
		this.thickness = thickness;
		this.fill = fill;
		this.roundedCorners = roundedCorners;
		this.color = color;
		this.fillColor = color;
		this.layer = layer;
		this.editable = editable;
//		this.myRect = new Rectangle2D.Float(Math.min(startX, stopX), Math.min(startY, stopY), Math.max(startX, stopX) - Math.min(startX, stopX), Math.max(startY, stopY) - Math.min(startY, stopY));
	}
	
	@Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + color;
//        result = prime * result + (editable ? 1231 : 1237);
        result = prime * result + (fill ? 1231 : 1237);
        result = prime * result + fillColor;
        result = prime * result + layer;
        result = prime * result + (roundedCorners ? 1231 : 1237);
        result = prime * result + java.lang.Float.floatToIntBits(startXValue);
        result = prime * result + java.lang.Float.floatToIntBits(startYValue);
        result = prime * result + java.lang.Float.floatToIntBits(stopXValue);
        result = prime * result + java.lang.Float.floatToIntBits(stopYValue);
        result = prime * result + thickness;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        GenericRectangle other = (GenericRectangle) obj;
        if (color != other.color)
            return false;
//        if (editable != other.editable)
//            return false;
        if (fill != other.fill)
            return false;
        if (fillColor != other.fillColor)
            return false;
        if (layer != other.layer)
            return false;
        if (roundedCorners != other.roundedCorners)
            return false;
        if (java.lang.Float.floatToIntBits(startXValue) != java.lang.Float.floatToIntBits(other.startXValue))
            return false;
        if (java.lang.Float.floatToIntBits(startYValue) != java.lang.Float.floatToIntBits(other.startYValue))
            return false;
        if (java.lang.Float.floatToIntBits(stopXValue) != java.lang.Float.floatToIntBits(other.stopXValue))
            return false;
        if (java.lang.Float.floatToIntBits(stopYValue) != java.lang.Float.floatToIntBits(other.stopYValue))
            return false;
        if (thickness != other.thickness)
            return false;
        return true;
    }

    public GenericRectangle(float startX, float startY, float stopX, float stopY, byte thickness, boolean fill, boolean roundedCorners, byte color, byte fillColor, byte layer, boolean editable) {
		this.startXValue = startX;
		this.startYValue = startY;
		this.stopXValue = stopX;
		this.stopYValue = stopY;
		this.thickness = thickness;
		this.fill = fill;
		this.roundedCorners = roundedCorners;
		this.color = color;
		this.fillColor = fillColor;
		this.layer = layer;
		this.editable = editable;
//		this.myRect = new Rectangle2D.Float(Math.min(startX, stopX), Math.min(startY, stopY), Math.max(startX, stopX) - Math.min(startX, stopX), Math.max(startY, stopY) - Math.min(startY, stopY));
	}

	public void setColor(byte color) {
		this.color=color;
	}
	
	public void setFillColor(byte color) {
		this.fillColor=color;
	}

	public float getStartXValue() {
		return startXValue;
	}

	public float getStartYValue() {
		return startYValue;
	}

//	public Rectangle2D getRectangle() {
//	    return myRect;
//	}
	
	public float getStopXValue() {
		return stopXValue;
	}

	public float getStopYValue() {
		return stopYValue;
	}

	public byte getThickness() {
		return thickness;
	}

	public boolean getFill() {
		return fill;
	}
	
	public boolean getEditable() {
	    return editable;
	}
	
	public boolean getRoundedCorners() {
		return roundedCorners;
	}
	
	public byte getColor() {
		return color;
	}

	public byte getFillColor() {
		return fillColor;
	}
	
	public byte getLayer() {
		return layer;
	}
	
	public static GenericRectangle[] addToArray(GenericRectangle rectangle, GenericRectangle[] array) {
		GenericRectangle[] newArray;
		
		newArray = new GenericRectangle[array.length+1];
		for (int i = 0; i < array.length; i++) {
			newArray[i] = array[i];
		}
		newArray[array.length] = rectangle;
		
		return newArray;
	}

    public void setEditable(boolean b) {
        this.editable = b;
    }

	// TODO remove if not currently being used
//	public static GenericRectangle[] removeFromArray(GenericRectangle[] array, int index) {
//		GenericRectangle[] newArray;
//		
//		newArray = new GenericRectangle[array.length+1];
//		for (int i = 0; i < array.length; i++) {
//			newArray[i] = array[i];
//		}
//		newArray[array.length] = rectangle;
//		
//		return null;
//	}
}
