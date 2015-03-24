package cnv.plots;

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

	public GenericRectangle(float startX, float startY, float stopX, float stopY, byte thickness, boolean fill, boolean roundedCorners, byte color, byte layer) {
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
	}
	
	public GenericRectangle(float startX, float startY, float stopX, float stopY, byte thickness, boolean fill, boolean roundedCorners, byte color, byte fillColor, byte layer) {
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
