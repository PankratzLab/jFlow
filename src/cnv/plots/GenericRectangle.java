package cnv.plots;

public class GenericRectangle {
	private float startX;
	private float startY;
	private float stopX;
	private float stopY;
	private byte thickness;
	private boolean fill;
	private boolean roundedCorners;
	private byte color;
	private byte layer;

	public GenericRectangle(float startX, float startY, float stopX, float stopY, byte thickness, boolean fill, boolean roundedCorners, byte color, byte layer) {
		this.startX = startX;
		this.startY = startY;
		this.stopX = stopX;
		this.stopY = stopY;
		this.thickness = thickness;
		this.fill = fill;
		this.roundedCorners = roundedCorners;
		this.color = color;
		this.layer = layer;
	}

	public void setColor(byte color) {
		this.color=color;
	}

	public float getStartX() {
		return startX;
	}

	public float getStartY() {
		return startY;
	}

	public float getStopX() {
		return stopX;
	}

	public float getStopY() {
		return stopY;
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

	public static GenericRectangle[] removeFromArray(GenericRectangle[] array, int index) {
		GenericRectangle[] newArray;
		
//		newArray = new GenericRectangle[array.length+1];
//		for (int i = 0; i < array.length; i++) {
//			newArray[i] = array[i];
//		}
//		newArray[array.length] = rectangle;
		
		return null;
	}
}
