package cnv.plots;

public class GenericLine {
	private float startX;
	private float startY;
	private float stopX;
	private float stopY;
	private byte thickness;
	private byte color;
	private byte layer;
	/** 0 - no direction, 1 - towards first point (start point), 2 - towards second point (end point) */
	private int direction;

	public GenericLine(float startX, float startY, float stopX, float stopY, byte thickness, byte color, byte layer, int direction) {
		this.startX = startX;
		this.startY = startY;
		this.stopX = stopX;
		this.stopY = stopY;
		this.thickness = thickness;
		this.color = color;
		this.layer = layer;
	}
	
	/*
	 * Stub constructor for backwards compatibility
	 */
	public GenericLine(float startX, float startY, float stopX, float stopY, byte thickness, byte color, byte layer) {
	    this(startX, startY, stopX, stopY, thickness, color, layer, 0);
	}
	
	/*
	 * Stub constructor for backwards compatibility
	 */
	public GenericLine(PlotPoint startPoint, PlotPoint endPoint, byte thickness, byte color, byte layer, boolean swapAxes) {
	    this(startPoint, endPoint, thickness, color, layer, swapAxes, 0);
	}
	
	/**
	 * Constructor for GenericLine from two {@link PlotPoint}
	 * 
	 * @param startPoint
	 * @param endPoint
	 * @param thickness
	 * @param color
	 * @param layer
	 * @param swapAxes
	 */
	public GenericLine(PlotPoint startPoint, PlotPoint endPoint, byte thickness, byte color, byte layer, boolean swapAxes, int direction) {
		this.startX = swapAxes ? startPoint.getRawY() : startPoint.getRawX();
		this.startY = swapAxes ? startPoint.getRawX() : startPoint.getRawY();
		this.stopX = swapAxes ? startPoint.getRawY() : endPoint.getRawX();
		this.stopY = swapAxes ? startPoint.getRawX() : endPoint.getRawY();
		this.thickness = thickness;
		this.color = color;
		this.layer = layer;
		this.direction = direction;
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

	public byte getColor() {
		return color;
	}

	public byte getLayer() {
		return layer;
	}
	
	public int getDirection() {
        return direction;
    }
	
}
