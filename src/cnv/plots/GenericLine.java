package cnv.plots;

public class GenericLine {
	private float startX;
	private float startY;
	private float stopX;
	private float stopY;
	private byte thickness;
	private byte color;
	private byte layer;

	public GenericLine(float startX, float startY, float stopX, float stopY, byte thickness, byte color, byte layer) {
		this.startX = startX;
		this.startY = startY;
		this.stopX = stopX;
		this.stopY = stopY;
		this.thickness = thickness;
		this.color = color;
		this.layer = layer;
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
	public GenericLine(PlotPoint startPoint, PlotPoint endPoint, byte thickness, byte color, byte layer, boolean swapAxes) {
		this.startX = swapAxes ? startPoint.getRawY() : startPoint.getRawX();
		this.startY = swapAxes ? startPoint.getRawX() : startPoint.getRawY();
		this.stopX = swapAxes ? startPoint.getRawY() : endPoint.getRawX();
		this.stopY = swapAxes ? startPoint.getRawX() : endPoint.getRawY();
		this.thickness = thickness;
		this.color = color;
		this.layer = layer;
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
}
