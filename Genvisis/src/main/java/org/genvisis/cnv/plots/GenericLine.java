package org.genvisis.cnv.plots;

public class GenericLine {
	private final float startX;
	private final float startY;
	private final float stopX;
	private final float stopY;
	private final byte thickness;
	private final byte color;
	private final byte layer;
	/**
	 * 0 - no direction, 1 - towards first point (start point), 2 - towards second point (end point)
	 */
	private int direction;
	private final boolean scale;

	public GenericLine(	float startX, float startY, float stopX, float stopY, byte thickness,
											byte color, byte layer, int direction) {
		this.startX = startX;
		this.startY = startY;
		this.stopX = stopX;
		this.stopY = stopY;
		this.thickness = thickness;
		this.color = color;
		this.layer = layer;
		scale = true;
	}

	public GenericLine(	float startX, float startY, float stopX, float stopY, byte thickness,
											byte color, byte layer, int direction, boolean shouldScale) {
		this.startX = startX;
		this.startY = startY;
		this.stopX = stopX;
		this.stopY = stopY;
		this.thickness = thickness;
		this.color = color;
		this.layer = layer;
		scale = shouldScale;
	}

	/*
	 * Stub constructor for backwards compatibility
	 */
	public GenericLine(	float startX, float startY, float stopX, float stopY, byte thickness,
											byte color, byte layer) {
		this(startX, startY, stopX, stopY, thickness, color, layer, 0);
	}

	/*
	 * Stub constructor for backwards compatibility
	 */
	public GenericLine(	PlotPoint startPoint, PlotPoint endPoint, byte thickness, byte color,
											byte layer, boolean swapAxes) {
		this(startPoint, endPoint, thickness, color, layer, swapAxes, 0, true);
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
	public GenericLine(	PlotPoint startPoint, PlotPoint endPoint, byte thickness, byte color,
											byte layer, boolean swapAxes, int direction, boolean shouldScale) {
		startX = swapAxes ? startPoint.getRawY() : startPoint.getRawX();
		startY = swapAxes ? startPoint.getRawX() : startPoint.getRawY();
		stopX = swapAxes ? startPoint.getRawY() : endPoint.getRawX();
		stopY = swapAxes ? startPoint.getRawX() : endPoint.getRawY();
		this.thickness = thickness;
		this.color = color;
		this.layer = layer;
		this.direction = direction;
		scale = shouldScale;
	}

	public boolean getShouldScale() {
		return scale;
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
