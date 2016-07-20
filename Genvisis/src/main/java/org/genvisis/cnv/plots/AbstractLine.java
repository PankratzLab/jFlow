package org.genvisis.cnv.plots;

public class AbstractLine {
	private float startX;
	private float startY;
	private float stopX;
	private float stopY;
	private byte thickness;
	private byte color;

	public AbstractLine(float startX, float startY, float stopX, float stopY, byte thickness, byte color) {
		this.startX = startX;
		this.startY = startY;
		this.stopX = stopX;
		this.stopY = stopY;
		this.thickness = thickness;
		this.color = color;
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
}