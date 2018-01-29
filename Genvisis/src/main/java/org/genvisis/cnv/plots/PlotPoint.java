package org.genvisis.cnv.plots;

import java.awt.Color;
import java.io.Serializable;

public class PlotPoint implements Serializable {
	private static final long serialVersionUID = 1L;

	public enum PointType {
		FILLED_CIRCLE,
		OPEN_CIRCLE,
		MISSING,
		NOT_A_NUMBER,
		FILLED_SQUARE,
		FILLED_TRIANGLE,
		OPEN_SQUARE;
	}

	public static final String MISSING_STR = "X";
	public static final String NAN_STR = "NaN";

	private final String id;
	private final float rawX;
	private final float rawY;
	private PointType type;
	private byte size;
	private Color tempColor;
	private byte color;
	private byte layer;
	private boolean highlight;
	private boolean visible;

	public PlotPoint(String id, PointType type, float rawX, float rawY, byte size, byte color,
									 byte layer) {
		this.id = id;
		this.type = type;
		this.rawX = rawX;
		this.rawY = rawY;
		this.size = size;
		this.color = color;
		this.layer = layer;
		highlight = false;
		visible = true;
	}

	/**
	 * Copy constructor. Don't depend on clone() to copy objects. Its messy!!!
	 *
	 * @param point a {@link PlotPoint} which is the source of copy
	 */
	public PlotPoint(PlotPoint point) {
		id = point.id;
		type = point.type;
		rawX = point.rawX;
		rawY = point.rawY;
		size = point.size;
		color = point.color;
		layer = point.layer;
		highlight = point.highlight;
		visible = point.visible;
	}

	public String getId() {
		return id;
	}

	public PointType getType() {
		return type;
	}

	public float getRawX() {
		return rawX;
	}

	public float getRawY() {
		return rawY;
	}

	public void setSize(byte size) {
		this.size = size;
	}

	public byte getSize() {
		return size;
	}

	public byte getColor() {
		return color;
	}

	public Color getTempColor() {
		return tempColor;
	}

	public byte getLayer() {
		return layer;
	}

	public void setType(PointType type) {
		this.type = type;
	}

	public void setLayer(byte layer) {
		this.layer = layer;
	}

	public void setHighlighted(boolean status) {
		highlight = status;
	}

	public void setVisible(boolean status) {
		visible = status;
	}

	public boolean isHighlighted() {
		return highlight;
	}

	public boolean isVisible() {
		return visible;
	}

	public void setTempColor(Color tempColor2) {
		tempColor = tempColor2;
	}

	public void setColor(byte color2) {
		color = color2;
	}
}
