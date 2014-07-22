package cnv.plots;

import java.io.Serializable;

public class PlotPoint implements Serializable {
	private static final long serialVersionUID = 1L;

	public static final byte FILLED_CIRCLE = 1;
	public static final byte OPEN_CIRCLE = 2;
	public static final byte MISSING = 3;
	public static final byte NOT_A_NUMBER = 4;
	public static final byte FILLED_SQUARE = 5;
	public static final byte FILLED_TRIANGLE = 6;
	public static final byte OPEN_SQUARE = 7;
	public static final String MISSING_STR = "X";
	public static final String NAN_STR = "NaN";

	private String id;
	private byte type;
	private float rawX;
	private float rawY;
	private byte size;
	private byte color;
	private byte layer;
	private boolean highlight;
	private boolean visible;

	public PlotPoint(String id, byte type, float rawX, float rawY, byte size, byte color, byte layer) {
		this.id = id;
		this.type = type;
		this.rawX = rawX;
		this.rawY = rawY;
		this.size = size;
		this.color = color;
		this.layer = layer;
		this.highlight = false;
		this.visible = true;
	}

	/**
	 * Copy constructor. Don't depend on clone() to copy objects. Its messy!!!
	 * 
	 * @param point
	 *            a {@link PlotPoint} which is the source of copy
	 */
	public PlotPoint(PlotPoint point) {
		this.id = point.id;
		this.type = point.type;
		this.rawX = point.rawX;
		this.rawY = point.rawY;
		this.size = point.size;
		this.color = point.color;
		this.layer = point.layer;
		this.highlight = point.highlight;
		this.visible = point.visible;
	}

	public String getId() {
		return id;
	}

	public byte getType() {
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

	public byte getLayer() {
		return layer;
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

	public boolean isVisble() {
		return visible;
	}
}
