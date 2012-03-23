package cnv.plots;

import java.io.Serializable;

public class PlotPoint implements Serializable {
	private static final long serialVersionUID = 1L;
	
	public static final byte FILLED_CIRCLE = 1;
	public static final byte OPEN_CIRCLE = 2;
	public static final byte MISSING = 3;
	public static final byte NOT_A_NUMBER = 4;
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

	public PlotPoint(String id, byte type, float rawX, float rawY, byte size, byte color, byte layer) {
		this.id = id;
		this.type = type;
		this.rawX = rawX;
		this.rawY = rawY;
		this.size = size;
		this.color = color;
		this.layer = layer;
		this.highlight = false;
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

	public void setHighlighted(boolean status) {
		highlight = status;
	}

	public boolean isHighlighted() {
		return highlight;
	}

}
