package org.genvisis.one.ben.fcs.gating;

import org.genvisis.one.ben.fcs.gating.Gate.RectangleGate;

public class GateDimension {
	String paramName;
	Gate owner;

	@Override
	public String toString() {
		return super.toString() + "--" + paramName;
	}

	public GateDimension(Gate gate, String param) {
		this.paramName = param;
		this.owner = gate;
	}

	public static class RectangleGateDimension extends GateDimension {
		public RectangleGateDimension(RectangleGate gate, String param) {
			super(gate, param);
		}

		public RectangleGateDimension(RectangleGate gate, String param, float min, float max) {
			super(gate, param);
			this.min = Math.min(min, max);
			this.max = Math.max(min, max);
		}

		private float min, max;

		public float getMin() {
			return min;
		}

		public float getMax() {
			return max;
		}

		public void setMin(float min2) {
			min = min2;
			owner.setChanged();
		}

		public void setMax(float max2) {
			max = max2;
			owner.setChanged();
		}
	}

	public String getParam() {
		return paramName;
	}

	public void setParam(String d) {
		this.paramName = d;
	}


}
