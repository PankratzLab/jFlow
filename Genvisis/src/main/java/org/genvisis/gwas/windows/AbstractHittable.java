package org.genvisis.gwas.windows;

import org.genvisis.cnv.util.Java6Helper;

/**
 * Abstract superclass for {@link Hittable} implementations. Covers basic identifying information
 * and the {@link #compareTo(Hittable)} implementation.
 */
public abstract class AbstractHittable implements Hittable {
	private String name;
	private byte chr;
	private int pos;

	public AbstractHittable(String name, byte chr, int pos) {
		this.name = name;
		this.chr = chr;
		this.pos = pos;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public byte getChr() {
		return chr;
	}

	@Override
	public int getPos() {
		return pos;
	}

	@Override
	public int compareTo(Hittable o) {
		int c = Java6Helper.compare(chr, o.getChr());
		if (c == 0) {
			c = Java6Helper.compare(pos, o.getPos());
		}
		if (c == 0) {
			c = name.compareTo(o.getName());
		}
		return c;
	}
}
