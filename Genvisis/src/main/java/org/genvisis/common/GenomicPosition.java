package org.genvisis.common;

import java.io.Serializable;

public class GenomicPosition implements Serializable, Comparable<GenomicPosition> {

	private static final long serialVersionUID = 1L;

	private final byte chr;
	private final int position;

	/**
	 * @param chr
	 * @param position
	 */
	public GenomicPosition(byte chr, int position) {
		super();
		this.chr = chr;
		this.position = position;
	}

	public byte getChr() {
		return chr;
	}

	public int getPosition() {
		return position;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + chr;
		result = prime * result + position;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		GenomicPosition other = (GenomicPosition) obj;
		if (chr != other.chr)
			return false;
		if (position != other.position)
			return false;
		return true;
	}

	@Override
	public int compareTo(GenomicPosition o) {
		int cmp = Byte.compare(chr, o.chr);
		if (cmp != 0)
			return cmp;
		cmp = Integer.compare(position, o.position);
		return cmp;
	}

	@Override
	public String toString() {
		return chr + ":" + position;
	}


}
