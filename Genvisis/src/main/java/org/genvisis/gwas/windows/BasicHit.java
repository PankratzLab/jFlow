/**
 * 
 */
package org.genvisis.gwas.windows;

import org.genvisis.cnv.util.Java6Helper;

/**
 * Has basic implementation for {@link Hittable}
 *
 */
public class BasicHit implements Hittable {

	private String name;
	private byte chr;
	private int pos;
	private double pVal;


	/**
	 * @param name usually marker name
	 * @param chr chromosome
	 * @param pos position
	 * @param pVal p-value of hit
	 */
	public BasicHit(String name, byte chr, int pos, double pVal) {
		super();
		this.name = name;
		this.chr = chr;
		this.pos = pos;
		this.pVal = pVal;
	}



	@Override
	public int compareTo(Hittable other) {
		int c = getChr() - other.getChr();
		if (c == 0) {
			c = Java6Helper.compare(getPos(), other.getPos());
		}
		return c;
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
	public double getPval() {
		return pVal;
	}



	@Override
	public String toString() {
		return "BasicHit [name=" + name + ", chr=" + chr + ", pos=" + pos + ", pVal=" + pVal + "]";
	}



}
