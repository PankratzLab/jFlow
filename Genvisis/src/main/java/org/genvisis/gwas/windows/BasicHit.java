/**
 * 
 */
package org.genvisis.gwas.windows;

/**
 * Has basic implementation for {@link Hittable}
 */
public class BasicHit extends AbstractHittable {

  private double pVal;

  /**
   * @param name usually marker name
   * @param chr chromosome
   * @param pos position
   * @param pVal p-value of hit
   */
  public BasicHit(String name, byte chr, int pos, double pVal) {
    super(name, chr, pos);
    this.pVal = pVal;
  }

  @Override
  public double getPval() {
    return pVal;
  }

  @Override
  public String toString() {
    return "BasicHit [name=" + getName() + ", chr=" + getChr() + ", pos=" + getPos() + ", pVal="
           + pVal + "]";
  }

}
