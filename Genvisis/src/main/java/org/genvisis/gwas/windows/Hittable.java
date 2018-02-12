/**
 * 
 */
package org.genvisis.gwas.windows;

/**
 * Classes implementing this interface can be used with the {@link GeneralHitWindowDetector}
 */
public interface Hittable extends Comparable<Hittable> {

  /**
   * @return the name of this hit (usually a marker name)
   */
  public String getName();

  /**
   * @return the chromosome of this hit
   */
  public byte getChr();

  /**
   * @return the bp of this hit
   */
  public int getPos();

  /**
   * @return the p- value for this hit
   */
  public double getPval();

}
