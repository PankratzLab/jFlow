/**
 * 
 */
package org.pankratzlab.utils.gwas.windows;

import java.util.List;

/**
 * @param <T>
 */
public class HitWindow<T extends Hittable> {

  private List<T> indexHits;
  private T firstHit;
  private T lastHit;
  private int window;
  private int numSigHits;
  private int numSuggestiveHits;
  private int numTotalHits;

  /**
   * @param window window size in bp
   * @param indexHits a list, since there can be ties
   * @param firstHit first hit
   * @param lastHit last hit
   * @param numSigHits number of sig hits
   * @param numSuggestiveHits number of suggestive hits
   * @param numTotalHits total hits in window
   */
  public HitWindow(int window, List<T> indexHits, T firstHit, T lastHit, int numSigHits,
                   int numSuggestiveHits, int numTotalHits) {
    super();
    this.window = window;
    this.indexHits = indexHits;
    this.firstHit = firstHit;
    this.lastHit = lastHit;
    this.numSigHits = numSigHits;
    this.numSuggestiveHits = numSuggestiveHits;
    this.numTotalHits = numTotalHits;
  }

  @Override
  public String toString() {
    return "HitWindow [indexHit=" + indexHits + ", firstHit=" + firstHit + ", lastHit=" + lastHit
           + ", window=" + window + ", numSigHits=" + numSigHits + ", numSuggestiveHits="
           + numSuggestiveHits + ", numTotalHits=" + numTotalHits + "]";
  }

  public List<T> getIndexHits() {
    return indexHits;
  }

  public T getFirstHit() {
    return firstHit;
  }

  public T getLastHit() {
    return lastHit;
  }

  public int getWindow() {
    return window;
  }

  public int getNumSigHits() {
    return numSigHits;
  }

  public int getNumSuggestiveHits() {
    return numSuggestiveHits;
  }

  public int getNumTotalHits() {
    return numTotalHits;
  }

}
