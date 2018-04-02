/**
 * 
 */
package org.genvisis.gwas.windows;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import org.genvisis.gwas.HitWindows;

/**
 * Looking to implement {@link HitWindows} functionality, but operate on more general pieces.
 *
 * @param <T> it's a hit
 */
public class GeneralHitWindowDetector<T extends Hittable> implements Iterator<HitWindow<T>> {

  private final List<T> hittables;
  private final int window; // in base pairs
  private final double windowExtensionThreshold; // pvalue to extend a region
  private final double indexThreshold; // pvalue to flag a region
  private int currentIndex;// where we are in the hittables
  private HitWindow<T> currentWindow;

  /**
   * @param hittablesInput search for {@link HitWindow}s within this list
   * @param window the bp window to search
   * @param windowExtensionThreshold pvalue threshold to extend window
   * @param indexThreshold pvalue to flag a region
   */
  public GeneralHitWindowDetector(List<T> hittablesInput, int window,
                                  double windowExtensionThreshold, double indexThreshold) {
    super();
    this.hittables = new ArrayList<>(hittablesInput);
    Collections.sort(hittables);// sorts by chr and pos
    this.window = window;
    this.windowExtensionThreshold = windowExtensionThreshold;
    this.indexThreshold = indexThreshold;
    this.currentIndex = 0;
  }

  /*
   * (non-Javadoc)
   * @see java.util.Iterator#hasNext()
   */
  @Override
  public boolean hasNext() {
    if (currentIndex >= hittables.size()) {
      return false;
    } else {
      return findNextWindow();
    }
  }

  /*
   * (non-Javadoc)
   * @see java.util.Iterator#next()
   */
  @Override
  public HitWindow<T> next() {
    if (currentWindow == null) {
      throw new NoSuchElementException("Internal error, null window");
    }
    return currentWindow;
  }

  @Override
  public void remove() {
    // Eclipse complains about this being missing. Maven doesn't seem to care. Java is dumb.
    throw new UnsupportedOperationException();
  }

  /**
   * This sets
   * 
   * @return true if a new window was found
   */
  private boolean findNextWindow() {
    currentWindow = null;
    for (int i = currentIndex; i < hittables.size(); i++) {
      if (hittables.get(i).getPval() < indexThreshold) {
        int startIndex = i;
        ArrayList<Integer> minIndices = new ArrayList<>();// in case of ties
        minIndices.add(i);
        double minPval = hittables.get(i).getPval();
        int offset = 0;
        int numSig = 1;
        int numSuggestive = 1;
        while (startIndex - offset - 1 >= 0
               && hittables.get(startIndex).getChr() == hittables.get(startIndex - offset
                                                                      - 1)
                                                                 .getChr()
               && hittables.get(startIndex).getPos()
                  - window * 2 <= hittables.get(startIndex - offset - 1).getPos()) { // *2
                                                                                                                                                                                                                                                       // required
                                                                                                                                                                                                                                                       // to ensure
                                                                                                                                                                                                                                                       // that
                                                                                                                                                                                                                                                       // there
                                                                                                                                                                                                                                                       // are no
                                                                                                                                                                                                                                                       // overlapping
                                                                                                                                                                                                                                                       // SNPs
                                                                                                                                                                                                                                                       // 500kb
                                                                                                                                                                                                                                                       // after
                                                                                                                                                                                                                                                       // last
                                                                                                                                                                                                                                                       // hit and
                                                                                                                                                                                                                                                       // 500kb
                                                                                                                                                                                                                                                       // before
                                                                                                                                                                                                                                                       // next
                                                                                                                                                                                                                                                       // hit is
                                                                                                                                                                                                                                                       // technically
                                                                                                                                                                                                                                                       // a 1M
                                                                                                                                                                                                                                                       // region
                                                                                                                                                                                                                                                       // that
                                                                                                                                                                                                                                                       // should
                                                                                                                                                                                                                                                       // be merged
          offset++;
          if (hittables.get(startIndex - offset).getPval() < windowExtensionThreshold) {
            startIndex -= offset;
            offset = 0;
            numSuggestive++;
          }
        }

        int stopIndex = i;
        offset = 0;

        while (stopIndex + offset + 1 < hittables.size()
               && hittables.get(stopIndex).getChr() == hittables.get(stopIndex + offset + 1)
                                                                .getChr()
               && hittables.get(stopIndex).getPos()
                  + window >= hittables.get(stopIndex + offset + 1).getPos()) { // don't want the
                                                                                                                                                                                                                                                            // 2* here,
                                                                                                                                                                                                                                                            // though
          offset++;
          if (hittables.get(stopIndex + offset).getPval() < indexThreshold) {
            numSig++;
          }
          if (hittables.get(stopIndex + offset).getPval() < windowExtensionThreshold) {
            stopIndex += offset;
            offset = 0;
            numSuggestive++;
          }

          if (hittables.get(stopIndex).getPval() < minPval) {
            minIndices = new ArrayList<>();
            minIndices.add(stopIndex);
            minPval = hittables.get(stopIndex).getPval();
          }
          if (!minIndices.contains(stopIndex) && hittables.get(stopIndex).getPval() == minPval) { // in
                                                                                                 // case
                                                                                                 // of
                                                                                                 // exact
                                                                                                 // ties,
                                                                                                 // like
                                                                                                 // in
                                                                                                 // the
                                                                                                 // coding test
            minIndices.add(stopIndex);

          }

        }
        ArrayList<T> indexHits = new ArrayList<>();
        for (Integer minIndex : minIndices) {
          indexHits.add(hittables.get(minIndex));
        }
        currentWindow = new HitWindow<>(window, indexHits, hittables.get(startIndex),
                                        hittables.get(stopIndex), numSig, numSuggestive,
                                        stopIndex - startIndex + 1);

        currentIndex = stopIndex + offset + 1;// +1 since this is not i
        return true;
      }
    }
    return false;
  }

}
