package org.genvisis.stats;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.NotPositiveException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.util.MathArrays;

/**
 * This class is extends the functionality of {@link ChiSquareTest} to use some of the more nuanced
 * functionality of {@link ChiSquaredDistribution} that was not exposed by Apache Commons
 */
public class ExpandedChiSquareTest extends ChiSquareTest {

  public ExpandedChiSquareTest() {
    super();
  }

  /**
   * Returns the <i>observed significance level</i>, or
   * <a href= "http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue"> p-value</a>,
   * associated with a <a href="http://www.itl.nist.gov/div898/handbook/prc/section4/prc45.htm">
   * chi-square test of independence</a> based on the input <code>counts</code> array, viewed as a
   * two-way table.
   * <p>
   * The rows of the 2-way table are <code>count[0], ... , count[count.length - 1] </code>
   * </p>
   * <p>
   * <strong>Preconditions</strong>:
   * <ul>
   * <li>All counts must be &ge; 0.</li>
   * <li>The count array must be rectangular (i.e. all count[i] subarrays must have the same
   * length).</li>
   * <li>The 2-way table represented by <code>counts</code> must have at least 2 columns and at
   * least 2 rows.</li></li>
   * </ul>
   * </p>
   * <p>
   * If any of the preconditions are not met, an <code>IllegalArgumentException</code> is thrown.
   * </p>
   *
   * @param counts array representation of 2-way table
   * @return p-value
   * @throws NullArgumentException if the array is null
   * @throws DimensionMismatchException if the array is not rectangular
   * @throws NotPositiveException if {@code counts} has negative entries
   * @throws MaxCountExceededException if an error occurs computing the p-value
   */
  public double chiSquareTestWithAccuracy(final long[][] counts,
                                          double inverseCumAccuracy) throws NullArgumentException,
                                                                     DimensionMismatchException,
                                                                     NotPositiveException,
                                                                     MaxCountExceededException {

    checkArray(counts);
    double df = ((double) counts.length - 1) * ((double) counts[0].length - 1);
    // pass a null rng to avoid unneeded overhead as we will not sample from this distribution
    final ChiSquaredDistribution distribution = new ChiSquaredDistribution(df, inverseCumAccuracy);
    return 1 - distribution.cumulativeProbability(chiSquare(counts));

  }

  /**
   * copied from {@link ChiSquareTest}
   */
  private void checkArray(final long[][] in) throws NullArgumentException,
                                             DimensionMismatchException, NotPositiveException {

    if (in.length < 2) {
      throw new DimensionMismatchException(in.length, 2);
    }

    if (in[0].length < 2) {
      throw new DimensionMismatchException(in[0].length, 2);
    }

    MathArrays.checkRectangular(in);
    MathArrays.checkNonNegative(in);

  }

}
