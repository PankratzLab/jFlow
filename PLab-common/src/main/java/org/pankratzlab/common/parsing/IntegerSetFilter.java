package org.pankratzlab.common.parsing;

import org.pankratzlab.common.stats.Maths.COMPARISON;

/**
 * A filter for comparing a value against a set of comparison values. Uses {@link COMPARISON}.
 */
public class IntegerSetFilter extends AbstractColumnFilter {

  FileColumn<Integer> valueColumn;
  COMPARISON comparison;
  int[] comparisonValues;

  public IntegerSetFilter(FileColumn<Integer> valueColumn, COMPARISON comparison,
                          int... comparisonValues) {
    super(valueColumn);
    this.valueColumn = valueColumn;
    this.comparison = comparison;
    this.comparisonValues = comparisonValues;
  }

  @Override
  public boolean filter(DataLine values) {
    if (values.hasValid(valueColumn)) {
      Integer valI = values.getUnsafe(valueColumn);
      for (int v : comparisonValues) {
        if (!comparison.check(valI, v)) return false;
      }
      return true;
    }
    return false;
  }
}
