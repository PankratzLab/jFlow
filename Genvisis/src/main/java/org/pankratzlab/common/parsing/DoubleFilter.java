package org.pankratzlab.common.parsing;

import org.pankratzlab.shared.stats.Maths.COMPARISON;

/**
 * A filter for comparing double values. Uses {@link COMPARISON}.
 */
public class DoubleFilter extends AbstractColumnFilter {

  FileColumn<Double> valueColumn;
  COMPARISON comparison;
  double comparisonValue;

  public DoubleFilter(FileColumn<Double> valueColumn, COMPARISON comparison,
                      double comparisonValue) {
    super(valueColumn);
    this.valueColumn = valueColumn;
    this.comparison = comparison;
    this.comparisonValue = comparisonValue;
  }

  @Override
  public boolean filter(DataLine values) {
    if (values.hasValid(valueColumn)) {
      Double valD = values.getUnsafe(valueColumn);
      return comparison.check(valD, comparisonValue);
    }
    return false;
  }
}
