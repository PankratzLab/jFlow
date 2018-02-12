package org.genvisis.gwas.parsing;

import org.genvisis.stats.Maths.COMPARISON;

/**
 * A filter for comparing integer values. Uses {@link COMPARISON}.
 */
public class IntegerFilter extends AbstractColumnFilter {

  FileColumn<Integer> valueColumn;
  COMPARISON comparison;
  int comparisonValue;

  public IntegerFilter(FileColumn<Integer> valueColumn, COMPARISON comparison,
                       int comparisonValue) {
    super(valueColumn);
    this.valueColumn = valueColumn;
    this.comparison = comparison;
    this.comparisonValue = comparisonValue;
  }

  @Override
  public boolean filter(DataLine values) {
    if (values.hasValid(valueColumn)) {
      Integer valI = values.getUnsafe(valueColumn);
      return comparison.check(valI, comparisonValue);
    }
    return false;
  }
}
