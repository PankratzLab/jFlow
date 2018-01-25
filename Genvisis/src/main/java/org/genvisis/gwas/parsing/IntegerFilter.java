package org.genvisis.gwas.results.files;

import java.util.Map;

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
	public boolean filter(Map<FileColumn<?>, String> values) {
		String valStr = values.get(valueColumn);
		int valI = Integer.parseInt(valStr);
		return comparison.check(valI, comparisonValue);
	}
}
