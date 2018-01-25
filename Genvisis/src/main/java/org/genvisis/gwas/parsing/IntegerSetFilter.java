package org.genvisis.gwas.results.files;

import java.util.Map;

import org.genvisis.stats.Maths.COMPARISON;

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
	public boolean filter(Map<FileColumn<?>, String> values) {
		String valStr = values.get(valueColumn);
		int valI = Integer.parseInt(valStr);
		for (int v : comparisonValues) {
			if (!comparison.check(valI, v))
				return false;
		}
		return true;
	}
}
