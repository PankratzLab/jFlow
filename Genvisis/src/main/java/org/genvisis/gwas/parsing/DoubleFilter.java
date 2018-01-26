package org.genvisis.gwas.parsing;

import java.util.Map;

import org.genvisis.stats.Maths.COMPARISON;

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
	public boolean filter(Map<FileColumn<?>, String> values) {
		String valStr = values.get(valueColumn);
		double valD = Double.parseDouble(valStr);
		return comparison.check(valD, comparisonValue);
	}
}
