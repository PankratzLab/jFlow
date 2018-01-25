package org.genvisis.gwas.parsing;

import java.util.Map;

public class DoubleWrapperColumn extends CachedFileColumn<Double> {
	private FileColumn<?> base;

	public DoubleWrapperColumn(FileColumn<?> base) {
		super(base.getName());
		this.base = base;
	}

	@Override
	public void initialize(Map<String, Integer> headerMap) {
		base.initialize(headerMap);
	}

	public String getBaseValue(String[] line) throws ParseFailureException {
		return base.getValue(line).toString();
	}

	@Override
	public Double calculateValue(String[] line) throws ParseFailureException {
		try {
			return Double.parseDouble(getBaseValue(line));
		} catch (NumberFormatException e) {
			throw new ParseFailureException(e);
		}
	}
}
