package org.genvisis.gwas.parsing;

import java.util.Map;

public class IntegerWrapperColumn extends CachedFileColumn<Integer> {
	private FileColumn<?> base;

	public IntegerWrapperColumn(FileColumn<?> base) {
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
	public Integer calculateValue(String[] line) throws ParseFailureException {
		try {
			return Integer.parseInt(getBaseValue(line));
		} catch (NumberFormatException e) {
			throw new ParseFailureException(e);
		}
	}
}
