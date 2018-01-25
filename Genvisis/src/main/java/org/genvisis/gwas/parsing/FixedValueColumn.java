package org.genvisis.gwas.results.files;

import java.util.Map;

public class FixedValueColumn extends AbstractFileColumn<String> {

	private String value;

	public FixedValueColumn(String name, String value) {
		super(name);
		this.value = value;
	}

	@Override
	public String getValue(String[] line) throws ParseFailureException {
		return value;
	}

	@Override
	public void initialize(Map<String, Integer> headerMap) {
		// no-op
	}

}
