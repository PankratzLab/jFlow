package org.genvisis.gwas.parsing;

import java.util.Map;

public class IndexedFileColumn extends AbstractFileColumn<String> {

	final int index;

	public IndexedFileColumn(String name, int inputIndex) {
		super(name);
		this.index = inputIndex;
	}

	@Override
	public void initialize(Map<String, Integer> headerMap) {
		if (!headerMap.values().contains(this.index)) {
			throw new IllegalStateException("Index " + this.index + " is not valid for header.");
		}
	}

	@Override
	public String getValue(String[] line) {
		return line[index];
	}

}
