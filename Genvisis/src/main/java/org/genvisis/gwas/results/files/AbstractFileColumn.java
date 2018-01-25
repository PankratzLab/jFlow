package org.genvisis.gwas.results.files;

public abstract class AbstractFileColumn<T> implements FileColumn<T> {
	private final String name;

	public AbstractFileColumn(String nm) {
		this.name = nm;
	}

	@Override
	public String getName() {
		return this.name;
	}

	// @Override
	// public abstract int hashCode();
	//
	// @Override
	// public abstract boolean equals(Object o);

}
