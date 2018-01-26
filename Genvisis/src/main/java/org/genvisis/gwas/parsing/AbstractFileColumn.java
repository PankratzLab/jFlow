package org.genvisis.gwas.parsing;

public abstract class AbstractFileColumn<T> implements FileColumn<T> {
	private final String name;

	public AbstractFileColumn(String nm) {
		this.name = nm;
	}

	@Override
	public String getName() {
		return this.name;
	}

	@Override
	// overridden to force subclasses to define their own hashCode method
	public abstract int hashCode();

	@Override
	// overridden to force subclasses to define their own equals method
	public abstract boolean equals(Object o);

}
