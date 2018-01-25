package org.genvisis.gwas.results.files;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public abstract class AbstractColumnFilter implements ColumnFilter {
	List<FileColumn<?>> columns;

	public AbstractColumnFilter() {
		this.columns = new ArrayList<>();
	}

	public AbstractColumnFilter(FileColumn<?>... columns) {
		this();
		for (FileColumn<?> fc : columns) {
			this.columns.add(fc);
		}
	}

	@Override
	public List<FileColumn<?>> getFilterColumns() {
		return Collections.unmodifiableList(columns);
	}

}
