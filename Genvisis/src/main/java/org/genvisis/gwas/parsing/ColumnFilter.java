package org.genvisis.gwas.parsing;

import java.util.List;
import java.util.Map;

/**
 * Interface for filtering lines of data. Subclasses should extend {@link AbstractColumnFilter}
 * rather than implementing this interface.
 */
public interface ColumnFilter {
	/**
	 * @return a list of {@link FileColumn} objects that this filter uses for filtering.
	 */
	public List<FileColumn<?>> getFilterColumns();

	/**
	 * 
	 * @param values Line of data.
	 * @return Whether the given line of data passes this filter.
	 */
	public boolean filter(Map<FileColumn<?>, String> values);
}
