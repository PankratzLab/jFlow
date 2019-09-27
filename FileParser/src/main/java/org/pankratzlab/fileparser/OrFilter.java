package org.pankratzlab.fileparser;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.ImmutableList;

public class OrFilter implements ColumnFilter {

  private final Collection<ColumnFilter> filters;

  /**
   * @param filters
   */
  public OrFilter(Collection<ColumnFilter> filters) {
    super();
    this.filters = Collections.unmodifiableCollection(filters);
  }

  /**
   * @param filters
   */
  public OrFilter(ColumnFilter... filters) {
    super();
    this.filters = Collections.unmodifiableList(Arrays.asList(filters));
  }

  @Override
  public List<FileColumn<?>> getFilterColumns() {
    return filters.stream().map(ColumnFilter::getFilterColumns)
                  .collect(ImmutableList.Builder<FileColumn<?>>::new, ImmutableList.Builder::addAll,
                           (b1, b2) -> b1.addAll(b2.build()))
                  .build();
  }

  @Override
  public boolean filter(DataLine values) {
    return filters.stream().anyMatch(f -> f.filter(values));
  }

}
