package org.pankratzlab.common.parsing;

public class ColumnFilters {

  private ColumnFilters() {}

  public static OrFilter or(ColumnFilter... filters) {
    return new OrFilter(filters);
  }

  public static ColumnFilter not(ColumnFilter filter) {
    return new AbstractColumnFilter(filter.getFilterColumns()) {

      @Override
      public boolean filter(DataLine values) {
        return !filter.filter(values);
      }
    };
  }

}
