package org.genvisis.gwas.parsing;

import java.util.ArrayList;
import java.util.List;
import com.google.common.collect.ImmutableList;

public abstract class AbstractColumnFilter implements ColumnFilter {

  private final List<FileColumn<?>> columns;

  public AbstractColumnFilter() {
    this.columns = new ArrayList<>();
  }

  public AbstractColumnFilter(FileColumn<?>... columns) {
    this.columns = ImmutableList.copyOf(columns);
  }

  public AbstractColumnFilter(Iterable<FileColumn<?>> columns) {
    this.columns = ImmutableList.copyOf(columns);
  }

  @Override
  public List<FileColumn<?>> getFilterColumns() {
    return columns;
  }

}
