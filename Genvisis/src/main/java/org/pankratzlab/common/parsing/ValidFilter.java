package org.pankratzlab.common.parsing;

import java.util.List;

public class ValidFilter extends AbstractColumnFilter {

  public ValidFilter(FileColumn<?>... columns) {
    super(columns);
  }

  public ValidFilter(List<FileColumn<?>> columns) {
    super(columns);
  }

  @Override
  public boolean filter(DataLine values) {
    return getFilterColumns().stream().allMatch(values::hasValid);
  }

}
