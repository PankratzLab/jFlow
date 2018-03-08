package org.genvisis.common.parsing;

public class StringWrapperColumn extends WrapperColumn<String, Object> {

  public StringWrapperColumn(FileColumn<? extends Object> base) {
    super(base);
  }

  public StringWrapperColumn(FileColumn<? extends Object> base, boolean dieOnParseFailure) {
    super(base, dieOnParseFailure);
  }

  @Override
  protected String calculateFromBase(Object base) {
    return base.toString();
  }

}
