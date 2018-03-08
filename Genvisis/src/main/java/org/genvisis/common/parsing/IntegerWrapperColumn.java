package org.genvisis.common.parsing;

public class IntegerWrapperColumn extends NumberWrapperColumn<Integer> {

  public IntegerWrapperColumn(FileColumn<?> base) {
    this(base, base.dieOnParseFailure());
  }

  public IntegerWrapperColumn(FileColumn<?> base, boolean dieOnParseFailure) {
    super(base, Integer::parseInt, dieOnParseFailure);
  }
}
