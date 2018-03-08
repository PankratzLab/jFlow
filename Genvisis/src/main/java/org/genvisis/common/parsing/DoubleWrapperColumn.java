package org.genvisis.common.parsing;

public class DoubleWrapperColumn extends NumberWrapperColumn<Double> {

  public <T> DoubleWrapperColumn(FileColumn<T> base) {
    this(base, base.dieOnParseFailure());
  }

  public <T> DoubleWrapperColumn(FileColumn<T> base, boolean dieOnParseFailure) {
    super(base, Double::parseDouble, dieOnParseFailure);
  }

}
