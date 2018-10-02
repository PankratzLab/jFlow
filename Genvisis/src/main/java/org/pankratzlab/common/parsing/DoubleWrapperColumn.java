package org.pankratzlab.common.parsing;

public class DoubleWrapperColumn extends NumberWrapperColumn<Double> {

  public <T> DoubleWrapperColumn(IndexedFileColumn<T> base) {
    this(base, base.dieOnParseFailure());
  }

  public <T> DoubleWrapperColumn(IndexedFileColumn<T> base, boolean dieOnParseFailure) {
    super(base, Double::parseDouble, dieOnParseFailure);
  }

}
