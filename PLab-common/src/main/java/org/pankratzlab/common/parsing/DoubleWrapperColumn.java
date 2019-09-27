package org.pankratzlab.common.parsing;

import org.pankratzlab.fileparser.IndexedFileColumn;
import org.pankratzlab.fileparser.NumberWrapperColumn;

public class DoubleWrapperColumn extends NumberWrapperColumn<Double> {

  public <T> DoubleWrapperColumn(IndexedFileColumn<T> base) {
    this(base, base.dieOnParseFailure());
  }

  public <T> DoubleWrapperColumn(IndexedFileColumn<T> base, boolean dieOnParseFailure) {
    super(base, Double::parseDouble, dieOnParseFailure);
  }

}
