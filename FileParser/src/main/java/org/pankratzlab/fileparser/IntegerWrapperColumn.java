package org.pankratzlab.fileparser;

public class IntegerWrapperColumn extends NumberWrapperColumn<Integer> {

  public IntegerWrapperColumn(IndexedFileColumn<?> base) {
    this(base, base.dieOnParseFailure());
  }

  public IntegerWrapperColumn(IndexedFileColumn<?> base, boolean dieOnParseFailure) {
    super(base, Integer::parseInt, dieOnParseFailure);
  }
}
