package org.pankratzlab.common.parsing;

import org.pankratzlab.common.ext;

public class RoundedDoubleWrapperColumn extends WrapperColumn<String, Double> {

  private final int numberOfSigDigits;

  public RoundedDoubleWrapperColumn(IndexedFileColumn<? extends Double> base, int sigDigits) {
    this(base, base.dieOnParseFailure(), sigDigits);
  }

  public RoundedDoubleWrapperColumn(IndexedFileColumn<? extends Double> base,
                                    boolean dieOnParseFailure, int sigDigits) {
    super(base, dieOnParseFailure);
    this.numberOfSigDigits = sigDigits;
  }

  @Override
  protected String calculateFromBase(Double base) throws ParseFailureException {
    return ext.formDeci(base.doubleValue(), numberOfSigDigits);
  }

}
