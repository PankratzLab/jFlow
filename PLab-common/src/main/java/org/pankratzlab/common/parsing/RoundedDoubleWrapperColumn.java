package org.pankratzlab.common.parsing;

import org.pankratzlab.common.ext;

public class RoundedDoubleWrapperColumn extends WrapperColumn<String, Double> {

  private final int numberOfSigDigits;
  private final boolean forceDigits;

  public RoundedDoubleWrapperColumn(String name, IndexedFileColumn<? extends Double> base,
                                    int sigDigits) {
    this(name, base, base.dieOnParseFailure(), sigDigits);
  }

  public RoundedDoubleWrapperColumn(String name, IndexedFileColumn<? extends Double> base,
                                    boolean dieOnParseFailure, int sigDigits) {
    this(name, base, dieOnParseFailure, sigDigits, true);
  }

  public RoundedDoubleWrapperColumn(String name, IndexedFileColumn<? extends Double> base,
                                    boolean dieOnParseFailure, int sigDigits, boolean forceDigits) {
    super(name, base, dieOnParseFailure);
    this.numberOfSigDigits = sigDigits;
    this.forceDigits = forceDigits;
  }

  @Override
  protected String calculateFromBase(Double base) throws ParseFailureException {
    return ext.formDeci(base.doubleValue(), numberOfSigDigits, forceDigits);
  }

}
