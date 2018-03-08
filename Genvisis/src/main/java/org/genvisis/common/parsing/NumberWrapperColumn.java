package org.genvisis.common.parsing;

import java.util.function.Function;

class NumberWrapperColumn<N extends Number> extends WrapperColumn<N, String> {

  private final Function<String, N> parseFunction;

  public NumberWrapperColumn(FileColumn<?> base, Function<String, N> parseFunction) {
    this(new StringWrapperColumn(base), parseFunction, base.dieOnParseFailure());
  }

  public NumberWrapperColumn(FileColumn<?> base, Function<String, N> parseFunction,
                             boolean dieOnParseFailure) {
    super(new StringWrapperColumn(base, dieOnParseFailure), dieOnParseFailure);
    this.parseFunction = parseFunction;
  }

  @Override
  protected N calculateFromBase(String base) throws ParseFailureException {
    try {
      return parseFunction.apply(base);
    } catch (NumberFormatException e) {
      throw new ParseFailureException(e);
    }
  }
}
