package org.genvisis.common.parsing;

public class IntegerWrapperColumn extends CachedFileColumn<Integer> {

  private FileColumn<?> base;

  public IntegerWrapperColumn(FileColumn<?> base) {
    this(base, base.dieOnParseFailure());
  }

  public IntegerWrapperColumn(FileColumn<?> base, boolean dieOnParseFailure) {
    super(base.getName(), dieOnParseFailure);
    this.base = base;
  }

  @Override
  public void initialize(FileParser parser) {
    base.initialize(parser);
  }

  public String getBaseValue(String[] line) throws ParseFailureException {
    return base.getValue(line).toString();
  }

  @Override
  public Integer calculateValue(String[] line) throws ParseFailureException {
    try {
      return Integer.parseInt(getBaseValue(line));
    } catch (NumberFormatException e) {
      throw new ParseFailureException(e);
    }
  }

}
