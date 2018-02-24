package org.genvisis.gwas.parsing;

public class DoubleWrapperColumn extends CachedFileColumn<Double> {

  private FileColumn<?> base;

  public DoubleWrapperColumn(FileColumn<?> base) {
    this(base, base.dieOnParseFailure());
  }

  public DoubleWrapperColumn(FileColumn<?> base, boolean dieOnParseFailure) {
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
  public Double calculateValue(String[] line) throws ParseFailureException {
    try {
      return Double.parseDouble(getBaseValue(line));
    } catch (NumberFormatException e) {
      throw new ParseFailureException(e);
    }
  }

}
