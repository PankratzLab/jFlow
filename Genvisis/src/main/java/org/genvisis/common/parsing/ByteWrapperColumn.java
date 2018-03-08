package org.genvisis.common.parsing;

public class ByteWrapperColumn extends CachedFileColumn<Byte> {

  private FileColumn<?> base;

  public ByteWrapperColumn(FileColumn<?> base, boolean dieOnParseFailure) {
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
  public Byte calculateValue(String[] line) throws ParseFailureException {
    try {
      return Byte.parseByte(getBaseValue(line));
    } catch (NumberFormatException e) {
      throw new ParseFailureException(e);
    }
  }

}
