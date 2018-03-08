package org.genvisis.common.parsing;

public class ArrayWrapperColumn extends CachedFileColumn<String[]> {

  private final FileColumn<String> base;
  private final String delimeter;

  /**
   * Construct an {@link ArrayWrapperColumn} that parses into a String[] based on delimiter
   * 
   * @param base underlying column to split into array
   * @param delimeter regex delimiter to split on
   */
  public ArrayWrapperColumn(FileColumn<String> base, String delimeter) {
    super(base.getName(), base.dieOnParseFailure());
    this.base = base;
    this.delimeter = delimeter;
  }

  @Override
  public String[] calculateValue(String[] line) throws ParseFailureException {
    return base.getValue(line).split(delimeter);
  }

  @Override
  public void initialize(FileParser parser) {
    base.initialize(parser);
  }

}
