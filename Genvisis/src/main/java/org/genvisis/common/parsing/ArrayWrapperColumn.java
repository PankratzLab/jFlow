package org.genvisis.common.parsing;

public class ArrayWrapperColumn extends WrapperColumn<String[], String> {

  private final String delimeter;

  /**
   * Construct an {@link ArrayWrapperColumn} that parses into a String[] based on delimiter
   * 
   * @param base underlying column to split into array
   * @param delimeter regex delimiter to split on
   */
  public ArrayWrapperColumn(IndexedFileColumn<String> base, String delimeter) {
    this(base, delimeter, base.dieOnParseFailure());
  }

  /**
   * Construct an {@link ArrayWrapperColumn} that parses into a String[] based on delimiter
   * 
   * @param base underlying column to split into array
   * @param delimeter regex delimiter to split on
   */
  public ArrayWrapperColumn(IndexedFileColumn<String> base, String delimeter,
                            boolean dieOnParseFailure) {
    super(base, dieOnParseFailure);
    this.delimeter = delimeter;
  }

  @Override
  protected String[] calculateFromBase(String base) {
    return base.split(delimeter);
  }

}
