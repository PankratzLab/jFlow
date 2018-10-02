package org.pankratzlab.common.parsing;

public abstract class CachedFileColumn<T> extends AbstractFileColumn<T> {

  public CachedFileColumn(String nm, boolean dieOnParseFailure) {
    super(nm, dieOnParseFailure);
  }

  private String[] currentLine;
  private T currentValue;

  @Override
  public T getValue(String[] line) throws ParseFailureException {
    if (line == currentLine) {
      return currentValue;
    }
    currentLine = line;
    return currentValue = calculateValue(line);
  }

  public abstract T calculateValue(String[] line) throws ParseFailureException;

  /**
   * Implementations of CachedFileColumn should be sure to call initialize on the column they're
   * wrapping.
   */
  @Override
  public abstract void initialize(FileParser parser);

}
