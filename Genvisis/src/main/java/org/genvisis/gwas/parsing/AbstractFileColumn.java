package org.genvisis.gwas.parsing;

public abstract class AbstractFileColumn<T> implements FileColumn<T> {

  private final String name;
  private final boolean dieOnParseFail;
  private FileParser parser;

  public AbstractFileColumn(String nm) {
    this(nm, false);
  }

  public AbstractFileColumn(String nm, boolean dieOnParseFailure) {
    this.name = nm;
    this.dieOnParseFail = dieOnParseFailure;
  }

  @Override
  public String getName() {
    return this.name;
  }

  @Override
  public boolean dieOnParseFailure() {
    return dieOnParseFail;
  }

  @Override
  public void initialize(FileParser parser) {
    if (this.parser != parser) {
      if (this.parser == null) {
        this.parser = parser;
      } else {
        throw new IllegalStateException("Initialize called on a new file for column " + getName()
                                        + ".  Current file: " + this.parser.getInputFile()
                                        + "; new file: " + parser.getInputFile());
      }
    }
  }

  @Override
  // overridden to force subclasses to define their own hashCode method
  public abstract int hashCode();

  @Override
  // overridden to force subclasses to define their own equals method
  public abstract boolean equals(Object o);

}
