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
  public int hashCode() {
    final int prime = 31;
    return prime + ((getName() == null) ? 0 : getName().hashCode());
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof FileColumn)) return false;
    FileColumn<?> other = (FileColumn<?>) obj;
    if (getName() == null) {
      if (other.getName() != null) return false;
    } else if (!getName().equals(other.getName())) return false;
    return true;
  }

}
