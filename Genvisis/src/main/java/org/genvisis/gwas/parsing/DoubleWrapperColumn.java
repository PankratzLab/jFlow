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

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((getName() == null) ? 0 : getName().hashCode());
    result = prime * result + base.hashCode();
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    DoubleWrapperColumn other = (DoubleWrapperColumn) obj;
    if (!base.equals(other.base)) return false;
    if (getName() == null) {
      if (other.getName() != null) return false;
    } else if (!getName().equals(other.getName())) return false;
    return true;
  }

}
