package org.genvisis.gwas.parsing;

import java.util.Map;

public class ByteWrapperColumn extends CachedFileColumn<Byte> {

  private FileColumn<?> base;

  public ByteWrapperColumn(FileColumn<?> base, boolean dieOnParseFailure) {
    super(base.getName(), dieOnParseFailure);
    this.base = base;
  }

  @Override
  public void initialize(Map<String, Integer> headerMap) {
    base.initialize(headerMap);
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
    ByteWrapperColumn other = (ByteWrapperColumn) obj;
    if (!base.equals(other.base)) return false;
    if (getName() == null) {
      if (other.getName() != null) return false;
    } else if (!getName().equals(other.getName())) return false;
    return true;
  }

}
