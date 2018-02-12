package org.genvisis.gwas.parsing;

import java.util.Map;

public class IndexedFileColumn extends AbstractFileColumn<String> {

  final int index;

  public IndexedFileColumn(String name, int inputIndex) {
    super(name);
    this.index = inputIndex;
  }

  @Override
  public void initialize(Map<String, Integer> headerMap) {
    if (!headerMap.values().contains(this.index)) {
      throw new IllegalStateException("Index " + this.index + " is not valid for header.");
    }
  }

  @Override
  public String getValue(String[] line) {
    return line[index];
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((getName() == null) ? 0 : getName().hashCode());
    result = prime * result + index;
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    IndexedFileColumn other = (IndexedFileColumn) obj;
    if (index != other.index) return false;
    if (getName() == null) {
      if (other.getName() != null) return false;
    } else if (!getName().equals(other.getName())) return false;
    return true;
  }

}
