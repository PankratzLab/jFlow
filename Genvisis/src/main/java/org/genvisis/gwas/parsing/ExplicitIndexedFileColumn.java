package org.genvisis.gwas.parsing;

import java.util.Map;

public class ExplicitIndexedFileColumn extends AbstractFileColumn<String> implements IndexedFileColumn<String> {

  final int index;

  public ExplicitIndexedFileColumn(String name, int inputIndex) {
    super(name);
    this.index = inputIndex;
  }

  @Override
  public void initialize(FileParser parser) {
    Map<String, Integer> headerMap = parser.getHeaderMap();
    if (!headerMap.values().contains(this.index)) {
      throw new IllegalStateException("Index " + this.index + " is not valid for header.");
    }
  }

  @Override
  public String getValue(String[] line) {
    return line[index];
  }

  @Override
  public int getIndex() {
    return index;
  }

}
