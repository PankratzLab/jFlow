package org.genvisis.gwas.parsing;

public class FixedValueColumn extends AbstractFileColumn<String> {

  private String value;

  public FixedValueColumn(String name, String value) {
    super(name);
    this.value = value;
  }

  @Override
  public String getValue(String[] line) throws ParseFailureException {
    return value;
  }

  @Override
  public void initialize(FileParser parser) {
    // no-op
  }

}
