package org.genvisis.cnv.prop;

import org.genvisis.cnv.filesys.Project;

public class StringProperty extends Property<String> {
  public StringProperty(Project proj, String name, String desc, String defVal) {
    super(proj, name, desc, defVal);
  }

  @Override
  public void parseValue(String valueStr) {
    setValue(valueStr);
  }
}