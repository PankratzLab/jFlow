package org.genvisis.cnv.prop;

import org.genvisis.cnv.filesys.Project;

public class BooleanProperty extends Property<Boolean> {
  public BooleanProperty(Project proj, String name, String desc, Boolean defVal) {
    super(proj, name, desc, defVal);
  }

  @Override
  public void parseValue(String valueStr) {
    Boolean newValue = valueStr.equals("") ? getDefaultValue() : Boolean.valueOf(valueStr);
    setValue(newValue);
  }
}