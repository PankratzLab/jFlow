package org.genvisis.cnv.prop;

import org.genvisis.cnv.filesys.Project;

public class EnumProperty<T extends Enum<T>> extends Property<T> {
  T[] enumValues;
  int defaultIndex;

  public EnumProperty(Project proj, String name, String description, int defaultIndex,
                      Class<T> opts) {
    super(proj, name, description, opts.getEnumConstants()[defaultIndex]);
    this.enumValues = opts.getEnumConstants();
    this.defaultIndex = defaultIndex;
  }

  public EnumProperty(Project proj, String name, String description, T defaultOpt,
                      Class<T> opts) {
    super(proj, name, description, defaultOpt);
    this.enumValues = opts.getEnumConstants();
    for (int i = 0; i < enumValues.length; i++) {
      if (defaultOpt == enumValues[i]) {
        this.defaultIndex = i;
        break;
      }
    }
  }

  @Override
  public void parseValue(String valueStr) {
    for (T val : enumValues) {
      if (val.toString().equals(valueStr)) {
        setValue(val);
      }
    }
  }
}