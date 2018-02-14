package org.genvisis.cnv.prop;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.COPY;
import org.genvisis.cnv.filesys.Project.GROUP;

public class EnumProperty<T extends Enum<T>> extends Property<T> {

  T[] enumValues;
  int defaultIndex;

  public EnumProperty(Project proj, String name, String description, GROUP group, boolean editable,
                      COPY copyOnCorrection, int defaultIndex, Class<T> opts) {
    super(proj, name, description, group, editable, copyOnCorrection,
          opts.getEnumConstants()[defaultIndex]);
    this.enumValues = opts.getEnumConstants();
    this.defaultIndex = defaultIndex;
  }

  public EnumProperty(Project proj, String name, String description, GROUP group, boolean editable,
                      COPY copyOnCorrection, T defaultOpt, Class<T> opts) {
    super(proj, name, description, group, editable, copyOnCorrection, defaultOpt);
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
