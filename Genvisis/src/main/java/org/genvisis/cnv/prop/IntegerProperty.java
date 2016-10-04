package org.genvisis.cnv.prop;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.GROUP;

public class IntegerProperty extends Property<Integer> {
  int currValue;
  int min, max;

  public IntegerProperty(Project proj, String name, String desc, GROUP group, boolean editable, int min, int max, int defValue) {
    super(proj, name, desc, group, editable, defValue);
    if (min > max || defValue < min || defValue > max || (max == min && defValue != max)) {
      throw new RuntimeException("Cannot initialize IntegerProperty with: min=" + min + ", max="
                                 + max + ", and default value=" + defValue);
    }
    this.min = min;
    this.max = max;
  }

  public int getMinValue() {
    return min;
  }

  public int getMaxValue() {
    return max;
  }
  
  @Override
  public void parseValue(String valueStr) {
    Integer newValue = valueStr.equals("") ? getDefaultValue() : Integer.valueOf(valueStr);
    setValue(newValue);
  }

  @Override
  public void setValue(Integer value) {
    if (value < min || value > max) {
      throw new RuntimeException("Error - values for property " + getName() + " must be within "
                                 + min + "-" + max + "; " + value + " is not valid");
    }
    super.setValue(value);
  }
}