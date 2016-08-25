package org.genvisis.cnv.prop;

import org.genvisis.cnv.filesys.Project;

public class DoubleProperty extends Property<Double> {
  double currValue;
  double min;
  double max;

  public DoubleProperty(Project proj, String name, String desc, double min, double max,
                        double defValue) {
    super(proj, name, desc, defValue);
    if (min > max || defValue < min || defValue > max || (max == min && defValue != max)) {
      throw new RuntimeException("Cannot initialize DoubleProperty['" + name + "'] with: min="
                                 + min + ", max=" + max + ", and default value=" + defValue);
    }
    this.min = min;
    this.max = max;
  }

  @Override
  public void parseValue(String valueStr) {
    Double newValue = valueStr.equals("") ? getDefaultValue() : Double.valueOf(valueStr);
    setValue(newValue);
  }

  @Override
  public void setValue(Double value) {
    if (value < min || value > max) {
      throw new RuntimeException("Error - values for property " + getName() + " must be within "
                                 + min + "-" + max + "; " + value + " is not valid");
    }
    super.setValue(value);
  }

  public double getMinValue() {
    return min;
  }

  public double getMaxValue() {
    return max;
  }
}