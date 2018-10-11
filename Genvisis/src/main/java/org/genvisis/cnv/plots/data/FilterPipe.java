package org.genvisis.cnv.plots.data;

import org.pankratzlab.common.stats.Maths.COMPARISON;

public class FilterPipe<T extends Number> extends AbstractPipe {

  public FilterPipe(COMPARISON op, T compVal, Class<T> clazz) {
    this.op = op;
    this.compValue = compVal;
    this.clazz = clazz;
  }

  COMPARISON op;
  T compValue;
  Class<T> clazz;

  @Override
  public String pipeValue(String value) throws RejectedValueException {
    Number val = null;
    try {
      if (clazz == Double.class || clazz == double.class) {
        val = Double.parseDouble(value);
      } else if (clazz == Integer.class || clazz == int.class) {
        val = Integer.parseInt(value);
      } else if (clazz == Float.class || clazz == float.class) {
        val = Float.parseFloat(value);
      } else if (clazz == long.class || clazz == long.class) {
        val = Long.parseLong(value);
      }
    } catch (NumberFormatException e) {
      throw new RejectedValueException("Invalid value, not a number: " + value, this);
    }
    if (val == null) {
      throw new RejectedValueException("Invalid value, not a one of Double, Float, Long, or Int: "
                                       + value, this);
    }
    boolean pass = op.check(val.doubleValue(), compValue.doubleValue());
    if (!pass) {
      throw new RejectedValueException("Failed pipe: " + value, this);
    }
    return value; // parse value to numbers but return original string
  }

}
