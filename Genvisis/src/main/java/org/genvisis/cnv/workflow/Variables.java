package org.genvisis.cnv.workflow;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class Variables {

  private Map<Requirement<?>, Object> lineValues;
  private static final Object FAIL = new Object();

  /**
   * Constructor
   * 
   * @param defaultFailValue
   */
  public Variables() {
    lineValues = new HashMap<>();
  }

  /**
   * Set value if valid, or set as a parse failure
   * 
   * @param fc Requirement
   * @param String Raw, unparsed variable value
   */
  public <T> void parseOrFail(Requirement<T> fc, String rawValue) {
    try {
      lineValues.put(fc, fc.parseValue(rawValue));
    } catch (Exception e) {
      lineValues.put(fc, FAIL);
    }
  }

  /**
   * Set value
   * 
   * @param <T>
   * @param fc Requirement
   * @param T Raw, unparsed variable value
   */
  public <T> void put(Requirement<T> fc, T value) {
    lineValues.put(fc, value);
  }

  /**
   * Is there <i>any</i> value for the given Requirement in this DataLine?<br />
   * (May be a parse failure value, but at least a value is present.)
   * 
   * @param fc
   * @return
   */
  public boolean has(Requirement<?> fc) {
    return lineValues.containsKey(fc) && lineValues.get(fc) != FAIL;
  }

  /**
   * Get the String representation of the value in this line for the given Requirement. Calls
   * {@link Object#toString()}, so if no value is present a {@link NullPointerException} will be
   * thrown.
   * 
   * @param fc
   * @return
   */
  public String getString(Requirement<?> fc) {
    return (has(fc) ? lineValues.get(fc) : null).toString();
  }

  /**
   * @param fc
   * @return
   */
  @SuppressWarnings("unchecked")
  public <T> T get(Requirement<T> fc) {
    return has(fc) ? (T) lineValues.get(fc) : fc.getDefaultValue();
  }

  public boolean parseFail(Requirement<?> fc) {
    return lineValues.get(fc) == FAIL;
  }

  public Set<Requirement<?>> keys() {
    return Collections.unmodifiableSet(lineValues.keySet());
  }

}
