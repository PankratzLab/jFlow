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
   * @param req Requirement
   * @param String Raw, unparsed variable value
   */
  public <T> void parseOrFail(Requirement<T> req, String rawValue) {
    try {
      lineValues.put(req, req.parseValue(rawValue));
    } catch (Exception e) {
      lineValues.put(req, FAIL);
    }
  }

  /**
   * Set value
   * 
   * @param <T>
   * @param req Requirement
   * @param T Raw, unparsed variable value
   */
  public <T> void put(Requirement<T> req, T value) {
    lineValues.put(req, value);
  }

  /**
   * Is there a <i>valid</i> value for the given Requirement in this Variables?<br />
   * 
   * @param req
   * @return
   */
  public boolean hasValid(Requirement<?> req) {
    return lineValues.containsKey(req) && lineValues.get(req) != FAIL;
  }

  /**
   * Get the String representation of the value in this line for the given Requirement or null if no
   * value is present.
   * 
   * @param req
   * @return
   */
  public String getString(Requirement<?> req) {
    return hasValid(req) ? lineValues.get(req).toString() : null;
  }

  /**
   * @param req
   * @return
   */
  @SuppressWarnings("unchecked")
  public <T> T get(Requirement<T> req) {
    return hasValid(req) ? (T) lineValues.get(req) : req.getDefaultValue();
  }

  public boolean parseFail(Requirement<?> req) {
    return lineValues.get(req) == FAIL;
  }

  public Set<Requirement<?>> keys() {
    return Collections.unmodifiableSet(lineValues.keySet());
  }

}
