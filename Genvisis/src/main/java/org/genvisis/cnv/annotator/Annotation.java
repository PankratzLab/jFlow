package org.genvisis.cnv.annotator;

/**
 * Wrapper class representing the name and value of an annotation. Holds on optional "suffix" as
 * well, which is used as a modifier on the base annotation name.
 */
public class Annotation {

  // Name of this annotation:
  private final String annoName;
  // Value of this annotation 
  private final String annoValue;

  /**
   * Make an annotation with default (empty) suffix
   */
  public Annotation(String annoName, String annoValue) {
    this(annoName, annoValue, "");
  }

  public Annotation(String annoName, String annoValue, String suffix) {
    this.annoName = annoName + suffix;
    this.annoValue = annoValue;
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((annoName == null) ? 0 : annoName.hashCode());
    result = prime * result + ((annoValue == null) ? 0 : annoValue.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) { // Note: equals method makes sure the value matches up to the key/hashcode.
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    Annotation other = (Annotation) obj;
    if (annoName == null) {
      if (other.annoName != null) return false;
    } else if (!annoName.equals(other.annoName)) return false;
    if (annoValue == null) {
      if (other.annoValue != null) return false;
    } else if (!annoValue.equals(other.annoValue)) return false;
    return true;
  }

  // Getter for annoName:
  public String getAnnotationName() {
    return this.annoName;
  }

  // Getter for annoValue:
  public String getAnnotationValue() {
    return this.annoValue;
  }
}
