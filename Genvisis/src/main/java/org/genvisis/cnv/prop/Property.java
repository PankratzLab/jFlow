package org.genvisis.cnv.prop;

import java.beans.PropertyChangeEvent;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.GROUP;

public abstract class Property<T> {
  private final Project myProj;
  private final String name;
  private final String desc;
  private final GROUP group;
  private final boolean editable;
  private final T defaultValue;
  private T value;

  public Property(Project proj, String name, String description, GROUP group, boolean editable, T defVal) {
    this.myProj = proj;
    this.name = name;
    this.desc = description;
    this.group = group;
    this.editable = editable;
    this.defaultValue = defVal;
    this.value = defaultValue;
  }

  public Project getProject() {
    return myProj;
  }

  public T getValue() {
    return value;
  }

  public void setValue(T value) {
    T oldValue = this.value;
    this.value = value;
    myProj.propertyChange(new PropertyChangeEvent(this, this.getName(), oldValue, this.value));
  }

  public String getName() {
    return name;
  }

  public String getDescription() {
    return desc;
  }
  
  public GROUP getGroup() {
    return group;
  }
  
  public boolean isEditable() {
    return editable;
  }

  public T getDefaultValue() {
    return defaultValue;
  }

  /**
   * Parse a String representation of whatever <code>type</code> of value this Property holds.<br />
   * All implementations of this method should use <code>setValue()</code> after parsing.
   * @param valueStr
   */
  public abstract void parseValue(String valueStr);

  public String getValueString() {
    return value.toString();
  }

  public String getDefaultValueString() {
    return defaultValue.toString();
  }

  @Override
  public String toString() {
    return getName() + "=" + getValueString();
  }
  
}