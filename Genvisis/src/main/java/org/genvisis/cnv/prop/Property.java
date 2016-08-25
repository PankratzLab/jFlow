package org.genvisis.cnv.prop;

import org.genvisis.cnv.filesys.Project;

public abstract class Property<T> {
  private final Project myProj;
  private final String name;
  private final String desc;
  private final T defaultValue;
  private T value;

  public Property(Project proj, String name, String description, T defVal) {
    this.myProj = proj;
    this.name = name;
    this.desc = description;
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
    this.value = value;
  }

  public String getName() {
    return name;
  }

  public String getDescription() {
    return desc;
  }

  public T getDefaultValue() {
    return defaultValue;
  }

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