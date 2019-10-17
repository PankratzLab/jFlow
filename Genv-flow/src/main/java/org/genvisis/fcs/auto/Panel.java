package org.genvisis.fcs.auto;

import java.util.HashSet;
import java.util.Set;

public class Panel {
  String name;
  Set<String> aliases;

  public Panel(String name, String... possibleAliases) {
    this.name = name;
    this.aliases = new HashSet<>();
    for (String s : possibleAliases) {
      aliases.add(s);
    }
  }

  public String getName() {
    return name;
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((aliases == null) ? 0 : aliases.hashCode());
    result = prime * result + ((name == null) ? 0 : name.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    Panel other = (Panel) obj;
    if (aliases == null) {
      if (other.aliases != null) return false;
    } else if (!aliases.equals(other.aliases)) return false;
    if (name == null) {
      if (other.name != null) return false;
    } else if (!name.equals(other.name)) return false;
    return true;
  }

}