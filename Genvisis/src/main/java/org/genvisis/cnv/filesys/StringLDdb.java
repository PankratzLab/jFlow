package org.genvisis.cnv.filesys;

import java.io.Serializable;
import java.util.Hashtable;
import java.util.List;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.ext;

public class StringLDdb implements Serializable {

  public static final long serialVersionUID = 1L;

  private final Hashtable<String, String> hash;
  private final Hashtable<String, String> missing;
  private final Hashtable<String, String> monomorphs;
  private boolean changed;

  public StringLDdb() {
    hash = new Hashtable<>();
    missing = new Hashtable<>();
    monomorphs = new Hashtable<>();
    changed = true;
  }

  public void add(String marker1, String marker2, float r2) {
    String str;

    str = ext.formDeci(r2, 3, true);
    if (str.startsWith("1")) {
      str = "1";
    } else {
      str = str.substring(2);
    }

    hash.put(marker1 + "\t" + marker2, str);
    changed = true;
  }

  public float get(String marker1, String marker2) {
    String str;

    str = hash.get(marker1 + "\t" + marker2);

    if (str == null) {
      str = hash.get(marker2 + "\t" + marker1);
    }

    if (str == null) {
      if (monomorphs.containsKey(marker1) && monomorphs.containsKey(marker2)) {
        return LDdatabase.MONOMORPH_BOTH;
      } else if (monomorphs.containsKey(marker1)) {
        return LDdatabase.MONOMORPH_FIRST;
      } else if (monomorphs.containsKey(marker2)) {
        return LDdatabase.MONOMORPH_SECOND;
      } else if (missing.containsKey(marker1) && missing.containsKey(marker2)) {
        return LDdatabase.MISSING_BOTH;
      } else if (missing.containsKey(marker1)) {
        return LDdatabase.MISSING_FIRST;
      } else if (missing.containsKey(marker2)) {
        return LDdatabase.MISSING_SECOND;
      } else {
        return LDdatabase.MISSING_INFO;
      }
    } else if (str.equals("1")) {
      return 1f;
    } else {
      return Float.parseFloat("0." + str);
    }
  }

  public void addMonomorphs(List<String> newMonomorphs) {
    for (int i = 0; i < newMonomorphs.size(); i++) {
      monomorphs.put(newMonomorphs.get(i), "");
      changed = true;
    }
  }

  public void addMissing(List<String> newMissings) {
    for (int i = 0; i < newMissings.size(); i++) {
      missing.put(newMissings.get(i), "");
      changed = true;
    }
  }

  public boolean getChanged() {
    return changed;
  }

  public void serialize(String root) {
    changed = false;
    SerializedFiles.writeSerial(this, root + ".slddb");
  }

  public static StringLDdb load(String root, boolean createIfAbsent) {
    if (Files.exists(root + ".slddb")) {
      return (StringLDdb) SerializedFiles.readSerial(root + ".slddb", false);
    } else if (createIfAbsent) {
      return new StringLDdb();
    } else {
      System.err.println("Error - '" + root + ".slddb" + "' does not exist");
      return null;
    }
  }
}
