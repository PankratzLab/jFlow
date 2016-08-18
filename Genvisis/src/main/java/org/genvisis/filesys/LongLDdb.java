package org.genvisis.filesys;

import java.io.Serializable;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;

public class LongLDdb implements Serializable {
  public static final long serialVersionUID = 1L;
  public static final long DUMMY_BASE = 999000000;

  private final Hashtable<Long, String> hash;
  private final Hashtable<String, Long> lookup;
  private final Hashtable<String, String> missing;
  private final Hashtable<String, String> monomorphs;
  private int count;
  private boolean changed;

  public LongLDdb() {
    hash = new Hashtable<Long, String>();
    lookup = new Hashtable<String, Long>();
    missing = new Hashtable<String, String>();
    monomorphs = new Hashtable<String, String>();
    count = 0;
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

    hash.put(convert(marker1, marker2), str);
    changed = true;
  }

  public float get(String marker1, String marker2) {
    String str;

    str = hash.get(convert(marker1, marker2));

    if (str == null) {
      str = hash.get(convert(marker2, marker1));
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

  public void addMonomorphs(Vector<String> newMonomorphs) {
    for (int i = 0; i < newMonomorphs.size(); i++) {
      monomorphs.put(newMonomorphs.elementAt(i), "");
      changed = true;
    }
  }

  public void addMissing(Vector<String> newMissings) {
    for (int i = 0; i < newMissings.size(); i++) {
      missing.put(newMissings.elementAt(i), "");
      changed = true;
    }
  }

  public boolean getChanged() {
    return changed;
  }

  public Long convert(String marker1, String marker2) {
    return convert(new String[] {marker1, marker2});
  }

  public Long convert(String[] markers) {
    long[] ls;

    ls = new long[2];
    for (int i = 0; i < 2; i++) {
      if (!markers[i].startsWith("rs") || markers[i].length() > 11) {
        if (lookup.containsKey(markers[i])) {
          ls[i] = lookup.get(markers[i]);
        } else {
          count++;
          ls[i] = DUMMY_BASE + count;
        }
      } else {
        ls[i] = Long.parseLong(markers[i].substring(2));
      }

    }

    return Long.valueOf(ls[0] * 1000000000 + ls[1]);
  }

  public void serialize(String root) {
    SerializedFiles.writeSerial(this, root + ".llddb");
  }

  public static LongLDdb load(String root, boolean jar, boolean createIfAbsent) {
    if (Files.exists(root + ".llddb", jar)) {
      return (LongLDdb) SerializedFiles.readSerial(root + ".llddb", jar, false);
    } else if (createIfAbsent) {
      return new LongLDdb();
    } else {
      System.err.println("Error - '" + root + ".llddb" + "' does not exist");
      return null;
    }
  }
}
