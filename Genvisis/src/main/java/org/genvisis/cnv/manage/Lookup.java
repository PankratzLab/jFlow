package org.genvisis.cnv.manage;

import java.util.Hashtable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;

public class Lookup {

  public static String lookup(String key, Project project, boolean lookupDNA) {
    SampleData sd = project.getSampleData(0, false);
    String[] lookup = sd.lookup(key);
    return lookup[lookupDNA ? 0 : 1];
  }

  /**
   * @param id Either <code>FID+"\t"+IID</code> or <code>IID</code> only
   * @param project
   * @return
   */
  public static String lookupDNA(String id, Project project) {
    SampleData sd = project.getSampleData(0, false);
    String[] lookup = sd.lookup(id);
    return lookup[0];
  }

  public static Hashtable<String, String> lookupDNABatch(java.util.Collection<String> ids,
      Project project) {
    SampleData sd = project.getSampleData(0, false);
    Hashtable<String, String> lookupTable = new Hashtable<String, String>();
    for (String key : ids) {
      lookupTable.put(key, sd.lookup(key)[0]);
    }
    return lookupTable;
  }

  public static Hashtable<String, String> lookupDNABatch(String[] ids, Project project) {
    SampleData sd = project.getSampleData(0, false);
    Hashtable<String, String> lookupTable = new Hashtable<String, String>();
    for (String key : ids) {
      lookupTable.put(key, sd.lookup(key)[0]);
    }
    return lookupTable;
  }

  public static void lookupFamilyDNAs(String fid, Project project) {
    throw new UnsupportedOperationException();
  }

  public static void lookupFamilyIDs(String fid, Project project) {
    throw new UnsupportedOperationException();
  }

  public static String lookupIDs(String dna, Project project) {
    SampleData sd = project.getSampleData(0, false);
    String[] lookup = sd.lookup(dna);
    return lookup[1];
  }

  public static Hashtable<String, String> lookupIDsBatch(java.util.Collection<String> dnas,
      Project project) {
    SampleData sd = project.getSampleData(0, false);
    Hashtable<String, String> lookupTable = new Hashtable<String, String>();
    for (String key : dnas) {
      lookupTable.put(key, sd.lookup(key)[1]);
    }
    return lookupTable;
  }

  public static Hashtable<String, String> lookupIDsBatch(String[] dnas, Project project) {
    SampleData sd = project.getSampleData(0, false);
    Hashtable<String, String> lookupTable = new Hashtable<String, String>();
    for (String key : dnas) {
      lookupTable.put(key, sd.lookup(key)[1]);
    }
    return lookupTable;
  }

  public static String lookupIIDOnly(String dna, Project project) {
    SampleData sd = project.getSampleData(0, false);
    String[] lookup = sd.lookup(dna);
    return lookup[2];
  }


  public static Hashtable<String, String> lookupIIDsOnlyBatch(java.util.Collection<String> dnas,
      Project project) {
    SampleData sd = project.getSampleData(0, false);
    Hashtable<String, String> lookupTable = new Hashtable<String, String>();
    for (String key : dnas) {
      lookupTable.put(key, sd.lookup(key)[2]);
    }
    return lookupTable;
  }

  public static Hashtable<String, String> lookupIIDsOnlyBatch(String[] dnas, Project project) {
    SampleData sd = project.getSampleData(0, false);
    Hashtable<String, String> lookupTable = new Hashtable<String, String>();
    for (String key : dnas) {
      lookupTable.put(key, sd.lookup(key)[2]);
    }
    return lookupTable;
  }

}
