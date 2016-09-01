package org.genvisis.cnv.manage;

import java.util.ArrayList;
import java.util.List;

import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.common.Logger;

public class TestResources {

  private static List<String> failures = new ArrayList<String>();

  private static void test(Resource value) {
    if (value.get() == null) {
      failures.add(value.getLocalPath());
    }
  }

  private static void testBinaries(Logger log) {
    for (Resource r : Resources.miniMac(log).getResources()) {
      test(r);
    }
    for (Resource r : Resources.shapeit(log).getResources()) {
      test(r);
    }
  }

  private static void testGenome(Logger log) {

    for (GENOME_BUILD build : GENOME_BUILD.values()) {
      for (Resource r : Resources.genome(build, log).getResources()) {
        test(r);
      }
    }
  }

  private static void testChr(Logger log) {
    for (GENOME_BUILD build : GENOME_BUILD.values()) {
      for (Resource r : Resources.genome(build, log).chr().getResources()) {
        test(r);
      }
    }
  }

  private static void testAffy(Logger log) {
    for (Resource r : Resources.affy(log).getResources()) {
      test(r);
    }
  }

  private static void testMitoCN(Logger log) {
    for (Resource r : Resources.mitoCN(log).getResources()) {
      test(r);
    }
  }

  public static void main(String[] args) {
    Logger log = new Logger();

    testBinaries(log);
    testGenome(log);
    testChr(log);
    testAffy(log);
    testMitoCN(log);
    System.out.println("--Test complete--");
    if (failures.isEmpty()) {
      System.out.println("---- No failures!");
    }
    else {
      System.out.println("The following resources failed:");
      for (String s : failures) {
        System.out.println(s);
      }
    }
  }
}
