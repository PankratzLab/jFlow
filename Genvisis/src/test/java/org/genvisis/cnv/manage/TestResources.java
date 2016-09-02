package org.genvisis.cnv.manage;

import java.util.ArrayList;
import java.util.List;

import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.common.Logger;

/**
 * Test class for {@link Resources}.
 */
public class TestResources {

  // Keep a list of all failures so we can summarize at the end
  private static List<String> failures = new ArrayList<String>();

  /**
   * Try to get the resource and record failure.
   *
   * @param value Resource to test
   */
  private static void test(Resource value) {
    if (value.get() == null) {
      failures.add(value.getLocalPath());
    }
  }

  /**
   * Helper method to test program downloads
   */
  private static void testBinaries(Logger log) {
    for (Resource r : Resources.miniMac(log).getResources()) {
      test(r);
    }
    for (Resource r : Resources.shapeit(log).getResources()) {
      test(r);
    }
  }

  /**
   * Helper method to test the {@link Resources#genome(GENOME_BUILD, Logger)} resources.
   */
  private static void testGenome(Logger log) {
    for (GENOME_BUILD build : GENOME_BUILD.values()) {
      for (Resource r : Resources.genome(build, log).getResources()) {
        test(r);
      }
    }
  }

  /**
   * Helper method to test the chromasome resources.
   */
  private static void testChr(Logger log) {
    for (GENOME_BUILD build : GENOME_BUILD.values()) {
      for (Resource r : Resources.genome(build, log).chr().getResources()) {
        test(r);
      }
    }
  }

  /**
   * Helper method to test {@link Resources#affy(Logger)} resources.
   */
  private static void testAffy(Logger log) {
    for (Resource r : Resources.affy(log).getResources()) {
      test(r);
    }
  }

  /**
   * Helper method to test {@link Resources#mitoCN(Logger)} resources.
   */
  private static void testMitoCN(Logger log) {
    for (Resource r : Resources.mitoCN(log).getResources()) {
      test(r);
    }
  }

  /**
   * Run all the tests and report the final outcome.
   */
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
