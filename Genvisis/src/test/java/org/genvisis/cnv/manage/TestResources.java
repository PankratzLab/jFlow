package org.genvisis.cnv.manage;

import java.util.ArrayList;
import java.util.List;

import org.genvisis.cnv.manage.Resources.CHROMASOME;
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
   * Helper method to run {@link #test(Resource)} on a collection of resources.
   */
  private static void test(List<Resource> resources) {
    for (Resource r : resources) {
      test(r);
    }
  }

  /**
   * Helper method to test program downloads
   */
  private static void testBinaries(Logger log) {
    test(Resources.miniMac(log).getResources());
    test(Resources.shapeit(log).getResources());
  }

  /**
   * Helper method to test the {@link Resources#genome(GENOME_BUILD, Logger)} resources.
   */
  private static void testGenome(Logger log) {
    for (GENOME_BUILD build : GENOME_BUILD.values()) {
      test(Resources.genome(build, log).getResources());
    }
  }

  /**
   * Helper method to test the chromasome resources.
   */
  private static void testChr(Logger log) {
    for (GENOME_BUILD build : GENOME_BUILD.values()) {
      for (CHROMASOME c : CHROMASOME.values()) {
        test(Resources.genome(build, log).chr(c).getResources());
      }
    }
  }

  /**
   * Helper method to test {@link Resources#affy(Logger)} resources.
   */
  private static void testAffy(Logger log) {
    test(Resources.affy(log).getResources());

    for (GENOME_BUILD build : GENOME_BUILD.values()) {
      test(Resources.affy(log).genome(build).getResources());
    }
  }

  /**
   * Helper method to test {@link Resources#mitoCN(Logger)} resources.
   */
  private static void testMitoCN(Logger log) {
    test(Resources.mitoCN(log).getResources());
  }

  /**
   * Helper method to test {@link Resources#path(Logger)} resources.
   */
  private static void testPathways(Logger log) {
    test(Resources.path(log).getResources());
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
    testPathways(log);
    System.out.println("--Test complete--");
    if (failures.isEmpty()) {
      System.out.println("---- No failures!");
    } else {
      System.out.println("The following resources failed:");
      for (String s : failures) {
        System.out.println(s);
      }
    }
  }
}
