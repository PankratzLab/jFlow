package org.genvisis;

import org.apache.commons.cli.ParseException;
import org.junit.Before;
import org.junit.Test;

import junit.framework.Assert;

/**
 * Tests for the {@link CLI} class.
 */
public class TestCLI {
  private CLI c;

  @Before
  public void setUp() {
    c = new CLI();
  }

  /**
   * Ensure type checking of arguments works as intended
   */
  @Test
  public void typeTest() {
    c.addArg("test", "An integer argument", true, CLI.Arg.NUMBER);

    // Try parsing a non-integer and ensure it fails
    boolean caught = false;
    try {
      c.parse(getClass(), "test=krakens");
    } catch (ParseException e) {
      caught = true;
    }
    Assert.assertTrue(caught);

    // Verify that parsing a numerical assignment does not fail
    try {
      c.parse(getClass(), "test=523");
    } catch (ParseException e) {
      Assert.fail();
    }
  }

  /**
   * Ensure the {@link CLI#defaultOptions()} includes help commands.
   */
  @Test
  public void helpTest() {
    // These calls should be successful even though we didn't add any options explicitly
    // as they are added in the default options
    for (String f : new String[] {"-h", "-help"}) {
      boolean caught = false;
      try {
        c.parse(getClass(), f);
      } catch (ParseException e) {
        caught = true;
      }
      Assert.assertTrue(caught);
    }
  }

  /**
   * Ensure parsing a flag works with the "-" symbol
   */
  @Test
  public void flagTest() {
    c.addFlag("testFlag", "this is a test flag", true);
    c.addFlag("skippedFlag", "test negative set");

    boolean caught = false;
    try {
      c.parse(getClass(), "-testFlag");
      Assert.assertFalse(c.has("skippedFlag"));
    } catch (ParseException e) {
      caught = true;
    }
    Assert.assertFalse(caught);

  }

  /**
   * Ensure:
   * <ul>
   * <li>parsing an argument works without the "-" symbol</li>
   * <li>the default value of an arg is used if not explicitly set</li>
   * </ul>
   */
  @Test
  public void argTest() {
    final String k1 = "test1";
    final String k2 = "test2";
    final String k3 = "test3";
    final String v1 = "I'm on the command line!";
    final String v2 = "Me too!";

    // Add two argument options, one with a default value and one without.
    c.addArg(k1, "This argument has a default value", v1, true);
    c.addArg(k2, "This argument does not have a default value", true);
    c.addArg(k3, "This argument is not required and does not have a default value");

    try {
      // Try parsing with just k2 set
      c.parse(getClass(), k2 + "=" + v2);

      // k1 should have a value since it had a default value
      Assert.assertEquals(v1, c.get(k1));

      // k2 should have the passed value
      Assert.assertEquals(v2, c.get(k2));

      // k3 should not be in the parsed output set
      Assert.assertFalse(c.has(k3));
      Assert.assertNull(c.get(k3));

    } catch (ParseException e) {
      Assert.fail(e.getMessage());
    }
  }
}
