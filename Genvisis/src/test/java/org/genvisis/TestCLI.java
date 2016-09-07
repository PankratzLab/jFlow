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
    c = new CLI(getClass());
  }

  /**
   * Ensure type checking of arguments works as intended
   */
  @Test
  public void typeTest() {
    c.addArg("test", "An integer argument", true, CLI.Arg.NUMBER);

    // Try parsing a non-integer and ensure it fails
    try {
      c.parse("test=krakens");
      Assert.fail("Parsing number as string did not fail");
    } catch (ParseException e) {
      // expected
    }

    // Verify that parsing a numerical assignment does not fail
    try {
      c.parse("test=523");
    } catch (ParseException e) {
      Assert.fail();
    }
  }

  /**
   * Ensure the default {@link CLI} includes help commands.
   */
  @Test
  public void helpTest() {
    // These calls should be successful even though we didn't add any options explicitly
    // as they are added in the default options
    for (String f : new String[] {"-h", "-help"}) {
      try {
        c.parse(f);
        Assert.fail("Parsing continued after printing help");
      } catch (ParseException e) {
        // expected
      }
    }

    // this parse is just to ensure the help message looks OK
    c.addArg("group1", "This arg is in a group");
    c.addArg("group2", "This arg is also in a group");
    c.addGroup("group1", "group2");
    c.addArg("required", "This arg is required", true);
    c.addArgWithDefault("default", "This arg has a default value", "this is the value");
    try {
      c.parse("-h");
    } catch (ParseException exc) {
      // expected
    }
  }

  /**
   * Ensure parsing a flag works with the "-" symbol
   */
  @Test
  public void flagTest() {
    c.addFlag("testFlag", "this is a test flag", true);
    c.addFlag("skippedFlag", "test negative set");

    try {
      c.parse("-testFlag");
      Assert.assertFalse(c.has("skippedFlag"));
    } catch (ParseException e) {
      Assert.fail(e.getMessage());
    }

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
    c.addArgWithDefault(k1, "This argument has a default value", v1);
    c.addArg(k2, "This argument does not have a default value", true);
    c.addArg(k3, "This argument is not required and does not have a default value");

    try {
      // Try parsing with just k2 set
      c.parse(k2 + "=" + v2);

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

  @Test
  public void groupTest() {
    final String k1 = "key1";
    final String k2 = "key2";
    final String k3 = "key3";
    c.addArg(k1, "Option 1");
    c.addFlag(k2, "Option 2");
    c.addArgWithDefault(k3, "Option 3", "default value");
    c.addGroup(k1, k2, k3);

    // This should be fine
    try {
      c.parse(k1 + "=hello");
      Assert.assertTrue(c.has(k1));
      // Even though a default value was given, it should have been removed due to addition to the group
      Assert.assertFalse(c.has(k3));
    } catch (ParseException exc) {
      Assert.fail(exc.getMessage());
    }

    try {
      c.parse("-" + k2);
      Assert.assertTrue(c.has(k2));
    } catch (ParseException exc) {
      Assert.fail(exc.getMessage());
    }


    try {
      c.parse("-" + k2, k1 + "=hello");
      Assert.fail("Parsing mutually exclusive options succeeded");
    } catch (ParseException exc) {
      // expected
    }

  }
}
