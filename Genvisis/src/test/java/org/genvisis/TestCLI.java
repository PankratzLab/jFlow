package org.genvisis;

import java.util.Map;

import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PatternOptionBuilder;
import org.junit.Test;

import junit.framework.Assert;

/**
 * Tests for the {@link CLI} class.
 */
public class TestCLI {

  /**
   * Ensure type checking of arguments works as intended
   */
  @Test
  public void testTypes() {
    Options options = CLI.defaultOptions();
    CLI.addArg(options, "test", "An integer argument", true, PatternOptionBuilder.NUMBER_VALUE);

    // Try parsing a non-integer and ensure it fails
    boolean caught = false;
    try {
      CLI.parse(getClass(), options, "test=krakens");
    } catch (ParseException e) {
      caught = true;
    }
    Assert.assertTrue(caught);

    // Verify that parsing a numerical assignment does not fail
    try {
      CLI.parse(getClass(), options, "test=523");
    } catch (ParseException e) {
      Assert.fail();
    }
  }

  /**
   * Ensure the {@link CLI#defaultOptions()} includes help commands.
   */
  @Test
  public void testHelp() {
    Options options = CLI.defaultOptions();
    // These calls should be successful even though we didn't add any options explicitly
    // as they are added in the default options
    for (String f : new String[] {"-h", "-help"}) {
      boolean caught = false;
      try {
        CLI.parse(getClass(), options, f);
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
  public void testFlags() {
    Options options = CLI.defaultOptions();
    CLI.addFlag(options, "testFlag", "this is a test flag", true);

    boolean caught = false;
    try {
      CLI.parse(getClass(), options, "-testFlag");
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
  public void testArgs() {
    Options options = CLI.defaultOptions();
    final String k1 = "test1";
    final String k2 = "test2";
    final String k3 = "test3";
    final String v1 = "I'm on the command line!";
    final String v2 = "Me too!";

    // Add two argument options, one with a default value and one without.
    CLI.addArg(options, k1, "This argument has a default value", v1, true);
    CLI.addArg(options, k2, "This argument does not have a default value", true);
    CLI.addArg(options, k3, "This argument is not required and does not have a default value");

    try {
      // Try parsing with just k2 set
      Map<String, String> parse = CLI.parse(getClass(), options, k2 + "=" + v2);

      // k1 should have a value since it had a default value
      Assert.assertEquals(v1, parse.get(k1));

      // k2 should have the passed value
      Assert.assertEquals(v2, parse.get(k2));

      // k3 should not be in the parsed output set
      Assert.assertNull(parse.get(k3));

    } catch (ParseException e) {
      Assert.fail(e.getMessage());
    }
  }
}
