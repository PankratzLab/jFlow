package org.pankratzlab.common;

import java.io.File;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;

/**
 * Static utility class for accessing bundled resources in a consistent manner.
 */
public final class Bundled {

  private Bundled() {
    // Prevent instantiation of utility class
  }

  /**
   * Load a resource using {@link Bundled} as the base class.
   *
   * @param path Path to the local resource
   * @return URL to the requested resource, or null if not found.
   * @see #get(Class, String)
   */
  public static URL get(String path) {
    return get(Bundled.class, path);
  }

  /**
   * Load a packaged resource using the specified base class.
   *
   * @param base Class to use as the base for loading the requested resource
   * @param path Path to the local resource
   * @return URL to the requested resource, or null if not found.
   */
  public static URL get(Class<?> base, String path) {
    return base.getClassLoader().getResource(path);
  }

  /**
   * Load a resource using {@link Bundled} as the base class.
   *
   * @param path Path to the local resource
   * @return A File reference to the loaded resource, or null if not found.
   * @see #getFile(Class, String)
   */
  public static File getFile(String path) {
    return getFile(Bundled.class, path);
  }

  /**
   * Load a packaged resource using the specified base class.
   *
   * @param base Class to use as the base for loading the requested resource
   * @param path Path to the local resource
   * @return A File reference to the loaded resource, or null if not found.
   */
  public static File getFile(Class<?> base, String path) {
    URL url = get(base, path);
    if (url == null) {
      return null;
    }
    File f = null;
    // From https://community.oracle.com/blogs/kohsuke/2007/04/25/how-convert-javaneturl-javaiofile
    try {
      f = new File(url.toURI());
    } catch (URISyntaxException e) {
      f = new File(url.getPath());
    }
    return f;
  }

  /**
   * Load a resource using {@link Bundled} as the base class.
   *
   * @param path Path to the local resource
   * @return InputStream to the requested resource, or null if not found.
   * @see #getStream(Class, String)
   */
  public static InputStream getStream(String path) {
    return getStream(Bundled.class, path);
  }

  /**
   * Load a packaged resource using the specified base class.
   *
   * @param base Class to use as the base for loading the requested resource
   * @param path Path to the local resource
   * @return InputStream to the requested resource, or null if not found.
   */
  public static InputStream getStream(Class<?> base, String path) {
    return base.getClassLoader().getResourceAsStream(path);
  }
}
