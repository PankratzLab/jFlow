package org.pankratzlab.common;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.util.function.Function;
import java.util.jar.Attributes;
import java.util.jar.JarFile;
import java.util.jar.Manifest;
import javax.annotation.Nullable;
import com.github.zafarkhaja.semver.Version;

/**
 * Describes the manifest of the jar containing the class that was used to launch this application.
 * Also provides the {@link #setLaunchClass(Class)} method to register the component that should be
 * used for version reporting.
 */
public class LauncherManifest {

  public static final String UNDETERMINED_VERSION = "0.0.0-unknown";
  public static final String BUILD_VERSION = "Implementation-Version";
  public static final String BUILD_TIMESTAMP = "Implementation-Date";
  public static final String COMPILE_TIME = "Compile-Time";
  public static final String BUILD_TYPE = "Build-Type";
  public static final String BUILD_AUTHOR = "Built-By";
  public static final String BUILD_COPYRIGHT = "Copyright";
  public static final String BUILD_COMMIT_SHA = "Implementation-Build";

  private static final Class<?> DEFAULT_LAUNCH_CLASS = LauncherManifest.class;

  private static Class<?> launchClass = null;

  /**
   * @return The class that should be used to determine version information for this application.
   */
  public static Class<?> getLaunchClass() {
    if (launchClass == null) {
      setLaunchClass(DEFAULT_LAUNCH_CLASS);
    }
    return launchClass;
  }

  /**
   * Set the "dominant" launch class, used to determine version information.
   *
   * @return If the requested class was successfully set. If not, it means the current entry point
   *         loaded a different class that called this method first.
   */
  public static boolean setLaunchClass(Class<?> c) {
    if (launchClass == null) {
      return lockedSetLaunchClass(c);
    }
    return launchClass == c;
  }

  /**
   * Synchronized helper method to ensure only the first request to set the launch class is
   * respected.
   */
  private static synchronized boolean lockedSetLaunchClass(Class<?> c) {
    if (launchClass == null) {
      launchClass = c;
    }
    return launchClass == c;
  }

  private Attributes attributes;

  public LauncherManifest() {
    this(new Attributes());
  }

  public LauncherManifest(Attributes attributes) {
    super();
    this.attributes = attributes;
  }

  public Attributes getAttributes() {
    return attributes;
  }

  public Version getVersion() {
    String v = null;
    if (attributes != null) {
      v = attributes.getValue(BUILD_VERSION);
    }

    v = v == null ? UNDETERMINED_VERSION : v;

    return Version.valueOf(v);
  }

  public String getBuiltBy() {
    return getAttribute(BUILD_AUTHOR);
  }

  public String getCopyright() {
    return getAttribute(BUILD_COPYRIGHT);
  }

  public String getBuildType() {
    return getAttribute(BUILD_TYPE);
  }

  public String getCompileTime() {
    return getAttribute(COMPILE_TIME);
  }

  public String getTimestamp() {
    return getAttribute(BUILD_TIMESTAMP);
  }

  public String getCommit() {
    return getAttribute(BUILD_COMMIT_SHA);
  }

  /**
   * Helper method to look up the {@link Attributes#getValue(String)} of a given key. Converts
   * {@code null} attributes and values to empty strings.
   */
  private String getAttribute(String key) {
    if (attributes == null) {
      return "";
    }

    String v = attributes.getValue(key);

    return v == null ? "" : v;
  }

  /**
   * @return The manifest for the jar file used to launch this application
   */
  public static LauncherManifest loadGenvisisManifest() {
    return validateAndLoadManifest(null);
  }

  /**
   * As {@link #loadManifest(File)}, but validates the file with the given method
   */
  public static LauncherManifest validateAndLoadManifest(@Nullable Function<File, File> fileValidator) {
    File file = getCurrentFile();
    if (fileValidator != null) {
      file = fileValidator.apply(file);
    }
    return loadManifest(file);
  }

  /**
   * @return A representation of the manifest of the given jar file
   */
  public static LauncherManifest loadManifest(File jarFile) {
    JarFile jar = null;

    LauncherManifest currentManifest = new LauncherManifest();
    try {

      if (jarFile != null && jarFile.exists()) {
        jar = new JarFile(jarFile);

        Manifest manifest = jar.getManifest();

        Attributes attributes = manifest.getMainAttributes();

        jar.close();
        currentManifest = new LauncherManifest(attributes);
      }
    } catch (IOException ioe) {

    }
    return currentManifest;
  }

  /**
   * @return Get the jar file that was used to launch this application
   */
  private static File getCurrentFile() {
    String jarPath = getLaunchClass().getProtectionDomain().getCodeSource().getLocation().getFile();
    try {
      jarPath = URLDecoder.decode(jarPath, StandardCharsets.UTF_8.name());
    } catch (UnsupportedEncodingException e) {
      e.printStackTrace();
    }

    File file = new File(jarPath);

    return file;
  }

}