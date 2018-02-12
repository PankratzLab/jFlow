package org.genvisis.common;

import com.github.zafarkhaja.semver.Version;

/**
 * Helper methods for the {@link Version} class
 */
public final class VersionHelper {

  private VersionHelper() {
    // Prevent instantiation of static utility class
  }

  /**
   * Helper method for finding previous semver versions. Accounts for SNAPSHOTs, such that the last
   * version of "2.0.1-SNAPSHOT" is "2.0.0". For non-SNAPSHOTs, the last release is just the given
   * version.
   *
   * @return The most recent release version based on the given version
   */
  public static Version lastRelease(Version v) {
    if ("SNAPSHOT".equals(v.getPreReleaseVersion())) {
      return decrement(v);
    }

    return v;
  }

  /**
   * @return Semver decrement of the given version.
   */
  public static Version decrement(Version v) {
    int major = v.getMajorVersion();
    int minor = v.getMinorVersion();
    int patch = v.getPatchVersion();

    if (patch > 0) {
      patch--;
    } else if (minor > 0) {
      minor--;
    } else if (major > 0) {
      major--;
    } else {
      throw new IllegalArgumentException("No previous version for: " + v.getNormalVersion());
    }

    return Version.forIntegers(major, minor, patch);
  }
}
