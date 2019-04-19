package org.genvisis.meta;

import java.io.Serializable;
import java.util.regex.Matcher;

import com.github.zafarkhaja.semver.Version;

/**
 * {@link Serializable} container for information about a particular version of genvisis.
 * <p>
 * NB: Sorting order puts newest version first.
 * </p>
 */
public class GenvisisVersion implements Comparable<GenvisisVersion>, Serializable {

  private static final long serialVersionUID = 1L;

  public static final String VERSION_PREFIX = "Genvisis v";

  /**
   * NB: must be transient since it's not serializable
   */
  private transient Version semVer;

  private boolean available;
  private final String version;

  /**
   * Constructs a new {@link GenvisisVersion} with the given {@link #isAvailable()} status.
   *
   * @see #GenvisisVersion(String)
   */
  public GenvisisVersion(String version, boolean available) {
    this(version);
    this.available = available;
  }

  /**
   * @param version The version of genvisis. This can be either a direct {@link Version}-compatible
   *          numerical representation, formatted per {@link Info#VERSION_PATTERN}, or per
   *          {@link #toString()} (starting with {@link #VERSION_PREFIX}
   */
  public GenvisisVersion(String version) {
    Matcher matcher = Info.VERSION_PATTERN.matcher(version);
    if (matcher.find()) {
      this.version = matcher.group(1);
    } else if (version.startsWith(VERSION_PREFIX)) {
      this.version = version.replace(VERSION_PREFIX, "");
    } else if (Info.LATEST_VER.equals(version) || Version.valueOf(version) != null) {
      this.version = version;
    } else {
      throw new IllegalArgumentException("Not a valid version: " + version);
    }
    available = false;
  }

  /**
   * Copying constructor. Copies version only.
   */
  public GenvisisVersion(GenvisisVersion other) {
    this.version = other.version;
  }

  /**
   * @return True if we know the location of this {@link GenvisisVersion}
   */
  public boolean isAvailable() {
    return available;
  }

  /**
   * @return The {@link Version} for this instance.
   */
  public Version version() {
    if (semVer == null) {
      semVer = Version.valueOf(version);
    }
    return semVer;
  }

  /**
   * @return {@link #toString()} with {@link Info#MENU_DISABLED} appended as appropriate
   */
  public String menuString() {
    return toString() + (isAvailable() ? "" : Info.MENU_DISABLED);
  }

  @Override
  public int compareTo(GenvisisVersion o) {
    return o.version().compareTo(version());
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((version == null) ? 0 : version.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    GenvisisVersion other = (GenvisisVersion) obj;
    if (version == null) {
      if (other.version != null) return false;
    } else if (!version.equals(other.version)) return false;
    return true;
  }

  @Override
  public String toString() {
    // Just use the lastest version string if applicable.
    return Info.LATEST_VER.equals(version) ? Info.LATEST_VER : (VERSION_PREFIX + version);
  }

}
