package org.genvisis.meta;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.github.zafarkhaja.semver.Version;

/**
 * Container for constants and information about the Genvisis environment.
 */
public class Info {

	public static final String JAR_PREFIX = "genvisis_v";

	/**
	 * Marker character for a menu entry that should be disabled
	 */
	public static final String MENU_DISABLED = "!";

	/**
	 * The prefix to a requested {@link Version} string when running from the CLI
	 */
	public static final String VERSION_CLI = "-ver=";

	/**
	 * Keyword to use the latest version
	 */
	public static final String LATEST_VER = "LATEST";

	/**
	 * Key for a boolean system environment variable to determine if we are running in a native
	 * (end-user app) environment.
	 */
	public static final String NATIVE_PROP_KEY = "Genvisis.isNative";

	/**
	 * Pattern to extract genvisis jar versions. Use {@link Matcher#group(int)} with {@code 1} for the
	 * version string.
	 */
	public static final Pattern VERSION_PATTERN = Pattern.compile(JAR_PREFIX
																																+ "([0-9]+[\\.[0-9]+]*)\\.jar");

	private static final String FILE_SEP = System.getProperty("file.separator");

	/**
	 * Path to directory with remote Genvisis versions
	 */
	public static final String REMOTE_URL = "http://www.genvisis.org/";

	/**
	 * Path to local Genvisis metadata storage
	 */
	public static final String GENVISIS_HOME = System.getProperty("user.home")
																						 + FILE_SEP + ".genvisis"
																						 + FILE_SEP;

	/**
	 * File on remote with list of release versions
	 */
	public static final String REMOTE_VERSION_FILE = REMOTE_URL + "releases.txt";

	/**
	 * Location (relative to {@link #GENVISIS_HOME} to store genvisis .jars
	 */
	public static final String LOCAL_VERSIONS_DIR = GENVISIS_HOME + ".versions" + FILE_SEP;

	/**
	 * Filename of the last requested version
	 */
	public static final String VERSION_TO_LOAD = GENVISIS_HOME + ".version.load.ser";

	/**
	 * List of all available versions
	 */
	public static final String VERSIONS_KNOWN_LIST = GENVISIS_HOME + ".version.avail.ser";
}
