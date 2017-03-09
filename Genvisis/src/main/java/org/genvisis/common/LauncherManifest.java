package org.genvisis.common;

import java.io.File;
import java.io.IOException;
import java.util.Date;
import java.util.Iterator;
import java.util.jar.Attributes;
import java.util.jar.JarFile;
import java.util.jar.Manifest;

import org.genvisis.cnv.Launch;

import com.github.zafarkhaja.semver.Version;

/**
 * Describes the manifest of the jar containing the class that was used to launch this application.
 */
public class LauncherManifest {

	public static final String UNDETERMINED_VERSION = "0.0.0-unknown";
	public static final String BUILD_VERSION = "Implementation-Version";
	public static final String BUILD_TIMESTAMP = "Implementation-Date";
	public static final String GIT_COMMIT_SHA = "Implementation-Build";
	private static final Class<?> DEFAULT_LAUNCH_CLASS = Launch.class;

	private static Class<?> launchClass = null;

	public static Class<?> getLaunchClass() {
		if (launchClass == null) {
			setLaunchClass(DEFAULT_LAUNCH_CLASS);
		}
		return launchClass;
	}

	public static boolean setLaunchClass(Class<?> c) {
		if (launchClass == null) {
			return lockedSetLaunchClass(c);
		}
		return launchClass == c;
	}

	private static synchronized boolean lockedSetLaunchClass(Class<?> c) {
		if (launchClass == null) {
			launchClass = c;
		}
		return launchClass == c;
	}

	private Attributes attributes;
	private Version version;
	private String compileTime;
	private String buildType;
	private String builtBy;
	private String copyright;

	public LauncherManifest() {
		this(new Attributes());
	}

	public LauncherManifest(Attributes attributes) {
		super();
		this.attributes = attributes;
	}

	private void populate() {

		version = Version.valueOf(UNDETERMINED_VERSION);
		compileTime = "";
		buildType = "";
		builtBy = "";
		if (attributes != null) {
			// FIXME use attributes.getValue(String)
			Iterator<Object> it = attributes.keySet().iterator();
			while (it.hasNext()) {
				java.util.jar.Attributes.Name key = (java.util.jar.Attributes.Name) it.next();
				String keyword = key.toString();
				if (keyword.equals(BUILD_VERSION)) {
					version = Version.valueOf((String) attributes.get(key));
				}
				if (keyword.equals("Compile-Time")) {
					compileTime = (String) attributes.get(key);
				}
				if (keyword.equals("Build-Type")) {
					buildType = (String) attributes.get(key);
				}
				if (keyword.equals("Built-By")) {
					builtBy = (String) attributes.get(key);
				}
				if (keyword.equals("Copyright")) {
					copyright = (String) attributes.get(key);
				}
			}
		}
	}

	public String getBuiltBy() {
		return builtBy;
	}

	public String getCopyright() {
		return copyright;
	}

	public Attributes getAttributes() {
		return attributes;
	}

	public String getBuildType() {
		return buildType;
	}

	public Version getVersion() {
		return version;
	}

	public String getCompileTime() {
		return compileTime;
	}

	public static LauncherManifest loadGenvisisManifest() {
		File file = getCurrentFile();
		return loadManifest(file);
	}

	public static String getGenvisisInfo() {
		try {
			LauncherManifest manifest = LauncherManifest.loadGenvisisManifest();// until it always works
			return "Genvisis, " + manifest.getVersion().getNormalVersion() + "\n"
						 + manifest.getCopyright()
						 + "\n\n" + (new Date());
		} catch (Exception e) {
			return "Genvisis, v0.0.0\n(c)2009-2015 Nathan Pankratz, GNU General Public License, v2\n\n"
						 + (new Date());
		}
	}

	// https://ant.apache.org/manual/Tasks/manifest.html
	public static LauncherManifest loadManifest(File file) {
		JarFile jar = null;

		LauncherManifest currentManifest = new LauncherManifest();
		try {

			if (file != null && file.exists()) {
				jar = new JarFile(file);

				Manifest manifest = jar.getManifest();

				Attributes attributes = manifest.getMainAttributes();

				jar.close();
				currentManifest = new LauncherManifest(attributes);
			}
		} catch (IOException ioe) {

		}
		currentManifest.populate();
		return currentManifest;
	}

	private static File getCurrentFile() {
		File file = new File(getLaunchClass().getProtectionDomain().getCodeSource()
																				 .getLocation().getFile());// get

		if (!file.exists() || !file.getAbsolutePath().endsWith(".jar")) {
			file = new File("../" + PSF.Java.GENVISIS);
		}
		return file;
	}

	public static void main(String[] args) {

		loadGenvisisManifest();

	}
}
