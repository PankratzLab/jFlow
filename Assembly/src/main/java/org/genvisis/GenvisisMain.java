package org.genvisis;

import java.awt.GraphicsEnvironment;
import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.util.jar.JarFile;
import java.util.jar.Manifest;

import javax.swing.JOptionPane;

import org.genvisis.cnv.Launch;
import org.genvisis.common.LauncherManifest;

/**
 * This is an application entry point that currently simply delegates to {@link Launch}. It was
 * originally added because https://javafx-maven-plugin.github.io/ was unhappy about building a jar
 * with no source files.
 * <p>
 * NB: This entry point validates the Java version before delegating the launch process. It is
 * intended to allow some level of backwards-compatibility with at least Java 6, thus this class
 * must be insulated from all classes with a compilation requirement more recent than Java 6.
 * </p>
 */
public final class GenvisisMain {

	private GenvisisMain() {
		// prevent instantiation of class
	}

	public static void main(String... args) {
		if (!javaCheck()) {
			System.exit(0);
		}

		LauncherManifest.setLaunchClass(GenvisisMain.class);
		Launch.main(args);
	}

	/**
	 * Check the minimum required Java version (in the manifest) and compare it to the runtime
	 * version.
	 *
	 * @return true iff the runtime version is equal to or greater than the manifest version.
	 */
	private static boolean javaCheck() {
		Manifest manifest = null;
		try {
			File jarFile = getManifestFile();

			if (jarFile != null && jarFile.exists()) {
				JarFile jar = new JarFile(jarFile);

				manifest = jar.getManifest();

				jar.close();
			}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}

		if (manifest == null) {
			return false;
		}
		String buildVersion = removeBuildNumber(manifest.getMainAttributes().getValue("Build-Jdk"));
		String runVersion = System.getProperty("java.version");
		double build = getVersionNumber(buildVersion);
		double run = getVersionNumber(runVersion);
		boolean passed = Double.compare(run, build) >= 0;

		if (!passed) {
			String error = "Error: this application requires a Java version of at least " + buildVersion
										 + ".\nDetected version: " + runVersion
										 + ".\nPlease use a standalone Genvisis application, or update your java version: https://java.com/download";
			System.err.println(error);
			if (!GraphicsEnvironment.isHeadless()) {
				error = "<html>Error: this application requires a Java version of at least " + buildVersion
								+ ".<br />Detected version: " + runVersion
								+ ".<br />Please either:"
								+ "<ul><li>Use a standalone Genvisis application (http://genvisis.org)</li>"
								+ "<li>update your java version (https://java.com/download)</li></ul></html>";
				JOptionPane.showMessageDialog(null, error, "Outdated Java version",
																			JOptionPane.ERROR_MESSAGE);
			}
		}
		return passed;
	}

	/**
	 * Helper method to remove the "_nnn" build numbers
	 */
	private static String removeBuildNumber(String value) {
		// Remove build numbers
		return value.split("_")[0];
	}

	/**
	 * @param versionString A Java JRE version, in the format "1.V.0_mmm"
	 * @return A fractional representation of the version for easy comparison.
	 */
	private static double getVersionNumber(String versionString) {
		// Strip non-numerics out of the version string and prepend a "0."
		// (e.g. 1.8.0_111 > 0.180111)

		String v = "0." + versionString.replaceAll("[^\\d]", "");
		return Double.parseDouble(v);
	}

	/**
	 * @return File reference to the manifest for this jar, or null if no manifest file can be
	 *         identified.
	 */
	private static File getManifestFile() {
		String jarPath = GenvisisMain.class.getProtectionDomain().getCodeSource()
																			 .getLocation().getFile();
		try {
			jarPath = URLDecoder.decode(jarPath, "UTF-8");
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		}

		File jarFile = new File(jarPath);

		if (!jarFile.exists() || !jarFile.getAbsolutePath().endsWith(".jar")) {
			return null;
		}
		return jarFile;
	}
}
