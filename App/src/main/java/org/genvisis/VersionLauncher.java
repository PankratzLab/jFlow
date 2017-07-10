package org.genvisis;

/**
 * Entry point for native applications. Determines available genvisis versions, whether they are
 * available remotely or locally, and loads the last requested version on startup (or the most
 * recent available version if no version is requested).
 */
public class VersionLauncher {

	public static void main(String... args) {

		//FIXME create standalone project above Genvisis that just has serialized version class and standardizes reading/writing?
		//FIXME In genvisis on startup, if serialized verison file exists, read it and create update menu
		/*
		 * 1. Find all local (.genvisis) + remote (genvisis.org) versions of genvisis
		 * 2. If no local versions, download latest remote version
		 * 3. Create .genvisis/versions.ser which reports available versions, which are locally available, and which is last requested
		 * 4. Add selected genvisis.jar to classpath
		 * 5. Inspect jar to determine main class.
		 * 6. Call main class via reflection
		 */
	}
}
