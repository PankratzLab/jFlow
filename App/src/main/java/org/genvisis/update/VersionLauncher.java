package org.genvisis.update;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.URI;
import java.net.URL;
import java.net.URLClassLoader;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.jar.JarFile;
import java.util.jar.Manifest;
import java.util.stream.Collectors;

import javax.swing.JOptionPane;
import javax.swing.ProgressMonitorInputStream;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;

import org.genvisis.meta.GenvisisVersion;
import org.genvisis.meta.Info;
import org.genvisis.meta.Persistence;

import com.github.zafarkhaja.semver.UnexpectedCharacterException;
import com.github.zafarkhaja.semver.Version;

/**
 * Entry point for native applications. Determines available genvisis versions, whether they are
 * available remotely or locally, and loads the last requested version on startup (or the most
 * recent available version if no version is requested).
 */
public final class VersionLauncher {

	/**
	 * Main class attribute in manifest
	 */
	private static final String MANIFEST_MAIN = "Main-Class";

	private static GenvisisVersion loaded = null;
	private static Map<GenvisisVersion, String> pathMap = new HashMap<>();

	private VersionLauncher() {
		// Disable construction of utility class
	}

	/**
	 * Identify the available app versions and launch either the last requested, or the most recent.
	 * 
	 * @param args
	 */
	private static void launchApp(String[] args) {
		System.setProperty(Info.NATIVE_PROP_KEY, "true");
		GenvisisVersion toLaunch = findLaunchVersion();

		if (toLaunch == null) {
			System.err.println("No Genvisis versions available.");
			return;
		}

		String jarPath = pathMap.get(toLaunch);

		// Check if the jar exists locally
		if (isRemote(jarPath)) {
			// If not, download from remote
			try {
				downloadAndRun(jarPath, toLaunch, args);
			} catch (IOException e) {
				System.err.println("Unable to download Genvisis from location: " + jarPath);
				return;
			}
		} else {
			launchVersion(jarPath, args);
		}
	}

	/**
	 * Launches the application at the specified version.
	 * <p>
	 * NB: this method should be called from the EDT
	 * </p>
	 *
	 * @param toLaunch Requested version to launch
	 * @param args
	 */
	private static void launchVersion(String localPath, String[] args) {
		// Since we are now running on the EDT, it is now safe to stop any Persistence requests
		Persistence.stop(VersionLauncher.class);
		try (JarFile jar = new JarFile(localPath)) {
			// Inspect jar manifest to find the main class name
			Manifest m = jar.getManifest();
			String mainClassName = m.getMainAttributes().getValue(MANIFEST_MAIN);

			// Add the target jar to the classpath
			// From: https://stackoverflow.com/a/60775/1027800
			URLClassLoader child = new URLClassLoader(new URL[] {new URL("file://" + localPath)},
																								VersionLauncher.class.getClassLoader());

			// Get a handle on the main class
			Class<?> mainClass = Class.forName(mainClassName, true, child);

			// Invoke main method of main class
			Method mainMethod = mainClass.getMethod("main", String[].class);
			SwingUtilities.invokeLater(() -> invokeSafely(mainMethod, args));
		} catch (FileNotFoundException e) {
			System.err.println("Genvisis application files not available.");
		} catch (NoSuchMethodException | SecurityException | IllegalArgumentException
						 | ClassNotFoundException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Run a main method and handle exceptions if necessary
	 */
	private static void invokeSafely(Method mainMethod, String[] args) {
		try {
			mainMethod.invoke(null, (Object) args);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * From https://stackoverflow.com/a/14413945/1027800
	 * <p>
	 * Just a simple method to download a genvisis.jar
	 * </p>
	 */
	private static void downloadAndRun(final String url,
																		 final GenvisisVersion toLaunch,
																		 String[] args) throws IOException {
		String localPath = Info.LOCAL_VERSIONS_DIR + Info.JAR_PREFIX
											 + toLaunch.version() + ".jar";

		// The app hasn't actually started yet, and can't start until the download is complete (since
		// we're downloading the jar we want to run).
		// SwingWorker.execute runs on a Daemon thread, so we need to do something to preserve the JVM
		Persistence.start(VersionLauncher.class);

		final URLConnection urlConn = new URL(url).openConnection();
		SwingWorker<Object, Integer> w = new SwingWorker<Object, Integer>() {
			@Override
			protected Object doInBackground() {
				try (FileOutputStream fos = new FileOutputStream(localPath);
						 ProgressMonitorInputStream monitorIS = new ProgressMonitorInputStream(null,
																																									 "Downloading "
																																												 + toLaunch.toString(),
																																									 new BufferedInputStream(urlConn.getInputStream()));) {
					byte[] buffer = new byte[4096]; // declare 4KB buffer
					int len;
					monitorIS.getProgressMonitor().setMillisToPopup(10);
					monitorIS.getProgressMonitor().setMaximum(urlConn.getContentLength());

					// while we have availble data, continue downloading and storing to local file
					while ((len = monitorIS.read(buffer)) > 0) {
						fos.write(buffer, 0, len);
					}
				} catch (IOException e) {
					File incompleteOutput = new File(localPath);
					if (incompleteOutput.exists()) {
						incompleteOutput.delete();
					}
				}
				return null;
			}

			@Override
			protected void done() {
				// NB: automatically queued on the EDT
				launchVersion(localPath, args);
			}
		};

		w.execute();
	}


	/**
	 * Check if there is a previously-requested app version file. If so, read and return the version
	 * from it. If not, create the file with the given fallback version.
	 *
	 * @param newestAvailable Default version to use if no previously-requested version available
	 * @return The version to load
	 */
	private static GenvisisVersion getRequestedVersion(GenvisisVersion newestAvailable,
																										 GenvisisVersion newestKnown) {
		// Check if cached version file exists
		File toLoad = new File(Info.VERSION_TO_LOAD);
		if (toLoad.exists()) {
			// If so, read contents
			GenvisisVersion cached = getLoadedVersion();

			// If set to use latest version, just return it immediately
			if (Info.LATEST_VER.equals(cached.toString())) {
				return newestAvailable;
			}

			// Check if the cached version is invalid, or if the user wants to update
			if (pathMap.containsKey(cached) && new File(pathMap.get(cached)).exists()
					&& !updateIfNew(newestAvailable, newestKnown)) {
				// Cached version was valid and user doesn't want to update
				return cached;
			}
		}

		// If the cached version did not exist or was invalid, use the fallback version
		write(Info.VERSION_TO_LOAD, newestAvailable);
		return newestAvailable;
	}

	/**
	 * @return The version to launch
	 */
	private static GenvisisVersion findLaunchVersion() {
		Set<GenvisisVersion> availableVersions = new HashSet<>();

		// Add the list of all local (cached in .genvisis) versions
		availableVersions.addAll(getLocalVersions());

		// Add the all remote versions (on genvisis.org) that are not available locally
		availableVersions.addAll(getRemoteVersions());

		Optional<GenvisisVersion> first = availableVersions.stream().sorted()
																											 .filter(GenvisisVersion::isAvailable)
																											 .findFirst();

		if (!first.isPresent()) {
			// No available versions
			return null;
		}

		GenvisisVersion firstAvailable = first.get();
		GenvisisVersion newestKnown = null;

		// Read the current versions known list, if it exists
		File versionsKnown = new File(Info.VERSIONS_KNOWN_LIST);
		if (versionsKnown.exists()) {
			try (BufferedInputStream fileIn = new BufferedInputStream(new FileInputStream(Info.VERSIONS_KNOWN_LIST));
					 ObjectInputStream in = new ObjectInputStream(fileIn)) {
				List<GenvisisVersion> knownVersions = (ArrayList<GenvisisVersion>) in.readObject();
				// Pop off the "Latest" version
				knownVersions.remove(0);

				// Cache the newest known version for update testing later
				newestKnown = knownVersions.get(0);
				// Combine the known + available versions list, stripping the "isAvailable" status from any
				// versions that are only "known"
				knownVersions.stream().filter(v -> !availableVersions.contains(v))
										 .forEach(v -> availableVersions.add(new GenvisisVersion(v)));
			} catch (IOException e) {
				e.printStackTrace();
			} catch (ClassNotFoundException e) {
				e.printStackTrace();
			}
		}

		// Sort the list of GenvisisVersions
		ArrayList<GenvisisVersion> allVersions = availableVersions.stream().sorted()
																															.collect(Collectors.toCollection(ArrayList::new));

		// Add a "latest version" entry, which is always available and allows the user to auto-update.
		allVersions.add(0, new GenvisisVersion(Info.LATEST_VER, true));

		// (over)write this list to the available versions file
		write(Info.VERSIONS_KNOWN_LIST, allVersions);

		return getRequestedVersion(firstAvailable, newestKnown);
	}

	/**
	 * Test if an update is available. If so, ask user if they want to update
	 */
	private static boolean updateIfNew(GenvisisVersion newestAvailable,
																		 GenvisisVersion newestKnown) {
		if (newestAvailable == null || newestKnown == null) {
			return false;
		}

		Version newVer = newestAvailable.version();
		Version prev = newestKnown.version();

		// If the newest available version is newer than the newest known, ask if user wants to update
		return newVer.greaterThan(prev)
					 && (JOptionPane.showConfirmDialog(null,
																						 "<html>New version: " + newVer.toString()
																									 + "<br />Would you like to update?</html>",
																						 "Update available",
																						 JOptionPane.YES_NO_OPTION) == JOptionPane.YES_OPTION);
	}

	/**
	 * @return List of all {@link GenvisisVersion}s available locally.
	 */
	private static Set<GenvisisVersion> getLocalVersions() {
		Set<GenvisisVersion> versions = new HashSet<>();
		File versionDir = new File(Info.LOCAL_VERSIONS_DIR);

		if (!versionDir.exists()) {
			versionDir.mkdir();
		}

		for (File jar : versionDir.listFiles()) {
			try {
				GenvisisVersion v = new GenvisisVersion(jar.getName(), true);
				versions.add(v);
				pathMap.put(v, jar.getAbsolutePath());
			} catch (UnexpectedCharacterException | IllegalArgumentException e) {
				// Ignore invalid files in this dir
			}
		}

		return versions;
	}

	/**
	 * @return List of all {@link GenvisisVersion}s available remotely.
	 */
	private static Set<GenvisisVersion> getRemoteVersions() {
		Set<GenvisisVersion> versions = new HashSet<>();

		// Parse remote version file
		try (BufferedReader reader = new BufferedReader(new InputStreamReader(new URL(Info.REMOTE_VERSION_FILE).openStream()));) {

			while (reader.ready()) {
				String version = reader.readLine();
				GenvisisVersion v = new GenvisisVersion(version, true);
				versions.add(v);
				// local versions take precedence
				if (!pathMap.containsKey(v)) {
					pathMap.put(v,
											Info.REMOTE_URL + Info.JAR_PREFIX + version + ".jar");
				}
			}
		} catch (IOException e) {
			System.err.println("Remote versions unavailable.");
		}

		return versions;
	}

	/**
	 * @return The active Genvisis version
	 */
	private static GenvisisVersion getLoadedVersion() {
		if (loaded == null) {
			readLoadedVersion();
		}
		return loaded;
	}

	/**
	 * Deserialize and cache the {@link #VERSION_TO_LOAD} file.
	 */
	private static synchronized void readLoadedVersion() {
		if (loaded != null) {
			return;
		}

		try (BufferedInputStream fileIn = new BufferedInputStream(new FileInputStream(Info.VERSION_TO_LOAD));
				 ObjectInputStream in = new ObjectInputStream(fileIn)) {
			loaded = (GenvisisVersion) in.readObject();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
	}

	/**
	 * If the first argument matches {@link Info#VERSION_CLI}, we write the value of that arg to the
	 * {@link Info#VERSION_TO_LOAD} file, setting it as the requested version.
	 *
	 * @param args
	 * @return
	 */
	private static String[] parseVersion(String[] args) {
		String[] rVal = args;
		if (rVal.length > 0 && rVal[0].startsWith(Info.VERSION_CLI)) {
			String v = args[0].replace(Info.VERSION_CLI, "");
			// Since we matched the version_cli flag, pop off the first arg
			rVal = Arrays.copyOfRange(rVal, 1, rVal.length);
			if (Info.LATEST_VER.equalsIgnoreCase(v)) {
				useLatest();
			} else {
				try {
					// If we have a valid semver version, write it to the version_to_load file
					GenvisisVersion toLaunch = new GenvisisVersion(v);
					write(Info.VERSION_TO_LOAD, toLaunch);
				} catch (Exception e) {
					System.out.println("Not a valid Genvisis version: " + v);
				}
			}
		}

		return rVal;
	}

	/**
	 * Clears any cached version requests, causing the most recent version to be used.
	 */
	private static void useLatest() {
		// Delete the version_to_load file and let it be auto-populated to the latest version
		File loadFile = new File(Info.VERSION_TO_LOAD);
		if (loadFile.exists()) {
			loadFile.delete();
		}
	}

	/**
	 * Serialize an object
	 */
	private static void write(String path, Object obj) {
		File f = new File(path);
		if (f.exists()) {
			f.delete();
		}
		try (BufferedOutputStream fileOut = new BufferedOutputStream(new FileOutputStream(path));
				 ObjectOutputStream out = new ObjectOutputStream(fileOut)) {
			out.writeObject(obj);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param path
	 * @return True if the given path is remote
	 */
	private static boolean isRemote(String path) {
		// copied from https://stackoverflow.com/a/25238203/1027800
		try {
			URI uri = new URI(path);
			return uri.getScheme().startsWith("http");
		} catch (Exception e) {
			return false;
		}
	}

	/**
	 * Entry point
	 *
	 * @param args
	 */
	public static void main(String... args) {
		String[] parsedArgs = parseVersion(args);
		launchApp(parsedArgs);
	}
}
