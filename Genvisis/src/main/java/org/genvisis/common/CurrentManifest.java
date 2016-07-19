package org.genvisis.common;

import java.io.File;
import java.io.IOException;
import java.util.Date;
import java.util.Iterator;
import java.util.jar.Attributes;
import java.util.jar.JarFile;

import org.genvisis.common.HttpUpdate.Version;

/**
 * Try to parse the attributes of currently running jar;
 *
 */
public class CurrentManifest {
	public static final String IMPLEMENTATION_VERSION = "Implementation-Version";
	private Attributes attributes;
	private Version version;
	private String compileTime;
	private String buildType;
	private String builtBy;
	private String copyright;

	public CurrentManifest() {

	}

	public CurrentManifest(Attributes attributes) {
		super();
		this.attributes = attributes;
	}

	private void populate() {

		this.version = new Version("v-1.-1.-1");
		this.compileTime = "";
		this.buildType = "";
		this.builtBy = "";
		if (attributes != null) {
			Iterator<Object> it = attributes.keySet().iterator();
			while (it.hasNext()) {
				java.util.jar.Attributes.Name key = (java.util.jar.Attributes.Name) it.next();
				String keyword = key.toString();
				if (keyword.equals(IMPLEMENTATION_VERSION)) {
					this.version = new Version((String) attributes.get(key));
				}
				if (keyword.equals("Compile-Time")) {
					this.compileTime = (String) attributes.get(key);
				}
				if (keyword.equals("Build-Type")) {
					this.buildType = (String) attributes.get(key);
				}
				if (keyword.equals("Built-By")) {
					this.builtBy = (String) attributes.get(key);
				}
				if (keyword.equals("Copyright")) {
					this.copyright = (String) attributes.get(key);
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

	public static CurrentManifest loadGenvisisManifest() {
		File file = getCurrentFile();
		return loadManifest(file);
	}

	public static String getGenvisisInfo() {
		try {
			CurrentManifest manifest = CurrentManifest.loadGenvisisManifest();// until it always works
			return "Genvisis, " + manifest.getVersion().getVersion() + "\n" + manifest.getCopyright() + "\n\n" + (new Date());
		} catch (Exception e) {
			return "Genvisis, v0.0.0\n(c)2009-2015 Nathan Pankratz, GNU General Public License, v2\n\n" + (new Date());
		}
	}

	// https://ant.apache.org/manual/Tasks/manifest.html
	public static CurrentManifest loadManifest(File file) {
		JarFile jar = null;

		CurrentManifest currentManifest = new CurrentManifest(new Attributes());
		try {

			if (file != null && file.exists()) {
				jar = new java.util.jar.JarFile(file);

				java.util.jar.Manifest manifest = jar.getManifest();

				Attributes attributes = manifest.getMainAttributes();

				jar.close();
				currentManifest = new CurrentManifest(attributes);
			}
		} catch (IOException ioe) {

		}
		currentManifest.populate();
		return currentManifest;
	}

	private static File getCurrentFile() {
		File file = new File(new CurrentManifest().getClass().getProtectionDomain().getCodeSource().getLocation().getFile());// get

		if (!file.exists() || !file.getAbsolutePath().endsWith(".jar")) {
			file = new File("../" + PSF.Java.GENVISIS);
		}
		return file;
	}

	public static void main(String[] args) {

		loadGenvisisManifest();

	}
}
