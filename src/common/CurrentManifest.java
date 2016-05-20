package common;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.jar.Attributes;
import java.util.jar.JarFile;

/**
 * Try to parse the attributes of currently running jar;
 *
 */
public class CurrentManifest {
	private Attributes attributes;
	private String version;
	private String compileTime;
	private String buildType;
	private String builtBy;

	public CurrentManifest() {

	}

	public CurrentManifest(Attributes attributes) {
		super();
		this.attributes = attributes;
	}

	private void populate() {

		this.version = "";
		this.compileTime = "";
		this.buildType = "";
		this.builtBy = "";
		if (attributes != null) {
			Iterator<Object> it = attributes.keySet().iterator();
			while (it.hasNext()) {
				java.util.jar.Attributes.Name key = (java.util.jar.Attributes.Name) it.next();
				String keyword = key.toString();
				if (keyword.equals("Implementation-Version")) {
					this.version = (String) attributes.get(key);
				}
				if (keyword.equals("Compile-Time")) {
					this.compileTime = (String) attributes.get(key);
				}
				if (keyword.equals("Build-Type")) {
					this.buildType = (String) attributes.get(key);
				}
				if (keyword.equals("Build-Creator")) {
					this.builtBy = (String) attributes.get(key);
				}
			}
		}
	}

	public Attributes getAttributes() {
		return attributes;
	}

	public String getBuildType() {
		return buildType;
	}

	public String getVersion() {
		return version;
	}

	public String getCompileTime() {
		return compileTime;
	}

	// https://ant.apache.org/manual/Tasks/manifest.html
	public static CurrentManifest loadManifest() {
		File file = null;
		JarFile jar = null;

		CurrentManifest currentManifest = new CurrentManifest(new Attributes());
		try {
			try {
				file = new File(new CurrentManifest().getClass().getProtectionDomain().getCodeSource().getLocation().getFile());// get
				jar = new java.util.jar.JarFile(file);
			} catch (FileNotFoundException fnfe) {// if running in eclipse will get this
				file = null;
			}
			if (file != null) {
				try {
					file = new File("../" + PSF.Java.GENVISIS);
					jar = new java.util.jar.JarFile(file);
				} catch (FileNotFoundException fnfe2) {
					file = null;
				}
			}
			if (file == null) {
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

	public static void main(String[] args) {

		loadManifest();

	}
}
