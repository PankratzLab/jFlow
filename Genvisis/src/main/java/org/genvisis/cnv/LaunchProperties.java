package org.genvisis.cnv;

import java.io.*;
import java.util.*;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class LaunchProperties extends Properties {
	public static final long serialVersionUID = 1L;

	public static final String DEFAULT_PROPERTIES_FILE = "launch.properties";
	public static final String PROJECTS_DIR = "PROJECTS_DIR";
	public static final String LAST_PROJECT_OPENED = "LAST_PROJECT_OPENED";
	public static final String DEBUG_PROJECT_FILENAME = "DEBUG_PROJECT_FILENAME";
	
	public String filename;

	public LaunchProperties(String filename) {
		InputStream is;

		this.filename = filename;
		
		try {
			is = new FileInputStream(filename);
			load(is);
			is.close();
        } catch (Exception e) {
        	System.err.println("Failed to load \""+filename+"\"");
        }
	}

	public String getDirectory() {
		String dir;

		dir = (String) get(PROJECTS_DIR);
		if (dir == null || !Files.exists(dir)) {
			new Logger().reportTimeWarning("Did not detect directory with projects (or was missing property). Defaulting to projects/");
			// System.err.println("Error - '" + filename + "' did not contain a property called \"" + PROJECTS_DIR + "\"");
			dir = "projects/";
			new File(dir).mkdirs();
		}

		dir = ext.verifyDirFormat(dir);

		return dir;
	}

	public String getFilename() {
		return filename;
	}

	public static String directoryOfLaunchProperties(String launchPropertiesFile) {
		String path = null;
		try {
			path = ext.parseDirectoryOfFile(new File(launchPropertiesFile).getCanonicalPath());
		} catch (IOException ioe) {
			path = "";
		}
		return path;
	}
	
	public void save() {
		FileOutputStream out;
		
		try {
			out = new FileOutputStream(filename);
			store(out, null);
			out.close();		
        } catch (Exception e) {
        	System.err.println("Failed to save \""+filename+"\"");
        }
	}
	
}
