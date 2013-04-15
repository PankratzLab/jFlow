package cnv;

import java.io.*;
import java.util.*;

import common.ext;

public class LaunchProperties extends Properties {
	public static final long serialVersionUID = 1L;

	public static final String DEFAULT_PROPERTIES_FILE = "launch.properties";
	public static final String PROJECTS_DIR = "PROJECTS_DIR";
	public static final String LAST_PROJECT_OPENED = "LAST_PROJECT_OPENED";
	
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
		
		dir = (String)get(PROJECTS_DIR);
		if (dir == null) {
			System.err.println("Error - '"+filename+"' did not contain a property called \""+PROJECTS_DIR+"\"");
			System.exit(1);
		}

		dir = ext.verifyDirFormat(dir);
			
		return dir;
	}
	
	public void save() {
		FileOutputStream out;
		
		try {
			out = new FileOutputStream(DEFAULT_PROPERTIES_FILE);
			store(out, null);
			out.close();		
        } catch (Exception e) {
        	System.err.println("Failed to save \""+DEFAULT_PROPERTIES_FILE+"\"");
        }
	}
	
}
