package common;

import java.io.IOException;
import java.util.Iterator;

public class Manifest {
	
//	https://ant.apache.org/manual/Tasks/manifest.html
	public static void manifest() throws IOException {
		java.io.File file = new java.io.File("/Users/Kitty/Downloads/genvisis.jar");// get
																					// can
																					// path
		java.util.jar.JarFile jar = new java.util.jar.JarFile(file);
		java.util.jar.Manifest manifest = jar.getManifest();

		String versionNumber = "";
		java.util.jar.Attributes attributes = manifest.getMainAttributes();
		if (attributes != null) {
			Iterator<Object> it = attributes.keySet().iterator();
			while (it.hasNext()) {
				java.util.jar.Attributes.Name key = (java.util.jar.Attributes.Name) it.next();
				String keyword = key.toString();
				if (keyword.equals("Implementation-Version") || keyword.equals("Bundle-Version")) {
					versionNumber = (String) attributes.get(key);
					break;
				}
			}
		}
		jar.close();

		System.out.println("Version: " + versionNumber); // "Version: 1.3.162"
	}

	public static void main(String[] args) {
		try {
			manifest();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
