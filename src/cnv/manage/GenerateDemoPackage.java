package cnv.manage;

import java.io.*;
//import java.util.*;

import cnv.filesys.Project;
import common.*;

public class GenerateDemoPackage {

	public static void makeUpScatter(Project proj, String filenameOfMarkerList, Logger log) {
		if (new File("demo/").exists()) {
			new File("demo/").renameTo(new File(Files.backup("demo", "", "")));
//			new File("projects/").renameTo(new File(Files.backup("demo/projects", "", "")));
			new File("projects/").delete();
		}
		new File("demo/data/").mkdirs();
		new File("demo/transposed/").mkdirs();
		new File("projects/").mkdirs();
		Files.copyFile(filenameOfMarkerList, "demo/data/test.txt");
		Files.copyFile(proj.getFilename(Project.SAMPLE_DATA_FILENAME), "demo/data/"+ext.removeDirectoryInfo(proj.getFilename(Project.SAMPLE_DATA_FILENAME)));
		
		
		
		
		
		proj.setJarStatus(true);
		proj.setProperty(Project.PROJECT_DIRECTORY, "demo/");
		proj.setProperty(Project.SOURCE_DIRECTORY, "demo/");
		proj.saveProperties("projects/"+ext.removeDirectoryInfo(proj.getPropertyFilename()));
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = cnv.Launch.getDefaultDebugProjectFile();
		String logfile = null;
		Logger log;
		Project proj;

		String usage = "\n" + "cnv.manage.GenerateDemoPackage requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(logfile);
			proj = new Project(filename, false);
			makeUpScatter(proj, proj.getFilename(Project.DISPLAY_MARKERS_FILENAME), log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
