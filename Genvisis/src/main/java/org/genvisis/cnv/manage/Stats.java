package org.genvisis.cnv.manage;

import org.genvisis.cnv.filesys.Project;

public class Stats {
	public static void runTheNumbers(Project proj) {
		System.out.println("Number of samples: "+proj.getSamples().length);
		System.out.println("Number of markers: "+proj.getMarkerNames().length);
		System.out.println("Number of markers in Lookup: "+proj.getMarkerLookup().getSize());
	}

	public static void main(String[] args) {
	    int numArgs = args.length;
	    String filename = null;

	    String usage = "\n"+
	    "cnv.manage.Stats requires 0-1 arguments\n"+
   		"   (1) project properties filename (i.e. proj="+org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"";
	    
	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
		try {
		    runTheNumbers(new Project(filename, false));
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
