package cnv.manage;

import java.io.*;
import java.util.*;

import cnv.filesys.*;
import common.*;

public class DumpSerializedObject {
	private static void dump(String filename, Logger log) {
		Object object;
		
		object = Files.readSerial(filename, false, log, false, false);
		log.report("Information on class:"+
				"\n"+"object.getClass().getName()="+object.getClass().getName()+
				"\n"+"object.getClass()="+object.getClass()+
				"\n"+"object.toString()="+object.toString());
		/**
		 * Classes we should support now:
		 * 
		 * outliers.ser
		 * AnnotationCollection
		 * Centroids
		 * ClusterFilterCollection
		 * MarkerSet
		 * CNVariant
		 * SNPMarkerSet
		 * 
		 */

		/**
		 * Classes we might consider supporting in the future:
		 * 
		 * BurdenMatrix
		 * GenotypeMatrix
		 * LDdatabase
		 * SerialStringMatrix
		 * ChromatinAccessibility
		 * SuperNovo
		 */

		if (object instanceof AnnotationCollection) {
			log.report("Detected an AnnotationCollection file");
//			AnnotationCollection.dump();
		}
		
		
    	if (filename.endsWith("outliers.ser")) { // Hashtable<String, Float>
    		Sample.dumpOutOfRangeValues(filename, ext.parseDirectoryOfFile(filename) + ext.rootOf(filename) + "_dump.xln", false);
    	}
		
		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String logfile = null;
		Logger log;

		String usage = "\n" + 
				"widgets.DumpSerializedObject requires 0-1 arguments\n" + 
				"   (1) filename (e.g. file=outliers.ser (default))\n" + "";

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
			if (logfile == null) {
				logfile = ext.rootOf(filename)+"_dump.log";
			}
			log = new Logger(logfile);
			dump(filename, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
