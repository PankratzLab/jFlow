package org.genvisis.seq;

//import java.io.*;
//import java.util.*;
//import common.*;

public class CompareRelatives {
	public static void compare(String dir, String filename) {
//		BufferedReader reader;
//        PrintWriter writer;
//        String[] line;
//        String temp, trav;
//        Hashtable<String,String> hash = new Hashtable<String,String>();
//        Vector<String> v = new Vector<String>();
//        int count;
//        long time;
//        
        
	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    String dir = "";
	    String filename = "CompareRelatives.dat";

	    String usage = "\n"+
	    "seq.CompareRelatives requires 0-1 arguments\n"+
	    "   (1) filename (i.e. file="+filename+" (default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("file=")) {
			    filename = args[i].split("=")[1];
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
		    compare(dir, filename);
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
