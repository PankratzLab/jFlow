package seq;

import java.io.*;
//import java.util.*;
import common.*;

public class FilterVariantLists {
	public static void filter(String dir, String chrs, String positions, boolean indelsOnly) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        int[][] locs;
        char[][] ops;
        String[] files;
        int loc;
        boolean passes;
        
        locs = new int[2][];
        ops = new char[2][];
        for (int i = 0; i<2; i++) {
            line = (i==0?chrs:positions).split("[\\s]+");
            locs[i] = new int[line.length];
            ops[i] = new char[line.length];
            for (int j = 0; j<line.length; j++) {
            	try {
            		locs[i][j] = Integer.parseInt(line[j].substring(1));
            		ops[i][j] = line[j].charAt(0);
            	} catch (Exception e) {
            		System.err.println("Error parsing location from "+line[i]);
            		e.printStackTrace();
            	}
            }
        }
        if (locs[0].length != locs[1].length) {
        	System.err.println("Error - unequal number of chrs and positions designations");
        	System.exit(1);
        }
        
        files = Files.list(dir, ".snplist", false);
        new File(dir+"filtered/").mkdirs();
        for (int i = 0; i<files.length; i++) {
    		try {
                reader = new BufferedReader(new FileReader(dir+files[i]));
                writer = new PrintWriter(new FileWriter(dir+"filtered/"+files[i]));
                while (reader.ready()) {
                	line = reader.readLine().trim().split("[\\s]+");
            		passes = true;
            		for (int j = 0; j<ops.length; j++) {
            			if (j==0) {
                			loc = Positions.chromosomeNumber(line[j].substring(line[j].startsWith("chr")?3:0));
            			} else {
            				loc = Integer.parseInt(line[j]);
            			}
            			for (int k = 0; k<ops[j].length; k++) {
                			switch (ops[j][k]) {
                            case '=':
                            	if (loc != locs[j][k]) {
                            		passes = false;
                            	}
    	                        break;
                            case '>':
                            	if (loc <= locs[j][k]) {
                            		passes = false;
                            	}
    	                        break;
                            case '<':
                            	if (loc >= locs[j][k]) {
                            		passes = false;
                            	}
    	                        break;
                            default:
                            	System.err.println("Error - unknown operator: "+ops[j][k]);
    	                        break;
                            }
                        }
                    }
            		if (indelsOnly && !line[2].equals("*")) {
            			passes = false;
            		}
            		if (passes) {
            			writer.println(Array.toStr(line));
            		}
                }
                reader.close();
                writer.close();
            } catch (FileNotFoundException fnfe) {
                System.err.println("Error: file \""+dir+files[i]+"\" not found in current directory");
                System.exit(1);
            } catch (IOException ioe) {
                System.err.println("Error reading file \""+dir+files[i]+"\"");
                System.exit(2);
            }
        }
	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\SequencingProjectWithCIDR\\SNPlists\\older_variant_files\\by_bed\\";
	    String chrs = "=2";
	    String positions = ">200000000";
	    boolean indelsOnly = true;
	    
	    String usage = "\n"+
	    "seq.FilterVariantLists requires 0-1 arguments\n"+
	    "   (1) filename (i.e. file="+dir+" (default))\n"+
	    "   (2) chrs (i.e. chrs="+chrs+" (default))\n"+
	    "   (3) positions (i.e. positions="+positions+" (default))\n"+
	    "   (4) keep indels only (i.e. --indeplsOnly (not the default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("file=")) {
			    dir = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("chr=")) {
			    chrs = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("positions=")) {
			    positions = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("positions=")) {
			    indelsOnly = ext.parseBooleanArg(args[i]);
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
		    filter(dir, chrs, positions, indelsOnly);
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
