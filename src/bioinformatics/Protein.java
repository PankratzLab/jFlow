package bioinformatics;

import java.io.*;

public class Protein {
	public static void parseAminoAcids(String filename) {
		BufferedReader reader;
        PrintWriter writer;
        String temp;
        int count;
        
		try {
	        reader = new BufferedReader(new FileReader(filename));
	        writer = new PrintWriter(new FileWriter(filename+"_aa.xln"));
	        count = 0;
	        while (reader.ready()) {
	        	temp = reader.readLine();
	        	if (temp.startsWith(">")) {
	        		System.out.println("Ignoring first line: "+temp);
	        	} else {
	        		for (int i = 0; i<temp.length(); i++) {
	        			count++;
	        			writer.println(count+"\t"+temp.charAt(i));
                    }
	        	}
	        }
	        reader.close();
            writer.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+filename+"\"");
	        System.exit(2);
        }
	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Reviews\\24 Rouleau random sequencing\\GCK5.txt";

	    String usage = "\n"+
	    "bioinformatics.Protein requires 0-1 arguments\n"+
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
		    parseAminoAcids(filename);
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
