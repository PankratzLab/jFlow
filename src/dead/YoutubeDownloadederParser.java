package dead;

import java.io.*;
import java.util.*;
import common.*;

public class YoutubeDownloadederParser {

	public static void parse(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		Vector<String> v = new Vector<String>();
		int count;
		long time;
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+".dat"));
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.contains("data-video-ids=")) {
					temp = temp.substring(temp.indexOf("data-video-ids=")+16);
					temp = temp.substring(0, temp.indexOf("\""));
					writer.println("http://www.youtube.com/watch?v="+temp);
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "YoutubeDownloadederParser.dat";

		String usage = "\n" + "dead.YoutubeDownloadederParser requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		filename = "C:\\Shows\\iPhone\\Caillou\\caillou.htm";
		try {
			parse(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
