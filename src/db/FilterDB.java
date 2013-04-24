package db;

import java.io.*;
import java.util.*;
import common.*;

public class FilterDB {
	public static void filter(String dbFilename, String filtersFilename, String outputFile, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String[] header;
		new Hashtable<String, String>();
		new Vector<String>();
		int count;
		String result;
//		Filter[] filters;
		Vector<Filter> filters;
		long time;

		time = new Date().getTime();
//		filters = new Filter[Files.countLines(filtersFilename, false)];
		filters = new Vector<Filter>(Files.countLines(filtersFilename, false));
		try {
			reader = new BufferedReader(new FileReader(filtersFilename));
//			for (int i = 0; i < filters.length; i++) {
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t", -1);
				if (!line[0].startsWith("#") && !line[0].startsWith("//")) {
					if (line.length != 2) {
						log.reportError("Malformed filter: "+Array.toStr(line, " "));
						log.reportError("     must have two values separated by a tab, where the first token is the filter and the second is the label to be used to describe it");
						reader.close();
						return;
					} else {
						try {
//							filters[i] = new Filter(line[0], line[1]);
							filters.add(new Filter(line[0], line[1]));
						} catch (Elision e) {
							log.reportException(e);
							reader.close();
							return;
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filtersFilename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filtersFilename + "\"");
			System.exit(2);
		}
		
		try {
			reader = new BufferedReader(new FileReader(dbFilename));
			header = reader.readLine().trim().split("[\\s]+");
//			for (int i = 0; i < filters.length; i++) {
			for (int i = 0; i < filters.size(); i++) {
				try {
//					filters[i].determineIndices(header);
					filters.elementAt(i).determineIndices(header);
				} catch (Elision e) {
					log.reportException(e);
					reader.close();
					return;
				}
			}

			writer = new PrintWriter(new FileWriter(outputFile));
			writer.println("Unit\tReason flagged\tnumFlags");
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				result = line[0] + "\t";
				count = 0;
//				for (int i = 0; i < filters.length; i++) {
				for (int i = 0; i < filters.size(); i++) {
//					if (filters[i].meetsCriteria(line, log)) {
					if (! filters.elementAt(i).meetsCriteria(line, log)) {
//						result += (count>0? ";" : "") + filters[i].getLabel();
						result += (count>0? ";" : "") + filters.elementAt(i).getLabel();
						count++;
					}
				}
				if (count > 0) {
					writer.println(result + "\t" + count);
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dbFilename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dbFilename + "\"");
			System.exit(2);
		}
		log.report("Program finished in " + (new Date().getTime() - time)/1000 + " sec. Marker filter results are now available at " + outputFile);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
//		String dbFilename = "db.txt";
//		String filterFilename = "filters.txt";
		String dbFilename = "D:/GEDI_exome/results/markerGenderChecks.xln";
		String filterFilename = "D:/GEDI_exome/qc/default_exclusion.criteria";
		String outputFilename = null;
		String logFilename = null;

		String usage = "\n" +
		"db.FilterDB requires 0-1 arguments\n" +
		"   (1) database filename (i.e. db=" + dbFilename+ " (default))\n" + 
		"   (2) filters filename (i.e. filters=" + dbFilename+ " (default))\n" + 
		"   (3) (optional) output filename (i.e. out=" + logFilename+ " (default))\n" + 
		"   (4) (optional) log filename (i.e. log=" + logFilename+ " (default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("db=")) {
				dbFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("filters=")) {
				filterFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outputFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logFilename = args[i].split("=")[1];
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
			if (logFilename == null) {
				logFilename = ext.parseDirectoryOfFile(filterFilename, false) + ext.rootOf(filterFilename) + ".log";
			}
			if (outputFilename == null) {
				outputFilename = ext.parseDirectoryOfFile(filterFilename, false) + ext.rootOf(filterFilename) + ".out";
			}
			filter(dbFilename, filterFilename, outputFilename, new Logger(logFilename));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
