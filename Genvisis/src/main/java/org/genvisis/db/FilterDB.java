package org.genvisis.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class FilterDB {
	public static void filter(String dbFilename, String filtersFilename, String outputFile,
														Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String[] header;
		new Hashtable<String, String>();
		new Vector<String>();
		int count;
		String result;
		Vector<Filter> filters;
		long time;

		time = new Date().getTime();
		filters = new Vector<Filter>();
		try {
			reader = new BufferedReader(new FileReader(filtersFilename));
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t", -1);
				if (!line[0].startsWith("#") && !line[0].startsWith("//")) {
					if (line.length != 2) {
						log.reportError("Malformed filter: " + Array.toStr(line, " "));
						log.reportError("     must have two values separated by a tab, where the first token is the filter and the second is the label to be used to describe it");
						reader.close();
						return;
					} else {
						try {
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
			new File(ext.parseDirectoryOfFile(outputFile)).mkdirs();
			writer = new PrintWriter(new FileWriter(outputFile));
		} catch (Exception e) {
			System.err.println("Error - could not write to file " + outputFile + "; is it locked?");
			return;
		}

		try {
			reader = new BufferedReader(new FileReader(dbFilename));
			header = reader.readLine().trim().split("[\\s]+");
			for (int i = 0; i < filters.size(); i++) {
				try {
					filters.elementAt(i).determineIndices(header);
				} catch (Elision e) {
					log.reportException(e);
					reader.close();
					writer.close();
					return;
				}
			}

			writer.println("Unit\tReasonFlagged\tnumFlags");
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				result = line[0] + "\t";
				count = 0;
				for (int i = 0; i < filters.size(); i++) {
					if (filters.elementAt(i).meetsCriteria(line, log)) {
						result += (count > 0 ? ";" : "") + filters.elementAt(i).getLabel(line);
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
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dbFilename + "\"");
			return;
		}
		log.report("Program finished in "	+ (new Date().getTime() - time) / 1000
								+ " sec. Marker filter results are now available at " + outputFile);
	}

	public static void fromParameters(String filename, Logger log) {
		Vector<String> params;
		String dbFilename, filtersFilename, output;
		String temp;

		params = Files.parseControlFile(filename, "FilterDB",
																		new String[] {"db=database.txt", "filters=filters.txt",
																									"out=output.txt",
																									"# the filters file will have two columns: allFilters  errorMessage. For example:",
																									"# CallRate>=0.95&CallRate<=0.99	Call freq 0.95 - 0.99",
																									"# pct_AA=1&CallRate<1	AA freq = 1 & call rate < 1",},
																		log);

		if (params != null) {
			dbFilename = null;
			filtersFilename = null;
			output = null;
			for (int i = 0; i < params.size(); i++) {
				temp = params.elementAt(i);

				if (temp.startsWith("db=")) {
					dbFilename = ext.parseStringArg(temp, null);
				} else if (temp.startsWith("filters=")) {
					filtersFilename = ext.parseStringArg(temp, null);
				} else if (temp.startsWith("out=")) {
					output = ext.parseStringArg(temp, null);
				}
			}

			if (!Files.exists(dbFilename)) {
				log.report("Error - could not find database file '" + dbFilename + "'; aborting");
				return;
			}

			if (!Files.exists(filtersFilename)) {
				log.report("Error - could not find filters file '" + filtersFilename + "'; aborting");
				return;
			}

			if (output == null) {
				output = ext.rootOf(dbFilename, false) + "_" + ext.rootOf(filtersFilename) + ".out";
			}

			filter(dbFilename, filtersFilename, output, log);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dbFilename = "db.txt";
		String filterFilename = "filters.txt";
		// String dbFilename = "D:/GEDI_exome/results/markerGenderChecks.xln";
		// String filterFilename = "D:/GEDI_exome/qc/default_exclusion.criteria";
		String outputFilename = null;
		String logFilename = null;

		String usage = "\n"	+ "db.FilterDB requires 0-1 arguments\n"
										+ "   (1) database filename (i.e. db=" + dbFilename + " (default))\n"
										+ "   (2) filters filename (i.e. filters=" + dbFilename + " (default))\n"
										+ "   (3) (optional) output filename (i.e. out=" + logFilename + " (default))\n"
										+ "   (4) (optional) log filename (i.e. log=" + logFilename + " (default))\n"
										+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("db=")) {
				dbFilename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("filters=")) {
				filterFilename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("out=")) {
				outputFilename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logFilename = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (logFilename == null) {
				logFilename = ext.parseDirectoryOfFile(filterFilename, false)	+ ext.rootOf(filterFilename)
											+ ".log";
			}
			if (outputFilename == null) {
				outputFilename = ext.parseDirectoryOfFile(filterFilename, false)
													+ ext.rootOf(filterFilename) + ".out";
			}
			filter(dbFilename, filterFilename, outputFilename, new Logger(logFilename));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
