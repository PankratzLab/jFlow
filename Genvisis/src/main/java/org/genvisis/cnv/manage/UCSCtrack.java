package org.genvisis.cnv.manage;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;

public class UCSCtrack {
	public static final String SAMPLE_DEMOGRAPHICS = "SampleData.txt";
	public static final String DEFAULT_SAMPLE_DEMOGRAPHIC_DIRECTORY = "data/";
	public static final String[][] NEEDS = new String[][] {{"CLASS=Original Phenotype",
																													"CLASS=Final Phenotype",
																													"CLASS=Phenotype", "CLASS=Affection",
																													"CLASS=Affected"}, {"COVAR=Age"}};

	public static void makeTrack(String filename, String outfile, Logger log) {
		BufferedReader reader;
		PrintWriter track;
		String[] line, header;
		String[] demo;
		Hashtable<String, String[]> hash;
		int[] indices;
		int idIndex;
		long time;

		hash = new Hashtable<String, String[]>();
		try {
			reader = Files.getReader(	SAMPLE_DEMOGRAPHICS,
																new String[] {"", DEFAULT_SAMPLE_DEMOGRAPHIC_DIRECTORY,
																							"../" + DEFAULT_SAMPLE_DEMOGRAPHIC_DIRECTORY});
			header = reader.readLine().trim().split("\t", -1);
			indices = ext.indexFactors(NEEDS, header, false, true, true, log, false);
			idIndex = ext.indexOfStr("IID", header);
			if (idIndex == -1) {
				throw new IOException("Column header 'IID' is required in sample database file in order to lookup affection status");
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				hash.put(line[idIndex], ArrayUtils.subArray(line, indices, "."));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError(fnfe.getMessage());
		} catch (IOException ioe) {
			log.reportError("Error reading Sample database file");
			log.reportException(ioe);
		}

		time = new Date().getTime();
		System.out.println("Generating " + outfile);
		try {
			reader = new BufferedReader(new FileReader(filename));
			if (!ext.checkHeader(	reader.readLine().trim().split("[\\s]+"), CNVariant.PLINK_CNV_HEADER,
														false)) {
				reader.close();
				return;
			}
			track = Files.getAppropriateWriter(filename.substring(0, filename.lastIndexOf("."))
																					+ ".bed.gz");
			track.println("track name=\""	+ ext.rootOf(filename)
										+ "\" description=\"CNV data\" visibility=2 itemRgb=\"On\"");
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				demo = hash.containsKey(line[1]) ? hash.get(line[1]) : ArrayUtils.stringArray(NEEDS.length, ".");
				track.print("chr"	+ line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[1]
										+ (demo[0].equals(".") || demo[1].equals(".")	? ""
																																	: (Integer.parseInt(demo[0]) == 2	? ";AOO="
																																																			+ demo[1]
																																																		: ";AAE="
																																																			+ demo[1]))
										+ "\t50\t.\t0\t0\t");
				// track.print("chr"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[1]+"\t50\t.\t0\t0\t");
				switch (demo[0].equals(".") ? 0 : Integer.parseInt(demo[0])) {
					// switch (2) {
					case 2:
						if (Integer.parseInt(line[5]) < 2) {
							track.println("255,0,0"); // deletions in red
						} else {
							track.println("255,100,100"); // duplications in pink
						}
						break;
					case 1:
						if (Integer.parseInt(line[5]) < 2) {
							track.println("50,100,50"); // deletions in dark green
						} else {
							track.println("100,200,50"); // duplications in light green
						}
						break;
					default:
						if (Integer.parseInt(line[5]) < 2) {
							track.println("0,0,0"); // deletions in black
						} else {
							track.println("125,125,125"); // duplications in grey
						}
						break;
				}
			}
			reader.close();
			track.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
		} catch (IOException ioe) {
			log.reportError("Error creating UCSC bed track from " + filename);
			log.reportException(ioe);
		}
		System.out.println("Finished in " + ext.getTimeElapsed(time));

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "penncnv.cnv";
		String outfile = null;

		String usage = "\\n"	+ "cnv.manage.UCSCtrack requires 0-1 arguments\n"
										+ "   (1) CNV filename to convert (i.e. file=" + filename + " (default))\n"
										+ "   (2) name of output file (i.e. out=[file].bed.gz (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("out=")) {
				outfile = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - don't know what to do with argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (outfile == null) {
				outfile = ext.rootOf(filename, false) + ".bed.gz";
			}
			makeTrack(filename, outfile, new Logger(ext.rootOf(filename, false) + ".log"));
		} catch (Exception e) {
			e.printStackTrace();
		}

		ext.waitForResponse();
	}

}
