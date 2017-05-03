// currently PD dependent
package org.genvisis.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

public class heritabilityEstimate {
	public static final String DEFAULT_DB = tools.CRF_DIR + "crf_db.dat";
	public static final String DEFAULT_TRAIT = "Depression";
	public static final String DEFAULT_ROOT = "solar";

	public static void estimate(String database, String trait, String[] covariates, boolean extended,
															boolean affonly, String root, boolean run) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line, phenoNames = null, data;
		int[] indices = null;
		Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
		Hashtable<String, String> affHash;
		Vector<String> vips = new Vector<String>(), v;
		Vector<String[]> members;
		String trav, prev;
		affHash = tools.getBestPDdx();

		try {
			reader = new BufferedReader(new FileReader(database));

			phenoNames = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
			v = ArrayUtils.toStringVector(covariates);
			v.insertElementAt(trait, 0);
			indices = ext.indexFactors(ArrayUtils.toStringArray(v), phenoNames, true, true);

			while (reader.ready()) {
				line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
				if (!line[indices[0]].equals(".")) {
					data = new String[indices.length];
					for (int i = 0; i < data.length; i++) {
						data[i] = line[indices[i]];
					}
					if (!affonly || tools.isAffected(affHash, line[1] + "\t" + line[2])) {
						hash.put(line[1] + "\t" + line[2], data);
						vips.add(line[1] + "\t" + line[2]);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + database + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + database + "\"");
			System.exit(2);
		}

		reader = tools.getNinfoReader(2);
		members = new Vector<String[]>();
		reader.readLine();
		prev = "";
		while (reader.ready()) {
			line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
			trav = line[0] + "\t" + line[1];
			if (!trav.equals(prev)) {
				members.add(new String[] {line[0], line[1], (line[4].equals(".") ? "0" : line[4]),
																	(line[5].equals(".") ? "0" : line[5]),
																	(line[2].toUpperCase().equals("M") ? "1"
																																		 : line[2].toUpperCase()
																																							.equals("F") ? "2"
																																													 : "0")});
			}
			prev = trav;
		}

		System.err.println("Error - conols.prune has been disabled");

		// consol.prune.procPeeps(members, consol.prune.findKeepers(members,
		// vips));

		new File(root).mkdirs();
		writer = Files.openAppropriateWriter(root + "/" + root + ".fam");
		writer.println("FAMID,ID,FA,MO,SEX");
		for (int i = 0; i < members.size(); i++) {
			line = members.elementAt(i);
			writer.println(ArrayUtils.toStr(line, ","));
		}
		writer.close();

		writer = Files.openAppropriateWriter(root + "/" + root + ".ptypes");
		writer.print("FAMID,ID");
		for (int indice : indices) {
			writer.print("," + phenoNames[indice]);
		}
		writer.println();
		for (int i = 0; i < members.size(); i++) {
			line = members.elementAt(i);
			data = hash.containsKey(line[0] + "\t" + line[1]) ? hash.get(line[0] + "\t" + line[1])
																												: ArrayUtils.stringArray(indices.length);
			for (int j = 0; j < data.length; j++) {
				if (data[j].equals(".")) {
					data[j] = "";
				}
			}
			writer.println(line[0] + "," + line[1] + "," + ArrayUtils.toStr(data, ","));
		}
		writer.close();

		writer = Files.openAppropriateWriter(root + "/" + root + ".ptypes");
		writer.print("FAMID,ID");
		for (int indice : indices) {
			writer.print("," + phenoNames[indice]);
		}
		writer.println();
		for (int i = 0; i < members.size(); i++) {
			line = members.elementAt(i);
			data = hash.containsKey(line[0] + "\t" + line[1]) ? hash.get(line[0] + "\t" + line[1])
																												: ArrayUtils.stringArray(indices.length);
			for (int j = 0; j < data.length; j++) {
				if (data[j].equals(".")) {
					data[j] = "";
				}
			}
			writer.println(line[0] + "," + line[1] + "," + ArrayUtils.toStr(data, ","));
		}
		writer.close();

		writer = Files.openAppropriateWriter(root + "/batch");
		writer.println("echo -e \"load ped " + root + ".fam\\nautomodel " + root
									 + ".ptypes Depression\\npolygenic -screen\\nquit\\n\" | solar > " + root
									 + ".log");
		writer.close();
		if (!Files.isWindows()) {
			try {
				Runtime.getRuntime().exec("chmod +x " + root + "/batch").waitFor();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	public static void procFile(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String trait, trav;
		String[] line;

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = Files.openAppropriateWriter(filename + ".batch");
			while (reader.ready()) {
				line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
				if (!line[0].startsWith("#")) {
					if (ext.indexOfStartsWith("root=", line, false) >= 0) {
						trav = line[ext.indexOfStartsWith("root=", line, false)].split("=", -1)[1];
					} else {
						trav = DEFAULT_ROOT;
					}
					if (ext.indexOfStartsWith("trait=", line, false) >= 0) {
						trait = line[ext.indexOfStartsWith("trait=", line, false)].split("=", -1)[1];
					} else {
						trait = DEFAULT_TRAIT;
					}
					main(ArrayUtils.addStrToArray("batch=", line));
					writer.println("cd " + trav);
					writer.println("./batch > batch.log");
					writer.println("cp " + trait + "/polygenic.out ../" + trav + "_" + trait
												 + "_polygenic.out");
					writer.println("cd ..");
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error processing file \"" + filename + "\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String database = DEFAULT_DB;
		String trait = DEFAULT_TRAIT;
		boolean extended = true;
		// String covariates = "";
		// String covariates = "BlessedFunctionality;Education;UPDRSliving;Hoehn&Yahr";
		// String covariates =
		// "BlessedFunctionality;Education;UPDRSliving;Hoehn&Yahr;MMSE;DurationFromAOO;UPDRSmotor;PIGD_score";
		// String covariates = "BlessedFunctionality;Education;UPDRSliving;MMSE";
		String covariates = "UPDRSmotor;PIGD_score;BlessedFunctionality;Education;UPDRSliving;Hoehn&Yahr;MMSE;DurationFromAOO";
		String root = DEFAULT_ROOT;
		String batch = "batch.heritabilities";
		boolean run = true;
		boolean affonly = true;

		String usage = "\n" + "park.heritabilityEstimate requires 0-5 arguments\n"
									 + "   (1) database filename (i.e. db=" + database + " (default)\n"
									 + "   (2) trait name (i.e. trait=" + trait + " (default)\n"
									 + "   (3) names of covariates separated by a semicoln (i.e. covars=" + trait
									 + " (default)\n" + "   (4) use extended family members (i.e. extended="
									 + extended + " (default)\n" + "   (5) affecteds only (i.e. affonly=" + affonly
									 + " (default)\n" + "   (6) root of output filenames (i.e. root=" + root
									 + " (default)\n" + "   (7) run solar if path exists (i.e. run (default)\n"
									 + " OR\n" + "   (1) batch file with the options listed above (i.e. batch="
									 + batch + " (default, if it exists)\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.equals("")) {
				System.err.println("Error - batch file contains a blank line");
				return;
			} else if (arg.startsWith("db=")) {
				database = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("trait=")) {
				trait = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("covars=")) {
				covariates = arg.split("=", -1)[1];
				numArgs--;
			} else if (arg.startsWith("extended=")) {
				extended = arg.split("=")[1].toLowerCase().equals("true");
				numArgs--;
			} else if (arg.startsWith("affonly=")) {
				affonly = arg.split("=")[1].toLowerCase().equals("true");
				numArgs--;
			} else if (arg.startsWith("root=")) {
				root = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("batch=")) {
				batch = arg.split("=", -1)[1];
				numArgs--;
			} else if (arg.equals("run")) {
				run = true;
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		try {
			if (new File(batch).exists()) {
				procFile(batch);
			} else {
				estimate(database, trait, covariates.split(";"), extended, affonly, root, run);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
