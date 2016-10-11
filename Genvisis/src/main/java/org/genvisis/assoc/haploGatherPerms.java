package org.genvisis.assoc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.ext;
import org.genvisis.stats.ProbDist;

public class haploGatherPerms {
	public haploGatherPerms(int numReps) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null, addons;
		String[] line;
		String temp;
		int count, df_probs;
		double[] globals = new double[numReps];
		int df = -1;
		int total, lt05, lt01, lt001, lt0001, lt00001;
		double average, sig;

		count = 0;
		df_probs = 0;
		addons = new PrintWriter(new FileWriter("gatherAdditionalPerms"));
		writer = new PrintWriter(new FileWriter("globals.xls"));

		for (int i = 1; i <= numReps; i++) {
			if (new File("bintests/" + i + "_bintest.out").exists()) {
				reader = new BufferedReader(new FileReader("bintests/" + i + "_bintest.out"));
				for (int j = 0; j < 6; j++) {
					reader.readLine();
				}
				line = reader.readLine().split("[\\s]+");
				if (line.length < 6) {
					System.err.println("Error - could not parse " + "bintests/" + i + "_bintest.out");
					globals[i - 1] = -1;
				} else {
					globals[i - 1] = Double.valueOf(line[2].substring(0, line[2].length() - 1)).doubleValue();
					if (globals[i - 1] < -200) {
						System.out.println("rep " + i);
					}
					writer.println(globals[i - 1]);
					if (df == -1) {
						df = Integer.valueOf(line[5].substring(0, line[5].length() - 1)).intValue();
					} else if (Integer.valueOf(line[5].substring(0, line[5].length() - 1)).intValue() < df) {
						System.err.print(".");
						globals[i - 1] = -1;
						specialPerm(i);
						addons.println("Splus7.0 BATCH sp_perm." + i + " sp_check." + i + ".out");
						addons.println("sleep 5");
						df_probs++;
					} else if (Integer.valueOf(line[5].substring(0, line[5].length() - 1)).intValue() > df) {
						df = Integer.valueOf(line[5].substring(0, line[5].length() - 1)).intValue();
						System.err.println("Uh oh - ned to redo reps 1-"	+ i
																+ ", cause they have less df than subsequent runs.");
					}
				}
				reader.close();
			} else {
				addons.println("Splus7.0 BATCH perm." + i + " check." + i + ".out");
				addons.println("sleep 5");
				count++;
				System.out.print(".");
			}
		}
		addons.close();
		writer.close();
		System.out.println();
		System.err.println();
		System.err.println("Warning - there are "	+ count
												+ " missing replicates. Run gatherAdditionalPerms to rerun just these replicates.");
		System.err.println("Warning - there are "	+ df_probs
												+ " replicates, that dropped haplotypes. Run gatherAdditionalPerms to rerun just these replicates.");

		total = lt05 = lt01 = lt001 = lt0001 = lt00001 = 0;
		average = 0;
		for (int i = 0; i < numReps; i++) {
			if (globals[i] >= 0) {
				average += globals[i];
				sig = ProbDist.ChiDist(globals[i], df);
				if (sig < 0.05) {
					lt05++;
				}
				if (sig < 0.01) {
					lt01++;
				}
				if (sig < 0.001) {
					lt001++;
				}
				if (sig < 0.0001) {
					lt0001++;
				}
				if (sig < 0.00001) {
					lt00001++;
				}
				total++;
			}
		}

		average = average / total;
		temp = "Mean stat: "	+ ext.formDeci(average, 3, true) + " with df=" + df + " (p="
						+ ProbDist.ChiDist(average, df) + ")\n";
		temp += "Alpha for p=0.05 is " + ext.formDeci(lt05 / (double) total, 5, true) + "\n";
		temp += "Alpha for p=0.01 is " + ext.formDeci(lt01 / (double) total, 5, true) + "\n";
		temp += "Alpha for p=0.001 is " + ext.formDeci(lt001 / (double) total, 5, true) + "\n";
		temp += "Alpha for p=0.0001 is " + ext.formDeci(lt0001 / (double) total, 5, true) + "\n";
		temp += "Alpha for p=0.00001 is " + ext.formDeci(lt00001 / (double) total, 5, true) + "\n";
		writer = new PrintWriter(new FileWriter("haploPerm_summary.out"));
		System.out.println(temp);
		writer.println(temp);
		writer.close();

	}

	public void specialPerm(int i) {
		PrintWriter writer = null;

		try {
			writer = new PrintWriter(new FileWriter("sp_perm." + i));
		} catch (IOException ioe) {
		}

		writer.println("library(haplo.stats)");
		writer.println("geno<-matrix(scan(\"pd4.pre\"), ncol = 13, byrow=T)[,6:13]");
		writer.println("affstat." + i + "<-read.table(\"aff." + i + "\", na.strings=\".\")[,2]");
		writer.println("all.em." + i + "<-haplo.em(geno = geno)");
		writer.println("sink(\"" + i + "_em.out\")");
		writer.println("");
		writer.println("print(all.em." + i + ")");
		writer.println("sink()");
		writer.println("sink(\"" + i + "_haplos.out\")");
		writer.println("");
		writer.println("summary(all.em." + i + ")");
		writer.println("sink()");
		// writer.println("bintest."+i+"<-haplo.score(affstat."+i+", geno,
		// trait.type=\"binomial\", skip.haplo=0.005)");
		writer.println("bintest."	+ i + "<-haplo.score(affstat." + i
										+ ", geno, trait.type=\"binomial\")");
		writer.println("sink(\"" + i + "_bintest.out\")");
		writer.println("");
		writer.println("print(bintest." + i + ")");
		writer.println("sink()");
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		int numReps = 10000;

		String usage = "\n"	+ "park.haploGatherPerms requires 1 arguments:\n"
										+ "   (1) the number of replicates to perform (i.e. reps=" + numReps
										+ " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("reps=")) {
				numReps = Integer.valueOf(arg.split("=")[1]).intValue();
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (args.length == 0) {
			System.err.println("Using defaults (reps=" + numReps + ")");
		}

		try {
			new haploGatherPerms(numReps);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
