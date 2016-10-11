package org.genvisis.assoc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.common.Array;

public class haploPerm {
	public haploPerm(int numReps) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		Vector<String> inds = new Vector<String>();
		Vector<String> affstat = new Vector<String>();
		int n;
		int[] randKey;

		if (!new File("aft.dat").exists()) {
			System.err.println("Error - could not find " + "aft.dat" + " in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader("aft.dat"));
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			inds.add(line[0]);
			affstat.add(line[1]);
		}
		reader.close();
		n = inds.size();

		for (int i = 1; i <= numReps; i++) {
			writer = new PrintWriter(new FileWriter("aff." + i));
			randKey = Array.random(n);
			for (int j = 0; j < n; j++) {
				writer.println(inds.elementAt(j) + "\t" + affstat.elementAt(randKey[j]));
			}
			writer.close();
			writer = new PrintWriter(new FileWriter("perm." + i));
			writer.println("library(haplo.stats)");
			writer.println("geno<-matrix(scan(\"first3.pre\"), ncol = 11, byrow=T)[,6:11]");
			// writer.println("geno<-matrix(scan(\"pd4.pre\"), ncol = 13,
			// byrow=T)[,6:13]");
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
											+ ", geno, trait.type=\"binomial\", skip.haplo=0.02)");
			// writer.println("bintest."+i+"<-haplo.score(affstat."+i+", geno,
			// trait.type=\"binomial\")");
			writer.println("sink(\"" + i + "_bintest.out\")");
			writer.println("");
			writer.println("print(bintest." + i + ")");
			writer.println("sink()");
			writer.close();
		}

		writer = new PrintWriter(new FileWriter("batchHaploPerm"));
		for (int i = 1; i <= numReps; i++) {
			writer.println("Splus7.0 BATCH perm." + i + " check." + i + ".out");
			writer.println("sleep 6");
		}
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		int numReps = 10;

		String usage = "\n"	+ "park.haploPerm requires 1 arguments:\n"
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
			new haploPerm(numReps);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
