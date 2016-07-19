package org.genvisis.link;

import java.io.*;

public class simRunAllegro {
	public simRunAllegro(int reps, boolean all) throws IOException {
		PrintWriter optfile = null;
		PrintWriter writer = new PrintWriter(new FileWriter("batch.1"));
		String chrome;

		for (int i = 1; i<=reps; i++) {
			for (int chr = 1; chr<=23; chr++) {
				chrome = (chr<10)?"0"+chr:""+chr;
				optfile = new PrintWriter(new FileWriter("temp-"+(all?chrome+"-":"")+i+".opt"));
				optfile.println("% Read input in LINKAGE style format:\n"+"PREFILE linkage-"+(all?chrome+"-":"")+i+".pre\n"+"DATFILE "+(all?"map"+chrome+".dat":"linkage.dat")+"\n\n"+"% Simulate stroke reconstruction pedigrees\n"+"MODEL mpt lin all equal output-"+(all?chrome+"-":"")+i+".prn trash\n\n"+"% Other options:\n"+"STEPS 10\n"+"MAXMEMORY 100");
				optfile.close();

				writer.println("/software/bin/allegro temp-"+(all?chrome+"-":"")+i+".opt");
			}
		}
		writer.close();
		try {
			Runtime.getRuntime().exec("chmod +x batch.1").waitFor();
			Runtime.getRuntime().exec("./batch.1");
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		int numReps = 1000;
		boolean allChrs = false;

		String usage = "\n"+"park.simRunAllegro requires 0-2 arguments\n"+"   (1) number of replicates (i.e. reps="+numReps+" (default)\n"+"   (2) all chromosomes (i.e. '-all' (optional))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("reps=")) {
				numReps = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].equals("-all")) {
				allChrs = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (numReps==1000) {
			System.out.println("Running "+numReps+" replicates.");
		}
		try {
			new simRunAllegro(numReps, allChrs);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
