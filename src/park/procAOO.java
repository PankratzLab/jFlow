package park;

import java.io.*;
import java.util.*;

import common.*;

public class procAOO {
	public static String MISSING_VALUE1 = "-99";

	public static String MISSING_VALUE2 = ".";

	public static int DECIMAL_POINTS = 2;

	public procAOO(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer, means, ghMeans, first, second, ghFirst, ghSecond, last, penultumate, ghLast, ghPenultumate, diff;
		String[] line;
		String trav, prev, val;
		Vector<String> values = new Vector<String>();
		String trait = filename.substring(0, filename.lastIndexOf("."));

		if (trait.equals("struct")) {
			trait = "AOO";
		}

		reader = new BufferedReader(new FileReader(filename));
		writer = new PrintWriter(new FileWriter("procd"+trait+".dat"));
		means = new PrintWriter(new FileWriter("mean"+trait+".dat"));
		ghMeans = new PrintWriter(new FileWriter("gh_mean"+trait+".dat"));
		first = new PrintWriter(new FileWriter("1st"+trait+".dat"));
		second = new PrintWriter(new FileWriter("2nd"+trait+".dat"));
		last = new PrintWriter(new FileWriter("last"+trait+".dat"));
		penultumate = new PrintWriter(new FileWriter("pu"+trait+".dat"));
		ghLast = new PrintWriter(new FileWriter("gh_last"+trait+".dat"));
		ghPenultumate = new PrintWriter(new FileWriter("gh_pu"+trait+".dat"));
		ghFirst = new PrintWriter(new FileWriter("gh_1st"+trait+".dat"));
		ghSecond = new PrintWriter(new FileWriter("gh_2nd"+trait+".dat"));
		diff = new PrintWriter(new FileWriter("diff"+trait+".dat"));
		writer.println("Fam\tMean\t1st\t2nd");

		prev = "";
		// reader.readLine(); // for gaw
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			trav = line[0];

			val = line[filename.startsWith("struct")?7:2];

			if ((!trav.equals(prev)&&!prev.equals(""))||!reader.ready()) {
				if (!reader.ready()&&!val.equals(MISSING_VALUE1)&&!val.equals(MISSING_VALUE2)) {
					values.add(val);
				}

				line = pertinent(values).split("[\\s]+");
				writer.println(prev+"\t"+pertinent(values));
				if (!line[0].equals(".")) {
					means.println(prev+"\t"+line[0]);
				}
				ghMeans.println(prev+"\t"+line[0]);

				if (!line[1].equals(".")) {
					first.println(prev+"\t"+line[1]);
				}
				ghFirst.println(prev+"\t"+line[1]);

				if (!line[2].equals(".")) {
					last.println(prev+"\t"+line[2]);
				}
				ghLast.println(prev+"\t"+line[2]);

				if (!line[3].equals(".")) {
					second.println(prev+"\t"+line[3]);
				}
				ghSecond.println(prev+"\t"+line[3]);

				if (!line[4].equals(".")) {
					penultumate.println(prev+"\t"+line[4]);
				}
				ghPenultumate.println(prev+"\t"+line[4]);

				if (!line[5].equals(".")) {
					diff.println(prev+"\t"+line[5]);
				}
				values.removeAllElements();
			}

			if (!val.equals(MISSING_VALUE1)&&!val.equals(MISSING_VALUE2)) {
				values.add(val);
			}

			prev = trav;
		}
		writer.close();
		means.close();
		ghMeans.close();
		first.close();
		second.close();
		ghFirst.close();
		ghSecond.close();
		last.close();
		penultumate.close();
		ghLast.close();
		ghPenultumate.close();
		diff.close();
	}

	public String pertinent(Vector<String> phenos) {
		// ext.formDeci 1, true for aoo --> 4, true for the others
		String str = "";
		int[] keys;
		int size = phenos.size();
		double[] realPhenos;
		double total = 0;
//		int count;

		if (size>0) {
			realPhenos = new double[size];
			for (int i = 0; i<size; i++) {
				realPhenos[i] = Double.valueOf(phenos.elementAt(i)).doubleValue();
				total += realPhenos[i];
			}
			str += ext.formDeci(total/(double)size, DECIMAL_POINTS, true)+"\t";

			keys = Sort.quicksort(realPhenos);
			str += ext.formDeci(realPhenos[keys[0]], 1)+"\t";
			str += ext.formDeci(realPhenos[keys[realPhenos.length-1]], 1)+"\t";

			if (size>1) {
				str += ext.formDeci(realPhenos[keys[1]], 1)+"\t";
				str += ext.formDeci(realPhenos[keys[realPhenos.length-2]], 1)+"\t";

				total = 0;
//				count = 0;
				for (int i = 0; i<keys.length-1; i++) {
					for (int j = i+1; j<keys.length; j++) {
						total += realPhenos[keys[j]]-realPhenos[keys[i]];
//						count++;
					}
				}
				str += ext.formDeci(total/(double)size, 1, true);

			} else {
				str += ".\t.\t.";
			}
		} else {
			str += ".\t.\t.\t.\t.\t.";
		}

		return str;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		// String filename = "duration_residuals.dat";
		String filename = "struct.dat";

		String usage = "\n"+"park.procAOO requires 0-1 arguments\n"+"   (1) the file to be processed (default: file="+filename+") or -all\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				if (!new File(filename).exists()) {
					System.err.println("Error - file '"+filename+"' does not exist");
					System.exit(2);
				}
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (args.length==0) {
			System.out.println("Warning - using defaults (processing data in "+filename);
		}
		try {
			new procAOO(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
