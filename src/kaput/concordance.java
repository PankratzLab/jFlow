package kaput;

import java.io.*;
import java.util.*;

public class concordance {
	public static int MMSE_CUTOFF_LTE = 25;

	public concordance(String filename, int numReps) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st;
		String trav, alleles, affStat, mmse;
		Hashtable<String,Vector<String>> hash = new Hashtable<String,Vector<String>>();
		Vector<String> families = new Vector<String>(), v;
		int group;
		int[] totals = new int[6], counts;
		String a1, a2;

		reader = new BufferedReader(new FileReader(filename));

		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			trav = st.nextToken();
			for (int i = 0; i<4; i++) {
				st.nextToken();
			}

			affStat = st.nextToken();
			alleles = st.nextToken();
			st.nextToken();
			st.nextToken(); // aoo
			st.nextToken();
			mmse = st.nextToken();
			// st.nextToken(); // add aae??

			if (alleles.equals("33")) {
				group = 1;
			} else if (alleles.equals("22")||alleles.equals("23")) {
				group = 2;
			} else if (alleles.equals("44")||alleles.equals("34")) {
				group = 3;
			} else {
				group = 0;
			}

			if (!mmse.equals(".")&&(group>0)&&affStat.equals("2")) {

				if (!hash.containsKey(trav)) {
					families.add(trav);
					hash.put(trav, v = new Vector<String>());
				} else {
					v = hash.get(trav);
				}
				v.add(group+"\t"+mmse);
			}
		}
		reader.close();

		// writer = new PrintWriter(new FileWriter("meanAPOE_MMSE_MCMC.out"));
		// writer.println("#aff4*\t#unaff4*\t#aff33\t#unaff33\t#aff2*\t#unaff2*");
		for (int rep = 0; rep<numReps; rep++) {
			counts = new int[6];
			writer = new PrintWriter(new FileWriter("count."+(rep+1)+".out"));
			for (int j = 0; j<families.size(); j++) {
				v = hash.get(families.elementAt(j));
				st = new StringTokenizer(v.elementAt(0));
				a1 = st.nextToken();
				a2 = st.nextToken();
				writer.println(a1+"\t"+a2);

				counts[2*(Integer.valueOf(a1).intValue()-1)+((Integer.valueOf(a2).intValue()<=MMSE_CUTOFF_LTE)?0:1)]++;
			}
			for (int j = 0; j<counts.length; j++) {
				System.out.print(counts[j]+"\t");
				totals[j] += counts[j];
			}
			System.out.println("");
			writer.close();
		}

	}

	public static void main(String[] args) throws IOException {
		String usage = "\n"+"park.concordance requires 2 arguments:\n"+"   (1) a pre-like file - FAMID INDID AFF/UNAFF group\n"+"       (for an APOE analysis this would be 2, 3 or 4)\n"+"   (2) number of replicates to perform\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			}
		}
		if (args.length!=2) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			new concordance(args[0], Integer.valueOf(args[1]).intValue());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
