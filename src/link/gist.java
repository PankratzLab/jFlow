// sets up for GIST, assumes 2 is the deleterious allele
package link;

import java.io.*;
import java.util.*;

import common.*;

public class gist {
	public static String TARGET_FREQ = "0.45";

	public static double[][][] GIST_MATRIX = { { {1.00, 1.00, 1.00}, {1.00, 0.50, 0.25}, {0.50, 0.50, 0.50}}, { {1.00, 0.50, 0.75}, {1.00, 0.00, 0.50}, {0.50, 0.00, 0.25}}, { {0.50, 0.50, 0.50}, {0.50, 0.00, 0.25}, {0.00, 0.00, 0.00}}};

	public gist(String mutation, int chr, int pos, String target) throws IOException {
		BufferedReader reader;
		PrintWriter writer;
		String[] line = null, alleles;
		String trav, prev, chrome = (chr<10)?"0"+chr:""+chr;
		Hashtable<String,String> gtypes = new Hashtable<String,String>(), npls = new Hashtable<String,String>();
		boolean done;
		Vector<String> founders;
		double c1, c2, c3, n;

		if (!new File(mutation).exists()) {
			System.err.println("Error - could not find "+mutation+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(mutation));
		reader.readLine();
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			gtypes.put(line[1]+"\t"+line[2], line[3]+"\t"+line[4]);
		}
		reader.close();

		if (!new File("chromf"+chrome+".lin.out").exists()) {
			System.err.println("Error - could not find "+"chromf"+chrome+".lin.out"+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader("chromf"+chrome+".lin.out"));
		reader.readLine();
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			if (Math.abs(pos-Double.valueOf(line[1]).doubleValue())<0.5) {
				npls.put(line[0], line[3]);
			}
		}
		reader.close();

		if (!new File("re_chrom"+chrome+".pre").exists()) {
			System.err.println("Error - could not find "+"re_chrom"+chrome+".pre"+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader("re_chrom"+chrome+".pre"));
		writer = new PrintWriter(new FileWriter("gist-"+chr+"@"+pos+"-"+target+".dat"));
		writer.println(TARGET_FREQ);

		founders = new Vector<String>();
		prev = "";
		done = false;
		c1 = c2 = c3 = n = 0;
		while (!done) {
			if (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				trav = line[0];
			} else {
				done = true;
				trav = "";
			}

			if (!trav.equals(prev)) {
				if (prev!="") {
					if (founders.size()!=2) {
						System.err.println("Error - there should only be sibpairs in file "+"re_chrom"+chrome+".pre"+"; problem encountered with family "+prev);
						System.exit(3);
					}
					if (n<2) {
						System.err.println("Family "+prev+" thrown out for lack of enough informative individuals");
						System.err.println(n+"\t"+c1+"\t"+c2+"\t"+c3);
					} else {
						if (!npls.containsKey(prev)) {
							System.err.println("How come there's no NPL score for family "+prev+"?");
						} else {
							writer.println(ext.formDeci(c1/n, 3)+" "+ext.formDeci(c2/n, 3)+" "+ext.formDeci(c3/n/2, 3)+" "+npls.get(prev));
						}
					}
				}

				founders.removeAllElements();
				c1 = c2 = c3 = n = 0;
				prev = trav;
			}

			if (line[2].equals("0")&&line[3].equals("0")) {
				if (!founders.contains(line[1])) {
					founders.add(line[1]);
				}
			} else {
				if (!founders.contains(line[2])) {
					founders.add(line[2]);
				}
				if (!founders.contains(line[3])) {
					founders.add(line[3]);
				}
			}

			if (line[5].equals("2")&&gtypes.containsKey(line[0]+"\t"+line[1])) {
				n++;
				alleles = (gtypes.get(line[0]+"\t"+line[1])).split("[\\s]+");
				if (alleles[0].equals(target)||alleles[1].equals(target)) {
					c1++;
				}
				if (alleles[0].equals(target)&&alleles[1].equals(target)) {
					c2++;
				}
				if (alleles[0].equals(target)) {
					c3++;
				}
				if (alleles[1].equals(target)) {
					c3++;
				}
			}

		}

		reader.close();
		writer.close();

	}

	public static void main(String[] args) throws IOException {
		String mutation = "BRI3.dat", target = "2";
		// String mutation = "BRI3_reversed.dat";
		int chr = 2, pos = 233;

		String usage = "\n"+"park.gist requires 3 arguments:\n"+"   (1) a chromosome#.dat-like file with genotypes of the putative mutation (i.e. mut="+mutation+" (default))\n"+"   (2) chromosome number (i.e. chr="+chr+" (default))\n"+"   (3) position of SNP (i.e. pos="+pos+" (default))\n"+"   (4) target allele (i.e. allele="+target+" (default))\n"+"";
		int numArgs = args.length;

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("mut=")) {
				mutation = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("chr=")) {
				chr = Integer.valueOf(args[i].split("=")[1]).intValue();
				numArgs--;
			} else if (args[i].startsWith("pos=")) {
				pos = Integer.valueOf(args[i].split("=")[1]).intValue();
				numArgs--;
			} else if (args[i].startsWith("allele=")) {
				target = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		System.out.println("Using counts of "+mutation+" on the estimates at "+pos+" of chromosome "+chr);
		try {
			new gist(mutation, chr, pos, target);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
