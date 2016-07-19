package kaput;

import java.io.*;
import java.util.*;
import stats.ProbDist;

public class computeLODs {

	public computeLODs(String[] args) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st;
		String temp, chrome, trav;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		int start, stop, n = 0, numLoci;
		double[][] dRay;
		double dub;

		if ((new File(args[0])).exists()) {
			reader = new BufferedReader(new FileReader(args[0]));
		} else {
			System.err.println("Error: "+args[0]+" could be found for questioning.");
			System.err.println("Correct usage : java collectibds diskey_file (Optional: chromosome_number, default all)");
			System.exit(0);
		}

		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			hash.put(st.nextToken()+st.nextToken(), "null");
			n++;
		}

		if (args.length>1) {
			start = stop = Integer.valueOf(args[1]).intValue();
		} else {
			start = 1;
			stop = 23;
		}

		for (int chromosome = start; chromosome<=stop; chromosome++) {
			chrome = (chromosome<10)?"0"+chromosome:""+chromosome;

			if ((new File("mibd"+chrome+".dat")).exists()) {
				reader = new BufferedReader(new FileReader("mibd"+chrome+".dat"));
			} else {
				System.err.println("Error: mibd"+chrome+".dat could be found for questioning.");
				System.exit(0);
			}

			try {
				writer = new PrintWriter(new FileWriter("chrom"+chrome+".out"));

				reader.readLine();
				numLoci = 0;
				do {
					temp = reader.readLine();
					st = new StringTokenizer(temp);
					numLoci++;
				} while (!st.nextToken().equals("-0.000"));

				reader.close();
				reader = new BufferedReader(new FileReader("mibd"+chrome+".dat"));

				dRay = new double[numLoci][n]; // create table size num loci x
				// num expected people
				n = 0; // determine actual num people
				temp = reader.readLine();
				while (reader.ready()) {
					st = new StringTokenizer(temp);
					st.nextToken();
					trav = st.nextToken()+st.nextToken();
					if (hash.containsKey(trav)) {
						// hash.remove(temp); // use if you want to know the
						// untyped individuals
						for (int i = 0; i<numLoci; i++) {
							st = new StringTokenizer(temp);
							st.nextToken();
							st.nextToken();
							st.nextToken();
							st.nextToken();
							dub = Double.valueOf(st.nextToken()).doubleValue()/2+Double.valueOf(st.nextToken()).doubleValue();
							// writer.println(park.ext.formDeci(dub, 5));
							dRay[i][n] = dub;
							if (reader.ready()) {
								temp = reader.readLine();
							}
						}
						n++;
					} else {
						for (int i = 0; i<numLoci; i++) {
							if (reader.ready()) {
								temp = reader.readLine();
							}
						}
					}
				}
				reader.close();

				for (int i = 225; i<235; i++) {
					System.out.println(calculatePVal(dRay[i], n));
					writer.println(Math.log(Math.exp(ProbDist.ChiDistReverse(2*calculatePVal(dRay[i], n), 1)/2))/Math.log(10));
				}
				writer.close();

				System.out.println("Completed chromosome "+chromosome);
			} catch (Exception e) {
				System.err.println("Correct usage : java collectibds diskey_file (Optional: chromosome_number, default all)");
				e.printStackTrace();
			}
		}
	}

	public double calculatePVal(double[] locus, int n) {
		// double sum, mean, variance, std, t, z, sigt=-999;
		// sum = 0;
		// for (int j=0; j<n; j++) {
		// sum += locus[j];
		// }
		// mean = sum/n;
		//
		// sum = 0;
		// for (int j=0; j<n; j++) {
		// sum += locus[j] * locus[j];
		// }
		// variance = (sum/n - mean*mean)*((double)n/((double)n-1));
		// std = Math.sqrt(variance);
		//
		// t = (mean-.5)/(std/Math.sqrt(n));
		// sigt = pt.stDist(df, Math.abs(t));
		//
		//
		// z = (mean-.5)/(Math.sqrt(.5*(1-.5)/n));
		// sigz = pt.stDist(325, Math.abs(z));
		//
		// p=1-2*sig_t;
		// pz=1-2*sig_z;
		// lod=log10(exp(cinv(p,1)/2));
		// lodz=log10(exp(cinv(pz,1)/2));
		//
		// return sigt;

		return -999;
	}

	public static void main(String[] args) throws IOException {
		if (args.length==0) {
			System.err.println("Correct usage : java collectibds diskey_file (Optional: chromosome_number, default all)");
		}
		try {
			new computeLODs(args);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
