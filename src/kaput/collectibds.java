package kaput;

import java.io.*;
import java.util.*;

import common.*;

public class collectibds {

	public collectibds(String[] args) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st;
		String temp, chrome;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		int chromosome, position, count = 0, n;
		double[] dRay;
		double dub, sum, mean, variance, std, t, z;

		if ((new File(args[0])).exists()) {
			reader = new BufferedReader(new FileReader(args[0]));
		} else {
			System.err.println("Error: "+args[0]+" could be found for questioning.");
			System.err.println("Correct usage : java collectibds diskey_file 2@170 12@1 21@42 ... etc.");
			System.exit(0);
		}

		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			hash.put(st.nextToken()+st.nextToken(), "null");
			count++;
		}

		for (int i = 1; i<args.length; i++) {
			chromosome = Integer.valueOf(args[i].substring(0, args[i].indexOf("@"))).intValue();
			position = Integer.valueOf(args[i].substring(args[i].indexOf("@")+1)).intValue();
			chrome = (chromosome<10)?"0"+chromosome:""+chromosome;

			if ((new File("mibd"+chrome+".dat")).exists()) {
				reader = new BufferedReader(new FileReader("mibd"+chrome+".dat"));
			} else {
				System.err.println("Error: mibd"+chrome+".dat could be found for questioning.");
				System.exit(0);
			}

			try {
				// reader.readLine(); // why was this here to begin with? I
				// don't know...
				writer = new PrintWriter(new FileWriter(chrome+"@"+position+".out"));

				dRay = new double[count];
				n = 0;
				while (reader.ready()) {
					st = new StringTokenizer(reader.readLine());
					if ((int)Double.valueOf(st.nextToken()).doubleValue()==position) {
						temp = st.nextToken()+st.nextToken();
						if (hash.containsKey(temp)) {
							// hash.remove(temp); // use if you want to know the
							// untyped individuals
							st.nextToken();
							dub = Double.valueOf(st.nextToken()).doubleValue()/2+Double.valueOf(st.nextToken()).doubleValue();
							writer.println(ext.formDeci(dub, 5));
							dRay[n] = dub;
							n++;
						}
					}
				}

				reader.close();
				System.out.println("Completed "+chromosome+" at "+position+"cM");

				sum = 0;
				for (int j = 0; j<n; j++) {
					sum += dRay[j];
				}
				mean = sum/n;
				sum = 0;
				for (int j = 0; j<n; j++) {
					sum += dRay[j]*dRay[j];
				}
				variance = (sum/n-mean*mean)*((double)n/((double)n-1));
				std = Math.sqrt(variance);
				t = (mean-.5)/(std/Math.sqrt(n));
				System.out.println("t = "+t);

				// sigt = pt.stDist(df, Math.abs(t));
				// System.out.println("p-value = "+sigt);

				z = (mean-.5)/(Math.sqrt(.5*(1-.5)/n));
				System.out.println("z = "+z);
				// sigz = pt.stDist(325, Math.abs(z));
				// System.out.println("p-value = "+sigz);

				// p=1-2*sig_t;
				// pz=1-2*sig_z;
				// lod=log10(exp(cinv(p,1)/2));
				// lodz=log10(exp(cinv(pz,1)/2));

			} catch (Exception e) {
				System.err.println("Correct usage : java collectibds diskey_file 2@170 12@1 21@42 ... etc.");
				e.printStackTrace();
			}
			writer.close();
		}
	}

	public static void main(String[] args) throws IOException {
		if (args.length==0) {
			System.err.println("Correct usage : java collectibds diskey_file 2@170 12@1 21@42 ... etc.");
		}
		try {
			new collectibds(args);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
