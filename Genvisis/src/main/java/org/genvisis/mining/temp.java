package org.genvisis.mining;

import java.io.*;

import org.genvisis.common.*;

public class temp {
	public temp(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		DoubleVector dv = new DoubleVector();

		double[] dist;

		try {
			reader = new BufferedReader(new FileReader("CRP.dat"));
			while (reader.ready()) {
				dv.add(Double.parseDouble(reader.readLine()));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+"CRP.dat"+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+"CRP.dat"+"\"");
			System.exit(2);
		}

		dist = dv.toArray();
		for (int k = 0; k<Transformations.NUM_TRANSFORMATIONS; k++) {
			System.err.println(Transformations.getLabel(k)+": "+ext.formDeci(Array.kurtosis(Transformations.transform(dist, k)), 4));
		}

		System.exit(1);

		try {
			reader = new BufferedReader(new FileReader("rank.txt"));
			while (reader.ready()) {
				dv.add(Double.parseDouble(reader.readLine()));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+"rank.txt"+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+"rank.txt"+"\"");
			System.exit(2);
		}

		double[] arr1 = dv.toArray();
		double[] arr2 = Transformations.rankTransform(arr1);

		for (int i = 0; i<arr1.length; i++) {
			System.out.println(arr1[i]+"\t"+arr2[i]);
		}

		System.exit(1);

		double[][] data = new double[250][];

		try {
			reader = new BufferedReader(new FileReader("kNNdata5.txt"));
			reader.readLine();
			for (int i = 0; i<data.length; i++) {
				line = reader.readLine().split("[\\s]+");
				data[i] = new double[] {Double.parseDouble(line[1]), Double.parseDouble(line[2])};
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+"kNNdata5.txt"+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+"kNNdata5.txt"+"\"");
			System.exit(2);
		}

		writer = new PrintWriter(new FileWriter("data5-std.xls"));
		data = Transformations.transform(data, Transformations.STANDARDIZE_RANGE);
		for (int i = 0; i<data.length; i++) {
			writer.println(data[i][0]+"\t"+data[i][1]);
		}
		writer.close();

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "temp.dat";

		String usage = "\n"+"park.temp requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default)\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			new temp(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
