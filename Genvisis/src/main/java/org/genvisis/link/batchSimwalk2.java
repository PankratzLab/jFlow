package org.genvisis.link;

import java.io.*;
import java.util.*;

public class batchSimwalk2 {
	@SuppressWarnings("resource")
	public batchSimwalk2(String models) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null, writer1 = null, writer2 = null, writer3 = null, writer4 = null;
		String chrome;
		Vector<double[]> modelParams = new Vector<double[]>();
		StringTokenizer st;
		double[] handle;

		reader = new BufferedReader(new FileReader(models));

		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			handle = new double[4];
			handle[0] = Double.valueOf(st.nextToken()).doubleValue();
			handle[1] = Double.valueOf(st.nextToken()).doubleValue();
			handle[2] = Double.valueOf(st.nextToken()).doubleValue();
			handle[3] = Double.valueOf(st.nextToken()).doubleValue();
			modelParams.add(handle);
		}
		reader.close();

		writer1 = new PrintWriter(new FileWriter("batch.1"));
		writer2 = new PrintWriter(new FileWriter("batch.2"));
		writer3 = new PrintWriter(new FileWriter("batch.3"));
		writer4 = new PrintWriter(new FileWriter("batch.4"));
		writer1.println("#/bin/sh");
		writer1.println();
		writer1.println("sleep 30");
		writer1.println();
		writer2.println("#/bin/sh");
		writer2.println();
		writer2.println("sleep 30");
		writer2.println();
		writer3.println("#/bin/sh");
		writer3.println();
		writer3.println("sleep 30");
		writer3.println();
		writer4.println("#/bin/sh");
		writer4.println();
		writer4.println("sleep 30");
		writer4.println();

		for (int i = 1; i<=22; i++) {
			chrome = (i<10)?"0"+i:""+i;

			if (i>=1&&i<=4) {
				writer = writer1;
			}
			if (i>=5&&i<=9) {
				writer = writer2;
			}
			if (i>=10&&i<=15) {
				writer = writer3;
			}
			if (i>=16&&i<=22) {
				writer = writer4;
			}
			if (writer==null) {
				System.err.println("Error: Not all chromosomes were told which file to be in.");
			}

			for (int j = 9; j<9+modelParams.size(); j++) {
				handle = modelParams.elementAt(j-9);
				writer.println("cd chr"+i);
				writer.println("mkdir Model"+(j+1));
				writer.println("cp re_chrom"+chrome+".pre Model"+(j+1)+"/pedin."+chrome);
				writer.println("cp map"+chrome+".dat Model"+(j+1)+"/datain.dat");
				writer.println("cp map."+chrome+" Model"+(j+1)+"/");
				writer.println("cd Model"+(j+1));
				writer.println("jcp makeMap4MLINK "+handle[0]+" "+handle[1]+" "+handle[2]+" "+handle[3]); // deprecated
				writer.println("mv datain.dat datain."+chrome);
				writer.println("dos2unix map."+chrome+" map."+chrome);
				writer.println("dos2unix pedin."+chrome+" pedin."+chrome);
				writer.println("dos2unix datain."+chrome+" datain."+chrome);
				writer.println("echo -e \"1\\n"+i+"\\n0\\n1\\n2\\n1\\n0\\n\" | mega2 > batchSimwalk2.log");
				writer.println("./loc."+chrome+".sh");
				writer.println("cd ../..");
				writer.println();
			}

			writer.println();
			writer.println();
		}

		writer1.close();
		writer2.close();
		writer3.close();
		writer4.close();
	}

	public static void main(String[] args) throws IOException {
		if (args.length!=1) {
			System.out.println("Expecting 1 argument: name of file containing dx models (Example: 0.00001 0.03 0.80 0.80)");
		} else {
			try {
				new batchSimwalk2(args[0]);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
