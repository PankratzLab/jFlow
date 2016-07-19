package link.bat;

import java.io.*;
import java.util.*;

import common.*;

public class gatherOSA {

	public gatherOSA(int chromosome) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st = null;
		String temp, chrome;
		boolean done;
		int increment;

		chrome = (chromosome<10)?"0"+chromosome:""+chromosome;

		writer = new PrintWriter(new FileWriter("chrom"+chrome+"a-surface.xls"));
		for (int i = 0; i<=254; i++) {
			writer.print("\t"+i);
		}
		writer.println();

		increment = 1;
		done = false;
		while (!done) {
			try {
				reader = new BufferedReader(new FileReader("chrom"+chrome+"/chrom"+chrome+"-a"+ext.formNum(increment, 4)+".out"));

				if (chromosome<23) {
					writer.print(increment);
					reader.readLine();
					temp = reader.readLine();
					while (!temp.equals("")) {
						st = new StringTokenizer(temp);
						st.nextToken();
						st.nextToken();
						st.nextToken();
						st.nextToken();
						writer.print("\t"+st.nextToken());

						temp = reader.readLine();
					}
				} else {
					writer.print(increment);
					do {
						temp = reader.readLine();
					} while (!temp.startsWith("Total"));
					temp = reader.readLine();
					while (!temp.equals("")) {
						st = new StringTokenizer(temp);
						st.nextToken();
						writer.print("\t"+st.nextToken());

						temp = reader.readLine();
					}
				}

				writer.println();

				reader.close();
			} catch (Exception e) {
				System.out.println("Ending with a total of "+(increment-1)+" increments.");
				done = true;
			}
			increment++;
		}
		writer.close();

		writer = new PrintWriter(new FileWriter("chrom"+chrome+"d-surface.xls"));
		for (int i = 0; i<=254; i++) {
			writer.print("\t"+i);
		}
		writer.println();

		increment = 1;
		done = false;
		while (!done) {
			try {
				reader = new BufferedReader(new FileReader("chrom"+chrome+"/chrom"+chrome+"-d"+ext.formNum(increment, 4)+".out"));

				if (chromosome<23) {
					writer.print(increment);
					reader.readLine();
					temp = reader.readLine();
					while (!temp.equals("")) {
						st = new StringTokenizer(temp);
						st.nextToken();
						st.nextToken();
						st.nextToken();
						st.nextToken();
						writer.print("\t"+st.nextToken());

						temp = reader.readLine();
					}
				} else {
					writer.print(increment);
					do {
						temp = reader.readLine();
					} while (!temp.startsWith("Total")&&reader.ready());
					if (!reader.ready()) {
						writer.print("\t0.25");
					} else {
						temp = reader.readLine();
						while (!temp.equals("")) {
							st = new StringTokenizer(temp);
							st.nextToken();
							writer.print("\t"+st.nextToken());

							temp = reader.readLine();
						}
					}
				}

				writer.println();

				reader.close();
			} catch (Exception e) {
				System.out.println("Ending with a total of "+(increment-1)+" increments.");
				done = true;
			}
			increment++;
		}
		writer.close();

	}

	public static void main(String[] args) {
		if (args.length!=1) {
			System.out.println("Expecting 1 optional argument: chromosome number");
		} else {
			try {
				new gatherOSA(Integer.valueOf(args[0]).intValue());
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
