package kaput;

import java.io.*;
import java.util.*;

public class reverseSolar {
	public reverseSolar() throws IOException {
		new reverseSolar(1, 22);
	}

	public reverseSolar(int start, int stop) throws IOException {
		BufferedReader genoReader = null, famReader = null;
		PrintWriter writer;
		StringTokenizer st;
		String chrome, FamID, IndID;
		int chromosome = 0;

		try {
			for (chromosome = start; chromosome<=stop; chromosome++) {
				chrome = (chromosome<10)?"0"+chromosome:""+chromosome;
				famReader = new BufferedReader(new FileReader("solar.fam"));
				genoReader = new BufferedReader(new FileReader("solar.gtypes."+chromosome));
				writer = new PrintWriter(new FileWriter("re_chrom"+chrome+".pre"));

				famReader.readLine();
				genoReader.readLine();

				while (famReader.ready()) {
					st = new StringTokenizer(famReader.readLine(), ",");
					FamID = st.nextToken();
					IndID = st.nextToken();
					writer.print(FamID+"\t"+IndID+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken());
					if (Integer.valueOf(IndID).intValue()<100) {
						writer.print("\t2");
					} else {
						writer.print("\t0");
					}

					st = new StringTokenizer(genoReader.readLine(), "/,");
					if (FamID.equals(st.nextToken())&&IndID.equals(st.nextToken())) {
						while (st.hasMoreTokens()) {
							writer.print("\t"+st.nextToken());
						}
						writer.println();
					} else {
						System.err.println("Error processing chromosome "+chromosome+" at individual "+FamID+"-"+IndID+".");
						System.exit(1);
					}
				}
				famReader.close();
				genoReader.close();
				writer.close();
			}
		} catch (Exception e) {
			System.err.println("Error processing chromosome "+chromosome+".");
			System.exit(1);
		}

	}

	public static void main(String[] args) throws IOException {
		if (args.length>2) {
			System.err.println("Error: maximum of 2 optional arguments - 1) chromosome to start with; and 2)chromosome to end with (doubly optional)");
		} else {
			try {
				if (args.length==0) {
					new reverseSolar();
				}
				if (args.length==1) {
					new reverseSolar(Integer.valueOf(args[0]).intValue(), Integer.valueOf(args[0]).intValue());
				}
				if (args.length==2) {
					new reverseSolar(Integer.valueOf(args[0]).intValue(), Integer.valueOf(args[1]).intValue());
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
