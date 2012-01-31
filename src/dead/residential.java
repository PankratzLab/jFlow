package dead;

import java.io.*;
import java.util.*;

public class residential {

	public residential(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st = null;
		String temp, trav, prev;
		int ruralYears = 0, urbanYears = 0, year1, year2;
		int missingData = 0;

		reader = new BufferedReader(new FileReader(filename));
		writer = new PrintWriter(new FileWriter(filename.substring(0, filename.length()-4)+"-db.dat"));

		writer.println("Family\tIndividual\tRural Years\tUrban Years\t%Rural\t#years");

		reader.readLine();
		prev = "";
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			trav = st.nextToken()+"-"+Integer.valueOf(st.nextToken()).intValue();
			if (!trav.equals(prev)) {
				if ((ruralYears+urbanYears)>0) {
					writer.println(prev.substring(0, prev.indexOf("-"))+"\t"+prev.substring(prev.indexOf("-")+1)+"\t"+ruralYears+"\t"+urbanYears+"\t"+(int)(100*ruralYears/(ruralYears+urbanYears))+"%"+"\t"+(ruralYears+urbanYears));
				} else {
					missingData++;
				}
				prev = trav;
				ruralYears = 0;
				urbanYears = 0;
			}

			if (st.nextToken().equals("Y")) {
				st.nextToken();
				year1 = Integer.valueOf(st.nextToken()).intValue();
				year2 = Integer.valueOf(st.nextToken()).intValue();
				if ((year2-year1)<0||(year2-year1)>100) {
					System.err.println("Flag individual "+trav+" for year beginning "+year1);
				}
				ruralYears += year2-year1+1;
			}

			if (st.hasMoreTokens()&&st.nextToken().equals("Y")) {
				year1 = Integer.valueOf(st.nextToken()).intValue();
				year2 = Integer.valueOf(st.nextToken()).intValue();
				if ((year2-year1)<0||(year2-year1)>100) {
					System.err.println("Flag individual "+trav+" for year beginning "+year1);
				}
				urbanYears += year2-year1+1;
			}

		}
		if ((ruralYears+urbanYears)>0) {
			writer.println(prev.substring(0, prev.indexOf("-"))+"\t"+prev.substring(prev.indexOf("-")+1)+"\t"+ruralYears+"\t"+urbanYears+"\t"+(int)(100*ruralYears/(ruralYears+urbanYears))+"%"+"\t"+(ruralYears+urbanYears));
		}
		reader.close();
		writer.close();

		reader = new BufferedReader(new FileReader(filename.substring(0, filename.length()-4)+"-db.dat"));
		writer = new PrintWriter(new FileWriter(filename.substring(0, filename.length()-4)+"-db-summary.dat"));

		// writer.println("Rurs\tUrbs\tConR\tDis\tConU");

		int urbs = 0, rurs = 0, totalUrbs = 0, totalRurs = 0, conR = 0, dis = 0, conU = 0;
		int majRurYears = 0, majCons = 0;

		reader.readLine();
		prev = "";
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			trav = st.nextToken();
			if (!trav.equals(prev)) {
				if ((urbs+rurs)>1) {
					if (majRurYears>=2) {
						majCons++;
						// writer.println(prev+"\t>15 years");
					}
					if (rurs>=2) {
						conR++;
						// writer.println(prev+"\t>25% of life");
					} else {}
					if (rurs>=1&&urbs>=1) {
						dis++;
					} else {}
					if (urbs>=2) {
						writer.println(prev+"\t<10% of life");
						conU++;
					} else {}
				}
				urbs = 0;
				rurs = 0;
				majRurYears = 0;
				prev = trav;
			}
			st.nextToken();
			if (Integer.valueOf(st.nextToken()).intValue()>=15) {
				majRurYears++;
			}
			st.nextToken();

			temp = st.nextToken();
			if (Integer.valueOf(temp.substring(0, temp.length()-1)).intValue()>25) {
				rurs++;
				totalRurs++;
			} else {
				urbs++;
				totalUrbs++;
			}

		}
		writer.println("Lived in Rural areas > 25%: "+totalRurs);
		writer.println("Lived in Rural areas < 25%: "+totalUrbs);
		writer.println("Concordant Rurals: "+conR);
		writer.println("Concordant 15 years or more Rurals: "+majCons);
		writer.println("Discordant families: "+dis);
		writer.println("Concordant Urbans: "+conU);
		writer.println();
		writer.println("Missing Data: "+missingData+" individuals");

		reader.close();
		writer.close();
	}

	public static void main(String[] args) {
		if (args.length!=1) {
			System.out.println("Expecting 1 argument: filename");
		} else {
			try {
				new residential(args[0]);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
