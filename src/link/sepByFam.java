package link;

import java.io.*;
import java.util.*;

public class sepByFam {

	public sepByFam(String[] args) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st;
		String chrome, fam;
		int chromosome, position;
		double pos;
		boolean parametric = true;

		for (int i = 0; i<args.length; i++) {
			chromosome = Integer.valueOf(args[i].substring(0, args[i].indexOf("@"))).intValue();
			position = Integer.valueOf(args[i].substring(args[i].indexOf("@")+1)).intValue();
			chrome = (chromosome<10)?"0"+chromosome:""+chromosome;

			writer = new PrintWriter(new FileWriter("lods"+chrome+"_"+position+".dat.xls"));
			// writer = new PrintWriter(new FileWriter("Chromosome "+chrome+" at
			// "+position+"cM.out.xls"));

			if ((new File("chromf"+chrome+".lin.out")).exists()) {
				reader = new BufferedReader(new FileReader("chromf"+chrome+".lin.out"));
				writer.println("FamID\tlinearLOD\tNPL");
				parametric = false;
			} else if ((new File("chromf"+chrome+".D.out")).exists()) {
				reader = new BufferedReader(new FileReader("chromf"+chrome+".D.out"));
				writer.println("FamID\tAD LOD");
			} else if ((new File("chromf"+chrome+".d.out")).exists()) {
				reader = new BufferedReader(new FileReader("chromf"+chrome+".d.out"));
				writer.println("FamID\tAD LOD");
			} else if ((new File("chromf"+chrome+".R.out")).exists()) {
				reader = new BufferedReader(new FileReader("chromf"+chrome+".R.out"));
				writer.println("FamID\tAR LOD");
			} else if ((new File("chromf"+chrome+".r.out")).exists()) {
				reader = new BufferedReader(new FileReader("chromf"+chrome+".r.out"));
				writer.println("FamID\tAR LOD");
			} else {
				if (!allFour(args)) {
					System.err.println("Error: Not chromf"+chrome+".lin.out, chromf"+chrome+".D.out or chromf"+chrome+".R.out could be found for questioning.");
					System.err.println("       Furthermore, could not find them in their respective struct/directories either.");
				}
				System.exit(0);
			}

			try {
				reader.readLine();
				while (reader.ready()) {
					st = new StringTokenizer(reader.readLine());
					fam = st.nextToken();
					pos = Double.valueOf(st.nextToken()).doubleValue();
					// st.nextToken(); // gets rid of lod score for NPL sorting
					if (Math.abs(position-pos)<0.5) {
						writer.println(fam+"\t"+st.nextToken()+(parametric?"":"\t"+st.nextToken()));
					}
				}

				reader.close();
				System.out.println("Completed chromosome "+chromosome+" at "+position+"cM");
			} catch (Exception e) {
				System.err.println("Correct usage : java sepByFam 2@170 12@1 21@42 ... etc.");
				e.printStackTrace();
			}
			writer.close();
		}
	}

	public boolean allFour(String[] args) throws IOException {
		BufferedReader[] reader = new BufferedReader[3];
		PrintWriter writer = null;
		StringTokenizer st;
		String temp, chrome, fam;
		int chromosome, position;
		double pos;
		boolean positive;

		for (int i = 0; i<args.length; i++) {
			chromosome = Integer.valueOf(args[i].substring(0, args[i].indexOf("@"))).intValue();
			position = Integer.valueOf(args[i].substring(args[i].indexOf("@")+1)).intValue();
			chrome = (chromosome<10)?"0"+chromosome:""+chromosome;

			writer = new PrintWriter(new FileWriter("lods"+chrome+"_"+position+".dat.xls"));

			if (!(new File("allegro/chromf"+chrome+".lin.out")).exists()||!(new File("parametric dominant/chromf"+chrome+".d.out")).exists()||!(new File("parametric recessive/chromf"+chrome+".r.out")).exists()) {
				writer.close();
				return false;
			}

			reader[0] = new BufferedReader(new FileReader("allegro/chromf"+chrome+".lin.out"));
			reader[1] = new BufferedReader(new FileReader("parametric dominant/chromf"+chrome+".d.out"));
			reader[2] = new BufferedReader(new FileReader("parametric recessive/chromf"+chrome+".r.out"));

			writer.println("Family\tlinearLOD\tNPL\tAD LOD\tAR LOD\t\tpositive");

			try {
				reader[0].readLine();
				reader[1].readLine();
				reader[2].readLine();
				while (reader[0].ready()) {
					positive = false;
					st = new StringTokenizer(reader[0].readLine());
					fam = st.nextToken();
					pos = Double.valueOf(st.nextToken()).doubleValue();

					if (Math.abs(position-pos)<0.5) {
						temp = st.nextToken();
						if (Double.valueOf(temp).doubleValue()>0) {
							positive = true;
						}
						writer.print(fam+"\t"+temp);
						temp = st.nextToken();
						if (Double.valueOf(temp).doubleValue()>0) {
							positive = true;
						}
						writer.print("\t"+temp);
					}

					st = new StringTokenizer(reader[1].readLine());
					fam = st.nextToken();
					pos = Double.valueOf(st.nextToken()).doubleValue();

					if (Math.abs(position-pos)<0.5) {
						temp = st.nextToken();
						if (Double.valueOf(temp).doubleValue()>0) {
							positive = true;
						}
						writer.print("\t"+temp);
					}

					st = new StringTokenizer(reader[2].readLine());
					fam = st.nextToken();
					pos = Double.valueOf(st.nextToken()).doubleValue();

					if (Math.abs(position-pos)<0.5) {
						temp = st.nextToken();
						if (Double.valueOf(temp).doubleValue()>0) {
							positive = true;
						}
						writer.print("\t"+temp);

						if (positive) {
							writer.println("\t\t1");
						} else {
							writer.println("\t\t0");
						}
					}

				}

				reader[0].close();
				reader[1].close();
				reader[2].close();
				System.out.println("Completed chromosome "+chromosome+" at "+position+"cM");
			} catch (Exception e) {
				System.err.println("Failed allFour routine...");
				System.err.println("Correct usage : java sepByFam 2@170 12@1 21@42 ... etc.");
				e.printStackTrace();
			}
			writer.close();
		}

		return true;
	}

	public static void main(String[] args) throws IOException {
		if (args.length==0) {
			System.err.println("Correct usage : java sepByFam 2@170 12@1 21@42 ... etc.");
		}
		try {
			new sepByFam(args);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
