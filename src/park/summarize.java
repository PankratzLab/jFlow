package park;

import java.io.*;
import java.util.*;

public class summarize {
	public int[] chrs = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

	public summarize(String[] arguments) throws IOException {
		BufferedReader reader = null, reader2 = null;
		PrintWriter writer = null;
		StringTokenizer st;
		String temp, chrome, name, dist;
		Vector<String> markerNames = new Vector<String>();
		Vector<String> markerDists = new Vector<String>();
		// double before, alpha, after;
		boolean[][] missing = new boolean[8][23];

		for (int i = 0; i<arguments.length; i++) {
			try {
				chrs = new int[1];
				chrs[0] = Integer.valueOf(arguments[i]).intValue();
			} catch (Exception e) {
				System.err.println("Error: invalid arguments (valid = command [chromosome #] [run]");
			}
		}

		for (int whichChr = 0; whichChr<chrs.length; whichChr++) {
			int i = chrs[whichChr];
			chrome = (i<10)?"0"+i:""+i;
			writer = new PrintWriter(new FileWriter("summary"+chrome+".out"));
			markerNames.removeAllElements();
			markerDists.removeAllElements();

			writer.println("Mapmaker\tno Dv\tDv");
			try {
				reader = new BufferedReader(new FileReader("mapmaker/chrom"+chrome+"-mls.out"));
				if (i<23) {
					reader2 = new BufferedReader(new FileReader("mapmaker/chrom"+chrome+"-Dv-mls.out"));
					reader.readLine();
					reader2.readLine();
					temp = reader.readLine();
					while (!temp.equals("")) {
						st = new StringTokenizer(temp);
						writer.print(st.nextToken()+"\t");
						st.nextToken();
						st.nextToken();
						st.nextToken();
						writer.print(st.nextToken()+"\t");
						st = new StringTokenizer(reader2.readLine());
						st.nextToken();
						st.nextToken();
						st.nextToken();
						st.nextToken();
						writer.println(st.nextToken());

						temp = reader.readLine();
					}
				} else {
					do {
						temp = reader.readLine();
					} while (!temp.startsWith("Total"));
					temp = reader.readLine();
					while (!temp.equals("")) {
						st = new StringTokenizer(temp);
						writer.print(st.nextToken()+"\t");
						temp = st.nextToken();
						writer.println(temp+"\t"+temp);

						temp = reader.readLine();
					}
				}

			} catch (FileNotFoundException fne) {
				writer.println("no Mapmaker data");
				missing[0][i-1] = true;
			} catch (Exception e) {
				writer.println("no Mapmaker data");
				e.printStackTrace();
				System.err.println("Error processing Mapmaker file for chromosome "+i);
			}

			writer.println();
			writer.println();
			writer.println("Merlin-sibpairs");
			try {
				reader = new BufferedReader(new FileReader("merlin-/chrom"+chrome+".out"));
				do {
					temp = reader.readLine();
				} while (!temp.startsWith("========"));
				reader.readLine();
				reader.readLine();
				reader.readLine();
				while (reader.ready()) {
					temp = reader.readLine();
					st = new StringTokenizer(temp);
					if (st.hasMoreTokens()) {
						dist = st.nextToken();
						writer.print(dist+"\t");
						st.nextToken();
						st.nextToken();
						st.nextToken();
						writer.println(st.nextToken());
					}
				}
			} catch (FileNotFoundException fne) {
				writer.println("no Merlin- data");
				missing[1][i-1] = true;
			} catch (Exception e) {
				writer.println("no Merlin- data");
				e.printStackTrace();
				System.err.println("Error processing Merlin- files for chromosome "+i);
			}

			writer.println();
			writer.println();
			writer.println("Merlin-extended");
			try {
				reader = new BufferedReader(new FileReader("merlin+/chrom"+chrome+".out"));
				do {
					temp = reader.readLine();
				} while (!temp.startsWith("========"));
				reader.readLine();
				reader.readLine();
				reader.readLine();
				while (reader.ready()) {
					temp = reader.readLine();
					st = new StringTokenizer(temp);
					if (st.hasMoreTokens()) {
						dist = st.nextToken();
						writer.print(dist+"\t");
						st.nextToken();
						st.nextToken();
						st.nextToken();
						writer.println(st.nextToken());
					}
				}
			} catch (FileNotFoundException fne) {
				writer.println("no Merlin+ data");
				missing[2][i-1] = true;
			} catch (Exception e) {
				writer.println("no Merlin+ data");
				e.printStackTrace();
				System.err.println("Error processing Merlin+ files for chromosome "+i);
			}

			writer.println();
			writer.println();
			writer.println("Allegro\tlinear");
			try {
				reader = new BufferedReader(new FileReader("allegro/chrom"+chrome+".lin.out"));
				reader.readLine();
				while (reader.ready()) {
					temp = reader.readLine();
					st = new StringTokenizer(temp);
					dist = st.nextToken();
					writer.println(dist+"\t"+st.nextToken());
					st.nextToken();
					st.nextToken();
					st.nextToken();
					name = st.nextToken();
					if (!name.equals("-")) {
						markerNames.add(name);
						markerDists.add(dist);
					}
				}
			} catch (FileNotFoundException fne) {
				writer.println("no Allegro data");
				missing[3][i-1] = true;
			} catch (Exception e) {
				writer.println("no Allegro data");
				e.printStackTrace();
				System.err.println("Error processing Allegro files for chromosome "+i);
			}

			writer.println();
			writer.println();
			// writer.println("Aspex2pt\t\tKosambi");
			// try {
			// reader = new BufferedReader(new
			// FileReader("aspex/chrom"+chrome+"-K-2pt.out"));
			//
			// reader.readLine();
			// for (int j=0; j<markerNames.size(); j++) {
			// st = new StringTokenizer(reader.readLine());
			// for (int k=0; k<9; k++) {
			// st.nextToken();
			// }
			// if (i<23) {
			// st.nextToken();
			// }
			// writer.println(markerNames.elementAt(j)+"\t"+markerDists.elementAt(j)+"\t"+st.nextToken());
			// }
			// if (markerNames.size()==0) {
			// writer.println("no Aspex data");
			// }
			// } catch (FileNotFoundException fne) {
			// writer.println("no Aspex data");
			// missing[4][i-1] = true;
			// } catch (Exception e) {
			// writer.println("no Aspex data");
			// e.printStackTrace();
			// System.err.println("Error processing Aspex 2pt files for
			// chromosome "+i);
			// }
			//
			// writer.println();
			// writer.println();
			// writer.println("AspexMpt\tKosambi");
			// try {
			// reader = new BufferedReader(new
			// FileReader("aspex/chrom"+chrome+"-K-mpt.out"));
			//
			// do {
			// temp = reader.readLine();
			// } while (!temp.startsWith("0.094"));
			//
			// while (reader.ready()) {
			// st = new StringTokenizer(reader.readLine());
			// writer.print(ext.formDeci(Double.valueOf(st.nextToken()).doubleValue()*100
			// -10, 1, true));
			// for (int j=0; j<(i<23?4:3); j++) {
			// st.nextToken();
			// }
			// writer.println("\t"+st.nextToken());
			// }
			// } catch (FileNotFoundException fne) {
			// writer.println("no Aspex data");
			// missing[5][i-1] = true;
			// } catch (Exception e) {
			// writer.println("no Aspex data");
			// e.printStackTrace();
			// System.err.println("Error processing Aspex mpt files for
			// chromosome "+i);
			// }
			//
			// writer.println();
			// writer.println();
			writer.println("autosomal dominant");
			try {
				reader = new BufferedReader(new FileReader("parametric dominant/chrom"+chrome+".d.out"));
				reader.readLine();
				while (reader.ready()) {
					temp = reader.readLine();
					st = new StringTokenizer(temp);
					writer.println(st.nextToken()+"\t"+st.nextToken());
				}
			} catch (FileNotFoundException fne) {
				writer.println("no autosomal dominant data");
				missing[6][i-1] = true;
			} catch (Exception e) {
				writer.println("no autosomal dominant data");
				e.printStackTrace();
				System.err.println("Error processing autosomal dominant files for chromosome "+i);
			}

			writer.println();
			writer.println();
			writer.println("autosomal dominant het");
			try {
				reader = new BufferedReader(new FileReader("parametric dominant/chrom"+chrome+".d.out"));
				reader.readLine();
				while (reader.ready()) {
					temp = reader.readLine();
					st = new StringTokenizer(temp);
					writer.print(st.nextToken());
					st.nextToken();
					st.nextToken();
					// before = Double.valueOf(st.nextToken()).doubleValue();
					// alpha = Double.valueOf(st.nextToken()).doubleValue();
					// before = Double.valueOf(st.nextToken()).doubleValue();

					writer.println("\t"+st.nextToken());
				}
			} catch (FileNotFoundException fne) {
				writer.println("no autosomal dominant data");
				missing[6][i-1] = true;
			} catch (Exception e) {
				writer.println("no autosomal dominant het data");
				e.printStackTrace();
				System.err.println("Error processing autosomal dominant files for chromosome "+i);
			}

			writer.println();
			writer.println();
			writer.println("autosomal recessive");
			try {
				reader = new BufferedReader(new FileReader("parametric recessive/chrom"+chrome+".r.out"));
				reader.readLine();
				while (reader.ready()) {
					temp = reader.readLine();
					st = new StringTokenizer(temp);
					writer.println(st.nextToken()+"\t"+st.nextToken());
				}
			} catch (FileNotFoundException fne) {
				writer.println("no autosomal recessive data");
				missing[7][i-1] = true;
			} catch (Exception e) {
				writer.println("no autosomal recessive data");
				e.printStackTrace();
				System.err.println("Error processing autosomal recessive files for chromosome "+i);
			}

			writer.println();
			writer.println();
			writer.println("autosomal recessive het");
			try {
				reader = new BufferedReader(new FileReader("parametric recessive/chrom"+chrome+".r.out"));
				reader.readLine();
				while (reader.ready()) {
					temp = reader.readLine();
					st = new StringTokenizer(temp);
					writer.print(st.nextToken());
					st.nextToken();
					st.nextToken();
					// before = Double.valueOf(st.nextToken()).doubleValue();
					// alpha = Double.valueOf(st.nextToken()).doubleValue();
					// before = Double.valueOf(st.nextToken()).doubleValue();

					writer.println("\t"+st.nextToken());
				}
			} catch (FileNotFoundException fne) {
				writer.println("no autosomal recessive data");
				missing[7][i-1] = true;
			} catch (Exception e) {
				writer.println("no autosomal recessive het data");
				e.printStackTrace();
				System.err.println("Error processing autosomal recessive files for chromosome "+i);
			}

			writer.println();
			writer.close();
		}

		System.err.println("  "+missed("Mapmaker", missing[0]));
		System.err.println("  "+missed("Merlin-", missing[1]));
		System.err.println("  "+missed("Merlin+", missing[2]));
		System.err.println("  "+missed("Allegro", missing[3]));
		// System.err.println(" "+missed("Aspex 2pt", missing[4]));
		// System.err.println(" "+missed("Aspex mpt", missing[5]));
		System.err.println("  "+missed("autosomal dominant", missing[6]));
		System.err.println("  "+missed("autosomal recessive", missing[7]));

	}

	public String missed(String method, boolean[] missing) throws IOException {
		String str;
		int count = 0;

		for (int i = 0; i<missing.length; i++) {
			if (missing[i]) {
				count++;
			}
		}
		if (count==0) {
			return "Complete data for "+method;
		}
		str = "Missing "+method+" data for chromosome"+(count==1?"":"s")+" ";
		for (int i = 0; i<missing.length; i++) {
			if (missing[i]) {
				str += (str.endsWith(" ")?"":",")+(i+1)+"";
				count = i;
				while (count<missing.length&&missing[count]) {
					count++;
				}
				if (count-1>i) {
					str += "-"+(count);
				}
				i = count;
			}
		}
		return str;
	}

	public static void main(String[] args) throws IOException {
		try {
			new summarize(args);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
