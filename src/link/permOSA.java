package link;

import java.io.*;
import java.util.*;

import common.*;

public class permOSA {
	public static int AROUND_LOCUS = 50; // 50 cM around each side of the

	// peak locus

	public permOSA(String filename, int chromosome, int numPerms, String direction) throws IOException {
		new permOSA(filename, chromosome, numPerms, direction, null);
	}

	public permOSA(String filename, int chromosome, int numPerms, String direction, String plug) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st = null;
		String temp, key, data, chrome, last, next;
		Vector<String> IDs = new Vector<String>();
		double[] covars;
		double trav, prev;
		int[] keys, randomKeys;
		int increment, numOfIncrements;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		double lod, maxLod;
		int maxIncrement, pos, maxPos;
		boolean problem;

		reader = new BufferedReader(new FileReader(filename));
		temp = reader.readLine();
		st = new StringTokenizer(temp);
		key = st.nextToken();
		data = st.nextToken();
		hash.put(key, data);
		IDs.add(key);
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			key = st.nextToken();
			data = st.nextToken();
			hash.put(key, data);
			IDs.add(key);
		}
		reader.close();

		covars = new double[IDs.size()];
		for (int i = 0; i<IDs.size(); i++) {
			covars[i] = Double.valueOf(hash.get(IDs.elementAt(i))).doubleValue();
		}
		keys = Sort.quicksort(covars);
		hash.clear();

		for (int perm = 1; perm<=numPerms; perm++) {
			if (plug!=null&&!(new File(plug).exists())) {
				System.exit(5);
			}
			problem = false;

			chrome = (chromosome<10)?"0"+chromosome:""+chromosome;
			randomKeys = random(IDs.size());
			if (!(new File("chrom"+chrome)).exists()) {
				(new File("chrom"+chrome)).mkdir();
			}

			reader = new BufferedReader(new FileReader("re_chrom"+chrome+".pre"));
			last = "";
			data = "";
			while (reader.ready()) {
				temp = reader.readLine();
				st = new StringTokenizer(temp);
				next = st.nextToken();
				if (!next.equals(last)) {
					if (!last.equals("")) {
						hash.put(last, data);
					}
					data = "";
				}
				data += temp+"\n";
				last = next;
			}
			hash.put(last, data);
			reader.close();

			increment = 0;
			prev = -99999.777;
			for (int i = 0; i<=IDs.size(); i++) {
				if (i!=IDs.size()) {
					if (direction.equals("a")) {
						trav = covars[keys[i]];
					} else if (direction.equals("d")) {
						trav = covars[keys[IDs.size()-i-1]];
					} else {
						trav = 0;
						System.err.println("Error: direction must either be \"a\" or \"d\"");
						System.exit(1);
					}
				} else {
					trav = -99999.777;
				}
				if (trav!=prev&&prev!=-99999.777) {
					increment++;
					writer = new PrintWriter(new FileWriter("chrom"+chrome+"/re_chrom"+chrome+"-"+direction+ext.formNum(increment, 4)+".pre"));
					for (int j = 0; j<i; j++) {
						if (hash.containsKey(IDs.elementAt(randomKeys[j]))) {
							writer.print(hash.get(IDs.elementAt(randomKeys[j])));
						} else {
							System.err.println("Error: could not find family "+IDs.elementAt(keys[j]));
						}
					}
					writer.close();
				}
				prev = trav;
			}
			numOfIncrements = increment;

			writer = new PrintWriter(new FileWriter("batch"));

			writer.println("#/bin/sh");
			writer.println();
			if ((new File("done"+chrome)).exists()) {
				writer.println("rm -r done"+chrome);
			}
			writer.println("cd chrom"+chrome);
			writer.println("cp ../map"+chrome+".dat .");
			writer.println("java -classpath /home/npankrat/" + common.PSF.Java.GENVISIS + " park.bat.dat2loc map"+chrome+".dat");
			for (int i = 1; i<=numOfIncrements; i++) {
				writer.println("echo -e \""+((chromosome==23)?"sex on\\n":"")+"pairs\\n3\\nload map"+chrome+".loc\\nprep re_chrom"+chrome+"-"+direction+ext.formNum(i, 4)+".pre\\nn\\nscan\\nestimate\\n"+((chromosome==23)?"":"n\\n")+"chrom"+chrome+"-"+direction+ext.formNum(i, 4)+".out\\nchrom"+chrome+"-"+direction+ext.formNum(i, 4)+"-share.ps\\nchrom"+chrome+"-"+direction+ext.formNum(i, 4)+"-mls.ps\\nquit\\n\" | /software/bin/sibs > /dev/null");
				writer.println();
			}
			writer.println("rm *.ps");
			writer.println("rm *.pre");
			writer.println("cd ..");
			writer.println("mv chrom"+chrome+" done"+chrome);
			writer.close();
			try {
				Process process = null;
				Runtime runtime = Runtime.getRuntime();
				process = runtime.exec("chmod +x batch");
				process.waitFor();
				Process process2 = null;
				Runtime runtime2 = Runtime.getRuntime();
				process2 = runtime2.exec("./batch");
				process2.waitFor();
			} catch (Exception e) {
				e.printStackTrace();
			}

			increment = 1;
			maxIncrement = -1;
			maxLod = 0;
			maxPos = -1;
			for (increment = 1; increment<=numOfIncrements; increment++) {
				try {
					reader = new BufferedReader(new FileReader("done"+chrome+"/chrom"+chrome+"-"+direction+ext.formNum(increment, 4)+".out"));
					pos = 0;

					if (chromosome<23) {
						reader.readLine();
						temp = reader.readLine();
						while (!temp.equals("")) {
							st = new StringTokenizer(temp);
							st.nextToken();
							st.nextToken();
							st.nextToken();
							st.nextToken();
							lod = Double.valueOf(st.nextToken()).doubleValue();
							if (lod>maxLod) {
								maxLod = lod;
								maxPos = pos;
								maxIncrement = increment;
							}

							temp = reader.readLine();
							pos++;
						}
					} else {
						do {
							temp = reader.readLine();
						} while (!temp.startsWith("Total"));
						temp = reader.readLine();
						while (!temp.equals("")) {
							st = new StringTokenizer(temp);
							st.nextToken();
							lod = Double.valueOf(st.nextToken()).doubleValue();
							if (lod>maxLod) {
								maxLod = lod;
								maxPos = pos;
								maxIncrement = increment;
							}

							temp = reader.readLine();
							pos++;
						}
					}

					reader.close();
				} catch (Exception e) {
					System.err.println("Ending with a total of "+(increment-1)+" increments.");
					e.printStackTrace();
					problem = true;
					break;
				}
			}

			if (!problem) {
				System.out.println(perm+"\t"+maxIncrement+"\t"+maxLod+"\t"+maxPos);
			}
		}

	}

	public int[] random(int size) throws IOException {
		IntVector source = new IntVector();
		int[] keys = new int[size];
		int temp;

		for (int i = 0; i<size; i++) {
			source.add(i);
		}

		temp = (int)(Math.random()*size);
		keys[0] = source.elementAt(temp);
		source.removeElement(temp);
		for (int i = size-1; i>0; i = i-1) {
			temp = source.elementAt((int)(Math.random()*i));
			keys[i] = temp;
			source.removeElement(temp);
		}

		return keys;
	}

	public static void main(String[] args) {
		if (args.length<4||args.length>5) {
			System.out.println("Expecting 4-5 arguments: filename, chromosome number, number of permutations, direction (a or d), and (optional) plug filename");
		} else {
			try {
				if (args.length==4) {
					new permOSA(args[0], Integer.valueOf(args[1]).intValue(), Integer.valueOf(args[2]).intValue(), args[3]);
				} else if (args.length==5) {
					new permOSA(args[0], Integer.valueOf(args[1]).intValue(), Integer.valueOf(args[2]).intValue(), args[3], args[4]);

				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
