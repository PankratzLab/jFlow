// program case-control "global haplotype analysis" (chi sq)

package org.genvisis.assoc;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import org.genvisis.common.*;

public class haplorParser {
	public static String[] TARGET_HAPLOTYPES = {"111X", "222X", "X1X2", "X2X1"};

	public class HapFam {
		public String id;

		public Hashtable<String,String> posts;

		public Hashtable<String,HapInd> members;

		public String maxConfig;

		public double maxPost;

		public HapFam(String famID) throws IOException {
			id = famID;
			posts = new Hashtable<String,String>();
			members = new Hashtable<String,HapInd>();
			maxConfig = "-";
			maxPost = 0;
		}
	}

	public class HapInd {
		public String id;

		public Hashtable<String,String> posts;

		public String maxHap;

		public double maxPost;

		public HapInd(String indID) throws IOException {
			id = indID;
			posts = new Hashtable<String,String>();
			maxHap = "-";
			maxPost = 0;
		}
	}

	public haplorParser(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		PrintWriter[] writers;
		String[] line, marks;
		String temp, trav, prev, maxPair;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		int numMarkers, numConfigs, numInds, count;
		Vector<String> fams = new Vector<String>();
		Vector<String> missingdata = new Vector<String>();
		String[] haplotypes;
		double[] hapFreqs;
		int[] hapKeys;
		String[] indIDs;
		double post, maxPost;
		Pattern delimiters = Pattern.compile("[\\s(),]+");
		Pattern comma = Pattern.compile(",");
		IntVector typed;
		HapFam hFam;
		HapInd hInd;
		Enumeration<String> enumer;
		boolean compare = (new File("f3_haplos.out")).exists()&&(new File("f3_em.out")).exists();
		boolean done;
		int offset;
		Vector<String> compHaps = new Vector<String>();
		Vector<String> comparisons = new Vector<String>();

		writers = new PrintWriter[TARGET_HAPLOTYPES.length];
		for (int i = 0; i<TARGET_HAPLOTYPES.length; i++) {
			writers[i] = new PrintWriter(new FileWriter("BRI3-haplotype-"+TARGET_HAPLOTYPES[i]+"-cases.dat"));
			writers[i].println("UniqueID\tFamID\tIndID\thaplotype-"+TARGET_HAPLOTYPES[i]);
		}

		if (compare) {
			reader = new BufferedReader(new FileReader("f3_em.out"));
			for (int i = 0; i<5; i++) {
				reader.readLine();
			}
			count = 1;
			done = false;
			while (!done) {
				line = reader.readLine().split("[\\s]+");
				offset = (line.length>1&&line[0].equals(""))?1:0;
				if (line[0+offset].equals(count+"")) {
					temp = line[1+offset];
					for (int i = 2+offset; i<line.length-1; i++) {
						temp += line[i];
					}
					compHaps.add(temp);
					count++;
				} else {
					done = true;
				}
			}
			reader.close();

			if (!new File("f3_aff.dat").exists()) {
				System.err.println("Error - could not find "+"f3_aff.dat"+" in current directory");
				System.exit(2);
			}
			reader = new BufferedReader(new FileReader("f3_aff.dat"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				comparisons.add(line[0]);
			}
			reader.close();

			reader = new BufferedReader(new FileReader("f3_haplos.out"));
			count = 1;
			done = false;
			trav = temp = "error";
			prev = "";
			maxPair = "";
			post = maxPost = -1;
			while (!done) {
				do {
					line = reader.readLine().split("[\\s]+");
					offset = (line.length>1&&line[0].equals(""))?1:0;
					if (line.length>1&&line[1].equals("Number")) {
						done = true;
						trav = "";
					}
				} while (!done&&!line[0+offset].equals(count+""));
				if (!done) {
					trav = comparisons.elementAt(Integer.valueOf(line[1+offset]).intValue()-1);
					temp = compHaps.elementAt(Integer.valueOf(line[2+offset]).intValue()-1)+"\t"+compHaps.elementAt(Integer.valueOf(line[3+offset]).intValue()-1);
					post = Double.valueOf(line[4+offset]).doubleValue();
				}
				if (trav.equals(prev)) {
					if (post>maxPost) {
						maxPost = post;
						maxPair = temp;
					}
				} else {
					hash.put(prev, maxPair+"\t"+maxPost);
					maxPost = post;
					maxPair = temp;
					prev = trav;
				}
				count++;
			}
			reader.close();
		}

		if (!new File(filename).exists()) {
			System.err.println("Error - could not find "+filename+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(filename));

		do {
			temp = reader.readLine();
		} while (reader.ready()&&!temp.startsWith("The haplotypes obtained "));

		line = reader.readLine().split("[\\s]+");
		fams.add(line[0]);
		line = reader.readLine().split("[\\s]+");
		numMarkers = line.length-1;
		reader.readLine();
		reader.readLine();
		temp = reader.readLine();

		do {
			line = reader.readLine().split("[\\s]+");
			if (!fams.contains(line[0])&&line.length>1) {
				fams.add(line[0]);
			}
			reader.readLine();
			reader.readLine();
			reader.readLine();
			marks = reader.readLine().split("[\\s]+");
			count = 0;
			for (int i = 1; i<marks.length; i++) {
				if (!marks[i].equals("0")) {
					count++;
				}
			}
			if (count<(double)numMarkers&&line.length>1&&Integer.valueOf(line[1]).intValue()<100) {
				System.out.print(".");
				missingdata.add(line[0]+"-"+line[1]);
			}
		} while (reader.ready()&&line.length>1);
		System.out.println();
		System.out.println("Warning - "+missingdata.size()+" individuals were excluded due to missing data");

		System.out.println("Determining haplotypes for "+fams.size()+" families");

		do {
			temp = reader.readLine();
		} while (reader.ready()&&!temp.startsWith("The total of "));

		writer = new PrintWriter(new FileWriter("haplor-key.out"));
		line = temp.split("[\\s]+");
		haplotypes = new String[Integer.valueOf(line[3]).intValue()];
		hapFreqs = new double[haplotypes.length];
		for (int i = 0; i<haplotypes.length; i++) {
			line = reader.readLine().split("[\\s]+");
			temp = "";
			for (int j = 0; j<numMarkers; j++) {
				temp += line[3+j];
			}
			haplotypes[i] = temp;
			hapFreqs[i] = Double.valueOf(line[3+numMarkers]).doubleValue();
			writer.println((i+1)+"\t"+haplotypes[i]+"\t"+hapFreqs[i]);
		}
		hapKeys = Sort.quicksort(haplotypes);
		writer.println("\n");
		for (int i = 0; i<haplotypes.length; i++) {
			writer.println((i+1)+"\t"+haplotypes[hapKeys[i]]+"\t"+hapFreqs[hapKeys[i]]);
		}
		writer.close();

		do {
			temp = reader.readLine();
		} while (reader.ready()&&!temp.startsWith("All possible haplotype configurations are:"));
		reader.readLine();

		writer = new PrintWriter(new FileWriter("haplotypes.dat"));
		for (int i = 0; i<fams.size(); i++) {
			// for (int i = 0; i<2; i++) {
			line = reader.readLine().split("[\\s]+");
			if (!fams.elementAt(i).equals(line[3])) {
				System.err.println("Error - failed family line up. Expecting "+fams.elementAt(i)+" found "+line[3]+".");
			}
			hFam = new HapFam(line[3]);
			// hash.put(line[3], hFam);
			numConfigs = Integer.valueOf(reader.readLine().split("[\\s]+")[9]).intValue();
			numInds = Integer.valueOf(reader.readLine().split("[\\s]+")[8]).intValue();
			reader.readLine();
			reader.readLine();

			line = reader.readLine().split("[\\s]+");
			indIDs = new String[numInds];
			typed = new IntVector();
			for (int j = 0; j<numInds; j++) {
				indIDs[j] = line[1+j];
				if (Integer.valueOf(indIDs[j]).intValue()<100) {
					typed.add(j);
					hFam.members.put(indIDs[j], new HapInd(indIDs[j]));
				}
			}

			for (int j = 0; j<numConfigs; j++) {
				line = delimiters.split(reader.readLine());
				post = Double.valueOf(line[line.length-1]).doubleValue();
				temp = "";
				for (int k = 0; k<numInds; k++) {
					if (typed.contains(k)) {
						trav = line[2+k*2]+"//"+line[2+k*2+1];
						temp += (temp.equals("")?"":",")+trav;
						hInd = hFam.members.get(indIDs[k]);
						if (hInd.posts.containsKey(trav)) {
							hInd.posts.put(trav, (Double.valueOf(hInd.posts.get(trav)).doubleValue()+post)+"");
						} else {
							hInd.posts.put(trav, post+"");
						}
					}
				}
				if (hFam.posts.containsKey(temp)) {
					hFam.posts.put(temp, (Double.valueOf(hFam.posts.get(temp)).doubleValue()+post)+"");
				} else {
					hFam.posts.put(temp, post+"");
				}
			}

			enumer = hFam.posts.keys();
			while (enumer.hasMoreElements()) {
				temp = enumer.nextElement();
				post = Double.valueOf(hFam.posts.get(temp)).doubleValue();
				if (post>hFam.maxPost) {
					hFam.maxConfig = temp;
					hFam.maxPost = post;
				}
			}
			line = comma.split(hFam.maxConfig);
			for (int j = 0; j<numInds; j++) {
				if (typed.contains(j)) {
					hInd = (HapInd)hFam.members.get(indIDs[j]);
					enumer = hInd.posts.keys();
					while (enumer.hasMoreElements()) {
						temp = enumer.nextElement();
						post = Double.valueOf(hInd.posts.get(temp)).doubleValue();
						if (post>hInd.maxPost) {
							hInd.maxHap = temp;
							hInd.maxPost = post;
						}
					}
					if (missingdata.contains(hFam.id+"-"+hInd.id)) {
						writer.print(hFam.id+"\t"+hInd.id+"\t.\t.\t.\t.\t.\t.");
						writer.println(compare?"\t.\t.\t.":"");
					} else {
						// writer.println(hFam.id+"\t"+hInd.id+"\t"+line[typed.indexOf(j)]+"\t"+ext.formDeci(hFam.maxPost,
						// 3,
						// true)+"\t"+hInd.maxHap+"\t"+ext.formDeci(hInd.maxPost,
						// 3, true));
						writer.print(hFam.id+"\t"+hInd.id+"\t"+translateHap(line[typed.indexOf(j)], haplotypes)+"\t"+ext.formDeci(hFam.maxPost, 3, true)+"\t"+translateHap(hInd.maxHap, haplotypes)+"\t"+ext.formDeci(hInd.maxPost, 3, true));
						temp = (Integer.valueOf(hFam.id).intValue()*1000+Integer.valueOf(hInd.id).intValue())+"";
						writer.println(compare?(hash.containsKey(temp)?"\t"+hash.get(temp):"\t.\t.\t."):"");
					}
					for (int k = 0; k<TARGET_HAPLOTYPES.length; k++) {
						writers[k].print(hFam.id+ext.formNum(hInd.id, 3)+"\t"+hFam.id+"\t"+hInd.id);
						if (missingdata.contains(hFam.id+"-"+hInd.id)) {
							writers[k].println("\t0\t0");
						} else {
							marks = translateHap(line[typed.indexOf(j)], haplotypes).split("[\\s]+");
							writers[k].println("\t"+compHap(marks[0], TARGET_HAPLOTYPES[k])+"\t"+compHap(marks[1], TARGET_HAPLOTYPES[k]));
						}
					}

					hFam.members.put(indIDs[j], new HapInd(indIDs[j]));
				}
			}

			reader.readLine();
		}
		reader.close();
		writer.close();

		for (int i = 0; i<TARGET_HAPLOTYPES.length; i++) {
			writers[i].close();
		}

	}

	public String translateHap(String hap_pair, String[] haps) throws IOException {
		String[] line = hap_pair.split("[/]+");
		String[] pair = new String[line.length];
		int[] keys;

		for (int i = 0; i<line.length; i++) {
			pair[i] = haps[Integer.valueOf(line[i]).intValue()-1];
		}
		keys = Sort.quicksort(pair);

		return pair[keys[0]]+"\t"+pair[keys[1]];
	}

	public String compHap(String hap, String target) throws IOException {
		for (int i = 0; i<(hap.length()<target.length()?hap.length():target.length()); i++) {
			if (target.charAt(i)!='X'&&hap.charAt(i)!=target.charAt(i)) {
				return "1";
			}
		}
		return "2";
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String filename = "pd_first3_all.out";
		String filename = "trythis.ped.out";

		String usage = "\n"+"park.haplorParser requires 1 arguments:\n"+"   (1) the name of the file to parse (i.e. file="+filename+" (default))\n"+"";

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
		if (args.length==0) {
			System.out.println("Using defaults (file="+filename+")");
		}

		try {
			new haplorParser(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
