package park;

import java.io.*;
import java.util.*;

import common.*;

import park.tools;

public class ParseRawSNPsOld48 {
	public static final String SNP_POSITIONS = "snp_positions.dat";

	public static void parse(String outfile) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		Vector<String[]> v = new Vector<String[]>();
		int count = 0, startpoint = 3, endpoint = -1, tally = 0;
		String[] header = null, famidPair;
		int[] keys, indIDs, snpPositions, snpOrder;
		SNP[] snps;

		while (new File("SNP"+(++count)+".dat").exists()) {
			try {
				reader = new BufferedReader(new FileReader("SNP"+count+".dat"));
				line = reader.readLine().split("\t", -1);
				if (count==1) {
					header = line;
					if (header[0].indexOf("SNPlex")==-1||header[1].indexOf("DNA")==-1||header[2].indexOf("Family")==-1) {
						System.err.println("Warning - first three headers are different than last time; might want to reconfigure");
					}
					for (int i = 0; i<header.length; i++) {
						if (header[i].toLowerCase().indexOf("failed")>=0) {
							endpoint = i;
						}
					}
					if (endpoint==-1) {
						endpoint = header.length;
					}
				} else {
					for (int i = 0; i<(header.length<line.length?header.length:line.length); i++) {
						if (!header[i].equals(line[i])) {
							System.err.println("Error - inconsistent header in file: "+"SNP"+count+".dat");
							System.err.println("        found '"+line[i]+"', expecting '"+header[i]+"'");
						}
					}
				}
				while (reader.ready()) {
					line = reader.readLine().split("\t", -1);
					if (!line[0].equals("")) {
						if (line[0].length()>2||line[0].equals(" ")) {
							System.err.println("errrr: "+line[0]);
						}
						v.add(line);
					} else {
						while (reader.ready()) {
							reader.readLine();
						}
						System.out.println("Parsed SNP"+count+".dat, adding "+(v.size()-tally));
						tally = v.size();
					}
				}
				reader.close();
			} catch (IOException ioe) {
				System.err.println("Error parsing "+"SNP"+count+".dat"+"");
				System.exit(3);
			}
		}
		System.out.println("    for a total of "+v.size()+" individuals");

		try {
			reader = new BufferedReader(new FileReader(SNP_POSITIONS));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				hash.put(line[0], line[1]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error - could not find "+SNP_POSITIONS+" in current directory");
			System.exit(2);
		} catch (IOException ioe) {
			System.err.println("Error parsing "+SNP_POSITIONS+"");
			System.exit(3);
		}

		if ((count = endpoint-startpoint)%2==1) {
			System.err.println("Warning - expecting even number of columns between FamID and #FailedSNPs");
			System.exit(1);
		}
		snps = new SNP[count /= 2];
		snpPositions = new int[snps.length];
		for (int i = 0; i<count; i++) {
			snps[i] = new SNP(header[startpoint+i*2]);
			if (hash.containsKey(snps[i].getName())) {
				snpPositions[i] = Integer.parseInt(hash.get(snps[i].getName()));
			} else {
				System.err.println("Error - marker '"+snps[i].getName()+"' was not listed in "+SNP_POSITIONS);
				System.exit(1);
			}
		}
		snpOrder = Sort.quicksort(snpPositions);

		indIDs = new int[v.size()];

		for (int i = 0; i<v.size(); i++) {
			line = v.elementAt(i);
			if (line[2].startsWith("*")) {
				line[2] = line[2].substring(1);
			}
			if (line[1].equals("2002PD0910")&&line[2].equals("7051-10")) {
				line[2] = "70541-10";
			}
			if (line[1].equals("2002PD0959")&&line[2].equals("7095-18")) {
				line[2] = "70395-18";
			}

			line[0] = line[1];
			famidPair = tools.getFamID(line[2]);
			line[1] = famidPair[0];
			line[2] = famidPair[1];
			indIDs[i] = Integer.parseInt(famidPair[0])*1000+Integer.parseInt(famidPair[1]);
			for (int j = 0; j<snps.length; j++) {
				snps[j].update(line[3+j*2+0], line[3+j*2+1]);
			}
		}

		writer = new PrintWriter(new FileWriter(outfile));
		writer.println("placeholder line");
		writer.print("DNA\tFamNo\tIndNo");
		for (int i = 0; i<snps.length; i++) {
			writer.print("\t"+snps[snpOrder[i]].getName()+"\t");
		}
		writer.println();

		keys = Sort.quicksort(indIDs);
		for (int i = 0; i<indIDs.length; i++) {
			line = v.elementAt(keys[i]);
			writer.print(line[0]+"\t"+line[1]+"\t"+line[2]);
			for (int j = 0; j<snps.length; j++) {
				for (int k = 0; k<2; k++) {
					writer.print("\t"+snps[snpOrder[j]].recode(line[3+snpOrder[j]*2+k]));
				}
			}
			writer.println();
		}
		writer.close();

		writer = new PrintWriter(new FileWriter(outfile+"-summary.xls"));
		for (int i = 0; i<snps.length; i++) {
			writer.println(snps[i]);
		}
		writer.close();

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String filename = "newSNPs.dat";
		String filename = "nextRound.dat";

		String usage = "\n"+"park.parseRawSNPs requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default)\n"+"";

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
			parse(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

class SNP {
	private String name;

	private String[] values;

	private int[] counts;

	public SNP(String name) {
		this.name = name;
		values = new String[] {"0"};
		counts = new int[] {0};
	}

	public void update(String allele1, String allele2) {
		String value;
		int index;

		for (int i = 0; i<2; i++) {
			value = i==0?allele1:allele2;
			index = ext.indexOfStr(value, values);
			if (index==-1) {
				values = Array.addStrToArray(value, values);
				counts = Array.addIntToArray(1, counts);
			} else {
				counts[index]++;
			}
		}
	}

	public String getName() {
		return name;
	}

	public String recode(String allele) {
		if (allele.equals("0")) {
			return "0";
		}
		if (allele.equals(values[1])) {
			return counts[1]<counts[2]?"1":"2";
		} else if (allele.equals(values[2])) {
			return counts[1]>counts[2]?"1":"2";
		} else {
			System.err.println("Error - allele '"+allele+"' is neither "+values[1]+" nor "+values[2]+" for marker "+name);
			return "0";
		}
	}

	public String toString() {
		int minor, major;

		while (values.length<3) {
			values = Array.addStrToArray("?", values);
			counts = Array.addIntToArray(0, counts);
		}
		minor = counts[1]<counts[2]?1:2;
		major = counts[1]<counts[2]?2:1;

		return name+"\t"+values[minor]+"\t"+ext.formDeci(ext.divide(counts[minor], counts[minor]+counts[major]), 4, true)+"\t"+values[major]+"\t"+ext.formDeci(ext.divide(counts[major], counts[minor]+counts[major]), 4, true)+"\t"+"Failed\t"+counts[0]+"\t"+ext.formDeci(ext.divide(counts[0], counts[0]+counts[1]+counts[2]), 4, true);
	}
}
