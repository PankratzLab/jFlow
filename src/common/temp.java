package common;

import java.io.*;
//import java.util.*;

import stats.Maths;

public class temp {

	public temp(String[] args) throws IOException {
		// BufferedReader reader = null;
		PrintWriter writer = null;
		// String temp, filename;
		// Hashtable hash = new Hashtable();
		// Vector<String> v = new Vector<String>();
		// String[] line, header;
		// double d;
		int count;

		// Bulls-eye leach (GND)
		String dir = "2007/09Jesse";
		String trav = dir.substring(7);
		count = 1;
		int fails = 0;
		new File(dir).mkdirs();
		while (fails<2) {
			if (Internat.downloadFile("http://www.bullz-eye.com/gnd/"+dir+"/"+trav+"-"+((count<10?"0":"")+count)+".jpg", dir+"/"+trav+"-"+((count<10?"0":"")+count)+".jpg")) {
				fails = 0;
			} else {
				fails++;
			}
			count++;
		}

		System.exit(1);

		System.out.println(Maths.nCr(3, 2));

		System.exit(1);

		CountVector cv = new CountVector();
		for (int i = 0; i<10000; i++) {
			count = 0;
			do {
				count++;
			} while (Math.random()<0.9998);
			cv.add(count+"");
		}
		String[] values = cv.getValues();
		int[] countz = cv.getCounts();
		writer = new PrintWriter(new FileWriter("prob-1.0"));
		for (int i = 0; i<values.length; i++) {
			writer.println(values[i]+"\t"+countz[i]);
		}
		writer.close();

		System.exit(1);

		// String db = "C:\\Download\\snp126";
		// String dir = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\PD-singleton\\";
		//
		// hash = HashVec.loadFileToHashString(dir+"poorly_localized.txt",
		// false);
		//
		// System.out.println(ext.getTime());
		// System.out.println("Searching database...");
		// try {
		// reader = new BufferedReader(new FileReader(db));
		// writer = new PrintWriter(new FileWriter(dir+"poor.out"));
		// reader.readLine();
		// while (reader.ready()) {
		// temp = reader.readLine();
		// line = temp.split("[\\s]+");
		// if (hash.containsKey(line[4])) {
		// writer.println(temp);
		// writer.flush();
		// }
		// }
		// reader.close();
		// writer.close();
		// } catch (FileNotFoundException fnfe) {
		// System.err.println("Error: file \""+db+"\" not found in current
		// directory");
		// System.exit(1);
		// } catch (IOException ioe) {
		// System.err.println("Error reading file \""+db+"\"");
		// System.exit(2);
		// }
		// System.out.println(ext.getTime());
		//
		//
		// System.exit(1);
		//
		// IntVector dv = new IntVector();
		// int[] woah = new int[] {1, 4, 2, 0, 1, -1, -5, 10};
		// for (int i = 0; i<woah.length; i++) {
		// System.out.println(dv.add(woah[i], true));
		// System.out.println(Array.toStr(dv.toArray()));
		// }
		//
		//
		// System.exit(1);
		//
		// filename = "C:\\Download\\Hs.data";
		// try {
		// reader = new BufferedReader(new FileReader(filename));
		// writer = new PrintWriter(new FileWriter(filename+".out"));
		// while (reader.ready()) {
		// temp = reader.readLine();
		// if (temp.indexOf("NM_001040139") >= 0) {
		// writer.println(temp);
		// }
		// }
		// writer.close();
		// reader.close();
		// } catch (FileNotFoundException fnfe) {
		// System.err.println("Error: file \""+filename+"\" not found in current
		// directory");
		// System.exit(1);
		// } catch (IOException ioe) {
		// System.err.println("Error reading file \""+filename+"\"");
		// System.exit(2);
		// }
		//
		// System.exit(1);
		//
		// filename = "C:\\Download\\seq_gene.md";
		// try {
		// reader = new BufferedReader(new
		// FileReader("C:\\Download\\temp.txt"));
		// while (reader.ready()) {
		// v.add(reader.readLine());
		// }
		// reader.close();
		// line = Array.toStringArray(v);
		//
		// reader = new BufferedReader(new FileReader(filename));
		// writer = new PrintWriter(new FileWriter(filename+".out"));
		// while (reader.ready()) {
		// temp = reader.readLine();
		// for (int i = 0; i<line.length; i++) {
		// if (temp.indexOf("GeneID:"+line[i]+"\t") >= 0) {
		// writer.println(temp);
		// }
		// }
		//
		// }
		// writer.close();
		// reader.close();
		// } catch (FileNotFoundException fnfe) {
		// System.err.println("Error: file \""+filename+"\" not found in current
		// directory");
		// System.exit(1);
		// } catch (IOException ioe) {
		// System.err.println("Error reading file \""+filename+"\"");
		// System.exit(2);
		// }
		//
		//
		//
		// System.exit(1);
		//
		// double[] huh = new double[] {1.1, 1.4, 1.5, 1.7, 1.9};
		// System.out.println(Array.toStr(huh));
		// System.out.println(java.util.Arrays.binarySearch(huh, 1.4));
		// System.out.println(-1*java.util.Arrays.binarySearch(huh, 1.6)-1);
		//
		// System.exit(1);
		//
		// filename = "GBA_Trend_master.dat";
		//
		// try {
		// reader = new BufferedReader(new FileReader(filename));
		// header = reader.readLine().split("[\\s]+");
		// reader.close();
		// writer = new PrintWriter(new FileWriter("permuted.xls"));
		// for (int i = 2; i<header.length; i++) {
		// reader = new BufferedReader(new
		// FileReader(header[i]+".dat-summary.out"));
		// reader.readLine();
		// reader.readLine();
		// reader.readLine();
		// line = reader.readLine().split("[\\s]+");
		// writer.println(header[i]+"\t"+line[1]+"\t"+line[4]+"\t"+line[5]+"\t"+line[8]+"\t"+line[10]);
		// reader.close();
		// }
		// writer.close();
		//
		// } catch (FileNotFoundException fnfe) {
		// System.err.println("Error: file \""+filename+"\" not found in current
		// directory");
		// System.exit(1);
		// } catch (IOException ioe) {
		// System.err.println("Error reading file \""+filename+"\"");
		// System.exit(2);
		// }
		//
		// System.exit(1);
		//
		// filename = "GBA_Trend_master.dat";
		//
		// try {
		// reader = new BufferedReader(new FileReader(filename));
		// header = reader.readLine().split("[\\s]+");
		// reader.close();
		// for (int i = 2; i<header.length; i++) {
		// reader = new BufferedReader(new FileReader(filename));
		// writer = new PrintWriter(new FileWriter(header[i]+".dat"));
		// System.out.println(" String filename = \""+header[i]+".dat"+"\";");
		// while (reader.ready()) {
		// line = reader.readLine().split("[\\s]+");
		// writer.println(line[0]+"\t"+line[i]+"\t"+line[1]);
		// }
		// reader.close();
		// writer.close();
		// }
		//
		// } catch (FileNotFoundException fnfe) {
		// System.err.println("Error: file \""+filename+"\" not found in current
		// directory");
		// System.exit(1);
		// } catch (IOException ioe) {
		// System.err.println("Error reading file \""+filename+"\"");
		// System.exit(2);
		// }
		//
		// System.exit(1);
		//
		// line = new String[] {"One", "Two", "Three", "Four"};
		// System.out.println("This is how we do it: "+ext.listWithCommas(line,
		// false)+".");
		// line = new String[] {"One", "Two", "Three"};
		// System.out.println("This is how we do it: "+ext.listWithCommas(line,
		// false)+".");
		// line = new String[] {"One", "Two"};
		// System.out.println("This is how we do it: "+ext.listWithCommas(line,
		// false)+".");
		// line = new String[] {"One"};
		// System.out.println("This is how we do it: "+ext.listWithCommas(line,
		// false)+".");
		// line = new String[] {};
		// System.out.println("This is how we do it: "+ext.listWithCommas(line,
		// false)+".");
		//
		// System.exit(1);
		//
		// double[] freqs = {0.286703601, 0.653739612, 0.059556787};
		// int[] counts;
		// int pick1, pick2;
		// double rand, chi;
		//
		// writer = new PrintWriter(new FileWriter("HWE_reps.xls"));
		// writer.println("pp\tpq\tpr\tqq\tqr\trr\tchi\t1df\t2df\t3df\t4df\t5df\6df");
		// for (int i = 0; i < 65520; i++) {
		// counts = new int[6];
		// for (int j = 0; j < 361; j++) {
		// rand = Math.random();
		// if (rand < freqs[0]) {
		// pick1 = 0;
		// } else if (rand > freqs[0] && rand < freqs[0]+freqs[1]) {
		// pick1 = 1;
		// } else {
		// pick1 = 2;
		// }
		// rand = Math.random();
		// if (rand < freqs[0]) {
		// pick2 = 0;
		// } else if (rand > freqs[0] && rand < freqs[0]+freqs[1]) {
		// pick2 = 1;
		// } else {
		// pick2 = 2;
		// }
		// counts[(pick1==0||pick2==0?0:1)+pick1+pick2]++;
		// }
		// chi = AlleleFreq.HWE(counts[0], counts[1], counts[2], counts[3],
		// counts[4], counts[5]);
		// writer.println(Array.toStr(counts)+"\t"+chi+"\t"+ProbDist.ChiDist(chi,
		// 1)+"\t"+ProbDist.ChiDist(chi, 2)+"\t"+ProbDist.ChiDist(chi,
		// 3)+"\t"+ProbDist.ChiDist(chi, 4)+"\t"+ProbDist.ChiDist(chi,
		// 5)+"\t"+ProbDist.ChiDist(chi, 6));
		// }
		// writer.close();
		//
		//
		//
		// System.exit(1);
		//
		// double[][] data;
		// double[] ci;
		//
		// data = new double[][] {{ 31, 68 },
		// { 95, 1448 }};
		// System.out.println("Aharon-Peretz 2004 - Ashkenazi GBA");
		// ci = ContingencyTable.oddsRatioCI(data);
		// System.out.println("The odds ratio is
		// "+ext.formDeci(ContingencyTable.oddsRatio(data), 4, true) +"
		// ("+ext.formDeci(ci[0], 4, true)+", "+ext.formDeci(ci[1], 4,
		// true)+")");
		// System.out.println("The Likelhood-Ratio chi-square statistic is
		// "+ext.formDeci(ContingencyTable.likelihoodRatioStatistic(data), 4,
		// true) +"
		// (p="+ProbDist.ChiDist(ContingencyTable.likelihoodRatioStatistic(data),
		// (data.length-1)*(data[0].length-1))+")");
		// System.out.println();
		//
		// data = new double[][] {{ 5, 83 },
		// { 1, 121 }};
		// System.out.println("Sato 2004 - Canadian");
		// ci = ContingencyTable.oddsRatioCI(data);
		// System.out.println("The odds ratio is
		// "+ext.formDeci(ContingencyTable.oddsRatio(data), 4, true) +"
		// ("+ext.formDeci(ci[0], 4, true)+", "+ext.formDeci(ci[1], 4,
		// true)+")");
		// System.out.println("The Likelhood-Ratio chi-square statistic is
		// "+ext.formDeci(ContingencyTable.likelihoodRatioStatistic(data), 4,
		// true) +"
		// (p="+ProbDist.ChiDist(ContingencyTable.likelihoodRatioStatistic(data),
		// (data.length-1)*(data[0].length-1))+")");
		// System.out.println();
		//
		// data = new double[][] {{ 12, 45 },
		// { 2, 42 }};
		// System.out.println("Lwin 2004 - US");
		// ci = ContingencyTable.oddsRatioCI(data);
		// System.out.println("The odds ratio is
		// "+ext.formDeci(ContingencyTable.oddsRatio(data), 4, true) +"
		// ("+ext.formDeci(ci[0], 4, true)+", "+ext.formDeci(ci[1], 4,
		// true)+")");
		// System.out.println("The Likelhood-Ratio chi-square statistic is
		// "+ext.formDeci(ContingencyTable.likelihoodRatioStatistic(data), 4,
		// true) +"
		// (p="+ProbDist.ChiDist(ContingencyTable.likelihoodRatioStatistic(data),
		// (data.length-1)*(data[0].length-1))+")");
		// System.out.println();
		//
		// data = new double[][] {{ 4, 29 },
		// { 1, 30 }};
		// System.out.println("Eblan 2006 - Venezualan");
		// ci = ContingencyTable.oddsRatioCI(data);
		// System.out.println("The odds ratio is
		// "+ext.formDeci(ContingencyTable.oddsRatio(data), 4, true) +"
		// ("+ext.formDeci(ci[0], 4, true)+", "+ext.formDeci(ci[1], 4,
		// true)+")");
		// System.out.println("The Likelhood-Ratio chi-square statistic is
		// "+ext.formDeci(ContingencyTable.likelihoodRatioStatistic(data), 4,
		// true) +"
		// (p="+ProbDist.ChiDist(ContingencyTable.likelihoodRatioStatistic(data),
		// (data.length-1)*(data[0].length-1))+")");
		// System.out.println();
		//
		// data = new double[][] {{ 4, 88 },
		// { 1, 91 }};
		// System.out.println("Ziegler 2007 - Taiwan");
		// ci = ContingencyTable.oddsRatioCI(data);
		// System.out.println("The odds ratio is
		// "+ext.formDeci(ContingencyTable.oddsRatio(data), 4, true) +"
		// ("+ext.formDeci(ci[0], 4, true)+", "+ext.formDeci(ci[1], 4,
		// true)+")");
		// System.out.println("The Likelhood-Ratio chi-square statistic is
		// "+ext.formDeci(ContingencyTable.likelihoodRatioStatistic(data), 4,
		// true) +"
		// (p="+ProbDist.ChiDist(ContingencyTable.likelihoodRatioStatistic(data),
		// (data.length-1)*(data[0].length-1))+")");
		// System.out.println();
		//
		// data = new double[][] {{ 7, 304 },
		// { 8, 466 }};
		// System.out.println("Toft 2006 - Norwegian");
		// ci = ContingencyTable.oddsRatioCI(data);
		// System.out.println("The odds ratio is
		// "+ext.formDeci(ContingencyTable.oddsRatio(data), 4, true) +"
		// ("+ext.formDeci(ci[0], 4, true)+", "+ext.formDeci(ci[1], 4,
		// true)+")");
		// System.out.println("The Likelhood-Ratio chi-square statistic is
		// "+ext.formDeci(ContingencyTable.likelihoodRatioStatistic(data), 4,
		// true) +"
		// (p="+ProbDist.ChiDist(ContingencyTable.likelihoodRatioStatistic(data),
		// (data.length-1)*(data[0].length-1))+")");
		// System.out.println();
		//
		// data = new double[][] {{ 63, 617 },
		// { 108, 2198 }};
		// System.out.println("Meta Analysis of all published studies");
		// ci = ContingencyTable.oddsRatioCI(data);
		// System.out.println("The odds ratio is
		// "+ext.formDeci(ContingencyTable.oddsRatio(data), 4, true) +"
		// ("+ext.formDeci(ci[0], 4, true)+", "+ext.formDeci(ci[1], 4,
		// true)+")");
		// System.out.println("The Likelhood-Ratio chi-square statistic is
		// "+ext.formDeci(ContingencyTable.likelihoodRatioStatistic(data), 4,
		// true) +"
		// (p="+ProbDist.ChiDist(ContingencyTable.likelihoodRatioStatistic(data),
		// (data.length-1)*(data[0].length-1))+")");
		// System.out.println();
		//
		//
		// data = new double[][] {{ 9, 759 },
		// { 3, 765 }};
		// System.out.println("Abou-Sleiman 2006 - PINK1");
		// ci = ContingencyTable.oddsRatioCI(data);
		// System.out.println("The odds ratio is
		// "+ext.formDeci(ContingencyTable.oddsRatio(data), 4, true) +"
		// ("+ext.formDeci(ci[0], 4, true)+", "+ext.formDeci(ci[1], 4,
		// true)+")");
		// System.out.println("The Likelhood-Ratio chi-square statistic is
		// "+ext.formDeci(ContingencyTable.likelihoodRatioStatistic(data), 4,
		// true) +"
		// (p="+ProbDist.ChiDist(ContingencyTable.likelihoodRatioStatistic(data),
		// (data.length-1)*(data[0].length-1))+")");
		// System.out.println();
		//
		// System.exit(1);
		//
		// try {
		// reader = new BufferedReader(new FileReader("freak2.dat"));
		// writer = new PrintWriter(new FileWriter("freak2.out"));
		// while (reader.ready()) {
		// temp = reader.readLine();
		// if (temp.indexOf("only") > 0) {
		// temp = reader.readLine();
		// }
		// temp = temp.substring(0, temp.length()-3);
		// line = temp.split("[\\s]+");
		// writer.println(line[line.length-1]+"\t"+reader.readLine());
		// }
		// reader.close();
		// writer.close();
		// } catch (FileNotFoundException fnfe) {
		// System.err.println("Error: file \""+"1_1_readin.txt"+"\" not found in
		// current directory");
		// System.exit(1);
		// } catch (IOException ioe) {
		// System.err.println("Error reading file \""+"1_1_readin.txt"+"\"");
		// System.exit(2);
		// }
		//
		// try {
		// reader = new BufferedReader(new FileReader("freak2.dat"));
		// writer = new PrintWriter(new FileWriter("freak2_errors.out"));
		// while (reader.ready()) {
		// temp = reader.readLine();
		// if (temp.indexOf("only") > 0) {
		// writer.println(temp);
		// }
		// }
		// reader.close();
		// writer.close();
		// } catch (FileNotFoundException fnfe) {
		// System.err.println("Error: file \""+"1_1_readin.txt"+"\" not found in
		// current directory");
		// System.exit(1);
		// } catch (IOException ioe) {
		// System.err.println("Error reading file \""+"1_1_readin.txt"+"\"");
		// System.exit(2);
		// }
		//
		//
		//
		// System.exit(1);
		//
		// // System.out.println("hi!");
		// // args = new String[] {"-pedfile", "rs2017143_001.pre", "-info",
		// "rs2017143_001.info", "-dprime", "-n"};
		// // edu.mit.wi.haploview.HVWrap.main(args);
		// // System.out.println("ho!");
		//
		// System.exit(1);
		//
		// filename = "Nichols_6Kv4_LinkageApril2007_FinalReport.txt";
		// String target = "rs748325";
		// try {
		// reader = new BufferedReader(new FileReader(filename));
		// writer = new PrintWriter(new FileWriter(target+".out"));
		// while (reader.ready()) {
		// temp = reader.readLine();
		// if (temp.indexOf(target) != -1) {
		// writer.println(temp);
		// }
		// }
		// reader.close();
		// writer.close();
		// } catch (FileNotFoundException fnfe) {
		// System.err.println("Error: file \""+filename+"\" not found in current
		// directory");
		// System.exit(1);
		// } catch (IOException ioe) {
		// System.err.println("Error reading file \""+filename+"\"");
		// System.exit(2);
		// }
		// System.exit(1);
		//
		// for (int i = 0; i<10; i++) {
		// System.out.println((int)(Math.random()*1000000));
		// }
		//
		//
		// System.exit(1);
		//
		// for (int i = 0; i<10; i++) {
		// new MarkerAnalysis("6K_dataset.dat", (i+1)+"_"+"1_rep.xls");
		// }
		//
		// System.exit(1);
		//
		// try {
		// reader = new BufferedReader(new FileReader("HWE.txt"));
		// writer = new PrintWriter(new FileWriter("HWE.xls"));
		// while (reader.ready()) {
		// temp = reader.readLine();
		// line = temp.split("[\\s]+");
		// writer.println(temp+"\t"+AlleleFreq.HWEsig(Double.parseDouble(line[0]),
		// Double.parseDouble(line[1]), Double.parseDouble(line[2])));
		// }
		// reader.close();
		// writer.close();
		// } catch (FileNotFoundException fnfe) {
		// System.err.println("Error: file \""+"HWE.txt"+"\" not found in
		// current directory");
		// System.exit(1);
		// } catch (IOException ioe) {
		// System.err.println("Error reading file \""+"HWE.txt"+"\"");
		// System.exit(2);
		// }
		//
		//
		// System.exit(1);
		//
		// // String filename = "Nichols_6Kv4_LinkageApril2007_FinalReport.txt";
		// // String target = "rs6907703";
		// // try {
		// // reader = new BufferedReader(new FileReader(filename));
		// // writer = new PrintWriter(new FileWriter(target+".out"));
		// // while (reader.ready()) {
		// // temp = reader.readLine();
		// // if (temp.indexOf(target) != -1) {
		// // writer.println(temp);
		// // }
		// // }
		// // reader.close();
		// // writer.close();
		// // } catch (FileNotFoundException fnfe) {
		// // System.err.println("Error: file \""+filename+"\" not found in
		// current directory");
		// // System.exit(1);
		// // } catch (IOException ioe) {
		// // System.err.println("Error reading file \""+filename+"\"");
		// // System.exit(2);
		// // }
		// // System.exit(1);
		//
		// int[] base = new int[] {1000, 2000, 3000};
		// int[] rep;
		// for (int i = 0; i<10; i++) {
		// rep = base.clone();
		// for (int j = 0; j<rep.length; j++) {
		// rep[j] += i;
		// }
		// System.out.println(Array.toStr(rep));
		// }
		//
		//
		// System.exit(1);
		//
		// for (int i = 0; i<255; i++) {
		// System.out.println(i+"\t"+(char)i);
		// }
		//
		//
		// System.exit(1);
		//
		// d = 1.0;
		// v.add(ext.formDeci(d, 5));
		// System.out.println(v.elementAt(0));
		// // int ii = Integer.parseInt(v.elementAt(0));
		//
		// System.exit(1);
		//
		// line = new String[] {"Alpha", "Squad", "Seven", "A", "Tek", "Jansen",
		// "Adventure"};
		// for (int i = 0; i < line.length; i++) {
		// HashVec.addIfAbsent(line[i], v, true);
		// }
		//
		// for (int i = 0; i < v.size(); i++) {
		// System.out.println(v.elementAt(i));
		// }
		// System.exit(1);
		//
		//
		// int[] keys;
		//
		// for (int j = 0; j < 20; j++) {
		// keys = Array.random(4);
		// for (int i = 0; i < keys.length; i++) {
		// System.out.print(keys[i] + "\t");
		// }
		// System.out.println("");
		// }
		//
		//
		// System.exit(1);
		//
		// String[] ya1 = new String[] {"Hello", "you"};
		// String[] ya2 = new String[2];
		// ya2[0] = "Hello";
		// ya2[1] = "you";
		//
		// if (ext.eqArrays(ya1, ya2)) {
		// System.out.println("worked");
		// } else {
		// System.out.println("please try again");
		// }

	}

	public static void main(String[] args) throws IOException {
		try {
			new temp(args);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

// test out perason correlation
// try {
// reader = new BufferedReader(new FileReader("pearson.txt"));
// while (reader.ready()) {
// line = reader.readLine().split("[\\s]+");
// v.add(new double[] {Double.parseDouble(line[0]),
// Double.parseDouble(line[1])});
// }
// reader.close();
// } catch (FileNotFoundException fnfe) {
// System.err.println("Error: file \"" + "pearson.txt" + "\" not found in
// current directory");
// System.exit(1);
// } catch (IOException ioe) {
// System.err.println("Error reading file \"" + "pearson.txt" + "\"");
// System.exit(2);
// }
// double[][] data = new double[2][v.size()];
// double[] ds;
// for (int i = 0; i < v.size(); i++) {
// ds = v.elementAt(i);
// data[0][i] = ds[0];
// data[1][i] = ds[1];
// }
//
// ds = Stats.PearsonCorrelation(data);
// System.out.println("r = "+ds[0]);
// System.out.println("sig = "+ds[1]);
//
