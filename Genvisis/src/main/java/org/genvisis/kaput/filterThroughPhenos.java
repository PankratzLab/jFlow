package org.genvisis.kaput;

import java.io.*;
import java.util.*;

public class filterThroughPhenos {
	public filterThroughPhenos(String structfile, String phenoFile, int affectedThreshold, int unaffectedThreshold, String unknownValue) throws IOException {
		BufferedReader reader = null;
		PrintWriter conAff = null, conUnaff = null, discord = null, unused = null, diskey = null;
		String temp, ID, trav, prev = "", rents = "shtoops";
		StringTokenizer st;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		Vector<String> affs = new Vector<String>();
		Vector<String> unfs = new Vector<String>();
		Vector<String> splats = new Vector<String>();
		Vector<String> disc = new Vector<String>();
		int pheno;

		reader = new BufferedReader(new FileReader(phenoFile));
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			ID = st.nextToken()+":"+st.nextToken();
			hash.put(ID, st.nextToken());
		}
		reader.close();

		reader = new BufferedReader(new FileReader(structfile));
		conAff = new PrintWriter(new FileWriter("concordant affected.out"));
		conUnaff = new PrintWriter(new FileWriter("concordant unaffected.out"));
		discord = new PrintWriter(new FileWriter("discordant.out"));
		unused = new PrintWriter(new FileWriter("unused individuals.out"));
		diskey = new PrintWriter(new FileWriter("diskey.dat"));

		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			trav = st.nextToken();
			if (trav.equals(prev)) {
				ID = st.nextToken();
				while (Integer.valueOf(ID).intValue()>100) {
					rents += "\n"+trav+"\t"+ID+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+unknownValue;
					st.nextToken();
					while (st.hasMoreTokens()) {
						rents += "\t"+st.nextToken();
					}
					st = new StringTokenizer(reader.readLine());
					st.nextToken();
					ID = st.nextToken();
				}
				if (!hash.containsKey(trav+":"+ID)) {
					temp = trav+"\t"+ID+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+unknownValue;
					st.nextToken();
					while (st.hasMoreTokens()) {
						temp += "\t"+st.nextToken();
					}
					splats.add(temp+" missing phenotype");
				} else {
					pheno = Integer.valueOf(hash.get(trav+":"+ID)).intValue();
					if (pheno>=affectedThreshold) {
						temp = trav+"\t"+ID+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+2;
						st.nextToken();
						while (st.hasMoreTokens()) {
							temp += "\t"+st.nextToken();
						}
						affs.add(temp);
						disc.add(temp);
					} else if (pheno>=unaffectedThreshold) {
						temp = trav+"\t"+ID+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+2;
						st.nextToken();
						while (st.hasMoreTokens()) {
							temp += "\t"+st.nextToken();
						}
						unfs.add(temp);
						disc.add(temp);
					} else if (unknownValue.equals(pheno+"")) {
						temp = trav+"\t"+ID+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+0;
						st.nextToken();
						while (st.hasMoreTokens()) {
							temp += "\t"+st.nextToken();
						}
						splats.add(temp+" crappy phenotype");
					} else {
						System.err.println("Error: Individual "+trav+"-"+ID+" has illegal phenotype: "+pheno);
					}
				}
			} else {
				if (affs.size()>=2) {
					conAff.println(rents);
					for (int i = 0; i<affs.size(); i++) {
						conAff.println(affs.elementAt(i));
					}
				}
				if (unfs.size()>=2) {
					conUnaff.println(rents);
					for (int i = 0; i<unfs.size(); i++) {
						conUnaff.println(unfs.elementAt(i));
					}
				}
				if (affs.size()>=1&&unfs.size()>=1) {
					discord.println(rents);
					for (int i = 0; i<disc.size(); i++) {
						discord.println(disc.elementAt(i));
					}
					String aff, unaff;

					for (int i = 0; i<affs.size(); i++) {
						st = new StringTokenizer(affs.elementAt(i));
						st.nextToken();
						aff = st.nextToken();
						for (int j = 0; j<unfs.size(); j++) {
							st = new StringTokenizer(unfs.elementAt(j));
							st.nextToken();
							unaff = st.nextToken();
							diskey.print(prev+"\t");
							if (Integer.valueOf(aff).intValue()<Integer.valueOf(unaff).intValue()) {
								for (int k = 0; k<disc.size(); k++) {
									st = new StringTokenizer(disc.elementAt(k));
									st.nextToken();
									if (st.nextToken().equals(aff)) {
										diskey.print((k+1));
									}
								}
								diskey.print("-");
								for (int k = 0; k<disc.size(); k++) {
									st = new StringTokenizer(disc.elementAt(k));
									st.nextToken();
									if (st.nextToken().equals(unaff)) {
										diskey.print((k+1));
									}
								}
							} else {
								for (int k = 0; k<disc.size(); k++) {
									st = new StringTokenizer(disc.elementAt(k));
									st.nextToken();
									if (st.nextToken().equals(unaff)) {
										diskey.print((k+1));
									}
								}
								diskey.print("-");
								for (int k = 0; k<disc.size(); k++) {
									st = new StringTokenizer(disc.elementAt(k));
									st.nextToken();
									if (st.nextToken().equals(aff)) {
										diskey.print((k+1));
									}
								}
							}
							diskey.println();
						}
					}
				}
				for (int i = 0; i<splats.size(); i++) {
					unused.println(splats.elementAt(i));
				}
				if (affs.size()<2&&unfs.size()<2&&(affs.size()==0||unfs.size()==0)) {
					for (int i = 0; i<affs.size(); i++) {
						unused.println(affs.elementAt(i)+" not enough pairs");
					}
					for (int i = 0; i<unfs.size(); i++) {
						unused.println(unfs.elementAt(i)+" not enough pairs");
					}
				}

				affs.removeAllElements();
				unfs.removeAllElements();
				splats.removeAllElements();
				disc.removeAllElements();
				prev = trav;
				rents = trav+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+unknownValue;
				st.nextToken();
				while (st.hasMoreTokens()) {
					rents += "\t"+st.nextToken();
				}
				st = new StringTokenizer(reader.readLine());
				rents += "\n"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken()+"\t"+unknownValue;
				st.nextToken();
				while (st.hasMoreTokens()) {
					rents += "\t"+st.nextToken();
				}
			}
		}

		reader.close();
		conAff.close();
		conUnaff.close();
		discord.close();
		diskey.close();
		unused.close();
	}

	public static void main(String[] args) throws IOException {
		if (args.length!=5) {
			System.out.println("Error: Expecting 5 arguments -  struct file, phenotype file (3 columns: FAM, ID and pheno), affected threshold (greater than or equal to),  unaffected threshold (greater than or equal to), and unknown value.");
		} else {
			try {
				new filterThroughPhenos(args[0], args[1], Integer.valueOf(args[2]).intValue(), Integer.valueOf(args[3]).intValue(), args[4]);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
