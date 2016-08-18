// When you write the compilation code (in crfdb?), make sure to put a reminder here

// Format 
package org.genvisis.park;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;
import org.genvisis.park.tools;

public class ParseRawSNPs {
	public static final String[] MISSING_CALL_CODES = {"Undetermined", "undetermined", "1"};
	public static final String[] MISSING_GENOTYPE_CODES = {"1", "*Undetermined", "U", ""};
	public static final String[] CONTROL_DNA_TYPES = {"Control", "Ceph", "From Core"};
	public static final String[] CONTROL_ID_TYPES = {"H2O", "Blank"};
	public static final String[][] REPLACEMENTS = { {"SNP VIC", "snp_vic"}, {"VIC", "snp_vic"}, {"undetermined", "Undetermined"}, };
	public static final String SNP_POSITIONS = "snp_positions.dat";

	// public static final String DEFAULT_DIR = "C:\\Documents and
	// Settings\\npankrat\\My Documents\\tWork\\SNCA Rep1\\extra SNPs\\";
	// public static final String DEFAULT_DIR = "C:\\Documents and
	// Settings\\npankrat\\My Documents\\_GIGYF2\\new data\\";
//	public static final String DEFAULT_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\mtDNA\\00src\\";
//	public static final String DEFAULT_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\mtDNA\\00src2\\";
//	public static final String DEFAULT_DIR = "C:\\Users\\npankrat\\Documents\\1_CRFdb\\GBA\\1st_RawData\\";
	public static final String DEFAULT_DIR = "C:\\Users\\npankrat\\Documents\\1_CRFdb\\GBA\\2nd_set\\";

	public static final int DEFAULT_VERBOSITY = 8;

	// public static final int DEFAULT_VERBOSITY = 4;

	public static void parse(String dir, int verbosity) {
		BufferedReader reader;
		PrintWriter writer = null;
		String[] line;
		CountHashHash checkCall = new CountHashHash();
		Hashtable<String,String> mergeDNAs = new Hashtable<String,String>();
		Hashtable<String,String> mergeIDs = new Hashtable<String,String>();
		Hashtable<String,String[]> alleleTranslation = new Hashtable<String,String[]>();
		Vector<String[]> v = new Vector<String[]>();
		int count = 0;
		String temp, genotype, snpName;
		int[] indices, counts, keys;
		CheckIDsAgainstDNAs checker = new CheckIDsAgainstDNAs();
		String[] values;
		String allele1, allele2;
		String[] ids;

		File[] files = new File(dir).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".snp");
			}
		});
		System.out.println("Found "+files.length+" files to parse");

		for (int i = 0; i<files.length; i++) {
			snpName = ext.rootOf(files[i].getName());
			try {
				v.clear();
				checkCall.clear();
				mergeDNAs.clear();
				mergeIDs.clear();
				reader = new BufferedReader(new FileReader(files[i]));
				System.out.print("Parsing "+snpName+"...");
				try {
					Thread.sleep(500);
				} catch (InterruptedException ie) {
				}
				indices = Array.intArray(4, -1);
				count = 0;
				if (verbosity>=5) {
					writer = new PrintWriter(new FileWriter(dir+files[i].getName()+"-discarded.out"));
				}
				while (reader.ready()) {
					temp = reader.readLine();
					line = temp.split("\t", -1);
					for (int j = 0; j < line.length; j++) {
						line[j] = ext.replaceAllWith(line[j], REPLACEMENTS);
					}
					count++;
					if (temp.contains("DNA")) {
						indices[0] = ext.indexOfStr("DNA", line, false, false);
						if (temp.contains("Fam")) {
							indices[1] = ext.indexOfStr("Fam", line, false, false);
							if (!line[indices[1]].contains("Ind")&&verbosity>=10) {
								System.err.println("Warning - Make sure column "+ext.getExcelColumn(indices[1])+" is the proper one to find the lab's Fam-Id");
							}
						} else if (verbosity>=10) {
							System.err.println("Error - found 'DNA' on line "+count+" but did not find 'Fam'");
						}
						if (temp.contains("Genotype")) {
							indices[3] = ext.indexOfStr("Genotype", line, false, false);
							if (temp.contains("Call")) {
								indices[2] = ext.indexOfStr("Call", line, false, false);
							} else if (verbosity>=10) {
								System.err.println("Error - found 'DNA' on line "+count+" but did not find 'Call'");
							}
						} else if (temp.contains("Result")) {
							indices[2] = ext.indexOfStr("Result", line, false, false);
							if (temp.contains("Call")) {
								indices[3] = ext.indexOfStr("Call", line, false, false);
							} else if (verbosity>=10) {
								System.err.println("Error - found 'DNA' on line "+count+" but did not find 'Call'");
							}
						} else if (verbosity>=10) {
							System.err.println("Error - found 'DNA' on line "+count+" but did not find 'Genotype' or 'Result'");
						}
						if (verbosity>=5) {
							writer.println(temp);
						}
					} else {
						if (("!"+temp).trim().split("\t").length>=Math.max(Array.max(indices)+1, 2)&&indices[0]!=-1&&!line[indices[1]].trim().equals("")&&(line[indices[0]].trim().equals("")||ext.indexOfStr(line[indices[0]], CONTROL_DNA_TYPES, false, false)==-1)&&ext.indexOfStr(line[indices[1]], CONTROL_ID_TYPES, false, false)==-1) {
							//	// doesn't pick up controls
							if (Array.min(indices)==-1) {
								System.err.println("Error - could not match up all indices for DNA, Fam-Id, Call, and Genotype before trying to add line:");
								System.err.println("    '"+temp+"'");
							}
							if (line[indices[0]].trim().equals("")&&line[indices[1]].startsWith("ND")) {
								line[indices[0]] = line[indices[1]];
							}
							if (line[indices[1]].startsWith("*")) {
								line[indices[1]] = line[indices[1]].substring(1);
							}
							if (line[indices[1]].endsWith("B")) {
								line[indices[1]] = line[indices[1]].substring(0, line[indices[1]].length()-1);
								System.err.println("Error - still found an individual with a 'B' after their name");
							}

							
							if (line[indices[3]].length()==3&&line[indices[3]].charAt(1)=='/') {
								// this messed up 84G
//								line[indices[3]] = line[indices[3]].charAt(0)+""+line[indices[2]];
								line[indices[3]] = line[indices[3]].charAt(0)+""+line[indices[3]].charAt(2);
							}

							checkCall.add(line[indices[2]], line[indices[3]]);
//							System.err.println(line[indices[2]]+"\t"+line[indices[3]]);

							if (ext.indexOfStr(line[indices[2]], MISSING_CALL_CODES)>=0&&ext.indexOfStr(line[indices[3]], MISSING_GENOTYPE_CODES)==-1) {
								System.err.println("Error - Call for "+line[indices[0]]+" is set to missing, but the genotype ("+line[indices[3]]+") is not; add new missing genotype code if necessary");
							}
							if (ext.indexOfStr(line[indices[2]], MISSING_CALL_CODES)==-1&&ext.indexOfStr(line[indices[3]], MISSING_GENOTYPE_CODES)>=0) {
								System.err.println("Error - Genotype for "+line[indices[0]]+" is set to missing, but the call ("+line[indices[2]]+") is not; add new missing call code if necessary");
							}

							if (!mergeDNAs.containsKey(line[indices[0]])) {
								mergeDNAs.put(line[indices[0]], line[indices[3]]);
							} else if (!mergeDNAs.get(line[indices[0]]).equals(line[indices[3]])) {
								if (ext.indexOfStr(mergeDNAs.get(line[indices[0]]), MISSING_GENOTYPE_CODES)>=0) {
									mergeDNAs.put(line[indices[0]], line[indices[3]]);
								} else if (ext.indexOfStr(line[indices[3]], MISSING_GENOTYPE_CODES)==-1) {
									System.err.println("Error - mismatched Genotype for "+line[indices[0]]+" (was "+mergeDNAs.get(line[indices[0]])+", but it is also "+line[indices[3]]+")");
								}
							}
							if (!mergeIDs.containsKey(line[indices[1]])) {
								mergeIDs.put(line[indices[1]], line[indices[3]]);
							} else if (!mergeIDs.get(line[indices[1]]).equals(line[indices[3]])) {
								if (ext.indexOfStr(mergeDNAs.get(line[indices[1]]), MISSING_GENOTYPE_CODES)>=0) {
									mergeDNAs.put(line[indices[1]], line[indices[3]]);
								} else if (ext.indexOfStr(line[indices[3]], MISSING_GENOTYPE_CODES)==-1) {
									System.err.println("Error - mismatched Genotype for "+line[indices[1]]+" (was "+mergeDNAs.get(line[indices[1]])+", but it is also "+line[indices[3]]+")");
								}
							}
							checker.checkPair(line[indices[1]], line[indices[0]], verbosity>=9).startsWith("\t");
							v.add(line);

						} else {
							if (verbosity>=5) {
								writer.println(temp);
							}
						}
					}
				}
				reader.close();
				if (verbosity>=5) {
					writer.close();
				}

				System.out.println("    for a total of "+v.size()+" individuals ("+mergeIDs.size()+" unique)");
				values = new String[v.size()];
				for (int j = 0; j<v.size(); j++) {
					line = v.elementAt(j);
					if (line[indices[0]].startsWith("ND")||line[indices[1]].startsWith("73008")) {
						values[j] = "B_"+line[indices[0]];
					} else if (line[indices[1]].startsWith("70")||line[indices[1]].startsWith("73")) {
						values[j] = "A_"+line[indices[1]];
					} else if (line[indices[1]].startsWith("71")||line[indices[1]].startsWith("95")) {
						values[j] = "C_"+line[indices[1]];
					} else if (line[indices[0]].contains("AD")) {
						values[j] = "D_"+line[indices[0]];
					} else {
						System.err.println("Error - unknown study pattern ("+line[indices[1]]+" / "+line[indices[0]]+")");
					}
				}
				keys = Sort.quicksort(values);

				line = checkCall.getKeys();
				if (line.length>4) {
					System.err.println("Error - There were more than 4 Call types: "+Array.toStr(line, ", "));
				}

				count = 0;
				allele1 = allele2 = "";
				for (int j = 0; j<line.length; j++) {
					values = checkCall.getValues(line[j]);
					counts = checkCall.getCounts(line[j]);
					if (values.length>1) {
						System.err.print("Error - Call '"+line[j]+"' was matched to multiple Genotypes (");
						for (int k = 0; k<values.length; k++) {
							System.err.print((k==0?"":", ")+"'"+values[k]+"' found "+counts[k]+" time"+(counts[k]==1?"":"s"));
						}
						System.err.println(")");
					} else if (ext.indexOfStr(line[j], MISSING_CALL_CODES)==-1) {
						for (int k = 0; k<2; k++) {
							if (allele1.equals("")) {
								allele1 = values[0].substring(k, k+1);
								count += counts[0];
							} else if (values[0].substring(k, k+1).equals(allele1)) {
								count += counts[0];
							} else if (allele2.equals("")) {
								allele2 = values[0].substring(k, k+1);
								count -= counts[0];
							} else if (values[0].substring(k, k+1).equals(allele2)) {
								count -= counts[0];
							} else {
								System.err.println("Error - more than two alleles are present for this SNP ("+allele1+", "+allele2+", and now we see "+values[0]+")");
							}
						}
					}
				}
				if (count<0) {
					temp = allele1;
					allele1 = allele2;
					allele2 = temp;
				}
				alleleTranslation.put(allele1+allele1, new String[] {"1", "1", "2", "0", "1", "0"});
				alleleTranslation.put(allele1+allele2, new String[] {"1", "2", "1", "1", "1", "1"});
				alleleTranslation.put(allele2+allele1, new String[] {"1", "2", "1", "1", "1", "1"});
				alleleTranslation.put(allele2+allele2, new String[] {"2", "2", "0", "2", "0", "1"});
				alleleTranslation.put(".", new String[] {".", ".", ".", ".", ".", "."});

				try {
					writer = new PrintWriter(new FileWriter(dir+snpName+".csv"));
					writer.println("DNA#,FamID,IndID,UniqueID,Genotype,Allele1,Allele2,CountAllele1,CountAllele2,CarrierAllele1,CarrierAllele2");
					for (int j = 0; j<v.size(); j++) {
						line = v.elementAt(keys[j]);
						genotype = mergeIDs.get(line[indices[1]]);
						ids = tools.getFamID(line[indices[1]]);
						if (!genotype.equals("used")) {
							mergeIDs.put(line[indices[1]], "used");
							writer.print(line[indices[0]]+","+ids[0]+","+ids[1]);
							if (line[indices[0]].startsWith("ND")||line[indices[0]].contains("AD")||ids[0].startsWith("71")||ids[0].startsWith("95")) {
								writer.print(","+line[indices[0]]);
							} else {
								writer.print(","+tools.getUniqueID(ids[0], ids[1]));
							}
							writer.print(","+genotype);
							if (ext.indexOfStr(genotype, MISSING_GENOTYPE_CODES)==-1) {
								writer.println(","+Array.toStr(alleleTranslation.get(genotype), ","));
							} else {
								writer.println(","+Array.toStr(alleleTranslation.get("."), ","));
							}
						}
					}
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing to "+dir+snpName+".csv");
				}
			} catch (IOException ioe) {
				System.err.println("Error parsing "+snpName);
				System.exit(3);
			}
			
			try {
				Thread.sleep(500);
			} catch (InterruptedException ie) {
			}
		}

		// try {
		// reader = new BufferedReader(new FileReader(SNP_POSITIONS));
		// while (reader.ready()) {
		// line = reader.readLine().split("[\\s]+");
		// hash.put(line[0], line[1]);
		// }
		// reader.close();
		// } catch (FileNotFoundException fnfe) {
		// System.err.println("Error - could not find "+SNP_POSITIONS+
		// " in current directory");
		// System.exit(2);
		// } catch (IOException ioe) {
		// System.err.println("Error parsing "+SNP_POSITIONS+"");
		// System.exit(3);
		// }
		//
		// if ((count = endpoint-startpoint) % 2 == 1) {
		// System.err.println("Warning - expecting even number of columns
		// between FamID and #FailedSNPs");
		// System.exit(1);
		// }
		// snps = new SNP[count /= 2];
		// snpPositions = new int[snps.length];
		// for (int i = 0; i<count; i++) {
		// snps[i] = new SNP(header[startpoint+i*2]);
		// if (hash.containsKey(snps[i].getName())) {
		// snpPositions[i] =
		// Integer.parseInt(hash.get(snps[i].getName()));
		// } else {
		// System.err.println("Error - marker '"+snps[i].getName()+"' was not
		// listed in "+SNP_POSITIONS);
		// System.exit(1);
		// }
		// }
		// snpOrder = Sort.quicksort(snpPositions);
		//
		// indIDs = new int[v.size()];
		//
		// for (int i = 0; i<v.size(); i++) {
		// line = v.elementAt(i);
		// if (line[2].startsWith("*")) {
		// line[2] = line[2].substring(1);
		// }
		// if (line[1].equals("2002PD0910") && line[2].equals("7051-10")) {
		// line[2] = "70541-10";
		// }
		// if (line[1].equals("2002PD0959") && line[2].equals("7095-18")) {
		// line[2] = "70395-18";
		// }
		//
		// line[0] = line[1];
		// famidPair = tools.getFamID(line[2]);
		// line[1] = famidPair[0];
		// line[2] = famidPair[1];
		// indIDs[i] =
		// Integer.parseInt(famidPair[0])*1000+Integer.parseInt(famidPair[1]);
		// for (int j = 0; j<snps.length; j++) {
		// snps[j].update(line[3+j*2+0], line[3+j*2+1]);
		// }
		// }
		//
		//
		// try {
		// writer = new PrintWriter(new FileWriter(outfile));
		// writer.println("placeholder line");
		// writer.print("DNA\tFamNo\tIndNo");
		// for (int i = 0; i<snps.length; i++) {
		// writer.print("\t"+snps[snpOrder[i]].getName()+"\t");
		// }
		// writer.println();
		//
		// keys = Sort.quicksort(indIDs);
		// for (int i = 0; i<indIDs.length; i++) {
		// line = v.elementAt(keys[i]);
		// writer.print(line[0]+"\t"+line[1]+"\t"+line[2]);
		// for (int j = 0; j<snps.length; j++) {
		// for (int k = 0; k<2; k++) {
		// writer.print("\t"+snps[snpOrder[j]].recode(line[3+snpOrder[j]*2+k]));
		// }
		// }
		// writer.println();
		// }
		// writer.close();
		// } catch (IOException ioe) {
		// System.err.println("Error writing to file");
		// }
		//
		//
		// try {
		// writer = new PrintWriter(new FileWriter(outfile+"-summary.xls"));
		// for (int i = 0; i<snps.length; i++) {
		// writer.println(snps[i]);
		// }
		// writer.close();
		// } catch (IOException ioe) {
		// System.err.println("Error writing to file");
		// }
		//

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String dir = DEFAULT_DIR;
		int verbosity = DEFAULT_VERBOSITY;

		String usage = "\n"+"park.ParseRawSNPs requires 0-1 arguments\n"+"   (1) dir (i.e. dir="+dir+" (default)\n"+"   (2) level of verbosity (i.e. verb="+verbosity+" (default)\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("verb=")) {
				verbosity = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			parse(dir, verbosity);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
