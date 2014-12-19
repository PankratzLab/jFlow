package cnv.analysis;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import common.Array;
import common.HashVec;
import common.ext;
import cnv.filesys.Project;
import cnv.var.CNVariant;
import stats.LeastSquares;

public class CNVBurdenIterator {
	
	static class Data {
		String fidiid;
		double age;
		double iq;
		boolean male;
		double pc1;
		double pc2;
		final ArrayList<CNVariant> cnvs = new ArrayList<CNVariant>();
		
		public int getCNVs(double minSizeKB) {
			int cnt = 0;
			for (CNVariant cnv : cnvs) {
				if (cnv.getSize() >= minSizeKB * 1000) {
					cnt++;
				}
			}
			return cnt;
		}
		
		public int getCNVs(int cn) {
			int cnt = 0;
			for (CNVariant cnv : cnvs) {
				if (cnv.getCN() == cn) {
					cnt++;
				}
			}
			return cnt;
		}
		
		public int getCNVs(double minSizeKB, int cn) {
			int cnt = 0;
			for (CNVariant cnv : cnvs) {
				if (cnv.getCN() == cn && cnv.getSize() >= minSizeKB * 1000) {
					cnt++;
				}
			}
			return cnt;
		}
		
	}

	private static void runProgram(String cnvFile, String famFile, String dataFile, String idsFile, boolean collapseCNVTypes, boolean binaryCNVcounts) {
		int males = 0, females = 0;
		HashSet<String> famIDs, tempIDs, usedIDs;
		Hashtable<String, String> dataTable;
		HashMap<String, Data> idData;
		
		famIDs = new HashSet<String>();
		if (famFile != null && !"".equals(famFile)) {
			famIDs.addAll(HashVec.loadFileToHashString(famFile, new int[]{0, 1}, null, false, "\t", false, false, false).keySet());
		}
		
		tempIDs = new HashSet<String>();
		if (idsFile != null && !"".equals(idsFile)) {
			tempIDs = HashVec.loadFileToHashSet(idsFile, false);
		}
		
		usedIDs = new HashSet<String>();
		usedIDs.addAll(famIDs);
		usedIDs.addAll(tempIDs);
		if (famIDs.size() > 0 && tempIDs.size() > 0) {
			usedIDs.retainAll(famIDs);
			usedIDs.retainAll(tempIDs);
		}
		
		dataTable = HashVec.loadFileToHashString(dataFile, new int[]{0, 1}, new int[]{2, 3, 4, 5, 6, 9}, false, "\t", true, false, true);
		
		idData = new HashMap<String, Data>();
		
		for (java.util.Map.Entry<String, String> entry : dataTable.entrySet()) {
			String[] dataLine = entry.getValue().split("\t");
			if (dataLine[5].equals("0") || 
					!ext.isValidDouble(dataLine[0]) || 
					!ext.isValidDouble(dataLine[1]) ||
					!ext.isValidDouble(dataLine[3]) || 
					!ext.isValidDouble(dataLine[4])) {
				continue;
			}
			if (!usedIDs.contains(entry.getKey())) {
				continue;
			}
			Data newData = new Data();
			newData.fidiid = entry.getKey();
			newData.iq = Double.parseDouble(dataLine[0]);
			newData.age = Double.parseDouble(dataLine[1]);
			newData.male = dataLine[2].equals("1");
			if (newData.male) {
				males++;
			} else {
				females++;
			}
			newData.pc1 = Double.parseDouble(dataLine[3]);
			newData.pc2 = Double.parseDouble(dataLine[4]);
			idData.put(newData.fidiid, newData);
		}
		
		Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(cnvFile, null, false);
		
		for (CNVariant cnv : cnvs) {
			Data indiv = idData.get(cnv.getFamilyID() + "\t" + cnv.getIndividualID());
			if (indiv != null) {
				indiv.cnvs.add(cnv);
			}
		}
		
		/* 
			for each sex : {BOTH, FEMALE, MALE}
				filter individuals in data
				
				for each CN : {0, 1, 3, 4+}
					for each length : {5 mb, 2.5 mb, 1 mb, 500 kb, 250 kb, 100kb}
						filter cnvs by individuals, CN, and length
						set as CNVBurden variable
						
						run leastsquares
		*/
		
		int[] CNV_FILTERS = collapseCNVTypes ? new int[]{1, 4} : new int[]{0, 1, 3, 4};
		double[] CNV_SIZES = new double[]{5000, 2500, 1000, 500, 250, 100};
		int[] VAR_SIZES = {idData.size(), males, females};
		
		double[][][] resultsB = new double[3][CNV_FILTERS.length][CNV_SIZES.length];
		double[][][] resultsP = new double[3][CNV_FILTERS.length][CNV_SIZES.length];
		double[][][] resultsR = new double[3][CNV_FILTERS.length][CNV_SIZES.length];
		
		
		double[] depVars = null;
		double[][] indepVars = new double[4][];
		String[] indepVarNames = new String[]{"Age", "PC1", "PC2", "CNVBurden"};
		
		
		for (int sex = 0; sex < 3; sex++) {
			for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
					depVars = new double[VAR_SIZES[sex]];
					indepVars = new double[VAR_SIZES[sex]][4];
					int cnt = 0;
					for (java.util.Map.Entry<String, Data> indiv : idData.entrySet()) {
						if (sex > 0) {
							if ((indiv.getValue().male && sex == 2) || (!indiv.getValue().male && sex == 1)) {
								continue;
							}
						}
						depVars[cnt] = indiv.getValue().iq;
						indepVars[cnt][0] = indiv.getValue().age;
						indepVars[cnt][1] = indiv.getValue().pc1;
						indepVars[cnt][2] = indiv.getValue().pc2;
						
						indepVars[cnt][3] = indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]) + (collapseCNVTypes ? indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]-1) : 0);
						if (binaryCNVcounts) {
							indepVars[cnt][3] = Math.min(1, indepVars[cnt][3]);
						}
						
						cnt++;
					}
					LeastSquares regression = new LeastSquares(depVars, indepVars, indepVarNames, false, true);
					
					resultsB[sex][cn][sz] = regression.getBetas().length > 4 ? regression.getBetas()[4] : Double.NaN;
					resultsP[sex][cn][sz] = regression.getSigs().length > 4 ? regression.getSigs()[4] : Double.NaN;
					resultsR[sex][cn][sz] = regression.getRsquare(); 
					
					regression.destroy();
				}
			}
		}
		
		String header = collapseCNVTypes ? "\tDeletions\tDuplications" : "\tCN=0\tCN=1\tCN=3\tCN=4";
		
		String outputFile = ext.rootOf(cnvFile, false) + ".burden";
		PrintWriter writer;
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(cnvFile, false) + ".counts.xln"));

			writer.print("FID\tIID");
			for (int sz = 0; sz < CNV_SIZES.length; sz++) {
				for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
					writer.print("\t" + CNV_SIZES[sz] + "_" + CNV_FILTERS[cn] + (collapseCNVTypes ? "/" + (CNV_FILTERS[cn] - 1) : ""));
				}
			}
			writer.println();
			
			for (java.util.Map.Entry<String, Data> indiv : idData.entrySet()) {
				String id = indiv.getKey();
				writer.print(id);
				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						String value = ""+(indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]) + (collapseCNVTypes ? indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn] - 1) : 0));
						writer.print("\t" + value);
					}
				}
				writer.println();
			}
			writer.close();
			
			writer = new PrintWriter(new FileWriter(outputFile));
			
			String[] labels = {"Both", "Males", "Females"};
			for (int sex = 0; sex < 3; sex++) {
				writer.println("Betas");
				writer.println(labels[sex] + header);
				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
					writer.print(CNV_SIZES[sz]);
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						writer.print("\t" + ext.formDeci(resultsB[sex][cn][sz], ext.getNumSigFig(resultsB[sex][cn][sz])));
					}
					writer.println();
				}
				writer.println("p-values");
				writer.println(labels[sex] + header);
				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
					writer.print(CNV_SIZES[sz]);
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						writer.print("\t" + ext.formDeci(resultsP[sex][cn][sz], ext.getNumSigFig(resultsP[sex][cn][sz])));
					}
					writer.println();
				}
				writer.println("Rsq");
				writer.println(labels[sex] + header);
				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
					writer.print(CNV_SIZES[sz]);
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						writer.print("\t" + ext.formDeci(resultsR[sex][cn][sz], ext.getNumSigFig(resultsR[sex][cn][sz])));
					}
					writer.println();
				}
			}
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public static void main(String[] args) {
		String cnvFile = "D:/SIDS and IQ/penncnv.cnv";
		String famFile = "D:/SIDS and IQ/penncnv.fam";
		String dataFile = "D:/SIDS and IQ/short_version.txt";
		String idsFile = "";//"D:/SIDS and IQ/ids_with_high_quality_CNV_data_from_blood_cluster.unrelateds_keep.dat";
		boolean collapse = false;
		boolean binaryCnts = false;
		
		
		String usage = "\n"+
				"one.CNVBurdenIterator requires 2-4+ arguments:\n" +
				"   (0) cnvFile " + 
				"   (1) dataFile " +
				"   (2) famFile " +
				"   (4) idsFile ";
		int numArgs = args.length;
		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				return;
			} else if (args[i].startsWith("cnvFile=")) {
				cnvFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("famFile=")) {
				famFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dataFile=")) {
				dataFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("idsFile=")) {
				idsFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-collapse")) {
				collapse = true;
				numArgs--;
			} else if (args[i].startsWith("-binary")) {
				binaryCnts = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			return;
		}
		
		try {
			runProgram(cnvFile, famFile, dataFile, idsFile, collapse, binaryCnts);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}