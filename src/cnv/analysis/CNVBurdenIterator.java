package cnv.analysis;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import stats.LeastSquares;
import stats.LogisticRegression;
import stats.RegressionModel;
import cnv.var.CNVariant;
import common.Array;
import common.Files;
import common.HashVec;
import common.ext;

public class CNVBurdenIterator {
	
	static abstract class CNVData {
		
		abstract String getID();
		abstract void setData(int ind, double data);
		abstract double getData(int ind);
		
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
	
	static class MappedData extends CNVData {
		String fidiid;
		HashMap<Integer, Double> data = new HashMap<Integer, Double>();
		
		@Override
		String getID() { return fidiid; }
		
		double getData(int ind) {
			Double d = data.get(ind);
//			if (d == null) {
//				System.out.println(fidiid + " -> " + ind);
//				return Double.NaN;
//			}
			return d.doubleValue(); 
		}
		
		public void setData(int ind, double data) {
			this.data.put(ind, data);
		}
		
		@Override
		public String toString() {
			return "{id: " + fidiid + "; data: " + data.toString() + "}";
		}
		
	}
	
//	static class EduData extends CNVData {
//		String fidiid;
//		boolean college;
//		boolean dropout;
//		double years;
//		double attain;
//	}
//	
//	static class Data extends CNVData {
//		String fidiid;
//		double age;
//		double iq;
//		boolean male;
//		double pc1;
//		double pc2;
//	}
	
	static int[] ID_COLS = {0, 1}; 
	static int[] DATA_COLS = {2, 3, 4, 5};
	static int EXCLUDE_COL = -1;
	
	HashSet<String> usedIDs;
	HashMap<String, MappedData> idData;
	boolean logistic = false;
	boolean collapseCNVTypes = false;
	boolean binaryCNVCounts = false;
	boolean stepSizes = false;
	int genderCol = -1;
	int males = 0, females = 0;
	int depInd;
	int[] indepInds;
	String depHdr;
	String[] indepHdrs;
	int[] CNV_FILTERS;
	double[] CNV_SIZES;
	int[] VAR_SIZES;
	double[][][] resultsB;
	double[][][] resultsP;
	double[][][] resultsR;
	double[][][] resultsSE;
	double[][][] resultsMean;
	
	private void loadIDs(String famFile, String idsFile) {
		HashSet<String> famIDs, tempIDs;
		
		famIDs = new HashSet<String>();
		if (famFile != null && !"".equals(famFile)) {
			famIDs.addAll(HashVec.loadFileToHashString(famFile, ID_COLS, null, false, "\t", false, false, false).keySet());
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
	}
	
	//	
	//	private static void runProgram(String cnvFile, String famFile, String dataFile, String idsFile, boolean collapseCNVTypes, boolean binaryCNVcounts) {
	//		int males = 0, females = 0;
	//		HashSet<String> famIDs, tempIDs, usedIDs;
	//		Hashtable<String, String> dataTable;
	//		HashMap<String, CNVData> idData;
	//		
	//		famIDs = new HashSet<String>();
	//		if (famFile != null && !"".equals(famFile)) {
	//			famIDs.addAll(HashVec.loadFileToHashString(famFile, ID_COLS, null, false, "\t", false, false, false).keySet());
	//		}
	//		
	//		tempIDs = new HashSet<String>();
	//		if (idsFile != null && !"".equals(idsFile)) {
	//			tempIDs = HashVec.loadFileToHashSet(idsFile, false);
	//		}
	//		
	//		usedIDs = new HashSet<String>();
	//		usedIDs.addAll(famIDs);
	//		usedIDs.addAll(tempIDs);
	//		if (famIDs.size() > 0 && tempIDs.size() > 0) {
	//			usedIDs.retainAll(famIDs);
	//			usedIDs.retainAll(tempIDs);
	//		}
	//		
	//		dataTable = HashVec.loadFileToHashString(dataFile, ID_COLS, DATA_COLS, false, "\t", true, false, true);
	//		idData = new HashMap<String, CNVData>();
	//		
	//		for (java.util.Map.Entry<String, String> entry : dataTable.entrySet()) {
	//			String[] dataLine = entry.getValue().split("\t");
	//			
	//			if (shouldSkip(dataLine)) {
	//				continue;
	//			}
	//			
	//			if (!usedIDs.contains(entry.getKey())) {
	//				continue;
	//			}
	//			
	//			// "short_version.xlsx"
	////			Data newData = new Data();
	////			newData.fidiid = entry.getKey();
	////			newData.iq = Double.parseDouble(dataLine[0]);
	////			newData.age = Double.parseDouble(dataLine[1]);
	////			newData.male = dataLine[2].equals("1");
	////			if (newData.male) {
	////				males++;
	////			} else {
	////				females++;
	////			}
	////			newData.pc1 = Double.parseDouble(dataLine[3]);
	////			newData.pc2 = Double.parseDouble(dataLine[4]);
	//			
	//			// "edutainment.txt"
	//			EduData newData = new EduData();
	//			newData.college = "1".equals(dataLine[0]);
	//			newData.dropout = "1".equals(dataLine[1]);
	//			newData.years = Integer.parseInt(dataLine[2]);
	//			newData.attain = Integer.parseInt(dataLine[3]);
	//			
	//			idData.put(newData.fidiid, newData);
	//		}
	//		
	//		Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(cnvFile, null, false);
	//		
	//		for (CNVariant cnv : cnvs) {
	//			CNVData indiv = idData.get(cnv.getFamilyID() + "\t" + cnv.getIndividualID());
	//			if (indiv != null) {
	//				indiv.cnvs.add(cnv);
	//			}
	//		}
	//		
	//		/* 
	//			for each sex : {BOTH, FEMALE, MALE}
	//				filter individuals in data
	//				
	//				for each CN : {0, 1, 3, 4+}
	//					for each length : {5 mb, 2.5 mb, 1 mb, 500 kb, 250 kb, 100kb}
	//						filter cnvs by individuals, CN, and length
	//						set as CNVBurden variable
	//						
	//						run leastsquares
	//		*/
	//		
	//		int[] CNV_FILTERS = collapseCNVTypes ? new int[]{1, 4} : new int[]{0, 1, 3, 4};
	//		double[] CNV_SIZES = new double[]{5000, 2500, 1000, 500, 250, 100};
	//		int[] VAR_SIZES = {idData.size(), males, females};
	//		
	//		double[][][] resultsB = new double[3][CNV_FILTERS.length][CNV_SIZES.length];
	//		double[][][] resultsP = new double[3][CNV_FILTERS.length][CNV_SIZES.length];
	//		double[][][] resultsR = new double[3][CNV_FILTERS.length][CNV_SIZES.length];
	//		double[][][] resultsSE = new double[3][CNV_FILTERS.length][CNV_SIZES.length];
	//
	//		double[] depVars = null;
	//		double[][] indepVars = new double[4][];
	//		String[] indepVarNames = new String[]{"Age", "PC1", "PC2", "CNVBurden"};
	//		
	//		for (int sex = 0; sex < 3; sex++) {
	//			for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
	//				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
	//					depVars = new double[VAR_SIZES[sex]];
	//					indepVars = new double[VAR_SIZES[sex]][4];
	//					int cnt = 0;
	//					ArrayList<Double> popIQs = new ArrayList<Double>();
	//					for (java.util.Map.Entry<String, Data> indiv : idData.entrySet()) {
	//						if (sex > 0) {
	//							if ((indiv.getValue().male && sex == 2) || (!indiv.getValue().male && sex == 1)) {
	//								continue;
	//							}
	//						}
	//						depVars[cnt] = indiv.getValue().iq;
	//						indepVars[cnt][0] = indiv.getValue().age;
	//						indepVars[cnt][1] = indiv.getValue().pc1;
	//						indepVars[cnt][2] = indiv.getValue().pc2;
	//						
	//						indepVars[cnt][3] = indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]) + (collapseCNVTypes ? indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]-1) : 0);
	//						if (binaryCNVcounts) {
	//							indepVars[cnt][3] = Math.min(1, indepVars[cnt][3]);
	//						}
	//						
	//						cnt++;
	//					}
	//					
	//					LeastSquares regression = new LeastSquares(depVars, indepVars, indepVarNames, false, true);
	//					
	//					resultsB[sex][cn][sz] = regression.getBetas().length > 4 ? regression.getBetas()[4] : Double.NaN;
	//					resultsP[sex][cn][sz] = regression.getSigs().length > 4 ? regression.getSigs()[4] : Double.NaN;
	//					resultsR[sex][cn][sz] = regression.getRsquare(); 
	//					resultsSE[sex][cn][sz] = regression.getSEofBs()[4];
	//					
	//					regression.destroy();
	//				}
	//			}
	//		}
	//		
	//		String header = collapseCNVTypes ? "\tDeletions\tDuplications" : "\tCN=0\tCN=1\tCN=3\tCN=4";
	//		StringBuilder header2 = new StringBuilder("Betas\t\t\t\t");
	//		if (!collapseCNVTypes) {
	//			header2 = header2.append("\t\t");
	//		}
	//		header2.append("p-values\t\t\t\t");
	//		if (!collapseCNVTypes) {
	//			header2 = header2.append("\t\t");
	//		}
	//		header2.append("Rsq\t\t\t\t");
	//		if (!collapseCNVTypes) {
	//			header2 = header2.append("\t\t");
	//		}
	//		header2.append("Counts\t\t\t\t");
	//		if (!collapseCNVTypes) {
	//			header2 = header2.append("\t\t");
	//		}
	//		header2.append("StdErr\t\t\t\t");
	//		if (!collapseCNVTypes) {
	//			header2 = header2.append("\t\t");
	//		}
	//		
	//		String outputFile = ext.rootOf(cnvFile, false) + ".burden";
	//		PrintWriter writer;
	//		try {
	//			int[][][] counts = new int[3][CNV_FILTERS.length][CNV_SIZES.length];
	//			for (java.util.Map.Entry<String, Data> indiv : idData.entrySet()) {
	////				String id = indiv.getKey();
	//				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
	//					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
	//						int val = indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]) + (collapseCNVTypes ? indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn] - 1) : 0);
	//						counts[0][cn][sz] += val;
	//						if (indiv.getValue().male) {
	//							counts[1][cn][sz] += val;
	//						} else {
	//							counts[2][cn][sz] += val;
	//						}
	//					}
	//				}
	//			}
	//			
	//			writer = new PrintWriter(new FileWriter(outputFile));
	//			
	//			String[] labels = {"Both", "Males", "Females"};
	//			
	//			for (int sex = 0; sex < 3; sex++) {
	//				writer.println(header2.toString());
	//				writer.println(labels[sex] + header + "\t\t" + labels[sex] + header + "\t\t" + labels[sex] + header + "\t\t" + labels[sex] + header + "\t\t" + labels[sex] + header + "\tn=" + VAR_SIZES[sex]);
	//				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
	//					writer.print(CNV_SIZES[sz]);
	//					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
	//						writer.print("\t" + ext.formDeci(resultsB[sex][cn][sz], ext.getNumSigFig(resultsB[sex][cn][sz])));
	//					}
	//					writer.print("\t\t");
	//					writer.print(CNV_SIZES[sz]);
	//					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
	//						writer.print("\t" + ext.formDeci(resultsP[sex][cn][sz], ext.getNumSigFig(resultsP[sex][cn][sz])));
	//					}
	//					writer.print("\t\t");
	//					writer.print(CNV_SIZES[sz]);
	//					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
	//						writer.print("\t" + ext.formDeci(resultsR[sex][cn][sz], ext.getNumSigFig(resultsR[sex][cn][sz])));
	//					}
	//					
	//					writer.print("\t\t");
	//					writer.print(CNV_SIZES[sz]);
	//					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
	//						writer.print("\t" + counts[sex][cn][sz]);
	//					}
	//					
	//					writer.print("\t\t");
	//					writer.print(CNV_SIZES[sz]);
	//					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
	//						writer.print("\t" + resultsSE[sex][cn][sz]);
	//					}
	//					
	//					writer.println();
	//					
	//				}
	//				writer.println();
	//			}
	//			
	//			writer.close();
	//
	//		} catch (IOException e) {
	//			e.printStackTrace();
	//		}
	//	}
		
		
	private static boolean shouldSkip(String[] dataLine) {
		for (String str : dataLine) {
			if (!ext.isValidDouble(str)) return true;
		}
		return false;
	}

	private void loadIDData(String dataFile, int[] idCols, int[] dataCols) {
		String[] hdr = Files.getHeaderOfFile(dataFile, null);
		depHdr = hdr[depInd]; 
		indepHdrs = new String[indepInds.length];
		for (int i = 0; i < indepInds.length; i++) {
			indepHdrs[i] = hdr[indepInds[i]];
		}
		int[] inclGenderDataCols;
		if (genderCol >= 0) {
			boolean found = false;
			for (int i : dataCols) {
				if (i == genderCol) {
					found = true;
					break;
				}
			}
			if (!found) {
				inclGenderDataCols = new int[dataCols.length + 1];
				for (int i = 0; i < dataCols.length; i++) {
					inclGenderDataCols[i] = dataCols[i];
				}
				inclGenderDataCols[dataCols.length] = genderCol;
			} else {
				inclGenderDataCols = dataCols;
			}
		} else {
			inclGenderDataCols = dataCols;
		}
		Hashtable<String, String> dataTable = HashVec.loadFileToHashString(dataFile, idCols, inclGenderDataCols, false, "\t", true, false, true);
		idData = new HashMap<String, MappedData>();
		for (java.util.Map.Entry<String, String> entry : dataTable.entrySet()) {
			String[] dataLine = entry.getValue().split("\t");
			
			if (shouldSkip(dataLine)) {
				continue;
			}
			
			if (!usedIDs.contains(entry.getKey())) {
				continue;
			}
			
			MappedData newData = new MappedData();
			newData.fidiid = entry.getKey();
			boolean setGender = false;
			for (int i = 0; i < dataCols.length; i++) {
				Double value = Double.parseDouble(dataLine[i]);
				newData.setData(dataCols[i], value);
				if (dataCols[i] == genderCol) {
					if (value.doubleValue() > 0.5) {
						males++;
					} else {
						females++;
					}
					setGender = true;
				}
			}
			if (!setGender && inclGenderDataCols.length - dataCols.length == 1) {
				Double gender = Double.parseDouble(dataLine[dataCols.length]);
				newData.setData(genderCol, gender);
				if (gender.doubleValue() > 0.5) {
					males++;
				} else {
					females++;
				}
			}

			idData.put(newData.fidiid, newData);
		}
		
	}
	
	private void loadCNVs(String cnvFile) {
		Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(cnvFile, null, false);
		
		for (CNVariant cnv : cnvs) {
			CNVData indiv = idData.get(cnv.getFamilyID() + "\t" + cnv.getIndividualID());
			if (indiv != null) {
				indiv.cnvs.add(cnv);
			}
		}
		
	}
	
	private void runCalc() {
		CNV_FILTERS = collapseCNVTypes ? new int[]{1, 4} : new int[]{0, 1, 3, 4};
		CNV_SIZES = new double[]{5000, 2500, 1000, 500, 250, 100, 0.001 };
		if (genderCol >= 0) {
			VAR_SIZES = new int[]{idData.size(), males, females}; 
		} else {
			VAR_SIZES = new int[]{idData.size()};
		}
		
		resultsB = new double[VAR_SIZES.length][CNV_FILTERS.length][CNV_SIZES.length];
		resultsP = new double[VAR_SIZES.length][CNV_FILTERS.length][CNV_SIZES.length];
		resultsR = new double[VAR_SIZES.length][CNV_FILTERS.length][CNV_SIZES.length];
		resultsSE = new double[VAR_SIZES.length][CNV_FILTERS.length][CNV_SIZES.length];
		resultsMean = new double[VAR_SIZES.length][CNV_FILTERS.length][CNV_SIZES.length];

		double[] depVars = null;
		double[][] indepVars = new double[indepInds.length][];
		String[] indepVarNames = new String[indepHdrs.length + 1];
		for (int i = 0; i < indepHdrs.length; i++) {
			indepVarNames[i] = indepHdrs[i];
		}
		indepVarNames[indepHdrs.length] = "CNVBurden";
		
		for (int sex = 0; sex < VAR_SIZES.length; sex++) {
			for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
					depVars = new double[VAR_SIZES[sex]];
					indepVars = new double[VAR_SIZES[sex]][indepInds.length + 1];
					int cnt = 0;
					for (java.util.Map.Entry<String, MappedData> indiv : idData.entrySet()) {
						if (sex > 0) {
							double doub = indiv.getValue().getData(genderCol); 
							if ((doub > 0.5 && sex == 2) || (doub < 0.5 && sex == 1)) {
								continue;
							}
						}
						
						depVars[cnt] = indiv.getValue().getData(depInd);
						
						if (stepSizes) {
							float larger = indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]) + (collapseCNVTypes ? indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]-1) : 0);
							float smaller = sz == 0 ? 0 : indiv.getValue().getCNVs(CNV_SIZES[sz - 1], CNV_FILTERS[cn]) + (collapseCNVTypes ? indiv.getValue().getCNVs(CNV_SIZES[sz - 1], CNV_FILTERS[cn]-1) : 0);
							indepVars[cnt][0] = larger - smaller;
						} else {
							indepVars[cnt][0] = indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]) + (collapseCNVTypes ? indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]-1) : 0);
						}
						if (binaryCNVCounts) {
							indepVars[cnt][0] = Math.min(1, indepVars[cnt][0]);
						}
						
						for (int i = 1; i < indepInds.length + 1; i++) {
							indepVars[cnt][i] = indiv.getValue().getData(indepInds[i - 1]);
						}
						
						cnt++;
					}
					
					RegressionModel model;
					if (logistic) {
						model = new LogisticRegression(depVars, indepVars, indepVarNames, false, true);
					} else {
						model = new LeastSquares(depVars, indepVars, indepVarNames, false, true);
					}

					int ind = -1;
					for (int i = 0; i < model.getVarNames().length; i++) {
						if ("Indep 1".equals(model.getVarNames()[i])) {
							ind = i;
							break;
						}
					}
					if (ind == -1) {
						resultsB[sex][cn][sz] = Double.NaN;
						resultsP[sex][cn][sz] = Double.NaN;
						resultsR[sex][cn][sz] = model.getRsquare(); 
						resultsSE[sex][cn][sz] = Double.NaN;
					} else {
						resultsB[sex][cn][sz] = model.getBetas()[ind];
						resultsP[sex][cn][sz] = model.getSigs()[ind];
						resultsR[sex][cn][sz] = model.getRsquare(); 
						resultsSE[sex][cn][sz] = model.getSEofBs()[ind];
					}
					
					model = null;
				}
			}
		}
	}
	
	private void writeOutput(String outFile) {
		String header = collapseCNVTypes ? "\tDeletions\tDuplications" : "\tCN=0\tCN=1\tCN=3\tCN=4";
		StringBuilder header2 = new StringBuilder("Betas\t\t\t\t");
		if (!collapseCNVTypes) {
			header2 = header2.append("\t\t");
		}
		header2.append("p-values\t\t\t\t");
		if (!collapseCNVTypes) {
			header2 = header2.append("\t\t");
		}
		header2.append("Rsq\t\t\t\t");
		if (!collapseCNVTypes) {
			header2 = header2.append("\t\t");
		}
		header2.append("Counts\t\t\t\t");
		if (!collapseCNVTypes) {
			header2 = header2.append("\t\t");
		}
		header2.append("StdErr\t\t\t\t");
		if (!collapseCNVTypes) {
			header2 = header2.append("\t\t");
		}
		header2.append("MeanIQ\t\t\t\t");
		if (!collapseCNVTypes) {
			header2 = header2.append("\t\t");
		}
		
		PrintWriter writer;
		try {
			int[][][] counts = new int[genderCol >= 0 ? 3 : 1][CNV_FILTERS.length][CNV_SIZES.length];
			for (java.util.Map.Entry<String, MappedData> indiv : idData.entrySet()) {
				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						int val = indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn]) + (collapseCNVTypes ? indiv.getValue().getCNVs(CNV_SIZES[sz], CNV_FILTERS[cn] - 1) : 0);
						if (stepSizes) {
							int smaller = sz == 0 ? 0 : indiv.getValue().getCNVs(CNV_SIZES[sz - 1], CNV_FILTERS[cn]) + (collapseCNVTypes ? indiv.getValue().getCNVs(CNV_SIZES[sz - 1], CNV_FILTERS[cn]-1) : 0);
							val -= smaller;
						}
						counts[0][cn][sz] += val > 0 ? 1 : 0;
						if (val > 0) {
							resultsMean[0][cn][sz] += indiv.getValue().getData(depInd);
						}
						if (genderCol >= 0) {
							if (indiv.getValue().getData(genderCol) > 0.5) {
								counts[1][cn][sz] += val > 0 ? 1 : 0;
								resultsMean[1][cn][sz] += val > 0 ? indiv.getValue().getData(depInd) : 0;
							} else {
								counts[2][cn][sz] += val > 0 ? 1 : 0;
								resultsMean[2][cn][sz] += val > 0 ? indiv.getValue().getData(depInd) : 0;
							}
						}
					}
				}
			}
			
			writer = new PrintWriter(new FileWriter(outFile));
			
			String[] labels = genderCol >= 0 ? new String[]{"BothSexes", "Males", "Females"} : new String[]{"BothSexes"};
			
			for (int sex = 0; sex < labels.length; sex++) {
				writer.println(header2.toString());
				writer.println(labels[sex] + header + "\t\t" + labels[sex] + header + "\t\t" + labels[sex] + header + "\t\t" + labels[sex] + header + "\t\t" + labels[sex] + header + "\t\t" + labels[sex] + header + "\tn=" + VAR_SIZES[sex]);
				for (int sz = 0; sz < CNV_SIZES.length; sz++) {
					writer.print(CNV_SIZES[sz]);
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						writer.print("\t" + ext.formDeci(resultsB[sex][cn][sz], ext.getNumSigFig(resultsB[sex][cn][sz])));
					}
					writer.print("\t\t");
					writer.print(CNV_SIZES[sz]);
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						writer.print("\t" + ext.formDeci(resultsP[sex][cn][sz], ext.getNumSigFig(resultsP[sex][cn][sz])));
					}
					writer.print("\t\t");
					writer.print(CNV_SIZES[sz]);
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						writer.print("\t" + ext.formDeci(resultsR[sex][cn][sz], ext.getNumSigFig(resultsR[sex][cn][sz])));
					}
					
					writer.print("\t\t");
					writer.print(CNV_SIZES[sz]);
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						writer.print("\t" + counts[sex][cn][sz]);
					}
					
					writer.print("\t\t");
					writer.print(CNV_SIZES[sz]);
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						writer.print("\t" + resultsSE[sex][cn][sz]);
					}
					
					writer.print("\t\t");
					writer.print(CNV_SIZES[sz]);
					for (int cn = 0; cn < CNV_FILTERS.length; cn++) {
						writer.print("\t" + (resultsMean[sex][cn][sz] / counts[sex][cn][sz]));
					}
					
					writer.println();
				}
				writer.println();
			}
			
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void runProgram(String cnvFile, String dataFile, String famFile, String idsFile, int[] idCols, int[] dataCols, int depCol, int genderCol, boolean logistic, boolean collapseCNVTypes, boolean binaryCNVcounts, boolean stepSizes) {
		this.logistic = logistic;
		this.genderCol = genderCol;
		this.collapseCNVTypes = collapseCNVTypes;
		this.binaryCNVCounts = binaryCNVcounts;
		this.stepSizes = stepSizes;
		depInd = depCol;
		indepInds = new int[dataCols.length-1];
		int ind = 0;
		for (int i : dataCols) {
			if (i == depCol) continue;
			indepInds[ind] = i;
			ind++;
		}
		loadIDs(famFile, idsFile);
		loadIDData(dataFile, idCols, dataCols);
		loadCNVs(cnvFile);
		runCalc();
		writeOutput(ext.rootOf(cnvFile, false) + ".burden");
	}
	
	public static void main(String[] args) {
		String cnvFile = "D:/SIDS and IQ/penncnv.cnv";
		String famFile = "D:/SIDS and IQ/penncnv.fam";
		String dataFile = "D:/SIDS and IQ/short_version.txt";
		String idsFile = "";//"D:/SIDS and IQ/ids_with_high_quality_CNV_data_from_blood_cluster.unrelateds_keep.dat";
		boolean collapse = false;
		boolean binaryCnts = false;
		boolean logisticRun = false;
		boolean stepwise = false;
		int[] dataCols = new int[]{2, 3, 5, 6}; 
		int[] idsCols = new int[]{0, 1};
		int depCol = 2;
		int genCol = 4;
		
		String usage = "\n"+
				"one.CNVBurdenIterator requires 2-4+ arguments:\n" +
				"   (0) cnvFile \n" + 
				"   (1) dataFile \n" +
				"   (2) famFile \n" +
				"   (4) idsFile \n" +
				"   (5) genderColumn \n" +
				"   (6) idColumns \n" + 
				"   (7) dataColumns \n" + 
				"   (8) dependColumn \n" + 
				"   (9) -collapse \n" + 
				"  (10) -binary \n" + 
				"  (11) -logistic \n" + 
				"  (12) -stepwise \n";
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
			} else if (args[i].startsWith("genderColumn=")) {
				genCol = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("idColumns=")) {
				String[] cols = args[i].split("=")[1].split(",");
				idsCols = new int[cols.length];
				for (int j = 0; j < cols.length; j++) {
					idsCols[j] = Integer.parseInt(cols[j]);
				}
				numArgs--;
			} else if (args[i].startsWith("dataColumns=")) {
				String[] cols = args[i].split("=")[1].split(",");
				dataCols = new int[cols.length];
				for (int j = 0; j < cols.length; j++) {
					dataCols[j] = Integer.parseInt(cols[j]);
				}
				numArgs--;
			} else if (args[i].startsWith("dependColumn=")) {
				depCol = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-collapse")) {
				collapse = true;
				numArgs--;
			} else if (args[i].startsWith("-binary")) {
				binaryCnts = true;
				numArgs--;
			} else if (args[i].startsWith("-logistic")) {
				logisticRun = true;
				numArgs--;
			} else if (args[i].startsWith("-stepwise")) {
				stepwise = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			return;
		}
		
		try {
			CNVBurdenIterator burden = new CNVBurdenIterator();
			burden.runProgram(cnvFile, dataFile, famFile, idsFile, idsCols, dataCols, depCol, genCol, logisticRun, collapse, binaryCnts, stepwise);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
 