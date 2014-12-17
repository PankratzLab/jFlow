package one;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import common.Array;
import common.HashVec;
import common.ext;
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
	
	 static class Results {
		double[] betas;
		double rsquared;
		double[] sigs;
		
		public Results(double[] b, double r, double[] sigs) {
			this.betas = b;
			this.rsquared = r;
			this.sigs = sigs;
		}
		
		@Override
		public String toString() {
			return rsquared + "\t" + Array.toStr(betas) + "\t" + Array.toStr(sigs);
		}
		
	}
	
	
	private static void runProgram() {
		String cnvFile = "D:/SIDS and IQ/penncnv.cnv";
		String famFile = "D:/SIDS and IQ/penncnv.fam";
		String dataFile = "D:/SIDS and IQ/short_version.txt";
		String idsFile = "";//"D:/SIDS and IQ/ids_with_high_quality_CNV_data_from_blood_cluster.unrelateds_keep.dat";
		int males = 0, females = 0;
		HashSet<String> famIDs, tempIDs, usedIDs;
		Hashtable<String, String> dataTable;
		HashMap<String, Data> idData;
		HashMap<String, Results> resultsMap;
		
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
			for each sex : {FEMALE, MALE, BOTH}
				filter individuals in data
			
				for each CN : {0, 1, 3, 4+}
					for each length : {5 mb, 2.5 mb, 1 mb, 500 kb, 250 kb, 100kb}
						filter cnvs by individuals, CN, and length
						set as CNVBurden variable
						
						run leastsquares
		*/
		
		resultsMap = new HashMap<String, CNVBurdenIterator.Results>();
		
		int[] CNV_FILTERS = new int[]{0, 1, 3, 4};
		double[] CNV_SIZES = new double[]{5000, 2500, 1000, 500, 250, 100};
		
		double[] depVars = null;
		double[][] indepVars = new double[4][];
		String[] indepVarNames = new String[]{"Age", "PC1", "PC2", "CNVBurden"};
		
		// both sexes
		for (int cn : CNV_FILTERS) {
			for (double sz : CNV_SIZES) {
				depVars = new double[idData.size()];
				indepVars = new double[idData.size()][4];
				int cnt = 0;
				for (java.util.Map.Entry<String, Data> indiv : idData.entrySet()) {
					depVars[cnt] = indiv.getValue().iq;
					indepVars[cnt][0] = indiv.getValue().age;
					indepVars[cnt][1] = indiv.getValue().pc1;
					indepVars[cnt][2] = indiv.getValue().pc2;
					indepVars[cnt][3] = indiv.getValue().getCNVs(sz, cn);
					cnt++;
				}
				LeastSquares lsBoth = new LeastSquares(depVars, indepVars, indepVarNames, false, true);
				Results res = new Results(lsBoth.getBetas(), lsBoth.getRsquare(), lsBoth.getSigs());
				resultsMap.put("both_cn" + cn + "_p" + sz, res);
				lsBoth.destroy();
			}
		}
		
		// males only
		for (int cn : CNV_FILTERS) {
			for (double sz : CNV_SIZES) {
				depVars = new double[males];
				indepVars = new double[males][4];
				int cnt = 0;
				for (java.util.Map.Entry<String, Data> indiv : idData.entrySet()) {
					if (!indiv.getValue().male) continue;
					depVars[cnt] = indiv.getValue().iq;
					indepVars[cnt][0] = indiv.getValue().age;
					indepVars[cnt][1] = indiv.getValue().pc1;
					indepVars[cnt][2] = indiv.getValue().pc2;
					indepVars[cnt][3] = indiv.getValue().getCNVs(sz, cn);
					cnt++;
				}
				LeastSquares lsMales = new LeastSquares(depVars, indepVars, indepVarNames, false, true);
				Results res = new Results(lsMales.getBetas(), lsMales.getRsquare(), lsMales.getSigs());
				resultsMap.put("male_cn" + cn + "_p" + sz, res);
				lsMales.destroy();
			}
		}
		
		// females only
		for (int cn : CNV_FILTERS) {
			for (double sz : CNV_SIZES) {
				depVars = new double[females];
				indepVars = new double[females][4];
				int cnt = 0;
				for (java.util.Map.Entry<String, Data> indiv : idData.entrySet()) {
					if (indiv.getValue().male) continue;
					depVars[cnt] = indiv.getValue().iq;
					indepVars[cnt][0] = indiv.getValue().age;
					indepVars[cnt][1] = indiv.getValue().pc1;
					indepVars[cnt][2] = indiv.getValue().pc2;
					indepVars[cnt][3] = indiv.getValue().getCNVs(sz, cn);
					cnt++;
				}
				LeastSquares lsFemales = new LeastSquares(depVars, indepVars, indepVarNames, false, true);
				Results res = new Results(lsFemales.getBetas(), lsFemales.getRsquare(), lsFemales.getSigs());
				resultsMap.put("female_cn" + cn + "_p" + sz, res);
				lsFemales.destroy();
			}
		}
		
		for (java.util.Map.Entry<String, Results> entry : resultsMap.entrySet()) {
			System.out.println(entry.getKey() + "\t" + entry.getValue().toString());
		}
		System.out.println();
	}
	
	public static void main(String[] args) {
		runProgram();
	}
	
}

