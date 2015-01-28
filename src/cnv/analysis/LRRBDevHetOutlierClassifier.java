package cnv.analysis;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import cnv.filesys.Project;
import cnv.var.SampleData;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

public class LRRBDevHetOutlierClassifier {
	
	static final class IndividualDatum {
		final String id;
		final double dataLRR;
		final double dataBDev;
		
		public IndividualDatum(String id, double lrr, double bdev) {
			this.id = id;
			this.dataLRR = lrr;
			this.dataBDev = bdev;
		}
	}
	
	static final class PopulationData {
		final double stdDevBDev;
		final double stdDevLRR;
		final double meanBDev;
		final double meanLRR;
		final double medianBDev;
		final double medianLRR;
		
		public PopulationData(double bdevMean, double lrrMean, double bdevSD, double lrrSD, double bdevMed, double lrrMed) {
			this.stdDevBDev = bdevSD;
			this.stdDevLRR = lrrSD;
			this.meanBDev = bdevMean;
			this.meanLRR = lrrMean;
			this.medianBDev = bdevMed;
			this.medianLRR = lrrMed;
		}
		
	}
	
	enum OutlierTest {
		
		SD_SCORE() {
			public boolean testIndiv(PopulationData popData, IndividualDatum indivDatum) {
				double sdCntBDev = (indivDatum.dataBDev - popData.meanBDev) / popData.stdDevBDev;
				double sdCntLRR = (indivDatum.dataLRR - popData.meanLRR) / popData.stdDevLRR;
				// 4.0 derived through examination of results from other measures - outliers classified by other measures all exceed 4.0 on this test. 
				return sdCntBDev > 4.0 || sdCntLRR > 4.0 || sdCntBDev < -4.0 || sdCntLRR < -4.0;
			}
			public int getMethodID() {
				return 1;
			}
		},
		/**
		 * Test if BDev exceeds 0.1 and LRR exceeds 0.1 or -0.1
		 */
		PT_1_TEST() {
			public boolean testIndiv(PopulationData popData, IndividualDatum indivDatum) {
//				System.out.print("\t" + indivDatum.id + "\t");
				if (indivDatum.dataBDev >= 0.1
						&& (indivDatum.dataLRR >= 0.1 || indivDatum.dataLRR <= -0.1)) {
					return true;
				}
				return false;
			}
			public int getMethodID() {
				return 2;
			}
		},
		/**
		 * Test if BDev or LRR exceed more than a certain number of SD's (default 4) from the mean
		 */
		EXTREME_SD() {
			private static final int EXTREME_SD_LIMIT = 4;
			public boolean testIndiv(PopulationData popData, IndividualDatum indivDatum) {
				if ((indivDatum.dataBDev > (popData.meanBDev + (popData.stdDevBDev * EXTREME_SD_LIMIT)))
						|| (indivDatum.dataBDev < (popData.meanBDev - (popData.stdDevBDev * EXTREME_SD_LIMIT)))
						|| (indivDatum.dataLRR > (popData.meanLRR + (popData.stdDevLRR * EXTREME_SD_LIMIT)))
						|| (indivDatum.dataLRR < (popData.meanLRR - (popData.stdDevLRR * EXTREME_SD_LIMIT)))) {
					return true;
				}
				return false;
			}
			public int getMethodID() {
				return 3;
			}
		},
		/**
		 * Test if BDev exceeds 0.1 and LRR is exceeds more than 2 SD's from the mean
		 */
		BD_PT1_LRRSD() {
			public boolean testIndiv(PopulationData popData, IndividualDatum indivDatum) {
				if (indivDatum.dataBDev >= 0.1
						&& (indivDatum.dataLRR >= (popData.meanLRR + (popData.stdDevLRR * 2)) || indivDatum.dataLRR <= (popData.meanLRR - (popData.stdDevLRR * 2)))) {
					return true;
				}
				return false;
			}
			public int getMethodID() {
				return 4;
			}
		};
		
		abstract boolean testIndiv(PopulationData popData, IndividualDatum indivDatum);
		abstract int getMethodID();
	}
	
	public LRRBDevHetOutlierClassifier(String file) {
		this.filename = file;
	}
	
	String filename;
	HashMap<String, IndividualDatum> dataMap;
	PopulationData populationData;
	private ArrayList<String> indivList;
	private ArrayList<String> outlierList;
	private HashMap<String, Integer> scoringMap;
	private HashSet<String> excludeList;
	
	private LRRBDevHetOutlierClassifier loadExcluded(String file, boolean project) {
		excludeList = new HashSet<String>();
		if (file == null) { return this; }
		
		if (project) {
			Project proj = new Project(file, false);
			SampleData sampleData = proj.getSampleData(0, false);
			for (String id : Array.subArray(proj.getSamples(), proj.getSamplesToExclude())) {
				for (String subID : sampleData.lookup(id)) {
					excludeList.add(subID);
				}
			}
		} else {
			excludeList.addAll(HashVec.loadFileToHashSet(file, false));
		}
		
		return this;
	}

	private LRRBDevHetOutlierClassifier loadData(int idCol, int lrrCol, int bdevCol) {
		dataMap = new HashMap<String, IndividualDatum>(10000);
		
		BufferedReader reader = Files.getReader(filename, false, true, false);
		if (reader != null) {
			try {
				String temp = reader.readLine();
				while ((temp = reader.readLine()) != null) {
					String[] line = temp.split("\t");
					if (!excludeList.contains(line[idCol])) {
						dataMap.put(line[idCol], new IndividualDatum(line[idCol], Double.parseDouble(line[lrrCol]), Double.parseDouble(line[bdevCol])));
					}
				}
				reader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			// TODO error logging?
		}
		
		return this;
	}
	
	private LRRBDevHetOutlierClassifier setPopData() {
		double[] bdev = new double[dataMap.size()];
		double[] lrr = new double[dataMap.size()];
		
		int index = 0;
		for (java.util.Map.Entry<String, IndividualDatum> entry : dataMap.entrySet()) {
			bdev[index] = entry.getValue().dataBDev;
			lrr[index] = entry.getValue().dataLRR;
			index++;
		}
		
		this.populationData = new PopulationData(Array.mean(bdev), Array.mean(lrr), Array.stdev(bdev), Array.stdev(lrr), Array.median(bdev), Array.median(lrr));
		
		return this;
	}
	
	private LRRBDevHetOutlierClassifier filterPotentialOutliers() {
		indivList = new ArrayList<String>();
		outlierList = new ArrayList<String>();
		scoringMap = new HashMap<String, Integer>();
		
		outer: for (java.util.Map.Entry<String, IndividualDatum> entry : dataMap.entrySet()) {
			for (OutlierTest test : OutlierTest.values()) {
				if (test.testIndiv(populationData, entry.getValue())) {
					outlierList.add(entry.getKey());
					scoringMap.put(entry.getKey(), test.getMethodID());
//					System.out.println("\t" + entry.getValue().id + "\t1"); 
					continue outer;
				}
			}
//			System.out.println();
			indivList.add(entry.getKey());
		}
		
		return this;
	}
	
	private LRRBDevHetOutlierClassifier writeResults(String outFile) {
		String fileRoot = outFile.startsWith("_") ? ext.rootOf(filename, false) : ext.parseDirectoryOfFile(filename);
		String file = fileRoot + (outFile == null ? "_outliers.txt" : outFile);
		
		PrintWriter writer = Files.getAppropriateWriter(file);
		writer.println("ID\tOutlier");
		for (String outlier : outlierList) {
			writer.print(outlier);
			writer.println("\t" + scoringMap.get(outlier));
		}
		for (String outlier : indivList) {
			writer.print(outlier);
			writer.println("\t0");
		}
		writer.flush();
		writer.close();
		
		return this;
	}
	
	private LRRBDevHetOutlierClassifier dispose() {
		this.filename = null;
		this.populationData = null;
		this.dataMap = null;
		this.indivList = null;
		this.outlierList = null;
		try {
			this.finalize();
		} catch (Throwable e) {
			// TODO Auto-generated catch block
		}
		return null;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "LRR_MEDIAN.xln";
		String proj = null;
		String exclFile = null;
		String outFile = null;
		int idCol = 0;
		int lrrCol = 1;
		int bdevCol = 2;
		
		String usage = "\n"
						+ "one.LRRBDevHetOutlierClassifier requires 4 arguments\n" 
						+ "   (1) filename (i.e. file=" + filename + " (default))\n"
						+ "   (2) Column index for IDs (i.e. id=" + idCol + " (default))\n"
						+ "   (3) Column index for LRR data (i.e. lrr=" + lrrCol + " (default))\n"
						+ "   (4) Column index for BDev data (i.e. bdev=" + bdevCol + " (default))\n"
						+ "       OPTIONAL:"
						+ "   (5a) (Optional) Use an existing Project to remove excluded data from analysis (i.e. proj=default.properties (not the default))\n"
						+ "       OR "
						+ "   (5b) (Optional) Excluded data from analysis based on a file of excluded sample IDs (i.e. excl=excluded.txt (not the default))\n"
						+ "\n";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("proj=")) {
				proj = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("excl=")) {
				exclFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("id=")) {
				idCol = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("lrr=")) {
				lrrCol = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("bdev=")) {
				bdevCol = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outFile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		String file = (proj == null ? exclFile : proj);
		boolean isProj = proj != null;
		String outfile = outFile == null ? "_" + idCol + "-" + lrrCol + "," + bdevCol + "_outliers.txt" : outFile;
		try {
			LRRBDevHetOutlierClassifier classifier = new LRRBDevHetOutlierClassifier(filename);
			classifier.loadExcluded(file, isProj).loadData(idCol, lrrCol, bdevCol).setPopData().filterPotentialOutliers().writeResults(outfile).dispose();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}


