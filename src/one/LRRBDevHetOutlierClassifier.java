package one;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import common.Array;
import common.Files;
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
				System.out.println("Scr_BD: " + sdCntBDev + " ;; Scr_LRR: " + sdCntLRR);
				return false;
			}
		},
		/**
		 * Test if BDev exceeds 0.1 and LRR exceeds 0.1 or -0.1
		 */
		PT_1_TEST() {
			public boolean testIndiv(PopulationData popData, IndividualDatum indivDatum) {
				if (indivDatum.dataBDev >= 0.1
						&& (indivDatum.dataLRR >= 0.1 || indivDatum.dataLRR <= -0.1)) {
					return true;
				}
				return false;
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
		};
		
		abstract boolean testIndiv(PopulationData popData, IndividualDatum indivDatum);
	}
	
	public LRRBDevHetOutlierClassifier(String file) {
		this.filename = file;
	}
	
	
	String filename;
	HashMap<String, IndividualDatum> dataMap;
	PopulationData populationData;
	
	private void loadData(int idCol, int lrrCol, int bdevCol) {
		dataMap = new HashMap<String, IndividualDatum>(10000);
		
		BufferedReader reader = Files.getReader(filename, false, true, false);
		if (reader != null) {
			try {
				String temp = reader.readLine();
				while ((temp = reader.readLine()) != null) {
					String[] line = temp.split("\t");
					dataMap.put(line[idCol], new IndividualDatum(line[idCol], Double.parseDouble(line[lrrCol]), Double.parseDouble(line[bdevCol])));
				}
				reader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			// TODO error logging?
		}
		
	}
	
	private void setPopData() {
		double[] bdev = new double[dataMap.size()];
		double[] lrr = new double[dataMap.size()];
		
		int index = 0;
		for (java.util.Map.Entry<String, IndividualDatum> entry : dataMap.entrySet()) {
			bdev[index] = entry.getValue().dataBDev;
			lrr[index] = entry.getValue().dataLRR;
			index++;
		}
		
		this.populationData = new PopulationData(Array.mean(bdev), Array.mean(lrr), Array.stdev(bdev), Array.stdev(lrr), Array.median(bdev), Array.median(lrr));
	}
	
	private ArrayList<String> getPotentialOutliers() {
		ArrayList<String> indivList = new ArrayList<String>();
		
		outer: for (java.util.Map.Entry<String, IndividualDatum> entry : dataMap.entrySet()) {
			for (OutlierTest test : OutlierTest.values()) {
				if (test.testIndiv(populationData, entry.getValue())) {
					indivList.add(entry.getKey());
					continue outer;
				}
			}
		}
		
		return indivList;
	}
	
	private void writeResults(ArrayList<String> potentialOutliers) {
		String fileRoot = ext.rootOf(filename, false);
		String file = fileRoot + "_outlier_IDs.txt";
		
		PrintWriter writer = Files.getAppropriateWriter(file);
		for (String outlier : potentialOutliers) {
			writer.println(outlier);
		}
		writer.flush();
		writer.close();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "LRR_MEDIAN.xln";

		int idCol = 0;
		int lrrCol = 1;
		int bdevCol = 2;
		
		String usage = "\n"
						+ "one.LRRBDevHetOutlierClassifier requires 4 arguments\n" 
						+ "   (1) filename (i.e. file=" + filename + " (default))\n"
						+ "   (2) Column index for IDs (i.e. id=" + idCol + " (default))\n"
						+ "   (3) Column index for LRR data (i.e. lrr=" + lrrCol + " (default))\n"
						+ "   (4) Column index for BDev data (i.e. bdev=" + bdevCol + " (default))\n";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
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
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			LRRBDevHetOutlierClassifier classifier = new LRRBDevHetOutlierClassifier(filename);
			classifier.loadData(idCol, lrrCol, bdevCol);
			classifier.setPopData();
			classifier.writeResults(classifier.getPotentialOutliers());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}


