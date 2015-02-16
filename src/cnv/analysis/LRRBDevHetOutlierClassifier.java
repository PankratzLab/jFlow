package cnv.analysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;

import cnv.filesys.Project;
import cnv.var.SampleData;
import common.Array;
import common.Files;
import common.HashVec;
import common.ext;

public class LRRBDevHetOutlierClassifier {
	
	static final class AnalysisData {
		private String lrrHdr;
		private String bdevHdr;
		private final int lrrCol;
		private final int bdevCol;
		
		final HashMap<String, IndividualData> dataMap = new HashMap<String, IndividualData>(10000);
		PopulationData populationData;

		HashSet<String> indivList = new HashSet<String>();
		HashSet<String> outlierList = new HashSet<String>();
		HashMap<String, Double> scoreMap;
		HashMap<String, Integer> validatedCodes;
		int popCount;
		private String ucscRegion;
		
		public AnalysisData(int lrrCol, int bdevCol) {
			this.lrrCol = lrrCol;
			this.bdevCol = bdevCol;
		}

		public String getLrrHdr() { return lrrHdr; }
		public void setLrrHdr(String lrrHdr) { this.lrrHdr = lrrHdr; }
		public String getBdevHdr() { return bdevHdr; }
		public void setBdevHdr(String bdevHdr) { this.bdevHdr = bdevHdr; }
		public int getLrrCol() { return lrrCol; }
		public int getBdevCol() { return bdevCol; }
		
		private void derivePopData() {
			double[] bdev = new double[dataMap.size()];
			double[] lrr = new double[dataMap.size()];
			
			int index = 0;
			for (java.util.Map.Entry<String, IndividualData> entry : dataMap.entrySet()) {
				bdev[index] = entry.getValue().dataBDev;
				lrr[index] = entry.getValue().dataLRR;
				index++;
			}
			
			this.populationData = new PopulationData(Array.mean(bdev), Array.mean(lrr), Array.stdev(bdev), Array.stdev(lrr), Array.median(bdev), Array.median(lrr));
		}

		public void addIndivData(String idStr, IndividualData individualData) {
			this.dataMap.put(idStr, individualData);
		}

		public void setUCSC(String ucsc) {
			this.ucscRegion = ucsc;
		}
		
		
	}
	
	static final class IndividualData {
		final String id;
		final double dataLRR;
		final double dataBDev;
		
		public IndividualData(String id, double lrr, double bdev) {
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
	
	
	enum OutlierClassifier {
		WEIGHTED_BDEV_SD_EUCLIDEAN_NORM() {
			@Override
			String getClassifierName() { return "StdDev BDev-Weighted Euclidean Norm"; }
			@Override
			void runClassifier(final AnalysisData analysis) {
				scoreEuclidean(analysis, this);
			}
		},
		WEIGHTED_LRR_SD_EUCLIDEAN_NORM() {
			@Override
			String getClassifierName() { return "StdDev LRR-Weighted Euclidean Norm"; }
			@Override
			void runClassifier(final AnalysisData analysis) {
				scoreEuclidean(analysis, this);
			}
		},
		UNWEIGHTED_SD_EUCLIDEAN_NORM() {
			@Override
			String getClassifierName() { return "StdDev Unweighted Euclidean Norm"; }
			@Override
			void runClassifier(final AnalysisData analysis) {
				scoreEuclidean(analysis, this);
			}
		};
		
		abstract void runClassifier(AnalysisData analysis);
		abstract String getClassifierName();
		
		private static void scoreEuclidean(final AnalysisData analysis, OutlierClassifier classifier) {
			analysis.scoreMap = new HashMap<String, Double>();
			final TreeSet<String> scoreList = new TreeSet<String>(new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					return analysis.scoreMap.get(o1).compareTo(analysis.scoreMap.get(o2));
				}
			});
			double[] scores = new double[analysis.dataMap.size()]; 
			
			int index = 0;
			for (java.util.Map.Entry<String, IndividualData> entry : analysis.dataMap.entrySet()) {
				IndividualData indivDatum = entry.getValue();
				double sdCntBDev = (indivDatum.dataBDev - analysis.populationData.meanBDev) / analysis.populationData.stdDevBDev;
				double sdCntLRR = (indivDatum.dataLRR - analysis.populationData.meanLRR) / analysis.populationData.stdDevLRR;
				
				double scr = 0.0;
				switch (classifier) {
					case WEIGHTED_BDEV_SD_EUCLIDEAN_NORM:
						scr = Math.sqrt((Math.abs(sdCntBDev * sdCntBDev * sdCntBDev)) + (sdCntLRR * sdCntLRR));
						break;
					case UNWEIGHTED_SD_EUCLIDEAN_NORM:
						scr = Math.sqrt((sdCntBDev * sdCntBDev) + (sdCntLRR * sdCntLRR));
						break;
					case WEIGHTED_LRR_SD_EUCLIDEAN_NORM:
						scr = Math.sqrt((Math.abs(sdCntLRR * sdCntLRR * sdCntLRR)) + (sdCntBDev * sdCntBDev));
						break;
				}
				
				analysis.scoreMap.put(indivDatum.id, Double.valueOf(scr));
				scoreList.add(indivDatum.id);
				scores[index] = scr;
				index++;
			}
			
			String indivLow = scoreList.pollFirst(); // remove
			String indivHigher = scoreList.first(); // look at
			double scrLow = analysis.scoreMap.get(indivLow);
			double scrHigher = analysis.scoreMap.get(indivHigher);
			while(scrHigher - scrLow < 1) {
				indivLow = scoreList.pollFirst();
				if (scoreList.size() == 0) break;
				indivHigher = scoreList.first();
				scrLow = analysis.scoreMap.get(indivLow);
				scrHigher = analysis.scoreMap.get(indivHigher);
			}
	
			analysis.outlierList = new HashSet<String>();
			analysis.outlierList.addAll(scoreList);
			
			analysis.indivList = new HashSet<String>();
			for (java.util.Map.Entry<String, IndividualData> entry : analysis.dataMap.entrySet()) {
				if (!analysis.outlierList.contains(entry.getKey())) {
					analysis.indivList.add(entry.getKey());
				}
			}
		}
	}
	
//	enum OutlierTest {
//		
////		SD_SCORE() {
////			
////			HashMap<String, Double> scoreMap = new HashMap<String, Double>();
////			TreeSet<String> scoreList = new TreeSet<String>(new Comparator<String>() {
////				@Override
////				public int compare(String o1, String o2) {
////					return scoreMap.get(o1).compareTo(scoreMap.get(o2));
////				}
////			});
////			
////			public boolean testIndiv(PopulationData popData, IndividualDatum indivDatum) {
////				double sdCntBDev = (indivDatum.dataBDev - popData.meanBDev) / popData.stdDevBDev;
////				double sdCntLRR = (indivDatum.dataLRR - popData.meanLRR) / popData.stdDevLRR;
////				
////				double scr = Math.sqrt( (2 * (sdCntBDev * sdCntBDev)) + sdCntLRR * sdCntLRR);
////				scoreMap.put(indivDatum.id, Double.valueOf(scr));
////				scoreList.add(indivDatum.id);
////				System.out.println(indivDatum.id + "\t" + scr);
////				
//				// 4.0 derived through examination of results from other measures - outliers classified by other measures all exceed 4.0 on this test. 
////				return sdCntBDev > 4.0 || sdCntLRR > 4.0 || sdCntBDev < -4.0 || sdCntLRR < -4.0;
////				
////				return scr > 8.0;
////			}
////			public int getMethodID() {
////				return 1;
////			}
////		},
//		/**
//		 * Test if BDev exceeds 0.1 and LRR exceeds 0.1 or -0.1
//		 */
//		PT_1_TEST() {
//			public boolean testIndiv(PopulationData popData, IndividualData indivDatum) {
////				System.out.print("\t" + indivDatum.id + "\t");
//				if (indivDatum.dataBDev >= 0.1
//						&& (indivDatum.dataLRR >= 0.1 || indivDatum.dataLRR <= -0.1)) {
//					return true;
//				}
//				return false;
//			}
//			public int getMethodID() {
//				return 1;
//			}
//		},
//		/**
//		 * Test if BDev or LRR exceed more than a certain number of SD's (default 4) from the mean
//		 */
//		EXTREME_SD() {
//			private static final int EXTREME_SD_LIMIT = 4;
//			public boolean testIndiv(PopulationData popData, IndividualData indivDatum) {
//				if ((indivDatum.dataBDev > (popData.meanBDev + (popData.stdDevBDev * EXTREME_SD_LIMIT)))
//						|| (indivDatum.dataBDev < (popData.meanBDev - (popData.stdDevBDev * EXTREME_SD_LIMIT)))
//						|| (indivDatum.dataLRR > (popData.meanLRR + (popData.stdDevLRR * EXTREME_SD_LIMIT)))
//						|| (indivDatum.dataLRR < (popData.meanLRR - (popData.stdDevLRR * EXTREME_SD_LIMIT)))) {
//					return true;
//				}
//				return false;
//			}
//			public int getMethodID() {
//				return 2;
//			}
//		},
//		/**
//		 * Test if BDev exceeds 0.1 and LRR is exceeds more than 2 SD's from the mean
//		 */
//		BD_PT1_LRRSD() {
//			public boolean testIndiv(PopulationData popData, IndividualData indivDatum) {
//				if (indivDatum.dataBDev >= 0.1
//						&& (indivDatum.dataLRR >= (popData.meanLRR + (popData.stdDevLRR * 2)) || indivDatum.dataLRR <= (popData.meanLRR - (popData.stdDevLRR * 2)))) {
//					return true;
//				}
//				return false;
//			}
//			public int getMethodID() {
//				return 3;
//			}
//		};
//		
//		abstract boolean testIndiv(PopulationData popData, IndividualData indivDatum);
//		abstract int getMethodID();
//	}
	
	public LRRBDevHetOutlierClassifier(String file, boolean isFileRoot) {
		this.filename = file;
		this.isFileRoot = isFileRoot;
	}
	
	String filename;
	boolean isFileRoot;
	private HashSet<String> excludeList;
	private HashSet<String> popList;
	ArrayList<AnalysisData> analyses = new ArrayList<LRRBDevHetOutlierClassifier.AnalysisData>();
	
	private LRRBDevHetOutlierClassifier loadExcluded(String file, boolean project) {
		excludeList = new HashSet<String>();
		popList = new HashSet<String>();
		if (file == null) { return this; }
		
		if (project) {
			Project proj = new Project(file, false);
			SampleData sampleData = proj.getSampleData(0, false);
			for (String id : Array.subArray(proj.getSamples(), proj.getSamplesToExclude())) {
				for (String subID : sampleData.lookup(id)) {
					excludeList.add(subID);
				}
			}
			for (String id : Array.subArray(proj.getSamples(), Array.booleanNegative(proj.getSamplesToExclude()))) {
				popList.add(sampleData.lookup(id)[0]);
			}
		} else {
			excludeList.addAll(HashVec.loadFileToHashSet(file, false));
		}
		
		return this;
	}
	
	private LRRBDevHetOutlierClassifier loadAllFilesAllData(int idCol, String lrrHdr, String bdevHdr) {
		String dir = ext.parseDirectoryOfFile(filename);
		String prefix = filename.substring(dir.length());
		String[] files = Files.list(dir, prefix, null, true, false);
		
		for (String file : files) {
			ArrayList<AnalysisData> fileData = new ArrayList<LRRBDevHetOutlierClassifier.AnalysisData>();
			BufferedReader reader = Files.getReader(dir + file, false, true, false);
			if (reader != null) {
				String temp;
				try {
					temp = reader.readLine();
					String[] hdr = temp.split("\t");
					
					int[] colSet = new int[]{-1, -1};
					for (int i = 0; i < hdr.length; i++) {
						if (hdr[i].contains(lrrHdr)) {
							if (colSet[0] != -1) {
								// TODO ERROR, no BDev column found between LRR columns
							} else {
								colSet[0] = i;
							}
						} else if (hdr[i].contains(bdevHdr)) {
							if (colSet[0] == -1) {
								// TODO ERROR, no LRR column found before BDev column
							} else if (colSet[1] != -1) {
								// TODO ERROR, additional BDev column found before new LRR column
							} else {
								colSet[1] = i;
								AnalysisData analysis = new AnalysisData(colSet[0], colSet[1]);
								analysis.setLrrHdr(hdr[colSet[0]]);
								analysis.setBdevHdr(hdr[colSet[1]]);
								String ucsc = "chr" + hdr[colSet[0]].split("chr")[1];
								analysis.setUCSC(ucsc);
								fileData.add(analysis);
								colSet = new int[]{-1, -1};
							}
						}
					}
					
					while ((temp = reader.readLine()) != null) {
						String[] line = temp.split("\t");
						if (!excludeList.contains(line[idCol])) {
							String id = line[idCol];
							for (AnalysisData analysis : fileData) {
								analysis.addIndivData(id, new IndividualData(id, Double.parseDouble(line[analysis.lrrCol]), Double.parseDouble(line[analysis.bdevCol])));
							}
						}
					}
					reader.close();
					
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				
				analyses.addAll(fileData);
			}
		}
		
		for (AnalysisData analysis : analyses) {
			analysis.popCount = popList.size();
			analysis.derivePopData();
		}
		
		return this;
	}
	
	private LRRBDevHetOutlierClassifier loadAllData(int idCol, String lrrHdr, String bdevHdr) {
		if (isFileRoot) {
			return loadAllFilesAllData(idCol, lrrHdr, bdevHdr);
		}
		BufferedReader reader = Files.getReader(filename, false, true, false);
		if (reader != null) {
			String temp;
			try {
				temp = reader.readLine();
				String[] hdr = temp.split("\t");
				
				int[] colSet = new int[]{-1, -1};
				for (int i = 0; i < hdr.length; i++) {
					if (hdr[i].contains(lrrHdr)) {
						if (colSet[0] != -1) {
							// TODO ERROR, no BDev column found between LRR columns
						} else {
							colSet[0] = i;
						}
					} else if (hdr[i].contains(bdevHdr)) {
						if (colSet[0] == -1) {
							// TODO ERROR, no LRR column found before BDev column
						} else if (colSet[1] != -1) {
							// TODO ERROR, additional BDev column found before new LRR column
						} else {
							colSet[1] = i;
							AnalysisData analysis = new AnalysisData(colSet[0], colSet[1]);
							analysis.setLrrHdr(hdr[colSet[0]]);
							analysis.setBdevHdr(hdr[colSet[1]]);
							String ucsc = "chr" + hdr[colSet[0]].split("chr")[1];
							analysis.setUCSC(ucsc);
							analyses.add(analysis);
							colSet = new int[]{-1, -1};
						}
					}
				}
				
				while ((temp = reader.readLine()) != null) {
					String[] line = temp.split("\t");
					if (!excludeList.contains(line[idCol])) {
						String id = line[idCol];
						for (AnalysisData analysis : analyses) {
							analysis.addIndivData(id, new IndividualData(id, Double.parseDouble(line[analysis.lrrCol]), Double.parseDouble(line[analysis.bdevCol])));
						}
					}
				}
				reader.close();
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		for (AnalysisData analysis : analyses) {
			analysis.popCount = popList.size();
			analysis.derivePopData();
		}
		
		return this;
	}
	
	private LRRBDevHetOutlierClassifier loadData(int idCol, int lrrCol, int bdevCol) {
		AnalysisData analysis = new AnalysisData(lrrCol, bdevCol);
		
		BufferedReader reader = Files.getReader(filename, false, true, false);
		if (reader != null) {
			try {
				String temp = reader.readLine();
				String[] hdr = temp.split("\t");
				analysis.setLrrHdr(hdr[lrrCol]);
				analysis.setBdevHdr(hdr[bdevCol]);
				String ucsc = "chr" + hdr[lrrCol].split("chr")[1];
				analysis.setUCSC(ucsc);
				
				int popCnt = 0;
				while ((temp = reader.readLine()) != null) {
					String[] line = temp.split("\t");
					if (!excludeList.contains(line[idCol])) {
						analysis.addIndivData(line[idCol], new IndividualData(line[idCol], Double.parseDouble(line[lrrCol]), Double.parseDouble(line[bdevCol])));
						popCnt++;
					}
				}
				reader.close();
				analysis.popCount = popCnt;
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			// TODO error logging?
		}
		
		analysis.derivePopData();
		analyses.add(analysis);
		return this;
	}
	
	private LRRBDevHetOutlierClassifier scoreOutliers() {
		for (final AnalysisData analysis : analyses) {
			OutlierClassifier.WEIGHTED_BDEV_SD_EUCLIDEAN_NORM.runClassifier(analysis);
//			OutlierClassifier.UNWEIGHTED_SD_EUCLIDEAN_NORM.runClassifier(analysis);
//			OutlierClassifier.WEIGHTED_LRR_SD_EUCLIDEAN_NORM.runClassifier(analysis);
		}
		
		return this;
	}
	
	/**
	 * Validation file must be in format: {ID ... UCSC ... 0/1 ... 0-5}
	 * 
	 * @param validationFile
	 * @return
	 * @throws IOException
	 */
	private LRRBDevHetOutlierClassifier validateResults(String validationFile) throws IOException {
		if (validationFile == null) return this;
		HashSet<String> locationsInValidationSet = new HashSet<String>();
		HashMap<String, HashSet<String>> outliersInRegionsMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashMap<String, Integer>> regionToIndivDetailMap = new HashMap<String, HashMap<String, Integer>>();
		
		BufferedReader reader = Files.getAppropriateReader(validationFile);
		String line = reader.readLine();
		while((line = reader.readLine()) != null) {
			String[] temp = line.split("\t");
			String id = temp[0];
			String ucsc = temp[1];
			String outlierStr = temp[2];
			String detailStr = temp[3];
			
			locationsInValidationSet.add(ucsc);
			
			HashSet<String> outliers = outliersInRegionsMap.get(ucsc);
			if (outliers == null) {
				outliers = new HashSet<String>();
				outliersInRegionsMap.put(ucsc, outliers);
			}
			if (outlierStr.equals("1")) outliers.add(id);
			
			HashMap<String, Integer> indivDetails = regionToIndivDetailMap.get(ucsc);
			if (indivDetails == null) {
				indivDetails = new HashMap<String, Integer>();
				regionToIndivDetailMap.put(ucsc, indivDetails);
			}
			indivDetails.put(id, Integer.parseInt(detailStr));
			
		}
		reader.close();

		int[] totalStats = {0, 0, 0, 0};
		// TP, FP, FN, TN
		HashMap<String, int[]> regionStats = new HashMap<String, int[]>();
		
		for (AnalysisData analysis : analyses) {
			int[] stats = {0, 0, 0, 0};
			
			HashSet<String> validOutliers = outliersInRegionsMap.get(analysis.ucscRegion);
			if (validOutliers != null) {
				for (String id : analysis.outlierList) {
					if (validOutliers.contains(id)) {
						stats[0]++;
						totalStats[0]++;
					} else {
						stats[1]++;
						totalStats[1]++;
					}
					validOutliers.remove(id);
				}
				// FN
				stats[2] = validOutliers.size();
				totalStats[2] += validOutliers.size();
			}
			// TN
			stats[3] = analysis.popCount - stats[0] - stats[1] - stats[2];
			totalStats[3] += stats[3];
			
			regionStats.put(analysis.ucscRegion, stats);
			
			HashMap<String, Integer> indivDetails = regionToIndivDetailMap.get(analysis.ucscRegion);
			if (indivDetails == null) continue;
			analysis.validatedCodes = new HashMap<String, Integer>();
			for (String id : analysis.dataMap.keySet()) {
				Integer code = indivDetails.get(id);
				if (code == null) {
					code = 0;
				}
				analysis.validatedCodes.put(id, code);
			}
		}
		
	
		PrintWriter writer = Files.getAppropriateWriter(ext.rootOf(filename, false) + "_validation.txt");
		String header = "UCSC\tTP\tFP\tFN\tTN\tScore";
		writer.println(header);
		for (java.util.Map.Entry<String, int[]> regionStat : regionStats.entrySet()) {
			int[] stats = regionStat.getValue();
			double score = 0.0;
			// balanced accuracy = (0.5*TP \ (TP + FN)) + (0.5*TN \ (TN + FP)) 
			score = ((stats[0] + stats[2]) == 0 ? 0 : ((0.5 * stats[0])/(stats[0] + stats[2]))) + ((0.5 * stats[3])/(stats[3] + stats[1]));
			writer.println(regionStat.getKey() + "\t" + stats[0] + "\t" + stats[1] + "\t" + stats[2] + "\t" + stats[3] + "\t" + score);
		}
		writer.println("\t\t" + totalStats[0] + "\t" + totalStats[1] + "\t" + totalStats[2] + "\t" + totalStats[3]);
		
		writer.flush();
		writer.close();
		
		return this;
	}

	private static void validateResults(String validationFile, String resultsFile, String outFile) throws IOException {
		HashSet<String> locationsInValidationSet = new HashSet<String>();
		HashMap<String, HashSet<String>> outliersInRegionsMap = new HashMap<String, HashSet<String>>();
		HashMap<String, HashMap<String, Integer>> regionToIndivDetailMap = new HashMap<String, HashMap<String, Integer>>();
		
		BufferedReader reader = Files.getAppropriateReader(validationFile);
		String line = reader.readLine();
		while((line = reader.readLine()) != null) {
			String[] temp = line.split("\t");
			String id = temp[0];
			String ucsc = temp[1];
			String outlierStr = temp[2];
			String detailStr = temp[3];
			
			locationsInValidationSet.add(ucsc);
			
			HashSet<String> outliers = outliersInRegionsMap.get(ucsc);
			if (outliers == null) {
				outliers = new HashSet<String>();
				outliersInRegionsMap.put(ucsc, outliers);
			}
			if (outlierStr.equals("1")) outliers.add(id);
			
			HashMap<String, Integer> indivDetails = regionToIndivDetailMap.get(ucsc);
			if (indivDetails == null) {
				indivDetails = new HashMap<String, Integer>();
				regionToIndivDetailMap.put(ucsc, indivDetails);
			}
			indivDetails.put(id, Integer.parseInt(detailStr));
			
		}
		reader.close();
		
		HashMap<String, HashSet<String>> calledOutliers = new HashMap<String, HashSet<String>>();
		
		reader = Files.getAppropriateReader(resultsFile);
		line = reader.readLine();
		int popCnt = Integer.valueOf(line.split("\t")[2].split("=")[1]);
		while((line = reader.readLine()) != null) {
			String[] temp = line.split("\t");
			String id = temp[0];
			String ucsc = temp[1];
			HashSet<String> regionOutliers = calledOutliers.get(ucsc);
			if (regionOutliers == null) {
				regionOutliers = new HashSet<String>();
				calledOutliers.put(ucsc, regionOutliers);
			}
			regionOutliers.add(id);
		}
		reader.close();
		
		
		// TP, FP, FN, TN
		HashMap<String, int[]> regionStats = new HashMap<String, int[]>();
		int[] totalStats = {0, 0, 0, 0};
		
		for (java.util.Map.Entry<String, HashSet<String>> regionOutliers : calledOutliers.entrySet()) {
			int[] stats = {0, 0, 0, 0};
			
			HashSet<String> validOutliers = outliersInRegionsMap.get(regionOutliers.getKey());
			if (validOutliers != null) {
				for (String id : regionOutliers.getValue()) {
					if (validOutliers.contains(id)) {
						stats[0]++;
						totalStats[0]++;
					} else {
						stats[1]++;
						totalStats[1]++;
					}
					validOutliers.remove(id);
				}
				// FN
				stats[2] = validOutliers.size();
				totalStats[2] += validOutliers.size();
			}
			// TN
			stats[3] = popCnt - stats[0] - stats[1] - stats[2];
			totalStats[3] += stats[3];
			
			regionStats.put(regionOutliers.getKey(), stats);
		}
	
		PrintWriter writer = Files.getAppropriateWriter(outFile);
		String header = "UCSC\tTP\tFP\tFN\tTN\tScore";
		writer.println(header);
		for (java.util.Map.Entry<String, int[]> regionStat : regionStats.entrySet()) {
			int[] stats = regionStat.getValue();
			double score = 0.0;
			// balanced accuracy = (0.5*TP \ (TP + FN)) + (0.5*TN \ (TN + FP)) 
			score = ((stats[0] + stats[2]) == 0 ? 0 : ((0.5 * stats[0])/(stats[0] + stats[2]))) + ((0.5 * stats[3])/(stats[3] + stats[1]));
			writer.println(regionStat.getKey() + "\t" + stats[0] + "\t" + stats[1] + "\t" + stats[2] + "\t" + stats[3] + "\t" + score);
		}
		writer.println("\t\t" + totalStats[0] + "\t" + totalStats[1] + "\t" + totalStats[2] + "\t" + totalStats[3]);
		
		writer.flush();
		writer.close();
	}

	private LRRBDevHetOutlierClassifier writeResults(String outFile, boolean writeOutliersOnly, boolean recodeOutliers, boolean writeScores, boolean writeHeader, boolean combine) {
		String fileRoot = ext.parseDirectoryOfFile(filename);
		String file = outFile == null ? ext.rootOf(filename, true) + "_outliers.txt" : outFile;
		PrintWriter writer = null;
		if (combine) {
			file = fileRoot + ext.replaceWithLinuxSafeCharacters(file, false);
			writer = Files.getAppropriateWriter(file);	
			if (writeOutliersOnly) {
				if (writeHeader) {
					StringBuilder hdrStr = new StringBuilder("ID\tUCSC");
					if (writeScores) {
						hdrStr.append("\tScores");
					}
					hdrStr.append("\tN=" + analyses.get(0).popCount);
					hdrStr.append("\n");
					writer.print(hdrStr.toString());
				}

				for (int i = 0; i < analyses.size(); i++) {
					for (String id : analyses.get(i).outlierList) {
						StringBuilder resultLine = new StringBuilder(id).append("\t").append(analyses.get(i).ucscRegion);
						if (writeScores) {
							resultLine.append("\t").append(analyses.get(i).scoreMap.get(id));
						}
						resultLine.append("\n");
						writer.print(resultLine.toString());
					}
				}
			} else {
				if (writeHeader) {
					StringBuilder hdrStr = new StringBuilder("ID");
					for (int i = 0; i < analyses.size(); i++) {
						hdrStr.append("\tOutliers_").append(analyses.get(i).ucscRegion);
						if (writeScores) {
							hdrStr.append("\tScores_").append(analyses.get(i).ucscRegion);
						}
					}
					hdrStr.append("\n");
					writer.print(hdrStr.toString());
				}

				for (String dnaID : popList) {
					StringBuilder resultLine = new StringBuilder(dnaID);
					for (int i = 0; i < analyses.size(); i++) {
						if (recodeOutliers && analyses.get(i).validatedCodes != null) {
							resultLine.append("\t").append(analyses.get(i).validatedCodes.get(dnaID));
						} else {
							resultLine.append(analyses.get(i).outlierList.contains(dnaID) ? "\t1" : "\t0");
						}
						if (writeScores) {
							resultLine.append("\t").append(analyses.get(i).scoreMap.get(dnaID));
						}
					}
					writer.println(resultLine.toString());
				}
			}

			writer.flush();
			writer.close();
		} else {
			for (AnalysisData analysis : analyses) {
				file = outFile == null ? analysis.ucscRegion + "_outliers.txt" : outFile;
				file = fileRoot + ext.replaceWithLinuxSafeCharacters(file, false);
				writer = Files.getAppropriateWriter(file);
				if (writeHeader) {
					writer.println("ID\tOutliers" + (writeScores ? "\tScore" : "") + "\tN=" + analyses.get(0).popCount);
				}

				for (String outlier : analysis.outlierList) {
					writer.print(outlier);
					if (recodeOutliers) {
						writer.print("\t"+analysis.validatedCodes.get(outlier));
					} else {
						writer.print("\t1");
					}
					if (writeScores) {
						writer.print("\t" + analysis.scoreMap.get(outlier));
					}
					writer.println();
				}
				if (!writeOutliersOnly) {
					for (String outlier : analysis.indivList) {
						writer.print(outlier);
						if (recodeOutliers) {
							writer.print("\t"+analysis.validatedCodes.get(outlier));
						} else {
							writer.print("\t0");
						}
						if (writeScores) {
							writer.print("\t" + analysis.scoreMap.get(outlier));
						}
						writer.println();
					}
				}
				
				writer.flush();
				writer.close();
			}
		}
		
		return this;
	}
	
	private LRRBDevHetOutlierClassifier dispose() {
		this.filename = null;
		this.excludeList = null;
		this.analyses = null;
		try {
			this.finalize();
		} catch (Throwable e) {
			// TODO Auto-generated catch block
		}
		return null;
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "LRR_MEDIAN_1.xln";
		String fileRoot = null;
		String validationFilename = null;
		String proj = null;
		String exclFile = null;
		String outFile = null;
		int idCol = 0;
		int lrrCol = 1;
		int bdevCol = 2;
		String lrrHdr = null;
		String bdevHdr = null;
		boolean outliersOnly = false;
		boolean writeScores = false;
		boolean noHeader = false;
		boolean splitResults = false;
		boolean runValidation = false;
		boolean validationOnly = false;
		
		String usage = "\n"
						+ "one.LRRBDevHetOutlierClassifier requires 2-4+ arguments\n" 
						+ "   (1a) filename (i.e. file=" + filename + " (default))\n"
						+ "     OR\n"
						+ "   (1b) root of multiple data files (i.e. fileRoot=" + fileRoot + " (default))" 
						+ "     AND\n"
						+ "   (2) Column index for IDs (i.e. id=" + idCol + " (default))\n"
						+ "   (3a) Column index for LRR data (i.e. lrr=" + lrrCol + " (default))\n"
						+ "   (4a) Column index for BDev data (i.e. bdev=" + bdevCol + " (default))\n"
						+ "     OR\n"
						+ "   (3b) Unique column header fragment for all columns of LRR data (i.e. lrrHdr=MEDIAN (not the default))"
						+ "   (4b) Unique column header fragment for all columns of BDev data (i.e. bdevHdr=BDeviation_Het (not the default))\n"	
						+ "       OPTIONAL:"
						+ "   (5a) (Optional) Use an existing Project to remove excluded data from analysis (i.e. proj=default.properties (not the default))\n"
						+ "     OR "
						+ "   (5b) (Optional) Excluded data from analysis based on a file of excluded sample IDs (i.e. excl=excluded.txt (not the default))\n"
						+ "   (6) filename of validation set (i.e. vfile=outliers.txt (not the default))\n"
						+ "   (7) Run validation and write results, do no write outliers or scores (i.e. -validateOnly (not the default))\n"
						+ "  OR  \n"
						+ "   (1) filename of outlier results (i.e. file=" + filename + " (default))\n"
						+ "   (2) filename of validation set (i.e. vfile=outliers.txt (not the default))\n"
						+ "   (3) run validation (i.e. -validate (not the default))\n"
						+ "\n";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("fileRoot=")) {
				fileRoot = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("vfile=")) {
				validationFilename = args[i].split("=")[1];
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
			} else if (args[i].startsWith("lrrHdr=")) {
				lrrHdr = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bdevHdr=")) {
				bdevHdr = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].equals("-validate")) {
				runValidation = true;
				numArgs--;
			} else if (args[i].equals("-validateOnly")) {
				validationOnly = true;
				numArgs--;
			} else if (args[i].equals("-scores")) {
				writeScores = true;
				numArgs--;
			} else if (args[i].equals("-outliers")) {
				outliersOnly = true;
				numArgs--;
			} else if (args[i].equals("-noHeader")) {
				noHeader = true;
				numArgs--;
			} else if (args[i].equals("-splitResults")) {
				splitResults = true;
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		if (runValidation) {
			try {
				validateResults(validationFilename, filename, ext.rootOf(filename, false) + "_validation.txt");
			} catch (IOException e) {
				e.printStackTrace();
			}
			return;
		}
		
		String file = (proj == null ? exclFile : proj);
		boolean isProj = proj != null;
		try {
			String f = fileRoot == null ? filename : fileRoot;
			boolean isRoot = fileRoot != null;
			LRRBDevHetOutlierClassifier classifier = new LRRBDevHetOutlierClassifier(f, isRoot).loadExcluded(file, isProj);
			if (lrrHdr != null && bdevHdr != null) {
				classifier.loadAllData(idCol, lrrHdr, bdevHdr);
			} else {
				classifier.loadData(idCol, lrrCol, bdevCol);
			}
			classifier.scoreOutliers().validateResults(validationFilename);
			if (!validationOnly) {
				classifier.writeResults(outFile, outliersOnly, validationFilename != null, writeScores, !noHeader, !splitResults);
			}
			classifier.dispose();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}




