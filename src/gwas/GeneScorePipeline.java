package gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import stats.Histogram;
import stats.ProbDist;
import stats.RegressionModel;
import bioinformatics.MapSNPsAndGenes;
import bioinformatics.Sequence;	
import common.Aliases;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Positions;
import common.ext;

public class GeneScorePipeline {
	
	private static float DEFAULT_INDEX_THRESHOLD = (float)0.00000005;
	private static int DEFAULT_WINDOW_MIN_SIZE_PER_SIDE = 500000;// 500kb each side is technically a 1M window until the next hit region, but we now take this into consideration in the main algorithm
	private static float DEFAULT_WINDOW_EXTENSION_THRESHOLD = (float)0.00000005; // (float)0.00001;
	private static String[] DEFAULT_ADDL_ANNOT_VAR_NAMES = new String[0];
	
	private static final String PLINK_FRQ_DIR = "D:/1000G/";
	private static final String TAG = "##";
	private static final String PLINK_FRQ_FILE_PATTERN = "chr" + TAG + "_eu_unrel.frq.xln";
	private static final int ACCEPTABLE_SEPARATION = 0;
	
	private static final String[][] LINKERS = {
		Aliases.MARKER_NAMES,
		Aliases.ALLELES[0],
		Aliases.EFFECTS 
	};
	
	private static final String REGRESSION_HEADER = "STUDY\tDATAFILE\tINDEX-THRESHOLD\tFACTOR\tBASE-R-SQR\tR-SQR\tSIG\tBETA\tSE\tNUM\t#DATASNPs\t#PLINKSNPs\t#HITSNPs\tB-F-SCORE\tINVCHI-SCORE\tP-VALUE";
	private String metaDir;
	
	private float[] indexThresholds = new float[]{DEFAULT_INDEX_THRESHOLD};
	private int[] windowMinSizePerSides = new int[]{DEFAULT_WINDOW_MIN_SIZE_PER_SIDE};
	private float[] windowExtensionThresholds = new float[]{DEFAULT_WINDOW_EXTENSION_THRESHOLD};
	
//	private int numThreads = 1;
//	private boolean runPlink = false;
//	private boolean runRegression = false;
//	private boolean writeHist = false;
	
	private ArrayList<String> dataFiles = new ArrayList<String>();
	private ArrayList<Study> studies = new ArrayList<GeneScorePipeline.Study>();
	private HashMap<String, Constraint> analysisConstraints = new HashMap<String, GeneScorePipeline.Constraint>();
	private HashMap<String, HashMap<String, Integer>> dataCounts = new HashMap<String, HashMap<String,Integer>>();
	
	private int bimChrIndex = 0;
	private int bimMkrIndex = 1;
	private int bimPosIndex = 3;
	private int bimA1Index = 4;
	private int bimA2Index = 5;

	private int hitsMkrIndex = 1;
	
	private class Study {
		String studyDir;
		String studyName;
		
		String plinkPref = "plink";
		
		ArrayList<String> phenoFiles = new ArrayList<String>();
		
		HashMap<String, PhenoData> phenoData = new HashMap<String, GeneScorePipeline.PhenoData>();
		
		//     constraint  ->  datafile   ->   phenofile
		HashMap<String, HashMap<String, HashMap<String, RegressionResult>>> regressions = new HashMap<String, HashMap<String, HashMap<String, RegressionResult>>>();
		//     constraint  ->  datafile   
		HashMap<String, HashMap<String, double[]>> scores = new HashMap<String, HashMap<String, double[]>>();
		//     constraint  ->  datafile   
		HashMap<String, HashMap<String, Integer>> hitSnpCounts = new HashMap<String, HashMap<String, Integer>>();
		//     constraint  ->  datafile
		HashMap<String, HashMap<String, Integer>> hitWindowCnts = new HashMap<String, HashMap<String,Integer>>();
//	     constraint  ->  datafile
		HashMap<String, HashMap<String, Integer>> dataCounts = new HashMap<String, HashMap<String,Integer>>();
	}
	
	private class RegressionResult {
		boolean logistic;
		double rsq;
		double baseRSq;
		double pval;
		double beta;
		double se;
		int num;
		public double stats;
	}
	
	private class Constraint {
		final float indexThreshold;
		final int windowMinSizePerSide;
		final float windowExtensionThreshold;
		public Constraint(float i, int m, float w) {
			this.indexThreshold = i;
			this.windowMinSizePerSide = m;
			this.windowExtensionThreshold = w;
		}
	}
	
	private class PhenoData {
		String phenoName;
		HashMap<String, PhenoIndiv> indivs = new HashMap<String, GeneScorePipeline.PhenoIndiv>();
		ArrayList<String> covars = new ArrayList<String>();
	}
	
	private class PhenoIndiv {
		String fid, iid;
		double depvar;
		HashMap<String, Double> covars = new HashMap<String, Double>();
	}
	
	private static String[] readPlinkFile(int chr, String file, HashMap<String, String> rsToFull) {
		String[] rsPos = new String[Positions.CHROMOSOME_LENGTHS_MAX[chr]];
		try {
			BufferedReader reader = Files.getAppropriateReader(file);
			String temp = reader.readLine();
			String[] line;
			while((temp = reader.readLine()) != null) {
				line = temp.split("[\\s]+");
				String[] parts = line[1].split(":");
				String rsOrChr = parts[0];
				String posStr = parts[1];
				// TODO incorporate allele examination
//				String a1 = parts[2];
//				String a2 = parts[3];
				String freqStr = line[4];
				int pos = Integer.parseInt(posStr);
				rsPos[pos] = (rsOrChr.startsWith("rs")) ? rsOrChr : line[1];
				rsToFull.put((rsOrChr.startsWith("rs")) ? rsOrChr : line[1], line[1] + "\t" + freqStr);
			}
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return rsPos;
	}
	
	private static HashMap<String, Double> get1000GFreq(HashMap<String, int[]> markerMap) {
//		BufferedReader reader = null;
		String currFile, /*temp,*/ plinkFile;
//		String[] line;
		HashMap<Integer, ArrayList<String>> mkrsByChr, nonRSByChr;
		ArrayList<String> chrMkrs, nonRS;
		HashMap<String, Double> plinkFreqs = new HashMap<String, Double>();
		
		mkrsByChr = new HashMap<Integer, ArrayList<String>>();
		nonRSByChr = new HashMap<Integer, ArrayList<String>>();
		for (java.util.Map.Entry<String, int[]> entry : markerMap.entrySet()) {
			chrMkrs = mkrsByChr.get(entry.getValue()[0]);
			nonRS = nonRSByChr.get(entry.getValue()[0]);
			if (chrMkrs == null) {
				chrMkrs = new ArrayList<String>();
				mkrsByChr.put(entry.getValue()[0], chrMkrs);
			}
			if (nonRS == null) {
				nonRS = new ArrayList<String>();
				nonRSByChr.put(entry.getValue()[0], nonRS);
			}
			chrMkrs.add(entry.getKey());
			if (!entry.getKey().startsWith("rs")) {
				nonRS.add(entry.getKey());
			}
		}
		
		plinkFile = PLINK_FRQ_DIR + PLINK_FRQ_FILE_PATTERN;
		for (Integer chr : mkrsByChr.keySet()) {
			// skip this chromosome if we don't have any markers for which we need freqencies
			if (mkrsByChr.get(chr) == null || mkrsByChr.get(chr).isEmpty()) continue; 
			currFile = plinkFile.replace(TAG, chr.toString());
			
			HashMap<String, String> rsToFullData = new HashMap<String, String>();
			String[] posRS = readPlinkFile(chr, currFile, rsToFullData);
			
			chrMkrs = mkrsByChr.get(chr);
			nonRS = nonRSByChr.get(chr);
			
			ArrayList<String> rsNotFound = new ArrayList<String>();
			
			for (String mkr : chrMkrs) {
				if (rsToFullData.containsKey(mkr)) {
					plinkFreqs.put(mkr, Double.valueOf(rsToFullData.get(mkr).split("\t")[1]));
				} else {
					rsNotFound.add(mkr);
					continue;
				}
			}
			
			nonRS.addAll(rsNotFound);
			
			for (String nonRSMkr : nonRS) {
				int pos = markerMap.get(nonRSMkr)[1];
				String mkrNm = pos > 0 && pos < posRS.length && null != posRS[pos] && !"".equals(posRS[pos]) ? posRS[pos] : "";
				plinkFreqs.put(nonRSMkr, "".equals(mkrNm) ? 0.0 : Double.valueOf(rsToFullData.get(mkrNm).split("\t")[1]));
			}
			
//			Hashtable<String, Vector<String>> plinkData = HashVec.loadFileToHashVec(currFile, 1, new int[]{4}, "\t", true, false);
//			
			
//			
//			ArrayList<String> found = new ArrayList<String>();
//			
//			for (String mkr : chrMkrs) {
//				if (plinkData.containsKey(mkr))
//			}
//			
//			try {
//				reader = Files.getAppropriateReader(currFile);
//				temp = reader.readLine();
//				while((temp = reader.readLine()) != null) {
//					line = temp.split("[\\s]+");
//					String rsOrPos = line[1];
//					String a1 = line[2];
//					String a2 = line[3];
//					String freqStr = line[4];
//					
//					if (rsOrPos.startsWith("rs")) {
//						if (mkrsByChr.get(chr).contains(rsOrPos.split(":")[0])) {
//							plinkFreqs.put(rsOrPos.split(":")[0], Double.valueOf(freqStr));
//							mkrsByChr.get(chr).remove(rsOrPos.split(":")[0]);
//						}
//					} else if (!nonRSByChr.get(chr).isEmpty()){
//						// check positions
//						int plinkPos = Integer.parseInt(rsOrPos.split(":")[1]);
//						for (String nonRSMkr : nonRSByChr.get(chr)) {
//							int mkrPos = markerMap.get(nonRSMkr)[1];
//							if (Math.abs(plinkPos - mkrPos) <= ACCEPTABLE_SEPARATION) {
//								plinkFreqs.put(nonRSMkr, Double.valueOf(freqStr));
//							}
//						}
//					}
//				}
//				reader.close();
//				reader = null;
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//				if (reader != null) {
//					try {
//						reader.close();
//					} catch (IOException e1) {
//						// TODO Auto-generated catch block
//						// e1.printStackTrace();
//					}
//					reader = null;
//				}
//			}
		}
		
		return plinkFreqs;
	}
	
	public static void preprocessDataFiles(String[] files) {
		BufferedReader reader;
		String temp, delimiter;
		String[] header, snps = null, data, line;
		String[][] factors;
		int[] indices;
		Logger log = new Logger();
		HashMap<String, int[]> markerMap;
		Hashtable<String, Vector<String>> fileData;
		HashMap<String, Double> freqs;
		
		factors = new String[][]{Aliases.MARKER_NAMES, Aliases.CHRS, Aliases.POSITIONS, Aliases.PVALUES, Aliases.ALLELE_FREQS};
		for (String filename : files) {
			try {
				reader = Files.getAppropriateReader(filename);
				temp = reader.readLine();
				delimiter = ext.determineDelimiter(temp);
				header = temp.trim().split(delimiter);
				indices = ext.indexFactors(factors, header, false, false, true, true, log, false);
				markerMap = new HashMap<String, int[]>();
				String errorMsg = "";
				if (indices[0] == -1) {
					errorMsg = "ERROR - no MarkerName column found";
					// ERROR - couldn't find MarkerName column!  COMPLETE FAIL
				}
				if (indices[3] == -1) {
					// NO PVALUES!  COMPLETE FAIL
					errorMsg = errorMsg.equals("") ? "ERROR - no P-Value column found" : errorMsg + "; no P-Value column found";
				}
				if (errorMsg.equals("")) {
					if (indices[1] == -1 || indices[2] == -1) {
						// No chromosomes/positions
						snps = HashVec.loadFileToStringArray(filename, false, true, new int[]{indices[0]}, false);//fileData.keySet().toArray(new String[fileData.size()]);
						Files.writeList(snps, ext.rootOf(filename, false) + ".snps");
						MapSNPsAndGenes.procSNPsToGenes(ext.parseDirectoryOfFile(filename), ext.rootOf(filename, true) + ".snps", MapSNPsAndGenes.DEFAULT_WIGGLE_ROOM, (byte) 37, log);
						data = Array.toStringArray(HashVec.loadFileToVec(ext.rootOf(filename, false)+"_positions.xln", true, false, false));
						for (int i = 0; i < data.length; i++) {
							line = data[i].trim().split("[\\s]+");
							markerMap.put(line[0], new int[]{Positions.chromosomeNumber(line[1]), Integer.parseInt(line[2])});
						}
					} else {
//						fileData = HashVec.loadFileToHashVec(filename, indices[0], new int[]{indices[1],  indices[2]}, "\t", true, false);
//						for (String key : fileData.keySet()) {
//							markerMap.put(key, new int[]{Positions.chromosomeNumber(fileData.get(key).get(0)), Integer.parseInt(fileData.get(key).get(1))});
//						}
//						snps = fileData.keySet().toArray(new String[fileData.size()]);
//						fileData = null;
					}
					if (indices[4] == -1) {
						// no frequencies
						if (markerMap.isEmpty()) {
							fileData = HashVec.loadFileToHashVec(filename, indices[0], new int[]{indices[1],  indices[2]}, "\t", true, false);
							for (String key : fileData.keySet()) {
								markerMap.put(key, new int[]{Positions.chromosomeNumber(fileData.get(key).get(0).split("\t")[0]), Integer.parseInt(fileData.get(key).get(0).split("\t")[1])});
							}
							fileData = null;
						}
						freqs = get1000GFreq(markerMap);
					} else {
						freqs = null;//new HashMap<String, Double>();
//						while ((temp = reader.readLine()) != null) {
//							line = temp.trim().split("[\\s]+");
//							freqs.put(line[indices[0]], Double.valueOf(line[indices[4]]));
//						}
					}
					
					StringBuilder newHeaderSB = new StringBuilder("SNP\tChr\tPos\tFreq\tP");
					for (int i = 0; i < header.length; i++) {
						if (i != indices[0] && i != indices[1] && i != indices[2] && i != indices[3] && i != indices[4]) {
							newHeaderSB.append("\t").append(header[i]);
						}
					}
					
					PrintWriter metaWriter = Files.getAppropriateWriter(ext.rootOf(filename, false) + ".meta");
					metaWriter.println(newHeaderSB.toString());
					while ((temp = reader.readLine()) != null) {
						line = temp.trim().split("[\\s]+");
						String snp = line[indices[0]];
						String chr = indices[1] == -1 ? "" + markerMap.get(snp)[0] : line[indices[1]];
						String pos = indices[2] == -1 ? "" + markerMap.get(snp)[1] : line[indices[2]];
						String pval = line[indices[3]];
						String freq = indices[4] == -1 ? "" + (freqs == null || freqs.isEmpty() || freqs.get(snp) == null ? 0.0 : freqs.get(snp)) : line[indices[4]];
						StringBuilder writeLineSB = new StringBuilder();
						writeLineSB.append(snp).append("\t")
									.append(chr).append("\t")
									.append(pos).append("\t")
									.append(freq).append("\t")
									.append(pval);
						for (int i = 0; i < line.length; i++) {
							if (i != indices[0] && i != indices[1] && i != indices[2] && i != indices[3] && i != indices[4]) {
								writeLineSB.append("\t").append(line[i]);
							}
						}
						metaWriter.println(writeLineSB.toString());
					}
					metaWriter.flush();
					metaWriter.close();
				} else {
					System.err.println(errorMsg);
					reader.close();
					continue;
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public GeneScorePipeline(String metaDir, /*int numThreads, boolean plink, boolean regression, boolean histogram,*/ float[] indexThresholds, int[] windowMins, float[] windowExtThresholds) {
		this.metaDir = metaDir;
//		this.numThreads = numThreads;
//		this.runPlink = plink;
//		this.runRegression = runPlink && regression;
//		this.writeHist = runPlink && histogram;
		
		this.indexThresholds = indexThresholds;
		this.windowMinSizePerSides = windowMins;
		this.windowExtensionThresholds = windowExtThresholds;
		setFilePrefices();
		loadStudyFolders();
		createAffectedPhenoFiles();
		loadPhenoFiles();
		// instantiate inner hashmaps:
		for (Study study : studies) {
			for (String pref : analysisConstraints.keySet()) {
				HashMap<String, HashMap<String, RegressionResult>> res = new HashMap<String, HashMap<String, RegressionResult>>();
				for (String dFile : dataFiles) {
					String dataFile = ext.rootOf(dFile, false);
					HashMap<String, RegressionResult> res2 = new HashMap<String, RegressionResult>();
					res.put(dataFile, res2);
				}
				study.regressions.put(pref, res);
				
				HashMap<String, Integer> cntMap = new HashMap<String, Integer>();
				HashMap<String, Integer> hitMap = new HashMap<String, Integer>();
				HashMap<String, Integer> cntMap2 = new HashMap<String, Integer>();
				study.hitWindowCnts.put(pref, cntMap);
				study.hitSnpCounts.put(pref, hitMap);
				study.dataCounts.put(pref, cntMap2);
			}
		}
		loadDataCounts();
		runMetaHitWindows();
	}
	
	private void setFilePrefices() {
		for (float i : indexThresholds) {
			for (int m : windowMinSizePerSides) {
				for (float w : windowExtensionThresholds) {
					StringBuilder prefixSB = new StringBuilder();
					prefixSB.append(ext.formSciNot(i, 4, false))
						.append("_").append(ext.formSciNot(m, 4, false))
						.append("_").append(ext.formSciNot(w, 4, false));
					analysisConstraints.put(prefixSB.toString(), new Constraint(i, m, w));
				}
			}
		}
	}
	
	private void loadStudyFolders() {
		File dir = new File(this.metaDir);
		File[] fs = dir.listFiles();
		for (File f : fs) {
			if (f.isDirectory()) {
				Study study = new Study();
				study.studyName = ext.rootOf(f.getAbsolutePath(), true);
				study.studyDir = f.getAbsolutePath() + "\\";
				for (File f1 : f.listFiles()) {
					if (f1.getName().endsWith(".pheno")) {
						study.phenoFiles.add(f1.getName());
					}
				}
				studies.add(study);
			} else if (f.getAbsolutePath().endsWith(".meta")) {
				dataFiles.add(f.getName());
			}
		}
	}
	
	private void createAffectedPhenoFiles() {
		studyLoop : for (Study study : studies) {
			String famFile = study.studyDir + "plink.fam";
			String affFile = study.studyDir + "AFFECTED.pheno";
			
			// if affected.pheno file already exists, skip
			if (study.phenoFiles.contains("AFFECTED.pheno")) {
				continue;
			}
			// if any of the data folders have been created, skip creating affected.pheno file
			for (String dFile : dataFiles) {
				String dataFile = ext.rootOf(dFile, false);
				if ((new File(study.studyDir + dataFile + "\\")).exists()) {
					continue studyLoop;
				}
			}
			
			BufferedReader reader;
			// 0 - fid
			// 1 - iid
			// 2
			// 3
			// 4 - sex
			// 5 - pheno
			ArrayList<String> fam = new ArrayList<String>();
			ArrayList<String> pheno = new ArrayList<String>();
			String temp;
			try {
				reader = Files.getAppropriateReader(famFile);
				while ((temp = reader.readLine()) != null) {
					String[] line = temp.split("[\\s]+");
					String affLine = line[0] + "\t" + line[1] + "\t" + (ext.isMissingValue(line[4]) ? "." : -1 * (Integer.parseInt(line[4]) - 2)) + "\t" + line[5];
					if (!ext.isMissingValue(line[5])) {
						if (ext.isValidDouble(line[5]) && Double.parseDouble(line[5]) != 0.0) {
							fam.add(affLine);
						}
						pheno.add(line[5]);
					}
				}
				
				String[] unique = Array.unique(pheno.toArray(new String[]{}));
				if (unique.length == 1) {
					System.out.println("Error - no variance in pheno data from .fam file for study '" + study.studyName + "'");
					continue;
				}
				
				String header = "FID\tIID\tPHENO\tMALE";
				PrintWriter writer = Files.getAppropriateWriter(affFile);
				writer.println(header);
				for (String line : fam) {
					writer.println(line);
				}
				writer.flush();
				writer.close();
				
				study.phenoFiles.add("AFFECTED.pheno");
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
	}
	
	private void loadPhenoFiles() {
		for (Study study : studies) {
			for (String pheno : study.phenoFiles) {
				PhenoData pd = new PhenoData();
				pd.phenoName = pheno;
				
				try {
					BufferedReader reader = Files.getAppropriateReader(study.studyDir + pheno);
					String[] header = reader.readLine().split("\t");
					// fid == header[0]
					// iid == header[1]
					// depvar = header[2];
					ArrayList<String> covars = new ArrayList<String>();
					for (int i = 3; i < header.length; i++) {
						covars.add(header[i]);
					}
					pd.covars.addAll(covars);
					String temp = reader.readLine();
					indiv: do {
						String[] line = temp.split("\t");
						
						if (!ext.isMissingValue(line[2])) {
							PhenoIndiv pi = new PhenoIndiv();
							pi.fid = line[0];
							pi.iid = line[1];
							pi.depvar = Double.parseDouble(line[2]);
							for (int i = 3; i < line.length; i++) {
								if (ext.isMissingValue(line[i])) {
									continue indiv;
								}
								pi.covars.put(header[i], Double.parseDouble(line[i]));
							}
							pd.indivs.put(pi.fid + "\t" + pi.iid, pi);
						}
						
					} while ((temp = reader.readLine()) != null);
					
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
				
				study.phenoData.put(pheno, pd);
			}
		}
	}
	
	private void loadDataCounts() {
		String countsFile = metaDir + "data.cnt";
		
		if ((new File(countsFile).exists())) {
			try {
				BufferedReader reader = Files.getAppropriateReader(countsFile);
				String line = null;
				while ((line = reader.readLine()) != null && !"".equals(line)) {
					String[] temp = line.split("\t");
					HashMap<String, Integer> dFileCnts = dataCounts.get(temp[0]);
					if (dFileCnts == null) {
						dFileCnts = new HashMap<String, Integer>();
						dataCounts.put(temp[0], dFileCnts);
					}
					dFileCnts.put(temp[1], Integer.parseInt(temp[2]));
				}
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		

		for (String dFile : dataFiles) {
			HashMap<String, Constraint> threshNeed = new HashMap<String, GeneScorePipeline.Constraint>();
			
			HashMap<String, Integer> cnts = dataCounts.get(dFile);
			if (cnts == null) {
				threshNeed.putAll(analysisConstraints);
				cnts = new HashMap<String, Integer>();
				dataCounts.put(dFile, cnts);
			} else {
				for (String pref : analysisConstraints.keySet()) {
					if (!cnts.containsKey(pref)) {
						threshNeed.put(pref, analysisConstraints.get(pref));
					}
				}
			}
			
			if (threshNeed.size() > 0) {
				try {
					BufferedReader reader = Files.getAppropriateReader(metaDir + dFile);
					String line = reader.readLine();
					String[] dataHdrs = line.split("[\\s]+");
					int[] indices = ext.indexFactors(Aliases.PVALUES, dataHdrs, false, false);
					int ind = -1;
					for (int i : indices) {
						if (i > 0) {
							ind = i;
							break;
						}
					}
					if (ind > 0) {
						while ((line = reader.readLine()) != null) {
							String[] parts = line.split("\t");
							if (!ext.isMissingValue(parts[ind]) && ext.isValidDouble(parts[ind])) {
								double pval = Double.parseDouble(parts[ind]);
								for (java.util.Map.Entry<String, Constraint> constraint : threshNeed.entrySet()) {
									if (pval < constraint.getValue().indexThreshold) {
										Integer cnt = dataCounts.get(dFile).get(constraint.getKey());
										if (cnt == null) {
											cnt = 0;
										}
										cnt = cnt + 1;
										dataCounts.get(dFile).put(constraint.getKey(), cnt);
									}
								}
							}
						}
					}
				} catch (NumberFormatException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		StringBuilder output = new StringBuilder();
		for (java.util.Map.Entry<String, HashMap<String, Integer>> entry : dataCounts.entrySet()) {
			for (java.util.Map.Entry<String, Integer> subEntry : entry.getValue().entrySet()) {
				output.append(entry.getKey()).append("\t").append(subEntry.getKey()).append("\t").append(subEntry.getValue()).append("\n");
			}
		}
		Files.write(output.toString(), countsFile);
		
	}
	
	private void runMetaHitWindows() {
		String[][] factors = new String[][]{Aliases.MARKER_NAMES, Aliases.EFFECTS, Aliases.ALLELE_FREQS, Aliases.PVALUES};
		
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
			for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
				String[][] results = HitWindows.determine(metaDir + dFile, filePrefix.getValue().indexThreshold, filePrefix.getValue().windowMinSizePerSide, filePrefix.getValue().windowExtensionThreshold, DEFAULT_ADDL_ANNOT_VAR_NAMES);
				System.out.println(ext.getTime()+"]\tFound " + results.length + " hit windows");
				
				HashSet<String> hitMkrSet = new HashSet<String>();
				for (String[] mkrLine : results) {
					hitMkrSet.add(mkrLine[1]);
				}
				
				// read betas and freqs for hitwindow markers
				HashMap<String, double[]> dataMarkers = new HashMap<String, double[]>();
				HashMap<String, double[]> bimSubsetMarkers = new HashMap<String, double[]>();
				try {
					BufferedReader reader = Files.getAppropriateReader(metaDir + dFile);
					String line = reader.readLine();
					String[] temp = line.split("[\\s]+");
					int[] indices = ext.indexFactors(factors, temp, false, false, true, true, new Logger(), false);
					while ((line = reader.readLine()) != null) {
						String mkr = line.split("[\\s]+")[indices[0]];
						if (hitMkrSet.contains(mkr)) {
							if (ext.isMissingValue(line.split("[\\s]+")[indices[1]]) 
									|| ext.isMissingValue(line.split("[\\s]+")[indices[2]])
									|| ext.isMissingValue(line.split("[\\s]+")[indices[3]])) {
								hitMkrSet.remove(mkr);
								continue;
							}
							dataMarkers.put(mkr, new double[]{Double.parseDouble(line.split("[\\s]+")[indices[1]]), Double.parseDouble(line.split("[\\s]+")[indices[2]]), Double.parseDouble(line.split("[\\s]+")[indices[3]])});
						}
					}
					reader.close();
					
					double dataScore1 = getBetaFreqScore(dataMarkers);
					double dataScore2 = getChiDistRevScore(dataMarkers);
					
					// cross-ref PLINK markers
					for (Study study : studies) {
						HashSet<String> bimMkrSet = new HashSet<String>();
						BufferedReader bimReader = Files.getAppropriateReader(study.studyDir + study.plinkPref + ".bim");
						line = bimReader.readLine();
						do {
							bimMkrSet.add(line.split("[\\s]+")[bimMkrIndex]);
						} while ((line = bimReader.readLine()) != null);
						bimReader.close();
						
						// pull betas and freqs for union markers
						for (String mkr : bimMkrSet) {
							if (hitMkrSet.contains(mkr)) {
								bimSubsetMarkers.put(mkr, dataMarkers.get(mkr));
							}
						}

						// apply equation and set value for overall results
						double bimScore1 = getBetaFreqScore(bimSubsetMarkers);
						double bimScore2 = getChiDistRevScore(bimSubsetMarkers);
						
						HashMap<String, double[]> fileMap = study.scores.get(filePrefix.getKey());
						if (fileMap == null) {
							fileMap = new HashMap<String, double[]>();
							study.scores.put(filePrefix.getKey(), fileMap);
						}
						fileMap.put(dataFile, new double[]{bimScore1 / dataScore1, bimScore2 / dataScore2});
					}
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	private double getBetaFreqScore(HashMap<String, double[]> markerMap) {
		double sum = 0.0;
		// (Beta^2 * 2 * MAF (1-MAF))
		for (java.util.Map.Entry<String, double[]> entry : markerMap.entrySet()) {
			double score = 2 * (entry.getValue()[0] * entry.getValue()[0]) * entry.getValue()[1] * (1 - entry.getValue()[1]);  
			sum += score;
		}
		return sum;
	}
	
	private double getChiDistRevScore(HashMap<String, double[]> markerMap) {
		double sum = 0.0;
		for (java.util.Map.Entry<String, double[]> entry : markerMap.entrySet()) {
			sum += ProbDist.ChiDistReverse(entry.getValue()[2], 1);
		}
		return sum;
	}
	
	public void runPipeline() {
		System.out.println(ext.getTime()+"]\tProcessing study data [" + studies.size() + " total]:");
//		if (numThreads == 1) {
//			for (String studyDir : studyFolders) {
//				processStudy(studyDir);
//			}
		for (Study study : studies) {
			processStudy(study);
		}
		writeResults();
//		} else {
//			ExecutorService server = Executors.newFixedThreadPool(numThreads);
//			
//		}
		System.out.println(ext.getTime()+"]\tProcessing Complete!");
	}
	
	private void createFolders(Study study) {
		for (String dataFile : dataFiles) {
			String dataFolder = study.studyDir + ext.rootOf(dataFile, true) + "\\";
			for (String constraints : analysisConstraints.keySet()) {
				String constraintFolder = dataFolder + constraints + "\\";
				File f = new File(constraintFolder); 
				if (!(f.exists())) {
					f.mkdirs();
				}
			}
		}
	}
	
	private void processStudy(Study study) {
		try {
			createFolders(study);
			crossFilterMarkerData(study);
			runHitWindows(study);
			extractHitMarkerData(study);
			runPlink(study);
//			writePlink(study);
			runRegression(study);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void crossFilterMarkerData(Study study) throws IOException {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
			for (java.util.Map.Entry<String, Constraint> constraintEntry : analysisConstraints.entrySet()) {
				String crossFilterFile = study.studyDir + dataFile + "\\" + constraintEntry.getKey() + "\\bimData.xln";
				if ((new File(crossFilterFile).exists())) {
					System.out.println(ext.getTime()+"]\tCross-filtered data file already exists! [ --> '" + crossFilterFile + "']");
					study.hitSnpCounts.get(constraintEntry.getKey()).put(dataFile, Files.countLines(crossFilterFile, true));
					continue;
				}
				System.out.println(ext.getTime()+"]\tCross-filtering data and .BIM files [ --> '" + crossFilterFile + "']");
				BufferedReader bimReader;
				BufferedReader dataReader;
				PrintWriter dataWriter;
				HashMap<String, String[]> mkrsBim;
				
				mkrsBim = new HashMap<String, String[]>();
				bimReader = Files.getAppropriateReader(study.studyDir + study.plinkPref + ".bim");
				String line = bimReader.readLine();
				int cntAmbig = 0;
				do {
					String[] parts = line.split("[\\s]+");
					String mkr = parts[bimMkrIndex];
					String a1 = parts[bimA1Index];
					String a2 = parts[bimA2Index];
					if (Sequence.validAllele(a1) && Sequence.validAllele(a2) && !a1.equals(Sequence.flip(a2))) {
						mkrsBim.put(mkr, parts);
					} else {
						cntAmbig++;
					}
				} while ((line = bimReader.readLine()) != null);
				bimReader.close();
				System.out.println(ext.getTime()+"]\tFound " + cntAmbig + " ambiguous markers (will be excluded)");
				dataReader = Files.getAppropriateReader(metaDir + dFile);
				dataWriter = new PrintWriter(crossFilterFile);
				
				int cnt = 0;
				String dataHdr = dataReader.readLine();
				String[] dataHdrs = dataHdr.split("[\\s]+");
				int[] indices = ext.indexFactors(Aliases.PVALUES, dataHdrs, false, false);
				int ind = -1;
				for (int i : indices) {
					if (i > 0) {
						ind = i;
						break;
					}
				}
				dataWriter.println("MarkerName\tChr\tPosition\t" + Array.toStr(Array.subArray(dataHdrs, 1))); //Allele1\tAllele2\tFreq.Allele1.HapMapCEU\tb\tSE\tp\tN 
				while((line = dataReader.readLine()) != null) {
					String[] parts = line.split("[\\s]+");
					if (mkrsBim.containsKey(parts[0]) && (ind == -1 || (ext.isValidDouble(parts[ind]) && Double.parseDouble(parts[ind]) < constraintEntry.getValue().indexThreshold))) {
						String[] bimParts = mkrsBim.get(parts[0]);
						dataWriter.print(parts[0]);
						dataWriter.print("\t");
						dataWriter.print(bimParts[bimChrIndex]);
						dataWriter.print("\t");
						dataWriter.print(bimParts[bimPosIndex]);
						dataWriter.print("\t");
						dataWriter.println(Array.toStr(Array.subArray(parts, 1)));
						cnt++;
					}
				}
				dataWriter.flush();
				dataWriter.close();
				dataReader.close();
				
				study.hitSnpCounts.get(constraintEntry.getKey()).put(dataFile, cnt);
			}
		}
	}
	
	private void runHitWindows(Study study) {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
			
			for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
				File prefDir = new File(study.studyDir + dataFile + "\\" + filePrefix.getKey() + "\\");
				String crossFilterFile = prefDir + "\\bimData.xln";
				String hitsFile = prefDir + "\\hits_" + filePrefix.getKey() + ".out";
				if ((new File(hitsFile)).exists()) {
					System.out.println(ext.getTime()+"]\tHit window analysis file already exists! [ --> '" + hitsFile + "']");
					study.hitWindowCnts.get(filePrefix.getKey()).put(dataFile, Files.countLines(hitsFile, true));
					continue;
				}
				System.out.println(ext.getTime()+"]\tRunning hit window analysis [ --> '" + hitsFile + "']");
				String[][] results = HitWindows.determine(crossFilterFile, filePrefix.getValue().indexThreshold, filePrefix.getValue().windowMinSizePerSide, filePrefix.getValue().windowExtensionThreshold, DEFAULT_ADDL_ANNOT_VAR_NAMES);
				System.out.println(ext.getTime()+"]\tFound " + results.length + " hit windows");
				Files.writeMatrix(results, hitsFile, "\t");
				study.hitWindowCnts.get(filePrefix.getKey()).put(dataFile, results.length);
			}
		}
	}
	
	private void extractHitMarkerData(Study study) {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);

			for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
				File prefDir = new File(study.studyDir + dataFile + "\\" + filePrefix.getKey() + "\\");

				String crossFilterFile = prefDir + "\\bimData.xln";
				String hitsFile = prefDir + "\\hits_" + filePrefix.getKey() + ".out";
				String mkrDataFile = prefDir + "\\subsetData_" + filePrefix.getKey() + ".xln";
				if ((new File(mkrDataFile)).exists()) {
					System.out.println(ext.getTime()+"]\tHit window marker data file already exists! [ --> '" + mkrDataFile + "']");
					continue;
				}
				System.out.println(ext.getTime()+"]\tExtracting data for hit window markers [ --> '" + mkrDataFile + "']");
				String[] hitMarkers = HashVec.loadFileToStringArray(hitsFile, true, new int[]{hitsMkrIndex}, false);
				HashSet<String> hitMrkSet = new HashSet<String>();
				for (String mkr : hitMarkers) {
					hitMrkSet.add(mkr);
				}
				
				String[] header = Files.getHeaderOfFile(crossFilterFile, null);
				int[] cols = ext.indexFactors(LINKERS, header, false, true, false, null, false);
				String[][] bimData = HashVec.loadFileToStringMatrix(crossFilterFile, true, cols, false);
				
				PrintWriter hitDataWriter = Files.getAppropriateWriter(mkrDataFile);
				hitDataWriter.println("MarkerName\tAllele1\tb");
				for (String[] markerData : bimData) {
					if (hitMrkSet.contains(markerData[0])) {
						hitDataWriter.println(Array.toStr(markerData));
					}
				}
				hitDataWriter.flush();
				hitDataWriter.close();
			}
		}
	}
	
	private void runPlink(Study study) {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
	
			for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
				File prefDir = new File(study.studyDir + dataFile + "\\" + filePrefix.getKey() + "\\");
				if (!prefDir.exists()) {
					System.out.println(ext.getTime()+"]\tError - no subfolder for '" + filePrefix.getKey() + "' analysis");
					continue;
				}
				if ((new File(prefDir + "\\plink.profile")).exists()) {
					System.out.println(ext.getTime()+"]\tPlink analysis results file already exists! [ --> '" + prefDir + "\\plink.profile" + "']");
					continue;
				}
				String mkrDataFile = prefDir + "\\subsetData_" + filePrefix.getKey() + ".xln";
				System.out.print(ext.getTime()+"]\tRunning plink command [ --> '");
				String cmd = "plink" + /*(plink2 ? "2" : "") +*/ " --noweb --bfile ../../" + study.plinkPref + " --score " + mkrDataFile;
				System.out.println(cmd + "']");
				/*boolean results = */CmdLine.run(cmd, prefDir.getAbsolutePath());
			}
		}
	}
	
	private void writePlink(Study study) {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
	
			for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
				File prefDir = new File(study.studyDir + dataFile + "\\" + filePrefix.getKey() + "\\");
				if (!prefDir.exists()) {
					System.out.println(ext.getTime()+"]\tError - no subfolder for '" + filePrefix + "' analysis");
					continue;
				}
				if ((new File(prefDir + "\\runPlink.sh")).exists()) {
					System.out.println(ext.getTime()+"]\tPlink analysis shell script already exists! [ --> '" + prefDir + "\\runPlink.sh" + "']");
					continue;
				}
				String mkrDataFile = prefDir + "\\subsetData_" + filePrefix.getKey() + ".xln";
				System.out.println(ext.getTime()+"]\tWriting plink command");
				String cmd = "plink" + /*(plink2 ? "2" : "") +*/ " --noweb --bfile ../../" + study.plinkPref + " --score " + mkrDataFile;
				Files.write(cmd, prefDir.getAbsolutePath() + "\\runPlink.sh");
			}
		}
	}
	
	private void runRegression(Study study) {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
		
			for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
				try {
					File prefDir = new File(study.studyDir + dataFile + "\\" + filePrefix.getKey() + "\\");
					
					String plinkFile = prefDir + "\\plink.profile";
					
					HashMap<String, double[]> plinkData = new HashMap<String, double[]>();
					
					BufferedReader plinkReader = Files.getAppropriateReader(plinkFile);
					String line = plinkReader.readLine();
					while((line = plinkReader.readLine()) != null) {
						String[] parts = line.split("[\\s]+");
						String pheno = parts[3];
						String score = parts[6];
						plinkData.put(parts[1] + "\t" + parts[2], new double[]{Double.parseDouble(pheno), Double.parseDouble(score)});
					}
					
					if (!(new File(prefDir + "\\scores.hist").exists())) {
						double[] scores = new double[plinkData.size()];
						int ind = 0;
						for (double[] data : plinkData.values()) {
							scores[ind] = data[1];
							ind++;
						}
						String[] unq = Array.unique(Array.toStringArray(scores));
						if (unq.length == 1) {
							System.out.println(ext.getTime()+"]\tError - no variance in scores for " + dataFile + " / " + filePrefix.getKey() + " -- no .hist file created");
						} else {
							Files.write((new Histogram(scores)).getSummary().trim(), prefDir + "\\scores.hist");
						}
					}
				
					for (int i = 0; i < study.phenoFiles.size(); i++) {
						
						PhenoData pd = study.phenoData.get(study.phenoFiles.get(i));
						ArrayList<Double> depData = new ArrayList<Double>();
						ArrayList<double[]> baselineIndeps = new ArrayList<double[]>();
						ArrayList<double[]> indepData = new ArrayList<double[]>();
						for (java.util.Map.Entry<String, PhenoIndiv> indiv : pd.indivs.entrySet()) {
							if (plinkData.containsKey(indiv.getKey())) {
								depData.add(pd.indivs.get(indiv.getKey()).depvar);
								double[] baseData = new double[pd.covars.size()];
								double[] covarData = new double[pd.covars.size() + 1];
								covarData[0] = plinkData.get(indiv.getKey())[1];
								for (int k = 1; k < pd.covars.size()+1; k++) {
									baseData[k - 1] = pd.indivs.get(indiv.getKey()).covars.get(pd.covars.get(k - 1));
									covarData[k] = pd.indivs.get(indiv.getKey()).covars.get(pd.covars.get(k - 1));
								}
								baselineIndeps.add(baseData);
								indepData.add(covarData);
							}
						}
						
						double[][] baseCovars = new double[baselineIndeps.size()][];
						double[][] covars = new double[indepData.size()][];
						for (int k = 0; k < covars.length; k++) {
							covars[k] = indepData.get(k);
							baseCovars[k] = baselineIndeps.get(k);
						}
						RegressionModel baseModel = RegressionModel.determineAppropriate(Array.toDoubleArray(depData), baseCovars, false, true);
						RegressionModel model = RegressionModel.determineAppropriate(Array.toDoubleArray(depData), covars, false, true);
						
						RegressionResult rr = new RegressionResult();
						if (model.analysisFailed()) {
							rr.baseRSq = baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare();
							rr.beta = Double.NaN;
							rr.se = Double.NaN;
							rr.rsq = Double.NaN;
							rr.pval = Double.NaN;
							rr.num = depData.size();
							rr.logistic = model.isLogistic();
						} else {
							int ind = -1;
							for (int l = 0; l < model.getVarNames().length; l++) {
								if ("Indep 1".equals(model.getVarNames()[l])) {
									ind = l;
									break;
								}
							}
							if (ind == -1) {
								rr.baseRSq = baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare();
								rr.beta = Double.NaN;
								rr.se = Double.NaN;
								rr.rsq = model.getRsquare();
								rr.pval = Double.NaN;
								rr.num = depData.size();
								rr.logistic = model.isLogistic();
								rr.stats = Double.NaN;
							} else {
								rr.baseRSq = baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare();
								rr.beta = model.getBetas()[ind];
								rr.se = model.getSEofBs()[ind];
								rr.rsq = model.getRsquare();
								rr.pval = model.getSigs()[ind];
								rr.num = depData.size();
								rr.logistic = model.isLogistic();
								rr.stats = model.getStats()[ind];
							}
						}
						study.regressions.get(filePrefix.getKey()).get(dataFile).put(pd.phenoName, rr);
						
						baseModel = null;
						model = null;
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	private void writeResults() {
		PrintWriter writer;
		
		String resFile = metaDir + "results.xln";
		System.out.println(ext.getTime()+"]\tWriting regression results... [ --> " + resFile + "]");
		writer = Files.getAppropriateWriter(resFile);
		writer.println(REGRESSION_HEADER);
		
		for (Study study : studies) {
			for (String dFile : dataFiles) {
				String dataFile = ext.rootOf(dFile, false);
				for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
					HashMap<String, RegressionResult> phenoResults = study.regressions.get(filePrefix.getKey()).get(dataFile);
					
					String resultPrefix = study.studyName + "\t" + dataFile + "\t" + ext.formSciNot(filePrefix.getValue().indexThreshold, 5, false) + "\t";
					for (String pheno : study.phenoFiles) {
						RegressionResult rr = phenoResults.get(pheno);
						
						String pvalExcl = rr.num == 0 ? "." : (rr.logistic ? "=NORMSDIST(" + Math.sqrt(rr.stats) + ")" : "=TDIST(" + Math.abs(rr.stats) + "," + rr.num + ",2)");
						
						StringBuilder sb = new StringBuilder(resultPrefix)
												.append(pheno).append("\t")
												.append(rr.baseRSq).append("\t")
												.append(rr.rsq).append("\t")
												.append(rr.pval).append("\t")
												.append(rr.beta).append("\t")
												.append(rr.se).append("\t")
												.append(rr.num).append("\t")
												.append(dataCounts.get(dFile).get(filePrefix.getKey())).append("\t")
												.append(study.hitSnpCounts.get(filePrefix.getKey()).get(dataFile)).append("\t")
												.append(study.hitWindowCnts.get(filePrefix.getKey()).get(dataFile)).append("\t")
												.append(study.scores.get(filePrefix.getKey()).get(dataFile)[0]).append("\t")
												.append(study.scores.get(filePrefix.getKey()).get(dataFile)[1]).append("\t")
												.append(pvalExcl);
						
						writer.println(sb.toString());
					}
				}
			}
		}
		
		writer.flush();
		writer.close();
	}
	
	
	
//	private void writeToForestInput() {
//		String resultsFile = metaDir + "regressions.out";
//		PrintWriter writer = Files.getAppropriateWriter(resultsFile);
//		StringBuilder header = new StringBuilder("Name\tbeta\tse");
//		for (Study study : studies) {
//			header.append("\tbeta.").append(study.studyName).append("\tse.").append(study.studyName);
//		}
//		writer.println(header);
//		
//		StringBuilder dataSB;
//		for (String constraint : analysisConstraints.keySet()) {
//			dataSB = new StringBuilder("GeneScore_").append(constraint);
//			dataSB.append("\t");
//			// TODO meta beta
//			dataSB.append("\t");
//			// TODO meta se
//			for (Study study : studies) {
//				dataSB.append("\t").append(study.regressionResults.get(constraint)[0]).append("\t").append(study.regressionResults.get(constraint)[1]);
//			}
//			writer.println(dataSB.toString());
//		}
//		
//		writer.flush();
//		writer.close();
//	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		
		String broot = null;
		
		float[] iT = new float[]{DEFAULT_INDEX_THRESHOLD};
		int[] mZ = new int[]{DEFAULT_WINDOW_MIN_SIZE_PER_SIDE};
		float[] wT = new float[]{DEFAULT_WINDOW_EXTENSION_THRESHOLD}; 
		
		boolean test = false;
		if (test) {
			preprocessDataFiles(new String[]{
					"D:/height/GeneScorePipeline/HeightScoring/ExtremeHeight.xln",
					"D:/height/GeneScorePipeline/HeightScoring/height_full.xln",
					"D:/height/GeneScorePipeline/HeightScoring/TannerSexCombined.xln"
			});
			return;
		}
		
//		int threads = 1;
		
		String usage =  "\n" + 
				"lab.MultiGeneScorePipeline requires 2+ arguments\n" + 
				"   (1) Metastudy directory root, containing subdirectories for each study (i.e. broot=C:/ (not the default))\n" +
				"       OPTIONAL:\n" + 
				"   (2) p-value threshold (or comma-delimited list) for index SNPs (i.e. indexThresh=" + DEFAULT_INDEX_THRESHOLD + " (default))\n" + 
				"   (3) minimum num bp per side of window (or comma delimited list) (i.e. minWinSize=" + DEFAULT_WINDOW_MIN_SIZE_PER_SIDE + " (default))\n" + 
				"   (4) p-value threshold to extend the window (or comma delimited list) (i.e. winThresh=" + DEFAULT_WINDOW_EXTENSION_THRESHOLD + " (default))\n" +
//				"   (8) Number of threads to use for computation (i.e. threads=" + threads + " (default))\n" + 
				"";
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("broot=")) {
				broot = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("indexThresh=")) {
				String[] lst = args[i].split("=")[1].split(",");
				int cntValid = 0;
				for (String poss : lst) {
					if (ext.isValidDouble(poss)) cntValid++;
				}
				iT = new float[cntValid];
				int ind = 0;
				for (String poss : lst) {
					if (ext.isValidDouble(poss)) {
						iT[ind] = Float.parseFloat(poss);
						ind++;
					}
				}
				numArgs--;
			} else if (args[i].startsWith("minWinSize=")) {
				String[] lst = args[i].split("=")[1].split(",");
				int cntValid = 0;
				for (String poss : lst) {
					if (ext.isValidDouble(poss)) cntValid++;
				}
				mZ = new int[cntValid];
				int ind = 0;
				for (String poss : lst) {
					if (ext.isValidInteger(poss)) {
						mZ[ind] = Integer.parseInt(poss);
						ind++;
					}
				}
				numArgs--;
			} else if (args[i].startsWith("winThresh=")) {
				String[] lst = args[i].split("=")[1].split(",");
				int cntValid = 0;
				for (String poss : lst) {
					if (ext.isValidDouble(poss)) cntValid++;
				}
				wT = new float[cntValid];
				int ind = 0;
				for (String poss : lst) {
					if (ext.isValidDouble(poss)) {
						wT[ind] = Float.parseFloat(poss);
						ind++;
					}
				}
				numArgs--;
			}else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0 || args.length == 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		File dir = new File(broot);
		if (!dir.isDirectory()) {
			System.err.println("Error - argument 'broot' must be a valid directory");
			System.exit(1);
		}
//		if (regress && !runPlink) {
//			System.err.println("Error - '-runPlink' option is required for '-regress' option");
//			System.exit(1);
//		}
//		if (writeHist && !runPlink) {
//			System.err.println("Error - '-runPlink' option is required for '-writeHist' option");
//			System.exit(1);
//		}
		
		GeneScorePipeline gsp = new GeneScorePipeline(broot, /*threads, runPlink, regress, writeHist,*/ iT, mZ, wT);
		gsp.runPipeline();
	}
	
}

