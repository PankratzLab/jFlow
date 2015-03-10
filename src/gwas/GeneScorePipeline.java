package gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import stats.RegressionModel;
import bioinformatics.Sequence;	
import common.Aliases;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.ext;

public class GeneScorePipeline {
	
	private static float DEFAULT_INDEX_THRESHOLD = (float)0.00000005;
	private static int DEFAULT_WINDOW_MIN_SIZE_PER_SIDE = 500000;// 500kb each side is technically a 1M window until the next hit region, but we now take this into consideration in the main algorithm
	private static float DEFAULT_WINDOW_EXTENSION_THRESHOLD = (float)0.00000005; // (float)0.00001;
	private static String[] DEFAULT_ADDL_ANNOT_VAR_NAMES = new String[0];
	
	private static final String[][] LINKERS = {
		Aliases.MARKER_NAMES,
		Aliases.ALLELES[0],
		Aliases.EFFECTS 
	};
	
	private static final String REGRESSION_HEADER = "STUDY\tDATAFILE\tINDEX-THRESHOLD\tFACTOR\tR-SQR\tSIG\tBETA\tSE\tNUM";
	private String metaDir;
	
	private float[] indexThresholds = new float[]{DEFAULT_INDEX_THRESHOLD};
	private int[] windowMinSizePerSides = new int[]{DEFAULT_WINDOW_MIN_SIZE_PER_SIDE};
	private float[] windowExtensionThresholds = new float[]{DEFAULT_WINDOW_EXTENSION_THRESHOLD};
	
//	private int numThreads = 1;
	private boolean runPlink = false;
	private boolean runRegression = false;
	
	private ArrayList<String> dataFiles = new ArrayList<String>();
	private ArrayList<Study> studies = new ArrayList<GeneScorePipeline.Study>();
	private HashMap<String, Constraints> analysisConstraints = new HashMap<String, GeneScorePipeline.Constraints>();
	
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
		HashMap<String, double[]> regressionResults = new HashMap<String, double[]>();
	}
	
	private class Constraints {
		final float indexThreshold;
		final int windowMinSizePerSide;
		final float windowExtensionThreshold;
//		final String phenoFile;
		public Constraints(float i, int m, float w/*, String pheno*/) {
			this.indexThreshold = i;
			this.windowMinSizePerSide = m;
			this.windowExtensionThreshold = w;
//			this.phenoFile = pheno;
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
	
	public GeneScorePipeline(String metaDir, /*String dataFile,*/ /*int numThreads,*/ boolean plink, boolean regression, float[] indexThresholds, int[] windowMins, float[] windowExtThresholds) {
		this.metaDir = metaDir;
//		this.dataFile = dataFile;
//		this.numThreads = numThreads;
		this.runPlink = plink;
		this.runRegression = runPlink && regression;
		
		this.indexThresholds = indexThresholds;
		this.windowMinSizePerSides = windowMins;
		this.windowExtensionThresholds = windowExtThresholds;
		setFilePrefices();
		loadStudyFolders();
		loadPhenoFiles();
	}
	
	private void setFilePrefices() {
		for (float i : indexThresholds) {
			for (int m : windowMinSizePerSides) {
				for (float w : windowExtensionThresholds) {
					StringBuilder prefixSB = new StringBuilder();
					prefixSB.append(ext.formSciNot(i, 4, false))
						.append("_").append(ext.formSciNot(m, 4, false))
						.append("_").append(ext.formSciNot(w, 4, false));
					analysisConstraints.put(prefixSB.toString(), new Constraints(i, m, w));
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
	
	public void runPipeline() {
		System.out.println("Processing study data[" + studies.size() + " total]:");
//		if (numThreads == 1) {
//			for (String studyDir : studyFolders) {
//				processStudy(studyDir);
//			}
			for (Study study : studies) {
				processStudy(study);
			}
			if (runRegression) {
				writeRegressionResults();
			}
//		} else {
//			ExecutorService server = Executors.newFixedThreadPool(numThreads);
//			
//		}
		System.out.println("Processing Complete!");
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
			if (runPlink) {
				runPlink(study);
			} else {
				writePlink(study);
			}
			if (runRegression) {
				runRegression(study);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void crossFilterMarkerData(Study study) throws IOException {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
			String crossFilterFile = study.studyDir + dataFile + "\\bimData.xln";
			if ((new File(crossFilterFile).exists())) {
				System.out.println("Cross-filtered data file already exists! [ --> '" + crossFilterFile + "']");
				continue;
			}
			System.out.println("Cross-filtering data and .BIM files [ --> '" + crossFilterFile + "']");
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
			System.out.println("Found " + cntAmbig + " ambiguous markers (will be excluded)");
			dataReader = Files.getAppropriateReader(metaDir + dFile);
			dataWriter = new PrintWriter(crossFilterFile);
			
			String[] dataHdrs = dataReader.readLine().split("[\\s]+");
			dataWriter.println("MarkerName\tChr\tPosition\t" + Array.toStr(Array.subArray(dataHdrs, 1))); //Allele1\tAllele2\tFreq.Allele1.HapMapCEU\tb\tSE\tp\tN
			while((line = dataReader.readLine()) != null) {
				String[] parts = line.split("[\\s]+");
				if (mkrsBim.containsKey(parts[0])) {
					String[] bimParts = mkrsBim.get(parts[0]);
					dataWriter.print(parts[0]);
					dataWriter.print("\t");
					dataWriter.print(bimParts[bimChrIndex]);
					dataWriter.print("\t");
					dataWriter.print(bimParts[bimPosIndex]);
					dataWriter.print("\t");
					dataWriter.println(Array.toStr(Array.subArray(parts, 1)));
				}
			}
			dataWriter.flush();
			dataWriter.close();
			dataReader.close();
		}
	}
	
	private void runHitWindows(Study study) {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
			String crossFilterFile = study.studyDir + dataFile + "\\bimData.xln";
			
			for (java.util.Map.Entry<String, Constraints> filePrefix : analysisConstraints.entrySet()) {
				File prefDir = new File(study.studyDir + dataFile + "\\" + filePrefix.getKey() + "\\");
				String hitsFile = prefDir + "\\hits_" + filePrefix.getKey() + ".out";
				if ((new File(hitsFile)).exists()) {
					System.out.println("Hit window analysis file already exists! [ --> '" + hitsFile + "']");
					continue;
				}
				System.out.println("Running hit window analysis [ --> '" + hitsFile + "']");
				String[][] results = HitWindows.determine(crossFilterFile, filePrefix.getValue().indexThreshold, filePrefix.getValue().windowMinSizePerSide, filePrefix.getValue().windowExtensionThreshold, DEFAULT_ADDL_ANNOT_VAR_NAMES);
				System.out.println("Found " + results.length + " hit windows");
				Files.writeMatrix(results, hitsFile, "\t");
			}
		}
	}
	
	private void extractHitMarkerData(Study study) {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
			String crossFilterFile = study.studyDir + dataFile + "\\bimData.xln";

			for (java.util.Map.Entry<String, Constraints> filePrefix : analysisConstraints.entrySet()) {
				File prefDir = new File(study.studyDir + dataFile + "\\" + filePrefix.getKey() + "\\");

				String hitsFile = prefDir + "\\hits_" + filePrefix.getKey() + ".out";
				String mkrDataFile = prefDir + "\\subsetData_" + filePrefix.getKey() + ".xln";
				if ((new File(mkrDataFile)).exists()) {
					System.out.println("Hit window marker data file already exists! [ --> '" + mkrDataFile + "']");
					continue;
				}
				System.out.println("Extracting data for hit window markers [ --> '" + mkrDataFile + "']");
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
	
			for (java.util.Map.Entry<String, Constraints> filePrefix : analysisConstraints.entrySet()) {
				File prefDir = new File(study.studyDir + dataFile + "\\" + filePrefix.getKey() + "\\");
	//			File prefDir = new File(study.studyDir + filePrefix.getKey() + "\\");
				if (!prefDir.exists()) {
					System.out.println("Error - no subfolder for '" + filePrefix.getKey() + "' analysis");
					continue;
				}
				if ((new File(prefDir + "\\plink.profile")).exists()) {
					System.out.println("Plink analysis results file already exists! [ --> '" + prefDir + "\\plink.profile" + "']");
					continue;
				}
				String mkrDataFile = prefDir + "\\subsetData_" + filePrefix.getKey() + ".xln";
				System.out.print("Running plink command [ --> '");
				String cmd = "plink" + /*(plink2 ? "2" : "") +*/ " --noweb --bfile ../../" + study.plinkPref + " --score " + mkrDataFile;
				System.out.println(cmd + "']");
				/*boolean results = */CmdLine.run(cmd, prefDir.getAbsolutePath());
			}
		}
	}
	
	private void writePlink(Study study) {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
	
			for (java.util.Map.Entry<String, Constraints> filePrefix : analysisConstraints.entrySet()) {
				File prefDir = new File(study.studyDir + dataFile + "\\" + filePrefix.getKey() + "\\");
//				File prefDir = new File(study.studyDir + filePrefix.getKey() + "\\");
				if (!prefDir.exists()) {
					System.out.println("Error - no subfolder for '" + filePrefix + "' analysis");
					continue;
				}
				if ((new File(prefDir + "\\runPlink.sh")).exists()) {
					System.out.println("Plink analysis shell script already exists! [ --> '" + prefDir + "\\runPlink.sh" + "']");
					continue;
				}
				String mkrDataFile = prefDir + "\\subsetData_" + filePrefix.getKey() + ".xln";
				System.out.println("Writing plink command");
				String cmd = "plink" + /*(plink2 ? "2" : "") +*/ " --noweb --bfile ../../" + study.plinkPref + " --score " + mkrDataFile;
				Files.write(cmd, prefDir.getAbsolutePath() + "\\runPlink.sh");
			}
		}
	}
	
	private void runRegression(Study study) {
		for (String dFile : dataFiles) {
			String dataFile = ext.rootOf(dFile, false);
		
			for (java.util.Map.Entry<String, Constraints> filePrefix : analysisConstraints.entrySet()) {
				try {
					File prefDir = new File(study.studyDir + dataFile + "\\" + filePrefix.getKey() + "\\");
//					File prefDir = new File(study.studyDir + filePrefix.getKey() + "\\");
//					if (!prefDir.exists()) {
//						System.out.println("Error - no subfolder for '" + filePrefix.getKey() + "' analysis");
//						continue;
//					}
					
					String plinkFile = prefDir + "\\plink.profile";
					
					HashMap<String, double[]> plinkData = new HashMap<String, double[]>();
					
					int cntTotal = 0, cntZero = 0;
					BufferedReader plinkReader = Files.getAppropriateReader(plinkFile);
					String line = plinkReader.readLine();
					while((line = plinkReader.readLine()) != null) {
						String[] parts = line.split("[\\s]+");
						String pheno = parts[3];
						String score = parts[6];
						plinkData.put(parts[1] + "\t" + parts[2], new double[]{Double.parseDouble(pheno), Double.parseDouble(score)});
						cntTotal++;
						if ("0".equals(pheno) || "0.0".equals(pheno) || (ext.isValidInteger(pheno) && Integer.parseInt(pheno) == 0)) {
							cntZero++;
						}
					}
				
					for (int i = 0; i < study.phenoFiles.size() + 1; i++) {
						
						if (i > 0) {
							PhenoData pd = study.phenoData.get(study.phenoFiles.get(i-1));
							ArrayList<Double> depData = new ArrayList<Double>();
							ArrayList<double[]> indepData = new ArrayList<double[]>();
	
							for (java.util.Map.Entry<String, PhenoIndiv> indiv : pd.indivs.entrySet()) {
								if (plinkData.containsKey(indiv.getKey())) {
									depData.add(pd.indivs.get(indiv.getKey()).depvar);
									double[] covarData = new double[pd.covars.size() + 1];
									covarData[0] = plinkData.get(indiv.getKey())[1];
									for (int k = 1; k < pd.covars.size()+1; k++) {
										covarData[k] = pd.indivs.get(indiv.getKey()).covars.get(pd.covars.get(k - 1));
									}
									indepData.add(covarData);
								}
							}
							
							
							double[][] covars = new double[indepData.size()][];
							for (int k = 0; k < covars.length; k++) {
								covars[k] = indepData.get(k);
							}
							RegressionModel model = RegressionModel.determineAppropriate(Array.toDoubleArray(depData), covars, false, true);
							
							if (model.analysisFailed()) {
								study.regressionResults.put(filePrefix.getKey() + dataFile + pd.phenoName, new double[]{Double.NaN, Double.NaN, Double.NaN, Double.NaN, depData.size()});
							} else {
								double beta = model.getBetas()[1];
								double se = model.getSEofBs()[1];
								double rsq = model.getRsquare();
								double sig = model.getSigs()[1];
								study.regressionResults.put(filePrefix.getKey() + dataFile + pd.phenoName, new double[]{rsq, sig, beta, se, depData.size()});
							}
							
							model = null;
						} else if (cntZero < cntTotal) {
							ArrayList<Double> depData = new ArrayList<Double>();
							ArrayList<Double> indepData = new ArrayList<Double>();
							
							for (java.util.Map.Entry<String, double[]> entry : plinkData.entrySet()) {
								if (entry.getValue()[0] != 0) {
									depData.add(entry.getValue()[0]);
									indepData.add(entry.getValue()[1]);
								}
							}
							
							RegressionModel model = RegressionModel.determineAppropriate(Array.toDoubleArray(depData), Array.toDoubleArray(indepData), false, true);
							
							if (model.analysisFailed()) {
								study.regressionResults.put(filePrefix.getKey() + dataFile, new double[]{Double.NaN, Double.NaN, Double.NaN, Double.NaN, depData.size()});
							} else {
								double beta = model.getBetas()[1];
								double se = model.getSEofBs()[1];
								double rsq = model.getRsquare();
								double sig = model.getSigs()[1];
								study.regressionResults.put(filePrefix.getKey() + dataFile, new double[]{rsq, sig, beta, se, depData.size()});
							}
							
							model = null;
						} else {
							System.out.println("Error - no variance in dependent variable!");
						}
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	private void writeRegressionResults() {
		PrintWriter writer;
		
		String resFile = metaDir + "results.xln";
		System.out.println("Writing regression results... [ --> " + resFile + "]");
		writer = Files.getAppropriateWriter(resFile);
		writer.println(REGRESSION_HEADER);
		
		for (Study study : studies) {
			for (String dFile : dataFiles) {
				String dataFile = ext.rootOf(dFile, false);
				for (java.util.Map.Entry<String, Constraints> filePrefix : analysisConstraints.entrySet()) {
					String resultPrefix = study.studyName + "\t" + dataFile + "\t" + ext.formSciNot(filePrefix.getValue().indexThreshold, 5, false) + "\t";
					writer.println(resultPrefix + "SCORE\t" + Array.toStr(study.regressionResults.get(filePrefix.getKey() + dataFile)));
					for (String pheno : study.phenoFiles) {
						double[] results = study.regressionResults.get(filePrefix.getKey() + dataFile + pheno);
						writer.println(resultPrefix + pheno + "\t" + Array.toStr(results));
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
		String dfile = "data.xln";
		
		float[] iT = new float[]{DEFAULT_INDEX_THRESHOLD};
		int[] mZ = new int[]{DEFAULT_WINDOW_MIN_SIZE_PER_SIDE};
		float[] wT = new float[]{DEFAULT_WINDOW_EXTENSION_THRESHOLD}; 
		
//		int threads = 1;
		boolean runPlink = false;
		boolean regress = false;
		
		
		String usage =  "\n" + 
				"lab.MultiGeneScorePipeline requires 2+ arguments\n" + 
				"   (1) Metastudy directory root, containing subdirectories for each study (i.e. broot=C:/ (not the default))\n" +
				"       OPTIONAL:\n" + 
				"   (2) p-value threshold (or comma-delimited list) for index SNPs (i.e. indexThresh=" + DEFAULT_INDEX_THRESHOLD + " (default))\n" + 
				"   (3) minimum num bp per side of window (or comma delimited list) (i.e. minWinSize=" + DEFAULT_WINDOW_MIN_SIZE_PER_SIDE + " (default))\n" + 
				"   (4) p-value threshold to extend the window (or comma delimited list) (i.e. winThresh=" + DEFAULT_WINDOW_EXTENSION_THRESHOLD + " (default))\n" +
				"   (5) run Plink's SNP scoring routine after completion (i.e. -runPlink (not the default))\n" +
				"   (6) run a regression analyis after completion (requires '-runPlink') (i.e. -regress (not the default))\n" +
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
			} else if (args[i].startsWith("-runPlink")) {
				runPlink = true;
				numArgs--;
			} else if (args[i].startsWith("-regress")) {
				regress = true;
				numArgs--;
			} else {
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
		if (regress && !runPlink) {
			System.err.println("Error - '-runPlink' option is required for '-regress' option");
			System.exit(1);
		}

		GeneScorePipeline gsp = new GeneScorePipeline(broot, /*dfile,*/ /*threads,*/ runPlink, regress, iT, mZ, wT);
		gsp.runPipeline();
	}
	
}
