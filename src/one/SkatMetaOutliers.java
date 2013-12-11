package one;

import gwas.PhenoPrep;

import java.io.*;
import java.util.*;
import common.*;
import stats.ProbDist;

public class SkatMetaOutliers {

	public static String calculateLambda(String dir, final String filenameExt, Logger log) {
		File[] fileNames;
		Vector<Double> pvals;
		Scanner reader;
		String[] line;
		int vecLength;
		Double median;

		if (log == null) {
			log = new Logger();
		}

		fileNames = new File(dir).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				boolean result = false;
				filename = ext.removeDirectoryInfo(filename);
				if (filename.endsWith(filenameExt)) {
					result = true;
				}
				return result;
			}
		});
		if (fileNames == null) {
			return null;
		}

		pvals = new Vector<Double>();
		for (int i = 0; i < fileNames.length; i++) {
			try {
				reader = new Scanner(new FileInputStream(fileNames[i]));
				reader.nextLine();
				while(reader.hasNext()) {
					line = reader.nextLine().split("[,\t]");
//					if (line[1].equals("NA")) {
//						System.out.println("");
//					}
					if (! line[1].equals("NA") && ! line[1].equals("N/A") && ! line[1].equals("Inf") && ! line[1].equals("") && ! (line[1]==null)) {
						pvals.add(Double.parseDouble(line[1]));
					}
				}
				reader.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		
		Collections.sort(pvals);
		vecLength = pvals.size();
		if ((vecLength % 2) == 0) {
			vecLength = vecLength/2;
			median = (pvals.elementAt(vecLength) + pvals.elementAt(vecLength - 1)) / 2;
		} else {
			median = pvals.elementAt((vecLength - 1)/2);
		}
		return ext.formDeci(ProbDist.ChiDistReverse(median, 1)/ProbDist.ChiDistReverse(0.50, 1), 4);
	}

	private static void summarizeAll(String dir, String idColName, String phenosCommaDelimited, String covarsCommaDelimited) {
		PrintWriter writer;
		String[] phenos, covars, transforms;
		Logger log;
		boolean winsorize, remove, makeResids, afterResids, inverseNormalize;
		String outFile, idFile;
		String[] rawData;
		double[] data;
		double mean, stdev, skewness, kurtosis;
		
		try {
			writer = new PrintWriter(new FileWriter(dir+"phenoSummary.xln"));
			writer.println("Trait\ttransform\twinsorize\tremoveOutliers\tmakeResiduals\tafterMakingResidualsDealWithOutliers\tN\tmean\tstdev\tskewness\tkurtosis\t'=SUM(ABS(SKEW)+ABS(KURT))");
		
			phenos = phenosCommaDelimited.split(",");
			covars = covarsCommaDelimited.split(",");
			
			transforms = new String[] {null, "ln", "sqrt"};
			
			log = new Logger(dir+"summarizeAll.log");
			
//			idFile = "EA_keeps.dat";
			idFile = null;
			inverseNormalize = false;
			for (int i = 0; i < phenos.length; i++) {
				for (int j = 0; j < transforms.length; j++) {
					for (int outlierMethods = 0; outlierMethods < 3; outlierMethods++) {
						if (outlierMethods == 0) {
							winsorize = false;
							remove = false;
						} else if (outlierMethods == 1) {
							winsorize = true;
							remove = false;
						} else {
							winsorize = false;
							remove = true;
						}
						for (int resids = 0; resids < (outlierMethods==0?1:3); resids++) {
							if (resids == 0) {
								makeResids = false;
								afterResids = false;
							} else if (resids == 1) {
								makeResids = true;
								afterResids = false;
							} else {
								makeResids = true;
								afterResids = true;
							}
							outFile = phenos[i];
							if (transforms[j] != null) {
								outFile += "_"+transforms[j];
							}
							if (winsorize) {
								outFile += "_win";
							}
							if (remove) {
								outFile += "_del";
							}
							if (makeResids) {
								if (afterResids) {
									outFile += "_afterResid";
								} else {
									outFile += "_beforeResid";
								}
							}
							System.out.println(outFile);
							outFile += ".csv";
							if (!Files.exists(dir+outFile)) {
								PhenoPrep.parse(dir, phenos[i]+".csv", idColName, phenos[i], transforms[j], 3.0, winsorize, remove, makeResids, afterResids, inverseNormalize, covars, idFile, null, outFile, log);
							}
							rawData = HashVec.loadFileToStringArray(dir+outFile, false, true, new int[] {1}, false, false, Files.determineDelimiter(dir+outFile, log));
							rawData = Array.removeFromArray(rawData, ext.MISSING_VALUES);
							data = Array.toDoubleArray(rawData);
							mean = Array.mean(data);
							stdev = Array.stdev(data);
							skewness = Array.skewness(data);
							kurtosis = Array.kurtosis(data);
							writer.println(phenos[i]+"\t"+transforms[j]+"\t"+winsorize+"\t"+remove+"\t"+makeResids+"\t"+afterResids+"\t"+data.length+"\t"+mean+"\t"+stdev+"\t"+skewness+"\t"+kurtosis+"\t"+(Math.abs(skewness)+Math.abs(kurtosis)));
						}
					}
				}
			}
			
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir+"phenoSummary.xln");
			e.printStackTrace();
		}
		
		
		
		
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int numArgs = args.length;
		Logger log;
		String dir;
//		String phenoFilename = "D:/SkatMeta/results_hemostasis/pheno_F7_del_3sd.csv";
		String phenoFilename = "D:/SkatMeta/test/pheno_F8_win_ln_chr1_beforeSave_results.csv";
		String[] processing = new String[] {"", "del", "winsorized"};
		String[] transform = new String[] {"", "ln", "sqrt"};
		String[] regression = new String[] {"", "AGE,SEX"};
		double[][] phenoTmp;
		double[] pheno;
		String idColName, phenos, covars;

		log = new Logger();
		
		dir = "D:/LITE/CHARGE-S/aric_wex_freeze3/testOutliers/";
		phenos = "Fibrinogen,F7,F8,vWF";
		covars = "V1AGE01,Sex,CenterF,CenterM";
		
		dir = "D:/ExomeChip/ARIC_primary/CompareTransformations/";
		phenos = "Hct,Hb,MCHC,MCV,RBC,MCH,RDW,WBC_TOTAL,WBC_NEUTRO,WBC_MONO,WBC_LYMPH,WBC_EOS,WBC_BASO";
		covars = "Age,Male";

		dir = "D:/LITE/CHARGE-S/aric_wex_freeze3/testOutliers/AAs/";
		idColName = "gwas_id";
		phenos = "Fibrinogen,F7,F8,vWF";
		covars = "V1AGE01,Sex";
		
		summarizeAll(dir, idColName, phenos, covars);
		System.exit(1);

		String usage = "\n" +
				"gwas.PhenoPrep requires 0-1 arguments\n" +
				"   (1) name of pheno file (i.e. pheno=" + phenoFilename + " (default))\n" + 
				"   (2) list of processing to apply to phenotype (i.e. processing=" + processing + " (default))\n" + 
				"   (3) list of transformation to apply to phenotype (i.e. transform=" + transform + " (default))\n" + 
				"   (4) list of regression to apply to phenotype (i.e. regression=" + regression + " (default))\n" + 
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("pheno=")) {
				phenoFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("processing=")) {
				processing = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("transform=")) {
				transform = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("regression=")) {
				regression = args[i].split("=")[1].split(",");
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		
		if (numArgs != 0) {
			log.reportError(usage);
			System.exit(1);
		}

//		double[] test = new double[] {0.1, 0.2, 0.3, 0.4, 0.5, .6, .7, .8, 1.9, 11};
//		log.report("Testing, testing ... \tskewness: " + Array.skewness(test) + "\tkurtosis: " + Array.kurtosis(test));


		dir = ext.parseDirectoryOfFile(phenoFilename) + ext.rootOf(phenoFilename) + "/results";
		String[] header = Files.getHeaderOfFile(phenoFilename, log);
		// just pheno: new int[] {1}
//		phenoTmp = Array.removeNaN(Matrix.toDoubleArrays(HashVec.loadFileToStringMatrix(phenoFilename, true, Array.subArray(Array.intArray(header.length), 1), false)));
		phenoTmp = Matrix.toDoubleArrays(HashVec.loadFileToStringMatrix(phenoFilename, true, Array.subArray(Array.intArray(header.length), 1), false));
		pheno = new double[phenoTmp.length];
		for (int i = 0; i < pheno.length; i++) {
			pheno[i] = phenoTmp[i][0];
		}
		log.report("lambda: " + calculateLambda(dir, ".csv", log) + "\tskewness: " + Array.skewness(pheno) + "\tkurtosis: " + Array.kurtosis(pheno));

//		for (int i = 0; i < processing.length; i ++) {
//			for (int j = 0; j < transform.length; j++) {
//				for (int k = 0; k < regression.length; k++) {
//					PhenoPrep.parse(dir, phenoFilename, phenoFilename, transform[j], null, processing[i], remove, makeResids, afterResids, covars, idFile, ext.rootOf(phenoFilename)+, null);
//					log.report(calculateLamda("D:/skatMeta/test", ".csv", log));
//				}
//			}
//		}

	}
}
