package org.genvisis.one;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.util.Collections;
import java.util.Scanner;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.stats.ProbDist;

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
			@Override
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
		for (File fileName : fileNames) {
			try {
				reader = new Scanner(new FileInputStream(fileName));
				reader.nextLine();
				while (reader.hasNext()) {
					line = reader.nextLine().split("[,\t]");
					// if (line[1].equals("NA")) {
					// System.out.println("");
					// }
					if (!line[1].equals("NA")	&& !line[1].equals("N/A") && !line[1].equals("Inf")
							&& !line[1].equals("") && !(line[1] == null)) {
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
			vecLength = vecLength / 2;
			median = (pvals.elementAt(vecLength) + pvals.elementAt(vecLength - 1)) / 2;
		} else {
			median = pvals.elementAt((vecLength - 1) / 2);
		}
		return ext.formDeci(ProbDist.ChiDistReverse(median, 1) / ProbDist.ChiDistReverse(0.50, 1), 4);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		Logger log;
		String dir;
		// String phenoFilename = "D:/SkatMeta/results_hemostasis/pheno_F7_del_3sd.csv";
		String phenoFilename = "D:/SkatMeta/test/pheno_F8_win_ln_chr1_beforeSave_results.csv";
		String[] processing = new String[] {"", "del", "winsorized"};
		String[] transform = new String[] {"", "ln", "sqrt"};
		String[] regression = new String[] {"", "AGE,SEX"};
		double[][] phenoTmp;
		double[] pheno;

		log = new Logger();

		String usage = "\n"	+ "gwas.PhenoPrep requires 0-1 arguments\n"
										+ "   (1) name of pheno file (i.e. pheno=" + phenoFilename + " (default))\n"
										+ "   (2) list of processing to apply to phenotype (i.e. processing="
										+ ArrayUtils.toStr(processing, ",") + " (default))\n"
										+ "   (3) list of transformation to apply to phenotype (i.e. transform="
										+ ArrayUtils.toStr(transform, ",") + " (default))\n"
										+ "   (4) list of regression to apply to phenotype (i.e. regression="
										+ ArrayUtils.toStr(regression, ",") + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("pheno=")) {
				phenoFilename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("processing=")) {
				processing = arg.split("=")[1].split(",");
				numArgs--;
			} else if (arg.startsWith("transform=")) {
				transform = arg.split("=")[1].split(",");
				numArgs--;
			} else if (arg.startsWith("regression=")) {
				regression = arg.split("=")[1].split(",");
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}

		if (numArgs != 0) {
			log.reportError(usage);
			System.exit(1);
		}

		// double[] test = new double[] {0.1, 0.2, 0.3, 0.4, 0.5, .6, .7, .8, 1.9, 11};
		// log.report("Testing, testing ... \tskewness: " + Array.skewness(test) + "\tkurtosis: " +
		// Array.kurtosis(test));


		dir = ext.parseDirectoryOfFile(phenoFilename) + ext.rootOf(phenoFilename) + "/results";
		String[] header = Files.getHeaderOfFile(phenoFilename, log);
		// just pheno: new int[] {1}
		// phenoTmp =
		// Array.removeNaN(Matrix.toDoubleArrays(HashVec.loadFileToStringMatrix(phenoFilename, true,
		// Array.subArray(Array.intArray(header.length), 1), false)));
		phenoTmp = Matrix.toDoubleArrays(HashVec.loadFileToStringMatrix(phenoFilename, true,
																																		ArrayUtils.subArray(	ArrayUtils.arrayOfIndices(header.length),
																																										1),
																																		false));
		pheno = new double[phenoTmp.length];
		for (int i = 0; i < pheno.length; i++) {
			pheno[i] = phenoTmp[i][0];
		}
		log.report("lambda: "	+ calculateLambda(dir, ".csv", log) + "\tskewness: "
								+ ArrayUtils.skewness(pheno) + "\tkurtosis: " + ArrayUtils.kurtosis(pheno));

		// for (int i = 0; i < processing.length; i ++) {
		// for (int j = 0; j < transform.length; j++) {
		// for (int k = 0; k < regression.length; k++) {
		// PhenoPrep.parse(dir, phenoFilename, phenoFilename, transform[j], null, processing[i], remove,
		// makeResids, afterResids, covars, idFile, ext.rootOf(phenoFilename)+, null);
		// log.report(calculateLamda("D:/skatMeta/test", ".csv", log));
		// }
		// }
		// }

	}
}
