package one;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.util.Collections;
import java.util.Scanner;
import java.util.Vector;

import stats.ProbDist;
import common.Array;
import common.Logger;
import common.ext;

import gwas.PhenoPrep;

public class SkatMetaOutliers {

	@SuppressWarnings("resource")
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

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int numArgs = args.length;
		Logger log;
		String dir;
		String phenoFilename = "";
		String[] processing = new String[] {"", "del", "winsorized"};
		String[] transform = new String[] {"", "ln", "sqrt"};
		String[] regression = new String[] {"", "AGE,SEX"};

		log = new Logger();

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

		dir = ext.parseDirectoryOfFile(phenoFilename);
		log.report(calculateLambda(dir, ".csv", log));
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
