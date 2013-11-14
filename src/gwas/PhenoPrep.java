package gwas;

import java.io.*;
import java.util.*;

import stats.LeastSquares;

import mining.Transformations;
import common.*;

public class PhenoPrep {
	private String[] finalHeader;
	private String[] finalIDs;
	private double[][] database;
	private Logger log;

	public static void parse(String dir, String filename, String pheno, String transform, double sdThreshold, boolean winsorize, boolean remove, boolean makeResids, boolean afterResids, boolean inverseNormalize, String[] covars, String idFile, String outFile, Logger log) {
		PhenoPrep prep;
		
		if (winsorize && remove) {
			log.reportError("Error - you have selected to both Winsorize and remove outliers for phenotype '"+pheno+"'; pick one or the other");
			return;
		}

		if (afterResids && !winsorize && !remove) {
			log.reportError("Error - you have selected the \"after residuals\" option for phenotype '"+pheno+"' but have not selected to Winsorize or remove outliers; aborting");
			return;
		}

		if (!makeResids && afterResids) {
			log.reportError("Error - you have selected to Winsorize or remove outliers with the \"after residuals\" option for phenotype '"+pheno+"' but have not selected the \"make residuals\" option; aborting");
			return;
		}
		
		if (makeResids && covars == null) {
			log.reportError("Error - you have selected to make residuals without specifying any covariates to regress out for phenotype '"+pheno+"'; aborting");
			return;
		}
		
		if (covars == null) {
			covars = new String[0];
		}
		
		
		prep = new PhenoPrep(dir+filename, idFile==null? null : dir+idFile, pheno, covars, log);
		
		if (transform != null) {
			prep.transform(transform);
		}

		if (!afterResids) {
			prep.dealWithOutliers(winsorize, remove, sdThreshold);
		}
		
		if (makeResids) {
			prep.computeResiduals();
		}
		
		if (afterResids) {
			prep.dealWithOutliers(winsorize, remove, sdThreshold);
		}
		
		if (inverseNormalize) {
			prep.inverseNormalize();
		}
		
		prep.writeFinalFile(dir+outFile);
	}
	
	public PhenoPrep(String filename, String idFile, String pheno, String[] covars, Logger log) {
		BufferedReader reader;
		String[] line;
		String temp;
		int[] indices;
		String[] idsWithDNA;
		String id;
		String delimiter;
		Vector<String> vIDs;
		Vector<double[]> vData;
		double[] data;
		boolean use;
		String[] header;

		this.log = log;
		
		temp = Files.getFirstNLinesOfFile(filename, 1, log)[0];
		delimiter = ext.determineDelimiter(temp);
		
		indices = ext.indexFactors(
				Array.insertStringAt(pheno, covars, 0),
				temp.split(delimiter, -1),
				false, log, true, true);
		
		if (idFile == null) {
			idsWithDNA = null;
		} else {
			idsWithDNA = HashVec.loadFileToStringArray(idFile, false, new int[] {0}, false);
		}

		header = null;
		vIDs = new Vector<String>();
		vData = new Vector<double[]>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			temp = reader.readLine();
			header = temp.split(delimiter, -1);
			finalHeader = Array.subArray(header, indices);
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.split(delimiter, -1);
				id = line[0];
				if (idsWithDNA == null || ext.indexOfStr(id, idsWithDNA) >= 0) {
					line = Array.subArray(line, indices);
					use = true;
					for (int i = 0; i < line.length; i++) {
						if (ext.isMissingValue(line[i])) {
							use = false;
						}
					}
					if (use) {
						vIDs.add(id);
						data = Array.toDoubleArray(line);
						if (data == null) {
							System.err.println("Error - failed to parse data for "+id);
						}
						vData.add(data);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
		
		finalIDs = Array.toStringArray(vIDs);
		database = Matrix.toDoubleArrays(vData);
	}
	
	public void transform(String transform) {
		double[] data;

		data = Matrix.extractColumn(database, 0);

		if (transform == null) {
			// do nothing
		} else if (transform.equalsIgnoreCase("ln")) {
			data = Transformations.naturalLogTransform(data);
		} else if (transform.equalsIgnoreCase("sqrt")) {
			data = Transformations.sqrtTransform(data);
		} else {
			System.err.println("Error - unkown transform: '"+transform+"'");
			return;
		}
		
		for (int i = 0; i < data.length; i++) {
			database[i][0] = data[i];
		}
		
	}

	public void dealWithOutliers(boolean winsorize, boolean remove, double sdThreshold) {
		double[] data;
		double mean, sd, lowerThreshold, upperThreshold;
		boolean[] rowsToUse;

		data = Matrix.extractColumn(database, 0);

		mean = Array.mean(data);
		sd = Array.stdev(data);
		lowerThreshold = mean-sdThreshold*sd;
		upperThreshold = mean+sdThreshold*sd;
		
		if (winsorize) {
			for (int i = 0; i < data.length; i++) {
				if (data[i] < lowerThreshold) {
					database[i][0] = lowerThreshold;
				} else if (data[i] > upperThreshold) {
					database[i][0] = upperThreshold;
				} else {
					database[i][0] = data[i];
				}
			}
			
		}

		if (remove) {
			rowsToUse = new boolean[data.length];

			for (int i = 0; i < data.length; i++) {
				rowsToUse[i] = true;
				if (data[i] < lowerThreshold || data[i] > upperThreshold) {
					rowsToUse[i] = false;
				}
			}
			
			finalIDs = Array.subArray(finalIDs, rowsToUse);
			database = Matrix.subset(database, rowsToUse);
		}
	}

	public void inverseNormalize() {
		Matrix.overwriteColumn(database, 0, Array.inverseNormalize(Matrix.extractColumn(database, 0)), log);
	}

	private void computeResiduals() {
		LeastSquares reg;
		double[] deps, resids;
		double[][] indeps;
		
		deps = Matrix.extractColumn(database, 0);
		indeps = Matrix.extractColumns(database, Array.subArray(Array.intArray(database[0].length), 1));
		
		reg = new LeastSquares(deps, indeps, null, false, true);
		
		if (reg.analysisFailed()) {
			log.reportError("Error performing the regression model; check for collinearity");
			System.exit(1);
		}
		
		resids = reg.getResiduals();
		
		if (deps.length != resids.length) {
			log.reportError("Error - lost a few rows in regression model; aborting");
			System.exit(1);
		}
		
		database = Matrix.toMatrix(resids);
		finalHeader = new String[] {finalHeader[0]};		
	}

	public void writeFinalFile(String filename) {
		PrintWriter writer;
		String delimiter;
		
		try {
			writer = new PrintWriter(new FileWriter(filename));
			delimiter = Files.suggestDelimiter(filename, log);
			writer.println("id"+delimiter+Array.toStr(finalHeader, delimiter));
			for (int i = 0; i < finalIDs.length; i++) {
				writer.println(finalIDs[i]+delimiter+Array.toStr(database[i], -1, -1, delimiter));
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + filename);
			e.printStackTrace();
		}
		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String dir;
		String filename;
		String logfile = null;
		Logger log;
		double sdThreshold = 3.0;
		String pheno;
		String transform = null;
		String[] covars = null;
		String idFile = null;
		String outFile;
		boolean winsorize = false;
		boolean remove = false;
		boolean makeResids = false;
		boolean afterResids = false;
		boolean inverseNormalize = false;
		
//		filename = "ARIC_Whites_WBC.csv";
//		pheno = "WBC";
//		outFile = "ARIC_EA_WBC.csv";
//		idFile = "keeps.txt";
//		transform = "ln";
//		covars = new String[] {"Age", "Male", "CenterF", "CenterJ", "CenterM"};
//		idFile = "EA_keeps.dat";
//		winsorize = true;

		dir = "D:/SkatMeta/results_hemostasis/";
		filename = "pheno_F7.csv";
		pheno = "F7";
		outFile = "pheno_F7_winsorize_3sd.csv";
//		covars = new String[] {"V1AGE01", "Sex"};
//		winsorize = false;
//		remove = false;
//		makeResids = false;
//		afterResids = false;

		String usage = "\n" +
				"gwas.PhenoPrep requires 0-1 arguments\n" +
				"   (0) name of directory (i.e. dir=" + dir + " (default))\n" + 
				"   (1) name of input file (i.e. file=" + filename + " (default))\n" + 
				"   (2) phenotype (i.e. pheno=" + pheno + " (default))\n" + 
				"   (3) covariates (i.e. covar=Age,Sex,Site1,Site2 (not the default))\n" + 
				"   (4) name of output file (i.e. out=" + outFile + " (default))\n" + 
				"   (5) transformation to apply to phenotype (i.e. transform=" + transform + " (default; current options are ln, sqrt, or null for none))\n" + 
				"   (6) name of file with IDs in genotype file (i.e. ids=" + idFile + " (default; set to null to include all rows with complete data))\n" + 
				"   (7) winsorize phenotype (i.e. winsorize=" + winsorize + " (default))\n" + 
				"   (8) remove outliers (i.e. remove=" + remove + " (default))\n" + 
				"   (9) threshold in standard deviation units at which to winsorize or remove outliers (i.e. sdThrehsold=" + sdThreshold + " (default))\n" + 
				"  (10) generate residuals instead of including covariates (i.e. makeResids=" + makeResids + " (default))\n" + 
				"  (11) winsorize/remove outliers after generating residuals (i.e. afterResids=" + afterResids + " (default))\n" + 
				"  (12) inverse quantile normalize the final phenotype (i.e. inverseNormalize=" + inverseNormalize + " (default))\n" + 
				"  (13) (optional) name of log file to write to (i.e. log=[pheno].log (default))\n" + 
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("file=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("pheno=")) {
				pheno = ext.parseStringArg(args[i], "missingPhenotype");
				numArgs--;
			} else if (args[i].startsWith("covar=")) {
				covars = null;
				if (ext.parseStringArg(args[i], null) != null) {
					covars = ext.parseStringArg(args[i], "").split(",");
				}
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outFile = ext.parseStringArg(args[i], "output_file.csv");
				numArgs--;
			} else if (args[i].startsWith("transform=")) {
				transform = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("ids=")) {
				idFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("winsorize=")) {
				winsorize = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("remove=")) {
				remove = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("sdThreshold=")) {
				sdThreshold = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("makeResids=")) {
				makeResids = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("afterResids=")) {
				afterResids = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("inverseNormalize=")) {
				inverseNormalize = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		
		if (logfile == null) {
			logfile = dir+ext.replaceWithLinuxSafeCharacters(pheno, true)+".log";
		}
		log = new Logger(logfile);
		
		if (numArgs != 0) {
			log.reportError(usage);
			System.exit(1);
		}
		
		try {
			parse(dir, filename, pheno, transform, sdThreshold, winsorize, remove, makeResids, afterResids, inverseNormalize, covars, idFile, outFile, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
