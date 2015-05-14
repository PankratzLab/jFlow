package gwas;

import java.io.*;
import java.util.*;

import stats.LeastSquares;
import mining.Transformations;
import common.*;

public class PhenoPrep {
	public static final String[] SUMMARY_INFO_HEADER = {"Race", "Trait", "meanTrait", "medianTrait", "stdevTrait", "minTrait", "maxTrait", "numFemales", "numMales", "meanAge", "medianAge", "stdevAge", "minAge", "maxAge", "numBelowLowerThrehsold", "numAboveUpperThrehsold"};
	public static final String[] NORMALIZATION_METHODS = {"none", "normalized", "normalizedSigned"};

	private String[] finalHeader;
	private String[] finalIDs;
	private double[][] database;
	private int numBelowLowerThreshold;
	private int numAboveUpperThreshold;
	private Logger log;

	public static void parse(String dir, String filename, String idColName, String[] phenos, String transform, double sdThreshold, boolean winsorize, boolean remove, boolean makeResids, boolean afterResids, boolean inverseNormalize, String covars, String idFile, boolean matchIdOrder, boolean plinkFormat, boolean pedFormat, boolean fastFormat, boolean excludeMissingValues, boolean variablesAllInOneFile, String extras, String[] outputs, boolean finalHeader, boolean addintercept, boolean sort, boolean zscore, boolean signZ, Logger log) {
		if (phenos == null) {
			log.reportError("Error - phenos is null");
			return;
		} else if (outputs == null) {
			log.reportError("Error - outputs is null");
			return;
		} else if (phenos.length != outputs.length) {
			log.reportError("Error - number of phenos is not equal to number of outputs");
			return;
		} else {
			for (int i = 0; i < phenos.length; i++) {
				parse(dir, filename, idColName, phenos[i], transform, sdThreshold, winsorize, remove, makeResids, afterResids, inverseNormalize, covars, idFile, matchIdOrder, plinkFormat, pedFormat, fastFormat, excludeMissingValues, variablesAllInOneFile, extras, outputs[i], finalHeader, addintercept, sort, zscore, signZ, log);
			}
		}
	}

	public static void parse(String dir, String filename, String idColName, String pheno, String transform, double sdThreshold, boolean winsorize, boolean remove, boolean makeResids, boolean afterResids, boolean inverseNormalize, String covarList, String idFile, boolean matchIdOrder, boolean plinkFormat, boolean pedFormat, boolean fastFormat, boolean excludeMissingValues, boolean variablesAllInOneFile, String extras, String outFile, boolean finalHeader, boolean addintercept, boolean sort, boolean zscore, boolean signZ, Logger log) {
		PhenoPrep prep;
		String[] covars;
		
		if (outFile == null) {
			outFile = pheno+"_out.csv";
			log.reportError("Warning - no output filename specified using [pheno]_out.csv ("+outFile+")");
		}

		log.report("Processing pheno: " + pheno + "\tout: " + outFile);

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
		
		if (makeResids && covarList == null) {
			log.reportError("Error - you have selected to make residuals without specifying any covariates to regress out for phenotype '"+pheno+"'; aborting");
			return;
		}
		
		if ((plinkFormat || fastFormat) && (idFile == null || !idFile.toLowerCase().endsWith(".fam"))) {
			log.reportError("Error - you have selected to make a PLINK or FAST formatted file with FID/IID, but have not provided a .fam file");
			return;
		}

		if (matchIdOrder && sort) {
			log.reportError("Error - you have selected both to match IDs order with another source and to sort IDs by ascending order");
			return;
		}

		// needed for Emmax, but not for anything else
//		if (variablesAllInOneFile) {
//			plinkFormat = true;
//		}

		if (covarList == null) {
			covars = new String[0];
		} else {
			covars = covarList.split(",");
		}

		prep = new PhenoPrep(dir+filename, idFile==null? null : dir+idFile, idColName, pheno, covars, log);
		
		if (prep.failed()) {
			log.report("Error - PhenoPrep failed for "+pheno);
			ext.waitForResponse();
			return;
		}
		
		if (transform != null && !transform.equals("none")) {
			if (!prep.transform(transform)) {
				return;
			}
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
		
		if (zscore) {
			prep.zscore(signZ);
		}
		
		if (extras != null) {
			prep.addExtraColumns(idColName, extras);
		}
		
		if (matchIdOrder) {
			if (idFile == null) {
				log.reportError("Error - match was selected, but no ID file was provided, skippping this step");
			} else {
				prep.matchIdOrder(idFile);
			}
		} else if (sort) {
			prep.sort();
		}

		prep.writeFinalFile(dir+outFile, plinkFormat, pedFormat, fastFormat, excludeMissingValues, variablesAllInOneFile, idFile, finalHeader);
		prep.summarizeCentralMoments(idFile);
	}
	
	private void summarizeCentralMoments(String idFile) {
		PrintWriter writer;
		boolean exists;
		double[] trait, ages;
		int[] males;
		
		exists = Files.exists("summary_stats.txt");
		try {
			writer = new PrintWriter(new FileWriter("summary_stats.txt", true));
			if (!exists) {
				writer.println(Array.toStr(SUMMARY_INFO_HEADER));
			}
			trait = Matrix.extractColumn(database, 0);
			if (ext.indexOfStr("Male", finalHeader) >= 0) {
				males = Array.toIntArray(Matrix.extractColumn(database, ext.indexOfStr("Male", finalHeader)));
			} else {
				males = null;
			}
			if (ext.indexOfStr("Age", finalHeader, false, false) >= 0) {
				ages = Matrix.extractColumn(database, ext.indexOfStr("Age", finalHeader, false, false));
			} else {
				ages = null;
			}
			writer.println((idFile == null?"All":ext.replaceAllWith(ext.rootOf(idFile), "_keeps", ""))+"\t"+finalHeader[0]+"\t"+ext.formDeci(Array.mean(trait), 4, false)+"\t"+ext.formDeci(Array.median(trait), 4, false)+"\t"+ext.formDeci(Array.stdev(trait), 4, false)+"\t"+ext.formDeci(Array.min(trait), 4, false)+"\t"+ext.formDeci(Array.max(trait), 4, false)
					+(males==null?"\t.\t.":"\t"+(males.length-Array.sum(males))+"\t"+Array.sum(males))
					+(ages==null?"\t.\t.\t.\t.":"\t"+ext.formDeci(Array.mean(ages), 4, false)+"\t"+ext.formDeci(Array.median(ages), 4, false)+"\t"+ext.formDeci(Array.stdev(ages), 4, false)+"\t"+ext.formDeci(Array.min(ages), 4, false)+"\t"+ext.formDeci(Array.max(ages), 4, false))
					+"\t"+(numBelowLowerThreshold < 0 ? "NA":numBelowLowerThreshold)+"\t"+(numAboveUpperThreshold < 0 ? "NA":numAboveUpperThreshold));
						
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + "summary_stats.txt");
			log.reportException(e);
		}
		
	}

	public PhenoPrep(String filename, String idFile, String idColName, String pheno, String[] covars, Logger log) {
		BufferedReader reader;
		String[] line;
		String temp;
		int[] indices;
		int idIndex;
		String[] idsWithDNA;
		String id;
		String delimiter;
		Vector<String> vIDs;
		Vector<double[]> vData;
		double[] data;
		boolean use;
		String[] header;

		this.log = log;
		numBelowLowerThreshold = -1;
		numAboveUpperThreshold = -1;
		
		if (idFile == null) {
			idsWithDNA = null;
		} else if (!Files.exists(idFile)) {
			log.reportError("ID file '"+idFile+"' does not exist; exiting");
			return;
		} else if (idFile.toLowerCase().endsWith(".fam")) {
			idsWithDNA = HashVec.loadFileToStringArray(idFile, false, new int[] { 1 }, false);
		} else {
			idsWithDNA = HashVec.loadFileToStringArray(idFile, false, new int[] { 0 }, false);
		}

		vIDs = new Vector<String>();
		vData = new Vector<double[]>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			temp = reader.readLine();
			delimiter = ext.determineDelimiter(temp);
			header = temp.split(delimiter, -1);
			idIndex = ext.indexOfStr(idColName, header, false, true, log, true);
			if (idIndex < 0) {
				log.reportError("Error - could not find specified id: "+idColName);
				reader.close();
				return;
			}
			indices = ext.indexFactors(Array.insertStringAt(pheno, covars, 0), header, false, log, true, true);
			finalHeader = Array.subArray(header, indices);
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.split(delimiter, -1);
				id = line[idIndex];
				if (idsWithDNA == null || ext.indexOfStr(id, idsWithDNA) >= 0) {
					line = Array.subArray(line, indices);
					use = true;
					for (int i = 0; i < line.length; i++) {
						if (ext.isMissingValue(line[i])) {
							use = false;
							break;
						}
					}
					if (use) {
						vIDs.add(id);
						data = Array.toDoubleArray(line);
						if (data == null) {
							log.reportError("Error - failed to parse data for "+id);
						}
						vData.add(data);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			log.reportException(fnfe);
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			log.reportException(ioe);
			return;
		}
		
		finalIDs = Array.toStringArray(vIDs);
		database = Matrix.toDoubleArrays(vData);

		if (finalIDs.length == 0) {
			log.reportError("Error - there are no indiviudals present in the final dataset"+(idFile != null?"; check the ids file to make sure the same set of IDs were used in both input files":""));
		}
	}
	
	public boolean failed() {
		return finalIDs == null;
	}
	
	public double[][] getDatabase() {
		return database;
	}

	public String[] getFinalIDs() {
		return finalIDs;
	}

	public boolean transform(String transform) {
		double[] data;
		int count;

		data = Matrix.extractColumn(database, 0);

		count = Array.countIf(Array.toStringArray(data), "0.0");
		if (count > 0 && (transform.equalsIgnoreCase("ln") || transform.equalsIgnoreCase("log10"))) {
			log.reportError("There "+(count == 1?"is one zero value":" are zero values")+", which will cause the "+transform+" transformation to fail; aborting");
			return false;
		}
		if (Array.min(data) < 0 && (transform.equalsIgnoreCase("ln") || transform.equalsIgnoreCase("log10") || transform.equalsIgnoreCase("sqrt"))) {
			log.reportError("Negative values will cause the "+transform+" transformation to fail; aborting");
			return false;
		}
		
		if (transform == null) {
			return true;
		} else if (transform.equalsIgnoreCase("ln")) {
			data = Transformations.naturalLogTransform(data);
		} else if (transform.equalsIgnoreCase("log10")) {
			data = Transformations.log10Transform(data);
		} else if (transform.equalsIgnoreCase("sqrt")) {
			data = Transformations.sqrtTransform(data);
		} else {
			log.reportError("Error - unknown transform: '"+transform+"'");
			return false;
		}
		
		for (int i = 0; i < data.length; i++) {
			database[i][0] = data[i];
		}
		
		return true;
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
		
		numBelowLowerThreshold = 0;
		numAboveUpperThreshold = 0;		
		if (winsorize) {
			for (int i = 0; i < data.length; i++) {
				if (data[i] < lowerThreshold) {
					numBelowLowerThreshold++;
					database[i][0] = lowerThreshold;
				} else if (data[i] > upperThreshold) {
					numAboveUpperThreshold++;
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
				if (data[i] < lowerThreshold) {
					numBelowLowerThreshold++;
					rowsToUse[i] = false;
				} else if (data[i] > upperThreshold) {
					numAboveUpperThreshold++;
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

	public void zscore(boolean signZ) {
		Matrix.overwriteColumn(database, 0, signZ?Array.normalizeSigned(Matrix.extractColumn(database, 0)):Array.normalize(Matrix.extractColumn(database, 0)), log);
	}

	public void computeResiduals() {
		LeastSquares reg;
		double[] deps, resids;
		double[][] indeps;
		
		deps = Matrix.extractColumn(database, 0);
		indeps = Matrix.extractColumns(database, Array.subArray(Array.intArray(database[0].length), 1));
		
		reg = new LeastSquares(deps, indeps, null, false, true);
		
		if (reg.analysisFailed()) {
			log.reportError("Error performing the regression model; check for collinearity if there are no other warnings above");
			ext.waitForResponse();
			System.exit(1);
		}
		
		resids = reg.getResiduals();
		
		if (deps.length != resids.length) {
			log.reportError("Error - lost a few rows in regression model; aborting");
			ext.waitForResponse();
			System.exit(1);
		}
		
		database = Matrix.toMatrix(resids);
		finalHeader = new String[] {finalHeader[0]};		
	}

	public void sort() {
		Arrays.sort(finalIDs);
	}

	public void addExtraColumns(String idColName, String extras) {
		String[] line;
		String temp;
		int[] indices;
		int idIndex;
		String delimiter;
		Vector<String> vIDs;
		Vector<double[]> vData;
		double[] data;
		boolean use;
		String[] header, newFinalHeader;
		boolean commaDelimitedFile;
		Hashtable<String, String> hash;

		temp = Files.getFirstNLinesOfFile(extras, 1, log)[0];
		delimiter = ext.determineDelimiter(temp);
		commaDelimitedFile = delimiter.equals(",");
		header = Files.getHeaderOfFile(extras, log);
		idIndex = ext.indexOfStr(idColName, header);
		if (idIndex == -1) {
			log.reportError("Error - extras file '"+extras+"' does not contain the same id linker ("+idColName+") as the main file; aborting all");
			System.exit(1);
		}
		indices = ext.indexFactors(Array.removeFromArray(header, idIndex), header, true, log, true, true);
		hash = HashVec.loadFileToHashString(extras, new int[] {idIndex}, indices, commaDelimitedFile, "\t", true, false, false);
		
		newFinalHeader = new String[finalHeader.length+indices.length] ;
		for (int i = 0; i < finalHeader.length; i++) {
			newFinalHeader[i] = finalHeader[i];
		}
		for (int i = 0; i < indices.length; i++) {
			newFinalHeader[finalHeader.length+i] = header[indices[i]];
		}

		vIDs = new Vector<String>();
		vData = new Vector<double[]>();
		for (int i = 0; i < finalIDs.length; i++) {
			if (hash.containsKey(finalIDs[i])) {
				line = hash.get(finalIDs[i]).split("[\\s]+");
				use = true;
				for (int j = 0; j < line.length; j++) {
					if (ext.isMissingValue(line[j])) {
						use = false;
					}
				}
				if (use) {
					vIDs.add(finalIDs[i]);
					data = new double[finalHeader.length+line.length];
					for (int j = 0; j < finalHeader.length; j++) {
						data[j] = database[i][j];
					}
					for (int j = 0; j < line.length; j++) {
						data[finalHeader.length+j] = Double.parseDouble(line[j]);
					}
					vData.add(data);
				}
			}
		}
		
		finalHeader = newFinalHeader;
		finalIDs = Array.toStringArray(vIDs);
		database = Matrix.toDoubleArrays(vData);
	}
	
	public void writeFinalFile(String filename, boolean plinkFormat, boolean pedFormat, boolean fastFormat, boolean excludeMissingValues, boolean variablesAllInOneFile, String idFile, boolean printFinalHeader)	{
		Hashtable<String, String> hash;
		PrintWriter writer;
		String delimiter;
		
		delimiter = Files.suggestDelimiter(filename, log);

		if (plinkFormat || pedFormat || fastFormat) {
			if (idFile == null || !idFile.toLowerCase().endsWith(".fam")) {
				log.reportError("Error - cannot export to plink format without an idFile.fam to lookup the FIDs");
				return;
			}
			
			if (pedFormat || fastFormat) {
				hash = HashVec.loadFileToHashString(idFile, new int[] {1}, new int[] {0, 1, 2, 3, 4}, false, delimiter, false, false, false);
			} else {
				hash = HashVec.loadFileToHashString(idFile, new int[] {1}, new int[] {0, 1}, false, delimiter, false, false, false);
			}
			if (variablesAllInOneFile || fastFormat) {
				try {
					writer = new PrintWriter(new FileWriter(filename));
					if (printFinalHeader) {
						if (fastFormat) {
							writer.println("#Fam_ID"+delimiter+"Ind_ID"+delimiter+"Dad_ID"+delimiter+"Mom_ID"+delimiter+"Sex"+delimiter+"Phenotype"+delimiter+Array.toStr(Array.subArray(finalHeader, 1), delimiter));
						} else {
							writer.println("FID"+delimiter+"IID"+delimiter+(pedFormat?"FA"+delimiter+"MO"+delimiter+"SEX"+delimiter:"")+ Array.toStr(finalHeader, delimiter));
						}
					}
					for (int i = 0; i < finalIDs.length; i++) {
						if (hash.containsKey(finalIDs[i])) {
							if (!excludeMissingValues || !Array.containsMissingValue(database[i])) {
								writer.println(hash.get(finalIDs[i]) + delimiter + Array.toStr(database[i], -1, -1, delimiter));
							}
						} else {
							log.report("Error - there was no record of " + finalIDs[i] + " in " + idFile + "; so no FID can be determined");
							writer.close();
							return;
						}
					}
					writer.close();
				} catch (Exception e) {
					log.reportError("Error writing to " + filename);
					e.printStackTrace();
				}
			} else {
				try {
					writer = new PrintWriter(new FileWriter(ext.addToRoot(filename, "_pheno")));
					if (printFinalHeader) {
						writer.println("FID"+delimiter+"IID"+delimiter+(pedFormat?"FA"+delimiter+"MO"+delimiter+"SEX"+delimiter:"")+finalHeader[0]);
					}
					for (int j = 0; j < finalIDs.length; j++) {
						if (hash.containsKey(finalIDs[j])) {
							if (!excludeMissingValues || !Array.containsMissingValue(database[j])) {
								writer.println(hash.get(finalIDs[j]) + delimiter + database[j][0]);
							}
						}
						else {
							log.report("Error - there was no record of " + finalIDs[j] + " in " + idFile + "; so no FID can be determined");
							writer.close();
							return;
						}
					}
					writer.close();
				} catch (Exception e) {
					log.reportError("Error writing to " + ext.addToRoot(filename, "_pheno"));
					e.printStackTrace();
				}
				if (finalHeader.length > 1) {
					try {
						writer = new PrintWriter(new FileWriter(ext.addToRoot(filename, "_covars")));
						if (printFinalHeader) {
							writer.println("FID"+delimiter+"IID"+delimiter + Array.toStr(Array.subArray(finalHeader, 1), delimiter));
						}
						for (int k = 0; k < finalIDs.length; k++) {
							if (hash.containsKey(finalIDs[k])) {
								if (!excludeMissingValues || !Array.containsMissingValue(database[k])) {
									writer.println(hash.get(finalIDs[k]).split(delimiter, -1)[0] + delimiter + finalIDs[k] + delimiter + Array.toStr(Array.subArray(database[k], 1), -1, -1, delimiter));
								}
							} else {
								log.report("Error - there was no record of " + finalIDs[k] + " in " + idFile + "; so no FID can be determined");
								writer.close();
								return;
							}
						}
						writer.close();
					}
					catch (Exception e) {
						log.reportError("Error writing to " + ext.addToRoot(filename, "_pheno"));
						e.printStackTrace();
					}
				}
			}
		} else {
			try {
				writer = new PrintWriter(new FileWriter(filename));
				if (printFinalHeader) {
					writer.println("id" + delimiter + Array.toStr(finalHeader, delimiter));
				}
				for (int m = 0; m < finalIDs.length; m++) {
					if (!excludeMissingValues || !Array.containsMissingValue(database[m])) {
						writer.println(finalIDs[m] + delimiter + Array.toStr(database[m], -1, -1, delimiter));
					}
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + filename);
				e.printStackTrace();
			}
		}
	}
	
	private void matchIdOrder(String idFile) {
		double[][] data;
		String[] newIDs;
		int[] idIndices;
		int index;
		
		if (idFile.endsWith(".fam")) {
			index = 1;
		} else {
			index = 0;
		}
		
		newIDs = HashVec.loadFileToStringArray(idFile, false, new int[] {index}, false);
		idIndices = new int[newIDs.length];
		for (int j = 0; j < newIDs.length; j++) {
			idIndices[j] = ext.indexOfStr(newIDs[j], finalIDs);
		}
		
		data = new double[newIDs.length][finalHeader.length];
		for (int j = 0; j < newIDs.length; j++) {
			if (idIndices[j] == -1) {
				data[j] = Array.doubleArray(finalHeader.length, Double.NaN);
			} else {
				data[j] = database[idIndices[j]];
			}
		}
		
		finalIDs = newIDs;
		database = data;
	}
	
	public static void fromParameters(String filename, Logger log) {
		Vector<String> params;

		params = Files.parseControlFile(filename, "PhenoPrep", new String[] {
				"dir=",
				"# name of input file",
				"file=input.txt",
				"# column name of the ID in the input file",
				"id=IID",
				"# phenotype column name",
				"pheno=outcomeVariable",
				"# covariate column names separated by a comma",
				"covar=Age,Sex,Site1,Site2",
				"# name of output file",
				"out=output.dat",
				"# transformation to apply to phenotype (current options are ln, log10, sqrt, or null for none)",
				"transform=none",
				"# name of file with IDs to use (e.g., that are in a genotype file); must be a plink .fam file if we are creating PLINK formatted files; otherwise, only the first column is used",
				"ids=plink.fam",
				"# winsorize phenotype (yes/no)",
				"winsorize=false",
				"# remove outliers (yes/no)",
				"remove=false",
				"# threshold in standard deviation units at which to winsorize or remove outliers",
				"sdThreshold=3.0",
				"# generate residuals instead of including covariates (yes/no)",
				"makeResids=false",
				"# winsorize/remove outliers after generating residuals (yes/no)",
				"afterResids=false",
				"# normalization of the final phenotype",
				"zscore=false",
				"# normalization of the final phenotype using sign-specific standard deviations",
				"signZ=false",
				"# inverse quantile normalize the final phenotype (e.g., after residuals are created if that is selected)",
				"inverseNormalize=false",
				"# name of file containing extra variables to include in final file but not in the outlier calculations; uncomment to use",
				"# (use this to add things like PCs, when you want to include all data in the SD/outlier calculation but only retain those you'll analyze)",
				"# extras=PrincipalComponentsFile.txt",
				"# match the order of the IDs in the idFile and the final file, using NA for missing data",
				"match=false",
				"# sort the IDs in the final file",
				"sort=false",
				"# output using FID and IID; FID is obtained from the ID file, which must have a .fam extension",
				"plinkFormat=false",
				"# use the .fam file to create .ped files instead (source must have a .fam extension)",
				"pedFormat=false",
				"# remove NaN values from the final file, even if matching, etc.",
				"excludeMissing=false",
				"# output using FID and IID, same as above, but have all variables in one file",
				"variablesAllInOneFile=false",
				"# whether to include a header row with the final file(s)",
				"finalHeader=true"
			}, log);

		if (params != null) {
			params.add("log=" + log.getFilename());
			main(Array.toStringArray(params));
		}
	}
	
	public static void summarizeFromParameters(String filename, Logger log) {
		Vector<String> params;

		params = Files.parseControlFile(filename, "bestTransformation", new String[] {
				"dir=",
				"# column name of the ID in the input file",
				"id=IID",
				"# phenotype names (requires a [phenoName].csv file as can be created by PhenoPrep)",
				"pheno=pheno1,pheno2,pheno3",
				"# covariate column names separated by a comma",
				"covar=Age,Sex,Site1,Site2",
				"# normalization of the final phenotype (0=none; 1=also normalization; 2=also normalization using sign-specific standard deviations)",
				"normalization=2",
			}, log);

		if (params != null) {
			params.add("-summarizeAll");
			params.add("log=" + log.getFilename());
			main(Array.toStringArray(params));
		}
	}
	
	public static void summarizeAll(String dir, String idColName, String phenosCommaDelimited, String covarsCommaDelimited, int normalization, String idFile) {
		PrintWriter writer;
		String[] phenos, transforms;
		Logger log;
		boolean winsorize, remove, makeResids, afterResids, inverseNormalize;
		String outFile;
		String[] rawData;
		double[] data;
		double mean, stdev, skewness, kurtosis;
		boolean normalize, normSigned;
		
		log = new Logger(dir+"summarizeAll.log");
		try {
			writer = new PrintWriter(new FileWriter(dir+"phenoSummary.xln"));
			writer.println("Trait\tshorthand\ttransform\twinsorize\tremoveOutliers\tmakeResiduals\tafterMakingResidualsDealWithOutliers\tnormalization\tN\tmean\tstdev\tskewness\tkurtosis\t'=SUM(ABS(SKEW)+ABS(KURT))");
		
			phenos = phenosCommaDelimited.split(",");
			
			transforms = new String[] {null, "ln", "sqrt"};
			
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
							for (int norm = 0; norm <= normalization; norm++) {
								normalize = false;
								normSigned = false;
								
								if (norm > 0) {
									normalize = true;
								}
								if (norm > 1) {
									normSigned = true;
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
								if (normalize) {
									outFile += "_"+NORMALIZATION_METHODS[norm];
								}
								log.report(outFile);
								outFile += ".csv";
								if (!Files.exists(dir+outFile)) {
									PhenoPrep.parse(dir, phenos[i]+".csv", idColName, phenos[i], transforms[j], 3.0, winsorize, remove, makeResids, afterResids, inverseNormalize, covarsCommaDelimited, idFile, false, false, false, false, true, true, null, outFile, true, false, false, normalize, normSigned, log);
								}
								if (Files.exists(dir+outFile)) {
									rawData = HashVec.loadFileToStringArray(dir+outFile, false, true, new int[] {1}, false, false, Files.determineDelimiter(dir+outFile, log));
									rawData = Array.removeFromArray(rawData, ext.MISSING_VALUES);
									data = Array.toDoubleArray(rawData);
									mean = Array.mean(data);
									stdev = Array.stdev(data);
									skewness = Array.skewness(data);
									kurtosis = Array.kurtosis(data);
									writer.println(phenos[i]+"\t"+ext.rootOf(outFile)+"\t"+transforms[j]+"\t"+winsorize+"\t"+remove+"\t"+makeResids+"\t"+afterResids+"\t"+(NORMALIZATION_METHODS[norm])+"\t"+data.length+"\t"+mean+"\t"+stdev+"\t"+skewness+"\t"+kurtosis+"\t"+(Math.abs(skewness)+Math.abs(kurtosis)));
								} else {
									writer.println(phenos[i]+"\t"+ext.rootOf(outFile)+"\t"+transforms[j]+"\t"+winsorize+"\t"+remove+"\t"+makeResids+"\t"+afterResids+"\t"+(NORMALIZATION_METHODS[norm])+"\tfailed\tfailed\tfailed\tfailed\tfailed\tfailed");
								}
							}
						}
					}
				}
			}
			
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + dir+"phenoSummary.xln");
			log.reportException(e);
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "";
		String filename = "pheno.csv";
		String logfile = null;
		Logger log;
		double sdThreshold = 3.0;
		String idColName = "id";
		String phenos = null;
		String transform = null;
		String covarsCommaDelimited = null;
		String idFile = null;
		String outFile = null;
		String[] outputs = null;
		boolean winsorize = false;
		boolean remove = false;
		boolean makeResids = false;
		boolean afterResids = false;
		boolean inverseNormalize = false;
		String extras = null;
		boolean matchIdOrder = false;
		boolean plinkFormat = false;
		boolean finalHeader = true;
		boolean variablesAllInOneFile = false;
		boolean summarizeAll = false;
		boolean addintercept = false;
		boolean sort = false;
		boolean zscore = false;
		boolean signZ = false;
		int normalization = 2;
		boolean pedFormat = false;
		boolean fastFormat = true;
		boolean excludeMissingValues = false;

//		dir = "";
//		filename = "N:/statgen/BOSS/phenotypes/PhenoPrep/taste/Taste_withOtherIDs.xln";
//		idColName = "GWAS_ID";
//		phenos = "AnyProp";
//		covarsCommaDelimited = "Age,Sex";
//		outFile = "N:/statgen/BOSS/phenotypes/PhenoPrep/taste/anyProp_1.dat";
//		finalHeader = true;
		
//		filename = "ARIC_Whites_WBC.csv";
//		pheno = "WBC";
//		outFile = "ARIC_EA_WBC.csv";
//		idFile = "keeps.txt";
//		transform = "ln";
//		covars = new String[] {"Age", "Male", "CenterF", "CenterJ", "CenterM"};
//		idFile = "EA_keeps.dat";
//		winsorize = true;

//		dir = "D:/SkatMeta/results_hemostasis/";
//		filename = "pheno_F7.csv";
//		pheno = "F7";
//		outFile = "pheno_F7_winsorize_3sd.csv";
		
//		covars = new String[] {"V1AGE01", "Sex"};
//		winsorize = false;
//		remove = false;
//		makeResids = false;
//		afterResids = false;
		

		
//		dir = "D:/LITE/CHARGE-S/aric_wex_freeze3/testOutliers/";
//		phenos = "Fibrinogen,F7,F8,vWF";
//		covarsCommaDelimited = "V1AGE01,Sex,CenterF,CenterM";
		
//		dir = "D:/ExomeChip/ARIC_primary/CompareTransformations/";
//		phenos = "Hct,Hb,MCHC,MCV,RBC,MCH,RDW,WBC_TOTAL,WBC_NEUTRO,WBC_MONO,WBC_LYMPH,WBC_EOS,WBC_BASO";
//		covarsCommaDelimited = "Age,Male";

//		dir = "D:/ExomeChip/APTT_and_ProteinC/";
////		filename = "coag_gwas_recoded2.txt";	// not relevant, data needs to be in separate files named [pheno].csv
//		idColName = "ID";
//		phenos = "ARIC_AA_APTT,ARIC_AA_ProteinC,ARIC_EA_APTT,ARIC_EA_ProteinC";
//		covarsCommaDelimited = "Age,Male";
//		
//		summarizeAll(dir, idColName, phenos, covarsCommaDelimited, null);
//		System.exit(1);

//		dir = "D:/LITE/phenotypes/1000G_FAST_runs/raw_PhenoScan/";
//		idColName = "ID";
//		phenos = "ARIC_AA_APTT,ARIC_AA_ProteinC,ARIC_EA_APTT,ARIC_EA_ProteinC";
//		covarsCommaDelimited = "Age,Male,PC1,PC2";
//		
//		summarizeAll(dir, idColName, phenos, covarsCommaDelimited, 0, null);
//		System.exit(1);


		String usage = "\n"+
				"gwas.PhenoPrep requires 0-1 arguments\n" +
				"	 (0) name of directory (i.e. dir=" + dir + " (default))\n" +
				"	 (1) name of input file (i.e. file=" + filename + " (default))\n" +
				"	 (2) id column name in input file (i.e. id=" + idColName + " (default))\n" +
				"	 (3) phenotype column name(s) (i.e. pheno=" + phenos + " (default; comma to delimit multiple phenos))\n" +
				"	 (4) covariate column name(s) (i.e. covar=Age,Sex,Site1,Site2 (not the default))\n" +
				"	 (5) name of file with IDs to use (e.g., that are in a genotype file) (i.e. ids=" + idFile + " (default; set to null to include all rows with complete data))\n" +
				"	 (6) name of output file (i.e. out=" + outFile + " (default))\n" +
				"	 (7) transformation to apply to phenotype (i.e. transform=" + transform + " (default; current options are ln, log10, sqrt, or null for none))\n" +
				"	 (8) winsorize phenotype (i.e. winsorize=" + winsorize + " (default))\n" +
				"	 (9) remove outliers (i.e. remove=" + remove + " (default))\n" +
				"	(10) threshold in standard deviation units at which to winsorize or remove outliers (i.e. sdThreshold=" + sdThreshold + " (default))\n" +
				"	(11) generate residuals instead of including covariates (i.e. makeResids=" + makeResids + " (default))\n" +
				"	(12) winsorize/remove outliers after generating residuals (i.e. afterResids=" + afterResids + " (default))\n" +
				"	(13) inverse quantile normalize the final phenotype (i.e. inverseNormalize=" + inverseNormalize + " (default))\n" +
				"	(14) name of file containing extra variables to include in final file but not in outlier calculations (i.e. extras=" + extras + " (default))\n" +
				"			 (use this to add things like PCs, when you want to include all data in the SD/outlier calculation but only retain those you'll analyze)\n" +
				"	(15) match the order of the IDs in the idFile and the final file, using NA for missing data (i.e. match=" + matchIdOrder + " (default))\n" +
				"	(16) output using FID and IID; FID is obtained from the ID file, which must have a .fam extension (i.e. plinkFormat=" + plinkFormat + " (default))\n" +
				"	(17) use PLINK FID and IID from .fam file, but have all variables in one file (i.e. variablesAllInOneFile=" + variablesAllInOneFile + " (default))\n" +
				"	(18) use the .fam file to create .ped files instead (i.e. pedFormat="+pedFormat+" (default))\n" +
				"	(19) use the .fam file to create FAST format .trait files instead (i.e. fastFormat="+fastFormat+" (default))\n" +
				"	(20) remove NaN values (i.e. excludeMissing="+excludeMissingValues+" (default))\n" +
				"	(21) include a header with the final file(s) (i.e. finalHeader=" + finalHeader + " (default))\n" +
				"	(22) add an intercept variable (value equals 1 constantly) as the 3rd column (i.e. addintercept=" + addintercept + " (default))\n" +
				"	(23) sort the output by the 1st column (i.e. sort=" + sort + " (default))\n" +
				"	(24) (optional) name of log file to write to (i.e. log=[pheno].log (default))\n" +
				"	(25) convert final phenotype into a z-score (i.e. zscore=" + zscore + " (default))\n" +
				"	(26) z-score uses positive-only (mirrored) and negative-only (mirrored) distributions to compute the standard deviation for the z-scores (i.e. signZ="+signZ+" (default))\n" +
				"  OR:\n" +
				"	 (6) run all possible combinations of transformations/outliers to assess normality (i.e. -summarizeAll (not the default))\n" +
				"	 (7) include normaliation transformations (i.e. normalization="+normalization+" (default; 0=none, 1=standard, 2=standard and sign-specific stdevs))\n" +
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
			} else if (args[i].startsWith("id=")) {
				idColName = ext.parseStringArg(args[i], "[firstColumn]");
				numArgs--;
			} else if (args[i].startsWith("pheno=")) {
				phenos = ext.parseStringArg(args[i], "missingPhenotype");
				numArgs--;
			} else if (args[i].startsWith("covar=")) {
				covarsCommaDelimited = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outFile = ext.parseStringArg(args[i], "output_file.csv");
				if (outFile.contains(",")) {
					outputs = outFile.split(",");
					outFile = null;
				}
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
			} else if (args[i].startsWith("extras=")) {
				extras = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("match=")) {
				matchIdOrder = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("plinkFormat=")) {
				plinkFormat = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("pedFormat=")) {
				pedFormat = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("fastFormat=")) {
				fastFormat = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("excludeMissing=")) {
				excludeMissingValues = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("variablesAllInOneFile=")) {
				variablesAllInOneFile = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("finalHeader=")) {
				finalHeader = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-summarizeAll")) {
				summarizeAll=true;
				numArgs--;
			} else if (args[i].startsWith("addintercept=")) {
				addintercept = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("sort=")) {
				sort = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("zscore=")) {
				zscore = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("signZ=")) {
				signZ = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("normalization=")) {
				normalization = ext.parseIntArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
				ext.waitForResponse();
			}
		}
		
		if (logfile == null) {
			logfile = dir+ext.replaceWithLinuxSafeCharacters(outFile==null?phenos+"_out":ext.rootOf(outFile), true)+".log";
		}
		log = new Logger(logfile);
		
		if (numArgs != 0) {
			log.reportError(usage);
			System.exit(1);
		}
		
		try {
			if (summarizeAll) {
				summarizeAll(dir, idColName, phenos, covarsCommaDelimited, normalization, idFile);
			} else if (phenos.contains(",")) {
				parse(dir, filename, idColName, phenos.split(","), transform, sdThreshold, winsorize, remove, makeResids, afterResids, inverseNormalize, covarsCommaDelimited, idFile, matchIdOrder, plinkFormat, pedFormat, fastFormat, excludeMissingValues, variablesAllInOneFile, extras, outputs, finalHeader, addintercept, sort, zscore, signZ, log);
			} else {
				parse(dir, filename, idColName, phenos, transform, sdThreshold, winsorize, remove, makeResids, afterResids, inverseNormalize, covarsCommaDelimited, idFile, matchIdOrder, plinkFormat, pedFormat, fastFormat, excludeMissingValues, variablesAllInOneFile, extras, outFile, finalHeader, addintercept, sort, zscore, signZ, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}



