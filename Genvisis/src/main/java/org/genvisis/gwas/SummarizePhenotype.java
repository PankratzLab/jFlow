package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class SummarizePhenotype {
	public static final String PARAMETERS_FILE_NAME = "parameters.txt";
	public static final String[][] FILETYPE_DELIMINATOR = {	{".csv", ","}, {".xln", "\t"},
																													{".txt", "\t"}};
	public static final String PARAMETERS_FILE_DELIMITER_FOR_DATAFILES = "[;|]";
	public static final String PARAMETERS_FILE_DELIMITER_FOR_VARIABLES_IN_MODEL = "[~=+*]";
	public static final String[] ID_NAMES = {	"IndividualID", "ID", "IID", "UID", "UniqueID", "IndID",
																						"Sample"};
	public static final String FILENAME_EXTIONSION_TO_OUTPUT_PHENOTYPE_DATA_WITH = ".xln";
	public static final String[][] DEFAULT_PARAMS =
																								{	{"Description", "statOperation", "columnLabel",
																									"columnCriteria", "parameter"},
																									{	"Trait - Minimum", "0 percentile", "trait",
																										"trait=ValidDouble", "-blank", "-sf=2"},
																									{	"Trait - 1st Quartile", "25 percentile", "trait",
																										"trait=ValidDouble", "-blank", "-sf=2"},
																									{	"Trait - Median", "50 percentile", "trait",
																										"trait=ValidDouble", "-blank", "-sf=2"},
																									{	"Trait - Mean", "mean", "trait",
																										"trait=ValidDouble", "-blank", "-sf=2"},
																									{	"Trait - 3rd Quartile", "75 percentile", "trait",
																										"trait=ValidDouble", "-blank", "-sf=2"},
																									{	"Trait - maximum", "100 percentile", "trait",
																										"trait=ValidDouble", "-blank", "-sf=2"},
																									{	"Trait- standard deviation", "stdev", "trait",
																										"trait=ValidDouble", "-blank", "-sf=2"},
																									{	"N", "count", "Female", "Female=ValidDouble",
																										"-blank", "-sf=0"},
																									{	"N - Female", "count", "Female", "Female=1",
																										"-blank", "-sf=0"},
																									{	"Age - Mean", "mean", "Age", "Age=ValidDouble",
																										"-blank", "-sf=2"},
																									{	"Age - standard deviation", "stdev", "Age",
																										"Age=ValidDouble", "-blank", "-sf=2",
																										"-stdev"},};

	public static String runDescriptiveStat(double[] data, String statOperation, boolean includeStdev,
																					boolean useBlankForNull, boolean showPercent, int sf,
																					Logger log) {
		double mean, stdev = Double.NaN;
		int percentile;
		String result;

		if (statOperation.equalsIgnoreCase("mean")) {
			mean = (Array.sum(data) / data.length);
			if (includeStdev) {
				stdev = Array.stdev(data);
			}
			result = (mean != 0	? (showPercent	? ext.formDeci(mean * 100, sf) + "%"
																					: ext.formDeci(mean, sf))
														+ (includeStdev	? " ("
																								+ (showPercent	? ext.formDeci(stdev * 100, sf) + "%"
																															: ext.formDeci(stdev, sf))
																							+ ")"
																						: "")
													: (useBlankForNull ? "" : "."));

		} else if (statOperation.equalsIgnoreCase("stdev")) {
			stdev = Array.stdev(data);
			result = (!Double.isNaN(stdev))	? (showPercent	? ext.formDeci(stdev * 100, sf) + "%"
																											: ext.formDeci(stdev, sf))
																			: (useBlankForNull ? "" : ".");

		} else if (statOperation.equalsIgnoreCase("count")) {
			result = Array.sum(data) > 0 ? data.length + "" : (useBlankForNull ? "" : ".");

		} else if (statOperation.contains("percentile")) {
			percentile = Integer.parseInt(statOperation.split(" ")[0].trim());
			Arrays.sort(data);
			if (percentile == 100) {
				result = (data.length == 0	? (useBlankForNull ? "" : "0")
																		: ext.formDeci(data[data.length - 1], sf));
			} else {
				result = (data.length == 0	? (useBlankForNull ? "" : "0")
																		: ext.formDeci(	data[(int) (((double) percentile / 100)
																																* data.length + .5)],
																										sf));
			}

		} else {
			log.reportError("Unrecognized statisitcal operation: " + statOperation);
			result = null;
		}

		return result;
	}


	public static String[] parseStatisticsFromParameters(	Vector<String[]> modelData,
																												String[] columnLabels,
																												String[][] parameters, Logger log) {
		boolean isStdev, isBlank, isPercent, isToInclude;
		int sf, numParameters, columnIndex, theSize;
		Vector<Integer> filters;
		Vector<String> filterText;
		Vector<Double> dataVec;
		double[] data;
		String[] line, result;

		numParameters = parameters.length;
		result = new String[numParameters];
		// data = new double[modelData.size()];
		// for (int i = 0; i < data.length; i ++) {
		// data[i] = Double.parseDouble(modelData.get(i));
		// }

		for (int i = 0; i < numParameters; i++) {
			isStdev = false;
			isBlank = false;
			isPercent = false;
			sf = 4;

			filters = new Vector<Integer>();
			filterText = new Vector<String>();
			for (int j = 3; j < parameters[i].length; j++) {
				if (parameters[i][j].equals("-stdev")) {
					isStdev = true;
				} else if (parameters[i][j].startsWith("-sf=")) {
					sf = ext.parseIntArg(parameters[i][j]);
				} else if (parameters[i][j].equals("-blank")) {
					isBlank = true;
				} else if (parameters[i][j].equals("-percent")) {
					isPercent = true;
				} else {
					line = parameters[i][j].split("=");
					if (line.length < 2 || line[1].equals("")) {
						log.reportError("Parameter contains invalid filter: "	+ parameters[i][j]
														+ ". Should include an '=' sign.");
						return null;
					}
					// if (line[0].equalsIgnoreCase("trait")) {
					// columnIndex = 1;
					// } else {
					columnIndex = ext.indexOfStr(line[0].trim(), columnLabels, false, true, log, true);
					if (columnIndex == -1) {
						// columnIndex = ext.indexOfStr(line[0].trim(), columnLabels, false, false, log, true);
						// if (columnIndex == -1) {
						log.reportError("Cannot find ANY match for the column specified by the parameter ("
														+ parameters[i][j] + ") in the data. Ignored this parameter.");
						return null;
						// } else {
						// log.reportError("Cannot find EXACT match for the column specified by the parameter ("
						// + parameters[i][j] + ") in the data. Instead, used the following column from the
						// data: " + columnLabels[columnIndex]);
						// }
					}
					// }
					line[1] = line[1].trim();
					if (!line[1].equalsIgnoreCase("validDouble")	&& !line[1].equalsIgnoreCase("validInteger")
							&& !ext.isValidInteger(line[1]) && !ext.isValidDouble(line[1])) {
						log.reportError("Unrecognized parameter: "	+ parameters[i][j]
														+ ". Expecting 'variableName = validDouble' or 'variableName = validInteger'.");
					}
					filters.add(columnIndex);
					filterText.add(line[1].trim());
				}
			}

			// if (parameters[i][2].equalsIgnoreCase("trait")) {
			// columnIndex = 1;
			// } else {
			columnIndex = ext.indexOfStr(parameters[i][2].trim(), columnLabels, false, true, log, true);
			if (columnIndex == -1) {
				// columnIndex = ext.indexOfStr(parameters[i][2].trim(), columnLabels, false, false, log,
				// true);
				// if (columnIndex == -1) {
				log.reportError("Cannot find ANY match for the column specified by the parameter ("
												+ parameters[i][2] + ") in the data. Ignored this parameter.");
				return null;
				// } else {
				// log.reportError("Cannot find EXACT match for the column specified by the parameter (" +
				// parameters[i][2] + ") in the data. Instead, used the following column from the data: " +
				// columnLabels[columnIndex]);
				// }
			}
			// }
			theSize = modelData.size();
			dataVec = new Vector<Double>(theSize);
			for (int j = 0; j < theSize; j++) {
				line = modelData.get(j);
				isToInclude = true;
				for (int k = 0; k < filters.size(); k++) {
					if (filterText.get(k).equalsIgnoreCase("validDouble")) {
						if (!ext.isValidDouble(line[filters.get(k)])) {
							isToInclude = false;
							break;
						}
					} else if (filterText.get(k).equalsIgnoreCase("validInteger")) {
						if (!ext.isValidInteger(line[filters.get(k)])) {
							isToInclude = false;
							break;
						}
					} else {
						if (!line[filters.get(k)].equalsIgnoreCase(filterText.get(k))) {
							isToInclude = false;
							break;
						}
					}
				}
				if (isToInclude) {
					dataVec.add(Double.parseDouble(line[columnIndex]));
				}
			}

			theSize = dataVec.size();
			data = new double[theSize];
			for (int j = 0; j < theSize; j++) {
				data[j] = dataVec.get(j);
			}

			result[i] = runDescriptiveStat(data, parameters[i][1], isStdev, isBlank, isPercent, sf, log);

		}
		return result;
	}

	public static void exportData(String filename, String[] header, Vector<String[]> data) {
		PrintWriter writer;
		int size;

		size = data.size();
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println(Array.toStr(header));
			for (int i = 0; i < size; i++) {
				writer.println(Array.toStr(data.get(i)));
			}
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static Vector<String[]> selectColumnsAndRows(Hashtable<String, String[]> traitDataTmp,
																											int[] columnsToKeep, String idListFileName,
																											Logger log) {
		BufferedReader reader;
		Vector<String[]> currentData;
		String[] temp2, temp3;
		String iid;

		currentData = new Vector<String[]>();
		try {
			reader = new BufferedReader(new FileReader(idListFileName));
			while (reader.ready()) {
				iid = reader.readLine();
				if (traitDataTmp.containsKey(iid)) {
					temp2 = new String[columnsToKeep.length];
					temp3 = traitDataTmp.get(iid);
					for (int i = 0; i < temp2.length; i++) {
						temp2[i] = temp3[columnsToKeep[i]];
					}
					currentData.add(temp2);
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return currentData;
	}

	public static Hashtable<String, String[]> loadData(String dataFileName, Logger log) {
		Hashtable<String, String[]> data;
		String delimiter;
		String[] line;
		BufferedReader reader;

		data = new Hashtable<String, String[]>();
		delimiter = null;
		for (String[] element : FILETYPE_DELIMINATOR) {
			if (dataFileName.toLowerCase().endsWith(element[0])) {
				delimiter = element[1];
			}
		}
		if (delimiter != null) {
			try {
				reader = new BufferedReader(new FileReader(dataFileName));
				line = reader.readLine().split(delimiter);
				if (ext.indexOfStr(line[0].trim(), ID_NAMES, false, true, log, true) == -1) {
					log.reportError("Cannot identify the 1st column of the following file as an ID '"
													+ dataFileName + "'");
				}
				for (String element : line) {
				}
				data.put("headerOfTheFile", line);

				while (reader.ready()) {
					line = reader.readLine().split(delimiter);
					data.put(line[0], line);
				}
				reader.close();

			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else {
			log.reportError("Unrecognized Extension in file name '" + dataFileName + "'");
		}

		return data;
	}


	/*
	 * Trim both the blank spaces and the blank tabs.
	 */
	public static String[] splitAndTrim(String line, String delimiter) {
		String[] tmp, result;
		int count;

		tmp = line.toLowerCase().split(delimiter);
		count = tmp.length;
		for (int i = 0; i < tmp.length; i++) {
			tmp[i] = tmp[i].trim();

			if (tmp[i].equals("")) {
				count--;
			}
		}

		if (count != tmp.length) {
			result = new String[count];
			count = 0;
			for (int i = 0; i < result.length; i++) {
				while (tmp[i + count].equals("")) {
					count++;
				}
				result[i] = tmp[i + count];
			}
		} else {
			result = tmp;
		}

		return result;
	}

	public static void parse2(String modelListFileName, String outFileName, Logger log) {
		String wkDir;
		String[] line = null;
		BufferedReader reader;
		Hashtable<String, String[][]> models;
		Vector<String> modelList;
		Hashtable<String, Hashtable<String, String[]>> data;
		Hashtable<String, String[]> currentDataTmp = null;
		String[] columnLabels = null;
		Vector<String[]> currentData;
		Vector<String[]> parametersVec;
		String[][] parameters;
		String[] variablesInModel = null;
		int index;
		int[] columnsToKeep = null;
		String[][] result;
		PrintWriter writer;
		String modelName, idListFileName;
		String[] dataFileNames;

		wkDir = ext.parseDirectoryOfFile(modelListFileName);
		try {
			models = new Hashtable<String, String[][]>();
			modelList = new Vector<String>();
			reader = new BufferedReader(new FileReader(modelListFileName));
			reader.readLine();
			while (reader.ready()) {
				line = splitAndTrim(reader.readLine(), "[\t]");
				if (!line[0].startsWith("#")) {
					if (line.length != 4) {
						log.reportError("Unrecognized format in " + modelListFileName + ": " + line.toString());
						System.exit(1);
					} else {
						modelName = line[0].trim().toLowerCase();
						models.put(	modelName,
												new String[][] {splitAndTrim(	line[1],
																											PARAMETERS_FILE_DELIMITER_FOR_DATAFILES),
																				splitAndTrim(	line[2],
																											PARAMETERS_FILE_DELIMITER_FOR_VARIABLES_IN_MODEL),
																				new String[] {line[3].toLowerCase()}});
						modelList.add(modelName);
					}
				}
			}
			reader.close();

			parametersVec = new Vector<String[]>();
			reader = new BufferedReader(new FileReader(wkDir + PARAMETERS_FILE_NAME));
			reader.readLine();
			while (reader.ready()) {
				parametersVec.add(reader.readLine().split("\t"));
			}
			reader.close();
			parameters = new String[parametersVec.size()][];
			for (int i = 0; i < parameters.length; i++) {
				parameters[i] = parametersVec.get(i);
			}

			data = new Hashtable<String, Hashtable<String, String[]>>();
			result = new String[modelList.size()][parametersVec.size()];
			for (int i = 0; i < result.length; i++) {
				modelName = modelList.get(i);
				log.report("Processing model: " + modelName);
				dataFileNames = models.get(modelName)[0];
				variablesInModel = models.get(modelName)[1];
				idListFileName = models.get(modelName)[2][0];
				currentData = new Vector<String[]>();

				for (int j = 0; j < dataFileNames.length; j++) {
					if (!data.containsKey(dataFileNames[j])) {
						data.put(dataFileNames[j], loadData(wkDir + dataFileNames[j], log));
					}
					currentDataTmp = data.get(dataFileNames[j]);
					columnLabels = currentDataTmp.get("headerOfTheFile");
					columnsToKeep = new int[variablesInModel.length + 1];
					for (int k = 1; k < columnsToKeep.length; k++) {
						index = ext.indexOfStr(variablesInModel[k - 1], columnLabels, false, true, log, true);
						if (index > -1) {
							columnsToKeep[k] = index;
						}
					}
					currentData.addAll(selectColumnsAndRows(currentDataTmp, columnsToKeep,
																									wkDir + idListFileName, log));
				}

				line = new String[variablesInModel.length + 1];
				line[0] = "iid";
				for (int j = 1; j < line.length; j++) {
					line[j] = variablesInModel[j - 1];
				}
				exportData(wkDir	+ modelName + FILENAME_EXTIONSION_TO_OUTPUT_PHENOTYPE_DATA_WITH, line,
										currentData);
				// result[i] = parseStatisticsFromParameters(currentData, variablesInModel, parameters,
				// log);
				result[i] = parseStatisticsFromParameters(currentData, line,
																									adjustParameters(parameters, line, log), log);
			}

			writer = new PrintWriter(new FileWriter(wkDir + outFileName));
			for (int i = 0; i < result.length; i++) {
				writer.print("\t" + modelList.get(i));
			}
			writer.println("");

			for (int i = 0; i < result[0].length; i++) {
				writer.print(parametersVec.get(i)[0]);
				for (String[] element : result) {
					if (element != null) {
						writer.print("\t" + element[i]);
					} else {
						writer.print("\t");
					}
				}
				writer.println("");
			}
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static String[][] adjustParameters(String[][] parameters, String[] columnLabels,
																						Logger log) {
		String[][] result;
		String[] line;
		int index;

		result = Arrays.copyOf(parameters, parameters.length);
		for (int i = 0; i < result.length; i++) {
			for (int j = 2; j < result[i].length; j++) {
				if (!result[i][j].startsWith("-")) {
					line = result[i][j].split("=");
					if (ext.indexOfStr(line[0], columnLabels, false, true, log, true) == -1) {
						if (line[0].toLowerCase().contains("trait")) {
							result[i][j] = result[i][j].replaceAll("trait", columnLabels[1]);
						} else {
							index = ext.indexOfStr(line[0], columnLabels, false, false, log, true);
							if (index > -1) {
								log.report("No exact match for parameter "	+ result[i][j]
														+ " in column labels. Instead, used " + columnLabels[index]);
							} else {
								log.reportError("No exact match for parameter "	+ result[i][j]
																+ " in the column labels. Program exit on error.");
								System.exit(0);
							}
							result[i][j] = result[i][j].replaceAll(line[0], columnLabels[index]);
						}
					}
				}
			}
		}

		return result;
	}

	public static void summarizeFromParameters(String filename, Logger log) {
		List<String> params;
		// String wkDir;
		String line;
		String[] defaults;
		String[][] parameters = null;
		String[] columnLabels;
		String phenoDataFileName = null, resultFileName = null;
		int index;


		defaults = new String[DEFAULT_PARAMS.length + 2];
		defaults[0] = "in=phenoDataFile.csv";
		defaults[1] = "out=statisticsResults.xln";
		for (int i = 0; i < DEFAULT_PARAMS.length; i++) {
			defaults[i + 2] = Array.toStr(DEFAULT_PARAMS[i]);
		}

		params = Files.parseControlFile(filename, "descriptive", defaults, log);

		if (params == null) {
			if (!Files.exists(PARAMETERS_FILE_NAME)) {
				Files.writeMatrix(DEFAULT_PARAMS, PARAMETERS_FILE_NAME, "\t");
			}
			log.report("Populated with default/example parameters");
			return;
		}


		// del = Files.determineDelimiter(filename, log);
		// del = ext.determineDelimiter(str);

		index = params.size();
		parameters = new String[index - 2][];
		for (int i = 0; i < index; i++) {
			line = params.get(i);
			if (!line.startsWith("#")) {
				if (line.startsWith("in=")) {
					phenoDataFileName = line.split("=")[1];
				} else if (line.startsWith("out=")) {
					resultFileName = line.split("=")[1];
					// } else if (line.startsWith("log=")) {
					// logFileName = line.split("=")[1];

					// TODO How to detect unrecognized lines?
				} else {
					parameters[i] = line.split("\t");
				}
			}
		}

		columnLabels = Files.getHeaderOfFile(filename, "\t", log);
		parameters = adjustParameters(parameters, columnLabels, log);

		Files.generateTables(	resultFileName, new String[] {phenoDataFileName},
													new String[] {ext.rootOf(phenoDataFileName)}, parameters, log);
	}


	// @SuppressWarnings("null")
	// public static void summarizeFromParameters(String filename, Logger log) {
	// Vector<String> params;
	// String wkDir;
	// String[] line, trav;
	// String[] defaults;
	// String[][] parameters = null;
	// BufferedReader reader;
	// String[] columnLabels;
	// int index;
	//
	// defaults = new String[DEFAULT_PARAMS.length+2];
	//// defaults[0] = "inputFile.csv";
	//// defaults[1] = "outputTable.xln";
	// for (int i = 0; i < DEFAULT_PARAMS.length; i++) {
	//// defaults[i+2] = Array.toStr(DEFAULT_PARAMS[i]);
	// defaults[i] = Array.toStr(DEFAULT_PARAMS[i]);
	// }
	//
	// wkDir = ext.parseDirectoryOfFile(filename);
	// params = Files.parseControlFile(wkDir + PARAMETERS_FILE_NAME, "Description", defaults, log);
	//
	// if (params == null) {
	// if (! Files.exists(PARAMETERS_FILE_NAME)) {
	// Files.writeMatrix(DEFAULT_PARAMS, PARAMETERS_FILE_NAME, "\t");
	// }
	// } else {
	//// params.add("log=" + log.getFilename());
	// }
	//
	//
	//// del = Files.determineDelimiter(filename, log);
	//// del = ext.determineDelimiter(str);
	//
	// try {
	// reader = new BufferedReader(new FileReader(filename));
	// columnLabels = reader.readLine().split("\t");
	// reader.close();
	//
	// parameters = new String[params.size()][];
	//// for (int i = 0; i < defaults.length; i++) {
	// for (int i = 0; i < parameters.length; i++) {
	// line = params.get(i).toLowerCase().split("\t");
	// for (int j = 2; j < line.length; j++) {
	// if (! line[j].startsWith("-")) {
	// trav = line[j].split("=");
	// if (ext.indexOfStr(trav[0], columnLabels, false, true, log, true) == -1) {
	// if (trav[0].toLowerCase().contains("trait")) {
	// line[j] = line[j].replaceAll("trait", columnLabels[1]);
	// } else {
	// index = ext.indexOfStr(trav[0], columnLabels, false, false, log, true);
	// if (index > -1) {
	// log.report("The parameter " + line[j] + " does not match any column label in the dataset.
	// Instead, used " + columnLabels[index]);
	// } else {
	// log.reportError("The parameter " + line[j] + " does not match any column label in the dataset.
	// Program exit on error.");
	// System.exit(0);
	// }
	// line[j] = line[j].replaceAll(trav[0], columnLabels[index]);
	// }
	// }
	// }
	// }
	// parameters[i] = line;
	// }
	//
	//// parseStatisticsFromParameters(modelData, columnLabels, parameters, log);
	// Files.generateTables(wkDir + "result_" + ext.rootOf(filename) + ".xln", new String[]
	// {filename}, new String[] {ext.rootOf(filename)}, parameters, log);
	// } catch (IOException e) {
	// e.printStackTrace();
	// }
	// }

	public static void filesFromParameters(String modelListFileName, Logger log) {
		List<String> params;

		params = Files.parseControlFile(modelListFileName, "phenotype",
																		new String[] {"#modelName\tdataFile\tvariableList\tIdListFile",
																									"CARDIA_blacks_Fibrinogen_y5\tBlacks.csv\tCL6FIBR_ln ~ EX3_AGE + Female\tdemoBlacks_Fibrinogen_y5.txt",
																									"CARDIA_whites_F7_y5\tWhites.csv\tCLEFACT7 ~ EX3_AGE + Female\tdemoWhites_F7_y5.txt"},
																		log);

		if (params == null) {
			if (!Files.exists(PARAMETERS_FILE_NAME)) {
				Files.writeMatrix(DEFAULT_PARAMS, PARAMETERS_FILE_NAME, "\t");
			}
		} else {
			params.add("log=" + log.getFilename());
			main(Array.toStringArray(params));
		}
	}

	// public static void filesFromParameters(String modelListFileName, Logger log) {
	// Vector<String> params;
	//
	// params = Files.parseControlFile(modelListFileName, "phenotype", new String[] {
	// "#modelName\tdataFile\tvariableList\tIdListFile",
	// "CARDIA_blacks_Fibrinogen_y5\tBlacks.csv\tCL6FIBR_ln ~ EX3_AGE +
	// Female\tdemoBlacks_Fibrinogen_y5.txt",
	// "CARDIA_whites_F7_y5\tWhites.csv\tCLEFACT7 ~ EX3_AGE + Female\tdemoWhites_F7_y5.txt"
	// }, log);
	//
	// if (params == null) {
	// if (!Files.exists(PARAMETERS_FILE_NAME)) {
	// Files.writeMatrix(DEFAULT_PARAMS, PARAMETERS_FILE_NAME, "\t");
	// }
	// } else {
	// params.add("log=" + log.getFilename());
	// main(Array.toStringArray(params));
	// }
	// }

	public static void main(String[] args) {
		int numArgs = args.length;
		// String filename = "pheno.csv";
		String modelList = "N:/statgen/ExomeChip/Hemostasis/CARDIA/models.txt";
		// String idList = null;
		String logfile = null;
		String outFile = "result.xln";
		Logger log;

		String usage = "\n"	+ "gwas.SummarizePhenotype requires 0-1 arguments\n"
										+ "   (1) name of model list file (i.e. models=" + modelList + " (default))\n"
										+ "   (2) name of output file (i.e. out=" + outFile + " (default))\n"
										+ "   (3) name of log file (i.e. log=mylog.txt (not the default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("models=")) {
				modelList = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("out=")) {
				outFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(logfile);
			log.report("Started summarizing genotypes.");
			parse2(modelList, outFile, log);
			// summarizeFromParameters("N:/statgen/ExomeChip/Hemostasis/CARDIA/cardia_blacks_f7_y5.xln",
			// log);
			log.report("Program finished.");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
