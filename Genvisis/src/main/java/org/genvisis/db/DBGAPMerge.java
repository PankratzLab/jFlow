package org.genvisis.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeSet;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class DBGAPMerge {

	Logger log;

	public static class FileSet {
		String dir;
		String dataFile;
		String dataDict;
		String varReport;

		private List<DataColumn> dataDefs;
		private HashMap<String, Integer> defMap;
		ArrayList<String> ids = new ArrayList<String>();
		HashMap<String, String[]> idDataMap = new HashMap<String, String[]>();

		public int getIndexOfDataColumn(String key) {
			return defMap.get(key);
		}

		public void setDataDefs(List<DataColumn> readDataDict) {
			dataDefs = readDataDict;
			defMap = new HashMap<String, Integer>();
			for (int i = 0; i < dataDefs.size(); i++) {
				DataColumn dc = dataDefs.get(i);
				String key = dc.varName + ";" + dc.varID + ";" + dc.table;;
				defMap.put(key, i);
			}
		}
	}

	public static class DataColumn {
		String source;
		String study;
		String table;
		// String consentGroup;
		String varID;
		String varName;
		String varDesc;
		HashMap<String, String> varValueDescriptions = new HashMap<String, String>();
		String varUnit;
		String comment;
	}

	private static final String DICT_SUFF = "data_dict.xml";
	private static final String VAR_RPT_SUFF = "var_report.xml";
	private static final String DATA_SUFF = "txt.gz";

	private static final String MISSING_DATA = ".";

	public ArrayList<FileSet> discoverFiles(String dir) {
		ArrayList<FileSet> fss = new ArrayList<DBGAPMerge.FileSet>();

		String[] allFiles = new File(dir).list();

		for (String gzFile : allFiles) {
			if (!gzFile.endsWith(DATA_SUFF)) {
				continue;
			}
			if (Files.countLines(dir + gzFile, 0) < 5) {
				System.err.println("Warning - file " + gzFile + " will be skipped due to lack of data.");
				continue;
			}
			String[] parts = gzFile.split("\\.");

			String key = parts[0] + "." + parts[1] + "." + parts[2] + "." + parts[3];

			String dict = null;
			String varRpt = null;
			for (String f : allFiles) {
				if (f.endsWith(".gz") || !f.startsWith(key)) {
					continue;
				}
				if (f.endsWith(DICT_SUFF)) {
					if (dict == null) {
						dict = f;
					} else {
						System.err.println("Error - duplicate data dictionary file detected: " + dict + " / "
															 + f);
					}
				} else if (f.endsWith(VAR_RPT_SUFF)) {
					if (varRpt == null) {
						varRpt = f;
					} else {
						System.err.println("Error - duplicate variable report file detected: " + varRpt + " / "
															 + f);
					}
				}
			}

			if (dict != null && varRpt != null) {
				FileSet fs = new FileSet();
				fs.dir = dir;
				fs.dataFile = gzFile;
				fs.dataDict = dict;
				fs.varReport = varRpt;
				fss.add(fs);
			} else {
				if (dict == null) {
					log.reportError("Error - no data dictionary file found for data file: " + gzFile);
				}
				if (varRpt == null) {
					log.reportError("Error - no variable report file found for data file: " + gzFile);
				}
			}
		}

		return fss;
	}

	public ArrayList<DataColumn> readDataDict(FileSet fs) throws ParserConfigurationException,
																											 SAXException, IOException {
		ArrayList<DataColumn> dict = new ArrayList<DBGAPMerge.DataColumn>();

		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		Document doc = builder.parse(new File(fs.dir + fs.dataDict));
		doc.getDocumentElement().normalize();

		Node tableNode = doc.getElementsByTagName("data_table").item(0);
		String table = tableNode.getAttributes().getNamedItem("id").getNodeValue();
		String study = tableNode.getAttributes().getNamedItem("study_id").getNodeValue();

		NodeList variables = doc.getElementsByTagName("variable");

		for (int i = 0, count = variables.getLength(); i < count; i++) {
			Node varNode = variables.item(i);

			DataColumn dc = new DataColumn();
			Node vnName = varNode.getAttributes().getNamedItem("id");
			dc.source = fs.dir + fs.dataFile;
			dc.varID = vnName.getNodeValue().trim();
			dc.table = table.trim();
			dc.study = study.trim();

			NodeList nlChildren = varNode.getChildNodes();
			for (int c = 0, countCh = nlChildren.getLength(); c < countCh; c++) {
				Node chNode = nlChildren.item(c);
				String chNm = chNode.getNodeName();
				if (chNm.equals("#text")) {
					continue;
				}
				String nodeNm = chNode.getNodeName();
				String nodeVal = chNode.getTextContent();
				if (nodeNm.equals("name")) {
					dc.varName = nodeVal;
				} else if (nodeNm.equals("description")) {
					dc.varDesc = nodeVal;
				} else if (nodeNm.equals("type")) {

				} else if (nodeNm.equals("value")) {
					NamedNodeMap nnm = chNode.getAttributes();
					Node cdNode = nnm.getNamedItem("code");
					if (cdNode != null) {
						dc.varValueDescriptions.put(cdNode.getNodeValue(), nodeVal);
					} else {
						System.out.println(ext.getTime() + "]\tValue without Code: " + nodeVal);
					}
				} else if (nodeNm.equals("unit")) {
					dc.varUnit = nodeVal;
				} else if (nodeNm.equals("comment")) {
					dc.comment = nodeVal.replaceAll("\\n", " - ");
				} else {
					log.reportError("Unknown node: " + nodeNm + " = " + nodeVal);
				}

			}
			dict.add(dc);
		}

		return dict;
	}

	public void readDataFile(FileSet fs) throws IOException {
		String line;
		String[] parts;
		BufferedReader reader = Files.getAppropriateReader(fs.dir + fs.dataFile);

		boolean firstAfterHeader = false;
		while ((line = reader.readLine()) != null) {
			if (line.startsWith("##")) {
				// TODO read header line
				firstAfterHeader = true;
			} else if (line.startsWith("#") || line.equals("")) {
				continue; // skip unimportant header content
			} else {
				if (firstAfterHeader) {
					firstAfterHeader = false;
					continue;
				}
				parts = line.split("\t", -1);
				fs.ids.add(parts[0]);
				while (parts.length - 1 < fs.dataDefs.size()) {
					parts = ArrayUtils.addStrToArray(MISSING_DATA, parts);
				}
				fs.idDataMap.put(parts[0], ArrayUtils.subArray(parts, 1));
			}
		}
		reader.close();
		reader = null;
		System.gc();
	}

	public void run(String[] dirs, String outFile, String outMap, Logger log) {
		this.log = log;
		ArrayList<FileSet> files = new ArrayList<DBGAPMerge.FileSet>();
		for (String dir : dirs) {
			files.addAll(discoverFiles(ext.verifyDirFormat(dir)));
		}

		for (FileSet fs : files) {
			try {
				fs.setDataDefs(readDataDict(fs));
			} catch (ParserConfigurationException e) {
				log.reportException(e);
			} catch (SAXException e) {
				log.reportException(e);
			} catch (IOException e) {
				log.reportException(e);
			}
			if (fs.dataDefs != null) {
				try {
					readDataFile(fs);
				} catch (IOException e) {
					log.reportException(e);
				}
			}
		}

		TreeSet<String> idSet = new TreeSet<String>();
		for (FileSet fs : files) {
			idSet.addAll(fs.ids);
		}

		ArrayList<String> dataColumnKeys = new ArrayList<String>();
		HashMap<String, ArrayList<FileSet>> dataColumnMap = new HashMap<String, ArrayList<FileSet>>();
		int total = 0;
		for (FileSet fs : files) {
			for (DataColumn dc : fs.dataDefs) {
				total++;
				String key = dc.varName + ";" + dc.varID + ";" + dc.table;
				if (!dataColumnKeys.contains(key)) {
					dataColumnKeys.add(key);
				}
				ArrayList<FileSet> fss = dataColumnMap.get(key);
				if (fss == null) {
					fss = new ArrayList<DBGAPMerge.FileSet>();
					dataColumnMap.put(key, fss);
				}
				fss.add(fs);
			}
		}

		System.out.println("Writing " + dataColumnKeys.size() + " data columns out of " + total);

		PrintWriter writer = Files.getAppropriateWriter(outFile);
		StringBuilder lineOut = null;

		lineOut = new StringBuilder("DBGAP_ID");
		for (String key : dataColumnKeys) {
			lineOut.append("\t").append(key);
		}
		writer.println(lineOut.toString());

		for (String id : idSet) {
			lineOut = new StringBuilder(id);

			for (String key : dataColumnKeys) {
				ArrayList<FileSet> sources = dataColumnMap.get(key);
				String val = MISSING_DATA;
				for (FileSet fs : sources) {
					String[] data = fs.idDataMap.get(id);
					if (data == null) {
						continue;
					}
					int ind = fs.getIndexOfDataColumn(key);
					if (MISSING_DATA.equals(val) && !"".equals(data[ind])) {
						val = data[ind];
					}
				}
				lineOut.append("\t").append(val);
			}

			writer.println(lineOut.toString());
		}
		writer.flush();
		writer.close();

		writer = Files.getAppropriateWriter(outMap);
		writer.println("Source\tStudy\tTable\tVariable\tVariableID\tfinalColumnName\tcustomColumnName\tDescription\tUnits\tVariableMapping\tComment");
		for (String key : dataColumnKeys) {
			lineOut = new StringBuilder();

			FileSet fs = dataColumnMap.get(key).get(0); // TODO fix selection of fileset?

			DataColumn dc = fs.dataDefs.get(fs.getIndexOfDataColumn(key));

			lineOut.append(dc.source).append("\t");
			lineOut.append(dc.study).append("\t");
			lineOut.append(dc.table).append("\t");
			lineOut.append(dc.varName).append("\t");
			lineOut.append(dc.varID).append("\t");
			lineOut.append(key).append("\t");
			lineOut.append(".").append("\t");
			lineOut.append(dc.varDesc).append("\t");
			lineOut.append(dc.varUnit == null ? "." : dc.varUnit).append("\t");
			if (dc.varValueDescriptions.isEmpty()) {
				lineOut.append(".").append("\t");
			} else {
				for (Entry<String, String> ent : dc.varValueDescriptions.entrySet()) {
					lineOut.append(ent.getKey()).append("=").append(ent.getValue()).append(";");
				}
			}
			if (dc.comment != null && !"".equals(dc.comment)) {
				lineOut.append("\t").append(dc.comment);
			} else {
				lineOut.append("\t").append(".");
			}

			writer.println(lineOut.toString());
		}
		writer.flush();
		writer.close();

	}

	public static void main(String[] args) {
		String out = "out";
		String outMap = "outMap";
		String dir = "dir";
		String logfile = "log";
		Logger log;

		CLI c = new CLI(DBGAPMerge.class);
		c.addArg(out, "Data Output Filename", true);
		c.addArg(outMap, "Map Output Filename", true);
		c.addArg(dir, "Input directory (or comma-separated list)", true);
		c.addArg(logfile, "Log file");

		c.parseWithExit(args);


		try {
			log = new Logger(c.get(logfile));
			(new DBGAPMerge()).run(c.get(dir).split(","), c.get(out), c.get(outMap), log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static class DBGapExtract {

		private static void extract(String varFile, String dataFile, String outputFile,
																String headIdent, Logger log) throws IOException {
			BufferedReader dataReader;
			PrintWriter writer;
			String[][] varData;
			String[] colsToLoad, parts;
			String line, delim;
			StringBuilder sb;

			varData = HashVec.loadFileToStringMatrix(varFile, true, new int[] {1, 2}); // ident
																																								 // col, repl
																																								 // col,
																																								 // ignore
																																								 // source
																																								 // and other
																																								 // cols
			colsToLoad = Matrix.extractColumn(varData, 0);

			dataReader = Files.getAppropriateReader(dataFile);
			line = dataReader.readLine();
			delim = ext.determineDelimiter(line);
			parts = line.split(delim, -1);
			int[] factors = ext.indexFactors(colsToLoad, parts, false);

			writer = Files.getAppropriateWriter(outputFile);
			sb = new StringBuilder();
			sb.append(headIdent);
			for (int i = 0; i < factors.length; i++) {
				if (factors[i] >= 0) {
					sb.append("\t").append(ext.isMissingValue(varData[i][1]) ? varData[i][0] : varData[i][1]);
				}
			}
			writer.println(sb.toString());

			while ((line = dataReader.readLine()) != null) {
				parts = line.split(delim, -1);
				sb = new StringBuilder();
				sb.append(parts[0]);
				for (int factor : factors) {
					if (factor >= 0) {
						sb.append("\t").append(parts[factor]);
					}
				}
				writer.println(sb.toString());
			}
			writer.flush();
			writer.close();
			dataReader.close();
		}

		public static void fromParameters(String filename, Logger log) {
			List<String> params;
			String[] args;

			params = Files.parseControlFile(filename, "dbgap",
																			new String[] {"variables=searchTerms.xln",
																										"dataFile=mergeOut.xln.gz",
																										"outputFile=searchedVariables.xln",
																										"head=IID",},
																			log);

			if (params != null) {
				args = params.toArray(new String[params.size()]);
				main(args);
			}
		}

		public static void main(String[] args) {
			int numArgs = args.length;
			String varFile = "searchTerms.xln";
			String dataFile = "mergeOut.xln.gz";
			String outputFile = "searchedVariables.xln";
			String headIdent = "IID";
			String logfile = null;
			Logger log;

			String usage = "\n"
										 + "one.ben.DBGAPExtract requires 4+ arguments\n"
										 + "   (1) Searched variables file (output from DBGapLookup) (i.e. variables="
										 + varFile
										 + " (default))\n"
										 + "   (2) Merged dbGap data file (data output from DBGapMerge ) (i.e. dataFile="
										 + dataFile + " (default))\n"
										 + "   (3) Extracted output filename (i.e. outputFile=" + outputFile
										 + " (default))\n" + "   (4) ID column name in output file (i.e. head="
										 + headIdent + " (default))\n" + "   (5) OPTIONAL: Log file (i.e. log="
										 + logfile + " (default))\n" + "";

			for (String arg : args) {
				if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
					System.err.println(usage);
					System.exit(1);
				} else if (arg.startsWith("variables=")) {
					varFile = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("dataFile=")) {
					dataFile = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("outputFile=")) {
					outputFile = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("head=")) {
					headIdent = arg.split("=")[1];
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
				extract(varFile, dataFile, outputFile, headIdent, log);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

	}

	public static class DBGapLookup {

		private static void search(String map, String[] search, int[] searchCols, int[] outputCols,
															 String output, String crf, Logger log) throws IOException {
			String line, delim, outputDelim;
			String[] parts;
			StringBuilder sb;

			BufferedReader reader = Files.getAppropriateReader(map);
			PrintWriter writer = Files.getAppropriateWriter(output);

			line = reader.readLine().trim();
			delim = ext.determineDelimiter(line);
			outputDelim = "\t";
			parts = line.split(delim, -1);

			sb = new StringBuilder();
			for (int outputCol : outputCols) {
				sb.append(parts[outputCol]).append(outputDelim);
			}
			for (int i = 0; i < search.length; i++) {
				sb.append(search[i]);
				if (i < search.length - 1) {
					sb.append(outputDelim);
				}
			}
			writer.println(sb.toString());

			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if ("".equals(line)) {
					continue;
				}
				parts = line.trim().split(delim, -1);

				int[] searchInds = ArrayUtils.intArray(search.length, -1);
				for (int s = 0; s < search.length; s++) {
					for (int col : searchCols) {
						searchInds[s] = Math.max(searchInds[s],
																		 parts[col].toLowerCase().indexOf(search[s].toLowerCase()));
					}
				}
				boolean foundAny = false;
				for (int s : searchInds) {
					if (s >= 0) {
						foundAny = true;
						break;
					}
				}
				if (foundAny) {
					sb = new StringBuilder();
					for (int outputCol : outputCols) {
						sb.append(parts[outputCol]).append(outputDelim);
					}
					for (int i = 0; i < searchInds.length; i++) {
						sb.append(searchInds[i] == -1 ? "." : searchInds[i]);
						if (i < search.length - 1) {
							sb.append(outputDelim);
						}
					}
					writer.println(sb.toString());
				}
			}
			reader.close();
			writer.flush();
			writer.close();
		}

		public static void fromParameters(String filename, Logger log) {
			List<String> params;
			String[] line, args;
			String mapFile;
			StringBuilder sb;

			params = Files.parseControlFile(filename,
																			"search",
																			new String[] {
																										"mergeMap.xln   searchCols=3,7,9,10 outputCols=0,5,6,3,7,9,10 out=searchTerms.xln  crf=dbgap.crf",
																										"# search terms below, one per line:"},
																			log);
			if (params != null) {
				line = params.remove(0).trim().split(PSF.Regex.GREEDY_WHITESPACE);

				mapFile = line[0];
				if (!Files.exists(mapFile)
						&& Files.exists(ext.verifyDirFormat(ext.parseDirectoryOfFile(filename)) + mapFile)) {
					mapFile = ext.verifyDirFormat(ext.parseDirectoryOfFile(filename)) + mapFile;
				}

				args = new String[line.length + (log.getFilename() == null ? 1 : 2)];
				args[0] = "map=" + mapFile;
				for (int i = 1; i < line.length; i++) {
					args[i] = line[i];
				}
				sb = new StringBuilder(params.get(0));
				for (int i = 1; i < params.size(); i++) {
					sb.append(",").append(params.get(i));
				}
				args[line.length] = "search=" + sb.toString();
				if (log.getFilename() != null) {
					args[line.length + 1] = "log=" + log.getFilename();
				}
				main(args);
			}
		}

		public static void main(String[] args) {
			int numArgs = args.length;

			String map = "mergeMap.xln";
			int[] searchCols = new int[] {3, 7, 9, 10, 11};
			int[] outputCols = new int[] {0, 5, 6, 3, 7, 9, 10, 11};
			String output = "searchTerms.xln";
			String crf = "dbgap.crf";
			String[] search = new String[] {"cancer", "height", "bmi"};

			String logfile = null;
			Logger log;

			String usage = "\n" + "one.ben.DBGAPMerge requires 6+ arguments\n"
										 + "   (1) dbGap Map File (from merged files) (i.e. map=" + map
										 + " (default))\n"
										 + "   (2) Column indices in which to search (i.e. searchCols="
										 + ArrayUtils.toStr(searchCols, ",") + " (default))\n"
										 + "   (3) Column indices from which to export data (i.e. outputCols="
										 + ArrayUtils.toStr(outputCols, ",") + " (default))\n"
										 + "   (4) Output File (i.e. out=" + output + " (default))\n"
										 + "   (5) Data export CRF filename (i.e. crf=" + crf + " (default))\n"
										 + "   (6) Search terms (i.e. search=" + ArrayUtils.toStr(search, ",")
										 + " (default))\n" + "   (7) OPTIONAL: Log file (i.e. log=" + logfile
										 + " (default))\n" + "";

			for (String arg : args) {
				if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
					System.err.println(usage);
					System.exit(1);
				} else if (arg.startsWith("map=")) {
					map = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("crf=")) {
					crf = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("out=")) {
					output = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("searchCols=")) {
					searchCols = ext.parseIntArrayArg(arg);
					numArgs--;
				} else if (arg.startsWith("outputCols=")) {
					outputCols = ext.parseIntArrayArg(arg);
					numArgs--;
				} else if (arg.startsWith("search=")) {
					search = ext.parseStringArrayArg(arg, search);
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
				search(map, search, searchCols, outputCols, output, crf, log);
				if (Files.exists(crf)) {
					log.reportError("Error - specified .CRF file already exists! -- " + crf);
				} else {
					Files.parseControlFile(crf, "dbgap",
																 new String[] {"variables=" + output, "dataFile=mergeOut.xln.gz",
																							 "outputFile=searchedVariables.xln", "head=IID",},
																 log);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

	}



}
