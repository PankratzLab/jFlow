package org.genvisis.cnv.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import javax.swing.JProgressBar;

import org.genvisis.cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import org.genvisis.cnv.manage.SourceFileParser;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class SourceFileHeaderData implements Serializable {
	/**
	 *
	 */
	private static final long serialVersionUID = -6906302109843776908L;

	String gsgtVersion = null;
	String processingDate = null;
	String content = null;
	int numSnps = -1;
	int totalSnps = -1;
	int numSamples = -1;
	int totalSamples = -1;
	int numFiles = -1;
	int currFile = -1;
	public int columnHeaderLineIndex = -1;

	public int colSampleIdent = -1;
	public int colSnpIdent = -1;
	public int colGenoAB1 = -1;
	public int colGenoAB2 = -1;
	public int colGeno1 = -1;
	public int colGeno2 = -1;
	public int colX = -1;
	public int colY = -1;
	public int colTheta = -1;
	public int colR = -1;
	public int colXRaw = -1;
	public int colYRaw = -1;
	public int colBAF = -1;
	public int colLRR = -1;
	public int colGC = -1;
	public String headerString = "";
	public String delimiter = "";
	public String[] cols = new String[0];

	// Order from ParseIllumina
	// 0 GC
	// 1 XRAW
	// 2 YRAW
	// 3 X
	// 4 Y
	// 5 Theta
	// 6 R
	// 7 B Allele Freq
	// 8 Log R Ratio

	private SourceFileHeaderData() {}

	public static SourceFileHeaderData parseHeader(String file, Logger log)	throws Elision,
																																					IOException {
		BufferedReader reader = Files.getAppropriateReader(file);
		String line = null;
		int lineCnt = 0;
		SourceFileHeaderData frhd = new SourceFileHeaderData();
		String delim = ",";
		while ((line = reader.readLine()) != null) {
			delim = ext.determineDelimiter(line, true); // TODO if file ends with .csv [or contains
																									// .csv?], can assume that delim is ','. Sim., if
																									// ends with .xln, can assume delim is '\t'
			if (",".equals(delim)) {
				delim = "[\\s]*,[\\s]*"; // ext.indexFactors doesn't call trim()
			}
			if ("[Data]".equals(line)	|| line.startsWith("rs") || line.toUpperCase().startsWith("SNP")
					|| ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line.split(delim), false, true,
															false, false)[0] != -1) {
				break;
			}
			String[] parts = line.trim().split(",");
			processLine(parts, frhd);
			lineCnt++;
		}
		if (!"[Data]".equals(line)
				&& !(line.startsWith("rs")	|| line.toUpperCase().startsWith("SNP")
							|| ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line.split(delim), false,
																	true, false, false)[0] != -1)) {
			log.reportError("Error - malformed or missing header.");
			throw new Elision(file);
		}
		if ("[Data]".equals(line)) {
			line = reader.readLine();
			delim = file.contains(".csv")	? "[\\s]*,[\\s]*"
																		: file.contains(".xln")	? "[\\s]*\t[\\s]*"
																														: ext.determineDelimiter(line, true);
			lineCnt++;
			if (!(line.startsWith("rs")	|| line.toUpperCase().startsWith("SNP")
						|| ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line.split(delim), false, true,
																false, false)[0] != -1)) {
				log.reportError("Error - malformed or missing header.  Header must start with 'rs' or 'SNP' or contain one of the following: "
												+ ArrayUtils.toStr(SourceFileParser.SNP_HEADER_OPTIONS[0]) + ".");
				throw new Elision(file);
			}
		}
		reader.close(); // done
		reader = null; // release
		String columnHeaders = line;
		delim = file.contains(".csv")	? SOURCE_FILE_DELIMITERS.COMMA.getDelimiter()
																	: file.contains(".xln")	? SOURCE_FILE_DELIMITERS.TAB.getDelimiter()
																													: SOURCE_FILE_DELIMITERS.getDelimiter(ext.determineDelimiter(	columnHeaders,
																																																												true))
																																									.getDelimiter();
		// if (",".equals(delim)) {
		// delim = "[\\s]*,[\\s]*";
		// }
		parseColumnsBestGuess(columnHeaders.split(delim), frhd);
		frhd.setSourceFileDelimiter(delim);
		frhd.setHeaderString(columnHeaders);
		frhd.columnHeaderLineIndex = lineCnt;
		return frhd;
	}

	private static final String[] SAMPLE_FIELD_ID = {"Sample ID", "Sample Name"};

	private static final String[] DATA_FIELDS_GC = {"GC Score", "GCscore", "confidence",
																									"confidenceScore"};
	private static final String[] DATA_FIELDS_XRAW = {"X Raw"};
	private static final String[] DATA_FIELDS_YRAW = {"Y Raw"};
	private static final String[] DATA_FIELDS_X = {	"X", "Xvalue", "Log Ratio", "intensity_1",
																									"Signal A"};
	private static final String[] DATA_FIELDS_Y = {	"Y", "Yvalue", "Strength", "intensity_2",
																									"Signal B"};
	private static final String[] DATA_FIELDS_THETA = {"Theta"};
	private static final String[] DATA_FIELDS_R = {"R"};
	private static final String[] DATA_FIELDS_BAF = {"B Allele Freq", "BAF"};
	private static final String[] DATA_FIELDS_LRR = {"Log R Ratio", "LRR"};
	private static final String[] GENOTYPE_FIELDS_A1_FOR = {"Allele1 - Forward", "Allele1",
																													"genotype1", "Allele1 - Top",
																													"Forward Strand Base Calls",
																													"Forced Call", "Forced Call Codes"};
	private static final String[] GENOTYPE_FIELDS_A2_FOR = {"Allele2 - Forward", "Allele B",
																													"genotype2", "Allele2 - Top",
																													"Forward Strand Base Calls",
																													"Forced Call", "Forced Call Codes"};
	private static final String[] GENOTYPE_FIELDS_A1_AB = {"Allele1 - AB", "Call Codes", "Call"};
	private static final String[] GENOTYPE_FIELDS_A2_AB = {"Allele2 - AB", "Call Codes", "Call"};

	private static final String[][] LOOKUP = {/* 0 */ SourceFileParser.SNP_HEADER_OPTIONS[0],
																						/* 1 */ SAMPLE_FIELD_ID, /* 2 */ DATA_FIELDS_GC,
																						/* 3 */ DATA_FIELDS_XRAW, /* 4 */ DATA_FIELDS_YRAW,
																						/* 5 */ DATA_FIELDS_X, /* 6 */ DATA_FIELDS_Y,
																						/* 7 */ DATA_FIELDS_THETA, /* 8 */ DATA_FIELDS_R,
																						/* 9 */ DATA_FIELDS_BAF, /* 10 */ DATA_FIELDS_LRR,
																						/* 11 */ GENOTYPE_FIELDS_A1_FOR,
																						/* 12 */ GENOTYPE_FIELDS_A2_FOR,
																						/* 13 */ GENOTYPE_FIELDS_A1_AB,
																						/* 14 */ GENOTYPE_FIELDS_A2_AB,};

	private static void parseColumnsBestGuess(String[] parts,
																						SourceFileHeaderData frhd) throws Elision {
		int[] indices = ext.indexFactors(LOOKUP, parts, false, true, false, false);
		if (indices[0] == -1) {
			throw new Elision("Error - missing SNP ID column");
		}
		frhd.cols = parts;
		frhd.colSnpIdent = indices[0];
		frhd.colSampleIdent = indices[1];
		// frhd.col_sampleIndex = indices[2];
		frhd.colGC = indices[2];
		frhd.colXRaw = indices[3];
		frhd.colYRaw = indices[4];
		frhd.colX = indices[5];
		frhd.colY = indices[6];
		frhd.colTheta = indices[7];
		frhd.colR = indices[8];
		frhd.colBAF = indices[9];
		frhd.colLRR = indices[10];
		frhd.colGeno1 = indices[11];
		frhd.colGeno2 = indices[12];
		frhd.colGenoAB1 = indices[13];
		frhd.colGenoAB2 = indices[14];
	}

	private static void processLine(String[] parts, SourceFileHeaderData frhd) {
		if ("[Header]".equals(parts[0])) {
			return;
		}
		if ("File".equals(parts[0])) {
			String[] fileParts = parts[parts.length - 1].split(" of ");
			if (ext.isValidInteger(fileParts[0])) {
				frhd.currFile = Integer.parseInt(fileParts[0]);
			} else {
				// TODO error
			}
			if (ext.isValidInteger(fileParts[1])) {
				frhd.numFiles = Integer.parseInt(fileParts[1]);
			} else {
				// TODO error
			}
		}
		if ("Total Samples".equals(parts[0])) {
			if (ext.isValidInteger(parts[parts.length - 1])) {
				frhd.totalSamples = Integer.parseInt(parts[parts.length - 1]);
			} else {
				// TODO error
			}
		}
		if ("Num Samples".equals(parts[0])) {
			if (ext.isValidInteger(parts[parts.length - 1])) {
				frhd.numSamples = Integer.parseInt(parts[parts.length - 1]);
			} else {
				// TODO error
			}
		}
		if ("Total SNPs".equals(parts[0])) {
			if (ext.isValidInteger(parts[parts.length - 1])) {
				frhd.totalSnps = Integer.parseInt(parts[parts.length - 1]);
			} else {
				// TODO error
			}
		}
		if ("Num SNPs".equals(parts[0])) {
			if (ext.isValidInteger(parts[parts.length - 1])) {
				frhd.numSnps = Integer.parseInt(parts[parts.length - 1]);
			} else {
				// TODO error
			}
		}
		if ("Content".equals(parts[0])) {
			frhd.content = parts[parts.length - 1];
		}
		if ("Processing Date".equals(parts[0])) {
			frhd.processingDate = parts[parts.length - 1];
		}
		if ("GSGT Version".equals(parts[0])) {
			frhd.gsgtVersion = parts[parts.length - 1];
		}
	}

	public static HashMap<String, SourceFileHeaderData> validate(	final String rawDir,
																																final String ext,
																																boolean fullValidation, Logger log,
																																JProgressBar progressBar) {
		String dir = rawDir.endsWith("/")
									|| rawDir.endsWith("\\")	? rawDir
																						: org.genvisis.common.ext.verifyDirFormat(rawDir);
		String[] possibleFiles = (new File(dir)).list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(ext);
			}
		});

		// No files found - not a valid project
		if (possibleFiles.length == 0) {
			log.report("Project validation failed: no project files found discovered.");
			return null;
		}

		if (progressBar != null) {
			progressBar.setVisible(true);
			progressBar.setMinimum(0);
			progressBar.setMaximum(possibleFiles.length);
			progressBar.setString(null);
			progressBar.setStringPainted(true);
		}

		boolean valid = false;
		HashMap<String, SourceFileHeaderData> headers = null;
		int progCnt = 0;
		try {
			headers = new HashMap<String, SourceFileHeaderData>();
			SourceFileHeaderData exemplar = null;
			for (String possFile : possibleFiles) {
				SourceFileHeaderData frhd;
				if (!fullValidation) {
					if (exemplar == null) {
						frhd = SourceFileHeaderData.parseHeader(dir + possFile, log);
						exemplar = frhd;
					} else {
						frhd = exemplar;
					}
				} else {
					frhd = SourceFileHeaderData.parseHeader(dir + possFile, log);
				}
				headers.put(possFile, frhd);
				if (progressBar != null) {
					progressBar.setValue(++progCnt);
				}
			}
			if (fullValidation) {
				if (progressBar != null) {
					progressBar.setIndeterminate(true);
					progressBar.setString("Verifying...");
				}
				String error = doFullValidation(headers, log);
				if (error != null) {
					throw new Elision(error);
				}
			}
			if (progressBar != null) {
				progressBar.setVisible(false);
			}
			valid = true;
		} catch (Elision e) {
			log.reportError(e.getMessage());
		} catch (IOException e) {
			log.reportException(e);
		}

		return valid ? headers : null;
	}

	public static String doFullValidation(HashMap<String, SourceFileHeaderData> headers, Logger log) {
		int cnt = headers.size();
		HashMap<Integer, ArrayList<String>> totSnpsSet = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> headerLineIndex = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> sampleID = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> snpIndex = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> genoAB1 = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> genoAB2 = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> genoForward1 = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> genoForward2 = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> x = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> y = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> theta = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> r = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> xRaw = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> yRaw = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> baf = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> lrr = new HashMap<Integer, ArrayList<String>>();
		HashMap<Integer, ArrayList<String>> gc = new HashMap<Integer, ArrayList<String>>();
		for (java.util.Map.Entry<String, SourceFileHeaderData> entry : headers.entrySet()) {
			SourceFileHeaderData headerData = entry.getValue();
			if (headerData.numFiles == -1) {
				headerData.numFiles = cnt;
			} else if (headerData.numFiles != cnt) {
				return "Number of Files listed in Source File {"	+ entry.getKey()
								+ "} does not equal the number of headers needing validation.  Please check source directory and extension and try again.";
			}
			ArrayList<String> files = totSnpsSet.get(headerData.totalSnps);
			if (files == null) {
				files = new ArrayList<String>();
				totSnpsSet.put(headerData.totalSnps, files);
			}
			files.add(entry.getKey());
			files = headerLineIndex.get(headerData.columnHeaderLineIndex);
			if (files == null) {
				files = new ArrayList<String>();
				headerLineIndex.put(headerData.columnHeaderLineIndex, files);
			}
			files.add(entry.getKey());
			files = sampleID.get(headerData.colSampleIdent);
			if (files == null) {
				files = new ArrayList<String>();
				sampleID.put(headerData.colSampleIdent, files);
			}
			files.add(entry.getKey());
			files = snpIndex.get(headerData.colSnpIdent);
			if (files == null) {
				files = new ArrayList<String>();
				snpIndex.put(headerData.colSnpIdent, files);
			}
			files.add(entry.getKey());
			files = genoAB1.get(headerData.colGenoAB1);
			if (files == null) {
				files = new ArrayList<String>();
				genoAB1.put(headerData.colGenoAB1, files);
			}
			files.add(entry.getKey());
			files = genoAB2.get(headerData.colGenoAB2);
			if (files == null) {
				files = new ArrayList<String>();
				genoAB2.put(headerData.colGenoAB2, files);
			}
			files.add(entry.getKey());
			files = genoForward1.get(headerData.colGeno1);
			if (files == null) {
				files = new ArrayList<String>();
				genoForward1.put(headerData.colGeno1, files);
			}
			files.add(entry.getKey());
			files = genoForward2.get(headerData.colGeno2);
			if (files == null) {
				files = new ArrayList<String>();
				genoForward2.put(headerData.colGeno2, files);
			}
			files.add(entry.getKey());
			files = x.get(headerData.colX);
			if (files == null) {
				files = new ArrayList<String>();
				x.put(headerData.colX, files);
			}
			files.add(entry.getKey());
			files = y.get(headerData.colY);
			if (files == null) {
				files = new ArrayList<String>();
				y.put(headerData.colY, files);
			}
			files.add(entry.getKey());
			files = theta.get(headerData.colTheta);
			if (files == null) {
				files = new ArrayList<String>();
				theta.put(headerData.colTheta, files);
			}
			files.add(entry.getKey());
			files = r.get(headerData.colR);
			if (files == null) {
				files = new ArrayList<String>();
				r.put(headerData.colR, files);
			}
			files.add(entry.getKey());
			files = xRaw.get(headerData.colXRaw);
			if (files == null) {
				files = new ArrayList<String>();
				xRaw.put(headerData.colXRaw, files);
			}
			files.add(entry.getKey());
			files = yRaw.get(headerData.colYRaw);
			if (files == null) {
				files = new ArrayList<String>();
				yRaw.put(headerData.colYRaw, files);
			}
			files.add(entry.getKey());
			files = baf.get(headerData.colBAF);
			if (files == null) {
				files = new ArrayList<String>();
				baf.put(headerData.colBAF, files);
			}
			files.add(entry.getKey());
			files = lrr.get(headerData.colLRR);
			if (files == null) {
				files = new ArrayList<String>();
				lrr.put(headerData.colLRR, files);
			}
			files.add(entry.getKey());
			files = gc.get(headerData.colGC);
			if (files == null) {
				files = new ArrayList<String>();
				gc.put(headerData.colGC, files);
			}
			files.add(entry.getKey());
		}
		int numErrors = 0;
		ArrayList<String> errorMsgs = new ArrayList<String>();
		
		String error;
		error = checkErrors(totSnpsSet, "Total SNPs");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(headerLineIndex, "Index of header line");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(sampleID, "Sample ID column index");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(snpIndex, "SNP column index");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(genoAB1, "AB Genotype column 1");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(genoAB2, "AB Genotype column 2");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(genoForward1, "Forward Genotype column 1");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(genoForward2, "Forward Genotype column 2");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(x, "X column");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(y, "Y column");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(theta, "Theta column");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(r, "R column");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(xRaw, "X Raw column");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(yRaw, "Y Raw column");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(baf, "B Allele Freq column");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(lrr, "Log R Ratio column");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		error = checkErrors(gc, "GC column");
		if (!"".equals(error)) {
			numErrors++;
			errorMsgs.add(error);
		}
		
		if (numErrors > 0 || !errorMsgs.isEmpty()) {
			if (!errorMsgs.isEmpty()) {
				for (String s : errorMsgs) {
					log.reportError(s);
				}
			}
			return "Found "	+ numErrors + " data or column index mismatches among source files.  Please check log for more details";
		}
		
		return null;
	}

	private static String checkErrors(HashMap<Integer, ArrayList<String>> valueMapping, String string) {
		if (valueMapping.size() == 1) {
			return "";
		}
		StringBuilder sb = new StringBuilder("Mismatch in ")
														.append(string)
														.append(" between ")
														.append(valueMapping.size())
														.append(" sets of values: {");
		
		ArrayList<Integer> values = new ArrayList<Integer>(valueMapping.keySet());
		Collections.sort(values);
		Collections.reverse(values);
		
		for (int i = 0, cnt = values.size(), cnt1 = cnt - 1; i < cnt; i++) {
			sb.append(values.get(i));
			if (i < cnt1) {
				sb.append(", ");
			}
		}
		
		sb.append("} with {");
		for (int i = 0, cnt = values.size(), cnt1 = cnt - 1; i < cnt; i++) {
			sb.append(valueMapping.get(values.get(i)).size());
			if (i < cnt1) {
				sb.append(", ");
			}
		}
		sb.append("} file(s) for each value, respectively.");
		if (valueMapping.size() == 2) {
			sb.append(" The file(s) in the second set are: ");
			ArrayList<String> files = valueMapping.get(values.get(values.size() - 1));
			for (int i = 0; i < files.size(); i++) {
				sb.append(files.get(i));
				if (i < files.size() - 1) {
					sb.append(", ");
				}
			}
			sb.append(".");
		}
		return sb.toString();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "D:/data/ny_registry/gedi_gwas/00src/";// null;
		String ext = ".csv.gz";
		String log = null;

		String usage = "\n"	+ "cnv.filesys.FinalReportHeaderData requires 2 argument\n"
										+ "   (1) Directory of FinalReport files (i.e. dir=" + dir + " (default))\n"
										+ "   (2) Extension of FinalReport files (i.e. ext=" + ext + " (default))\n";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("ext=")) {
				ext = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				log = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0 || dir == null) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			validate(dir, ext, true, log == null ? new Logger() : new Logger(log), null);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public String getHeaderString() {
		return headerString;
	}

	public void setHeaderString(String headerString) {
		this.headerString = headerString;
	}

	private void setSourceFileDelimiter(String delim) {
		delimiter = delim;
	}

	public String getSourceFileDelimiter() {
		return delimiter;
	}
}
