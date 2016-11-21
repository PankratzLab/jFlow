package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.filesys.SerialHash;
import org.genvisis.parse.GenParser;
import org.genvisis.stats.Maths;
import org.genvisis.stats.Maths.OPERATOR;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class MarkerQC {
	public static final String DEFAULT_FILENAME = "thresholds.properties";

	public static final String DIR_KEY = "dir";
	public static final String FILE_KEY = "file";
	public static final String FILE_DEFAULT = "markerQC.xln";
	public static final String MARKERS_KEY = "markers";
	public static final String MARKERS_DEFAULT = "freq.frq";

	public static final String[] FINAL_HEADER = {	"SNP", "CHR", "MAF", "F_MISS", "P", "HETERO p-value",
																								"miss.hap min p-value", "P_MISS"};
	public static final String[] HWE_HEADER = {"CHR", "SNP", "TEST", "A1", "A2", "GENO", "O(HET)",
																							"E(HET)", "P"};
	// public static final String[] THRESHOLDS = {"snp", "chr", "maf", "f_miss", "hwe", "hetero",
	// "minmishap", "p_miss"};


	private static final String[] FRQ_HEADER = {"CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"};
	private static final String[] REF_FREQ_HEADER = {"SNP", "P"};
	private static final String[] LMISS_HEADER = {"CHR", "SNP", "N_MISS", "N_GENO", "F_MISS"};
	private static final String[] MISHAP_HEADER = {	"SNP", "HAPLOTYPE", "F_0", "F_1", "M_H1", "M_H2",
																									"CHISQ", "P", "FLANKING"};
	private static final String[] MISSTEST_HEADER = {"CHR", "SNP", "F_MISS_A", "F_MISS_U", "P"};
	private static final String[] GENDER_ASSOC_HEADER = {	"CHR", "SNP", "BP", "A1", "F_A", "F_U", "A2",
																												"CHISQ", "P", "OR"};

	private enum METRIC {
		// FILE("file", "markerQC.xln"),
		// MARKERS("markers", "freq.frq", "1", "header"),
		// CHR("chr", "freq.frq", lessThan(1), "0:1"),
		FRQ("maf", "freq.frq", FRQ_HEADER, lessThan(0.01), "1", "4=maf"),
		REF_FREQ_MISMATCH("ref_freq", "refFreqMismatch.dat", REF_FREQ_HEADER, lessThan(0.0001), "0",
											"1=ref_freq"),
		CALLRATE("callrate", "missing.lmiss", LMISS_HEADER, lessThan(0.98), "1", "$#1-4=callrate"),
		HWE("hwe", "hardy.hwe", HWE_HEADER, lessThan(0.00001), "!2!ALL", "!2!AFF", "1", "8==hwe"),
		// allow for the possibility of a quantitative trait, which lists the test as ALL(QT);
		// otherwise it should only accept the UNAFF test
		MISHAP_HETERO("mishap_hetero", "mishap.missing.hap", MISHAP_HEADER, lessThan(0.0001),
									"!1=HETERO", "0", "7=mishap_hetero"),
		MISHAP_MIN("mishap_min", "mishap.missing.hap", MISHAP_HEADER, lessThan(0.0001), "0", "$MIN7"),
		P_MISS("p_miss", "test.missing.missing", MISSTEST_HEADER, lessThan(0.0001), "1", "4=p_miss"),
		P_GENDER(	"p_gender", "gender.assoc", GENDER_ASSOC_HEADER, lessThan(1E-7), "!0<24", "1",
							"8=p_gender"),
		P_GENDER_MISS("p_gender_miss", "gender.missing", MISSTEST_HEADER, lessThan(0.0001), "!0<24",
									"1", "4=p_gender_miss");

		private static final Map<String, METRIC> KEY_MAP;
		static {
			Map<String, METRIC> keyMap = Maps.newHashMap();
			for (METRIC metric : values()) {
				keyMap.put(metric.getKey(), metric);
			}
			KEY_MAP = Collections.unmodifiableMap(keyMap);
		}

		private String key;
		private String defaultFile;
		private String[] fileHeader;
		private String defaultArg;
		private List<String> lookupArgs;

		/**
		 * @param key
		 * @param defaultFile
		 * @param fileHeader
		 * @param defaultArg
		 */
		private METRIC(	String key, String defaultFile, String[] fileHeader, String defaultArg,
										String... lookupArgs) {
			this.key = key;
			this.defaultFile = defaultFile;
			this.fileHeader = fileHeader;
			this.defaultArg = defaultArg;
			this.lookupArgs = ImmutableList.copyOf(lookupArgs);
		}

		private static String lessThan(Number threshold) {
			return OPERATOR.LESS_THAN.getSymbol() + threshold.toString();
		}

		public String getKey() {
			return key;
		}

		public String getDefaultFile() {
			return defaultFile;
		}

		public String[] getFileHeader() {
			return fileHeader;
		}

		public String getDefaultArg() {
			return defaultArg;
		}

		public String getLookupString(String dir, String filename) {
			return generateLookupString(dir, filename, lookupArgs);
		}

		public String getExample() {
			return generateExampleCode(key, defaultFile, defaultArg);
		}

		public static METRIC getMetric(String key) {
			return KEY_MAP.get(key);
		}

	}

	private static String generateExampleCode(String key, String file, String... args) {
		List<String> exampleCode =  Lists.newArrayList(key + "=" + file);
		for (String arg : args) {
			exampleCode.add(arg);
		}
		return Joiner.on(',').join(exampleCode);
	}
	
	private static String generateLookupString(String dir, String filename, Iterable<String> args) {
		Collection<String> lookupTokens = Lists.newArrayList();
		lookupTokens.add("\"" + dir + filename + "\"");
		for (String arg : args) {
			lookupTokens.add(arg);
		}
		return Joiner.on(' ').join(lookupTokens);
	}

	public static int parse(String dir, String[][] params, Logger log, boolean kill) {
		log.report(ext.getTime());
		long time = new Date().getTime();
		try {
			List<String> fileLookups = Lists.newArrayList();
			List<String[]> headers = Lists.newArrayList();
			for (int i = 3; i < params.length; i++) {
				String key = params[i][0];
				String file = params[i][1];
				METRIC metric = METRIC.getMetric(key);
				if (metric != null) {
					switch (metric) {
						case MISHAP_MIN:
							try {
								prepareMishapMin(dir, file);;
							} catch (FileNotFoundException fnfe) {
								log.reportError("Error: file \"" + file + "\" not found in current directory");
								log.reportException(fnfe);
								if (kill) {
									System.exit(1);
								} else {
									return 1;
								}
							} catch (IOException ioe) {
								log.reportError("Error reading file \"" + file + "\"");
								log.reportException(ioe);
								if (kill) {
									System.exit(2);
								} else {
									return 2;
								}
							}
							break;
						case P_MISS:
							if (!Files.exists(dir + file) || new File(dir + file).length() == 0) {
								continue;
							}
							break;
						default:
							break;
					}
					fileLookups.add(metric.getLookupString(dir, file));
					headers.add(metric.fileHeader);
				} else {
					log.report("Using user-defined file: " + file);
					fileLookups.add(generateLookupString(	dir, file,
																								ImmutableList.of((params[i].length > 3	? params[i][3]
																																												: "0 1")
																																	+ "=" + key))); // allows for
																																									// additional
																																									// user-defined
																																									// files
					headers.add(null);

				}
			}
			String[] markerNames = HashVec.loadFileToStringArray(dir	+ params[1][1],
																														params[1].length > 3 && params[1][3].equals("header"),
																														new int[] {Integer.parseInt(params[1][2])},
																														false);
			log.report("Found " + markerNames.length + " markers to parse in " + params[1][1]);
			Files.writeIterable(fileLookups, dir + "whatGoesIn.out");
			// Files.combine(markerNames, Array.toStringArray(v), Matrix.toStringArrays(headers),
			// "Marker", ".", dir+params[0][1], log, true, true, false);
			Files.combineWithLessMemory(markerNames, Array.toStringArray(fileLookups),
																	Matrix.toStringArrays(headers), "Marker", ".", dir + params[0][1],
																	log, true, true, false, false);
			log.report("Finished in " + ext.getTimeElapsed(time));
		} catch (Exception e) {
			log.reportException(e);
			if (kill) {
				System.exit(2);
			} else {
				return 2;
			}
		}
		return 0;
	}

	private static void prepareMishapMin(String dir, String file)	throws FileNotFoundException,
																																IOException {
		BufferedReader reader = new BufferedReader(new FileReader(dir + file));
		ext.checkHeader(reader.readLine().trim().split("[\\s]+"), METRIC.MISHAP_MIN.getFileHeader(),
										true);
		String prev = "";
		List<String> mishaps = Lists.newArrayList();
		Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
		boolean done = false;
		while (!done) {
			String[] line;
			if (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
			} else {
				done = true;
				line = Array.stringArray(8, "");
			}
			if (!line[0].equals(prev)) {
				double min = 2;
				for (int j = 0; j < mishaps.size(); j++) {
					double current = Double.parseDouble(mishaps.get(j));
					if (current < min) {
						min = current;
					}
				}
				hash.put(prev, new String[] {mishaps.isEmpty() ? "." : Double.toString(min)});
				mishaps.clear();
			}
			if (line.length >= 7 && !line[7].equals("NA")) {
				mishaps.add(line[7]);
			}
			prev = line[0];
		}
		hash.put("!colNames", new String[] {METRIC.MISHAP_MIN.getKey()});
		reader.close();
		SerialHash.createSerializedStringArrayHash(dir
																									+ GenParser.parseSerialFilename(new String[] {file,
																																															"0", "$MIN7"}),
																								hash);
	}

	public static int testThresholds(String dir, String[][] params, Logger log, boolean kill) {
		BufferedReader reader;
		PrintWriter[] writers;
		PrintWriter writer;
		String[] line, ops;
		double[] thresholds;
		int[][] counts;
		int[] indices;
		int count, fail, maxSize;
		String reason, simpleReason;

		try {
			if (Files.exists(dir + params[0][1], false)) {
				log.report("Using existing file: " + dir + params[0][1]);
			} else {
				log.report("Generating new source file: " + dir + params[0][1]);
				parse(dir, params, log, kill);
			}

			ops = new String[params.length - 3];
			thresholds = new double[params.length - 3];
			for (int i = 0; i < params.length - 3; i++) {
				if (!params[i + 3][0].equals("dir")) {
					for (String element : Maths.OPERATORS) {
						// try {
						if (params[i + 3][2].startsWith(element)) {
							ops[i] = element;
						}
						// } catch (Exception e) {
						// System.err.println("Error: "+Array.toStr(params[i+3]));
						// }
					}
					if (ops[i] == null) {
						log.reportError("Error - invalid operator for "	+ params[i + 3][0] + " ('"
														+ params[i + 3][2] + "')");
						System.exit(1);
					}
					try {
						thresholds[i] = Double.parseDouble(params[i + 3][2].substring(ops[i].length()));
					} catch (NumberFormatException nfe) {
						log.reportError("Error - threshold for "	+ params[i + 3][0] + " ('"
														+ params[i + 3][2].substring(1) + "') is not a valid number");
						if (kill) {
							System.exit(1);
						} else {
							return 1;
						}
					}
				}
			}

			count = 0;
			indices = null;
			counts = null;
			try {
				reader = new BufferedReader(new FileReader(dir + params[0][1]));
				writers = new PrintWriter[4];
				writers[0] = new PrintWriter(new FileWriter(dir + params[2][1] + "_drops.dat"));
				writers[1] = new PrintWriter(new FileWriter(dir + params[2][1] + "_singles.out"));
				writers[2] = new PrintWriter(new FileWriter(dir + params[2][1] + "_annotated.xln"));
				writers[3] = new PrintWriter(new FileWriter(dir + params[2][1] + "_allAnnotations.out"));
				line = reader.readLine().trim().split("\\t");
				indices = ext.indexFactors(	Array.subArray(Matrix.extractColumn(params, 0), 3), line, false,
																		log, true, false);
				counts = new int[3][indices.length]; // all, primary, only
				while (reader.ready()) {
					line = reader.readLine().trim().split("\\t");
					fail = 0;
					reason = "";
					simpleReason = "";
					for (int lap = 0; lap < 2; lap++) {
						for (int i = 0; i < indices.length; i++) {
							if (indices[i] >= 0	&& !ext.isMissingValue(line[indices[i]])
									&& Maths.op(Double.parseDouble(line[indices[i]]), thresholds[i], ops[i])) {
								if (lap == 0) {
									counts[0][i]++;
									if (fail == 0) {
										counts[1][i]++;
									}
									reason += (reason.length() == 0 ? "" : "; ")	+ params[i + 3][0] + "="
														+ line[indices[i]];
									simpleReason += (simpleReason.length() == 0 ? "" : ";")	+ params[i + 3][0]
																	+ params[i + 3][2];
									fail++;
								} else if (fail == 1) {
									counts[2][i]++;
									writers[1].println(line[0] + "\t" + reason);
								}
							}
						}
					}
					if (fail > 0) {
						writers[0].println(line[0]);
						writers[3].println(line[0] + "\t" + simpleReason);
					}
					writers[2].println(line[0] + "\t" + reason);
					count++;
				}
				reader.close();
				for (PrintWriter writer2 : writers) {
					writer2.close();
				}
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + params[0][1] + "\" not found in current directory");
				log.reportException(fnfe);
				if (kill) {
					System.exit(1);
				} else {
					return 1;
				}
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + params[0][1] + "\"");
				log.reportException(ioe);
				if (kill) {
					System.exit(2);
				} else {
					return 2;
				}
			}

			maxSize = 0;
			for (int i = 0; i < indices.length; i++) {
				maxSize = Math.max(maxSize, (params[i + 3][0] + params[i + 3][2] + ":").length());
			}

			writer = new PrintWriter(new FileWriter(dir + params[2][1] + ".out"));
			writer.println("Number of markers failed: "	+ Array.sum(counts[1]) + " of " + count + " ("
											+ ext.formDeci((double) Array.sum(counts[1]) / (double) count * 100, 2)
											+ "%)");
			writer.println();

			writer.println("Number of total markers failed for each filter: ");
			for (int i = 0; i < indices.length; i++) {
				writer.println(ext.formStr(params[i + 3][0] + params[i + 3][2] + ":", maxSize + 2, true)
													+ "\t" + counts[0][i] + " ("
												+ ext.formDeci((double) counts[0][i] / (double) count * 100, 2) + "%)");
			}
			writer.println();

			writer.println("Number of additional markers failed for each consecutive filter: ");
			for (int i = 0; i < indices.length; i++) {
				writer.println(ext.formStr(params[i + 3][0] + params[i + 3][2] + ":", maxSize + 2, true)
													+ "\t" + counts[1][i] + " ("
												+ ext.formDeci((double) counts[1][i] / (double) count * 100, 2) + "%)");
			}
			writer.println();

			writer.println("Number of markers failed for only a single reason: "	+ Array.sum(counts[2])
											+ " (" + ext.formDeci((double) Array.sum(counts[2]) / (double) count * 100, 2)
											+ "%)");
			for (int i = 0; i < indices.length; i++) {
				writer.println(ext.formStr(params[i + 3][0] + params[i + 3][2] + ":", maxSize + 2, true)
													+ "\t" + counts[2][i] + " ("
												+ ext.formDeci((double) counts[2][i] / (double) count * 100, 2) + "%)");
			}
			writer.println();

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + params[2][1] + ".out");
			log.reportException(e);
			if (kill) {
				System.exit(1);
			} else {
				return 1;
			}
		}
		return 0;
	}

	public static int parseParameters(String filename, Logger log, boolean kill) {
		String[] line, record;
		Vector<String[]> v;
		String key;
		String[] file, markers;
		Vector<String> paramV;
		String dir;

		dir = "";
		List<String> sampleCode = Lists.newArrayList();
		sampleCode.add(generateExampleCode(FILE_KEY, FILE_DEFAULT));
		sampleCode.add(generateExampleCode(MARKERS_KEY, MARKERS_DEFAULT));
		METRIC[] metrics = METRIC.values();
		for (METRIC metric : metrics) {
			sampleCode.add(metric.getExample());
		}
		paramV = Files.parseControlFile(filename, "miss", sampleCode.toArray(new String[sampleCode.size()]), log);
		if (paramV != null) {
			file = null;
			markers = null;
			v = new Vector<String[]>();
			for (int i = 0; i < paramV.size(); i++) {
				line = paramV.elementAt(i).trim().split("=");
				key = line[0];
				line = line[1].trim().split(",");
				record = new String[line.length + 1];
				record[0] = key;
				for (int j = 0; j < line.length; j++) {
					// FIXME This is unclear and potentially unpredictable
					// Is this an attempt to replace ":" with " " unless the ":" is a part of a filepath/url?
					record[j + 1] = ext.replaceAllWith(line[j], ":", "QWERTY");
					record[j + 1] = ext.replaceAllWith(record[j + 1], "QWERTY/", ":/");
					record[j + 1] = ext.replaceAllWith(record[j + 1], "QWERTY", " ");
				}
				if (key.equals(FILE_KEY)) {
					file = record;
					file[0] = dir + file[0];
				} else if (key.equals(DIR_KEY)) {
					dir = line[0];
					System.out.println("Running out of: " + dir);
				} else if (key.equals(MARKERS_KEY)) {
					if (record.length < 3) {
						markers = new String[] {key, line[0], "0"};
					} else {
						markers = record;
					}
				} else {
					v.add(record);
				}
			}
			if (file == null) {
				log.reportError("Error - no file listed (i.e. file=markerQC.xln)");
				if (kill) {
					System.exit(1);
				} else {
					return 1;
				}
			} else {
				v.insertElementAt(file, 0);
			}

			if (markers == null) {
				log.reportError("Error - need to define the marker list/order, even if it's just markers=freq.frq,1");
				if (kill) {
					System.exit(1);
				} else {
					return 1;
				}
			} else {
				v.insertElementAt(markers, 1);
			}

			v.insertElementAt(new String[] {"root", ext.rootOf(filename)}, 2);
			int testCode = testThresholds(dir, Matrix.toStringArrays(v), log, kill);
			if (testCode != 0) {
				if (kill) {
					System.exit(testCode);
				} else {
					return testCode;
				}
			}
		}
		return 0;
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = DEFAULT_FILENAME;
		Logger log;

		String usage = "\n"	+ "gwas.MarkerQC requires 0-1 arguments\n"
										+ "   (0) properties file (i.e. file=" + filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		filename = "miss.crf";

		log = new Logger(ext.rootOf(filename, false) + ".log");
		try {
			parseParameters(filename, log, true);
		} catch (Exception e) {
			log.reportException(e);
		}
	}
}
