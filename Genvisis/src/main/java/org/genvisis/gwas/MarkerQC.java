package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Date;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.filesys.SerialHash;
import org.genvisis.parse.GenParser;
import org.genvisis.stats.Maths;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class MarkerQC {
	public static final String DEFAULT_FILENAME = "thresholds.properties";
	public static final String COMMAND = "miss";
	public static final String[] FRQ_HEADER = {"CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"};
	public static final String[] LMISS_HEADER = {"CHR", "SNP", "N_MISS", "N_GENO", "F_MISS"};
	public static final String[] HWE_HEADER = {"CHR", "SNP", "TEST", "A1", "A2", "GENO", "O(HET)",
																						 "E(HET)", "P"};
	public static final String[] MISSHAP_HEADER = {"SNP", "HAPLOTYPE", "F_0", "F_1", "M_H1", "M_H2",
																								 "CHISQ", "P", "FLANKING"};
	public static final String[] MISSTEST_HEADER = {"CHR", "SNP", "F_MISS_A", "F_MISS_U", "P"};
	public static final String[] GENDER_ASSOC_HEADER = {"CHR", "SNP", "BP", "A1", "F_A", "F_U", "A2",
																											"CHISQ", "P", "OR"};
	public static final String[] FINAL_HEADER = {"SNP", "CHR", "MAF", "F_MISS", "P", "HETERO p-value",
																							 "miss.hap min p-value", "P_MISS"};
	// public static final String[] THRESHOLDS = {"snp", "chr", "maf", "f_miss", "hwe", "hetero",
	// "minmishap", "p_miss"};

	public static final Map<QC_METRIC, String> DEFAULT_METRIC_THRESHOLDS;
	public static final Map<QC_METRIC, String> DEFAULT_METRIC_FILENAMES;

	public static final String DEFAULT_ILLUMINA_CALLRATE_THRESHOLD = "<0.98";
	public static final String DEFAULT_AFFY_CALLRATE_THRESHOLD = "<0.95";

	private static final List<String> DEFAULT_CRF_START_LINES = ImmutableList.of("file=markerQC.xln",
																																							 "markers=freq.frq,1,header");

	static {
		Map<QC_METRIC, String> defaultMetricThresholds = Maps.newHashMap();
		defaultMetricThresholds.put(QC_METRIC.CHR, "<1");
		defaultMetricThresholds.put(QC_METRIC.MAF, "<0.01");
		defaultMetricThresholds.put(QC_METRIC.CALLRATE, DEFAULT_ILLUMINA_CALLRATE_THRESHOLD);
		defaultMetricThresholds.put(QC_METRIC.HWE, "<0.00001");
		defaultMetricThresholds.put(QC_METRIC.MISHAP_HETERO, "<0.0001");
		defaultMetricThresholds.put(QC_METRIC.MISHAP_MIN, "<0.0001");
		defaultMetricThresholds.put(QC_METRIC.P_MISS, "<0.0001");
		defaultMetricThresholds.put(QC_METRIC.P_GENDER, "<1E-7");
		defaultMetricThresholds.put(QC_METRIC.P_GENDER_MISS, "<0.0001");
		DEFAULT_METRIC_THRESHOLDS = Collections.unmodifiableMap(defaultMetricThresholds);


		Map<QC_METRIC, String> defaultMetricFilenames = Maps.newEnumMap(QC_METRIC.class);
		defaultMetricFilenames.put(QC_METRIC.CHR, "freq.frq");
		defaultMetricFilenames.put(QC_METRIC.MAF, "freq.frq");
		defaultMetricFilenames.put(QC_METRIC.CALLRATE, "missing.lmiss");
		defaultMetricFilenames.put(QC_METRIC.HWE, "hardy.hwe");
		defaultMetricFilenames.put(QC_METRIC.MISHAP_HETERO, "mishap.missing.hap");
		defaultMetricFilenames.put(QC_METRIC.MISHAP_MIN, "mishap.missing.hap");
		defaultMetricFilenames.put(QC_METRIC.P_MISS, "test.missing.missing");
		defaultMetricFilenames.put(QC_METRIC.P_GENDER, "gender.assoc");
		defaultMetricFilenames.put(QC_METRIC.P_GENDER_MISS, "gender.missing");
		DEFAULT_METRIC_FILENAMES = Collections.unmodifiableMap(defaultMetricFilenames);
	}

	public static enum QC_METRIC {
		CHR("Chromosome Number"),
		MAF("Minor Allele Frequency"),
		CALLRATE("Marker Callrate"),
		HWE("p-value for Hardy-Weinberg Equilibrium"),
		MISHAP_HETERO("p-value for heterozygous flanking mishaps"),
		MISHAP_MIN("p-value for all mishaps"),
		P_MISS("p-value for case-control association with missingess"),
		P_GENDER("p-value for sex association"),
		P_GENDER_MISS("p-value for sex association with missigness");

		private final String key;
		private final String description;

		QC_METRIC(String description) {
			this(description, null);
		}

		QC_METRIC(String description, String key) {
			this.description = description;
			this.key = key != null ? key : this.name().toLowerCase();
		}

		public String getKey() {
			return key;
		}

		public String getDescription() {
			return description;
		}

		public String getCLIDescription() {
			String userDescription = getUserDescription();
			return userDescription.substring(0, 1).toLowerCase() + userDescription.substring(1);
		}

		public String getUserDescription() {
			StringBuilder descriptionBuilder = new StringBuilder();
			descriptionBuilder.append("Threshold to reject ").append(getDescription())
												.append(", using: ");
			descriptionBuilder.append(Joiner.on(", ").join(Maths.OPERATORS));
			descriptionBuilder.append(" (<0 to not filter)");
			return descriptionBuilder.toString();
		}

	}


	/**
	 * 
	 * @param metricThresholds Map from {@link QC_METRIC} to desired threshold, all of the keys will
	 *        be included
	 * @param metricFilenames Map from QC_METRIC to filenames, keys not in {@code metricThresholds}
	 *        will be ignored
	 * @return lines to use in CRF
	 */
	private static List<String> generateCRFLines(Map<QC_METRIC, String> metricThresholds,
																							 Map<QC_METRIC, String> metricFilenames) {
		List<String> crfLines = Lists.newArrayList(DEFAULT_CRF_START_LINES);
		for (Map.Entry<QC_METRIC, String> metricThreshold : metricThresholds.entrySet()) {
			QC_METRIC metric = metricThreshold.getKey();
			String threshold = metricThreshold.getValue();
			crfLines.add(metric.getKey() + "=" + metricFilenames.get(metric) + "," + threshold);
		}
		return crfLines;
	}

	/**
	 * 
	 * {@link #generateCRF(String, String, Map, Map)} using the default thresholds and filenames
	 */
	public static void generateCRF(String targetDir, String outfile) {
		generateCRF(targetDir, outfile, DEFAULT_METRIC_THRESHOLDS);
	}

	/**
	 * {@link #generateCRF(String, String, Map, Map)} using the default filenames
	 */
	public static void generateCRF(String targetDir, String outfile,
																 Map<QC_METRIC, String> metricThresholds) {
		generateCRF(targetDir, outfile, metricThresholds, DEFAULT_METRIC_FILENAMES);
	}


	/**
	 * Generates an appropriate {@link MarkerQC} CRF
	 * 
	 * @param targetDir directory with files to parse
	 * @param outfile the filename of the CRF to generate
	 * @param metricThresholds thresholds to use for each desired {@link QC_METRIC}
	 * @param metricFilenames filenames to use for each desired {@code QC_METRIC}, any not in
	 *        {@code metricThresholds} are ignored
	 */
	public static void generateCRF(String targetDir, String outfile,
																 Map<QC_METRIC, String> metricThresholds,
																 Map<QC_METRIC, String> metricFilenames) {
		List<String> lines = generateCRFLines(metricThresholds, metricFilenames);
		lines.add(0, COMMAND);
		lines.add(1, "dir=" + targetDir);
		Files.writeIterable(lines, outfile);
	}

	public static int parse(String dir, String[][] params, Logger log, boolean kill) {
		BufferedReader reader;
		String[] line;
		Hashtable<String, String[]> hash;
		Vector<String> v;
		Vector<String[]> headers;
		Vector<String> mishaps;
		String[] markerNames;
		String prev, trav, min;
		long time;
		boolean done;

		log.report(ext.getTime());
		time = new Date().getTime();
		try {
			v = new Vector<String>();
			headers = new Vector<String[]>();
			for (int i = 3; i < params.length; i++) {
				if (params[i][0].equals("maf")) {
					v.add("\"" + dir + params[i][1] + "\" 1 4=" + params[i][0]);
					headers.add(FRQ_HEADER);
				} else if (params[i][0].equals("callrate")) {
					v.add("\"" + dir + params[i][1] + "\" 1 $#1-4=" + params[i][0]);
					headers.add(LMISS_HEADER);
				} else if (params[i][0].equals("hwe")) {
					v.add("\"" + dir + params[i][1] + "\" !2!ALL !2!AFF 1 8=" + params[i][0]); // allow for
																																										 // the
																																										 // possibility
																																										 // of a
																																										 // quantitative
																																										 // trait,
																																										 // which
																																										 // lists the
																																										 // test as
																																										 // ALL(QT);
																																										 // otherwise
																																										 // it
																																										 // should only
																																										 // accept the
																																										 // UNAFF test
					headers.add(HWE_HEADER);
				} else if (params[i][0].equals("mishap_hetero")) {
					v.add("\"" + dir + params[i][1] + "\" !1=HETERO 0 7=" + params[i][0]);
					headers.add(MISSHAP_HEADER);
				} else if (params[i][0].equals("mishap_min")) {
					try {
						reader = new BufferedReader(new FileReader(dir + params[i][1]));
						ext.checkHeader(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE),
														MISSHAP_HEADER, true);
						prev = "";
						mishaps = new Vector<String>();
						hash = new Hashtable<String, String[]>();
						done = false;
						while (!done) {
							if (reader.ready()) {
								line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
							} else {
								done = true;
								line = ArrayUtils.stringArray(8, "");
							}
							if (!line[0].equals(prev)) {
								min = "2";
								for (int j = 0; j < mishaps.size(); j++) {
									trav = mishaps.elementAt(j);
									if (Double.parseDouble(trav) < Double.parseDouble(min)) {
										min = trav;
									}
								}
								hash.put(prev, new String[] {mishaps.size() > 0 ? min : "."});
								mishaps.clear();
							}
							if (line.length >= 7 && !line[7].equals("NA")) {
								mishaps.add(line[7]);
							}
							prev = line[0];
						}
						hash.put("!colNames", new String[] {params[i][0]});
						reader.close();
						SerialHash.createSerializedStringArrayHash(dir
																											 + GenParser.parseSerialFilename(new String[] {params[i][1],
																																																		 "0",
																																																		 "$MIN7"}),
																											 hash);
					} catch (FileNotFoundException fnfe) {
						log.reportError("Error: file \"" + params[i][1] + "\" not found in current directory");
						log.reportException(fnfe);
						if (kill) {
							System.exit(1);
						} else {
							return 1;
						}
					} catch (IOException ioe) {
						log.reportError("Error reading file \"" + params[i][1] + "\"");
						log.reportException(ioe);
						if (kill) {
							System.exit(2);
						} else {
							return 2;
						}
					}
					v.add("\"" + dir + params[i][1] + "\" 0 $MIN7");
					headers.add(MISSHAP_HEADER);
				} else if (params[i][0].equals("p_miss")) {
					if (Files.exists(dir + params[i][1]) && new File(dir + params[i][1]).length() > 0) {
						v.add("\"" + dir + params[i][1] + "\" 1 4=" + params[i][0]);
						headers.add(MISSTEST_HEADER);
					}
				} else if (params[i][0].equals("p_gender")) {
					v.add("\"" + dir + params[i][1] + "\" !0<24 1 8=" + params[i][0]);
					headers.add(GENDER_ASSOC_HEADER);
				} else if (params[i][0].equals("p_gender_miss")) {
					v.add("\"" + dir + params[i][1] + "\" !0<24 1 4=" + params[i][0]);
					headers.add(MISSTEST_HEADER);
				} else {
					log.report("Using user-defined file: " + params[i][1]);
					v.add("\"" + dir + params[i][1] + "\" " + (params[i].length > 3 ? params[i][3] : "0 1")
								+ "=" + params[i][0]); // allows for additional user-defined files
					headers.add(null);
				}
			}
			markerNames = HashVec.loadFileToStringArray(dir + params[1][1],
																									params[1].length > 3 && params[1][3].equals("header"),
																									new int[] {Integer.parseInt(params[1][2])},
																									false);
			log.report("Found " + markerNames.length + " markers to parse in " + params[1][1]);
			Files.writeArray(ArrayUtils.toStringArray(v), dir + "whatGoesIn.out");
			// Files.combine(markerNames, Array.toStringArray(v), Matrix.toStringArrays(headers),
			// "Marker", ".", dir+params[0][1], log, true, true, false);
			Files.combineWithLessMemory(markerNames, ArrayUtils.toStringArray(v),
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
			if (Files.exists(dir + params[0][1])) {
				log.report("Using existing file: " + dir + params[0][1]);
			} else {
				log.report("Generating new source file: " + dir + params[0][1]);
				parse(dir, params, log, kill);
			}

			ops = new String[params.length - 3];
			thresholds = new double[params.length - 3];
			for (int i = 0; i < params.length - 3; i++) {
				if (!params[i + 3][0].equals("dir")) {
					Maths.OPERATOR op = findOperator(params[i + 3][2]);
					if (op == null) {
						log.reportError("Error - invalid operator for " + params[i + 3][0] + " ('"
														+ params[i + 3][2] + "')");
						System.exit(1);
					} else {
						ops[i] = op.getSymbol();
					}
					try {
						thresholds[i] = Double.parseDouble(params[i + 3][2].substring(ops[i].length()));
					} catch (NumberFormatException nfe) {
						log.reportError("Error - threshold for " + params[i + 3][0] + " ('"
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
				writers[0] = Files.openAppropriateWriter(dir + params[2][1] + "_drops.dat");
				writers[1] = Files.openAppropriateWriter(dir + params[2][1] + "_singles.out");
				writers[2] = Files.openAppropriateWriter(dir + params[2][1] + "_annotated.xln");
				writers[3] = Files.openAppropriateWriter(dir + params[2][1] + "_allAnnotations.out");
				line = reader.readLine().trim().split("\\t");
				indices = ext.indexFactors(ArrayUtils.subArray(Matrix.extractColumn(params, 0), 3), line,
																	 false, log, true, false);
				counts = new int[3][indices.length]; // all, primary, only
				while (reader.ready()) {
					line = reader.readLine().trim().split("\\t");
					fail = 0;
					reason = "";
					simpleReason = "";
					for (int lap = 0; lap < 2; lap++) {
						for (int i = 0; i < indices.length; i++) {
							if (indices[i] >= 0 && !ext.isMissingValue(line[indices[i]])
									&& Maths.op(Double.parseDouble(line[indices[i]]), thresholds[i], ops[i])) {
								if (lap == 0) {
									counts[0][i]++;
									if (fail == 0) {
										counts[1][i]++;
									}
									reason += (reason.length() == 0 ? "" : "; ") + params[i + 3][0] + "="
														+ line[indices[i]];
									simpleReason += (simpleReason.length() == 0 ? "" : ";") + params[i + 3][0]
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

			writer = Files.openAppropriateWriter(dir + params[2][1] + ".out");
			writer.println("Number of markers failed: " + ArrayUtils.sum(counts[1]) + " of " + count
										 + " ("
										 + ext.formDeci((double) ArrayUtils.sum(counts[1]) / (double) count * 100, 2)
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

			writer.println("Number of markers failed for only a single reason: "
										 + ArrayUtils.sum(counts[2]) + " ("
										 + ext.formDeci((double) ArrayUtils.sum(counts[2]) / (double) count * 100, 2)
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

	public static Maths.OPERATOR findOperator(String thresholdString) {
		for (Maths.OPERATOR op : Maths.OPERATOR.values()) {
			if (thresholdString.startsWith(op.getSymbol()))
				return op;
		}
		return null;
	}

	public static int parseParameters(String filename, Logger log, boolean kill) {
		String[] line, record;
		Vector<String[]> v;
		String key;
		String[] file, markers;
		List<String> paramV;
		String dir;

		dir = "";
		paramV = Files.parseControlFile(filename, COMMAND,
																		ArrayUtils.toStringArray(generateCRFLines(DEFAULT_METRIC_THRESHOLDS,
																																							DEFAULT_METRIC_FILENAMES)),
																		log);
		if (paramV != null) {
			file = null;
			markers = null;
			v = new Vector<String[]>();
			for (int i = 0; i < paramV.size(); i++) {
				line = paramV.get(i).trim().split("=");
				key = line[0];
				line = line[1].trim().split(",");
				record = new String[line.length + 1];
				record[0] = key;
				for (int j = 0; j < line.length; j++) {
					record[j + 1] = ext.replaceAllWith(line[j], ":", "QWERTY");
					record[j + 1] = ext.replaceAllWith(record[j + 1], "QWERTY/", ":/");
					record[j + 1] = ext.replaceAllWith(record[j + 1], "QWERTY", " ");
				}
				if (key.equals("file")) {
					file = record;
					file[0] = dir + file[0];
				} else if (key.equals("dir")) {
					dir = line[0];
					System.out.println("Running out of: " + dir);
				} else if (key.equals("markers")) {
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

		String usage = "\n" + "gwas.MarkerQC requires 0-1 arguments\n"
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
