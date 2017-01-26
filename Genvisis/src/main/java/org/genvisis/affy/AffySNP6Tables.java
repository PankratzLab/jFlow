package org.genvisis.affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.SourceFileParser;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class AffySNP6Tables {
	public static final String[][] SNP_HEADER_OPTIONS = {{"SNP Name", "rsID", "Probe Set ID",
																												"probeset_id"}};
	StringBuilder[] indFiles;
	private final String outputDirectory;
	private String callFile;
	private String confFile;
	private final String sigFile;
	private final int maxWriters;
	private PrintWriter[] writers;
	private final Logger log;

	public AffySNP6Tables(String outputDirectory, String sigFile, int maxWriters, Logger log) {
		this.outputDirectory = outputDirectory;
		this.sigFile = sigFile;
		this.maxWriters = maxWriters;
		writers = null;
		this.log = log;
	}

	public AffySNP6Tables(String outputDirectory, String callFile, String confFile, String sigFile,
												int maxWriters, Logger log) {
		this.outputDirectory = outputDirectory;
		this.callFile = callFile;
		this.confFile = confFile;
		this.sigFile = sigFile;
		this.maxWriters = maxWriters;
		writers = null;
		this.log = log;
	}

	public String parseCall(String callCode) {
		String call = "NC";
		if (callCode.equals("0")) {
			call = "AA";
		} else if (callCode.equals("1")) {
			call = "AB";
		} else if (callCode.equals("2")) {
			call = "BB";
		} else if (callCode.equals("-1")) {
			call = "NC";
		} else {
			log.reportError("unknown call code: " + callCode);
		}

		return call;
	}

	public double power2(String signal) {
		return (Math.pow(2, Double.parseDouble(signal)));
	}

	public void parseSNPLine(String[] calls, String[] confs, String[] sigA, String[] sigB) {
		for (int j = 1; j < calls.length; j++) {
			indFiles[j - 1].append(calls[0]	+ "\t" + parseCall(calls[j]) + "\t" + confs[j] + "\t"
															+ Double.toString(power2(sigA[j])) + "\t"
															+ Double.toString(power2(sigB[j])) + "\t" + parseCall(calls[j])
															+ "\n");
		}
	}

	// setting CN probeset calls to NC (-1), confidence to 0, and sigB to 0;
	public void parseCNLine(String[] sigA) {
		for (int j = 1; j < sigA.length; j++) {
			indFiles[j - 1].append(sigA[0]	+ "\tNC\t0\t" + Double.toString(power2(sigA[j])) + "\t"
															+ Double.toString(power2(sigA[j])) + "\tNC\n");
		}
	}

	public static PrintWriter getWriter(String filename, boolean append, Logger log) {
		PrintWriter writer = null;
		try {

			writer = new PrintWriter(new FileWriter(filename, append));
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \""	+ filename
													+ "\" could not be written to (it's probably open)");
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		return writer;
	}

	public boolean allExist(String[] header, String outputDirectory) {
		for (int i = 1; i < header.length; i++) {
			if (!Files.exists(outputDirectory + ext.rootOf(header[i]) + ".txt")) {
				return false;
			}
		}
		return true;
	}

	@SuppressWarnings("resource")
	public void printIt(String[] header, int chunkCount) {
		boolean append = true;
		long time = new Date().getTime();
		for (int j = 1; j < header.length; j++) {
			PrintWriter writer = writers == null ? null : writers[j - 1];
			if (chunkCount == 0) {
				new File(outputDirectory).mkdirs();
				append = false;
				writer = getWriter(outputDirectory + ext.rootOf(header[j]) + ".txt", append, log);
				writer.print("Probe Set ID\tCall Codes\tConfidence\tSignal A\tSignal B\tForced Call Codes\n");
				writer.print(indFiles[j - 1]);
				if (header.length < maxWriters) {
					if (writers == null) {
						log.reportTimeInfo("Initializing static writers since number of samples is less than "
																+ maxWriters);
						writers = new PrintWriter[header.length - 1];
					}
					writers[j - 1] = writer;
				} else {
					writer.close();
				}
				indFiles[j - 1] = new StringBuilder();
			} else {
				if (writer == null) {
					writer = getWriter(outputDirectory + ext.rootOf(header[j]) + ".txt", append, log);
				}
				writer.print(indFiles[j - 1]);
				if (writers == null) {
					writer.close();
				}
				indFiles[j - 1] = new StringBuilder();
			}
		}
		if (writers == null) {
			log.reportTimeInfo("Printing took " + ext.getTimeElapsed(time));
		}
	}

	public void parseCNTable(int numLinesBuffer) {
		BufferedReader sigReader;
		String[] sigALine;
		int chunkCount = 0;
		int lineCount = 0;
		String delimiter = "\t";
		String sigHeader;
		try {

			sigReader = Files.getAppropriateReader(sigFile);
			sigHeader = getHeader(sigReader, delimiter, sigFile);
			String[] header = sigHeader.split(delimiter, -1);
			if (!allExist(header, outputDirectory)) {

				int numFiles = header.length - 1;
				indFiles = new StringBuilder[numFiles];
				for (int i = 0; i < indFiles.length; i++) {
					indFiles[i] = new StringBuilder();
				}
				while (sigReader.ready()) {
					do {
						sigALine = sigReader.readLine().trim().split(delimiter, -1);
					} while (!sigALine[0].substring(0, 3).equals("CN_"));
					if (sigALine[0].substring(0, 3).equals("CN_")) {
						parseCNLine(sigALine);
						lineCount++;
						if (lineCount >= numLinesBuffer || writers != null) {
							if (writers == null) {
								log.reportTimeInfo("Parsed " + chunkCount * numLinesBuffer + " lines");
							} else if (chunkCount % 10000 == 0) {
								log.reportTimeInfo("Static writers: parsed " + chunkCount + " lines");
							}
							// log.reportTimeInfo(ext.getTime() + " Free memory: " + ((float) 100 *
							// Runtime.getRuntime().freeMemory() / Runtime.getRuntime().totalMemory()) + "%");
							printIt(header, chunkCount);
							lineCount = 0;
							chunkCount++;
						}
					} else {
						log.reportError("This Should Not Happen");
					}
				}
				if (lineCount < numLinesBuffer) {
					printIt(header, chunkCount);
					int numLinesTotal = chunkCount * numLinesBuffer + lineCount + 1;
					if (writers == null) {
						log.reportTimeInfo("Parsed a total of " + numLinesTotal + " CN probesets");
					}
				}
				sigReader.close();
			} else {
				log.reportTimeInfo("All files exist in "	+ outputDirectory
														+ ", skipping parsing. Please delete these files if you would like to re-parse");

			}
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: one of the input tables was not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + "\"");
		}
		if (writers != null) {
			for (PrintWriter writer : writers) {
				writer.close();
			}
		}
	}

	public void parseSNPTables(int numLinesBuffer) {
		BufferedReader callReader, confReader, sigReader;
		String callHeader, confHeader, sigHeader;
		String[] callLine, confLine, sigALine, sigBLine;
		int chunkCount = 0;
		int lineCount = 0;
		String delimiter = "\t";
		// int totalLines = 0;

		try {
			callReader = Files.getAppropriateReader(callFile);
			confReader = Files.getAppropriateReader(confFile);
			sigReader = Files.getAppropriateReader(sigFile);
			callHeader = getHeader(callReader, delimiter, callFile);
			confHeader = getHeader(confReader, delimiter, confFile);
			sigHeader = getHeader(sigReader, delimiter, sigFile);

			if (callHeader.equals(confHeader) && callHeader.equals(sigHeader)) {

				String[] header = callHeader.split(delimiter, -1);
				if (!allExist(header, outputDirectory)) {

					int numFiles = header.length - 1;
					indFiles = new StringBuilder[numFiles];
					for (int i = 0; i < indFiles.length; i++) {
						indFiles[i] = new StringBuilder();
					}
					callLine = callReader.readLine().trim().split(delimiter, -1);
					confLine = confReader.readLine().trim().split(delimiter, -1);
					sigALine = allignFiles(sigReader, callLine, confLine, delimiter);
					sigBLine = sigReader.readLine().trim().split(delimiter, -1);

					while (callReader.ready()) {
						// totalLines++;
						if (callLine[0].equals(confLine[0])	&& sigALine[0].equals(callLine[0] + "-A")
								&& sigBLine[0].equals(callLine[0] + "-B")) {
							parseSNPLine(callLine, confLine, sigALine, sigBLine);
							lineCount++;
							if (lineCount >= numLinesBuffer || writers != null) {
								printIt(header, chunkCount);
								lineCount = 0;
								chunkCount++;
								if (writers == null) {
									log.reportTimeInfo("Parsed " + chunkCount * numLinesBuffer + " lines");
								} else if (chunkCount % 10000 == 0) {
									log.reportTimeInfo("Static writers: parsed " + chunkCount + " lines");
								}
							}
						} else if (!callLine[0].equals(confLine[0])	|| !sigALine[0].equals(callLine[0] + "-A")
												|| !sigBLine[0].equals(callLine[0] + "-B")) {
							log.reportError("Error: probeset identifier mismatch between calls/confidence/signal files ");
							System.exit(1);
						} else if (!sigReader.ready()) {
							log.reportError("Error: probeset identifier discordance between calls/confidence/signal files");
							return;
						} else {
							log.reportError("This Should Not Happen");
							System.exit(1);
						}
						callLine = callReader.readLine().trim().split(delimiter, -1);
						confLine = confReader.readLine().trim().split(delimiter, -1);
						sigALine = sigReader.readLine().trim().split(delimiter, -1);
						sigBLine = sigReader.readLine().trim().split(delimiter, -1);
					}
					if (lineCount < numLinesBuffer) {
						parseLastBatch(callLine, confLine, sigALine, sigBLine, chunkCount, header);
						int numLinesTotal = chunkCount * numLinesBuffer + lineCount + 1;
						if (writers == null) {
							log.reportTimeInfo("Parsed a total of " + numLinesTotal + " SNP probesets");
						}
					}
					callReader.close();
					confReader.close();
					sigReader.close();
				} else {
					log.reportTimeInfo("All files exist in "	+ outputDirectory
															+ ", skipping parsing. Please delete these files if you would like to re-parse");
				}
			} else {
				log.reportTimeInfo("table headers do not match ");
			}
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: one of the input tables was not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + "\"");
		}
		if (writers != null) {
			for (PrintWriter writer : writers) {
				writer.close();
			}
		}
	}

	private void parseLastBatch(String[] callLine, String[] confLine, String[] sigALine,
															String[] sigBLine, int chunkCount, String[] header) {
		parseSNPLine(callLine, confLine, sigALine, sigBLine); // Last Line
		printIt(header, chunkCount);
	}

	private String[] allignFiles(	BufferedReader sigReader, String[] callLine, String[] confLine,
																String delimiter) throws IOException {
		String[] sigALine;

		do {
			sigALine = sigReader.readLine().trim().split(delimiter, -1);
		} while (callLine[0].equals(confLine[0]) && !sigALine[0].equals(callLine[0] + "-A"));
		return sigALine;
	}

	public String getHeader(BufferedReader reader, String delimiter, String tableName) {
		String[] line;
		String header;

		try {
			do {
				line = reader.readLine().trim().split(delimiter, -1);
			} while (reader.ready() && (ext.indexFactors(	SNP_HEADER_OPTIONS, line, false, true, false,
																										false)[0] == -1));
			header = ArrayUtils.toStr(line);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + tableName + "\"");
			return "Error reading file \"" + tableName + "\"";
		}
		return header;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj;
		String filename = null;;
		boolean SNP = false;
		boolean CN = false;
		boolean merge = false;
		boolean create = false;
		boolean comboMergeCreate = false;
		int numLinesBuffer = 100;
		int numThreads = 2;
		String calls = "C:/data/AFFYtable/00src/SNP_/birdseed-v2.calls.txt";
		String conf = "C:/data/AFFYtable/00src/SNP_/birdseed-v2.confidences.txt";
		String sig = "C:/data/AFFYtable/00src/SNP_/quant-norm.pm-only.med-polish.expr.summary.txt";
		String output = "C:/data/AFFYtable/00src/SNP_/";
		String inputFileNameExt = ".txt";
		String affyResultsDir = "";
		String commonSubFolderPattern = "";
		int maxwriters = 5000;

		for (String arg : args) {

			if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			}
			if (arg.startsWith("affyResultsDir=")) {
				affyResultsDir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("threads=")) {
				numThreads = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("-SNP_")) {
				SNP = true;
				numArgs--;
			} else if (arg.startsWith("-CN_")) {
				CN = true;
				numArgs--;
			} else if (arg.startsWith("numLinesBuffer=")) {
				numLinesBuffer = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("-merge")) {
				merge = true;
				numArgs--;
			} else if (arg.startsWith("-create")) {
				create = true;
				numArgs--;
			} else if (arg.startsWith("-combo")) {
				comboMergeCreate = true;
				numArgs--;
			} else if (arg.startsWith("calls=")) {
				calls = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("conf=")) {
				conf = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("sig=")) {
				sig = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("out=")) {
				output = arg.split("=")[1];
				numArgs--;
			}

		}
		if (numArgs != 0) {
			// TODO usage
			System.err.println("TODO usuage");
			System.exit(1);
		}
		proj = new Project();
		if (filename != null) {
			proj = new Project(filename, false);
		}
		try {
			if (SNP) {
				AffySNP6Tables AS6T =
														new AffySNP6Tables(output, calls, conf, sig, maxwriters, new Logger());
				AS6T.parseSNPTables(numLinesBuffer);
			}
			if (CN) {
				AffySNP6Tables AS6T = new AffySNP6Tables(output, sig, maxwriters, new Logger());
				AS6T.parseCNTable(numLinesBuffer);
			}
			if (merge && !comboMergeCreate) {
				MergeChp.combineChpFiles(	affyResultsDir, numThreads, commonSubFolderPattern,
																	inputFileNameExt, output, new Logger());

			}
			if (create && !comboMergeCreate) {
				SourceFileParser.createFiles(proj, numThreads);
			} else if (comboMergeCreate) {
				MergeChp.combineChpFiles(	affyResultsDir, numThreads, commonSubFolderPattern,
																	inputFileNameExt, output, new Logger());
				SourceFileParser.createFiles(proj, numThreads);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
