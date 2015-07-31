package affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Date;

import common.Array;
import common.Files;
import common.ext;
import cnv.filesys.Project;
import cnv.manage.ParseAffySNP6;

public class AffySNP6Tables {
	public static final String[][] SNP_HEADER_OPTIONS = { { "SNP Name", "rsID", "Probe Set ID", "probeset_id" } };
	String[] indFiles;
	private String outputDirectory;
	private String callFile;
	private String confFile;
	private String sigFile;

	public AffySNP6Tables(String outputDirectory, String sigFile) {
		this.outputDirectory = outputDirectory;
		this.sigFile = sigFile;
	}

	public AffySNP6Tables(String outputDirectory, String callFile, String confFile, String sigFile) {
		this.outputDirectory = outputDirectory;
		this.callFile = callFile;
		this.confFile = confFile;
		this.sigFile = sigFile;
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
			System.err.println("unknown call code: " + callCode);
		}

		return call;
	}

	public double power2(String signal) {
		return (Math.pow(2, Double.parseDouble(signal)));
	}

	public void parseSNPLine(String[] calls, String[] confs, String[] sigA, String[] sigB) {
		for (int j = 1; j < calls.length; j++) {
			indFiles[j - 1] += calls[0] + "\t" + parseCall(calls[j]) + "\t" + confs[j] + "\t" + Double.toString(power2(sigA[j])) + "\t" + Double.toString(power2(sigB[j])) + "\t" + parseCall(calls[j]) + "\n";
		}
	}

	// setting CN probeset calls to NC (-1), confidence to 0, and sigB to 0;
	public void parseCNLine(String[] sigA) {
		for (int j = 1; j < sigA.length; j++) {
			indFiles[j - 1] += sigA[0] + "\tNC\t0\t" + Double.toString(power2(sigA[j])) + "\t" + Double.toString(power2(sigA[j])) + "\tNC\n";
		}
	}

	public static PrintWriter getWriter(String filename, boolean append) {
		PrintWriter writer = null;
		try {

			writer = new PrintWriter(new FileWriter(filename, append));
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" could not be written to (it's probably open)");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		return writer;
	}

	public void mkdir(String outputDirectory) {
		File file = new File(outputDirectory);
		if (!file.exists()) {
			if (file.mkdirs()) {
				System.out.println("Directory " + outputDirectory + " is created!");
			} else {
				System.out.println("Failed to create  " + outputDirectory + "!");
			}
		}
	}

	public void printIt(String[] header, int chunkCount) {
		PrintWriter writer = null;
		boolean append = true;
		long time = new Date().getTime();
		for (int j = 1; j < header.length; j++) {
			if (chunkCount == 0) {
				mkdir(outputDirectory);
				append = false;
				writer = getWriter(outputDirectory + header[j] + ".IND.txt", append);
				writer.print("Probe Set ID\tCall Codes\tConfidence\tSignal A\tSignal B\tForced Call Codes\n");
				writer.print(indFiles[j - 1]);
				writer.close();
				indFiles[j - 1] = "";
			} else {
				writer = getWriter(outputDirectory + header[j] + ".IND.txt", append);
				writer.print(indFiles[j - 1]);
				writer.close();
				indFiles[j - 1] = "";
			}
		}
		System.out.println("Printing took " + ext.getTimeElapsed(time));
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
			int numFiles = header.length - 1;
			indFiles = new String[numFiles];
			Arrays.fill(indFiles, "");
			while (sigReader.ready()) {
				do {
					sigALine = sigReader.readLine().trim().split(delimiter, -1);
				} while (!sigALine[0].substring(0, 3).equals("CN_"));
				if (sigALine[0].substring(0, 3).equals("CN_")) {
					parseCNLine(sigALine);
					lineCount++;
					if (lineCount >= numLinesBuffer) {
						System.out.println("Parsed " + chunkCount * numLinesBuffer + " lines");
						System.out.println(ext.getTime() + " Free memory: " + ((float) 100 * Runtime.getRuntime().freeMemory() / Runtime.getRuntime().totalMemory()) + "%");
						printIt(header, chunkCount);
						lineCount = 0;
						chunkCount++;
					}
				} else {
					System.err.println("This Should Not Happen");
				}
			}
			if (lineCount < numLinesBuffer) {
				printIt(header, chunkCount);
				int numLinesTotal = chunkCount * numLinesBuffer + lineCount + 1;
				System.out.println("Parsed a total of " + numLinesTotal + " CN probesets");
			}
			sigReader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: one of the input tables was not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + "\"");
		}
	}

	public void parseSNPTables(int numLinesBuffer) {
		BufferedReader callReader, confReader, sigReader;
		String callHeader, confHeader, sigHeader;
		String[] callLine, confLine, sigALine, sigBLine;
		int chunkCount = 0;
		int lineCount = 0;
		String delimiter = "\t";
		int totalLines = 0;

		try {
			callReader = Files.getAppropriateReader(callFile);
			confReader = Files.getAppropriateReader(confFile);
			sigReader = Files.getAppropriateReader(sigFile);
			callHeader = getHeader(callReader, delimiter, callFile);
			confHeader = getHeader(confReader, delimiter, confFile);
			sigHeader = getHeader(sigReader, delimiter, sigFile);

			if (callHeader.equals(confHeader) && callHeader.equals(sigHeader)) {

				String[] header = callHeader.split(delimiter, -1);
				int numFiles = header.length - 1;
				indFiles = new String[numFiles];
				Arrays.fill(indFiles, "");
				callLine = callReader.readLine().trim().split(delimiter, -1);
				confLine = confReader.readLine().trim().split(delimiter, -1);
				sigALine = allignFiles(sigReader, callLine, confLine, delimiter);
				sigBLine = sigReader.readLine().trim().split(delimiter, -1);

				while (callReader.ready()) {
					totalLines++;
					if (callLine[0].equals(confLine[0]) && sigALine[0].equals(callLine[0] + "-A") && sigBLine[0].equals(callLine[0] + "-B")) {
						parseSNPLine(callLine, confLine, sigALine, sigBLine);
						lineCount++;
						if (lineCount >= numLinesBuffer) {
							printIt(header, chunkCount);
							lineCount = 0;
							chunkCount++;
							System.out.println("Parsed " + chunkCount * numLinesBuffer + " marker lines");
						}
					} else if (!callLine[0].equals(confLine[0]) || !sigALine[0].equals(callLine[0] + "-A") || !sigBLine[0].equals(callLine[0] + "-B")) {
						System.err.println("Error: probeset identifier mismatch between calls/confidence/signal files ");
						System.exit(1);
					} else if (!sigReader.ready()) {
						System.err.println("Error: probeset identifier discordance between calls/confidence/signal files");
						return;
					} else {
						System.err.println("This Should Not Happen");
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
					System.out.println("Parsed a total of " + numLinesTotal + " SNP probesets");
				}
				callReader.close();
				confReader.close();
				sigReader.close();
			} else {
				System.out.println("table headers do not match ");
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: one of the input tables was not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + "\"");
		}
	}

	private void parseLastBatch(String[] callLine, String[] confLine, String[] sigALine, String[] sigBLine, int chunkCount, String[] header) {
		parseSNPLine(callLine, confLine, sigALine, sigBLine); // Last Line
		printIt(header, chunkCount);
	}

	private String[] allignFiles(BufferedReader sigReader, String[] callLine, String[] confLine, String delimiter) throws IOException {
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
			} while (reader.ready() && (ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1));
			header = Array.toStr(line);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + tableName + "\"");
			return "Error reading file \"" + tableName + "\"";
		}
		return header;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj;
		String filename = null;
		;
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

		for (int i = 0; i < args.length; i++) {

			if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
			if (args[i].startsWith("affyResultsDir=")) {
				affyResultsDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-SNP_")) {
				SNP = true;
				numArgs--;
			} else if (args[i].startsWith("-CN_")) {
				CN = true;
				numArgs--;
			} else if (args[i].startsWith("numLinesBuffer=")) {
				numLinesBuffer = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-merge")) {
				merge = true;
				numArgs--;
			} else if (args[i].startsWith("-create")) {
				create = true;
				numArgs--;
			} else if (args[i].startsWith("-combo")) {
				comboMergeCreate = true;
				numArgs--;
			} else if (args[i].startsWith("calls=")) {
				calls = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("conf=")) {
				conf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("sig=")) {
				sig = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				output = args[i].split("=")[1];
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
				AffySNP6Tables AS6T = new AffySNP6Tables(output, calls, conf, sig);
				AS6T.parseSNPTables(numLinesBuffer);
			}
			if (CN) {
				AffySNP6Tables AS6T = new AffySNP6Tables(output, sig);
				AS6T.parseCNTable(numLinesBuffer);
			}
			if (merge && !comboMergeCreate) {
				MergeChp.combineChpFiles(affyResultsDir, numThreads, commonSubFolderPattern, inputFileNameExt, output);

			}
			if (create && !comboMergeCreate) {
				ParseAffySNP6.createFiles(proj, numThreads);
			} else if (comboMergeCreate) {
				MergeChp.combineChpFiles(affyResultsDir, numThreads, commonSubFolderPattern, inputFileNameExt, output);
				ParseAffySNP6.createFiles(proj, numThreads);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
