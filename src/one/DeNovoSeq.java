package one;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TimeZone;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.Positions;
import common.Sort;
import common.ext;

public class DeNovoSeq {
	public static final String[] SAMPLE_SUFFIX = new String[] {"C", "D", "M"};
	public static final char[] BASES = new char[] {'A', 'T', 'G', 'C'};
	public static final char[] BASES_WITH_N = new char[] {'A', 'T', 'G', 'C', 'N'};
	public static final String[] BASES_WITH_N_INDEL = new String[] {"A", "T", "G", "C", "N", "Ins", "Del"};
	public static final byte INDEX_OF_N = 4;
	public static final byte INDEX_OF_INS = 5;
	public static final byte INDEX_OF_DEL = 6;
	public static final String[] ALLELE_COUNTS_ARRAY_STRUCT = new String[] {"A", "T", "G", "C", "N", "Ins", "Del", "totalReads", "numAllelesStrictCount", "numAllelesLooseCount"};
	public static final byte INDEX_OF_TOTAL_READS = 7;
	public static final byte INDEX_OF_NUM_ALLELES_STRICT = 8;
	public static final byte INDEX_OF_NUM_ALLELES_LOOSE = 9;
	public static final byte MARKERDATA_NUMSAMPLES_START = 0;
	public static final byte MARKERDATA_NUMSAMPLES_LEN = 4;
	public static final int DEFAULT_PHRED_SCORE_FOR_DELETION = 30;
	public static final double THRESHOLD_PHRED_SCORE_FOR_INS_DEL = .10;
	public static final int MAX_ALLELE_COUNT_TREATED_AS_ZERO = 2;
	public static final int MIN_ALLELE_COUNT_FOR_DENOVO_MUTATION = 5;
	public static final double MIN_ALLELE_FREQ_FOR_DENOVO_MUTATION = .20;
	public static final int MIN_READ_DEPTH = 10;
	public static final int WINDOW_SIZE_FOR_NEARBY_INDEL = 60;
	public static final int WINDOW_SIZE_FOR_NEARBY_VARIANCE = 30;
	public static final double DISCOUNT_FOR_NEARBY_INDEL = .80;
	public static final double DISCOUNT_FOR_ON_INDEL_SITE = .25;
	public static final double DISCOUNT_FOR_NEARBY_VARIANCE = .95;
	public static final double DISCOUNT_FOR_ON_VARIANCE_SITE = .70;
	public static final double DISCOUNT_FOR_AT_VARIANCE = .60;
	public static final double DISCOUNT_FOR_N = .50;
	public static final double DISCOUNT_FOR_3ALLELES_LOOSE = .25;
	public static final double DISCOUNT_FOR_3ALLELES_STRICT = .50;
	public static final double DISCOUNT_FOR_NEARBY_DENOVOMUTATION = .90;
	public static final int REGION_LEN_AT_A_TIME = 100000;

	public static void generateScriptForSamtools(String fileNameOfDeNovoPointMutationCandidateList, String bamFileDir, String scriptDir) {
		BufferedReader candidateList;
		String[] line;
		String commands;
		Vector<String[]> iterationsVec;
		String[][] iterations;

		commands = "samtools view [%0]C.bam  [%1]:[%2]-[%3] > " + "[%0]C_[%1]_[%3].txt\n"
				+ "samtools view [%0]D.bam  [%1]:[%2]-[%3] > " + "[%0]D_[%1]_[%3].txt\n"
				+ "samtools view [%0]M.bam  [%1]:[%2]-[%3] > " + "[%0]M_[%1]_[%3].txt";
		try {
			candidateList = new BufferedReader(new FileReader(fileNameOfDeNovoPointMutationCandidateList));
			candidateList.readLine();
			iterationsVec = new Vector<String[]>(); 
			while(candidateList.ready()) {
				line = candidateList.readLine().split("\t");
				iterationsVec.add(new String[] {line[0].substring(0, line[0].length()-1), line[1], line[2], line[3]});
			}
			candidateList.close();
			iterations = new String[iterationsVec.size()][4];
			for (int i = 0; i < iterations.length; i++) {
				iterations[i] = iterationsVec.elementAt(i); 
			}
			Files.qsub(scriptDir + "getReads_", bamFileDir, 16, commands, iterations, -1, -1);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void screenDeNovoPointMutation(String fileNameOfDeNovoPointMutationCandidateList, String bamFileDir, int thresholdFor3Alleles) {
		BufferedReader candidateList;
		PrintWriter writer;
		String[] line;
		int[][] alleleCounts;
		String notes;
		int score;

		try {
			candidateList = new BufferedReader(new FileReader(fileNameOfDeNovoPointMutationCandidateList));
			candidateList.readLine();
			writer = new PrintWriter(ext.parseDirectoryOfFile(fileNameOfDeNovoPointMutationCandidateList) + ext.rootOf(fileNameOfDeNovoPointMutationCandidateList) + "_screened.txt");
			writer.println("ID\tchrs\tposition\tlookup\tSarver\tREF\tALT\tMendelianLikelihood\tMendelianPP\tMendelianGT\tSnpCode\tCode\tDeNovoLikelihood\tDeNovoPP\tActualDenovo\tConf\tPankratzScore\tNotes");
			while(candidateList.ready()) {
				line = candidateList.readLine().split("\t");
				//look for the chr location in bcfFile
				alleleCounts = new int[3][];
				for (int i = 0; i < alleleCounts.length; i++) {
					alleleCounts[i] = getAlleleCounts(bamFileDir + line[0].substring(0, line[0].length()-1) + SAMPLE_SUFFIX[i] + "_" + line[1] + "_" + line[3] + ".txt", Integer.parseInt(line[1].split("chr")[1]), Integer.parseInt(line[3]), thresholdFor3Alleles);
				}
				score = getScoreForDeNovoMutation(alleleCounts, line[9].charAt(0), line[8].charAt(0), thresholdFor3Alleles);
				notes = getNotesForDeNovoMutation(alleleCounts, line[9].charAt(0), thresholdFor3Alleles);
				writer.println(line[0] + "\t" + line[1] + "\t" + line[3] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\t" + line[9] + "\t" + line[10] + "\t" + line[11] + "\t" + line[12] + "\t" + line[13] + "\t" + line[14] + "\t" + line[15] + "\t" + line[16] + "\t" + line[17] + "\t" + line[18] + "\t" + score + "\t" + notes);
			}
			candidateList.close();
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static int[] getAlleleCounts(String readsFileFullPath, int chr, int position, int thresholdFor3Alleles) {
		int[] result;
		BufferedReader bamFile;
		String[] line;
//		char[] bases;
		int index, startPositionOfTheRead, lengthOfCurrentSegmentOfCurrentRead, positionOnThisRead;
		String[][] operatorsOperatorIndicesAndSplit;

		result = new int[] {0, 0, 0, 0, 0, 0, 0, 0};	//A, T, G, C, N, Ins, Del, Total, #Alleles
		try {
			bamFile =  new BufferedReader(new FileReader(readsFileFullPath));
			while(bamFile.ready()) {
				positionOnThisRead = position;
				line = bamFile.readLine().split("\t");
				if (line[2].equalsIgnoreCase("chr" + chr) && Integer.parseInt(line[3])<=position && (Integer.parseInt(line[3]) + Math.abs(Integer.parseInt(line[8]))) >=position && Integer.parseInt(line[1]) < 512) {
					startPositionOfTheRead = Integer.parseInt(line[3]);
					index = 0;
					operatorsOperatorIndicesAndSplit = ext.getOperatorsOperatorIndicesAndSplit(line[5], "DIMNPSH");
					for (int i = 0; i < operatorsOperatorIndicesAndSplit[0].length; i++) {
						lengthOfCurrentSegmentOfCurrentRead = Integer.parseInt(operatorsOperatorIndicesAndSplit[2][i]);
						if ((startPositionOfTheRead + index + lengthOfCurrentSegmentOfCurrentRead) > positionOnThisRead) {
							if (operatorsOperatorIndicesAndSplit[0][i].equals("I") || operatorsOperatorIndicesAndSplit[0][i].equals("S") || operatorsOperatorIndicesAndSplit[0][i].equals("H") || operatorsOperatorIndicesAndSplit[0][i].equals("N")) {
								positionOnThisRead += lengthOfCurrentSegmentOfCurrentRead;
								index += lengthOfCurrentSegmentOfCurrentRead;
							} else if (operatorsOperatorIndicesAndSplit[0][i].equals("D") || operatorsOperatorIndicesAndSplit[0][i].equals("P")) {
								result[6] ++;
								break;
							} else {
								index += (positionOnThisRead - startPositionOfTheRead - index);
//								index = (positionOnThisRead - startPositionOfTheRead);
								index = ext.indexOfChar(line[9].charAt(index), BASES_WITH_N);
								if (index >= 0) {
									result[index] ++;
								} else {
									System.out.println("Error - unrecognized base " + line[9].charAt(index) + " in the following file and read:\n" + readsFileFullPath + "\nRead ID: " + line[0]);
								}
								break;
							}
						} else {
							if (operatorsOperatorIndicesAndSplit[0][i].equals("I") || operatorsOperatorIndicesAndSplit[0][i].equals("S") || operatorsOperatorIndicesAndSplit[0][i].equals("H") || operatorsOperatorIndicesAndSplit[0][i].equals("N")) {
								positionOnThisRead += lengthOfCurrentSegmentOfCurrentRead;
								index += lengthOfCurrentSegmentOfCurrentRead;
							} else if (operatorsOperatorIndicesAndSplit[0][i].equals("D") || operatorsOperatorIndicesAndSplit[0][i].equals("P")) {
								startPositionOfTheRead += lengthOfCurrentSegmentOfCurrentRead;
							} else {
								index += lengthOfCurrentSegmentOfCurrentRead;
							}
						}
					}
				}
			}
			bamFile.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		for (int i = 0; i < 4; i++) {
			result[7] += result[i];
			if (result[i] > thresholdFor3Alleles) {
				result[8] ++;
			}
		}

		return result;
	}

//	public static boolean[] isThreeAlleles(int[][] alleleCounts, int threshold) {
//		byte count;
//		boolean[] result;
//
//		result = new boolean[alleleCounts.length];
//		for (int i = 0; i < result.length; i++) {
//			count = 0;
//			for (int j = 0; j < alleleCounts[i].length; j++) {
//				if (alleleCounts[i][j] > threshold) {
//					count ++;
//				}
//			}
//			if(count > 2) {
//				result[i] = true;
//			}
//		}
//		return result;
//	}

//	public static boolean isDeNovoMutation(int[][] alleleCounts, int threshold) {
//		for (int j = 0; j < alleleCounts[0].length; j++) {
//			if (alleleCounts[0][j] > threshold && alleleCounts[1][j] <= threshold && alleleCounts[2][j] <= threshold) {
//				return true;
//			}
//		}
//		return false;
//	}

	/**
	 * Annotate the markers with allele frequencies, allele counts, and deletion and insertion counts.
	 * @param alleleCounts
	 * @param alternativeAllele
	 * @param thresholdFor3Alleles
	 * @return
	 */
	public static String getNotesForDeNovoMutation(int[][] alleleCounts, char alternativeAllele, int thresholdFor3Alleles) {
		byte index;
		index = (byte) ext.indexOfChar(alternativeAllele, BASES);

		return	((alleleCounts[1][index] > thresholdFor3Alleles || alleleCounts[2][index] > thresholdFor3Alleles)? ext.formDeci(100 * alleleCounts[0][index] / (double) alleleCounts[0][5], 1) + "%," + alleleCounts[1][index] + "," + alleleCounts[2][index] + "; " : "")
			+	"C" + (alleleCounts[0][6]>2? "(" + alleleCounts[0][6]+ " alleles)" : "") + ":" + alleleCounts[0][0] + "," + alleleCounts[0][1] + "," + alleleCounts[0][2] + "," + alleleCounts[0][3] + ",(" + alleleCounts[0][4] + "," + alleleCounts[0][5] + ")"
			+ "; D" + (alleleCounts[1][6]>2? "(" + alleleCounts[1][6] + "alleles)" : "") + ":" + alleleCounts[1][0] + "," + alleleCounts[1][1] + "," + alleleCounts[1][2] + "," + alleleCounts[1][3] + ",(" + alleleCounts[1][4] + "," + alleleCounts[1][5] + ")"
			+ "; M" + (alleleCounts[2][6]>2? "(" + alleleCounts[2][6] + "alleles)" : "") + ":" + alleleCounts[2][0] + "," + alleleCounts[2][1] + "," + alleleCounts[2][2] + "," + alleleCounts[2][3] + ",(" + alleleCounts[2][4] + "," + alleleCounts[2][5] + ")";
	}

	/**
	 * A decimal value ranging 0 through 1 to indicate how likely the marker is a De Novo point mutation, with 1 being very likely and 0 very unlikely 
	 * @param alleleCounts
	 * @param alternativeAllele
	 * @param referencedAllele
	 * @param thresholdFor3Alleles
	 * @return
	 */
	//This is functioning but naive. Has been replaced by getDenovoMarkerCandidateScores(...)
	public static int getScoreForDeNovoMutation(int[][] alleleCounts, char alternativeAllele, char referencedAllele, int thresholdFor3Alleles) {
		int DenovoMutationScore;
		int index;

		index = ext.indexOfChar(alternativeAllele, BASES);
		if (alleleCounts[1][index] > thresholdFor3Alleles || alleleCounts[2][index] > thresholdFor3Alleles || alleleCounts[0][8] > 2) {
			DenovoMutationScore = 0;
		} else {
			DenovoMutationScore = 1;
		}
		return	DenovoMutationScore;
	}

	public static boolean isInsertionDeletionOtherMutationsNearby() {
		return false;
	}

	public static void processGenomeOfAllTriosInDir(String bamDir, String refFastaFilename, String bedFilename, String outputDir, String scriptDir, int numThreads, Logger log) {
		String[] bamFilenames;
		String[][] bamFilenamesByTrios;
		String command;
		String trioId;
		Vector<String> qsubFilesVec;
		PrintWriter writer;

		bamFilenames = Files.list(bamDir, ".bam", false);
		bamFilenamesByTrios = groupNamesByTrios(bamFilenames);
		qsubFilesVec = new Vector<String>(bamFilenamesByTrios.length);
		command = "cd " + outputDir + "\njcp one.DeNovoSeq bed=" + bedFilename + " outdir=" + outputDir + " bamdir=" + bamDir + " reffasta=" + refFastaFilename + " numthreads=" + numThreads + " bamset=";
		for (int i = 0; i < bamFilenamesByTrios.length; i++) {
//			processGenomeOfOneTrio(bamDir, bamFilenames, refFastaFilename, bedFilename, outputDir, numThreads, log);
			trioId = getRootOf(bamFilenamesByTrios[i]);
			Files.qsub(scriptDir + trioId + ".qsub", command + bamFilenamesByTrios[i][0] + "," + bamFilenamesByTrios[i][1] + "," + bamFilenamesByTrios[i][2], 2000, 36, 1);
			qsubFilesVec.add(scriptDir + trioId + ".qsub");
		}

		try {
			writer = new PrintWriter(new FileOutputStream(scriptDir + "master_run"));
			writer.println("cd " + scriptDir);
			for (int i = 0; i < qsubFilesVec.size(); i++) {
				writer.println("qsub " + qsubFilesVec.elementAt(i));
			}
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		log.report("Scripts for all trios in " + bamDir + " are ready at " + scriptDir);
	}

	public static void processGenomeOfOneTrio(String bamDir, String[] bamFilenamesOfTheTrio, String refFastaFilename, String bedFilename, String outputDir, int numThreads, Logger log) {
		BufferedReader reader;
		String[] line;
		String chr, prevChr;
		int start;
		int stop;
		Vector<String> chrs = null;
		Vector<Integer> starts = null;
		Vector<Integer> stops = null;
		int loop;
		PrintWriter writer;
		String trioId;
		String outFileName;
		ExecutorService executor = null;
		Runnable worker = null;
		int processId;
		long timer;
		SimpleDateFormat timeFormat;

		if (log == null) {
			log = new Logger();
		}

		timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
     	timer = new Date().getTime();
		trioId = getRootOf(bamFilenamesOfTheTrio);
		outFileName = outputDir + trioId + "_denovoSeq_" + ext.rootOf(bedFilename) + ".txt";
		processId = 0;
		prevChr = "start";
		try {
			writer = new PrintWriter(outFileName);
			writer.println("id\tchr\tpos\tlookup\tref\talt\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore\t1\t2\t4\t5\t7\t8");
			reader = Files.getAppropriateReader(bedFilename);
			reader.readLine();
			reader.readLine();
			if (numThreads > 1) {
				executor = Executors.newFixedThreadPool(numThreads);
				chrs = new Vector<String>();
				starts = new Vector<Integer>();
				stops = new Vector<Integer>();
			}
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				chr = line[0].substring(3);
				if (!chr.equals(prevChr) && numThreads == 1) {
					log.report(ext.getTime() + "\t" + "Starting chr" + chr + "; identified " + Files.countLines(outFileName, true) + " possible de novo events so far");
					prevChr = chr;
				}
				if (line.length < 2 || line[1] == null || line[1].equals("")) {
					start = 1;
					stop = Positions.CHROMOSOME_LENGTHS_MAX[ext.indexOfStr(chr, Positions.CHR_CODES)];
					log.report("No start position is specified. Using default region: chr" + chr + ":" + start + "-" + stop);
				} else {
					start = Integer.parseInt(line[1]);
					stop = Integer.parseInt(line[2]);
				}
				while (start <= stop) {
					loop = Math.min(start + REGION_LEN_AT_A_TIME, stop);
					if (numThreads > 1) {
						chrs.add(chr);
						starts.add(start);
						stops.add(loop);
					} else {
						processRegion(bamDir, bamFilenamesOfTheTrio, trioId, refFastaFilename, chr, start, loop, writer, null);
					}
					start += (REGION_LEN_AT_A_TIME + 1);
				}
			}
			if (numThreads > 1) {
				for (int i = 0; i < numThreads; i++) {
					worker = new WorkerThread(bamDir, bamFilenamesOfTheTrio, trioId, refFastaFilename, chrs, starts, stops, writer, null, processId);
					executor.execute(worker);
				}
				executor.awaitTermination(7, TimeUnit.DAYS);
				executor.shutdown();
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + bedFilename + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + bedFilename + "\"");
			return;
		} catch (InterruptedException e) {
		}

		log.report("DeNovo mutation result is ready at: " + outFileName + "\nTotal time used " + timeFormat.format(new Date().getTime() - timer));
	}

	public static class WorkerThread implements Runnable {
		private String bamDir;
		private String[] bamFilenames;
		private String trioId;
		private String refFastaFilename;
		private Vector<String> chr;
		private Vector<Integer> start;
		private Vector<Integer> stop;
		private PrintWriter writer;
		private String outAlleleCountsFileName;
		private int threadId;

		public WorkerThread(String bamDir, String[] bamFilenames, String trioId, String refFastaFilename, Vector<String> chr, Vector<Integer> start, Vector<Integer> stop, PrintWriter writer, String outAlleleCountsFileName, int threadId) {
			super();
			this.bamDir = bamDir;
			this.bamFilenames = bamFilenames;
			this.trioId = trioId;
			this.refFastaFilename = refFastaFilename;
			this.chr = chr;
			this.start = start;
			this.stop = stop;
			this.writer = writer;
			this.outAlleleCountsFileName = outAlleleCountsFileName;
			this.threadId = threadId;
		}

		@Override
		public void run() {
			int loop;
			loop = chr.size();
			for (int i = 0; i < loop; i++) {
				processRegion(bamDir, bamFilenames, trioId, refFastaFilename, chr.elementAt(i), start.elementAt(i), stop.elementAt(i), writer, outAlleleCountsFileName);
			}
		}
	}

//	public static byte[][] processRegion(String dir, String bamFilename, byte chr, int start, int stop) {
//		PrintStream stream;
//		similar to piping or redirecting in linux
//		BufferedReader reader = new BufferedReader();
//		CmdLine.run("samtools view "+bamFilename+" chr"+chr+":"+start+"-"+stop, dir, stream);

	public static void processRegion(String bamDir, String[] bamFilenames, String refFastaFilename, String chr, int start, int stop, String outputDir, boolean isToOutputAlleleCounts) {
		PrintWriter writer;
		String trioId;
		String outFileName;
//		PrintWriter writer2;
		String outAlleleCountsFileName;
		long timer;
		SimpleDateFormat timeFormat;

		trioId = getRootOf(bamFilenames);
		outFileName = outputDir + trioId + "_denovoMutations_chr" + chr + "_" + start + "_" + stop + ".txt";
		if (isToOutputAlleleCounts) {
			outAlleleCountsFileName = outputDir + trioId + "_alleleCounts_chr" + chr + "_" + start + "_" + stop + ".txt";
		} else {
			outAlleleCountsFileName = null;
		}
        timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
       	timer = new Date().getTime();
		try {
			writer = new PrintWriter(outFileName);
//			writer.println("id\tchr\tpos\tlookup\tsarver\tref\talt\tmendelianLikelihood\tmendelianPP\tmendelianGT\tsnpCode\tcode\tdeNovoLikelihood\tdeNovoPP\tactualDenovo\tconf\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tchildQuality\tdadQuality\tmomQuality");
			writer.println("id\tchr\tpos\tlookup\tref\talt\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore\t1\t2\t4\t5\t7\t8");
			processRegion(bamDir, bamFilenames, trioId, refFastaFilename, chr, start, stop, writer, outAlleleCountsFileName);
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		timer = new Date().getTime() - timer;
		System.out.println("processRegion result is ready at: " + outFileName + (isToOutputAlleleCounts? "\nAllele counts result is ready at:" + outAlleleCountsFileName : "") + "\nTotal time used " + timeFormat.format(timer));
	}

	public static void processRegion(String bamDir, String[] bamFilenames, String trioId, String refFastaFilename, String chr, int start, int stop, PrintWriter writer, String outAlleleCountsFileName) {
		Process p;
//		ProcessBuilder ps;
		BufferedReader reader;
		BufferedReader error;
		String line;
		Vector<String[]> bamContentVec;
		int numLines;
		int numMarkers;
		int[][][] alleleCounts = null;
		int[][][] phredScores = null;
		int[][][] mappingScores = null;
		String refAlleles;
		byte[] denovoMutationScores;
		String[] denovoMarkerNotes;

        try {
            numMarkers = stop - start + 1;
			alleleCounts = new int[SAMPLE_SUFFIX.length][numMarkers][ALLELE_COUNTS_ARRAY_STRUCT.length];
			phredScores = new int[SAMPLE_SUFFIX.length][numMarkers][BASES_WITH_N_INDEL.length];
			mappingScores = new int[SAMPLE_SUFFIX.length][numMarkers][BASES_WITH_N_INDEL.length];

			for (int i = 0; i < bamFilenames.length; i++) {
				bamContentVec = new Vector<String[]>();

				p = Runtime.getRuntime().exec("samtools view " + bamFilenames[i] + " chr" + chr + ":" + start + "-" + stop, null, new File(bamDir));

//				ps=new ProcessBuilder("samtools", "view", bamFilenames[i], "chr" + chr + ":" + start + "-" + stop);
//		        ps.redirectErrorStream(true);
//		        ps.directory(new File(bamDir));
//		        p = ps.start();  

//		        p.waitFor();
				reader = new BufferedReader(new InputStreamReader(p.getInputStream()));

		        while ((line = reader.readLine()) != null) {
		        	bamContentVec.add(line.split("[\t]"));
		        }
		        p.waitFor();
		        reader.close();
		
				numLines = bamContentVec.size();
		        for (int j = 0; j < numLines; j++) {
		        	getAlleleCountsPhredScoresMappingScores(bamContentVec.elementAt(j), chr, start, stop, 0, alleleCounts[i], phredScores[i], mappingScores[i]);
				};

//		        error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
//		        while ((line = error.readLine()) != null) {
//		            System.out.println(line);
//		        }
//		        p.waitFor();
//				error.close();

				displayErrorStream(p);
			}
			if (outAlleleCountsFileName != null) {
				saveAlleleCountsToFile(outAlleleCountsFileName, alleleCounts, start, true);
			}

			p = Runtime.getRuntime().exec("samtools faidx " + refFastaFilename + " chr" + chr + ":" + start + "-" + stop, null, new File(ext.parseDirectoryOfFile(refFastaFilename)));

//			ps=new ProcessBuilder("samtools", "faidx", refFastaFilename, "chr" + chr + ":" + start + "-" + stop);
//	        ps.redirectErrorStream(true);
//	        ps.directory(new File(ext.parseDirectoryOfFile(refFastaFilename)));
//	        p = ps.start();  

//	        p.waitFor();
	        reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
	        line = reader.readLine();
        	refAlleles = "";
	        while ((line = reader.readLine()) != null) {
	        	refAlleles += line;
	        }
	        p.waitFor();
	        reader.close();

//	        error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
//	        while ((line = error.readLine()) != null) {
//	            System.out.println(line);
//	        }
//			error.close();
	        
	        displayErrorStream(p);

			denovoMutationScores = new byte[numMarkers];
			denovoMarkerNotes = new String[numMarkers];
			getDenovoMutationScores(alleleCounts, phredScores, mappingScores, denovoMutationScores, denovoMarkerNotes);
			for (int i = 0; i < numMarkers; i++) {
				if(denovoMutationScores[i] > 0) {
					exportInfoForPosition(writer, trioId, chr, start + i, refAlleles.charAt(i), denovoMutationScores[i], alleleCounts, phredScores, mappingScores, i, denovoMarkerNotes, "Phred score proportions < " + THRESHOLD_PHRED_SCORE_FOR_INS_DEL);
				}
			}
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        } catch (InterruptedException e) {
			e.printStackTrace();
		}
        
	}

	//	resultAlleleCounts = new int[markers][6];	//A, T, G, C, Ind, Del
	//	resultPhredScores = new int[markers][6];	//A, T, G, C, Ind, Del
	//	resultMappingScores = new int[markers];
	public static void getAlleleCountsPhredScoresMappingScores(String[] aSingleLineOfBamFile, String chr, int start, int stop, int thresholdFor3Alleles, int[][] output1_AlleleCounts, int[][] output2_PhredScores, int[][] output3_MappingScores) {
		int readPointer;
		int outputArrayPointer;
		int currentPosition;
		int lengthOfCurrentSegment;
		int outputArrayLength;
		int loop;
		int indexInBases;
		String[][] readSegments;
		int currentMappingScore;

		if (aSingleLineOfBamFile[2].equalsIgnoreCase("chr" + chr) && Integer.parseInt(aSingleLineOfBamFile[1]) < 512 && ! aSingleLineOfBamFile[9].equalsIgnoreCase("*")) {
			currentPosition = Integer.parseInt(aSingleLineOfBamFile[3]);
			currentMappingScore = Integer.parseInt(aSingleLineOfBamFile[4]);
			readSegments = ext.getOperatorsOperatorIndicesAndSplit(aSingleLineOfBamFile[5], "MIDNSHP=X");
			readPointer = 0;
			outputArrayPointer = currentPosition - start;
			outputArrayLength = stop - start + 1;

			for (int i = 0; (outputArrayPointer < outputArrayLength) && (i < readSegments[0].length); i++) {
				lengthOfCurrentSegment = Integer.parseInt(readSegments[2][i]);
				if (readSegments[0][i].equals("M") || readSegments[0][i].equals("=") || readSegments[0][i].equals("X")) {	// present in read sequence, and also in reference sequence.
					if ((outputArrayPointer + lengthOfCurrentSegment) <= 0) {
						currentPosition += lengthOfCurrentSegment;
						outputArrayPointer += lengthOfCurrentSegment;
						readPointer += lengthOfCurrentSegment;
					} else {
						if (outputArrayPointer < 0) {
							currentPosition -= outputArrayPointer;
							readPointer -= outputArrayPointer;
							lengthOfCurrentSegment += outputArrayPointer;
							outputArrayPointer = 0;
						}

						loop = outputArrayPointer + Math.min(lengthOfCurrentSegment, outputArrayLength - outputArrayPointer);
						while (outputArrayPointer < loop) {
							indexInBases = ext.indexOfChar(aSingleLineOfBamFile[9].charAt(readPointer), BASES_WITH_N);
							if (indexInBases >= 0) {
								output1_AlleleCounts[outputArrayPointer][indexInBases] ++;
								output2_PhredScores[outputArrayPointer][indexInBases] += convertToPhredScore(aSingleLineOfBamFile[10].charAt(readPointer));
								output3_MappingScores[outputArrayPointer][indexInBases] += currentMappingScore;
							} else {
								System.err.println("Error - unrecognized base (" + aSingleLineOfBamFile[9].charAt(readPointer) + ") at chr" + chr + ":" + start + outputArrayPointer + " Read ID: " + aSingleLineOfBamFile[0]);
							}
							readPointer ++;
							outputArrayPointer ++;
						}

						currentPosition += lengthOfCurrentSegment;
					}

				} else if (readSegments[0][i].equals("I")) {	// present in read sequence, but NOT in reference sequence.
					if (outputArrayPointer >= 0) {
						output1_AlleleCounts[outputArrayPointer][INDEX_OF_INS] ++;
						output2_PhredScores[outputArrayPointer][INDEX_OF_INS] += convertToPhredScore(aSingleLineOfBamFile[10].charAt(readPointer));
						output3_MappingScores[outputArrayPointer][INDEX_OF_INS] += currentMappingScore;
					}
					readPointer += lengthOfCurrentSegment;

				} else if (readSegments[0][i].equals("D")) {	// NOT present in read sequence, but does in reference sequence.
					if ((outputArrayPointer + lengthOfCurrentSegment) <= 0) {
						currentPosition += lengthOfCurrentSegment;
						outputArrayPointer += lengthOfCurrentSegment;
					} else {
						if (outputArrayPointer < 0) {
							currentPosition -= outputArrayPointer;
//							readPointer -= outputArrayPointer;
							lengthOfCurrentSegment += outputArrayPointer;
							outputArrayPointer = 0;
						}

						loop = outputArrayPointer + Math.min(lengthOfCurrentSegment, outputArrayLength - outputArrayPointer);
						while (outputArrayPointer < loop) {
							output1_AlleleCounts[outputArrayPointer][INDEX_OF_DEL] ++;
							output2_PhredScores[outputArrayPointer][INDEX_OF_DEL] += DEFAULT_PHRED_SCORE_FOR_DELETION;
							output3_MappingScores[outputArrayPointer][INDEX_OF_DEL] += currentMappingScore;
							outputArrayPointer ++;
						}

						currentPosition += lengthOfCurrentSegment;
					}

				} else if (readSegments[0][i].equals("N")) {	// NOT present in read sequence, but does in reference sequence. Similar to D.
					currentPosition += lengthOfCurrentSegment;
					outputArrayPointer += lengthOfCurrentSegment;

				} else if (readSegments[0][i].equals("S")) {	// present in read sequence, and also in reference sequence, but do not match. Lower case letters in read sequence.
					if (i == 0) {
						readPointer += lengthOfCurrentSegment;
					} else {
						currentPosition += lengthOfCurrentSegment;
						outputArrayPointer += lengthOfCurrentSegment;
						readPointer += lengthOfCurrentSegment;
					}

				} else if (readSegments[0][i].equals("H")) {	// NOT present in read sequence, but in reference sequence. Similar to D.
					currentPosition += lengthOfCurrentSegment;
					outputArrayPointer += lengthOfCurrentSegment;

				} else if (readSegments[0][i].equals("P")) {	// present in read sequence, but NOT in reference sequence. Similar to I.

				} else {
					System.err.println("Unrecognized CIGAR string: " + readSegments[0][i]);
				}
			}
		}
	}

	public static int convertToPhredScore(char asciiString) {
		if(asciiString == '*') {
			return -1;
		} else {
			return asciiString - 33;
		}
	}

	public static int[] convertToPhredScore(String asciiString) {
		int[] result;
		int stringLength;

		if(asciiString.equals("*")) {
			result = null;
		} else {
			stringLength = asciiString.length();
			result = new int[stringLength];
			for (int i = 0; i < stringLength; i++) {
				result[i] = asciiString.charAt(i) - 33;
			}
		}

		return result;
	}

	public static String getForwardGenotypes(int[][][] alleleCounts, int currentMarkerIndex) {
		String result;
		int[] tempAlleleCounts;
		int[] tempIndices;
		int tempIndex1;
		int tempIndex2;
		int tempIndex3;

		tempAlleleCounts = new int[BASES.length]; 
		result = "";
		for (int i = 0; i < alleleCounts.length; i++) {
			tempIndices = Sort.quicksort(alleleCounts[i][currentMarkerIndex]);
			tempIndex1 = tempIndices.length - 1;
			while (tempIndex1 >= 0) {
				if (tempIndices[tempIndex1] >= BASES.length) {
					tempIndex1 --;
				} else {
					break;
				}
			}
			tempIndex2 = tempIndex1 - 1;
			while (tempIndex2 >= 0) {
				if (tempIndices[tempIndex2] >= BASES.length) {
					tempIndex2 --;
				} else {
					break;
				}
			}
			if (alleleCounts[i][currentMarkerIndex][tempIndices[tempIndex2]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
				tempIndex2 = tempIndex1;
			}
			
			if (tempIndices[tempIndex1] > tempIndices[tempIndex2]) {
				tempIndex3 = tempIndex2;
				tempIndex2 = tempIndex1;
				tempIndex1 = tempIndex3;
			}
			result += BASES[tempIndices[tempIndex1]] + "" + BASES[tempIndices[tempIndex2]] + (i == (alleleCounts.length - 1)? "" : ",");
		}
		return result;
	}

	public static String getNotes(String id, byte chr, int start, int stop, int[][] output1_AlleleCounts, int[][] output2_PhredScores, int[][] output3_MappingScores) {
		String notes;
		notes = "";
		
		return notes;
	}

	public static void getDenovoMutationScores(int[][][] alleleCounts, int[][][] phredScores, int[][][] mappingScores, byte[] output1DenovoMutationScores, String[] output2DenovoMutationNotes) {
		int loop;

		if (output1DenovoMutationScores.length != alleleCounts[0].length) {
			System.err.println("The length of the array outputDenovoMutationCandidateScores (" + output1DenovoMutationScores.length + ") is not consistent with that of the array alleleCounts (" + alleleCounts[0].length + ")");
		} else {
			getDenovoMutationScoresFromAlleleCounts(alleleCounts, output1DenovoMutationScores, output2DenovoMutationNotes);
			adjDenovoMutationScoresForNearbyIndels(alleleCounts, output1DenovoMutationScores, output2DenovoMutationNotes);
			adjDenovoMutationScoresForNearbyVariances(alleleCounts, output1DenovoMutationScores, output2DenovoMutationNotes);
			adjDenovoMutationScoresForNearbyDeNovoMutations(output1DenovoMutationScores, output2DenovoMutationNotes);
		}
	}

	public static void getDenovoMutationScoresFromAlleleCounts(int[][][] alleleCounts, byte[] output1DenovoMutationScores, String[] output2DenovoMutationNotes) {
		for (int i = 0; i < output1DenovoMutationScores.length; i++) {
			for (int j = 0; j < BASES.length; j++) {
				if (alleleCounts[0][i][j] > 0 && alleleCounts[1][i][j] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && alleleCounts[2][i][j] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
					for (int k = 0; k < alleleCounts.length; k++) {
						alleleCounts[k][i][INDEX_OF_TOTAL_READS] = alleleCounts[k][i][0] + alleleCounts[k][i][1] + alleleCounts[k][i][2] + alleleCounts[k][i][3] + alleleCounts[k][i][4] + alleleCounts[k][i][5];
					}

					if (alleleCounts[1][i][INDEX_OF_TOTAL_READS] > MIN_READ_DEPTH && alleleCounts[2][i][INDEX_OF_TOTAL_READS] > MIN_READ_DEPTH) {
						if (alleleCounts[0][i][j] > MIN_ALLELE_COUNT_FOR_DENOVO_MUTATION && alleleCounts[1][i][j] == 0 && alleleCounts[2][i][j] == 0) {
							output1DenovoMutationScores[i] = 100;

						} else if (alleleCounts[0][i][j] > MAX_ALLELE_COUNT_TREATED_AS_ZERO && alleleCounts[1][i][j] == 0 && alleleCounts[2][i][j] == 0) {
							output1DenovoMutationScores[i] = 80;

						} else if (alleleCounts[0][i][j] > (MIN_ALLELE_COUNT_FOR_DENOVO_MUTATION + 3 * (alleleCounts[1][i][j] + alleleCounts[2][i][j]))) {
							output1DenovoMutationScores[i] = 70;

						} else {
							output1DenovoMutationScores[i] = 25;
						}
		
						adjDenovoMutationScoresForMismatches(alleleCounts, i, output1DenovoMutationScores, output2DenovoMutationNotes);
						adjDenovoMutationScoresForThreeAlleles(alleleCounts, i, output1DenovoMutationScores, output2DenovoMutationNotes);
		
						break;
					}
				}
			}
		}
	}

	public static void adjDenovoMutationScoresForMismatches(int[][][] alleleCounts, int indexCurrentMarker, byte[] output1DenovoMutationScores, String[] output2DenovoMutationNotes) {
		String note;

		note = "";
		for (int j = 0; j < SAMPLE_SUFFIX.length; j++) {
			if (alleleCounts[j][indexCurrentMarker][INDEX_OF_N] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
				if (note.equals("")) {
					note = SAMPLE_SUFFIX[j];
				} else {
					note += "," + SAMPLE_SUFFIX[j];
				}
			}
		}

		if (! note.equals("")) {
			output1DenovoMutationScores[indexCurrentMarker] = (byte) (output1DenovoMutationScores[indexCurrentMarker] * DISCOUNT_FOR_N);

			if (output2DenovoMutationNotes[indexCurrentMarker] == null) {
				output2DenovoMutationNotes[indexCurrentMarker] = note + " mismatch(es)";
			} else {
				output2DenovoMutationNotes[indexCurrentMarker] += (";" + note + " mismatch(es)");
			}
		}
	}

	public static void adjDenovoMutationScoresForThreeAlleles(int[][][] alleleCounts, int indexCurrentMarker, byte[] output1DenovoMutationScores, String[] output2DenovoMutationNotes) {
		String note;
		double discountRate;

		note = "";
		discountRate = 1.0;
		updateNumAlleles(alleleCounts, indexCurrentMarker);
		for(int i = 0; i < SAMPLE_SUFFIX.length; i++) {
			if(alleleCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] > 2) {
				if (note.equals("")) {
					note = SAMPLE_SUFFIX[i];
				} else {
					note += "," + SAMPLE_SUFFIX[i];
				}

				if (alleleCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] > 2 || alleleCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] > 3) {
					discountRate = DISCOUNT_FOR_3ALLELES_LOOSE;
				} else {
					discountRate = DISCOUNT_FOR_3ALLELES_STRICT;
				}
			}
		}

		if (! note.equals("")) {
			output1DenovoMutationScores[indexCurrentMarker] *= discountRate;
			if (output2DenovoMutationNotes[indexCurrentMarker] == null) {
				output2DenovoMutationNotes[indexCurrentMarker] = note + ":3 alleles";
			} else {
				output2DenovoMutationNotes[indexCurrentMarker] += (note + ":3 alleles");
			}
		}
	}

	public static void updateNumAlleles(int[][][] alleleCounts, int indexCurrentMarker) {
		for(int i = 0; i < SAMPLE_SUFFIX.length; i++) {
			for (int j = 0; j < BASES.length; j++) {
				if (alleleCounts[i][indexCurrentMarker][j] > 0) {
					alleleCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] ++;
				} else if (alleleCounts[i][indexCurrentMarker][j] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
					alleleCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] ++;
				}
			}
		}
	}


	public static void adjDenovoMutationScoresForNearbyIndels(int[][][] alleleCounts, byte[] output1DenovoMutationScores, String[] output2DenovoMutationNotes) {
		byte[] indicesOfInDels;
		String note;
		int loop;
		double discountDifference;

		note = "";
		discountDifference = DISCOUNT_FOR_NEARBY_INDEL - DISCOUNT_FOR_ON_INDEL_SITE;
		indicesOfInDels = new byte[] {INDEX_OF_INS, INDEX_OF_DEL};
		for (int i = 0; i < output1DenovoMutationScores.length; i++) {
			for (int j = 0; j < indicesOfInDels.length; j++) {
				note = ""; 
				for (int k = 0; k < alleleCounts.length; k++) {
					if (alleleCounts[k][i][indicesOfInDels[j]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
						if (note.equals("")) {
							note = SAMPLE_SUFFIX[k];
						} else {
							note += "," + SAMPLE_SUFFIX[k];
						}
					}
	
				}
				if (! note.equals("")) {
					loop = Math.min(3, i);
					for (int k = 0; k < loop; k++) {
						if (output1DenovoMutationScores[i - k] > 0) {
							output1DenovoMutationScores[i - k] *= DISCOUNT_FOR_ON_INDEL_SITE;
							if (output2DenovoMutationNotes[i - k] == null) {
								output2DenovoMutationNotes[i - k] = note + ":" + BASES_WITH_N_INDEL[indicesOfInDels[j]] + " on site";
							} else if (output2DenovoMutationNotes[i - k].contains(BASES_WITH_N_INDEL[indicesOfInDels[j]] + " on site")) {
							} else {
								output2DenovoMutationNotes[i - k] += "; " + note + ":" + BASES_WITH_N_INDEL[indicesOfInDels[j]] + " on site";
							}
						}
					}
					loop = Math.min(WINDOW_SIZE_FOR_NEARBY_INDEL, i);
					for (int k = 3; k <= loop; k++) {
						if (output1DenovoMutationScores[i - k] > 0) {
							output1DenovoMutationScores[i - k] *= (DISCOUNT_FOR_ON_INDEL_SITE + Math.pow((k - 3) / (WINDOW_SIZE_FOR_NEARBY_INDEL - 3), 3) * discountDifference);
							if (output2DenovoMutationNotes[i - k] == null) {
								output2DenovoMutationNotes[i - k] = note + ":" + BASES_WITH_N_INDEL[indicesOfInDels[j]] + " nearby";
							} else if (output2DenovoMutationNotes[i - k].contains(BASES_WITH_N_INDEL[indicesOfInDels[j]] + " nearby")) {
							} else {
								output2DenovoMutationNotes[i - k] += "; " + note + ":" + BASES_WITH_N_INDEL[indicesOfInDels[j]] + " nearby";
							}
						}
					}
				}
			}
		}
	}

	public static void adjDenovoMutationScoresForNearbyVariances(int[][][] alleleCounts, byte[] output1DenovoMutationScores, String[] output2DenovoMutationNotes) {
		String note;
		int loop;

		for (int i = 0; i < output1DenovoMutationScores.length; i++) {
			if(alleleCounts[0][i][INDEX_OF_NUM_ALLELES_STRICT] == 0) {
				updateNumAlleles(alleleCounts, i);
			}

			if (output1DenovoMutationScores[i] > 0 && (alleleCounts[1][i][INDEX_OF_NUM_ALLELES_STRICT] == 2 || alleleCounts[2][i][INDEX_OF_NUM_ALLELES_STRICT] == 2)) {
				output1DenovoMutationScores[i] *= DISCOUNT_FOR_ON_VARIANCE_SITE;
				if (output2DenovoMutationNotes[i] == null) {
					output2DenovoMutationNotes[i] = "on variance site";
				} else {
					output2DenovoMutationNotes[i] += "; on variance site";
				}
			}

			note = "";
			for (int j = 0; j < alleleCounts.length; j++) {
				if (alleleCounts[j][i][INDEX_OF_NUM_ALLELES_LOOSE] == 2) {
					if (note.equals("")) {
						note = SAMPLE_SUFFIX[j];
					} else {
						note += "," + SAMPLE_SUFFIX[j];
					}
				}
			}
			if (! note.equals("")) {
				loop = Math.min(WINDOW_SIZE_FOR_NEARBY_VARIANCE, i);
				for (int j = 1; j <= loop; j++) {
					if (output1DenovoMutationScores[i - j] > 0) {
						output1DenovoMutationScores[i - j] *= DISCOUNT_FOR_NEARBY_VARIANCE;
						if (output2DenovoMutationNotes[i - j] == null) {
							output2DenovoMutationNotes[i - j] = "variance(s) nearby";
						} else if (output2DenovoMutationNotes[i - j].contains("variance(s) nearby")) {
						} else {
							output2DenovoMutationNotes[i - j] += "; variance(s) nearby";
						}
					}
				}
			}
		}
	}

	public static void adjDenovoMutationScoresForNearbyDeNovoMutations(byte[] denovoMutationScores, String[] outputDenovoMutationNotes) {
		int loop;

		for (int i = 0; i < denovoMutationScores.length; i++) {
			if (denovoMutationScores[i] >= 80) {
				loop = Math.min(WINDOW_SIZE_FOR_NEARBY_INDEL, i);
				for (int j = 1; j <= loop; j++) {
					if (denovoMutationScores[i - j] > 0) {
						denovoMutationScores[i - j] = (byte) (denovoMutationScores[i - j] * DISCOUNT_FOR_NEARBY_DENOVOMUTATION);
						if (outputDenovoMutationNotes[i - j] == null) {
							outputDenovoMutationNotes[i - j] = "DeNovo mutation nearby";
						} else if (outputDenovoMutationNotes[i - j].contains("DeNovo mutation nearby")) {
						} else {
							outputDenovoMutationNotes[i - j] += "; DeNovo mutation nearby";
						}
					}
				}
			}
		}
	}

	public static boolean[] isThreeAlleles(int[][][] alleleCounts, int indexCurrentMarker, int thresholdFor3Alleles) {
		boolean[] isThreeAlleles;
		int DenovoMutationScore;

		isThreeAlleles = new boolean[alleleCounts.length];
		for(int i = 0; i < alleleCounts.length; i++) {
			DenovoMutationScore = 0;
			for (int j = 0; j < BASES.length; j++) {
				if (alleleCounts[i][indexCurrentMarker][j] > thresholdFor3Alleles) {
					DenovoMutationScore ++;
				}
			}
			if(DenovoMutationScore > 2) {
				isThreeAlleles[i] = true;
			}
		}

		return	isThreeAlleles;
	}

	public static void exportInfoForPosition(PrintWriter writer, String id, String chr, int pos, char ref, byte call, int[][][] alleleCounts, int[][][] phredScores, int[][][] mappingScores, int indexOfCurrentMarker, String[] denovoMarkerNotes, String flag) {
		double totalPhredScores;
		double phredScoreProportion;
		String[][] phredScoreProportions;
		String forwardGenotypes;

		phredScoreProportions = new String[SAMPLE_SUFFIX.length][BASES.length];
		for (int i = 0; i < phredScoreProportions.length; i++) {
			totalPhredScores = phredScores[i][indexOfCurrentMarker][0] + phredScores[i][indexOfCurrentMarker][1] + phredScores[i][indexOfCurrentMarker][2] + phredScores[i][indexOfCurrentMarker][3];
			for (int j = 0; j < phredScoreProportions[i].length; j++) {
				phredScoreProportion = phredScores[i][indexOfCurrentMarker][j] / totalPhredScores;
				phredScoreProportions[i][j] = (phredScoreProportion==0? "-" : ext.formDeci(phredScoreProportion, 3));
			}
		}
		forwardGenotypes = getForwardGenotypes(alleleCounts, indexOfCurrentMarker);
		writer.println(id
				+ "\t" + chr
				+ "\t" + pos
				+ "\tchr" + chr + ":" + pos
				+ "\t" + ref
				+ "\t" + (forwardGenotypes.charAt(4) == forwardGenotypes.charAt(3)? forwardGenotypes.charAt(4) + "" : forwardGenotypes.charAt(7) + "")	//TODO alt allele
				+ "\t" + (call==100? "1.00" : "." + call)
				+ "\t(" + (alleleCounts[0][indexOfCurrentMarker][0]==0? "-" : alleleCounts[0][indexOfCurrentMarker][0]) + "," + (alleleCounts[0][indexOfCurrentMarker][1]==0? "-" : alleleCounts[0][indexOfCurrentMarker][1]) + "," + (alleleCounts[0][indexOfCurrentMarker][2]==0? "-" : alleleCounts[0][indexOfCurrentMarker][2]) + "," + (alleleCounts[0][indexOfCurrentMarker][3]==0? "-" : alleleCounts[0][indexOfCurrentMarker][3]) + "/" + (alleleCounts[0][indexOfCurrentMarker][4]==0? "-" : alleleCounts[0][indexOfCurrentMarker][4]) + "," + (alleleCounts[0][indexOfCurrentMarker][5]==0? "-" : alleleCounts[0][indexOfCurrentMarker][5]) + "," + (alleleCounts[0][indexOfCurrentMarker][6]==0? "-" : alleleCounts[0][indexOfCurrentMarker][6]) + ") (" + (alleleCounts[1][indexOfCurrentMarker][0]==0? "-" : alleleCounts[1][indexOfCurrentMarker][0]) + "," + (alleleCounts[1][indexOfCurrentMarker][1]==0? "-" : alleleCounts[1][indexOfCurrentMarker][1]) + "," + (alleleCounts[1][indexOfCurrentMarker][2]==0? "-" : alleleCounts[1][indexOfCurrentMarker][2]) + "," + (alleleCounts[1][indexOfCurrentMarker][3]==0? "-" : alleleCounts[1][indexOfCurrentMarker][3]) + "/" + (alleleCounts[1][indexOfCurrentMarker][4]==0? "-" : alleleCounts[1][indexOfCurrentMarker][4]) + "," + (alleleCounts[1][indexOfCurrentMarker][5]==0? "-" : alleleCounts[1][indexOfCurrentMarker][5]) + "," + (alleleCounts[1][indexOfCurrentMarker][6]==0? "-" : alleleCounts[1][indexOfCurrentMarker][6]) + ") (" + (alleleCounts[2][indexOfCurrentMarker][0]==0? "-" : alleleCounts[2][indexOfCurrentMarker][0]) + "," + (alleleCounts[2][indexOfCurrentMarker][1]==0? "-" : alleleCounts[2][indexOfCurrentMarker][1]) + "," + (alleleCounts[2][indexOfCurrentMarker][2]==0? "-" : alleleCounts[2][indexOfCurrentMarker][2]) + "," + (alleleCounts[2][indexOfCurrentMarker][3]==0? "-" : alleleCounts[2][indexOfCurrentMarker][3]) + "/" + (alleleCounts[2][indexOfCurrentMarker][4]==0? "-" : alleleCounts[2][indexOfCurrentMarker][4]) + "," + (alleleCounts[2][indexOfCurrentMarker][5]==0? "-" : alleleCounts[2][indexOfCurrentMarker][5]) + "," + (alleleCounts[2][indexOfCurrentMarker][6]==0? "-" : alleleCounts[2][indexOfCurrentMarker][6]) + ")" + (denovoMarkerNotes[indexOfCurrentMarker]==null? "" : (" " + denovoMarkerNotes[indexOfCurrentMarker]))
				+ "\t" + forwardGenotypes
				+ "\t" + flag
				+ "\t" + alleleCounts[0][indexOfCurrentMarker][INDEX_OF_TOTAL_READS]
				+ "\t" + alleleCounts[1][indexOfCurrentMarker][INDEX_OF_TOTAL_READS]
				+ "\t" + alleleCounts[2][indexOfCurrentMarker][INDEX_OF_TOTAL_READS]
				+ "\t(" + phredScoreProportions[0][0] + "," + phredScoreProportions[0][1] + "," + phredScoreProportions[0][2] + "," + phredScoreProportions[0][3] + ") (" + phredScoreProportions[1][0] + "," + phredScoreProportions[1][1] + "," + phredScoreProportions[1][2] + "," + phredScoreProportions[1][3] + ") (" + phredScoreProportions[2][0] + "," + phredScoreProportions[2][1] + "," + phredScoreProportions[2][2] + "," + phredScoreProportions[2][3] + ")"
				+ "\t" + (alleleCounts[0][indexOfCurrentMarker][INDEX_OF_TOTAL_READS] == 0? "" : ((mappingScores[0][indexOfCurrentMarker][0] + mappingScores[0][indexOfCurrentMarker][1] + mappingScores[0][indexOfCurrentMarker][2] + mappingScores[0][indexOfCurrentMarker][3]) / alleleCounts[0][indexOfCurrentMarker][INDEX_OF_TOTAL_READS]))
				+ "\t" + (alleleCounts[1][indexOfCurrentMarker][INDEX_OF_TOTAL_READS] == 0? "" : ((mappingScores[1][indexOfCurrentMarker][0] + mappingScores[1][indexOfCurrentMarker][1] + mappingScores[1][indexOfCurrentMarker][2] + mappingScores[1][indexOfCurrentMarker][3]) / alleleCounts[1][indexOfCurrentMarker][INDEX_OF_TOTAL_READS]))
				+ "\t" + (alleleCounts[2][indexOfCurrentMarker][INDEX_OF_TOTAL_READS] == 0? "" : ((mappingScores[2][indexOfCurrentMarker][0] + mappingScores[2][indexOfCurrentMarker][1] + mappingScores[2][indexOfCurrentMarker][2] + mappingScores[2][indexOfCurrentMarker][3]) / alleleCounts[2][indexOfCurrentMarker][INDEX_OF_TOTAL_READS]))
				+ "\t" + forwardGenotypes.charAt(0) + "\t" + forwardGenotypes.charAt(1) + "\t" + forwardGenotypes.charAt(3) + "\t" + forwardGenotypes.charAt(4) + "\t" + forwardGenotypes.charAt(6) + "\t" + forwardGenotypes.charAt(7)
				);
		writer.flush();
	}

	public static void displayErrorStream(Process p) {
		BufferedReader error;
        String line;

        error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
		try {
			while ((line = error.readLine()) != null) {
				System.out.println(line);
			}
			p.waitFor();
			error.close();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	public static void saveAlleleCountsToFile(String filename, int[][][] alleleCounts, int startPosition, boolean toIgnorePositionsWithZeroCount) {
		PrintWriter writer;
		String line;
		int total;

        try {
			writer = new PrintWriter(new FileWriter(filename));
			line = "position";
			for (int i = 0; i < SAMPLE_SUFFIX.length; i++) {
				line += ("\t" + SAMPLE_SUFFIX[i]);
			}
			writer.println(line);
			for (int i = 0; i < alleleCounts[0].length; i++) {
				line = startPosition + "";
				total = 0;
				for (int j = 0; j < alleleCounts.length; j++) {
					line += "\t(";
					for (int k = 0; k < BASES.length; k++) {
						line += ((k == 0? "" : ",") + (alleleCounts[j][i][k] == 0? "-" : alleleCounts[j][i][k]));
						total += alleleCounts[j][i][k];
					}
					line += "/";
					for (int k = BASES.length; k < BASES_WITH_N_INDEL.length; k++) {
						line += ((k == 0? "" : ",") + (alleleCounts[j][i][k] == 0? "-" : alleleCounts[j][i][k]));
						total += alleleCounts[j][i][k];
					}
					line += ")";
				}
				if (! (toIgnorePositionsWithZeroCount && total == 0)) {
					writer.println(line);
				}
				startPosition ++;
			}
			writer.close();
		} catch (Exception e) {
			System.out.println("Error writing to " + filename);
			System.out.println(e);
		}
	}

	/**
	 * Get the common root of all the filenames, assuming the filenames following the format of "*_*_*"
	 * @param filenames
	 * @return
	 */
	public static String getRootOf(String[] filenames) {
		String[][] roots;
		String commonRoot = null;
		int maxLength;
		int index1 = 0;
		int index2 = 0;
		boolean found;
		int loop;
		String str;

		roots = new String[filenames.length][];
		maxLength = 0;
		for (int i = 0; i < roots.length; i++) {
			roots[i] = ext.rootOf(filenames[i]).split("_");
			if (roots[i].length > maxLength) {
				maxLength = roots[i].length;
			}
		}

		found = false;
		for (int i = 0; i < maxLength; i++) {
			for (int j = 0; j < roots.length; j++) {
				if (roots[j].length > i) {
					commonRoot = roots[j][i];
					index1 = j;
					index2 = i;
					break;
				}
			}
			for (int j = index1; j < roots.length; j++) {
				if(! roots[j][i].equalsIgnoreCase(commonRoot)) {
					found = true;
					break;
				}
			}
			if (found) {
				break;
			}
		}

		if (found) {
			found = false;
			maxLength = 0;
			loop = commonRoot.length();
			for (int i = 0; i <= loop; i++) {
				for (int j = 0; j < roots.length; j++) {
					if (! roots[j][index2].substring(0, i).equalsIgnoreCase(commonRoot.substring(0, i))) {
						if (i > 0) {
							commonRoot = commonRoot.substring(0, i - 1);
						} else if (index2 > 0) { 
							commonRoot = roots[index1][index2 - 1] + commonRoot.substring(0, 0);
						}
						found = true;
						break;
					}
				}
				if (found) {
					break;
				}
			}
		} else {
			commonRoot = ext.rootOf(filenames[0]);
		}

		return commonRoot;
	}

	public static String[][] groupNamesByTrios(String[] filenames) {
		boolean[] isVisited;
		String currentTrio;
		Vector<String[]> result;
		String[][] result2;
		int location;
		int[] indices;

		isVisited = new boolean[filenames.length];
		indices = new int[3];
		for (int i = 0; i < indices.length; i++) {
			indices[i] = -1;
		}
		result = new Vector<String[]>(filenames.length/3);
		for (int i = 0; i < filenames.length; i++) {
			if (! isVisited[i] && filenames[i].contains("C_L")) {
				location = filenames[i].indexOf("C_L");
				currentTrio = filenames[i].substring(0, location);
				indices[0] = i;
				for (int j = 0; j < filenames.length; j++) {
					if (! isVisited[j] && i != j && filenames[j].substring(0, location).equals(currentTrio)) {
						isVisited[j] = true;

						if (filenames[j].contains("D_L")) {
							indices[1] = j;
						} else if (filenames[j].contains("M_L")) {
							indices[2] = j;
						} else {
							System.out.println("These file names do not make sense: " + filenames[j] + "\n" + filenames[indices[0]]);
						}
						
						if (indices[0] >= 0 && indices[1] >= 0 && indices[2] >= 0) {
							result.add(new String[] {filenames[indices[0]], filenames[indices[1]], filenames[indices[2]]});
							for (int k = 0; k < indices.length; k++) {
								indices[k] = -1;
							}
							break;
						}
					}
				}
			}
		}

		result2 = new String[result.size()][];
		for (int i = 0; i < result2.length; i++) {
			result2[i] = result.elementAt(i);
		}
		return result2;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String[] commands;
		String fileNameOfDeNovoPointMutationCandidateList;
		String[] bamFilenamesOfTheTrio;
		String bamDir;
		String scriptDir;
		String outputDir;
		boolean isToAnnotate;
		int numThreads;
		String bamFilenamesForHelpMenu;
		String chr;
		int start;
		int stop;
		String refFastaFilename;
		String bedFilename;
		Logger log;

		fileNameOfDeNovoPointMutationCandidateList = "N:/statgen/OS_Logan/IGV_validation/results_test1.txt";
		bamDir = "/home/spectorl/shared/exome_processing/bam/";
//		bamFilenamesOfTheTrio = new String[] {"rrd_F10639C_L008.bam", "rrd_F10639D_L007.bam", "rrd_F10639M_L007.bam"};
		bamFilenamesOfTheTrio = null;
//		bamDir = "D:/bams/";
//		bamFilenamesTheTrio = new String[] {"F10639C_chr17_39346502_39346622.txt", "F10639D_chr17_39346502_39346622.txt", "F10639M_chr17_39346502_39346622.txt"};
		refFastaFilename = "/home/spectorl/xuz2/hg19_canonical.fa";
		scriptDir = "/home/spectorl/xuz2/scripts/";
		outputDir = "/home/spectorl/xuz2/outputs/";
		isToAnnotate = false;
		numThreads = 1;
		chr = "";
		start = 39346600;
		stop = 39346622;
		bedFilename = "/home/spectorl/xuz2/outputs/S04380219_Regions.bed";

		bamFilenamesForHelpMenu = (bamFilenamesOfTheTrio == null || bamFilenamesOfTheTrio[0] == null)? "" : bamFilenamesOfTheTrio[0];
		for (int i = 1; bamFilenamesOfTheTrio != null && bamFilenamesOfTheTrio[i] != null && i < bamFilenamesOfTheTrio.length; i++) {
			bamFilenamesForHelpMenu += "," + bamFilenamesOfTheTrio[i];
		}

		commands = new String[] {"-annotation", "candidate=", "bamdir=", "scriptdir=", "outdir=", "bamset=", "reffasta=", "chr=", "start=", "stop=", "bed=", "numthreads="};
		String usage = "\n"
				+ "To annotate a list of candidate markers:"
				+ "   (1) command for annotatation (i.e. " + commands[0] + " (default))\n"
				+ "   (2) full path of the candidate list file (i.e. " + commands[1] + fileNameOfDeNovoPointMutationCandidateList + " (default))\n"
				+ "   (3) directory of the bam files (i.e. " + commands[2] + bamDir + " (default))\n"
				+ "   (4) directory for the output script file to extract information from bam files (i.e. " + commands[3] + scriptDir + " (default))\n"
				+ "To process a region:"
				+ "   (1) chr of the region (i.e. " + commands[7] + chr + " (default))\n"
				+ "   (2) start position of the region (i.e. " + commands[8] + start + " (default))\n"
				+ "   (3) stop position of the region (i.e. " + commands[9] + stop + " (default))\n"
				+ "   (4) directory of the bam files (i.e. " + commands[2] + bamDir + " (default))\n"
				+ "   (4) names of bam/sam files of the trio, if just for a trio (i.e. " + commands[5] + bamFilenamesForHelpMenu + " (default))\n"
				+ "   (5) full path of the reference Fasta file (i.e. " + commands[6] + refFastaFilename + " (default))\n"
				+ "   (6) directory for the output files (i.e. " + commands[4] + outputDir + " (default))\n"
				+ "To process genome:"
				+ "   (1) full path of the bed file (i.e. " + commands[10] + bedFilename + " (default))\n"
				+ "   (2) directory of the bam files (i.e. " + commands[2] + bamDir + " (default))\n"
				+ "   (3) names the bam/sam files of the trio, if just for a trio (i.e. " + commands[5] + bamFilenamesForHelpMenu + " (default))\n"
				+ "   (4) full path of the reference Fasta file (i.e. " + commands[7] + refFastaFilename + " (default))\n"
				+ "   (5) directory for the output files (i.e. " + commands[4] + outputDir + " (default))\n"
				+ "   (6) number of threads (i.e. " + commands[11] + numThreads + " (default))\n"
				+ "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith(commands[0])) {
				isToAnnotate = true;
				numArgs--;
			} else if (args[i].startsWith(commands[1])) {
				fileNameOfDeNovoPointMutationCandidateList = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[2])) {
				bamDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[3])) {
				scriptDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[4])) {
				outputDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[5])) {
				bamFilenamesOfTheTrio = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith(commands[6])) {
				refFastaFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[7])) {
				chr = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[8])) {
				start = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith(commands[9])) {
				stop = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith(commands[10])) {
				if (args[i].split("=").length < 2) {
					bedFilename = null;
				} else {
					bedFilename = args[i].split("=")[1];
				}
				numArgs--;
			} else if (args[i].startsWith(commands[11])) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		log = new Logger(outputDir + "DeNovoSeq_" + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log");
		
		if (isToAnnotate) {
//			getAlleleCounts("D:/bams/F10639C_chr11_89531764.txt", 11, 89531764, 3);	//Testing only. Feature already included in screenDeNovoPointMutation()
//			getAlleleCounts("D:/bams/F10639C_chr5_140558628.txt", 5, 140558628, 3);
			generateScriptForSamtools(fileNameOfDeNovoPointMutationCandidateList, bamFilenamesOfTheTrio[0], scriptDir);
			screenDeNovoPointMutation(fileNameOfDeNovoPointMutationCandidateList, bamFilenamesOfTheTrio[0], 1);
		} else if (chr != null && !chr.equals("") && start > 0 && stop > 0 && stop >= start) {
			processRegion(bamDir, bamFilenamesOfTheTrio, refFastaFilename, chr, start, stop, outputDir, false);
		} else if (bedFilename != null && new File(bedFilename).exists()) {
			if (bamFilenamesOfTheTrio == null || bamFilenamesOfTheTrio.length < 2) {
				processGenomeOfAllTriosInDir(bamDir, refFastaFilename, bedFilename, outputDir, scriptDir, numThreads, log);
			} else {
				processGenomeOfOneTrio(bamDir, bamFilenamesOfTheTrio, refFastaFilename, bedFilename, outputDir, numThreads, log);
			}
		} else {
			log.reportError("No command executed, due to:\n1) " + commands[0] + " is not specified; or 2) chr, start and stop are not complete; or 3) .bed file is not found");
		}
	}

}
