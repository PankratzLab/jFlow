package one;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
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
import common.Sort;
import common.ext;

public class DeNovoSeq {
	public static final String[] SAMPLE_SUFFIX = new String[] {"C", "D", "M"};
	public static final char[] BASES = new char[] {'A', 'T', 'G', 'C'};
	public static final char[] BASES_WITH_N = new char[] {'A', 'T', 'G', 'C', 'N'};
	public static final String[] BASES_WITH_N_INS_DEL = new String[] {"A", "T", "G", "C", "N", "Ins", "Del"};
	public static final byte INDEX_OF_N = 4;
	public static final byte INDEX_OF_INS = 5;
	public static final byte INDEX_OF_DEL = 6;
	public static final byte MARKERDATA_NUMSAMPLES_START = 0;
	public static final byte MARKERDATA_NUMSAMPLES_LEN = 4;
	public static final int DEFAULT_PHRED_SCORE_FOR_DELETION = 30;
	public static final double THRESHOLD_PHRED_SCORE_FOR_INS_DEL = .10;
	public static final int MAX_ALLELE_COUNT_TREATED_AS_ZERO = 2;
	public static final int MIN_ALLELE_COUNT_FOR_DENOVO_MUTATION = 5;
	public static final double MIN_ALLELE_FREQ_FOR_DENOVO_MUTATION = .20;
	public static final int MIN_READ_DEPTH = 10;
	public static final int NEIGHBOR = 60;
	public static final double DISCOUNT_FOR_NEIGHBOR_MUTATION_OR_INDEL = .75;
	public static final double DISCOUNT_FOR_N = .50;
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

//	public static void processGenome(String bamDir, String[] bamFilenames, String refFastaFilename, String bedFilename, String outputDir, Logger log) {
//		BufferedReader reader;
//		String[] line;
//		byte chr;
//		int start;
//		int stop;
//		PrintWriter writer;
//		String trioId;
//		String outFileName;
//		long timer;
//		SimpleDateFormat timeFormat;
//
//		if (log == null) {
//			log = new Logger();
//		}
//
//		trioId = getRootOf(bamFilenames);
//		outFileName = outputDir + trioId + "_denovoMutations_" + ext.rootOf(bedFilename) + ".txt";
//        timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
//		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
//       	timer = new Date().getTime();
//		try {
//			writer = new PrintWriter(outFileName);
////			writer.println("id\tchr\tpos\tlookup\tsarver\tref\talt\tmendelianLikelihood\tmendelianPP\tmendelianGT\tsnpCode\tcode\tdeNovoLikelihood\tdeNovoPP\tactualDenovo\tconf\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tchildQuality\tdadQuality\tmomQuality");
//			writer.println("id\tchr\tpos\tlookup\tref\talt\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore\t1\t2\t4\t5\t7\t8");
//			reader = Files.getAppropriateReader(bedFilename);
//			reader.readLine();
//			reader.readLine();
//			while (reader.ready()) {
//				line = reader.readLine().split("\t");
//				chr = Byte.parseByte(line[0].substring(3));
//				start = Integer.parseInt(line[1]);
//				stop = Integer.parseInt(line[2]);
//				while (start + REGION_LEN_AT_A_TIME <= stop) {
//					processRegion(bamDir, bamFilenames, trioId, refFastaFilename, chr, start, start + REGION_LEN_AT_A_TIME, writer, null);
//					start = start + REGION_LEN_AT_A_TIME + 1;
//				}
//				if (stop > start) {
//					processRegion(bamDir, bamFilenames, trioId, refFastaFilename, chr, start, stop, writer, null);
//				}
//			}
//			reader.close();
//			writer.close();
//		} catch (FileNotFoundException fnfe) {
//			log.reportError("Error: file \"" + bedFilename + "\" not found in current directory");
//			return;
//		} catch (IOException ioe) {
//			log.reportError("Error reading file \"" + bedFilename + "\"");
//			return;
//		}
//	}

	public static void processGenomeOfAllTriosInDir(String bamDir, String refFastaFilename, String bedFilename, String outputDir, int numThreads, Logger log) {
		String[] bamFilenames;
		String[][] bamFilenamesByTrios;
//		String trioId;
//		String outFileName;
//		PrintWriter writer;

		bamFilenames = Files.list(bamDir, ".bam", false);
		bamFilenamesByTrios = getNameByTrios(bamFilenames);
		for (int i = 0; i < bamFilenamesByTrios.length; i++) {
			processGenomeOfOneTrio(bamDir, bamFilenames, refFastaFilename, bedFilename, outputDir, numThreads, log);
		}

//		trioId = getRootOf(bamFilenamesByTrios[0]);
//		outFileName = outputDir + trioId + "_denovoMutations_genome" + ext.rootOf(bedFilename) + ".txt";
//		try {
//			writer = new PrintWriter(outFileName);
////			writer.println("id\tchr\tpos\tlookup\tsarver\tref\talt\tmendelianLikelihood\tmendelianPP\tmendelianGT\tsnpCode\tcode\tdeNovoLikelihood\tdeNovoPP\tactualDenovo\tconf\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tchildQuality\tdadQuality\tmomQuality");
//			writer.println("id\tchr\tpos\tlookup\tref\talt\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore\t1\t2\t4\t5\t7\t8");
//			for (int i = 0; i < bamFilenamesByTrios.length; i++) {
//				processGenomeOfOneTrio(bamDir, bamFilenames, refFastaFilename, bedFilename, writer, numThreads, log);
//			}
//			writer.close();
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
	}

//	public static void processGenomeOfOneTrio(String bamDir, String[] bamFilenamesOfTheTrio, String refFastaFilename, String bedFilename, String outputDir, int numThreads, Logger log) {
//		String trioId;
//		String outFileName;
//		PrintWriter writer;
//		long timer;
//		SimpleDateFormat timeFormat;
//
//		trioId = getRootOf(bamFilenamesOfTheTrio);
//		outFileName = outputDir + trioId + "_denovoSeq_" + ext.rootOf(bedFilename) + ".txt";
//        timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
//		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
//       	timer = new Date().getTime();
//		try {
//			writer = new PrintWriter(outFileName);
////			writer.println("id\tchr\tpos\tlookup\tsarver\tref\talt\tmendelianLikelihood\tmendelianPP\tmendelianGT\tsnpCode\tcode\tdeNovoLikelihood\tdeNovoPP\tactualDenovo\tconf\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tchildQuality\tdadQuality\tmomQuality");
//			writer.println("id\tchr\tpos\tlookup\tref\talt\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore\t1\t2\t4\t5\t7\t8");
//			processGenomeOfOneTrio(bamDir, bamFilenamesOfTheTrio, refFastaFilename, bedFilename, writer, numThreads, log);
//			writer.close();
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
//
//		log.report("DeNovo mutation result is ready at: " + outFileName + "\nTotal time used " + timeFormat.format(new Date().getTime() - timer));
//	}

	public static void processGenomeOfOneTrio(String bamDir, String[] bamFilenamesOfTheTrio, String refFastaFilename, String bedFilename, String outputDir, int numThreads, Logger log) {
		BufferedReader reader;
		String[] line;
		String chr, prevChr;
		int start;
		int stop;
		int loop;
		PrintWriter writer;
		String trioId;
		String outFileName;
		ExecutorService executor = null;
		int processId;

		if (log == null) {
			log = new Logger();
		}

		trioId = getRootOf(bamFilenamesOfTheTrio);
		outFileName = outputDir + trioId + "_denovoMutations_genome" + ".txt";
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
			}
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				chr = line[0].substring(3);
				if (!chr.equals(prevChr) && numThreads == 1) {
					System.out.println(ext.getTime()+"\t"+"Starting chr"+chr+"; identified "+Files.countLines(outFileName, true)+" possible de novo events so far");
					prevChr = chr;
				}
				start = Integer.parseInt(line[1]);
				stop = Integer.parseInt(line[2]);
				while (start <= stop) {
					loop = Math.min(start + REGION_LEN_AT_A_TIME, stop);
					if (numThreads > 1) {
						Runnable worker = new WorkerThread(bamDir, bamFilenamesOfTheTrio, trioId, refFastaFilename, chr, start, loop, writer, null, processId);
						executor.execute(worker);
						processId ++;
					} else {
						System.out.println(chr + "\t" + start + "\t" + loop);
						processRegion(bamDir, bamFilenamesOfTheTrio, trioId, refFastaFilename, chr, start, loop, writer, null);
					}
					start += (REGION_LEN_AT_A_TIME + 1);
				}
			}
			if (numThreads > 1) {
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
	}

//	public static void processGenomeMultiThread(String bamDir, String[] bamFilenames, String refFastaFilename, String bedFilename, String outputDir, int numThreads, Logger log) {
//		BufferedReader reader;
//		String[] line;
//		byte chr;
//		int start;
//		int stop;
//		PrintWriter writer;
//		String trioId;
//		String outFileName;
//		long timer;
//		SimpleDateFormat timeFormat;
//		int processId;
//
//		if (log == null) {
//			log = new Logger();
//		}
//
//		trioId = getRootOf(bamFilenames);
//		outFileName = outputDir + trioId + "_denovoMutations_genome" + ".txt";
//		processId = 0;
//        timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
//		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
//       	timer = new Date().getTime();
//		try {
//			writer = new PrintWriter(outFileName);
////			writer.println("id\tchr\tpos\tlookup\tsarver\tref\talt\tmendelianLikelihood\tmendelianPP\tmendelianGT\tsnpCode\tcode\tdeNovoLikelihood\tdeNovoPP\tactualDenovo\tconf\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tchildQuality\tdadQuality\tmomQuality");
//			writer.println("id\tchr\tpos\tlookup\tref\talt\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore\t1\t2\t4\t5\t7\t8");
//			reader = Files.getAppropriateReader(bedFilename);
//			reader.readLine();
//			reader.readLine();
//			ExecutorService executor = Executors.newFixedThreadPool(numThreads);
//			while (reader.ready()) {
//				line = reader.readLine().split("\t");
//				chr = Byte.parseByte(line[0].substring(3));
//				start = Integer.parseInt(line[1]);
//				stop = Integer.parseInt(line[2]);
//				while (start + REGION_LEN_AT_A_TIME <= stop) {
//					Runnable worker = new WorkerThread(bamDir, bamFilenames, trioId, refFastaFilename, chr, start, start + REGION_LEN_AT_A_TIME, writer, null, processId);
//					executor.execute(worker);
//					start = start + REGION_LEN_AT_A_TIME + 1;
//					processId ++;
//				}
//				if (stop > start) {
//					Runnable worker = new WorkerThread(bamDir, bamFilenames, trioId, refFastaFilename, chr, start, stop, writer, null, processId);
//					executor.execute(worker);
//					processId ++;
//				}
//			}
//			reader.close();
//			writer.close();
//			executor.awaitTermination(7, TimeUnit.DAYS);
//			executor.shutdown();
//		} catch (FileNotFoundException fnfe) {
//			log.reportError("Error: file \"" + bedFilename + "\" not found in current directory");
//			return;
//		} catch (IOException ioe) {
//			log.reportError("Error reading file \"" + bedFilename + "\"");
//			return;
//		} catch (InterruptedException e) {
//		}
//	}

	public static class WorkerThread implements Runnable {
		private String bamDir;
		private String[] bamFilenames;
		private String trioId;
		private String refFastaFilename;
		private String chr;
		private int start;
		private int stop;
		private PrintWriter writer;
		private String outAlleleCountsFileName;
		private int threadId;

		public WorkerThread(String bamDir, String[] bamFilenames, String trioId, String refFastaFilename, String chr, int start, int stop, PrintWriter writer, String outAlleleCountsFileName, int threadId) {
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
			processRegion(bamDir, bamFilenames, trioId, refFastaFilename, chr, start, stop, writer, outAlleleCountsFileName);
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
		PrintWriter writer2;
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
		int[][] depths;
/*		int[] depthes;
		double totalPhredScores;
		String[][] phredScoresProportions;
*/		String refAlleles;
		byte[] denovoMutationScores;
		String[] denovoMarkerNotes;

        try {
            numMarkers = stop - start + 1;
			alleleCounts = new int[SAMPLE_SUFFIX.length][numMarkers][BASES_WITH_N_INS_DEL.length];
			phredScores = new int[SAMPLE_SUFFIX.length][numMarkers][BASES_WITH_N_INS_DEL.length];
			mappingScores = new int[SAMPLE_SUFFIX.length][numMarkers][BASES_WITH_N_INS_DEL.length];

			for (int i = 0; i < bamFilenames.length; i++) {
				bamContentVec = new Vector<String[]>();

				p = Runtime.getRuntime().exec("samtools view " + bamFilenames[i] + " chr" + chr + ":" + start + "-" + stop, null, new File(bamDir));

//				ps=new ProcessBuilder("samtools", "view", bamFilenames[i], "chr" + chr + ":" + start + "-" + stop);
//		        ps.redirectErrorStream(true);
//		        ps.directory(new File(bamDir));
//		        p = ps.start();  

//		        p.waitFor();
				reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
/*				reader = new BufferedReader(new FileReader(bamDir + bamFilenames[i]));
*/
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
				saveToFile(outAlleleCountsFileName, alleleCounts, start, true);
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

	        depths = new int[SAMPLE_SUFFIX.length][numMarkers];
			denovoMutationScores = new byte[numMarkers];
			denovoMarkerNotes = new String[numMarkers];
			getDenovoMutationScores(alleleCounts, phredScores, mappingScores, depths, denovoMutationScores, denovoMarkerNotes);
/*			depthes = new int[NUMSAMPLES_IN_A_TRIO];
			phredScoresProportions = new double[NUMSAMPLES_IN_A_TRIO][BASES.length];
*/			for (int i = 0; i < numMarkers; i++) {
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
/*	public static void getAlleleCountsPhredScoresMappingScores(String[] aSingleLineOfBamFile, int chr, int start, int stop, int thresholdFor3Alleles, int[][] output1_AlleleCounts, int[][] output2_PhredScores, int[][] output3_MappingScores) {
//		char[] bases;
		int pointerOfCurrentSearchOnTheRead, index2, startPositionOfTheRead, lengthOfCurrentSegmentOfCurrentRead, positionToLookForOnCurrentRead, numMarkers, currentPosition;
		String[][] operatorsOperatorIndicesAndSplit;
		int lengthOfTheRead;
		int currentMappingScore;

		numMarkers = stop - start + 1;
//		lengthOfTheRead = Math.abs(Integer.parseInt(aSingleLineOfBamFile[8]));
		lengthOfTheRead = aSingleLineOfBamFile[9].length();
		currentMappingScore = Integer.parseInt(aSingleLineOfBamFile[4]);
		currentPosition = start;
		for (int i = 0; i < numMarkers; i++) {
			positionToLookForOnCurrentRead = currentPosition;
			if (aSingleLineOfBamFile[2].equalsIgnoreCase("chr" + chr) && Integer.parseInt(aSingleLineOfBamFile[3])<=start && (Integer.parseInt(aSingleLineOfBamFile[3]) + lengthOfTheRead) >=start && Integer.parseInt(aSingleLineOfBamFile[1]) < 512) {
				startPositionOfTheRead = Integer.parseInt(aSingleLineOfBamFile[3]);
				pointerOfCurrentSearchOnTheRead = 0;
				operatorsOperatorIndicesAndSplit = ext.getOperatorsOperatorIndicesAndSplit(aSingleLineOfBamFile[5], "DIMNPSH");
				for (int j = 0; j < operatorsOperatorIndicesAndSplit[0].length; j++) {
					lengthOfCurrentSegmentOfCurrentRead = Integer.parseInt(operatorsOperatorIndicesAndSplit[2][j]);
					if ((startPositionOfTheRead + pointerOfCurrentSearchOnTheRead + lengthOfCurrentSegmentOfCurrentRead) > positionToLookForOnCurrentRead) {
						if (operatorsOperatorIndicesAndSplit[0][j].equals("I")) {
							positionToLookForOnCurrentRead += lengthOfCurrentSegmentOfCurrentRead;
							pointerOfCurrentSearchOnTheRead += lengthOfCurrentSegmentOfCurrentRead;
							if ((positionToLookForOnCurrentRead - (startPositionOfTheRead + pointerOfCurrentSearchOnTheRead)) <= 1) {
								output1_AlleleCounts[i][4] ++;
								output2_PhredScores[i][4] += DEFAULT_PHRED_SCORE_FOR_INS_DEL;
								output3_MappingScores[i][4] += currentMappingScore; 
							}
						} else if (operatorsOperatorIndicesAndSplit[0][j].equals("S") || operatorsOperatorIndicesAndSplit[0][j].equals("H") || operatorsOperatorIndicesAndSplit[0][j].equals("N")) {
							positionToLookForOnCurrentRead += lengthOfCurrentSegmentOfCurrentRead;
							pointerOfCurrentSearchOnTheRead += lengthOfCurrentSegmentOfCurrentRead;
						} else if (operatorsOperatorIndicesAndSplit[0][j].equals("D")) {
							output1_AlleleCounts[i][5] ++;
							output2_PhredScores[i][5] += DEFAULT_PHRED_SCORE_FOR_INS_DEL;
							output3_MappingScores[i][5] += currentMappingScore; 
							break;
						} else if (operatorsOperatorIndicesAndSplit[0][j].equals("P")) {
							break;
						} else {
							pointerOfCurrentSearchOnTheRead += (positionToLookForOnCurrentRead - startPositionOfTheRead - pointerOfCurrentSearchOnTheRead);
//								index = (positionOnThisRead - startPositionOfTheRead);
							index2 = ext.indexOfChar(aSingleLineOfBamFile[9].charAt(pointerOfCurrentSearchOnTheRead), BASES);
							if (index2 >= 0) {
								output1_AlleleCounts[i][index2] ++;
								output2_PhredScores[i][index2] += convertToPhredScore(aSingleLineOfBamFile[10].charAt(pointerOfCurrentSearchOnTheRead));
								output3_MappingScores[i][index2] += currentMappingScore; 
							} else {
								System.out.println("Error - unrecognized base (neither A, T, G, nor C) at chr" + chr + ":" + positionToLookForOnCurrentRead + " Read ID: " + aSingleLineOfBamFile[0]);
							}
							break;
						}
					} else {
						if (operatorsOperatorIndicesAndSplit[0][j].equals("I") || operatorsOperatorIndicesAndSplit[0][j].equals("S") || operatorsOperatorIndicesAndSplit[0][j].equals("H") || operatorsOperatorIndicesAndSplit[0][j].equals("N")) {
							positionToLookForOnCurrentRead += lengthOfCurrentSegmentOfCurrentRead;
							pointerOfCurrentSearchOnTheRead += lengthOfCurrentSegmentOfCurrentRead;
						} else if (operatorsOperatorIndicesAndSplit[0][j].equals("D") || operatorsOperatorIndicesAndSplit[0][j].equals("P")) {
							startPositionOfTheRead += lengthOfCurrentSegmentOfCurrentRead;
						} else {
							pointerOfCurrentSearchOnTheRead += lengthOfCurrentSegmentOfCurrentRead;
						}
					}
				}
			}
			currentPosition ++;
		}
	}
*/
/*	public static void getAlleleCountsPhredScoresMappingScores(String[] aSingleLineOfBamFile, int chr, int start, int stop, int thresholdFor3Alleles, int[][] output1_AlleleCounts, int[][] output2_PhredScores, int[][] output3_MappingScores) {
		int readPointer;
		int outputArrayPointer;
		int positionPointer;
		int lengthOfCurrentSegment;
		int indexInBases;
		String[][] readSegments;
		int currentMappingScore;
		boolean isFirstPosition;

		if (aSingleLineOfBamFile[2].equalsIgnoreCase("chr" + chr) && Integer.parseInt(aSingleLineOfBamFile[1]) < 512) {
			positionPointer = Integer.parseInt(aSingleLineOfBamFile[3]);
			currentMappingScore = Integer.parseInt(aSingleLineOfBamFile[4]);
			readSegments = ext.getOperatorsOperatorIndicesAndSplit(aSingleLineOfBamFile[5], "DIMNPSH");
			readPointer = 0;
			outputArrayPointer = Math.max(0, positionPointer - start);
			isFirstPosition = true;
			for (int i = 0; i < readSegments[0].length; i++) {
				lengthOfCurrentSegment = Integer.parseInt(readSegments[2][i]);
				if (readSegments[0][i].equals("I")) {
					if (Math.abs(start + outputArrayPointer - positionPointer) <= 1) {
						output1_AlleleCounts[outputArrayPointer][INDEX_OF_INS] ++;
						output2_PhredScores[outputArrayPointer][INDEX_OF_INS] += convertToPhredScore(aSingleLineOfBamFile[10].charAt(readPointer));	//TODO
						output3_MappingScores[outputArrayPointer][INDEX_OF_INS] += currentMappingScore;
					}
					readPointer += lengthOfCurrentSegment;
				} else if (readSegments[0][i].equals("S") || readSegments[0][i].equals("H") || readSegments[0][i].equals("N")) {
					readPointer += lengthOfCurrentSegment;
				} else if (readSegments[0][i].equals("D")) {
					if ((start + outputArrayPointer - positionPointer) <= lengthOfCurrentSegment) {
						while ((start + outputArrayPointer) <= stop && lengthOfCurrentSegment > 0) {
							output1_AlleleCounts[outputArrayPointer][INDEX_OF_DEL] ++;
							output2_PhredScores[outputArrayPointer][INDEX_OF_DEL] += DEFAULT_PHRED_SCORE_FOR_DELETION;
							output3_MappingScores[outputArrayPointer][INDEX_OF_DEL] += currentMappingScore;
							positionPointer ++;
							outputArrayPointer ++;
							lengthOfCurrentSegment --;
						}
					} else {
						positionPointer += lengthOfCurrentSegment;
					}
				} else if (readSegments[0][i].equals("P")) {
					positionPointer += lengthOfCurrentSegment;
				} else if (readSegments[0][i].equals("M")) {
					if ((start + outputArrayPointer - positionPointer) < lengthOfCurrentSegment) {
						if (isFirstPosition) {
							readPointer += (start + outputArrayPointer - positionPointer);
							lengthOfCurrentSegment -= (start + outputArrayPointer - positionPointer);
							positionPointer = start + outputArrayPointer;
							isFirstPosition = false;
						}
						while ((start + outputArrayPointer) <= stop && lengthOfCurrentSegment > 0) {
							indexInBases = ext.indexOfChar(aSingleLineOfBamFile[9].charAt(readPointer), BASES_WITH_N);
//							if (positionPointer == 39346603) {
//								System.out.println(aSingleLineOfBamFile[0] + "\t" + aSingleLineOfBamFile[3] + "\t" + aSingleLineOfBamFile[5] + "\t" + aSingleLineOfBamFile[9] + "\t" + readPointer + "\t" + outputArrayPointer + "\t" + (outputArrayPointer + start) + "\t" + BASES[index2]);
//							}
							if (indexInBases >= 0) {
								output1_AlleleCounts[outputArrayPointer][indexInBases] ++;
								output2_PhredScores[outputArrayPointer][indexInBases] += convertToPhredScore(aSingleLineOfBamFile[10].charAt(readPointer));
								output3_MappingScores[outputArrayPointer][indexInBases] += currentMappingScore;
							} else {
								System.err.println("Error - unrecognized base (" + aSingleLineOfBamFile[9].charAt(readPointer) + ") at chr" + chr + ":" + start + outputArrayPointer + " Read ID: " + aSingleLineOfBamFile[0]);
							}
							positionPointer ++;
							readPointer ++;
							outputArrayPointer ++;
							lengthOfCurrentSegment --;
						}
					} else {
						positionPointer += lengthOfCurrentSegment;
						readPointer += lengthOfCurrentSegment;
					}
				} else {
					System.err.println("Unrecognized CIGAR string: " + readSegments[0][i]);
				}
				
				if(start + outputArrayPointer > stop) {
					break;
				}
			}
		}

//		if (readSegments[0][i].equals("I") || readSegments[0][i].equals("S") || readSegments[0][i].equals("H") || readSegments[0][i].equals("N")) {
//		} else if (readSegments[0][i].equals("D") || readSegments[0][i].equals("P")) {
//		} else {
//		}
	}
*/
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

	public static void getDenovoMutationScores(int[][][] alleleCounts, int[][][] phredScores, int[][][] mappingScores, int[][] output1ReadDepths, byte[] output2DenovoMutationScores, String[] output3DenovoMutationNotes) {
		boolean[] isThreeAlleles;
		String output3DenovoMutationNote;
		int loop;
		byte[] indicesOfInDels;

		if (output2DenovoMutationScores.length != alleleCounts[0].length) {
			System.err.println("The length of the array outputDenovoMutationCandidateScores (" + output2DenovoMutationScores.length + ") is not consistent with that of the array alleleCounts (" + alleleCounts[0].length + ")");
		} else {
			for (int i = 0; i < output2DenovoMutationScores.length; i++) {
				for (int j = 0; j < alleleCounts.length; j++) {
					output1ReadDepths[j][i] = alleleCounts[j][i][0] + alleleCounts[j][i][1] + alleleCounts[j][i][2] + alleleCounts[j][i][3] + alleleCounts[j][i][4] + alleleCounts[j][i][5];
				}
				if (output1ReadDepths[0][i] > MIN_READ_DEPTH && output1ReadDepths[1][i] > MIN_READ_DEPTH && output1ReadDepths[2][i] > MIN_READ_DEPTH) {
					for (int j = 0; j < BASES.length; j++) {
						if ((alleleCounts[0][i][j] > MIN_ALLELE_COUNT_FOR_DENOVO_MUTATION && alleleCounts[1][i][j] == 0 && alleleCounts[2][i][j] == 0)
								|| (alleleCounts[1][i][j] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && alleleCounts[2][i][j] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && (alleleCounts[1][i][j] + alleleCounts[2][i][j]) > 0 && alleleCounts[0][i][j] > (MIN_ALLELE_COUNT_FOR_DENOVO_MUTATION + 2 * (alleleCounts[1][i][j] + alleleCounts[2][i][j])))) {
							output2DenovoMutationScores[i] = 100;
							break;
						}
					}

					if (output2DenovoMutationScores[i] > 0) {
						output3DenovoMutationNote = "";
						for (int j = 0; j < SAMPLE_SUFFIX.length; j++) {
							if (alleleCounts[j][i][INDEX_OF_N] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
								if (output3DenovoMutationNote.equals("")) {
									output3DenovoMutationNote = SAMPLE_SUFFIX[j];
								} else {
									output3DenovoMutationNote += "," + SAMPLE_SUFFIX[j];
								}
							}
						}
						if (! output3DenovoMutationNote.equals("")) {
							output2DenovoMutationScores[i] = (byte) (output2DenovoMutationScores[i] * DISCOUNT_FOR_N);
							if (output3DenovoMutationNotes[i] == null) {
								output3DenovoMutationNotes[i] = output3DenovoMutationNote + " missing values";
							} else {
								output3DenovoMutationNotes[i] += (output3DenovoMutationNote + " missing values");
							}
						}

						isThreeAlleles = isThreeAlleles(alleleCounts, i, MAX_ALLELE_COUNT_TREATED_AS_ZERO);
						output3DenovoMutationNote = "";
						for (int l = 0; l < SAMPLE_SUFFIX.length; l++) {
							if(isThreeAlleles[l]) {
								output2DenovoMutationScores[i] = 25;
								if (output3DenovoMutationNote.equals("")) {
									output3DenovoMutationNote = SAMPLE_SUFFIX[l];
								} else {
									output3DenovoMutationNote += "," + SAMPLE_SUFFIX[l];
								}
							}
						}
						if (! output3DenovoMutationNote.equals("")) {
							output2DenovoMutationScores[i] = (byte) (output2DenovoMutationScores[i] * .25);
							if (output3DenovoMutationNotes[i] == null) {
								output3DenovoMutationNotes[i] = output3DenovoMutationNote + ":3 alleles";
							} else {
								output3DenovoMutationNotes[i] += (output3DenovoMutationNote + ":3 alleles");
							}
						} else {
							loop = Math.min(NEIGHBOR, i);
							for (int l = 1; l <= loop; l++) {
								if (output2DenovoMutationScores[i - l] > 0) {
									output2DenovoMutationScores[i - l] = (byte) (output2DenovoMutationScores[i - l] * DISCOUNT_FOR_NEIGHBOR_MUTATION_OR_INDEL);
									if (output3DenovoMutationNotes[i - l] == null) {
										output3DenovoMutationNotes[i - l] = "DeNovo mutation nearby";
									} else if (output3DenovoMutationNotes[i - l].contains("DeNovo mutation nearby")) {
									} else {
										output3DenovoMutationNotes[i - l] += "; DeNovo mutation nearby";
									}
								}
							}
						}
					}
				}

				indicesOfInDels = new byte[] {INDEX_OF_INS, INDEX_OF_DEL};
				for (int j = 0; j < indicesOfInDels.length; j++) {
					output3DenovoMutationNote = ""; 
					for (int m = 0; m < alleleCounts.length; m++) {
						if (alleleCounts[m][i][indicesOfInDels[j]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
							if (output3DenovoMutationNote.equals("")) {
								output3DenovoMutationNote = SAMPLE_SUFFIX[m];
							} else {
								output3DenovoMutationNote += "," + SAMPLE_SUFFIX[m];
							}
						}

					}
					if (! output3DenovoMutationNote.equals("")) {
						loop = Math.min(NEIGHBOR, i);
						for (int l = 1; l <= loop; l++) {
							if (output2DenovoMutationScores[i - l] > 0) {
								output2DenovoMutationScores[i - l] = (byte) (output2DenovoMutationScores[i - l] * DISCOUNT_FOR_NEIGHBOR_MUTATION_OR_INDEL);
								if (output3DenovoMutationNotes[i - l] == null) {
									output3DenovoMutationNotes[i - l] = output3DenovoMutationNote + ":" + BASES_WITH_N_INS_DEL[indicesOfInDels[j]] + " nearby";
								} else if (output3DenovoMutationNotes[i - l].contains(BASES_WITH_N_INS_DEL[indicesOfInDels[j]] + " nearby")) {
								} else {
									output3DenovoMutationNotes[i - l] += "; " + output3DenovoMutationNote + ":" + BASES_WITH_N_INS_DEL[indicesOfInDels[j]] + " nearby";
								}
							}
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
		int[] depths;
		String[][] phredScoreProportions;
		String forwardGenotypes;

		phredScoreProportions = new String[SAMPLE_SUFFIX.length][BASES.length];
		depths = new int[SAMPLE_SUFFIX.length];
		for (int i = 0; i < phredScoreProportions.length; i++) {
			depths[i] = alleleCounts[i][indexOfCurrentMarker][0] + alleleCounts[i][indexOfCurrentMarker][1] + alleleCounts[i][indexOfCurrentMarker][2] + alleleCounts[i][indexOfCurrentMarker][3];
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
				+ "\t" + depths[0]
				+ "\t" + depths[1]
				+ "\t" + depths[2]
				+ "\t(" + phredScoreProportions[0][0] + "," + phredScoreProportions[0][1] + "," + phredScoreProportions[0][2] + "," + phredScoreProportions[0][3] + ") (" + phredScoreProportions[1][0] + "," + phredScoreProportions[1][1] + "," + phredScoreProportions[1][2] + "," + phredScoreProportions[1][3] + ") (" + phredScoreProportions[2][0] + "," + phredScoreProportions[2][1] + "," + phredScoreProportions[2][2] + "," + phredScoreProportions[2][3] + ")"
				+ "\t" + (depths[0] == 0? "" : ((mappingScores[0][indexOfCurrentMarker][0] + mappingScores[0][indexOfCurrentMarker][1] + mappingScores[0][indexOfCurrentMarker][2] + mappingScores[0][indexOfCurrentMarker][3]) / depths[0]))
				+ "\t" + (depths[1] == 0? "" : ((mappingScores[1][indexOfCurrentMarker][0] + mappingScores[1][indexOfCurrentMarker][1] + mappingScores[1][indexOfCurrentMarker][2] + mappingScores[1][indexOfCurrentMarker][3]) / depths[1]))
				+ "\t" + (depths[2] == 0? "" : ((mappingScores[2][indexOfCurrentMarker][0] + mappingScores[2][indexOfCurrentMarker][1] + mappingScores[2][indexOfCurrentMarker][2] + mappingScores[2][indexOfCurrentMarker][3]) / depths[2]))
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

	public static void saveToFile(String filename, int[][][] alleleCounts, int startPosition, boolean toIgnorePositionsWithZeroCount) {
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
					for (int k = 0; k < BASES_WITH_N_INS_DEL.length; k++) {
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

	public static String[][] getNameByTrios(String[] filenames) {

		return null;
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
		boolean isToProcessRegion;
		boolean isToProcessGenome;
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
		bamFilenamesOfTheTrio = new String[] {"rrd_F10639C_L008.bam", "rrd_F10639D_L007.bam", "rrd_F10639M_L007.bam"};
//		bamFilenamesOfTheTrio = new String[] {"rrd_F10008C_L005.bam", "rrd_F10008D_L004.bam", "rrd_F10008M_L003.bam"};
//		bamDir = "D:/bams/";
//		bamFilenamesTheTrio = new String[] {"F10639C_chr17_39346502_39346622.txt", "F10639D_chr17_39346502_39346622.txt", "F10639M_chr17_39346502_39346622.txt"};
		refFastaFilename = "/home/spectorl/xuz2/hg19_canonical.fa";
		scriptDir = "D:/logan/DeNovos/scripts/";
		outputDir = "/home/spectorl/xuz2/outputs/";
		isToAnnotate = false;
		isToProcessRegion = false;
		isToProcessGenome = true;
		numThreads = 1;
		chr = "17";
		start = 39346600;
		stop = 39346622;
		bedFilename = "/home/spectorl/xuz2/outputs/S04380219_Regions.bed";

		bamFilenamesForHelpMenu = (bamFilenamesOfTheTrio == null || bamFilenamesOfTheTrio[0] == null)? "" : bamFilenamesOfTheTrio[0];
		for (int i = 1; i < bamFilenamesOfTheTrio.length; i++) {
			bamFilenamesForHelpMenu += "," + bamFilenamesOfTheTrio[i];
		}

		commands = new String[] {"-annotation", "candidate=", "bamdir=", "scriptdir=", "outputdir=", "-processregion", "bamset=", "reffasta=", "chr=", "start=", "stop=", "-processgenome", "bed=", "numthreads="};
		String usage = "\n"
				+ "To annotate a list of candidate markers:"
				+ "   (1) command for annotatation (i.e. " + commands[0] + " (default))\n"
				+ "   (2) full path of the candidate list file (i.e. " + commands[1] + fileNameOfDeNovoPointMutationCandidateList + " (default))\n"
				+ "   (3) directory of the bam files (i.e. " + commands[2] + bamDir + " (default))\n"
				+ "   (4) directory for the output script file to extract information from bam files (i.e. " + commands[3] + scriptDir + " (default))\n"
				+ "To process a region:"
				+ "   (1) command for processing a region (i.e. " + commands[5] + " (default))\n"
				+ "   (2) directory of the bam files (i.e. " + commands[6] + bamFilenamesForHelpMenu + " (default))\n"
				+ "   (3) full path of the reference Fasta file (i.e. " + commands[7] + refFastaFilename + " (default))\n"
				+ "   (4) chr of the region (i.e. " + commands[8] + chr + " (default))\n"
				+ "   (5) start position of the region (i.e. " + commands[9] + start + " (default))\n"
				+ "   (6) stop position of the region (i.e. " + commands[10] + stop + " (default))\n"
				+ "   (7) directory for the output files (i.e. " + commands[4] + outputDir + " (default))\n"
				+ "To process genome:"
				+ "   (1) command for processing a region (i.e. " + commands[11] + " (default))\n"
				+ "   (2) directory of the bam files (i.e. " + commands[6] + bamFilenamesForHelpMenu + " (default))\n"
				+ "   (3) full path of the reference Fasta file (i.e. " + commands[7] + refFastaFilename + " (default))\n"
				+ "   (4) full path of the bed file (i.e. " + commands[12] + bedFilename + " (default))\n"
				+ "   (5) directory for the output files (i.e. " + commands[4] + outputDir + " (default))\n"
				+ "   (6) number of threads (i.e. " + commands[13] + numThreads + " (default))\n"
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
				isToProcessRegion = true;
				numArgs--;
			} else if (args[i].startsWith(commands[6])) {
				bamFilenamesOfTheTrio = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith(commands[7])) {
				refFastaFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[8])) {
				chr = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[9])) {
				start = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith(commands[10])) {
				stop = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith(commands[11])) {
				isToProcessGenome = true;
				numArgs--;
			} else if (args[i].startsWith(commands[12])) {
				bedFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[13])) {
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
		} else if (isToProcessRegion) {
			processRegion(bamDir, bamFilenamesOfTheTrio, refFastaFilename, chr, start, stop, outputDir, false);
		} else if (isToProcessGenome) {
			if (bamFilenamesOfTheTrio.length < 2 || bamFilenamesOfTheTrio == null) {
				processGenomeOfAllTriosInDir(bamDir, refFastaFilename, bedFilename, outputDir, numThreads, log);
			} else {
				processGenomeOfOneTrio(bamDir, bamFilenamesOfTheTrio, refFastaFilename, bedFilename, outputDir, numThreads, log);
			}
		}
	}

}
