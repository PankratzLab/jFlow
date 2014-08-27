package one;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TimeZone;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import bioinformatics.Samtools;
import bioinformatics.SeattleSeq;
import common.Array;
import common.Files;
import common.Logger;
import common.Matrix;
import common.Positions;
import common.Sort;
import common.ext;

public class SuperNovo {
	public static final String[] SAMPLE_SUFFIX = new String[] {"C", "D", "M"};
	public static final char[] BASES = new char[] {'A', 'T', 'G', 'C'};
	public static final char[] BASES_WITH_N = new char[] {'A', 'T', 'G', 'C', 'N'};
	public static final String[] BASES_WITH_N_INDEL = new String[] {"A", "T", "G", "C", "N", "Ins", "Del"};
	public static final String[] READS_COUNTS_ARRAY_STRUCT = new String[] {"A", "T", "G", "C", "N", "Ins", "Del", "numAllelesStrictCount", "numAllelesLooseCount"};
	public static final byte INDEX_OF_N = 4;
	public static final byte INDEX_OF_INS = 5;
	public static final byte INDEX_OF_DEL = 6;
	public static final byte INDEX_OF_NUM_ALLELES_STRICT = 7;
	public static final byte INDEX_OF_NUM_ALLELES_LOOSE = 8;
	public static final String[] DENOVO_MUTATION_ARRAY_STRUCT = new String[] {"i", "score", "C_readDepth", "D_readDepth", "M_readDepth", "C_allele1", "C_allele2", "D_allele1", "D_allele2", "M_allele1", "M_allele2"};
	public static final String[] DENOVO_MUTATION_NOTES_ARRAY_STRUCT = new String[] {"Ins", "Del", "Var", "3Allelic", "DeNovo"};
	public static final int MIN_AVERAGE_MAPPING_SCORE = 20;
	public static final int DEFAULT_PHRED_SCORE_FOR_DELETION = 30;
	public static final double THRESHOLD_PHRED_SCORE_FOR_INS_DEL = .10;
	public static final int MAX_ALLELE_COUNT_TREATED_AS_ZERO = 2;
	public static final int MIN_ALLELE_COUNT_FOR_DENOVO_MUTATION = 5;
	public static final double MIN_ALLELE_FREQ_FOR_DENOVO_MUTATION = .20;
	public static final int MIN_READ_DEPTH = 8;
	public static final int WINDOW_SIZE_FOR_NEARBY_INDEL = 60;
	public static final int WINDOW_SIZE_FOR_NEARBY_VARIANCE = 30;
	public static final int WINDOW_SIZE_FOR_NEARBY_DENOVO = 30;
	public static final double DISCOUNT_FOR_NEARBY_INDEL = .80;
	public static final double DISCOUNT_FOR_ON_INDEL_SITE = .50;
	public static final double DISCOUNT_FOR_NEARBY_VARIANCE = .95;
	public static final double DISCOUNT_FOR_ON_VARIANCE_SITE = .70;
	public static final double DISCOUNT_FOR_AT_VARIANCE = .60;
	public static final double DISCOUNT_FOR_N = .50;
	public static final double DISCOUNT_FOR_3ALLELES_LOOSE = .25;
	public static final double DISCOUNT_FOR_3ALLELES_STRICT = .50;
	public static final double DISCOUNT_FOR_NEARBY_DENOVOMUTATION = .90;
	public static final int REGION_LEN_AT_A_TIME_DEFAULT = 100000;
//	public static final String OUTPUT_FILE_HEADER = "id\tchr\tpos\tlookup\tsarver\tref\talt\tmendelianLikelihood\tmendelianPP\tmendelianGT\tsnpCode\tcode\tdeNovoLikelihood\tdeNovoPP\tactualDenovo\tconf\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tchildQuality\tdadQuality\tmomQuality";
//	public static final String OUTPUT_FILE_HEADER = "id\tchr\tpos\tlookup\tref\talt\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore\t1\t2\t4\t5\t7\t8";
	public static final String OUTPUT_FILE_HEADER = "id\tchr\tpos\tlookup\tref\talt\tcall\treadsCounts\tnumNearbyInDelVars\tposNearbyInDelVars\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore";
//	public static final String OUTPUT_FILE_HEADER2 = "id\tchr\tpos\tlookup\tref\talt\tcall\tcA\tcT\tcG\tcC\tcN\tcI\tcD\tdA\tdT\tdG\tdC\tdN\tdI\tdD\tmA\tmT\tmG\tmC\tmN\tmI\tmD\tcIns\tdIns\tmIns\tcDel\tdDel\tmDel\tcVar\tdVar\tmVar\tposNearbyInDelVars\tdeNovoGT\tflag\tcDepth\tdDepth\tmDepth\tcAPhred\tcTPhred\tcGPhred\tcCPhred\tdAPhred\tdTPhred\tdGPhred\tdCPhred\tmAPhred\tmTPhred\tmGPhred\tmCPhred\tcMap\tdMap\tmMap";
	public static final String OUTPUT_FILE_HEADER2 = "id\tchr\tpos\tlookup\tref\talt\tcall\t#Reads_A_Child\t#Reads_T_Child\t#Reads_G_Child\t#Reads_C_Child\t#Reads_N_Child\t#Reads_I_Child\t#Reads_D_Child\t#Reads_A_Dad\t#Reads_T_Dad\t#Reads_G_Dad\t#Reads_C_Dad\t#Reads_N_Dad\t#Reads_I_Dad\t#Reads_D_Dad\t#Reads_A_Mom\t#Reads_T_Mom\t#Reads_G_Mom\t#Reads_C_Mom\t#Reads_N_Mom\t#Reads_I_Mom\t#Reads_D_Mom\t#NearbySites_I_Child\t#NearbySites_I_Dad\t#NearbySites_I_Mom\t#NearbySites_D_Child\t#NearbySites_D_Dad\t#NearbySites_D_Mom\t#NearbySites_V_Child\t#NearbySites_V_Dad\t#NearbySites_V_Mom\t#NearbySites_3_Child\t#NearbySites_3_Dad\t#NearbySites_3_Mom\t#NearbySites_DNV\tposNearbyInDelVars\tdeNovoGT\tflag\tcDepth\tdDepth\tmDepth\tcAPhred\tcTPhred\tcGPhred\tcCPhred\tdAPhred\tdTPhred\tdGPhred\tdCPhred\tmAPhred\tmTPhred\tmGPhred\tmCPhred\tcMap\tdMap\tmMap";
	public static final String PARSE_RESULT_HEADER = "chr\tpos\tlookup\tnumHits\tid\tref\talt\tcall\tnote\tfwdGenotype\treadDepth_C\treadDepth_D\treadDepth_M\tPhred_A_Child\tPhred_T_Child\tPhred_G_Child\tPhred_C_Child\tPhred_A_Dad\tPhred_T_Dad\tPhred_G_Dad\tPhred_C_Dad\tPhred_A_Mom\tPhred_T_Mom\tPhred_G_Mom\tPhred_C_Mom\tMappingQuality_C\tMappingQuality_D\tMappingQuality_M\taltAllele%_C\taltAllele%_D\taltAllele%_M";
	public static final String PARSE_RESULT_HEADER2 = "chr\tpos\tlookup\tnumTrios";


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
		int[][] readsCounts;
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
				readsCounts = new int[3][];
				for (int i = 0; i < readsCounts.length; i++) {
					readsCounts[i] = getAlleleCounts(bamFileDir + line[0].substring(0, line[0].length()-1) + SAMPLE_SUFFIX[i] + "_" + line[1] + "_" + line[3] + ".txt", Integer.parseInt(line[1].split("chr")[1]), Integer.parseInt(line[3]), thresholdFor3Alleles);
				}
				score = getScoreForDeNovoMutation(readsCounts, line[9].charAt(0), line[8].charAt(0), thresholdFor3Alleles);
				notes = getNotesForDeNovoMutation(readsCounts, line[9].charAt(0), thresholdFor3Alleles);
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

	/**
	 * Annotate the markers with allele frequencies, allele counts, and deletion and insertion counts.
	 * @param readsCounts
	 * @param alternativeAllele
	 * @param thresholdFor3Alleles
	 * @return
	 */
	public static String getNotesForDeNovoMutation(int[][] readsCounts, char alternativeAllele, int thresholdFor3Alleles) {
		byte index;
		index = (byte) ext.indexOfChar(alternativeAllele, BASES);

		return	((readsCounts[1][index] > thresholdFor3Alleles || readsCounts[2][index] > thresholdFor3Alleles)? ext.formDeci(100 * readsCounts[0][index] / (double) readsCounts[0][5], 1) + "%," + readsCounts[1][index] + "," + readsCounts[2][index] + "; " : "")
			+	"C" + (readsCounts[0][6]>2? "(" + readsCounts[0][6]+ " alleles)" : "") + ":" + readsCounts[0][0] + "," + readsCounts[0][1] + "," + readsCounts[0][2] + "," + readsCounts[0][3] + ",(" + readsCounts[0][4] + "," + readsCounts[0][5] + ")"
			+ "; D" + (readsCounts[1][6]>2? "(" + readsCounts[1][6] + "alleles)" : "") + ":" + readsCounts[1][0] + "," + readsCounts[1][1] + "," + readsCounts[1][2] + "," + readsCounts[1][3] + ",(" + readsCounts[1][4] + "," + readsCounts[1][5] + ")"
			+ "; M" + (readsCounts[2][6]>2? "(" + readsCounts[2][6] + "alleles)" : "") + ":" + readsCounts[2][0] + "," + readsCounts[2][1] + "," + readsCounts[2][2] + "," + readsCounts[2][3] + ",(" + readsCounts[2][4] + "," + readsCounts[2][5] + ")";
	}

	/**
	 * A decimal value ranging 0 through 1 to indicate how likely the marker is a De Novo point mutation, with 1 being very likely and 0 very unlikely 
	 * @param readsCounts
	 * @param alternativeAllele
	 * @param referencedAllele
	 * @param thresholdFor3Alleles
	 * @return
	 */
	//This is functioning but naive. Has been replaced by getDenovoMarkerCandidateScores(...)
	public static int getScoreForDeNovoMutation(int[][] readsCounts, char alternativeAllele, char referencedAllele, int thresholdFor3Alleles) {
		int DenovoMutationScore;
		int index;

		index = ext.indexOfChar(alternativeAllele, BASES);
		if (readsCounts[1][index] > thresholdFor3Alleles || readsCounts[2][index] > thresholdFor3Alleles || readsCounts[0][8] > 2) {
			DenovoMutationScore = 0;
		} else {
			DenovoMutationScore = 1;
		}
		return	DenovoMutationScore;
	}

	public static void processGenomeOfAllTriosInDir(String bamDir, String refFastaFilename, String bedFilename, String outputDir, String scriptDir, String fullPathToTrioNameList, int regionLegnthATime, int numThreads, Logger log) {
		String[][] bamFilenamesByTrios;
		String command;
		String trioId;
		Vector<String> qsubFilesVec;
		PrintWriter writer;

		if (fullPathToTrioNameList == null) {
			bamFilenamesByTrios = groupNamesByTrios(Files.list(bamDir, ".bam", false));
		} else {
			bamFilenamesByTrios = loadNamesFromList(fullPathToTrioNameList);
		}
		qsubFilesVec = new Vector<String>(bamFilenamesByTrios.length);
		command = "cd " + outputDir + "\njcp one.SuperNovo bed=" + bedFilename + " outdir=" + outputDir + " bamdir=" + bamDir + " reffasta=" + refFastaFilename + " regionlengthatime=" + regionLegnthATime + " numthreads=" + numThreads;
		for (int i = 0; i < bamFilenamesByTrios.length; i++) {
//			processGenomeOfOneTrio(bamDir, bamFilenames, refFastaFilename, bedFilename, outputDir, numThreads, log);
			trioId = bamFilenamesByTrios[i][0];
			if (trioId == null || trioId.equals("")) {
				trioId = getRootOf(bamFilenamesByTrios[i], true);
			}
			Files.qsub(scriptDir + trioId + ".qsub", command + " trioid=" + trioId + " bamset=" + bamFilenamesByTrios[i][1] + "," + bamFilenamesByTrios[i][2] + "," + bamFilenamesByTrios[i][3], 2000, 36, 1);
			qsubFilesVec.add(scriptDir + trioId + ".qsub");
		}

		try {
			writer = new PrintWriter(new FileOutputStream(scriptDir + "master_run.sh"));
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

	public static void processGenomeOfOneTrio(String bamDir, String trioId, String[] bamFilenamesOfTheTrio, String refFastaFilename, String bedFilename, String outputDir, int regionLengthATime, int numThreads, Logger log) {
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
//		String trioId;
		String outFileName;
		ExecutorService executor = null;
		int numElementsATime;
		int startElement;
		long timer;
		SimpleDateFormat timeFormat;

		if (log == null) {
			log = new Logger();
		}

		if (regionLengthATime <= 0) {
			regionLengthATime = REGION_LEN_AT_A_TIME_DEFAULT;
		}
		timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
     	timer = new Date().getTime();
     	if (trioId == null) {
     		trioId = getRootOf(bamFilenamesOfTheTrio, false);
     	}
		outFileName = outputDir + trioId + "_superNovo_" + ext.rootOf(bedFilename) + ".txt";
		prevChr = "start";
		try {
			writer = new PrintWriter(outFileName);
			writer.println(OUTPUT_FILE_HEADER2);
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
				if (! line[0].equals("") && line[0] != null && line[0].toLowerCase().startsWith("chr")) {
					chr = line[0].substring(3).trim();
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
						loop = Math.min(start + regionLengthATime - 1, stop);
						if (numThreads > 1) {
							chrs.add(chr);
							starts.add(start);
							stops.add(loop);
						} else {
							processRegion(bamDir, bamFilenamesOfTheTrio, trioId, refFastaFilename, chr, start, loop, writer, null, log);
						}
						start += regionLengthATime;
					}
				} else {
					log.reportError("Warning: unrecognized chr number '" + line[0] + "' in bed file " + bedFilename);
				}
			}
			if (numThreads > 1) {
				if (chrs.size() < numThreads) {
					numThreads = chrs.size();
				}
				numElementsATime = (int) Math.ceil((double) chrs.size() / numThreads);
				startElement = 0;
				for (int i = 0; i < numThreads; i++) {
					if ((i + 1) == numThreads) {
						numElementsATime = chrs.size() % numElementsATime;
					}
					executor.execute(new WorkerThread(bamDir, bamFilenamesOfTheTrio, trioId, refFastaFilename, chrs, starts, stops, startElement, numElementsATime, writer, null, i, log));
					startElement += numElementsATime;
				}
				executor.shutdown();
				executor.awaitTermination(2, TimeUnit.DAYS);
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
		private String[] chrs;
		private int[] startPositions;
		private int[] stopPositions;
		private PrintWriter writer;
		private String outAlleleCountsFileName;
		private int threadId;
		private Logger log;

		public WorkerThread(String bamDir, String[] bamFilenames, String trioId, String refFastaFilename, Vector<String> chrs, Vector<Integer> starts, Vector<Integer> stops, int startIndexOfTheVectors, int numElementsOfTheVectors, PrintWriter writer, String outAlleleCountsFileName, int threadId, Logger log) {
			super();
			this.bamDir = bamDir;
			this.bamFilenames = bamFilenames;
			this.trioId = trioId;
			this.refFastaFilename = refFastaFilename;
			this.chrs = new String[numElementsOfTheVectors];
			this.startPositions = new int[numElementsOfTheVectors];
			this.stopPositions = new int[numElementsOfTheVectors];
			for (int i = 0; i < numElementsOfTheVectors; i++) {
				this.chrs[i] = chrs.elementAt(i + startIndexOfTheVectors);
				this.startPositions[i] = starts.elementAt(i + startIndexOfTheVectors);
				this.stopPositions[i] = stops.elementAt(i + startIndexOfTheVectors);
			}
			this.writer = writer;
			this.outAlleleCountsFileName = outAlleleCountsFileName;
			this.threadId = threadId;
			this.log = log;
		}

		@Override
		public void run() {
			for (int i = 0; i < this.chrs.length; i++) {
				processRegion(bamDir, bamFilenames, trioId, refFastaFilename, this.chrs[i], this.startPositions[i], this.stopPositions[i], writer, outAlleleCountsFileName, log);
			}
		}
		
		public int getThredId() {
			return this.threadId;
		}
	}

//	public static byte[][] processRegion(String dir, String bamFilename, byte chr, int start, int stop) {
//		PrintStream stream;
//		similar to piping or redirecting in linux
//		BufferedReader reader = new BufferedReader();
//		CmdLine.run("samtools view "+bamFilename+" chr"+chr+":"+start+"-"+stop, dir, stream);

	public static void processRegion(String bamDir, String[] bamFilenames, String refFastaFilename, String chr, int start, int stop, String outputDir, boolean isToOutputAlleleCounts, Logger log) {
		PrintWriter writer;
		String trioId;
		String outFileName;
		String outAlleleCountsFileName;
		long timer;
		SimpleDateFormat timeFormat;

		trioId = getRootOf(bamFilenames, false);
		outFileName = outputDir + trioId + "_denovoMutations_chr" + chr + "_" + start + "_" + stop + ".txt";
		if (isToOutputAlleleCounts) {
			outAlleleCountsFileName = outputDir + trioId + "_readsCounts_chr" + chr + "_" + start + "_" + stop + ".txt";
		} else {
			outAlleleCountsFileName = null;
		}
        timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
       	timer = new Date().getTime();
		try {
			writer = new PrintWriter(outFileName);
			writer.println(OUTPUT_FILE_HEADER2);
			processRegion(bamDir, bamFilenames, trioId, refFastaFilename, chr, start, stop, writer, outAlleleCountsFileName, log);
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		timer = new Date().getTime() - timer;
		System.out.println("processRegion result is ready at: " + outFileName + (isToOutputAlleleCounts? "\nAllele counts result is ready at:" + outAlleleCountsFileName : "") + "\nTotal time used " + timeFormat.format(timer));
	}

	public static void processRegion(String bamDir, String[] bamFilenames, String trioId, String refFastaFilename, String chr, int startPosition, int stopPosition, PrintWriter writer, String outAlleleCountsFileName, Logger log) {
		Process p;
//		ProcessBuilder ps;
		BufferedReader reader;
//		BufferedReader error;
		String line;
		Vector<String[]> bamContentVec;
		int numLines;
		int numMarkersPlusWindow;
		int window;
		int[][][] readsCounts = null;
		int[][][] phredScores = null;
		int[][][] mappingScores = null;
		String[] refAlleles;
		Vector<int[]> denovoMutations;
		Vector<Vector<Integer>[][]> denovoMutationNotes;
		int startPosAdjForWindow;
		int stopPosAdjForWindow;
		SimpleDateFormat timeFormat;
		long timer;
		String status;

        try {
        	window = Math.max(WINDOW_SIZE_FOR_NEARBY_INDEL, WINDOW_SIZE_FOR_NEARBY_VARIANCE);
        	startPosAdjForWindow = Math.max(0, startPosition - window);
        	stopPosAdjForWindow = Math.min(Positions.CHROMOSOME_LENGTHS_MAX[ext.indexOfStr(chr, Positions.CHR_CODES, false, true)], stopPosition + window);
    		numMarkersPlusWindow = stopPosAdjForWindow - startPosAdjForWindow + 1;
			readsCounts = new int[SAMPLE_SUFFIX.length][numMarkersPlusWindow][READS_COUNTS_ARRAY_STRUCT.length];
			phredScores = new int[SAMPLE_SUFFIX.length][numMarkersPlusWindow][BASES_WITH_N_INDEL.length];
			mappingScores = new int[SAMPLE_SUFFIX.length][numMarkersPlusWindow][BASES_WITH_N_INDEL.length];
	        timeFormat = new SimpleDateFormat("mm:ss.SSS");
			timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
			status = trioId + "\tchr" + chr + ":" + startPosAdjForWindow + "-" + stopPosAdjForWindow + "\t(" + startPosition + "-" + stopPosition + ")";

			timer = new Date().getTime();
			for (int i = 0; i < bamFilenames.length; i++) {
				bamContentVec = new Vector<String[]>();

				p = Runtime.getRuntime().exec("samtools view " + bamFilenames[i] + " chr" + chr + ":" + startPosAdjForWindow + "-" + stopPosAdjForWindow, null, new File(bamDir));

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
		        	cigarDecoding(bamContentVec.elementAt(j), chr, startPosAdjForWindow, stopPosAdjForWindow, 0, readsCounts[i], phredScores[i], mappingScores[i]);
				};
				status += ("\t" + numLines);

//		        error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
//		        while ((line = error.readLine()) != null) {
//		            System.out.println(line);
//		        }
//		        p.waitFor();
//				error.close();

//				displayErrorStream(p);
			}
			if (outAlleleCountsFileName != null) {
				saveAlleleCountsToFile(outAlleleCountsFileName, readsCounts, startPosAdjForWindow, true);
			}
			status += ("\t" + (timeFormat.format(new Date().getTime() - timer)));
			log.report(status);

			denovoMutations = new Vector<int[]>(10);
			filterByAlleleCount(readsCounts, mappingScores, startPosition - startPosAdjForWindow, stopPosition - startPosAdjForWindow, denovoMutations);
			denovoMutationNotes = new Vector<Vector<Integer>[][]>(denovoMutations.size());
			for (int j = 0; j < denovoMutations.size(); j++) {
				denovoMutationNotes.add(new Vector[DENOVO_MUTATION_NOTES_ARRAY_STRUCT.length][SAMPLE_SUFFIX.length]);
			}
			filterByNearbyInDelVars(readsCounts, denovoMutations, denovoMutationNotes);
			filterByNearbyDeNovos(readsCounts, denovoMutations, denovoMutationNotes);
			refAlleles = getRefFromFasta(refFastaFilename, denovoMutations, chr, startPosAdjForWindow);
			exportResult(writer, trioId, chr, startPosAdjForWindow, startPosition, refAlleles, readsCounts, phredScores, mappingScores, denovoMutations, denovoMutationNotes, "Phred score proportions < " + THRESHOLD_PHRED_SCORE_FOR_INS_DEL, (byte) 2);
			writer.flush();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(-1);
        } catch (InterruptedException e) {
			e.printStackTrace();
		}
        
	}

	public static void cigarDecoding(String[] aSingleLineOfBamFile, String chr, int startPos, int stopPos, int thresholdFor3Alleles, int[][] output1_ReadsCounts, int[][] output2_PhredScores, int[][] output3_MappingScores) {
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
			outputArrayPointer = currentPosition - startPos;
			outputArrayLength = stopPos - startPos + 1;

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
								output1_ReadsCounts[outputArrayPointer][indexInBases] ++;
								output2_PhredScores[outputArrayPointer][indexInBases] += convertToPhredScore(aSingleLineOfBamFile[10].charAt(readPointer));
								output3_MappingScores[outputArrayPointer][indexInBases] += currentMappingScore;
							} else {
								System.err.println("Error - unrecognized base (" + aSingleLineOfBamFile[9].charAt(readPointer) + ") at chr" + chr + ":" + startPos + outputArrayPointer + " Read ID: " + aSingleLineOfBamFile[0]);
							}
							readPointer ++;
							outputArrayPointer ++;
						}

						currentPosition += lengthOfCurrentSegment;
					}

				} else if (readSegments[0][i].equals("I")) {	// present in read sequence, but NOT in reference sequence.
					if (outputArrayPointer >= 0) {
						output1_ReadsCounts[outputArrayPointer][INDEX_OF_INS] ++;
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
							lengthOfCurrentSegment += outputArrayPointer;
							outputArrayPointer = 0;
						}

						loop = outputArrayPointer + Math.min(lengthOfCurrentSegment, outputArrayLength - outputArrayPointer);
						while (outputArrayPointer < loop) {
							output1_ReadsCounts[outputArrayPointer][INDEX_OF_DEL] ++;
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

	public static void filterByAlleleCount(int[][][] readsCounts, int[][][] mappingScores, int inputArrayStartIndex, int inputArrayStopIndex, Vector<int[]> output1DenovoMutations) {
		int[][] orderedIndices;
		int minRead;
		int[] temp;

		orderedIndices = new int[readsCounts.length][];
		minRead = MIN_READ_DEPTH / 2;
		for (int i = inputArrayStartIndex; i <= inputArrayStopIndex; i++) {
			for (int j = 0; j < orderedIndices.length; j++) {
				orderedIndices[j] = Sort.quicksort(new int[] {readsCounts[j][i][0], readsCounts[j][i][1], readsCounts[j][i][2], readsCounts[j][i][3]}, Sort.DESCENDING);
			}
			updateNumAlleles(readsCounts, orderedIndices, i);

			if (isReadDepth(readsCounts, orderedIndices, i, minRead)
					&& isAlleleCounts(readsCounts, orderedIndices, i)
//					&& is3Allelic(readsCounts, orderedIndices, i)
//					&& isMappingQuality(mappingScores, readsCounts, orderedIndices, i)
//					&& isInDel(readsCounts, orderedIndices, i)
//					&& isMismatch(readsCounts, orderedIndices, i)
					) {

				temp = new int[DENOVO_MUTATION_ARRAY_STRUCT.length];
				temp[0] = i;
				if ((readsCounts[1][i][orderedIndices[0][0]] == 0 && readsCounts[2][i][orderedIndices[0][0]] == 0) || (readsCounts[1][i][orderedIndices[0][1]] == 0 && readsCounts[2][i][orderedIndices[0][1]] == 0)) {
					temp[1] = 100;
				} else {
					temp[1] = 80;
				}
				for (int j = 0; j < readsCounts.length; j++) {
					for (int k = 0; k < 7; k++) {
						temp[j + 2] += readsCounts[j][i][k];
					}
					temp[2 * j + 5] = orderedIndices[j][0];
					if (readsCounts[j][i][orderedIndices[j][1]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
						temp[2 * j + 6] = orderedIndices[j][1];
					} else {
						temp[2 * j + 6] = -1;
					}
				}
				output1DenovoMutations.add(temp);
			}
		}
	}

	public static boolean isReadDepth(int[][][] readsCounts, int[][] orderedIndices, int i, int minRead) {
		return (   readsCounts[0][i][orderedIndices[0][0]] >= minRead
				&& readsCounts[1][i][orderedIndices[1][0]] >= minRead
				&& readsCounts[2][i][orderedIndices[2][0]] >= minRead
				&& (readsCounts[0][i][orderedIndices[0][0]] + readsCounts[0][i][orderedIndices[0][1]]) >= MIN_READ_DEPTH
				&& (readsCounts[1][i][orderedIndices[1][0]] + readsCounts[1][i][orderedIndices[1][1]]) >= MIN_READ_DEPTH
				&& (readsCounts[2][i][orderedIndices[2][0]] + readsCounts[2][i][orderedIndices[2][1]]) >= MIN_READ_DEPTH
//				&& readsCounts[0][i][INDEX_OF_TOTAL_READS] >= MIN_READ_DEPTH
//				&& readsCounts[1][i][INDEX_OF_TOTAL_READS] >= MIN_READ_DEPTH
//				&& readsCounts[2][i][INDEX_OF_TOTAL_READS] >= MIN_READ_DEPTH
				);
	}

	public static boolean isAlleleCounts(int[][][] readsCounts, int[][] orderedIndices, int i) {
		return (	(readsCounts[0][i][orderedIndices[0][1]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
						&& readsCounts[1][i][orderedIndices[0][1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
						&& readsCounts[2][i][orderedIndices[0][1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
//						&& readsCounts[0][i][orderedIndices[0][1]] > .18 * readsCounts[0][i][orderedIndices[0][0]]
//						&& readsCounts[0][i][orderedIndices[0][1]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 * (readsCounts[1][i][orderedIndices[0][1]] + readsCounts[2][i][orderedIndices[0][1]]))
						&& readsCounts[0][i][orderedIndices[0][1]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + readsCounts[1][i][orderedIndices[0][1]] + readsCounts[2][i][orderedIndices[0][1]])
//						&& (readsCounts[1][i][orderedIndices[0][1]] * 19) < (readsCounts[1][i][orderedIndices[0][0]] + readsCounts[1][i][orderedIndices[0][2]] + readsCounts[1][i][orderedIndices[0][3]])
//						&& (readsCounts[2][i][orderedIndices[0][1]] * 19) < (readsCounts[2][i][orderedIndices[0][0]] + readsCounts[2][i][orderedIndices[0][2]] + readsCounts[2][i][orderedIndices[0][3]])
				) || (readsCounts[0][i][orderedIndices[0][0]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
						&& readsCounts[1][i][orderedIndices[0][0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
						&& readsCounts[2][i][orderedIndices[0][0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
//						&& readsCounts[0][i][orderedIndices[0][0]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 * (readsCounts[1][i][orderedIndices[0][0]] + readsCounts[2][i][orderedIndices[0][0]]))
						&& readsCounts[0][i][orderedIndices[0][0]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + readsCounts[1][i][orderedIndices[0][0]] + readsCounts[2][i][orderedIndices[0][0]])
//						&& (readsCounts[1][i][orderedIndices[0][0]] * 19) < (readsCounts[1][i][orderedIndices[0][1]] + readsCounts[1][i][orderedIndices[0][2]] + readsCounts[1][i][orderedIndices[0][3]])
//						&& (readsCounts[2][i][orderedIndices[0][0]] * 19) < (readsCounts[2][i][orderedIndices[0][1]] + readsCounts[2][i][orderedIndices[0][2]] + readsCounts[2][i][orderedIndices[0][3]])
				) || (readsCounts[0][i][INDEX_OF_INS] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
						&& readsCounts[1][i][INDEX_OF_INS] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
						&& readsCounts[2][i][INDEX_OF_INS] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
//						&& readsCounts[0][i][INDEX_OF_INS] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 * (readsCounts[1][i][INDEX_OF_INS] + readsCounts[2][i][INDEX_OF_INS]))
						&& readsCounts[0][i][INDEX_OF_INS] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + readsCounts[1][i][INDEX_OF_INS] + readsCounts[2][i][INDEX_OF_INS])
				) || (readsCounts[0][i][INDEX_OF_DEL] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
						&& readsCounts[1][i][INDEX_OF_DEL] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
						&& readsCounts[2][i][INDEX_OF_DEL] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
//						&& readsCounts[0][i][INDEX_OF_DEL] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 * (readsCounts[1][i][INDEX_OF_DEL] + readsCounts[2][i][INDEX_OF_DEL]))
						&& readsCounts[0][i][INDEX_OF_DEL] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + readsCounts[1][i][INDEX_OF_DEL] + readsCounts[2][i][INDEX_OF_DEL])
				));
	}

	public static boolean isMismatch(int[][][] readsCounts, int[][] orderedIndices, int i) {
		return (   (readsCounts[0][i][INDEX_OF_N] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO || (readsCounts[0][i][orderedIndices[0][0]] + readsCounts[0][i][orderedIndices[0][1]]) >= 20 * readsCounts[0][i][INDEX_OF_N])
				&& (readsCounts[1][i][INDEX_OF_N] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO || (readsCounts[1][i][orderedIndices[0][0]] + readsCounts[1][i][orderedIndices[0][1]]) >= 20 * readsCounts[1][i][INDEX_OF_N])
				&& (readsCounts[2][i][INDEX_OF_N] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO || (readsCounts[2][i][orderedIndices[0][0]] + readsCounts[2][i][orderedIndices[0][1]]) >= 20 * readsCounts[2][i][INDEX_OF_N]));
	}

	public static boolean isInDel(int[][][] readsCounts, int[][] orderedIndices, int i) {
		int[] numInsDelMismatch;
		
		numInsDelMismatch = new int[] {readsCounts[0][i][INDEX_OF_INS] + readsCounts[0][i][INDEX_OF_DEL], readsCounts[1][i][INDEX_OF_INS] + readsCounts[1][i][INDEX_OF_DEL], readsCounts[2][i][INDEX_OF_INS] + readsCounts[2][i][INDEX_OF_DEL]};
		return (   (numInsDelMismatch[0] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO || (readsCounts[0][i][orderedIndices[0][0]] + readsCounts[0][i][orderedIndices[0][1]]) >= 20 * numInsDelMismatch[0])
				&& (numInsDelMismatch[1] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO || (readsCounts[1][i][orderedIndices[0][0]] + readsCounts[1][i][orderedIndices[0][1]]) >= 20 * numInsDelMismatch[1])
				&& (numInsDelMismatch[2] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO || (readsCounts[2][i][orderedIndices[0][0]] + readsCounts[2][i][orderedIndices[0][1]]) >= 20 * numInsDelMismatch[2]));
	}

	public static boolean is3Allelic(int[][][] readsCounts, int[][] orderedIndices, int i) {
		return (   (readsCounts[0][i][orderedIndices[0][2]] + readsCounts[0][i][orderedIndices[0][3]]) <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
				&& (readsCounts[1][i][orderedIndices[1][2]] + readsCounts[1][i][orderedIndices[1][3]]) <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
				&& (readsCounts[2][i][orderedIndices[2][2]] + readsCounts[2][i][orderedIndices[2][3]]) <= MAX_ALLELE_COUNT_TREATED_AS_ZERO);
	}

	public static boolean isMappingQuality(int[][][] mappingScores, int[][][] readsCounts, int[][] orderedIndices, int i) {
		return (   mappingScores[0][i][orderedIndices[0][0]] >= (MIN_AVERAGE_MAPPING_SCORE * readsCounts[0][i][orderedIndices[0][0]])
				&& mappingScores[1][i][orderedIndices[1][0]] >= (MIN_AVERAGE_MAPPING_SCORE * readsCounts[0][i][orderedIndices[1][0]])
				&& mappingScores[2][i][orderedIndices[2][0]] >= (MIN_AVERAGE_MAPPING_SCORE * readsCounts[0][i][orderedIndices[2][0]])
				&& mappingScores[0][i][orderedIndices[0][1]] >= (MIN_AVERAGE_MAPPING_SCORE * readsCounts[0][i][orderedIndices[0][1]])
				&& mappingScores[1][i][orderedIndices[1][1]] >= (MIN_AVERAGE_MAPPING_SCORE * readsCounts[0][i][orderedIndices[1][1]])
				&& mappingScores[2][i][orderedIndices[2][1]] >= (MIN_AVERAGE_MAPPING_SCORE * readsCounts[0][i][orderedIndices[2][1]]));
	}

//	public static void adjDenovoMutationScoresForMismatches(int[][][] readsCounts, int inputArrayMarkerIndex, byte[] output1DenovoMutationScores, String[] output2DenovoMutationNotes) {
//		String note;
//
//		note = "";
//		for (int j = 0; j < SAMPLE_SUFFIX.length; j++) {
//			if (readsCounts[j][inputArrayMarkerIndex][INDEX_OF_N] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
//				if (note.equals("")) {
//					note = SAMPLE_SUFFIX[j];
//				} else {
//					note += "," + SAMPLE_SUFFIX[j];
//				}
//			}
//		}
//
//		if (! note.equals("")) {
//			output1DenovoMutationScores[inputArrayMarkerIndex] = (byte) (output1DenovoMutationScores[inputArrayMarkerIndex] * DISCOUNT_FOR_N);
//
//			if (output2DenovoMutationNotes[inputArrayMarkerIndex] == null) {
//				output2DenovoMutationNotes[inputArrayMarkerIndex] = note + " mismatch(es)";
//			} else {
//				output2DenovoMutationNotes[inputArrayMarkerIndex] += (";" + note + " mismatch(es)");
//			}
//		}
//	}

//	public static void adjDenovoMutationScoresForThreeAlleles(int[][][] readsCounts, int indexCurrentMarker, byte[] output1DenovoMutationScores, String[] output2DenovoMutationNotes) {
//		String note;
//		double discountRate;
//
//		note = "";
//		discountRate = 1.0;
//		updateNumAlleles(readsCounts, indexCurrentMarker);
//		for(int i = 0; i < SAMPLE_SUFFIX.length; i++) {
//			if(readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] > 2) {
//				if (note.equals("")) {
//					note = SAMPLE_SUFFIX[i];
//				} else {
//					note += "," + SAMPLE_SUFFIX[i];
//				}
//
//				if (readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] > 2 || readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] > 3) {
//					discountRate = DISCOUNT_FOR_3ALLELES_LOOSE;
//				} else {
//					discountRate = DISCOUNT_FOR_3ALLELES_STRICT;
//				}
//			}
//		}
//
//		if (! note.equals("")) {
//			output1DenovoMutationScores[indexCurrentMarker] *= discountRate;
//			if (output2DenovoMutationNotes[indexCurrentMarker] == null) {
//				output2DenovoMutationNotes[indexCurrentMarker] = note + ":3 alleles";
//			} else {
//				output2DenovoMutationNotes[indexCurrentMarker] += (note + ":3 alleles");
//			}
//		}
//	}

	// (#A=10, #T=3, #G=2, #C=1) is counted as 3 allelic, though #G<=2 and #C<=2
	public static void updateNumAlleles(int[][][] readsCounts, int[][] orderedIndices, int indexCurrentMarker) {
		int sum;
		for(int i = 0; i < SAMPLE_SUFFIX.length; i++) {
//			for (int j = orderedIndices[0].length - 1; j >= 0; j--) {
//				if (readsCounts[i][indexCurrentMarker][orderedIndices[i][j]] > 0) {
//					if (readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] == 0) {
//						if (readsCounts[i][indexCurrentMarker][orderedIndices[i][j]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
//							readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] = j + 1;
//							readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] = j + 1;
//							break;
//						} else {
//							readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] = j + 1;
//						}
//					} else {
//						if (readsCounts[i][indexCurrentMarker][orderedIndices[i][j]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO || (readsCounts[i][indexCurrentMarker][orderedIndices[i][j]] + readsCounts[i][indexCurrentMarker][orderedIndices[i][j+1]]) > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
//							readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] = j + 1;
//							break;
//						}
//					}
//				} 
//			}

			readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] = 4;
			readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] = 0;
			for (int j = 0; j < orderedIndices[0].length; j++) {
				if (readsCounts[i][indexCurrentMarker][orderedIndices[i][j]] == 0) {
					readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] = j;
					break;
				}
			}
			sum = 0;
			for (int j = readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] - 1; j >= 0; j--) {
				sum += readsCounts[i][indexCurrentMarker][orderedIndices[i][j]];
				if (sum > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
					readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] = j + 1;
					break;
				}
			}
		}
	}

	public static void filterByNearbyInDelVars(int[][][] readsCounts, Vector<int[]> denovoMutations, Vector<Vector<Integer>[][]> output2DenovoMutationNotes) {
		byte[] indicesOfInDels;
		double discountDifference;
		int index;
		int loop2;
		int loop1;
		Vector<Integer>[][] temp;
	
		discountDifference = DISCOUNT_FOR_NEARBY_INDEL - DISCOUNT_FOR_ON_INDEL_SITE;
		indicesOfInDels = new byte[] {INDEX_OF_INS, INDEX_OF_DEL, INDEX_OF_NUM_ALLELES_LOOSE};
		loop1 = denovoMutations.size();
		for (int i = 0; i < loop1; i++) {
			index = denovoMutations.elementAt(i)[0];
			temp = output2DenovoMutationNotes.elementAt(i);
			loop2 = Math.min(readsCounts[0].length, index + WINDOW_SIZE_FOR_NEARBY_INDEL);
			for (int j = Math.max(0, index - WINDOW_SIZE_FOR_NEARBY_INDEL); j < loop2; j++) {
				for (int k = 0; k < readsCounts.length; k++) {
					for (int l = 0; l < DENOVO_MUTATION_NOTES_ARRAY_STRUCT.length; l++) {
						if ((l < 2 && readsCounts[k][j][indicesOfInDels[l]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) || (l == 2 && readsCounts[k][j][indicesOfInDels[l]] == 2) || (l == 3 && readsCounts[k][j][INDEX_OF_NUM_ALLELES_LOOSE] > 2)) {
							if (temp[l][k] == null) {
								temp[l][k] = new Vector<Integer>();
							}
							temp[l][k].addElement(j - index);
						}
					}
				}
			}
		}
	}

	public static void filterByNearbyDeNovos(int[][][] readsCounts, Vector<int[]> denovoMutations, Vector<Vector<Integer>[][]> output2DenovoMutationNotes) {
		int index2;
		int loop1;
		Vector<Integer>[][] temp;
		int[] currentDenovoMuations;
		int[] orderedIndices;

		loop1 = denovoMutations.size();
		for (int i = 0; i < loop1; i++) {
			currentDenovoMuations = denovoMutations.elementAt(i);
			temp = output2DenovoMutationNotes.elementAt(i);
			for (int j = 0; j < loop1; j++) {
				if (i != j) {
					index2 = denovoMutations.elementAt(j)[0];
					if (Math.abs(currentDenovoMuations[0] - index2) <= WINDOW_SIZE_FOR_NEARBY_DENOVO && currentDenovoMuations[3] > 10 && currentDenovoMuations[4] > 10) {
						orderedIndices = Sort.quicksort(new int[] {readsCounts[0][index2][0], readsCounts[0][index2][1], readsCounts[0][index2][2], readsCounts[0][index2][3]}, Sort.DESCENDING);
						if ((readsCounts[0][index2][orderedIndices[0]] >= 5 && readsCounts[1][index2][orderedIndices[0]] <= 1 && readsCounts[2][index2][orderedIndices[0]] <= 1) || (readsCounts[0][index2][orderedIndices[1]] >= 5 && readsCounts[1][index2][orderedIndices[1]] <= 1 && readsCounts[2][index2][orderedIndices[1]] <= 1)) {
							if (temp[4][0] == null) {
								temp[4][0] = new Vector<Integer>();
							}
							temp[4][0].addElement(j - index2);
						}
					}
				}
			}
		}
	}

//	public static boolean[] isThreeAlleles(int[][][] readsCounts, int indexCurrentMarker, int thresholdFor3Alleles) {
//		boolean[] isThreeAlleles;
//		int DenovoMutationScore;
//
//		isThreeAlleles = new boolean[readsCounts.length];
//		for(int i = 0; i < readsCounts.length; i++) {
//			DenovoMutationScore = 0;
//			for (int j = 0; j < BASES.length; j++) {
//				if (readsCounts[i][indexCurrentMarker][j] > thresholdFor3Alleles) {
//					DenovoMutationScore ++;
//				}
//			}
//			if(DenovoMutationScore > 2) {
//				isThreeAlleles[i] = true;
//			}
//		}
//
//		return	isThreeAlleles;
//	}

	public static byte adjDeNovoMutationScore (byte oldScore, double discountRate) {
		oldScore *= discountRate;
		if (oldScore == 0) {
			oldScore = 1;
		}
		return oldScore;
	}

	public static double getFisherExcat (int[][] contingencyTable) {
		return	getFactorial(contingencyTable[0][0] + contingencyTable[0][1]);
	}

	public static int getFactorial (int a) {
		int result;

		if (a == 0) {
			result = 0;
		} else {
			result = 1;
			for (int i = 2; i <= a; i++) {
				result *= i;
			}
		}

		return	result;
	}

	public static String[] getRefFromFasta(String refFastaFilename, Vector<int[]> denovoMutations, String chr, int startPosAdjForWindow) {
		int loop;
		int[] temp;
		Process p;
		BufferedReader reader;
//		BufferedReader error;
		String line;
		int pos;
		String[] refs;

		loop = denovoMutations.size();
		refs = new String[loop];
		try {
			for (int i = 0; i < loop; i++) {
				temp = denovoMutations.elementAt(i);
				pos = temp[0] + startPosAdjForWindow;
				p = Runtime.getRuntime().exec("samtools faidx " + refFastaFilename + " chr" + chr + ":" + pos + "-" + pos, null, new File(ext.parseDirectoryOfFile(refFastaFilename)));
	
//				ps=new ProcessBuilder("samtools", "faidx", refFastaFilename, "chr" + chr + ":" + start + "-" + stop);
//				ps.redirectErrorStream(true);
//				ps.directory(new File(ext.parseDirectoryOfFile(refFastaFilename)));
//				p = ps.start();  
	
//				p.waitFor();
				reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
				line = reader.readLine();
				while ((line = reader.readLine()) != null) {
					refs[i] = line;
				}
				p.waitFor();
				reader.close();
	
//				error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
//				while ((line = error.readLine()) != null) {
//					System.out.println(line);
//				}
//				error.close();

				displayErrorStream(p);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		return refs;
	}

	public static void exportResult(PrintWriter writer, String id, String chr, int startPosAdjForWindow, int startPos, String[] refsForDeNovoMutations, int[][][] readsCounts, int[][][] phredScores, int[][][] mappingScores, Vector<int[]> denovoMutations, Vector<Vector<Integer>[][]> denovoMutationNotes, String flag, byte format) {
		int[] temp;
		int loop;
		int l;

		loop = denovoMutations.size();
		for (int i = 0; i < loop; i++) {
			temp = denovoMutations.elementAt(i);
			l = temp[0];
			writer.println(id
					+ "\t" + chr
					+ "\t" + (startPosAdjForWindow + l)
					+ "\tchr" + chr + ":" + (startPosAdjForWindow + l)
					+ "\t" + refsForDeNovoMutations[i]
					+ "\t" + formatAltAllele(temp, refsForDeNovoMutations[i])
					+ "\t" + ext.formDeci(temp[1] / (double)100, 2)
					+ "\t" + formatAlleleCounts(readsCounts, l, format)
					+ "\t" + formatNotes(denovoMutationNotes, i, format)
					+ "\t" + formatForwardGenotypes(temp, refsForDeNovoMutations[i])
					+ "\t" + flag
					+ "\t" + temp[2]
					+ "\t" + temp[3]
					+ "\t" + temp[4]
					+ "\t" + formatPhred(phredScores, l, format)
					+ "\t" + formatMapping(mappingScores, temp, l, format) 
					);
		}
//		writer.flush();
	}

	public static String formatAlleleCounts(int[][][] readsCounts, int indexOfCurrentMarker, byte format) {
		String result;
		
		result = "";
		if (format == 2) {
			for (int i = 0; i < readsCounts.length; i++) {
				for (int j = 0; j < BASES_WITH_N_INDEL.length; j++) {
					result += ((i==0 && j == 0? "" : "\t") + (readsCounts[i][indexOfCurrentMarker][j]==0? "" : readsCounts[i][indexOfCurrentMarker][j]));
				}
			}
		} else {
			for (int i = 0; i < readsCounts.length; i++) {
				result += (i==0? "(": " (");
				for (int j = 0; j < BASES.length; j++) {
					result += ((j == 0? "" : ",") + (readsCounts[i][indexOfCurrentMarker][j]==0? "-" : readsCounts[i][indexOfCurrentMarker][j]));
				}
				result += "/";
				for (int j = BASES.length; j < BASES_WITH_N_INDEL.length; j++) {
					result += ((j == BASES.length? "" : ",") + (readsCounts[i][indexOfCurrentMarker][j]==0? "-" : readsCounts[i][indexOfCurrentMarker][j]));
				}
				result += ")";
			}
		}

		return result;
//		return (readsCounts[0][indexOfCurrentMarker][0]==0? "-" : readsCounts[0][indexOfCurrentMarker][0]) + "," + (readsCounts[0][indexOfCurrentMarker][1]==0? "-" : readsCounts[0][indexOfCurrentMarker][1]) + "," + (readsCounts[0][indexOfCurrentMarker][2]==0? "-" : readsCounts[0][indexOfCurrentMarker][2]) + "," + (readsCounts[0][indexOfCurrentMarker][3]==0? "-" : readsCounts[0][indexOfCurrentMarker][3]) + "/" + (readsCounts[0][indexOfCurrentMarker][4]==0? "-" : readsCounts[0][indexOfCurrentMarker][4]) + "," + (readsCounts[0][indexOfCurrentMarker][5]==0? "-" : readsCounts[0][indexOfCurrentMarker][5]) + "," + (readsCounts[0][indexOfCurrentMarker][6]==0? "-" : readsCounts[0][indexOfCurrentMarker][6]) + ") (" + (readsCounts[1][indexOfCurrentMarker][0]==0? "-" : readsCounts[1][indexOfCurrentMarker][0]) + "," + (readsCounts[1][indexOfCurrentMarker][1]==0? "-" : readsCounts[1][indexOfCurrentMarker][1]) + "," + (readsCounts[1][indexOfCurrentMarker][2]==0? "-" : readsCounts[1][indexOfCurrentMarker][2]) + "," + (readsCounts[1][indexOfCurrentMarker][3]==0? "-" : readsCounts[1][indexOfCurrentMarker][3]) + "/" + (readsCounts[1][indexOfCurrentMarker][4]==0? "-" : readsCounts[1][indexOfCurrentMarker][4]) + "," + (readsCounts[1][indexOfCurrentMarker][5]==0? "-" : readsCounts[1][indexOfCurrentMarker][5]) + "," + (readsCounts[1][indexOfCurrentMarker][6]==0? "-" : readsCounts[1][indexOfCurrentMarker][6]) + ") (" + (readsCounts[2][indexOfCurrentMarker][0]==0? "-" : readsCounts[2][indexOfCurrentMarker][0]) + "," + (readsCounts[2][indexOfCurrentMarker][1]==0? "-" : readsCounts[2][indexOfCurrentMarker][1]) + "," + (readsCounts[2][indexOfCurrentMarker][2]==0? "-" : readsCounts[2][indexOfCurrentMarker][2]) + "," + (readsCounts[2][indexOfCurrentMarker][3]==0? "-" : readsCounts[2][indexOfCurrentMarker][3]) + "/" + (readsCounts[2][indexOfCurrentMarker][4]==0? "-" : readsCounts[2][indexOfCurrentMarker][4]) + "," + (readsCounts[2][indexOfCurrentMarker][5]==0? "-" : readsCounts[2][indexOfCurrentMarker][5]) + "," + (readsCounts[2][indexOfCurrentMarker][6]==0? "-" : readsCounts[2][indexOfCurrentMarker][6]) + ")";
	}

	public static String formatNotes(Vector<Vector<Integer>[][]> denovoMarkerNotes, int indexOfCurrentElement, byte format) {
		String result2;
		String result1;
		Vector<Integer>[][] temp;
		int loop;
		boolean isFirstInTheSection;
		boolean isFirstSection;
		
		result1 = "";
		result2 = "";
		temp = denovoMarkerNotes.elementAt(indexOfCurrentElement);
		isFirstSection = true;
		if (format == 2) {
			for (int i = 0; i < temp.length; i++) {
				isFirstInTheSection = true;
				for (int j = 0; j < temp[0].length; j++) {
					if (temp[i][j] != null) {
						if (isFirstInTheSection) {
							result2 += ((isFirstSection? "" : " ") + DENOVO_MUTATION_NOTES_ARRAY_STRUCT[i] + "(");
							isFirstInTheSection = false;
							isFirstSection = false;
						} else {
							result2 += " ";
						}
						result2 += (SAMPLE_SUFFIX[j] + ":" + temp[i][j].elementAt(0));
						loop = temp[i][j].size();
						for (int k = 1; k < loop; k++) {
							result2 += ("," + temp[i][j].elementAt(k));
						}
						result1 += ((i==0 && j==0? "" : "\t") + temp[i][j].size());
					} else {
						result1 += (i==0 && j==0? "" : (i == 4 && j > 0? "" : "\t"));
					}
				}
				if (! isFirstInTheSection) {
					result2 += ")";
				}
			}
		} else {
			for (int i = 0; i < temp.length; i++) {
				isFirstInTheSection = true;
				result1 += ((result1.equals("")? "": " ") + "(");
				for (int j = 0; j < temp[0].length; j++) {
					if (temp[i][j] != null) {
						if (isFirstInTheSection) {
							result2 += (DENOVO_MUTATION_NOTES_ARRAY_STRUCT[i] + "(");
							isFirstInTheSection = false;
						} else {
							result2 += " ";
						}
						result2 += (SAMPLE_SUFFIX[j] + ":" + temp[i][j].elementAt(0));
						loop = temp[i][j].size();
						for (int k = 1; k < loop; k++) {
							result2 += ("," + temp[i][j].elementAt(k));
						}
						result1 += ((j==0? "" : ",") + temp[i][j].size());
					} else {
						result1 += ((j==0? "" : ",") + (i == 4 && j > 0? "" : "-"));
					}
				}
				if (! isFirstInTheSection) {
					result2 += ")";
				}
				result1 += ")";
			}
		}

		return result1 + "\t" + result2;
	}

	public static String formatAltAllele(int[] temp, String ref) {
		String result;

		if (temp[6] < 0) {
			result = "";
		} else {
			result = (BASES[temp[6]] + "");
			if (result.equalsIgnoreCase(ref)) {
				result = (BASES[temp[5]] + "");
			}
		}

		return result;
	}

	public static String formatForwardGenotypes(int[] temp, String ref) {
		String allele1;
		String result;

		result = "";
		for (int j = 0; j < SAMPLE_SUFFIX.length; j++) {
			result += (j == 0? "" : ",");
			allele1 = BASES[temp[5 + 2 * j]] + "";
			if (temp[6 + 2 * j] < 0) {
				result += (allele1 + allele1);
			} else {
				if (allele1.equalsIgnoreCase(ref)) {
					result += (allele1 + BASES[temp[6 + 2 * j]]);
				} else {
					result += (BASES[temp[6 + 2 * j]] + allele1);
				}
			}
		}

		return result;
	}

	public static String formatPhred(int[][][] phredScores, int indexOfCurrentMarker, byte format) {
		double[][] phredScoreProportions;
		double totalPhredScores;
		String result;

		phredScoreProportions = new double[SAMPLE_SUFFIX.length][BASES.length];
		for (int j = 0; j < phredScoreProportions.length; j++) {
			totalPhredScores = phredScores[j][indexOfCurrentMarker][0] + phredScores[j][indexOfCurrentMarker][1] + phredScores[j][indexOfCurrentMarker][2] + phredScores[j][indexOfCurrentMarker][3];
			for (int k = 0; k < phredScoreProportions[j].length; k++) {
				phredScoreProportions[j][k] = phredScores[j][indexOfCurrentMarker][k] / totalPhredScores;
			}
		}

		result = "";
		if (format == 2) {
			for (int i = 0; i < phredScoreProportions.length; i++) {
				for (int j = 0; j < phredScoreProportions[i].length; j++) {
					result += ((i==0 && j==0? "" : "\t") + (phredScoreProportions[i][j] == 0? "" : ext.formDeci(phredScoreProportions[i][j], 3)));
				}
			}
		} else {
			for (int i = 0; i < phredScoreProportions.length; i++) {
				result += ((i==0? "": " ") + "(");
				for (int j = 0; j < phredScoreProportions[i].length; j++) {
					result += ((j==0? "" : ",") + (phredScoreProportions[i][j] == 0? "-" : ext.formDeci(phredScoreProportions[i][j], 3)));
				}
				result += ")";
			}
		}

		return result;
	}

	public static String formatMapping(int[][][] mappingScores, int[] temp, int indexOfCurrentMarker, byte format) {
		int sum;
		String result;

		result = "";
		if (format == 2) {
			for (int i = 0; i < mappingScores.length; i++) {
				if (temp[i + 2] == 0) {
					result += (i==0? "" : "\t");
				} else {
					sum = 0;
					for (int j = 0; j < BASES.length; j++) {
						sum += mappingScores[i][indexOfCurrentMarker][j];
					}
					result += ((i==0? "" : "\t") + (sum/temp[i + 2]));
				}
			}
		} else {
			for (int i = 0; i < mappingScores.length; i++) {
				result += (i==0? "(" : ",");
				if (temp[i + 2] == 0) {
					result += "-";
				} else {
					sum = 0;
					for (int j = 0; j < BASES.length; j++) {
						sum += mappingScores[i][indexOfCurrentMarker][j];
					}
					result += (sum/temp[i + 2]);
				}
			}
			result += ")";
		}

		return result;
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

	public static void parseResults(String resultDir, double callScoreThreshold, Logger log) {
		String[] resultFilenames;
		String[] trioIds;
		BufferedReader reader;
		String[] line;
		HashMap<String, HashMap<String, String[]>> result;

		resultFilenames = Files.list(resultDir, ".txt", false);
		trioIds = new String[resultFilenames.length];
		result = new HashMap<String, HashMap<String, String[]>>();
		for (int i = 0; i < trioIds.length; i++) {
			trioIds[i] = resultFilenames[i].split("_")[0];
	        try {
	        	reader = new BufferedReader(new FileReader(resultDir + resultFilenames[i]));
	        	reader.readLine();
	        	while(reader.ready()) {
	        		line = reader.readLine().split("\t");
	        		if (Double.parseDouble(line[6]) >= callScoreThreshold) {
	        			if (! result.containsKey(line[3])) {
	        				result.put(line[3], new HashMap<String, String[]>());
	        			}
	        			result.get(line[3]).put(line[0], new String[] {line[7], line[8]});
	        		}
	        	}
				reader.close();
			} catch (Exception e) {
				System.out.println("Error reading from " + resultDir + resultFilenames[i]);
				System.out.println(e);
			}
		}

		saveParsedResults(resultDir + "SuperNovo_summary.txt", trioIds, result, (byte) 2, PARSE_RESULT_HEADER2, log);
	}

	public static void parseResults(String resultDir, double callScoreThreshold, byte outFormat, Logger log) {
		String[] resultFilenames;
		String[] trioIds;
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Object[] keys;
		String filename;
		HashMap<String, HashMap<String, String[]>> result;
		HashMap<String, String[]> temp;
		Set<Entry<String, String[]>> keySet;
		String key;
		String[] tmp1;
		String[] tmp2;
		int[][] readsCounts;
		int[] orderedIndices;
		int numTmp;
		double[] altProportion;
		boolean isThreeAllelesOrInDel;

		resultFilenames = Files.list(resultDir, ".txt", false);
		trioIds = new String[resultFilenames.length];
		result = new HashMap<String, HashMap<String, String[]>>();
		readsCounts = new int[3][7];
		for (int i = 0; i < trioIds.length; i++) {
			trioIds[i] = resultFilenames[i].split("_")[0];
	        try {
	        	reader = new BufferedReader(new FileReader(resultDir + resultFilenames[i]));
	        	reader.readLine();
	        	while(reader.ready()) {
	        		line = reader.readLine().split("\t");
	        		if (Double.parseDouble(line[6]) >= callScoreThreshold) {
	        			tmp1 = line[7].split(" ");
	        			for (int j = 0; j < 3; j++) {
		        			tmp2 = tmp1[j].substring(1).split("[,/()]");
	        				for (int k = 0; k < 7; k++) {
	        					if (tmp2[k].equals("-")) {
	        						readsCounts[j][k] = 0;
	        					} else {
	        						readsCounts[j][k] = Integer.parseInt(tmp2[k]);
	        					}
							}
						}

	        			orderedIndices = Sort.quicksort(new int[] {readsCounts[0][0], readsCounts[0][1], readsCounts[0][2], readsCounts[0][3]}, Sort.DESCENDING);
	        			numTmp = readsCounts[0][4] + readsCounts[0][5] + readsCounts[0][6];
	        			if (Integer.parseInt(line[14]) > 20 && Integer.parseInt(line[15]) > 20 && Integer.parseInt(line[16]) > 20
//	        					&& ((readsCounts[0][orderedIndices[1]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO && readsCounts[1][orderedIndices[1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && readsCounts[2][orderedIndices[1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && readsCounts[0][orderedIndices[1]] > .18 * readsCounts[0][orderedIndices[0]] && readsCounts[0][orderedIndices[1]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 * (readsCounts[1][orderedIndices[1]] + readsCounts[2][orderedIndices[1]])) && (readsCounts[1][orderedIndices[1]] * 19 ) < (readsCounts[1][orderedIndices[0]] + readsCounts[1][orderedIndices[2]] + readsCounts[1][orderedIndices[3]]) && (readsCounts[2][orderedIndices[1]] * 19 ) < (readsCounts[2][orderedIndices[0]] + readsCounts[2][orderedIndices[2]] + readsCounts[2][orderedIndices[3]]))
//	        							|| (readsCounts[0][orderedIndices[0]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO && readsCounts[1][orderedIndices[0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && readsCounts[2][orderedIndices[0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO &&  readsCounts[0][orderedIndices[0]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 * (readsCounts[1][orderedIndices[0]] + readsCounts[2][orderedIndices[0]]))) && (readsCounts[1][orderedIndices[0]] * 19 ) < (readsCounts[1][orderedIndices[1]] + readsCounts[1][orderedIndices[2]] + readsCounts[1][orderedIndices[3]]) && (readsCounts[2][orderedIndices[0]] * 19 ) < (readsCounts[2][orderedIndices[1]] + readsCounts[2][orderedIndices[2]] + readsCounts[2][orderedIndices[3]]))
	        					&& ((readsCounts[0][orderedIndices[1]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO && readsCounts[1][orderedIndices[1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && readsCounts[2][orderedIndices[1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && readsCounts[0][orderedIndices[1]] > .18 * readsCounts[0][orderedIndices[0]] && readsCounts[0][orderedIndices[1]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 * (readsCounts[1][orderedIndices[1]] + readsCounts[2][orderedIndices[1]])))
	        							|| (readsCounts[0][orderedIndices[0]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO && readsCounts[1][orderedIndices[0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && readsCounts[2][orderedIndices[0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO &&  readsCounts[0][orderedIndices[0]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 * (readsCounts[1][orderedIndices[0]] + readsCounts[2][orderedIndices[0]]))))
	        					&& (readsCounts[0][orderedIndices[2]] + readsCounts[0][orderedIndices[3]]) <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
	        					&& (numTmp <= MAX_ALLELE_COUNT_TREATED_AS_ZERO || (readsCounts[0][orderedIndices[0]] + readsCounts[0][orderedIndices[1]]) >= 20 * numTmp)) {

	        				isThreeAllelesOrInDel = false;
		        			altProportion = new double[3];
		        			for (int j = 0; j < altProportion.length; j++) {
			        			orderedIndices = Sort.quicksort(new int[] {readsCounts[j][0], readsCounts[j][1], readsCounts[j][2], readsCounts[j][3]}, Sort.DESCENDING);
			        			if ((readsCounts[j][orderedIndices[2]] + readsCounts[j][orderedIndices[3]]) > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
			        				isThreeAllelesOrInDel = true;
			        				break;
			        			}
			        			numTmp = readsCounts[j][4] + readsCounts[j][5] + readsCounts[j][6];
			        			if (numTmp > MAX_ALLELE_COUNT_TREATED_AS_ZERO && (readsCounts[j][orderedIndices[0]] + readsCounts[j][orderedIndices[1]]) < 20 * numTmp) {
			        				isThreeAllelesOrInDel = true;
			        				break;
			        			}

			        			numTmp = ext.indexOfChar(line[4].toUpperCase().charAt(0), BASES);
			        			if (orderedIndices[0] != numTmp) {
			        				altProportion[j] = readsCounts[j][orderedIndices[0]] / (double) (readsCounts[j][numTmp] + readsCounts[j][orderedIndices[0]]);
			        			} else {
			        				altProportion[j] = readsCounts[j][orderedIndices[1]] / (double) (readsCounts[j][numTmp] + readsCounts[j][orderedIndices[1]]);
			        			}
							}

		        			if (! isThreeAllelesOrInDel) {
			        			if (! result.containsKey(line[3])) {
			        				result.put(line[3], new HashMap<String, String[]>());
			        			}
			        			result.get(line[3]).put(line[0], new String[] {line[4], line[6], line[7], line[8], line[10], line[11], line[12], line[13], line[14], line[15], line[16], ext.formDeci(altProportion[0], 3), ext.formDeci(altProportion[1], 3), ext.formDeci(altProportion[2], 3)});
		        			}
	        			}
	        		}
	        	}
				reader.close();
			} catch (Exception e) {
				log.report("Error reading from " + resultDir + resultFilenames[i]);
				e.printStackTrace();
			}
		}

		filename = resultDir + "SuperNovo_summary.txt";
//		filename = "N:/statgen/OS_Logan/SuperNovo/SuperNovo_summary.txt";
		try {
			keys = result.keySet().toArray();
			Arrays.sort(keys);
			writer = new PrintWriter(new FileWriter(filename));
			if (outFormat == 1) {
				writer.println("id\tchr\tposition\tlookup\tnumHits\tref\talt\tcall\tnote\tfwdGenotype\treadDepth_C\treadDepth_D\treadDepth_M\tPhredScores\tMappingQuality_C\tMappingQuality_D\tMappingQuality_M\taltAllele%_C\taltAllele%_D\taltAllele%_M");
				for (int i = 0; i < keys.length; i++) {
					temp = result.get(keys[i]);
					keySet = temp.entrySet();
					numTmp = temp.size();
					for (Entry <String, String[]> entry : keySet) {
					    line = entry.getValue();
					    key = (String) keys[i];
						writer.println(entry.getKey() + "\t" + key.substring(3, key.indexOf(":")) + "\t" + key.substring(key.indexOf(":") + 1) + "\t" + key + "\t" + numTmp + "\t" + line[0] + "\t" + (line[3].substring(0, 1).equalsIgnoreCase(line[0])? line[3].charAt(1) : line[3].charAt(0)) + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\t" + line[9] + "\t" + line[10] + "\t" + line[11] + "\t" + line[12] + "\t" + line[13]);
					}
				}
			} else if (outFormat == 2) {
				writer.print("chr\tposition\tlookup\tnumTrios");
				for (int i = 0; i < trioIds.length; i++) {
					writer.print("\t" + trioIds[i]);
				}
				writer.println();
				for (int i = 0; i < keys.length; i++) {
					temp = result.get(keys[i]);
					line = ((String) keys[i]).split(":");
					writer.print(line[0] + "\t" + line[1] + "\t" + keys[i] + "\t" + temp.size());
					for (int j = 0; j < trioIds.length; j++) {
						if (temp.containsKey(trioIds[j])) {
							line = temp.get(trioIds[j]);
							writer.print("\t" + line[0] + "; " + line[1]);
						} else {
							writer.print("\t");
						}
					}
					writer.println();
				}
			}

			writer.close();
			log.report("Summerized result is available at: " + filename);

		} catch (Exception e) {
			log.reportError("Error writing to " + filename);
			e.printStackTrace();
		}
	}

	public static void parseResults(String resultDir, String annotationDir, String miniBamDir, String fullPathToTrioNameList, byte resultFormat, double callScoreThreshold, byte outFormat, Logger log) {
		String[] resultFilenames;
		String[] trioIds;
		BufferedReader reader;
		String[] line;
		HashMap<String, HashMap<String, String[]>> result;
		int[][][] readsCounts;
		int[][][] phredScores;
		int[][][] mappingScores;
		int[][] orderedIndices;
		int index;
		double[] altProportion;
		String[] annotation;
		Hashtable<String, String[]> annotationHash;
		HashSet<String> miniBamHash;
		Vector<String> annotationNeedsVars, annotationNeedsInDels, miniBamNeeds;
		String markerName;
		String[] resultElements;
		String miniBamLink;
		String[] seatleSeekHeader;
		
		annotationNeedsVars = new Vector<String>();
		annotationNeedsInDels =  new Vector<String>();
		miniBamNeeds = new Vector<String>();
		annotationHash = SeattleSeq.loadAllAnnotationInDir(annotationDir, log);
		miniBamHash = Samtools.listFiles(miniBamDir, log);
		resultFilenames = Files.list(resultDir, ".txt", false);
		trioIds = new String[resultFilenames.length];
		result = new HashMap<String, HashMap<String, String[]>>();
		orderedIndices = new int[SAMPLE_SUFFIX.length][];
		seatleSeekHeader = Matrix.extractColumn(SeattleSeq.RELEVANTS, 0);
		for (int i = 0; i < trioIds.length; i++) {
			trioIds[i] = resultFilenames[i].split("_")[0];
	        try {
	        	reader = new BufferedReader(new FileReader(resultDir + resultFilenames[i]));
	        	reader.readLine();
	        	while(reader.ready()) {
	        		line = reader.readLine().split("\t");
	        		if (Double.parseDouble(line[6]) >= callScoreThreshold) {
						readsCounts = new int[SAMPLE_SUFFIX.length][1][BASES_WITH_N_INDEL.length];
						phredScores = new int[SAMPLE_SUFFIX.length][1][BASES_WITH_N_INDEL.length];
						mappingScores = new int[SAMPLE_SUFFIX.length][1][BASES_WITH_N_INDEL.length];
	        			parseLine(line, resultFormat, readsCounts, phredScores, mappingScores);

	        			for (int j = 0; j < SAMPLE_SUFFIX.length; j++) {
	        				orderedIndices[j] = Sort.quicksort(new int[] {readsCounts[j][0][0], readsCounts[j][0][1], readsCounts[j][0][2], readsCounts[j][0][3]}, Sort.DESCENDING);
						}
	        			if (isMappingQuality(mappingScores, readsCounts, orderedIndices, 0)
	        					&& isAlleleCounts(readsCounts, orderedIndices, 0)
	        					&& is3Allelic(readsCounts, orderedIndices, 0)
//	        					&& isInDel(readsCounts, orderedIndices, 0)
	        					) {

		        			markerName = "chr" + line[1] + ":" + line[2] + "_" + line[4].toUpperCase() + "_" + line[5].toUpperCase();
		        			if (annotationHash.containsKey(markerName)) {
		        				annotation = annotationHash.get(markerName);
		        			} else {
		        				annotation = Array.stringArray(seatleSeekHeader.length);
		        				if (readsCounts[0][0][orderedIndices[0][1]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
		        					annotationNeedsVars.add(line[1] + "\t" + line[2] + "\t" + line[5] + "\t" + line[5]);	//chr + pos + alt + alt
		        				} else {
		        					annotationNeedsInDels.add(line[1] + "\t" + line[2] + "\t" + line[5] + "\t" + line[5]);
		        				}
		        			}

		        			if (miniBamHash.contains(line[0] + "\t" + line[1] + "\t" + line[2])) {	//trioId + chr + pos
//		        				miniBamLink = "=HYPERLINK(something)";
		        				miniBamLink = "exists";
		        			} else {
		        				miniBamLink = "missing";
		        				miniBamNeeds.add(line[0] + "\t" + line[1] + "\t" + line[2]);	//trioId + chr + pos
		        			}

		        			altProportion = new double[3];
		        			for (int j = 0; j < altProportion.length; j++) {
			        			index = ext.indexOfChar(line[4].toUpperCase().charAt(0), BASES);
			        			if (orderedIndices[j][0] == index || index < 0) {
			        				altProportion[j] = readsCounts[j][0][orderedIndices[j][1]] / (double) (readsCounts[j][0][orderedIndices[j][0]] + readsCounts[j][0][orderedIndices[j][1]]);
			        			} else {
			        				altProportion[j] = readsCounts[j][0][orderedIndices[j][0]] / (double) (readsCounts[j][0][index] + readsCounts[j][0][orderedIndices[j][0]]);
			        			}
							}

		        			if (! result.containsKey(line[3])) {
		        				result.put(line[3], new HashMap<String, String[]>());
		        			}

		        			resultElements = new String[] {
		        					line[4], //ref
		        					line[5], //alt
		        					line[6], //position
		        					line[41], //posNearbyInDelVars
		        					line[42], //deNovoGT
		        					line[44], //cDepth
		        					line[45], //dDepth
		        					line[46], //mDepth
		        					line[47], //cAPhred
		        					line[48], //cTPhred
		        					line[49], //cGPhred
		        					line[50], //cCPhred
		        					line[51], //dAPhred
		        					line[52], //dTPhred
		        					line[53], //dGPhred
		        					line[54], //dCPhred
		        					line[55], //mAPhred
		        					line[56], //mTPhred
		        					line[57], //mGPhred	
		        					line[58], //mCPhred
		        					line[59], //cMap
		        					line[60], //dMap
		        					line[61], //mMap
		        					ext.formDeci(altProportion[0], 3), 
		        					ext.formDeci(altProportion[1], 3), 
		        					ext.formDeci(altProportion[2], 3)
		        					};
		        			
		        			resultElements = Array.concatAll(resultElements, annotation);

		        			result.get(line[3]).put(line[0], resultElements);//	        			}
	        			}
	        		}
	        	}
				reader.close();
			} catch (Exception e) {
				log.report("Error reading from " + resultDir + resultFilenames[i]);
				e.printStackTrace();
			}
		}

		if (annotationNeedsVars.size() > 0) {
			Files.writeList(Array.toStringArray(annotationNeedsVars), annotationDir + "seattleSeq_input_" + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".txt");
		}
		if (annotationNeedsInDels.size() > 0) {
			Files.writeList(Array.toStringArray(annotationNeedsInDels), annotationDir + "seattleSeq_input_InDels_" + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".txt");
		}
		if (miniBamNeeds.size() > 0) {
			Files.writeList(Array.toStringArray(miniBamNeeds), miniBamDir + "miniBamNeeds_" + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".txt");
			Samtools.writerToFile(miniBamDir + "generateMiniBams" + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".sh", miniBamNeeds, loadNamesFromList(fullPathToTrioNameList), 50, log);
		}

		saveParsedResults(resultDir + "SuperNovo_summary.txt", trioIds, result, outFormat, PARSE_RESULT_HEADER + "\t" + Array.toStr(seatleSeekHeader), log);
	}

	public static void saveParsedResults(String resultsFilename, String[] trioIds, HashMap<String, HashMap<String, String[]>> result, byte format, String header, Logger log) {
		PrintWriter writer;
		String[] line;
		Object[] keys;
		HashMap<String, String[]> temp;
		Set<Entry<String, String[]>> keySet;
		String key;
		int numTmp;

		try {
			keys = result.keySet().toArray();
			Arrays.sort(keys);
			writer = new PrintWriter(new FileWriter(resultsFilename));
			if (format == 1) {
//				writer.println("id\tchr\tposition\tlookup\tnumHits\tref\talt\tcall\tnote\tfwdGenotype\treadDepth_C\treadDepth_D\treadDepth_M\tPhredScores\tMappingQuality_C\tMappingQuality_D\tMappingQuality_M\taltAllele%_C\taltAllele%_D\taltAllele%_M");
				writer.println(header);
				for (int i = 0; i < keys.length; i++) {
					temp = result.get(keys[i]);
					keySet = temp.entrySet();
					numTmp = temp.size();
					for (Entry <String, String[]> entry : keySet) {
					    line = entry.getValue();
					    key = (String) keys[i];
					    writer.print(key.substring(3, key.indexOf(":")) + "\t" + key.substring(key.indexOf(":") + 1) + "\t" + key + "\t" + numTmp + "\t" + entry.getKey());
					    for (int j = 0; j < line.length; j++) {
					    	writer.print("\t" + line[j]);
						}
						writer.println();
					}
				}
			} else if (format == 2) {
//				writer.print("chr\tposition\tlookup\tnumTrios");
				writer.print(header);
				for (int i = 0; i < trioIds.length; i++) {
					writer.print("\t" + trioIds[i]);
				}
				writer.println();
				for (int i = 0; i < keys.length; i++) {
					temp = result.get(keys[i]);
					line = ((String) keys[i]).split(":");
					writer.print(line[0] + "\t" + line[1] + "\t" + keys[i] + "\t" + temp.size());
					for (int j = 0; j < trioIds.length; j++) {
						if (temp.containsKey(trioIds[j])) {
							line = temp.get(trioIds[j]);
							writer.print("\t" + line[0] + "; " + line[1]);
						} else {
							writer.print("\t");
						}
					}
					writer.println();
				}
			}

			writer.close();
			log.report("Summerized result is available at: " + resultsFilename);

		} catch (Exception e) {
			log.reportError("Error writing to " + resultsFilename);
			e.printStackTrace();
		}
	}

	public static void parseLine(String[] line, byte resultFormat, int[][][] readsCounts, int[][][] phredScores, int[][][] mappingScores) {
		String[] tmp1;
		String[] tmp2;
		int index;
		int score;

		if (resultFormat == 1) {
			tmp1 = line[7].split(" ");
			for (int j = 0; j < 3; j++) {
				tmp2 = tmp1[j].substring(1).split("[,/()]");
				for (int k = 0; k < 7; k++) {
					if (tmp2[k].equals("-")) {
						readsCounts[j][0][k] = 0;
					} else {
						readsCounts[j][0][k] = Integer.parseInt(tmp2[k]);
					}
				}
			}
		} else {
			index = 7;
			for (int i = 0; i < SAMPLE_SUFFIX.length; i++) {
				for (int j = 0; j < BASES_WITH_N_INDEL.length; j++) {
					if (line[index] != null && ! line[index].isEmpty()) {
						readsCounts[i][0][j] = Integer.parseInt(line[index]);
					}
					index ++;
				}
			}

			index = 47;
			for (int i = 0; i < SAMPLE_SUFFIX.length; i++) {
				for (int j = 0; j < BASES_WITH_N_INDEL.length; j++) {
//					phredScores[i][0][j] = Double.parseDouble(line[index]);
					index ++;
				}
			}

			index = 59;
			for (int i = 0; i < SAMPLE_SUFFIX.length; i++) {
				if (line[index] != null && ! line[index].isEmpty()) {
					score = Integer.parseInt(line[index]);
					for (int j = 0; j < mappingScores[i][0].length; j++) {
						mappingScores[i][0][j] = score * readsCounts[i][0][j];
					}
				}
				index ++;
			}
		}
	}

	public static void saveAlleleCountsToFile(String filename, int[][][] readsCounts, int startPosition, boolean toIgnorePositionsWithZeroCount) {
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
			for (int i = 0; i < readsCounts[0].length; i++) {
				line = startPosition + "";
				total = 0;
				for (int j = 0; j < readsCounts.length; j++) {
					line += "\t(";
					for (int k = 0; k < BASES.length; k++) {
						line += ((k == 0? "" : ",") + (readsCounts[j][i][k] == 0? "-" : readsCounts[j][i][k]));
						total += readsCounts[j][i][k];
					}
					line += "/";
					for (int k = BASES.length; k < BASES_WITH_N_INDEL.length; k++) {
						line += ((k == 0? "" : ",") + (readsCounts[j][i][k] == 0? "-" : readsCounts[j][i][k]));
						total += readsCounts[j][i][k];
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
	public static String getRootOf(String[] filenames, boolean isToExcludeFirstColumn) {
		String[][] roots;
		String commonRoot = null;
		int maxLength;
		int index1 = 0;
		int index2 = 0;
		boolean found;
		int loop;
		String[] tmp;

		if (isToExcludeFirstColumn) {
			tmp = new String[filenames.length - 1];
			for (int i = 0; i < tmp.length; i++) {
				tmp[i] = filenames[i + 1];
			}
			filenames = tmp;
		}

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

	public static String[][] loadNamesFromList(String fullPathToTheList) {
		BufferedReader reader;
		String currentTrio;
		Vector<String[]> result;
		String[][] result2;
		int location;

		result = new Vector<String[]>();
		try {
			reader = new BufferedReader(new FileReader(fullPathToTheList));
			while(reader.ready()) {
				result.add(reader.readLine().split("\t"));
			}
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		result2 = new String[result.size()][4];
		for (int i = 0; i < result2.length; i++) {
			result2[i] = result.elementAt(i);
		}

		return result2;
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

	// Not yet finished.
	public static String[][] groupNamesByTrios2(String[] filenames) {
		boolean[] isVisited;
		String currentTrio;
		Vector<String[]> result;
		String[][] result2;
		int location;
		int[] indices;
		String[] nameSegments;
		String temp;
		int numMatches;

		isVisited = new boolean[filenames.length];
		indices = new int[3];
		for (int i = 0; i < indices.length; i++) {
			indices[i] = -1;
		}
		result = new Vector<String[]>(filenames.length/3);
		nameSegments = filenames[0].split("_");
		numMatches = 0;
		for (int i = 0; i < nameSegments.length; i++) {
			for (int j = 0; j < filenames.length; j++) {
				if (filenames[j].split("_")[i].equalsIgnoreCase(nameSegments[i])) {
					numMatches ++;
				}
			}
			if (numMatches == filenames.length) {
				
			}
		}
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

	public static void getScriptToExtractSmallerPiecesBamSamFiles(String outputfilenames) {

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String[] commands;
		String fileNameOfDeNovoPointMutationCandidateList;
		String[] bamFilenamesOfTheTrio;
		String bamDir;
		String scriptDir;
		String outputDir, annotationDir, miniBamDir;
		boolean isToAnnotate;
		int regionLegnthATime;
		int numThreads;
		String bamFilenamesForHelpMenu;
		String chr;
		int start;
		int stop;
		String refFastaFilename;
		String bedFilename;
		boolean isParseResult;
		String fullPathToTrioList;
		String trioId;
		Logger log;

		fileNameOfDeNovoPointMutationCandidateList = "N:/statgen/OS_Logan/IGV_validation/results_test1.txt";
		bamDir = "/home/spectorl/shared/exome_processing/bam/";
		bamFilenamesOfTheTrio = new String[] {"rrd_F10639C_L008.bam", "rrd_F10639D_L007.bam", "rrd_F10639M_L007.bam"};
		if (bamFilenamesOfTheTrio != null && !bamFilenamesOfTheTrio[0].equals("")) {
			trioId = getRootOf(bamFilenamesOfTheTrio, false);
		} else {
			trioId = null;
		}
//		bamDir = "D:/bams/";
//		bamFilenamesTheTrio = new String[] {"F10639C_chr17_39346502_39346622.txt", "F10639D_chr17_39346502_39346622.txt", "F10639M_chr17_39346502_39346622.txt"};
		fullPathToTrioList = "/home/spectorl/xuz2/lists/trios_list.txt";
		refFastaFilename = "/home/spectorl/xuz2/ref_fasta/hg19_canonical.fa";
		scriptDir = "/home/spectorl/xuz2/scripts/";
		outputDir = "/home/spectorl/xuz2/outputs/";
		isToAnnotate = false;
		regionLegnthATime = -1;
		numThreads = 1;
		chr = "";
		start = 39346600;
		stop = 39346622;
		bedFilename = "/home/spectorl/xuz2/outputs/WholeGenome.bed";
		annotationDir = "/home/spectorl/xuz2/outputs/SeattleSeqAnnotation/";
		miniBamDir = "/home/spectorl/xuz2/outputs/mini_bams/";
		isParseResult = false;

//		bamDir = "D:/logan/DeNovos/bams/";
//		bedFilename = "D:/logan/DeNovos/outputs/WholeGenome.bed";

//		isParseResult = true;
//		bedFilename = null;
//		outputDir = "N:/statgen/OS_Logan/SuperNovo/rawOutput/";
//		annotationDir = "N:/statgen/OS_Logan/SuperNovo/SeattleSeqAnnotation/";
//		miniBamDir = "D:/logan/DeNovos/mini_bams/";
//		fullPathToTrioList = "N:/statgen/OS_Logan/SuperNovo/triolist_rrd.txt";

		bamFilenamesForHelpMenu = (bamFilenamesOfTheTrio == null || bamFilenamesOfTheTrio[0] == null)? "" : bamFilenamesOfTheTrio[0];
		for (int i = 1; bamFilenamesOfTheTrio != null && i < bamFilenamesOfTheTrio.length && bamFilenamesOfTheTrio[i] != null; i++) {
			bamFilenamesForHelpMenu += "," + bamFilenamesOfTheTrio[i];
		}

		commands = new String[] {"-annotation", "candidate=", "bamdir=", "scriptdir=", "outdir=", "bamset=", "reffasta=", "chr=", "start=", "stop=", "bed=", "numthreads=", "-parseresult", "trioid=", "triolistfile=", "regionlengthatime="};
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
				+ "   (5) names of bam/sam files of the trio, if just for a trio (i.e. " + commands[5] + bamFilenamesForHelpMenu + " (default))\n"
				+ "   (6) label of the trio, to be used in name of the output file (i.e. " + commands[13] + trioId + " (default))\n"
				+ "   (7) full path of the reference Fasta file (i.e. " + commands[6] + refFastaFilename + " (default))\n"
				+ "   (8) directory for the output files (i.e. " + commands[4] + outputDir + " (default))\n"
				+ "To process genome:"
				+ "   (1) full path of the bed file (i.e. " + commands[10] + bedFilename + " (default))\n"
				+ "   (2) directory of the bam files (i.e. " + commands[2] + bamDir + " (default))\n"
				+ "   (3) names the bam/sam files of the trio, if just for a trio (i.e. " + commands[5] + bamFilenamesForHelpMenu + " (default))\n"
				+ "   (4) full path of the reference Fasta file (i.e. " + commands[7] + refFastaFilename + " (default))\n"
				+ "   (5) directory for the output files (i.e. " + commands[4] + outputDir + " (default))\n"
				+ "   (6) full path to the trio list file (i.e. " + commands[14] + fullPathToTrioList + " (default))\n"
				+ "   (7) number of threads (i.e. " + commands[11] + numThreads + " (default))\n"
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
				if (args[i].split("=").length < 2) {
					bamFilenamesOfTheTrio = null;
				} else {
					bamFilenamesOfTheTrio = args[i].split("=")[1].split(",");
				}
				numArgs--;
			} else if (args[i].startsWith(commands[6])) {
				refFastaFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[7])) {
				if (args[i].split("=").length < 2) {
					chr = null;
				} else {
					chr = args[i].split("=")[1];
				}
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
			} else if (args[i].startsWith(commands[12])) {
				isParseResult = true;
				numArgs--;
			} else if (args[i].startsWith(commands[13])) {
				trioId = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[14])) {
				fullPathToTrioList = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(commands[15])) {
				regionLegnthATime = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		log = new Logger(outputDir + "SuperNovo_" + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log");
		
		if (isToAnnotate) {
//			getAlleleCounts("D:/bams/F10639C_chr11_89531764.txt", 11, 89531764, 3);	//Testing only. Feature already included in screenDeNovoPointMutation()
//			getAlleleCounts("D:/bams/F10639C_chr5_140558628.txt", 5, 140558628, 3);
			generateScriptForSamtools(fileNameOfDeNovoPointMutationCandidateList, bamFilenamesOfTheTrio[0], scriptDir);
			screenDeNovoPointMutation(fileNameOfDeNovoPointMutationCandidateList, bamFilenamesOfTheTrio[0], 1);
		} else if (chr != null && !chr.equals("") && start > 0 && stop > 0 && stop >= start) {
			processRegion(bamDir, bamFilenamesOfTheTrio, refFastaFilename, chr, start, stop, outputDir, false, log);
		} else if (bedFilename != null && new File(bedFilename).exists()) {
			if (bamFilenamesOfTheTrio == null || bamFilenamesOfTheTrio.length < 2) {
				processGenomeOfAllTriosInDir(bamDir, refFastaFilename, bedFilename, outputDir, scriptDir, fullPathToTrioList, regionLegnthATime, numThreads, log);
			} else {
				processGenomeOfOneTrio(bamDir, trioId, bamFilenamesOfTheTrio, refFastaFilename, bedFilename, outputDir, regionLegnthATime, numThreads, log);
			}
		} else if (isParseResult) {
			parseResults(outputDir, annotationDir, miniBamDir, fullPathToTrioList, (byte) 2, 0, (byte) 1, log);
		} else {
			log.reportError("No command executed, due to:\n1) " + commands[0] + " is not specified; or 2) chr, start and stop are not complete; or 3) .bed file is not found");
		}
	}
}
