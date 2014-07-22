package one;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Vector;

import common.Array;
import common.CmdLine;
import common.Files;
import common.Sort;
import common.ext;

public class DeNovoSNP {
	public static final String[] SAMPLE_SUFFIX = new String[] {"C", "D", "M"};
	public static final char[] BASES = new char[] {'A', 'T', 'G', 'C'};
	public static final String[] BASES_AND_INDELS = new String[] {"A", "T", "G", "C", "Ins", "Del"};
	public static final byte MARKERDATA_NUMSAMPLES_START = 0;
	public static final byte MARKERDATA_NUMSAMPLES_LEN = 4;
	public static final int DEFAULT_PHRED_SCORE_FOR_DELETION = 30;
	public static final double THRESHOLD_PHRED_SCORE_FOR_INS_DEL = .10;
	public static final int MAX_ALLELE_COUNT_TREATED_AS_ZERO = 2;
	public static final int MIN_READ_DEPTH = 10;

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

		result = new int[] {0, 0, 0, 0, 0, 0, 0};	//A, T, G, C, Del, Total, #Alleles
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
								result[4] ++;
								break;
							} else {
								index += (positionOnThisRead - startPositionOfTheRead - index);
//								index = (positionOnThisRead - startPositionOfTheRead);
								index = ext.indexOfChar(line[9].charAt(index), BASES);
								if (index >= 0) {
									result[index] ++;
								} else {
									System.out.println("Error - unrecognized base (neither A, T, G, or C) in the following file and read:\n" + readsFileFullPath + "\nRead ID: " + line[0]);
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
			result[5] += result[i];
			if (result[i] > thresholdFor3Alleles) {
				result[6] ++;
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
		if (alleleCounts[1][index] > thresholdFor3Alleles || alleleCounts[2][index] > thresholdFor3Alleles || alleleCounts[0][6] > 2) {
			DenovoMutationScore = 0;
		} else {
			DenovoMutationScore = 1;
		}
		return	DenovoMutationScore;
	}

	public static boolean isInsertionDeletionOtherMutationsNearby() {
		return false;
	}

//	public static byte[][] processRegion(String dir, String bamFilename, byte chr, int start, int stop) {
//		PrintStream stream;
//		similar to piping or redirecting in linux
//		BufferedReader reader = new BufferedReader();
//		CmdLine.run("samtools view "+bamFilename+" chr"+chr+":"+start+"-"+stop, dir, stream);

	public static void processRegion(String bamDir, String[] bamFilenames, String refFastaFilename, byte chr, int start, int stop, String outputDir) {
		Process p;
//		ProcessBuilder ps;
		BufferedReader reader;
		BufferedReader error;
		String line;
		String outFileName;
		String trioId;
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
*/		PrintWriter writer;
		String refAlleles;
		byte[] denovoMarkerCandidateScores;
		String[] denovoMarkerNotes;

        try {
            numMarkers = stop - start + 1;
			alleleCounts = new int[SAMPLE_SUFFIX.length][numMarkers][BASES_AND_INDELS.length];
			phredScores = new int[SAMPLE_SUFFIX.length][numMarkers][BASES_AND_INDELS.length];
			mappingScores = new int[SAMPLE_SUFFIX.length][numMarkers][BASES_AND_INDELS.length];

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
				System.out.println(numLines);
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
			denovoMarkerCandidateScores = new byte[numMarkers];
			denovoMarkerNotes = new String[numMarkers];
			getDenovoMarkerCandidateScores(alleleCounts, phredScores, mappingScores, depths, denovoMarkerCandidateScores, denovoMarkerNotes);
			trioId = rootOf(bamFilenames);
			outFileName = outputDir + trioId + "_denovoSnp.txt";
/*			depthes = new int[NUMSAMPLES_IN_A_TRIO];
			phredScoresProportions = new double[NUMSAMPLES_IN_A_TRIO][BASES.length];
*/			writer = new PrintWriter(outFileName);
//			writer.println("id\tchr\tpos\tlookup\tsarver\tref\talt\tmendelianLikelihood\tmendelianPP\tmendelianGT\tsnpCode\tcode\tdeNovoLikelihood\tdeNovoPP\tactualDenovo\tconf\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tchildQuality\tdadQuality\tmomQuality");
			writer.println("id\tchr\tpos\tlookup\tref\talt\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore\t1\t2\t4\t5\t7\t8");
			for (int i = 0; i < numMarkers; i++) {
				if(denovoMarkerCandidateScores[i] > 0) {
					exportInfoForPosition(writer, trioId, chr, start + i, refAlleles.charAt(i), denovoMarkerCandidateScores[i], alleleCounts, phredScores, mappingScores, i, denovoMarkerNotes, "Phred score proportions < " + THRESHOLD_PHRED_SCORE_FOR_INS_DEL);
				}
			}

/*			for (int i = 0; i < numMarkers; i++) {
				for (int j = 0; j < depthes.length; j++) {
					depthes[j] = alleleCounts[j][i][0] + alleleCounts[j][i][1] + alleleCounts[j][i][2] + alleleCounts[j][i][3];
					totalPhredScores = phredScores[j][i][0] + phredScores[j][i][1] + phredScores[j][i][2] + phredScores[j][i][3];
					for (int k = 0; k < BASES.length; k++) {
						phredScoresProportions[j][k] = phredScores[j][i][k] / totalPhredScores;
					}
				}
				for (int j = 0; j < BASES.length; j++) {
//					if (alleleCounts[0][i][j] > MAX_ALLELE_COUNT_TREATED_AS_ZERO && alleleCounts[1][i][j] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && alleleCounts[2][i][j] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && phredScoresProportions[0][j] < THRESHOLD_PHRED_SCORE_FOR_INS_DEL) {
					if (alleleCounts[0][i][j] > MAX_ALLELE_COUNT_TREATED_AS_ZERO && alleleCounts[1][i][j] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && alleleCounts[2][i][j] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
						forwardGenotypes = getForwardGenotypes(alleleCounts, i);
						exportInfoForPosition(writer, trioId, chr, start + i, refAlleles.charAt(i) + "", BASES[j] + "", .75 + "",
								"(" + (alleleCounts[0][i][0]==0? "-" : alleleCounts[0][i][0]) + "," + (alleleCounts[0][i][1]==0? "-" : alleleCounts[0][i][1]) + "," + (alleleCounts[0][i][2]==0? "-" : alleleCounts[0][i][2]) + "," + (alleleCounts[0][i][3]==0? "-" : alleleCounts[0][i][3]) + "," + (alleleCounts[0][i][4]==0? "-" : alleleCounts[0][i][4]) + "," + (alleleCounts[0][i][5]==0? "-" : alleleCounts[0][i][5]) + ") (" + (alleleCounts[1][i][0]==0? "-" : alleleCounts[1][i][0]) + "," + (alleleCounts[1][i][1]==0? "-" : alleleCounts[1][i][1]) + "," + (alleleCounts[1][i][2]==0? "-" : alleleCounts[1][i][2]) + "," + (alleleCounts[1][i][3]==0? "-" : alleleCounts[1][i][3]) + "," + (alleleCounts[1][i][4]==0? "-" : alleleCounts[1][i][4]) + "," + (alleleCounts[1][i][5]==0? "-" : alleleCounts[1][i][5]) + ") (" + (alleleCounts[2][i][0]==0? "-" : alleleCounts[2][i][0]) + "," + (alleleCounts[2][i][1]==0? "-" : alleleCounts[2][i][1]) + "," + (alleleCounts[2][i][2]==0? "-" : alleleCounts[2][i][2]) + "," + (alleleCounts[2][i][3]==0? "-" : alleleCounts[2][i][3]) + "," + (alleleCounts[2][i][4]==0? "-" : alleleCounts[2][i][4]) + "," + (alleleCounts[2][i][5]==0? "-" : alleleCounts[2][i][5]) + ")",
								forwardGenotypes, "Phred score proportions < " + THRESHOLD_PHRED_SCORE_FOR_INS_DEL, depthes[0], depthes[1], depthes[2],
//								phredScores[0][i][j] / totalPhredScores[0], phredScores[1][i][j] / totalPhredScores[1], phredScores[2][i][j] / totalPhredScores[2],
								"(" + (phredScoresProportions[0][0]==0? "-" : ext.formDeci(phredScoresProportions[0][0], 3)) +"," + (phredScoresProportions[0][1]==0? "-" : ext.formDeci(phredScoresProportions[0][1], 3)) + "," + (phredScoresProportions[0][2]==0? "-" : ext.formDeci(phredScoresProportions[0][2], 3)) + "," + (phredScoresProportions[0][3]==0? "-" : ext.formDeci(phredScoresProportions[0][3], 3)) + ") (" + (phredScoresProportions[1][0]==0? "-" : ext.formDeci(phredScoresProportions[1][0], 3)) + "," + (phredScoresProportions[1][1]==0? "-" : ext.formDeci(phredScoresProportions[1][1], 3)) + "," + (phredScoresProportions[1][2]==0? "-" : ext.formDeci(phredScoresProportions[1][2], 3)) + "," + (phredScoresProportions[1][3]==0? "-" : ext.formDeci(phredScoresProportions[1][3], 3)) + ") (" + (phredScoresProportions[2][0]==0? "-" : ext.formDeci(phredScoresProportions[2][0], 3)) + "," + (phredScoresProportions[2][1]==0? "-" : ext.formDeci(phredScoresProportions[2][1], 3)) + "," + (phredScoresProportions[2][2]==0? "-" : ext.formDeci(phredScoresProportions[2][2], 3)) + "," + (phredScoresProportions[2][3]==0? "-" : ext.formDeci(phredScoresProportions[2][3], 3)) + ")",
								(mappingScores[0][i][0] + mappingScores[0][i][1] + mappingScores[0][i][2] + mappingScores[0][i][3]) / depthes[0], (mappingScores[1][i][0] + mappingScores[1][i][1] + mappingScores[1][i][2] + mappingScores[1][i][3]) / depthes[1], (mappingScores[2][i][0] + mappingScores[2][i][1] + mappingScores[2][i][2] + mappingScores[2][i][3]) / depthes[2],
								forwardGenotypes.charAt(0), forwardGenotypes.charAt(1), forwardGenotypes.charAt(3), forwardGenotypes.charAt(4), forwardGenotypes.charAt(6), forwardGenotypes.charAt(7));
						break;
					}
				}
			}
*/			writer.close();
			System.out.println("processRegion result is ready at: " + outFileName);

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
	public static void getAlleleCountsPhredScoresMappingScores(String[] aSingleLineOfBamFile, int chr, int start, int stop, int thresholdFor3Alleles, int[][] output1_AlleleCounts, int[][] output2_PhredScores, int[][] output3_MappingScores) {
//		char[] bases;
		int readPointer;
		int outputArrayPointer;
//		int positionOffset;
		int index2, positionPointer, lengthOfCurrentSegment, positionToLookFor;
		String[][] readSegments;
		int currentMappingScore;

		currentMappingScore = Integer.parseInt(aSingleLineOfBamFile[4]);
		positionToLookFor = start;
		readSegments = ext.getOperatorsOperatorIndicesAndSplit(aSingleLineOfBamFile[5], "DIMNPSH");
		positionPointer = Integer.parseInt(aSingleLineOfBamFile[3]);
		readPointer = 0;
//		positionOffset = 0;
		outputArrayPointer = Math.max(0, positionPointer - start);
		for (int i = 0; i < readSegments[0].length; i++) {
			lengthOfCurrentSegment = Integer.parseInt(readSegments[2][i]);
			if (readSegments[0][i].equals("I")) {
				if ((positionToLookFor - positionPointer) == 1) {
					output1_AlleleCounts[outputArrayPointer][4] ++;
					output2_PhredScores[outputArrayPointer][4] += convertToPhredScore(aSingleLineOfBamFile[10].charAt(readPointer));	//TODO
					output3_MappingScores[outputArrayPointer][4] += currentMappingScore;
				}
				readPointer += lengthOfCurrentSegment;
			} else if (readSegments[0][i].equals("S") || readSegments[0][i].equals("H") || readSegments[0][i].equals("N")) {
				readPointer += lengthOfCurrentSegment;
			} else if (readSegments[0][i].equals("D")) {
				if ((positionToLookFor - positionPointer) <= lengthOfCurrentSegment) {
					while (positionToLookFor <= stop && lengthOfCurrentSegment > 0) {
						output1_AlleleCounts[outputArrayPointer][5] ++;
						output2_PhredScores[outputArrayPointer][5] += DEFAULT_PHRED_SCORE_FOR_DELETION;
						output3_MappingScores[outputArrayPointer][5] += currentMappingScore;
						positionPointer ++;
						outputArrayPointer ++;
						positionToLookFor ++;
						lengthOfCurrentSegment --;
					}
				} else {
					positionPointer += lengthOfCurrentSegment;
				}
			} else if (readSegments[0][i].equals("P")) {
				positionPointer += lengthOfCurrentSegment;
			} else if (readSegments[0][i].equals("M")) {
//				lengthAlreadySearchedOnCurrentRead += lengthOfCurrentSegmentOfCurrentRead;
				if ((positionToLookFor - positionPointer) <= lengthOfCurrentSegment) {
					while (positionToLookFor <= stop && lengthOfCurrentSegment > 0) {
						index2 = ext.indexOfChar(aSingleLineOfBamFile[9].charAt(readPointer), BASES);
						if (index2 >= 0) {
							output1_AlleleCounts[outputArrayPointer][index2] ++;
							output2_PhredScores[outputArrayPointer][index2] += convertToPhredScore(aSingleLineOfBamFile[10].charAt(readPointer));
							output3_MappingScores[outputArrayPointer][index2] += currentMappingScore; 
						} else {
							System.err.println("Error - unrecognized base (neither A, T, G, nor C) at chr" + chr + ":" + positionToLookFor + " Read ID: " + aSingleLineOfBamFile[0]);
						}
						positionPointer ++;
						readPointer ++;
						outputArrayPointer ++;
						positionToLookFor ++;
						lengthOfCurrentSegment --;
					}
				} else {
					positionPointer += lengthOfCurrentSegment;
					readPointer += lengthOfCurrentSegment;
				}
			} else {
				System.err.println("Unrecognized CIGAR string: " + readSegments[0][i]);
			}
		}

//				if (readSegments[0][i].equals("I") || readSegments[0][i].equals("S") || readSegments[0][i].equals("H") || readSegments[0][i].equals("N")) {
////					positionToLookForOnCurrentRead += lengthOfCurrentSegmentOfCurrentRead;
////					lengthAlreadySearchedOnCurrentRead += lengthOfCurrentSegmentOfCurrentRead;
//				} else if (readSegments[0][i].equals("D") || readSegments[0][i].equals("P")) {
////					startPositionOfCurrentRead -= lengthOfCurrentSegmentOfCurrentRead;
//					lengthAlreadySearchedOnCurrentRead += lengthOfCurrentSegmentOfCurrentRead;
//				} else {
//					lengthAlreadySearchedOnCurrentRead += lengthOfCurrentSegmentOfCurrentRead;
//				}
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

	public static void getDenovoMarkerCandidateScores(int[][][] alleleCounts, int[][][] phredScores, int[][][] mappingScores, int[][] output1ReadDepths, byte[] output2DenovoMutationCandidateScores, String[] output3DenovoMutationNotes) {
		boolean[] isThreeAlleles;
		String output3DenovoMutationNote;

		if (output2DenovoMutationCandidateScores.length != alleleCounts[0].length) {
			System.err.println("The length of the array outputDenovoMutationCandidateScores (" + output2DenovoMutationCandidateScores.length + ") is not consistent with that of the array alleleCounts (" + alleleCounts[0].length + ")");
		} else {
			for (int i = 0; i < output2DenovoMutationCandidateScores.length; i++) {
				for (int j = 0; j < alleleCounts.length; j++) {
					output1ReadDepths[j][i] = alleleCounts[j][i][0] + alleleCounts[j][i][1] + alleleCounts[j][i][2] + alleleCounts[j][i][3] + alleleCounts[j][i][4] + alleleCounts[j][i][5];
				}
//				System.out.println(39346502 + i + ":\t" + output1ReadDepths[0][i] + "\t" + output1ReadDepths[1][i] + "\t" + output1ReadDepths[2][i]);
				if (output1ReadDepths[0][i] > MIN_READ_DEPTH && output1ReadDepths[1][i] > MIN_READ_DEPTH && output1ReadDepths[2][i] > MIN_READ_DEPTH) {
					for (int j = 0; j < BASES.length; j++) {
						if (alleleCounts[0][i][j] > MAX_ALLELE_COUNT_TREATED_AS_ZERO && alleleCounts[1][i][j] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO && alleleCounts[2][i][j] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
							isThreeAlleles = isThreeAlleles(alleleCounts, i, MAX_ALLELE_COUNT_TREATED_AS_ZERO);
							output3DenovoMutationNote = "";
							for (int l = 0; l < SAMPLE_SUFFIX.length; l++) {
								if(isThreeAlleles[l]) {
									output2DenovoMutationCandidateScores[i] = 25;
									if (output3DenovoMutationNote.equals("")) {
										output3DenovoMutationNote = SAMPLE_SUFFIX[l];
									} else {
										output3DenovoMutationNote += "," + SAMPLE_SUFFIX[l];
									}
								}
							}
							if (! output3DenovoMutationNote.equals("")) {
								output3DenovoMutationNotes[i] = output3DenovoMutationNote + ":3 alleles";
								output2DenovoMutationCandidateScores[i] = 25;
							} else {
								output2DenovoMutationCandidateScores[i] = 100;
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

	public static void exportInfoForPosition(PrintWriter writer, String id, byte chr, int pos, char ref, byte call, int[][][] alleleCounts, int[][][] phredScores, int[][][] mappingScores, int indexOfCurrentMarker, String[] denovoMarkerNotes, String flag) {
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
				+ "\t(" + (alleleCounts[0][indexOfCurrentMarker][0]==0? "-" : alleleCounts[0][indexOfCurrentMarker][0]) + "," + (alleleCounts[0][indexOfCurrentMarker][1]==0? "-" : alleleCounts[0][indexOfCurrentMarker][1]) + "," + (alleleCounts[0][indexOfCurrentMarker][2]==0? "-" : alleleCounts[0][indexOfCurrentMarker][2]) + "," + (alleleCounts[0][indexOfCurrentMarker][3]==0? "-" : alleleCounts[0][indexOfCurrentMarker][3]) + "," + (alleleCounts[0][indexOfCurrentMarker][4]==0? "-" : alleleCounts[0][indexOfCurrentMarker][4]) + "," + (alleleCounts[0][indexOfCurrentMarker][5]==0? "-" : alleleCounts[0][indexOfCurrentMarker][5]) + ") (" + (alleleCounts[1][indexOfCurrentMarker][0]==0? "-" : alleleCounts[1][indexOfCurrentMarker][0]) + "," + (alleleCounts[1][indexOfCurrentMarker][1]==0? "-" : alleleCounts[1][indexOfCurrentMarker][1]) + "," + (alleleCounts[1][indexOfCurrentMarker][2]==0? "-" : alleleCounts[1][indexOfCurrentMarker][2]) + "," + (alleleCounts[1][indexOfCurrentMarker][3]==0? "-" : alleleCounts[1][indexOfCurrentMarker][3]) + "," + (alleleCounts[1][indexOfCurrentMarker][4]==0? "-" : alleleCounts[1][indexOfCurrentMarker][4]) + "," + (alleleCounts[1][indexOfCurrentMarker][5]==0? "-" : alleleCounts[1][indexOfCurrentMarker][5]) + ") (" + (alleleCounts[2][indexOfCurrentMarker][0]==0? "-" : alleleCounts[2][indexOfCurrentMarker][0]) + "," + (alleleCounts[2][indexOfCurrentMarker][1]==0? "-" : alleleCounts[2][indexOfCurrentMarker][1]) + "," + (alleleCounts[2][indexOfCurrentMarker][2]==0? "-" : alleleCounts[2][indexOfCurrentMarker][2]) + "," + (alleleCounts[2][indexOfCurrentMarker][3]==0? "-" : alleleCounts[2][indexOfCurrentMarker][3]) + "," + (alleleCounts[2][indexOfCurrentMarker][4]==0? "-" : alleleCounts[2][indexOfCurrentMarker][4]) + "," + (alleleCounts[2][indexOfCurrentMarker][5]==0? "-" : alleleCounts[2][indexOfCurrentMarker][5]) + ")" + (denovoMarkerNotes[indexOfCurrentMarker]==null? "" : (" " + denovoMarkerNotes[indexOfCurrentMarker]))
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

	public static String rootOf(String[] filenames) {
		String[] roots;
		String commomRoot;

		roots = new String[filenames.length];
		for (int i = 0; i < roots.length; i++) {
			roots[i] = ext.rootOf(filenames[i]);
		}
		
//		return commonRoot;
		return "F10639";
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String[] commands;
		String fileNameOfDeNovoPointMutationCandidateList;
		String[] bamFilenamesTheTrio;
		String bamDir;
		String scriptDir;
		String outputDir;
		boolean isToAnnotate;
		boolean isToProcessRegion;
		String bamFilenamesForHelpMenu;
		byte chr;
		int start;
		int stop;
		String refFastaFilename;

		fileNameOfDeNovoPointMutationCandidateList = "N:/statgen/OS_Logan/IGV_validation/results_test1.txt";
		bamDir = "/home/spectorl/shared/exome_processing/bam/";
		bamFilenamesTheTrio = new String[] {"rrd_F10639C_L008.bam", "rrd_F10639D_L007.bam", "rrd_F10639M_L007.bam"};
//		bamDir = "D:/bams/";
//		bamFilenamesTheTrio = new String[] {"F10639C_chr17_39346502_39346622.txt", "F10639D_chr17_39346502_39346622.txt", "F10639M_chr17_39346502_39346622.txt"};
		refFastaFilename = "/home/spectorl/xuz2/hg19_canonical.fa";
		scriptDir = "D:/logan/DeNovos/scripts/";
		outputDir = "/home/spectorl/xuz2/outputs/";
		isToAnnotate = false;
		isToProcessRegion = true;
		chr = (byte) 17;
		start = 39346622;
		stop = 39346622;
//		start = 39345622;
//		stop = 39346722;

		bamFilenamesForHelpMenu = "";
		for (int i = 0; i < bamFilenamesTheTrio.length; i++) {
			bamFilenamesForHelpMenu += bamFilenamesTheTrio[i] + ";";
		}
		bamFilenamesForHelpMenu.substring(0, bamFilenamesForHelpMenu.length() - 1);

		commands = new String[] {"-annotation", "candidate=", "bamdir=", "scriptdir=", "outputdir=", "-processregion", "bamset=", "reffasta=", "chr=", "start=", "stop="};
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
				bamFilenamesTheTrio = args[i].split("=")[1].split(";");
				numArgs--;
			} else if (args[i].startsWith(commands[7])) {
				chr = Byte.parseByte(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith(commands[8])) {
				chr = Byte.parseByte(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith(commands[9])) {
				start = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith(commands[10])) {
				stop = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		if (isToAnnotate) {
//			getAlleleCounts("D:/bams/F10639C_chr11_89531764.txt", 11, 89531764, 3);	//Testing only. Feature already included in screenDeNovoPointMutation()
//			getAlleleCounts("D:/bams/F10639C_chr5_140558628.txt", 5, 140558628, 3);
			generateScriptForSamtools(fileNameOfDeNovoPointMutationCandidateList, bamFilenamesTheTrio[0], scriptDir);
			screenDeNovoPointMutation(fileNameOfDeNovoPointMutationCandidateList, bamFilenamesTheTrio[0], 1);
		} else if (isToProcessRegion) {
			processRegion(bamDir, bamFilenamesTheTrio, refFastaFilename, chr, start, stop, outputDir);
		}
	}

}
