package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.EOFException;
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

import org.genvisis.common.Array;
import org.genvisis.common.CountHash;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.Positions;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.filesys.SegmentLists;

import com.google.common.base.Strings;

public class SuperNovo {
  public static final String[] SAMPLE_SUFFIX = new String[] {"C", "D", "M"};
  public static final char[] BASES = new char[] {'A', 'T', 'G', 'C'};
  public static final char[] BASES_WITH_N = new char[] {'A', 'T', 'G', 'C', 'N'};
  public static final String[] BASES_WITH_N_INDEL = new String[] {"A", "T", "G", "C", "N", "Ins",
                                                                  "Del"};
  public static final String[] READS_COUNTS_ARRAY_STRUCT = new String[] {"A", "T", "G", "C", "N",
                                                                         "Ins", "Del",
                                                                         "numAllelesStrictCount",
                                                                         "numAllelesLooseCount"};
  public static final byte INDEX_OF_N = 4;
  public static final byte INDEX_OF_INS = 5;
  public static final byte INDEX_OF_DEL = 6;
  public static final byte INDEX_OF_NUM_ALLELES_STRICT = 7;
  public static final byte INDEX_OF_NUM_ALLELES_LOOSE = 8;
  public static final String[] DENOVO_MUTATION_ARRAY_STRUCT =
                                                            new String[] {"i", "score",
                                                                          "C_readDepth",
                                                                          "D_readDepth",
                                                                          "M_readDepth",
                                                                          "C_allele1", "C_allele2",
                                                                          "D_allele1", "D_allele2",
                                                                          "M_allele1", "M_allele2"};
  public static final String[] DENOVO_MUTATION_NOTES_ARRAY_STRUCT = new String[] {"Ins", "Del",
                                                                                  "Var", "3Allelic",
                                                                                  "DeNovo"};
  public static final int MIN_READ_DEPTH = 8;
  public static final int MIN_MAPPING_SCORE = 45;
  public static final int MIN_PHRED = 26;
  public static final int DEFAULT_PHRED_SCORE_FOR_DELETION = 30;
  public static final double THRESHOLD_PHRED_SCORE_FOR_INS_DEL = .10;
  public static final int MAX_ALLELE_COUNT_TREATED_AS_ZERO = 2;
  public static final int MIN_ALLELE_COUNT_FOR_DENOVO_MUTATION = 5;
  public static final double MIN_ALLELE_FREQ_FOR_DENOVO_MUTATION = .20;
  public static final int MAX_NUM_VAR_INS_DEL_DENOVO_NEARBY = 5;
  public static final int REGION_LEN_AT_A_TIME_DEFAULT = 100000;
  public static final int WINDOW_SIZE_FOR_NEARBY_INDEL = 60;
  public static final int WINDOW_SIZE_FOR_NEARBY_VARIANCE = 30;
  public static final int WINDOW_SIZE_FOR_NEARBY_DENOVO = 30;
  public static final int WINDOW_SIZE_FOR_HAPLOTYPE_CHECK = 102;
  public static final double DISCOUNT_FOR_NEARBY_INDEL = .80;
  public static final double DISCOUNT_FOR_ON_INDEL_SITE = .50;
  public static final double DISCOUNT_FOR_NEARBY_VARIANCE = .95;
  public static final double DISCOUNT_FOR_ON_VARIANCE_SITE = .70;
  public static final double DISCOUNT_FOR_AT_VARIANCE = .60;
  public static final double DISCOUNT_FOR_N = .50;
  public static final double DISCOUNT_FOR_3ALLELES_LOOSE = .25;
  public static final double DISCOUNT_FOR_3ALLELES_STRICT = .50;
  public static final double DISCOUNT_FOR_NEARBY_DENOVOMUTATION = .90;
  // public static final String OUTPUT_FILE_HEADER =
  // "id\tchr\tpos\tlookup\tsarver\tref\talt\tmendelianLikelihood\tmendelianPP\tmendelianGT\tsnpCode\tcode\tdeNovoLikelihood\tdeNovoPP\tactualDenovo\tconf\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tchildQuality\tdadQuality\tmomQuality";
  // public static final String OUTPUT_FILE_HEADER =
  // "id\tchr\tpos\tlookup\tref\talt\tcall\tnote\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore\t1\t2\t4\t5\t7\t8";
  public static final String OUTPUT_FILE_HEADER =
                                                "id\tchr\tpos\tlookup\tref\talt\tcall\treadsCounts\tnumNearbyInDelVars\tposNearbyInDelVars\tdeNovoGT\tflag\tchildDepth\tdadDepth\tmomDepth\tPhredScores\tchildMappingScore\tdadMappingScore\tmomMappingScore";
  // public static final String OUTPUT_FILE_HEADER2 =
  // "id\tchr\tpos\tlookup\tref\talt\tcall\tcA\tcT\tcG\tcC\tcN\tcI\tcD\tdA\tdT\tdG\tdC\tdN\tdI\tdD\tmA\tmT\tmG\tmC\tmN\tmI\tmD\tcIns\tdIns\tmIns\tcDel\tdDel\tmDel\tcVar\tdVar\tmVar\tposNearbyInDelVars\tdeNovoGT\tflag\tcDepth\tdDepth\tmDepth\tcAPhred\tcTPhred\tcGPhred\tcCPhred\tdAPhred\tdTPhred\tdGPhred\tdCPhred\tmAPhred\tmTPhred\tmGPhred\tmCPhred\tcMap\tdMap\tmMap";
  public static final String OUTPUT_FILE_HEADER2 =
                                                 "id\tchr\tpos\tlookup\tref\talt\tcall\t#Reads_A_Child\t#Reads_T_Child\t#Reads_G_Child\t#Reads_C_Child\t#Reads_N_Child\t#Reads_I_Child\t#Reads_D_Child\t#Reads_A_Dad\t#Reads_T_Dad\t#Reads_G_Dad\t#Reads_C_Dad\t#Reads_N_Dad\t#Reads_I_Dad\t#Reads_D_Dad\t#Reads_A_Mom\t#Reads_T_Mom\t#Reads_G_Mom\t#Reads_C_Mom\t#Reads_N_Mom\t#Reads_I_Mom\t#Reads_D_Mom\t#Alleles_Child\t#Alleles_Dad\t#Alleles_Mom\t#NearbySites_I_Child\t#NearbySites_I_Dad\t#NearbySites_I_Mom\t#NearbySites_D_Child\t#NearbySites_D_Dad\t#NearbySites_D_Mom\t#NearbySites_Var_Child\t#NearbySites_Var_Dad\t#NearbySites_Var_Mom\t#NearbySites_3_Child\t#NearbySites_3_Dad\t#NearbySites_3_Mom\t#NearbySites_DeNovoVar\tposNearbyInDelVars\tdeNovoGT\tflag\tDepth_Child\tDepth_Dad\tDepth_Mom\t%Phred_A_Child\t%Phred_T_Child\t%Phred_G_Child\t%Phred_C_Child\t%Phred_A_Dad\t%Phred_T_Dad\t%Phred_G_Dad\t%Phred_C_Dad\t%Phred_A_Mom\t%Phred_T_Mom\t%Phred_G_Mom\t%Phred_C_Mom\tPhred_A_Child\tPhred_T_Child\tPhred_G_Child\tPhred_C_Child\tPhred_N_Child\tPhred_I_Child\tPhred_D_Child\tPhred_A_Dad\tPhred_T_Dad\tPhred_G_Dad\tPhred_C_Dad\tPhred_N_Dad\tPhred_I_Dad\tPhred_D_Dad\tPhred_A_Mom\tPhred_T_Mom\tPhred_G_Mom\tPhred_C_Mom\tPhred_N_Mom\tPhred_I_Mom\tPhred_D_Mom\tMapping_Child\tMapping_Dad\tMapping_Mom";
  public static final String PARSE_RESULT_HEADER =
                                                 "chr\tpos\tlookup\tgeneCount_Pos\tgeneCount_Trio\tgeneCount_diff\tid\tnumVars\tnumDels\tnumIns\tnumDenovo\tgeneList\thasBam\tref\talt\tcall\tnote\tfwdGenotype\treadDepth_C\treadDepth_D\treadDepth_M\tPhred_A_Child\tPhred_T_Child\tPhred_G_Child\tPhred_C_Child\tPhred_A_Dad\tPhred_T_Dad\tPhred_G_Dad\tPhred_C_Dad\tPhred_A_Mom\tPhred_T_Mom\tPhred_G_Mom\tPhred_C_Mom\tMappingQuality_C\tMappingQuality_D\tMappingQuality_M\taltAllele%_C\taltAllele%_D\taltAllele%_M";
  public static final String PARSE_RESULT_HEADER_FORMAT2 = "chr\tpos\tlookup\tnumTrios";
  public static final String READ_COUNTS_HEADER_FORMAT1 =
                                                        "chr\tpos\t#Reads_byAllele_C\t#Reads_byAllele_D\t#Reads_byAllele_M\t#Reads_total_C\t#Reads_total_D\t#Reads_total_M\tMap_C\tMap_D\tMap_M\tPhred_C\tPhred_D\tPhred_M\tHaploCount_C\tHaploCount_D\tHaploCount_M";
  public static final String READ_COUNTS_HEADER_FORMAT2 =
                                                        "chr\tpos\t#Reads_total_C\t#Reads_total_D\t#Reads_total_M\tMap_C\tMap_D\tMap_M\tPhred_C\tPhred_D\tPhred_M";
  public static final String[] COLUMNS_NEEDED_FROM_PHASE1_OUTPUT =
                                                                 {"ref", "alt", "call",
                                                                  "posNearbyInDelVars", "deNovoGT",
                                                                  "Depth_Child", "Depth_Dad",
                                                                  "Depth_Mom", "Phred_A_Child",
                                                                  "Phred_T_Child", "Phred_G_Child",
                                                                  "Phred_C_Child", "Phred_A_Dad",
                                                                  "Phred_T_Dad", "Phred_G_Dad",
                                                                  "Phred_C_Dad", "Phred_A_Mom",
                                                                  "Phred_T_Mom", "Phred_G_Mom",
                                                                  "Phred_C_Mom", "Mapping_Child",
                                                                  "Mapping_Dad", "Mapping_Mom"};
  public static final String DEFAULT_OUTPUT_FILENAME = "SuperNovo_summary.xln";
  public static final String DEFAULT_READ_COUNTS_FILENAME_SUFFIX = ".txt.gz";
  public static final String INSERTION_DELIMITER_LEFT = " < ";
  public static final String INSERTIONS_DELIMITER_RIGHT = " > ";

  public static void generateScriptForSamtools(String fileNameOfDeNovoPointMutationCandidateList,
                                               String bamFileDir, String scriptDir) {
    BufferedReader candidateList;
    String[] line;
    String commands;
    Vector<String[]> iterationsVec;
    String[][] iterations;

    commands = "samtools view [%0]C.bam  [%1]:[%2]-[%3] > " + "[%0]C_[%1]_[%3].txt\n"
               + "samtools view [%0]D.bam  [%1]:[%2]-[%3] > " + "[%0]D_[%1]_[%3].txt\n"
               + "samtools view [%0]M.bam  [%1]:[%2]-[%3] > " + "[%0]M_[%1]_[%3].txt";
    try {
      candidateList =
                    new BufferedReader(new FileReader(fileNameOfDeNovoPointMutationCandidateList));
      candidateList.readLine();
      iterationsVec = new Vector<String[]>();
      while (candidateList.ready()) {
        line = candidateList.readLine().split("\t");
        iterationsVec.add(new String[] {line[0].substring(0, line[0].length() - 1), line[1],
                                        line[2], line[3]});
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

  public static void screenDeNovoPointMutation(String fileNameOfDeNovoPointMutationCandidateList,
                                               String bamFileDir, int thresholdFor3Alleles) {
    BufferedReader candidateList;
    PrintWriter writer;
    String[] line;
    int[][] readsCounts;
    String notes;
    int score;

    try {
      candidateList =
                    new BufferedReader(new FileReader(fileNameOfDeNovoPointMutationCandidateList));
      candidateList.readLine();
      writer = new PrintWriter(ext.parseDirectoryOfFile(fileNameOfDeNovoPointMutationCandidateList)
                               + ext.rootOf(fileNameOfDeNovoPointMutationCandidateList)
                               + "_screened.txt");
      writer.println("ID\tchrs\tposition\tlookup\tSarver\tREF\tALT\tMendelianLikelihood\tMendelianPP\tMendelianGT\tSnpCode\tCode\tDeNovoLikelihood\tDeNovoPP\tActualDenovo\tConf\tPankratzScore\tNotes");
      while (candidateList.ready()) {
        line = candidateList.readLine().split("\t");
        // look for the chr location in bcfFile
        readsCounts = new int[3][];
        for (int i = 0; i < readsCounts.length; i++) {
          readsCounts[i] = getAlleleCounts(bamFileDir + line[0].substring(0, line[0].length() - 1)
                                           + SAMPLE_SUFFIX[i] + "_" + line[1] + "_" + line[3]
                                           + ".txt", Integer.parseInt(line[1].split("chr")[1]),
                                           Integer.parseInt(line[3]), thresholdFor3Alleles);
        }
        score = getScoreForDeNovoMutation(readsCounts, line[9].charAt(0), line[8].charAt(0),
                                          thresholdFor3Alleles);
        notes = getNotesForDeNovoMutation(readsCounts, line[9].charAt(0), thresholdFor3Alleles);
        writer.println(line[0] + "\t" + line[1] + "\t" + line[3] + "\t" + line[6] + "\t" + line[7]
                       + "\t" + line[8] + "\t" + line[9] + "\t" + line[10] + "\t" + line[11] + "\t"
                       + line[12] + "\t" + line[13] + "\t" + line[14] + "\t" + line[15] + "\t"
                       + line[16] + "\t" + line[17] + "\t" + line[18] + "\t" + score + "\t"
                       + notes);
      }
      candidateList.close();
      writer.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static int[] getAlleleCounts(String readsFileFullPath, int chr, int position,
                                      int thresholdFor3Alleles) {
    int[] result;
    BufferedReader bamFile;
    String[] line;
    // char[] bases;
    int index, startPositionOfTheRead, lengthOfCurrentSegmentOfCurrentRead, positionOnThisRead;
    String[][] operatorsOperatorIndicesAndSplit;

    result = new int[] {0, 0, 0, 0, 0, 0, 0, 0}; // A, T, G, C, N, Ins, Del, Total, #Alleles
    try {
      bamFile = new BufferedReader(new FileReader(readsFileFullPath));
      while (bamFile.ready()) {
        positionOnThisRead = position;
        line = bamFile.readLine().split("\t");
        if (line[2].equalsIgnoreCase("chr" + chr) && Integer.parseInt(line[3]) <= position
            && (Integer.parseInt(line[3]) + Math.abs(Integer.parseInt(line[8]))) >= position
            && Integer.parseInt(line[1]) < 512) {
          startPositionOfTheRead = Integer.parseInt(line[3]);
          index = 0;
          operatorsOperatorIndicesAndSplit = ext.getOperatorsOperatorIndicesAndSplit(line[5],
                                                                                     "DIMNPSH");
          for (int i = 0; i < operatorsOperatorIndicesAndSplit[0].length; i++) {
            lengthOfCurrentSegmentOfCurrentRead =
                                                Integer.parseInt(operatorsOperatorIndicesAndSplit[2][i]);
            if ((startPositionOfTheRead + index
                 + lengthOfCurrentSegmentOfCurrentRead) > positionOnThisRead) {
              if (operatorsOperatorIndicesAndSplit[0][i].equals("I")
                  || operatorsOperatorIndicesAndSplit[0][i].equals("S")
                  || operatorsOperatorIndicesAndSplit[0][i].equals("H")
                  || operatorsOperatorIndicesAndSplit[0][i].equals("N")) {
                positionOnThisRead += lengthOfCurrentSegmentOfCurrentRead;
                index += lengthOfCurrentSegmentOfCurrentRead;
              } else if (operatorsOperatorIndicesAndSplit[0][i].equals("D")
                         || operatorsOperatorIndicesAndSplit[0][i].equals("P")) {
                result[6]++;
                break;
              } else {
                index += (positionOnThisRead - startPositionOfTheRead - index);
                // index = (positionOnThisRead - startPositionOfTheRead);
                index = ext.indexOfChar(line[9].charAt(index), BASES_WITH_N);
                if (index >= 0) {
                  result[index]++;
                } else {
                  System.out.println("Error - unrecognized base " + line[9].charAt(index)
                                     + " in the following file and read:\n" + readsFileFullPath
                                     + "\nRead ID: " + line[0]);
                }
                break;
              }
            } else {
              if (operatorsOperatorIndicesAndSplit[0][i].equals("I")
                  || operatorsOperatorIndicesAndSplit[0][i].equals("S")
                  || operatorsOperatorIndicesAndSplit[0][i].equals("H")
                  || operatorsOperatorIndicesAndSplit[0][i].equals("N")) {
                positionOnThisRead += lengthOfCurrentSegmentOfCurrentRead;
                index += lengthOfCurrentSegmentOfCurrentRead;
              } else if (operatorsOperatorIndicesAndSplit[0][i].equals("D")
                         || operatorsOperatorIndicesAndSplit[0][i].equals("P")) {
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

    // FIXME: was referencing result[8] which is guaranteed to be AIOOB
    // presumably this step is totaling ATGC counts and counting them as alleles if the count
    // exceeds the threshold?
    for (int i = 0; i < 4; i++) {
      result[6] += result[i];
      if (result[i] > thresholdFor3Alleles) {
        result[7]++;
      }
    }

    return result;
  }

  /**
   * Annotate the markers with allele frequencies, allele counts, and deletion and insertion counts.
   *
   * @param readsCounts
   * @param alternativeAllele
   * @param thresholdFor3Alleles
   * @return
   */
  public static String getNotesForDeNovoMutation(int[][] readsCounts, char alternativeAllele,
                                                 int thresholdFor3Alleles) {
    byte index;
    index = (byte) ext.indexOfChar(alternativeAllele, BASES);

    return ((readsCounts[1][index] > thresholdFor3Alleles
             || readsCounts[2][index] > thresholdFor3Alleles) ? ext.formDeci(100
                                                                             * readsCounts[0][index]
                                                                             / (double) readsCounts[0][5],
                                                                             1)
                                                                + "%," + readsCounts[1][index] + ","
                                                                + readsCounts[2][index] + "; "
                                                              : "")
           + "C" + (readsCounts[0][6] > 2 ? "(" + readsCounts[0][6] + " alleles)" : "") + ":"
           + readsCounts[0][0] + "," + readsCounts[0][1] + "," + readsCounts[0][2] + ","
           + readsCounts[0][3] + ",(" + readsCounts[0][4] + "," + readsCounts[0][5] + ")" + "; D"
           + (readsCounts[1][6] > 2 ? "(" + readsCounts[1][6] + "alleles)" : "") + ":"
           + readsCounts[1][0] + "," + readsCounts[1][1] + "," + readsCounts[1][2] + ","
           + readsCounts[1][3] + ",(" + readsCounts[1][4] + "," + readsCounts[1][5] + ")" + "; M"
           + (readsCounts[2][6] > 2 ? "(" + readsCounts[2][6] + "alleles)" : "") + ":"
           + readsCounts[2][0] + "," + readsCounts[2][1] + "," + readsCounts[2][2] + ","
           + readsCounts[2][3] + ",(" + readsCounts[2][4] + "," + readsCounts[2][5] + ")";
  }

  /**
   * A decimal value ranging 0 through 1 to indicate how likely the marker is a De Novo point
   * mutation, with 1 being very likely and 0 very unlikely
   *
   * @param readsCounts
   * @param alternativeAllele
   * @param referencedAllele
   * @param thresholdFor3Alleles
   * @return
   */
  // This is functioning but naive. Has been replaced by getDenovoMarkerCandidateScores(...)
  public static int getScoreForDeNovoMutation(int[][] readsCounts, char alternativeAllele,
                                              char referencedAllele, int thresholdFor3Alleles) {
    int DenovoMutationScore;
    int index;

    index = ext.indexOfChar(alternativeAllele, BASES);
    if (readsCounts[1][index] > thresholdFor3Alleles || readsCounts[2][index] > thresholdFor3Alleles
        || readsCounts[0][8] > 2) {
      DenovoMutationScore = 0;
    } else {
      DenovoMutationScore = 1;
    }
    return DenovoMutationScore;
  }

  public static void generateScriptsAllTriosInADirectoryToScanSamFilesForDenovoMutation(String bamDir,
                                                                                        String refFastaFilename,
                                                                                        String bimFilename,
                                                                                        String outputDir,
                                                                                        String scriptDir,
                                                                                        String fullPathToTrioNameList,
                                                                                        String pathAndRootToOutputReadCounts,
                                                                                        int regionLegnthATime,
                                                                                        int numThreads,
                                                                                        Logger log) {
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
    command = "cd " + outputDir + "\n" + Files.getRunString() + " one.SuperNovo -denovo bim="
              + bimFilename + " outdir=" + outputDir + " bamdir=" + bamDir + " reffasta="
              + refFastaFilename + " regionlengthatime=" + regionLegnthATime + " numthreads="
              + numThreads + " -denovo";
    for (String[] bamFilenamesByTrio : bamFilenamesByTrios) {
      // processGenomeOfOneTrio(bamDir, bamFilenames, refFastaFilename, bedFilename, outputDir,
      // numThreads, log);
      trioId = bamFilenamesByTrio[0];
      if (trioId == null || trioId.equals("")) {
        trioId = getRootOf(bamFilenamesByTrio, true);
      }
      Files.qsub(scriptDir + trioId + ".qsub",
                 command + " trioid=" + trioId + " bamset=" + bamFilenamesByTrio[1] + ","
                                               + bamFilenamesByTrio[2] + "," + bamFilenamesByTrio[3]
                                               + (pathAndRootToOutputReadCounts == null ? ""
                                                                                        : (" outputreadcounts="
                                                                                           + pathAndRootToOutputReadCounts
                                                                                           + trioId
                                                                                           + DEFAULT_READ_COUNTS_FILENAME_SUFFIX)),
                 2000, 11, 1);
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

  public static void scanMultipleRegionsInSamFilesOfATrioForDenovoMutations(String trioId,
                                                                            String[] bamFilenamesOfTheTrio,
                                                                            String refFastaFilename,
                                                                            String bimFilename,
                                                                            String outputDir,
                                                                            String fullpathToSaveAlleleCounts,
                                                                            int regionLengthATime,
                                                                            int numThreads,
                                                                            Logger log) {
    BufferedReader reader;
    String[] line;
    String chr;
    // String prevChr;
    int start;
    int stop;
    Vector<String> chrs = null;
    Vector<Integer> starts = null;
    Vector<Integer> stops = null;
    int loop;
    PrintWriter writer;
    PrintWriter alleleCountWriter;
    // String trioId;
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
    outFileName = outputDir + trioId + "_superNovo_" + ext.rootOf(bimFilename) + ".txt";
    // prevChr = "start";
    try {
      reader = Files.getAppropriateReader(bimFilename);
      chrs = new Vector<String>();
      starts = new Vector<Integer>();
      stops = new Vector<Integer>();
      while (reader.ready()) {
        line = reader.readLine().split("\t");
        if (!line[0].equals("") && line[0] != null && line[0].toLowerCase().startsWith("chr")) {
          chr = line[0].substring(3).trim();
          // if (!chr.equals(prevChr) && numThreads == 1) {
          // log.report(ext.getTime() + "\t" + "Starting chr" + chr + "; identified " +
          // Files.countLines(outFileName, true) + " possible de novo events so far");
          // prevChr = chr;
          // }
          if (line.length < 2 || line[1] == null || line[1].equals("")) {
            start = 1;
            stop = Positions.CHROMOSOME_LENGTHS_MAX[ext.indexOfStr(chr, Positions.CHR_CODES)];
            log.report("No start position is specified. Using default region: chr" + chr + ":"
                       + start + "-" + stop);
          } else {
            start = Integer.parseInt(line[1]);
            stop = Integer.parseInt(line[2]);
          }
          while (start <= stop) {
            loop = Math.min(start + regionLengthATime - 1, stop);
            // if (numThreads > 1) {
            chrs.add(chr);
            starts.add(start);
            stops.add(loop);
            // } else {
            // processRegion(bamDir, bamFilenamesOfTheTrio, trioId, refFastaFilename, chr, start,
            // loop, writer, null, log);
            // }
            start += regionLengthATime;
          }
        } else {
          log.reportError("Warning: unrecognized chr number '" + line[0] + "' in the file '"
                          + bimFilename + "'. Skipped the region.");
        }
      }
      reader.close();

      if (fullpathToSaveAlleleCounts != null) {
        alleleCountWriter = Files.getAppropriateWriter(fullpathToSaveAlleleCounts);
        alleleCountWriter.println(READ_COUNTS_HEADER_FORMAT2);
      } else {
        alleleCountWriter = null;
      }
      writer = new PrintWriter(outFileName);
      writer.println(OUTPUT_FILE_HEADER2);
      if (numThreads > 1) {
        if (chrs.size() < numThreads) {
          numThreads = chrs.size();
        }
        executor = Executors.newFixedThreadPool(numThreads);
        numElementsATime = (int) Math.ceil((double) chrs.size() / numThreads);
        startElement = 0;
        for (int i = 0; i < numThreads; i++) {
          if ((i + 1) == numThreads) {
            numElementsATime = chrs.size() % numElementsATime;
          }
          executor.execute(new WorkerThread(bamFilenamesOfTheTrio, trioId, refFastaFilename, chrs,
                                            starts, stops, startElement, numElementsATime, writer,
                                            alleleCountWriter, i, log));
          startElement += numElementsATime;
        }
        executor.shutdown();
        executor.awaitTermination(1, TimeUnit.HOURS);
      } else {
        loop = chrs.size();
        for (int i = 0; i < loop; i++) {
          scanARegionInSamFilesOfATrioForDenovoMutations(bamFilenamesOfTheTrio, trioId,
                                                         refFastaFilename, chrs.elementAt(i),
                                                         starts.elementAt(i), stops.elementAt(i),
                                                         writer, alleleCountWriter, (byte) 2, log);
        }
      }
      writer.close();
      if (fullpathToSaveAlleleCounts != null) {
        alleleCountWriter.close();
      }
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + bimFilename + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + bimFilename + "\"");
      return;
    } catch (InterruptedException e) {
    }

    log.report("DeNovo mutation result is ready at: " + outFileName + "\nTotal time used "
               + timeFormat.format(new Date().getTime() - timer));
  }

  public static class WorkerThread implements Runnable {
    private final String[] bamFilenames;
    private final String trioId;
    private final String refFastaFilename;
    private final String[] chrs;
    private final int[] startPositions;
    private final int[] stopPositions;
    private final PrintWriter writer;
    private final PrintWriter fullpathToOutputAlleleCounts;
    private final int threadId;
    private final Logger log;

    public WorkerThread(String[] bamFilenames, String trioId, String refFastaFilename,
                        Vector<String> chrs, Vector<Integer> starts, Vector<Integer> stops,
                        int startIndexOfTheVectors, int numElementsOfTheVectors, PrintWriter writer,
                        PrintWriter fullpathToOutputAlleleCounts, int threadId, Logger log) {
      super();
      this.bamFilenames = bamFilenames;
      this.trioId = trioId;
      this.refFastaFilename = refFastaFilename;
      this.chrs = new String[numElementsOfTheVectors];
      startPositions = new int[numElementsOfTheVectors];
      stopPositions = new int[numElementsOfTheVectors];
      for (int i = 0; i < numElementsOfTheVectors; i++) {
        this.chrs[i] = chrs.elementAt(i + startIndexOfTheVectors);
        startPositions[i] = starts.elementAt(i + startIndexOfTheVectors);
        stopPositions[i] = stops.elementAt(i + startIndexOfTheVectors);
      }
      this.writer = writer;
      this.fullpathToOutputAlleleCounts = fullpathToOutputAlleleCounts;
      this.threadId = threadId;
      this.log = log;
    }

    @Override
    public void run() {
      for (int i = 0; i < chrs.length; i++) {
        scanARegionInSamFilesOfATrioForDenovoMutations(bamFilenames, trioId, refFastaFilename,
                                                       chrs[i], startPositions[i], stopPositions[i],
                                                       writer, fullpathToOutputAlleleCounts,
                                                       (byte) 2, log);
      }
    }

    public int getThredId() {
      return threadId;
    }
  }

  // public static byte[][] processRegion(String dir, String bamFilename, byte chr, int start, int
  // stop) {
  // PrintStream stream;
  // similar to piping or redirecting in linux
  // BufferedReader reader = new BufferedReader();
  // CmdLine.run("samtools view "+bamFilename+" chr"+chr+":"+start+"-"+stop, dir, stream);

  public static void scanARegionInSamFilesOfATrioForDeNovoMutations(String[] bamFilenames,
                                                                    String refFastaFilename,
                                                                    String chr, int start, int stop,
                                                                    String outputDir,
                                                                    String fullpathToOutputAlleleCounts,
                                                                    Logger log) {
    PrintWriter writer;
    PrintWriter alleleCountWriter;
    String trioId;
    String outFileName;
    long timer;
    SimpleDateFormat timeFormat;

    trioId = getRootOf(bamFilenames, false);
    outFileName = outputDir + trioId + "_denovoMutations_chr" + chr + "_" + start + "_" + stop
                  + ".txt";
    timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
    timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
    timer = new Date().getTime();
    try {
      if (fullpathToOutputAlleleCounts != null) {
        alleleCountWriter = Files.getAppropriateWriter(fullpathToOutputAlleleCounts);
        alleleCountWriter.println(READ_COUNTS_HEADER_FORMAT2);
      } else {
        alleleCountWriter = null;
      }
      writer = new PrintWriter(outFileName);
      writer.println(OUTPUT_FILE_HEADER2);
      scanARegionInSamFilesOfATrioForDenovoMutations(bamFilenames, trioId, refFastaFilename, chr,
                                                     start, stop, writer, alleleCountWriter,
                                                     (byte) 2, log);
      writer.close();
      if (fullpathToOutputAlleleCounts != null) {
        alleleCountWriter.close();
      }
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }

    timer = new Date().getTime() - timer;
    System.out.println("processRegion result is ready at: " + outFileName
                       + (fullpathToOutputAlleleCounts != null ? "\nAllele counts result is ready at:"
                                                                 + fullpathToOutputAlleleCounts
                                                               : "")
                       + "\nTotal time used " + timeFormat.format(timer));
  }

  public static void scanARegionInSamFilesOfATrioForDenovoMutations(String[] samFilenames,
                                                                    String trioId,
                                                                    String refFastaFilename,
                                                                    String chr, int startPosition,
                                                                    int stopPosition,
                                                                    PrintWriter writer,
                                                                    PrintWriter fullpathToSaveReadCounts,
                                                                    byte readCountsFileFormat,
                                                                    Logger log) {
    Vector<String[]> samContentVec;
    int numLines;
    int numMarkersPlusWindow;
    int window;
    int[][][] readsCounts = null;
    int[][][] phredScores = null;
    int[][][] mappingScores = null;
    String[] refAlleles;
    Vector<int[]> denovoMutations;
    Vector<Vector<Integer>[][]> denovoMutationNotes;
    Hashtable<String, String> altAllelesForInsertions_Child;
    int startPosAdjForWindow;
    int stopPosAdjForWindow;
    long timer;
    String status;

    window = Math.max(WINDOW_SIZE_FOR_NEARBY_INDEL, WINDOW_SIZE_FOR_NEARBY_VARIANCE);
    startPosAdjForWindow = Math.max(0, startPosition - window);
    stopPosAdjForWindow = Math.min(
                                   Positions.CHROMOSOME_LENGTHS_MAX[ext.indexOfStr(chr,
                                                                                   Positions.CHR_CODES,
                                                                                   false, true)],
                                   stopPosition + window);
    numMarkersPlusWindow = stopPosAdjForWindow - startPosAdjForWindow + 1;
    readsCounts =
                new int[SAMPLE_SUFFIX.length][numMarkersPlusWindow][READS_COUNTS_ARRAY_STRUCT.length];
    phredScores = new int[SAMPLE_SUFFIX.length][numMarkersPlusWindow][BASES_WITH_N_INDEL.length];
    mappingScores = new int[SAMPLE_SUFFIX.length][numMarkersPlusWindow][BASES_WITH_N_INDEL.length];
    altAllelesForInsertions_Child = new Hashtable<String, String>();
    status = trioId + "\tchr" + chr + ":" + startPosAdjForWindow + "-" + stopPosAdjForWindow + "\t("
             + startPosition + "-" + stopPosition + ")";

    timer = new Date().getTime();
    for (int i = 0; i < samFilenames.length; i++) {
      samContentVec = loadBamByPipeline(samFilenames[i], chr, startPosAdjForWindow,
                                        stopPosAdjForWindow);
      numLines = samContentVec.size();
      for (int j = 0; j < numLines; j++) {
        cigarDecoding(samContentVec.elementAt(j), chr, startPosAdjForWindow, stopPosAdjForWindow, 0,
                      readsCounts[i], phredScores[i], mappingScores[i],
                      (i == 0 ? altAllelesForInsertions_Child : null));
      } ;
      status += ("\t" + numLines);
    }
    if (fullpathToSaveReadCounts != null) {
      saveAlleleCountsToFile(fullpathToSaveReadCounts, chr, startPosAdjForWindow, readsCounts,
                             phredScores, mappingScores, readCountsFileFormat, true);
    }
    status += ("\t" + (new Date().getTime() - timer));
    timer = new Date().getTime();

    denovoMutations = new Vector<int[]>(10);
    getFilteredResults(readsCounts, mappingScores, startPosition - startPosAdjForWindow,
                       stopPosition - startPosAdjForWindow, denovoMutations);
    denovoMutationNotes = new Vector<Vector<Integer>[][]>(denovoMutations.size());
    for (int j = 0; j < denovoMutations.size(); j++) {
      denovoMutationNotes.add(new Vector[DENOVO_MUTATION_NOTES_ARRAY_STRUCT.length][SAMPLE_SUFFIX.length]);
    }
    getNearbyInDelVars(readsCounts, denovoMutations, denovoMutationNotes);
    getNearbyDeNovos(readsCounts, denovoMutations, denovoMutationNotes);
    refAlleles = getRefFromFasta(refFastaFilename, denovoMutations, chr, startPosAdjForWindow);
    exportResult(writer, trioId, chr, startPosAdjForWindow, startPosition, refAlleles,
                 altAllelesForInsertions_Child, readsCounts, phredScores, mappingScores,
                 denovoMutations, denovoMutationNotes,
                 "Phred score proportions < " + THRESHOLD_PHRED_SCORE_FOR_INS_DEL, (byte) 2, log);
    writer.flush();

    status += ("\t" + (new Date().getTime() - timer));
    log.report(status);
  }

  public static void cigarDecoding(String[] aSingleLineOfBamFile, String chr, int startPos,
                                   int stopPos, int thresholdFor3Alleles,
                                   int[][] output1_ReadsCounts, int[][] output2_PhredScores,
                                   int[][] output3_MappingScores,
                                   Hashtable<String, String> altAllelesForInDels) {
    int readPointer;
    int outputArrayPointer;
    int currentPosition;
    int lengthOfCurrentSegment;
    int outputArrayLength;
    int loop;
    int indexInBases;
    String[][] readSegments;
    int currentMappingScore;
    String[] tmp1;
    String tmp2;
    String tmp3;
    boolean found;

    if (aSingleLineOfBamFile[2].equalsIgnoreCase("chr" + chr)
        && Integer.parseInt(aSingleLineOfBamFile[1]) < 512
        && !aSingleLineOfBamFile[9].equalsIgnoreCase("*")) {
      currentPosition = Integer.parseInt(aSingleLineOfBamFile[3]);
      currentMappingScore = Integer.parseInt(aSingleLineOfBamFile[4]);
      readSegments = ext.getOperatorsOperatorIndicesAndSplit(aSingleLineOfBamFile[5], "MIDNSHP=X");
      readPointer = 0;
      outputArrayPointer = currentPosition - startPos;
      outputArrayLength = stopPos - startPos + 1;

      for (int i = 0; (outputArrayPointer < outputArrayLength)
                      && (i < readSegments[0].length); i++) {
        lengthOfCurrentSegment = Integer.parseInt(readSegments[2][i]);
        if (readSegments[0][i].equals("M") || readSegments[0][i].equals("=")
            || readSegments[0][i].equals("X")) { // present in read sequence, and also in reference
                                                 // sequence.
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

            loop = outputArrayPointer
                   + Math.min(lengthOfCurrentSegment, outputArrayLength - outputArrayPointer);
            while (outputArrayPointer < loop) {
              indexInBases = ext.indexOfChar(aSingleLineOfBamFile[9].charAt(readPointer),
                                             BASES_WITH_N);
              if (indexInBases >= 0) {
                output1_ReadsCounts[outputArrayPointer][indexInBases]++;
                output2_PhredScores[outputArrayPointer][indexInBases] +=
                                                                      convertToPhredScore(aSingleLineOfBamFile[10].charAt(readPointer));
                output3_MappingScores[outputArrayPointer][indexInBases] += currentMappingScore;
              } else {
                System.err.println("Error - unrecognized base ("
                                   + aSingleLineOfBamFile[9].charAt(readPointer) + ") at chr" + chr
                                   + ":" + startPos + outputArrayPointer + " Read ID: "
                                   + aSingleLineOfBamFile[0]);
              }
              readPointer++;
              outputArrayPointer++;
            }

            currentPosition += lengthOfCurrentSegment;
          }

        } else if (readSegments[0][i].equals("I")) { // present in read sequence, but NOT in
                                                     // reference sequence.
          if (outputArrayPointer >= 0) {
            output1_ReadsCounts[outputArrayPointer][INDEX_OF_INS]++;
            output2_PhredScores[outputArrayPointer][INDEX_OF_INS] +=
                                                                  convertToPhredScore(aSingleLineOfBamFile[10].charAt(readPointer));
            output3_MappingScores[outputArrayPointer][INDEX_OF_INS] += currentMappingScore;
            if (altAllelesForInDels != null) {
              tmp3 = "+" + aSingleLineOfBamFile[9].substring(readPointer,
                                                             readPointer + lengthOfCurrentSegment);
              if (altAllelesForInDels.containsKey(currentPosition + "")) {
                tmp2 = altAllelesForInDels.get(currentPosition + "");
                if (!tmp2.equalsIgnoreCase(tmp3)) {
                  found = false;
                  tmp1 = tmp2.split(",");
                  for (String element : tmp1) {
                    if (element.equalsIgnoreCase(tmp3)) {
                      found = true;
                      break;
                    }
                  }
                  if (!found) {
                    altAllelesForInDels.put(currentPosition + "", tmp2 + "," + tmp3);
                  }
                }
              } else {
                altAllelesForInDels.put(currentPosition + "", tmp3);
              }
            }
          }
          readPointer += lengthOfCurrentSegment;

        } else if (readSegments[0][i].equals("D")) { // NOT present in read sequence, but does in
                                                     // reference sequence.
          if ((outputArrayPointer + lengthOfCurrentSegment) <= 0) {
            currentPosition += lengthOfCurrentSegment;
            outputArrayPointer += lengthOfCurrentSegment;
          } else {
            if (outputArrayPointer < 0) {
              currentPosition -= outputArrayPointer;
              lengthOfCurrentSegment += outputArrayPointer;
              outputArrayPointer = 0;
            }

            loop = outputArrayPointer
                   + Math.min(lengthOfCurrentSegment, outputArrayLength - outputArrayPointer);
            while (outputArrayPointer < loop) {
              output1_ReadsCounts[outputArrayPointer][INDEX_OF_DEL]++;
              output2_PhredScores[outputArrayPointer][INDEX_OF_DEL] +=
                                                                    DEFAULT_PHRED_SCORE_FOR_DELETION;
              output3_MappingScores[outputArrayPointer][INDEX_OF_DEL] += currentMappingScore;
              outputArrayPointer++;
            }

            currentPosition += lengthOfCurrentSegment;
          }

        } else if (readSegments[0][i].equals("N")) { // NOT present in read sequence, but does in
                                                     // reference sequence. Similar to D.
          currentPosition += lengthOfCurrentSegment;
          outputArrayPointer += lengthOfCurrentSegment;

        } else if (readSegments[0][i].equals("S")) { // present in read sequence, and also in
                                                     // reference sequence, but do not match. Lower
                                                     // case letters in read sequence.
          if (i == 0) {
            readPointer += lengthOfCurrentSegment;
          } else {
            currentPosition += lengthOfCurrentSegment;
            outputArrayPointer += lengthOfCurrentSegment;
            readPointer += lengthOfCurrentSegment;
          }

        } else if (readSegments[0][i].equals("H")) { // NOT present in read sequence, but in
                                                     // reference sequence. Similar to D.
          currentPosition += lengthOfCurrentSegment;
          outputArrayPointer += lengthOfCurrentSegment;

        } else if (readSegments[0][i].equals("P")) { // present in read sequence, but NOT in
                                                     // reference sequence. Similar to I.

        } else {
          System.err.println("Unrecognized CIGAR string: " + readSegments[0][i]);
        }
      }
    }
  }

  public static int convertToPhredScore(char asciiString) {
    if (asciiString == '*') {
      return -1;
    } else {
      return asciiString - 33;
    }
  }

  public static int[] convertToPhredScore(String asciiString) {
    int[] result;
    int stringLength;

    if (asciiString.equals("*")) {
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

  public static void getFilteredResults(int[][][] readsCounts, int[][][] mappingScores,
                                        int inputArrayStartIndex, int inputArrayStopIndex,
                                        Vector<int[]> output1DenovoMutations) {
    int[][] orderedIndices;
    int[] temp;

    orderedIndices = new int[readsCounts.length][];
    for (int i = inputArrayStartIndex; i <= inputArrayStopIndex; i++) {
      for (int j = 0; j < orderedIndices.length; j++) {
        orderedIndices[j] = Sort.quicksort(
                                           new int[] {readsCounts[j][i][0], readsCounts[j][i][1],
                                                      readsCounts[j][i][2], readsCounts[j][i][3]},
                                           Sort.DESCENDING);
      }
      updateNumAlleles(readsCounts, orderedIndices, i);

      if (filterVariants(readsCounts, orderedIndices, i)
      // if (filterReadDepth(readsCounts, orderedIndices, i, minRead)
      // && filterReadsCount(readsCounts, orderedIndices, i)

      // && filter3Allelic(readsCounts, orderedIndices, i)
      // && filterMappingQuality(mappingScores, readsCounts, orderedIndices, i)
      // && filterInDel(readsCounts, orderedIndices, i)
      // && filterMismatch(readsCounts, orderedIndices, i)
      ) {

        temp = new int[DENOVO_MUTATION_ARRAY_STRUCT.length];
        temp[0] = i;
        if ((readsCounts[1][i][orderedIndices[0][0]] == 0
             && readsCounts[2][i][orderedIndices[0][0]] == 0)
            || (readsCounts[1][i][orderedIndices[0][1]] == 0
                && readsCounts[2][i][orderedIndices[0][1]] == 0)) {
          temp[1] = 100;
        } else {
          temp[1] = 80;
        }
        for (int j = 0; j < readsCounts.length; j++) {
          for (int k = 0; k < 7; k++) {
            temp[j + 2] += readsCounts[j][i][k];
          }
          for (int k = 0; k < 2; k++) {
            if (readsCounts[j][i][orderedIndices[j][k]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
              temp[2 * j + 5 + k] = orderedIndices[j][k];
            } else {
              temp[2 * j + 5 + k] = -1;
            }
          }
        }
        output1DenovoMutations.add(temp);
      }
    }
  }

  public static boolean filterReadDepth(int[][][] readsCounts, int[][] orderedIndices, int i,
                                        int minRead) {
    return (readsCounts[0][i][orderedIndices[0][0]] >= minRead
            && readsCounts[1][i][orderedIndices[1][0]] >= minRead
            && readsCounts[2][i][orderedIndices[2][0]] >= minRead
            && (readsCounts[0][i][orderedIndices[0][0]]
                + readsCounts[0][i][orderedIndices[0][1]]) >= MIN_READ_DEPTH
            && (readsCounts[1][i][orderedIndices[1][0]]
                + readsCounts[1][i][orderedIndices[1][1]]) >= MIN_READ_DEPTH
            && (readsCounts[2][i][orderedIndices[2][0]]
                + readsCounts[2][i][orderedIndices[2][1]]) >= MIN_READ_DEPTH
    // && readsCounts[0][i][INDEX_OF_TOTAL_READS] >= MIN_READ_DEPTH
    // && readsCounts[1][i][INDEX_OF_TOTAL_READS] >= MIN_READ_DEPTH
    // && readsCounts[2][i][INDEX_OF_TOTAL_READS] >= MIN_READ_DEPTH
    );
  }

  public static boolean filterReadsCount(int[][][] readsCounts, int[][] orderedIndices, int i) {
    return ((readsCounts[0][i][orderedIndices[0][1]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
             && readsCounts[1][i][orderedIndices[0][1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
             && readsCounts[2][i][orderedIndices[0][1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
             // && readsCounts[0][i][orderedIndices[0][1]] > .18 *
             // readsCounts[0][i][orderedIndices[0][0]]
             // && readsCounts[0][i][orderedIndices[0][1]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5
             // * (readsCounts[1][i][orderedIndices[0][1]] +
             // readsCounts[2][i][orderedIndices[0][1]]))
             && readsCounts[0][i][orderedIndices[0][1]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO
                                                           + readsCounts[1][i][orderedIndices[0][1]]
                                                           + readsCounts[2][i][orderedIndices[0][1]])
    // && (readsCounts[1][i][orderedIndices[0][1]] * 19) < (readsCounts[1][i][orderedIndices[0][0]]
    // + readsCounts[1][i][orderedIndices[0][2]] + readsCounts[1][i][orderedIndices[0][3]])
    // && (readsCounts[2][i][orderedIndices[0][1]] * 19) < (readsCounts[2][i][orderedIndices[0][0]]
    // + readsCounts[2][i][orderedIndices[0][2]] + readsCounts[2][i][orderedIndices[0][3]])
    ) || (readsCounts[0][i][orderedIndices[0][0]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
          && readsCounts[1][i][orderedIndices[0][0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
          && readsCounts[2][i][orderedIndices[0][0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
          // && readsCounts[0][i][orderedIndices[0][0]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 *
          // (readsCounts[1][i][orderedIndices[0][0]] + readsCounts[2][i][orderedIndices[0][0]]))
          && readsCounts[0][i][orderedIndices[0][0]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO
                                                        + readsCounts[1][i][orderedIndices[0][0]]
                                                        + readsCounts[2][i][orderedIndices[0][0]])
    // && (readsCounts[1][i][orderedIndices[0][0]] * 19) < (readsCounts[1][i][orderedIndices[0][1]]
    // + readsCounts[1][i][orderedIndices[0][2]] + readsCounts[1][i][orderedIndices[0][3]])
    // && (readsCounts[2][i][orderedIndices[0][0]] * 19) < (readsCounts[2][i][orderedIndices[0][1]]
    // + readsCounts[2][i][orderedIndices[0][2]] + readsCounts[2][i][orderedIndices[0][3]])
    ) || (readsCounts[0][i][INDEX_OF_INS] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
          && readsCounts[1][i][INDEX_OF_INS] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
          && readsCounts[2][i][INDEX_OF_INS] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
          // && readsCounts[0][i][INDEX_OF_INS] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 *
          // (readsCounts[1][i][INDEX_OF_INS] + readsCounts[2][i][INDEX_OF_INS]))
          && readsCounts[0][i][INDEX_OF_INS] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO
                                                + readsCounts[1][i][INDEX_OF_INS]
                                                + readsCounts[2][i][INDEX_OF_INS]))
            || (readsCounts[0][i][INDEX_OF_DEL] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
                && readsCounts[1][i][INDEX_OF_DEL] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                && readsCounts[2][i][INDEX_OF_DEL] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                // && readsCounts[0][i][INDEX_OF_DEL] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 *
                // (readsCounts[1][i][INDEX_OF_DEL] + readsCounts[2][i][INDEX_OF_DEL]))
                && readsCounts[0][i][INDEX_OF_DEL] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO
                                                      + readsCounts[1][i][INDEX_OF_DEL]
                                                      + readsCounts[2][i][INDEX_OF_DEL])));
  }

  public static boolean filterVariants(int[][][] readsCounts, int[][] orderedIndices, int i) {
    return (readsCounts[0][i][orderedIndices[0][1]] >= MAX_ALLELE_COUNT_TREATED_AS_ZERO
            // || readsCounts[1][i][orderedIndices[1][1]] >= MAX_ALLELE_COUNT_TREATED_AS_ZERO
            // || readsCounts[2][i][orderedIndices[2][1]] >= MAX_ALLELE_COUNT_TREATED_AS_ZERO
            && (readsCounts[0][i][orderedIndices[0][0]]
                + readsCounts[0][i][orderedIndices[0][1]]) >= MIN_READ_DEPTH);
  }

  public static boolean filterMismatch(int[][][] readsCounts, int[][] orderedIndices, int i) {
    return ((readsCounts[0][i][INDEX_OF_N] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
             || (readsCounts[0][i][orderedIndices[0][0]]
                 + readsCounts[0][i][orderedIndices[0][1]]) >= 20 * readsCounts[0][i][INDEX_OF_N])
            && (readsCounts[1][i][INDEX_OF_N] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                || (readsCounts[1][i][orderedIndices[0][0]]
                    + readsCounts[1][i][orderedIndices[0][1]]) >= 20
                                                                  * readsCounts[1][i][INDEX_OF_N])
            && (readsCounts[2][i][INDEX_OF_N] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                || (readsCounts[2][i][orderedIndices[0][0]]
                    + readsCounts[2][i][orderedIndices[0][1]]) >= 20
                                                                  * readsCounts[2][i][INDEX_OF_N]));
  }

  public static boolean filterInDel(int[][][] readsCounts, int[][] orderedIndices, int i) {
    int[] numInsDelMismatch;

    numInsDelMismatch =
                      new int[] {readsCounts[0][i][INDEX_OF_INS] + readsCounts[0][i][INDEX_OF_DEL],
                                 readsCounts[1][i][INDEX_OF_INS] + readsCounts[1][i][INDEX_OF_DEL],
                                 readsCounts[2][i][INDEX_OF_INS] + readsCounts[2][i][INDEX_OF_DEL]};
    return ((numInsDelMismatch[0] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
             || (readsCounts[0][i][orderedIndices[0][0]]
                 + readsCounts[0][i][orderedIndices[0][1]]) >= 20 * numInsDelMismatch[0])
            && (numInsDelMismatch[1] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                || (readsCounts[1][i][orderedIndices[0][0]]
                    + readsCounts[1][i][orderedIndices[0][1]]) >= 20 * numInsDelMismatch[1])
            && (numInsDelMismatch[2] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                || (readsCounts[2][i][orderedIndices[0][0]]
                    + readsCounts[2][i][orderedIndices[0][1]]) >= 20 * numInsDelMismatch[2]));
  }

  public static boolean filterThreeAllelic(int[][][] readsCounts, int[][] orderedIndices, int i) {
    return ((readsCounts[0][i][orderedIndices[0][2]]
             + readsCounts[0][i][orderedIndices[0][3]]) <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
            && (readsCounts[1][i][orderedIndices[1][2]]
                + readsCounts[1][i][orderedIndices[1][3]]) <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
            && (readsCounts[2][i][orderedIndices[2][2]]
                + readsCounts[2][i][orderedIndices[2][3]]) <= MAX_ALLELE_COUNT_TREATED_AS_ZERO);
  }

  public static boolean filterChildAltAlleleCount(int[][][] readsCounts, int[][] orderedIndices,
                                                  int i, int altAlleleCountThreshold) {
    return ((readsCounts[0][i][orderedIndices[0][1]]) >= altAlleleCountThreshold
    // && (readsCounts[1][i][orderedIndices[1][2]] + readsCounts[1][i][orderedIndices[1][3]]) <=
    // altAlleleCountThreshold
    // && (readsCounts[2][i][orderedIndices[2][2]] + readsCounts[2][i][orderedIndices[2][3]]) <=
    // altAlleleCountThreshold
    );
  }

  public static boolean filterMappingQuality(int[][][] mappingScores, int[][][] readsCounts,
                                             int[][] orderedIndices, int i) {
    return (mappingScores[0][i][orderedIndices[0][0]] >= (MIN_MAPPING_SCORE
                                                          * readsCounts[0][i][orderedIndices[0][0]])
            && mappingScores[1][i][orderedIndices[1][0]] >= (MIN_MAPPING_SCORE
                                                             * readsCounts[0][i][orderedIndices[1][0]])
            && mappingScores[2][i][orderedIndices[2][0]] >= (MIN_MAPPING_SCORE
                                                             * readsCounts[0][i][orderedIndices[2][0]])
            && mappingScores[0][i][orderedIndices[0][1]] >= (MIN_MAPPING_SCORE
                                                             * readsCounts[0][i][orderedIndices[0][1]])
            && mappingScores[1][i][orderedIndices[1][1]] >= (MIN_MAPPING_SCORE
                                                             * readsCounts[0][i][orderedIndices[1][1]])
            && mappingScores[2][i][orderedIndices[2][1]] >= (MIN_MAPPING_SCORE
                                                             * readsCounts[0][i][orderedIndices[2][1]]));
  }

  public static boolean filterPhred(int[][][] phredScores, int[][][] readsCounts,
                                    int[][] orderedIndices, int i) {
    return (phredScores[0][i][orderedIndices[0][0]] >= (MIN_PHRED
                                                        * readsCounts[0][i][orderedIndices[0][0]])
            && phredScores[1][i][orderedIndices[1][0]] >= (MIN_PHRED
                                                           * readsCounts[0][i][orderedIndices[1][0]])
            && phredScores[2][i][orderedIndices[2][0]] >= (MIN_PHRED
                                                           * readsCounts[0][i][orderedIndices[2][0]])
            && phredScores[0][i][orderedIndices[0][1]] >= (MIN_PHRED
                                                           * readsCounts[0][i][orderedIndices[0][1]])
            && phredScores[1][i][orderedIndices[1][1]] >= (MIN_PHRED
                                                           * readsCounts[0][i][orderedIndices[1][1]])
            && phredScores[2][i][orderedIndices[2][1]] >= (MIN_PHRED
                                                           * readsCounts[0][i][orderedIndices[2][1]]));
  }

  public static boolean filterNearbyVarInDelDenovo(int[] numVarInDelDenovoNearby) {
    return ((numVarInDelDenovoNearby[0] + numVarInDelDenovoNearby[1] + numVarInDelDenovoNearby[2]
             + numVarInDelDenovoNearby[3]) <= MAX_NUM_VAR_INS_DEL_DENOVO_NEARBY);
  }

  // public static void adjDenovoMutationScoresForMismatches(int[][][] readsCounts, int
  // inputArrayMarkerIndex, byte[] output1DenovoMutationScores, String[] output2DenovoMutationNotes)
  // {
  // String note;
  //
  // note = "";
  // for (int j = 0; j < SAMPLE_SUFFIX.length; j++) {
  // if (readsCounts[j][inputArrayMarkerIndex][INDEX_OF_N] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
  // if (note.equals("")) {
  // note = SAMPLE_SUFFIX[j];
  // } else {
  // note += "," + SAMPLE_SUFFIX[j];
  // }
  // }
  // }
  //
  // if (! note.equals("")) {
  // output1DenovoMutationScores[inputArrayMarkerIndex] = (byte)
  // (output1DenovoMutationScores[inputArrayMarkerIndex] * DISCOUNT_FOR_N);
  //
  // if (output2DenovoMutationNotes[inputArrayMarkerIndex] == null) {
  // output2DenovoMutationNotes[inputArrayMarkerIndex] = note + " mismatch(es)";
  // } else {
  // output2DenovoMutationNotes[inputArrayMarkerIndex] += (";" + note + " mismatch(es)");
  // }
  // }
  // }

  // public static void adjDenovoMutationScoresForThreeAlleles(int[][][] readsCounts, int
  // indexCurrentMarker, byte[] output1DenovoMutationScores, String[] output2DenovoMutationNotes) {
  // String note;
  // double discountRate;
  //
  // note = "";
  // discountRate = 1.0;
  // updateNumAlleles(readsCounts, indexCurrentMarker);
  // for(int i = 0; i < SAMPLE_SUFFIX.length; i++) {
  // if(readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] > 2) {
  // if (note.equals("")) {
  // note = SAMPLE_SUFFIX[i];
  // } else {
  // note += "," + SAMPLE_SUFFIX[i];
  // }
  //
  // if (readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] > 2 ||
  // readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] > 3) {
  // discountRate = DISCOUNT_FOR_3ALLELES_LOOSE;
  // } else {
  // discountRate = DISCOUNT_FOR_3ALLELES_STRICT;
  // }
  // }
  // }
  //
  // if (! note.equals("")) {
  // output1DenovoMutationScores[indexCurrentMarker] *= discountRate;
  // if (output2DenovoMutationNotes[indexCurrentMarker] == null) {
  // output2DenovoMutationNotes[indexCurrentMarker] = note + ":3 alleles";
  // } else {
  // output2DenovoMutationNotes[indexCurrentMarker] += (note + ":3 alleles");
  // }
  // }
  // }

  // (#A=10, #T=3, #G=2, #C=1) is counted as 3 allelic, though #G<=2 and #C<=2
  public static void updateNumAlleles(int[][][] readsCounts, int[][] orderedIndices,
                                      int indexCurrentMarker) {
    int sum;
    for (int i = 0; i < SAMPLE_SUFFIX.length; i++) {
      // for (int j = orderedIndices[0].length - 1; j >= 0; j--) {
      // if (readsCounts[i][indexCurrentMarker][orderedIndices[i][j]] > 0) {
      // if (readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] == 0) {
      // if (readsCounts[i][indexCurrentMarker][orderedIndices[i][j]] >
      // MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
      // readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] = j + 1;
      // readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] = j + 1;
      // break;
      // } else {
      // readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] = j + 1;
      // }
      // } else {
      // if (readsCounts[i][indexCurrentMarker][orderedIndices[i][j]] >
      // MAX_ALLELE_COUNT_TREATED_AS_ZERO ||
      // (readsCounts[i][indexCurrentMarker][orderedIndices[i][j]] +
      // readsCounts[i][indexCurrentMarker][orderedIndices[i][j+1]]) >
      // MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
      // readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] = j + 1;
      // break;
      // }
      // }
      // }
      // }

      readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] = 4;
      readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] = 0;
      for (int j = 0; j < orderedIndices[0].length; j++) {
        if (readsCounts[i][indexCurrentMarker][orderedIndices[i][j]] == 0) {
          readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] = j;
          break;
        }
      }
      sum = 0;
      for (int j =
                 readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_STRICT] - 1; j >= 0; j--) {
        sum += readsCounts[i][indexCurrentMarker][orderedIndices[i][j]];
        if (sum > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
          readsCounts[i][indexCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE] = j + 1;
          break;
        }
      }
    }
  }

  public static void getNearbyInDelVars(int[][][] readsCounts, Vector<int[]> denovoMutations,
                                        Vector<Vector<Integer>[][]> outputDenovoMutationNotes) {
    byte[] indicesOfInDels;
    int index;
    int loop2;
    int loop1;
    Vector<Integer>[][] temp;

    indicesOfInDels = new byte[] {INDEX_OF_INS, INDEX_OF_DEL, INDEX_OF_NUM_ALLELES_LOOSE};
    loop1 = denovoMutations.size();
    for (int i = 0; i < loop1; i++) {
      index = denovoMutations.elementAt(i)[0];
      temp = outputDenovoMutationNotes.elementAt(i);
      loop2 = Math.min(readsCounts[0].length, index + WINDOW_SIZE_FOR_NEARBY_INDEL);
      for (int j = Math.max(0, index - WINDOW_SIZE_FOR_NEARBY_INDEL); j < loop2; j++) {
        for (int k = 0; k < readsCounts.length; k++) {
          for (int l = 0; l < DENOVO_MUTATION_NOTES_ARRAY_STRUCT.length; l++) {
            if ((l < 2 && readsCounts[k][j][indicesOfInDels[l]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO)
                || (l == 2 && readsCounts[k][j][indicesOfInDels[l]] == 2)
                || (l == 3 && readsCounts[k][j][INDEX_OF_NUM_ALLELES_LOOSE] > 2)) {
              if (temp[l][k] == null) {
                temp[l][k] = new Vector<Integer>();
              }
              temp[l][k].add(j - index);
            }
          }
        }
      }
    }
  }

  public static void getNearbyDeNovos(int[][][] readsCounts, Vector<int[]> denovoMutations,
                                      Vector<Vector<Integer>[][]> outputDenovoMutationNotes) {
    int index2;
    int loop1;
    Vector<Integer>[][] temp;
    int[] currentDenovoMuations;
    int[] orderedIndices;

    loop1 = denovoMutations.size();
    for (int i = 0; i < loop1; i++) {
      currentDenovoMuations = denovoMutations.elementAt(i);
      temp = outputDenovoMutationNotes.elementAt(i);
      for (int j = 0; j < loop1; j++) {
        if (i != j) {
          index2 = denovoMutations.elementAt(j)[0];
          if (Math.abs(currentDenovoMuations[0] - index2) <= WINDOW_SIZE_FOR_NEARBY_DENOVO
              && currentDenovoMuations[3] > 10 && currentDenovoMuations[4] > 10) {
            orderedIndices = Sort.quicksort(new int[] {readsCounts[0][index2][0],
                                                       readsCounts[0][index2][1],
                                                       readsCounts[0][index2][2],
                                                       readsCounts[0][index2][3]},
                                            Sort.DESCENDING);
            if ((readsCounts[0][index2][orderedIndices[0]] >= 5
                 && readsCounts[1][index2][orderedIndices[0]] <= 1
                 && readsCounts[2][index2][orderedIndices[0]] <= 1)
                || (readsCounts[0][index2][orderedIndices[1]] >= 5
                    && readsCounts[1][index2][orderedIndices[1]] <= 1
                    && readsCounts[2][index2][orderedIndices[1]] <= 1)) {
              if (temp[4][0] == null) {
                temp[4][0] = new Vector<Integer>();
              }
              temp[4][0].add(currentDenovoMuations[0] - index2); // temp[indexOf_DenovoCounts][indexOfChild].add(distanceToTheNearbyDenovo)
            }
          }
        }
      }
    }
  }

  // public static boolean[] isThreeAlleles(int[][][] readsCounts, int indexCurrentMarker, int
  // thresholdFor3Alleles) {
  // boolean[] isThreeAlleles;
  // int DenovoMutationScore;
  //
  // isThreeAlleles = new boolean[readsCounts.length];
  // for(int i = 0; i < readsCounts.length; i++) {
  // DenovoMutationScore = 0;
  // for (int j = 0; j < BASES.length; j++) {
  // if (readsCounts[i][indexCurrentMarker][j] > thresholdFor3Alleles) {
  // DenovoMutationScore ++;
  // }
  // }
  // if(DenovoMutationScore > 2) {
  // isThreeAlleles[i] = true;
  // }
  // }
  //
  // return isThreeAlleles;
  // }

  public static byte adjDeNovoMutationScore(byte oldScore, double discountRate) {
    oldScore *= discountRate;
    if (oldScore == 0) {
      oldScore = 1;
    }
    return oldScore;
  }

  public static double getFisherExcat(int[][] contingencyTable) {
    return getFactorial(contingencyTable[0][0] + contingencyTable[0][1]);
  }

  public static int getFactorial(int a) {
    int result;

    if (a == 0) {
      result = 0;
    } else {
      result = 1;
      for (int i = 2; i <= a; i++) {
        result *= i;
      }
    }

    return result;
  }

  public static String[] getRefFromFasta(String refFastaFilename, Vector<int[]> denovoMutations,
                                         String chr, int startPosAdjForWindow) {
    int loop;
    int[] temp;
    Process p;
    BufferedReader reader;
    // BufferedReader error;
    String line;
    int pos;
    String[] refs;

    loop = denovoMutations.size();
    refs = new String[loop];
    try {
      for (int i = 0; i < loop; i++) {
        temp = denovoMutations.elementAt(i);
        pos = temp[0] + startPosAdjForWindow;
        p =
          Runtime.getRuntime()
                 .exec("samtools faidx " + refFastaFilename + " chr" + chr + ":" + pos + "-" + pos,
                       null, new File(ext.parseDirectoryOfFile(refFastaFilename)));

        // ps=new ProcessBuilder("samtools", "faidx", refFastaFilename, "chr" + chr + ":" + start +
        // "-" + stop);
        // ps.redirectErrorStream(true);
        // ps.directory(new File(ext.parseDirectoryOfFile(refFastaFilename)));
        // p = ps.start();

        // p.waitFor();
        reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
        line = reader.readLine();
        while ((line = reader.readLine()) != null) {
          refs[i] = line;
        }
        p.waitFor();
        reader.close();

        // error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
        // while ((line = error.readLine()) != null) {
        // System.out.println(line);
        // }
        // error.close();

        displayErrorStream(p);
      }
    } catch (IOException e) {
      e.printStackTrace();
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    return refs;
  }

  public static void exportResult(PrintWriter writer, String id, String chr,
                                  int startPosAdjForWindow, int startPos,
                                  String[] refsForDeNovoMutations,
                                  Hashtable<String, String> altAllelesForInsertionForChild,
                                  int[][][] readsCounts, int[][][] phredScores,
                                  int[][][] mappingScores, Vector<int[]> denovoMutations,
                                  Vector<Vector<Integer>[][]> denovoMutationNotes, String flag,
                                  byte format, Logger log) {
    int[] deNovoMutation;
    int loop;
    int l;
    int currentPos;
    String altAlleleForInsertionForChild;

    loop = denovoMutations.size();
    for (int i = 0; i < loop; i++) {
      deNovoMutation = denovoMutations.elementAt(i);
      l = deNovoMutation[0];
      currentPos = startPosAdjForWindow + l;
      altAlleleForInsertionForChild = altAllelesForInsertionForChild.get(currentPos + "");
      writer.println(id + "\t" + chr + "\t" + currentPos + "\tchr" + chr + ":" + currentPos + "\t"
                     + refsForDeNovoMutations[i] + "\t"
                     + formatAltAllele(deNovoMutation, readsCounts, l, refsForDeNovoMutations[i],
                                       altAlleleForInsertionForChild, log)
                     + "\t" + ext.formDeci(deNovoMutation[1] / (double) 100, 2) + "\t"
                     + formatAlleleCounts(readsCounts, l, altAlleleForInsertionForChild, format)
                     + "\t" + formatNotes(denovoMutationNotes, i, format) + "\t"
                     + formatForwardGenotypes(deNovoMutation, refsForDeNovoMutations[i]) + "\t"
                     + flag + "\t" + deNovoMutation[2] + "\t" + deNovoMutation[3] + "\t"
                     + deNovoMutation[4] + "\t" + formatPhred(phredScores, readsCounts, l, format)
                     + "\t" + formatMapping(mappingScores, deNovoMutation, l, format));
    }
    // writer.flush();
  }

  public static String formatAlleleCounts(int[][][] readsCounts, int indexOfCurrentMarker,
                                          String altAlleleForInsertionForChild, byte format) {
    String result;
    int numAlleles;

    result = "";
    if (format == 2) {
      for (int i = 0; i < readsCounts.length; i++) {
        for (int j = 0; j < BASES_WITH_N_INDEL.length; j++) {
          result +=
                 ((i == 0 && j == 0 ? "" : "\t")
                  + (readsCounts[i][indexOfCurrentMarker][j] == 0 ? ""
                                                                  : readsCounts[i][indexOfCurrentMarker][j]));
        }
      }
      for (int i = 0; i < readsCounts.length; i++) {
        numAlleles = readsCounts[i][indexOfCurrentMarker][INDEX_OF_NUM_ALLELES_LOOSE];
        if (readsCounts[i][indexOfCurrentMarker][INDEX_OF_INS] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
          if (i == 0) {
            if (altAlleleForInsertionForChild != null) {
              numAlleles += altAlleleForInsertionForChild.split(",").length;
            }
          } else {
            numAlleles += 1;
          }
        }
        if (readsCounts[i][indexOfCurrentMarker][INDEX_OF_DEL] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
          numAlleles += 1;
        }
        result += ("\t" + (numAlleles == 0 ? "" : numAlleles));
      }
    } else {
      for (int i = 0; i < readsCounts.length; i++) {
        result += (i == 0 ? "(" : " (");
        for (int j = 0; j < BASES.length; j++) {
          result +=
                 ((j == 0 ? "" : ",")
                  + (readsCounts[i][indexOfCurrentMarker][j] == 0 ? "-"
                                                                  : readsCounts[i][indexOfCurrentMarker][j]));
        }
        result += "/";
        for (int j = BASES.length; j < BASES_WITH_N_INDEL.length; j++) {
          result +=
                 ((j == BASES.length ? "" : ",")
                  + (readsCounts[i][indexOfCurrentMarker][j] == 0 ? "-"
                                                                  : readsCounts[i][indexOfCurrentMarker][j]));
        }
        result += ")";
      }
    }

    return result;
    // return (readsCounts[0][indexOfCurrentMarker][0]==0? "-" :
    // readsCounts[0][indexOfCurrentMarker][0]) + "," + (readsCounts[0][indexOfCurrentMarker][1]==0?
    // "-" : readsCounts[0][indexOfCurrentMarker][1]) + "," +
    // (readsCounts[0][indexOfCurrentMarker][2]==0? "-" : readsCounts[0][indexOfCurrentMarker][2]) +
    // "," + (readsCounts[0][indexOfCurrentMarker][3]==0? "-" :
    // readsCounts[0][indexOfCurrentMarker][3]) + "/" + (readsCounts[0][indexOfCurrentMarker][4]==0?
    // "-" : readsCounts[0][indexOfCurrentMarker][4]) + "," +
    // (readsCounts[0][indexOfCurrentMarker][5]==0? "-" : readsCounts[0][indexOfCurrentMarker][5]) +
    // "," + (readsCounts[0][indexOfCurrentMarker][6]==0? "-" :
    // readsCounts[0][indexOfCurrentMarker][6]) + ") (" +
    // (readsCounts[1][indexOfCurrentMarker][0]==0? "-" : readsCounts[1][indexOfCurrentMarker][0]) +
    // "," + (readsCounts[1][indexOfCurrentMarker][1]==0? "-" :
    // readsCounts[1][indexOfCurrentMarker][1]) + "," + (readsCounts[1][indexOfCurrentMarker][2]==0?
    // "-" : readsCounts[1][indexOfCurrentMarker][2]) + "," +
    // (readsCounts[1][indexOfCurrentMarker][3]==0? "-" : readsCounts[1][indexOfCurrentMarker][3]) +
    // "/" + (readsCounts[1][indexOfCurrentMarker][4]==0? "-" :
    // readsCounts[1][indexOfCurrentMarker][4]) + "," + (readsCounts[1][indexOfCurrentMarker][5]==0?
    // "-" : readsCounts[1][indexOfCurrentMarker][5]) + "," +
    // (readsCounts[1][indexOfCurrentMarker][6]==0? "-" : readsCounts[1][indexOfCurrentMarker][6]) +
    // ") (" + (readsCounts[2][indexOfCurrentMarker][0]==0? "-" :
    // readsCounts[2][indexOfCurrentMarker][0]) + "," + (readsCounts[2][indexOfCurrentMarker][1]==0?
    // "-" : readsCounts[2][indexOfCurrentMarker][1]) + "," +
    // (readsCounts[2][indexOfCurrentMarker][2]==0? "-" : readsCounts[2][indexOfCurrentMarker][2]) +
    // "," + (readsCounts[2][indexOfCurrentMarker][3]==0? "-" :
    // readsCounts[2][indexOfCurrentMarker][3]) + "/" + (readsCounts[2][indexOfCurrentMarker][4]==0?
    // "-" : readsCounts[2][indexOfCurrentMarker][4]) + "," +
    // (readsCounts[2][indexOfCurrentMarker][5]==0? "-" : readsCounts[2][indexOfCurrentMarker][5]) +
    // "," + (readsCounts[2][indexOfCurrentMarker][6]==0? "-" :
    // readsCounts[2][indexOfCurrentMarker][6]) + ")";
  }

  public static String formatNotes(Vector<Vector<Integer>[][]> denovoMarkerNotes,
                                   int indexOfCurrentElement, byte format) {
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
              result2 +=
                      ((isFirstSection ? "" : " ") + DENOVO_MUTATION_NOTES_ARRAY_STRUCT[i] + "(");
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
            result1 += ((i == 0 && j == 0 ? "" : "\t") + temp[i][j].size());
          } else {
            result1 += (i == 0 && j == 0 ? "" : (i == 4 && j > 0 ? "" : "\t"));
          }
        }
        if (!isFirstInTheSection) {
          result2 += ")";
        }
      }
    } else {
      for (int i = 0; i < temp.length; i++) {
        isFirstInTheSection = true;
        result1 += ((result1.equals("") ? "" : " ") + "(");
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
            result1 += ((j == 0 ? "" : ",") + temp[i][j].size());
          } else {
            result1 += ((j == 0 ? "" : ",") + (i == 4 && j > 0 ? "" : "-"));
          }
        }
        if (!isFirstInTheSection) {
          result2 += ")";
        }
        result1 += ")";
      }
    }

    return result1 + "\t" + result2;
  }


  /*
   * Definition of alt allele: alt = The one of the 2 most frequent alleles that is not equal to
   * ref; alt = "", if the most frequent allele has a read count <= 2 or (the most frequent allele
   * == ref && the 2nd most frequent allele has a read count <= 2)
   */
  public static String formatAltAllele(int[] deNovoMutations, int[][][] readsCounts,
                                       int indexCurrentMarker, String ref,
                                       String altAllelesForInsertions, Logger log) {
    String result;

    if (deNovoMutations[5] < 0) { // "C_allele1", see DENOVO_MUTATION_ARRAY_STRUCT. ==-1 if
                                  // !(readsCounts[j][][orderedIndices[j][1]] >
                                  // MAX_ALLELE_COUNT_TREATED_AS_ZERO)
      result = "";
    } else {
      result = BASES[deNovoMutations[5]] + "";
      if (result.equalsIgnoreCase(ref)) {
        if (deNovoMutations[6] < 0) { // "C_allele2", see DENOVO_MUTATION_ARRAY_STRUCT. ==-1 if
                                      // !(readsCounts[j][][orderedIndices[j][1]] >
                                      // MAX_ALLELE_COUNT_TREATED_AS_ZERO)
          result = "";
        } else {
          result = BASES[deNovoMutations[6]] + "";
        }
      }
    }

    if (altAllelesForInsertions != null) {
      if (!result.equals("")) {
        result += ",";
      }
      result += altAllelesForInsertions;
    }

    if (readsCounts[0][indexCurrentMarker][INDEX_OF_DEL] > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
      if (!result.equals("")) {
        result += ",";
      }
      result += "-" + ref;
    }

    return result;
  }

  public static String formatForwardGenotypes(int[] temp, String ref) {
    String allele1;
    String result;

    result = "";
    for (int i = 0; i < SAMPLE_SUFFIX.length; i++) {
      result += (i == 0 ? "" : ",");
      if (temp[5 + 2 * i] == -1) {
        allele1 = "";
      } else {
        allele1 = BASES[temp[5 + 2 * i]] + "";
      }
      if (temp[6 + 2 * i] < 0) {
        result += (allele1 + allele1);
      } else {
        if (allele1.equalsIgnoreCase(ref)) {
          result += (allele1 + BASES[temp[6 + 2 * i]]);
        } else {
          result += (BASES[temp[6 + 2 * i]] + allele1);
        }
      }
    }

    return result;
  }

  public static String formatPhred(int[][][] phredScores, int[][][] readsCounts,
                                   int indexOfCurrentMarker, byte format) {
    double[][] phredScoreProportions;
    double totalPhredScores;
    String result;

    phredScoreProportions = new double[SAMPLE_SUFFIX.length][BASES.length];
    for (int j = 0; j < phredScoreProportions.length; j++) {
      totalPhredScores = phredScores[j][indexOfCurrentMarker][0]
                         + phredScores[j][indexOfCurrentMarker][1]
                         + phredScores[j][indexOfCurrentMarker][2]
                         + phredScores[j][indexOfCurrentMarker][3];
      for (int k = 0; k < phredScoreProportions[j].length; k++) {
        phredScoreProportions[j][k] = phredScores[j][indexOfCurrentMarker][k] / totalPhredScores;
      }
    }

    result = "";
    if (format == 2) {
      for (int i = 0; i < phredScoreProportions.length; i++) {
        for (int j = 0; j < phredScoreProportions[i].length; j++) {
          result += ((i == 0 && j == 0 ? "" : "\t")
                     + (phredScoreProportions[i][j] == 0 ? ""
                                                         : ext.formDeci(phredScoreProportions[i][j],
                                                                        3)));
        }
      }

      for (int i = 0; i < phredScores.length; i++) {
        for (int j = 0; j < phredScores[i][indexOfCurrentMarker].length; j++) {

          result +=
                 ("\t"
                  + (phredScores[i][indexOfCurrentMarker][j] == 0 ? ""
                                                                  : (phredScores[i][indexOfCurrentMarker][j]
                                                                     / readsCounts[i][indexOfCurrentMarker][j])));
        }
      }

    } else {
      for (int i = 0; i < phredScoreProportions.length; i++) {
        result += ((i == 0 ? "" : " ") + "(");
        for (int j = 0; j < phredScoreProportions[i].length; j++) {
          result += ((j == 0 ? "" : ",") + (phredScoreProportions[i][j] == 0 ? "-"
                                                                             : ext.formDeci(phredScoreProportions[i][j],
                                                                                            3)));
        }
        result += ")";
      }
    }

    return result;
  }

  public static String formatMapping(int[][][] mappingScores, int[] denovoMutation,
                                     int indexOfCurrentMarker, byte format) {
    int sum;
    String result;

    result = "";
    if (format == 2) {
      for (int i = 0; i < mappingScores.length; i++) {
        if (denovoMutation[i + 2] == 0) {
          result += (i == 0 ? "" : "\t");
        } else {
          sum = 0;
          for (int j = 0; j < BASES.length; j++) {
            sum += mappingScores[i][indexOfCurrentMarker][j];
          }
          result += ((i == 0 ? "" : "\t") + (sum / denovoMutation[i + 2]));
        }
      }
    } else {
      for (int i = 0; i < mappingScores.length; i++) {
        result += (i == 0 ? "(" : ",");
        if (denovoMutation[i + 2] == 0) {
          result += "-";
        } else {
          sum = 0;
          for (int j = 0; j < BASES.length; j++) {
            sum += mappingScores[i][indexOfCurrentMarker][j];
          }
          result += (sum / denovoMutation[i + 2]);
        }
      }
      result += ")";
    }

    return result;
  }

  public static Vector<String[]> loadBamByPipeline(String fullpathToSamFile, String chr, int begin,
                                                   int end) {
    Process p;
    // ProcessBuilder ps;
    BufferedReader reader;
    // BufferedReader error;
    String line;
    Vector<String[]> samContentVec;

    samContentVec = new Vector<String[]>();

    try {
      p = Runtime.getRuntime().exec("samtools view " + ext.rootOf(fullpathToSamFile) + " chr" + chr
                                    + ":" + begin + "-" + end, null,
                                    new File(ext.parseDirectoryOfFile(fullpathToSamFile)));

      // ps=new ProcessBuilder("samtools", "view", bamFilenames[i], "chr" + chr + ":" + start + "-"
      // + stop);
      // ps.redirectErrorStream(true);
      // ps.directory(new File(bamDir));
      // p = ps.start();

      // p.waitFor();
      reader = new BufferedReader(new InputStreamReader(p.getInputStream()));

      while ((line = reader.readLine()) != null) {
        samContentVec.add(line.split("[\t]"));
      }
      p.waitFor();
      reader.close();

      displayErrorStream(p);
    } catch (IOException e) {
      e.printStackTrace();
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    return samContentVec;
  }

  public static Vector<String[]> loadBamByFileReader(String fullpathToSamFile, String chr,
                                                     int begin, int end) {
    BufferedReader reader;
    String line;
    Vector<String[]> samContentVec;

    samContentVec = new Vector<String[]>();
    try {
      reader = new BufferedReader(new FileReader(fullpathToSamFile));
      while ((line = reader.readLine()) != null) {
        samContentVec.add(line.split("[\t]"));
      }
      reader.close();
    } catch (IOException e) {
      e.printStackTrace();
    }

    return samContentVec;
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
    Hashtable<String, Hashtable<String, String[]>> result;

    resultFilenames = Files.list(resultDir, ".txt", false);
    trioIds = new String[resultFilenames.length];
    result = new Hashtable<String, Hashtable<String, String[]>>();
    for (int i = 0; i < trioIds.length; i++) {
      trioIds[i] = resultFilenames[i].split("_")[0];
      try {
        reader = new BufferedReader(new FileReader(resultDir + resultFilenames[i]));
        reader.readLine();
        while (reader.ready()) {
          line = reader.readLine().split("\t");
          if (Double.parseDouble(line[6]) >= callScoreThreshold) {
            if (!result.containsKey(line[3])) {
              result.put(line[3], new Hashtable<String, String[]>());
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

    saveParsedResults(resultDir + DEFAULT_OUTPUT_FILENAME, trioIds, result, null, null, (byte) 2,
                      PARSE_RESULT_HEADER_FORMAT2, log);
  }

  public static void parseResults(String resultDir, double callScoreThreshold, byte outFormat,
                                  Logger log) {
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
        while (reader.ready()) {
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

            orderedIndices = Sort.quicksort(
                                            new int[] {readsCounts[0][0], readsCounts[0][1],
                                                       readsCounts[0][2], readsCounts[0][3]},
                                            Sort.DESCENDING);
            numTmp = readsCounts[0][4] + readsCounts[0][5] + readsCounts[0][6];
            if (Integer.parseInt(line[14]) > 20 && Integer.parseInt(line[15]) > 20
                && Integer.parseInt(line[16]) > 20
                // && ((readsCounts[0][orderedIndices[1]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO &&
                // readsCounts[1][orderedIndices[1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO &&
                // readsCounts[2][orderedIndices[1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO &&
                // readsCounts[0][orderedIndices[1]] > .18 * readsCounts[0][orderedIndices[0]] &&
                // readsCounts[0][orderedIndices[1]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 *
                // (readsCounts[1][orderedIndices[1]] + readsCounts[2][orderedIndices[1]])) &&
                // (readsCounts[1][orderedIndices[1]] * 19 ) < (readsCounts[1][orderedIndices[0]] +
                // readsCounts[1][orderedIndices[2]] + readsCounts[1][orderedIndices[3]]) &&
                // (readsCounts[2][orderedIndices[1]] * 19 ) < (readsCounts[2][orderedIndices[0]] +
                // readsCounts[2][orderedIndices[2]] + readsCounts[2][orderedIndices[3]]))
                // || (readsCounts[0][orderedIndices[0]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO &&
                // readsCounts[1][orderedIndices[0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO &&
                // readsCounts[2][orderedIndices[0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO &&
                // readsCounts[0][orderedIndices[0]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5 *
                // (readsCounts[1][orderedIndices[0]] + readsCounts[2][orderedIndices[0]]))) &&
                // (readsCounts[1][orderedIndices[0]] * 19 ) < (readsCounts[1][orderedIndices[1]] +
                // readsCounts[1][orderedIndices[2]] + readsCounts[1][orderedIndices[3]]) &&
                // (readsCounts[2][orderedIndices[0]] * 19 ) < (readsCounts[2][orderedIndices[1]] +
                // readsCounts[2][orderedIndices[2]] + readsCounts[2][orderedIndices[3]]))
                && ((readsCounts[0][orderedIndices[1]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
                     && readsCounts[1][orderedIndices[1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                     && readsCounts[2][orderedIndices[1]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                     && readsCounts[0][orderedIndices[1]] > .18 * readsCounts[0][orderedIndices[0]]
                     && readsCounts[0][orderedIndices[1]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO + 2.5
                                                                                                * (readsCounts[1][orderedIndices[1]]
                                                                                                   + readsCounts[2][orderedIndices[1]])))
                    || (readsCounts[0][orderedIndices[0]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
                        && readsCounts[1][orderedIndices[0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                        && readsCounts[2][orderedIndices[0]] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                        && readsCounts[0][orderedIndices[0]] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO
                                                                + 2.5
                                                                  * (readsCounts[1][orderedIndices[0]]
                                                                     + readsCounts[2][orderedIndices[0]]))))
                && (readsCounts[0][orderedIndices[2]]
                    + readsCounts[0][orderedIndices[3]]) <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                && (numTmp <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                    || (readsCounts[0][orderedIndices[0]]
                        + readsCounts[0][orderedIndices[1]]) >= 20 * numTmp)) {

              isThreeAllelesOrInDel = false;
              altProportion = new double[3];
              for (int j = 0; j < altProportion.length; j++) {
                orderedIndices = Sort.quicksort(
                                                new int[] {readsCounts[j][0], readsCounts[j][1],
                                                           readsCounts[j][2], readsCounts[j][3]},
                                                Sort.DESCENDING);
                if ((readsCounts[j][orderedIndices[2]]
                     + readsCounts[j][orderedIndices[3]]) > MAX_ALLELE_COUNT_TREATED_AS_ZERO) {
                  isThreeAllelesOrInDel = true;
                  break;
                }
                numTmp = readsCounts[j][4] + readsCounts[j][5] + readsCounts[j][6];
                if (numTmp > MAX_ALLELE_COUNT_TREATED_AS_ZERO
                    && (readsCounts[j][orderedIndices[0]]
                        + readsCounts[j][orderedIndices[1]]) < 20 * numTmp) {
                  isThreeAllelesOrInDel = true;
                  break;
                }

                numTmp = ext.indexOfChar(line[4].toUpperCase().charAt(0), BASES);
                if (orderedIndices[0] != numTmp) {
                  altProportion[j] = readsCounts[j][orderedIndices[0]]
                                     / (double) (readsCounts[j][numTmp]
                                                 + readsCounts[j][orderedIndices[0]]);
                } else {
                  altProportion[j] = readsCounts[j][orderedIndices[1]]
                                     / (double) (readsCounts[j][numTmp]
                                                 + readsCounts[j][orderedIndices[1]]);
                }
              }

              if (!isThreeAllelesOrInDel) {
                if (!result.containsKey(line[3])) {
                  result.put(line[3], new HashMap<String, String[]>());
                }
                result.get(line[3])
                      .put(line[0],
                           new String[] {line[4], line[6], line[7], line[8], line[10], line[11],
                                         line[12], line[13], line[14], line[15], line[16],
                                         ext.formDeci(altProportion[0], 3),
                                         ext.formDeci(altProportion[1], 3),
                                         ext.formDeci(altProportion[2], 3)});
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

    filename = resultDir + DEFAULT_OUTPUT_FILENAME;
    // filename = "N:/statgen/OS_Logan/SuperNovo/"+DEFAULT_OUTPUT_FILENAME;
    try {
      keys = result.keySet().toArray();
      Arrays.sort(keys);
      writer = new PrintWriter(new FileWriter(filename));
      if (outFormat == 1) {
        writer.println("id\tchr\tposition\tlookup\tnumHits\tref\talt\tcall\tnote\tfwdGenotype\treadDepth_C\treadDepth_D\treadDepth_M\tPhredScores\tMappingQuality_C\tMappingQuality_D\tMappingQuality_M\taltAllele%_C\taltAllele%_D\taltAllele%_M");
        for (Object key2 : keys) {
          temp = result.get(key2);
          keySet = temp.entrySet();
          numTmp = temp.size();
          for (Entry<String, String[]> entry : keySet) {
            line = entry.getValue();
            key = (String) key2;
            writer.println(entry.getKey() + "\t" + key.substring(3, key.indexOf(":")) + "\t"
                           + key.substring(key.indexOf(":") + 1) + "\t" + key + "\t" + numTmp + "\t"
                           + line[0] + "\t"
                           + (line[3].substring(0, 1).equalsIgnoreCase(line[0]) ? line[3].charAt(1)
                                                                                : line[3].charAt(0))
                           + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4]
                           + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8]
                           + "\t" + line[9] + "\t" + line[10] + "\t" + line[11] + "\t" + line[12]
                           + "\t" + line[13]);
          }
        }
      } else if (outFormat == 2) {
        writer.print("chr\tposition\tlookup\tnumTrios");
        for (String trioId : trioIds) {
          writer.print("\t" + trioId);
        }
        writer.println();
        for (Object key2 : keys) {
          temp = result.get(key2);
          line = ((String) key2).split(":");
          writer.print(line[0] + "\t" + line[1] + "\t" + key2 + "\t" + temp.size());
          for (String trioId : trioIds) {
            if (temp.containsKey(trioId)) {
              line = temp.get(trioId);
              writer.print("\t" + line[0] + "; " + line[1]);
            } else {
              writer.print("\t");
            }
          }
          writer.println();
        }
      }

      writer.close();
      log.report("Summarized result is available at: " + filename);

    } catch (Exception e) {
      log.reportError("Error writing to " + filename);
      e.printStackTrace();
    }
  }

  public static void parseResults(String resultDir, String annotationDir, String miniSamScriptsDir,
                                  String miniSamDir, String fullPathToTrioNameList,
                                  byte resultFormat, double callScoreThreshold, byte outFormat,
                                  Logger log) {
    int index;
    int[] indices;
    int[][] orderedIndices;
    int[][][] readsCounts, phredScores, mappingScores;
    double[] pctDenovoAllele;
    String markerName, miniBamLink = null, miniSamScriptSubDir, trioId;
    String[] alleles, resultFilenames, annotation = null, resultElements, seatleSeekHeader, header,
        genes, line;
    BufferedReader reader;
    Vector<String> annotationNeedsVars, annotationNeedsInDels;
    CountHash geneCounts;
    CountHash geneTrioCounts;
    Set<String> geneSet;
    HashSet<String> miniSamHash;
    Hashtable<String, Vector<String>> miniSamNeeded;
    Hashtable<String, Hashtable<String, String[]>> result;
    Hashtable<String, String[]> annotationHash;
    HashSet<String> geneTrioCount;
    long time;

    time = new Date().getTime();
    if (!Files.exists(annotationDir)) {
      log.report("Creating new annotation directory in: " + annotationDir);
      new File(annotationDir).mkdirs();
    }
    annotationNeedsVars = new Vector<String>();
    annotationNeedsInDels = new Vector<String>();
    // miniBamNeeds = new Vector<String>();
    miniSamHash = Samtools.listFiles(miniSamScriptsDir, log);
    annotationHash = SeattleSeq.loadAllAnnotationInDir(annotationDir, log);
    miniSamNeeded = new Hashtable<String, Vector<String>>();
    resultFilenames = Files.list(resultDir, ".txt", false);
    result = new Hashtable<String, Hashtable<String, String[]>>();
    orderedIndices = new int[SAMPLE_SUFFIX.length][];
    seatleSeekHeader = Matrix.extractColumn(SeattleSeq.RELEVANTS, 0);
    geneCounts = new CountHash();
    geneTrioCounts = new CountHash();
    geneTrioCount = new HashSet<String>();
    for (String resultFilename : resultFilenames) {
      trioId = resultFilename.split("_")[0];
      log.report(ext.getTime() + "\tProcessing " + trioId);

      try {
        reader = new BufferedReader(new FileReader(resultDir + resultFilename));
        // header = reader.readLine().trim().split("[\\s]+");
        header = reader.readLine().trim().split("\t");
        indices = ext.indexFactors(COLUMNS_NEEDED_FROM_PHASE1_OUTPUT, header, false, true);
        while (reader.ready()) {
          line = reader.readLine().split("\t");
          if (Double.parseDouble(line[6]) >= callScoreThreshold) {
            readsCounts = new int[SAMPLE_SUFFIX.length][1][BASES_WITH_N_INDEL.length];
            phredScores = new int[SAMPLE_SUFFIX.length][1][BASES_WITH_N_INDEL.length];
            mappingScores = new int[SAMPLE_SUFFIX.length][1][BASES_WITH_N_INDEL.length];
            parseLine(line, resultFormat, readsCounts, phredScores, mappingScores);
            getNumVarsInDelsDenovos(line[44]);

            for (int j = 0; j < SAMPLE_SUFFIX.length; j++) {
              orderedIndices[j] = Sort.quicksort(new int[] {readsCounts[j][0][0],
                                                            readsCounts[j][0][1],
                                                            readsCounts[j][0][2],
                                                            readsCounts[j][0][3]},
                                                 Sort.DESCENDING);
            }
            if (true
            // if (filterMappingQuality(mappingScores, readsCounts, orderedIndices, 0)
            // && filter3Allelic(readsCounts, orderedIndices, 0)
            // && filterPhred(phredScores, readsCounts, orderedIndices, 0)
            // && filterNearbyVarInDelDenovo(numVarInDelDenovos)

            // && filterInDel(readsCounts, orderedIndices, 0)
            ) {

              // line[5] = formatAltAllele(new int[] {0, 1, 2, 3, 4, orderedIndices[0][0],
              // orderedIndices[0][1]}, readsCounts, 0, line[4], (line[5].startsWith("+")? line[5] :
              // null), log);
              if (line[5] == null || line[5].equals("")) {
                alleles = new String[] {line[4]};
              } else {
                alleles = line[5].split(",");
              }

              for (String allele : alleles) {
                if (allele.startsWith("-")) {
                  markerName = "chr" + line[1] + ":" + (Integer.parseInt(line[2]) - 1) + "_"
                               + line[4].toUpperCase() + "_" + allele.toUpperCase();
                } else {
                  markerName = "chr" + line[1] + ":" + line[2] + "_" + line[4].toUpperCase() + "_"
                               + allele.toUpperCase();
                }

                if (annotationHash != null && annotationHash.containsKey(markerName)) {
                  annotation = annotationHash.get(markerName);
                  // if (annotation[0].equals("1")) {
                  if (true) {
                    if (allele.startsWith("+")) {
                      index = 5; // TODO if there are two strings for this insertion, say "+T" and
                                 // "+TGC", then cannot tell how many reads are for the former and
                                 // how many the latter. Only have the number of the reads for the
                                 // combination of the two.
                    } else if (allele.startsWith("-")) {
                      index = 6;
                    } else {
                      index = ext.indexOfChar(allele.toUpperCase().charAt(0), BASES); // TODO maybe
                                                                                      // add a
                                                                                      // column in
                                                                                      // phase1's
                                                                                      // output:
                                                                                      // DeNovo
                                                                                      // Allele
                      if (!(readsCounts[0][0][index] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
                            && readsCounts[1][0][index] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                            && readsCounts[2][0][index] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                            && (readsCounts[0][0][index] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO
                                                            + readsCounts[1][0][index] + readsCounts[2][0][index])))) {
                        index = ext.indexOfChar(line[4].toUpperCase().charAt(0), BASES);
                      }
                    }

                    pctDenovoAllele = new double[SAMPLE_SUFFIX.length];
                    for (int k = 0; k < pctDenovoAllele.length; k++) {
                      pctDenovoAllele[k] = readsCounts[k][0][index]
                                           / (double) (readsCounts[k][0][0] + readsCounts[k][0][1]
                                                       + readsCounts[k][0][2] + readsCounts[k][0][3]
                                                       + readsCounts[k][0][5]
                                                       + readsCounts[k][0][6]);
                    }

                    if (pctDenovoAllele[0] >= .2) {
                      genes = annotation[7].split(",");
                      for (int k = 0; k < genes.length; k++) {
                        if (!genes[k].equals("none")) {
                          geneCounts.add(genes[k]);
                          if (!geneTrioCount.contains(genes[k] + "\t" + line[0])) { // genes[k] \t
                                                                                    // trioId
                            geneTrioCounts.add(genes[k]);
                            geneTrioCount.add(genes[k] + "\t" + line[0]); // genes[k] \t trioId
                          }
                        }

                        // TODO "none"?
                        miniSamScriptSubDir = miniSamScriptsDir + genes[k] + "/";
                        if (!new File(miniSamScriptSubDir).exists()) {
                          new File(miniSamScriptSubDir).mkdir();
                        }

                        // if (new File(miniSamSubDir).exists()) {
                        // miniSamHash = Samtools.listFiles(miniSamSubDir, log);
                        // TODO is the new region contained in the older region?
                        if (miniSamHash.contains(line[0] + "\t" + line[1] + "\t" + line[2])) { // trioId
                                                                                               // +
                                                                                               // "\t"
                                                                                               // +
                                                                                               // chr
                                                                                               // +
                                                                                               // "\t"
                                                                                               // +
                                                                                               // pos
                          // if (! new File(miniSamSubDir + line[0] + "_" + line[1] + "_" + line[2]
                          // + ".xml").exists()) {
                          // Files.write(Samtools.getIgvXmlScript(miniSamDir, line[1], line[2],
                          // miniSamFilenamesOfOneTrio), miniSamSubDir + line[0] + "_" + line[1] +
                          // "_" + line[2] + ".xml");
                          // }
                          if (!new File(miniSamScriptSubDir + line[0] + "_chr" + line[1] + "_"
                                        + line[2] + ".bat").exists()) {
                            Files.write(Samtools.getIgvLaunchScript(miniSamScriptSubDir + line[0]
                                                                    + "_" + line[1] + "_" + line[2]
                                                                    + ".xml"),
                                        miniSamScriptSubDir + line[0] + "_" + line[1] + "_"
                                                                               + line[2] + ".bat");
                          }
                          if (!new File(miniSamScriptSubDir + line[0] + "_chr" + line[1] + "_"
                                        + line[2] + ".xml").exists()) {
                            Files.write(Samtools.getIgvXmlScript(miniSamDir, line[1], line[2],
                                                                 new String[] {line[0],
                                                                               line[0] + "_chr"
                                                                                        + line[1]
                                                                                        + "_"
                                                                                        + line[2]
                                                                                        + "_C.bam",
                                                                               line[0] + "_chr" + line[1]
                                                                                                    + "_"
                                                                                                    + line[2]
                                                                                                    + "_D.bam",
                                                                               line[0] + "_chr" + line[1] + "_" + line[2]
                                                                                                                + "_M.bam"}),
                                        miniSamScriptSubDir + line[0] + "_chr" + line[1] + "_" + line[2] + ".xml");
                          }
                          miniBamLink = "=HYPERLINK(\"" + miniSamScriptSubDir + line[0] + "_chr"
                                        + line[1] + "_" + line[2] + ".bat\",\"IGV\")";
                          // miniBamLink = miniSamSubDir;
                        } else {
                          miniBamLink = "";
                          if (!miniSamNeeded.containsKey(genes[k])) {
                            miniSamNeeded.put(genes[k], new Vector<String>());
                          }
                          miniSamNeeded.get(genes[k])
                                       .add(line[0] + "\t" + line[1] + "\t" + line[2]); // trioId +
                                                                                        // chr + pos
                        }
                        // } else {
                        // miniBamLink = "";
                        // if (! miniSamNeeded.containsKey(genes[k])) {
                        // miniSamNeeded.put(genes[k], new Vector<String>());
                        // }
                        // miniSamNeeded.get(genes[k]).add(line[0] + "\t" + line[1] + "\t" +
                        // line[2]); //trioId + chr + pos
                        // }
                      }

                      resultElements = Array.subArray(line, indices);
                      resultElements = Array.concatAll(new String[] {annotation[7], miniBamLink},
                                                       resultElements);
                      resultElements = Array.concatAll(resultElements,
                                                       new String[] {ext.formDeci(pctDenovoAllele[0],
                                                                                  3),
                                                                     ext.formDeci(pctDenovoAllele[1],
                                                                                  3),
                                                                     ext.formDeci(pctDenovoAllele[2],
                                                                                  3)});
                      resultElements = Array.concatAll(resultElements, annotation);

                      if (!result.containsKey(line[3])) {
                        result.put(line[3], new Hashtable<String, String[]>());
                      }
                      result.get(line[3]).put(line[0] + "_" + allele, resultElements);

                    }
                  }
                } else {
                  // System.out.println("missing "+markerName);
                  // if (alleles[j].startsWith("+")) {
                  // System.out.println("missing "+markerName);
                  // }
                  annotation = Array.stringArray(seatleSeekHeader.length);
                  // this is capturing a subset of the indels
                  // if (readsCounts[0][0][orderedIndices[0][1]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO)
                  // {
                  // annotationNeedsVars.add(line[1] + "\t" + line[2] + "\t" + allelesOfInDels[j] +
                  // "\t" + allelesOfInDels[j]); //chr + pos + alt + alt
                  // } else {
                  if (allele.startsWith("+")) {
                    annotationNeedsInDels.add(line[1] + "\t" + line[2] + "\t" + line[2] + "\t"
                                              + allele.toUpperCase());
                  } else if (allele.startsWith("-")) {
                    annotationNeedsInDels.add(line[1] + "\t" + (Integer.parseInt(line[2]) - 1)
                                              + "\t"
                                              + (Integer.parseInt(line[2]) + allele.length() - 2)
                                              + "\t" + allele.toUpperCase());
                  } else {
                    annotationNeedsVars.add(line[1] + "\t" + line[2] + "\t" + allele.toUpperCase()
                                            + "\t" + allele.toUpperCase()); // chr + pos + alt + alt
                  }
                  // }
                }

              }
            }
          }
        }
        reader.close();

      } catch (Exception e) {
        log.report("Error reading from " + resultDir + resultFilename);
        e.printStackTrace();
      }
    }

    if (annotationNeedsVars.size() > 0) {
      Files.writeArray(Array.toStringArray(annotationNeedsVars),
                      annotationDir + "seattleSeq_input_" + new SimpleDateFormat("yyyyMMdd_hhmmss").format(new Date())
                                                                + ".txt");
    }
    if (annotationNeedsInDels.size() > 0) {
      Files.writeArray(Array.toStringArray(annotationNeedsInDels),
                      annotationDir + "seattleSeq_input_InDels_" + new SimpleDateFormat("yyyyMMdd_hhmmss").format(new Date())
                                                                  + ".txt");
    }
    if (miniSamNeeded.size() > 0) {
      // Files.writeList(Array.toStringArray(miniSamNeeded), miniSamDir + "miniBamNeeds_" + new
      // SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".txt");
      geneSet = miniSamNeeded.keySet();
      for (String gene : geneSet) {
        if (!new File(miniSamScriptsDir + gene + "/").exists()) {
          new File(miniSamScriptsDir + gene + "/").mkdir();
        }
        // Samtools.saveScriptsGeneratingMiniSams(miniSamScriptsDir + gene + "/generateMiniSams" +
        // new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".sh", miniSamDir,
        // miniSamNeeded.get(gene), loadNamesFromList(fullPathToTrioNameList), 200, log);
        Samtools.saveScriptsGeneratingMiniSamsOfOneGene(miniSamScriptsDir + gene
                                                        + "/generateMiniSams_" + gene + ".sh",
                                                        "/home/spectorl/xuz2/mini_bams/",
                                                        miniSamNeeded.get(gene),
                                                        loadNamesFromList(fullPathToTrioNameList),
                                                        200, log);
        Samtools.saveScriptsLaunchingIgv(miniSamScriptsDir + gene + "/", miniSamDir,
                                         miniSamNeeded.get(gene),
                                         loadNamesFromList(fullPathToTrioNameList), 50000, log);
      }
      Samtools.saveScriptsGeneratingMiniSamsAllGenesAtOnce(miniSamScriptsDir + "generateMiniSams_"
                                                           + new SimpleDateFormat("yyyyMMdd_hhmmss").format(new Date())
                                                           + ".sh",
                                                           "/home/spectorl/xuz2/mini_bams/",
                                                           miniSamNeeded,
                                                           loadNamesFromList(fullPathToTrioNameList),
                                                           200, log);
    }

    saveParsedResults(resultDir + DEFAULT_OUTPUT_FILENAME, null, result, geneCounts, geneTrioCounts,
                      outFormat, PARSE_RESULT_HEADER + "\tBAD\t" + Array.toStr(seatleSeekHeader),
                      log);
    log.report(ext.getTime() + "\tFinished in " + ext.getTimeElapsed(time));
  }

  public static void parseResults2(String resultDir, String seatleSeqDir,
                                   String fullpathToAdditionalComments, String miniSamScriptsDir,
                                   String miniSamDir, String fullPathToTrioNameList,
                                   byte resultFormat, double callScoreThreshold, byte outFormat,
                                   Logger log) {
    int index;
    int[] indices;
    int[][] orderedIndices;
    int[][][] readsCounts, phredScores, mappingScores;
    double[] pctDenovoAllele;
    String markerName, miniBamLink = null, miniSamScriptSubDir, trioId, tmp;
    String[] alleles, resultFilenames, annotation = null, resultElements, seatleSeekHeader, header,
        genes, line, temp;
    BufferedReader reader;
    Vector<String> annotationNeedsVars, annotationNeedsInDels;
    Hashtable<String, String[]> additionalComment;
    CountHash geneCounts;
    CountHash geneTrioCounts;
    Set<String> geneSet;
    HashSet<String> miniSamHash;
    Hashtable<String, Vector<String>> miniSamNeeded;
    Hashtable<String, Hashtable<String, String[]>> result;
    Hashtable<String, String[]> annotationHash;
    Hashtable<String, Hashtable<String, String[]>> additionalComments = null;
    HashSet<String> geneTrioCount;
    long time;

    time = new Date().getTime();
    if (!Files.exists(seatleSeqDir)) {
      log.report("Creating new annotation directory in: " + seatleSeqDir);
      new File(seatleSeqDir).mkdirs();
    }
    annotationNeedsVars = new Vector<String>();
    annotationNeedsInDels = new Vector<String>();
    // miniBamNeeds = new Vector<String>();
    miniSamHash = Samtools.listFiles(miniSamScriptsDir, log);
    annotationHash = SeattleSeq.loadAllAnnotationInDir(seatleSeqDir, log);
    additionalComments = loadAdditionalComments(fullpathToAdditionalComments);
    miniSamNeeded = new Hashtable<String, Vector<String>>();
    resultFilenames = Files.list(resultDir, ".txt", false);
    result = new Hashtable<String, Hashtable<String, String[]>>();
    orderedIndices = new int[SAMPLE_SUFFIX.length][];
    seatleSeekHeader = Matrix.extractColumn(SeattleSeq.RELEVANTS, 0);
    geneCounts = new CountHash();
    geneTrioCounts = new CountHash();
    geneTrioCount = new HashSet<String>();
    for (String resultFilename : resultFilenames) {
      trioId = resultFilename.split("_")[0];
      log.report(ext.getTime() + "\tProcessing " + trioId);

      try {
        reader = new BufferedReader(new FileReader(resultDir + resultFilename));
        // header = reader.readLine().trim().split("[\\s]+");
        header = reader.readLine().trim().split("\t");
        indices = ext.indexFactors(COLUMNS_NEEDED_FROM_PHASE1_OUTPUT, header, false, true);
        while (reader.ready()) {
          line = reader.readLine().split("\t");
          if (Double.parseDouble(line[6]) >= callScoreThreshold) {
            readsCounts = new int[SAMPLE_SUFFIX.length][1][BASES_WITH_N_INDEL.length];
            phredScores = new int[SAMPLE_SUFFIX.length][1][BASES_WITH_N_INDEL.length];
            mappingScores = new int[SAMPLE_SUFFIX.length][1][BASES_WITH_N_INDEL.length];
            parseLine(line, resultFormat, readsCounts, phredScores, mappingScores);
            getNumVarsInDelsDenovos(line[44]);

            for (int j = 0; j < SAMPLE_SUFFIX.length; j++) {
              orderedIndices[j] = Sort.quicksort(new int[] {readsCounts[j][0][0],
                                                            readsCounts[j][0][1],
                                                            readsCounts[j][0][2],
                                                            readsCounts[j][0][3]},
                                                 Sort.DESCENDING);
            }

            if (true
                // && filterMappingQuality(mappingScores, readsCounts, orderedIndices, 0)
                // && filter3Allelic(readsCounts, orderedIndices, 0)
                // && filterPhred(phredScores, readsCounts, orderedIndices, 0)
                // && filterNearbyVarInDelDenovo(numVarInDelDenovos)
                && filterChildAltAlleleCount(readsCounts, orderedIndices, 0, 5)
            // && filterInDel(readsCounts, orderedIndices, 0)
            ) {

              // line[5] = formatAltAllele(new int[] {0, 1, 2, 3, 4, orderedIndices[0][0],
              // orderedIndices[0][1]}, readsCounts, 0, line[4], (line[5].startsWith("+")? line[5] :
              // null), log);
              if (line[5] == null || line[5].equals("")) {
                alleles = new String[] {line[4]};
              } else {
                alleles = line[5].split(",");
              }

              for (String allele : alleles) {
                if (allele.startsWith("-")) {
                  markerName = "chr" + line[1] + ":" + (Integer.parseInt(line[2]) - 1) + "_"
                               + line[4].toUpperCase() + "_" + allele.toUpperCase();
                } else {
                  markerName = "chr" + line[1] + ":" + line[2] + "_" + line[4].toUpperCase() + "_"
                               + allele.toUpperCase();
                }

                if (annotationHash != null && annotationHash.containsKey(markerName)) {
                  annotation = annotationHash.get(markerName);
                  // if (annotation[0].equals("1")) {
                  if (true) {
                    if (allele.startsWith("+")) {
                      index = 5; // TODO if there are two strings for this insertion, say "+T" and
                                 // "+TGC", then cannot tell how many reads are for the former and
                                 // how many the latter. Only have the number of the reads for the
                                 // combination of the two.
                    } else if (allele.startsWith("-")) {
                      index = 6;
                    } else {
                      index = ext.indexOfChar(allele.toUpperCase().charAt(0), BASES); // TODO maybe
                                                                                      // add a
                                                                                      // column in
                                                                                      // phase1's
                                                                                      // output:
                                                                                      // DeNovo
                                                                                      // Allele
                      if (!(readsCounts[0][0][index] > MAX_ALLELE_COUNT_TREATED_AS_ZERO
                            && readsCounts[1][0][index] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                            && readsCounts[2][0][index] <= MAX_ALLELE_COUNT_TREATED_AS_ZERO
                            && (readsCounts[0][0][index] > (MAX_ALLELE_COUNT_TREATED_AS_ZERO
                                                            + readsCounts[1][0][index] + readsCounts[2][0][index])))) {
                        index = ext.indexOfChar(line[4].toUpperCase().charAt(0), BASES);
                      }
                    }

                    pctDenovoAllele = new double[SAMPLE_SUFFIX.length];
                    for (int k = 0; k < pctDenovoAllele.length; k++) {
                      pctDenovoAllele[k] = readsCounts[k][0][index]
                                           / (double) (readsCounts[k][0][0] + readsCounts[k][0][1]
                                                       + readsCounts[k][0][2] + readsCounts[k][0][3]
                                                       + readsCounts[k][0][5]
                                                       + readsCounts[k][0][6]);
                    }

                    if (pctDenovoAllele[0] >= .2) {
                      genes = annotation[7].split(",");
                      for (int k = 0; k < genes.length; k++) {
                        if (!genes[k].equals("none")) {
                          geneCounts.add(genes[k]);
                          if (!geneTrioCount.contains(genes[k] + "\t" + line[0])) { // genes[k] \t
                                                                                    // trioId
                            geneTrioCounts.add(genes[k]);
                            geneTrioCount.add(genes[k] + "\t" + line[0]); // genes[k] \t trioId
                          }
                        }

                        // TODO "none"?
                        miniSamScriptSubDir = miniSamScriptsDir + genes[k] + "/";
                        if (!new File(miniSamScriptSubDir).exists()) {
                          new File(miniSamScriptSubDir).mkdir();
                        }

                        // if (new File(miniSamSubDir).exists()) {
                        // miniSamHash = Samtools.listFiles(miniSamSubDir, log);
                        // TODO is the new region contained in the older region?
                        if (miniSamHash.contains(line[0] + "\t" + line[1] + "\t" + line[2])) { // trioId
                                                                                               // +
                                                                                               // "\t"
                                                                                               // +
                                                                                               // chr
                                                                                               // +
                                                                                               // "\t"
                                                                                               // +
                                                                                               // pos
                          // if (! new File(miniSamSubDir + line[0] + "_" + line[1] + "_" + line[2]
                          // + ".xml").exists()) {
                          // Files.write(Samtools.getIgvXmlScript(miniSamDir, line[1], line[2],
                          // miniSamFilenamesOfOneTrio), miniSamSubDir + line[0] + "_" + line[1] +
                          // "_" + line[2] + ".xml");
                          // }
                          if (!new File(miniSamScriptSubDir + line[0] + "_chr" + line[1] + "_"
                                        + line[2] + ".bat").exists()) {
                            Files.write(Samtools.getIgvLaunchScript(miniSamScriptSubDir + line[0]
                                                                    + "_" + line[1] + "_" + line[2]
                                                                    + ".xml"),
                                        miniSamScriptSubDir + line[0] + "_" + line[1] + "_"
                                                                               + line[2] + ".bat");
                          }
                          if (!new File(miniSamScriptSubDir + line[0] + "_chr" + line[1] + "_"
                                        + line[2] + ".xml").exists()) {
                            Files.write(Samtools.getIgvXmlScript(miniSamDir, line[1], line[2],
                                                                 new String[] {line[0],
                                                                               line[0] + "_chr"
                                                                                        + line[1]
                                                                                        + "_"
                                                                                        + line[2]
                                                                                        + "_C.bam",
                                                                               line[0] + "_chr" + line[1]
                                                                                                    + "_"
                                                                                                    + line[2]
                                                                                                    + "_D.bam",
                                                                               line[0] + "_chr" + line[1] + "_" + line[2]
                                                                                                                + "_M.bam"}),
                                        miniSamScriptSubDir + line[0] + "_chr" + line[1] + "_" + line[2] + ".xml");
                          }
                          miniBamLink = "=HYPERLINK(\"" + miniSamScriptSubDir + line[0] + "_chr"
                                        + line[1] + "_" + line[2] + ".bat\",\"IGV\")";
                          // miniBamLink = miniSamSubDir;
                        } else {
                          miniBamLink = "";
                          if (!miniSamNeeded.containsKey(genes[k])) {
                            miniSamNeeded.put(genes[k], new Vector<String>());
                          }
                          miniSamNeeded.get(genes[k])
                                       .add(line[0] + "\t" + line[1] + "\t" + line[2]); // trioId +
                                                                                        // chr + pos
                        }
                        // } else {
                        // miniBamLink = "";
                        // if (! miniSamNeeded.containsKey(genes[k])) {
                        // miniSamNeeded.put(genes[k], new Vector<String>());
                        // }
                        // miniSamNeeded.get(genes[k]).add(line[0] + "\t" + line[1] + "\t" +
                        // line[2]); //trioId + chr + pos
                        // }
                      }

                      additionalComment = additionalComments.get(line[1] + ":" + line[2]);
                      // additionalComment = additionalComments.get(line[3].replaceAll("chr", ""));

                      resultElements = Array.subArray(line, indices);
                      resultElements = Array.concatAll(new String[] {annotation[7], miniBamLink},
                                                       resultElements);
                      resultElements = Array.concatAll(resultElements,
                                                       new String[] {ext.formDeci(pctDenovoAllele[0],
                                                                                  3),
                                                                     ext.formDeci(pctDenovoAllele[1],
                                                                                  3),
                                                                     ext.formDeci(pctDenovoAllele[2],
                                                                                  3)});
                      resultElements = Array.concatAll(resultElements, annotation);
                      if (additionalComment != null
                          && additionalComment.containsKey(line[45].substring(0, 2))) {
                        resultElements =
                                       Array.concatAll(resultElements,
                                                       additionalComment.get(line[45].substring(0,
                                                                                                2)));
                      } else {
                        resultElements =
                                       Array.concatAll(resultElements,
                                                       new String[] {".", ".", ".", ".", ".", "."});
                      }
                      tmp = "";
                      if (additionalComment != null) {
                        for (String allelePair : additionalComment.keySet()) {
                          if (!allelePair.equals(line[45].substring(0, 2))) {
                            if (!tmp.equals("")) {
                              tmp += " /";
                            }
                            tmp += allelePair;
                            temp = additionalComment.get(allelePair);
                            for (String element : temp) {
                              tmp += ("; " + element);
                            }
                          }
                        }
                      }
                      resultElements = Array.concatAll(resultElements, new String[] {tmp});

                      if (!result.containsKey(line[3])) {
                        result.put(line[3], new Hashtable<String, String[]>());
                      }
                      result.get(line[3]).put(line[0] + "_" + allele, resultElements);

                    }
                  }
                } else {
                  // System.out.println("missing "+markerName);
                  // if (alleles[j].startsWith("+")) {
                  // System.out.println("missing "+markerName);
                  // }
                  annotation = Array.stringArray(seatleSeekHeader.length);
                  // this is capturing a subset of the indels
                  // if (readsCounts[0][0][orderedIndices[0][1]] > MAX_ALLELE_COUNT_TREATED_AS_ZERO)
                  // {
                  // annotationNeedsVars.add(line[1] + "\t" + line[2] + "\t" + allelesOfInDels[j] +
                  // "\t" + allelesOfInDels[j]); //chr + pos + alt + alt
                  // } else {
                  if (allele.startsWith("+")) {
                    annotationNeedsInDels.add(line[1] + "\t" + line[2] + "\t" + line[2] + "\t"
                                              + allele.toUpperCase());
                  } else if (allele.startsWith("-")) {
                    annotationNeedsInDels.add(line[1] + "\t" + (Integer.parseInt(line[2]) - 1)
                                              + "\t"
                                              + (Integer.parseInt(line[2]) + allele.length() - 2)
                                              + "\t" + allele.toUpperCase());
                  } else {
                    annotationNeedsVars.add(line[1] + "\t" + line[2] + "\t" + allele.toUpperCase()
                                            + "\t" + allele.toUpperCase()); // chr + pos + alt + alt
                  }
                  // }
                }

              }
            }
          }
        }
        reader.close();

      } catch (Exception e) {
        log.report("Error reading from " + resultDir + resultFilename);
        e.printStackTrace();
      }
    }

    if (annotationNeedsVars.size() > 0) {
      Files.writeArray(Array.toStringArray(annotationNeedsVars),
                      seatleSeqDir + "seattleSeq_input_" + new SimpleDateFormat("yyyyMMdd_hhmmss").format(new Date())
                                                                + ".txt");
    }
    if (annotationNeedsInDels.size() > 0) {
      Files.writeArray(Array.toStringArray(annotationNeedsInDels),
                      seatleSeqDir + "seattleSeq_input_InDels_" + new SimpleDateFormat("yyyyMMdd_hhmmss").format(new Date())
                                                                  + ".txt");
    }
    if (miniSamNeeded.size() > 0) {
      // Files.writeList(Array.toStringArray(miniSamNeeded), miniSamDir + "miniBamNeeds_" + new
      // SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".txt");
      geneSet = miniSamNeeded.keySet();
      for (String gene : geneSet) {
        if (!new File(miniSamScriptsDir + gene + "/").exists()) {
          new File(miniSamScriptsDir + gene + "/").mkdir();
        }
        // Samtools.saveScriptsGeneratingMiniSams(miniSamScriptsDir + gene + "/generateMiniSams" +
        // new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".sh", miniSamDir,
        // miniSamNeeded.get(gene), loadNamesFromList(fullPathToTrioNameList), 200, log);
        Samtools.saveScriptsGeneratingMiniSamsOfOneGene(miniSamScriptsDir + gene
                                                        + "/generateMiniSams_" + gene + ".sh",
                                                        "/home/spectorl/xuz2/mini_bams/",
                                                        miniSamNeeded.get(gene),
                                                        loadNamesFromList(fullPathToTrioNameList),
                                                        200, log);
        Samtools.saveScriptsLaunchingIgv(miniSamScriptsDir + gene + "/", miniSamDir,
                                         miniSamNeeded.get(gene),
                                         loadNamesFromList(fullPathToTrioNameList), 50000, log);
      }
      Samtools.saveScriptsGeneratingMiniSamsAllGenesAtOnce(miniSamScriptsDir + "generateMiniSams_"
                                                           + new SimpleDateFormat("yyyyMMdd_hhmmss").format(new Date())
                                                           + ".sh",
                                                           "/home/spectorl/xuz2/mini_bams/",
                                                           miniSamNeeded,
                                                           loadNamesFromList(fullPathToTrioNameList),
                                                           200, log);
    }

    saveParsedResults(resultDir + DEFAULT_OUTPUT_FILENAME, null, result, geneCounts, geneTrioCounts,
                      outFormat,
                      PARSE_RESULT_HEADER + "\tBAD\t" + Array.toStr(seatleSeekHeader)
                                 + "\tSKATgene\tfunc_region\tMAF_whites\tn_whites\tMAF_blacks\tn_blacks\totherAllelePairs",
                      log);
    log.report(ext.getTime() + "\tFinished in " + ext.getTimeElapsed(time));
  }

  public static Hashtable<String, Hashtable<String, String[]>> loadAdditionalComments(String filename) {
    BufferedReader reader;
    String[] line;
    Hashtable<String, Hashtable<String, String[]>> result;
    Hashtable<String, String[]> allelePairs;
    int i;

    result = new Hashtable<String, Hashtable<String, String[]>>();
    try {
      reader = Files.getAppropriateReader(filename);
      reader.readLine();
      while (reader.ready()) {
        line = reader.readLine().split("\t");
        if (result.containsKey(line[0])) {
          allelePairs = result.get(line[0]);
        } else {
          allelePairs = new Hashtable<String, String[]>();
          result.put(line[0], allelePairs);
        }
        i = 3;
        while (line.length > i) {
          if (line.length > 9) {
            System.out.println();
          }
          // allelePairs.put(line[1] + "\t" + line[2], new String[] {line[i], line[1 + i], line[2 +
          // i], line[3 + i], line[4 + i], line[5 + i]});
          allelePairs.put(line[1] + line[2], new String[] {line[i], line[1 + i], line[2 + i],
                                                           line[3 + i], line[4 + i], line[5 + i]});
          i += 6;
        }
      }
      reader.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }

    return result;
  }

  public static void saveParsedResults(String resultsFilename, String[] trioIds,
                                       Hashtable<String, Hashtable<String, String[]>> result,
                                       CountHash geneCounts, CountHash geneTrioCounts, byte format,
                                       String header, Logger log) {
    PrintWriter writer;
    String[] line;
    Object[] keys;
    Hashtable<String, String[]> temp;
    Set<Entry<String, String[]>> keySet;
    String key;
    int geneCount, geneTrioCount;
    int[] numVarInDelDenovoNearby;

    try {
      keys = result.keySet().toArray();
      Arrays.sort(keys);
      // int[] order = Sort.orderTwoLayers(chrs, positions);
      // keys = Sort.putInOrder(keys, order);
      writer = new PrintWriter(new FileWriter(resultsFilename));
      if (format == 1) {
        // writer.println("id\tchr\tposition\tlookup\tnumHits\tref\talt\tcall\tnote\tfwdGenotype\treadDepth_C\treadDepth_D\treadDepth_M\tPhredScores\tMappingQuality_C\tMappingQuality_D\tMappingQuality_M\taltAllele%_C\taltAllele%_D\taltAllele%_M");
        writer.println(header);
        for (Object key2 : keys) {
          keySet = result.get(key2).entrySet();
          for (Entry<String, String[]> entry : keySet) {
            key = (String) key2;
            line = entry.getValue();
            if (geneCounts == null) {
              geneCount = 0;
            } else {
              geneCount = geneCounts.getCount(line[35].split(",")[0]);
            }
            if (geneTrioCounts == null) {
              geneTrioCount = 0;
            } else {
              geneTrioCount = geneTrioCounts.getCount(line[35].split(",")[0]);
            }

            // writer.print(key.substring(3, key.indexOf(":")) + "\t" +
            // key.substring(key.indexOf(":") + 1) + "\t" + key.substring(0, key.indexOf("_")) +
            // "\t" + (geneCount == 0? "." : geneCount) + "\t" + (geneTrioCount == 0? "." :
            // geneTrioCount) + "\t" + ((geneCount - geneTrioCount) == 0? "." : (geneCount -
            // geneTrioCount)) + "\t" + entry.getKey());
            writer.print(key.substring(3, key.indexOf(":")) + "\t"
                         + key.substring(key.indexOf(":") + 1) + "\t" + key + "\t"
                         + (geneCount == 0 ? "." : geneCount) + "\t"
                         + (geneTrioCount == 0 ? "." : geneTrioCount) + "\t"
                         + ((geneCount - geneTrioCount) == 0 ? "." : (geneCount - geneTrioCount))
                         + "\t" + entry.getKey().substring(0, entry.getKey().lastIndexOf("_")));
            numVarInDelDenovoNearby = getNumVarsInDelsDenovos(line[5]);
            for (int element : numVarInDelDenovoNearby) {
              writer.print("\t" + (element == 0 ? "" : element));
            }
            for (String element : line) {
              writer.print("\t" + element);
            }
            writer.println();
          }
        }
      } else if (format == 2) {
        // writer.print("chr\tposition\tlookup\tnumTrios");
        writer.print(header);
        for (String trioId : trioIds) {
          writer.print("\t" + trioId);
        }
        writer.println();
        for (Object key2 : keys) {
          temp = result.get(key2);
          line = ((String) key2).split(":");
          writer.print(line[0] + "\t" + line[1] + "\t" + key2 + "\t" + temp.size());
          for (String trioId : trioIds) {
            if (temp.containsKey(trioId)) {
              line = temp.get(trioId);
              writer.print("\t" + line[0] + "; " + line[1]);
            } else {
              writer.print("\t");
            }
          }
          writer.println();
        }
      }

      writer.close();
      log.report(ext.getTime() + "\tSummarized result is available at: " + resultsFilename);

    } catch (Exception e) {
      log.reportError("Error writing to " + resultsFilename);
      e.printStackTrace();
    }
  }

  public static int[] getNumVarsInDelsDenovos(String line) {
    String[] tmp1;
    String[] tmp2;
    String[] tmp3;
    int index, tmp;
    HashSet[] counter;

    counter = new HashSet[4];
    for (int i = 0; i < counter.length; i++) {
      counter[i] = new HashSet();
    }
    tmp1 = line.split("\\)");
    for (int i = 0; i < tmp1.length; i++) {
      tmp1[i] = tmp1[i].trim();
      if (tmp1[i].startsWith("Var")) {
        index = 0;
      } else if (tmp1[i].startsWith("Del")) {
        index = 1;
      } else if (tmp1[i].startsWith("Ins")) {
        index = 2;
      } else if (tmp1[i].startsWith("DeNovo")) {
        index = 3;
      } else {
        index = -1;
      }

      if (index != -1) {
        tmp2 = tmp1[i].substring(tmp1[i].indexOf("(") + 1).split(" ");
        for (String element : tmp2) {
          tmp3 = element.substring(element.indexOf(":") + 1).split(",");
          for (String element2 : tmp3) {
            tmp = Integer.parseInt(element2);
            if (tmp <= 50 && tmp >= -50 && tmp != 0) {
              counter[index].add(element2);
            }
          }
        }
      }
    }

    // return counter[0] + "\t" + counter[1] + "\t" + counter[2] + "\t" + counter[3];
    return new int[] {counter[0].size(), counter[1].size(), counter[2].size(), counter[3].size()};
  }

  public static void parseLine(String[] line, byte resultFormat, int[][][] readsCounts,
                               int[][][] phredScores, int[][][] mappingScores) {
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
          if (line[index] != null && !line[index].isEmpty()) {
            readsCounts[i][0][j] = Integer.parseInt(line[index]);
          }
          index++;
        }
      }

      index = 62;
      for (int i = 0; i < SAMPLE_SUFFIX.length; i++) {
        for (int j = 0; j < BASES_WITH_N_INDEL.length; j++) {
          phredScores[i][0][j] = (line[index] == null
                                  || line[index].equals("") ? 0
                                                            : (Integer.parseInt(line[index])
                                                               * readsCounts[i][0][j]));
          index++;
        }
      }

      index = 83;
      for (int i = 0; i < SAMPLE_SUFFIX.length; i++) {
        if (line[index] != null && !line[index].isEmpty() && !line[index].equals(".")) {
          score = Integer.parseInt(line[index]);
          for (int j = 0; j < mappingScores[i][0].length; j++) {
            mappingScores[i][0][j] = score * readsCounts[i][0][j];
          }
        }
        index++;
      }
    }
  }

  public static void saveAlleleCountsToFile(PrintWriter writer, String chr, int startPosition,
                                            int[][][] readsCounts, int[][][] phredScores,
                                            int[][][] mappingScores, byte format,
                                            boolean toIgnorePositionsWithZeroCount) {
    String line;
    int[] totalReads;
    int[] totalMappingScores;
    int[] totalPhred;
    int total;

    totalReads = new int[readsCounts.length];
    totalMappingScores = new int[readsCounts.length];
    totalPhred = new int[readsCounts.length];
    try {
      for (int i = 0; i < readsCounts[0].length; i++) {
        line = chr + "\t" + startPosition;

        if (format == 1) {
          for (int[][] readsCount : readsCounts) {
            line += "\t(";
            for (int k = 0; k < BASES.length; k++) {
              line += ((k == 0 ? "" : ",") + (readsCount[i][k] == 0 ? "-" : readsCount[i][k]));
            }
            line += "/";
            for (int k = BASES.length; k < BASES_WITH_N_INDEL.length; k++) {
              line += ((k == 0 ? "" : ",") + (readsCount[i][k] == 0 ? "-" : readsCount[i][k]));
            }
            line += ")";
          }
        }

        for (int j = 0; j < readsCounts.length; j++) {
          totalReads[j] = 0;
          totalMappingScores[j] = 0;
          totalPhred[j] = 0;
          for (int k = 0; k < BASES.length; k++) {
            totalReads[j] += readsCounts[j][i][k];
            totalMappingScores[j] += mappingScores[j][i][k];
            totalPhred[j] += phredScores[j][i][k];
          }
          for (int k = BASES.length; k < BASES_WITH_N_INDEL.length; k++) {
            totalReads[j] += readsCounts[j][i][k];
            totalMappingScores[j] += mappingScores[j][i][k];
            totalPhred[j] += phredScores[j][i][k];
          }
        }

        total = 0;
        for (int totalRead : totalReads) {
          total += totalRead;
        }

        if (!(toIgnorePositionsWithZeroCount && total == 0)) {
          for (int totalRead : totalReads) {
            line += "\t" + (totalRead == 0 ? "" : totalRead);
          }

          for (int j = 0; j < totalReads.length; j++) {
            if (totalReads[j] != 0) {
              line += "\t" + (totalMappingScores[j] / totalReads[j]);
            } else {
              line += "\t";
            }
          }

          for (int j = 0; j < totalReads.length; j++) {
            if (totalReads[j] != 0) {
              line += "\t" + (totalPhred[j] / totalReads[j]);
            } else {
              line += "\t";
            }
          }

          writer.println(line);
        }
        startPosition++;
      }
    } catch (Exception e) {
      System.out.println(e);
    }
  }

  /**
   * Get the common root of all the filenames, assuming the filenames following the format of
   * "*_*_*"
   *
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
        if (!roots[j][i].equalsIgnoreCase(commonRoot)) {
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
          if (!roots[j][index2].substring(0, i).equalsIgnoreCase(commonRoot.substring(0, i))) {
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
    Vector<String[]> result;
    String[][] result2;
    result = new Vector<String[]>();
    try {
      reader = new BufferedReader(new FileReader(fullPathToTheList));
      while (reader.ready()) {
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
    result = new Vector<String[]>(filenames.length / 3);
    for (int i = 0; i < filenames.length; i++) {
      if (!isVisited[i] && filenames[i].contains("C_L")) {
        location = filenames[i].indexOf("C_L");
        currentTrio = filenames[i].substring(0, location);
        indices[0] = i;
        for (int j = 0; j < filenames.length; j++) {
          if (!isVisited[j] && i != j && filenames[j].substring(0, location).equals(currentTrio)) {
            isVisited[j] = true;

            if (filenames[j].contains("D_L")) {
              indices[1] = j;
            } else if (filenames[j].contains("M_L")) {
              indices[2] = j;
            } else {
              System.out.println("These file names do not make sense: " + filenames[j] + "\n"
                                 + filenames[indices[0]]);
            }

            if (indices[0] >= 0 && indices[1] >= 0 && indices[2] >= 0) {
              result.add(new String[] {filenames[indices[0]], filenames[indices[1]],
                                       filenames[indices[2]]});
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
    int numMatches;

    isVisited = new boolean[filenames.length];
    indices = new int[3];
    for (int i = 0; i < indices.length; i++) {
      indices[i] = -1;
    }
    result = new Vector<String[]>(filenames.length / 3);
    nameSegments = filenames[0].split("_");
    numMatches = 0;
    for (int i = 0; i < nameSegments.length; i++) {
      for (String filename : filenames) {
        if (filename.split("_")[i].equalsIgnoreCase(nameSegments[i])) {
          numMatches++;
        }
      }
      if (numMatches == filenames.length) {

      }
    }
    for (int i = 0; i < filenames.length; i++) {
      if (!isVisited[i] && filenames[i].contains("C_L")) {
        location = filenames[i].indexOf("C_L");
        currentTrio = filenames[i].substring(0, location);
        indices[0] = i;
        for (int j = 0; j < filenames.length; j++) {
          if (!isVisited[j] && i != j && filenames[j].substring(0, location).equals(currentTrio)) {
            isVisited[j] = true;

            if (filenames[j].contains("D_L")) {
              indices[1] = j;
            } else if (filenames[j].contains("M_L")) {
              indices[2] = j;
            } else {
              System.out.println("These file names do not make sense: " + filenames[j] + "\n"
                                 + filenames[indices[0]]);
            }

            if (indices[0] >= 0 && indices[1] >= 0 && indices[2] >= 0) {
              result.add(new String[] {filenames[indices[0]], filenames[indices[1]],
                                       filenames[indices[2]]});
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

  public static void getHaplotypes(String fullPathToParsedResult, String fullPathToTrioNameList,
                                   String bamDir, boolean toOutputHaplotypeStrings, Logger log) {
    String dirAndRoot = null;
    PrintWriter writer;
    BufferedReader reader;
    String[][] bamFilenamesByTrios;
    String[] line;
    int x;

    bamFilenamesByTrios = loadNamesFromList(fullPathToTrioNameList);
    if (toOutputHaplotypeStrings) {
      dirAndRoot = ext.parseDirectoryOfFile(fullPathToParsedResult)
                   + ext.rootOf(fullPathToParsedResult);
    }
    try {
      writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(fullPathToParsedResult)
                                              + ext.rootOf(fullPathToParsedResult)
                                              + "_haplotypeCount.txt"));
      reader = Files.getAppropriateReader(fullPathToParsedResult);
      reader.readLine();
      while (reader.ready()) {
        line = reader.readLine().split("\t");
        x = -1;
        for (int i = 0; i < bamFilenamesByTrios.length; i++) {
          if (line[6].equals(bamFilenamesByTrios[i][0])) {
            x = i;
            break;
          }
        }
        if (x != -1) {
          writer.println(line[2] + "\t" + line[6] + "\t"
                         + getHaplotypes(new String[] {bamDir + bamFilenamesByTrios[x][1],
                                                       bamDir + bamFilenamesByTrios[x][2],
                                                       bamDir + bamFilenamesByTrios[x][3]},
                                         line[6], line[0], Integer.parseInt(line[1]),
                                         Integer.parseInt(line[1]),
                                         (toOutputHaplotypeStrings ? (dirAndRoot + "_" + line[6])
                                                                   : null),
                                         log));
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + fullPathToParsedResult
                      + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + fullPathToParsedResult + "\"");
      return;
    }
  }

  public static int getHaplotypes(String[] fullpathsToSamFiles, String trioId, String chr,
                                  int beginPosition, int endPosition, String haplotypeFilename,
                                  Logger log) {
    Vector<String[]> samContentVec;
    int numLines, window, beginPositionExt, endPositionExt;
    char[] ref;
    String[] tmp;
    Hashtable<Integer, String> insertionStringOfCurrentRead, insertionPhredOfCurrentRead;
    Hashtable<Byte, Hashtable<Integer, String>> haplotypeDb;
    Hashtable<Byte, Vector<int[]>> readCounts;

    window = WINDOW_SIZE_FOR_HAPLOTYPE_CHECK;
    beginPositionExt = Math.max(0, beginPosition - window);
    endPositionExt = Math.min(
                              Positions.CHROMOSOME_LENGTHS_MAX[ext.indexOfStr(chr,
                                                                              Positions.CHR_CODES,
                                                                              false, true)],
                              endPosition + window);
    haplotypeDb = new Hashtable<Byte, Hashtable<Integer, String>>();
    readCounts = new Hashtable<Byte, Vector<int[]>>();
    ref = new char[endPositionExt - beginPositionExt + 1];

    for (int i = 0; i < fullpathsToSamFiles.length; i++) {
      if (!new File(fullpathsToSamFiles[i]).exists()) {
        return -1;
      }

      // samContentVec = loadBamByPipeline(samDir, samFilenames[i], chr, beginPositionExt,
      // endPositionExt);
      samContentVec = loadBamByFileReader(fullpathsToSamFiles[i], chr, beginPositionExt,
                                          endPositionExt);
      numLines = samContentVec.size();
      for (int j = 0; j < numLines; j++) {
        insertionStringOfCurrentRead = new Hashtable<Integer, String>();
        insertionPhredOfCurrentRead = new Hashtable<Integer, String>();
        tmp = cigarToReFormatedString(samContentVec.elementAt(j), beginPositionExt, endPositionExt,
                                      insertionStringOfCurrentRead, insertionPhredOfCurrentRead);
        if (tmp != null && (beginPositionExt + Integer.parseInt(tmp[0])) <= endPosition
            && (beginPositionExt + Integer.parseInt(tmp[0]) + tmp[1].length()) >= beginPosition) { // TODO
          updateRefAndAlt(tmp[1], ref, Integer.parseInt(tmp[0]), insertionStringOfCurrentRead,
                          haplotypeDb, readCounts);
          // log.report(haplotypeDb.size() + "\t" + tmp[0] + "\t" + tmp[1] + "\t" +
          // formatHaplotype(haplotypeDb, beginPositionExt));
          // log.report(haplotypeDb.size() + "\t" + tmp[0] + "\t" + tmp[1]);
          // log.report(samContentVec.elementAt(j)[0]
          // + "\t" + samContentVec.elementAt(j)[1]
          // + "\t" + samContentVec.elementAt(j)[2]
          // + "\t" + samContentVec.elementAt(j)[3]
          // + "\t" + samContentVec.elementAt(j)[4]
          // + "\t" + samContentVec.elementAt(j)[5]
          // + "\t" + samContentVec.elementAt(j)[6]
          // + "\t" + samContentVec.elementAt(j)[7]
          // + "\t" + samContentVec.elementAt(j)[8]
          // + "\t" + samContentVec.elementAt(j)[9]
          // + "\t" + samContentVec.elementAt(j)[10]);
        }
      }
    }

    if (haplotypeFilename != null) {
      exportHaplotypes(haplotypeDb, readCounts, ref, beginPositionExt, haplotypeFilename);
    }
    return (haplotypeDb.size() + 1);
  }

  /*
   * Original, working version
   */
  /*
   * public static int getHaplotypes(String samDir, String[] samFilenames, String trioId, String
   * chr, int beginPosition, int endPosition, String haplotypeFilename, Logger log) { Process p; //
   * ProcessBuilder ps; BufferedReader reader; // BufferedReader error; String line;
   * Vector<String[]> samContentVec; int numLines, numMarkersPlusWindow, window, begin, end,
   * beginPositionExt, endPositionExt; // int[][][] readsCounts = null; // int[][][] phredScores =
   * null; // int[][][] mappingScores = null; char[] ref; // Vector<String> haplotypeStrings = null;
   * String[] tmp; // Hashtable<Integer, String> altAllelesForInsertions_Child; Hashtable<Integer,
   * String> insertionStringOfCurrentRead, insertionPhredOfCurrentRead; Hashtable<Byte,
   * Hashtable<Integer, String>> haplotypeDb; Hashtable<Byte, Vector<int[]>> readCounts;
   *
   * window = WINDOW_SIZE_FOR_HAPLOTYPE_CHECK; beginPositionExt = Math.max(0, beginPosition -
   * window); endPositionExt = Math.min(Positions.CHROMOSOME_LENGTHS_MAX[ext.indexOfStr(chr,
   * Positions.CHR_CODES, false, true)], endPosition + window); numMarkersPlusWindow =
   * endPositionExt - beginPositionExt + 1; // readsCounts = new
   * int[SAMPLE_SUFFIX.length][numMarkersPlusWindow][READS_COUNTS_ARRAY_STRUCT.length]; //
   * phredScores = new int[SAMPLE_SUFFIX.length][numMarkersPlusWindow][BASES_WITH_N_INDEL.length];
   * // mappingScores = new
   * int[SAMPLE_SUFFIX.length][numMarkersPlusWindow][BASES_WITH_N_INDEL.length]; //
   * altAllelesForInsertions_Child = new Hashtable<Integer, String>(); haplotypeDb = new
   * Hashtable<Byte, Hashtable<Integer, String>>(); readCounts = new Hashtable<Byte,
   * Vector<int[]>>(); ref = new char[endPositionExt - beginPositionExt + 1]; // haplotypeStrings =
   * new Vector<String>(); // haplotypeStrings.add("");
   *
   * if (! Files.exists(samDir, samFilenames)) { return -1; }
   *
   * for (int i = 0; i < samFilenames.length; i++) { // samContentVec = loadBamByPipeline(samDir,
   * samFilenames[i], chr, beginPositionExt, endPositionExt); samContentVec =
   * loadBamByFileReader(samDir, samFilenames[i], chr, beginPositionExt, endPositionExt); numLines =
   * samContentVec.size(); for (int j = 0; j < numLines; j++) { insertionStringOfCurrentRead = new
   * Hashtable<Integer, String>(); insertionPhredOfCurrentRead = new Hashtable<Integer, String>();
   * tmp = cigarToReFormatedString(samContentVec.elementAt(j), beginPositionExt, endPositionExt,
   * insertionStringOfCurrentRead, insertionPhredOfCurrentRead); if (tmp != null) {
   * updateRefAndAlt(tmp[1], ref, Integer.parseInt(tmp[0]), insertionStringOfCurrentRead,
   * haplotypeDb, readCounts); // log.report(haplotypeDb.size() + "\t" + tmp[0] + "\t" + tmp[1] +
   * "\t" + formatHaplotype(haplotypeDb, beginPositionExt)); // log.report(haplotypeDb.size() + "\t"
   * + tmp[0] + "\t" + tmp[1]); // log.report(samContentVec.elementAt(j)[0] // + "\t" +
   * samContentVec.elementAt(j)[1] // + "\t" + samContentVec.elementAt(j)[2] // + "\t" +
   * samContentVec.elementAt(j)[3] // + "\t" + samContentVec.elementAt(j)[4] // + "\t" +
   * samContentVec.elementAt(j)[5] // + "\t" + samContentVec.elementAt(j)[6] // + "\t" +
   * samContentVec.elementAt(j)[7] // + "\t" + samContentVec.elementAt(j)[8] // + "\t" +
   * samContentVec.elementAt(j)[9] // + "\t" + samContentVec.elementAt(j)[10]); } } }
   *
   * if (haplotypeFilename != null) { exportHaplotypes(haplotypeDb, readCounts, ref,
   * beginPositionExt, haplotypeFilename); } return (haplotypeDb.size() + 1); }
   */
  public static String[] cigarToReFormatedString(String[] aSingleLineOfBamFile, int beginPos,
                                                 int endPos,
                                                 Hashtable<Integer, String> insertionStrings,
                                                 Hashtable<Integer, String> insertionPhreds) {
    int readPointer, outputStringPointer, currentPosition, lengthOfCurrentSegment,
        outputStringStopLocation, loop, beginPosOfOutputString = -1;
    String[][] readSegments;
    String reFormatedSeqString = null, reFormatedPhredString = null;

    if (Integer.parseInt(aSingleLineOfBamFile[1]) < 512
        && !aSingleLineOfBamFile[9].equalsIgnoreCase("*")) {
      reFormatedSeqString = "";
      reFormatedPhredString = "";
      currentPosition = Integer.parseInt(aSingleLineOfBamFile[3]);
      readSegments = ext.getOperatorsOperatorIndicesAndSplit(aSingleLineOfBamFile[5], "MIDNSHP=X");
      readPointer = 0;
      outputStringPointer = currentPosition - beginPos;
      outputStringStopLocation = endPos - beginPos + 1;
      beginPosOfOutputString = Math.max(outputStringPointer, 0);

      for (int i = 0; (outputStringPointer < outputStringStopLocation)
                      && (i < readSegments[0].length); i++) {
        lengthOfCurrentSegment = Integer.parseInt(readSegments[2][i]);
        if (readSegments[0][i].equals("M") || readSegments[0][i].equals("=")
            || readSegments[0][i].equals("X")) { // present in read sequence, and also in reference
                                                 // sequence.
          if ((outputStringPointer + lengthOfCurrentSegment) <= 0) {
            currentPosition += lengthOfCurrentSegment;
            outputStringPointer += lengthOfCurrentSegment;
            readPointer += lengthOfCurrentSegment;
          } else {
            if (outputStringPointer < 0) {
              currentPosition -= outputStringPointer;
              readPointer -= outputStringPointer;
              lengthOfCurrentSegment += outputStringPointer;
              outputStringPointer = 0;
              // if (beginPosOfOutputString < 0) {
              // beginPosOfOutputString = 0;
              // }
            }

            loop = Math.min(lengthOfCurrentSegment, outputStringStopLocation - outputStringPointer);
            outputStringPointer += loop;
            loop += readPointer;
            reFormatedSeqString += aSingleLineOfBamFile[9].substring(readPointer, loop);
            reFormatedPhredString += aSingleLineOfBamFile[10].substring(readPointer, loop);
            readPointer = loop;
            // currentPosition += lengthOfCurrentSegment;
            currentPosition += loop;
          }

        } else if (readSegments[0][i].equals("I")) { // present in read sequence, but NOT in
                                                     // reference sequence.
          if (outputStringPointer > 0) {
            loop = readPointer + lengthOfCurrentSegment;
            // reFormatedSeqString += (INSERTION_DELIMITER_LEFT +
            // aSingleLineOfBamFile[9].substring(readPointer, loop) + INSERTIONS_DELIMITER_RIGHT);
            // reFormatedPhredString += (INSERTION_DELIMITER_LEFT +
            // aSingleLineOfBamFile[10].substring(readPointer, loop) + INSERTIONS_DELIMITER_RIGHT);
            insertionStrings.put(outputStringPointer - 1,
                                 "+" + aSingleLineOfBamFile[9].substring(readPointer, loop));
            insertionPhreds.put(outputStringPointer - 1,
                                "+" + aSingleLineOfBamFile[10].substring(readPointer, loop));

            // output1_ReadsCounts[outputStringPointer][INDEX_OF_INS] ++;
            // output2_PhredScores[outputStringPointer][INDEX_OF_INS] +=
            // convertToPhredScore(aSingleLineOfBamFile[10].charAt(readPointer));
            // output3_MappingScores[outputStringPointer][INDEX_OF_INS] += currentMappingScore;
          }
          readPointer += lengthOfCurrentSegment;

        } else if (readSegments[0][i].equals("D")) { // NOT present in read sequence, but does in
                                                     // reference sequence.
          if ((outputStringPointer + lengthOfCurrentSegment) <= 0) {
            currentPosition += lengthOfCurrentSegment;
            outputStringPointer += lengthOfCurrentSegment;
          } else {
            if (outputStringPointer < 0) {
              currentPosition -= outputStringPointer;
              lengthOfCurrentSegment += outputStringPointer;
              outputStringPointer = 0;
            }

            loop = outputStringPointer + Math.min(lengthOfCurrentSegment,
                                                  outputStringStopLocation - outputStringPointer);
            while (outputStringPointer < loop) {
              reFormatedSeqString += "-";
              reFormatedPhredString += "-";
              // output1_ReadsCounts[outputStringPointer][INDEX_OF_DEL] ++;
              // output2_PhredScores[outputStringPointer][INDEX_OF_DEL] +=
              // DEFAULT_PHRED_SCORE_FOR_DELETION;
              // output3_MappingScores[outputStringPointer][INDEX_OF_DEL] += currentMappingScore;

              outputStringPointer++;
            }

            currentPosition += lengthOfCurrentSegment;
          }

        } else if (readSegments[0][i].equals("N")) { // NOT present in read sequence, but does in
                                                     // reference sequence. Similar to D.
          currentPosition += lengthOfCurrentSegment;
          outputStringPointer += lengthOfCurrentSegment;

        } else if (readSegments[0][i].equals("S")) { // present in read sequence, and also in
                                                     // reference sequence, but do not match. Lower
                                                     // case letters in read sequence.
          if (i == 0) {
            readPointer += lengthOfCurrentSegment;
          } else {
            currentPosition += lengthOfCurrentSegment;
            outputStringPointer += lengthOfCurrentSegment;
            readPointer += lengthOfCurrentSegment;
          }

        } else if (readSegments[0][i].equals("H")) { // NOT present in read sequence, but in
                                                     // reference sequence. Similar to D.
          currentPosition += lengthOfCurrentSegment;
          outputStringPointer += lengthOfCurrentSegment;

        } else if (readSegments[0][i].equals("P")) { // present in read sequence, but NOT in
                                                     // reference sequence. Similar to I.

        } else {
          System.err.println("Unrecognized CIGAR string: " + readSegments[0][i]);
        }
      }
    }

    if (reFormatedSeqString == null || reFormatedSeqString.equals("")) {
      return null;
    } else {
      return new String[] {beginPosOfOutputString + "", reFormatedSeqString, reFormatedPhredString,
                           aSingleLineOfBamFile[4]};
    }
  }

  /**
   *
   * @param reFormattedCigarString
   * @param begin
   * @param end inclusive
   * @param altAlleleDb "alleleId, Hash<position, alleleString>"
   */
  public static void updateRefAndAlt(String reFormattedCigarString, char[] ref, int begin,
                                     Hashtable<Integer, String> varsInCurrentString,
                                     Hashtable<Byte, Hashtable<Integer, String>> haplotypeDb,
                                     Hashtable<Byte, Vector<int[]>> readCounts) {
    // boolean isTheSame;
    byte haplotypeId;
    int end;
    reFormattedCigarString = reFormattedCigarString.toUpperCase();
    end = begin + reFormattedCigarString.length();
    haplotypeId = 0;
    // varsInCurrentString = new Hashtable<Integer, String>();

    compareWithRefAndGetListOfVars(reFormattedCigarString, ref, begin, varsInCurrentString);

    if (varsInCurrentString.size() > 0) {
      haplotypeId = compareWithHaplotypesAndGetHaplotypeId(begin, end, varsInCurrentString,
                                                           haplotypeDb, readCounts);
    }

    addToReadCounts(begin, reFormattedCigarString, haplotypeId, readCounts);
  }

  public static void compareWithRefAndGetListOfVars(String reFormattedCigarString, char[] ref,
                                                    int begin,
                                                    Hashtable<Integer, String> varsInCurrentString) {
    int end;

    end = begin + reFormattedCigarString.length();
    for (int i = begin; i < end; i++) {
      if (ref[i] == 0) {
        ref[i] = reFormattedCigarString.charAt(i - begin);
      } else if (ref[i] != reFormattedCigarString.charAt(i - begin)) {
        if (varsInCurrentString.containsKey(i)) {
          varsInCurrentString.put(i,
                                  varsInCurrentString.get(i)
                                     + (reFormattedCigarString.substring(i - begin,
                                                                         i - begin + 1)));
        } else {
          varsInCurrentString.put(i, reFormattedCigarString.substring(i - begin, i - begin + 1));
        }
      }
    }
  }

  public static byte compareWithHaplotypesAndGetHaplotypeId(int begin, int end,
                                                            Hashtable<Integer, String> varsInCurrentString,
                                                            Hashtable<Byte, Hashtable<Integer, String>> haplotypeDb,
                                                            Hashtable<Byte, Vector<int[]>> readCounts) {
    byte haplotypeId, haplotypeIdTmp;
    int[] varsNotInHaplotypeDbSorted, numMatches, sortedIndices;
    Hashtable<Integer, String> haplotypeDbTmp;
    Vector<Integer>[] varsNotInHaplotypeDbTmp;

    if (haplotypeDb.size() == 0) {
      haplotypeId = 1;
      haplotypeDb.put(haplotypeId, varsInCurrentString);
    } else {
      numMatches = new int[haplotypeDb.size()];
      varsNotInHaplotypeDbTmp = new Vector[haplotypeDb.size()];
      for (byte i = 0; i < haplotypeDb.size(); i++) {
        haplotypeDbTmp = haplotypeDb.get((byte) (i + 1));
        for (int j : varsInCurrentString.keySet()) {
          if (haplotypeDbTmp.containsKey(j)) {
            if (haplotypeDbTmp.get(j).equals(varsInCurrentString.get(j))) {
              numMatches[i]++;
            } else {
              numMatches[i] = -1;
              varsNotInHaplotypeDbTmp[i] = null;
              break;
            }
          } else {
            if (varsNotInHaplotypeDbTmp[i] == null) {
              varsNotInHaplotypeDbTmp[i] = new Vector<Integer>();
            }
            varsNotInHaplotypeDbTmp[i].add(j);
          }
        }
        for (int j : haplotypeDbTmp.keySet()) {
          if (j >= begin && j < end) {
            if (varsInCurrentString.containsKey(j)) {
              if (!haplotypeDbTmp.get(j).equals(varsInCurrentString.get(j))) {
                numMatches[i] = -1;
                varsNotInHaplotypeDbTmp[i] = null;
                break;
              }
            } else {
              if (varsNotInHaplotypeDbTmp[i] == null) {
                varsNotInHaplotypeDbTmp[i] = new Vector<Integer>();
              }
              varsNotInHaplotypeDbTmp[i].add(j);
            }
          }
        }
      }
      sortedIndices = Sort.quicksort(numMatches, Sort.DESCENDING);

      if (numMatches[sortedIndices[0]] == varsInCurrentString.size()) {
        haplotypeId = (byte) sortedIndices[0];
      } else {
        haplotypeId = -1;
        for (int i = 0; i < sortedIndices.length && numMatches[sortedIndices[i]] != -1; i++) {
          haplotypeIdTmp = (byte) (sortedIndices[i] + 1);
          haplotypeDbTmp = haplotypeDb.get(haplotypeIdTmp);

          if (varsNotInHaplotypeDbTmp[sortedIndices[i]].size() > 0) {
            varsNotInHaplotypeDbSorted = new int[varsNotInHaplotypeDbTmp[sortedIndices[i]].size()];
            for (int j = 0; j < varsNotInHaplotypeDbSorted.length; j++) {
              varsNotInHaplotypeDbSorted[j] =
                                            varsNotInHaplotypeDbTmp[sortedIndices[i]].elementAt(j);
            }
            Arrays.sort(varsNotInHaplotypeDbSorted);
            if (getNumReadsCoveringThisRegion(haplotypeIdTmp, varsNotInHaplotypeDbSorted[0],
                                              varsNotInHaplotypeDbSorted[varsNotInHaplotypeDbSorted.length
                                                                         - 1] - varsNotInHaplotypeDbSorted[0]
                                                                                             + 1,
                                              readCounts) < 1
                && haplotypeId == -1) {
              for (int element : varsNotInHaplotypeDbSorted) {
                if (varsInCurrentString.containsKey(element)) {
                  haplotypeDbTmp.put(element, varsInCurrentString.get(element));
                } else {
                  haplotypeDbTmp.put(element, haplotypeDbTmp.get(element));
                }
              }
              haplotypeId = haplotypeIdTmp;
              break;
            }
          }
        }

        if (haplotypeId == -1) {
          haplotypeId = (byte) (haplotypeDb.size() + 1);
          haplotypeDb.put(haplotypeId, varsInCurrentString);
        }
      }
    }

    return haplotypeId;
  }

  // TODO how to tell which ones are different, which ones are blank?
  public static void addToReadCounts(int begin, String reformattedCigarString, byte alleleId,
                                     Hashtable<Byte, Vector<int[]>> readCounts) {
    Vector<int[]> currentReadCount;
    int[] currentReadCountElement;
    int length;
    boolean found;

    length = reformattedCigarString.length();
    if (readCounts.containsKey(alleleId)) {
      found = false;
      currentReadCount = readCounts.get(alleleId);
      for (int i = 0; i < currentReadCount.size(); i++) {
        currentReadCountElement = currentReadCount.elementAt(i);
        if (currentReadCountElement[0] == begin && currentReadCountElement[1] == length) {
          currentReadCountElement[2]++;
          found = true;
          break;
        }
      }
      if (!found) {
        currentReadCount.add(new int[] {begin, length, 1});
      }
    } else {
      currentReadCount = new Vector<int[]>();
      currentReadCount.add(new int[] {begin, length, 1});
      readCounts.put(alleleId, currentReadCount);
    }
  }

  public static int getNumReadsCoveringThisRegion(byte alleleId, int begin, int length,
                                                  Hashtable<Byte, Vector<int[]>> readCounts) {
    Vector<int[]> currentReadCount;
    int[] currentReadCountElement;
    int end, counter;

    counter = 0;
    end = begin + length - 1;
    currentReadCount = readCounts.get(alleleId);
    for (int i = 0; i < currentReadCount.size(); i++) {
      currentReadCountElement = currentReadCount.elementAt(i);
      if ((begin >= currentReadCountElement[0]
           && begin <= (currentReadCountElement[0] + currentReadCountElement[1] - 1))
          || (end >= currentReadCountElement[0]
              && end <= (currentReadCountElement[0] + currentReadCountElement[1] - 1))) {
        counter++;
      }
    }

    return counter;
  }

  private static void exportHaplotypes(Hashtable<Byte, Hashtable<Integer, String>> altAlleles,
                                       Hashtable<Byte, Vector<int[]>> readCounts, char[] ref,
                                       int begin, String outputFileFullPath) {
    PrintWriter writer;
    HashSet<Integer> positions;
    Hashtable<Integer, String> currentAltAllele;
    Object[] tmp;
    String allele;
    Logger log;

    positions = new HashSet<Integer>();
    for (byte i = 1; i <= altAlleles.size(); i++) {
      currentAltAllele = altAlleles.get(i);
      tmp = currentAltAllele.keySet().toArray();
      for (Object element : tmp) {
        positions.add((Integer) element);
      }
    }
    tmp = positions.toArray();
    Arrays.sort(tmp);
    try {
      writer = new PrintWriter(new FileWriter(outputFileFullPath));
      writer.print("alleleId");
      for (Object element : tmp) {
        writer.print("\t" + (begin + (Integer) element));
      }
      writer.println();
      writer.print("ref");
      for (Object element : tmp) {
        writer.print("\t" + ref[(Integer) element]);
      }
      writer.println();
      for (byte i = 1; i <= altAlleles.size(); i++) {
        writer.print(i);
        currentAltAllele = altAlleles.get(i);
        for (Object element : tmp) {
          writer.print("\t");
          allele = currentAltAllele.get(element);
          if (allele != null) {
            writer.print(allele);
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      log = new Logger();
      log.reportError("Error writing to " + outputFileFullPath);
      log.reportException(e);
    }
  }

  public static void getCoverages(String dirReadCountsFiles, String suffixsOfReadCountsFiles,
                                  String bedfile, String thresholdsForReadCounts,
                                  String thresholdsForMapping, Logger log) {
    Segment[][] segments;
    String[] filenames;
    String trioIdOfTheReadCountsFile, root, result, dir;

    dir = ext.parseDirectoryOfFile(dirReadCountsFiles);
    root = ext.rootOf(dirReadCountsFiles);
    filenames = Files.list(dir, root, suffixsOfReadCountsFiles, true, false);
    segments = getSegments(bedfile, log);
    result = getHeader(thresholdsForMapping);
    for (String filename : filenames) {
      trioIdOfTheReadCountsFile = filename.substring(root.length()).split("[_.]")[0];
      result += assessCoverage(dir + filename, trioIdOfTheReadCountsFile, segments,
                               thresholdsForReadCounts, thresholdsForMapping, false, log);
      log.report("");
    }

    Files.write(result, dir + "coverage.txt");
    log.report("Output file is ready at " + dir + "coverage.txt");
  }

  public static String getCoverages(String fullpathToReadCountsFile,
                                    String trioIdOfTheReadCountsFile, String bedfile,
                                    String thresholdsForReadCounts, String thresholdsForMapping,
                                    boolean addHeaderInOutputFile, Logger log) {
    return assessCoverage(fullpathToReadCountsFile, trioIdOfTheReadCountsFile,
                          getSegments(bedfile, log), thresholdsForReadCounts, thresholdsForMapping,
                          false, log);
  }

  public static String assessCoverage(String fullpathToReadCountsFile,
                                      String trioIdOfTheReadCountsFile, Segment[][] segments,
                                      String thresholdsForReadCounts, String thresholdsForMapping,
                                      boolean addHeaderInOutputFile, Logger log) {
    BufferedReader reader;
    String[] line;
    String temp, logFilename = null;
    int count;

    byte chr, prevChr;
    int[] readCounts, mappingScores, thresholdsForReadCount, thresholdForMapping;
    int[][] runningCounts;
    int onTargetReads, numQualify;
    boolean allQualify, isLogFromCaller = true;
    long timer;
    SimpleDateFormat timeFormat;

    if (log == null) {
      logFilename = ext.parseDirectoryOfFile(fullpathToReadCountsFile) + "coverage_"
                    + ext.rootOf(fullpathToReadCountsFile) + ".out";
      log = new Logger(logFilename);
      isLogFromCaller = false;
    }

    timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
    timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
    timer = new Date().getTime();

    try {
      thresholdsForReadCount = Array.toIntArray(thresholdsForReadCounts.split(","));
      thresholdForMapping = Array.toIntArray(thresholdsForMapping.split(","));
    } catch (Exception e) {
      log.reportError("Error - failed to parse coverage thresholds: " + thresholdsForReadCounts);
      log.reportException(e);
      return null;
    }
    onTargetReads = 0;
    runningCounts = new int[thresholdsForReadCount.length][4 * thresholdForMapping.length];

    prevChr = -1;
    temp = "not even started";
    try {
      reader = Files.getAppropriateReader(fullpathToReadCountsFile);
      do {
        temp = reader.readLine();
      } while (!temp.startsWith(READ_COUNTS_HEADER_FORMAT2));

      if (!ext.checkHeader(temp.split("[\\s]+"), READ_COUNTS_HEADER_FORMAT2.split("\t"), false)) {
        log.reportError("Problem with header for file: " + fullpathToReadCountsFile);
        reader.close();
        return null;
      }

      count = 0;
      while (reader.ready()) {
        if (count % 10000000 == 0) {
          log.report(ext.getTimestampForFilename() + "\t" + count);
        }
        try {
          temp = reader.readLine();
        } catch (EOFException eofe) {
          log.report("Truncated at line: " + count);
          reader.close();
          continue;
        }
        line = temp.trim().split("\t");
        if (line.length < 5) {
          log.report("Truncated at line: " + count);
          continue;
        }
        chr = Positions.chromosomeNumber(line[0], log);
        if (chr != prevChr) {
          log.report(ext.getTimestampForFilename() + "\tStarting chromosome "
                     + Positions.CHR_CODES[chr]);
          prevChr = chr;
        }
        if (segments == null
            || (segments[chr] != null
                && Segment.overlapsAny(new Segment(chr, Integer.parseInt(line[1]),
                                                   Integer.parseInt(line[1])),
                                       segments[chr]))) {
          readCounts = new int[] {(line[2].equals("") ? 0 : Integer.parseInt(line[2])),
                                  (line[3].equals("") ? 0 : Integer.parseInt(line[3])),
                                  (line[4].equals("") ? 0 : Integer.parseInt(line[4]))};
          mappingScores = new int[] {(line[5].equals("") ? 0 : Integer.parseInt(line[5])),
                                     (line[6].equals("") ? 0 : Integer.parseInt(line[6])),
                                     (line[7].equals("") ? 0 : Integer.parseInt(line[7]))};
          for (int i = 0; i < thresholdsForReadCount.length; i++) {
            for (int j = 0; j < thresholdForMapping.length; j++) {
              allQualify = true;
              numQualify = 0;
              for (int k = 0; k < 3; k++) {
                if (readCounts[k] > thresholdsForReadCount[i]
                    && mappingScores[k] >= thresholdForMapping[j]) {
                  runningCounts[i][k + j * 4]++;
                  numQualify++;
                } else {
                  allQualify = false;
                }
              }
              if (allQualify) {
                runningCounts[i][3 + j * 4]++;
              }
              if (numQualify == 0) {
                break;
              }
            }
          }
          onTargetReads++;
        }
        count++;
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + fullpathToReadCountsFile
                         + "\" not found in current directory");
      System.exit(1);
    } catch (Exception e) {
      System.err.println("Error reading file \"" + fullpathToReadCountsFile + "\" at line:");
      System.err.println(temp);
      e.printStackTrace();
    }

    if (trioIdOfTheReadCountsFile == null) {
      trioIdOfTheReadCountsFile = ext.rootOf(fullpathToReadCountsFile).split("_")[0];
    }
    if (addHeaderInOutputFile) {
      temp = getHeader(thresholdsForMapping);
    } else {
      temp = "";
    }
    for (int i = 0; i < thresholdsForReadCount.length; i++) {
      // temp += ((i == 0? "" : "\n") + trioIdOfTheReadCountsFile + "\t>=" +
      // thresholdsForReadCount[i]);
      temp += ("\n" + trioIdOfTheReadCountsFile + "\t>=" + thresholdsForReadCount[i]);
      for (int j = 0; j < thresholdForMapping.length; j++) {
        for (int k = 0; k < 4; k++) {
          temp += ("\t"
                   + ext.formPercent((double) runningCounts[i][k + j * 4] / (double) onTargetReads,
                                     1));
        }
      }
    }
    log.report(getHeader(thresholdsForMapping) + temp);

    if (!isLogFromCaller) {
      log.report("Result is also saved at " + logFilename + ".out");
    }

    log.report("Time used for trio " + trioIdOfTheReadCountsFile + "\t"
               + timeFormat.format(new Date().getTime() - timer));
    return temp;
  }

  public static Segment[][] getSegments(String bedfile, Logger log) {
    SegmentLists segmentLists;
    Segment[][] segments;

    if (log == null) {
      log = new Logger();
    }

    if (bedfile == null) {
      segments = null;
    } else if (Files.exists(bedfile + ".ser")) {
      log.report("Reading in preserialized " + bedfile + ".ser");
      segments = SegmentLists.load(bedfile + ".ser", false).getLists();
    } else {
      log.report("Importing " + bedfile);
      segmentLists = SegmentLists.parseSegmentList(bedfile, 0, 1, 2, true);
      segmentLists.serialize(bedfile + ".ser");
      segments = segmentLists.getLists();
    }

    return segments;
  }

  public static String getHeader(String thresholdsForMapping) {
    String temp;
    int[] thresholdForMapping;

    thresholdForMapping = Array.toIntArray(thresholdsForMapping.split(","));
    temp = "TrioId\t#Reads_Threshold";
    for (int element : thresholdForMapping) {
      temp += "\tMapping>=" + element + "_Child\tMapping>=" + element + "_Dad\tMapping>=" + element
              + "_Mom\tMapping>=" + element + "_Total";
    }

    return temp;
  }

  public static void test() {
    loadAdditionalComments("N:/statgen/OS_Logan/SuperNovo/charge_fibrinogen_mafs_and_counts.xln.gz");

    // try {
    // BufferedReader reader = new BufferedReader(new
    // FileReader("D:/logan/DeNovos/bams/test2.txt"));
    // String line;
    // String[] results;
    // Hashtable<Integer, String> tmp1 = new Hashtable<Integer, String>(), tmp2 = new
    // Hashtable<Integer, String>();
    // while (reader.ready()) {
    // line = reader.readLine();
    //// results = cigarToReFormatedString(line.split("\t"), 39346524, 39346623, null);
    //// results = cigarToReFormatedString(line.split("\t"), 39346596, 39346600, tmp1, tmp2); //for
    // Insertion at the 1st place
    // results = cigarToReFormatedString(line.split("\t"), 105815200 - 60, 105815204, tmp1, tmp2);
    // for (int i = 0; results != null && i < results.length; i++) {
    // System.out.println(results[i]);
    // }
    // }
    // } catch (FileNotFoundException e) {
    // e.printStackTrace();
    // } catch (IOException e) {
    // e.printStackTrace();
    // }

    // try {
    // BufferedReader reader = new BufferedReader(new
    // FileReader("D:/logan/DeNovos/tests/haplotypes.log"));
    // String[] line;
    // char[] ref = new char[300];
    // Hashtable<Integer, String> varsInCurrentString;
    // Hashtable<Byte, Hashtable<Integer, String>> haplotypeDb = new Hashtable<Byte,
    // Hashtable<Integer, String>>();
    // Hashtable<Byte, Vector<int[]>> readCounts = new Hashtable<Byte, Vector<int[]>>();
    // while (reader.ready()) {
    // varsInCurrentString = new Hashtable<Integer, String>();
    // line = reader.readLine().split("\t");
    // if (line[2].endsWith("AGGGCTTCCT")) {
    // System.out.println();
    // }
    // updateRefAndAlt(line[2], ref, Integer.parseInt(line[1]), varsInCurrentString, haplotypeDb,
    // readCounts);
    // }
    // reader.close();
    // } catch (FileNotFoundException e) {
    // e.printStackTrace();
    // } catch (IOException e) {
    // e.printStackTrace();
    // }

    // Hashtable<Integer, Vector<int[]>> test = new Hashtable<Integer, Vector<int[]>>();
    // test.put(53700, new Vector<int[]>());
    // test.get(53700).add(new int[] {53, 1, 9});
    // test.get(53700).add(new int[] {93, 0, 13});
    // test.put(63700, new Vector<int[]>());
    // test.get(63700).add(new int[] {66, 1, 69});
    // test.get(63700).add(new int[] {66, 0, 63});
    // int[] tmp = test.get(53700).elementAt(0);
    // for (int i = 0; i < tmp.length; i++) {
    // System.out.println(tmp[i]);
    // }
    // tmp = test.get(63700).elementAt(1);
    // for (int i = 0; i < tmp.length; i++) {
    // System.out.println(tmp[i]);
    // }

    // Hashtable<Byte, Hashtable<Integer, String>> altAlleles = new Hashtable<Byte,
    // Hashtable<Integer, String>>();
    // char[] ref = new char[10];
    // Hashtable<Byte, Vector<int[]>> readCounts = new Hashtable<Byte, Vector<int[]>>();
    // addOneReadToAlleleStringsAndCounts("AAAAAA", ref, 2, altAlleles, readCounts);
    // addOneReadToAlleleStringsAndCounts("BBB", ref, 0, altAlleles, readCounts);
    // addOneReadToAlleleStringsAndCounts("C", ref, 2, altAlleles, readCounts);
    // addOneReadToAlleleStringsAndCounts("C", ref, 2, altAlleles, readCounts);
    // addOneReadToAlleleStringsAndCounts("C", ref, 2, altAlleles, readCounts);
    // addOneReadToAlleleStringsAndCounts("BCA", ref, 1, altAlleles, readCounts);
    // addOneReadToAlleleStringsAndCounts("D", ref, 1, altAlleles, readCounts);
    // addOneReadToAlleleStringsAndCounts("D", ref, 9, altAlleles, readCounts);
    // addOneReadToAlleleStringsAndCounts("DC", ref, 1, altAlleles, readCounts);
    // addOneReadToAlleleStringsAndCounts("BAA", ref, 1, altAlleles, readCounts);
  }

  public static void fromParameters(String filename, Logger log) {
    Vector<String> params;

    params = Files.parseControlFile(filename, "SuperNovo",
                                    new String[] {"file=snps.txt",
                                                  "# set the program to parse results",
                                                  "-parseresult",
                                                  "# input: directory of the 1st phase results",
                                                  "outputdir=",
                                                  "# input: full path to the trio list file (leave as null, if no mini bam files are needed)",
                                                  "triolistfile=",
                                                  "# directory to output the SeattleSeq .needs files",
                                                  "seattleseqdir=",
                                                  "# directory to output the script files for generating the mini bam files",
                                                  "minibamdir="},
                                    log);

    if (params != null) {
      params.add("log=" + log.getFilename());
      main(Array.toStringArray(params));
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String[] fullpathsToSamFilenamesOfTheTrio;
    String commandIsToScanSamFilesForDenovoMutations, commandIsToFilterAndSummarizeResults,
        commandIsToHaploidCount, commandIsToOutputHaplotypeStrings, commandAnnotation,
        commandOutDir, commandTrioList, commandDirScripts, commandDirSummary, commandCandidate,
        commandBamDir, commandBamSet, commandRefFasta, commandChr, commandStart, commandStop,
        commandBim, commandNumThreads, commandTrioId, commandSeattleSeqDir, commandMiniBamScriptDir,
        commandOutputReadCounts, commandIsToAssessCoverage, commandReadCountsDir,
        commandReadCountsThreshold, commandMappingThreshold, commandHaploTypeOut, commandMiniBamDir,
        commandIsToTestExperiementalProgram;
    String helpMenu, fullPathVariantCandidatesList, fullPathRefFasta, fullPathBim,
        fullPathReadCounts, fullPathTrioList, fullPathToSummaryOutputFile,
        fullPathToOutputHaplotypeStrings;
    String dirSam, dirScript, dirOutputsFromTheScanningOfSamFilesForDenovoMutations, dirSeattleSeq,
        dirMiniSam, dirMiniSamScripts, dirReadCounts;
    String trioId, thresholdsForReadCounts, thresholdsForMapping;
    String chr;
    int regionLegnthATime, numThreads;
    int begin, end;
    boolean isToAnnotate, isToSummarizeResults, isToAssessCoverage,
        isToScanSamFilesForDenovoMutations, isGetHaplotypes, isToTest, isToOutputHaplotypeStrings;
    Logger log;

    fullPathVariantCandidatesList = "N:/statgen/OS_Logan/IGV_validation/results_test1.txt";
    dirSam = "/home/spectorl/shared/exome_processing/bam/";
    // dirSam = "/home/spectorl/shared/exome/project126/all_101/";
    fullpathsToSamFilenamesOfTheTrio =
                                     new String[] {"/home/spectorl/shared/exome_processing/bam/rrd_F10639C_L008.bam",
                                                   "/home/spectorl/shared/exome_processing/bam/rrd_F10639D_L007.bam",
                                                   "/home/spectorl/shared/exome_processing/bam/rrd_F10639M_L007.bam"};
    if (fullpathsToSamFilenamesOfTheTrio != null
        && !fullpathsToSamFilenamesOfTheTrio[0].equals("")) {
      trioId = getRootOf(fullpathsToSamFilenamesOfTheTrio, false);
    } else {
      trioId = null;
    }
    fullPathTrioList = "/home/spectorl/xuz2/lists/trios_list.txt";
    fullPathRefFasta = "/home/spectorl/xuz2/ref_fasta/hg19_canonical.fa";
    fullPathBim = "/home/spectorl/xuz2/outputs/wholegenome.bim";
    fullPathReadCounts = "/home/spectorl/xuz2/readcounts/wholegenome_rrd/readCount_F10639.txt.gz";
    fullPathToSummaryOutputFile =
                                "/home/spectorl/xuz2/outputs/wholegenome_rrd_46/SuperNovo_summary.xln";
    fullPathToOutputHaplotypeStrings = "/home/spectorl/xuz2/outputs/tests/haplotypes.txt";
    dirScript = "/home/spectorl/xuz2/scripts/";
    dirOutputsFromTheScanningOfSamFilesForDenovoMutations = "/home/spectorl/xuz2/outputs/";
    dirReadCounts = "/home/spectorl/xuz2/readcounts/wholegenome_rrd/readcounts_";
    dirSeattleSeq = "/home/spectorl/xuz2/outputs/SeattleSeqAnnotation/";
    dirMiniSam = "/home/spectorl/xuz2/outputs/mini_bams/";
    dirMiniSamScripts = "/home/spectorl/xuz2/outputs/mini_bams/";
    regionLegnthATime = -1;
    numThreads = 1;
    chr = "";
    begin = 39346600;
    end = 39346622;
    thresholdsForReadCounts = "8,10,12";
    thresholdsForMapping = "0,10,15,20,25,40,50";
    isToAnnotate = false;
    isToScanSamFilesForDenovoMutations = false;
    isToSummarizeResults = false;
    isToAssessCoverage = false;
    isGetHaplotypes = false;
    isToTest = false;
    isToOutputHaplotypeStrings = false;

    commandIsToScanSamFilesForDenovoMutations = "-denovo";
    commandIsToFilterAndSummarizeResults = "-summarizeresults";
    commandIsToHaploidCount = "-haploidcount";
    commandIsToAssessCoverage = "-coverage";
    commandIsToTestExperiementalProgram = "-test";
    commandIsToOutputHaplotypeStrings = "outhaplotypes=";
    commandAnnotation = "-annotation";
    commandTrioList = "triolistfile=";
    commandOutDir = "outdir=";
    commandDirScripts = "scriptdir=";
    commandDirSummary = "parsedresult=";
    commandBamDir = "bamdir=";
    commandBamSet = "bamset=";
    commandRefFasta = "reffasta=";
    commandChr = "chr=";
    commandStart = "start=";
    commandStop = "stop=";
    commandBim = "bim=";
    commandNumThreads = "numthreads=";
    commandTrioId = "trioid=";
    commandSeattleSeqDir = "seattleseqdir=";
    commandMiniBamScriptDir = "minibamscriptdir=";
    commandMiniBamDir = "minibamdir=";
    commandOutputReadCounts = "outputreadcounts=";
    commandReadCountsDir = "readcountsdir=";
    commandReadCountsThreshold = "thresholdreadcount=";
    commandMappingThreshold = "thresholdmapping=";
    commandHaploTypeOut = "haplotypeout=";

    commandCandidate = "candidate=";

    helpMenu =
             (fullpathsToSamFilenamesOfTheTrio == null
              || fullpathsToSamFilenamesOfTheTrio[0] == null) ? ""
                                                              : fullpathsToSamFilenamesOfTheTrio[0];
    for (int i = 1; fullpathsToSamFilenamesOfTheTrio != null
                    && i < fullpathsToSamFilenamesOfTheTrio.length
                    && fullpathsToSamFilenamesOfTheTrio[i] != null; i++) {
      helpMenu += "," + fullpathsToSamFilenamesOfTheTrio[i];
    }

    String usage =
                 "SuperNovo analyzes DNA sequencing data of family studies for De Novo variants through 2 steps (phases):"
                   + "\nPhase 1: Scan the .bam/.sam file set of each family (trio) for candidate SNPs of De Novo variants, and output a text file for each family (trio);"
                   + "\nPhase 2: Scan all the output files from phase 1 for all the families (trios) in the project, add gene function annotations, and output a text summary file."
                   + "\n" + "\nNote:"
                   + "\n   .sam/.bam file names in the following commands and the trio list file must follow the order of Child, Dad, and Mom."
                   + "\nTo scan one region in .sam or .bam files of a trio for De Novo Mutations (phase 1 scan):"
                   + "\n   (1) " + commandIsToScanSamFilesForDenovoMutations + " (not default))"
                   + "\n   (2) full paths to the .sam (or .bam) files of the trio (i.e. "
                   + commandBamSet + helpMenu + " (not default))"
                   + "\n   (3) label of the trio, to be used in the output file for identification (i.e. "
                   + commandTrioId + trioId + " (default))"
                   + "\n   (4) chromosome number of the region (i.e. " + commandChr + chr
                   + " (not default))" + "\n   (5) beginning position of the region (i.e. "
                   + commandStart + begin + " (not default))"
                   + "\n   (6) ending position of the region (i.e. " + commandStop + end
                   + " (not default))" + "\n   (7) full path to the reference Fasta file (i.e. "
                   + commandRefFasta + fullPathRefFasta + " (default))"
                   + "\n   (8) directory for the output files (i.e. " + commandOutDir
                   + dirOutputsFromTheScanningOfSamFilesForDenovoMutations + " (default))"
                   + "\n   (9) to output read counts, full path to the read counts output file (i.e. "
                   + commandOutputReadCounts + fullPathReadCounts + " (not default))" + "\n"
                   + "\nTo scan multiple regions in .sam or .bam files of a trio for De Novo Mutations (phase 1 scan):"
                   + "\n   (1) " + commandIsToScanSamFilesForDenovoMutations + " (not default))"
                   + "\n   (2) full paths to the .sam (or .bam) files of the trio (i.e. "
                   + commandBamSet + helpMenu + " (not default))"
                   + "\n   (3) label of the trio, to be used in the output file for identification (i.e. "
                   + commandTrioId + trioId + " (default))"
                   + "\n   (4) full path of the .bim file (i.e. " + commandBim + fullPathBim
                   + " (default))" + "\n   (5) full path of the reference Fasta file (i.e. "
                   + commandRefFasta + fullPathRefFasta + " (default))"
                   + "\n   (6) directory for the output files (i.e. " + commandOutDir
                   + dirOutputsFromTheScanningOfSamFilesForDenovoMutations + " (default))"
                   + "\n   (7) to output read count files, directory and prefix for the read counts files (i.e. "
                   + commandReadCountsDir + dirReadCounts + " (default)"
                   + "\n   (8) maximum length that a region is to be broken into, for the scan of .sam/.bam files at a time (i.e. "
                   + regionLegnthATime + "\n   (9) number of threads (i.e. " + commandNumThreads
                   + numThreads + " (default))" + "\n"
                   + "\nTo generate scripts for all the trios in a directory, to scan multiple regions in .sam or .bam files for each one of them at a time for De Novo Mutations (phase 1 scan):"
                   + "\n   (1) " + commandIsToScanSamFilesForDenovoMutations + " (not default))"
                   + "\n   (2) directory of the .sam files (i.e. " + commandBamDir + dirSam
                   + " (default))" + "\n   (3) full path to the trio list file (i.e. "
                   + commandTrioList + fullPathTrioList + " (default))"
                   + "\n   (4) full path of the .bim file (i.e. " + commandBim + fullPathBim
                   + " (default))" + "\n   (5) full path of the reference Fasta file (i.e. "
                   + commandRefFasta + fullPathRefFasta + " (default))"
                   + "\n   (6) directory for the output files (i.e. " + commandOutDir
                   + dirOutputsFromTheScanningOfSamFilesForDenovoMutations + " (default))"
                   + "\n   (7) directory for the script files (i.e. " + commandDirScripts
                   + dirScript
                   + "\n   (8) to output read count files, directory and prefix for the read counts files (i.e. "
                   + commandReadCountsDir + dirReadCounts + " (default)"
                   + "\n   (9) maximum length that a region is to be broken into, for the scan of .sam/.bam files at a time (i.e. "
                   + regionLegnthATime + "\n   (10) number of threads (i.e. " + commandNumThreads
                   + numThreads + " (default))" + "\n"
                   + "\nTo filter and summarize the outputs of the .sam files' scanning of all trios (phase 2 scan):"
                   + "\nYou might need to repeat this command for several times, after each of which you'll use the updated lists to download the annotations from SeattleSeq and output of this summary will be based on the updated SeattleSeq annotations before the run."
                   + "\n   (1) command for filter and summarize results (i.e. "
                   + commandIsToFilterAndSummarizeResults + " (default))"
                   + "\n   (2) directory for the output files of the scan of .sam files (phase 1 scan) (i.e. "
                   + commandOutDir + dirOutputsFromTheScanningOfSamFilesForDenovoMutations
                   + " (default))"
                   + "\n   (3) full path to the trio list file (to be used as the source for trio IDs) (i.e. "
                   + commandTrioList + fullPathTrioList
                   + " (default), can be left as null if no mini bam script is needed)"
                   + "\n   (4) directory of the SeattleSeq annotation files (existing annotation files will be incorporated into the summary, while new list of the additional files needed will be generated) (i.e. "
                   + commandSeattleSeqDir + dirSeattleSeq + " (default))"
                   + "\n   (5) directory to output the scripts for generating the mini .bam files (to call samtools and extract regions from the .sam/.bam files) (i.e. "
                   + commandMiniBamScriptDir + dirMiniSamScripts + " (default))"
                   + "\n   (6) directory to output the mini .bam files (i.e. " + commandMiniBamDir
                   + dirMiniSam + " (default))" + "\n" + "\nTo assess coverage of a specific trio:"
                   + "\n   (1) command for assess coverage (i.e. " + commandIsToAssessCoverage
                   + " (not default))"
                   + "\n   (2) full path to the read counts file of the trio (i.e. "
                   + commandOutputReadCounts + fullPathReadCounts + " (default)"
                   + "\n   (3) full path of the bed file (i.e. " + commandBim + fullPathBim
                   + " (default))" + "\n   (4) thresholds of read counts (i.e. "
                   + commandReadCountsThreshold + thresholdsForReadCounts + " (default))"
                   + "\n   (5) thresholds of mapping quality scores (i.e. "
                   + commandMappingThreshold + thresholdsForMapping + " (default))" + "\n"
                   + "\nTo assess coverage of all trios in a directory:"
                   + "\n   (1) command for assess coverage (i.e. " + commandIsToAssessCoverage
                   + " (not default))" + "\n   (2) dir the read counts files (i.e. "
                   + commandReadCountsDir + dirReadCounts + " (default)"
                   + "\n   (3) full path of the bed file (i.e. " + commandBim + fullPathBim
                   + " (default))" + "\n   (4) thresholds of read counts (i.e. "
                   + commandReadCountsThreshold + thresholdsForReadCounts + " (default))"
                   + "\n   (5) thresholds of mapping quality scores (i.e. "
                   + commandMappingThreshold + thresholdsForMapping + " (default)" + "\n"
                   + "\nTo get haploid counts for a region of a trio (basing on output of the phase 2 parsing):" // TODO
                                                                                                                 // to
                                                                                                                 // do
                                                                                                                 // haploid
                                                                                                                 // count
                                                                                                                 // for
                                                                                                                 // each
                                                                                                                 // individual,
                                                                                                                 // rather
                                                                                                                 // than
                                                                                                                 // trio
                   + "\n   (1) command for haplotypes (i.e. " + commandIsToHaploidCount
                   + " (not default))"
                   + "\n   (2) label of the trio, to be used in the name of the output file (i.e. "
                   + commandTrioId + trioId + " (default))"
                   + "\n   (3) full paths to the .sam (or .bam) files of the trio (i.e. "
                   + commandBamSet + helpMenu + " (not default))"
                   + "\n   (4) chr of the region (i.e. " + commandChr + chr + " (default))"
                   + "\n   (5) start position of the region (i.e. " + commandStart + begin
                   + " (default))" + "\n   (6) stop position of the region (i.e. " + commandStop
                   + end + " (default))"
                   + "\n   (7) full path to output the haploytype strings, if needed (i.e. "
                   + commandHaploTypeOut + fullPathToOutputHaplotypeStrings + " (default))" + "\n"
                   + "\nTo get haploid counts of all trios appeared in an output file of the filter and summary (output of the phase 2 scan):"
                   + "\n   (1) command for haplotypes (i.e. " + commandIsToHaploidCount
                   + " (not default))"
                   + "\n   (2) full path to the output file of the filter and summary (output of the phase 2 scan) (i.e. "
                   + commandTrioList + fullPathToSummaryOutputFile + " (default))"
                   + "\n   (3) full path to the trio list file (i.e. " + commandTrioList
                   + fullPathTrioList + " (default))" + "\n   (4) directory of the bam files (i.e. "
                   + commandBamDir + dirSam + " (default))"
                   + "\n   (5) whether to output the haplotype strings (i.e. "
                   + commandIsToOutputHaplotypeStrings + isToOutputHaplotypeStrings + " (default))"
                   + "\n" + "\n(Legacy, obsoleted code): To annotate a list of candidate markers:"
                   + "\n   (1) command for annotatation (i.e. " + commandAnnotation + " (default))"
                   + "\n   (2) full path of the candidate list file (i.e. " + commandCandidate
                   + fullPathVariantCandidatesList + " (default))"
                   + "\n   (3) directory of the bam files (i.e. " + commandBamDir + dirSam
                   + " (default))"
                   + "\n   (4) directory for the output script file to extract information from bam files (i.e. "
                   + commandDirScripts + dirScript + " (default))" + "";

    fullpathsToSamFilenamesOfTheTrio = null;
    chr = null;
    begin = -1;
    end = -1;
    dirReadCounts = null;
    fullPathReadCounts = null;

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith(commandIsToScanSamFilesForDenovoMutations)) {
        isToScanSamFilesForDenovoMutations = true;
        numArgs--;
      } else if (arg.startsWith(commandAnnotation)) {
        isToAnnotate = true;
        numArgs--;
      } else if (arg.startsWith(commandCandidate)) {
        fullPathVariantCandidatesList = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandBamDir)) {
        dirSam = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandDirScripts)) {
        dirScript = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandOutDir)) {
        dirOutputsFromTheScanningOfSamFilesForDenovoMutations = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandBamSet)) {
        if (arg.split("=").length < 2) {
          fullpathsToSamFilenamesOfTheTrio = null;
        } else {
          fullpathsToSamFilenamesOfTheTrio = arg.split("=")[1].split(",");
        }
        numArgs--;
      } else if (arg.startsWith(commandRefFasta)) {
        fullPathRefFasta = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandChr)) {
        if (arg.split("=").length < 2) {
          chr = null;
        } else {
          chr = arg.split("=")[1];
        }
        numArgs--;
      } else if (arg.startsWith(commandStart)) {
        begin = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith(commandStop)) {
        end = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith(commandBim)) {
        if (arg.split("=").length < 2) {
          fullPathBim = null;
        } else {
          fullPathBim = arg.split("=")[1];
        }
        numArgs--;
      } else if (arg.startsWith(commandNumThreads)) {
        numThreads = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith(commandIsToFilterAndSummarizeResults)) {
        isToSummarizeResults = true;
        numArgs--;
      } else if (arg.startsWith(commandTrioId)) {
        trioId = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandTrioList)) {
        fullPathTrioList = arg.split("=")[1];
        numArgs--;
        // } else if (args[i].startsWith(commandRegionLength)) {
        // regionLegnthATime = Integer.parseInt(args[i].split("=")[1]);
        // numArgs--;
      } else if (arg.startsWith(commandSeattleSeqDir)) {
        dirSeattleSeq = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandMiniBamScriptDir)) {
        dirMiniSamScripts = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandOutputReadCounts)) {
        fullPathReadCounts = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandIsToAssessCoverage)) {
        isToAssessCoverage = true;
        numArgs--;
      } else if (arg.startsWith(commandReadCountsDir)) {
        dirReadCounts = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandReadCountsThreshold)) {
        thresholdsForReadCounts = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandMappingThreshold)) {
        thresholdsForMapping = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandIsToTestExperiementalProgram)) {
        isToTest = true;
        numArgs--;
      } else if (arg.startsWith(commandIsToHaploidCount)) {
        isGetHaplotypes = true;
        numArgs--;
      } else if (arg.startsWith(commandHaploTypeOut)) {
        fullPathToOutputHaplotypeStrings = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandIsToOutputHaplotypeStrings)) {
        isToOutputHaplotypeStrings = Boolean.parseBoolean(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith(commandDirSummary)) {
        fullPathToSummaryOutputFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(commandMiniBamDir)) {
        dirMiniSam = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    // dirBam = "D:/logan/DeNovos/mini_bams/";
    // fullPathBed = "D:/logan/DeNovos/beds/regions_chr17test.bed";
    // fullPathRefFasta = "D:/logan/DeNovos/hg19_canonical.fa";
    // dirDenovoVars = "D:/logan/DeNovos/outputs/tests/";
    // bamFilenamesOfTheTrio = new String[] {"F10639_chr17_38000000_C.bam",
    // "F10639_chr17_38000000_D.bam", "F10639_chr17_38000000_M.bam"};
    // trioId = "F10639_mini";

    isToSummarizeResults = true;
    // dirDenovoVars = "D:/logan/DeNovos/phase1outputs_denovos_rrd42_rrd46_dedup4/";
    dirOutputsFromTheScanningOfSamFilesForDenovoMutations =
                                                          "D:/logan/DeNovos/phase1outputs_inherited/";
    dirSeattleSeq = "N:/statgen/OS_Logan/SuperNovo/SeattleSeqAnnotation/";
    dirMiniSamScripts = "D:/logan/DeNovos/mini_bam_scripts/";
    dirMiniSam = "D:/logan/DeNovos/mini_bams/";
    fullPathTrioList = "N:/statgen/OS_Logan/SuperNovo/beds/triolist_rrd42_rrd46-6_dedup4.txt";

    // isGetHaplotypes = true;
    // fullPathToParsedResult =
    // "/home/spectorl/xuz2/outputs/wholegenome_rrd_46/SuperNovo_summary.xln";
    // fullPathTrioList = "/home/spectorl/xuz2/lists/triolist_rrd_46_dedup.txt";
    // isToOutputHaplotypeStrings = false;

    // isGetHaplotypes = true;
    // chr = 11 + "";
    // dirSam = "/home/spectorl/shared/exome/project126/all_101/";
    // bamFilenamesOfTheTrio = new String[] {"dedup_F10476C.bam", "dedup_F10476D.bam",
    // "dedup_F10476M.bam"};
    // trioId = "dedup_F10476";
    // begin = 130079255;
    // end = 130079459;
    // fullPathToOutputHaplotypeStrings =
    // "~/outputs/tests/dedup_F10476_chr11_130079357_haplotypes.txt";

    // isGetHaplotypes = true;
    // chr = 17 + "";
    // dirSam = "D:/logan/DeNovos/mini_bams/";
    // bamFilenamesOfTheTrio = new String[] {"F10503_chr17_7577102_C.sam",
    // "F10503_chr17_7577102_D.sam", "F10503_chr17_7577102_M.sam"};
    // trioId = "dedup_F10476";
    // begin = 7577102;
    // end = 7577102;
    // fullPathToOutputHaplotypeStrings = "D:/logan/DeNovos/tests/haplotypes.txt";

    // isToTest = true;

    if (isToScanSamFilesForDenovoMutations) {
      log = new Logger(dirOutputsFromTheScanningOfSamFilesForDenovoMutations + "SuperNovo_"
                       + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log");
      if (chr != null && !chr.equals("") && begin > 0 && end > 0 && end >= begin
          && (fullPathBim == null || fullPathBim.equals(""))) {
        log =
            new Logger(dirOutputsFromTheScanningOfSamFilesForDenovoMutations + "SuperNovo_"
                       + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log");
        scanARegionInSamFilesOfATrioForDeNovoMutations(fullpathsToSamFilenamesOfTheTrio,
                                                       fullPathRefFasta, chr, begin, end,
                                                       dirOutputsFromTheScanningOfSamFilesForDenovoMutations,
                                                       fullPathReadCounts, log);

      } else if (!Strings.isNullOrEmpty(fullPathBim) && new File(fullPathBim).exists()
                 && Strings.isNullOrEmpty(chr) && Strings.isNullOrEmpty(dirSam)) {
        scanMultipleRegionsInSamFilesOfATrioForDenovoMutations(trioId,
                                                               fullpathsToSamFilenamesOfTheTrio,
                                                               fullPathRefFasta, fullPathBim,
                                                               dirOutputsFromTheScanningOfSamFilesForDenovoMutations,
                                                               fullPathReadCounts,
                                                               regionLegnthATime, numThreads, log);

      } else if (!Strings.isNullOrEmpty(dirSam) && fullpathsToSamFilenamesOfTheTrio == null) {
        generateScriptsAllTriosInADirectoryToScanSamFilesForDenovoMutation(dirSam, fullPathRefFasta,
                                                                           fullPathBim,
                                                                           dirOutputsFromTheScanningOfSamFilesForDenovoMutations,
                                                                           dirScript,
                                                                           fullPathTrioList,
                                                                           dirReadCounts,
                                                                           regionLegnthATime,
                                                                           numThreads, log);

      } else {
        log.reportError("Error - for " + isToScanSamFilesForDenovoMutations
                        + ", due to some of the following:\n 1) chr and .bim cannot be specified at the same time;\n 2) "
                        + commandBamDir + " and "
                        + "\n 3) begin is less than end, or any one of them being negative;\n 4) file does not exist: "
                        + fullPathBim);

      }
    } else if (isToSummarizeResults) {
      log = new Logger(dirOutputsFromTheScanningOfSamFilesForDenovoMutations + "SuperNovo_"
                       + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log");
      // parseResults(dirOutputsFromTheScanningOfSamFilesForDenovoMutations, dirSeattleSeq,
      // dirMiniSamScripts, dirMiniSam, fullPathTrioList, (byte) 2, 0, (byte) 1, log);
      parseResults2(dirOutputsFromTheScanningOfSamFilesForDenovoMutations, dirSeattleSeq,
                    "N:/statgen/OS_Logan/SuperNovo/charge_fibrinogen_mafs_and_counts.xln.gz",
                    dirMiniSamScripts, dirMiniSam, fullPathTrioList, (byte) 2, 0, (byte) 1, log);

    } else if (isToAssessCoverage) {
      if ((dirReadCounts != null && !dirReadCounts.equals("")) || fullPathReadCounts == null
          || fullPathReadCounts.equals("")) {
        log =
            new Logger(dirReadCounts + "SuperNovo_"
                       + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log");
        getCoverages(dirReadCounts, DEFAULT_READ_COUNTS_FILENAME_SUFFIX, fullPathBim,
                     thresholdsForReadCounts, thresholdsForMapping, log);

      } else {
        log =
            new Logger(ext.parseDirectoryOfFile(fullPathReadCounts) + "SuperNovo_"
                       + new SimpleDateFormat("yyyy.MM.dd_hh.mm.ssa").format(new Date()) + ".log");
        getCoverages(fullPathReadCounts, null, fullPathBim, thresholdsForReadCounts,
                     thresholdsForMapping, true, log);

      }
    } else if (isGetHaplotypes) {
      if (chr == null) {
        log = new Logger(ext.parseDirectoryOfFile(fullPathToSummaryOutputFile) + "haplotypes.log");
        getHaplotypes(fullPathToSummaryOutputFile, fullPathTrioList, dirSam,
                      isToOutputHaplotypeStrings, log);

      } else {
        if (fullPathToOutputHaplotypeStrings == null) {
          log = new Logger();
        } else {
          log = new Logger(ext.parseDirectoryOfFile(fullPathToOutputHaplotypeStrings)
                           + "haplotypes.log");
        }
        log.report("Number of haplotypes: "
                   + getHaplotypes(fullpathsToSamFilenamesOfTheTrio, trioId, chr, begin, end,
                                   fullPathToOutputHaplotypeStrings, log));
      }
    } else if (isToTest) {
      log = new Logger();
      test();

    } else if (isToAnnotate) {
      // getAlleleCounts("D:/bams/F10639C_chr11_89531764.txt", 11, 89531764, 3); //Testing only.
      // Feature already included in screenDeNovoPointMutation()
      // getAlleleCounts("D:/bams/F10639C_chr5_140558628.txt", 5, 140558628, 3);
      generateScriptForSamtools(fullPathVariantCandidatesList, fullpathsToSamFilenamesOfTheTrio[0],
                                dirScript);
      screenDeNovoPointMutation(fullPathVariantCandidatesList, fullpathsToSamFilenamesOfTheTrio[0],
                                1);

    } else {
      System.out.println("\nNo command executed, due to none of the following is specified: "
                         + commandIsToScanSamFilesForDenovoMutations + ", "
                         + commandIsToFilterAndSummarizeResults + ", " + commandIsToAssessCoverage
                         + ", or " + commandAnnotation + ".\n" + usage);
    }
  }
}
