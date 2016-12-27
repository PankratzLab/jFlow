package org.genvisis.seq.qc.contamination;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.Blast.BlastResultsSummary;

/**
 * @author lane0212 Class for assessing taxonomic contamination from .fastq files
 */
public class BlastContamination {
	private final WorkerTrain<BlastResultsSummary[]> train;
	private final Hashtable<String, Integer> taxonCounts;

	public BlastContamination(String currentDB, WorkerTrain<BlastResultsSummary[]> train,
	                          Logger log) {
		super();
		this.train = train;
		taxonCounts = new Hashtable<String, Integer>();
	}

	public void runContam() {
		while (train.hasNext()) {
			BlastResultsSummary[] summaries = train.next();
			for (BlastResultsSummary summarie : summaries) {
				Set<String> hits = summarie.getHitCounts().keySet();
				for (String hit : hits) {
					if (!taxonCounts.containsKey(hit)) {
						taxonCounts.put(hit, 0);
					}
					int cur = taxonCounts.get(hit);
					taxonCounts.put(hit, cur + summarie.getHitCounts().get(hit));
				}
			}
		}
	}

	public Hashtable<String, Integer> getTaxonCounts() {
		return taxonCounts;
	}

	// public static void batch(String)

	public static void setup(String dbDir, String fastaqDir, String fileOfFastqs, int numReads,
	                         int numReadsPerThread, int blastWordSize, int reportWordSize,
	                         int numSampThreads, String outputDir, String outputRoot, int numBatches,
	                         int memoryInMB, int wallTimeInHours, Logger log) {
		String[] fastaDbs = Files.toFullPaths(Files.list(dbDir, "nt.", "nsq", true, false), dbDir);
		String[] fastaqs = null;
		if (fastaDbs.length == 0) {
			log.reportError("Did not find any \"nt\" database files listed in " + dbDir);
			return;
		} else {
			log.reportTimeInfo("Found  " + fastaDbs.length + " database files to search");
		}
		if (fastaqDir == null && fileOfFastqs == null) {
			log.reportError("A file or directory with .fastq files must be provided");
			return;
		}
		if (fastaqDir != null && !fastaqDir.equals("")) {
			fastaqs = Files.listFullPaths(fastaqDir, ".fastq", false);
			if (fastaqs.length == 0) {
				log.reportError("Did not find any .fastq files listed in " + fastaqDir);
				return;
			} else {
				log.reportTimeInfo("Found  " + fastaqs.length + " database files to search");
			}
		} else {
			if (Files.exists(fileOfFastqs)) {
				log.reportTimeInfo("Attempting to load fastq files listed in " + fileOfFastqs);
				fastaqs = HashVec.loadFileToStringArray(fileOfFastqs, false, new int[] {0}, true);
				if (fastaqs.length == 0) {
					log.reportError("Did not find any .fastq files listed in " + fileOfFastqs);
					return;
				} else {
					log.reportTimeInfo("Found  " + fastaqs.length + " .fastq files");
				}
			} else {
				log.reportFileNotFound(fileOfFastqs);
			}
		}
		new File(outputDir).mkdirs();

		if (numBatches > 1) {
			batchIt(dbDir, numReads, numReadsPerThread, blastWordSize, reportWordSize, numSampThreads,
			        outputDir, outputRoot, numBatches, memoryInMB, wallTimeInHours, log, fastaqs);

		} else {
			String outputFile = outputDir + outputRoot;
			Logger newLog = new Logger(ext.rootOf(outputFile, false) + ".log");
			runContams(fastaDbs, numReads, fastaqs, numSampThreads, numReadsPerThread, blastWordSize,
			           reportWordSize, outputFile, newLog);
		}

	}

	private static void runContams(String[] fastaDbs, int numReads, String[] fastaqs,
	                               int numSampThreads, int numReadsPerThread, int blastWordSize,
	                               int reportWordSize, String outputFile, Logger log) {
		ArrayList<Hashtable<String, Integer>> popCounts = new ArrayList<Hashtable<String, Integer>>();
		for (int i = 0; i < fastaqs.length; i++) {
			Hashtable<String, Integer> allCounts = new Hashtable<String, Integer>();
			for (int j = 0; j < fastaDbs.length; j++) {
				log.reportTimeInfo("Currently processing file " + ext.removeDirectoryInfo(fastaqs[i]) + " ("
				                   + i + " of " + fastaqs.length + ") .fastqs; (" + j + " of "
				                   + fastaDbs.length + ") dbs");
				Hashtable<String, Integer> curCounts = runContam(ext.rootOf(fastaDbs[j], false), numReads,
				                                                 fastaqs[i], numSampThreads,
				                                                 numReadsPerThread, blastWordSize,
				                                                 reportWordSize, log).getTaxonCounts();
				allCounts = ext.addHashCounts(allCounts, curCounts);
			}
			popCounts.add(allCounts);
		}
		ArrayList<String> allTaxa = new ArrayList<String>();
		HashSet<String> tmpUniq = new HashSet<String>();
		for (int i = 0; i < popCounts.size(); i++) {
			tmpUniq.addAll(popCounts.get(i).keySet());
		}
		allTaxa.addAll(tmpUniq);

		String[] allTaxaA = allTaxa.toArray(new String[allTaxa.size()]);
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
			writer.println("Taxa\t" + Array.toStr(fastaqs));
			for (String element : allTaxaA) {
				writer.print(element);
				for (int j = 0; j < fastaqs.length; j++) {
					if (popCounts.get(j).containsKey(element)) {
						writer.print("\t" + popCounts.get(j).get(element));
					} else {
						writer.print("\t" + 0);
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outputFile);
			log.reportException(e);
		}

	}

	private static BlastContamination runContam(String fastaDb, int numReads, String fastaq,
	                                            int numSampThreads, int numReadsPerThread,
	                                            int blastWordSize, int reportWordSize, Logger log) {
		BlastFastaq bFastaq = new BlastFastaq(fastaq, numReads, fastaDb, blastWordSize, reportWordSize,
		                                      numReadsPerThread, log);
		WorkerTrain<BlastResultsSummary[]> train = new WorkerTrain<BlastResultsSummary[]>(bFastaq,
		                                                                                  numSampThreads,
		                                                                                  10, log);
		BlastContamination blastContamination = new BlastContamination(fastaDb, train, new Logger());
		blastContamination.runContam();
		return blastContamination;
	}

	private static void batchIt(String dbDir, int numReads, int numReadsPerThread, int blastWordSize,
	                            int reportWordSize, int numSampThreads, String outputDir,
	                            String outputRoot, int numBatches, int memoryInMB,
	                            int wallTimeInHours, Logger log, String[] fastaqs) {
		List<String[]> splits = Array.splitUpArray(fastaqs, numBatches, log);
		String command = Array.toStr(PSF.Load.getAllModules(), "\n");
		String[][] batches = new String[splits.size()][1];
		for (int i = 0; i < batches.length; i++) {
			batches[i][0] = "batch_" + i + "_" + outputRoot;
			Files.writeArray(splits.get(i), outputDir + batches[i][0] + ".txt");
		}
		command += "\njava -Xmx" + memoryInMB + "m -jar " + PSF.Java.GENVISIS
		           + " seq.qc.contamination.BlastContamination" + PSF.Ext.SPACE;
		command += FASTQ_LIST_FILE + outputDir + "[%0].txt" + PSF.Ext.SPACE;
		command += DB_DIR + dbDir + PSF.Ext.SPACE;
		command += NUM_READS + numReads + PSF.Ext.SPACE;
		command += NUM_READS_THREAD + numReadsPerThread + PSF.Ext.SPACE;
		command += BLAST_WORD_SIZE + blastWordSize + PSF.Ext.SPACE;
		command += REPORT_WORD_SIZE + reportWordSize + PSF.Ext.SPACE;
		command += OUTPUT_DIR + outputDir + PSF.Ext.SPACE;
		command += OUTPUT_ROOT + "[%0].results.txt" + PSF.Ext.SPACE;
		command += PSF.Ext.NUM_THREADS_COMMAND + numSampThreads + PSF.Ext.SPACE;
		command += FASTQ_DIR + PSF.Ext.SPACE;

		// command += PSF.Ext.MEMORY_MB + memoryInMB + PSF.Ext.SPACE;
		// command += PSF.Ext.WALLTIME_HRS + wallTimeInHours + PSF.Ext.SPACE;
		Files.qsub("Contam" + outputRoot, command, batches, memoryInMB, wallTimeInHours,
		           numSampThreads);
	}

	private static final String FASTQ_LIST_FILE = "fastqFile=";
	private static final String FASTQ_DIR = "fastqDir=";

	private static final String DB_DIR = "dbDir=";
	private static final String NUM_READS = "numReads=";
	private static final String NUM_READS_THREAD = "numReadsPerThread=";
	private static final String BLAST_WORD_SIZE = "blastWordSize=";
	private static final String REPORT_WORD_SIZE = "reportWordSize=";
	private static final String OUTPUT_DIR = "outDir=";
	private static final String OUTPUT_ROOT = "outRoot=";
	private static final String NUM_BATCHES = "numBatches=";

	public static void main(String[] args) {
		int numArgs = args.length;
		String dbDir = "/home/ntDB/";
		String fastaqDir = null;
		int numReads = 100000;
		int blastWordSize = 100;
		int reportWordSize = 100;
		int numReadsPerThread = 10000;
		int numSampThreads = 10;
		int wallTimeInHours = 16;
		int memoryInMB = 22000;
		String fileOfFastqs = null;
		String outputDir = "contamination/";
		String outputRoot = "contam";
		int numBatches = 1;
		String usage = "\n" + "seq.qc.contamination.BlastContamination requires 0-1 arguments\n";
		usage += "   (1) full path to directory containing nt* database files (i.e. " + DB_DIR + dbDir
		         + " (default))\n" + "";
		usage += "   (2) full path to directory containing *.fastq files (i.e. " + FASTQ_DIR
		         + " (no default))\n" + "";
		usage += "   (3) number of reads to blast from each *.fastq file (i.e. " + NUM_READS + numReads
		         + " (default))\n" + "";
		usage += "   (4) number of reads to blast on each thread (i.e. " + NUM_READS_THREAD
		         + numReadsPerThread + " (default))\n" + "";
		usage += "   (5) blast word size to intialize search (i.e. " + BLAST_WORD_SIZE + blastWordSize
		         + " (default))\n" + "";
		usage += "   (6) blast word size to report (i.e. " + REPORT_WORD_SIZE + reportWordSize
		         + " (default))\n" + "";
		usage += "   (7) full path to output directory (i.e. " + OUTPUT_DIR + outputDir
		         + " (default))\n" + "";
		usage += "   (8) root output (relative to output directory) (i.e. " + OUTPUT_ROOT + outputRoot
		         + " (default))\n" + "";
		usage += "   (9) number of batches to break up the job into (i.e. " + NUM_BATCHES + numBatches
		         + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(10, numSampThreads);
		usage += PSF.Ext.getWallTimeCommand(11, wallTimeInHours);
		usage += PSF.Ext.getMemoryMbCommand(12, memoryInMB);
		usage += "   (13) full path to a file listing full paths of fastqs to use (i.e. "
		         + FASTQ_LIST_FILE + fileOfFastqs + " (default))\n" + "";
		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith(DB_DIR)) {
				dbDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(FASTQ_DIR)) {
				fastaqDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(OUTPUT_ROOT)) {
				outputRoot = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(FASTQ_LIST_FILE)) {
				fileOfFastqs = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(OUTPUT_DIR)) {
				outputDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(NUM_READS)) {
				numReads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(NUM_READS)) {
				numReads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(NUM_READS_THREAD)) {
				numReadsPerThread = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(BLAST_WORD_SIZE)) {
				blastWordSize = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(REPORT_WORD_SIZE)) {
				reportWordSize = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(NUM_BATCHES)) {
				numBatches = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numSampThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.MEMORY_MB)) {
				memoryInMB = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.WALLTIME_HRS)) {
				wallTimeInHours = ext.parseIntArg(arg);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			setup(dbDir, fastaqDir, fileOfFastqs, numReads, numReadsPerThread, blastWordSize,
			      reportWordSize, numSampThreads, outputDir, outputRoot, numBatches, memoryInMB,
			      wallTimeInHours, new Logger());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
// public static void test() {
// String dbDir = "/panfs/roc/rissdb/blast/current/";
// Logger log = new Logger();
// String[] fastaDbs = Files.toFullPaths(Files.list(dbDir, "nt.", "nsq", true, false), dbDir);
// String fastaqDir =
// "/panfs/roc/data_release/2/bmgc/tsaim/hiseq/141126_SN261_0580_AC5R7WACXX/Project_Tsai_Project_021/";
// String[] fastaqs = Files.listFullPaths(fastaqDir, ".fastq", false);
// String outputDir = "/home/tsaim/shared/Project_Tsai_Project_021/contamination/";
// log.reportTimeInfo("Found " + fastaqs.length + " fastq files in " + fastaqDir);
// ArrayList<Hashtable<String, Integer>> popCounts = new ArrayList<Hashtable<String, Integer>>();
// // fastaqs = Array.subArray(fastaqs, 0, 1);
// for (int i = 0; i < fastaqs.length; i++) {
// Hashtable<String, Integer> allCounts = new Hashtable<String, Integer>();
// for (int j = 0; j < fastaDbs.length; j++) {
// log.reportTimeInfo(i + " of " + fastaqs.length + ".fastqs, " + j + " of " + fastaDbs.length + "
// dbs");
// Hashtable<String, Integer> curCounts = runContam(ext.rootOf(fastaDbs[j], false), fastaqs[i], 10,
// log).getTaxonCounts();
// allCounts = ext.addHashCounts(allCounts, curCounts);
// }
// popCounts.add(allCounts);
// }
//
// ArrayList<String> allTaxa = new ArrayList<String>();
// for (int i = 0; i < popCounts.size(); i++) {
// allTaxa.addAll(popCounts.get(i).keySet());
// }
// String output = outputDir + "taxonCounts.txt";
// String[] allTaxaA = allTaxa.toArray(new String[allTaxa.size()]);
// try {
// PrintWriter writer = new PrintWriter(new FileWriter(output));
// writer.println("Taxa\t" + Array.toStr(fastaqs));
// for (int i = 0; i < allTaxaA.length; i++) {
// writer.print(allTaxaA[i]);
// for (int j = 0; j < fastaqs.length; j++) {
// if (popCounts.get(j).containsKey(allTaxaA[i])) {
// writer.print("\t" + popCounts.get(j).get(allTaxaA[i]));
// } else {
// writer.print("\t" + 0);
// }
// }
// writer.println();
// }
// writer.close();
// } catch (Exception e) {
// log.reportError("Error writing to " + output);
// log.reportException(e);
// }
// log.reportTimeInfo("Found " + fastaDbs.length + " fasta databases to search");
//
// }

