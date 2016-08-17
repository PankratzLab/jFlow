package org.genvisis.seq.analysis;



import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * Note we currently use bwa mem -M ref.fa read1.fq read2.fq > aln-pe.sam to do the allignments
 */
public class BWA {
	public static final String INDEX_COMMAND = "index";
	public static final String ALIGN_COMMAND = "aln";
	public static final String ALIGN_COMBINED_COMMAND = "sampe";
	public static final String BWA_MEM_COMMAND = "mem";
	public static final String BWA_COMMAND = "bwa";

  public static void main(String[] args) {
    String dir = "/home/pankrat2/shared/testGATK/fastq/";
    // String fq1 = dir + "F10004C_ATCACG_L006_R1_001.fastq";
    // String fq2 = dir + "F10004C_ATCACG_L006_R2_001.fastq";
    String fq1 = dir + "F10004M_ATCACG_L005_R1_001.fastq";
    String fq2 = dir + "F10004M_ATCACG_L005_R2_001.fastq";
    String ref = "/home/pankrat2/lanej/bin/ref/hg19_canonical.fa";
    test(fq1, fq2, ref);
  }

  public static void test(String fq1, String fq2, String ref) {
    BWA bwa = new BWA(null, true, true, new Logger("/home/pankrat2/shared/testGATK/testing.log"));
    // bwa.indexReferenceGenome(ref);
    String refIndex = ref;

    String finalSame = ext.rootOf(fq2, false) + ".test.sam";

    bwa.bwaMEM(refIndex, fq1, fq2, finalSame, "test", 8, null);

  }

  private String bwaLocation;

  private final boolean fail, verbose, overwriteExisting;

  private final Logger log;

  // TODO
  // public boolean indexReferenceGenome(String fullPathToReferenceFasta) {
  // boolean success = false;
  // String command = bwaLocation + " " + INDEX_COMMAND + " " + fullPathToReferenceFasta;
  // String[] expectedOutput = new String{fullPathToReferenceFasta+"."
  // }

	public BWA(String bwaLocation, boolean overwriteExisting, boolean verbose, Logger log) {
		super();
		this.bwaLocation = bwaLocation;
		this.log = log;
		this.overwriteExisting = overwriteExisting;
		this.verbose = verbose;
		fail = !validBWA();
	}


  public boolean bwaMEM(String fullPathToReferenceIndexedFasta, String fullPathToForwardReadFQ,
                        String fullPathToRevereseReadFQ, String fullPathToOutputFile,
                        String readGroup, int numThreads, Logger altLog) {
    String[] inputFiles = new String[] {fullPathToReferenceIndexedFasta, fullPathToForwardReadFQ,
                                        fullPathToRevereseReadFQ};
    String[] outputFiles = new String[] {fullPathToOutputFile};
    String[] commandArray = new String[] {bwaLocation, BWA_MEM_COMMAND, "-M", "-R", readGroup,
                                          numThreads > 1 ? " -t " + numThreads : "",
                                          fullPathToReferenceIndexedFasta, fullPathToForwardReadFQ,
                                          fullPathToRevereseReadFQ, ">", fullPathToOutputFile};
    if (!fail) {
      return CmdLine.runCommandWithFileChecks(commandArray, "", inputFiles, outputFiles,
                                              verbose, overwriteExisting, true,
                                              altLog == null ? log : altLog);
    } else {
      return false;
    }
  }

  public String getBwaLocation() {
    return bwaLocation;
  }

	public boolean isFail() {
		return fail;
	}

	private boolean validBWA() {

		if (bwaLocation != null && !bwaLocation.equals("") && !bwaLocation.equals(BWA_COMMAND)) {
			if (Files.exists(bwaLocation)) {
				if (verbose) {
					log.report(ext.getTime() + " Info - using bwa located at " + bwaLocation);
				}
				return true;
			} else {
				log.reportError("Error - invalid bwa location, could not find" + bwaLocation);
				return false;
			}
		} else {
			if (CmdLine.run(BWA_COMMAND, "")) {
				if (verbose) {
					log.report(ext.getTime() + " Info - using bwa set by the system's path varaible");
				}
				bwaLocation = BWA_COMMAND;
				return true;
			} else {
				log.reportError("Error - a path to bwa was not supplied and bwa was not detected on the system's path");
				return false;
			}
		}
	}
}
//
// public boolean combinePairedReads(String fullPathToReferenceIndexedFasta, String
// fullPathToForwardReadFQ, String fullPathToRevereseReadFG, String fullPathToForwardReadSAI, String
// fullPathToRevereseReadSAI, String fullPathToOutputFile) {
// boolean success = false;
// String[] inputFiles = new String[] { fullPathToReferenceIndexedFasta, fullPathToForwardReadFQ,
// fullPathToRevereseReadFG, fullPathToForwardReadSAI, fullPathToRevereseReadSAI };
// if (verbose) {
// log.report(ext.getTime() + " Info - attempting to allign " + Array.toStr(inputFiles, "\n"));
// }
// if (!Files.exists(fullPathToOutputFile) || overwriteExisting) {
// if (!fail) {
// if (checkAllExist(inputFiles, log)) {
// if (CmdLine.run(bwaLocation + " " + ALIGN_COMBINED_COMMAND + " " + fullPathToForwardReadSAI + " "
// + fullPathToRevereseReadSAI + " " + fullPathToForwardReadFQ + " " + fullPathToRevereseReadFG,
// "")) {
// if (!Files.exists(fullPathToOutputFile)) {
// log.reportError("Error - the seemed to run, but the output file " + fullPathToOutputFile + " was
// not found");
// } else {
// if (verbose) {
// log.report(ext.getTime() + " Info - finished alligning " + Array.toStr(inputFiles, "\n"));
// }
// success = true;
// }
// } else {
// log.reportError("Error - the command " + bwaLocation + " " + ALIGN_COMBINED_COMMAND + " " +
// fullPathToForwardReadSAI + " " + fullPathToRevereseReadSAI + " " + fullPathToForwardReadFQ + " "
// + fullPathToRevereseReadFG + " has failed");
// }
// } else {
// log.reportError("Error - could not find all neccesary input files");
// }
// } else {
// log.reportError("Error - cannot allign combined forward and reverse reads, failed on a previous
// step");
// }
// } else {
// log.report(ext.getTime() + " Warning - the output file " + fullPathToOutputFile + " already
// exists and the overwrite option was not flagged, skipping");
// success = true;
// }
// return success;
// }
//
// public boolean allignSeparateReads(String fullPathToReferenceIndexedFasta, String
// fullPathToForwardOrReverseReadFQ, String fullPathToOutputFile) {
// boolean success = false;
// String[] inputFiles = new String[] { fullPathToReferenceIndexedFasta,
// fullPathToForwardOrReverseReadFQ };
// if (verbose) {
// log.report(ext.getTime() + " Info - attempting to allign " + Array.toStr(inputFiles, "\n"));
// }
// if (!Files.exists(fullPathToOutputFile) || overwriteExisting) {
// if (!fail) {
// if (checkAllExist(inputFiles, log)) {
// if (CmdLine.run(bwaLocation + " " + ALIGN_COMMAND + " " + fullPathToReferenceIndexedFasta + " " +
// fullPathToForwardOrReverseReadFQ + " > " + fullPathToOutputFile, "")) {
// if (!Files.exists(fullPathToOutputFile)) {
// log.reportError("Error - the command ran successfully, but the output file " +
// fullPathToOutputFile + " was not found");
// } else {
// if (verbose) {
// log.report(ext.getTime() + " Info - finished alligning " + fullPathToReferenceIndexedFasta + "
// and " + fullPathToForwardOrReverseReadFQ);
// }
// success = true;// finally
// }
// } else {
// log.reportError("Error - the command " + bwaLocation + " " + ALIGN_COMMAND + " " +
// fullPathToReferenceIndexedFasta + " " + fullPathToForwardOrReverseReadFQ + " > " +
// fullPathToOutputFile + " has failed");
// }
// } else {
// log.reportError("Error - could not find all neccesary input files");
// }
// } else {
// log.reportError("Error - cannot allign " + fullPathToReferenceIndexedFasta + " and " +
// fullPathToForwardOrReverseReadFQ + ", failed on a previous step");
// }
// } else {
// log.report(ext.getTime() + " Warning - the output file " + fullPathToOutputFile + " already
// exists and the overwrite option was not flagged, skipping");
// success = true;
// }
// return success;
// }

// if (verbose) {
// log.report(ext.getTime() + " Info - attempting to allign \n" + Array.toStr(inputFiles, "\n"));
// }
// if (numThreads > 10) {
//
// }
// if (!Files.exists(fullPathToOutputFile) || overwriteExisting) {
// if (!fail) {
// if (checkAllExist(inputFiles, log)) {
// try {
// PrintStream printStream = new PrintStream(log.getFilename() + "." + BWA_MEM_COMMAND);
// if (verbose) {
// log.report(ext.getTime() + " Info - attempting to run command" + command);
// }
// if (CmdLine.run(command, "", printStream)) {
// if (!Files.exists(fullPathToOutputFile)) {
// log.reportError("Error - the command " + command + "\nseemed to run, but the output file " +
// fullPathToOutputFile + " was not found");
// } else {
// if (verbose) {
// log.report(ext.getTime() + " Info - finished alligning " + Array.toStr(inputFiles, "\n"));
// }
// success = true;
// }
// } else {
// log.reportError("Error - the command " + command + " has failed");
// }
// printStream.close();
// } catch (FileNotFoundException e) {
// log.reportError("Error - could not generate log stream for bwa mem");
// log.reportException(e);
// }
// } else {
// log.reportError("Error - could not find all neccesary input files");
// }
// } else {
// log.reportError("Error - cannot allign combined forward and reverse reads, failed on a previous
// step");
// }
// } else {
// log.report(ext.getTime() + " Warning - the output file " + fullPathToOutputFile + " already
// exists and the overwrite option was not flagged, skipping");
// success = true;
// }
// return success;
