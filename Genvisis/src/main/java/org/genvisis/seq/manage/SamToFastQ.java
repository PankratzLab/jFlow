/**
 * 
 */
package org.genvisis.seq.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.common.ext;

/**
 *
 * 
 * Prepares bam input for re-genotyping by handling bam-> fastq conversion.
 *
 */
public class SamToFastQ {
	private SamToFastQ() {

	}

	/**
	 * @param inputBam bam to convert
	 * @param samToFastQLoc path to picard.jar, or path to SamToFastq.jar
	 * @param r1 output fastq for
	 * @param r2
	 * @param log
	 * @return
	 */
	private static boolean convertToFasta(String inputBam, String samToFastQLoc, String r1, String r2,
																				Logger log) {
		String[] inputs = new String[] {inputBam};
		String[] outputs = new String[] {r1, r2};
		ArrayList<String> command = new ArrayList<>();
		command.add("java");
		command.add("-jar");
		command.add(samToFastQLoc);
		if (samToFastQLoc.endsWith("picard.jar")) {
			command.add("SamToFastq");
		}
		command.add("I=" + inputBam);
		command.add("F=" + r1);
		command.add("F2=" + r2);

		return CmdLine.runCommandWithFileChecks(ArrayUtils.toStringArray(command), "", inputs, outputs,
																						true, false, false, log);
	}

	private static void prepBams(String bams, final String outDir, final String tag,
															 final String samToFastQ, int numThreads) {
		new File(outDir).mkdirs();
		final Logger log = new Logger(outDir + "samToFastq.log");
		final String[] bamsFiles = Files.isDirectory(bams) ? Files.listFullPaths(bams, ".bam")
																											 : HashVec.loadFileToStringArray(bams, false,
																																											 new int[] {0},
																																											 true);
		log.reportTimeInfo("Found " + bamsFiles.length + " bams from input " + bams);

		Producer<Boolean> prepProducer = new Producer<Boolean>() {
			private int index = 0;

			@Override
			public boolean hasNext() {

				return index < bamsFiles.length;
			}

			@Override
			public Callable<Boolean> next() {
				final String bamFile = bamsFiles[index];

				Callable<Boolean> callable = new Callable<Boolean>() {

					@Override
					public Boolean call() throws Exception {
						String sampleName = null;
						try {
							if (BamOps.getHeader(bamFile, log).getReadGroups().size() != 1) {
								throw new IllegalArgumentException("This method currently supports bam to .fastq conversion for 1 and only 1 PE readgroups");
							}
							sampleName = BamOps.getSampleName(bamFile, log);
							String rootOut = outDir + (tag == null ? "" : tag + "-") + sampleName
															 + "_UUUUU-UUUUU_L001_.fastq";
							String r1 = ext.addToRoot(rootOut, "R1_001") + ".gz";
							String r2 = ext.addToRoot(rootOut, "R2_001") + ".gz";
							boolean success = convertToFasta(bamFile, samToFastQ, r1, r2, log);

							if (!success) {
								log.reportError("Could not parse " + bamFile + ", removing any output");
								new File(r1).delete();
								new File(r2).delete();
							}
							return success;
						} catch (Exception e) {
							log.reportError("Could not process " + bamFile);
							log.reportException(e);
						}
						return false;
					}
				};
				index++;
				return callable;
			}

			@Override
			public void shutdown() {
				//
			}

			@Override
			public void remove() {
				//
			}
		};
		WorkerTrain<Boolean> train = new WorkerTrain<>(prepProducer, numThreads, 10, log);
		while (train.hasNext()) {
			train.next();
		}
	}



	public static void main(String[] args) {
		CLI c = new CLI(SamToFastQ.class);
		c.addArgWithDefault("bams",
												"file listing bam files to analyze, one per line - or a directory of bams",
												"bams.txt");
		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "out/");
		c.addArgWithDefault("samToFastq", "full path to SamToFastq.jar", "SamToFastq.jar");
		c.addArgWithDefault("tag", "custom ID tag to add to files", null);

		c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, CLI.EXAMPLE_THREADS);
		c.parseWithExit(args);

		prepBams(c.get("bams"), c.get(CLI.ARG_OUTDIR), c.get("tag"), c.get("samToFastq"),
						 c.getI(CLI.ARG_THREADS));

	}

}
