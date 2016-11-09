/**
 * 
 */
package org.genvisis.seq.analysis.mtdna;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.CLI;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.BamOps;

/**
 * @author Kitty
 * 
 *         Prepares bam input for genotyping
 *
 */
public class MtDNAGenotypePrep {
	private MtDNAGenotypePrep() {

	}

	private static boolean convertToFasta(String inputBam, String samToFastQLoc, String r1, String r2,
																				Logger log) {
		String[] inputs = new String[] {inputBam};
		String[] outputs = new String[] {r1, r2};
		ArrayList<String> command = new ArrayList<String>();
		command.add("java");
		command.add("-jar");
		command.add(samToFastQLoc);
		command.add("I=");
		command.add(inputBam);
		command.add("F=");
		command.add(r1);
		command.add("F2=");
		command.add(r2);

		return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, outputs, true,
																						false, false, log);
	}

	private static void prepBams(	String bams, final String outDir, final String tag,
																final String samToFastQ, int numThreads) {
		new File(outDir).mkdirs();
		final Logger log = new Logger(outDir + "mtdnaPrep.log");
		final String[] bamsFiles = HashVec.loadFileToStringArray(bams, false, new int[] {0}, true);
		log.reportTimeInfo("Found " + bamsFiles.length + " bams");

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
							sampleName = BamOps.getSampleName(bamFile, log);
							String rootOut = outDir+ (tag == null ? "" : tag + "-") + sampleName
																+ "_UUUUU-UUUUU_L001_.fastq";
							String r1 = ext.addToRoot(rootOut, "R1_001");
							String r2 = ext.addToRoot(rootOut, "R2_001");
							boolean success = convertToFasta(bamFile, samToFastQ, r1, r2, log);

							if (!success) {
								log.reportTimeError("Could not parse " + bamFile + ", removing any output");
								new File(r1).delete();
								new File(r2).delete();
							}
							return success;
						} catch (Exception e) {
							log.reportTimeError("Could not process " + bamFile);
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
		};
		WorkerTrain<Boolean> train = new WorkerTrain<Boolean>(prepProducer, numThreads, 10, log);
		while (train.hasNext()) {
			train.next();
		}
	}

	public static void main(String[] args) {
		CLI c = new CLI(MtDNAGenotypePrep.class);
		c.addArgWithDefault("bams", "file listing bam files to analyze, one per line", "bams.txt");
		c.addArgWithDefault("outDir", "output directory", "out/");
		c.addArgWithDefault("samToFastq", "full path to SamToFastq.jar", "SamToFastq.jar");
		c.addArgWithDefault("tag", "custom ID tag to add to files", null);

		c.addArgWithDefault("threads", "number of threads", "24");
		c.parseWithExit(args);

		prepBams(c.get("bams"), c.get("outDir"), c.get("tag"), c.get("samToFastq"), c.getI("threads"));

	}

}
