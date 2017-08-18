package org.genvisis.one.kir;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.common.ext;
import org.genvisis.seq.telomere.Computel;

/**
 *
 *
 * Testing run of Ping and setup
 *
 * https://github.com/Hollenbach-lab/PING http://hollenbachlab.ucsf.edu/ping
 *
 *
 */
public class Ping {

	private static final String CWD = "CURRENT_WD";
	private static final String SAMP_LOC = "SAMPLE_LOC";
	private static final String OUT_DIR = "OUT_DIR";
	private static final String BOWTIE_THREADS = "MANY_THREADS";

	private static final String PING_SCRIPT = "source(\"" + CWD + "PING_run_all.R\")\r\n" + "setwd(\"" + CWD + "\")\r\n"
			+ "  sample.location = \"" + SAMP_LOC + "\"\r\n " + "fastq.pattern.1 = \"r1.fastq\"\r\n"
			+ "  fastq.pattern.2 = \"r2.fastq\"\r\n " + "bowtie.threads = " + BOWTIE_THREADS + "\r\n "
			+ " results.directory = \"" + OUT_DIR + "\"\r\n"
			+ "ping_run_all(sample.location = sample.location,fastq.pattern.1 = fastq.pattern.1,fastq.pattern.2 = fastq.pattern.2,bowtie.threads = bowtie.threads,results.directory = results.directory)\r\n";

	private static final class ConvProducer implements Producer<Boolean> {
		private final String[] bams;
		private final String outDir;
		private final String samToFastQLoc;
		private final Logger log;
		int index = 0;

		private ConvProducer(String[] bams, String outDir, String samToFastQLoc, Logger log) {
			this.bams = bams;
			this.outDir = outDir;
			this.samToFastQLoc = samToFastQLoc;
			this.log = log;
		}

		@Override
		public boolean hasNext() {
			return index < bams.length;
		}

		@Override
		public Callable<Boolean> next() {
			final String bam = bams[index];
			final String r1 = outDir + ext.rootOf(bam, true) + "_r1.fastq";
			final String r2 = outDir + ext.rootOf(bam, true) + "_r2.fastq";
			index++;
			return () -> convertToFasta(bam, samToFastQLoc, r1, r2, log);

		}

		@Override
		public void shutdown() {
			//

		}

		private static boolean convertToFasta(String inputBam, String samToFastQLoc, String r1, String r2, Logger log) {
			String[] inputs = new String[] { inputBam };
			String[] outputs = new String[] { r1, r2 };
			ArrayList<String> command = new ArrayList<>();
			command.add("java");
			command.add("-jar");
			command.add(samToFastQLoc);
			command.add("I=");
			command.add(inputBam);
			command.add("F=");
			command.add(r1);
			command.add("F2=");
			command.add(r2);

			return CmdLine.runCommandWithFileChecks(ArrayUtils.toStringArray(command), "", inputs, outputs, true, false,
					false, log);
		}
	}

	private Ping() {

	}

	private static void convertBams(String bamDir, String samToFastQLoc, String outDir, int numthreads, Logger log) {
		new File(outDir).mkdirs();
		final String[] bams = Files.listFullPaths(bamDir, ".bam");
		Producer<Boolean> convProducer = new ConvProducer(bams, outDir, samToFastQLoc, log);
		WorkerTrain<Boolean> train = new WorkerTrain<>(convProducer, numthreads, 10, log);
		int index = 0;
		log.reportTimeInfo("found " + bams.length + " bams to convert");
		while (train.hasNext()) {
			if (!train.next()) {
				log.reportError("Could not convert " + bams[index] + " to fastq");
			}
			log.reportTimeInfo("finished converting " + bams[index]);
			index++;
		}
	}

	private static void runPING(String bamDir, String samToFastQLoc, String pingLoc, String outDir, int numthreads) {
		Logger log = new Logger(outDir + "Ping.log");
		new File(outDir).mkdirs();
		convertBams(bamDir, samToFastQLoc, getFastQLoc(outDir), numthreads, log);
		try {
			File outPing = new File(getPingLoc(outDir));
			Files.copyRecursive(new File(pingLoc), outPing);
			String pingScript = PING_SCRIPT.replaceAll(CWD, getPingLoc(outDir)).replaceAll(SAMP_LOC, getFastQLoc(outDir))
					.replaceAll(OUT_DIR, outDir + "results/").replaceAll(BOWTIE_THREADS, Integer.toString(numthreads));

			String pingR = outDir + ext.getTimestampForFilename() + "ping.Rscript";
			CmdLine.prepareBatchForCommandLine(new String[] { pingScript }, pingR, true, log);
			CmdLine.runCommandWithFileChecks(new String[] { "Rscript", pingR }, "", new String[] { pingR }, null, true,
					false, false, log);

		} catch (IOException e) {
			log.reportError("could not copy Ping to correct location");
			log.reportException(e);
		}

	}

	private static String getPingLoc(String outDir) {
		return outDir + "ping/";
	}

	private static String getFastQLoc(String outDir) {
		return outDir + "fastq/";
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		CLI c = new CLI(Computel.class);
		c.addArgWithDefault("bams", "directory of bams to analyze", "bams/");

		final String ping = "ping";
		c.addArgWithDefault(ping, "full PING directory (as git clone ideally)", ping);

		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "out/");
		c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, "24");
		c.addArgWithDefault("samToFastQ", "samToFastQ jar location", "samToFastQ/");

		c.parseWithExit(args);

		runPING(c.get("bams"), c.get("samToFastQ"), c.get("ping"), c.get(CLI.ARG_OUTDIR), c.getI(CLI.ARG_THREADS));

	}

}
