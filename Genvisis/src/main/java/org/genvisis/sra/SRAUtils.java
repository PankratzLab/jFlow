package org.genvisis.sra;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;

/**
 * @author lane0212 Converts .sra files to .bam format for future processing.
 */
public class SRAUtils {
	private static final String SAM_DUMP = "sam-dump.2.6.3";
	private static final String SRA_EXT = ".sra";

	private SRAUtils() {

	}

	/**
	 * Some info about the conversion
	 *
	 */
	public static class SRAConversionResult {
		private String outputBam;
		private boolean valid;

		/**
		 * @param outputBam
		 *            the bam that was writing to
		 * @param valid
		 *            whether the conversion was successful
		 * @param log
		 */
		public SRAConversionResult(String outputBam, boolean valid) {
			super();
			this.outputBam = outputBam;
			this.valid = valid;
		}

		public String getOutputBam() {
			return outputBam;
		}

		public boolean isValid() {
			return valid;
		}

	}

	private static class SRABamWorker implements Callable<SRAConversionResult> {
		private String inputSra;
		private String outputBam;
		private Logger log;

		public SRABamWorker(String inputSra, String outputBam, Logger log) {
			super();
			this.inputSra = inputSra;
			this.outputBam = outputBam;
			this.log = log;
		}

		@Override
		public SRAConversionResult call() throws Exception {
			boolean valid = dumpSra(inputSra, outputBam, log);
			return new SRAConversionResult(outputBam, valid);
		}

		/**
		 * @param inputSra
		 *            full path
		 * @param outputBam
		 *            full path
		 * @param log
		 * @return
		 */
		private static boolean dumpSra(String inputSra, String outputBam, Logger log) {
			String[] inputs = new String[] { inputSra };
			String[] outputs = new String[] { outputBam };
			ArrayList<String> command = new ArrayList<String>();
			command.add("cd " + ext.parseDirectoryOfFile(inputSra) + "\n");
			command.add(SAM_DUMP);
			command.add("-u");// output un-mapped reads as well
			command.add(ext.rootOf(inputSra, true));
			command.add("|");
			command.add("samtools");
			command.add("view");
			command.add("-bS");// convert to bam
			command.add("-");// pipe input
			command.add(">");
			command.add(outputBam);

			String[] bat = CmdLine.prepareBatchForCommandLine(Array.toStringArray(command), outputBam + ".bat", true,
					log);
			return CmdLine.runCommandWithFileChecks(bat, "", inputs, outputs, true, false, false, log);
		}

	}

	private static class SRABamProducer extends AbstractProducer<SRAConversionResult> {
		private String[] inputSras;
		private String outDir;
		private Logger log;
		private int index;

		public SRABamProducer(String[] inputSras, String outDir, Logger log) {
			super();
			this.inputSras = inputSras;
			this.outDir = outDir;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < inputSras.length;
		}

		@Override
		public Callable<SRAConversionResult> next() {
			SRABamWorker worker = new SRABamWorker(inputSras[index], outDir + ext.rootOf(inputSras[index]) + ".bam",
					log);
			index++;
			return worker;
		}
	}

	// sam-dump.2.6.3 SRR1737697 |samtools view -bS -

	/**
	 * @param sraDir
	 *            your directory that S
	 * @param outDir
	 * @param threads
	 * @return
	 */
	public static List<SRAConversionResult> run(String sraDir, String outDir, int threads) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "sraConv.log");
		String[] sraFiles = Files.listFullPaths(sraDir, SRA_EXT, false);
		log.reportTimeInfo("Found " + sraFiles.length + " " + SRA_EXT + " files");
		String bamDir = outDir + "bams/";
		new File(bamDir).mkdirs();

		SRABamProducer producer = new SRABamProducer(sraFiles, bamDir, log);
		WorkerTrain<SRAConversionResult> train = new WorkerTrain<SRAUtils.SRAConversionResult>(producer, threads, 10,
				log);
		ArrayList<SRAConversionResult> results = new ArrayList<SRAUtils.SRAConversionResult>();
		while (train.hasNext()) {
			results.add(train.next());
		}
		return results;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String sraDir = "/scratch.global/lanej/aric_raw/sra/";
		String outDir = "/scratch.global/lanej/aric_raw/";
		int threads = 24;

		String usage = "\n" + " SRAUtils requires 0-1 arguments\n" + "   (1) SRA directory (i.e. sraDir=" + sraDir
				+ " (default))\n" + "   (2) out directory (i.e. outDir=" + outDir + " (default))\n"
				+ PSF.Ext.getNumThreadsCommand(3, threads) + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("sraDir")) {
				sraDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outDir")) {
				outDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				threads = ext.parseIntArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			run(sraDir, outDir, threads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
