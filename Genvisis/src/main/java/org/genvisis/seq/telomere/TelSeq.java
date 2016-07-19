package org.genvisis.seq.telomere;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.filesys.LocusSet;
import org.genvisis.seq.analysis.MitoSeqCN;
import org.genvisis.seq.manage.BEDFileReader;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.BEDFileReader.BEDFeatureSeg;
import org.genvisis.seq.telomere.SRAUtils.SRAConversionResult;

//http://goggable.areteh.co:3000/RotBlauer/IntallingTelSeq for telseq install instructions
public class TelSeq {

	private enum TYPE {
		BASE,
		BED,
		BUFFERED_BED;
	}

	private static class TelSeqResult {
		private String output;
		private Ran ran;
		private String sample;
		private int readSizeUsed;
		private TYPE type;

		public TelSeqResult(String output, Ran ran, String sample, TYPE type, int readSizeUsed) {
			super();
			this.output = output;
			this.ran = ran;
			this.sample = sample;
			this.type = type;
			this.readSizeUsed = readSizeUsed;
		}

	}

	private static class TelSeqWorker implements Callable<TelSeqResult> {
		private String inputBam;
		private ArrayList<String> additionalArgs;
		private String outputDir;
		private TYPE type;

		private Logger log;

		public TelSeqWorker(String inputBam, ArrayList<String> additionalArgs, String outputDir, TYPE type, Logger log) {
			super();
			this.inputBam = inputBam;
			this.additionalArgs = additionalArgs;
			this.outputDir = outputDir;
			this.type = type;
			this.log = log;
		}

		@Override
		public TelSeqResult call() throws Exception {
			String out = outputDir + ext.rootOf(inputBam) + ".telseq";
			int readSize = BamOps.estimateReadSize(inputBam, 20000, log);
			log.reportTimeInfo("Estimated readsize for " + inputBam + " to be " + readSize);
			Ran ran = telSeqIt(inputBam, out, readSize, additionalArgs, log);
			String sampleName = "NA";
			try {
				sampleName = BamOps.getSampleName(inputBam);
			} catch (Exception e) {
				log.reportTimeError("Could not get sample name from " + inputBam);
			}

			return new TelSeqResult(out, ran, sampleName, type, readSize);
		}
	}

	private static class TelSeqProducer implements Producer<TelSeqResult> {
		private String[] inputBams;
		private ArrayList<String> additionalArgs;
		private String outputDir;
		private TYPE type;
		private Logger log;
		private int index;

		public TelSeqProducer(String[] inputBams, ArrayList<String> additionalArgs, String outputDir, TYPE type, Logger log) {
			super();
			this.inputBams = inputBams;
			this.additionalArgs = additionalArgs;
			this.outputDir = outputDir;
			this.type = type;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < inputBams.length;
		}

		@Override
		public Callable<TelSeqResult> next() {
			TelSeqWorker worker = new TelSeqWorker(inputBams[index], additionalArgs, outputDir, type, log);
			index++;
			return worker;
		}

		@Override
		public void shutdown() {

		}

	}

	private static class Ran {
		private boolean valid;
		private ArrayList<String> command;

		public Ran(boolean valid, ArrayList<String> command) {
			super();
			this.valid = valid;
			this.command = command;
		}

	}

	private static Ran telSeqIt(String inputBam, String output, int readSize, ArrayList<String> additionalArgs, Logger log) {
		String[] outputs = new String[] { output };
		String[] input = new String[] { inputBam };
		ArrayList<String> command = new ArrayList<String>();
		command.add("telseq");
		command.add("-o");
		command.add(output);
		command.add("-r");
		command.add(readSize + "");

		command.addAll(additionalArgs);

		command.add(inputBam);

		boolean valid = CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs, true, false,
				false, log);
		return new Ran(valid, command);
	}

	/**
	 * Likely tmp method for dumping SRA format first
	 */
	public static void runTelSeqSRA(String sraDir, String outDir, String optionalBed, String referenceGenomeFasta,
			int threads) {
		Logger log = new Logger(outDir + ".telseq.log");

		ArrayList<SRAConversionResult> conv = SRAUtils.run(sraDir, outDir, threads);

		ArrayList<String> bamst = new ArrayList<String>();
		for (int i = 0; i < conv.size(); i++) {
			if (conv.get(i).isValid() && new File(conv.get(i).getOutputBam()).getTotalSpace() > 0) {
				bamst.add(conv.get(i).getOutputBam());

			} else {
				log.reportTimeWarning(conv.get(i).getOutputBam() + " must have failed, skipping");
			}
		}
		log.reportTimeInfo("Found " + bamst.size() + " bams to analyze");

		runTelSeq(Array.toStringArray(bamst), outDir, optionalBed, referenceGenomeFasta, threads, false, false, log);

		// log.reportTimeError("John remember to add back in SRA");
//		String[] bams = Files.listFullPaths(outDir + "bams/", ".bam", false);
//		log.reportTimeInfo("Found " + bams.length + " bams");
//		runTelSeq(bams, outDir, optionalBed, referenceGenomeFasta, threads, false, false, log);

	}

	/**
	 * @param bams
	 *            bams to run
	 * @param outDir
	 *            where to put things
	 * @param optionalBed
	 *            for WES
	 * @param referenceGenomeFasta
	 *            will likely be removed later
	 * @param threads
	 * @param log
	 * 
	 */
	private static void runTelSeq(String[] bams, String outDir, String optionalBed, String referenceGenomeFasta,
			int threads, boolean onlyExome, boolean chr, Logger log) {
		if (log == null) {
			log = new Logger(outDir + ".telseq.log");

		}
		String telseqDir = outDir + "telseq/";
		new File(telseqDir).mkdirs();
		log.reportTimeInfo("Assuming telseq is on system path");
		ArrayList<TelSeqResult> results = new ArrayList<TelSeq.TelSeqResult>();
		ArrayList<String> argPopulator = new ArrayList<String>();
		argPopulator.add("-m");// doesn't look like telseq handles RGs properly
		String baseDir = telseqDir + "base/";
		new File(baseDir).mkdirs();
		// TODO, do either or with optional bed, currently testing
		if (!onlyExome) {
			runType(threads, log, bams, results, argPopulator, baseDir, TYPE.BASE);
		}
		if (optionalBed != null) {
			if (Files.exists(optionalBed)) {
				BEDFileReader reader = new BEDFileReader(optionalBed, false);

				LocusSet<BEDFeatureSeg> segs = reader.loadAll(log);
				reader.close();

				String dirBed = telseqDir + ext.rootOf(optionalBed) + "/";
				String baseBed = dirBed + "base_" + ext.rootOf(optionalBed) + ".bed";
				ArrayList<String> argPopulatorBed = new ArrayList<String>();
				argPopulatorBed.addAll(argPopulator);
				argPopulatorBed.add("-e");
				argPopulatorBed.add(baseBed);
				log.reportTimeInfo("writing bed to " + baseBed);
				new File(dirBed).mkdirs();

				segs.writeSegmentRegions(baseBed, !chr, log);

				runType(threads, log, bams, results, argPopulatorBed, dirBed, TYPE.BED);
				if (!onlyExome) {

					String buffDir = telseqDir + "buff_20KB_" + ext.rootOf(optionalBed) + "/";
					new File(buffDir).mkdirs();
					String buffBed = buffDir + "buff_20KB_" + ext.rootOf(optionalBed) + ".bed";
					log.reportTimeInfo("writing bed to " + buffBed);
					segs.getBufferedSegmentSet(20000).writeSegmentRegions(buffBed, !chr, log);

					ArrayList<String> argPopulatorBuffBed = new ArrayList<String>();
					argPopulatorBuffBed.addAll(argPopulator);
					argPopulatorBuffBed.add("-e");
					argPopulatorBuffBed.add(buffBed);

					runType(threads, log, bams, results, argPopulatorBuffBed, buffDir, TYPE.BUFFERED_BED);
				}

			} else {
				log.reportFileNotFound(optionalBed);
			}
		}

		// summarize
		String finalOut = telseqDir + "telseq.summary.txt";
		String[] telHeader = Files.getHeaderOfFile(results.get(0).output, log);

		ArrayList<String> result = new ArrayList<String>();
		result.add("BAM\t" + Array.toStr(telHeader) + "\tType\tSampleName\tReadSize");
		for (TelSeqResult telSeqResult : results) {
			if (Files.exists(telSeqResult.output)) {
				String[][] data = HashVec.loadFileToStringMatrix(telSeqResult.output, true, null, false);
				for (int i = 0; i < data.length; i++) {
					result.add(ext.rootOf(telSeqResult.output) + "\t" + Array.toStr(data[i]) + "\t"
							+ telSeqResult.type + "\t" + telSeqResult.sample + "\t" + telSeqResult.readSizeUsed);
				}
			}
		}
		Files.writeArrayList(result, finalOut);

		// can kill this later... going to do mtDNA CN

		if (optionalBed != null) {

			String mitoDir = outDir + "mitoCN/";
			new File(mitoDir).mkdirs();
			String baseBed = mitoDir + "base_" + ext.rootOf(optionalBed) + ".bed";

			BEDFileReader reader = new BEDFileReader(optionalBed, false);
			LocusSet<BEDFeatureSeg> segs = reader.loadAll(log);
			reader.close();
			segs.writeSegmentRegions(baseBed, true, log);

			String bamsToMito = mitoDir + "bams.txt";
			Files.writeList(bams, bamsToMito);
			MitoSeqCN.run(bamsToMito, mitoDir, baseBed, referenceGenomeFasta, chr, threads);
		}
	}

	private static void runType(int threads, Logger log, String[] bams, ArrayList<TelSeqResult> results,
			ArrayList<String> argPopulator, String baseDir, TYPE type) {
		TelSeqProducer producer = new TelSeqProducer(bams, argPopulator, baseDir, type, log);
		WorkerTrain<TelSeqResult> train = new WorkerTrain<TelSeq.TelSeqResult>(producer, threads, 100, log);
		while (train.hasNext()) {
			results.add(train.next());
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String sraDir = "/scratch.global/lanej/aric_raw/sra/";
		String outDir = "/scratch.global/lanej/aric_raw/";
		String captureBed = "/home/pankrat2/public/bin/ref/VCRome_2_1_hg19_capture_targets.bed";
		String refGenome = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
		String bamList = null;
		boolean onlyExome = false;
		boolean chr = false;
		// String captureBed = null;

		int threads = 24;

		String usage = "\n" + "telomere.SRAUtils requires 0-1 arguments\n" + "   (1) SRA directory (i.e. sraDir="
				+ sraDir + " (default))\n" + "   (2) out directory (i.e. outDir=" + outDir + " (default))\n"
				+ "   (3) capture bed (i.e. bed=" + outDir + " (default))\n" +

				PSF.Ext.getNumThreadsCommand(3, threads) + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("sraDir")) {
				sraDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bed=")) {
				captureBed = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outDir=")) {
				outDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bams=")) {
				bamList = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-exome")) {
				onlyExome = true;
				numArgs--;
			} else if (args[i].startsWith("-chr")) {
				// chr = true;
				System.out.println("-chr doesn't work, ignoring");
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
			if (bamList != null) {
				runTelSeq(HashVec.loadFileToStringArray(bamList, false, new int[] { 0 }, true), outDir, captureBed, refGenome, threads, onlyExome, chr, null);
			} else {
				runTelSeqSRA(sraDir, outDir, captureBed, refGenome, threads);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
