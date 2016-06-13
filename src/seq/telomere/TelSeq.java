package seq.telomere;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import seq.manage.BEDFileReader;
import seq.manage.BEDFileReader.BEDFeatureSeg;
import seq.telomere.SRAUtils.SRAConversionResult;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.PSF;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;
import filesys.LocusSet;

//http://goggable.areteh.co:3000/RotBlauer/IntallingTelSeq for telseq install instructions
public class TelSeq {

	private static class TelSeqResult {
		private String output;
		private Ran ran;

		public TelSeqResult(String output, Ran ran) {
			super();
			this.output = output;
			this.ran = ran;
		}

	}

	private static class TelSeqWorker implements Callable<TelSeqResult> {
		private String inputBam;
		private ArrayList<String> additionalArgs;
		private String outputDir;
		private Logger log;

		public TelSeqWorker(String inputBam, ArrayList<String> additionalArgs, String outputDir, Logger log) {
			super();
			this.inputBam = inputBam;
			this.additionalArgs = additionalArgs;
			this.outputDir = outputDir;
			this.log = log;
		}

		@Override
		public TelSeqResult call() throws Exception {
			String out = outputDir + ext.rootOf(inputBam) + ".telseq";
			Ran ran = telSeqIt(inputBam, out, additionalArgs, log);

			return new TelSeqResult(out, ran);
		}
	}

	private static class TelSeqProducer implements Producer<TelSeqResult> {
		private String[] inputBams;
		private ArrayList<String> additionalArgs;
		private String outputDir;
		private Logger log;
		private int index;

		public TelSeqProducer(String[] inputBams, ArrayList<String> additionalArgs, String outputDir, Logger log) {
			super();
			this.inputBams = inputBams;
			this.additionalArgs = additionalArgs;
			this.outputDir = outputDir;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < inputBams.length;
		}

		@Override
		public Callable<TelSeqResult> next() {
			TelSeqWorker worker = new TelSeqWorker(inputBams[index], additionalArgs, outputDir, log);
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

	private static Ran telSeqIt(String inputBam, String output, ArrayList<String> additionalArgs, Logger log) {
		String[] outputs = new String[] { output };
		String[] input = new String[] { inputBam };
		ArrayList<String> command = new ArrayList<String>();
		command.add("telseq");
		command.add("-o");
		command.add(output);
		command.addAll(additionalArgs);

		command.add(inputBam);

		boolean valid = CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", input, outputs, true, false, false, log);
		return new Ran(valid, command);
	}

	public static void run(String sraDir, String outDir, String optionalBed, int threads) {
		ArrayList<SRAConversionResult> conv = SRAUtils.run(sraDir, outDir, threads);
		String telseqDir = outDir + "telseq/";
		new File(telseqDir).mkdirs();
		Logger log = new Logger(telseqDir + ".telseq.log");
		log.reportTimeInfo("Assuming telseq is on system path");
		String[] bams = new String[conv.size()];
		for (int i = 0; i < bams.length; i++) {

		}

		ArrayList<TelSeqResult> results = new ArrayList<TelSeq.TelSeqResult>();
		ArrayList<String> argPopulator = new ArrayList<String>();
		String baseDir = telseqDir + "base/";
		new File(baseDir).mkdirs();
		runType(threads, log, bams, results, argPopulator, baseDir);

		System.exit(1);
		if (optionalBed != null) {
			if (Files.exists(optionalBed)) {
				ArrayList<String> argPopulatorBed = new ArrayList<String>();
				argPopulatorBed.addAll(argPopulator);
				argPopulatorBed.add("-e");
				argPopulatorBed.add(optionalBed);
				String dirBed = telseqDir + ext.rootOf(optionalBed) + "/";
				new File(dirBed).mkdirs();
				runType(threads, log, bams, results, argPopulatorBed, dirBed);

				ArrayList<String> argPopulatorBuffBed = new ArrayList<String>();
				argPopulatorBuffBed.addAll(argPopulator);

				String buffDir = telseqDir + "buff_20KB" + ext.rootOf(optionalBed) + "/";
				new File(buffDir).mkdirs();

				BEDFileReader reader = new BEDFileReader(optionalBed, false);

				LocusSet<BEDFeatureSeg> segs = reader.loadAll(log);

			} else {
				log.reportFileNotFound(optionalBed);
			}
		}

	}

	private static void runType(int threads, Logger log, String[] bams, ArrayList<TelSeqResult> results, ArrayList<String> argPopulator, String baseDir) {
		TelSeqProducer producer = new TelSeqProducer(bams, argPopulator, baseDir, log);
		WorkerTrain<TelSeqResult> train = new WorkerTrain<TelSeq.TelSeqResult>(producer, threads, 10, log);
		while (train.hasNext()) {
			results.add(train.next());

		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String sraDir = "/scratch.global/lanej/aric_raw/sra/";
		String outDir = "/scratch.global/lanej/aric_raw/";
		String captureBed = "/home/pankrat2/public/bin/ref/VCRome_2_1_hg19_capture_targets.bed";
		int threads = 24;

		String usage = "\n" +
				"telomere.SRAUtils requires 0-1 arguments\n" +
				"   (1) SRA directory (i.e. sraDir=" + sraDir + " (default))\n" +
				"   (2) out directory (i.e. outDir=" + outDir + " (default))\n" +
				"   (3) capture bed (i.e. bed=" + outDir + " (default))\n" +

				PSF.Ext.getNumThreadsCommand(3, threads) +
				"";

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
			run(sraDir, outDir, captureBed, threads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
