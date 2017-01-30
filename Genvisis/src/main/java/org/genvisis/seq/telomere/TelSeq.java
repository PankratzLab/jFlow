package org.genvisis.seq.telomere;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.manage.BEDFileReader;
import org.genvisis.seq.manage.BEDFileReader.BEDFeatureSeg;
import org.genvisis.seq.manage.BamOps;

/**
 * @author Kitty Class for assisting in generating telomere length estimates with TelSeq
 */
public class TelSeq {

	private static final String[] TELSEQ_REPORT = new String[] {"ReadGroup", "Library", "Sample",
																															"Total", "Mapped", "Duplicates",
																															"LENGTH_ESTIMATE"};

	private enum TYPE {
											BASE, BED, BUFFERED_BED;
	}

	private static class TelSeqResult {
		private final String output;
		private final String sample;
		private final int readSizeUsed;
		private final TYPE type;

		public TelSeqResult(String output, Ran ran, String sample, TYPE type, int readSizeUsed) {
			super();
			this.output = output;
			this.sample = sample;
			this.type = type;
			this.readSizeUsed = readSizeUsed;
		}

	}

	private static class TelSeqWorker implements Callable<TelSeqResult> {
		private final String inputBam;
		private final List<String> additionalArgs;
		private final String outputDir;
		private final TYPE type;

		private final Logger log;

		public TelSeqWorker(String inputBam, List<String> additionalArgs, String outputDir, TYPE type,
												Logger log) {
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
			int readSize = BamOps.estimateReadSize(inputBam, log);
			log.reportTimeInfo("Estimated readsize for " + inputBam + " to be " + readSize);
			Ran ran = telSeqIt(inputBam, out, readSize, additionalArgs, log);
			String sampleName = "NA";
			try {
				sampleName = BamOps.getSampleName(inputBam, log);
			} catch (Exception e) {
				log.reportError("Could not get sample name from " + inputBam);
				log.reportException(e);
			}

			return new TelSeqResult(out, ran, sampleName, type, readSize);
		}

		private static Ran telSeqIt(String inputBam, String output, int readSize,
																List<String> additionalArgs, Logger log) {
			String[] outputs = new String[] {output};
			String[] input = new String[] {inputBam};
			ArrayList<String> command = new ArrayList<String>();
			command.add("telseq");
			command.add("-o");
			command.add(output);
			command.add("-r");
			command.add(Integer.toString(readSize));

			command.addAll(additionalArgs);

			command.add(inputBam);

			boolean valid = CmdLine.runCommandWithFileChecks(	ArrayUtils.toStringArray(command), "", input,
																												outputs, true, false, false, log);
			return new Ran(valid, command);
		}
	}

	private static class TelSeqProducer extends AbstractProducer<TelSeqResult> {
		private final String[] inputBams;
		private final List<String> additionalArgs;
		private final String outputDir;
		private final TYPE type;
		private final Logger log;
		private int index;

		public TelSeqProducer(String[] inputBams, List<String> additionalArgs, String outputDir,
													TYPE type, Logger log) {
			super();
			this.inputBams = inputBams;
			this.additionalArgs = additionalArgs;
			this.outputDir = outputDir;
			this.type = type;
			this.log = log;
			index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < inputBams.length;
		}

		@Override
		public Callable<TelSeqResult> next() {
			TelSeqWorker worker =
													new TelSeqWorker(inputBams[index], additionalArgs, outputDir, type, log);
			index++;
			return worker;
		}

	}

	private static class Ran {
		public Ran(boolean valid, List<String> command) {
			super();
		}

	}

	/**
	 * @param bams bams to run
	 * @param outDir where to put things
	 * @param captureBed for {@link ASSAY_TYPE#WXS}
	 * @param threads number of threads
	 * @param aType see {@link ASSAY_TYPE}
	 * @param aName {@link ASSEMBLY_NAME}
	 * @param captureBufferSize number of base pairs to buffer the caputure bed file
	 */
	public static String runTelSeq(	String[] bams, String outDir, String captureBed, int threads,
																	ASSAY_TYPE aType, ASSEMBLY_NAME aName, int captureBufferSize,
																	Logger log) {

		log.reportTimeInfo("Assuming telseq is on system path");
		ArrayList<TelSeqResult> results = new ArrayList<TelSeq.TelSeqResult>();
		ArrayList<String> argPopulator = new ArrayList<String>();
		argPopulator.add("-m");// doesn't look like telseq handles RGs properly

		switch (aType) {
			case WGS:
				runType(threads, log, bams, results, argPopulator, outDir, TYPE.BASE);
				break;
			case WXS:
				processWXS(	bams, outDir, captureBed, threads, aName, captureBufferSize, log, results,
										argPopulator);
				break;
			default:
				break;

		}

		// summarize
		String finalOut = outDir + "telseq.summary.txt";

		String[] telHeader = Files.getHeaderOfFile(results.get(0).output, log);

		// Only interested in these columns currently, do not know what to do
		// with TEL* and GC*
		int[] indices = ext.indexFactors(TELSEQ_REPORT, telHeader, true, false);
		if (ArrayUtils.countIf(indices, -1) > 0) {
			throw new IllegalStateException("Missing proper heading for " + results.get(0).output);
		}

		ArrayList<String> result = new ArrayList<String>();
		result.add("BAM\t" + ArrayUtils.toStr(TELSEQ_REPORT) + "\tType\tSampleName\tReadSize");
		for (TelSeqResult telSeqResult : results) {
			if (Files.exists(telSeqResult.output)) {
				String[][] data = HashVec.loadFileToStringMatrix(telSeqResult.output, true, null, false);
				for (String[] element : data) {
					result.add(ext.rootOf(telSeqResult.output)	+ "\t"
											+ ArrayUtils.toStr(ArrayUtils.subArray(element, indices)) + "\t" + telSeqResult.type
											+ "\t" + telSeqResult.sample + "\t" + telSeqResult.readSizeUsed);
				}
			}
		}
		Files.writeIterable(result, finalOut);
		return finalOut;
	}

	private static void processWXS(	String[] bams, String outDir, String captureBed, int threads,
																	ASSEMBLY_NAME aName, int captureBufferSize, Logger log,
																	ArrayList<TelSeqResult> results, ArrayList<String> argPopulator) {
		if (Files.exists(captureBed)) {

			BEDFileReader reader = new BEDFileReader(captureBed, false);

			LocusSet<BEDFeatureSeg> segs = reader.loadAll(log);
			reader.close();

			String buffDir = outDir + "buff_" + captureBufferSize + "_" + ext.rootOf(captureBed) + "/";
			new File(buffDir).mkdirs();
			String buffBed = buffDir	+ "buff_" + captureBufferSize + "bp_" + ext.rootOf(captureBed)
												+ ".bed";
			log.reportTimeInfo("writing bed to " + buffBed);
			segs.getBufferedSegmentSet(captureBufferSize)
					.writeSegmentRegions(buffBed, aName == ASSEMBLY_NAME.GRCH37, log);

			ArrayList<String> argPopulatorBuffBed = new ArrayList<String>();
			argPopulatorBuffBed.addAll(argPopulator);
			argPopulatorBuffBed.add("-e");
			argPopulatorBuffBed.add(buffBed);
			runType(threads, log, bams, results, argPopulatorBuffBed, buffDir, TYPE.BUFFERED_BED);
		} else {
			log.reportFileNotFound(captureBed);
		}
	}

	private static void runType(int threads, Logger log, String[] bams,
															ArrayList<TelSeqResult> results, ArrayList<String> argPopulator,
															String baseDir, TYPE type) {
		TelSeqProducer producer = new TelSeqProducer(bams, argPopulator, baseDir, type, log);
		WorkerTrain<TelSeqResult> train = new WorkerTrain<TelSeq.TelSeqResult>(	producer, threads, 100,
																																						log);
		while (train.hasNext()) {
			results.add(train.next());
		}
	}
}
