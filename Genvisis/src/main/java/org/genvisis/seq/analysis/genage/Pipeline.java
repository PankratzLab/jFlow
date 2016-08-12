package org.genvisis.seq.analysis.genage;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.analysis.MitoSeqCN;
import org.genvisis.seq.telomere.TelSeq;

/**
 * Going to be the pipeline of execution for a single input bam file
 *
 */
public class Pipeline {

	private Pipeline() {

	}

	private static final String MITO_DIR = "mtDNACN/";
	private static final String COMPUTEL_DIR = "computel/";
	private static final String TELSEQ_DIR = "telseq/";
	private static final String GENVISIS_DIR = "genvisis/";

	private abstract static class PipelinePart implements Callable<PipelinePart> {
		private List<String> output;

		private PipelinePart() {
		}

		protected void setOutput(List<String> output) {
			this.output = output;
		}

	}

	private static class MitoPipeResult extends PipelinePart {
		private String bamFile;
		private String rootOutDir;
		private String captureBed;
		private String referenceGenomeFasta;
		private ASSEMBLY_NAME aName;
		private ASSAY_TYPE aType;
		private int numthreads;
		private Logger log;

		private MitoPipeResult(String bamFile, String rootOutDir, String captureBed, String referenceGenomeFasta,
				ASSEMBLY_NAME aName, ASSAY_TYPE aType, int numthreads, Logger log) {
			super();
			this.bamFile = bamFile;
			this.rootOutDir = rootOutDir;
			this.captureBed = captureBed;
			this.referenceGenomeFasta = referenceGenomeFasta;
			this.aName = aName;
			this.aType = aType;
			this.numthreads = numthreads;
			this.log = log;
		}

		@Override
		public PipelinePart call() throws Exception {
			String mitoDir = rootOutDir + MITO_DIR + ext.rootOf(bamFile) + "/";
			String bamList = mitoDir + "bam.list.txt";
			new File(mitoDir).mkdirs();
			Files.write(bamFile, bamList);
			String result = MitoSeqCN.run(bamList, rootOutDir, aType == ASSAY_TYPE.WGS ? null : captureBed,
					referenceGenomeFasta, aName, numthreads, log);
			ArrayList<String> output = new ArrayList<String>();
			output.add(result);
			setOutput(output);
			return this;
		}

	}

	private static class TelSeqResult extends PipelinePart {

		private String bam;
		private String rootOutDir;
		private String captureBed;
		private ASSEMBLY_NAME aName;
		private ASSAY_TYPE aType;
		private int numthreads;
		private int captureBufferSize;
		private Logger log;

		private TelSeqResult(String bam, String rootOutDir, String captureBed, ASSEMBLY_NAME aName, ASSAY_TYPE aType,
				int numthreads, int captureBufferSize, Logger log) {
			super();
			this.bam = bam;
			this.rootOutDir = rootOutDir;
			this.captureBed = captureBed;
			this.aName = aName;
			this.aType = aType;
			this.numthreads = numthreads;
			this.captureBufferSize = captureBufferSize;
			this.log = log;
		}

		@Override
		public PipelinePart call() throws Exception {
			String telSeqDir = rootOutDir + TELSEQ_DIR + ext.rootOf(bam);
			new File(telSeqDir).mkdir();
			String result = TelSeq.runTelSeq(new String[] { bam }, telSeqDir, captureBed, numthreads, aType, aName,
					captureBufferSize, log);
			ArrayList<String> output = new ArrayList<String>();
			output.add(result);
			setOutput(output);
			return this;
		}

	}

	public static List<PipelinePart> pipeline(String inputBam, String rootOutDir, String referenceGenome,
			String captureBed, ASSAY_TYPE atype, ASSEMBLY_NAME aName, int numThreads, Logger log) {
		if (!Files.exists(inputBam)) {
			throw new IllegalArgumentException("Bam file " + inputBam + " must exist");
		}

		if (!Files.exists(referenceGenome)) {
			throw new IllegalArgumentException("Reference Genome " + inputBam + " must exist");
		} else {
			log.reportTimeWarning("Assuming " + referenceGenome + " matches assembly type " + aName);
		}
		if (atype == ASSAY_TYPE.WXS && (!Files.exists(captureBed))) {
			throw new IllegalArgumentException(captureBed + " must exist");
		}

		WorkerHive<PipelinePart> hive = new WorkerHive<Pipeline.PipelinePart>(numThreads, 10, log);
		// mtDNA CN
		hive.addCallable(new MitoPipeResult(inputBam, rootOutDir, captureBed, referenceGenome, aName, atype, 1, log));

		hive.addCallable(new TelSeqResult(inputBam, rootOutDir, captureBed, aName, atype, 1, 100, log));
		hive.execute(true);

		return hive.getResults();

	}

}
