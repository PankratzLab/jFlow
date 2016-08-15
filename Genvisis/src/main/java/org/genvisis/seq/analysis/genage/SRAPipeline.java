package org.genvisis.seq.analysis.genage;

import java.io.File;
import java.util.concurrent.Callable;

import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.sra.SRAUtils;
import org.genvisis.sra.SRAUtils.SRABamWorker;
import org.genvisis.sra.SRAUtils.SRAConversionResult;

/**
 * more specific version of {@link Pipeline} that starts with a single SRA file
 *
 */
public class SRAPipeline implements Callable<Boolean> {

	private String inputSRA;
	private String rootOutDir;
	private String referenceGenome;
	private String captureBed;
	private ASSAY_TYPE atype;
	private ASSEMBLY_NAME aName;
	private int numThreads;
	private Logger log;

	/**
	 * @param inputSRA
	 *            the input sra file in appropriate sra-toolkit directory
	 * @param rootOutDir
	 *            the output directory for the analysis
	 * @param referenceGenome
	 *            proper reference genome
	 * @param captureBed
	 *            the capture bed, only utilized with {@link ASSAY_TYPE#WXS}
	 * @param atype
	 *            {@link ASSAY_TYPE} of the sample
	 * @param aName
	 *            {@link ASSEMBLY_NAME} for the sample
	 * @param numThreads
	 *            number of threads for the pipeline branches
	 * @param log
	 */
	public SRAPipeline(String inputSRA, String rootOutDir, String referenceGenome, String captureBed, ASSAY_TYPE atype,
			ASSEMBLY_NAME aName, int numThreads, Logger log) {
		super();
		this.inputSRA = inputSRA;
		this.rootOutDir = rootOutDir;
		this.referenceGenome = referenceGenome;
		this.captureBed = captureBed;
		this.atype = atype;
		this.aName = aName;
		this.numThreads = numThreads;
		this.log = log;
	}

	@Override
	public Boolean call() throws Exception {
		String bamDir = rootOutDir + "bams/";
		new File(bamDir).mkdirs();
		String bam = bamDir + ext.rootOf(inputSRA) + ".bam";
		WorkerHive<SRAConversionResult> hive = new WorkerHive<SRAUtils.SRAConversionResult>(1, 10, log);
		hive.addCallable(new SRABamWorker(inputSRA, bam, log));
		hive.execute(true);

		Pipeline.pipeline(bam, rootOutDir, referenceGenome, captureBed, atype, aName, numThreads, log);

		return true;
	}

}
