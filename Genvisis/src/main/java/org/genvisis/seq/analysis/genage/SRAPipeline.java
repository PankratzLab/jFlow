package org.genvisis.seq.analysis.genage;

import java.io.File;
import java.util.Map;
import java.util.concurrent.Callable;

import org.apache.commons.cli.Options;
import org.genvisis.CLI;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.telomere.Computel;
import org.genvisis.sra.SRARunTable;
import org.genvisis.sra.SRASample;
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

	/**
	 * This will be a bit "reversed" in the final pipeline version...This method
	 * is for processing many pre-downloaded files
	 */
	private static void runAll(String sraDir, String sraRunTableFile, int numThreads, Logger log) {
		String[] sraFiles = Files.list(sraDir, ".sra", false);
		SRARunTable srRunTable = SRARunTable.load(sraRunTableFile, log);
		log.reportTimeInfo("Found " + sraFiles.length + " sra files in " + sraDir);
		WorkerHive<SRAPipeline> hive = new WorkerHive<SRAPipeline>(numThreads, 10, log);
		for (int i = 0; i < sraFiles.length; i++) {
			SRASample sample = srRunTable.get(ext.rootOf(sraFiles[i]));
			

		}
	}

	public static void main(String[] args) {
		String sraDir = "sra/";
		String outDir = "out/";

		Options options = CLI.defaultOptions();
		final String sra = "sraDir";
		CLI.addArg(options, sra, "directory with .sra files", sraDir);

		final String outdir = "outDir";
		CLI.addArg(options, outdir, "the output directory for results", outDir);

		
		
		Map<String, String> parsed = CLI.parseWithExit(Computel.class, options, args);
	}

}
