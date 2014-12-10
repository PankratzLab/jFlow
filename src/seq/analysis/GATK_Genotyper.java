package seq.analysis;

import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import common.Files;
import common.Logger;
import common.ext;

public class GATK_Genotyper {
	private GATK gatk;
	private GATK.SingleSampleHaplotypeCaller[] siSampleHaplotypeCallers;
	private boolean fail, verbose;
	private int numBetweenSampleThreads;
	private Logger log;

	public GATK_Genotyper(GATK gatk, int numBetweenSampleThreads, boolean verbose, Logger log) {
		super();
		this.gatk = gatk;
		this.numBetweenSampleThreads = numBetweenSampleThreads;
		this.log = log;
	}

	public boolean isFail() {
		return fail;
	}

	public void setFail(boolean fail) {
		this.fail = fail;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	public void runSingleSampleAllSites(String[] inputBams) {

		if (!isFail() && Files.checkAllFiles("", inputBams, verbose, log)) {
			if (inputBams != null) {
				this.siSampleHaplotypeCallers = new GATK.SingleSampleHaplotypeCaller[inputBams.length];
				ExecutorService executor = Executors.newFixedThreadPool(numBetweenSampleThreads);
				Hashtable<String, Future<GATK.SingleSampleHaplotypeCaller>> tmpResults = new Hashtable<String, Future<GATK.SingleSampleHaplotypeCaller>>();
				for (int i = 0; i < inputBams.length; i++) {
					Logger altLog = new Logger(ext.rootOf(inputBams[i], false) + ".HC_ERC.log");
					tmpResults.put(i + "", executor.submit(new WorkerSingleSampleAllSites(gatk, inputBams[i], ext.rootOf(inputBams[i]), altLog)));
				}
				for (int i = 0; i < siSampleHaplotypeCallers.length; i++) {
					try {
						siSampleHaplotypeCallers[i] = tmpResults.get(i + "").get();
						if (siSampleHaplotypeCallers[i].isFail() && !isFail()) {
							log.reportError("Error - failed single sample haplotype calling for " + siSampleHaplotypeCallers[i].getInputBam());
							setFail(true);
						}
					} catch (InterruptedException e) {
						log.reportError("Error - when running GATK single sample haplotype calling on internal index " + i);
						log.reportException(e);
						setFail(true);
					} catch (ExecutionException e) {
						log.reportError("Error - when running GATK single sample haplotype calling on internal index " + i);
						log.reportException(e);
						setFail(true);
					}
				}
				executor.shutdown();
				try {
					executor.awaitTermination(10, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					log.reportException(e);
				}
			} else {
				// TODO better check
			}
		}

	}

	private static class WorkerSingleSampleAllSites implements Callable<GATK.SingleSampleHaplotypeCaller> {
		private GATK GATK;
		private String inputBam, baseId;
		private Logger altLog;

		public WorkerSingleSampleAllSites(seq.analysis.GATK gATK, String inputBam, String baseId, Logger altLog) {
			super();
			this.GATK = gATK;
			this.inputBam = inputBam;
			this.baseId = baseId;
			this.altLog = altLog;
		}

		@Override
		public GATK.SingleSampleHaplotypeCaller call() {
			return GATK.haplotypeCallABam(baseId, inputBam, altLog);
		}
	}

}
