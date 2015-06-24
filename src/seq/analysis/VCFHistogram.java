package seq.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Set;
import java.util.concurrent.Callable;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import common.Array;
import common.Files;
import common.Logger;
import common.PSF;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;
import seq.manage.ReferenceGenome;
import seq.manage.VCFOps;
import seq.manage.VCOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import stats.Histogram.DynamicHistogram;
import stats.Rscript.COLUMNS_MULTIPLOT;
import stats.Rscript.PLOT_DEVICE;
import stats.Rscript.RScatter;
import stats.Rscript.RScatters;
import stats.Rscript.SCATTER_TYPE;

/**
 * @author lane0212 Currently used to make histograms of the differences between cases and controls
 */
public class VCFHistogram implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static final String[] METRICS_TRACKED = new String[] { "Avg_GQ", "Avg_DP", "GC_100bp" };

	public enum HIST_WALKER {
		/**
		 * Two histograms for each metric, one for CASE only, and one for Case with Control
		 */
		CASE_V_CONTROL_PRESENT;

	}

	private String vcfFile;
	private VcfPopulation vpop;
	private HIST_WALKER walker;
	private Logger log;

	private DynamicHistogram[][] histograms;
	private String[][] histTitles;

	public VCFHistogram(String vcfFile, String vpopFile, HIST_WALKER walker, Logger log) {
		super();
		this.vcfFile = vcfFile;
		VCFOps.verifyIndex(vcfFile, log);
		this.vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.CASE_CONTROL, log);
		vpop.report();
		this.walker = walker;
		this.log = log;
		initHists();
	}

	public void writeSerial(String fileName) {
		Files.writeSerial(this, fileName, true);
	}

	public static VCFHistogram readSerial(String filename, Logger log) {
		return (VCFHistogram) Files.readSerial(filename, false, log, false, true);
	}

	public RScatters dumpAndPlot(String dir, String root) {
		ArrayList<RScatter> rScatters = new ArrayList<RScatter>();

		for (int i = 0; i < histograms[0].length; i++) {
			ArrayList<String> tmpTitles = new ArrayList<String>();
			ArrayList<DynamicHistogram> tmpHists = new ArrayList<DynamicHistogram>();
			String output = dir + root;
			for (int j = 0; j < histograms.length; j++) {
				tmpTitles.add(histTitles[j][i] + "_Proportion");
				tmpHists.add(histograms[j][i]);
				output += "_" + histTitles[j][i];

			}
			output += "_" + METRICS_TRACKED[i] + ".txt";

			String[] aTmpTitles = tmpTitles.toArray(new String[tmpTitles.size()]);
			DynamicHistogram.dumpToSameFile(tmpHists.toArray(new DynamicHistogram[tmpHists.size()]), aTmpTitles, output, true, log);
			RScatter rScatter = new RScatter(output, output + ".rscript", ext.rootOf(output), output + ".pdf", "Bin", aTmpTitles, SCATTER_TYPE.POINT, log);
			rScatter.setOverWriteExisting(true);
			rScatter.setyLabel("Proportion");
			rScatter.setxLabel(METRICS_TRACKED[i]);
			rScatter.setTitle(root);
			rScatter.execute();
			rScatters.add(rScatter);
		}
		RScatters rScattersAll = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]), dir + root + "full.rscript", dir + root + "full.pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, log);
		rScattersAll.execute();
		return rScattersAll;
	}

	public void populateHists(ReferenceGenome referenceGenome, double maf) {
		VCFFileReader reader = new VCFFileReader(vcfFile, false);
		Set<String> cases = vpop.getSubPop().get(VcfPopulation.CASE);
		Set<String> controls = vpop.getSubPop().get(VcfPopulation.CONTROL);
		int count = 0;
		int caseUniqueCount = 0;
		int caseControlShared = 0;
		for (VariantContext vc : reader) {
			count++;
			if (count % 50000 == 0) {
				log.reportTimeInfo("Scanned " + count + " variants " + caseUniqueCount + " were unique to cases " + caseControlShared + " were shared");
//				reader.close();
//				return;
			}
			VariantContext vcCase = VCOps.getSubset(vc, cases, VC_SUBSET_TYPE.SUBSET_STRICT, true);
			VariantContext vcControls = VCOps.getSubset(vc, controls, VC_SUBSET_TYPE.SUBSET_STRICT, true);

			double mafCase = VCOps.getMAF(vcCase, null);
			int mafControl = VCOps.getNumWithAlt(vcControls);

			if (mafCase > maf) {
				switch (walker) {
				case CASE_V_CONTROL_PRESENT:
					GenotypesContext genotypesContextCase = vcCase.getGenotypes();
					double avgDP = 0;
					double avgGQ = 0;
					double avgGC = -1;
					for (Genotype g : genotypesContextCase) {
						avgDP += g.getDP();
						avgGQ += g.getGQ();
					}
					avgDP /= vcCase.getSampleNames().size();
					avgGQ /= vcCase.getSampleNames().size();
					avgGC = referenceGenome.getGCContentFor(vc);
					avgGC = avgGC * 100.0;
					if (mafControl > 0) {
						caseControlShared++;

						histograms[0][0].addDataPointToHistogram(avgGQ);
						histograms[0][1].addDataPointToHistogram(avgDP);
						histograms[0][2].addDataPointToHistogram(avgGC);

						// }
					} else {
						caseUniqueCount++;
						// for (Genotype g : genotypesContextCase) {
						histograms[1][0].addDataPointToHistogram(avgGQ);
						histograms[1][1].addDataPointToHistogram(avgDP);
						histograms[1][2].addDataPointToHistogram(avgGC);
						// }
					}
					break;
				default:
					log.reportTimeError("Invalid walker " + walker);
					break;
				}
			}

		}
		reader.close();
	}

	private void initHists() {
		switch (walker) {
		case CASE_V_CONTROL_PRESENT:
			this.histograms = new DynamicHistogram[2][METRICS_TRACKED.length];
			this.histTitles = new String[2][3];
			for (int i = 0; i < histograms.length; i++) {
				for (int j = 0; j < histograms[i].length; j++) {
					histTitles[i][j] = i == 0 ? "SHARED_Variants" : "UNIQ_Variants";

					if (j == 2) {
						histograms[i][j] = new DynamicHistogram(0, 100, 0);// GC

					} else if (j == 0) {
						histograms[i][j] = new DynamicHistogram(0, 99, 0);// GQ,DP

					} else {
						histograms[i][j] = new DynamicHistogram(0, 150, 0);// GQ,DP
					}
				}
			}
			break;
		default:
			log.reportTimeError("Invalid walker " + walker);
			break;
		}
	}

	private static String[] divideToNewCaseControlStatus(VcfPopulation vcfPopulation, String dir, Logger log) {
		ArrayList<String> newVpops = new ArrayList<String>();
		for (String superPop : vcfPopulation.getSuperPop().keySet()) {
			for (String superPopComp : vcfPopulation.getSuperPop().keySet()) {
				if (!superPop.equals(superPopComp)) {
					Set<String> newCases = vcfPopulation.getSuperPop().get(superPop);
					Set<String> newControls = vcfPopulation.getSuperPop().get(superPopComp);
					String newFile = dir + superPop + "_V_" + superPopComp + ".vpop";
					try {
						PrintWriter writer = new PrintWriter(new FileWriter(newFile));
						writer.println(Array.toStr(VcfPopulation.HEADER));
						for (String newCase : newCases) {
							writer.println(newCase + "\t" + VcfPopulation.CASE + "\t" + superPop);
						}
						for (String newControl : newControls) {
							writer.println(newControl + "\t" + VcfPopulation.CONTROL + "\t" + superPop);
						}
						writer.close();
						VcfPopulation.load(newFile, POPULATION_TYPE.CASE_CONTROL, log).report();
						newVpops.add(newFile);

					} catch (Exception e) {
						log.reportError("Error writing to " + newFile);
						log.reportException(e);
					}
				}
			}
		}

		return newVpops.toArray(new String[newVpops.size()]);

	}

	private static class HistInit {
		private String vpop;
		private double maf;

		public HistInit(String vpop, double maf) {
			super();
			this.vpop = vpop;
			this.maf = maf;
		}

		public String getVpop() {
			return vpop;
		}

		public double getMaf() {
			return maf;
		}

	}

	private static class HistProducer implements Producer<RScatters> {
		private HistInit[] histInits;
		private String vcf, outputDir;
		private String referenceGenomeFasta;
		private Logger log;
		private int index = 0;

		public HistProducer(HistInit[] histInits, String vcf, String outputDir, String referenceGenomeFasta, Logger log) {
			super();
			this.histInits = histInits;
			this.vcf = vcf;
			this.referenceGenomeFasta = referenceGenomeFasta;
			this.outputDir = outputDir;
			this.log = log;
		}

		@Override
		public boolean hasNext() {
			return index < histInits.length;
		}

		@Override
		public Callable<RScatters> next() {
			final String currentRoot = ext.rootOf(histInits[index].getVpop());
			HistWorker worker = new HistWorker(histInits[index], vcf, referenceGenomeFasta, outputDir, currentRoot, log);
			index++;
			return worker;
		}

		@Override
		public void remove() {

		}

		@Override
		public void shutdown() {

		}

	}

	private static class HistWorker implements Callable<RScatters> {
		private HistInit histInits;
		private String vcf;
		private String referenceGenomeFasta;
		private String outputDir;
		private String outputRoot;
		private Logger log;

		public HistWorker(HistInit histInits, String vcf, String referenceGenomeFasta, String outputDir, String outputRoot, Logger log) {
			super();
			this.histInits = histInits;
			this.vcf = vcf;
			this.referenceGenomeFasta = referenceGenomeFasta;
			this.outputDir = outputDir;
			this.outputRoot = outputRoot;
			this.log = log;
		}

		@Override
		public RScatters call() throws Exception {
			return createHistogram(vcf, histInits.getVpop(), outputDir, outputRoot, referenceGenomeFasta, histInits.getMaf(), log);
		}
	}

	private static RScatters createHistogram(String vcf, String vpop, String outputDir, String outputRoot, String referenceGenomeFasta, double maf, Logger log) {

		outputDir = outputDir == null ? ext.parseDirectoryOfFile(vpop) : outputDir;
		outputRoot += ".MAF_" + ext.formDeci(maf, 2);
		System.out.println(outputRoot);
		new File(outputDir).mkdirs();
		String ser = outputDir + outputRoot + ".ser";
		VCFHistogram histogram = null;
		if (Files.exists(ser)) {
			log.reportTimeInfo("Loading serialized data " + ser);
			histogram = readSerial(ser, log);
		} else {
			ReferenceGenome referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);
			referenceGenome.setDefaultBuffer(100);
			histogram = new VCFHistogram(vcf, vpop, HIST_WALKER.CASE_V_CONTROL_PRESENT, log);
			histogram.populateHists(referenceGenome, maf);
			histogram.writeSerial(ser);
		}
		return histogram.dumpAndPlot(outputDir, outputRoot);
	}

	public static void createHistograms(String vcf, String vpop, String outputDir, String outputRoot, String referenceGenomeFasta, double[] mafs, int numthreads, Logger log) {
		VcfPopulation vcfPopulation = VcfPopulation.load(vpop, POPULATION_TYPE.ANY, log);

		
		String[] allIters = divideToNewCaseControlStatus(vcfPopulation, outputDir, log);
		HistInit[] histInits = new HistInit[mafs.length * allIters.length];
		int index = 0;
		for (int i = 0; i < mafs.length; i++) {
			for (int j = 0; j < allIters.length; j++) {
				histInits[index] = new HistInit(allIters[j], mafs[i]);
				System.out.println(histInits[index].getMaf()+"\t"+histInits[index].getVpop());
				index++;
			}
		}
		log.reportTimeInfo(histInits.length +" total comparisions");

		HistProducer producer = new HistProducer(histInits, vcf, outputDir, referenceGenomeFasta, log);
		WorkerTrain<RScatters> train = new WorkerTrain<RScatters>(producer, numthreads, numthreads, log);
		while (train.hasNext()) {
			train.next();
		}
	}	

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = "aVcf.vcf";
		String vpop = "vpop.vpop";
		String outputDir = null;
		String outputRoot = "hists";
		String referenceGenomeFasta = "ref.fa";
		String logfile = null;
		Logger log;
		int numthreads = 2;
		double[] mafs = new double[] {  .05,0 };

		String usage = "\n" + "seq.analysis.VCFHistogram requires 0-1 arguments\n";
		usage += "   (1) vcf (i.e. vcf=" + vcf + " (default))\n" + "";
		usage += "   (2) vpop (i.e. vpop=" + vcf + " (default))\n" + "";
		usage += PSF.Ext.getOutputDirCommand(3, null);
		usage += "   (4) output root (i.e. root=" + outputRoot + " (default))\n" + "";
		usage += "   (5) referenece genome (i.e. ref=" + referenceGenomeFasta + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(6, numthreads);

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcf=")) {
				vcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("vpop=")) {
				vpop = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
				outputDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("root=")) {
				outputRoot = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ref=")) {
				referenceGenomeFasta = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
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
			log = new Logger(logfile);
			createHistograms(vcf, vpop, outputDir, outputRoot, referenceGenomeFasta, mafs, numthreads, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
