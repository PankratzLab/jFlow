package org.genvisis.seq.qc.contamination;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.PennCNV;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.stats.LeastSquares;
import org.genvisis.stats.LeastSquares.LS_TYPE;
import org.genvisis.stats.RegressionModel;

import com.google.common.primitives.Doubles;

public class BAFContamination {
	private final double[] sampBAF;
	private final double[] mafs, callRate;
	private final byte[] sampGenotypes;
	private final double minFreq, minCallRate;
	private final boolean fail, verbose;
	private static double MIN_MAF = 0.05;
	private static double MIN_CALL_RATE = 0.98;

	private final Project proj;
	private final Logger log;

	public BAFContamination(Project proj, final double[] sampBAF, final byte[] sampGenotypes,
													final double[] maf, final double[] callRate, double minFreq,
													double minCallRate, boolean verbose, Logger log) {
		super();
		this.proj = proj;
		this.sampBAF = sampBAF;
		this.minFreq = minFreq;
		this.callRate = callRate;
		this.sampGenotypes = sampGenotypes;
		mafs = maf;
		this.verbose = verbose;
		this.minCallRate = minCallRate;
		this.log = log;

		fail = !verify();
		System.out.println(fail + "|" + this.verbose);
	}

	public BAFContaminationResults getContamination() {
		byte[] chrs = proj.getMarkerSet().getChrs();
		int subIndex = Array.indexOfFirstMaxByte(chrs, (byte) 23);
		// mafs = getMafs(mafs);
		boolean[] use = getPFBsToUse(	getMafs(mafs), callRate, sampGenotypes, minFreq, minCallRate,
																	subIndex);
		double[][] indeps = new double[Array.booleanArraySum(use)][2];

		log.reportTimeInfo("Found "	+ Array.booleanArraySum(use)
												+ " homozygous autosomal markers passing freq threashold of " + minFreq
												+ " and callrate " + minCallRate);
		int index = 0;
		for (int i = 0; i < indeps.length; i++) {
			if (use[i]) {
				double pfbA = 1 - mafs[i];

				indeps[index][0] = sampGenotypes[i] == 2 ? mafs[i] : -1 * pfbA;
				// System.out.println( sampGenotypes[i]+"\t"+sampBAF[i]+"\t"+mafs[i]+"\t"+indeps[index][0]);
				// indeps[index][0] = mafs[i];

				indeps[index][1] = sampGenotypes[i] == 2 ? 0 : 1;
				index++;
			}
		}

		RegressionModel model = new LeastSquares(	Array.subArray(sampBAF, use), indeps, null, false,
																							verbose, LS_TYPE.REGULAR);
		BAFContaminationResults results = new BAFContaminationResults(model.getOverallSig(),
																																	model.getRsquare(),
																																	model.getBetas());
		results.setsAlleleDeviation(new StandardAlleleDeviation(Array.subArray(sampGenotypes, use),
																														Array.subArray(sampBAF, use)));
		// System.out.println(Array.toStr(results.getBetas()));
		return results;
	}

	private static boolean[] getPFBsToUse(double[] maf, double[] callRate, byte[] sampGenotypes,
																				double minFreq, double minCallRate, int sexStart) {
		boolean[] use = new boolean[maf.length];
		for (int i = 0; i < use.length; i++) {
			if (maf[i] >= minFreq	&& (sampGenotypes[i] == 0 || sampGenotypes[i] == 2) && i < sexStart
					&& callRate[i] >= minCallRate) {
				use[i] = true;
			}
		}
		return use;
	}

	private static double[] getMafs(double[] pfb) {
		double[] minpfb = new double[pfb.length];
		for (int i = 0; i < minpfb.length; i++) {
			minpfb[i] = Math.min(pfb[i], 1 - pfb[i]);
		}
		return minpfb;
	}

	// private static double[] getBAFDiffs(double[] pfb, double[] bafs) {
	// double[] diffs = new double[pfb.length];
	// for (int i = 0; i < diffs.length; i++) {
	// diffs[i] = Math.abs(pfb[i] - bafs[i]);
	// }
	// return diffs;
	// }

	private boolean verify() {
		boolean verfied = true;
		if (sampBAF.length != mafs.length) {
			log.reportError("Mismatched array sizes for sample baf and pfb");
			verfied = false;
		}

		return verfied;
	}

	public static class BAFContaminationResults {
		private final double pval;
		private final double rsquared;
		private final double[] betas;
		private StandardAlleleDeviation sAlleleDeviation;

		public BAFContaminationResults(double pval, double rsquared, double[] betas) {
			super();
			this.pval = pval;
			this.rsquared = rsquared;
			this.betas = betas;
		}

		public double getPval() {
			return pval;
		}

		public double getRsquared() {
			return rsquared;
		}

		public double[] getBetas() {
			return betas;
		}

		public StandardAlleleDeviation getsAlleleDeviation() {
			return sAlleleDeviation;
		}

		public void setsAlleleDeviation(StandardAlleleDeviation sAlleleDeviation) {
			this.sAlleleDeviation = sAlleleDeviation;
		}

	}

	private static class StandardAlleleDeviation {
		private final double AA_stdev;
		private final double BB_stdev;

		// private double AB_stdev;

		public StandardAlleleDeviation(byte[] genotypes, double[] bafs) {
			AA_stdev = getStdev((byte) 0, genotypes, bafs);
			// this.AB_stdev = getStdev((byte) 1, genotypes, bafs);
			BB_stdev = getStdev((byte) 2, genotypes, bafs);
		}

		private static double getStdev(byte geno, byte[] genotypes, double[] bafs) {
			ArrayList<Double> tmp = new ArrayList<Double>();
			for (int i = 0; i < bafs.length; i++) {
				if (genotypes[i] == geno) {
					tmp.add(geno == 2 ? bafs[i] : 1 - bafs[i]);
				}
			}
			double[] tmp2 = Doubles.toArray(tmp);
			double cv = Array.stdev(tmp2, true) / Array.mean(tmp2, true);
			System.out.println(geno + "\t" + cv);

			return cv;
		}

		public double getAA_stdev() {
			return AA_stdev;
		}

		public double getBB_stdev() {
			return BB_stdev;
		}

		// public double getAB_stdev() {
		// return AB_stdev;
		// }

	}

	private static class ContaminationProducer extends AbstractProducer<BAFContaminationResults> {
		private final Project proj;
		private final double[] pfb, callRate;
		private final String[] samples;
		private int index;

		public ContaminationProducer(Project proj, double[] pfb, double[] callRate, String[] samples) {
			super();
			this.proj = proj;
			this.pfb = pfb;
			this.samples = samples;
			index = 0;
			this.callRate = callRate;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<BAFContaminationResults> next() {
			BAFContaminationWorker worker =
																		new BAFContaminationWorker(proj, samples[index], pfb, callRate);
			index++;
			return worker;
		}
	}

	private static class BAFContaminationWorker implements Callable<BAFContaminationResults> {
		private final Project proj;
		private final String sample;
		private final double[] pfbs;
		private final double[] callRate;

		public BAFContaminationWorker(Project proj, String sample, double[] pfbs, double[] callRate) {
			super();
			this.proj = proj;
			this.sample = sample;
			this.pfbs = pfbs;
			this.callRate = callRate;
		}

		@Override
		public BAFContaminationResults call() throws Exception {
			Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
			ClusterFilterCollection clusterFilterCollection = new ClusterFilterCollection();
			byte[] genos = samp.getAB_GenotypesAfterFilters(proj.getMarkerNames(),
																											clusterFilterCollection, (float) 0.20);
			BAFContamination bafContamination = new BAFContamination(	proj,
																																Array.toDoubleArray(samp.getBAFs()),
																																genos, pfbs, callRate, MIN_MAF,
																																MIN_CALL_RATE, true, proj.getLog());
			BAFContaminationResults results = bafContamination.getContamination();
			return results;
		}

	}

	private static ExtProjectDataParser getParser(Project proj, String file, String dataKeyColumnName,
																								String[] numericColumns) throws FileNotFoundException {
		ExtProjectDataParser.ProjectDataParserBuilder builder =
																													new ExtProjectDataParser.ProjectDataParserBuilder();
		builder.sampleBased(false);
		builder.treatAllNumeric(false);
		builder.requireAll(true);
		builder.verbose(true);
		builder.dataKeyColumnName(dataKeyColumnName);
		builder.numericDataTitles(numericColumns);

		// builder.headerFlags(Array.concatAll(new String[]{dataKeyColumnName}, numericColumns));
		return builder.build(proj, file);
	}

	public static void test(Project proj, int numThreads) {
		// Project proj = new
		// Project("C:/workspace/Genvisis/projects/testVCFNORMALIZED_GC_CORRECTED.properties", false);
		// String[] testSamps = new String[] { "PT21_03_AGGCAGAA-ACTGCATA",
		// "PT292_03_TAAGGCGA-CTAAGCCT", "PT280_03_TAAGGCGA-AAGGAGTA", "PT99_03_AGGCAGAA-TAGATCGC",
		// "CONTROL_1_TAGGCATG-CTAAGCCT", "CONTROL_2_CTCTCTAC-TAGATCGC" };
		String[] testSamps = proj.getSamples();
		String pfb = proj.PROJECT_DIRECTORY.getValue() + "custom.pfb";
		String mafs = proj.PROJECT_DIRECTORY.getValue() + "mafs.txt";
		if (!Files.exists(mafs)) {
			proj.getLog().reportTimeInfo("Creating " + mafs + " for contamination detection");
			MAF.computeMAF(proj);
		}
		if (!Files.exists(pfb)) {
			proj.getLog().reportTimeInfo("Creating " + pfb + " for contamination detection");
			PennCNV.populationBAF(proj);
		}
		if (!Files.exists(proj.SAMPLE_QC_FILENAME.getValue())) {
			proj.getLog().reportTimeInfo("Creating "	+ proj.SAMPLE_QC_FILENAME.getValue()
																		+ " for contamination detection");
			LrrSd.init(proj, null, null, numThreads);
		}

		// builder.numericDataTitles(new String[] { "PFB" });
		// builder.numericDataTitles(new String[] { "MAF", "CALLRATE" });
		// builder.headerFlags(new String[] { "NAME", "MAF", "CALLRATE" });

		// public static SampleQC loadSampleQC(Project proj, String sampleColumnName, String[]
		// qcTitlesToLoad, boolean generateSampleQC, boolean gcCorrectedLrrSd, String duplicatesSetFile)
		// {

		SampleQC sampleQC = SampleQC.loadSampleQC(proj, LrrSd.SAMPLE_COLUMN,
																							new String[] {"LRR_SD", "BAF1585_SD", "AB_callrate",},
																							true, true, null);

		try {
			ExtProjectDataParser parserPfb = getParser(proj, pfb, "Name", new String[] {"PFB"});
			ExtProjectDataParser parserMaf = getParser(	proj, mafs, "NAME",
																									new String[] {"MAF", "CALLRATE"});
			parserPfb.determineIndicesFromTitles();
			parserMaf.determineIndicesFromTitles();
			parserPfb.loadData();
			parserMaf.loadData();

			double[] pfbs = parserPfb.getNumericData()[0];
			double[] callRates = parserMaf.getNumericData()[1];
			proj.getLog().reportTimeInfo("found " + pfbs.length + "pfbs");
			String output = proj.PROJECT_DIRECTORY.getValue() + "contamination.txt";
			ContaminationProducer producer = new ContaminationProducer(	proj, pfbs, callRates,
																																	proj.getSamples());
			WorkerTrain<BAFContaminationResults> train = new WorkerTrain<BAFContamination.BAFContaminationResults>(	producer,
																																																							numThreads,
																																																							1,
																																																							proj.getLog());

			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println("Sample\tLRR_SD\tBAF1585_SD\tAB_callrate\tAA_STDEV\tBB_STDEV\tContamBeta");
				int index = 0;
				while (train.hasNext()) {
					BAFContaminationResults results = train.next();
					// System.out.println(results.getBetas().length);
					String result = testSamps[index]	+ "\t" + sampleQC.getDataFor("LRR_SD")[index] + "\t"
													+ sampleQC.getDataFor("BAF1585_SD")[index] + "\t"
													+ sampleQC.getDataFor("AB_callrate")[index];
					result += "\t"	+ results.getsAlleleDeviation().getAA_stdev() + "\t"
										+ results.getsAlleleDeviation().getBB_stdev();
					result += "\t"	+ sampleQC.getDataFor("AB_callrate")[index] + "\t"
										+ Array.toStr(results.getBetas());
					writer.println(result);
					index++;
					writer.flush();
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + output);
				proj.getLog().reportException(e);
			}

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static void main(String[] args) {
		// int numArgs = args.length;
		Project proj = new Project(null, false);
		test(proj, 4);
	}

}
