package org.genvisis.cnv.analysis.pca;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.pca.CorrectionIterator.ITERATION_TYPE;
import org.genvisis.cnv.analysis.pca.CorrectionIterator.MODEL_BUILDER_TYPE;
import org.genvisis.cnv.analysis.pca.CorrectionIterator.ORDER_TYPE;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.cnv.util.Java6Helper;
import org.genvisis.common.AlleleFreq;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.SciStringComparator;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.gwas.windows.BasicHit;
import org.genvisis.gwas.windows.GeneralHitWindowDetector;
import org.genvisis.seq.manage.StrandOps;
import org.genvisis.seq.manage.StrandOps.CONFIG;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.stats.Correlation;
import org.genvisis.stats.LeastSquares;
import org.genvisis.stats.RegressionModel;
import org.genvisis.stats.Rscript.COLUMNS_MULTIPLOT;
import org.genvisis.stats.Rscript.PLOT_DEVICE;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.RScatters;
import org.genvisis.stats.Rscript.Restrictions;
import org.genvisis.stats.Rscript.SCATTER_TYPE;
import org.genvisis.stats.Rscript.VertLine;

import com.google.common.primitives.Doubles;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Going to handle the beta-optimization via correlating effects to a previous meta-analysis
 *
 */
public class BetaOptimizer {
	//

	private static final Segment HG_19_HLA = new Segment((byte) 6, 25000000, 40000000);
	public static final String[] BETA_HEADER = new String[] {"rsID", "ref", "alt", "beta", "beta",
																													 "p"};

	private static final String SUB_DIR = "_eval/typed/";

	private static class BetaCorrelationResult {
		private final int comparisonIndex;
		private final double correlPearsonSigned;
		private final double pPearsonSigned;
		private final double correlSpearmanSigned;
		private final double pSpearmanSigned;

		private final double correlPearsonUnSigned;
		private final double pPearsonUnSigned;
		private final double correlSpearmanUnSigned;
		private final double pSpearmanUnSigned;
		private final int numMarkers;

		public BetaCorrelationResult(int comparisonIndex, double correlPearsonSigned,
																 double pPearsonSigned, double correlSpearmanSigned,
																 double pSpearmanSigned, double correlPearsonUnSigned,
																 double pPearsonUnSigned, double correlSpearmanUnSigned,
																 double pSpearmanUnSigned, int numMarkers) {
			super();
			this.comparisonIndex = comparisonIndex;
			this.correlPearsonSigned = correlPearsonSigned;
			this.pPearsonSigned = pPearsonSigned;
			this.correlSpearmanSigned = correlSpearmanSigned;
			this.pSpearmanSigned = pSpearmanSigned;
			this.correlPearsonUnSigned = correlPearsonUnSigned;
			this.pPearsonUnSigned = pPearsonUnSigned;
			this.correlSpearmanUnSigned = correlSpearmanUnSigned;
			this.pSpearmanUnSigned = pSpearmanUnSigned;
			this.numMarkers = numMarkers;
		}

		private static String[] getHeader() {
			ArrayList<String> header = new ArrayList<String>();
			header.add("PC");
			header.add("correl_pearsonSigned");
			header.add("pval_pearsonSigned");
			header.add("correl_pearsonUnSigned");
			header.add("pval_pearsonUnSigned");
			header.add("correl_spearmanSigned");
			header.add("pval_spearmanSigned");
			header.add("correl_spearmanUnSigned");
			header.add("pval_spearmanUnSigned");
			header.add("numberOfMarkers");
			return ArrayUtils.toStringArray(header);
		}

		private String[] getSummary() {
			ArrayList<String> summary = new ArrayList<String>();
			summary.add(comparisonIndex + "");
			summary.add(correlPearsonSigned + "");
			summary.add(pPearsonSigned + "");
			summary.add(correlPearsonUnSigned + "");
			summary.add(pPearsonUnSigned + "");
			summary.add(correlSpearmanSigned + "");
			summary.add(pSpearmanSigned + "");
			summary.add(correlSpearmanUnSigned + "");
			summary.add(pSpearmanUnSigned + "");
			summary.add(numMarkers + "");
			return ArrayUtils.toStringArray(summary);
		}
	}

	private static double[] inverseNormalizeTo(double[] data, boolean[] samples) {
		double[] newData = ArrayUtils.doubleArray(data.length, Double.NaN);
		double[] tmp = ArrayUtils.inverseNormalize(ArrayUtils.subArray(data, samples));
		int retrieve = 0;
		for (int i = 0; i < newData.length; i++) {
			if (samples[i]) {
				newData[i] = tmp[retrieve];
				retrieve++;
			} else {
				newData[i] = data[i];
			}
		}
		return newData;
	}

	private static class BetaWorker implements Callable<BetaCorrelationResult[]> {
		private final byte[][] genotypes;
		private final double[] data;
		private final ArrayList<MetaBeta> metaBetas;
		private final int comparisonIndex;
		private final boolean[] sampleDef;
		private final double maf;
		private final Logger log;

		public BetaWorker(byte[][] genotypes, boolean[] sampleDef, double[] data, int comparisonIndex,
											ArrayList<MetaBeta> metaBetas, double maf, Logger log) {
			super();
			this.genotypes = genotypes;
			this.data = data;
			this.metaBetas = metaBetas;
			this.comparisonIndex = comparisonIndex;
			this.sampleDef = sampleDef;
			this.maf = maf;
			this.log = log;
		}

		@Override
		public BetaCorrelationResult[] call() throws Exception {
			BetaCorrelationResult result = getResult(data);
			BetaCorrelationResult resultInv = getResult(inverseNormalizeTo(data, sampleDef));
			return new BetaCorrelationResult[] {result, resultInv};
		}

		private BetaCorrelationResult getResult(double[] data) {
			ArrayList<Double> dataBetas = new ArrayList<Double>();
			ArrayList<Double> meta = new ArrayList<Double>();
			int countMafRemoved = 0;
			for (int i = 0; i < genotypes.length; i++) {
				if (i % 2000 == 0) {
					log.reportTimeInfo("computing betas for genotype " + (i + 1) + " of " + genotypes.length
														 + " for index " + comparisonIndex + " on "
														 + Thread.currentThread().getName());
					log.reportTimeInfo("removed " + countMafRemoved + " markers for maf filter of " + maf);
				}
				ArrayList<Double> tmpGeno = new ArrayList<Double>();
				ArrayList<Double> tmpData = new ArrayList<Double>();
				for (int j = 0; j < genotypes[i].length; j++) {
					if (genotypes[i][j] >= 0 && !Double.isNaN(data[j])) {
						tmpGeno.add((double) genotypes[i][j]);
						tmpData.add(data[j]);

					}
				}
				if (tmpGeno.size() > 0) {
					// whew, is this inefficient or what
					double[] genos = Doubles.toArray(tmpGeno);
					int[] tmpInt = ArrayUtils.toIntArray(genos);
					double mafMark = AlleleFreq.calcMAF(ArrayUtils.countIf(tmpInt, 0),
																							ArrayUtils.countIf(tmpInt, 1),
																							ArrayUtils.countIf(tmpInt, 2));
					if (mafMark >= maf) {
						RegressionModel model = new LeastSquares(Doubles.toArray(tmpData), genos);
						if (!model.analysisFailed()) {
							dataBetas.add(model.getBetas()[1]);
							meta.add(metaBetas.get(i).getBeta());
						}
					} else {
						countMafRemoved++;
					}
				}

			}
			double[] correlB = Doubles.toArray(dataBetas);
			double[] correlM = Doubles.toArray(meta);
			double[] pearsonSigned = ArrayUtils.doubleArray(2, Double.NaN);
			double[] spearmanSigned = ArrayUtils.doubleArray(2, Double.NaN);
			double[] pearsonUnSigned = ArrayUtils.doubleArray(2, Double.NaN);
			double[] spearmanUnSigned = ArrayUtils.doubleArray(2, Double.NaN);
			if (correlB.length > 2) {
				pearsonSigned = Correlation.Pearson(correlB, correlM);
				spearmanSigned = Correlation.Spearman(new double[][] {correlB, correlM});
				pearsonUnSigned = Correlation.Pearson(ArrayUtils.abs(correlB), ArrayUtils.abs(correlM));
				spearmanUnSigned = Correlation.Spearman(new double[][] {ArrayUtils.abs(correlB),
																																ArrayUtils.abs(correlM)});
			}
			BetaCorrelationResult result = new BetaCorrelationResult(comparisonIndex, pearsonSigned[0],
																															 pearsonSigned[1], spearmanSigned[0],
																															 spearmanSigned[1],
																															 pearsonUnSigned[0],
																															 pearsonUnSigned[1],
																															 spearmanUnSigned[0],
																															 spearmanUnSigned[1],
																															 dataBetas.size());
			return result;
		}

	}

	private static class BetaProducer extends AbstractProducer<BetaCorrelationResult[]> {
		private final byte[][] genotypes;
		private final ArrayList<MetaBeta> metaBetas;
		private final ExtProjectDataParser parser;
		private final boolean[] sampleDef;
		private int index;
		private final int max;
		private final double maf;
		private final Logger log;

		public BetaProducer(byte[][] genotypes, boolean[] sampleDef, ArrayList<MetaBeta> metaBetas,
												ExtProjectDataParser parser, int max, double maf, Logger log) {
			super();
			this.genotypes = genotypes;
			this.metaBetas = metaBetas;
			this.parser = parser;
			this.sampleDef = sampleDef;
			this.log = log;
			index = 0;
			this.max = max;
			this.maf = maf;
		}

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return index < parser.getNumericData().length && index < max;
		}

		@Override
		public Callable<BetaCorrelationResult[]> next() {
			double[] data = parser.getNumericData()[index];
			BetaWorker worker = new BetaWorker(genotypes, sampleDef, data, index, metaBetas, maf, log);
			index++;
			return worker;
		}
	}

	private static void analyzeAll(Project proj, String pcFile, String samplesToBuildModels,
																 MarkerDetailSet markerSet, ABLookup abLookup, String dbsnpVCF,
																 String[] namesToQuery, String outpuDir, String[] betas,
																 double[] pvals, double markerCallRate, int maxPCs, int numthreads,
																 String usedInPCFile, int pvalRefineCutoff, double minPval,
																 Logger log) {
		analyze(proj, pcFile, samplesToBuildModels, markerSet, abLookup, dbsnpVCF, namesToQuery,
						outpuDir, Arrays.asList(betas), pvals, markerCallRate, maxPCs, numthreads, usedInPCFile,
						pvalRefineCutoff, minPval, log);

	}

	public enum RUN_TYPES {
		ALL_PC_SAMPS, RACE_IN, RACE_OUT;
	}


	private static FilterGenoResult filterGenos(RUN_TYPES runtype, boolean[] samplesForModels,
																							byte[][] genos, boolean[] sampsPCs,
																							double markerCallRate) {
		boolean[] sampleDef = ArrayUtils.booleanArray(samplesForModels.length, false);
		byte[][] analysisGenos = new byte[genos.length][genos[0].length];
		int markersRemoved = 0;
		int numSamps = 0;
		boolean[] genoDef = ArrayUtils.booleanArray(genos.length, true);

		switch (runtype) {
			case RACE_OUT:// typcially non-whites
				if (ArrayUtils.booleanArraySum(samplesForModels) > 0) {
					for (int j = 0; j < genos.length; j++) {
						int numMissing = 0;
						int numSampsHere = 0;
						for (int j2 = 0; j2 < genos[j].length; j2++) {
							if (!samplesForModels[j2] && sampsPCs[j2]) {
								numSampsHere++;
								sampleDef[j2] = true;
								if (genos[j][j2] < 0) {
									numMissing++;
								}
							}
							analysisGenos[j][j2] = ((!samplesForModels[j2] && sampsPCs[j2]) ? genos[j][j2] : -1);// set
																																																	 // to
																																																	 // missing
						}
						double cr = (double) numMissing / numSampsHere;
						if ((1 - cr) < markerCallRate) {
							markersRemoved++;
							analysisGenos[j] = ArrayUtils.byteArray(analysisGenos[j].length, (byte) -1);
							genoDef[j] = false;
							// log.reportTimeInfo("Removing " +
							// current.get(j).getMarkerRsFormat().getMarkerName() + " for callrate
							// " + markerCallRate + "(" + cr + ")");
						}
						numSamps = numSampsHere;
					}
				}

				break;
			case RACE_IN:// typically whites
				for (int j = 0; j < genos.length; j++) {
					int numMissing = 0;
					int numSampsHere = 0;

					for (int j2 = 0; j2 < genos[j].length; j2++) {
						if (samplesForModels[j2] && sampsPCs[j2]) {
							numSampsHere++;
							sampleDef[j2] = true;
							if (genos[j][j2] < 0) {
								numMissing++;
							}
						}
						analysisGenos[j][j2] = ((samplesForModels[j2] && sampsPCs[j2]) ? genos[j][j2] : -1);// set
																																																// to
																																																// missing
					}
					numSamps = numSampsHere;
					double cr = (double) numMissing / numSamps;
					if ((1 - cr) < markerCallRate) {
						markersRemoved++;
						// log.reportTimeInfo("Removing " +
						// current.get(j).getMarkerRsFormat().getMarkerName() + " for callrate "
						// + markerCallRate + "(" + cr + ")");
						analysisGenos[j] = ArrayUtils.byteArray(analysisGenos[j].length, (byte) -1);
						genoDef[j] = false;

					}
				}
				break;
			case ALL_PC_SAMPS:
				for (int j = 0; j < genos.length; j++) {
					int numMissing = 0;
					int numSampsHere = 0;

					for (int j2 = 0; j2 < genos[j].length; j2++) {
						if (sampsPCs[j2]) {
							numSampsHere++;
							sampleDef[j2] = true;
							if (genos[j][j2] < 0) {
								numMissing++;
							}
						}
						analysisGenos[j][j2] = (sampsPCs[j2] ? genos[j][j2] : -1);// set to
																																			// missing
					}
					numSamps = numSampsHere;
					double cr = (double) numMissing / numSamps;
					if ((1 - cr) < markerCallRate) {
						markersRemoved++;
						// log.reportTimeInfo("Removing " +
						// current.get(j).getMarkerRsFormat().getMarkerName() + " for callrate "
						// + markerCallRate + "(" + cr + ")");
						analysisGenos[j] = ArrayUtils.byteArray(analysisGenos[j].length, (byte) -1);
						genoDef[j] = false;

					}
				}
			default:
				break;
		}
		return new FilterGenoResult(sampleDef, analysisGenos, markersRemoved, numSamps, genoDef);
	}


	private static class FilterGenoResult {
		private boolean[] sampleDef;
		private boolean[] genoDef;
		private byte[][] analysisGenos;
		private int markersRemoved;
		private int numSamps;

		public FilterGenoResult(boolean[] sampleDef, byte[][] analysisGenos, int markersRemoved,
														int numSamps, boolean[] genoDef) {
			super();
			this.sampleDef = sampleDef;
			this.analysisGenos = analysisGenos;
			this.markersRemoved = markersRemoved;
			this.numSamps = numSamps;
			this.genoDef = genoDef;
		}


	}

	private static ArrayList<MetaBeta> getWindowPruned(FilterGenoResult filterGenoResultPrimary,
																										 ArrayList<MetaBeta> metaPrimary, double pval,
																										 Logger log) {
		ArrayList<MetaBeta> filtered = new ArrayList<BetaOptimizer.MetaBeta>();
		if (filterGenoResultPrimary.genoDef.length != metaPrimary.size()) {
			throw new IllegalStateException("Mismatched meta and genotype info, should have "
																			+ filterGenoResultPrimary.genoDef.length + " but got "
																			+ metaPrimary.size());
		}

		for (int i = 0; i < filterGenoResultPrimary.genoDef.length; i++) {
			if (filterGenoResultPrimary.genoDef[i]) {
				filtered.add(metaPrimary.get(i));
			}
		}
		log.reportTimeInfo(filtered.size() + " markers remaining post QC");
		GeneralHitWindowDetector<MetaBeta> hitWindowDetector = new GeneralHitWindowDetector<MetaBeta>(filtered,
																																																	125000,
																																																	pval,
																																																	pval);
		ArrayList<MetaBeta> windows = new ArrayList<BetaOptimizer.MetaBeta>();
		while (hitWindowDetector.hasNext()) {
			List<MetaBeta> index = hitWindowDetector.next().getIndexHits();
			for (MetaBeta in : index) {
				windows.add(in);
			}
		}
		log.reportTimeInfo(windows.size() + " markers remaining post window selection");
		return windows;
	}


	private static void analyze(Project proj, String pcFile, String singleRaceSamples,
															MarkerDetailSet markerSet, ABLookup abLookup, String dbsnpVCF,
															String[] namesToQuery, String outpuDir, List<String> betaFiles,
															double[] pvals, double markerCallRate, int maxPCs, int numthreads,
															String usedInPCFile, int pvalRefineCutoff, double minPval,
															Logger log) {

		String subDir = ext.rootOf(pcFile, false) + SUB_DIR;

		// ArrayList<RScatter> rScatters = new ArrayList<RScatter>();

		for (String betaFile : betaFiles) {
			String bigSummaryOut = outpuDir + ext.rootOf(betaFile) + "_beta_summary.txt";

			if (!Files.exists(bigSummaryOut)) {

				try {
					PrintWriter writer = new PrintWriter(new FileWriter(bigSummaryOut));
					writer.println(ArrayUtils.toStr(BetaCorrelationResult.getHeader()) + "\t"
												 + ArrayUtils.toStr(ArrayUtils.tagOn(BetaCorrelationResult.getHeader(),
																														 "inv_", null))
												 + "\tBetaFile\tMethod\tpvalCutoff\tnumSamples\tmethodKey\tmarkerCallRateThreshold");
					ArrayList<MetaBeta> metaBetas = prep(proj, markerSet, abLookup, dbsnpVCF, namesToQuery,
																							 outpuDir, betaFile, ArrayUtils.max(pvals), log);
					if (pvalRefineCutoff > 0 && metaBetas.size() > pvalRefineCutoff) {
						ArrayList<Double> tmpPvals = new ArrayList<Double>();

						for (double pval : pvals) {
							tmpPvals.add(pval);
						}
						double seed = ArrayUtils.min(pvals);
						ArrayList<MetaBeta> tmp = filter(metaBetas, seed);
						while (tmp.size() > pvalRefineCutoff && seed > minPval) {
							seed = seed / 10;
							tmp = filter(metaBetas, seed);
							if (tmp.size() > pvalRefineCutoff) {
								tmpPvals.add(seed);
							}
						}
						pvals = Doubles.toArray(tmpPvals);
					}

					boolean[] samplesForModels = ArrayUtils.booleanArray(proj.getSamples().length, false);
					String[] pcSamps = HashVec.loadFileToStringArray(usedInPCFile, false, false,
																													 new int[] {0}, false, true, "\t");
					int[] indicesPC = ext.indexLargeFactors(pcSamps, proj.getSamples(), true, proj.getLog(),
																									true, false);
					boolean[] sampsPCs = ArrayUtils.booleanArray(proj.getSamples().length, false);

					for (int i = 0; i < indicesPC.length; i++) {
						sampsPCs[indicesPC[i]] = true;
					}

					if (singleRaceSamples != null && Files.exists(singleRaceSamples)) {
						log.reportTimeInfo("Loading samples from " + singleRaceSamples);
						String[] sampsForMods = HashVec.loadFileToStringArray(singleRaceSamples, false,
																																	new int[] {0}, true);
						log.reportTimeInfo("Loaded " + sampsForMods.length + " samples from "
															 + singleRaceSamples);

						int[] indices = ext.indexLargeFactors(sampsForMods, proj.getSamples(), true,
																									proj.getLog(), true, false);
						for (int i = 0; i < indices.length; i++) {
							samplesForModels[indices[i]] = true;
						}
					} else {
						log.reportTimeWarning("A file of unrelated and single race samples was not provided for optimization");
					}
					for (double pval : pvals) {
						ArrayList<MetaBeta> current = filter(metaBetas, pval);
						if (current.size() > 2) {
							byte[][] genos = loadGenos(proj, markerSet, current);

							for (RUN_TYPES runtype : RUN_TYPES.values()) {
								for (ORDER_TYPE oType : ORDER_TYPE.values()) {
									String file = subDir + MODEL_BUILDER_TYPE.WITH_QC_BUILDERS + "_" + oType + "_"
																+ ITERATION_TYPE.WITHOUT_INDEPS + "_finalSummary.estimates.txt.gz";
									if (Files.exists(file)) {
										FilterGenoResult filterGenoResultPrimary = filterGenos(runtype,
																																					 samplesForModels, genos,
																																					 sampsPCs,
																																					 markerCallRate);

										if (filterGenoResultPrimary.markersRemoved > 0) {
											log.reportTimeWarning("Removed " + filterGenoResultPrimary.markersRemoved
																						+ " markers for pval cut-off " + pval
																						+ " at callrate threshold " + markerCallRate);
										}
										log.reportTimeInfo("numsamps =" + filterGenoResultPrimary.numSamps);
										if (ArrayUtils.booleanArraySum(filterGenoResultPrimary.sampleDef) != filterGenoResultPrimary.numSamps) {
											throw new IllegalArgumentException("Mismatched sample definitions");
										}
										ArrayList<MetaBeta> filtered = getWindowPruned(filterGenoResultPrimary, current,
																																	 pval, log);

										String rootDump = outpuDir + MODEL_BUILDER_TYPE.WITH_QC_BUILDERS + "_" + oType
																			+ "_" + ITERATION_TYPE.WITHOUT_INDEPS + "_" + runtype
																			+ "_pval_" + pval + "_" + ext.rootOf(betaFile)
																			+ "_markersUsed";
										dumpToFile(current, rootDump + "_preFilter.txt", log);
										dumpToFile(filtered, rootDump + "_postFilter.txt", log);

										if (filtered.size() > pvalRefineCutoff) {

											byte[][] genosSecondary = loadGenos(proj, markerSet, filtered);
											FilterGenoResult filterGenoResultSecondary = filterGenos(runtype,
																																							 samplesForModels,
																																							 genosSecondary,
																																							 sampsPCs,
																																							 markerCallRate);

											if (filterGenoResultSecondary.numSamps > 0) {
												ProjectDataParserBuilder builder = new ProjectDataParserBuilder();
												builder.sampleBased(true);
												builder.requireAll(true);
												builder.dataKeyColumnIndex(0);
												builder.treatAllNumeric(true);
												try {
													ExtProjectDataParser parser = builder.build(proj, file);
													parser.determineIndicesFromTitles();
													parser.loadData();
													BetaProducer producer = new BetaProducer(filterGenoResultSecondary.analysisGenos,
																																	 filterGenoResultSecondary.sampleDef,
																																	 current, parser, maxPCs, 0, log);
													WorkerTrain<BetaCorrelationResult[]> train = new WorkerTrain<BetaOptimizer.BetaCorrelationResult[]>(producer,
																																																															numthreads,
																																																															10,
																																																															log);
													while (train.hasNext()) {
														BetaCorrelationResult[] results = train.next();
														String method = "_" + oType;
														switch (runtype) {
															case ALL_PC_SAMPS:
																method = "AllPCASamps_" + method;
																break;
															case RACE_IN:
																method = ext.rootOf(singleRaceSamples) + method;
																break;
															case RACE_OUT:
																method = "not_" + ext.rootOf(singleRaceSamples) + method;
																break;
															default:
																break;
														}
														writer.println(ArrayUtils.toStr(results[0].getSummary()) + "\t"
																					 + ArrayUtils.toStr(results[1].getSummary()) + "\t"
																					 + betaFile + "\t" + method + "\t" + pval + "\t"
																					 + filterGenoResultSecondary.numSamps + "\t" + method
																					 + "_" + pval + "\t" + markerCallRate);
													}
												} catch (FileNotFoundException e) {
													log.reportException(e);
												}
											}

										}
									} else {
										log.reportTimeWarning("Expecting file " + file + ", and did not find it");
									}
								}
							}
						}
						// }
					}
					writer.close();

				} catch (Exception e) {
					log.reportError("Error writing to " + bigSummaryOut);
					log.reportException(e);
				}
			}
			String rootOutBetas = outpuDir + ext.rootOf(betaFile) + "_summaryPlots";
			String rootOutInvBetas = outpuDir + ext.rootOf(betaFile) + "_Inv_summaryPlots";

			ArrayList<RScatter> rScattersBetas = new ArrayList<RScatter>();
			ArrayList<RScatter> rScattersInvBetas = new ArrayList<RScatter>();

			String[] summHeader = Files.getHeaderOfFile(bigSummaryOut, log);
			Hashtable<String, Hashtable<String, Vector<String>>> info = HashVec.loadFileToHashHashVec(bigSummaryOut,
																																																ext.indexOfStr("pvalCutoff",
																																																							 summHeader),
																																																ext.indexOfStr("Method",
																																																							 summHeader),
																																																new int[] {ext.indexOfStr("numberOfMarkers",
																																																													summHeader),
																																																					 ext.indexOfStr("numSamples",
																																																													summHeader)},
																																																false,
																																																false);

			for (double pval : pvals) {
				String pvalKey = pval + "".replaceAll("\\.", "_") + ".pval";
				VertLine[] vertLines = new VertLine[] {new VertLine(15)};
				Restrictions rPval = new Restrictions(new String[] {"pvalCutoff"}, new double[] {pval},
																							new String[] {"=="}, null);
				Hashtable<String, Vector<String>> currentInfo = info.get(pval + "");

				String[][] altLegends = null;
				if (currentInfo != null) {
					String[] tmpMethods = new String[currentInfo.size()];
					Map<String, String> legendMap = new HashMap<String, String>();
					int ind = 0;
					for (String method : currentInfo.keySet()) {
						String[] inf = currentInfo.get(method).get(0).split("\t");
						String legend = method + " (m=" + inf[0] + ", s=" + inf[1] + ")";

						tmpMethods[ind] = method;
						legendMap.put(method, legend);


						ind++;
					}
					Arrays.sort(tmpMethods, new SciStringComparator());
					altLegends = new String[2][currentInfo.size()];

					for (int i = 0; i < tmpMethods.length; i++) {
						altLegends[0][i] = tmpMethods[i];
						altLegends[1][i] = legendMap.get(tmpMethods[i]);
					}
				}

				String methodProotUnsigned = rootOutBetas + "pvalue_method_pearsonSigned" + pvalKey;

				RScatter methodPU = new RScatter(bigSummaryOut, methodProotUnsigned + ".rscript",
																				 ext.rootOf(methodProotUnsigned),
																				 methodProotUnsigned + ".jpeg", "PC",
																				 new String[] {"correl_pearsonSigned"}, "Method",
																				 SCATTER_TYPE.POINT, log);
				methodPU.setOverWriteExisting(true);
				methodPU.setRestrictions(new Restrictions[] {rPval});// rNat
				methodPU.setxRange(new double[] {0, 150});
				methodPU.setyRange(new double[] {-1, 1});
				methodPU.setyLabel("correlation (pearson) of betas with p < " + pval);
				methodPU.setTitle("ABS Betas from " + ext.rootOf(betaFile) + "\nMarker Call rate >= "
													+ markerCallRate);
				methodPU.setxLabel("PC");
				methodPU.setVertLines(vertLines);
				methodPU.setAltLegendTitles(altLegends);
				methodPU.setLegendTitle("Method");
				methodPU.execute();
				rScattersBetas.add(methodPU);

				String methodProotsignednat = rootOutBetas + "pvalue_method_spearmanSignedNat" + pvalKey;

				RScatter methodPSNat = new RScatter(bigSummaryOut, methodProotsignednat + ".rscript",
																						ext.rootOf(methodProotsignednat),
																						methodProotsignednat + ".jpeg", "PC",
																						new String[] {"correl_spearmanSigned"}, "Method",
																						SCATTER_TYPE.POINT, log);
				methodPSNat.setOverWriteExisting(true);
				methodPSNat.setRestrictions(new Restrictions[] {rPval});// , rNat
				methodPSNat.setxRange(new double[] {0, 150});
				methodPSNat.setyRange(new double[] {-1, 1});
				methodPSNat.setyLabel("correlation (spearman) of betas with p < " + pval);
				methodPSNat.setTitle("Betas from " + ext.rootOf(betaFile) + "\nMarker Call rate >= "
														 + markerCallRate);
				methodPSNat.setxLabel("PC");
				methodPSNat.setVertLines(vertLines);
				methodPSNat.execute();
				methodPSNat.setAltLegendTitles(altLegends);
				methodPSNat.setLegendTitle("Method");
				rScattersBetas.add(methodPSNat);

				String methodProotUnsignedInv = rootOutBetas + "pvalue_method_pearsonSigned_Inv" + pvalKey;

				RScatter methodPUInv = new RScatter(bigSummaryOut, methodProotUnsignedInv + ".rscript",
																						ext.rootOf(methodProotUnsignedInv),
																						methodProotUnsignedInv + ".jpeg", "PC",
																						new String[] {"inv_correl_pearsonSigned"}, "Method",
																						SCATTER_TYPE.POINT, log);
				methodPUInv.setOverWriteExisting(true);
				methodPUInv.setRestrictions(new Restrictions[] {rPval});// rNat
				methodPUInv.setxRange(new double[] {0, 150});
				methodPUInv.setyRange(new double[] {-1, 1});
				methodPUInv.setyLabel("correlation (pearson) of betas with p < " + pval);
				methodPUInv.setTitle("ABS Betas from " + ext.rootOf(betaFile) + "\nMarker Call rate >= "
														 + markerCallRate + "\nInverse normalized mtDNA CN estimate");
				methodPUInv.setxLabel("PC");
				methodPUInv.setVertLines(vertLines);
				methodPUInv.setAltLegendTitles(altLegends);
				methodPUInv.setLegendTitle("Method");
				methodPUInv.execute();
				rScattersInvBetas.add(methodPUInv);

				String methodProotsignednatInv = rootOutBetas + "pvalue_method_spearmanSignedNat_Inv"
																				 + pvalKey;

				RScatter methodPSNatInv = new RScatter(bigSummaryOut, methodProotsignednatInv + ".rscript",
																							 ext.rootOf(methodProotsignednatInv),
																							 methodProotsignednatInv + ".jpeg", "PC",
																							 new String[] {"inv_correl_spearmanSigned"}, "Method",
																							 SCATTER_TYPE.POINT, log);
				methodPSNatInv.setOverWriteExisting(true);
				methodPSNatInv.setRestrictions(new Restrictions[] {rPval});// , rNat
				methodPSNatInv.setxRange(new double[] {0, 150});
				methodPSNatInv.setyRange(new double[] {-1, 1});
				methodPSNatInv.setyLabel("correlation (spearman) of betas with p < " + pval);
				methodPSNatInv.setTitle("Betas from " + ext.rootOf(betaFile) + "\nMarker Call rate >= "
																+ markerCallRate + "\nInverse normalized mtDNA CN estimate");
				methodPSNatInv.setxLabel("PC");
				methodPSNatInv.setVertLines(vertLines);
				methodPSNatInv.execute();
				methodPSNatInv.setAltLegendTitles(altLegends);
				methodPSNatInv.setLegendTitle("Method");
				rScattersInvBetas.add(methodPSNatInv);

			}
			RScatters rscScattersBetas = new RScatters(rScattersBetas, rootOutBetas + ".rscript",
																								 rootOutBetas + ".pdf",
																								 COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1,
																								 PLOT_DEVICE.PDF, log);
			rscScattersBetas.execute();

			RScatters rscScattersInvBetas = new RScatters(rScattersInvBetas, rootOutInvBetas + ".rscript",
																										rootOutInvBetas + ".pdf",
																										COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1,
																										PLOT_DEVICE.PDF, log);
			rscScattersInvBetas.execute();

		}
		// String rootOut = outpuDir + "finalPlots";
		// RScatters rscScatters = new RScatters(rScatters, rootOut + ".rscript", rootOut + ".pdf",
		// COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, log);
		// rscScatters.execute();

	}

	private static ArrayList<MetaBeta> filter(ArrayList<MetaBeta> metaBetas, double pval) {
		ArrayList<MetaBeta> filt = new ArrayList<MetaBeta>();
		for (MetaBeta m : metaBetas) {
			if (m.getPval() < pval) {
				filt.add(m);
			}
		}
		return filt;
	}

	private static byte[][] loadGenos(Project proj, MarkerDetailSet markerSet,
																		ArrayList<MetaBeta> metaBetas) {
		byte[][] genos = new byte[metaBetas.size()][];
		String[] markerNames = new String[metaBetas.size()];

		for (int i = 0; i < genos.length; i++) {
			markerNames[i] = metaBetas.get(i).getMarkerRsFormat().getMarkerName();
		}
		MDL mdl = new MDL(proj, markerSet, markerNames, 2, 100);
		int index = 0;
		proj.getLog().reportTimeInfo("Loading genotypes...");
		while (mdl.hasNext()) {
			genos[index] = mdl.next().getAbGenotypes();
			index++;
		}
		return genos;
	}

	private static void dumpToFile(ArrayList<MetaBeta> metaBetas, String output, Logger log) {
		StringBuilder builder = new StringBuilder();
		builder.append("rsID\tbeta\tp\tchr\tpos\tmarkerName\n");
		for (MetaBeta metaBeta : metaBetas) {
			builder.append(metaBeta.getMarkerRsFormat().rs + "\t" + metaBeta.beta + "\t"
										 + metaBeta.getPval() + "\t" + metaBeta.getChr() + "\t" + metaBeta.getPos()
										 + "\t" + metaBeta.getMarkerRsFormat().markerName + "\n");

		}
		Files.write(builder.toString(), output);
	}

	private static ArrayList<MetaBeta> prep(Project proj, MarkerSetInfo markerSet, ABLookup abLookup,
																					String dbsnpVCF, String[] namesToQuery, String outpuDir,
																					String betaFile, double minPval, Logger log) {
		new File(outpuDir).mkdirs();
		String outSer = outpuDir + "rsIdLookup.ser";
		ArrayList<MarkerRsFormat> markerRsFormats = null;
		if (Files.exists(outSer)) {
			try {
				log.reportTimeInfo("Trying to load " + outSer);
				markerRsFormats = MarkerRsFormat.readSerial(outSer, log);
				log.reportTimeInfo("Loaded " + outSer);
			} catch (Exception e) {

			}
		}
		if (markerRsFormats == null) {
			markerRsFormats = mapToRsIds(proj, abLookup, dbsnpVCF, namesToQuery, outSer, log);
		}
		ArrayList<MetaBeta> metaBetas = loadBetas(markerRsFormats, betaFile, minPval, markerSet, log);
		log.reportTimeInfo("Loaded " + metaBetas.size()
											 + " valid rsIds, having valid betas and  pval less than " + minPval);
		return metaBetas;

	}

	private static boolean useConfig(CONFIG config) {
		switch (config) {

			case STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:
			case STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
			case STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
			case STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
				return true;
			case STRAND_CONFIG_DIFFERENT_ALLELES:
			case STRAND_CONFIG_BOTH_NULL:
			case STRAND_CONFIG_SPECIAL_CASE:
			case STRAND_CONFIG_UNKNOWN:
			case STRAND_CONFIG_AMBIGOUS:
				return false;

			default:
				return false;

		}
	}

	private static ArrayList<MetaBeta> loadBetas(ArrayList<MarkerRsFormat> markerRsFormats,
																							 String betaFile, double minPval,
																							 MarkerSetInfo markerSet,
																							 Logger log) {
		String[] header = Files.getHeaderOfFile(betaFile, log);
		int[] indices = ext.indexFactors(BETA_HEADER, header, false, false);
		int ambi = 0;

		if (ArrayUtils.countIf(indices, -1) > 0) {
			log.reportError("Did not detect proper header in " + betaFile + ", requires "
											+ ArrayUtils.toStr(BETA_HEADER));
			return null;
		} else {
			ArrayList<MetaBeta> metaBetas = new ArrayList<BetaOptimizer.MetaBeta>();
			Hashtable<String, Integer> index = new Hashtable<String, Integer>();
			for (int i = 0; i < markerRsFormats.size(); i++) {
				index.put(markerRsFormats.get(i).getRs(), i);
			}
			HashSet<String> added = new HashSet<String>();
			try {
				BufferedReader reader = Files.getAppropriateReader(betaFile);
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					String rsId = line[indices[0]];
					if (index.containsKey(rsId)) {
						try {
							double beta = Double.parseDouble(line[indices[3]]);
							double betaSE = Double.parseDouble(line[indices[4]]);
							double p = Double.parseDouble(line[indices[5]]);

							if (Java6Helper.isFinite(beta) && Java6Helper.isFinite(p) && p < minPval) {
								MarkerRsFormat current = markerRsFormats.get(index.get(rsId));
								String[] betaAlleles = new String[] {line[indices[1]].toUpperCase(),
																										 line[indices[2]].toUpperCase()};
								String[] markerAlleles = new String[] {current.getMarkerAlleles()[0],
																											 current.getMarkerAlleles()[1]};
								CONFIG strandConfig = StrandOps.determineStrandConfig(markerAlleles, betaAlleles);
								if (strandConfig == CONFIG.STRAND_CONFIG_AMBIGOUS) {
									ambi++;
								}
								if (useConfig(strandConfig) && current.isValidMatch()) {

									switch (strandConfig) {
										case STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:
										case STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
											beta *= -1;
											break;
										case STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
										case STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
											break;
										default:
											throw new IllegalArgumentException("Should not have strand config "
																												 + strandConfig + " in current");
									}
									if (!added.contains(current.getMarkerName())) {
										Segment markerSeg = new Segment(markerSet.getChrs()[current.getProjectIndex()],
																										markerSet.getPositions()[current.getProjectIndex()],
																										markerSet.getPositions()[current.getProjectIndex()] + 1);
										if (!markerSeg.overlaps(HG_19_HLA)) {
											MetaBeta me = new MetaBeta(current,
																								 markerSet.getChrs()[current.getProjectIndex()],
																								 markerSet.getPositions()[current.getProjectIndex()],
																								 beta, betaSE, p);
											metaBetas.add(me);
										}
									}
									added.add(current.getMarkerName());
								}
							} else if (!Java6Helper.isFinite(beta) || !Java6Helper.isFinite(p)) {
								log.reportTimeWarning("Invalid number on line " + ArrayUtils.toStr(line));
								log.reportTimeWarning(line[indices[3]] + "\t" + line[indices[4]]);
							}
						} catch (NumberFormatException nfe) {
							log.reportTimeWarning("Invalid number on line " + ArrayUtils.toStr(line));
						}
					}

				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + betaFile + "\" not found in current directory");
				return null;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + betaFile + "\"");
				return null;
			}
			if (ambi > 0) {
				log.reportTimeWarning("Removed " + ambi + " ambiguous markers");
			}
			return metaBetas;

		}

	}

	public enum SITE_TYPE {
		BIALLELIC, TRIALLELIC, UNKNOWN;
	}

	/**
	 * Map a projects markers to rsIds by chr,pos,and alleles
	 *
	 */
	public static ArrayList<MarkerRsFormat> mapToRsIds(Project proj, ABLookup abLookup,
																										 String dbsnpVCF, String[] namesToQuery,
																										 String outSer, Logger log) {
		MarkerSetInfo markerSet = proj.getMarkerSet();
		String[] markerNames = markerSet.getMarkerNames();
		int[] indices = ext.indexLargeFactors(namesToQuery, markerNames, true, log, true, false);
		int[] posIndices = indices;
		if (proj.GENOME_BUILD_VERSION.getValue() != GENOME_BUILD.HG19) {
			if (proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6
					|| proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6_CN) {
				proj.getLog().reportTimeInfo("Attempting to use " + GENOME_BUILD.HG19.getBuild()
																		 + " positions for rsID lookup");
				Resource affyhg19 = Resources.affy(log).genome(GENOME_BUILD.HG19).getMarkerPositions();
				String tmpSer = ext.rootOf(outSer, false) + "hg19.positions.ser";

				if (!Files.exists(tmpSer)) {
					Markers.orderMarkers(null, affyhg19.get(), tmpSer, proj.getLog());

				}
				markerSet = MarkerSet.load(tmpSer, false);
				posIndices = ext.indexLargeFactors(namesToQuery, markerSet.getMarkerNames(), true, log,
																					 true, false);

			} else {
				throw new IllegalArgumentException("Genome version must be "
																					 + GENOME_BUILD.HG19.getBuild());
			}
		}

		Segment[] segs = new Segment[indices.length];
		for (int i = 0; i < posIndices.length; i++) {
			if (posIndices[i] >= 0) {
				segs[i] = new Segment(markerSet.getChrs()[posIndices[i]],
															markerSet.getPositions()[posIndices[i]],
															markerSet.getPositions()[posIndices[i]] + 1);
			} else {
				segs[i] = new Segment((byte) 0, 0, 0 + 1);
			}
		}

		VCFFileReader reader = new VCFFileReader(new File(dbsnpVCF), true);
		ArrayList<MarkerRsFormat> markerRsFormats = new ArrayList<MarkerRsFormat>();
		String outTxt = ext.rootOf(outSer, false) + ".txt";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outTxt));
			writer.println("markerName\trsID\tref\talt\tA\tB");
			log.reportTimeInfo("Attempting to look up " + namesToQuery.length + " markers as rsIds from "
												 + dbsnpVCF);
			for (int i = 0; i < segs.length; i++) {
				Segment current = segs[i];
				String[] allelesMarker = new String[] {abLookup.getLookup()[indices[i]][0]
																							 + "".toUpperCase(),
																							 abLookup
																											 .getLookup()[indices[i]][1]
																																	 + "".toUpperCase()};
				if (i % 10000 == 0 && i > 0) {
					log.reportTimeInfo("queried " + (i + 1) + " markers");
				}
				MarkerRsFormat markerRsFormat = new MarkerRsFormat(namesToQuery[i], current.getStart(),
																													 indices[i], "NA",
																													 new String[] {"N", "N"}, allelesMarker,
																													 CONFIG.STRAND_CONFIG_UNKNOWN,
																													 SITE_TYPE.UNKNOWN);

				if (ArrayUtils.countIf(allelesMarker, "N") == 0) {
					CloseableIterator<VariantContext> vcIter = reader.query(Positions.getChromosomeUCSC(current.getChr(),
																																															false),
																																	current.getStart() - 2,
																																	current.getStop() + 2);
					boolean foundGood = false;
					while (!foundGood && vcIter.hasNext()) {
						VariantContext vc = vcIter.next();

						if (!vc.isPointEvent() || current.getStart() == vc.getStart()
								|| Integer.parseInt(VCOps.getAnnotationsFor(new String[] {"RSPOS"}, vc,
																														"-1")[0]) == current.getStart()) {
							String[][] allelesVC = new String[][] {{"N", "N"}};
							if (vc.isPointEvent()) {
								allelesVC = new String[vc.getAlternateAlleles().size()][];// tri allele possibility
								for (int j = 0; j < allelesVC.length; j++) {
									allelesVC[j] = new String[] {vc.getReference().getBaseString(),
																							 vc.getAlternateAllele(j).getBaseString()};
								}
							} else if (vc.isIndel()) {
								if (vc.getReference().getBaseString().length() > vc.getAlternateAllele(0)
																																	 .length()) {
									allelesVC = new String[][] {{"I", "D"}};
								} else {
									allelesVC = new String[][] {{"D", "I"}};
								}
							}
							for (String[] element : allelesVC) {
								if (!foundGood) {
									String[] tmpMarkerAlleles = new String[] {allelesMarker[0], allelesMarker[1]};
									String[] tmpVcAllel = new String[] {element[0], element[1]};

									CONFIG config = StrandOps.determineStrandConfig(tmpMarkerAlleles, tmpVcAllel);
									markerRsFormat = new MarkerRsFormat(namesToQuery[i], current.getStart(),
																											indices[i], vc.getID(), element,
																											allelesMarker, config,
																											vc.isBiallelic() ? SITE_TYPE.BIALLELIC
																																			 : SITE_TYPE.TRIALLELIC);
									switch (config) {

										case STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:
										case STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
										case STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
										case STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
											markerRsFormat.setValidMatch(true);
											foundGood = true;
											break;
										case STRAND_CONFIG_BOTH_NULL:
										case STRAND_CONFIG_DIFFERENT_ALLELES:
										case STRAND_CONFIG_SPECIAL_CASE:
										default:
											break;
									}
								}
							}
						}

					}
				}
				markerRsFormats.add(markerRsFormat);
				writer.println(markerRsFormat.getMarkerName() + "\t" + markerRsFormat.getRs() + "\t"
											 + ArrayUtils.toStr(markerRsFormat.getDbSnpAlleles()) + "\t"
											 + ArrayUtils.toStr(markerRsFormat.getMarkerAlleles()));
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outTxt);
			log.reportException(e);
		}
		reader.close();
		MarkerRsFormat.writeSerial(markerRsFormats, outSer);
		return markerRsFormats;
	}

	public static class MarkerRsFormat implements Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 1L;
		private final int projectIndex;
		private final String markerName;
		private final String rs;
		private final String[] dbSnpAlleles;
		private final String[] markerAlleles;
		private final int posMarker;
		private final CONFIG config;
		private boolean validMatch;
		private final SITE_TYPE type;

		private MarkerRsFormat(String markerName, int posMarker, int projectIndex, String rs,
													 String[] refs, String[] markers, CONFIG config, SITE_TYPE type) {
			super();
			this.markerName = markerName;
			this.posMarker = posMarker;
			this.projectIndex = projectIndex;
			this.rs = rs;
			dbSnpAlleles = refs;
			markerAlleles = markers;
			this.config = config;
			validMatch = false;
			this.type = type;
		}

		public int getPosMarker() {
			return posMarker;
		}

		public int getProjectIndex() {
			return projectIndex;
		}

		public boolean flippedStrand() {
			return config == CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND
						 || config == CONFIG.STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND;
		}

		public SITE_TYPE getType() {
			return type;
		}

		public String getRs() {
			return rs;
		}

		public CONFIG getConfig() {
			return config;
		}

		private void setValidMatch(boolean validMatch) {
			this.validMatch = validMatch;
		}

		public boolean isValidMatch() {
			return validMatch;
		}

		public String getMarkerName() {
			return markerName;
		}

		public String[] getDbSnpAlleles() {
			return dbSnpAlleles;
		}

		public String[] getMarkerAlleles() {
			return markerAlleles;
		}

		boolean flipBetas() {
			if (!validMatch) {
				throw new IllegalArgumentException("Did not have valid rs id match, this method should not be used");
			}
			return config == CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND
						 || config == CONFIG.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
		}

		public static void writeSerial(ArrayList<MarkerRsFormat> markerRsFormats, String fileName) {
			SerializedFiles.writeSerial(markerRsFormats, fileName, true);
		}

		@SuppressWarnings("unchecked")
		public static ArrayList<MarkerRsFormat> readSerial(String fileName, Logger log) {
			return (ArrayList<MarkerRsFormat>) SerializedFiles.readSerial(fileName, false, log, false,
																																		true);
		}

	}

	/**
	 * Storing the betas of variants from a meta-analysis
	 *
	 */
	private static class MetaBeta extends BasicHit {
		private final MarkerRsFormat markerRsFormat;
		private final double beta;
		private final double betaSE; // beta standard error


		public MetaBeta(MarkerRsFormat markerRsFormat, byte chr, int pos, double beta, double betaSE,
										double p) {
			super(markerRsFormat.getMarkerName(), chr, pos, p);
			this.markerRsFormat = markerRsFormat;
			this.beta = beta;
			this.betaSE = betaSE;
		}

		private MarkerRsFormat getMarkerRsFormat() {
			return markerRsFormat;
		}

		private double getBeta() {
			return beta;
		}
	}

	public static void optimize(Project proj, String pcFile, String outDir, String betaLoc,
															String unRelatedFile, String pcSamps, double[] pvals, int maxPCs,
															double markerCallRate, int pvalRefineCutoff, double minPval,
															int numthreads) throws IllegalStateException {
		new File(outDir).mkdirs();

		proj.AB_LOOKUP_FILENAME.setValue(outDir + "AB_LookupBeta.dat");
		ABLookup abLookup = new ABLookup();
		if (!Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
			if (Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
				Files.copyFile(proj.AB_LOOKUP_FILENAME.getValue(), outDir + "AB_LookupBeta.dat");
			} else {
				if (proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6
						|| proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6_CN) {
					Resource abAffy = Resources.affy(proj.getLog()).getABLookup();
					if (abAffy.isAvailable()) {
						String tmpAB = abAffy.get();
						Files.copyFile(tmpAB, proj.AB_LOOKUP_FILENAME.getValue());
					} else {
						throw new IllegalStateException("Could not retrieve required AB lookup for affymetrix array, halting");
					}
				}
				if (!Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
					abLookup.parseFromAnnotationVCF(proj);
					abLookup.writeToFile(proj.AB_LOOKUP_FILENAME.getValue(), proj.getLog());
				}
				if (!Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
					throw new IllegalArgumentException("Could not create required AB Lookup file, halting");
				}
			}
		}
		MarkerDetailSet markerSet = proj.getMarkerSet();
		abLookup = new ABLookup(markerSet.getMarkerNames(), proj.AB_LOOKUP_FILENAME.getValue(), true,
														true, proj.getLog());
		Resource dbsnp = Resources.genome(GENOME_BUILD.HG19, proj.getLog()).getDBSNP();// TODO, need hg
																																									 // 18
		// db snp
		if (!dbsnp.isAvailable()) {
			throw new IllegalStateException("Could not retrieve required dbSNP vcf");
		}
		String[] betaFiles = null;
		if (Files.isDirectory(betaLoc)) {
			proj.getLog().reportTimeInfo("Searching " + betaLoc + " for files ending with " + ".beta");
			betaFiles = Files.list(betaLoc, "", ".beta", true, false, true);
			proj.getLog().reportTimeInfo("found " + betaFiles.length + " beta files in " + betaLoc);
		} else {
			betaFiles = new String[] {betaLoc};
		}
		analyzeAll(proj, pcFile, unRelatedFile, markerSet, abLookup, dbsnp.get(),
							 proj.getNonCNMarkers(), outDir, betaFiles, pvals, markerCallRate, maxPCs, numthreads,
							 pcSamps, pvalRefineCutoff, minPval, proj.getLog());
	}

	public static void main(String[] args) {
		if (args.length == 2 && args[1].equals("test")) {
			maisn(args);
		} else {
			int numArgs = args.length;
			String filename = "BetaOptimizer.dat";
			String out = "betaOpt/";
			double[] pvals = new double[] {.05, .01, .001, .0001};
			double markerCallRate = .96;
			String unRelatedFile = null;
			String pcFile = null;
			int numthreads = 24;
			int maxPCs = 120;
			String pcSamps = null;
			String betaDir = null;
			String usage = "\n" + "one.JL.BetaOptimizer requires 0-1 arguments\n"
										 + "   (1) project filename (i.e. proj=" + filename + " (default))\n"
										 + "   (2) output directory (relative to project directory) (i.e. out=" + out
										 + " (default))\n"
										 + "   (3) pcFile(relative to project directory) (i.e. pcFile= (no default))\n"
										 + "   (4) full path to a file of unrelated, and single race samples  (i.e. unrelated= (no default))\n"
										 + "   (5) maximum number of pcs to optimze to (i.e. maxPCs=" + maxPCs
										 + " ( default))\n"
										 + "   (6) comma delimited list of pvalue thresholds (i.e. pvals="
										 + ArrayUtils.toStr(ArrayUtils.toStringArray(pvals), ",") + " ( default))\n"
										 + PSF.Ext.getNumThreadsCommand(7, numthreads)
										 + "   (8) full path to a file of samples used to generate the pcs (i.e. pcSamps= (no default))\n"
										 + "   (9) full path a directory of .beta files (i.e. betaDir= (no default))\n"
										 + "   (10) call rate for markers (i.e. cr=" + markerCallRate + " (default))\n"
										 +

										 "";

			for (String arg : args) {
				if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
					System.err.println(usage);
					System.exit(1);
				} else if (arg.startsWith("proj=")) {
					filename = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("pcFile=")) {
					pcFile = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("unrelated=")) {
					unRelatedFile = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("betaDir=")) {
					betaDir = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("out=")) {
					out = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith("pvals=")) {
					pvals = ArrayUtils.toDoubleArray(arg.split("=")[1].split(","));
					numArgs--;
				} else if (arg.startsWith("maxPCs=")) {
					maxPCs = ext.parseIntArg(arg);
					numArgs--;
				} else if (arg.startsWith("cr=")) {
					markerCallRate = ext.parseDoubleArg(arg);
					numArgs--;
				} else if (arg.startsWith("pcSamps=")) {
					pcSamps = arg.split("=")[1];
					numArgs--;
				} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
					numthreads = ext.parseIntArg(arg);
					numArgs--;
				} else {
					System.err.println("Error - invalid argument: " + arg);
				}
			}
			if (numArgs != 0) {
				System.err.println(usage);
				System.exit(1);
			}
			try {

				Project proj = new Project(filename, false);
				optimize(proj, pcFile, proj.PROJECT_DIRECTORY.getValue() + out, betaDir, unRelatedFile,
								 pcSamps, pvals, maxPCs, markerCallRate, 25, -1, numthreads);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	public static void maisn(String[] args) {
		Project proj = new Project("/home/pankrat2/lanej/projects/Aric_gw6.properties", false);
		String out = proj.PROJECT_DIRECTORY.getValue() + "betaOpt/";
		new File(out).mkdirs();
		proj.AB_LOOKUP_FILENAME.setValue(out + "AB_LookupBeta.dat");
		ABLookup abLookup = new ABLookup();
		if (!Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
			abLookup.parseFromAnnotationVCF(proj);
			abLookup.writeToFile(proj.AB_LOOKUP_FILENAME.getValue(), proj.getLog());
		}

		MarkerDetailSet markerSet = proj.getMarkerSet();
		abLookup = new ABLookup(markerSet.getMarkerNames(), proj.AB_LOOKUP_FILENAME.getValue(), true,
														true, proj.getLog());
		Resource dbsnp = Resources.genome(GENOME_BUILD.HG19, proj.getLog()).getDBSNP();
		String betaFileDir = "/home/pankrat2/shared/MitoPipeLineResources/betas/" + args[0] + "/";

		String[] betaFiles = Files.list(betaFileDir, "", ".beta", true, false, true);
		proj.getLog().reportTimeInfo("found " + betaFiles.length + " beta files in " + betaFileDir);
		String pcFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/ohw_ws_20_ALL1000PCs_gc_corrected_OnTheFly_SampLRR_Recomp_LRR_035CR_096.PCs.extrapolated.txt";
		String toUseFile = out + "Whites.txt";
		String usedInPCFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/ohw_ws_20_ALL1000PCs_gc_corrected_OnTheFly_SampLRR_Recomp_LRR_035CR_096.samples.USED_PC.txt";
		double[] pvals = new double[] {.05, .01, .001, .0001};
		double markerCallRate = .96;
		int numthreads = 24;
		proj.getLog().reportTimeInfo("Using " + numthreads + " of "
																 + Runtime.getRuntime().availableProcessors() + " available cores");
		int maxPCs = 120;
		analyzeAll(proj, pcFile, toUseFile, markerSet, abLookup, dbsnp.get(), proj.getNonCNMarkers(),
							 out, betaFiles, pvals, markerCallRate, maxPCs, numthreads, usedInPCFile, 25, -1,
							 proj.getLog());
	}
}
