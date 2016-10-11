package org.genvisis.cnv.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.MosaicismQuant.ComputeParams.MosaicParamsBuilder;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.hmm.CNVCaller;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Numbers;
import org.genvisis.common.PSF;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.Sort;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.mining.Calcfc;
import org.genvisis.mining.Cobyla;
import org.genvisis.mining.CobylaExitStatus;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.stats.Rscript.COLUMNS_MULTIPLOT;
import org.genvisis.stats.Rscript.PLOT_DEVICE;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.RScatters;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import be.ac.ulg.montefiore.run.distributions.GaussianMixtureDistribution;

/**
 * @author lane0212 Class to quantify mosiacism at a particular locus, based on procedure from
 *         http://www.ncbi.nlm.nih.gov/pubmed/22277120
 */
public class MosaicismQuant implements Calcfc {
	private static final double MIN_BAF = 0.15;
	private static final double MAX_BAF = 0.85;
	private final Project proj;
	private final SampleMosiac sampleMosiac;
	private final Segment evalSegment;
	private final MarkerSet markerSet;
	private final int[][] indicesByChr;
	private final MOSAIC_TYPE type;

	// private ComputeParams computeParams;

	public MosaicismQuant(Project proj, SampleMosiac sampleMosiac, MOSAIC_TYPE type,
												ComputeParams computeParams, Segment evalSegment, MarkerSet markerSet,
												int[][] indicesByChr) {
		super();
		this.proj = proj;
		this.sampleMosiac = new SampleMosiac(sampleMosiac);
		this.evalSegment = evalSegment;
		this.markerSet = markerSet;
		this.indicesByChr = indicesByChr;
		this.type = type;
		// this.computeParams = computeParams;
		// this.proj.getLog().reportTimeInfo("MOSAIC TYPE: " + type);
	}

	public enum MOSAIC_TYPE {
														MONOSOMY_DISOMY, TRISOMY_DISOMY;
	}

	public Project getProj() {
		return proj;
	}

	@Override
	public double Compute(int n, int m, double[] x, double[] con) {
		// if (x[0] - computeParams.getMinF() < 0) {
		// con[0] = -100000;
		// }
		// if (computeParams.getMaxF() - x[0] < 0) {
		// con[1] = -100000;
		// }
		// if (x[1] - computeParams.getMinShift() < 0) {
		// con[2] = -100000;
		// }
		// if (computeParams.getMaxShift() - x[1] < 0) {
		// con[3] = -100000;
		// }
		// if (x[2] - computeParams.getMinNullD() < 0) {
		// con[4] = 0;
		// }
		// if (computeParams.getMaxNullD() - x[2] < 0) {
		// con[5] = 0;
		// }
		return getValue(x);
	}

	private double getValue(double[] x) {
		double f = x[0] / 100;
		double shift = x[1];
		double delta = getDelta(f, type);
		int count = sampleMosiac.getCdf().getVals().length;
		int halfCount = getHalfCount(count);
		double halfSigma = getHalfSigma(halfCount, shift);
		double[] currentCase = sampleMosiac.getCdf().getValsInOrder();
		double[] currentControl = sampleMosiac.getSmoothedControl().getValsInOrder();
		double diffAverage = Array.mean(currentCase) - Array.mean(currentControl);
		double[] shiftControl = getShift(currentControl, diffAverage, delta, halfSigma);
		int[] shiftControlRank = Sort.quicksort(shiftControl);
		double[] resids = getResid(shiftControl, shiftControlRank, currentCase);
		return getSumResids(resids);
	}

	/**
	 *
	 * Loads current sample and controls, develops cdfs
	 */
	public void prepareMosiacQuant(int numThreads, double minBaf, double maxBaf) {
		sampleMosiac.load(numThreads);
		sampleMosiac.developCDFs(markerSet, indicesByChr, evalSegment, minBaf, maxBaf);

	}

	private static double getSumResids(double[] resids) {
		double sum = 0;
		for (double resid : resids) {
			sum += Math.sqrt(Math.pow(resid, 2));
		}
		return sum;

	}

	private static double[] getResid(	double[] shiftControl, int[] shiftControlRank,
																		double[] currentCase) {
		double[] resid = new double[shiftControl.length];
		for (int i = 0; i < resid.length; i++) {
			resid[i] = currentCase[i] - shiftControl[shiftControlRank[i]];
		}
		return resid;
	}

	private static double[] getShift(	double[] currentControl, double diffAverage, double delta,
																		double halfSigma) {
		double[] shift = new double[currentControl.length];
		for (int i = 0; i < currentControl.length; i++) {
			double val = diffAverage + currentControl[i];
			if (i - halfSigma > 0) {
				val += delta * -1;
			} else {
				val += delta;
			}
			shift[i] = val;
		}
		return shift;
	}

	// private static double getHalfDelta(double f) {
	// return 1 * (0.5 - (1 / (2 - f)));
	// }

	private static double getHalfSigma(int halfCount, double shift) {
		return halfCount - shift;
	}

	private static int getHalfCount(int count) {
		return (int) (count + 0.5) / 2;
	}

	// private static int getHalfNull(double nullD) {
	// return (int) (nullD + 0.5) / 2;
	// }

	public static double getDisomyF(double delta) {
		if (Double.isNaN(delta)) {
			return 0;
		}
		if (delta == .5) {
			return 1;
		}
		double f = (4 * delta);
		double div = 2 * delta + 1;
		return f / div;
	}

	private static double getDelta(double f, MOSAIC_TYPE type) {
		double delta = Double.NaN;
		switch (type) {
			case MONOSOMY_DISOMY:
				delta = -1 * (0.5 - (1 / (2 - f)));
				break;
			case TRISOMY_DISOMY:
				delta = 0.5 - (1 / (2 + f));
				break;
			default:
				break;

		}
		return delta;
	}

	public SampleMosiac getSampleMosiac() {
		return sampleMosiac;
	}

	public static SampleMosiac prep(Project proj, String sampleName, String qcMetric,
																	int numControls) {
		if (numControls > proj.getSamples().length - 1) {
			proj.getLog()
					.reportTimeWarning("Only have "	+ proj.getSamples().length
															+ " total samples, changing number of controls to "
															+ (proj.getSamples().length - 1));
			numControls = proj.getSamples().length - 1;
		}
		SampleMosiacBase[] controls = new SampleMosiacBase[numControls];
		int sampleIndex = ext.indexOfStr(sampleName, proj.getSamples());
		SampleQC sampleQC = SampleQC.loadSampleQC(proj, LrrSd.SAMPLE_COLUMN, new String[] {qcMetric},
																							true, false, null);
		double[] data = sampleQC.getDataFor(qcMetric);
		int[] sorted = Sort.quicksort(data);
		int sampleSortIndex = ext.indexOfStr(sampleIndex + "", Array.toStringArray(sorted));
		int minIndex = Math.max(0, sampleSortIndex - numControls);
		if (sampleSortIndex == minIndex) {
			minIndex++;
		}
		int maxIndex = Math.min(proj.getSamples().length, sampleSortIndex + numControls);
		if (sampleSortIndex == maxIndex) {
			maxIndex--;
		}
		ArrayList<Double> controlData = new ArrayList<Double>();
		ArrayList<String> controlDataNames = new ArrayList<String>();

		double sampData = data[sorted[sampleSortIndex]];
		for (int i = minIndex; i < maxIndex; i++) {
			if (i != sampleSortIndex) {
				controlData.add(data[sorted[i]]);
				controlDataNames.add(proj.getSamples()[sorted[i]]);
			}
		}
		if (minIndex == maxIndex) {
			controlData.add(data[sorted[minIndex]]);
			controlDataNames.add(proj.getSamples()[sorted[minIndex]]);
		}

		double[] distanceToSample = Array.InplaceAbs(Array.InplaceSub(Doubles.toArray(controlData),
																																	sampData));
		int[] minDists = Sort.quicksort(distanceToSample);

		for (int i = 0; i < controls.length; i++) {
			int index = minDists[i];
			SampleMosiacBase control = new SampleMosiacBase(proj, controlDataNames.get(index), qcMetric,
																											controlData.get(index),
																											distanceToSample[index], true);
			controls[i] = control;
		}
		return new SampleMosiac(proj, sampleName, qcMetric, sampData, controls);
	}

	public static class SampleMosiac extends SampleMosiacBase {

		private SampleMosiacBase[] controls;
		private CDF smoothedControl;
		private boolean loaded;

		public SampleMosiac(SampleMosiac sampleMosiac) {
			super(sampleMosiac);
			controls = sampleMosiac.controls;
			smoothedControl = sampleMosiac.smoothedControl;
			loaded = sampleMosiac.loaded;
		}

		public CDF getSmoothedControl() {
			return smoothedControl;
		}

		public SampleMosiac(Project proj, String sampleName, String qcMetric, double qcValue,
												SampleMosiacBase[] controls) {
			super(proj, sampleName, qcMetric, qcValue, qcValue, false);
			this.controls = controls;
			loaded = false;
		}

		// public SampleMosiacBase[] getControls() {
		// return controls;
		// }

		public void load(int numThreads) {
			if (!loaded) {
				load();
				WorkerHive<SampleMosiacBase> hive =
																					new WorkerHive<MosaicismQuant.SampleMosiacBase>(numThreads,
																																													10,
																																													getProj().getLog());
				for (SampleMosiacBase control2 : controls) {
					hive.addCallable(control2);
				}
				hive.execute(true);
				ArrayList<SampleMosiacBase> results = hive.getResults();
				controls = results.toArray(new SampleMosiacBase[results.size()]);
				loaded = true;
			}

		}

		public void developCDFs(MarkerSet markerSet, int[][] indicesByChr, Segment seg, double minBaf,
														double maxBaf) {

			developCDF(markerSet, indicesByChr, seg, minBaf, maxBaf, -1);

			double[] smoothed = new double[getCdf().getVals().length];
			for (int i = 0; i < controls.length; i++) {
				controls[i].developCDF(	markerSet, indicesByChr, seg, minBaf, maxBaf,
																getCdf().getVals().length);
				if (i > 0) {// smoothed random
					CDF tmp = controls[i].getCdf();
					for (int j = 0; j < smoothed.length; j++) {
						smoothed[j] += tmp.getVals()[tmp.getOrder()[j]];
					}
				}
			}
			CDF tmpRandomControl = controls[0].getCdf();
			for (int i = 0; i < smoothed.length; i++) {
				double tmp = smoothed[i] + tmpRandomControl.getVals()[tmpRandomControl.getOrder()[i]];
				smoothed[i] = tmp / controls.length;
			}
			smoothed = CNVCaller.adjustBaf(smoothed, minBaf, maxBaf, false, proj.getLog());
			smoothedControl = new CDF(new BafSelection(smoothed, null));// many different indices for
																																	// controls
		}

		// private void plotCDFs(String output) {
		// try {
		// PrintWriter writer = new PrintWriter(new FileWriter(output));
		// writer.println("CDF\tSample\tSmoothedControl");
		// for (int i = 0; i < smoothedControl.getVals().length; i++) {
		// double cdf = (double) i / smoothedControl.getVals().length;
		// writer.println(cdf + "\t" + getCdf().getVals()[getCdf().getOrder()[i]] + "\t" +
		// smoothedControl.getVals()[smoothedControl.getOrder()[i]]);
		// }
		// writer.close();
		// RScatter rScatter = new RScatter(output, output + ".rscript",
		// ext.removeDirectoryInfo(output), output + ".jpeg", "CDF", new String[] { "Sample",
		// "SmoothedControl" }, SCATTER_TYPE.POINT, proj.getLog());
		// rScatter.setOverWriteExisting(true);
		// rScatter.execute();
		// } catch (Exception e) {
		// proj.getLog().reportError("Error writing to " + output);
		// proj.getLog().reportException(e);
		// }
		// }
	}

	private static class BafSelection {
		private final double[] bafs;
		private final int[] projectIndices;

		public BafSelection(double[] bafs, int[] projectIndices) {
			super();
			this.bafs = bafs;
			this.projectIndices = projectIndices;
		}

		public double[] getBafs() {
			return bafs;
		}

		public int[] getProjectIndices() {
			return projectIndices;
		}

	}

	private static class SampleMosiacBase implements Callable<SampleMosiacBase> {
		protected Project proj;
		protected Sample samp;
		protected CDF cdf;
		protected String sampleName;
		private final String qcMetric;
		private final double qcValue;
		private final double qcDistance;
		private final boolean control;

		public SampleMosiacBase(SampleMosiacBase s) {
			proj = s.proj;
			samp = s.samp;
			cdf = s.cdf;
			sampleName = s.sampleName;
			qcMetric = s.qcMetric;
			qcValue = s.qcValue;
			qcDistance = s.qcDistance;
			control = s.control;
		}

		protected Sample getSamp() {
			return samp;
		}

		public SampleMosiacBase(Project proj, String sampleName, String qcMetric, double qcValue,
														double qcDistance, boolean control) {
			super();
			this.proj = proj;
			this.sampleName = sampleName;
			this.qcMetric = qcMetric;
			this.qcValue = qcValue;
			this.qcDistance = qcDistance;
			this.control = control;
		}

		//
		// protected String getSampleName() {
		// return sampleName;
		// }

		protected Project getProj() {
			return proj;
		}

		protected CDF getCdf() {
			return cdf;
		}

		// protected String getQcMetric() {
		// return qcMetric;
		// }
		//
		// protected double getQcValue() {
		// return qcValue;
		// }
		//
		// protected double getQcDistance() {
		// return qcDistance;
		// }

		protected void load() {
			samp = proj.getFullSampleFromRandomAccessFile(sampleName);
		}

		protected void developCDF(MarkerSet markerSet, int[][] indicesByChr, Segment seg, double minBaf,
															double maxBaf, int numControlForce) {
			BafSelection bafSelection = selectBafs(	markerSet, indicesByChr, seg, minBaf, maxBaf,
																							numControlForce);
			double[] bafs = CNVCaller.adjustBaf(bafSelection.getBafs(), minBaf, maxBaf, false,
																					proj.getLog());
			cdf = new CDF(new BafSelection(bafs, bafSelection.getProjectIndices()));
		}

		private BafSelection selectBafs(MarkerSet markerSet, int[][] indicesByChr, Segment seg,
																		double minBaf, double maxBaf, int numControlForce) {
			String[] markersInSeg = markerSet.getMarkersIn(seg, indicesByChr);
			if (markersInSeg.length < 1) {
				return new BafSelection(new double[0], new int[0]);
			}
			int[] indicesInSeg = ext.indexLargeFactors(	markersInSeg, markerSet.getMarkerNames(), true,
																									proj.getLog(), true, false);
			double[] bafs = Array.toDoubleArray(samp.getBAFs());
			ArrayList<Integer> bafIndicesToUse = new ArrayList<Integer>();

			for (int index : indicesInSeg) {
				if (useBaf(minBaf, maxBaf, bafs[index])) {
					bafIndicesToUse.add(index);
				}
			}
			if (!control) {
				int[] projectIndices = Ints.toArray(bafIndicesToUse);
				double[] bafsSelection = Array.subArray(bafs, projectIndices);
				return new BafSelection(bafsSelection, projectIndices);
			} else {
				if (bafIndicesToUse.size() > numControlForce) {
					Collections.shuffle(bafIndicesToUse);
					ArrayList<Integer> bafIndicesToUseTmp = new ArrayList<Integer>();
					for (int i = 0; i < numControlForce; i++) {
						bafIndicesToUseTmp.add(bafIndicesToUse.get(i));
					}
					int[] projectIndices = Ints.toArray(bafIndicesToUseTmp);
					double[] bafsSelection = Array.subArray(bafs, projectIndices);
					return new BafSelection(bafsSelection, projectIndices);
				} else {
					boolean forward = true;
					boolean canGoForward = true;
					boolean canGoBackward = true;
					int nextIndexForward = indicesInSeg[indicesInSeg.length - 1];
					int nextIndexBackward = indicesInSeg[0];

					while (bafIndicesToUse.size() < numControlForce && (canGoBackward || canGoForward)) {
						if (forward) {
							if (nextIndexForward < markerSet.getMarkerNames().length) {
								if (useBaf(minBaf, maxBaf, bafs[nextIndexForward])) {
									bafIndicesToUse.add(nextIndexForward);
								}
							} else {
								canGoForward = false;
							}
							nextIndexForward++;
							forward = false;
						} else {
							if (nextIndexBackward >= 0) {
								if (useBaf(minBaf, maxBaf, bafs[nextIndexBackward])) {
									bafIndicesToUse.add(nextIndexBackward);
								}
							} else {
								canGoBackward = false;
							}
							nextIndexBackward--;
							forward = true;
						}
					}
					if (bafIndicesToUse.size() != numControlForce) {
						proj.getLog()
								.reportTimeWarning("Not enough baf values between "	+ minBaf + " and " + maxBaf
																		+ " for sample " + sampleName
																		+ ", imputting with random re-sampling");
						Collections.shuffle(bafIndicesToUse);
						ArrayList<Integer> bafIndicesToUseTmp = new ArrayList<Integer>();
						bafIndicesToUseTmp.addAll(bafIndicesToUse);
						int add = 0;
						while (bafIndicesToUseTmp.size() < numControlForce) {
							bafIndicesToUseTmp.add(bafIndicesToUse.get(add));
							add++;
							if (add >= bafIndicesToUse.size()) {
								add = 0;
								Collections.shuffle(bafIndicesToUse);
							}
						}
						int[] projectIndices = Ints.toArray(bafIndicesToUseTmp);
						double[] bafsSelection = Array.subArray(bafs, projectIndices);
						return new BafSelection(bafsSelection, projectIndices);
					} else {
						int[] projectIndices = Ints.toArray(bafIndicesToUse);
						double[] bafsSelection = Array.subArray(bafs, projectIndices);
						return new BafSelection(bafsSelection, projectIndices);
					}

				}
			}
		}

		private boolean useBaf(double minBaf, double maxBaf, double baf) {
			return Numbers.isFinite(baf) && baf > minBaf && baf < maxBaf;
		}

		@Override
		public SampleMosiacBase call() throws Exception {
			load();
			return this;
		}
	}

	private static class CDF {
		BafSelection bafSelection;
		private final int[] order;

		public CDF(BafSelection bafSelection) {
			super();
			this.bafSelection = bafSelection;
			order = Sort.quicksort(bafSelection.getBafs());
		}

		// public BafSelection getBafSelection() {
		// return bafSelection;
		// }

		public double[] getValsInOrder() {
			double[] orderedVals = new double[bafSelection.getBafs().length];
			for (int i = 0; i < orderedVals.length; i++) {
				orderedVals[i] = bafSelection.getBafs()[order[i]];
			}
			return orderedVals;
		}

		public double[] getVals() {
			return bafSelection.getBafs();
		}

		public int[] getOrder() {
			return order;
		}

	}

	public static class MosaicQuantResults implements Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 1L;
		private final String sample;
		private final double[] fs;
		private final double[] shifts;
		private final int[] numMarkers;
		private final int numErrors;

		// private MOSAIC_TYPE type;

		public MosaicQuantResults(String sample, double[] fs, double[] shifts, int[] numMarkers,
															int numErrors, MOSAIC_TYPE type) {
			super();
			this.sample = sample;
			this.fs = fs;
			this.shifts = shifts;
			this.numMarkers = numMarkers;
			this.numErrors = numErrors;
			// this.type = type;
		}

		public String getSample() {
			return sample;
		}

		public double[] getFs() {
			return fs;
		}

		public double[] getShifts() {
			return shifts;
		}

		public int[] getNumMarkers() {
			return numMarkers;
		}

		public int getNumErrors() {
			return numErrors;
		}

	}

	public static class MosaicQuantWorker implements Callable<MosaicQuantResults[]> {
		private final LocusSet<Segment> bins;
		private final Project proj;
		private final String sample;
		private final int numControls;
		private final MOSAIC_TYPE[] types;

		public MosaicQuantWorker(	Segment[] bins, Project proj, String sample, MOSAIC_TYPE[] types,
															int numControls) {
			super();
			this.bins = new LocusSet<Segment>(bins, true, proj.getLog()) {

				/**
				 *
				 */
				private static final long serialVersionUID = 1L;

			};
			this.proj = proj;
			this.sample = sample;
			this.numControls = numControls;
			this.types = types;
		}

		public MosaicQuantWorker(	LocusSet<Segment> bins, Project proj, String sample,
															MOSAIC_TYPE[] types, int numControls) {
			super();
			this.bins = bins;
			this.proj = proj;
			this.sample = sample;
			this.numControls = numControls;
			this.types = types;
		}

		@Override
		public MosaicQuantResults[] call() throws Exception {
			SampleMosiac sampleMosiac = prep(proj, sample, "BAF1585_SD", numControls);
			sampleMosiac.load();
			System.out.println(bins.getLoci()[0].getUCSClocation());
			double[][] fs = new double[types.length][bins.getLoci().length];
			double[][] shifts = new double[types.length][bins.getLoci().length];
			int[][] numMarkers = new int[types.length][bins.getLoci().length];

			int numErrors = 0;
			MosaicParamsBuilder builder = new MosaicParamsBuilder();
			ComputeParams params = builder.build();
			MarkerSet markerSet = proj.getMarkerSet();
			int[][] indices = markerSet.getIndicesByChr();
			MosaicQuantResults[] results = new MosaicQuantResults[types.length];
			for (int i = 0; i < bins.getLoci().length; i++) {
				if ((i + 1) % 10 == 0 || bins.getLoci().length < 100) {
					proj.getLog().reportTimeInfo("Currently on "	+ bins.getLoci()[i].getUCSClocation()
																				+ " for sample " + sample);
				}
				for (int j = 0; j < types.length; j++) {
					MosaicismQuant mosiacismQuant =
																				new MosaicismQuant(	proj, sampleMosiac, types[j], params,
																														bins.getLoci()[i], markerSet, indices);
					mosiacismQuant.prepareMosiacQuant(1, MIN_BAF, MAX_BAF);
					numMarkers[j][i] = mosiacismQuant.getSampleMosiac().getCdf().getVals().length;
					// System.out.println("Num Markers\t" + numMarkers[j][i]);
					if (numMarkers[j][i] > 0) {
						double[] x = params.getX();
						CobylaExitStatus cobylaExitStatus = Cobyla.FindMinimum(	mosiacismQuant,
																																		params.getX().length,
																																		params.getCon().length, x,
																																		params.getX_Bounds(), 100, .001,
																																		0, 50000);
						switch (cobylaExitStatus) {
							case DivergingRoundingErrors:
								proj.getLog()
										.reportTimeWarning("Could not Estimate "	+ bins.getLoci()[i].getUCSClocation()
																				+ " due to " + cobylaExitStatus);
								fs[j][i] = Double.NaN;
								shifts[j][i] = Double.NaN;
								numErrors++;
								break;
							case MaxIterationsReached:
								proj.getLog()
										.reportTimeWarning("Could not Estimate "	+ bins.getLoci()[i].getUCSClocation()
																				+ " due to " + cobylaExitStatus);
								fs[j][i] = Double.NaN;
								shifts[j][i] = Double.NaN;
								numErrors++;
								break;
							case Normal:
								fs[j][i] = x[0];
								shifts[j][i] = x[1];
								break;
							default:
								break;

						}
					} else {
						fs[j][i] = Double.NaN;
						shifts[j][i] = Double.NaN;
					}
				}

			}
			for (int i = 0; i < results.length; i++) {
				results[i] = new MosaicQuantResults(sample, fs[i], shifts[i], numMarkers[i], numErrors,
																						types[i]);
			}
			return results;
		}
	}

	public static class MosaicQuantProducer extends AbstractProducer<MosaicQuantResults[]> {

		private final Project proj;
		private final String[] samples;
		private final LocusSet<Segment> bins;
		private int index;
		private final int numControls;
		private final MOSAIC_TYPE[] types;

		public MosaicQuantProducer(	Project proj, String[] samples, LocusSet<Segment> bins,
																MOSAIC_TYPE[] types, int numControls) {
			super();
			this.proj = proj;
			this.samples = samples;
			this.bins = bins;
			this.numControls = numControls;
			this.types = types;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<MosaicQuantResults[]> next() {
			MosaicQuantWorker worker = new MosaicQuantWorker(	bins, proj, samples[index], types,
																												numControls);
			index++;
			return worker;
		}
	}

	public static class FullMosiacResults implements Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 1L;
		private final MosaicQuantResults[] sampleMosaicQuantResults;
		private final LocusSet<Segment> set;

		public FullMosiacResults(MosaicQuantResults[] sampleMosaicQuantResults, LocusSet<Segment> set) {
			super();
			this.sampleMosaicQuantResults = sampleMosaicQuantResults;
			this.set = set;
		}

		public MosaicQuantResults[] getSampleMosaicQuantResults() {
			return sampleMosaicQuantResults;
		}

		public LocusSet<Segment> getSet() {
			return set;
		}

		private void writeSerial(String fileName) {
			SerializedFiles.writeSerial(this, fileName, true);
		}

		private static FullMosiacResults readSerial(String fileName, Logger log) {
			return (FullMosiacResults) SerializedFiles.readSerial(fileName, false, log, false, true);
		}

	}

	public static FullMosiacResults quantMosaic(Project proj, int bpWindow, int numControls,
																							int numThreads) {
		System.out.println("Sorry, reformat for multitype");
		System.exit(1);
		String outDir = proj.PROJECT_DIRECTORY.getValue() + "MosaicResults/";
		FullMosiacResults fullMosiacResults = null;
		new File(outDir).mkdirs();
		String rootOut = outDir + "MosaicResults_window_" + bpWindow + "_controls_" + numControls;
		String out = rootOut + ".txt";
		String ser = rootOut + ".ser";
		if (!Files.exists(ser)) {
			ReferenceGenome referenceGenome =
																			new ReferenceGenome(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(),
																													proj.getLog());

			LocusSet<Segment> set = referenceGenome.getBins(bpWindow);
			MosaicQuantProducer mProducer = new MosaicQuantProducer(proj, proj.getSamples(), set,
																															MOSAIC_TYPE.values(), numControls);
			WorkerTrain<MosaicQuantResults[]> train =
																							new WorkerTrain<MosaicismQuant.MosaicQuantResults[]>(	mProducer,
																																																		numThreads,
																																																		1,
																																																		proj.getLog());

			try {
				PrintWriter writer = new PrintWriter(new FileWriter(out));
				writer.print("Sample");
				for (int i = 0; i < set.getLoci().length; i++) {
					writer.print("\t" + set.getLoci()[i].getUCSClocation());
				}
				writer.println();
				MosaicQuantResults[] sampleMosaicQuantResults =
																											new MosaicQuantResults[proj.getSamples().length];
				int index = 0;
				while (train.hasNext()) {

					MosaicQuantResults[] results = train.next();
					for (MosaicQuantResults result : results) {
						sampleMosaicQuantResults[index] = result;

						// writer.println(results.getSample() + "\t" + Array.toStr(results.getFs()));

					}
					index++;
				}

				writer.close();
				fullMosiacResults = new FullMosiacResults(sampleMosaicQuantResults, set);
				fullMosiacResults.writeSerial(ser);
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + out);
				proj.getLog().reportException(e);
			}
		} else {
			proj.getLog().reportTimeInfo("Loading pre-computed results from " + ser);
			fullMosiacResults = FullMosiacResults.readSerial(ser, proj.getLog());
		}
		return fullMosiacResults;

	}

	// select case marker/indices check
	// select controls check
	// select control markers check
	// adjust case bafs check
	// adjust control bafs check
	// create case cdf check
	// create control cdf check
	// Non-linear fitting

	public static void test(Project proj) {

		SampleMosiac sampleMosiac = prep(proj, "7330686030_R02C01", "BAF1585_SD", 5);
		// SampleMosiac sampleMosiac = prep(proj, "7355066051_R03C01", "BAF1585_SD", 5);

		MarkerSet markerSet = proj.getMarkerSet();
		int[][] indices = markerSet.getIndicesByChr();
		// ReferenceGenome referenceGenome = new
		// ReferenceGenome(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(), proj.getLog());
		// LocusSet<Segment> set = referenceGenome.getBins(500000);
		// Segment segTest = new Segment("chr17:42,963,198-78,940,173");
		// Segment segTest = new Segment("chr17:42,963,198-78,940,173");
		Segment segTest = new Segment((byte) 11, 0, Integer.MAX_VALUE);
		int[] currentIndices = ext.indexLargeFactors(	markerSet.getMarkersIn(segTest, indices),
																									markerSet.getMarkerNames(), true, proj.getLog(),
																									true, false);
		MosaicParamsBuilder builder = new MosaicParamsBuilder();
		ComputeParams params = builder.build();
		MosaicismQuant mosiacismQuant =
																	new MosaicismQuant(	proj, sampleMosiac, MOSAIC_TYPE.TRISOMY_DISOMY,
																											params, segTest, markerSet, indices);
		mosiacismQuant.prepareMosiacQuant(5, MIN_BAF, MAX_BAF);

		String testDir = proj.PROJECT_DIRECTORY.getValue() + "TestMosaic/";
		new File(testDir).mkdirs();
		String out = testDir	+ ext.replaceWithLinuxSafeCharacters(segTest.getUCSClocation(), true)
									+ ".txt";
		int[] autosomal = proj.getAutosomalMarkerIndices();
		double[] val0_33 = Array.getValuesBetween(Array.toDoubleArray(Array.subArray(	mosiacismQuant.getSampleMosiac()
																																																.getSamp()
																																																.getBAFs(),
																																									autosomal)),
																							0, .33);
		double[] val33_66 = Array.getValuesBetween(
																								Array.toDoubleArray(Array.subArray(	mosiacismQuant.getSampleMosiac()
																																																	.getSamp()
																																																	.getBAFs(),
																																										autosomal)),
																								.33, .66);
		double[] val66_33 = Array.getValuesBetween(
																								Array.toDoubleArray(Array.subArray(	mosiacismQuant.getSampleMosiac()
																																																	.getSamp()
																																																	.getBAFs(),
																																										autosomal)),
																								.66, 1);
		double[] vars = new double[3];
		double[] means = new double[3];
		vars[0] = Math.pow(Array.stdev(val0_33), 2);
		vars[1] = Math.pow(Array.stdev(val33_66), 2);
		vars[2] = Math.pow(Array.stdev(val66_33), 2);
		means[0] = Array.mean(val0_33);
		means[1] = Array.mean(val33_66);
		means[2] = Array.mean(val66_33);
		// int[] totals = new int[] { val0_33.length, val33_66.length, val66_33.length };
		double[] blanks = new double[3];
		Arrays.fill(blanks, 1);
		double[] props = new double[3];
		for (int i = 0; i < props.length; i++) {
			props[i] = (double) 1 / 3;
		}
		GaussianMixtureDistribution gd = new GaussianMixtureDistribution(means, vars, props);

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(out));

			int[] mas = new int[] {10, 20, 50, 100, 250, 500, 1000};
			String[] MAtitles = Array.tagOn(Array.toStringArray(mas), "MA", null);
			writer.println("Position\tBaf\tRandBaf\tPval\tBaseLinePval\tDistMember\t"
											+ Array.toStr(MAtitles));
			double[] pvals = new double[currentIndices.length];
			double[] rand = new double[currentIndices.length];
			double baseLine = 0;
			int[] distMember = new int[currentIndices.length];
			for (int j = 0; j < gd.distributions().length; j++) {
				baseLine += gd.distributions()[j].probability(means[j] + 2 * Math.sqrt(vars[j]))
										* Math.sqrt(vars[j]);
			}

			for (int i = 0; i < currentIndices.length; i++) {
				double baf = mosiacismQuant.getSampleMosiac().getSamp().getBAFs()[currentIndices[i]];
				pvals[i] = 0;
				for (int j = 0; j < gd.distributions().length; j++) {
					if (j == 0 || j == 2) {
						if (j == 0 && Math.abs(baf - means[j]) < 2 * Math.sqrt(vars[j])) {
							pvals[i] = Double.NaN;
							System.out.println("0"	+ "\t" + baf + "\t" + (baf - means[j]) + "\t"
																	+ (means[j] + Math.sqrt(vars[j])));

						} else if (Math.abs(baf - means[j]) < 2 * Math.sqrt(vars[j])) {
							pvals[i] = Double.NaN;

							System.out.println("2"	+ "\t" + baf + "\t" + (baf - means[j]) + "\t"
																	+ (means[j] + Math.sqrt(vars[j])));

						}

					}
					double test = !Numbers.isFinite(baf) ? 0 : (double) (baf);
					System.out.println(test);
					double tmp = gd.distributions()[j].probability(test) * Math.sqrt(vars[j]);

					if (tmp > pvals[i] && !Double.isNaN(pvals[i])) {
						pvals[i] = tmp;
						distMember[i] = j;
					}
				}

				System.out.println(pvals[i]);
				rand[i] = gd.generate();
				if (rand[i] < 0) {
					rand[i] = 0;
				}
				if (rand[i] > 1) {
					rand[i] = 1;
				}
			}

			double[][] madatas = new double[mas.length][];
			for (int i = 0; i < madatas.length; i++) {
				madatas[i] = Array.movingAverageForward(mas[i], pvals, true);
			}

			for (int i = 0; i < currentIndices.length; i++) {
				double baf = mosiacismQuant.getSampleMosiac().getSamp().getBAFs()[currentIndices[i]];

				writer.print(markerSet.getPositions()[currentIndices[i]]	+ "\t"
											+ (!Numbers.isFinite(baf) ? 0 : baf) + "\t" + rand[i] + "\t" + pvals[i] + "\t"
											+ baseLine + "\t" + distMember[i]);
				for (double[] madata : madatas) {
					writer.print("\t" + madata[i]);
				}
				writer.println();
			}
			writer.close();
			String outBase = out + ".base";
			ArrayList<RScatter> rd = new ArrayList<RScatter>();
			RScatter rsScatter = new RScatter(out, outBase + ".rscript", ext.removeDirectoryInfo(outBase),
																				outBase + ".jpeg", "Position",
																				new String[] {"Baf", "Pval", "BaseLinePval"},
																				SCATTER_TYPE.POINT, proj.getLog());
			rsScatter.setOverWriteExisting(true);
			rsScatter.execute();
			rd.add(rsScatter);
			String outMA = out + ".ma";
			RScatter rsScatterMa = new RScatter(out, outMA + ".rscript", ext.removeDirectoryInfo(outMA),
																					outMA + ".jpeg", "Position",
																					Array.concatAll(new String[] {"Baf", "BaseLinePval"},
																													MAtitles),
																					SCATTER_TYPE.POINT, proj.getLog());
			rsScatterMa.setOverWriteExisting(true);
			rsScatterMa.execute();
			rd.add(rsScatterMa);
			String outMA50 = out + ".ma50";

			RScatter rsScatter50 = new RScatter(out, outMA50 + ".rscript",
																					ext.removeDirectoryInfo(outMA50), outMA50 + ".jpeg",
																					"Position",
																					new String[] {"Baf", "Pval", "MA50", "BaseLinePval"},
																					SCATTER_TYPE.POINT, proj.getLog());
			rsScatter50.setOverWriteExisting(true);
			rsScatter50.execute();
			rd.add(rsScatter50);
			String outRand = out + ".rand";

			RScatter rsScatterRand = new RScatter(out, outRand + ".rscript",
																						ext.removeDirectoryInfo(outRand), outRand + ".jpeg",
																						"Position", new String[] {"Baf", "RandBaf"},
																						SCATTER_TYPE.POINT, proj.getLog());
			rsScatterRand.setOverWriteExisting(true);
			// rsScatterRand.setTitle("Sample SD=" + sd + ", mean=" + mean);
			rsScatterRand.execute();
			rd.add(rsScatterRand);

			RScatters rScatters = new RScatters(rd.toArray(new RScatter[rd.size()]), out + ".rscript",
																					out + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1,
																					PLOT_DEVICE.PDF, proj.getLog());
			rScatters.execute();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + out);
			proj.getLog().reportException(e);
		}

		System.exit(1);

		// double[][] dubs = new double[vals.length][1];
		// Point[] pointsSample = new Point[vals.length];
		// for (int i = 0; i < dubs.length; i++) {
		// pointsSample[i] = new Point(vals[i]);
		// dubs[i] = new double[] { vals[i], 1 };
		// }
		//
		// GaussianDistribution gDistribution = new GaussianDistribution(.5, Math.pow(.04, 2));
		// int index = 0;
		// Point[] points = new Point[10000];
		// while (index < points.length) {
		// double d = gDistribution.generate();
		// if (d > .15 && d < .85) {
		// points[index] = new Point(d);
		// index++;
		// }
		// }
		// double max = -1 * Double.MAX_VALUE;
		// for (int i = 30; i < 70; i++) {
		// for (int j = 30; j < 70; j++) {
		// for (int j2 = 0; j2 < 10; j2++) {
		// UnivariateGaussianMixtureModel ugmmInit = ExpectationMaximization.initialize(pointsSample,
		// new double[] { (double) i / 100, (double) j / 100 }, new double[] { (double) j2 / 1000,
		// (double) j2 / 1000 });
		//
		// UnivariateGaussianMixtureModel ugmmFinal = ExpectationMaximization.run(pointsSample,
		// ugmmInit);
		// double lods = ExpectationMaximization.logLikelihood(pointsSample, ugmmFinal);
		// if (lods > max) {
		// max = lods;
		// System.out.println(ugmmFinal.toString());
		//
		// System.out.println(lods + "\t" + i + "\t" + j + "\t" + j2);
		//
		// }
		// }
		//
		// }
		// }

		System.exit(1);

		// MultivariateNormalMixtureExpectationMaximization mle = new
		// MultivariateNormalMixtureExpectationMaximization(dubs);
		// MixtureMultivariateNormalDistribution initialMix =
		// MultivariateNormalMixtureExpectationMaximization.estimate(dubs, 2);
		// mle.fit(initialMix, 500, 1E-5);
		// MixtureMultivariateNormalDistribution finalMix = mle.getFittedModel();
		// List<Pair<Double, MultivariateNormalDistribution>> dists = finalMix.getComponents();
		// for (Pair<Double, MultivariateNormalDistribution> pair : dists) {
		// System.out.println(Array.toStr(pair.getValue().getMeans()));
		// System.out.println(Array.toStr(pair.getValue().getStandardDeviations()));
		//
		// }
		// mle.getFittedModel().getComponents();
		// System.out.println(mle.getLogLikelihood());
		// // mle.fit(new Mi);
		System.exit(1);
		// String testDir = proj.PROJECT_DIRECTORY.getValue() + "TestMosaic/";
		// new File(testDir).mkdirs();
		// String out = testDir + ext.replaceWithLinuxSafeCharacters(segTest.getUCSClocation(), true);
		// mosiacismQuant.getSampleMosiac().plotCDFs(out);
		// double[] x = params.getX();
		// CobylaExitStatus cobylaExitStatus = Cobyla.FindMinimum(mosiacismQuant, params.getX().length,
		// params.getCon().length, x, params.getX_Bounds(), 1, .001, 0, 50000);

		// for (int i = 0; i < set.getLoci().length; i++) {
		//
		// proj.getLog().reportTimeInfo(set.getLoci()[i].getUCSClocation());
		// try {
		//
		//
		// String testDir = proj.PROJECT_DIRECTORY.getValue() + "TestMosaic/";
		// new File(testDir).mkdirs();
		// String out = testDir + ext.replaceWithLinuxSafeCharacters(segTest.getUCSClocation(), true);
		// mosiacismQuant.getSampleMosiac().plotCDFs(out);
		//
		// System.out.println(mosiacismQuant.Compute(1, 1, params.getX(), params.getCon()));
		// double[] x = params.getX();
		// CobylaExitStatus cobylaExitStatus = Cobyla.FindMinimum(mosiacismQuant, params.getX().length,
		// params.getCon().length, x, params.getX_Bounds(), 1, .001, 0, 50000);
		// System.out.println(Array.toStr(x));
		// System.out.println(cobylaExitStatus);
		// BruteForce bruteForce = new BruteForce(mosiacismQuant, params);
		// } catch (Exception e) {
		// proj.getLog().reportException(e);
		// }
		// double[][] results = bruteForce.force();
		// }

		//
		//
		//
		// try {
		// PrintWriter writer = new PrintWriter(new FileWriter(testDir +
		// ext.replaceWithLinuxSafeCharacters(segTest.getUCSClocation(), true)+"force.txt"));
		// writer.println("F\tShift\tVAL");
		// for (int i = 0; i < results.length; i++) {
		// writer.println(Array.toStr(results[i]));
		// }
		// writer.close();
		// } catch (Exception e) {
		// proj.getLog().reportError("Error writing to " + testDir + "force.txt");
		// proj.getLog().reportException(e);
		// }
		System.exit(1);
	}

	// private static class BruteForce {
	// private MosaicismQuant mosiacismQuant;
	// private ComputeParams computeParams;
	//
	// public BruteForce(MosaicismQuant mosiacismQuant, ComputeParams computeParams) {
	// super();
	// this.mosiacismQuant = mosiacismQuant;
	// this.computeParams = computeParams;
	// }
	//
	// public double[][] force() {
	// double f = computeParams.getMinF();
	// double shift = computeParams.getMinShift();
	// ArrayList<Double> fs = new ArrayList<Double>();
	// ArrayList<Double> shifts = new ArrayList<Double>();
	//
	// while (f < computeParams.getMaxF()) {
	// f = f + .1;
	// fs.add(f);
	// }
	// while (shift < computeParams.getMaxShift()) {
	// shift = shift + 20;
	// shifts.add(shift);
	// }
	// double[][] results = new double[shifts.size() * fs.size()][3];
	// int index = 0;
	// for (int i = 0; i < fs.size(); i++) {
	// for (int j = 0; j < shifts.size(); j++) {
	// double[] params = new double[] { fs.get(i), shifts.get(j), 0 };
	// double result = mosiacismQuant.Compute(1, 1, params, null);
	// double[] tmp = new double[] { fs.get(i), shifts.get(j), result };
	// results[index] = tmp;
	// index++;
	// }
	// }
	// return results;
	// }
	//
	// }

	// private static BruteForceResult{
	//
	// }

	public static class ComputeParams {

		private final double f;
		private final double minF;
		private final double maxF;
		private final double shift;
		private final double minShift;
		private final double maxShift;
		private final double nullD;
		private final double minNullD;
		private final double maxNullD;

		private double[] getX() {
			return new double[] {f == 0 ? f + .1 : f, shift, nullD};
		}

		private double[][] getX_Bounds() {
			double[][] xbounds = new double[3][];
			xbounds[0] = new double[] {minF, maxF};
			xbounds[1] = new double[] {minShift, maxShift};
			xbounds[2] = new double[] {minNullD, maxNullD};
			return xbounds;
		}

		private double[] getCon() {
			return new double[] {minF, maxF, minShift, maxShift, minNullD, maxNullD};
		}

		public double getMinF() {
			return minF;
		}

		public double getMaxF() {
			return maxF;
		}

		public double getMinShift() {
			return minShift;
		}

		public double getMaxShift() {
			return maxShift;
		}

		public double getMinNullD() {
			return minNullD;
		}

		public double getMaxNullD() {
			return maxNullD;
		}

		public static class MosaicParamsBuilder {
			private double f = .5;
			private double minF = 0;
			private double maxF = 100;
			private double shift = 0;
			private double minShift = -100;
			private double maxShift = 300;
			private double nullD = 0;
			private double minNullD = 0;
			private double maxNullD = 100;

			public MosaicParamsBuilder f(double f) {
				this.f = f;
				return this;
			}

			public MosaicParamsBuilder minF(double minF) {
				this.minF = minF;
				return this;
			}

			public MosaicParamsBuilder maxF(double maxF) {
				this.maxF = maxF;
				return this;
			}

			public MosaicParamsBuilder shift(double shift) {
				this.shift = shift;
				return this;
			}

			public MosaicParamsBuilder minShift(double minShift) {
				this.minShift = minShift;
				return this;
			}

			public MosaicParamsBuilder maxShift(double maxShift) {
				this.maxShift = maxShift;
				return this;
			}

			public MosaicParamsBuilder nullD(double nullD) {
				this.nullD = nullD;
				return this;
			}

			public MosaicParamsBuilder minNullD(double minNullD) {
				this.minNullD = minNullD;
				return this;
			}

			public MosaicParamsBuilder maxNullD(double maxNullD) {
				this.maxNullD = maxNullD;
				return this;
			}

			public ComputeParams build() {
				return new ComputeParams(this);
			}
		}

		private ComputeParams(MosaicParamsBuilder builder) {
			f = builder.f;
			minF = builder.minF;
			maxF = builder.maxF;
			shift = builder.shift;
			minShift = builder.minShift;
			maxShift = builder.maxShift;
			nullD = builder.nullD;
			minNullD = builder.minNullD;
			maxNullD = builder.maxNullD;
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		int bpWindow = 1000000;
		int numControls = 5;
		int numThreads = 24;

		String usage = "\n" + "cnv.hmm.CircDuCnv requires 0-1 arguments\n";
		usage += "   (1) filename (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) window (i.e. window=" + bpWindow + " (default))\n" + "";
		usage += "   (3) numControls (i.e. numControls=" + numControls + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(4, numThreads);

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("window=")) {
				bpWindow = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("numControls=")) {
				numControls = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(arg);
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
			// test(proj);
			quantMosaic(proj, bpWindow, numControls, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
