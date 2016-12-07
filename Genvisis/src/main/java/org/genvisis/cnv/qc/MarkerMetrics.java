package org.genvisis.cnv.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.AnnotationCollection;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.qc.MendelErrors.MendelErrorCheck;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.AlleleFreq;
import org.genvisis.common.Array;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.db.FilterDB;
import org.genvisis.stats.LeastSquares;
import org.genvisis.stats.LogisticRegression;
import org.genvisis.stats.RegressionModel;
import org.genvisis.stats.Ttest;

import com.google.common.primitives.Booleans;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Floats;

public class MarkerMetrics {
	public static final String[] FULL_QC_HEADER = {	"MarkerName", "Chr", "CallRate", "meanTheta_AA",
																									"meanTheta_AB", "meanTheta_BB", "diffTheta_AB-AA",
																									"diffTheta_BB-AB", "sdTheta_AA", "sdTheta_AB",
																									"sdTheta_BB", "meanR_AA", "meanR_AB", "meanR_BB",
																									"num_AA", "num_AB", "num_BB", "pct_AA", "pct_AB",
																									"pct_BB", "MAF", "HetEx", "num_NaNs", "LRR_SEX_z",
																									"LRR_SD", "LRR_num_NaNs", "MendelianErrors"};
	public static final String[] LRR_VARIANCE_HEADER = {"MarkerName", "Chr", "Position", "SD_LRR",
																											"MeanAbsLRR", "SD_BAF1585", "MeanAbsBAF1585"};
	private static final String[] MENDEL_HEADER = {"FID", "KID", "CHR", "SNP", "CODE", "ERROR"};


	public static final String DEFAULT_REVIEW_CRITERIA = "cnv/qc/default_review.criteria";
	public static final String DEFAULT_EXCLUSION_CRITERIA = "cnv/qc/default_exclusion.criteria";
	public static final String DEFAULT_COMBINED_CRITERIA = "cnv/qc/default_combined.criteria";
	public static final String DEFAULT_MENDEL_FILE_SUFFIX = ".mendel";

	/**
	 * @param proj
	 * @param samplesToExclude these samples will not be included in the QC computation
	 * @param markersToInclude compute qc over the markers in this file only
	 * @param numThreads
	 */
	public static void fullQC(Project proj, boolean[] samplesToExclude, String markersToInclude,
														boolean checkMendel, int numThreads) {
		String[] markerNames;
		String finalQcFile = proj.MARKER_METRICS_FILENAME.getValue(true, false);

		if (markersToInclude != null) {
			if (Files.exists(markersToInclude)) {
				markerNames = HashVec.loadFileToStringArray(markersToInclude, false, new int[] {0}, false);
			} else if (Files.exists(proj.PROJECT_DIRECTORY.getValue(false, true) + markersToInclude)) {
				markerNames = HashVec.loadFileToStringArray(proj.PROJECT_DIRECTORY.getValue(false, true) + markersToInclude, false, new int[] {0}, false);
			} else {
				proj.getLog().reportTimeWarning("Markers to QC file not found, using all markers in project: " + markersToInclude);
				markerNames = proj.getMarkerNames();
			}
		} else {
			markerNames = proj.getMarkerNames();
		}
		proj.verifyAndGenerateOutliers(false);

		if (numThreads <= 1) {
			fullQC(proj, samplesToExclude, markerNames, finalQcFile, checkMendel);
		} else {
			WorkerHive<Boolean> hive = new WorkerHive<Boolean>(numThreads, 10, proj.getLog());
			ArrayList<String[]> batches = Array.splitUpArray(markerNames, numThreads, proj.getLog());
			String[] tmpQc = new String[batches.size()];
			String[] tmpMendel = new String[batches.size()];
			for (int i = 0; i < batches.size(); i++) {
				String tmp = ext.addToRoot(finalQcFile, "tmp" + i);
				hive.addCallable(new MarkerMetricsWorker(	proj, samplesToExclude, batches.get(i), tmp,
																									checkMendel));
				tmpQc[i] = tmp;
				if (checkMendel) {
					tmpMendel[i] = ext.rootOf(tmp, false) + DEFAULT_MENDEL_FILE_SUFFIX;
				}
			}

			hive.execute(true);
			ArrayList<Boolean> complete = hive.getResults();
			if (Array.booleanArraySum(Booleans.toArray(complete)) == complete.size()) {
				Files.cat(tmpQc, finalQcFile, new int[0], proj.getLog());
				if (Files.exists(finalQcFile) && Files.countLines(finalQcFile, 1) == markerNames.length) {
					for (String element : tmpQc) {
						new File(element).delete();
					}
				} else {
					proj.getLog()
							.reportError("Could not collapse temporary marker files to " + finalQcFile);
				}

				if (checkMendel) {
					Files.cat(tmpMendel, ext.rootOf(finalQcFile, false) + DEFAULT_MENDEL_FILE_SUFFIX,
										new int[0], proj.getLog());
					if (Files.exists(ext.rootOf(finalQcFile, false) + DEFAULT_MENDEL_FILE_SUFFIX)) {
						for (String element : tmpMendel) {
							new File(element).delete();
						}
					} else {
						proj.getLog()
								.reportError("Could not collapse temporary mendel files to "
																	+ ext.rootOf(finalQcFile, false) + DEFAULT_MENDEL_FILE_SUFFIX);
					}
				}
			} else {
				proj.getLog().reportError("Could not complete marker QC");
			}
		}
	}


	private static void fullQC(	Project proj, boolean[] samplesToExclude, String[] markerNames,
															String fullPathToOutput, boolean checkMendel) {
		PrintWriter writer, mendelWriter = null;
		String[] samples;
		float[] thetas, rs, lrrs;
		MarkerData markerData;
		byte[] abGenotypes;
		String markerName;
		ClusterFilterCollection clusterFilterCollection;
		float gcThreshold, lrrsd;
		long time;
		// MarkerDataLoader markerDataLoader;
		String eol, mendelLine;
		int[] counts, sexes;
		double[] sumTheta, sumR, meanTheta, sdTheta;
		double temp, lrrSexZ;
		int numNaNs, numLRRNaNs, count;
		ArrayList<Float> aLRR;
		Logger log;
		boolean[] toInclude;
		MendelErrorCheck[] mecArr;

		log = proj.getLog();

		if (Files.isWindows()) {
			eol = "\r\n";
		} else {
			eol = "\n";
		}

		samples = proj.getSamples();
		clusterFilterCollection = proj.getClusterFilterCollection();
		gcThreshold = proj.getProperty(proj.GC_THRESHOLD).floatValue();
		sexes = getSexes(proj, samples);
		Pedigree pedigree = proj.loadPedigree();

		try {
			writer = new PrintWriter(new FileWriter(fullPathToOutput));
			writer.println(Array.toStr(FULL_QC_HEADER));

			if (pedigree != null && checkMendel) {
				mendelWriter = new PrintWriter(new FileWriter(ext.rootOf(fullPathToOutput, false)
																											+ DEFAULT_MENDEL_FILE_SUFFIX));
				mendelWriter.println(Array.toStr(MENDEL_HEADER));
			}

			// markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj,
			// markerNames);
			MDL mdl = new MDL(proj, proj.getMarkerSet(), markerNames, 2, 100);
			time = new Date().getTime();
			// for (int i = 0; i < markerNames.length; i++) {
			int index = 0;
			while (mdl.hasNext()) {
				markerData = mdl.next();
				index++;
				if (index % 1000 == 0) {
					log.report(ext.getTime() + "\tMarker " + index + " of " + markerNames.length);
				}

				markerName = markerData.getMarkerName();
				thetas = markerData.getThetas();
				rs = markerData.getRs();
				abGenotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerName,
																														gcThreshold, log);
				lrrs = markerData.getLRRs();
				aLRR = new ArrayList<Float>(samples.length);

				numLRRNaNs = 0;
				numNaNs = 0;
				counts = new int[4];
				sumTheta = new double[counts.length];
				sumR = new double[counts.length];
				for (int j = 0; j < samples.length; j++) {
					if (samplesToExclude == null || !samplesToExclude[j]) {
						counts[abGenotypes[j] + 1]++;
						sumTheta[abGenotypes[j] + 1] += thetas == null ? 0 : thetas[j];
						sumR[abGenotypes[j] + 1] += rs == null ? 0 : rs[j];
						if (thetas != null && Float.isNaN(thetas[j])) {
							numNaNs++;
						}
						if (lrrs != null) {
							if (Float.isNaN(lrrs[j])) {
								numLRRNaNs++;
							} else {
								aLRR.add(lrrs[j]);
							}
						}
					}
				}

				meanTheta = new double[counts.length];
				for (int j = 1; j < meanTheta.length; j++) {
					meanTheta[j] = sumTheta[j] / counts[j];
				}
				sdTheta = new double[counts.length];
				for (int j = 0; j < samples.length; j++) {
					if (samplesToExclude == null || !samplesToExclude[j]) {
						temp = thetas == null ? 0 : (thetas[j] - meanTheta[abGenotypes[j] + 1]);
						sdTheta[abGenotypes[j] + 1] += temp * temp;
					}
				}
				for (int j = 1; j < sdTheta.length; j++) {
					if (counts[j] == 0) {
						sdTheta[j] = Double.NaN;
					} else {
						sdTheta[j] = Math.sqrt(sdTheta[j] / (counts[j] - 1));
					}
				}

				if (lrrs != null && aLRR.size() > 0) {
					lrrsd = Array.stdev(Floats.toArray(aLRR), true);
					lrrSexZ = getSexZscore(sexes, lrrs, samplesToExclude, log);
				} else {
					lrrsd = Float.NaN;
					lrrSexZ = Float.NaN;
				}

				String mecCnt = ".";
				if (pedigree != null && checkMendel) {
					toInclude = samplesToExclude == null	? Array.booleanArray(samples.length, true)
																								: Array.booleanNegative(samplesToExclude);
					mecArr = Pedigree.PedigreeUtils.checkMendelErrors(pedigree, markerData, toInclude, null,
																														clusterFilterCollection, gcThreshold,
																														log);
					count = 0;
					for (int i = 0; i < mecArr.length; i++) {
						if (!toInclude[i]) {
							continue;
						}
						if (mecArr[i].getErrorCode() != -1) {
							count++;

							mendelLine = pedigree.getFID(i)	+ "\t" + pedigree.getIID(i) + "\t"
														+ markerData.getChr() + "\t" + markerName + "\t"
														+ mecArr[i].getErrorCode() + "\t" + mecArr[i].getError() + eol;

							mendelWriter.print(mendelLine);
						}
					}
					mecCnt = "" + count;
				}


				String line = markerName	+ "\t" + markerData.getChr() + "\t"
											+ (1 - ((float) counts[0] / (counts[0] + counts[1] + counts[2] + counts[3])))
											+ "\t" + meanTheta[1] + "\t" + meanTheta[2] + "\t" + meanTheta[3] + "\t"
											+ (meanTheta[2] - meanTheta[1]) + "\t" + (meanTheta[3] - meanTheta[2]) + "\t"
											+ sdTheta[1] + "\t" + sdTheta[2] + "\t" + sdTheta[3] + "\t"
											+ (sumR[1] / counts[1]) + "\t" + (sumR[2] / counts[2]) + "\t"
											+ (sumR[3] / counts[3]) + "\t" + counts[1] + "\t" + counts[2] + "\t"
											+ counts[3] + "\t"
											+ ((float) counts[1] / (counts[0] + counts[1] + counts[2] + counts[3])) + "\t"
											+ ((float) counts[2] / (counts[0] + counts[1] + counts[2] + counts[3])) + "\t"
											+ ((float) counts[3] / (counts[0] + counts[1] + counts[2] + counts[3])) + "\t"
											+ (float) (counts[1] < counts[3]	? (counts[1] + counts[2])
																												: (counts[2] + counts[3]))
												/ (counts[0] + counts[1] + 2 * counts[2] + counts[3])
											+ "\t" + AlleleFreq.HetExcess(counts[1], counts[2], counts[3])[0] + "\t"
											+ numNaNs + "\t" + lrrSexZ + "\t" + lrrsd + "\t" + numLRRNaNs + "\t" + mecCnt;
				writer.println(line);

				// if (line.length() > 25000) {
				// writer.print(line);
				// writer.flush();
				// line = "";
				// }
				// markerDataLoader.releaseIndex(i);
			}
			mdl.shutdown();
			if (pedigree != null && checkMendel) {
				mendelWriter.flush();
				mendelWriter.close();
			}
			// writer.print(line);
			writer.close();
			log.report("Finished analyzing " + markerNames.length + " in " + ext.getTimeElapsed(time));
		} catch (Exception e) {
			log.reportError("Error writing marker metrics to " + fullPathToOutput);
			log.reportException(e);
		}
	}

	/**
	 * Retrieves Sex coding (1=male, 2=female) for all samples
	 * <p>
	 *
	 * @author John Lane
	 */
	private static int[] getSexes(Project proj, String[] samples) {
		int[] sexes = new int[samples.length];
		SampleData sampleData = proj.getSampleData(2, false);
		for (int i = 0; i < samples.length; i++) {
			sexes[i] = sampleData.getSexForIndividual(samples[i]);
		}
		return sexes;
	}

	/**
	 * Computes z-score to compare female and male intensity data means
	 * <p>
	 * Sex coding must be 1=male, 2=female
	 * <p>
	 * boolean[] samplesToExclude can be null
	 *
	 * @author John Lane
	 */
	public static double getSexZscore(int[] sexes, float[] independantData,
																		boolean[] samplesToExclude, Logger log) {
		double zscore = Double.NaN;
		DoubleVector[] values = new DoubleVector[3];
		for (int s = 0; s < sexes.length; s++) {
			if (!Double.isNaN(independantData[s])	&& sexes != null && (sexes[s] == 1 || sexes[s] == 2)
					&& (samplesToExclude == null || !samplesToExclude[s])) {
				if (values[sexes[s]] == null) {
					values[sexes[s]] = new DoubleVector();
				}
				values[sexes[s]].add((double) independantData[s]);
			}
		}
		if (values[1] != null && values[2] != null) {
			double[] maleValues = Doubles.toArray(values[1]);
			zscore = (Array.mean(maleValues) - Array.mean(Doubles.toArray(values[2])))
								/ Array.stdev(maleValues, false);
		}
		return zscore;
	}

	public static void lrrVariance(	Project proj, boolean[] samplesToInclude,
																	String markersToInclude) {
		PrintWriter writer;
		float[] lrrs, bafs;
		boolean[] useBAFs;
		MarkerData markerData;
		String markerName;
		long time;
		MarkerDataLoader markerDataLoader;
		String[] markerNames;
		String line, eol;
		Logger log;

		log = proj.getLog();
		if (Files.isWindows()) {
			eol = "\r\n";
		} else {
			eol = "\n";
		}

		try {
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()
																							+ "lrrVariance.xln"));
			writer.println(Array.toStr(LRR_VARIANCE_HEADER));

			if (markersToInclude != null) {
				markerNames =
										HashVec.loadFileToStringArray(proj.PROJECT_DIRECTORY.getValue(false, true)
																									+ markersToInclude, false, new int[] {0}, false);
			} else {
				markerNames = proj.getMarkerNames();
			}
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
			line = "";
			time = new Date().getTime();
			for (int i = 0; i < markerNames.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);

				markerName = markerData.getMarkerName();
				lrrs = markerData.getLRRs();
				bafs = markerData.getBAFs();

				line += markerName + "\t" + markerData.getChr() + "\t" + markerData.getPosition();

				if (lrrs == null) {
					log.reportError("Error - null lrr array for marker " + markerName);
				}
				lrrs = Array.removeNaN(Array.subArray(lrrs, samplesToInclude));
				line += "\t" + Array.stdev(lrrs, true);
				for (int j = 0; j < lrrs.length; j++) {
					lrrs[j] = Math.abs(lrrs[j]);
				}
				line += "\t" + Array.mean(lrrs);

				useBAFs = Arrays.copyOf(samplesToInclude, samplesToInclude.length);
				for (int j = 0; j < bafs.length; j++) {
					if (bafs[j] < 0.15 || bafs[j] > 0.85) {
						useBAFs[j] = false;
					}
				}
				bafs = Array.removeNaN(Array.subArray(bafs, useBAFs));
				line += "\t" + Array.stdev(bafs, true);
				for (int j = 0; j < bafs.length; j++) {
					bafs[j] = Math.abs(bafs[j] - 0.50f);
				}
				line += "\t" + Array.mean(bafs);

				line += eol;

				if (line.length() > 25000) {
					writer.print(line);
					writer.flush();
					line = "";
				}
				markerDataLoader.releaseIndex(i);
			}
			writer.print(line);
			writer.close();
			log.report("Finished analyzing " + markerNames.length + " in " + ext.getTimeElapsed(time));
		} catch (Exception e) {
			log.reportError("Error writing results");
			log.reportException(e);
		}
	}

	public static void separationOfSexes(Project proj, String subset) {
		PrintWriter writer;
		String[] samples;
		float[] xs, ys;
		MarkerData markerData;
		SampleData sampleData;
		int[] sexes;
		byte[] abGenotypes;
		String markerName;
		ClusterFilterCollection clusterFilterCollection;
		float gcThreshold;
		long time;
		DoubleVector[][][] values; // x/y, genotype, sex
		double femaleComp;
		double[][] zScores, tVals, counts;
		double[] maleValues, femaleValues;
		double zMean, tMean, zMin, tMin, count;
		double maf, countBs;
		MarkerDataLoader markerDataLoader;
		String[] markerList;
		String line, eol;
		Logger log;

		if (Files.isWindows()) {
			eol = "\r\n";
		} else {
			eol = "\n";
		}

		log = proj.getLog();
		sampleData = proj.getSampleData(2, false);
		samples = proj.getSamples();
		sexes = new int[samples.length];
		for (int i = 0; i < samples.length; i++) {
			sexes[i] = Math.max(0, sampleData.getSexForIndividual(samples[i]));
		}

		clusterFilterCollection = proj.getClusterFilterCollection();
		// gcThreshold = Float.parseFloat(proj.getProperty(proj.GC_THRESHOLD));
		gcThreshold = proj.getProperty(proj.GC_THRESHOLD).floatValue();

		try {
			writer = new PrintWriter(new FileWriter(proj.RESULTS_DIRECTORY.getValue(false, true)
																							+ "markerGenderChecks.xln"));
			writer.println("SNP\tX_zAA\tY_zBB\tX_tAA\tY_tBB\tMean_Zs\tMean_Ts\tMin_Z\tMin_T\tMAF");

			if (subset != null) {
				markerList = HashVec.loadFileToStringArray(proj.PROJECT_DIRECTORY.getValue(false, true)
																										+ subset, false, new int[] {0}, false);
			} else {
				markerList = proj.getMarkerNames();
			}
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerList);
			time = new Date().getTime();
			line = "";
			for (int i = 0; i < markerList.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);

				markerName = markerData.getMarkerName();
				xs = markerData.getXs();
				ys = markerData.getYs();
				abGenotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerName,
																														gcThreshold, log);

				values = new DoubleVector[3][4][3]; // x/y/r, genotype, sex
				for (int s = 0; s < samples.length; s++) {
					if (ext.isValidDouble(xs[s] + "")) {
						if (values[0][abGenotypes[s] + 1][sexes[s]] == null) {
							values[0][abGenotypes[s] + 1][sexes[s]] = new DoubleVector();
							values[1][abGenotypes[s] + 1][sexes[s]] = new DoubleVector();
							values[2][abGenotypes[s] + 1][sexes[s]] = new DoubleVector();
						}
						values[0][abGenotypes[s] + 1][sexes[s]].add((double) xs[s]);
						values[1][abGenotypes[s] + 1][sexes[s]].add((double) ys[s]);
						values[2][abGenotypes[s] + 1][sexes[s]].add((double) xs[s] + ys[s]);
					}
				}

				zScores = new double[3][3];
				tVals = new double[3][3];
				counts = new double[3][3];
				for (int k = 0; k < values.length; k++) {
					for (int ab = 1; ab < values[k].length; ab++) {
						if (values[k][ab][1] != null) {
							maleValues = Doubles.toArray(values[k][ab][1]);
							counts[k][ab - 1] += maleValues.length;
							// log.report(markerName+"\t"+FullSample.AB_PAIRS[ab-1]+"/"+(k==0?"X":"Y")+"\t"+maleValues.length+"
							// males\t"+Array.mean(maleValues));

							// use all the females with a missing genotype for comparison, if any females have
							// this genotype, only use that group if they outnumber those with a missing genotype
							femaleValues = null;
							if (values[k][ab][2] != null) {
								femaleValues = Doubles.toArray(values[k][ab][2]);
								counts[k][ab - 1] += femaleValues.length;
							}
							if (values[k][0][2] != null
									&& (femaleValues == null || values[k][0][2].size() > femaleValues.length)) {
								femaleValues = Doubles.toArray(values[k][0][2]);
							}

							if (femaleValues == null) {
								femaleComp = 0;
								tVals[k][ab - 1] = Double.NaN;
							} else {
								// log.report(markerName+"\t"+FullSample.AB_PAIRS[ab-1]+"/"+(k==0?"X":"Y")+"\t"+femaleValues.length+"
								// females\t"+Array.mean(femaleValues));
								femaleComp = Array.mean(femaleValues);
								tVals[k][ab - 1] = Math.abs(new Ttest(maleValues, femaleValues).getT());
							}
							zScores[k][ab - 1] = (Array.mean(maleValues) - femaleComp)
																		/ Array.stdev(maleValues, false);
						} else {
							zScores[k][ab - 1] = Double.NaN;
							tVals[k][ab - 1] = Double.NaN;
						}
					}
				}

				// int[][] interestedPairs = new int[][] {{0,0}, {1,2}, {2,1}};
				int[][] interestedPairs = new int[][] {{0, 0}, {1, 2}};
				zMean = count = 0;
				zMin = Double.MAX_VALUE;
				line += markerName;
				for (int[] interestedPair : interestedPairs) {
					if (Double.isNaN(zScores[interestedPair[0]][interestedPair[1]])) {
						line += "\t.";
					} else {
						line += "\t" + zScores[interestedPair[0]][interestedPair[1]];
						zMean += zScores[interestedPair[0]][interestedPair[1]]
											* counts[interestedPair[0]][interestedPair[1]];
						zMin = Math.min(zMin, zScores[interestedPair[0]][interestedPair[1]]);
						count += counts[interestedPair[0]][interestedPair[1]];
					}
				}
				zMean /= count;


				tMean = count = 0;
				tMin = Double.MAX_VALUE;
				for (int[] interestedPair : interestedPairs) {
					if (Double.isNaN(tVals[interestedPair[0]][interestedPair[1]])) {
						line += "\t.";
					} else {
						line += "\t" + tVals[interestedPair[0]][interestedPair[1]];
						tMean += tVals[interestedPair[0]][interestedPair[1]]
											* counts[interestedPair[0]][interestedPair[1]];
						tMin = Math.min(tMin, tVals[interestedPair[0]][interestedPair[1]]);
						count += counts[interestedPair[0]][interestedPair[1]];
					}
				}
				tMean /= count;

				line += "\t" + zMean + "\t" + tMean + "\t" + zMin + "\t" + tMin;

				count = countBs = 0;
				for (byte abGenotype : abGenotypes) {
					if (abGenotype >= 0) {
						countBs += abGenotype;
						count++;
					}
				}
				maf = countBs / (2 * count);
				if (maf > 0.5) {
					maf = 1 - maf;
				}
				line += "\t" + maf + eol;
				if (line.length() > 25000) {
					writer.print(line);
					writer.flush();
					line = "";
				}
			}
			writer.print(line);
			log.report("Finished analyzing " + markerList.length + " in " + ext.getTimeElapsed(time));

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing results");
			log.reportException(e);
		}

	}

	public static void filterMetrics(Project proj) {
		String markerMetricsFilename, reviewCriteriaFilename, exclusionCriteriaFilename,
				combinedCriteriaFilename;
		Logger log;

		log = proj.getLog();
		markerMetricsFilename = proj.MARKER_METRICS_FILENAME.getValue(false, false);
		if (!Files.exists(markerMetricsFilename)) {
			log.reportError("Error - marker metrics file not found at " + markerMetricsFilename);
			return;
		}

		reviewCriteriaFilename = proj.MARKER_REVIEW_CRITERIA_FILENAME.getValue(false, false);
		if (Files.exists(reviewCriteriaFilename)) {
			log.report("Using " + reviewCriteriaFilename + " for the review criteria");
		} else {
			log.report("Could not find "	+ reviewCriteriaFilename
									+ ", so generating from default parameters");
			Files.copyFileFromJar(DEFAULT_REVIEW_CRITERIA, reviewCriteriaFilename);
		}

		exclusionCriteriaFilename = proj.MARKER_EXCLUSION_CRITERIA_FILENAME.getValue(false, false);
		if (Files.exists(exclusionCriteriaFilename)) {
			log.report("Using " + exclusionCriteriaFilename + " for the review criteria");
		} else {
			log.report("Could not find "	+ reviewCriteriaFilename
									+ ", so generating from default parameters");
			Files.copyFileFromJar(DEFAULT_EXCLUSION_CRITERIA, exclusionCriteriaFilename);
		}

		combinedCriteriaFilename = proj.MARKER_COMBINED_CRITERIA_FILENAME.getValue(false, false);
		if (Files.exists(combinedCriteriaFilename)) {
			log.report("Using " + combinedCriteriaFilename + " for the review criteria");
		} else {
			log.report("Could not find "	+ combinedCriteriaFilename
									+ ", so generating from default parameters");
			Files.copyFileFromJar(DEFAULT_COMBINED_CRITERIA, combinedCriteriaFilename);
		}

		FilterDB.filter(markerMetricsFilename, reviewCriteriaFilename,
										proj.RESULTS_DIRECTORY.getValue(false, true) + "markersToReview.out", log);
		FilterDB.filter(markerMetricsFilename, exclusionCriteriaFilename,
										proj.RESULTS_DIRECTORY.getValue(false, true) + "markersToExclude.out", log);
		FilterDB.filter(markerMetricsFilename, combinedCriteriaFilename,
										proj.RESULTS_DIRECTORY.getValue(false, true) + "markersToReviewCombined.out",
										log);
	}

	public static void tallyFlaggedReviewedChangedAndDropped(	Project proj,
																														boolean checkForDeletedMarkers) {
		BufferedReader reader;
		PrintWriter writer, writerMissed;
		String[] line;
		Hashtable<String, Vector<String>> warningHash;
		Hashtable<String, String> flaggedMarkers, droppedMarkers, annotatedMarkers, reclusteredMarkers;
		HashSet<String> allOtherMarkers;
		Vector<String> v;
		String[] filenames;
		String dir;
		String[] header, expectedHeader, warnings, markerNames;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		ClusterFilterCollection clusterFilterCollection;
		AnnotationCollection annotationCollection;
		byte[] genotypes;
		boolean zeroedOut;
		float gcThreshold;
		String[] warningKeys;
		int numAnnotated, numReclustered, numDropped;
		String warning;
		boolean problem;
		String missedOutputFile;
		boolean annotated;
		char[] annotationKeys;
		Hashtable<String, HashSet<String>> warningHashHash;
		int[] warningCounts;
		String[] markersWithAnnotation;
		boolean[] shouldBeExcluded;
		Logger log;

		log = proj.getLog();
		dir = proj.PROJECT_DIRECTORY.getValue();
		// gcThreshold = proj.getFloat(proj.GC_THRESHOLD);
		gcThreshold = proj.GC_THRESHOLD.getValue().floatValue();
		shouldBeExcluded = proj.getSamplesToExclude();

		filenames = new String[] {"results/markersToExclude.out", "results/markersToReview.out"};
		filenames = new String[] {"results/markersToBoth.out"};
		filenames = new String[] {"results/markersToReviewOrFlaggedByQC.out"};
		filenames = new String[] {"results/markersToReview.out", "results/markersFlaggedByQC.out",
															"results/markersToReviewOrFlaggedByQC.out",
															"results/unrelatedWhiteFlags.out"};
		filenames = new String[] {"results/unrelatedWhiteFlags.out"};
		filenames = new String[] {"results/markersToReviewOrUnrelatedWhiteFlags.out"};
		filenames = new String[] {"results/markersToReview.out", "results/unrelatedWhiteFlags.out",
															"results/markersToReviewOrUnrelatedWhiteFlags.out"};
		filenames = new String[] {"results/unrelatedWhiteFlags_woCallRate.out"};
		filenames = new String[] {"results/markersToReviewCombined.out"};



		expectedHeader = new String[] {"Unit", "ReasonFlagged"};

		problem = false;
		for (String filename : filenames) {
			if (Files.exists(dir + filename)) {
				header = Files.getHeaderOfFile(dir + filename, log);
				if (!ext.checkHeader(header, expectedHeader, new int[] {0, 1}, false, log, false)) {
					problem = true;
				}
			} else {
				log.reportError("Error - could not find " + dir + filename);
				problem = true;
			}
		}
		if (problem) {
			return;
		}

		clusterFilterCollection = proj.getClusterFilterCollection();
		annotationCollection = proj.getAnnotationCollection();

		flaggedMarkers = new Hashtable<String, String>();
		for (String filename : filenames) {
			v = new Vector<String>();
			warningHash = new Hashtable<String, Vector<String>>();
			try {
				reader = Files.getAppropriateReader(dir + filename);
				if (reader == null) {
					log.reportError("Error - could not find " + dir + filename);
				}
				header = reader.readLine().trim().split("\t");
				ext.checkHeader(header, expectedHeader, new int[] {0, 1}, false, log, false);
				while (reader.ready()) {
					line = reader.readLine().trim().split("\t");
					if (!flaggedMarkers.containsKey(line[0])) {
						flaggedMarkers.put(line[0], "flaggedButNotAnnotated");
					}
					flaggedMarkers.put(line[0], flaggedMarkers.get(line[0]) + ";" + line[1]);
					warnings = line[1].split(";");
					for (String warning2 : warnings) {
						warning = warning2.substring(0, warning2.indexOf(" (") > 0	? warning2.indexOf(" (")
																																				: warning2.length());
						if (!warningHash.containsKey(warning)) {
							warningHash.put(warning, new Vector<String>());
						}
						warningHash.get(warning).add(line[0]);
						HashVec.addIfAbsent(line[0], v);
					}

				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + dir + filename + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + dir + filename + "\"");
				return;
			}

			markerNames = Array.toStringArray(v);
			log.report(markerNames.length + " markers met at least one criterion in " + filename);
			if (checkForDeletedMarkers) {
				markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(	proj,
																																										markerNames);
			} else {
				markerDataLoader = null;
			}
			annotatedMarkers = new Hashtable<String, String>();
			reclusteredMarkers = new Hashtable<String, String>();
			droppedMarkers = new Hashtable<String, String>();
			for (int j = 0; j < markerNames.length; j++) {
				if (j % 1000 == 0) {
					log.report((j + 1) + " of " + markerNames.length);
				}
				if (checkForDeletedMarkers) {
					markerData = markerDataLoader.requestMarkerData(j);
					zeroedOut = true;
					genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerNames[j],
																														gcThreshold, log);
					for (int k = 0; k < genotypes.length; k++) {
						if (!shouldBeExcluded[k] && genotypes[k] != -1) {
							zeroedOut = false;
						}
					}
				} else {
					zeroedOut = false;
				}
				if (zeroedOut) {
					droppedMarkers.put(markerNames[j], "");
				}
				if (annotationCollection != null
						&& annotationCollection.markerHasAnyAnnotation(markerNames[j])) {
					annotatedMarkers.put(markerNames[j], "");
				}
				if (clusterFilterCollection != null
						&& clusterFilterCollection.getClusterFilters(markerNames[j]) != null) {
					reclusteredMarkers.put(markerNames[j], "");
				}
				if (checkForDeletedMarkers) {
					markerDataLoader.releaseIndex(j);
				}
			}

			warningKeys = HashVec.getKeys(warningHash);
			try {
				writer = new PrintWriter(new FileWriter(dir + ext.rootOf(filename, false) + "_counts.out"));
				writer.println("WarningCriterion\tnumMarkers\tnumAnnotated\tnumReclustered\tnumDropped");
				for (String warningKey : warningKeys) {
					v = warningHash.get(warningKey);
					numReclustered = numAnnotated = numDropped = 0;
					for (int k = 0; k < v.size(); k++) {
						if (reclusteredMarkers.containsKey(v.elementAt(k))) {
							numReclustered++;
						}
						if (annotatedMarkers.containsKey(v.elementAt(k))) {
							numAnnotated++;
						}
						if (droppedMarkers.containsKey(v.elementAt(k))) {
							numDropped++;
						}
					}
					writer.println(warningKey	+ "\t" + v.size() + "\t" + numAnnotated + "\t" + numReclustered
													+ "\t" + numDropped);
				}
				numReclustered = numAnnotated = numDropped = 0;
				for (String markerName : markerNames) {
					if (reclusteredMarkers.containsKey(markerName)) {
						numReclustered++;
					}
					if (annotatedMarkers.containsKey(markerName)) {
						numAnnotated++;
					}
					if (droppedMarkers.containsKey(markerName)) {
						numDropped++;
					}
				}
				writer.println("Any criteria\t"	+ markerNames.length + "\t" + numAnnotated + "\t"
												+ numReclustered + "\t" + numDropped);

				missedOutputFile = dir + filename + "_missed.out";
				writerMissed = new PrintWriter(new FileWriter(missedOutputFile));
				allOtherMarkers = HashVec.loadToHashSet(proj.getMarkerNames());
				for (int j = 0; j < markerNames.length; j++) {
					if (!annotatedMarkers.containsKey(markerNames[j])) {
						writerMissed.println(markerNames[j] + "\t" + flaggedMarkers.get(markerNames[j]));
					}
					allOtherMarkers.remove(markerNames[j]);
				}
				numReclustered = numAnnotated = numDropped = 0;
				Set<String> toRemove = new HashSet<String>();
				writer.print("Everything else\t" + allOtherMarkers.size());

				for (String markerName : allOtherMarkers) {
					if (annotationCollection != null
							&& annotationCollection.markerHasAnyAnnotation(markerName)) {
						numAnnotated++;
					}
					if (clusterFilterCollection != null
							&& clusterFilterCollection.getClusterFilters(markerName) != null) {
						numReclustered++;
					} else {
						toRemove.add(markerName);
					}
				}

				allOtherMarkers.removeAll(toRemove);

				markerNames = allOtherMarkers.toArray(new String[allOtherMarkers.size()]);
				if (checkForDeletedMarkers) {
					markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(	proj,
																																											markerNames);
				}
				for (int j = 0; j < markerNames.length; j++) {
					if (checkForDeletedMarkers) {
						markerData = markerDataLoader.requestMarkerData(j);
						zeroedOut = true;
						genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection,
																															markerNames[j], gcThreshold, log);
						for (int k = 0; k < genotypes.length; k++) {
							if (!shouldBeExcluded[k] && genotypes[k] != -1) {
								zeroedOut = false;
							}
						}
					} else {
						zeroedOut = false;
					}
					if (annotationCollection != null
							&& annotationCollection.markerHasAnyAnnotation(markerNames[j])) {
						annotated = true;
					} else {
						annotated = false;
					}

					if (zeroedOut) {
						if (!annotated) {
							writerMissed.println(markerNames[j] + "\tdroppedButUnannotated");
						}
						numDropped++;
					} else {
						if (!annotated) {
							writerMissed.println(markerNames[j] + "\treclusteredButUnannotated");
						}
					}
					if (checkForDeletedMarkers) {
						markerDataLoader.releaseIndex(j);
					}
				}
				writer.println("\t" + numAnnotated + "\t" + numReclustered + "\t" + numDropped);
				writer.close();
				writerMissed.close();
				if (new File(missedOutputFile).length() == 0) {
					new File(missedOutputFile).delete();
				}
			} catch (Exception e) {
				log.reportError("Error writing to " + dir + ext.rootOf(filename, false) + "_counts.out");
				log.reportException(e);
			}

			warningHashHash = new Hashtable<String, HashSet<String>>();
			warningCounts = new int[warningKeys.length];
			for (int k = 0; k < warningKeys.length; k++) {
				String key = warningKeys[k];
				warningHashHash.put(key, HashVec.loadToHashSet(Array.toStringArray(warningHash.get(key))));
				warningCounts[k] = warningHashHash.get(key).size();
			}
			try {
				writer = new PrintWriter(new FileWriter(dir + ext.rootOf(filename, false) + "_matrix.out"));
				writer.println("\tN=\t" + Array.toStr(warningCounts));
				writer.println("Annotation\tCounts\t" + Array.toStr(warningKeys));

				annotationKeys = annotationCollection.getKeys();
				for (char annotationKey : annotationKeys) {
					markersWithAnnotation = annotationCollection.getMarkerLists(annotationKey);

					warningCounts = new int[warningKeys.length];
					for (String element : markersWithAnnotation) {
						for (int w = 0; w < warningKeys.length; w++) {
							if (warningHashHash.get(warningKeys[w]).contains(element)) {
								warningCounts[w]++;
							}
						}
					}
					writer.println(annotationCollection.getDescriptionForComment(annotationKey, false, false)
														+ "\t" + annotationCollection.getMarkerLists(annotationKey).length + "\t"
													+ Array.toStr(warningCounts));
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + dir + ext.rootOf(filename, false) + "_matrix.out");
				log.reportException(e);
			}
		}

		annotationReports(proj, checkForDeletedMarkers,
											proj.RESULTS_DIRECTORY.getValue(false, true) + "annotationReport.xln");
	}

	public static void tallyClusterFilters(	Project proj, boolean[] samplesToInclude,
																					String markersSubset) {
		PrintWriter writer;
		String[] markerNames;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		ClusterFilterCollection clusterFilterCollection;
		byte[] genotypesBefore, genotypesAfter;
		float gcThreshold;
		String filename;
		int numGenotypesAffected, numNonMissingBefore, numNonMissingAfter;
		Logger log;

		log = proj.getLog();
		// gcThreshold = proj.getFloat(proj.GC_THRESHOLD);
		gcThreshold = proj.GC_THRESHOLD.getValue().floatValue();
		clusterFilterCollection = proj.getClusterFilterCollection();
		if (samplesToInclude == null) {
			samplesToInclude = proj.getSamplesToInclude(null);
		}

		filename = proj.RESULTS_DIRECTORY.getValue(false, true) + "reclusteredMarkers.xln";
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("Marker\tChr\tPosition\tnumClusterFilters\tnumGenotypesAffected\tproportionGenotypesAffected\tCallrateBefore\tCallrateAfter\tCallrateChange\tmafBefore\tmafAfter");

			markerNames = clusterFilterCollection.getMarkerNames();
			log.report("Found " + markerNames.length + " markers with cluster filters");
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);

			for (int i = 0; i < markerNames.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);
				numGenotypesAffected = numNonMissingBefore = numNonMissingAfter = 0;
				genotypesBefore = markerData.getAbGenotypes();
				genotypesAfter = markerData.getAbGenotypesAfterFilters(	clusterFilterCollection,
																																markerNames[i], gcThreshold, log);
				for (int k = 0; k < genotypesBefore.length; k++) {
					if (samplesToInclude[k]) {
						if (genotypesBefore[k] != genotypesAfter[k]) {
							numGenotypesAffected++;
						}
						if (genotypesBefore[k] != -1) {
							numNonMissingBefore++;
						}
						if (genotypesAfter[k] != -1) {
							numNonMissingAfter++;
						}
					}
				}
				writer.println(markerNames[i]	+ "\t" + markerData.getChr() + "\t" + markerData.getPosition()
												+ "\t" + clusterFilterCollection.getClusterFilters(markerNames[i]).size()
												+ "\t" + numGenotypesAffected + "\t"
												+ ((double) numGenotypesAffected / (double) genotypesBefore.length) + "\t"
												+ ((double) numNonMissingBefore / (double) genotypesBefore.length) + "\t"
												+ ((double) numNonMissingAfter / (double) genotypesBefore.length) + "\t"
												+ (((double) numNonMissingAfter / (double) genotypesBefore.length)
														- ((double) numNonMissingBefore / (double) genotypesBefore.length))
												+ "\t" + (Array.mean(Array.removeAllValues(genotypesBefore, (byte) -1)) / 2)
												+ "\t"
												+ (Array.mean(Array.removeAllValues(genotypesAfter, (byte) -1)) / 2));
				markerDataLoader.releaseIndex(i);
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + filename);
			log.reportException(e);
		}
		log.report("Output successfully written to " + filename);
	}

	public static void annotationReports(	Project proj, boolean checkForDeletedMarkers,
																				String outputFile) {
		PrintWriter writer, writerMissed;
		Hashtable<String, String> reclusteredMarkers, droppedMarkers;
		HashSet<String> allOtherMarkers;
		String[] markerNames, markersWithAnnotation;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		ClusterFilterCollection clusterFilterCollection;
		AnnotationCollection annotationCollection;
		byte[] genotypes;
		boolean zeroedOut;
		float gcThreshold;
		int numReclustered, numDropped;
		char[] annotationKeys;
		String annotation;
		String missedOutputFile;
		boolean[] shouldBeExcluded;
		Logger log;

		log = proj.getLog();
		// gcThreshold = proj.getFloat(proj.GC_THRESHOLD);
		gcThreshold = proj.GC_THRESHOLD.getValue().floatValue();
		shouldBeExcluded = proj.getSamplesToExclude();

		clusterFilterCollection = proj.getClusterFilterCollection();
		annotationCollection = proj.getAnnotationCollection();

		markerNames = annotationCollection.getMarkerLists();
		log.report(markerNames.length + " markers have an annotation");
		if (checkForDeletedMarkers) {
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
		} else {
			markerDataLoader = null;
		}
		reclusteredMarkers = new Hashtable<String, String>();
		droppedMarkers = new Hashtable<String, String>();
		for (int j = 0; j < markerNames.length; j++) {
			if (j % 1000 == 0) {
				log.report((j + 1) + " of " + markerNames.length + " annotated markers");
			}
			if (checkForDeletedMarkers) {
				markerData = markerDataLoader.requestMarkerData(j);
				zeroedOut = true;
				genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerNames[j],
																													gcThreshold, log);
				for (int k = 0; k < genotypes.length; k++) {
					if (!shouldBeExcluded[k] && genotypes[k] != -1) {
						zeroedOut = false;
					}
				}
			} else {
				zeroedOut = false;
			}
			if (zeroedOut) {
				droppedMarkers.put(markerNames[j], "");
			}
			if (clusterFilterCollection != null
					&& clusterFilterCollection.getClusterFilters(markerNames[j]) != null) {
				reclusteredMarkers.put(markerNames[j], "");
			}
			if (checkForDeletedMarkers) {
				markerDataLoader.releaseIndex(j);
			}
		}

		annotationKeys = annotationCollection.getKeys();
		try {
			writer = new PrintWriter(new FileWriter(outputFile));
			writer.println("Annotation\tnumMarkers\tnumReclustered\tnumDropped");
			for (char annotationKey : annotationKeys) {
				annotation = annotationCollection.getDescriptionForComment(annotationKey, false, false);
				markersWithAnnotation = annotationCollection.getMarkerLists(annotationKey);
				numReclustered = numDropped = 0;
				for (String element : markersWithAnnotation) {
					if (reclusteredMarkers.containsKey(element)) {
						numReclustered++;
					}
					if (droppedMarkers.containsKey(element)) {
						numDropped++;
					}
				}
				writer.println(annotation	+ "\t" + markersWithAnnotation.length + "\t" + numReclustered
												+ "\t" + numDropped);
			}
			writer.println("Any annotation\t"	+ markerNames.length + "\t" + reclusteredMarkers.size()
											+ "\t" + droppedMarkers.size());
			Files.writeArray(HashVec.getKeys(droppedMarkers), proj.RESULTS_DIRECTORY.getValue(false, true)
																												+ "markers_that_were_dropped.out");

			allOtherMarkers = HashVec.loadToHashSet(proj.getMarkerNames());
			for (String markerName : markerNames) {
				allOtherMarkers.remove(markerName);
			}
			numReclustered = numDropped = 0;
			writer.print("Everything else\t" + allOtherMarkers.size());
			Files.writeArray(markerNames, proj.RESULTS_DIRECTORY.getValue(false, true)
																		+ "markers_not_yet_annotated.out");

			Set<String> toRemove = new HashSet<String>();
			for (String markerName : allOtherMarkers) {
				if (clusterFilterCollection != null
						&& clusterFilterCollection.getClusterFilters(markerName) != null) {
					numReclustered++;
				} else {
					toRemove.add(markerName);
				}
			}
			allOtherMarkers.removeAll(toRemove);

			markerNames = allOtherMarkers.toArray(new String[allOtherMarkers.size()]);
			if (checkForDeletedMarkers) {
				markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(	proj,
																																										markerNames);
			} else {
				markerDataLoader = null;
			}
			missedOutputFile = ext.addToRoot(outputFile, "_missed");
			writerMissed = new PrintWriter(new FileWriter(missedOutputFile));
			for (int j = 0; j < markerNames.length; j++) {
				if (checkForDeletedMarkers) {
					markerData = markerDataLoader.requestMarkerData(j);
					zeroedOut = true;
					genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerNames[j],
																														gcThreshold, log);
					for (int k = 0; k < genotypes.length; k++) {
						if (!shouldBeExcluded[k] && genotypes[k] != -1) {
							zeroedOut = false;
						}
					}
				} else {
					zeroedOut = false;
				}
				if (zeroedOut) {
					writerMissed.println(markerNames[j] + "\tdropped");
					numDropped++;
				} else {
					writerMissed.println(markerNames[j] + "\treclustered");
				}
				if (checkForDeletedMarkers) {
					markerDataLoader.releaseIndex(j);
				}
			}
			writer.println("\t" + numReclustered + "\t" + numDropped);
			writer.close();
			writerMissed.close();
			if (new File(missedOutputFile).length() == 0) {
				new File(missedOutputFile).delete();
			}
		} catch (Exception e) {
			log.reportError("Error writing to " + outputFile);
			log.reportException(e);
		}
	}

	public static void regress(Project proj, String phenotype) {
		PrintWriter writer;
		String trav;
		Hashtable<String, String> hash;
		String[] samples, header;
		double[] deps, indeps;
		String filename;
		int[] indices;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		String[] markerNames;
		RegressionModel model;
		boolean binary;
		Logger log;

		log = proj.getLog();
		// filename = proj.getFilename(proj.SAMPLE_DATA_FILENAME);
		filename = proj.SAMPLE_DATA_FILENAME.getValue();
		header = Files.getHeaderOfFile(filename, log);
		indices = ext.indexFactors(new String[] {"DNA", phenotype}, header, false, true);
		hash = HashVec.loadFileToHashString(filename, new int[] {indices[0]}, new int[] {indices[1]},
																				filename.endsWith(".csv"), null, true,
																				proj.JAR_STATUS.getValue(), false);

		samples = proj.getSamples();
		deps = new double[samples.length];
		for (int i = 0; i < samples.length; i++) {
			trav = hash.get(samples[i]);
			deps[i] = trav == null || ext.isMissingValue(trav) ? Double.NaN : Double.parseDouble(trav);
		}
		binary = RegressionModel.isBinaryTrait(Array.toStringArray(deps), log);

		markerNames = proj.getMarkerNames();
		try {
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()
																								+ ext.replaceWithLinuxSafeCharacters(phenotype, true)
																							+ "_regress.xln"));
			writer.println("MarkerName\tBAF_p");
			// markerDataLoader = new MarkerDataLoader(proj, markerNames, -1, log);
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
			for (int i = 0; i < markerNames.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);
				indeps = Array.toDoubleArray(markerData.getBAFs());
				if (binary) {
					model = new LogisticRegression(deps, indeps);
				} else {
					model = new LeastSquares(deps, indeps);
				}
				writer.println(markerNames[i] + "\t" + model.getSigs()[1]);
				if (i < 5) {
					log.report(model.getSummary());
				}
				markerDataLoader.releaseIndex(i);
			}

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to "	+ proj.PROJECT_DIRECTORY.getValue()
											+ ext.replaceWithLinuxSafeCharacters(phenotype, true) + "_regress.xln");
			log.reportException(e);
		}
	}

	/**
	 * @author lane0212 Computes Qc metrics for a subset of markers
	 */
	private static class MarkerMetricsWorker implements Callable<Boolean> {
		private final Project proj;
		private final boolean[] samplesToExclude;
		private final String[] markerNames;
		private final String fullPathToOutput;
		private final boolean checkMendel;

		public MarkerMetricsWorker(	Project proj, boolean[] samplesToExclude, String[] markerNames,
																String fullPathToOutput, boolean checkMendel) {
			super();
			this.proj = proj;
			this.samplesToExclude = samplesToExclude;
			this.markerNames = markerNames;
			this.fullPathToOutput = fullPathToOutput;
			this.checkMendel = checkMendel;
		}

		@Override
		public Boolean call() throws Exception {

			fullQC(proj, samplesToExclude, markerNames, fullPathToOutput, checkMendel);
			if (Files.exists(fullPathToOutput)
					&& Files.countLines(fullPathToOutput, 1) == markerNames.length) {
				return true;
			} else {
				if (!Files.exists(fullPathToOutput)) {
					proj.getLog().reportError("Could not compute marker metrics on "
																				+ Thread.currentThread().toString());
					proj.getLog().reportFileNotFound(fullPathToOutput);
				} else {
					proj.getLog()
							.reportError("Found "	+ Files.countLines(fullPathToOutput, 1) + " markers in "
																+ fullPathToOutput + " but should have found "
																+ markerNames.length);
				}
				return false;
			}
		}
	}

	public static ExtProjectDataParser developParser(Project proj, String markerMetricsFile) {
		ProjectDataParserBuilder builder = new ProjectDataParserBuilder();
		builder.sampleBased(false);
		builder.requireAll(false);
		builder.dataKeyColumnName(FULL_QC_HEADER[0]);
		builder.numericDataTitles(FULL_QC_HEADER);
		builder.setInvalidNumericToNaN(true);

		try {
			ExtProjectDataParser parser = builder.build(proj, markerMetricsFile);
			parser.determineIndicesFromTitles();
			parser.loadData();
			return parser;
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String markersSubset = null;
		String samples = null;
		String logfile = null;
		Project proj = null;
		String filename = null;
		boolean sexSeparation = false;
		boolean fullQC = false;
		boolean filter = false;
		boolean lrrVariance = false;
		boolean tally = false;
		boolean checkForDeletedMarkers = true;
		String pheno = null;
		boolean countFilters = false;
		int numThreads = 1;
		boolean checkMendel = true;

		String usage = "\n"	+ "cnv.qc.MarkerMetrics requires 0-1 arguments\n"
										+ "   (1) project properties filename (i.e. proj="
										+ org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
										+ "   (2) filename of subset of samples to include (i.e. samples=" + samples
										+ " (default; if null, uses all samples except those marked in the \"Excluded\" column in SampleData.txt))\n"
										+ "   (3) filename of subset of markers to include / otherwise all markers (i.e. markers="
										+ markersSubset + " (default))\n" + PSF.Ext.getNumThreadsCommand(4, numThreads)
										+

										"  AND\n"
										+ "   (4) look at intensity separation between males and females (i.e. -separation (not the default))\n"
										+ "  OR\n"
										+ "   (4) perform full list of checks on marker quality (i.e. -fullQC (not the default))\n"
										+ "   (5) do not check for mendelian errors when performing full qc (i.e. -skipMendel (not the default))\n"
										+ "  OR\n"
										+ "   (4) filter markers based on filter criteria (i.e. -filter (not the default))\n"
										+ "  OR\n"
										+ "   (4) check variance of LRR to help determine regions of instability (i.e. -lrrVar (not the default))\n"
										+ "  OR\n"
										+ "   (4) tally the number of reviewed markers that were changed or dropped (i.e. -tally (not the default))\n"
										+ "   (5) check for deleted markers (i.e. checkForDeleted="
										+ checkForDeletedMarkers + " (default))\n" + "  OR\n"
										+ "   (3) list which markers were adjusted using a cluster filter and how many genotypes changed class (i.e. -countFilters (not the default))\n"
										+ "  OR\n"
										+ "   (2) variable name in SampleData.txt to use as the outcome variable in the regression analyses (i.e. pheno=Class=ExamplePheno (not the default))\n"
										+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("markers=")) {
				markersSubset = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("samples=")) {
				samples = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-separation")) {
				sexSeparation = true;
				numArgs--;
			} else if (arg.startsWith("-skipMendel")) {
				checkMendel = false;
				numArgs--;
			} else if (arg.startsWith("-fullQC")) {
				fullQC = true;
				numArgs--;
			} else if (arg.startsWith("-lrrVar")) {
				lrrVariance = true;
				numArgs--;
			} else if (arg.startsWith("-filter")) {
				filter = true;
				numArgs--;
			} else if (arg.startsWith("-tally")) {
				tally = true;
				numArgs--;
			} else if (arg.startsWith("-countFilters")) {
				countFilters = true;
				numArgs--;
			} else if (arg.startsWith("checkForDeleted=")) {
				checkForDeletedMarkers = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("pheno=")) {
				pheno = arg.substring(6);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
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
			// subset = "data/test.txt";
			// filename = "C:/workspace/Genvisis/projects/gedi_exome.properties";
			// filename = "/home/npankrat/projects/GEDI_exomeRAF.properties";
			// filename = "/home/npankrat/projects/exome_chip_win.properties";
			// fullQC = true;
			// filter = true;

			// filename = "C:/workspace/Genvisis/projects/SingaporeReplication.properties";
			// fullQC = true;
			// filter = true;

			// markersSubset = "";

			// filename = "/home/npankrat/projects/SingaporeReplication.properties";
			// filename = "/home/npankrat/projects/GEDI_exomeRAF.properties";
			// filename = "/home/npankrat/projects/BOSS.properties";
			// filename = "/home/npankrat/projects/SOL_Metabochip.properties";
			// filename = "/home/npankrat/projects/WinterHillsCombo.properties";
			// fullQC = true;
			// markersSubset = "nans.txt";
			// tally = true;
			// checkForDeletedMarkers = true;
			// pheno = "Class=BAF_Outliers";
			// countFilters = true;

			proj = new Project(filename, logfile, false);

			if (sexSeparation) {
				separationOfSexes(proj, markersSubset);
			}
			if (fullQC) {
				fullQC(	proj,
								(samples == null)	? proj.getSamplesToExclude()
																	: Array.booleanNegative(proj.getSamplesToInclude(samples)),
								markersSubset, checkMendel, numThreads);
			}
			if (lrrVariance) {
				lrrVariance(proj, proj.getSamplesToInclude(samples), markersSubset);
			}
			if (filter) {
				filterMetrics(proj);
			}
			if (tally) {
				tallyFlaggedReviewedChangedAndDropped(proj, checkForDeletedMarkers);
			}
			if (countFilters) {
				tallyClusterFilters(proj, proj.getSamplesToInclude(samples), markersSubset);
			}
			if (pheno != null) {
				regress(proj, pheno);
			}
		} catch (Exception e) {
			if (proj != null) {
				proj.getLog().reportException(e);
			} else {
				e.printStackTrace();
			}
		}
	}
}
