// -Xms1024M -Xmx1024M
package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GraphicsEnvironment;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import javax.swing.JOptionPane;

import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.SerialFloatArray;

public class QQPlot {
	public static final long serialVersionUID = 1L;

	public static final String[] DEFAULT_FILES = {
																								// "C:\\Documents and Settings\\npankrat\\My
																								// Documents\\Downloads\\allelic.txt"
																								// "C:\\Documents and Settings\\npankrat\\My
																								// Documents\\LOAD\\QQplots\\dominant.txt"
																								"C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\QQplots\\no_covariates.txt,0",
																								// "C:\\Documents and Settings\\npankrat\\My
																								// Documents\\LOAD\\QQplots\\4PCs.txt",
																								"C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\QQplots\\E2_E4.txt,0",
																							 // "C:\\Documents and Settings\\npankrat\\My
																							 // Documents\\LOAD\\QQplots\\4PCs_E2_E4.txt",
																							 // "C:\\Documents and Settings\\npankrat\\My
																							 // Documents\\LOAD\\QQplots\\E4_binary.txt"
	};

	public static final Color[] COLOR_SCHEME = new Color[] {Color.BLACK, Color.GRAY,
																													new Color(33, 31, 53), // dark dark
																													new Color(23, 58, 172), // dark blue
																													new Color(201, 30, 10), // deep red
																													new Color(140, 20, 180), // deep purple
																													new Color(33, 87, 0), // dark green
																													new Color(55, 129, 252), // light blue
																													new Color(94, 88, 214), // light purple
																													new Color(189, 243, 61), // light green
																													new Color(217, 109, 194), // pink
																													new Color(0, 0, 128), // ALL KINDS OF
																																								// BLUES
																													new Color(100, 149, 237),
																													new Color(72, 61, 139),
																													new Color(106, 90, 205),
																													new Color(123, 104, 238),
																													new Color(132, 112, 255),
																													new Color(0, 0, 205),
																													new Color(65, 105, 225),
																													new Color(0, 0, 255),
																													new Color(30, 144, 255),
																													new Color(0, 191, 255),
																													new Color(135, 206, 250),
																													new Color(135, 206, 250),
																													new Color(70, 130, 180),
																													new Color(176, 196, 222),
																													new Color(173, 216, 230),
																													new Color(176, 224, 230),
																													new Color(175, 238, 238),
																													new Color(0, 206, 209),
																													new Color(72, 209, 204),
																													new Color(64, 224, 208),
																													new Color(0, 255, 255),
																													new Color(224, 255, 255),
																													new Color(255, 192, 0), // yellowy orange
																													new Color(227, 108, 9), // halloween
																																									// orange
	};

	private final QQPanel qqPanel;
	private final Logger log;

	/**
	 *
	 * @param labels labels for each set of pvals
	 * @param pvals array of pvals for each label in labels
	 * @param log10 use -log10 of p values
	 * @param rotated display QQ rotated
	 * @param symmetric force symmetric axes
	 * @param maxValue maximum -log10 pval to display
	 * @param log
	 */
	public QQPlot(String[] labels, double[][] pvals, boolean[][] pvalsToUseForLambdaCalc,
								boolean log10, boolean rotated,
								boolean symmetric, float maxValue, Logger log) {
		this.log = log;
		log.report("Loading data for " + ext.listWithCommas(labels));

		qqPanel = new QQPanel(labels, pvals, pvalsToUseForLambdaCalc, log10, rotated, maxValue,
													COLOR_SCHEME, log);
		qqPanel.setSymmetricAxes(symmetric);
		if (!log10) {
			qqPanel.setForcePlotXmax(1);
			qqPanel.setForcePlotYmax(1);
		}

		qqPanel.setTrackedSize((int) qqPanel.getSize().getWidth(), (int) qqPanel.getSize().getHeight());
	}

	public void screenCap(String outFile) {

		int descrHeight = qqPanel.getPvals().length * 35;
		qqPanel.setSize(new Dimension(2000, 1440 + descrHeight));
		qqPanel.validate();
		qqPanel.createImage();
		qqPanel.screenCapture(outFile);

	}


	protected QQPanel getQqPanel() {
		return qqPanel;
	}

	protected Logger getLog() {
		return log;
	}

	/**
	 *
	 * @param filenames files to load pvals from
	 * @param plotLabel label for the plot
	 * @param displayQuantiles display quantiles instead of -log10 pvals
	 * @param displayStandardQQ display standard QQ plot with -log10 pvals
	 * @param displayRotatedQQ display rotated QQ plot with -log10 pvals
	 * @param maxToPlot -log10(p) at which to start truncating
	 * @param symmetric force symmetric axes
	 * @param maxValue maximum -log10 p-value to plot
	 * @param log
	 */
	public static QQPlot loadPvals(String[] filenames, String plotLabel, boolean displayQuantiles,
																 boolean displayStandardQQ, boolean displayRotatedQQ,
																 double maxToPlot, double mafLwrBnd, boolean symmetric,
																 float maxValue, Logger log) {
		BufferedReader reader;
		String[] labels;
		double[][] pvals;
		boolean[][] usePvalForLambda;
		int count;
		String fullLine;
		String[] lineParts;
		String trav;
		boolean error, header;
		int[] cols;
		int[] mafCols;
		double minPval;
		String delimiter;
		String temp;
		int invalids;

		if (filenames == null || filenames.length == 0) {
			reportQQError("There are no files selected for viewing in QQ plot; please add the names of the files you want to view after QQ_FILENAMES= in the project's .properties file (also be sure to uncomment the property by removing the hash symbol (\"#\"))",
										log);
			return null;
		}

		if (maxToPlot < 0) {
			minPval = -1;
		} else {
			minPval = 1 / Math.pow(10, maxToPlot);
		}

		error = false;
		labels = new String[filenames.length];
		pvals = new double[filenames.length][];
		usePvalForLambda = new boolean[filenames.length][];
		cols = ArrayUtils.intArray(filenames.length, 0);
		mafCols = ArrayUtils.intArray(filenames.length, -1);

		for (int i = 0; i < filenames.length; i++) {
			if (filenames[i].indexOf("=") > 0) {
				labels[i] = filenames[i].substring(filenames[i].lastIndexOf("=") + 1);
				filenames[i] = filenames[i].substring(0, filenames[i].lastIndexOf("="));
			} else {
				labels[i] = ext.rootOf(filenames[i]);
			}
			if (filenames[i].indexOf(",") > 0) {
				try {
					String[] inds = filenames[i].substring(filenames[i].indexOf(',') + 1).split(",");
					if (inds.length == 1) {
						cols[i] = Integer.parseInt(inds[0]);
						filenames[i] = filenames[i].substring(0, filenames[i].lastIndexOf(","));
						labels[i] += " '" + Files.getHeaderOfFile(filenames[i], log)[cols[i]] + "'";
						if (mafLwrBnd > 0) {
							String[] v = Files.getHeaderOfFile(filenames[i], log);
							int[] poss = ext.indexFactors(Aliases.ALLELE_FREQS, v, false);
							for (int p : poss) {
								if (p >= 0) {
									mafCols[i] = p;
									break;
								}
							}
						}
					} else if (inds.length == 2) {
						cols[i] = Integer.parseInt(inds[0]);
						if (mafLwrBnd > 0) {
							mafCols[i] = Integer.parseInt(inds[1]);
						}
						filenames[i] = filenames[i].substring(0, filenames[i].lastIndexOf(","));
						labels[i] += " '" + Files.getHeaderOfFile(filenames[i], log)[cols[i]] + "'";
					}
				} catch (Exception e) {
				}
			} else {
				String[] v = Files.getHeaderOfFile(filenames[i], log);
				int[] poss = ext.indexFactors(Aliases.PVALUES, v, false);
				for (int p : poss) {
					if (p >= 0) {
						cols[i] = p;
						labels[i] += " '" + v[cols[i]] + "'";
						break;
					}
				}
				if (mafLwrBnd > 0) {
					poss = ext.indexFactors(Aliases.ALLELE_FREQS, v, false);
					for (int p : poss) {
						if (p >= 0) {
							mafCols[i] = p;
							break;
						}
					}
				}
			}
			log.report("Loading " + labels[i]);
			if (!Files.exists(filenames[i])) {
				reportQQError("Error - file " + filenames[i] + " does not exist", log);
				return null;
			}
			delimiter = Files.determineDelimiter(filenames[i], log);
			try {
				reader = Files.getReader(filenames[i], true, true);
				count = 0;
				temp = reader.readLine();
				try {
					if (delimiter.equals(",")) {
						trav = ext.splitCommasIntelligently(temp, true, log)[cols[i]];
					} else {
						trav = temp.trim().split(delimiter, -1)[cols[i]];
					}
				} catch (Exception e) {
					log.reportError("Error - could not parse " + filenames[i] + " completely:");
					log.reportError(temp);
					log.reportException(e);
					return null;
				}
				invalids = 0;
				try {
					if (!ext.isMissingValue(trav)) {
						if (Double.parseDouble(trav) <= 0) {
							reportQQError("Error - one of the p-values in file " + filenames[i]
														+ " is near zero (" + trav + ") for line:\n"
														+ ext.replaceAllWith(temp, delimiter, "  "), log);
							invalids++;
						}
						header = false;
						count++;
					} else {
						header = false;
					}
				} catch (NumberFormatException nfe) {
					header = true;
				}
				while (reader.ready()) {
					temp = reader.readLine();
					if (delimiter.equals(",")) {
						trav = ext.splitCommasIntelligently(temp, true, log)[cols[i]];
					} else {
						trav = temp.trim().split(delimiter, -1)[cols[i]];
					}
					if (!ext.isMissingValue(trav)) {
						try {
							if (Double.parseDouble(trav) <= 0) {
								if (invalids < 3) {
									reportQQError("Error - one of the p-values in file " + filenames[i]
																+ " is near zero (" + trav + ") for line:\n"
																+ ext.replaceAllWith(temp, delimiter, "  "), log);
								}
								invalids++;
							} else {
								count++;
							}
						} catch (NumberFormatException nfe) {
							if (invalids < 3) {
								reportQQError("Error - one of the p-values in file " + filenames[i]
															+ " is not a number (" + trav + ") for line:\n"
															+ ext.replaceAllWith(temp, delimiter, "  "), log);
							}
							invalids++;
						}
					}
				}
				reader.close();
				if (invalids > 2) {
					reportQQError("There were " + invalids
												+ " total markers that had an invalid p-value for file " + filenames[i],
												log);
				}

				reader = Files.getReader(filenames[i], true, true);
				pvals[i] = new double[count];
				usePvalForLambda[i] = ArrayUtils.booleanArray(count, true);
				count = 0;
				if (header) {
					reader.readLine();
				}
				while (reader.ready()) {
					fullLine = reader.readLine();
					if (delimiter.equals(",")) {
						lineParts = ext.splitCommasIntelligently(fullLine, true, log);
					} else {
						lineParts = fullLine.trim().split(delimiter, -1);
					}
					trav = lineParts[cols[i]];
					if (!ext.isMissingValue(trav)) {
						try {
							pvals[i][count] = Double.parseDouble(trav);
							if (pvals[i][count] > 0) {
								if (pvals[i][count] < minPval) {
									pvals[i][count] = minPval;
								}
							}
						} catch (NumberFormatException nfe) {
						}
						if (mafLwrBnd > 0) {
							if (mafCols[i] >= 0) {
								trav = lineParts[mafCols[i]];
								if (!ext.isMissingValue(trav)) {
									try {
										usePvalForLambda[i][count] = Double.parseDouble(trav) > mafLwrBnd;
									} catch (NumberFormatException nfe) {
										usePvalForLambda[i][count] = false;
									}
								}
							}
						}
						count++;
					}
				}
				reader.close();
				if (count != pvals[i].length) {
					reportQQError("Error - mismatched number of values: " + count + " of " + pvals[i].length
												+ " were valid", log);
					return null;
				}

				Arrays.sort(pvals[i]);
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error - missing file: \"" + filenames[i] + "\"");
				error = false;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + filenames[i] + "\"");
				error = false;
			}
		}

		if (error) {
			return null;
		}

		if (displayQuantiles) {
			return new QQPlot(labels, pvals, usePvalForLambda, false, false, false, maxValue, log);
		}
		if (displayStandardQQ) {
			return new QQPlot(labels, pvals, usePvalForLambda, true, false, symmetric, maxValue, log);
		}
		if (displayRotatedQQ) {
			return new QQPlot(labels, pvals, usePvalForLambda, true, true, false, maxValue, log);
		}
		return null;
	}

	private static void reportQQError(String error, Logger log) {
		if (!GraphicsEnvironment.isHeadless()) {
			JOptionPane.showMessageDialog(null, error, "Error", JOptionPane.ERROR_MESSAGE);
		} else {
			log.reportError(error);
		}
	}

	public static void computeCI(String dir, String prefix, int max, Logger log) {
		PrintWriter writer;
		int count, length;
		float[] array;
		float[][] dists;

		if (max == -1) {
			max = Integer.MAX_VALUE;
		}
		count = 1;
		while (new File(dir + prefix + "." + count + ".results").exists() && count < max) {
			count++;
		}
		count--;
		log.report("Found " + (count == 0 ? "no" : count + "") + " replicates to process");
		dists = new float[0][0];
		length = -1;
		for (int i = 1; i <= count; i++) {
			if (i % 10 == 0) {
				log.report(i + "");
			}
			array = SerialFloatArray.load(dir + prefix + "." + i + ".results").getArray();
			Arrays.sort(array);
			if (i == 1) {
				length = array.length;
				dists = new float[length][count];
			} else if (array.length != length) {
				log.reportError("Error - mismatched number of p-values at rep " + i);
			}
			for (int j = 0; j < array.length; j++) {
				dists[j][i - 1] = array[j];
				if (j == array.length - 1 && Float.isNaN(dists[j][i - 1])) {
					log.reportError("Error - rep " + i + " has NaNs");
				}
			}
		}

		try {
			writer = Files.openAppropriateWriter(dir + "CIs.xln");
			writer.println("2.5%ile\t5%ile\t95%ile\t97.5%ile");
			for (float[] dist : dists) {
				writer.println(ArrayUtils.toStr(ArrayUtils.quants(dist, new double[] {0.025, 0.05, 0.95,
																																							0.975})));
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + "CIs.xln");
			log.reportException(e);
		}
	}

	public static void main(String[] args) {
		QQPlotFrame.main(args);
	}
}
