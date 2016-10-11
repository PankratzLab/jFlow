package org.genvisis.gwas;

import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class MetaAnalysis {
	public static String[] inverseVarianceWeighting(String[] lines, Logger log) {
		String[] line;
		Vector<double[]> v;
		int countInvalids;
		int numColumns;
		boolean parentheses, sdConversion;
		double[] result;
		int n, totalN;

		totalN = 0;
		parentheses = false;
		sdConversion = false;
		countInvalids = 0;
		lines = ext.getClipboard().trim().split("\\n");
		v = new Vector<double[]>();
		numColumns = lines[0].trim().split("[\\s]+").length;
		for (int i = 0; i < lines.length; i++) {
			line = lines[i].trim().split("[\\s]+");
			if (line.length == 1 && line[0].equals("")) {
				// do nothing
			} else {
				if (line.length < 2) {
					log.reportError("Line # "	+ (i + 1)
													+ " did not have at least two columns with beta and standard error data: "
													+ Array.toStr(line, " / "));
				} else if (line.length != numColumns) {
					log.reportError("Line # "	+ (i + 1) + " had a different number of columns (n="
													+ line.length + ") than the previous columns (n=" + numColumns + "): "
													+ Array.toStr(line, " / "));
				}
				if (line[1].startsWith("(") && line[1].endsWith(")")) {
					line[1] = line[1].substring(1, line[1].length() - 1);
					parentheses = true;
					if (line.length > 2) {
						log.reportError("Assuming the third column is an integer N to convert the SD to stderr");
						n = Integer.parseInt(line[2]);
						line[1] = (Double.parseDouble(line[1]) / Math.sqrt(n)) + "";
						totalN += n;
						sdConversion = true;
					}
				}
				if (ext.isValidDouble(line[0]) && ext.isValidDouble(line[1])) {
					v.add(new double[] {Double.parseDouble(line[0]), Double.parseDouble(line[1])});
				} else {
					if (countInvalids < 5) {
						log.reportError("Line # "	+ (i + 1) + " had an invalid double: "
														+ Array.toStr(line, " / "));
					}
					countInvalids++;
				}
			}
		}
		if (countInvalids >= 5) {
			log.reportError("There were " + countInvalids + " invalid doubles in the data");
		}

		result = inverseVarianceWeighting(Matrix.toDoubleArrays(v), log);
		if (sdConversion) {
			result[1] *= Math.sqrt(totalN);
		}
		if (parentheses) {
			return new String[] {ext.formDeci(result[0], 9)	+ " (" + ext.formDeci(result[1], 9) + ")"
														+ (sdConversion ? "\t" + totalN : "")};
		} else {
			return Array.toStringArray(result);
		}
	}

	public static double[] inverseVarianceWeighting(double[][] betaStderrs, Logger log) {
		double squaredWeight, weightedSum, sumInverseWeights;

		weightedSum = sumInverseWeights = 0;
		for (int i = 0; i < betaStderrs.length; i++) {
			if (Double.isNaN(betaStderrs[i][0]) || Double.isNaN(betaStderrs[i][1])) {
				log.reportError("Entry # "	+ (i + 1) + " had an invalid beta or stderr: "
												+ Array.toStr(Array.toStringArray(betaStderrs[i]), " / "));
			} else {
				squaredWeight = betaStderrs[i][1] * betaStderrs[i][1];
				weightedSum += betaStderrs[i][0] / squaredWeight;
				sumInverseWeights += 1 / squaredWeight;
			}
		}

		return new double[] {weightedSum / sumInverseWeights, Math.sqrt(1 / sumInverseWeights)};
	}
}
