package org.genvisis.stats;

import org.genvisis.common.Array;
import org.genvisis.mining.Transformations;

public class Correlation {
	public static double[] Pearson(double[] data1, double[] data2) {
		return Pearson(new double[][] {data1, data2});
	}

	public static double[] Pearson(double[][] data) {
		double[][] zScores;
		double[] results = new double[] {-999, -999};
		double t;
		int count;
		
		if (data.length!=2) {
			System.err.println("Error - PearsonCorrelation requires a 2xN table to calculate a coefficient");
		} else if (data[0].length!=data[1].length) {
			System.err.println("Error - PearsonCorrelation requires the 2 arrays to have equal N");
		} else {
			zScores = new double[][] {Array.normalize(data[0]), Array.normalize(data[1])};
			results[0] = 0;
			count = 0;
			for (int i = 0; i<data[0].length; i++) {
				results[0] += zScores[0][i]*zScores[1][i];
				count++;
			}
			results[0] /= count-1;
			t = results[0]*Math.sqrt(data[0].length-2)/Math.sqrt(1-results[0]*results[0]);
			results[1] = ProbDist.TDist(t, data[0].length-2);
		}
	
		return results;
	}

	public static double[] Spearman(double[][] data) {
		return Pearson(new double[][] {Transformations.rankTransform(data[0]), Transformations.rankTransform(data[1])});
	}
}
