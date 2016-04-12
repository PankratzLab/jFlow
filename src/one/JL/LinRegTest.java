package one.JL;

import java.util.ArrayList;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.apache.log4j.helpers.Loader;

import common.Logger;

import stats.LeastSquares;
import stats.LeastSquares.LS_TYPE;

public class LinRegTest {

	private static void run() {

		Logger log = new Logger();
		int numSamps = 10000;
		int numVar = 100;
		double[] y = new double[numSamps];
		double[][] x2 = new double[numSamps][numVar];
		for (int i = 0; i < 10000; i++) {
			y[i] = Math.random();
			for (int j = 0; j < x2[i].length; j++) {
				x2[i][j] = Math.random();
			}

		}

		long time = System.currentTimeMillis();
		LeastSquares ls = new LeastSquares(y, x2, true, true, LS_TYPE.REGULAR);
		log.reportTimeElapsed(LS_TYPE.REGULAR.toString(), time);
		System.out.println(ls.getRsquare());
		time = System.currentTimeMillis();
		LeastSquares lssvd = new LeastSquares(y, x2, true, true, LS_TYPE.SVD);
		log.reportTimeElapsed(LS_TYPE.SVD.toString(), time);
		System.out.println(lssvd.getRsquare());
		time = System.currentTimeMillis();
		OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
		ols.newSampleData(y, x2);
		
		// ols.estimateRegressionParameters();
		System.out.println(ols.calculateRSquared());
		log.reportTimeElapsed(LS_TYPE.OLS.toString(), time);

		for (int i = 0; i < ls.getBetas().length; i++) {
			if (i < 10) {
				System.out.println(LS_TYPE.REGULAR + ": " + ls.getBetas()[i] + "\t" + LS_TYPE.SVD + ": " + lssvd.getBetas()[i] + "\tOLS: " + ols.estimateRegressionParameters()[i]);
			}
		}

	}

	public static void main(String[] args) {
		run();

	}
}
