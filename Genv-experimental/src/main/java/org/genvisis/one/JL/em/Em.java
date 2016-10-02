package org.genvisis.one.JL.em;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.MixtureMultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization;
import org.apache.commons.math3.util.Pair;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

import com.google.common.primitives.Doubles;

import edu.stanford.facs.logicle.Logicle;

/**
 * Explore em for flow
 *
 */
public class Em {
	private Em() {

	}

	private static Logicle getBiexScale(double max) {
		double w = 2; // linear decades // W=1 is weird, W=3 -> "scale() didn't
						// converge"
		double a = 0; // negative decades
		double decCnt = count(max);

		return new Logicle(max, Math.min(w, decCnt / 2), decCnt, a);
	}

	private static double count(double max) {
		double dec = max;
		double decCnt = 0;
		while (dec > 1) {
			dec = dec / 10d;
			decCnt++;
		}
		return decCnt;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String export = "/Users/Kitty/temp/emflow/export.xln";
		String[][] data = HashVec.loadFileToStringMatrix(export, false, new int[] { 1, 2 }, false);
		Logger log = new Logger();
		log.reportTimeInfo(data.length + " by " + data[0].length);

		ArrayList<String> out = new ArrayList<String>();
		ArrayList<Double> xt = new ArrayList<Double>();
		ArrayList<Double> yt = new ArrayList<Double>();

		for (int i = 0; i < data.length; i++) {
			xt.add(Double.parseDouble(data[i][0]));
			yt.add(Double.parseDouble(data[i][1]));
		}

		double[] x = Doubles.toArray(xt);
		double[] y = Doubles.toArray(yt);

		Logicle mehx = getBiexScale(Array.max(x));
		Logicle mehy = getBiexScale(Array.max(y));

		for (int i = 0; i < y.length; i++) {
			x[i] = mehx.scale(x[i]);
			y[i] = mehy.scale(y[i]);
			out.add(x[i] + "\t" + y[i]);

		}

		double[][] m = new double[][] { x, y };
		MultivariateNormalMixtureExpectationMaximization mle = new MultivariateNormalMixtureExpectationMaximization(
				Matrix.transpose(m));
		MixtureMultivariateNormalDistribution initialMix = MultivariateNormalMixtureExpectationMaximization
				.estimate(Matrix.transpose(m), 3);
		mle.fit(initialMix, 500, 1E-5);
		MixtureMultivariateNormalDistribution finalMix = mle.getFittedModel();
		List<Pair<Double, MultivariateNormalDistribution>> dists = finalMix.getComponents();
		for (Pair<Double, MultivariateNormalDistribution> pair : dists) {
			log.reportTimeInfo("Means: " + Array.toStr(pair.getValue().getMeans()));
			log.reportTimeInfo("SDs: " + Array.toStr(pair.getValue().getStandardDeviations()));
		}

		Files.writeIterable(out, ext.parseDirectoryOfFile(export) + "export.biexp.txt");

	}
}
