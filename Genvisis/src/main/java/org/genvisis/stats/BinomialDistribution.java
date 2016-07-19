package org.genvisis.stats;

import java.io.*;

import org.genvisis.common.ext;

public class BinomialDistribution {
	public static final double DEFAULT_ALPHA = 0.05;
	public static final int DEFAULT_SIGFIGS = 3;
	public static final int NUM_STEPS = 10;
	public static final int CI_SYMMETRIC = 0;
	public static final int CI_SCORE_TEST = 1;
	public static final int CI_LOGLIKELIHOOD = 2;

	public static double fact(int num) {
		double d = 1;

		for (int i = 2; i<=num; i++) {
			d *= i;
		}

		return d;
	}

	public static double factUntil(int num, int stop) {
		double d = 1;

		for (int i = stop+1; i<=num; i++) {
			d *= i;
		}

		return d;
	}

	public static double probability(int y, int n, double pi) {
		double multiplier;
		if (y>n-y) {
			multiplier = factUntil(n, y)/fact(n-y);
		} else {
			multiplier = factUntil(n, n-y)/fact(y);
		}

		return multiplier*Math.pow(pi, y)*Math.pow(1-pi, n-y);
	}

	public static double probabilityGTE(int y, int n, double pi) {
		double p = 0;

		for (int i = y; i<=n; i++) {
			p += probability(i, n, pi);
		}

		return p;
	}

	public static double probabilityLTE(int y, int n, double pi) {
		double p = 0;

		for (int i = y; i>=0; i--) {
			p += probability(i, n, pi);
		}

		return p;
	}

	public static double midProbabilityGTE(int y, int n, double pi) {
		double p = 0;

		for (int i = y; i<=n; i++) {
			p += probability(i, n, pi);
			if (i==y) {
				p /= 2;
			}
		}

		return p;
	}

	public static double midProbabilityLTE(int y, int n, double pi) {
		double p = 0;

		for (int i = y; i>=0; i--) {
			p += probability(i, n, pi);
			if (i==y) {
				p /= 2;
			}
		}

		return p;
	}

	public static double mean(int n, double pi) {
		return (double)n*pi;
	}

	public static double stdev(int n, double pi) {
		return Math.sqrt((double)n*pi*(1-pi));
	}

	public static double stderr(double pi, int n) {
		return Math.sqrt(pi*(1-pi)/(double)n);
	}

	public static double scoreTest(int observed, int n, double pi0) {
		return ((double)observed/(double)n-pi0)/stderr(pi0, n);
	}

	public static double[] confidenceInterval(int observed, int n, int type) {
		return confidenceInterval(observed, n, DEFAULT_ALPHA, DEFAULT_SIGFIGS, type);
	}

	public static double[] confidenceInterval(int observed, int n, double alpha, int sigfigs, int type) {
		double pi, piHat, z, se;
		double low, high, convergence_limit, step, value, bestVal, target;
		int bestPos;
		double[] ci = new double[2];

		piHat = (double)observed/(double)n;
		z = ProbDist.NormDistReverse(alpha);
		se = stderr(piHat, n);

		if (type==CI_SYMMETRIC) {
			ci = new double[] {piHat-z*se, piHat+z*se};
			ci[0] = Double.parseDouble(ext.formDeci(ci[0], sigfigs));
			ci[1] = Double.parseDouble(ext.formDeci(ci[1], sigfigs));
		} else {
			switch (type) {
			case CI_SCORE_TEST:
				target = z;
				break;
			case CI_LOGLIKELIHOOD:
				target = z*z;
				break;
			default:
				System.err.println("Error - Unknown confidence interval type ("+type+"); using score test");
				type = CI_SCORE_TEST;
				target = z;
				break;
			}

			convergence_limit = Math.pow(0.1, sigfigs+1);
			for (int i = 0; i<2; i++) {
				low = i==0?0:piHat;
				high = i==0?piHat:1;

				while (high-low>convergence_limit) {
					// System.out.println(ext.formDeci(low, sigfigs) +"\t"+
					// ext.formDeci(high, sigfigs));
					bestPos = -1;
					bestVal = Double.POSITIVE_INFINITY;

					step = (high-low)/NUM_STEPS;
					for (int j = 0; j<NUM_STEPS; j++) {
						pi = low+j*step;
						switch (type) {
						case CI_LOGLIKELIHOOD:
							value = loglikelihoodTest(observed, n, pi);
							break;
						case CI_SCORE_TEST:
						default:
							value = scoreTest(observed, n, pi);
							break;
						}
						if (Math.abs(value-(i==0||type!=CI_SCORE_TEST?1:-1)*target)<Math.abs(bestVal-(i==0||type!=CI_SCORE_TEST?1:-1)*target)) {
							bestPos = j;
							bestVal = value;
						}

					}
					if (bestPos==-1) {
						System.err.println("Error - yo, you got serious problems of non-convergence");
						return new double[] {-1, -1};
					} else {
						low = low+(bestPos-1)*step;
						high = high-(NUM_STEPS-bestPos-1)*step;
					}
				}

				ci[i] = Double.parseDouble(ext.formDeci((high+low)/2, sigfigs));
			}
		}

		return ci;
	}

	public static double waldTest(int observed, int n, double pi0) {
		double piHat = (double)observed/(double)n;

		return Math.pow((piHat-pi0)/stderr(piHat, n), 2);
	}

	public static double loglikelihoodTest(int observed, int n, double pi0) {
		double piHat = (double)observed/(double)n;

		return -2*Math.log(probability(observed, n, pi0)/probability(observed, n, piHat));
	}

	public static void demo() {
		int obs, n;
		double piHat, pi0;
		double[] ci;

		// obs = 400;
		// n = 893;
		// pi0 = 0.5;
		// obs = 5;
		// n = 25;
		// pi0 = 0.5;
		obs = 150;
		n = 250;
		pi0 = 0.5;

		piHat = (double)obs/(double)n;

		System.out.println("The mean is "+mean(n, piHat)+" and the standard deviation is "+ext.formDeci(stdev(n, piHat), 3));
		System.out.println("The probability of "+obs+" out of "+n+" trials, given piHat = "+ext.formDeci(piHat, 5)+" is: "+ext.formDeci(probability(obs, n, piHat), DEFAULT_SIGFIGS));
		System.out.println();
		System.out.println("With n = "+n+" and pi0 = "+pi0+" the standard error would be "+ext.formDeci(stderr(pi0, n), 4));
		System.out.println("With "+obs+" observed, the score test would yield a z of "+ext.formDeci(scoreTest(obs, n, pi0), DEFAULT_SIGFIGS));
		System.out.println("This yields a two-sided p-value of "+ext.prettyP(ProbDist.NormDist(scoreTest(obs, n, pi0))));
		ci = confidenceInterval(obs, n, CI_SCORE_TEST);
		System.out.println("The confidence interval using the duality of the score test is ("+ci[0]+", "+ci[1]+")");
		System.out.println();
		System.out.println("The Wald test would yield a chi-square of "+waldTest(obs, n, pi0));
		System.out.println("This yields a two-sided p-value of "+ext.prettyP(ProbDist.ChiDist(waldTest(obs, n, pi0), 1)));
		ci = confidenceInterval(obs, n, CI_SYMMETRIC);
		System.out.println("The symmetric confidence interval (i.e. Wald test) is ("+ci[0]+", "+ci[1]+")");
		System.out.println();
		if (Double.isNaN(loglikelihoodTest(obs, n, pi0))) {
			System.err.println();
			System.err.println("Error - the counts involved are too large to calculate exact probabilities/log-likelihoof ratios");
			System.err.println("        on the upsaide, the Score/Wald tests should be sufficiently accurate");
		} else {
			System.out.println("The log-likelihood test would yield a chi-square of "+ext.formDeci(loglikelihoodTest(obs, n, pi0), DEFAULT_SIGFIGS));
			System.out.println("This yields a two-sided p-value of "+ext.prettyP(ProbDist.ChiDist(loglikelihoodTest(obs, n, pi0), 1)));
			ci = confidenceInterval(obs, n, CI_LOGLIKELIHOOD);
			System.out.println("The confidence interval using the duality of the loglikelihood test is ("+ci[0]+", "+ci[1]+")");
			System.out.println();
			System.out.println("Exact probability LTE "+obs+" equals "+ext.formDeci(probabilityLTE(obs, n, pi0), 5));
			System.out.println("Exact probability GTE "+obs+" equals "+ext.formDeci(probabilityGTE(obs, n, pi0), 5));
			System.out.println();
			System.out.println("Mid probability LTE "+obs+" equals "+ext.formDeci(midProbabilityLTE(obs, n, pi0), 5));
			System.out.println("Mid probability GTE "+obs+" equals "+ext.formDeci(midProbabilityGTE(obs, n, pi0), 5));
			System.out.println();
		}
		// System.out.println("y\tP(y)\tp-value\tmid p-value");
		// for (int i = 0; i<=10; i++) {
		// System.out.println(i+"\t"+ext.formDeci(probability(i, 10, 0.5), 3,
		// true)+"\t"+ext.formDeci(probabilityGTE(i, 10, 0.5), 3,
		// true)+"\t"+ext.formDeci(midProbabilityGTE(i, 10, 0.5), 3, true));
		// }
	}

	public static void main(String[] args) throws IOException {
		demo();
	}
}
