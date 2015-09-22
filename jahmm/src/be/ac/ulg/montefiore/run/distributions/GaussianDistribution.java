/* jahmm package - v0.6.1 */

/*
  *  Copyright (c) 2004-2006, Jean-Marc Francois.
 *
 *  This file is part of Jahmm.
 *  Jahmm is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Jahmm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Jahmm; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */

package be.ac.ulg.montefiore.run.distributions;

import java.util.*;


/**
 * This class implements a Gaussian distribution.
 */
public class GaussianDistribution
implements RandomDistribution {
	
	private double mean;
	private double deviation;
	private double variance;
	private final static Random randomGenerator = new Random();
	
	
	/**
	 * Creates a new pseudo-random, Gaussian distribution with zero mean
	 * and unitary variance.
	 */
	public GaussianDistribution()
	{
		this(0., 1.);
	}
	
	
	/**
	 * Creates a new pseudo-random, Gaussian distribution.
	 *
	 * @param mean The mean value of the generated numbers.
	 * @param variance The variance of the generated numbers.
	 */
	public GaussianDistribution(double mean, double variance)
	{
		if (variance <= 0.)
			throw new IllegalArgumentException("Variance must be positive");
		
		this.mean = mean;
		this.variance = variance;
		this.deviation = Math.sqrt(variance);
	}
	
	
	/**
	 * Returns this distribution's mean value.
	 *
	 * @return This distribution's mean value.
	 */
	public double mean()
	{
		return mean;
	}
	
	
	/**
	 * Returns this distribution's variance.
	 *
	 * @return This distribution's variance.
	 */
	public double variance()
	{
		return variance;
	}
	
	
	public double generate()
	{
		return randomGenerator.nextGaussian() * deviation + mean;
	}
	
	
	public double probability(double n) {
		double expArg = -.5 * (n - mean) * (n - mean) / variance;
		return Math.pow(2. * Math.PI * variance, -.5) * Math.exp(expArg);
	}

	/**
	 * 
	 * @param n
	 * @return the cdf of this value using taylor aproximation, custom genvisis code.
	 * 
	 * Nabbed from http://introcs.cs.princeton.edu/java/22library/Gaussian.java.html
	 */
	public double cdf(double n) {
		double z = (n - mean) / Math.sqrt(variance);
		if (z < -8.0)
			return 0.0;
		if (z > 8.0)
			return 1.0;
		double sum = 0.0, term = z;
		for (int i = 3; sum + term != sum; i += 2) {
			sum = sum + term;
			term = term * z * z / i;
		}
		return 0.5 + sum * probability(n);
	}
	
	private static final long serialVersionUID = 9127329839769283975L;
	
	public static void main(String[] args) {
		// can use http://www.solvemymath.com/online_math_calculator/statistics/cdf_calculator.php to compare
		GaussianDistribution distribution = new GaussianDistribution(2, 20);
		System.out.println(distribution.probability(.5));
		System.out.println(distribution.cdf(.5));
		System.out.println(distribution.probability(0));
		System.out.println(distribution.cdf(0));

	}
}
