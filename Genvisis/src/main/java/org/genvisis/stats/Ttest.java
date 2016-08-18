package org.genvisis.stats;

import java.io.*;
import java.util.*;

import org.genvisis.common.Array;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.ext;

import com.google.common.primitives.Doubles;

public class Ttest {
	private double meanDiff;
	private double stdev;
	private double[] t;
	private double[] p;
	private int df;
	private double v;
	private double[] ci05;
	private double[] ci025;
	private int testType;
	private double Ftest;
	private double Fprob;
	private boolean equal;

	public Ttest(double[][] data) { // paired t-test
		double[] diffs = new double[data.length];
		int n = data.length;

		t = new double[1];
		p = new double[1];
		ci05 = new double[1];
		ci025 = new double[1];
		meanDiff = 0;
		stdev = 0;

		for (int i = 0; i<n; i++) {
			diffs[i] = data[i][0]-data[i][1];
			meanDiff += diffs[i];
		}
		meanDiff /= n;

		for (int i = 0; i<n; i++) {
			stdev += (diffs[i]-meanDiff)*(diffs[i]-meanDiff);
		}
		stdev /= n-1;
		stdev = Math.sqrt(stdev);

		ci05[0] = ProbDist.TDistReverse(0.1, n-1)*stdev/Math.sqrt(n); // distribution
		// assumes
		// a
		// 2-tailed
		// test,
		// workaround
		ci025[0] = ProbDist.TDistReverse(0.05, n-1)*stdev/Math.sqrt(n);

		t[0] = (meanDiff-0)/(stdev/Math.sqrt(n));
		df = n-1;
		p[0] = ProbDist.TDist(Math.abs(t[0]), df);

		System.out.println("Running a paired t-test");

		testType = 1;
	}

	public Ttest(double[] data1, double[] data2) {
		this(Array.mean(data1), Array.stdev(data1), data1.length, Array.mean(data2), Array.stdev(data2), data2.length);

		Ftest = LevenesTest(new double[][] {data1, data2});
		Fprob = ProbDist.FDist(Ftest, 1, df);
		equal = Fprob>0.05;
	}

	public Ttest(int[] groupings, double[] values) {
		this(splitOut(groupings, values, Array.min(groupings)), splitOut(groupings, values, Array.min(groupings)+1));
	}

	public static double[] splitOut(int[] groupings, double[] values, int group) {
		DoubleVector dv = new DoubleVector();

		for (int i = 0; i<groupings.length; i++) {
			if (groupings[i]==group) {
				dv.add(values[i]);
			}
		}

		return Doubles.toArray(dv);
	}

	public Ttest(double x1Hat, double s1, int n1, double x2Hat, double s2, int n2) { // independent
																						// sample
																						// t-test
		this(x1Hat, s1, (double)n1, x2Hat, s2, (double)n2);
	}

	public Ttest(double x1Hat, double s1, double n1, double x2Hat, double s2, double n2) { // independent
																							// sample
																							// t-test
		double pooledVariance, se;

		t = new double[2];
		p = new double[2];
		ci05 = new double[2];
		ci025 = new double[2];

		meanDiff = x1Hat-x2Hat;

		// with equal variances
		pooledVariance = ((n1-1)*s1*s1+(n2-1)*s2*s2)/(n1+n2-2);
		stdev = Math.sqrt(pooledVariance);
		se = Math.sqrt(pooledVariance*(1/(double)n1+1/(double)n2));
		t[0] = meanDiff/se;
		df = (int)(n1+n2-2);
		p[0] = ProbDist.TDist(Math.abs(t[0]), df);

		ci05[0] = ProbDist.TDistReverse(0.1, n1+n2-2)*se; // distribution
		// assumes a
		// 2-tailed
		// test,
		// workaround
		ci025[0] = ProbDist.TDistReverse(0.05, n1+n2-2)*se;

		// with unequal variances
		t[1] = meanDiff/Math.sqrt((s1*s1/(double)n1+s2*s2/(double)n2));
		v = Math.pow(s1*s1/n1+s2*s2/n2, 2)/(Math.pow(s1*s1/n1, 2)/(n1-1)+Math.pow(s2*s2/n2, 2)/(n2-1));
		p[1] = ProbDist.TDist(t[1], Math.floor(v));

		ci05[1] = ProbDist.TDistReverse(0.1, n1+n2-2)*Math.sqrt((s1*s1/(double)n1+s2*s2/(double)n2)); // distribution
		// assumes
		// a
		// 2-tailed
		// test,
		// workaround
		ci025[1] = ProbDist.TDistReverse(0.05, n1+n2-2)*Math.sqrt((s1*s1/(double)n1+s2*s2/(double)n2));

		Ftest = -1;
		Fprob = -999;
		equal = false;

		testType = 2;
	}

	public String getReport() {
		String str = "";

		if (testType==1) {
			str += "Difference in means: "+ext.formDeci(meanDiff, 2, true)+"\n"+(testType<3?"Standard Deviation: "+ext.formDeci(stdev, 2, true)+"\n":"")+"t: "+ext.formDeci(t[0], 2, true)+"\n"+"df: "+df+"\n"+"one-tailed test: p="+ext.formDeci(p[0]/2, 4)+" ("+ext.formDeci(meanDiff-ci05[0], 2)+" OR "+ext.formDeci(meanDiff+ci05[0], 2)+")"+"\n"+"two-tailed test: p="+ext.formDeci(p[0], 4)+" ("+ext.formDeci(meanDiff-ci025[0], 2)+", "+ext.formDeci(meanDiff+ci025[0], 2)+")"+"\n";
		} else {
			str += "Difference in means: "+ext.formDeci(meanDiff, 2, true)+"\n"+"Variances\tT\tDF\t\tp-value\n"+"Equal\t"+ext.formDeci(t[0], 3, true)+"\t"+df+"\t\t"+ext.formDeci(p[0], 3)+"\n"+"Unqual\t"+ext.formDeci(t[1], 3, true)+"\t"+ext.formDeci(v, 3, true)+"\t"+ext.formDeci(p[1], 3)+"\n"+"\n"+"For H0: Variances are equal. "+ext.formDeci(Ftest, 3, true)+", p="+ext.formDeci(Fprob, 3, true);
		}

		return str;
	}

	public double getT() {
		if (testType==1||equal) {
			return t[0];
		} else {
			return t[1];
		}
	}

	public double getPvalue() {
		if (testType==1||equal) {
			return p[0];
		} else {
			return p[1];
		}
	}

	public double getStdev() {
		return stdev;
	}

	public double getDiff() {
		return meanDiff;
	}

	public boolean hasEqualVariances() {
		return equal;
	}

	public static double LevenesTest(double[][] data) {
		int k = data.length;
		int N = 0;
		int[] Ns = new int[k];
		double[] means = new double[k];
		double[][] Zijs = new double[k][];
		double[] ZiMeans = new double[k];
		double ZiMeansMean = 0;
		double[][] ZiMeansSq = new double[k][];
		double denomin = 0;
		double[] ZiMeansMeanSqN = new double[k];
		double numer = 0;

		for (int i = 0; i<k; i++) {
			Ns[i] = data[i].length;
			N += Ns[i];
			means[i] = Array.mean(data[i]);
		}

		for (int i = 0; i<k; i++) {
			Zijs[i] = new double[Ns[i]];
			for (int j = 0; j<Ns[i]; j++) {
				Zijs[i][j] = Math.abs(data[i][j]-means[i]);
			}
			ZiMeans[i] = Array.mean(Zijs[i]);
			ZiMeansMean += Array.sum(Zijs[i]);
		}
		ZiMeansMean /= N;

		for (int i = 0; i<k; i++) {
			ZiMeansSq[i] = new double[Ns[i]];
			for (int j = 0; j<Ns[i]; j++) {
				ZiMeansSq[i][j] = Math.pow(Zijs[i][j]-ZiMeans[i], 2);
			}
			denomin += Array.sum(ZiMeansSq[i]);
		}

		for (int i = 0; i<k; i++) {
			ZiMeansMeanSqN[i] = Ns[i]*Math.pow(ZiMeans[i]-ZiMeansMean, 2);
		}
		numer = (N-k)*Array.sum(ZiMeansMeanSqN);

		return numer/denomin;
	}

	public static void doPairedTtest(String filename) {
		BufferedReader reader;
		Vector<double[]> v = new Vector<double[]>();
		String[] line;
		double[][] data;

		try {
			reader = new BufferedReader(new FileReader(filename));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				v.add(new double[] {Double.parseDouble(line[0]), Double.parseDouble(line[1])});
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		data = new double[v.size()][];
		for (int i = 0; i<data.length; i++) {
			data[i] = v.elementAt(i);
		}

		Ttest test = new Ttest(data);
		System.out.println(test.getReport());
	}

	public static void doIndependentTtest(String filename) {
		BufferedReader reader;
		DoubleVector dv = new DoubleVector();
		String[] line = new String[] {" "};
		double[][] data = new double[2][];
		try {
			reader = new BufferedReader(new FileReader(filename));
			for (int i = 0; i<2; i++) {
				while (reader.ready()&&!line[0].equals("")) {
					line = reader.readLine().split("[\\s]+");
					if (!line[0].equals("")) {
						dv.add(Double.parseDouble(line[0]));
					}
				}
				line[0] = " ";
				data[i] = Doubles.toArray(dv);
				dv.removeAllElements();
				if (data[i].length==0) {
					System.err.println("Error in input file; separate samples by a blank line.");
					System.exit(1);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		Ttest test = new Ttest(data[0], data[1]);
		System.out.println(test.getReport());
	}

	public static void main(String[] args) throws IOException {
		// System.err.println("Ttest is usually called from within a program");
		// doPairedTtest("pairedData1.dat");
		doIndependentTtest("indepData2.dat");
		// System.out.println(new Ttest(new double[0], new
		// double[0]).getReport());
	}
}
