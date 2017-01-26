package org.genvisis.stats;

import java.math.BigInteger;
import java.util.List;
import java.util.Map;

import org.genvisis.common.ArrayUtils;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;

public class Maths {
	
	public static final List<String> OPERATORS;

	static {
		ImmutableList.Builder<String> builder = ImmutableList.builder();
		for (OPERATOR operator : OPERATOR.values()) {
			builder.add(operator.getSymbol());
		}
		OPERATORS = builder.build();
	}

	private Maths() {
		// Prevent construction
	}

	
	public enum OPERATOR {
		// Valid operators that will be searched for in the following order:
		LESS_THAN_OR_EQUAL("<="),
		LESS_THAN("<"),
		GREATER_THAN_OR_EQUAL(">="),
		GREATER_THAN(">"),
		EQUAL("="),
		NOT("!");
		
		private static final Map<String, OPERATOR> SYMBOL_MAP;
		static {
			ImmutableMap.Builder<String, OPERATOR> builder = ImmutableMap.builder();
			for (OPERATOR operator : OPERATOR.values()) {
				builder.put(operator.getSymbol(), operator);
			}
			SYMBOL_MAP = builder.build();
		}
		
		private final String symbol;
		
		OPERATOR(String symbol) {
			this.symbol = symbol;
		}

		public String getSymbol() {
			return symbol;
		}
		
		public static OPERATOR forSymbol(String symbol) {
			return SYMBOL_MAP.get(symbol);
		}
		
		public boolean check(double num1, double num2) {
			switch (this) {
				case LESS_THAN_OR_EQUAL:
					return num1 <= num2;
				case LESS_THAN:
					return num1 < num2;
				case GREATER_THAN_OR_EQUAL:
					return num1 >= num2;
				case GREATER_THAN:
					return num1 > num2;
				case EQUAL:
					return num1 == num2;
				case NOT:
					return num1 != num2;
				default:
					System.err.println("Cannot perform check, invalid " + getClass().getName() + ":" + this.name());
					return false;
			}
		}

	}
	
	public static double limit(double d, double min, double max) {
		return d < min ? min : (d > max ? max : d);
	}

	private static BigInteger fact(int n) {
		BigInteger fact = BigInteger.ONE;
		for (int i = n; i > 1; i--) {
			fact = fact.multiply(BigInteger.valueOf(n));
		}
		return fact;
	}

	public static int nCr(int n, int r) {
		return Integer.parseInt(fact(n).divide(fact(r).multiply(fact(n - r))).toString());
	}

	public static double roundDouble(double d, int numPlaces) {
		double rounder = Math.pow(10, numPlaces);
		return Math.round(d * rounder) / rounder;
	}

	public static int[][] getIndicesForAllCombinations(int n, int r) {
		int[][] indices;
		int total, j, k;

		if (r > n) {
			System.err.println("Error - r cannot be less than n when generating combinations");
			System.exit(2);
		}
		if (n < 1) {
			System.err.println("Error - n cannot be less than 1 when generating combinations");
			System.exit(3);
		}

		total = Maths.nCr(n, r);
		indices = new int[total][];
		indices[0] = ArrayUtils.arrayOfIndices(r);

		for (int i = 1; i < total; i++) {
			indices[i] = indices[i - 1].clone();
			j = r - 1;
			while (indices[i][j] == n - r + j) {
				j--;
			}
			indices[i][j] = indices[i][j] + 1;
			for (k = j + 1; k < r; k++) {
				indices[i][k] = indices[i][j] + k - j;
			}
		}

		return indices;
	}

	public static double log2(double num) {
		return Math.log(num) / Math.log(2);
	}

	public static double min(double num1, double num2) {
		if (Double.isNaN(num1) && Double.isNaN(num2)) {
			return Double.NaN;
		} else if (Double.isNaN(num1)) {
			return num2;
		} else if (Double.isNaN(num2)) {
			return num1;
		} else {
			return Math.min(num1, num2);
		}
	}

	public static double max(double num1, double num2) {
		if (Double.isNaN(num1) && Double.isNaN(num2)) {
			return Double.NaN;
		} else if (Double.isNaN(num1)) {
			return num2;
		} else if (Double.isNaN(num2)) {
			return num1;
		} else {
			return Math.max(num1, num2);
		}
	}

	public static float min(float num1, float num2) {
		if (Float.isNaN(num1) && Float.isNaN(num2)) {
			return Float.NaN;
		} else if (Float.isNaN(num1)) {
			return num2;
		} else if (Float.isNaN(num2)) {
			return num1;
		} else {
			return Math.min(num1, num2);
		}
	}

	public static float max(float num1, float num2) {
		if (Float.isNaN(num1) && Float.isNaN(num2)) {
			return Float.NaN;
		} else if (Float.isNaN(num1)) {
			return num2;
		} else if (Float.isNaN(num2)) {
			return num1;
		} else {
			return Math.max(num1, num2);
		}
	}

	public static double logit(double x) {
		return Math.exp(x) / (1 + Math.exp(x));
	}

	public static double[] slopeAndIntercept(double x1, double y1, double x2, double y2) {
		double slope;
		double intercept;

		slope = (y2 - y1) / (x2 - x1);
		intercept = y1 - (slope * x1);

		return new double[] {slope, intercept};
	}

	public static String[] parseOp(String str) {
		for (String element : OPERATORS) {
			if (str.indexOf(element) >= 0) {
				return new String[] {	str.substring(0, str.indexOf(element)), element,
															str.substring(str.indexOf(element) + element.length())};
			}
		}

		return null;
	}

	public static boolean op(double num1, double num2, String operatorSymbol) {
		OPERATOR operator = OPERATOR.forSymbol(operatorSymbol);
		if (operator == null) {
			System.err.println("Cannot perform operation, " + OPERATOR.class.getName() + " for symbol " + operatorSymbol + " does not exist");
			return false;
		}
		return operator.check(num1, num2);
	}
	public static void main(String[] args) {
		// System.out.println(nCr(11, 5));
		// System.out.println(Array.toStr(slopeAndIntercept(432, 234, 132, 324)));
		System.out.println(ArrayUtils.toStr(slopeAndIntercept(188, 900, 154, 1000)));

	}

	public static double[][] cartesianToPolar(double[] x, double[] y) {
		double[][] polar = new double[x.length][2];
		for (int i = 0; i < x.length; i++) {
			polar[i][0] = Math.atan2(y[i], x[i]);
			polar[i][1] = Math.hypot(x[i], y[i]);
		}
		return polar;
	}

	public static double[][] polarToCartesian(double[][] polar) {
		double[] x = new double[polar.length];
		double[] y = new double[polar.length];

		for (int i = 0; i < polar.length; i++) {
			double theta = polar[i][0];
			double r = polar[i][1];
			x[i] = r * Math.cos(theta);
			y[i] = r * Math.sin(theta);
		}

		return new double[][] {x, y};
	}
}
