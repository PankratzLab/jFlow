package stats;

import java.math.BigInteger;

import common.Array;

public class Maths {
	public static final String[] OPERATORS = {"<=", "<", ">=", ">", "=", "!"};
	
	public static double limit(double d, double min, double max) {
		return d<min?min:(d>max?max:d);
	}

	private static BigInteger fact(int n) {
		BigInteger fact = BigInteger.ONE;
		for (int i = n; i>1; i--) {
			fact = fact.multiply(new BigInteger(i+""));
		}
		return fact;
	}

	public static int nCr(int n, int r) {
		return Integer.parseInt(fact(n).divide(fact(r).multiply(fact(n-r))).toString());
	}

	public static int[][] getIndicesForAllCombinations(int n, int r) {
		int[][] indices;
		int total, j, k;

		if (r>n) {
			System.err.println("Error - r cannot be less than n when generating combinations");
			System.exit(2);
		}
		if (n<1) {
			System.err.println("Error - n cannot be less than 1 when generating combinations");
			System.exit(3);
		}

		total = Maths.nCr(n, r);
		indices = new int[total][];
		indices[0] = Array.intArray(r);

		for (int i = 1; i<total; i++) {
			indices[i] = indices[i-1].clone();
			j = r-1;
			while (indices[i][j]==n-r+j) {
				j--;
			}
			indices[i][j] = indices[i][j]+1;
			for (k = j+1; k<r; k++) {
				indices[i][k] = indices[i][j]+k-j;
			}
		}

		return indices;
	}
	
	public static double log2(double num) {
		return (Math.log(num)/Math.log(2));
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
		return Math.exp(x)/(1+Math.exp(x));
	}
	
	public static double[] slopeAndIntercept(double x1, double y1, double x2, double y2) {
		double slope;
		double intercept;
		
		slope = (y2-y1)/(x2-x1);
		intercept = y1-(slope*x1);
		
		return new double[] {slope, intercept};		
	}
	
	public static String[] parseOp(String str) {
		for (int i = 0; i<OPERATORS.length; i++) {
			if (str.indexOf(OPERATORS[i]) >= 0) {
				return new String[] {str.substring(0, str.indexOf(OPERATORS[i])), OPERATORS[i], str.substring(str.indexOf(OPERATORS[i])+OPERATORS[i].length())};
			}
        }
		
		return null;
	}
	
	public static boolean op(double num1, double num2, String operator) {
		if (operator.equals("<=") && num1 <= num2) {
			return true;
		} else if (operator.equals("<") && num1 < num2) {
			return true;
		} else if (operator.equals(">=") && num1 >= num2) {
			return true;
		} else if (operator.equals(">") && num1 > num2) {
			return true;
		} else if (operator.equals("=") && num1 == num2) {
			return true;
		} else if (operator.equals("!") && num1 != num2) {
			return true;
		} 
		return false;
	}

	public static void main(String[] args) {
//		System.out.println(nCr(11, 5));
//		System.out.println(Array.toStr(slopeAndIntercept(432, 234, 132, 324)));
		System.out.println(Array.toStr(slopeAndIntercept(188, 900, 154, 1000)));

	}
}
