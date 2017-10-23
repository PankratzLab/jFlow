package org.genvisis.stats;

import org.genvisis.common.ext;

public class FishersExact2by2Calculator {
	private double[] logFactorial;
	private int a;
	private int b;
	private int c;
	private int d;
	private int n;
	private int tmp;
	
	public FishersExact2by2Calculator() {
		logFactorial = null;

		this.a = 0;
    	this.b = 0;
    	this.c = 0;
    	this.d = 0;
	}
	
    public void loadMatrix(int[][] matrix) {
    	loadCounts(matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]);
    }

    public void loadCounts(int a, int b, int c, int d) {
    	this.a = a;
    	this.b = b;
    	this.c = c;
    	this.d = d;
    }
    
    private double testComputation() {
        return Math.exp(
        	(logFactorial[a + b] +
        	logFactorial[c + d] +
        	logFactorial[a + c] +
        	logFactorial[b + d])
        	-
        	(logFactorial[a + b + c + d] +
        	logFactorial[a] +
        	logFactorial[b] +
        	logFactorial[c] +
        	logFactorial[d])
        );
    }
    
    public double getPvalue() {
    	return getPvalue(false);
    }
    
    public double getPvalue(boolean oneTailed) {
        n=a+b+c+d;

        if (logFactorial == null || logFactorial.length < n+1) {
            logFactorial = new double[n+1];
            logFactorial[0] = 0.0;
            for (int i = 1; i <= n; i++) {
                logFactorial[i] = logFactorial[i-1] + Math.log(i);
            }
        }
        
        if (a * d > b * c) {
        	tmp = a;
        	a = b;
        	b = tmp;

        	tmp = c;
        	c = d;
        	d = tmp;
        }
        if (a > d) {
        	tmp = a;
        	a = d;
        	d = tmp;
        }
        if (b > c) {
        	tmp = b;
        	b = c;
        	c = tmp;
        }

        int original_a = a;
        double sumP = 0;

        double p = testComputation();
        double p_1 = p;

        while (a >= 0) {
            sumP += p;
            if (a == 0) {
            	break;
            }
            a--;
            b++;
            c++;
            d--;
            p = testComputation();
        }
        if (oneTailed) {
        	return Math.min(sumP, 1);
        }

        a = b;
        b = 0;
        c = c - a;
        d = d + a;
        p = testComputation();

        while (p < p_1) {
            if (a == original_a) {
            	break;
            }
            sumP += p;
            a--;
            b++;
            c++;
            d--;
            p = testComputation();
        }
    	return Math.min(sumP, 1);
    }

    public static void main(String[] args) {
		double chiP, fishP;
		int[][][] testMatrices = new int[][][] {
				{{0, 20}, {0, 25}},
				{{2, 20}, {3, 25}},
				{{0, 20}, {5, 20}},
				{{0, 2100}, {0, 2500}},
				{{0, 2100}, {5, 2500}},
				{{0, 2100}, {10, 2500}},
				{{210, 2100}, {250, 2500}},
		};
		
		FishersExact2by2Calculator calc = new FishersExact2by2Calculator();
		
		System.out.println("chi^2 p\tFisher's Exact p");
		for (int i = 0; i < testMatrices.length; i++) {
			chiP = ProbDist.ChiDist(ContingencyTable.ChiSquare(testMatrices[i], false), 1);
			calc.loadMatrix(testMatrices[i]);
			fishP = calc.getPvalue();
			System.out.println(ext.prettyP(chiP)+"\t"+fishP);
		}			
	}
}
