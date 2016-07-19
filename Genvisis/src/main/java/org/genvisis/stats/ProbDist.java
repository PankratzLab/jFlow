package org.genvisis.stats;

public class ProbDist {
	public static double TDist(double t, double df) {
		double pi = Math.PI;
		double pj2 = pi/2;

		t = Math.abs(t);
		double rt = t/Math.sqrt(df);
		double fk = Math.atan(rt);
		if (df==1) {
			return 1-fk/pj2;
		}
		double ek = Math.sin(fk);
		double dk = Math.cos(fk);
		if ((df%2)==1) {
			return 1-(fk+ek*dk*StatCom(dk*dk, 2, df-3, -1))/pj2;
		} else {
			return 1-ek*StatCom(dk*dk, 1, df-3, -1);
		}
	}

	public static double TDistReverse(double p, double df) {
		double v = 0.5, dv = 0.5, t = 0;

		while (dv>1e-9) {
			t = 1/v-1;
			dv = dv/2;
			if (TDist(t, df)>p) {
				v = v-dv;
			} else {
				v = v+dv;
			}
		}

		return t;
	}

	public static double StatCom(double q, int i, double j, int b) {
		double zz = 1;
		double z = zz;
		int k = i;

		while (k<=j) {
			zz = zz*q*k/(k-b);
			z = z+zz;
			k = k+2;
		}

		return z;
	}

	public static double FDist(double f, int df1, int df2) {
		double x = (double)df2/((double)df1*f+(double)df2);
		if ((df1%2)==0) {
			return StatCom(1-x, df2, df1+df2-4, df2-2)*Math.pow(x, (double)df2/2);
		}
		if ((df2%2)==0) {
			return 1-StatCom(x, df1, df1+df2-4, df1-2)*Math.pow(1-x, (double)df1/2);
		}
		double th = Math.atan(Math.sqrt(df1*f/df2));
		double a = th/(Math.PI/2);
		double sth = Math.sin(th);
		double cth = Math.cos(th);
		if (df2>1) {
			a = a+sth*cth*StatCom(cth*cth, 2, df2-3, -1)/(Math.PI/2);
		}
		if (df1==1) {
			return 1-a;
		}
		double c = 4*StatCom(sth*sth, df2+1, df1+df2-4, df2-2)*sth*Math.pow(cth, df2)/Math.PI;
		if (df2==1) {
			return 1-a+c/2;
		}
		int k = 2;
		while (k<=(df2-1)/2) {
			c = c*k/(k-0.5);
			k = k+1;
		}
		return 1-a+c;
	}
	
	public static double ChiDist(double x, int n) {
		if (Double.isNaN(x)) {
			return Double.NaN;
		}
		
		if (x>1000|n>1000) { // bitwise OR - prevent short-circuiting and can be faster
			double q = NormDist((Math.pow(x/n, 1.0/3.0)+2.0/(9.0*n)-1)/Math.sqrt(2.0/(9.0*n)))/2.0;
			if (x>n) {
				return q;
			} else {
				return 1-q;
			}
		}
		double p = Math.exp(-0.5*x);
		if ((n%2)==1) {
			p = p*Math.sqrt(2*x/Math.PI);
		}
		int k = n;
		while (k>=2) {
			p = p*x/k;
			k = k-2;
		}
		double t = p;
		int a = n;
		while (t>1e-15*p) {
			a = a+2;
			t = t*x/a;
			p = p+t;
		}
		return 1-p;
	}

	// 2-sided test
	public static double NormDist(double z) {
		double q = z*z;
		if (Math.abs(z)>7) {
			return (1-1/q+3/(q*q))*Math.exp(-q/2)/(Math.abs(z)*Math.sqrt(Math.PI/2));
		} else {
			return ChiDist(q, 1);
		}
	}

	// 2-sided test
	public static double NormDistReverse(double p) {
		double v = 0.5, dv = 0.5, z = 0;
		while (dv>1e-6) {
			z = 1/v-1;
			dv = dv/2;
			if (NormDist(z)>p) {
				v = v-dv;
			} else {
				v = v+dv;
			}
		}
		return z;
	}

	public static double ChiDistReverse(double p, int n) {
		double v = 0.5, dv = 0.5, x = 0;
		while (dv>1e-10) {
			x = 1/v-1;
			dv = dv/2;
			if (ChiDist(x, n)>p) {
				v = v-dv;
			} else {
				v = v+dv;
			}
		}
		return x;
	}

	public static double FDistReverse(double p, int n1, int n2) {
		double v = 0.5, dv = 0.5, f = 0;
		while (dv>1e-10) {
			f = 1/v-1;
			dv = dv/2;
			if (FDist(f, n1, n2)>p) {
				v = v-dv;
			} else {
				v = v+dv;
			}
		}
		return f;
	}
}
