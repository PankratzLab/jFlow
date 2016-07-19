// master in 2.5 sec, Flat in 12.0 sec
package org.genvisis.stats;

import java.io.*;
import java.util.Date;

import org.genvisis.common.*;

public class FisherTemp {
	public static void yap() {
		int[][] matrix = {{0, 1, 2}, {0, 0, 0}, {2, 0, 1}};
//		int[][] matrix = {{1, 0, 2}, {0, 0, 0}, {0, 2, 1}};
//		int[][] matrix = {{14, 8, 20}, {13, 7, 5}, {40, 23, 9}};
//		int[][] matrix = {{14, 89, 120}, {13, 77, 58}, {40, 123, 91}};
		
		boolean equalSizedColumns;
		long time;
		int[][] work;
		int[] r, c;
//		int n;
		PrintWriter writer;
//		int numR;
//		int numC;
		int mRi;
		int mCi;
//		int trav;
		int c1;
		int c2;
		int poss;

		
//		n = 0;
		r = new int[matrix.length];
		c = new int[matrix[0].length];
		work = new int[r.length][c.length];
		
		for (int i = 0; i<r.length; i++) {
			for (int j = 0; j<c.length; j++) {
				work[i][j] = matrix[i][j];
            }
        }

		for (int i = 0; i<r.length; i++) {
			for (int j = 0; j<c.length; j++) {
				r[i] += work[i][j];
				c[j] += work[i][j];
//				n += work[i][j];
            }
        }
		
		equalSizedColumns = true;
		for (int i = 0; i<c.length-1; i++) {
			if (c[i] != c[i+1]) {
				equalSizedColumns = false;
			}
        }
		
		if(equalSizedColumns) {
			for (int i = 0; i<r.length; i++) {
				for (int j = 0; j<c.length; j++) {
					work[i][j] = matrix[j][i];
	            }
	        }

			r = new int[matrix.length];
			c = new int[matrix[0].length];
			for (int i = 0; i<r.length; i++) {
				for (int j = 0; j<c.length; j++) {
					r[i] += work[i][j];
					c[j] += work[i][j];
	            }
	        }
		}
		 
		c = Sort.putInOrder(c);
		
		try {
            time = new Date().getTime();
			writer = new PrintWriter(new FileWriter("master.out"));
			for (int i = 0; i<matrix.length; i++) {
				writer.println(Array.toStr(matrix[i]));
            }
			writer.println();
			writer.println("r:\t"+Array.toStr(r));
			writer.println("c:\t"+Array.toStr(c));
			writer.println();
			for(int o=0; o<=c[0]; o++) {
				for(int p=0; p<=c[0]-o;p++) {
					for(int q=0; q<=c[1];q++) {
						for(int s=0; s<=c[1]-q;s++) {
							work[0][0] = o;
							work[1][0] = p;
							work[2][0] = c[0]-work[0][0]-work[1][0];
							work[0][1] = q;
							work[1][1] = s;
							work[2][1] = c[1]-work[0][1]-work[1][1];
							work[0][2] = r[0]-work[0][0]-work[0][1];
							work[1][2] = r[1]-work[1][0]-work[1][1];
							work[2][2] = r[2]-work[2][0]-work[2][1];

//							tot = 0;
//							for (int i = 0; i<r.length; i++) {
//								for (int j = 0; j<c.length; j++) {
//									tot += Math.abs(work[i][j]);
//	                            }
//	                        }
//							if(tot==n) {
								for (int i = 0; i<r.length; i++) {
									writer.println(Array.toStr(work[i]));
	                            }
								writer.println();
//							}
						}
					}
				}
			}
			writer.close();
			System.out.println("master in "+ext.getTimeElapsed(time));
        } catch (Exception e) {
	        System.err.println("Error writing to master");
	        e.printStackTrace();
        }
		
		
		try {
            time = new Date().getTime();
			writer = new PrintWriter(new FileWriter("checks.out"));
			for (int i = 0; i<matrix.length; i++) {
				writer.println(Array.toStr(matrix[i]));
            }
			writer.println();
			writer.println("r:\t"+Array.toStr(r));
			writer.println("c:\t"+Array.toStr(c));
			writer.println();
//			numR = r.length;
//			numC = c.length;
			mRi = r.length-1;
			mCi = c.length-1;
			
//			trav = 0;
			c1 = (int)((c[0]+2)*((c[0]+1)/2.0));
			c2 = (int)((c[1]+2)*((c[1]+1)/2.0));
			poss = c1*c2;
			for (int rep = 0; rep<poss; rep++) {
				work = FisherMod.mod(rep, r, c);

				for (int j = 0; j<mCi; j++) {
					work[mRi][j] = c[j];
					for (int i = 0; i<mRi; i++) {
						work[mRi][j] -= work[i][j];
	                }
	            }
				for (int i = 0; i<=mRi; i++) {
					work[i][mCi] = r[i];
					for (int j = 0; j<mRi; j++) {
						work[i][mCi] -= work[i][j];
	                }
	            }

				for (int i = 0; i<r.length; i++) {
					writer.println(Array.toStr(work[i]));
	            }
				writer.println();
			}					
			
			
			writer.close();
			System.out.println("checks in "+ext.getTimeElapsed(time));
        } catch (Exception e) {
	        System.err.println("Error writing to checks");
	        e.printStackTrace();
        }
        
        try {
        	BufferedReader reader = new BufferedReader(new FileReader("checks.out"));
        	BufferedReader reader2 = new BufferedReader(new FileReader("master.out"));
	        while (reader.ready()) {
	        	if (!reader.readLine().equals(reader2.readLine())) {
	        		System.err.println("FAILED!!!");
	        		reader.close();
	        		reader2.close();
	        		return;
	        	}
	        }
	        reader.close();
    		reader2.close();
	        if (reader2.ready()) {
	        	System.out.println("incomplete");
	        } else {
	        	System.out.println("passed");
	        }
	        System.out.println();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+"\"");
	        System.exit(2);
        }
    }
	

	public static void main(String[] args) {
		for (int i = 0; i<1; i++) {
			yap();
        }
	}
}
