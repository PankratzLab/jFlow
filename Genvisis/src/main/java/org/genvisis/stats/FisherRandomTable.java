// check if any of the marginals is zero, collapse if need be
package org.genvisis.stats;

import org.genvisis.common.Array;

public class FisherRandomTable {
	public static int[][] rcont(int[] r, int[] c) {

//			****************************************************************************80
//			  Purpose:
//			    RCONT generates a random two-way table with given marginal totals.
//			  Discussion:
//			    Each time the program is called, another table will be randomly
//			    generated.
//			    Note that it should be the case that the sum of the row totals
//			    is equal to the sum of the column totals.  However, this program
//			    does not check for that condition.
//			  Author:
//			    Original FORTRAN77 version by James Boyett.
//			    C++ version by John Burkardt.
//			  Reference:
//			    James Boyett,
//			    Algorithm AS 144:
//			    Random R x C Tables with Given Row and Column Totals,
//			    Applied Statistics,
//			    Volume 28, Number 3, pages 329-332, 1979.
//			  Parameters:
//			    Input, int NROW, the number of rows in the observed matrix.
//			    Input, int NCOL, the number of columns in the observed matrix.
//			    Input, int NROWT[NROW], the row totals of the observed matrix.
//			    Input, int NCOLT[NCOL], the column totals of the observed matrix.
//			    Input/output, int NSUBT[NCOL], used by RCONT for partial column sums.
//			    Must not be changed by the calling program.
//			    Output, int MATRIX[NROW*NCOL], the randomly generated matrix.
//			    Input/output, bool *KEY, should be set to FALSE by the user before
//			    the initial call.  RCONT will reset it to TRUE, and it should be left
//			    at that value for subsequent calls in which the same values of NROW,
//			    NCOL, NROWT and NCOLT are being used.
//			    Output, int *IFAULT, fault indicator.
//			    0, no error occured.
//			    1, NROW <= 0.
//			    2, NCOL <= 1.
//			    3, some entry of NROWT is less than 0.
//			    4, some entry of NCOLT is less than 0.

		int i;
		int ii;
		int j;
		int k;
		int limit;
		int[] nnvect;
		int noct;
		int ntemp;
		int ntotal;
		int[] nvect = null;
		int[][] matrix;
		int[] nsubt;
		int ncol, nrow;
		
		nrow = r.length;
		ncol = c.length;
		
		nsubt = new int[c.length];

		nsubt[0] = c[0];
		for ( j = 1; j < ncol; j++ ) {
			nsubt[j] = nsubt[j-1] + c[j];
		}
		ntotal = nsubt[ncol-1];

		nvect = new int[ntotal];
		for ( i = 0; i < ntotal; i++ ) {
			nvect[i] = i + 1;
		}

//			  Initialize vector to be permuted.
		nnvect = new int[ntotal];
		for ( i = 0; i < ntotal; i++ ) {
			nnvect[i] = nvect[i];
		}

//			  Permute vector.
		ntemp = ntotal;
		for ( i = 0; i < ntotal; i++ ) {
			noct = ( int ) ( Math.random() * ( double ) ( ntemp ) + 1.0 );
			nvect[i] = nnvect[noct-1];
			nnvect[noct-1] = nnvect[ntemp-1];
			ntemp--;
		}

		matrix = new int[nrow][ncol];

		ii = 0;
		for ( i = 0; i < nrow; i++ ) {
			limit = r[i];

			for ( k = 0; k < limit; k++ ) {
				for ( j = 0; j < ncol; j++ ) {
			        if ( nvect[ii] <= nsubt[j] ) {
			          ii = ii + 1;
			          matrix[i][j]++;
			          break;
			        }
				}
			}
		}

//		delete [] nnvect;
		return matrix;
	}
	
	public static int[][] rcont2 (int[] nrowt, int[] ncolt, int ntotal, double[] fact) {

//			****************************************************************************80
//			  Purpose:
//			    RCONT2 constructs a random two-way contingency table with given sums.
//			  Discussion:
//			    It is possible to specify row and column sum vectors which
//			    correspond to no table at all.  As far as I can see, this routine does
//			    not detect such a case.
//			  Author:
//			    Original FORTRAN77 version by WM Patefield.
//			    C++ version by John Burkardt.
//			  Reference:
//			    WM Patefield,
//			    Algorithm AS 159:
//			    An Efficient Method of Generating RXC Tables with
//			    Given Row and Column Totals,
//			    Applied Statistics,
//			    Volume 30, Number 1, 1981, pages 91-97.
//			  Parameters:
			//
//			    Input, int NROW, NCOL, the number of rows and columns 
//			    in the table.  NROW and NCOL must each be at least 2.
			//
//			    Input, int NROWT[NROW], NCOLT[NCOL], the row and column 
//			    sums.  Each entry must be positive.
			//
//			    Input/output, bool *KEY, a flag that indicates whether data has
//			    been initialized for this problem.  Set KEY = .FALSE. before the first
//			    call.
			//
//			    Input/output, int *SEED, a seed for the random number generator.
			//
//			    Output, int MATRIX[NROW*NCOL], the matrix.
			//
//			    Output, int *IERROR, an error flag, which is returned 
//			    as 0 if no error occurred.
			//
		int nrow, ncol;
		int matrix[];
		boolean done1;
		boolean done2;
		int i;
		int ia;
		int iap;
		int ib;
		int ic;
		int id;
		int idp;
		int ie;
		int igp;
		int ihp;
		int ii;
		int iip;
		int j;
		int jc;
		int[] jwork;
		int l;
		boolean lsm;
		boolean lsp;
		int m;
		int nll;
		int nlm;
		int nlmp;
		int nrowtl;
		double r;
		double sumprb;
		double x;
		double y;

		nrow = nrowt.length;
		ncol = ncolt.length;
		matrix = new int[nrow*ncol];
		done2 = false;
		ib = -1;

		jwork = new int[ncol];
		for ( i = 0; i < ncol - 1; i++ ) {
			jwork[i] = ncolt[i];
		}

		jc = ntotal;

		for (l = 0; l < nrow - 1; l++ ) {
			nrowtl = nrowt[l];
			ia = nrowtl;
			ic = jc;
			jc = jc - nrowtl;

			for ( m = 0; m < ncol - 1; m++ ) {
				id = jwork[m];
				ie = ic;
				ic = ic - id;
				ib = ie - ia;
				ii = ib - id;

//		Test for zero entries in matrix.
				if ( ie == 0 ) {
					ia = 0;
					for ( j = m; j < ncol; j++ ) {
						matrix[l+j*nrow] = 0;
					}
					break;
				}

//		Generate a pseudo-random number.
				r = Math.random();

//		Compute the conditional expected value of MATRIX(L,M).

				done1 = false;
				for ( ; ; ) {
					nlm = ( int ) ( ( double ) ( ia * id ) / ( double ) ( ie ) + 0.5 );
					iap = ia + 1;
					idp = id + 1;
					igp = idp - nlm;
					ihp = iap - nlm;
					nlmp = nlm + 1;
					iip = ii + nlmp;
					x = Math.exp( fact[iap-1] + fact[ib] + fact[ic] + fact[idp-1] - 
							fact[ie] - fact[nlmp-1] - fact[igp-1] - fact[ihp-1] - fact[iip-1] );

					if ( r <= x ) {
						break;
					}

					sumprb = x;
					y = x;
					nll = nlm;
					lsp = false;
					lsm = false;

//		 Increment entry in row L, column M.
					while ( !lsp ) {
						j = ( id - nlm ) * ( ia - nlm );

						if ( j == 0 ) {
							lsp = true;
						} else {
							nlm = nlm + 1;
							x = x * ( double ) ( j ) / ( double ) ( nlm * ( ii + nlm ) );
							sumprb = sumprb + x;

							if ( r <= sumprb ) {
								done1 = true;
								break;
							}
						}

						done2 = false;

						while ( !lsm ) {
//		 Decrement the entry in row L, column M.
							j = nll * ( ii + nll );

							if ( j == 0 ) {
								lsm = true;
								break;
							}

							nll = nll - 1;
							y = y * ( double ) ( j ) / ( double ) ( ( id - nll ) * ( ia - nll ) );
							sumprb = sumprb + y;

							if ( r <= sumprb ) {
								nlm = nll;
								done2 = true;
								break;
							}

							if ( !lsp ) {
								break;
							}
						}

						if ( done2 ) {
							break;
						}
					}

					if ( done1 ) {
						break;
					}

					if ( done2 ) {
						break;
					}

					r = Math.random();
					r = sumprb * r;
				}

				matrix[l+m*nrow] = nlm;
				ia = ia - nlm;
				jwork[m] = jwork[m] - nlm;
			}
			matrix[l+(ncol-1)*nrow] = ia;
		}

//		 Compute the last row.
		for ( j = 0; j < ncol - 1; j++ ) {
			matrix[nrow-1+j*nrow] = jwork[j];
		}
		matrix[nrow-1+(ncol-1)*nrow] = ib - matrix[nrow-1+(ncol-2)*nrow];

		int[][] actualMatrix;
		actualMatrix = new int[nrow][ncol];
		int count = 0;
		for (i = 0; i<nrow; i++) {
			for (j = 0; j<ncol; j++) {
				actualMatrix[i][j] = matrix[count];
				count++;
            }
        }
		 
		return actualMatrix;
	}

	public static void main(String[] args) {
		int[][] matrix;
		
		for (int rep = 0; rep<100; rep++) {
			matrix = rcont(new int[] {3,0,3}, new int[] {1,2,3});
			for (int i = 0; i<matrix.length; i++) {
				System.out.println(Array.toStr(matrix[i]));
	        }
        }
	}
}
