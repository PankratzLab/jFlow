package org.genvisis.stats;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Date;

import org.genvisis.common.Array;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class FishersExact {
	public static double calc(int[][] matrix) {
		return calc(matrix, 0, true);
	}

	public static double calc(int[][] matrix, int permutations, boolean verbose) {
		double num, denom, pex, totprob;
		boolean equalSizedColumns;
		double[] facts, lookup;
		int[][] work;
		int[] r, c;
		int n;

		n = 0;
		r = new int[matrix.length];
		c = new int[matrix[0].length];
		work = new int[r.length][c.length];

		for (int i = 0; i < r.length - 1; i++) {
			if (matrix[i].length != matrix[i + 1].length) {
				System.err.println("Error - different number of columns in row " + (i + 2));
				return Double.NEGATIVE_INFINITY;
			}
		}

		for (int i = 0; i < r.length; i++) {
			for (int j = 0; j < c.length; j++) {
				work[i][j] = matrix[i][j];
			}
		}

		for (int i = 0; i < r.length; i++) {
			for (int j = 0; j < c.length; j++) {
				r[i] += work[i][j];
				c[j] += work[i][j];
				n += work[i][j];
			}
		}

		equalSizedColumns = true;
		for (int i = 0; i < c.length - 1; i++) {
			if (c[i] != c[i + 1]) {
				equalSizedColumns = false;
			}
		}

		if (equalSizedColumns) {
			work = new int[c.length][r.length];
			for (int i = 0; i < r.length; i++) {
				for (int j = 0; j < c.length; j++) {
					work[j][i] = matrix[i][j];
				}
			}

			r = new int[work.length];
			c = new int[work[0].length];
			for (int i = 0; i < r.length; i++) {
				for (int j = 0; j < c.length; j++) {
					r[i] += work[i][j];
					c[j] += work[i][j];
				}
			}
		}

		lookup = new double[n + 1];
		for (int i = 1; i < lookup.length; i++) {
			lookup[i] = lookup[i - 1] + Math.log(i);
		}

		c = Sort.putInOrder(c);

		facts = new double[r.length + c.length];
		for (int i = 0; i < r.length; i++) {
			facts[i] = lookup[r[i]];
		}
		for (int j = 0; j < c.length; j++) {
			facts[r.length + j] = lookup[c[j]];
		}
		num = Array.sum(Sort.putInOrder(facts));

		facts = new double[r.length * c.length + 1];
		for (int i = 0; i < r.length; i++) {
			for (int j = 0; j < c.length; j++) {
				facts[c.length * i + j + 1] = lookup[work[i][j]];
			}
		}
		facts[0] = lookup[n];
		denom = Array.sum(Sort.putInOrder(facts));
		pex = Math.exp(num - denom);

		totprob = 0;
		switch (r.length * 10 + c.length) {
			case 22:
				for (int a1 = 0; a1 <= c[0]; a1++) {
					work[0][0] = a1;
					work[1][0] = c[0] - work[0][0];
					work[0][1] = r[0] - work[0][0];
					work[1][1] = r[1] - work[1][0];

					if (checkN(work, n)) {
						totprob += workUp(work, facts, num, pex, lookup);
					}
				}
				break;
			case 23:
				for (int a1 = 0; a1 <= c[0]; a1++) {
					work[0][0] = a1;
					work[1][0] = c[0] - work[0][0];

					for (int b1 = 0; b1 <= c[1]; b1++) {
						work[0][1] = b1;
						work[1][1] = c[1] - work[0][1];

						work[0][2] = r[0] - work[0][0] - work[0][1];
						work[1][2] = r[1] - work[1][0] - work[1][1];

						if (checkN(work, n)) {
							totprob += workUp(work, facts, num, pex, lookup);
						}
					}
				}
				break;
			case 32:
				for (int a1 = 0; a1 <= c[0]; a1++) {
					work[0][0] = a1;
					for (int a2 = 0; a2 <= c[0] - a1; a2++) {
						work[1][0] = a2;
						work[2][0] = c[0] - work[0][0] - work[1][0];

						work[0][1] = r[0] - work[0][0];
						work[1][1] = r[1] - work[1][0];
						work[2][1] = r[2] - work[2][0];

						if (checkN(work, n)) {
							totprob += workUp(work, facts, num, pex, lookup);
						}
					}
				}
				break;
			case 33:
				for (int a1 = 0; a1 <= c[0]; a1++) {
					work[0][0] = a1;
					for (int a2 = 0; a2 <= c[0] - a1; a2++) {
						work[1][0] = a2;
						work[2][0] = c[0] - work[0][0] - work[1][0];
						for (int b1 = 0; b1 <= c[1]; b1++) {
							work[0][1] = b1;

							for (int b2 = 0; b2 <= c[1] - b1; b2++) {
								work[1][1] = b2;
								work[2][1] = c[1] - work[0][1] - work[1][1];

								work[0][2] = r[0] - work[0][0] - work[0][1];
								work[1][2] = r[1] - work[1][0] - work[1][1];
								work[2][2] = r[2] - work[2][0] - work[2][1];

								if (checkN(work, n)) {
									totprob += workUp(work, facts, num, pex, lookup);
								}
							}
						}
					}
				}
				break;
			case 43:
				for (int a1 = 0; a1 <= c[0]; a1++) {
					work[0][0] = a1;
					for (int a2 = 0; a2 <= c[0] - a1; a2++) {
						work[1][0] = a2;
						for (int a3 = 0; a3 <= c[0] - a1 - a2; a3++) {
							work[2][0] = a3;
							work[3][0] = c[0] - work[0][0] - work[1][0] - work[2][0];
							for (int b1 = 0; b1 <= c[1]; b1++) {
								work[0][1] = b1;
								for (int b2 = 0; b2 <= c[1] - b1; b2++) {
									work[1][1] = b2;
									for (int b3 = 0; b3 <= c[1] - b1 - b2; b3++) {
										work[2][1] = b3;
										work[3][1] = c[1] - work[0][1] - work[1][1] - work[2][1];

										work[0][2] = r[0] - work[0][0] - work[0][1];
										work[1][2] = r[1] - work[1][0] - work[1][1];
										work[2][2] = r[2] - work[2][0] - work[2][1];
										work[3][2] = r[3] - work[3][0] - work[3][1];

										if (checkN(work, n)) {
											totprob += workUp(work, facts, num, pex, lookup);
										}
									}
								}
							}
						}
					}
				}
				break;
			case 53:
				for (int a1 = 0; a1 <= c[0]; a1++) {
					work[0][0] = a1;
					for (int a2 = 0; a2 <= c[0] - a1; a2++) {
						work[1][0] = a2;
						for (int a3 = 0; a3 <= c[0] - a1 - a2; a3++) {
							work[2][0] = a3;
							for (int a4 = 0; a4 <= c[0] - a1 - a2 - a3; a4++) {
								work[3][0] = a4;
								work[4][0] = c[0] - work[0][0] - work[1][0] - work[2][0] - work[3][0];
								for (int b1 = 0; b1 <= c[1]; b1++) {
									work[0][1] = b1;
									for (int b2 = 0; b2 <= c[1] - b1; b2++) {
										work[1][1] = b2;
										for (int b3 = 0; b3 <= c[1] - b1 - b2; b3++) {
											work[2][1] = b3;
											for (int b4 = 0; b4 <= c[1] - b1 - b2 - b3; b4++) {
												work[3][1] = b4;
												work[4][1] = c[1] - work[0][1] - work[1][1] - work[2][1] - work[3][1];

												work[0][2] = r[0] - work[0][0] - work[0][1];
												work[1][2] = r[1] - work[1][0] - work[1][1];
												work[2][2] = r[2] - work[2][0] - work[2][1];
												work[3][2] = r[3] - work[3][0] - work[3][1];
												work[4][2] = r[4] - work[4][0] - work[4][1];

												if (checkN(work, n)) {
													totprob += workUp(work, facts, num, pex, lookup);
												}
											}
										}
									}
								}
							}
						}
					}
				}
				break;
			case 63:
				for (int a1 = 0; a1 <= c[0]; a1++) {
					work[0][0] = a1;
					for (int a2 = 0; a2 <= c[0] - a1; a2++) {
						work[1][0] = a2;
						for (int a3 = 0; a3 <= c[0] - a1 - a2; a3++) {
							work[2][0] = a3;
							for (int a4 = 0; a4 <= c[0] - a1 - a2 - a3; a4++) {
								work[3][0] = a4;
								for (int a5 = 0; a5 <= c[0] - a1 - a2 - a3 - a4; a5++) {
									work[4][0] = a5;
									work[5][0] = c[0]	- work[0][0] - work[1][0] - work[2][0] - work[3][0]
																- work[4][0];
									for (int b1 = 0; b1 <= c[1]; b1++) {
										work[0][1] = b1;
										for (int b2 = 0; b2 <= c[1] - b1; b2++) {
											work[1][1] = b2;
											for (int b3 = 0; b3 <= c[1] - b1 - b2; b3++) {
												work[2][1] = b3;
												for (int b4 = 0; b4 <= c[1] - b1 - b2 - b3; b4++) {
													work[3][1] = b4;
													for (int b5 = 0; b5 <= c[1] - b1 - b2 - b3 - b4; b5++) {
														work[4][1] = b5;
														work[5][1] = c[1]	- work[0][1] - work[1][1] - work[2][1] - work[3][1]
																					- work[4][1];

														work[0][2] = r[0] - work[0][0] - work[0][1];
														work[1][2] = r[1] - work[1][0] - work[1][1];
														work[2][2] = r[2] - work[2][0] - work[2][1];
														work[3][2] = r[3] - work[3][0] - work[3][1];
														work[4][2] = r[4] - work[4][0] - work[4][1];
														work[5][2] = r[5] - work[5][0] - work[5][1];

														if (checkN(work, n)) {
															totprob += workUp(work, facts, num, pex, lookup);
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
				break;
		}

		if (totprob >= 1) {
			totprob = 1.0;
		}

		return totprob;
	}

	public static boolean checkN(int[][] work, int n) {
		int tot = 0;

		for (int[] element : work) {
			for (int element2 : element) {
				tot += Math.abs(element2);
			}
		}

		return tot == n;
	}

	public static double workUp(int[][] work, double[] facts, double num, double pex,
															double[] lookup) {
		double denom, temp;

		for (int i = 0; i < work.length; i++) {
			for (int j = 0; j < work[0].length; j++) {
				facts[work[0].length * i + j + 1] = lookup[work[i][j]];
			}
		}
		denom = Array.sumInOrder(facts, Sort.quicksort(facts)); // 14 secs
		temp = Math.exp(num - denom);
		if (temp <= pex) {
			return temp;
		}

		return 0;
	}

	public static void main(String[] args) {
		int[][] matrix;
		double[] probs = {
											// 1.7685658194048082E-7,
											// 1.2346368281696975E-4,
											// 0.06782864732521461,
											0.005591506952405223, 0.027958060654154904, 0.0019136629015811507,
											0.008846883022512499, 0.009796226182780747, 0.007435141888923397,
											0.9887934356789729, 0.9887934356789694,};
		int[][][] tests = {{{3, 85}, {0, 183}}
				// {{0,1,14}, {3,4,38}, {0,1,33}, {0,1,47}, {3,1,18}},
				// {{10,15,236}, {2,45,362}, {13,37,485}, {2,56,459}, {7,133,1076}},
				// {{26,70,52}, {21,119,114}},
				// {{1, 2, 6}, {8, 7, 3}},
				// {{14, 19, 12}, {13, 27, 8}, {40, 22, 9}}, // 0.005591506952405223 time=7.37 sec, then 2.0
				// sec, now 0.72 sec and then 0.62
				// {{8, 5, 3}, {7, 8, 9}, {1, 2, 10}}, // 0.027958060654154904
				// {{25, 5, 5}, {7, 5, 12}, {15, 15, 9}}, // 0.0019136629015811507
				// {{15, 15, 5}, {7, 25, 2}, {15, 35, 19}}, // 0.008846883022512499
				// {{0, 2, 4}, {2, 4, 0}, {4, 0, 2}}, // 0.009796226182780747
				// {{0, 1, 5}, {2, 4, 0}, {4, 1, 1}}, // 0.007435141888923397
				// {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, // 0.9887934356789729
				// {{9, 8, 7}, {6, 5, 4}, {3, 2, 1}}, // 0.9887934356789694
				// {{13,132,78}, {9,50,89}, {12,92,150}},
				// {{14, 89, 120}, {13, 77, 58}, {40, 123, 91}},
				// {{1, 2, 10}, {8, 5, 3}, {7, 8, 9}},
		};

		PrintWriter writer;
		try {
			writer = new PrintWriter(new FileWriter("data4.xln"));
			matrix = tests[0];
			for (int i = 0; i < matrix.length; i++) {
				for (int j = 0; j < matrix[i].length; j++) {
					for (int k = 0; k < matrix[i][j]; k++) {
						writer.println(i + "\t" + j);
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + "data.xln");
			e.printStackTrace();
		}

		long time;
		double d;
		// System.out.println("Generic FishersExact");
		for (int i = 0; i < tests.length; i++) {
			time = new Date().getTime();
			matrix = tests[i];
			for (int j = 0; j < matrix.length; j++) {
				System.out.print((j == 0 ? "" : " / ") + Array.toStr(matrix[j], ","));
			}
			d = calc(matrix, 0, true);
			System.out.print("\tExact: "	+ d + "\t"
												+ (Math.abs(d - probs[i]) < 0.000000000000001	? "checks out"
																																			: "FAILED!!!"));
			System.out.print(" in " + ext.getTimeElapsed(time));
			System.out.println();
			System.out.println("Chi Sq: "
													+ ProbDist.ChiDist(	ContingencyTable.ChiSquare(matrix, false),
																							(matrix.length - 1) * (matrix[0].length - 1)));
		}
	}
}
