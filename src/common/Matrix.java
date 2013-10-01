package common;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Vector;

public class Matrix {
	/**
     * Transposes an array of doubles
     * 
     * @param a
     *            original matrix of doubles
     * @return transposed array
     */
	public static double[][] transpose(double[][] a) {
		double m[][] = new double[a[0].length][a.length];

		for (int i = 0; i<a.length; i++) {
			for (int j = 0; j<a[i].length; j++) {
				m[j][i] = a[i][j];
			}
		}

		return m;
	}

	/**
     * Creates an array of double arrays of given size and initializes each
     * element with the given value
     * 
     * @param numRows
     *            number of rows for the of array
     * @param numColumns
     *            number of columns for the of array
     * @param initValue
     *            initial value of each element
     * @return array of arrays of numbers initialized to the given value
     */
    public static double[][] doubleMatrix(int numRows, int numColumns, double initValue) {
    	double[][] arr = new double[numRows][numColumns];
    	for (int i = 0; i<arr.length; i++) {
    		for (int j = 0; j<arr[i].length; j++) {
    			arr[i][j] = initValue;
    		}
    	}
    	return arr;
    }

    /**
     * Creates a matrix of integers of given size and initializes each element with the given value
     * 
     * @param numRows
     *            number of rows for the of array
     * @param numColumns
     *            number of columns for the of array
     * @param initValue
     *            initial value of each element
     * @return array of arrays of integers initialized to the given value
     */
    public static int[][] intMatrix(int numRows, int numColumns, int initValue) {
    	int[][] arrs = new int[numRows][numColumns];
    	for (int i = 0; i<arrs.length; i++) {
    		for (int j = 0; j<arrs[i].length; j++) {
    			arrs[i][j] = initValue;
    		}
    	}
    	return arrs;
    }

	/**
     * Creates a matrix of bytes of given size and initializes each element with the given value
     * 
     * @param numRows
     *            number of rows for the of array
     * @param numColumns
     *            number of columns for the of array
     * @param initValue
     *            initial value of each element
     * @return array of arrays of bytes initialized to the given value
     */
    public static byte[][] byteMatrix(int numRows, int numColumns, int initValue) {
    	byte[][] arrs = new byte[numRows][numColumns];
    	
    	if (initValue > Byte.MAX_VALUE || initValue < Byte.MIN_VALUE) {
    		System.err.println("Error - '"+initValue+"' is an invalid initValue for a byte matrix");
    		System.exit(1);
    	}
    	
    	for (int i = 0; i<arrs.length; i++) {
    		for (int j = 0; j<arrs[i].length; j++) {
    			arrs[i][j] = (byte)initValue;
    		}
    	}
    	return arrs;
    }

	/**
     * Creates a matrix of booleans of given size and initializes each element with the given value
     * 
     * @param numRows
     *            number of rows for the matrix
     * @param numColumns
     *            number of columns for the matrix
     * @param initValue
     *            initial value of each element
     * @return array of arrays of booleans initialized to the given value
     */
    public static boolean[][] booleanMatrix(int numRows, int numColumns, boolean initValue) {
    	boolean[][] arrs = new boolean[numRows][numColumns];
    	
    	for (int i = 0; i<arrs.length; i++) {
    		for (int j = 0; j<arrs[i].length; j++) {
    			arrs[i][j] = initValue;
    		}
    	}
    	return arrs;
    }

	/**
     * Creates a matrix of String of given size and initializes each element with the given value
     * 
     * @param numRows
     *            number of rows for the matrix
     * @param numColumns
     *            number of columns for the matrix
     * @param initValue
     *            initial value of each element
     * @return array of arrays of String initialized to the given value
     */
    public static String[][] stringMatrix(int numRows, int numColumns, String initValue) {
    	String[][] arrs = new String[numRows][numColumns];
    	
    	for (int i = 0; i<arrs.length; i++) {
    		for (int j = 0; j<arrs[i].length; j++) {
    			arrs[i][j] = initValue;
    		}
    	}
    	return arrs;
    }

	/**
     * Creates an array of integer arrays from the contents of a matrix of
     * numbers
     * 
     * @param matrix
     *            matrix of numbers to be converted
     * @return matrix of the converted numbers
     */
    public static int[][] toIntArraysRounded(double[][] matrix) {
    	int[][] mat = new int[matrix.length][];
    
    	for (int i = 0; i<matrix.length; i++) {
    		mat[i] = new int[matrix[i].length];
    		for (int j = 0; j<matrix[i].length; j++) {
    			mat[i][j] = (int)Math.round(matrix[i][j]);
    		}
    	}
    
    	return mat;
    }

	/**
     * Creates an array of and array of double arrays of given size and
     * initializes each element with the given value
     * 
     * @param length
     *            size of first dimension
     * @param breath
     *            size of second dimension
     * @param width
     *            size of third dimension
     * @param initValue
     *            initial value of each element
     * @return array of arrays of numbers initialized to the given value
     */
    public static double[][][] tripleArrays(int length, int breadth, int width, double initValue) {
    	double[][][] arr = new double[length][breadth][width];
    	for (int i = 0; i<arr.length; i++) {
    		for (int j = 0; j<arr[i].length; j++) {
    			for (int k = 0; k<arr[i][j].length; k++) {
    				arr[i][j][k] = initValue;
    			}
    		}
    	}
    	return arr;
    }

	/**
     * Creates an array of and array of double arrays of given size and
     * initializes each element with the given value
     * 
     * @param length
     *            size of first dimension
     * @param breath
     *            size of second dimension
     * @param width
     *            size of third dimension
     * @param initValue
     *            initial value of each element
     * @return array of arrays of numbers initialized to the given value
     */
    public static byte[][][] tripleArrays(int length, int breadth, int width, byte initValue) {
    	byte[][][] arrs = new byte[length][breadth][width];
    	for (int i = 0; i<arrs.length; i++) {
    		for (int j = 0; j<arrs[i].length; j++) {
    			for (int k = 0; k<arrs[i][j].length; k++) {
    				arrs[i][j][k] = initValue;
    			}
    		}
    	}
    	return arrs;
    }

    /**
     * Calculates the sum of a matrix
     * 
     * @param matrix
     *            a matrix of numbers
     * @return sum of the matrix
     */
    public static double sum(double[][] matrix) {
    	double sum = 0;
    
    	for (int i = 0; i<matrix.length; i++) {
    		sum += Array.sum(matrix[i]);
    	}
    
    	return sum;
    }

	/**
     * Calculates the sum of a matrix
     * 
     * @param matrix
     *            a matrix of integers
     * @return sum of the matrix
     */
    public static int sum(int[][] matrix) {
    	int sum = 0;
    
    	for (int i = 0; i<matrix.length; i++) {
    		sum += Array.sum(matrix[i]);
    	}
    
    	return sum;
    }
    
	/**
     * Creates an array of an array of Strings and copies the contents of a
     * vector into it
     * 
     * @param v
     *            vector of String arrays
     * @return an array of an array of Strings from the vector
     */
    public static String[][] toStringArrays(Vector<String[]> v) {
    	String[][] arrays = new String[v.size()][];
    	for (int i = 0; i<v.size(); i++) {
    		arrays[i] = v.elementAt(i);
    	}
    	return arrays;
    }

	/**
     * Creates a matrix of doubles and copies the contents of a vector of double
     * arrays into it
     * 
     * @param v
     *            a vector of double arrays
     * @return a matrix of doubles copied from a vector of double arrays
     */
    public static double[][] toDoubleArrays(Vector<double[]> v) {
    	double[][] matrix = new double[v.size()][];
    	for (int i = 0; i<v.size(); i++) {
    		matrix[i] = v.elementAt(i);
    	}
    	return matrix;
    }

	/**
     * Creates a matrix of doubles and copies the contents of a matrix of integers into it
     * 
     * @param v a matrix of integers
     * @return a matrix of doubles copied from a matrix of integers
     */
    public static double[][] toDoubleArrays(int[][] iMatrix) {
    	double[][] matrix = new double[iMatrix.length][];
    
    	for (int i = 0; i<iMatrix.length; i++) {
    		matrix[i] = Array.toDoubleArray(iMatrix[i]);
    	}
    	
    	return matrix;
    }

	/**
     * Creates a matrix of doubles and copies the contents of a matrix of String into it
     * 
     * @param v a matrix of String
     * @return a matrix of doubles copied from a matrix of String
     */
    public static double[][] toDoubleArrays(String[][] sMatrix) {
    	double[][] matrix = new double[sMatrix.length][];
    
    	for (int i = 0; i<sMatrix.length; i++) {
    		matrix[i] = Array.toDoubleArray(sMatrix[i]);
    	}
    	
    	return matrix;
    }

	/**
     * Creates a matrix of ints and copies the contents of a matrix of Strings into it
     * 
     * @param v a matrix of integers
     * @return a matrix of doubles copied from a matrix of integers
     */
    public static int[][] toIntArrays(String[][] sMatrix) {
    	int[][] matrix = new int[sMatrix.length][];
    
    	for (int i = 0; i<sMatrix.length; i++) {
    		matrix[i] = Array.toIntArray(sMatrix[i]);
    	}
    	
    	return matrix;
    }

	/**
     * Creates a matrix of float and copies the contents of a vector of float
     * arrays into it
     * 
     * @param v
     *            a vector of float arrays
     * @return a matrix of floats copied from a vector of float arrays
     */
    public static float[][] toFloatArrays(Vector<float[]> v) {
    	float[][] matrix = new float[v.size()][];
    	for (int i = 0; i<v.size(); i++) {
    		matrix[i] = v.elementAt(i);
    	}
    	return matrix;
    }

	/**
     * Creates a matrix of short and copies the contents of a vector of short
     * arrays into it
     * 
     * @param v
     *            a vector of short arrays
     * @return a matrix of shorts copied from a vector of short arrays
     */
    public static short[][] toShortArrays(Vector<short[]> v) {
    	short[][] matrix = new short[v.size()][];
    	for (int i = 0; i<v.size(); i++) {
    		matrix[i] = v.elementAt(i);
    	}
    	return matrix;
    }

	/**
     * Creates a matrix of int and copies the contents of a vector of int
     * arrays into it
     * 
     * @param v
     *            a vector of int arrays
     * @return a matrix of ints copied from a vector of int arrays
     */
	public static int[][] toMatrix(Vector<int[]> vs) {
    	int[][] matrix = new int[vs.size()][];
    
    	for (int i = 0; i<vs.size(); i++) {
    		matrix[i] = vs.elementAt(i);
    	}
    
    	return matrix;
    }
	
	/**
     * Creates a matrix of double and copies the contents of an array of double into the first column
     * 
     * @param ds
     *            an array of double
     * @return a matrix of double copied from an array of double
     */
	public static double[][] toMatrix(double[] ds) {
    	double[][] matrix = new double[ds.length][1];
    
    	for (int i = 0; i<ds.length; i++) {
    		matrix[i][0] = ds[i];
    	}
    
    	return matrix;
    }
	
	/**
     * Creates a matrix of String and copies the contents of an array of String into the first column
     * 
     * @param ds
     *            an array of double
     * @return a matrix of double copied from an array of double
     */
	public static String[][] toMatrix(String[] ds) {
    	String[][] matrix = new String[ds.length][1];
    
    	for (int i = 0; i<ds.length; i++) {
    		matrix[i][0] = ds[i];
    	}
    
    	return matrix;
    }
	
	/**
	 * Extracts a column from a matrix
	 * 
	 * @param data
	 *            two-dimensional matrix
	 * @param col
	 *            number of column to be extracted
	 * @return array of the elements in specified column
	 */
	public static double[] extractColumn(double[][] data, int col) {
		double[] array = new double[data.length];

		for (int i = 0; i<array.length; i++) {
			array[i] = data[i][col];
		}

		return array;
	}

	/**
	 * Extracts a column from a matrix
	 * 
	 * @param data
	 *            two-dimensional matrix
	 * @param col
	 *            number of column to be extracted
	 * @return array of the elements in specified column
	 */
	public static float[] extractColumn(float[][] data, int col) {
		float[] array = new float[data.length];

		for (int i = 0; i<array.length; i++) {
			array[i] = data[i][col];
		}

		return array;
	}

	/**
	 * Extracts a column from a matrix
	 * 
	 * @param data
	 *            two-dimensional matrix
	 * @param col
	 *            number of column to be extracted
	 * @return array of the elements in specified column
	 */
	public static int[] extractColumn(int[][] data, int col) {
		int[] array = new int[data.length];

		for (int i = 0; i<array.length; i++) {
			array[i] = data[i][col];
		}

		return array;
	}

	/**
	 * Extracts a column from a matrix
	 * 
	 * @param data
	 *            two-dimensional matrix
	 * @param col
	 *            number of column to be extracted
	 * @return array of the elements in specified column
	 */
	public static String[] extractColumn(String[][] data, int col) {
		String[] array = new String[data.length];

		for (int i = 0; i<array.length; i++) {
			if (data[i].length - 1 < col) {
				System.err.println("Error - trying to extract column index "+(col)+" from a row that doesn't have "+(col+1)+" columns (row index "+i+" only has "+(data[i].length-1)+" columns)");
			}
			array[i] = data[i][col];
		}

		return array;
	}

	/**
	 * Extracts certain columns from a matrix
	 * 
	 * @param data
	 *            two-dimensional matrix
	 * @param cols
	 *            the indices of the columns to be included in the resulting array in their final order
	 * @return array containing the elements of the matrix collapsed into a single String and in the specified order
	 */
	public static String[] extractColumns(String[][] data, int[] cols, String delimiter) {
		String[] array = new String[data.length];

		for (int i = 0; i<array.length; i++) {
			array[i] = "";
			for (int j = 0; j < cols.length; j++) {
				if (data[i].length - 1 < cols[j]) {
					System.err.println("Error - trying to extract column index "+(cols[j])+" from a row that doesn't have "+(cols[i]+1)+" columns (row index "+i+" only has "+(data[i].length-1)+" columns)");
				}
				array[i] += (j==0?"":delimiter)+data[i][cols[j]];
			}
		}

		return array;
	}
	
	/**
	 * Creates a new matrix from an old matrix using the specified columns in the specified order
	 * 
	 * @param data
	 *            two-dimensional matrix
	 * @param cols
	 *            the indices of the columns to be included in the final matrix in their final order
	 * @return the rearranged matrix in specified order
	 */
	public static String[][] extractColumns(String[][] data, int[] cols) {
		String[][] matrix = new String[data.length][cols.length];

		for (int i = 0; i<matrix.length; i++) {
			for (int j = 0; j < cols.length; j++) {
				if (data[i].length - 1 < cols[j]) {
					System.err.println("Error - trying to extract column index "+(cols[j])+" from a row that doesn't have "+(cols[i]+1)+" columns (row index "+i+" only has "+(data[i].length-1)+" columns)");
				}
				matrix[i][j] = data[i][cols[j]];
			}
		}

		return matrix;
	}
	
	/**
	 * Creates a new matrix from an old matrix using the specified columns in the specified order
	 * 
	 * @param data
	 *            two-dimensional matrix
	 * @param cols
	 *            the indices of the columns to be included in the final matrix in their final order
	 * @return the rearranged matrix in specified order
	 */
	public static double[][] extractColumns(double[][] data, int[] cols) {
		double[][] matrix = new double[data.length][cols.length];

		for (int i = 0; i<matrix.length; i++) {
			for (int j = 0; j < cols.length; j++) {
				if (data[i].length - 1 < cols[j]) {
					System.err.println("Error - trying to extract column index "+(cols[j])+" from a row that doesn't have "+(cols[i]+1)+" columns (row index "+i+" only has "+(data[i].length-1)+" columns)");
				}
				matrix[i][j] = data[i][cols[j]];
			}
		}

		return matrix;
	}
	
	public static int getSize(String[][] matrix) {
    	int count = 0;
    
    	for (int i = 0; i<matrix.length; i++) {
    		count += matrix[i].length;
    	}
    
    	return count;
    }

	/**
     * Returns an array formed by lining up the rows of a matrix
     * 
     * @param matrix
     *            a matrix of numbers
     * @return array of all the numbers
     */
    public static double[] collapseMatrix(double[][] matrix) {
    	double[] array;
    	int size = 0, count;
    
    	for (int i = 0; i<matrix.length; i++) {
    		size += matrix[i].length;
    	}
    
    	array = new double[size];
    
    	count = 0;
    	for (int i = 0; i<matrix.length; i++) {
    		for (int j = 0; j<matrix[i].length; j++) {
    			array[count++] += matrix[i][j];
    		}
    	}
    
    	return array;
    }

	/**
     * Returns an array formed by lining up the rows of a matrix
     * 
     * @param matrix
     *            a matrix of integers
     * @return array of all the integers
     */
    public static int[] collapseMatrix(int[][] matrix) {
    	int[] array;
    	int size = 0, count;
    
    	for (int i = 0; i<matrix.length; i++) {
    		size += matrix[i].length;
    	}
    
    	array = new int[size];
    
    	count = 0;
    	for (int i = 0; i<matrix.length; i++) {
    		for (int j = 0; j<matrix[i].length; j++) {
    			array[count++] += matrix[i][j];
    		}
    	}
    
    	return array;
    }
    
    public static int indexInArray(int i, int j, int maxJ) {
		return i*maxJ+j;
	}

    public static int[] indicesInMatrix(int indexInArray, int maxJ) {
		return new int[] {indexInArray/maxJ,  indexInArray%maxJ};
	}
    
    public static String toStr(String[][] matrix) {
    	String temp;
    	
    	temp = "";
    	for (int i = 0; i<matrix.length; i++) {
    		for (int j = 0; j<matrix[i].length; j++) {
    			temp += (j==0?"":"\t")+matrix[i][j];
            }
    		temp += "\n";
        }
		
    	return temp;
	}
    
    public static void writeToFile(String[][] matrix, String filename) {
        PrintWriter writer;

        try {
	        writer = new PrintWriter(new FileWriter(filename));
	    	for (int i = 0; i<matrix.length; i++) {
	    		for (int j = 0; j<matrix[i].length; j++) {
	    	        writer.print((j==0?"":"\t")+matrix[i][j]);
	            }
	    		writer.println();
	        }

	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+filename);
	        e.printStackTrace();
        }
		
	}
    
	/**
	 * Removes NaN values from a set of arrays
	 * 
	 * @param matrix
	 *            a matrix of doubles
	 * @return scrubbed matrix
	 */
	public static double[][] removeRowsWithNaN(double[][] matrix) {
		double[][] newMatrix;
		boolean[] use;
		int count;

		for (int i = 1; i<matrix.length; i++) {
			if (matrix[i].length != matrix[0].length) {
				System.err.println("Error - matrix does not have a consistent number of elements");
				return null;
			}
		}
		
		use = new boolean[matrix[0].length];
		for (int i = 0; i<use.length; i++) {
			use[i] = true;
			for (int j = 0; j<matrix.length; j++) {
				if ((matrix[j][i]+"").equals("NaN")) {
					use[i] = false;					
				}
            }
        }

		count = Array.booleanArraySum(use);
		if (count < use.length) {
			newMatrix = new double[matrix.length][count];
			count = 0;
			for (int i = 0; i<use.length; i++) {
				if (use[i]) {
					for (int j = 0; j<matrix.length; j++) {
						newMatrix[j][count] = matrix[j][i]; 
	                }
					count++;
				}
	        }
			
			return newMatrix;
		} else {
			return matrix;
		}
	}
	
	/**
	 * Only keeps specified columns from a matrix
	 * 
	 * @param matrix
	 *            a matrix of String
	 * @return pruned matrix
	 */
	public static String[][] prune(String[][] matrix, int[] rowsToKeep, int[] columnsToKeep, Logger log) {
		if (matrix == null || matrix.length == 0) {
			log.reportError("Error - null/empty matrix; can't be pruned");
			return null;
		}

		return prune(matrix, rowsToKeep==null?null:Array.indicesToBooleanArray(rowsToKeep, matrix.length), columnsToKeep==null?null:Array.indicesToBooleanArray(columnsToKeep, matrix[0].length), log);
	}
	
	/**
	 * Only keeps specified columns from a matrix
	 * 
	 * @param matrix
	 *            a matrix of String
	 * @return pruned matrix
	 */
	public static String[][] prune(String[][] matrix, boolean[] rowsToKeep, boolean[] columnsToKeep, Logger log) {
		String[][] newMatrix;
		int row, col;
		
		if (matrix == null) {
			log.reportError("Error - null matrix; can't be pruned");
			return null;
		}

		if (matrix.length == 0) {
			log.reportError("Error - empty matrix; can't be pruned");
			return new String[0][0];
		}
		
		row = 0;
		newMatrix = new String[rowsToKeep == null?matrix.length:Array.booleanArraySum(rowsToKeep)][columnsToKeep == null?matrix[0].length:Array.booleanArraySum(columnsToKeep)];
		for (int i = 0; i < matrix.length; i++) {
			if (rowsToKeep == null || rowsToKeep[i]) {
				if (columnsToKeep == null) {
					newMatrix[row] = matrix[i];
				} else {
					col = 0;
					for (int j = 0; j < matrix[i].length; j++) {
						if (columnsToKeep[j]) {
							newMatrix[row][col] = matrix[i][j];
							col++;
						}
					}
				}
				row++;
			}
		}
		
		return newMatrix;
	}
	
	/**
	 * Removes rows/columns with maginal counts of zero
	 * 
	 * @param matrix
	 *            a matrix of ints
	 * @return pruned matrix
	 */
	public static int[][] prune(int[][] matrix) {
		boolean[] use;
		int[][] temp, newMatrix;
		int[] r, c;
		int count;
		boolean pruned;
		
		r = new int[matrix.length];
		c = new int[matrix[0].length];
		pruned = false;
		
		for (int i = 0; i<r.length; i++) {
			for (int j = 0; j<c.length; j++) {
				r[i] += matrix[i][j];
				c[j] += matrix[i][j];
            }
        }
		
		newMatrix = matrix;
		if (Array.min(c) == 0) {
			use = new boolean[c.length];
			for (int j = 0; j<c.length; j++) {
				use[j] = c[j]!=0;
            }
			temp = new int[r.length][Array.booleanArraySum(use)];
			count = 0;
			for (int j = 0; j<c.length; j++) {
				if (use[j]) {
					for (int i = 0; i<r.length; i++) {
						temp[i][count] = newMatrix[i][j];
		            }
					count++;
				}
            }
			newMatrix = null;
			newMatrix = temp;
			pruned = true;
		}
		
		if (Array.min(r) == 0) {
			use = new boolean[r.length];
			for (int i = 0; i<r.length; i++) {
				use[i] = r[i]!=0;
            }
			temp = new int[Array.booleanArraySum(use)][newMatrix[0].length];
			count = 0;
			for (int i = 0; i<r.length; i++) {
				if (use[i]) {
					for (int j = 0; j<newMatrix[0].length; j++) {
						temp[count][j] = newMatrix[i][j];
		            }
					count++;
				}
            }
			newMatrix = null;
			newMatrix = temp;
			pruned = true;
		}

		if (pruned) {
			return newMatrix;
		} else {
			return matrix;
		}
	}
	
	public static boolean equals(int[][] m1, int[][] m2) {
		if (m1.length != m2.length) {
			return false;
		}
		for (int i = 0; i < m1.length; i++) {
			if (m1[i].length != m2[i].length) {
				return false;
			}
			for (int j = 0; j < m1[i].length; j++) {
				if (m1[i][j] != m2[i][j]) {
					return false;
				}
			}
		}
		return true;		
	}
	
	public static int[][] clone(int[][] matrix) {
		int[][] clone;
		
		clone = new int[matrix.length][];
		for (int i = 0; i < matrix.length; i++) {
			clone[i] = new int[matrix[i].length];
			for (int j = 0; j < matrix[i].length; j++) {
				clone[i][j] = matrix[i][j];
			}
		}
		
		return clone;		
	}

	public static String[][] clone(String[][] matrix) {
		String[][] clone;
		
		clone = new String[matrix.length][];
		for (int i = 0; i < matrix.length; i++) {
			clone[i] = new String[matrix[i].length];
			for (int j = 0; j < matrix[i].length; j++) {
				clone[i][j] = matrix[i][j];
			}
		}

		return clone;
	}

	public static int[][] putInOrder(int[][] matrix, int[] order) {
		int[][] newMatrix;
		
		newMatrix = new int[matrix.length][];
		for (int i = 0; i < matrix.length; i++) {
			newMatrix[i] = new int[matrix[i].length];
			for (int j = 0; j < matrix[i].length; j++) {
				newMatrix[i][j] = matrix[order[i]][j];
			}
		}
		
		return newMatrix;		
	}
	
	/**
	 * Return the minimum in a matrix of integers
	 * 
	 * @param matrix
	 *            matrix of integers
	 * @return the minimum
	 */
	public static int min(int[][] matrix) {
		int min;

		min = matrix[0][0];
		for (int i = 0; i<matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				min = Math.min(matrix[i][j], min);
			}
		}
		return min;
	}
	
	// TODO needs documentation
	public static double[][] computeRowProportions(int[][] counts) {
		double[][] rowProportions;
		
		rowProportions = new double[counts.length][];
		for (int i = 0; i<counts.length; i++) {
			rowProportions[i] = Array.getProportions(counts[i]);
		}
		
		return rowProportions;
	}
	
	public static String[][] addRow(String[][] matrix, String[] row) {
		String[][] newMatrix;
		
		newMatrix = new String[matrix.length+1][];
		for (int i = 0; i<matrix.length; i++) {
			newMatrix[i] = matrix[i];
		}
		newMatrix[matrix.length] = row;
		
		return newMatrix;
	}
	
	public static void main(String[] args) {
		int[][] matrix = {{1, 2, 4}, {1, 4, 0}, {0, 1, 0}};
		
		for (int i = 0; i<matrix.length; i++) {
			System.out.println(Array.toStr(matrix[i]));
        }
		System.out.println();
		matrix = prune(matrix);
		System.out.println();
		for (int i = 0; i<matrix.length; i++) {
			System.out.println(Array.toStr(matrix[i]));
        }
		System.out.println();
		
    }

	public static double[][] subset(double[][] matrix, boolean[] rowsToKeep) {
		double[][] subset;
		int count;
		
		count = 0;
		subset = new double[Array.booleanArraySum(rowsToKeep)][];
		for (int i = 0; i < matrix.length; i++) {
			if (rowsToKeep[i]) {
				subset[count++] = matrix[i];
			}
		}
		return null;
	}
}
