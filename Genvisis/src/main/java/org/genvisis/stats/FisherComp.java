// -Xms1024M -Xmx1024M
package org.genvisis.stats;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Matrix;

public class FisherComp {
	public static void comp(String filename, int permutations) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Vector<String> v;
		int[][] matrix, prunedMatrix;
		long time;

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename + "_comp.xln"));
			writer.println(reader.readLine() + "\tChiSq\tFishersExact\ttime\tcollapsed");
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t", -1);
				v = new Vector<String>();
				for (int i = 1; i <= 3; i++) {
					if (!line[i].equals("")) {
						v.add(line[i]);
					}
				}
				matrix = new int[v.size()][];
				for (int i = 0; i < v.size(); i++) {
					matrix[i] = Array.toIntArray(v.elementAt(i).split("/"));
				}
				prunedMatrix = Matrix.prune(matrix);
				time = new Date().getTime();
				writer.println(Array.toStr(line)	+ "\t"
												+ ProbDist.ChiDist(	ContingencyTable.ChiSquare(prunedMatrix, false),
																						(prunedMatrix.length - 1) * (prunedMatrix[0].length
																																					- 1))
												+ "\t" + FishersExact.calc(prunedMatrix, 0, false) + "\t"
												+ (new Date().getTime() - time) + "\t"
												+ (matrix.length == prunedMatrix.length
														&& matrix[0].length == prunedMatrix[0].length ? 0 : 1));
				writer.flush();
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		// String filename = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\Consortium\\00src\\Miami\\fishFood1.dat";
		String filename =
										"C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\00src\\Miami\\homogeneityTests_NaNs.xln";
		int perms = 10000;

		String usage = "\n"	+ "stats.FisherComp requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n"
										+ "   (2) number of permutations (i.e. perms=" + perms + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("perms=")) {
				perms = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			comp(filename, perms);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
