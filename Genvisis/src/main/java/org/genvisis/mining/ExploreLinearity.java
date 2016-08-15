package org.genvisis.mining;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;
import org.genvisis.stats.Correlation;
import org.genvisis.stats.LeastSquares;

import com.google.common.primitives.Doubles;

public class ExploreLinearity {
	private double[][][] correlations;

	private double[][][] R_squares;

	private String depLabel = null;

	private String[] indepLabels = null;

	private int preferredDepTrans = -1;

	public ExploreLinearity(double[] deps, double[][] indeps) {
		int M = indeps[0].length;
		double[][][] transIndeps = new double[M][Transformations.NUM_TRANSFORMATIONS][deps.length];
		double[][] transDeps = new double[Transformations.NUM_TRANSFORMATIONS][deps.length];

		for (int i = 0; i<M; i++) {
			for (int j = 0; j<Transformations.NUM_TRANSFORMATIONS; j++) {
				transIndeps[i][j] = Transformations.transform(Matrix.extractColumn(indeps, i), j);
			}
		}

		correlations = new double[M][Transformations.NUM_TRANSFORMATIONS][Transformations.NUM_TRANSFORMATIONS];
		R_squares = new double[M][Transformations.NUM_TRANSFORMATIONS][Transformations.NUM_TRANSFORMATIONS];
		for (int k = 0; k<Transformations.NUM_TRANSFORMATIONS; k++) {
			transDeps[k] = Transformations.transform(deps, k);
			for (int i = 0; i<M; i++) {
				for (int j = 0; j<Transformations.NUM_TRANSFORMATIONS; j++) {
					correlations[i][j][k] = Correlation.Pearson(transDeps[k], transIndeps[i][j])[0];
					R_squares[i][j][k] = new LeastSquares(transDeps[k], LeastSquares.Transpose(new double[][] {transIndeps[i][j]})).getRsquare();
				}
			}
		}

	}

	public void setLabels(String dep, String[] indeps) {
		this.depLabel = dep;
		this.indepLabels = indeps;
	}

	public void setPreferredDepTrans(int value) {
		this.preferredDepTrans = value;
	}

	public void report(String filename, int measure) {
		PrintWriter writer = null;

		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (int i = 0; i<correlations.length; i++) {
				writer.println((depLabel==null?"Dep":depLabel)+" versus "+(indepLabels==null?"Indep "+(i+1):indepLabels[i]));
				writer.println();
				if (preferredDepTrans==-1) {
					writer.println(Transformations.getLabel(-1)+"\tTransformation of "+(depLabel==null?"Dep":depLabel));
					writer.print(Transformations.getLabel(-1));
					for (int j = 0; j<Transformations.NUM_TRANSFORMATIONS; j++) {
						writer.print(ext.formStr("     "+(j+1), 8, true));
					}
					writer.println();
				}
				for (int j = 0; j<Transformations.NUM_TRANSFORMATIONS; j++) {
					writer.print((preferredDepTrans==-1?(j+1)+" ":"")+Transformations.getLabel(j));
					for (int k = (preferredDepTrans==-1?0:preferredDepTrans); k<(preferredDepTrans==-1?Transformations.NUM_TRANSFORMATIONS:preferredDepTrans+1); k++) {
						writer.print(ext.formStr(measure==1?(correlations[i][j][k]<0?"":" ")+ext.formDeci(correlations[i][j][k], 3, true):ext.formDeci(R_squares[i][j][k], 3, true), 8, true));
					}
					writer.println();
				}
				writer.println();
				writer.println();
			}

			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error writing file \""+filename+"\"");
			System.exit(2);
		}
	}

	public static void exploreFromFile(String filename, boolean lastNotFirst) {
		BufferedReader reader = null;
		String[] line;
		Vector<double[]> v = new Vector<double[]>();
		double[][] data;
		double[] dataline;
		DoubleVector dv = new DoubleVector();
		ExploreLinearity el;
		String depLabel;
		String[] indepLabels;

		try {
			reader = new BufferedReader(new FileReader(filename));
			line = reader.readLine().split("[\\s]+");
			depLabel = line[lastNotFirst?line.length-1:0];
			indepLabels = new String[line.length-1];
			for (int i = 0; i<indepLabels.length; i++) {
				indepLabels[i] = line[lastNotFirst?i:i+1];
			}
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				dv.add(Double.parseDouble(line[lastNotFirst?line.length-1:0]));
				dataline = new double[line.length-1];
				for (int i = 0; i<dataline.length; i++) {
					dataline[i] = Double.parseDouble(line[lastNotFirst?i:i+1]);
				}
				v.add(dataline);
			}
			reader.close();

			dataline = Doubles.toArray(dv);
			data = new double[v.size()][];
			for (int i = 0; i<v.size(); i++) {
				data[i] = v.elementAt(i);
			}

			el = new ExploreLinearity(dataline, data);
			el.setLabels(depLabel, indepLabels);
			el.setPreferredDepTrans(2);
			el.report(filename+"-linearity-pearson.out", 1);
			el.report(filename+"-linearity-rsquare.out", 2);

		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "Boston.xls";

		String usage = "\n"+"park.ExploreLinearity requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default)\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			exploreFromFile(filename, true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
