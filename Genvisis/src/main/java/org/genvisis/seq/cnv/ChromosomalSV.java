package org.genvisis.seq.cnv;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.filesys.MarkerSet.PreparedMarkerSet;
import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.seq.manage.BamImport.NGS_MARKER_TYPE;

/**
 * Looking at large ( all or half) of chromosome structural variants
 *
 */
public class ChromosomalSV {

	public static void run(Project proj, int numthreads) {
		Logger log = proj.getLog();
		PreparedMarkerSet preparedMarkerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());
		String[] samples = proj.getSamples();
		ChrProducer producer = new ChrProducer(proj, samples, preparedMarkerSet);
		WorkerTrain<ChrResult[][]> train = new WorkerTrain<ChrResult[][]>(producer, numthreads, 10, proj.getLog());
		int index = 0;
		ChrResult[][][] allResults = new ChrResult[samples.length][][];
		Hashtable<String, ArrayList<Double>> summaryMedian = new Hashtable<String, ArrayList<Double>>();

		while (train.hasNext()) {
			allResults[index] = train.next();
//			for (int i = 0; i < allResults[index].length; i++) {
//				if (!summaryMedian.containsKey(i + "")) {
//					summaryMedian.put(i + "", new ArrayList<Double>());
//				}
//				if (!Double.isNaN(allResults[index][i].getMedian())) {
//					summaryMedian.get(i + "").add(allResults[index][i].getMedian());
//				}
//			}
			index++;
			proj.getLog().reportTimeInfo(index + "");
		}

		double[] allMedians = new double[preparedMarkerSet.getChrs().length];
		for (int i = 0; i < allMedians.length; i++) {
			allMedians[i] = Double.NaN;
			if (summaryMedian.containsKey(i + "") && summaryMedian.get(i + "").size() > 2) {
				ArrayList<Double> tmp = summaryMedian.get(i + "");
				double[] d = Array.removeNaN(Array.toDoubleArray(tmp));
				allMedians[i] = Array.median(d);
			}
		}

		String outDir = proj.PROJECT_DIRECTORY.getValue() + "chrSV/";
		new File(outDir).mkdirs();
		String outFile = outDir + "chr.svs.txt";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			writer.println(Array.toStr(new String[] { "Sample", "Chr", "Median", "TYPE" }));

			for (int i = 0; i < allResults.length; i++) {
				for (int j = 0; j < allResults[i].length; j++) {
					for (int j2 = 0; j2 < allResults[i][j].length; j2++) {
						double median = allResults[i][j][j2].getMedian();
						if (!Double.isNaN(median)) {
							writer.println(samples[i] + "\t" + j + "\t" + median + "\t" + allResults[i][j][j2].getType());
						}
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outFile);
			log.reportException(e);
		}
	}

	enum TYPE {
		LEFT, RIGHT, ALL;
	}

	private static class ChrResult {

		//private byte chr;
		private double median;
		//private double mad;
		//private double stDev;
		private TYPE type;

		private ChrResult(byte chr, double median, double mad, double stDev, TYPE type) {
			super();
			//this.chr = chr;
			this.median = median;
			//this.mad = mad;
			//this.stDev = stDev;
			this.type = type;
		}

		public double getMedian() {
			return median;
		}

		public TYPE getType() {
			return type;
		}

	}

	private static class ChrProducer implements Producer<ChrResult[][]> {
		private Project proj;
		private String[] samples;
		private PreparedMarkerSet preparedMarkerSet;
		private int index;

		public ChrProducer(Project proj, String[] samples, PreparedMarkerSet preparedMarkerSet) {
			super();
			this.proj = proj;
			this.samples = samples;
			this.preparedMarkerSet = preparedMarkerSet;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<ChrResult[][]> next() {
			String samp = samples[index];
			index++;
			return new ChrWorker(proj, samp, preparedMarkerSet);
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static class ChrWorker implements Callable<ChrResult[][]> {
		private Project proj;
		private String sample;
		private PreparedMarkerSet preparedMarkerSet;

		public ChrWorker(Project proj, String sample, PreparedMarkerSet preparedMarkerSet) {
			super();
			this.proj = proj;
			this.sample = sample;
			this.preparedMarkerSet = preparedMarkerSet;
		}

		@Override
		public ChrResult[][] call() throws Exception {
			return tallyChrs(proj, sample, preparedMarkerSet);
		}

	}

	private static ChrResult[][] tallyChrs(Project proj, String sample, PreparedMarkerSet preparedMarkerSet) {
		Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
		int[] pos = preparedMarkerSet.getPositions();
		int[][] boundaries = Positions.determineCentromereBoundariesFromMarkerSet(preparedMarkerSet.getChrs(), preparedMarkerSet.getPositions(), 37, proj.getLog());
		double[] lrrs = Array.toDoubleArray(samp.getLRRs());
		ChrResult[][] results = new ChrResult[preparedMarkerSet.getIndicesByChr().length][3];
		String[] markerNames = preparedMarkerSet.getMarkerNames();
		for (int i = 0; i < preparedMarkerSet.getIndicesByChr().length; i++) {
			int[] indices = preparedMarkerSet.getIndicesByChr()[i];
			ArrayList<Integer> noOffTarget = new ArrayList<Integer>();
			for (int j = 0; j < indices.length; j++) {
				if (markerNames[indices[j]].contains(NGS_MARKER_TYPE.OFF_TARGET.getFlag())) {
					noOffTarget.add(indices[j]);
				}
			}
			indices=Array.toIntArray(noOffTarget);
			ChrResult chrResultAll = getResult(lrrs, i, indices, TYPE.ALL);
			results[i][0] = chrResultAll;
			int[] bounds = boundaries[i];
			ArrayList<Integer> left = new ArrayList<Integer>();
			ArrayList<Integer> right = new ArrayList<Integer>();

			for (int j = 0; j < indices.length; j++) {
				if (pos[indices[j]] <= bounds[0]) {
					left.add(indices[j]);
				} else if (pos[indices[j]] >= bounds[1]) {
					right.add(indices[j]);
				}
			}
			ChrResult chrResultLeft = getResult(lrrs, i, Array.toIntArray(left), TYPE.LEFT);
			ChrResult chrResultRight = getResult(lrrs, i, Array.toIntArray(right), TYPE.RIGHT);
			results[i][1] = chrResultLeft;
			results[i][2] = chrResultRight;

		}
		return results;
	}

	private static ChrResult getResult(double[] lrrs, int i, int[] indices, TYPE type) {
		byte chr = (byte) i;
		double median = Double.NaN;
		double mad = Double.NaN;
		double stDev = Double.NaN;
		if (indices.length > 0) {
			double[] subLrr = Array.removeNaN(Array.subArray(lrrs, indices));
			if (subLrr.length > 2) {
				median = Array.median(subLrr);
				mad = Array.mad(subLrr);
				stDev = Array.stdev(subLrr);
			}
		}
		ChrResult chrResult = new ChrResult(chr, median, mad, stDev, type);
		return chrResult;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		int numthreads = 7;

		String usage = "\n" +
				"seq.cnv.ChromosomalSV requires 0-1 arguments\n" +
				"   (1) proj (i.e. proj=" + filename + " (default))\n" +
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Project proj = new Project(filename, false);
		run(proj, numthreads);

	}

}
