package cnv.filesys;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;

import stats.Correlation;
import common.Array;

public class MarkerData implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String[][] TYPES = { {"X Raw", "Y Raw"}, {"X", "Y"}, {"Theta", "R"}, {"B Allele Freq", "Log R Ratio"}};

	private String markerName;
	private byte chr;
	private int position; 
	private long fingerprint;
	private float[] gcs;
	private float[] xRaws;
	private float[] yRaws;
	private float[] xs;
	private float[] ys;
	private float[] thetas;
	private float[] rs;
	private float[] lrrs;
	private float[] bafs;
	private byte[] abGenotypes;
	private String[] alleleMappings;

	public MarkerData(String markerName, byte chr, int position, long fingerprint, float[] gcs, float[] xRaws, float[] yRaws, float[] xs, float[] ys, float[] thetas, float[] rs, float[] bafs, float[] lrrs, byte[] abGenotypes, String[] alleleMappings) {
		this.markerName = markerName;
		this.chr = chr;
		this.position = position;
		this.fingerprint = fingerprint;
		this.gcs = gcs;
		this.xRaws = xRaws;
		this.yRaws = yRaws;
		this.xs = xs;
		this.ys = ys;
		this.thetas = thetas;
		this.rs = rs;
		this.bafs = bafs;
		this.lrrs = lrrs;
		this.abGenotypes = abGenotypes;
		this.alleleMappings = alleleMappings;
	}

	public float[][] getDatapoints(int type) {
		switch (type) {
		case 0:
			return new float[][] {xRaws, yRaws};
		case 1:
			return new float[][] {xs, ys};
		case 2:
			return new float[][] {thetas, rs};
		case 3:
			return new float[][] {bafs, lrrs};
		default:
			System.err.println("Error - invalid plot type");
			return null;
		}
	}

	public String getMarkerName() {
		return markerName;
	}

	public byte getChr() {
		return chr;
	}

	public int getPosition() {
		return position;
	}
	
	public long getFingerprint() {
		return fingerprint;
	}

	public float[] getGCs() {
		return gcs;
	}

	public float[] getX_Raws() {
		return xRaws;
	}

	public float[] getY_Raws() {
		return yRaws;
	}

	public float[] getXs() {
		return xs;
	}

	public float[] getYs() {
		return ys;
	}

	public float[] getThetas() {
		if (thetas == null) {
			thetas = new float[xs.length];
			for (int i = 0; i<xs.length; i++) {
				thetas[i] = Centroids.calcTheta(xs[i], ys[i]);
            }
		}
		return thetas;
	}

	public float[] getRs() {
		if (rs == null) {
			rs = new float[xs.length];
			for (int i = 0; i<xs.length; i++) {
				rs[i] = Centroids.calcR(xs[i], ys[i]);
            }
		}
		return rs;
	}

	public void recompute(float[][] centroids) {
		float[] thetas, rs;

		thetas = getThetas();
		rs = getRs();
		
		bafs = new float[xs.length];
		lrrs = new float[xs.length];
		for (int i = 0; i<xs.length; i++) {
			bafs[i] = Centroids.calcBAF(thetas[i], centroids);
			lrrs[i] = Centroids.calcLRR(thetas[i], rs[i], centroids);
        }
	}

	public float[] getBAFs() {
		return bafs;
	}

	public float[] getLRRs() {
		return lrrs;
	}

	public byte[] getAB_Genotypes() {
		return abGenotypes;
	}

	public byte[] getAB_GenotypesAfterFilters(ArrayList<MarkerFilter> filters, float gcThreshold) {
		byte[] newGenotypes;
		
		// apply filters from abGenotypes -> newGenotypes
		newGenotypes = abGenotypes;
		
		
		return newGenotypes;
	}

	public String[] getAlleleMappings() {
		return alleleMappings;
	}
	
	public double[] compareLRRs(float[][] centroids) {
		double[] originalLRRs, compLRRs;
		double error;
		int count;
		
		error = 0;
		count = 0;
		originalLRRs = new double[lrrs.length];
		compLRRs = new double[lrrs.length];
		for (int i = 0; i<lrrs.length; i++) {
			if (!Float.isNaN(lrrs[i])) {
				originalLRRs[count] = lrrs[i];
				compLRRs[count] = Centroids.calcLRR(Centroids.calcTheta(xs[i], ys[i]), Centroids.calcR(xs[i], ys[i]), centroids);
				if (Double.isNaN(compLRRs[count])) {
					System.err.println("Error - compLRR is invalid ("+compLRRs[count]+") where oriLRR is not ("+lrrs[count]+")");
				} else {
					error += Math.abs(compLRRs[count]-originalLRRs[count]);
					count++;
					if (Double.isNaN(error) || Double.isInfinite(error)) {
						System.err.println("Started with index "+i+", compLRR of '"+compLRRs[count]+"', and oriLRR of '"+originalLRRs[count]+"'");
						return new double[] {-999,-999};
					}
				}
			}
        }
		
//		System.out.println("error="+error+" count="+count);
		return new double[] {Correlation.Pearson(Array.subArray(originalLRRs, 0, count), Array.subArray(compLRRs, 0, count))[0], error/count};
	}
	
	public void writeToFile(String[] samples, String filename) {
		PrintWriter writer;

		if (samples.length!=lrrs.length) {
			System.err.println("Error - Number of samples (n="+samples.length+") does not match up with the number of LRRs/BAFs/etc (n="+lrrs.length+")");
			System.exit(1);
		}
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("Sample\tGC Score\tRaw X\tRaw Y\tX\tY\tTheta\tR\tLRR\tBAF\tGenotypes");
			for (int i = 0; i<samples.length; i++) {
				writer.println(samples[i]+"\t"+gcs[i]+"\t"+xRaws[i]+"\t"+yRaws[i]+"\t"+xs[i]+"\t"+ys[i]+"\t"+thetas[i]+"\t"+rs[i]+"\t"+lrrs[i]+"\t"+bafs[i]+"\t"+alleleMappings[abGenotypes[i]+1]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+filename);
			e.printStackTrace();
		}
	}
}
