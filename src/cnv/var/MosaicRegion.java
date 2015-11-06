package cnv.var;

import java.io.Serializable;
import java.util.ArrayList;

import common.Array;

public class MosaicRegion extends CNVariant implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static final String[] ADD_HEADER = new String[] { "bpWeightedScore", "nearestStatScore", "pdfScore", "delta", "f", "customF", "numFMarkers", "numCustomFMarkers", "beastScore", "beastHeight", "beastLength" };
	private double nearestStateScore;
	private double bpWeightedScore;
	private double pdfScore;
	private double delta;
	private double f;
	private double customF;
	private int numFMarkers;
	private int numCustomFMarkers;
	private int beastLength;
	private double beastHeight;
	private double beastScore;

	public MosaicRegion(CNVariant cnv, double bpWeightedScore, double nearestStateScore, double pdfScore, double delta, double f, double customF) {
		super(cnv);
		this.bpWeightedScore = bpWeightedScore;
		this.nearestStateScore = nearestStateScore;
		this.pdfScore = pdfScore;
		this.delta = delta;
		this.f = f;
		this.customF = customF;
		this.numFMarkers = 0;
		this.numCustomFMarkers = 0;
		this.beastScore = 0;
	}

	public MosaicRegion(CNVariant cnv, MosaicRegion another) {
		super(cnv);
		this.bpWeightedScore = another.bpWeightedScore;
		this.nearestStateScore = another.nearestStateScore;
		this.pdfScore = another.pdfScore;
		this.delta = another.delta;
		this.f = another.f;
		this.customF = another.customF;
	}

	public void setF(double f) {
		this.f = f;
	}
	
	

	public void setBeastLength(int beastLength) {
		this.beastLength = beastLength;
	}

	public void setBeastHeight(double beastHeight) {
		this.beastHeight = beastHeight;
	}

	public void setBeastScore(double beastScore) {
		this.beastScore = beastScore;
	}

	public void setNumFMarkers(int numFMarkers) {
		this.numFMarkers = numFMarkers;
	}

	public void setNumCustomFMarkers(int numCustomFMarkers) {
		this.numCustomFMarkers = numCustomFMarkers;
	}

	@Override
	public String toAnalysisString() {
		ArrayList<String> tmp = new ArrayList<String>();
		tmp.add(bpWeightedScore + "");
		tmp.add(nearestStateScore + "");
		tmp.add(pdfScore + "");
		tmp.add(delta + "");
		tmp.add(f + "");
		tmp.add(customF + "");
		tmp.add(numFMarkers + "");
		tmp.add(numCustomFMarkers + "");
		tmp.add(beastScore + "");
		tmp.add(beastHeight + "");
		tmp.add(beastLength + "");

		String[] s = Array.concatAll(toPlinkFormat().split("\t"), Array.toStringArray(tmp));
		return Array.toStr(s);
	}


	public double getPdfScore() {
		return pdfScore;
	}

	public double getNearestStateScore() {
		return nearestStateScore;
	}

	public double getBpWeightedScore() {
		return bpWeightedScore;
	}

}
